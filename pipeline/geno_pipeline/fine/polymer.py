"""L'étage 3D : un polymère de chromatine sous Langevin, avec des liaisons qui bougent.

Champ de force, cinq termes, aucun ajustable après coup :

1. **squelette** — liaison harmonique entre monomères consécutifs, longueur au
   repos `sigma` (un diamètre de monomère) ;
2. **rigidité** — terme angulaire `k·(1 − cos θ)`, qui donne à la fibre une
   longueur de persistance de l'ordre de quelques monomères ;
3. **volume exclu mou** — répulsion tronquée à `trunc` kT au contact plein. Elle
   est **volontairement franchissable** : à l'échelle de temps de l'extrusion, la
   topoisomérase II laisse passer les brins, et un volume exclu dur figerait des
   enlacements que la cellule défait. C'est le choix de polychrom, et il se paie
   — la topologie de la région n'est pas une prédiction du modèle ;
4. **confinement** — mur sphérique mou au rayon que le modèle de noyau alloue à
   cette quantité de chromatine (§ `region.blob_radius`) ;
5. **cohésines** — une liaison par emplacement, déplacée d'un cran à chaque mise
   à jour.

**Pourquoi une sphère et pas une boîte périodique.** Une boîte périodique n'a pas
de paroi, ce qui est séduisant, mais elle remplace l'environnement du segment par
ses propres images : les statistiques à grande échelle saturent alors à la taille
de la boîte, qui vaut ici 13,8 diamètres pour 8,6 de rayon de blob — l'artefact
tombe en plein dans la fenêtre où le raccord se mesure. La sphère a une paroi, ce
qui est un défaut qu'on énonce, mais elle fixe *exactement* le volume que la
semaine 6 accorde à ces mégabases, et c'est cette quantité-là que le raccord
compare.

**Unités.** Longueur = un diamètre de monomère `sigma`, énergie = kT à 300 K.
OpenMM lit des nanomètres et des kJ/mol ; on lui donne `sigma = 1 nm` et on
convertit à la sortie. Le facteur de conversion est `sigma_nm`, fixé par la loi
de la semaine 6 et pas par cette semaine-ci.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

KT = 2.494339          # kJ/mol à 300 K
MASS = 100.0           # uma, convention polychrom
TIMESTEP_FS = 70.0
FRICTION_PER_PS = 0.05


@dataclass(frozen=True)
class Field:
    """Les constantes du champ de force, toutes en unités réduites."""

    trunc: float = 3.0            # coût kT d'un recouvrement complet
    stiffness: float = 1.5        # kT, terme angulaire
    wiggle: float = 0.05          # jeu de la liaison de squelette, en sigma
    lef_wiggle: float = 0.2       # jeu de la liaison de cohésine
    wall: float = 30.0            # kT/sigma, raideur du mur sphérique


def bond_catalogue(legs) -> tuple[np.ndarray, dict[tuple[int, int], int]]:
    """Toutes les paires qu'une cohésine tiendra **un jour** dans cette trajectoire.

    OpenMM refuse de changer les *particules* d'une liaison après coup — seules
    les longueurs et les raideurs sont modifiables dans un contexte vivant. Une
    liaison par emplacement de cohésine, qu'on déplacerait au fil de l'extrusion,
    ne marche donc pas : `updateParametersInContext` lève « The set of particles
    in a bond has changed ».

    La parade est de déclarer d'emblée **toutes** les paires du parcours, avec une
    raideur nulle, et de n'allumer que celles qui sont tenues à l'instant présent.
    Le coût est nul : quelques milliers de liaisons à raideur zéro ne pèsent rien
    devant les termes non liés, et rien ne bouge dans le champ de force sinon deux
    scalaires par cohésine et par mise à jour.
    """
    seen: dict[tuple[int, int], int] = {}
    for frame in legs:
        for i, j in frame:
            if i < 0 or j < 0 or i == j:
                continue
            key = (int(min(i, j)), int(max(i, j)))
            if key not in seen:
                seen[key] = len(seen)
    pairs = np.array(sorted(seen, key=seen.get), dtype=np.int64).reshape(-1, 2)
    return pairs, seen


def _build(n: int, confine: float | None, pairs: np.ndarray, field: Field):
    import openmm as mm

    system = mm.System()
    for _ in range(n):
        system.addParticle(MASS)

    backbone = mm.HarmonicBondForce()
    k_bb = 2.0 * KT / field.wiggle**2      # E = ½k(r−r0)² = kT(r−r0)²/wiggle²
    for i in range(n - 1):
        backbone.addBond(i, i + 1, 1.0, k_bb)
    system.addForce(backbone)

    if field.stiffness > 0.0:
        angle = mm.CustomAngleForce("kang*(1 - cos(theta - 3.141592653589793))")
        angle.addGlobalParameter("kang", field.stiffness * KT)
        for i in range(n - 2):
            angle.addAngle(i, i + 1, i + 2, [])
        system.addForce(angle)

    # Répulsion molle, monotone, nulle et de dérivée nulle au contact. Une forme
    # en (1 − (r/rc)²)² serait plus douce mais sa force s'annule aussi en r = 0 :
    # deux monomères parfaitement superposés ne se sépareraient jamais.
    rep = mm.CustomNonbondedForce("trunc*(1 - r/rc)^2")
    rep.addGlobalParameter("trunc", field.trunc * KT)
    rep.addGlobalParameter("rc", 1.0)
    rep.setCutoffDistance(1.0)
    rep.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    for _ in range(n):
        rep.addParticle([])
    system.addForce(rep)

    if confine is not None:
        # Mur linéaire lissé : force de rappel constante au-delà du rayon, pas de
        # discontinuité à la traversée. Un mur harmonique deviendrait d'autant plus
        # raide qu'on s'en éloigne et rendrait le pas de temps instable au moment
        # précis où le polymère se replie.
        wall = mm.CustomExternalForce(
            "step(rr - Rc) * kw * (sqrt((rr - Rc)^2 + 0.01) - 0.1);"
            "rr = sqrt(x*x + y*y + z*z + 0.0001)"
        )
        wall.addGlobalParameter("kw", field.wall * KT)
        wall.addGlobalParameter("Rc", confine)
        for i in range(n):
            wall.addParticle(i, [])
        system.addForce(wall)

    lef = mm.HarmonicBondForce()
    for i, j in pairs:
        lef.addBond(int(i), int(j), 1.0, 0.0)   # déclarée, éteinte
    system.addForce(lef)

    return system, lef


def start(n: int, confine: float | None, rng: np.random.Generator) -> np.ndarray:
    """Conformation de départ : marche aléatoire de pas unité, réfléchie sur la paroi.

    La tentation est d'écrire une marche libre puis de la **comprimer** d'un
    facteur d'échelle jusqu'à ce qu'elle tienne dans la sphère. Le rayon de
    giration paraît alors correct dès le premier bloc — mesuré, 7,02 sigma au pas
    2 000 et 7,05 au pas 40 000, une courbe parfaitement plate — et on en conclut
    que l'équilibrage est immédiat. C'est faux : la compression divise *toutes*
    les distances par le même facteur, donc elle écrase la structure interne
    exactement autant que la structure globale, et `Rg` ne peut pas le voir.

    Une marche réfléchie garde des pas de longueur unité à toutes les échelles.
    Elle n'est pas pour autant un fondu équilibré ; ce qu'on en attend, c'est
    qu'elle ne mente pas sur le point de départ. La convergence, elle, se mesure
    sur `R(s)` aux grandes séparations (§ `equilibration`).
    """
    steps = rng.normal(size=(n, 3))
    steps /= np.linalg.norm(steps, axis=1, keepdims=True)
    if confine is None:
        x = np.cumsum(steps, axis=0)
        return x - x.mean(axis=0)

    x = np.empty((n, 3))
    x[0] = rng.normal(size=3) * confine / 4.0
    for i in range(1, n):
        cand = x[i - 1] + steps[i]
        for _ in range(12):
            if np.linalg.norm(cand) <= confine:
                break
            v = rng.normal(size=3)
            cand = x[i - 1] + v / np.linalg.norm(v)
        else:                                  # coincé dans un coin : on rentre
            cand = x[i - 1] * (1.0 - 1.0 / max(np.linalg.norm(x[i - 1]), 1e-9))
        x[i] = cand
    return x - x.mean(axis=0)


def _switch(lef, pairs, index, frame, live: set[int], k_lef: float) -> set[int]:
    """Allume les paires tenues maintenant, éteint celles qui ne le sont plus.

    Ne touche qu'aux liaisons dont l'état **change**. Reparcourir le catalogue
    entier à chaque instantané coûterait quelques millions d'appels Python pour un
    résultat identique.
    """
    want = {
        index[(int(min(i, j)), int(max(i, j)))]
        for i, j in frame
        if i >= 0 and j >= 0 and i != j
    }
    for b in live - want:
        lef.setBondParameters(b, int(pairs[b, 0]), int(pairs[b, 1]), 1.0, 0.0)
    for b in want - live:
        lef.setBondParameters(b, int(pairs[b, 0]), int(pairs[b, 1]), 1.0, k_lef)
    return want


def simulate(
    n: int,
    legs: tuple[np.ndarray, ...] | list[np.ndarray],
    *,
    confine: float | None,
    field: Field = Field(),
    md_per_step: int = 200,
    relax: int = 20_000,
    seed: int = 0,
    threads: int = 1,
    progress=None,
) -> np.ndarray:
    """Une trajectoire de Langevin dont les liaisons de cohésine suivent `legs`.

    Renvoie `(len(legs), n, 3)` en unités de `sigma`.

    `relax` est le nombre de pas de mise en place, *avant* le premier instantané.
    Le choisir demande de regarder `R(s)` aux grandes séparations et pas `Rg`,
    pour la raison expliquée dans `start` : le rodage tourne **avec** les
    cohésines du premier instantané, sans quoi la chaîne relaxerait nue puis
    encaisserait d'un coup une vingtaine de liaisons.

    L'aléa entre dans deux endroits seulement : la conformation initiale et le
    thermostat, tous deux par `seed`. Deux exécutions de même graine et de même
    `legs` sont identiques à l'arithmétique flottante près.
    """
    import openmm as mm
    import openmm.unit as u

    rng = np.random.default_rng(seed)
    pairs, index = bond_catalogue(legs)
    system, lef = _build(n, confine, pairs, field)
    k_lef = 2.0 * KT / field.lef_wiggle**2

    integrator = mm.LangevinMiddleIntegrator(
        300.0 * u.kelvin,
        FRICTION_PER_PS / u.picosecond,
        TIMESTEP_FS * u.femtosecond,
    )
    integrator.setRandomNumberSeed(int(seed) % (2**31 - 1) + 1)
    platform = mm.Platform.getPlatformByName("CPU")
    context = mm.Context(system, integrator, platform, {"Threads": str(max(1, threads))})

    context.setPositions(start(n, confine, rng))

    live: set[int] = set()
    if len(legs):
        live = _switch(lef, pairs, index, legs[0], live, k_lef)
        lef.updateParametersInContext(context)
    mm.LocalEnergyMinimizer.minimize(context, 1.0, 5_000)
    if relax:
        integrator.step(relax)

    out = np.empty((len(legs), n, 3), dtype=np.float32)
    for t, frame in enumerate(legs):
        live = _switch(lef, pairs, index, frame, live, k_lef)
        lef.updateParametersInContext(context)
        integrator.step(md_per_step)
        state = context.getState(getPositions=True, enforcePeriodicBox=False)
        out[t] = state.getPositions(asNumpy=True).value_in_unit(u.nanometer)
        if progress is not None:
            progress(t + 1, len(legs))

    if not np.isfinite(out).all():
        raise RuntimeError("la trajectoire a divergé — pas de temps ou raideur à revoir")
    return out


def equilibration(
    n: int,
    *,
    confine: float | None,
    seps: tuple[int, ...] = (10, 100, 400, 1000),
    field: Field = Field(),
    blocks: int = 20,
    per_block: int = 5_000,
    legs: tuple[np.ndarray, ...] | list[np.ndarray] | None = None,
    seed: int = 0,
    threads: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Courbe de mise en place : `Rg` et `R(s)` bloc par bloc.

    Renvoie `(rg, rs)` de formes `(blocks,)` et `(blocks, len(seps))`. Sert à
    choisir `relax` **en le regardant**, et à choisir la bonne observable pour ça :
    `Rg` est aveugle à la structure interne (§ `start`), alors que `R(s)` aux
    grandes séparations est justement le mode le plus lent. Démarrer les
    instantanés avant sa convergence laisserait la conformation initiale décider
    du raccord de fin de semaine.

    `legs` permet de mesurer la mise en place **avec** les cohésines, qui
    réorganisent la chaîne localement et changent donc la réponse. Sans lui, la
    chaîne relaxe nue.
    """
    frames = (
        [np.full((1, 2), -1, dtype=np.int64)] * blocks
        if legs is None
        else [legs[min(i, len(legs) - 1)] for i in range(blocks)]
    )
    traj = simulate(
        n,
        frames,
        confine=confine,
        field=field,
        md_per_step=per_block,
        relax=0,
        seed=seed,
        threads=threads,
    )
    rg = np.array([float(np.sqrt(((f - f.mean(0)) ** 2).sum(1).mean())) for f in traj])
    rs = np.array(
        [[float(np.linalg.norm(f[s:] - f[:-s], axis=1).mean()) for s in seps] for f in traj]
    )
    return rg, rs
