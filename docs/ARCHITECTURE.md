# Geno — Architecture technique

Décisions arrêtées avant le code. Ce document est révisé en semaine 4, après le spike de
rendu, puis figé.

---

## 1. Structure du dépôt

```
geno/
├── pipeline/              # Python — données → structures 3D
│   ├── rules/             # Snakemake
│   ├── geno_pipeline/
│   │   ├── hic/           # cooler, cooltools : balancing, eigs, insulation, boucles
│   │   ├── recon/         # contact→distance, recuit, ensembles
│   │   ├── polymer/       # polychrom / OpenMM : extrusion de boucles
│   │   └── export/        # écriture .g3d
│   └── envs/              # conteneurs
├── packages/
│   ├── core/              # TS — modèle de données, index génomique, mapping coord↔3D
│   ├── g3d/               # TS — lecture du format, streaming, octree
│   ├── viewer/            # TS — moteur Three.js, LOD, picking GPU
│   └── annotate/          # TS — requêtes d'annotation, composition de la fiche
├── apps/
│   └── web/               # l'application, les deux modes
└── docs/
```

Monorepo pnpm. Le pipeline Python et le front TypeScript ne se parlent que par le format
`.g3d` et des fichiers d'annotation indexés — jamais par un couplage de code.

---

## 2. Les deux modes : un noyau, deux coquilles

Le mode Découverte **est** le mode Recherche, avec un préréglage curaté et une couche
narrative par-dessus. Pas de fork, pas de duplication.

```
packages/viewer  ─┐
packages/core    ─┼─→  apps/web  ─┬─→ /recherche   (feature flags : import, params, export, ensemble)
packages/annotate ┘               └─→ /decouverte  (preset figé + narration + glossaire)
```

Un contexte `GenoMode` porte les flags. Règle : **aucun `if (mode === ...)` dans
`packages/*`.** Le branchement vit uniquement dans `apps/web`. Si un package a besoin de
connaître le mode, c'est que l'abstraction est fausse.

---

## 3. Format `.g3d`

Binaire, multi-résolution, pensé pour la requête HTTP Range. Un fichier = une structure (ou
un médoïde d'ensemble) pour une lignée cellulaire et un assemblage.

```
┌─ header (JSON, longueur préfixée) ───────────────────────────┐
│  assembly: "GRCh38" | "CHM13v2.0"                            │
│  cell_type, ploidy, n_structures                             │
│  levels: [{ bin_size, n_beads, offset, length, octree_off }] │
│  provenance: { accessions[], pipeline_version, sha256 }      │
│  evidence: "measured"|"simulated"|"deterministic"|"predicted" │
└──────────────────────────────────────────────────────────────┘
┌─ par niveau ─────────────────────────────────────────────────┐
│  positions   Float32Array  [x,y,z] × n_beads                 │
│  variability Float32Array  écart-type sur l'ensemble          │
│  chrom_index Uint32Array   bille → (chrom, start, haplotype)  │
│  octree      nœuds → plages d'octets                         │
└──────────────────────────────────────────────────────────────┘
```

Le champ `evidence` est structurel, pas cosmétique : c'est lui qui pilote le code couleur
exigé par le principe n°2 de la feuille de route. `predicted` est réservé aux positions issues
d'un modèle séquence → contact (Akita, Orca, C.Origami) — hors périmètre v1, mais le champ le
prévoit pour qu'une telle position ne puisse jamais se faire passer pour une mesure.
Voir [`DATA_SOURCES.md` § 7](DATA_SOURCES.md).

---

## 4. Budget de niveaux de détail

Génome diploïde = 6,2 Gb. Positions = 3 × float32 = 12 octets par bille.

| Niveau | Taille de bin | Billes (diploïde) | Positions | Chargement |
|--------|---------------|-------------------|-----------|-----------|
| Noyau | 1 Mb | 6 200 | 74 ko | global, immédiat |
| Territoire | 250 kb | 24 800 | 298 ko | global |
| Compartiment | 100 kb | 62 000 | 744 ko | global |
| TAD | 25 kb | 248 000 | 3,0 Mo | par octant visible |
| Boucle | 10 kb | 620 000 | 7,4 Mo | par octant visible |
| Fibre | 1 kb | 6,2 M | — | fenêtre visible seulement |
| Séquence | 1 pb | — | — | géométrie procédurale |

Ce tableau est de l'**arithmétique d'octets, vérifiée**. Ce n'est pas encore une mesure de
rendu : la colonne qui manque est le temps par frame à chaque palier, sur trois cibles
matérielles. Le spike de la semaine 4 est écrit et sa correction est prouvée
([`VALIDATION.md` § S4](VALIDATION.md)), mais l'environnement de développement ne rend que via
SwiftShader, un rasteriseur logiciel — ses images par seconde ne disent rien d'un GPU.

`pnpm build` produit `dist/spike.html` ; l'ouvrir sur GPU desktop, iGPU portable et téléphone
rend les trois tableaux qui manquent. Le go/no-go de la semaine 4 attend ces chiffres.

Une fois le budget arrêté, il est **tenu par des tests de non-régression** : un commit qui le
fait dépasser casse la CI.

---

## 5. Rendu

**Imposteurs de sphères, pas de maillages.** Un quad orienté caméra, le fragment shader
calcule la normale de sphère et la profondeur. C'est ce que font Mol* et NGL. Une vraie
géométrie de sphère, même grossière à 80 triangles, ferait 80 millions de triangles pour un
million de billes ; un imposteur en fait deux millions.

Le point critique est que `gl_FragDepth` vienne du **point d'impact réel**, pas du quad. Sans
ça, deux billes qui s'interpénètrent se découpent selon l'arête du quad au lieu de la courbe
d'intersection. Implémenté et testé dans `packages/viewer` — le test discriminant est que le
résultat ne dépend pas de l'ordre de dessin.

**Picking GPU.** Rendu des identifiants dans un attachement `R32UI`, en même temps que la
couleur (cible multiple : une seule passe géométrique, ce qui à un million d'instances n'est
pas une nuance). Lecture d'un pixel au clic, coût O(1) quel que soit le nombre d'objets — le
raycasting CPU de Three.js est inutilisable à cette échelle. Vérifié exact parmi 250 005
instances, identifiants sur toute la plage 32 bits.

Les shaders sont écrits en GLSL ES 3.0 brut dans `packages/viewer/src/impostors.ts` et se
transposent tels quels dans un `RawShaderMaterial` de Three.js.

**Fibre et boucles.** Géométrie tubulaire générée depuis la polyligne, tessellation adaptative
selon la distance caméra. Les ancres de boucles sont des arcs explicites, pas des artefacts
de la polyligne.

**Lisibilité du volume.** Occlusion ambiante en espace écran : sans elle, un nuage de billes
sphériques est illisible. Profondeur de champ réservée aux vues macro.

**Transitions entre niveaux.** Fondu croisé sur un recouvrement de distance caméra, avec
préchargement du niveau suivant. L'utilisateur ne doit jamais voir un niveau apparaître.

---

## 6. Mapping coordonnée ↔ position 3D

Les deux sens sont nécessaires, et les deux doivent être rapides.

- **Génomique → 3D** : arbre d'intervalles sur `(chrom, start, end)` → indices de billes →
  décalages dans le `Float32Array`. Sert à la recherche, aux signets, au brossage depuis
  `igv.js`.
- **3D → génomique** : l'identifiant rendu par le picking GPU indexe directement
  `chrom_index`. Sert au clic.

Les deux haplotypes sont distincts : une coordonnée génomique correspond à **deux** positions
3D. L'interface doit le montrer, pas le masquer — c'est une information biologique réelle
(expression allélique, inactivation du X, empreinte parentale).

---

## 7. Vues liées

Trois vues, une seule sélection :

```
        ┌──────────────┐
        │  sélection   │   (chrom, start, end, haplotype)
        └──────┬───────┘
     ┌─────────┼─────────┐
  3D viewer  igv.js   HiGlass
```

Un store unique émet la sélection ; chaque vue s'y abonne et la publie. Pas de
synchronisation deux-à-deux — c'est la source des boucles infinies.

---

## 8. Socle 1D — index d'intervalles et mesures

Un intervalle `[s, e)` chevauche une requête `[qs, qe)` si et seulement si `s < qe` **et**
`e > qs`. La première condition se résout par dichotomie sur les débuts triés. La seconde est
le piège : les intervalles sont triés par début, pas par fin, donc un candidat peut se trouver
arbitrairement loin en arrière. Un balayage naïf devient O(n) dès qu'un gène long traîne en
tête de chromosome — et il y en a : DMD fait 2,2 Mb, CNTNAP2 2,3 Mb.

D'où un **maximum de fin par bloc** de 512 entrées : un bloc dont le maximum est ≤ `qs` ne
peut rien contenir et se saute d'un seul test. Le filtrage des blocs est une opération numpy
vectorisée, pas une boucle Python.

Mesures sur **1 000 000 d'intervalles** (dont 0,2 % de très longs), 2 000 requêtes par ligne :

| Fenêtre | Features | Index seul (méd / p99) | Index + attributs (méd / p99) |
|---------|---------:|------------------------|-------------------------------|
| 3 kb — un gène | 42 | 0,028 / 0,074 ms | 0,204 / 0,472 ms |
| 30 kb — un gène + flancs | 87 | 0,028 / 0,069 ms | 0,338 / **0,819 ms** |
| 300 kb — un TAD | 530 | 0,031 / 0,082 ms | 1,535 / 6,090 ms |
| 3 Mb — un compartiment | 5 094 | 0,044 / 0,110 ms | 17,8 / 42,6 ms |

**Cible de la semaine 2 tenue** : p99 à 0,819 ms à l'échelle d'un locus, attributs compris.
L'index seul reste sous 0,11 ms à toutes les échelles, y compris celle d'un compartiment.

### Ce que ça dit pour la semaine 9

La latence ne suit pas la taille du magasin, elle suit le **nombre de features ramenées** :
environ 3 µs chacune, dominés par le parsing JSON des attributs. Conséquence directe sur le
format `.g3d` : les attributs ne doivent pas y être du JSON par enregistrement, mais des
**colonnes** — tableaux numpy de codes catégoriels plus une table de chaînes. On lit alors les
colonnes utiles et on ne parse rien.

Ce n'est pas un problème pour la v1 de la CLI : une requête qui ramène un compartiment entier
est un export, pas une interaction, et la fiche de la semaine 13 clique **une** bille.

---

## 9. Ce qui ne rentre pas dans le navigateur

La reconstruction (recuit simulé sur ~10 000 billes, × 200 structures) et la simulation
polymère restent côté pipeline. Le mode Recherche déclenche un travail serveur et récupère
un `.g3d`. On n'essaiera pas de porter OpenMM en WASM.

Un ensemble complet (200 × 620 000 billes) pèse ~1,5 Go : jamais chargé entier. On streame le
médoïde plus un scalaire de variabilité par bille ; l'ensemble ne se charge qu'à la demande,
en résolution grossière.

---

## 10. Modèle de noyau — billes TAD, territoires, périphérie

### La bille

Une bille est un TAD. Encore faut-il dire lequel : Dixon 2012 en compte ~2 200 sur le génome
haploïde, taille moyenne ~880 kb ; Rao 2014 en compte 9 274, médiane 185 kb. Ce ne sont pas
deux mesures du même objet, ce sont deux définitions. Le modèle de noyau travaille à
**750 kb par bille**, donc à l'échelle Dixon, et le génome diploïde tient en ~8 100 billes.

La chromatine ayant une densité locale à peu près constante, le volume d'une bille suit sa
longueur et son rayon la racine cubique : `r ∝ L^(1/3)`. La constante vient d'une seule
quantité physique, la **fraction volumique** `phi` occupée par les billes — 0,30 par défaut,
dans la fourchette 12–52 % mesurée par ChromEMT (Ou 2017).

### Le rayon nucléaire n'est qu'une unité

Le nombre de billes et la fraction volumique déterminent le rapport `r/R` :

```
N · (4/3)π r³ = phi · (4/3)π R³     ⟹     r / R = (phi / N)^(1/3)
```

Une fois `phi` et `N` fixés, toute la géométrie est en unités de `R`. Changer `R` change les
micromètres affichés et rien d'autre. Les deux seuls paramètres physiques du modèle sont donc
`phi` et la résolution.

### Trois façons de dire « à la lamina », et elles ne coïncident pas

Une bille **touche** l'enveloppe si l'écart entre sa surface et celle-ci est sous un demi-rayon.
Cette bande est plus fine que l'écart entre deux couches empilées (~1,63 r), donc elle ne
contient qu'une monocouche. Deux bornes l'encadrent :

- **part à densité uniforme** — la part de volume de la coquille dans la boule accessible aux
  centres, soit `[(R−r)³ − (R−1,5r)³] / (R−r)³`, environ `1,5 · r/(R−r)` ;
- **borne d'empilement** — au plus `eta · 4(R−r)²/r²` billes peuvent toucher à la fois, avec
  `eta ≤ 0,9069` la densité hexagonale.

Les deux sont **linéaires en r**, donc dépendent de la résolution du modèle :

| Résolution | Billes | Rayon | Densité uniforme | Borne d'empilement | Niveau `.g3d` |
|-----------:|-------:|------:|-----------------:|-------------------:|--------------:|
| 3 Mb | 2 021 | 264,7 nm | 8,2 % | 57 % | 24 ko |
| **750 kb** | **8 083** | **166,8 nm** | **5,1 %** | **38 %** | **97 ko** |
| 250 kb | 24 248 | 115,6 nm | 3,5 % | 27 % | 291 ko |
| 100 kb | 60 621 | 85,2 nm | 2,6 % | 20 % | 727 ko |
| 10 kb | 606 208 | 39,5 nm | 1,2 % | 9 % | 7,3 Mo |

Conséquence directe : **« 35 % du génome est en LAD » et « 35 % du génome touche la lamina »
ne sont pas la même phrase**, et la seconde ne découle pas de la première. Atteindre la borne
d'empilement suppose un noyau à croûte dense et intérieur creux ; rester à densité uniforme
plafonne à quelques pour cent. Une « fraction de LADs à la lamina » n'est donc pas comparable
entre deux modèles de granularité différente sans dire laquelle.

La dernière colonne donne le poids du niveau dans le `.g3d` de la semaine 9, à 12 octets par
bille (§ 4) : à 750 kb le niveau fait **97 ko**, largement sous le budget de 500 ko de la
semaine 10. La résolution du modèle de noyau tombe entre les paliers « Noyau » (1 Mb) et
« Territoire » (250 kb) du § 4 ; c'est l'échelle TAD de Dixon qui l'a fixée, pas le budget.

Sur les totaux : le caryotype 46,XX de GM12878 fait **6,06 Gb** d'assemblage primaire, deux
exemplaires de chr1–22 plus deux X. Le « 6,2 Gb » du § 4 est la valeur ronde usuelle, scaffolds
et chrY compris.

### Hiérarchie des contraintes

Trois conditions dures, une préférence molle, et l'ordre compte :

| | Contrainte | Nature |
|---|---|---|
| 1 | volume exclu, `d ≥ r_i + r_j` | dure |
| 2 | longueur de liaison **maximale**, `d ≤ stretch · (r_i + r_j)` | dure, unilatérale |
| 3 | confinement, `\|x\| ≤ R − r` | dure |
| 4 | rappel radial vers un rayon cible ordonné par contenu LAD | molle |

La liaison est unilatérale : elle empêche la chaîne de casser, pas les billes de se rapprocher.
Une longueur *imposée* se bat contre le volume exclu dans les replis serrés, et c'est le volume
exclu qui cède — or c'est lui la condition.

Le rappel radial vise un rayon par **rang** de contenu LAD, pas par valeur : la bille de rang
`q` vise `q^(1/3)`, la loi des rayons d'une sphère uniforme. Le biais trie sans tasser. Le
polissage final le coupe entièrement.

### Topologie : ce que le recuit ne peut pas faire

Une relaxation sous contraintes ne fait **jamais** se croiser deux chaînes. La topologie du
noyau est donc celle de son initialisation, et la territorialité est une **entrée** du modèle,
pas un résultat. Ce qui la justifie est la mitose : les chromosomes se décondensent là où la
télophase les a laissés. Le seul énoncé vérifiable est que le recuit la conserve.

Le solveur fait croître les rayons pendant le recuit plutôt que de partir à taille pleine
(Lubachevsky–Stillinger). D'où une condition d'étanchéité : à l'échelle `s`, deux billes liées
s'excluent au rayon `s·r` mais restent séparées d'au plus `stretch·2r`, donc une troisième
bille passe entre elles dès que **`s ≤ stretch / 2`**. En dessous du seuil une chaîne traverse
une liaison, et le défaut est difficile à défaire. Le défaut n'est pas construit : il est
mesuré (`docs/VALIDATION.md` § S6).

---

## 11. Ensembles — `ensemble.zarr`

Le principe n° 1 du projet dit qu'une structure Hi-C unique est un artefact statistique. Un
ensemble en est la conséquence dans le format : **un génome, N repliements**.

```
ensemble.zarr/
  .zattrs           caryotype, assemblage, provenance des longueurs, lad_source,
                    lad_seed, first_seed, réglages de construction, avertissement
  beads/            copy_id, start, end, radius, lad, acrocentric     ← écrits UNE fois
  coords    (N, n, 3) float32, un chunk par structure
  quality/          seed, max_overlap, bond_stretch, outside, shakes  (N,)
  done      (N,) bool
```

Le découpage en billes et la piste LAD ne portent **pas** d'indice de structure. Ce n'est pas
une économie de place — 200 copies du génome pèseraient 10 Mo — c'est ce qui rend l'objet
interprétable : si chaque structure portait le sien, la variabilité mesurée mélangerait la
variabilité de repliement et celle du génome, sans moyen de les séparer après coup. La
génération le vérifie à chaque structure plutôt que de s'y fier (§ `ensemble/generate`).

`done` est écrit **après** les coordonnées : une génération interrompue ne prétend jamais
contenir une structure à moitié écrite, et `--resume` reprend où elle s'est arrêtée.

### Ce qu'un ensemble permet de mesurer, et ce qu'il interdit

Deux noyaux recuits séparément ne partagent **aucun repère** : ni orientation — rien ne
distingue un axe dans une sphère — ni placement des territoires, puisque chr7:a atterrit
ailleurs à chaque tirage. Une variance par bille en x, y, z est donc un nombre sans objet, et
l'alignement de Procruste ne la sauve pas. C'est mesuré, pas supposé : après alignement
optimal, **88 % de l'écart au hasard subsiste** ([`VALIDATION.md` § S7](VALIDATION.md)).

Restent les grandeurs invariantes par rotation, et elles seules :

| Grandeur | Ce qu'elle porte |
|---|---|
| position radiale par bille | profondeur nucléaire — la quantité comparable entre structures |
| distances par paires | carte de contacts, P(s), fraction trans |
| par copie | rayon de giration, profondeur du centroïde |

Conséquence pour le viewer (semaine 10 et suivantes) : ce qui se streame d'un ensemble n'est
pas « la » structure plus un nuage de points, c'est **le médoïde plus un scalaire de
variabilité par bille** — et ce scalaire est une variabilité *de profondeur*, pas de position.
L'interface doit dire laquelle.

## 12. Échelle fine — extrusion de boucles sous le TAD

La semaine 6 descend jusqu'à la bille de 750 kb, et pas plus bas. En dessous, le modèle de
noyau n'a rien à dire : ses billes sont des sphères sans intérieur. La semaine 8 ouvre cet
intérieur sur une région de quelques mégabases, et le fait avec un mécanisme plutôt qu'une
interpolation — **l'extrusion de boucles** (Fudenberg 2016), qui est à ce jour la seule
hypothèse qui explique d'un coup les domaines, les points d'angle et la règle de convergence
des motifs CTCF.

```
geno_pipeline/fine/
  region.py      la région, ses monomères, ses sites CTCF orientés, la vérité plantée
  extrusion.py   modèle 1D : où sont les pieds des cohésines
  polymer.py     dynamique de Langevin (OpenMM), liaisons pilotées par l'extrusion
  observe.py     carte de contacts, insulation, points d'angle, P(s), raccord R(s)
  run.py         assemblage, réplicats parallèles, sauvegarde `.npz` + sidecar JSON
```

### La règle de convergence, traduite en deux tableaux

Tout le pouvoir explicatif du modèle tient dans une asymétrie :

| Motif | Arrête une jambe qui va… |
|---|---|
| **forward (+)** | vers la **gauche** |
| **reverse (−)** | vers la **droite** |

Un cohésine chargé entre un `+` en `i` et un `−` en `j > i` bloque donc ses deux jambes, et la
paire `(i, j)` accumule des contacts : c'est un point d'angle. Une paire **divergente** (`−`
puis `+`) n'arrête rien. Inverser ces deux lignes produirait des boucles exactement là où Rao
2014 n'en voit pas — c'est une prédiction binaire, et c'est ce qui rend le modèle réfutable.

### Le rayon d'un monomère n'est pas un paramètre de la semaine 8

Il sort de la loi de la semaine 6 (§ 10) : volume de chromatine proportionnel aux paires de
bases, constante fixée par la fraction volumique `phi`. À 2 kb par monomère, `r = 23,1 nm`,
donc `sigma = 46,2 nm` de diamètre. Le confinement suit : `n` monomères occupent la sphère de
rayon `R = r·(n/phi)^(1/3)`, c'est-à-dire **exactement le volume que le modèle de noyau alloue
à cette quantité de chromatine**.

C'est ce qui rend le raccord falsifiable. Si la semaine 8 avait le droit de choisir la taille
de ses monomères, on pourrait toujours faire coïncider les deux modèles et la coïncidence ne
dirait rien.

### Champ de force, et ce que chaque terme coûte

| Terme | Forme | Valeur | Pourquoi |
|---|---|---|---|
| squelette | harmonique | `wiggle` 0,05 σ | la chaîne ne casse pas |
| rigidité | `k(1 − cos θ)` | 1,5 kT | longueur de persistance de quelques monomères |
| volume exclu | `trunc·(1 − r/σ)²` | 3 kT | **franchissable** — voir ci-dessous |
| confinement | mur linéaire lissé | 30 kT/σ | densité du noyau, pas une forme |
| cohésines | harmonique | `wiggle` 0,2 σ | plus molle que le squelette, pour encaisser les mises à jour |

Le volume exclu est volontairement mou : à l'échelle de temps de l'extrusion, la
topoisomérase II laisse passer les brins, et un volume exclu dur figerait des enlacements que
la cellule défait. Le prix est explicite — **la topologie de la région n'est pas une
prédiction du modèle**.

### Une sphère, pas une boîte périodique

Une boîte périodique n'a pas de paroi, ce qui serait préférable, mais elle remplace
l'environnement du segment par ses propres images : `R(s)` y sature à la taille de la boîte.
Pour 2 000 monomères à `phi = 0,30`, la boîte fait 13,8 σ de côté pour un blob de 9,4 σ de
rayon — l'artefact tombe en plein dans la fenêtre où le raccord se mesure. La sphère a une
paroi, c'est un défaut qu'on énonce, mais elle fixe le bon volume.

### OpenMM ne déplace pas une liaison

`updateParametersInContext` change les longueurs et les raideurs, jamais les particules d'une
liaison. Une liaison par cohésine, qu'on déplacerait au fil de l'extrusion, lève « The set of
particles in a bond has changed ». D'où `bond_catalogue` : **toutes** les paires du parcours
sont déclarées d'emblée à raideur nulle, et seules celles tenues à l'instant présent sont
allumées. Quelques milliers de liaisons éteintes ne pèsent rien devant les termes non liés.

Corollaire de format : les instantanés d'extrusion gardent **un emplacement par cohésine**,
`-1` quand elle est décrochée, au lieu de la liste compacte des seules cohésines chargées. Si
l'identité des emplacements changeait d'un instantané à l'autre, une liaison sauterait d'une
paire à une autre sans rapport et la dynamique encaisserait un choc là où l'extrusion n'avance
que d'un cran.

### Le raccord avec le noyau entier

Les deux modèles se recouvrent entre 750 kb — une bille de la semaine 6 — et l'étendue de la
région fine, **strictement exclue** : sur `n` monomères, la séparation `n` n'existe pas. Sur
cette fenêtre les deux prétendent à la même grandeur, `R(s)`, la distance spatiale moyenne
entre deux morceaux de chromatine séparés de `s` paires de bases, et ils ne partagent aucun
paramètre ajusté. Le rapport des deux courbes est le critère de fin de semaine
([`VALIDATION.md` § S8](VALIDATION.md)).
