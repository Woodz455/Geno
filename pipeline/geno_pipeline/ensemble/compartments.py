"""Compartiments A/B du modèle, lus dans sa propre carte de contacts.

Aucun compartiment n'est donné au modèle de noyau : ses seuls ingrédients sont la chaîne,
le volume exclu, le confinement et un terme radial ordonné par la piste LAD (§ S6). Si un
damier A/B apparaît dans sa carte de contacts, il en **émerge** — et c'est exactement ce
qu'on mesure ici, avec les appelants validés de la semaine 3 (`hic.features`) appliqués à
la carte qu'une expérience Hi-C aurait tirée de l'ensemble.

Trois choix, chacun pour une raison :

- **Carte haploïde.** Le Hi-C ne distingue pas les homologues : les contacts des deux copies
  d'un chromosome tombent dans les mêmes bins, contacts entre homologues compris. La carte
  agrège donc les copies par chromosome, X compris (Xa et Xi), comme le ferait le séquençage.
- **Orientation par la piste LAD.** Le signe du premier vecteur propre est arbitraire ; en
  vrai, on l'oriente par le contenu GC ou la densité génique, qu'on ne peut pas lire ici. On
  oriente donc A du côté pauvre en LAD. Conséquence : la corrélation négative entre E1 et
  LAD est en partie **imposée par le signe** ; ce qui ne l'est pas, c'est son amplitude, et
  la force du damier.
- **Xi sans compartiment.** L'X inactif n'a pas de compartiments A/B (il se replie en deux
  superdomaines, Rao 2014) : ses billes reçoivent « absent », pas l'E1 de l'Xa. Le modèle,
  lui, ne replie pas l'Xi autrement que l'Xa — c'est dit, pas corrigé.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .stats import contact_pairs


@dataclass(frozen=True)
class Compartments:
    """E1 par bille du modèle, et ce qui dit s'il veut dire quelque chose."""

    e1: np.ndarray            # (n_billes,) diploïde ; NaN pour l'Xi et les bins non appelés
    bins: pd.DataFrame        # bins haploïdes : chrom, start, end, E1, lad, radial
    strength: float           # (AA + BB) / (AB + BA), quintiles extrêmes, O/E cis — croisée
    strength_in_sample: float # la même, E1 et carte tirés des mêmes structures : biaisée
    saddle: np.ndarray        # (k, k) O/E moyen par paire de groupes de E1, croisé
    split_half_r: float       # E1 de la première moitié des structures contre la seconde
    r_lad: float              # Pearson E1 ~ fraction LAD, sur les bins appelés
    r_radial: float           # Pearson E1 ~ position radiale moyenne
    fraction_a: float
    called: int               # bins avec un E1
    cutoff: float
    n_structures: int

    def __str__(self) -> str:
        return (
            f"{self.called:,} bins appelés sur {len(self.bins):,} · A {self.fraction_a:.0%} · "
            f"force du damier {self.strength:.2f} (croisée ; {self.strength_in_sample:.2f} sur "
            f"l'échantillon même) · moitié contre moitié r = {self.split_half_r:+.3f} · "
            f"E1 ~ LAD r = {self.r_lad:+.2f} · E1 ~ rayon r = {self.r_radial:+.2f}"
        )


def haploid_bins(labels: tuple[str, ...], copy_id: np.ndarray, start: np.ndarray,
                 end: np.ndarray) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    """Bin haploïde de chaque bille, la table des bins, et le masque des billes de l'Xi.

    Les copies d'un même chromosome sont découpées de la même façon — on le **vérifie** au
    lieu de le supposer : un découpage différent d'un homologue à l'autre ferait tomber
    dans un même bin des billes qui ne couvrent pas la même séquence.
    """
    chrom_of = [lab.split(":", 1)[0] for lab in labels]
    order: list[str] = list(dict.fromkeys(chrom_of))
    bin_of = np.full(len(copy_id), -1, dtype=np.int64)
    rows = []
    base = 0
    for chrom in order:
        copies = [c for c, ch in enumerate(chrom_of) if ch == chrom]
        ref = copies[0]
        m_ref = copy_id == ref
        s_ref, e_ref = start[m_ref], end[m_ref]
        for c in copies:
            m = copy_id == c
            if not (np.array_equal(start[m], s_ref) and np.array_equal(end[m], e_ref)):
                raise ValueError(f"{labels[c]} n'est pas découpée comme {labels[ref]}")
            bin_of[m] = base + np.arange(m.sum())
        rows.append(pd.DataFrame({"chrom": chrom, "start": s_ref, "end": e_ref}))
        base += int(m_ref.sum())
    xi = np.array([labels[c].endswith(":Xi") for c in copy_id])
    return bin_of, pd.concat(rows, ignore_index=True), xi


def count_contacts(coords: np.ndarray, radius: np.ndarray, bin_of: np.ndarray,
                   n_bins: int, cutoff: float) -> pd.DataFrame:
    """Pixels (bin1 ≤ bin2, compte) agrégés sur toutes les structures données."""
    keys = []
    for x in coords:
        pairs = contact_pairs(x, radius, cutoff)
        if len(pairs) == 0:
            continue
        a, b = bin_of[pairs[:, 0]], bin_of[pairs[:, 1]]
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        keys.append(lo * n_bins + hi)
    if not keys:
        return pd.DataFrame({"bin1_id": [], "bin2_id": [], "count": []})
    uniq, counts = np.unique(np.concatenate(keys), return_counts=True)
    return pd.DataFrame({"bin1_id": uniq // n_bins, "bin2_id": uniq % n_bins,
                         "count": counts.astype(np.int64)})


def _call(bins: pd.DataFrame, pixels: pd.DataFrame, phasing: pd.DataFrame, workdir: Path):
    """E1 et matrices cis équilibrées, par les appelants de la semaine 3."""
    import cooler

    from ..hic.features import balance, compartments
    from ..hic.synthetic import write_cool

    path = write_cool(workdir / "modele.cool", bins, pixels)
    clr = cooler.Cooler(str(path))
    balance(clr)
    e1 = compartments(clr, phasing_track=phasing)["E1"].to_numpy()
    return e1, clr


def _saddle(clr, e1: np.ndarray, groups: int = 5) -> np.ndarray:
    """O/E cis moyen entre groupes de E1 (quantiles). Le coin AA/BB contre AB dit la force."""
    ok = np.isfinite(e1)
    edges = np.quantile(e1[ok], np.linspace(0, 1, groups + 1))
    g = np.full(len(e1), -1)
    g[ok] = np.clip(np.searchsorted(edges, e1[ok], side="right") - 1, 0, groups - 1)
    total = np.zeros((groups, groups))
    n = np.zeros((groups, groups))
    offset = 0
    for chrom in clr.chromnames:
        m = clr.matrix(balance=True).fetch(chrom)
        k = len(m)
        gi = g[offset:offset + k]
        offset += k
        oe = np.full_like(m, np.nan)
        for d in range(2, k):        # les deux premières diagonales sont ignorées par ICE
            diag = np.diagonal(m, d)
            mean = np.nanmean(diag) if np.isfinite(diag).any() else np.nan
            if mean > 0:
                idx = np.arange(k - d)
                oe[idx, idx + d] = diag / mean
        i, j = np.nonzero(np.isfinite(oe))
        keep = (gi[i] >= 0) & (gi[j] >= 0)
        np.add.at(total, (gi[i[keep]], gi[j[keep]]), oe[i[keep], j[keep]])
        np.add.at(n, (gi[i[keep]], gi[j[keep]]), 1)
    s = total / np.maximum(n, 1)
    return (s + s.T) / 2


def compute(coords_um: np.ndarray, radius_um: np.ndarray, copy_id: np.ndarray,
            labels: tuple[str, ...], start: np.ndarray, end: np.ndarray, lad: np.ndarray,
            nuclear_radius: float, *, cutoff: float = 1.5) -> Compartments:
    """Compartiments de l'ensemble, et les trois mesures qui disent ce qu'ils valent."""
    bin_of, bins, xi = haploid_bins(labels, copy_id, start, end)
    n_bins = len(bins)

    # Piste LAD et rayon moyen par bin haploïde (moyenne sur les copies, Xi exclu du rayon :
    # c'est de l'Xa qu'on attend un comportement de compartiment).
    def per_bin(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        s = np.bincount(bin_of[mask], weights=values[mask], minlength=n_bins)
        c = np.bincount(bin_of[mask], minlength=n_bins)
        return s / np.maximum(c, 1)

    bins["lad"] = per_bin(lad.astype(np.float64), np.ones(len(lad), bool))
    radial = (np.linalg.norm(coords_um, axis=2) / nuclear_radius).mean(axis=0)
    bins["radial"] = per_bin(radial, ~xi)
    phasing = bins[["chrom", "start", "end"]].assign(anti_lad=1.0 - bins["lad"])

    half = len(coords_um) // 2
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "a").mkdir()
        (tmp / "b").mkdir()
        (tmp / "all").mkdir()
        pix = count_contacts(coords_um, radius_um, bin_of, n_bins, cutoff)
        e1, clr = _call(bins[["chrom", "start", "end"]], pix, phasing, tmp / "all")
        e1a, clr_a = _call(bins[["chrom", "start", "end"]],
                           count_contacts(coords_um[:half], radius_um, bin_of, n_bins, cutoff),
                           phasing, tmp / "a")
        e1b, clr_b = _call(bins[["chrom", "start", "end"]],
                           count_contacts(coords_um[half:], radius_um, bin_of, n_bins, cutoff),
                           phasing, tmp / "b")
        # Le saddle mesuré sur la carte même dont E1 est tiré est biaisé vers le haut :
        # le vecteur propre est, par construction, la partition qui fait le meilleur
        # damier de **cette** carte, bruit compris. Sur un ensemble sans aucun
        # compartiment, le test jouet le mesure à 2,8 au lieu de 1. On croise donc : E1
        # d'une moitié des structures, carte de l'autre, dans les deux sens.
        in_sample = _saddle(clr, e1)
        saddle = (_saddle(clr_b, e1a) + _saddle(clr_a, e1b)) / 2

    bins["E1"] = e1
    ok = np.isfinite(e1)
    both = ok & np.isfinite(e1a) & np.isfinite(e1b)
    per_bead = e1[bin_of].copy()
    per_bead[xi] = np.nan
    k = saddle.shape[0] - 1
    return Compartments(
        e1=per_bead,
        bins=bins,
        strength=float((saddle[0, 0] + saddle[k, k]) / (saddle[0, k] + saddle[k, 0])),
        strength_in_sample=float((in_sample[0, 0] + in_sample[k, k])
                                 / (in_sample[0, k] + in_sample[k, 0])),
        saddle=saddle,
        split_half_r=float(np.corrcoef(e1a[both], e1b[both])[0, 1]),
        r_lad=float(np.corrcoef(e1[ok], bins["lad"].to_numpy()[ok])[0, 1]),
        r_radial=float(np.corrcoef(e1[ok], bins["radial"].to_numpy()[ok])[0, 1]),
        fraction_a=float((e1[ok] > 0).mean()),
        called=int(ok.sum()),
        cutoff=cutoff,
        n_structures=len(coords_um),
    )
