# Semaine 21 — critères de validation, fixés avant toute donnée

Ce document est commité **avant** le téléchargement des données de comparaison et avant toute
mesure. Une modification ultérieure sera datée ici et dite dans `VALIDATION.md` § S21, jamais faite
en silence.

## Ce qu'on valide

Le modèle **tel quel** (décision de l'utilisateur, 2026-11-24) :
- le noyau diploïde de 8 082 billes de 750 kb ;
- son ensemble de 200 repliements (`pipeline/data/ensemble/gm12878.zarr`).

Ce noyau n'a été construit à partir d'**aucune donnée de conformation mesurée**. Ses contraintes
sont le volume exclu, les liaisons, le confinement et un rappel radial ordonné par une **piste LAD
synthétique** (`nucleus/beads.py`). On s'attend donc à des échecs au-delà de ce que la distance le
long du fil explique. Ce rapport les chiffre.

Les régions fines (2 kb, ACTB et TP53) ne sont **pas** validées : aucune des données ci-dessous ne
les mesure à cette échelle dans GM12878, et leurs CTCF sont plantés.

## Les données

| Jeu | Lignée | Assemblage | Source |
|---|---|---|---|
| SPRITE, matrices intra et inter à 1 Mb (Quinodoz et al. 2018) | **GM12878** | hg19, levé en hg38 | GEO GSE114242 |
| Traçage multiplexé, chr21 à 50 kb, et 1 041 loci sur le génome (Su et al. 2020) | **IMR90**, pas GM12878 | hg38 | Zenodo 10.5281/zenodo.3928890 |
| pA-DamID LMNB1 (4DN, 4DNFINXDTU5Z) | **K562**, pas GM12878 | hg38 | déjà verrouillé (S13) |
| FISH de territoires en lymphoblastes (Croft et al. 1999 ; Boyle et al. 2001) | lymphoblastes | — | **texte** des articles |

Écartés, avec leur raison :
- **GAM** : pas de jeu humain en lymphoblastes ; Beagrie 2017 est murin ;
- **Bintu 2018** : 2,5 Mb de chr21, soit 3 billes du modèle, ce qui ne laisse rien à corréler.

## Ce que le modèle prédit, et comment on le lit

- **Proximité** de deux régions : la **distance médiane** entre leurs billes, sur les 200
  repliements et sur toutes les paires de copies (les deux homologues confondus, comme dans des
  données non phasées). C'est la sortie native du modèle, dense à toute séparation.
- **Contact** : la définition de la S7, des centres à moins de 1,5 × (r_i + r_j)
  (`ensemble/stats.contact_pairs`).
- **Position radiale** : la distance au centre du noyau sur le rayon du noyau (5 µm), moyennée sur
  les repliements et les copies.
- **Affectation** : une bille va au bin de 1 Mb qui contient son milieu ; un locus de traçage va à
  la bille qui le contient.

## Les trois références

- **N1, la distance seule** : un prédicteur qui ne connaît que la séparation génomique, la loi
  moyenne du modèle à cette séparation. Une corrélation brute du modèle ne compte que si elle
  **bat** N1. Sinon, ce qu'elle mesure est la distance génomique, pas la 3D.
- **N2, les permutations** : 1 000 permutations des identités des billes **à l'intérieur de chaque
  chromosome**, graine fixée. « Accord » exige ρ au-dessus du 95ᵉ centile des permutations.
- **Le témoin positif** : la première moitié des repliements contre la seconde, mesurée comme des
  données. Il doit donner ρ ≥ 0,9 ; sinon l'outil est cassé et le rapport invalide.

## Les critères

ρ est le coefficient de Spearman. « Corrigé de la distance » veut dire que chaque valeur est
divisée par la moyenne de sa diagonale (observé/attendu), des deux côtés.

| # | Données | Mesure | Accord si |
|---|---|---|---|
| C1a | SPRITE GM12878, intra | ρ brut entre proximité du modèle et contacts SPRITE, toutes paires intra à 1 Mb | ρ(modèle) > ρ(N1) |
| C1b | SPRITE GM12878, intra | ρ corrigé de la distance | ρ ≥ 0,3 et > permutations p95 |
| C1c | SPRITE GM12878, intra | compartiments : signe de E1 (orienté par la densité de gènes GENCODE), modèle contre SPRITE | ≥ 70 % des bins d'accord |
| C1d | SPRITE GM12878, intra | pente de P(s) entre 1 et 10 Mb, contacts du modèle contre SPRITE | écart ≤ 0,2 |
| C2 | SPRITE GM12878, inter | enrichissement par paire de chromosomes (observé / attendu par les marges, 253 paires), contacts du modèle contre SPRITE | ρ ≥ 0,3 et > permutations p95 |
| C3a | Su 2020 chr21 (IMR90) | exposant ν de R(s) ∝ s^ν entre 0,75 et 10 Mb | écart ≤ 0,1 |
| C3b | Su 2020 chr21 (IMR90) | ρ corrigé de la distance, distances médianes par paire de billes | ρ ≥ 0,3 et > permutations p95 |
| C3c | Su 2020 chr21 (IMR90) | échelle absolue : distance médiane modèle / traçage à 1 Mb et à 5 Mb (±20 %) | les deux rapports entre 1/1,5 et 1,5 |
| C4a | Su 2020, génome (IMR90) | ρ entre positions radiales, par locus | ρ ≥ 0,3 et > permutations p95 |
| C4b | Su 2020, génome (IMR90) | ρ entre positions radiales moyennes, par chromosome (23 valeurs) | ρ ≥ 0,5 |
| C4c | Su 2020, génome (IMR90) | ρ entre distances médianes de loci de chromosomes différents | ρ ≥ 0,3 et > permutations p95 |
| C5 | pA-DamID LMNB1 (K562) | ρ par bille entre position radiale et signal lamina moyen | ρ ≤ −0,3 et au-delà des permutations |
| C6 | FISH, lymphoblastes (texte) | chr19 plus central que chr18, en position radiale moyenne | vrai |

Précisions pour le traçage :
- **Position radiale d'un locus** : sa distance au centroïde de tous les loci détectés de la
  cellule, divisée par le 95ᵉ centile de ces distances dans la cellule, puis la médiane sur les
  cellules. Les noyaux d'IMR90 sont aplatis, ce qui sera dit.
- **Distances** : la médiane sur les cellules de toutes les paires de spots détectés, copies
  confondues.

**Verdicts**
- « **Accord** » : le seuil est franchi.
- « **Partiel** » : du bon signe et au-delà des permutations, mais sous le seuil.
- « **Échec** » : tout le reste.
- Une lignée étrangère est **étiquetée** à chaque ligne, sans seuil adouci.

## Ce qui fait la semaine

« Fait quand le rapport existe, avec ses corrélations et ses échecs nommés » (feuille de route).
La semaine est faite quand :
- les 13 critères sont mesurés ;
- les témoins tiennent : témoin positif ≥ 0,9, permutations centrées sur 0 ;
- chaque échec est dit avec sa cause probable.

Elle ne demande pas que le modèle réussisse.
