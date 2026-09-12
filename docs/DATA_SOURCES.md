# Geno — Catalogue des sources de données

Toute donnée entrant dans le projet figure ici avec son accession, sa licence et son
empreinte. Livrable de la semaine 1 ; les accessions marquées *à figer en S1* doivent être
remplacées par l'identifiant exact avant tout téléchargement automatisé.

**Règle** : rien n'entre dans le dépôt ou dans un livrable public sans une ligne de licence
explicite dans ce tableau.

> **Décision arrêtée le 2026-09-12** — **GM12878** est la lignée de référence du projet.
> Référentiel primaire hg38, secondaire T2T-CHM13v2.0.

Le manifeste exécutable vit dans [`pipeline/data/manifest.tsv`](../pipeline/data/manifest.tsv) ;
ce document en est la justification scientifique.

---

## 1. Référentiels

| Jeu | Accession | Apport | Licence |
|-----|-----------|--------|---------|
| GRCh38 / hg38 | `GCA_000001405.15` | Référentiel primaire, coordonnées de travail | Libre |
| T2T-CHM13v2.0 | `GCA_009914755.4` | Centromères, ADNr, bras acrocentriques — absents de hg38 | Domaine public |
| Chaînes de conversion hg38 ↔ CHM13 | UCSC goldenPath | Double affichage des coordonnées | Libre |

CHM13 est une môle hydatiforme : essentiellement haploïde et homozygote. Le chrY de la v2.0
provient de HG002. À garder en tête quand on parle d'haplotypes.

## 2. Conformation (3D)

| Jeu | Accession | Apport | Licence |
|-----|-----------|--------|---------|
| Hi-C in situ GM12878 | GEO `GSE63525` (Rao et al. 2014, *Cell* 159:1665) | **Carte de référence.** La plus profonde publiée ; compartiments, TADs, boucles | Libre, citer |
| Micro-C H1-hESC | 4DN Data Portal — *à figer en S1* (Krietenstein et al. 2020, *Mol Cell*) | Échelle fine, sous le TAD | 4DN — libre après publication |
| Micro-C HFFc6 | 4DN Data Portal — *à figer en S1* (Akgol Oksuz et al. 2021, *Nat Methods*) | Second type cellulaire, contrôle | 4DN |
| SPRITE | Quinodoz et al. 2018, *Cell* | Contacts d'ordre supérieur — **validation orthogonale** | Libre, citer |
| GAM | Beagrie et al. 2017, *Nature* | Contacts sans ligature — validation orthogonale | Libre, citer |

GM12878 est retenue comme lignée de référence : caryotype quasi normal, profondeur de
séquençage inégalée, et lignée féminine — ce qui donne le X inactif et son superdomaine
comme démonstration native.

## 3. Repères nucléaires

| Jeu | Accession | Apport | Licence |
|-----|-----------|--------|---------|
| LADs (DamID lamin B1) | Guelen et al. 2008, *Nature* ; 4DN — *à figer en S1* | Contrainte de périphérie, indispensable à la reconstruction | Libre, citer |
| TSA-seq SON (speckles) | Chen et al. 2018, *J Cell Biol* ; 4DN | Distance aux speckles nucléaires | Libre, citer |
| NADs (nucléole) | *à figer en S1* | Contrainte périnucléolaire | À vérifier |
| Repli-seq | ENCODE / 4DN — *à figer en S1* | Timing de réplication, corrèle avec A/B | ENCODE |

## 4. Annotation 1D

| Jeu | Source | Apport | Licence |
|-----|--------|--------|---------|
| GENCODE | gencodegenes.org (version figée en S1) | Gènes, transcrits, biotypes | Libre |
| cCREs / SCREEN | ENCODE | Éléments cis-régulateurs candidats | ENCODE Data Use Policy — libre, citer |
| ChromHMM 18 états | Roadmap / ENCODE | État chromatinien | Libre, citer |
| CTCF ChIP-seq GM12878 | ENCODE — *à figer en S1* | Sites CTCF **avec orientation du motif** | ENCODE |
| Cohésine (RAD21, SMC3) | ENCODE — *à figer en S1* | Ancres de boucles | ENCODE |
| Cytobandes, gaps, centromères | UCSC goldenPath | Repères caryotypiques | Libre |

L'orientation du motif CTCF n'est pas un détail : la règle des motifs convergents est ce qui
rend la simulation d'extrusion de boucles prédictive plutôt que décorative.

## 5. Fonction et clinique

| Jeu | Source | Apport | Licence |
|-----|--------|--------|---------|
| ClinVar | NCBI | Variants et significations cliniques | Domaine public |
| GWAS Catalog | EMBL-EBI | Associations génotype–phénotype | Libre, citer |
| gnomAD v4 | Broad Institute | Contrainte : pLI, LOEUF, fréquences alléliques | Libre |
| GTEx | GTEx Portal | Expression par tissu | Accès libre aux données résumées |
| MONDO | OBO Foundry | Ontologie de maladies | CC-BY 4.0 |
| HPO | Monarch / Jax | Phénotypes | Libre, avec licence |
| Orphanet | Orphanet | Maladies rares | CC-BY 4.0 |

### Attention licence

| Ressource | Statut | Ce qu'on fait |
|-----------|--------|---------------|
| **OMIM** | **Non redistribuable.** Le téléchargement en masse exige une licence. | On n'embarque **aucun contenu OMIM**. On affiche le numéro MIM et un lien sortant. La description de maladie vient de MONDO / Orphanet. |
| **DECIPHER** | Conditions d'usage spécifiques, certaines données exigent un enregistrement | Hors périmètre v1. À réévaluer si besoin, avec lecture des conditions. |
| **4DN** | Libre après publication, mais certains jeux sont sous embargo | Vérifier le statut de chaque jeu au moment du téléchargement, pas une fois pour toutes. |

---

## 6. Vérification

Chaque entrée du catalogue produit une ligne dans `pipeline/data/manifest.tsv` :

```
nom	url	sha256	octets	licence	date_acces	version_pipeline
```

`make data` télécharge, calcule les sha256 et échoue si l'un d'eux diverge du manifeste. Une
structure 3D qui ne peut pas remonter à un manifeste vérifié n'est pas publiable.

Trois garanties, vérifiées par `make selftest` (13 assertions, sans réseau, via des URL
`file://`) :

- une source altérée en amont est **rejetée**, et le fichier corrompu n'est pas conservé ;
- une accession non résolue fait **échouer** le fetch au lieu d'être devinée ;
- le verrouillage d'une empreinte est une étape **explicite et séparée** — un téléchargement
  ne réécrit jamais un sha256 de lui-même.

---

## 7. Ce qu'on accepte comme source

### Sources primaires

Le dépôt d'origine, celui qui porte l'accession citée dans la publication : **4DN**,
**GEO/SRA**, **ENCODE**, **UCSC**, **EBI**, **NCBI**. Rien d'autre.

### Ce qu'on refuse comme source primaire

Toute copie ré-hébergée dont on ne peut pas remonter la chaîne de traitement — **HuggingFace
Hub**, Kaggle, un Zenodo tiers, un miroir de laboratoire.

La raison est technique, pas dogmatique. Pour du Hi-C, les choix de traitement — aligneur,
filtrage des duplicats et des self-ligations, correction ICE, taille de bin — **changent les
appels de compartiments et de TADs**. Un `.mcool` ré-hébergé sans ses paramètres est un
fichier dont on ne sait pas ce qu'il mesure. Le principe n° 3 (chaque bille remonte à une
accession et une version de pipeline) devient alors invérifiable, et tout ce que le modèle
affirme devient incontrôlable.

Une copie ré-hébergée reste utile pour **dégrossir** — regarder à quoi ressemble un jeu avant
de lancer un téléchargement de 40 Go. Elle ne produit jamais une position dans le modèle livré.

### Recoupements admis

**3D Genome Browser** (Yue lab), **3DIV**, **Nucleome Browser**, tilesets publics HiGlass :
excellents pour comparer nos appels à ceux d'autres équipes, jamais comme source de positions.

### Où HuggingFace a une vraie valeur — en v2, et pour des modèles

Pas pour les données. Pour les **modèles séquence → contact**, qui prédisent une carte de
contacts à partir de la séquence : **Akita** (Fudenberg 2020), **Orca** (Zhou 2022),
**C.Origami** (Tan 2023) ; et plus loin les modèles de séquence génomique (Enformer, Borzoi,
Nucleotide Transformer, DNABERT-2, HyenaDNA, Evo).

Pourquoi ça nous concerne précisément : la **semaine 15** prévoit un mode variant structural
qui, tel qu'il est écrit, applique une **règle** — on supprime une frontière de TAD, on montre
le détournement d'enhancer attendu. Un modèle séquence → contact le rendrait **prédictif** :
on entre le variant, il sort la carte de contacts prédite, on reconstruit dessus. C'est le
seul endroit du projet où ces modèles apportent quelque chose que les données seules ne
donnent pas.

**Contrainte** : une position issue d'un modèle prédictif ne porte jamais l'étiquette
« mesuré ». Elle relève d'un quatrième statut, `predicted`, distinct de `measured`,
`simulated` et `deterministic` — voir le champ `evidence` dans
[`ARCHITECTURE.md`](ARCHITECTURE.md).

---

## 8. État au 2026-09-12

`make status` : **28 entrées** — 15 `unverified`, 13 `unresolved`.

Aucune n'est encore `ok`, et ce n'est pas un oubli. L'environnement d'exécution de cette
session applique une politique d'egress fermée : `hgdownload.soe.ucsc.edu`,
`ftp.ebi.ac.uk`, `data.4dnucleome.org` et `ftp.ncbi.nlm.nih.gov` répondent tous **403** au
tunnel. `make data-core` échoue donc correctement sur les six entrées, sans rien laisser
derrière lui.

Deux façons de débloquer :

1. autoriser ces quatre hôtes dans la politique réseau de l'environnement ;
2. exécuter `make data && make data-lock` sur une machine qui a l'accès, et rapatrier le
   manifeste verrouillé — c'est précisément pour ça que les empreintes sont versionnées et
   pas les données.

Les 13 entrées `unresolved` attendent une accession exacte. Elles ne sont pas inventées : une
accession plausible mais fausse empoisonnerait la provenance plus sûrement qu'une case vide.
