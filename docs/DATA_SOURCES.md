# Geno — Catalogue des sources de données

Toute donnée entrant dans le projet figure ici avec son accession, sa licence et son
empreinte. Livrable de la semaine 1 ; les accessions marquées *à figer en S1* doivent être
remplacées par l'identifiant exact avant tout téléchargement automatisé.

**Règle** : rien n'entre dans le dépôt ou dans un livrable public sans une ligne de licence
explicite dans ce tableau.

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
