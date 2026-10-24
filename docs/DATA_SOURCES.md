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

## 8. État au 2026-10-19 : les sources de la fiche (semaine 13)

`make status` : **38 entrées** — **33 `ok`**, 1 `unverified`, 4 `unresolved`. La semaine 13 a résolu, téléchargé et verrouillé les 24 sources de la fiche
d'information. Le md5 de chaque fichier a été confronté à celui que **son dépôt publie**, partout
où il en publie un ; ailleurs, l'empreinte inscrite est la nôtre, prise au premier
téléchargement, et c'est dit.

| Entrée | Fichier | Contrôle indépendant |
|---|---|---|
| `encode.ccres` | ENCODE ENCSR820WFY / ENCFF590IMH — cCRE SCREEN en GM12878, GRCh38 | md5 = ENCODE |
| `chromhmm.gm12878` | ENCSR988QYW / ENCFF338RIC — ChromHMM **18 états**, GRCh38 | md5 = ENCODE |
| `ctcf.peaks.gm12878` | ENCSR000DZN / ENCFF960ZGP — pics IDR optimaux | md5 = ENCODE |
| `ctcf.motif.orientation` | JASPAR MA0139.1 (jaspar.genereg.net ; le miroir elixir coupe le tunnel) | — |
| `repliseq.gm12878` + `.S1` … `.G2` | UW, 6 phases, « percentage normalized signal », **hg19** | md5 = ENCODE (×6) |
| `rao2014.gm12878.subcompartments`, `.domains`, `.loops` | GEO GSE63525, **hg19** ; boucles : la liste HiCCUPS *avec motifs* | — |
| `chain.hg19.to.hg38` | Ensembl GRCh37 → GRCh38 | — |
| `lads.damid.lmnb1` | 4DNFINXDTU5Z — pA-DamID LMNB1, **K562**, 5 kb, moyenne des réplicats | md5 = 4DN |
| `tsaseq.son` | 4DNFIVZSO9RI — TSA-seq v2 SON, **K562**, 20 kb | md5 = 4DN |
| `nads.hffc6.nucleolus`, `.control` | 4DNFIW2JQD83, 4DNFI2S373O8 — NAD-seq **HFFc6**, RPKM | md5 = 4DN (×2) |
| `clinvar.vcf` | ClinVar GRCh38, version **20261019** | md5 = `.md5` publié par NCBI |
| `gwas.catalog` | GWAS Catalog, `associations_ontology-annotated-full.zip` du 2026-10-05 | — |
| `gnomad.constraint` | gnomAD v4.1, `constraint_metrics.tsv` | — |
| `gtex.median.tpm` | GTEx v10, TPM médian par tissu | — |
| `mondo.obo`, `hpo.obo` | ontologies MONDO et HPO | — |

**Trois remarques qui changent la lecture.**

- **Aucune mesure de lamina, de speckles ni de nucléole n'existe en GM12878** dans 4DN : DamID,
  pA-DamID et TSA-seq y couvrent K562, HCT116, H1, HFFc6, RPE — pas les lymphoblastes. Décision
  (question posée) : la fiche montre la mesure de la lignée la plus proche, **étiquetée** «
  mesuré dans K562, pas GM12878 » (NAD : HFFc6), à côté de la valeur du modèle, étiquetée «
  simulé ». Le timing de réplication, lui, est mesuré en GM12878.
- **Rao 2014 et le Repli-seq UW sont publiés en hg19.** Ils sont levés en hg38 par la chaîne
  Ensembl, et ce qui ne se lève pas est perdu et compté : 47 sous-compartiments sur 3 849, 38
  domaines sur 9 274, 22 boucles sur 9 448, 828 bins de 10 kb sur 211 913.
- **UCSC ne passe plus le TLS** depuis le 2026-10-19 (« unable to get local issuer certificate »
  — il répondait 200 la veille). Pas de contournement : tout ce dont la semaine 13 avait besoin
  était déjà local (cytobandes, chaînes, `.2bit`) ou existe ailleurs (chaîne GRCh37 → GRCh38
  chez Ensembl). Ce qui manque à cause de cela : la cytobande T2T — la fiche donne la position
  T2T et la bande hg38.

**La chaîne hg38 → T2T, contrôlée sans UCSC.** Pour 10 000 positions hg38 tirées au hasard, le
21-mer centré est comparé dans `hg38.2bit` et dans `hs1.2bit` à la position levée : 9 773
positions se lèvent, 9 558 de celles-là ont le même 21-mer (97,8 %) et 9 764 la même base
(99,91 %) — ce qu'on attend de deux assemblages qui diffèrent d'environ une SNV par kb. Témoin,
position levée décalée de 1 kb : 0 21-mer, 25,7 % de bases, le hasard. (`geno fiche check`, VALIDATION § S13.)

**Licences.** ClinVar : domaine public. GWAS Catalog : libre, citation demandée. gnomAD : libre.
GTEx : données résumées en accès ouvert. MONDO : CC-BY 4.0, cité dans la fiche. HPO : licence
propre, citation demandée. ENCODE et 4DN : politique d'usage de chaque consortium, citée par
accession. **OMIM** : ClinVar en cite des numéros ; la fiche les garde comme identifiants et
n'affiche aucun texte OMIM — les noms de maladies viennent de MONDO, de HPO, ou du libellé
ClinVar.

**Restent `unresolved`, hors de la fiche** : le mcool de 4DN, les deux Micro-C, RAD21. Le `.hic`
de Rao 2014 (plusieurs dizaines de Go) reste `unverified` : la fiche n'en a pas besoin, et
`fetch.sh` accepte désormais `ONLY=` pour ne pas le tirer avec le reste de son tier.

## 8 ter. Semaine 14 : aucune source nouvelle, deux produits dérivés

Les vues liées n'ajoutent rien au manifeste. Ce qu'igv.js et HiGlass lisent est **écrit** par
`make linked`, depuis ce qui est déjà verrouillé, et chaque fichier porte l'empreinte de ce dont il
vient (`data/g3d/igv/index.json`, champ `provenance`) :

| Fichier | Vient de | Évidence |
|---|---|---|
| `igv/genes.bed.gz` | `gencode.annotation.gtf` (v50) | annotation |
| `igv/ctcf.bed.gz` | `ctcf.peaks.gm12878`, brin = motif MA0139.1 scanné sur `hg38.2bit` (§ 8) | mesuré |
| `igv/loops.bedpe.gz` | `rao2014.gm12878.loops`, levées par `chain.hg19.to.hg38` (22 perdues sur 9 448) | mesuré |
| `igv/compartments.bedgraph.gz` | l'E1 du `.g3d` (carte de contacts du modèle) | simulé |
| `igv/cytoband.txt.gz`, `igv/chrom.sizes` | `hg38.cytoband`, `hg38.chrom.sizes` | référence |
| `gm12878.contacts` | l'ensemble de 200 repliements (semaine 7), les régions fines (semaine 8) | **simulé** |
| séquence d'igv | `hg38.2bit`, servi tel quel par plages | référence |

**Pas de Hi-C mesuré dans HiGlass.** Le `.hic` de Rao 2014 (`rao2014.gm12878.hic`) reste
`unverified` : plusieurs dizaines de Go, jamais téléchargé, donc jamais verrouillé ; le `.mcool` de
4DN est `unresolved`. Les montrer est le rattrapage S1–S8, pas la semaine 14.

## 8 quinquies. Semaine 16 : aucune source nouvelle, un index dérivé

`gm12878.search` est **écrit** par `make search` depuis ce qui est déjà verrouillé, et porte
l'empreinte de chaque source dans son en-tête (`provenance`) :

| Table | Vient de | Ce qui en est gardé |
|---|---|---|
| gènes | `gencode.annotation.gtf` (v50), toutes les lignes `gene` | nom, ENSG, chrom, début, fin, brin, biotype |
| rsID | `clinvar.vcf` (champ `RS`, toutes classes) | rsID, chrom, POS − 1 |
| rsID | `gwas.catalog` (SNPS = un seul rsID, position GRCh38 du catalogue) | rsID, chrom, CHR_POS − 1 |
| bornes | `hg38.chrom.sizes`, chromosomes primaires | validation des positions |

**Pas dbSNP**, par décision de l'utilisateur : un rsID sans variant ClinVar ni association GWAS
n'est pas dans l'index, et la page le dit. Ajouter dbSNP serait une source nouvelle (des dizaines
de Go) à verrouiller.

## 8 quater. Semaine 15 : une source lue, pas téléchargée

Aucune entrée nouvelle au manifeste. Le locus EPHA4 se construit depuis ce qui est déjà verrouillé
— GENCODE v50, les pics CTCF d'ENCODE GM12878 orientés par MA0139.1 (§ 8), les boucles de Rao 2014
levées en hg38 —, et les couches cliniques depuis le fichier de fiches (ClinVar, catalogue GWAS,
gnomAD v4.1 constraint, § 8).

**Lupiáñez et al., Cell 161:1012–1025 (2015)**, doi:10.1016/j.cell.2015.04.004 : le **texte**
seulement, lu par l'API BioC de PMC, pour contraindre la reconstruction (longueurs des délétions,
distance enhancers → promoteur, témoin DelBS). La **Table S1** (points de cassure des patients) n'a
pas été lue : PMC répond par un reCAPTCHA, Europe PMC par 403/500, les fichiers supplémentaires ne
sont pas exposés, et les délétions ne sont ni dans ClinVar ni dans dbVar. Les copies hébergées
ailleurs ne sont pas acceptées comme source (§ 7). Toute coordonnée qui en dépend est donc
étiquetée **reconstruite depuis le texte** dans `epha4.variant.json` (`evidence`), et les mesures
portent sur toute la famille de points de cassure que le texte autorise. Les CTCF sont ceux d'une
lignée **lymphoblastoïde**, pas d'un bourgeon de membre.

## 8 bis. État au 2026-10-18 : le réseau s'est ouvert

Jusqu'à la semaine 11, l'environnement d'exécution appliquait une politique d'egress fermée :
`hgdownload.soe.ucsc.edu`, `ftp.ebi.ac.uk`, `data.4dnucleome.org` et `ftp.ncbi.nlm.nih.gov`
répondaient **403** au tunnel, et tout ce qui suit dans ce document était écrit sans un octet
récupéré. Le 2026-10-18, au début de la semaine 12, ils répondent **200** — et ENCODE, 4DN, EBI,
NCBI, l'API REST d'Ensembl et celle d'UCSC aussi.

`make status` : **28 entrées** — 9 `ok`, 7 `unverified`, 12 `unresolved`.

| Entrée | Statut | Contrôle indépendant |
|---|---|---|
| `hg38.chrom.sizes`, `hg38.cytoband`, `hg38.gap`, `hs1.chrom.sizes`, chaînes hg38 ↔ hs1 | `ok` | — |
| `hg38.2bit` (835 Mo) | `ok` | md5 `dcc3ea27…` = `md5sum.txt` publié par UCSC |
| `hs1.2bit` (812 Mo) | `ok` | md5 `68202ff4…` = `md5sum.txt` publié par UCSC |
| `gencode.annotation.gtf` — **v50, sous-ensemble basic**, release figée | `ok` | md5 `11e77cf1…` = `MD5SUMS` publié par EBI |

Le téléchargement de `hg38.2bit` a été coupé deux fois par le tunnel (« connection reset »,
`ws_closed_mid_exchange` dans l'état du proxy), puis repris par plage (`curl -C -`) : c'est
l'empreinte, confrontée à celle du dépôt, qui dit que le fichier est entier, pas l'absence
d'erreur.

**Ce que la semaine 12 en fait, et seulement ça.** La séquence des deux régions fines est lue
dans `hg38.2bit` et écrite dans le `.g3d` ; la fenêtre autour de TP53 c.818 est identique,
base pour base, à celle que rend l'API UCSC (`AAACACGCACCT` pour chr17:7 673 797–7 673 808).
TP53 et ACTB sont situés par les CDS de leur transcrit MANE dans GENCODE v50 ; c.818 de TP53
retombe sur chr17:7 673 802, la position de rs28934576 (R273H) dans dbSNP via Ensembl.

**La table interne des longueurs de chromosomes est revérifiée** : ses 24 valeurs concordent à
la base près avec `hg38.chrom.sizes`. Le caryotype lisait le fichier à `data/core/`, où rien
n'arrive jamais — `make data-core` dépose dans `data/raw/core/` ; le chemin est corrigé, et les
modèles construits désormais portent la provenance du fichier au lieu de `builtin`, pour des
longueurs identiques. Les modèles déjà écrits (semaines 6 à 8) ne sont pas reconstruits : ils
n'en changeraient pas d'une base.

**Ce qui reste à faire, et qui est maintenant faisable.** Les 12 entrées `unresolved` attendent
toujours une accession exacte (4DN, ENCODE, JASPAR) — elles ne sont pas inventées ; les 7
`unverified` (Hi-C de Rao 2014, ClinVar, GWAS, MONDO) n'ont pas encore été récupérées. Les
positions du modèle, elles, restent toutes **simulées** : le noyau sur une piste LAD synthétique,
les régions fines sur des sites CTCF **plantés**. Remplacer cela — DamID publié, CTCF d'ENCODE
orientés par le motif sur la vraie séquence, Hi-C de Rao 2014 en validation — est le rattrapage
des semaines 1 à 8, décidé comme une étape à part, après la semaine 12.

### Ce qui tourne quand même, et à quel prix

Le modèle de noyau de la semaine 6 a besoin des longueurs de chromosomes. Une longueur de
chromosome est une propriété déterministe de l'assemblage, pas une mesure, mais elle reste une
donnée et la règle du projet est qu'une donnée remonte à un fichier empreint. D'où une
hiérarchie explicite :

1. si `data/raw/core/hg38.chrom.sizes` existe — donc récupéré et verrouillé par `make data-core` —
   c'est lui qui fait foi ;
2. sinon, une table interne à `geno_pipeline/nucleus/karyotype.py`, et **tout ce qui en découle
   est estampillé `builtin`**.

L'estampille n'est pas décorative : elle sort dans l'en-tête de `make nucleus` et dans le
`.json` frère de chaque structure écrite.

Même logique pour la piste LAD, qui n'a toujours aucune entrée résolue : le modèle en fabrique
une, synthétique, calibrée sur les statistiques publiées de taille et de couverture, et
l'estampille `lad_source: "synthetic"`. Ce qui s'en déduit est géométrique ; la corrélation avec
le DamID publié attend que l'accession soit figée.
