# Geno — Feuille de route 24 semaines

**Objectif** : un modèle 3D navigable du génome humain, construit sur des données de
conformation réelles, où **chaque point cliqué répond de lui-même** : où suis-je sur le
génome, dans quelle structure nucléaire, quel gène, quelle régulation, quelle pathologie.

Deux publics servis par un seul moteur : **Mode Recherche** (bio-informatique, données
propres, export publication) et **Mode Découverte** (grand public, parcours guidé).

---

## 1. Ce que nous construisons : les 7 échelles

Le cœur du projet est un **zoom continu sur 4 ordres de grandeur**, du noyau (10 µm) à la
paire de bases (2 nm). Dit autrement : 2 mètres d'ADN diploïde repliés dans une sphère de
10 µm, et il faut pouvoir descendre jusqu'à la base.

| # | Échelle | Objet biologique | Résolution génomique | Origine des positions |
|---|---------|------------------|---------------------|----------------------|
| 0 | ~10 µm | Noyau, lamina, nucléole, speckles | — | TSA-seq, DamID, imagerie |
| 1 | 1–2 µm | Territoires chromosomiques | chromosome entier | Hi-C 1 Mb + contraintes radiales |
| 2 | 0,5–1 µm | Compartiments A/B (A1, A2, B1–B4) | 100 kb – 1 Mb | Vecteur propre Hi-C, SNIPER/Calder |
| 3 | 100–500 nm | TADs | 10–100 kb | Insulation score, Arrowhead |
| 4 | 50–200 nm | Boucles CTCF, contacts enhancer–promoteur | 5–25 kb | HiCCUPS, Mustache, Micro-C |
| 5 | 10–30 nm | Fibre de chromatine, nucléosomes | 200 pb – 5 kb | Micro-C + simulation polymère |
| 6 | 2 nm | Double hélice B | 1 pb | Séquence de référence, géométrie procédurale |

Les niveaux 0–4 viennent de **mesures**. Le niveau 5 vient d'une **simulation contrainte par
les mesures**. Le niveau 6 est de la **géométrie déterministe**. L'interface doit dire
laquelle des trois on regarde, à tout moment. C'est non négociable (§4).

---

## 2. Les deux modes

Un seul noyau de code, deux coquilles. Pas deux applications.

**Mode Recherche**
- Import de ses propres `.mcool` / `.hic`
- Choix des paramètres de reconstruction, exploration de l'ensemble de structures
- Statistiques : position radiale, fréquence de contact, eigenvector de compartiment
- Export figures (SVG/PNG haute résolution) et structures (glTF, PDB, 3DG)
- Manifeste de session JSON → reproductibilité bit-à-bit

**Mode Découverte**
- Parcours narratif « du noyau à une base en 7 étapes »
- Échelles comparées, glossaire au survol, aucun jargon sans définition
- Haltes commentées : le locus HBB, l'enhancer FTO→IRX3, la frontière TAD de HOXD et les
  malformations des membres, le superdomaine Xist du X inactif
- Accessibilité : navigation clavier, palettes sûres pour daltoniens, `prefers-reduced-motion`
- FR / EN

Techniquement : même `packages/core`, un contexte `mode`, des feature flags. Le mode
Découverte est le mode Recherche avec un préréglage curaté et une couche narrative.

---

## 3. Pile technique

| Couche | Choix | Pourquoi |
|--------|-------|----------|
| Pipeline données | Python — `cooler`, `cooltools`, `bioframe`, `polychrom`/OpenMM | Standard de fait du domaine (consortium Open2C) |
| Reconstruction | Chrom3D-like (recuit simulé) + Pastis pour comparaison | Chrom3D intègre les contraintes de lamina et sort un noyau diploïde complet |
| Orchestration | Snakemake + conteneur | Chaque figure regénérable depuis les accessions brutes |
| Format 3D | `.g3d` binaire multi-résolution + index d'intervalles | Requêtes HTTP Range, « Google Maps du génome » |
| Rendu | Three.js / WebGL2, imposteurs de sphères, picking GPU | Millions de billes à 60 fps ; techniques éprouvées par Mol*/NGL |
| Vues liées | `igv.js` (piste linéaire) + HiGlass (matrice de contacts) | La 3D seule désoriente ; il faut les trois vues synchronisées |
| App | TypeScript, Vite, monorepo pnpm | — |

---

## 4. Principes non négociables

1. **Une structure Hi-C unique est un artefact statistique.** Le Hi-C moyenne des millions de
   cellules. On produit et on affiche un **ensemble** (N ≈ 200), avec la variabilité par bille
   visible. Jamais une structure présentée comme « le » génome.
2. **Mesuré ≠ simulé ≠ interpolé.** Un code couleur permanent distingue les trois.
3. **Provenance traçable.** Chaque bille remonte à un accession + un hash de données + une
   version de pipeline.
4. **Validation contre données orthogonales.** DNA-FISH, chromatin tracing (ORCA), DamID
   lamin B1, SPRITE/GAM. Un modèle non validé est une illustration, pas un résultat.
5. **Coordonnées doubles.** hg38 *et* T2T-CHM13v2.0. Les centromères, l'ADNr et les bras
   acrocentriques n'existent pas dans hg38 — on les marque explicitement comme modélisés
   par contrainte.

---

## 5. Les 6 phases

| Phase | Semaines | Titre | Fin de phase = |
|-------|----------|-------|----------------|
| 0 | 1–4 | Socle données & décisions | Toute donnée téléchargeable et vérifiée par `make data` |
| 1 | 5–8 | Reconstruction 3D | Un ensemble de noyaux diploïdes validé quantitativement |
| 2 | 9–12 | Moteur de rendu multi-échelle | Zoom continu noyau → paire de bases, sans couture |
| 3 | 13–16 | Couche d'annotation & vues liées | Clic sur un point → fiche complète, 3 vues synchronisées |
| 4 | 17–20 | Les deux modes | Recherche et Découverte livrés |
| 5 | 21–24 | Validation, perf, publication | Déployé, DOI, méthodes publiables |

---

## 6. Semaine par semaine

### Phase 0 — Socle données & décisions (S1–S4)

#### Semaine 1 — Cadrage scientifique et lignée de référence
- **Objectif** : figer le référentiel et la lignée cellulaire de référence.
- **Décisions** : hg38 primaire, T2T-CHM13v2.0 secondaire. **GM12878** comme lignée de
  référence (Rao 2014, GSE63525 — la carte Hi-C in situ la plus profonde publiée,
  caryotype quasi normal, femelle → offre le X inactif en démonstration).
  **H1-hESC / HFFc6 Micro-C** pour l'échelle fine.
- **Livrable** : `docs/DATA_SOURCES.md` complété — accession, taille, licence, checksum
  pour chaque jeu. Piège licence : OMIM n'est pas redistribuable → on passe par
  MONDO / HPO / Orphanet et on ne garde d'OMIM que les identifiants liés.
- **Fait quand** : `make data` télécharge tout et vérifie les sha256.

#### Semaine 2 — Socle de données 1D
- **Objectif** : pouvoir interroger n'importe quel intervalle génomique en moins d'une milliseconde.
- **Travaux** : GENCODE (GTF → binaire indexé), cytobandes, centromères, télomères, gaps ;
  cCREs ENCODE/SCREEN ; états ChromHMM ; pics CTCF **avec orientation du motif** (essentiel
  pour l'extrusion de boucles).
- **Livrable** : magasin d'intervalles unifié + CLI.
- **Fait quand** : `geno query chr7:5,527,000-5,530,600` renvoie le locus ACTB complet.

#### Semaine 3 — Socle de données 3D (Hi-C)
- **Objectif** : extraire les features de conformation des cartes brutes.
- **Travaux** : chargement `.mcool`, balancing ICE, décomposition en vecteurs propres →
  compartiments A/B (orientation fixée par densité génique), insulation score → frontières
  de TAD, HiCCUPS + Mustache → boucles. Résolutions 1 Mb / 250 kb / 100 kb / 25 kb / 10 kb / 5 kb / 1 kb.
- **Fait quand** : nos appels de compartiments et de boucles recoupent ceux publiés par
  Rao 2014 à plus de 80 %. Écart mesuré et expliqué, pas ignoré.

#### Semaine 4 — Spike technique & architecture
- **Objectif** : savoir ce que le navigateur encaisse *avant* d'écrire le moteur.
- **Travaux** : prototype Three.js — 1 M d'imposteurs de sphères, mesure du framerate et de
  la latence de picking GPU sur 3 cibles (GPU desktop, iGPU portable, téléphone).
- **Livrable** : budget de billes par niveau de LOD, écrit et chiffré. Squelette du monorepo.
  ADR sur la stratégie deux-modes.
- **Fait quand** : go/no-go signé sur le budget de rendu.

---

### Phase 1 — Reconstruction 3D (S5–S8)

#### Semaine 5 — Contact → distance
- **Objectif** : la conversion la plus critique du projet. `d ∝ f^(-α)`, on balaie α.
- **Travaux** : implémenter et comparer loi de puissance, approche Poisson (Pastis), MDS
  métrique (ShRec3D).
- **Fait quand** : le modèle reproduit des faits connus indépendants — chr19 (dense en gènes)
  plus central que chr18, taille des territoires cohérente avec la FISH publiée.

#### Semaine 6 — Noyau diploïde complet
- **Objectif** : première structure entière.
- **Travaux** : modèle bille-par-TAD, 46 chromosomes, contraintes LAD → périphérie, rayon
  nucléaire, volume exclu, recuit simulé.
- **Fait quand** : un noyau diploïde, ~6 000–10 000 billes TAD, sans interpénétration.
- **Fait.** `make nucleus` — 8 082 billes de 750 kb, chevauchement maximal sous 0,5 % du
  contact sur six graines, aucune bille hors du noyau. Résultat de fond : « en LAD » n'est pas
  « à la lamina », et l'écart est une loi d'échelle ([`VALIDATION.md` § S6](VALIDATION.md)).

#### Semaine 7 — Ensembles, pas une structure
- **Objectif** : appliquer le principe n°1 (§4).
- **Travaux** : N = 200 structures depuis des graines différentes ; variabilité par bille,
  structure médoïde, distributions radiales.
- **Livrable** : `ensemble.zarr` + rapport de validation chiffré (corrélation position radiale
  modélisée vs. LADs DamID publiés).
- **Fait quand** : le rapport existe avec ses nombres, réussite ou échec.

#### Semaine 8 — Échelle fine : polymère et extrusion de boucles
- **Objectif** : descendre sous le TAD.
- **Travaux** : sur une région choisie de 2–4 Mb, simulation polymère d'extrusion de boucles
  (polychrom/OpenMM) conditionnée par les sites CTCF orientés et le Micro-C. Résolution 1–5 kb.
- **Fait quand** : le segment fin se raccorde géométriquement au modèle de noyau entier.

---

### Phase 2 — Moteur de rendu multi-échelle (S9–S12)

#### Semaine 9 — Format et streaming
- **Objectif** : ne jamais charger ce qu'on ne regarde pas.
- **Travaux** : format `.g3d` — en-tête (assemblage, niveaux, type cellulaire, hash de
  provenance) + tableaux de positions par LOD + index d'intervalles. Découpage octree,
  requêtes HTTP Range.
- **Fait quand** : niveau noyau < 500 ko, premier rendu < 1,5 s en 4G.

#### Semaine 10 — Niveaux macro
- Noyau, territoires, compartiments. Imposteurs de sphères, occlusion ambiante écran-espace
  (sans elle, un nuage de billes est illisible). Caméra orbitale + transition « plongée ».
- **Fait quand** : on descend du noyau aux compartiments en plongée continue, à 60 fps sur la
  cible matérielle médiane.

#### Semaine 11 — Niveaux TAD, boucle, fibre
- Géométrie tubulaire depuis la polyligne, tessellation adaptative. Ancres de boucles rendues
  comme arcs de liaison. Nucléosomes à l'échelle 10 nm.
- **Note scientifique** : on représente la fibre de 10 nm en chaîne désordonnée, pas en fibre
  régulière de 30 nm — la ChromEMT en cellule (Ou et al. 2017) n'observe pas cette dernière
  *in situ*.
- **Fait quand** : une boucle CTCF se lit comme telle, ancres comprises, à l'échelle où elle
  existe.

#### Semaine 12 — Niveau séquence : la double hélice procédurale
- **Objectif** : la démonstration qui fait le projet.
- **Travaux** : génération de la géométrie B-ADN (10,5 pb/tour, pas de 3,6 nm, rise 0,34 nm,
  grand et petit sillon) depuis la séquence réelle de la fenêtre visible uniquement.
- **Fait quand** : un zoom continu enregistré, du noyau jusqu'à une base précise de TP53,
  sans coupure ni chargement visible.

---

### Phase 3 — Annotation & vues liées (S13–S16)

#### Semaine 13 — Picking GPU et fiche d'information
- **Objectif** : « des détails et informations sur les endroits pointés » — le cœur de la demande.
- **Travaux** : buffer d'identifiants entier (RGBA32UI), lecture d'un pixel au clic → O(1)
  quel que soit le nombre de billes. Mapping inverse position 3D → coordonnée génomique.
- **La fiche affiche** :
  - Coordonnées : chr, position hg38 **et** T2T, cytobande
  - Contexte structural : compartiment (A/B + sous-compartiment), TAD et distance à la
    frontière la plus proche, ancres de boucles, position radiale normalisée (0 = centre,
    1 = périphérie), statut LAD/NAD, distance aux speckles (TSA-seq), timing de réplication
  - Gènes : GENCODE chevauchants ou voisins, brin, biotype, expression GTEx
  - Régulation : cCREs, état ChromHMM, sites CTCF avec orientation
  - Clinique : variants ClinVar, maladies (MONDO/HPO), hits GWAS, contrainte gnomAD (pLI, LOEUF)
  - Séquence : au zoom maximal, les bases réelles de la référence
- **Fait quand** : un clic sur n'importe quelle bille produit la fiche complète en moins de 100 ms.

#### Semaine 14 — Vues liées
- `igv.js` et HiGlass montés à côté de la 3D. Une sélection unique publiée par un store —
  pas de synchronisation deux-à-deux, source connue des boucles infinies. Brossage :
  sélection dans le navigateur linéaire → surbrillance en 3D, et réciproquement.
- **Fait quand** : une sélection faite dans l'une des trois vues apparaît dans les deux
  autres, sans rebond.

#### Semaine 15 — Couche clinique et fonctionnelle
- Surimpression ClinVar / GWAS / contrainte gnomAD sur la structure.
- **Mode variant structural** : supprimer une frontière de TAD et voir le détournement
  d'enhancer prédit — la démonstration de Lupiáñez 2015 (EPHA4, malformations des membres)
  et Franke 2016 (néo-TAD, SOX9), rejouable en direct.
- **Fait quand** : la suppression d'une frontière rejoue en direct le cas EPHA4 des
  malformations des membres.

#### Semaine 16 — Recherche et navigation
- Recherche « BRCA1 », « chr17:43,044,295 », « rs334 » → vol de caméra vers la cible.
- Signets, URL partageables encodant caméra + sélection.
- **Fait quand** : les trois formes de requête amènent la caméra au bon endroit, et l'URL
  rejoue la vue à l'identique.

---

### Phase 4 — Les deux modes (S17–S20)

#### Semaines 17–18 — Mode Recherche
- Import `.mcool` / `.hic` de l'utilisateur, reconstruction déclenchable, panneau de paramètres,
  explorateur d'ensemble, statistiques, exports (PNG/SVG publication, glTF, PDB, 3DG),
  manifeste de session reproductible.
- **Contrainte mémoire** : un ensemble complet (200 × 620 000 billes) pèse ~1,5 Go → on
  streame le médoïde + un scalaire de variabilité par bille ; l'ensemble complet ne se charge
  qu'à la demande et en résolution grossière.
- **Fait quand** : un `.mcool` importé ressort en structure, et le manifeste de session la
  reproduit à l'identique.

#### Semaines 19–20 — Mode Découverte
- Parcours guidé, comparaisons d'échelles, glossaire, haltes commentées, accessibilité
  complète, i18n FR/EN.
- **Fait quand** : un lecteur sans formation en génomique termine le parcours et sait
  expliquer ce qu'est un TAD.

---

### Phase 5 — Validation, performance, publication (S21–S24)

#### Semaine 21 — Validation scientifique
- Confrontation aux données orthogonales : distances DNA-FISH publiées, chromatin tracing
  (Bintu 2018, Su 2020), DamID lamin B1, contacts SPRITE et GAM.
- **Livrable** : rapport de validation quantitatif. C'est ce qui sépare un outil d'une jolie image.
- **Fait quand** : le rapport existe, avec ses corrélations **et ses échecs nommés**.

#### Semaine 22 — Performance et compatibilité
- Budget de frame, mémoire, mobile, chemin de repli WebGPU, particularités Safari.
- Tests de non-régression de performance automatisés.
- **Fait quand** : les budgets tiennent sur les trois cibles, et un commit qui les dépasse
  casse la CI.

#### Semaine 23 — Documentation, reproductibilité, tests
- Pipeline Snakemake conteneurisé ; toute figure regénérable depuis les accessions brutes.
- Documentation d'API + document de méthodes énonçant **toutes** les limites sans les adoucir.
- **Fait quand** : toute figure du projet se regénère depuis les accessions brutes, par une
  seule commande.

#### Semaine 24 — Publication et livraison
- Déploiement, DOI Zenodo, section méthodes prête pour préprint, vidéo de démonstration,
  feuille de route v2.
- **Fait quand** : l'outil est en ligne, citable, et la v2 est écrite.

---

## 7. Jalons de démonstration

| Fin de | On peut montrer |
|--------|-----------------|
| S4 | Les données sont là, vérifiées ; on sait ce que le GPU encaisse |
| S8 | Un noyau diploïde qui tourne, calculé depuis du Hi-C réel |
| S12 | **Le zoom continu du noyau à la paire de bases** |
| S16 | Clic sur n'importe quel point → fiche complète, trois vues synchronisées |
| S20 | Les deux applications, utilisables par leurs publics respectifs |
| S24 | Déployé, validé, citable |

---

## 8. Risques et parades

| # | Risque | Parade |
|---|--------|--------|
| R1 | Le Hi-C moyenne les cellules → structure trompeuse | Ensembles + variabilité visible. Décidé S7, jamais renégocié |
| R2 | Plafond de résolution (le 1 kb exige des milliards de lectures) | Sous le TAD, c'est de la simulation contrainte, étiquetée comme telle |
| R3 | Centromères, ADNr, bras acrocentriques absents de hg38 | Bascule T2T-CHM13v2.0, régions marquées « modélisées par contrainte » |
| R4 | Mémoire et perf navigateur | Budget LOD figé en S4, tenu par des tests automatisés |
| R5 | Dérive de périmètre entre les deux modes | Un seul noyau + feature flags, ADR en S4 |
| R6 | Licences (OMIM non redistribuable) | MONDO / HPO / Orphanet ; d'OMIM, seulement les identifiants liés |
| R7 | Le zoom continu est le point dur technique | Attaqué en S12, après le spike S4 — pas en fin de projet |

---

## 9. Ce dont j'ai besoin de toi

- **S1** : valider GM12878 comme lignée de référence, ou en proposer une autre selon ton usage.
- **S4** : arbitrer le budget de rendu si les cibles basses ne suivent pas.
- **S8** : choisir la région de 2–4 Mb pour la simulation fine (une région qui t'intéresse
  scientifiquement vaut mieux qu'une région choisie au hasard).
- **S15** : la liste des loci cliniques prioritaires.
- **S19** : le ton du mode Découverte — musée, cours, ou documentaire.
