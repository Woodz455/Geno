# Geno

Un modèle 3D navigable du génome humain, construit sur des données de conformation réelles,
où **chaque point cliqué répond de lui-même** : où suis-je sur le génome, dans quelle
structure nucléaire, quel gène, quelle régulation, quelle pathologie.

Deux mètres d'ADN diploïde repliés dans une sphère de dix micromètres — et il faut pouvoir
descendre jusqu'à la base.

## État

**Phase de cadrage.** Le plan est écrit, le code ne l'est pas encore.

## Ce qu'on construit

Un zoom continu sur quatre ordres de grandeur, du noyau (10 µm) à la paire de bases (2 nm),
en sept niveaux : noyau et repères nucléaires, territoires chromosomiques, compartiments A/B,
TADs, boucles CTCF, fibre de chromatine, double hélice.

Les niveaux macro viennent de **mesures** (Hi-C, Micro-C, DamID, TSA-seq). L'échelle fine
vient d'une **simulation contrainte par les mesures**. Le niveau séquence est de la
**géométrie déterministe**. L'interface dit toujours laquelle des trois on regarde.

Deux modes servis par un seul moteur :

- **Recherche** — import de ses propres `.mcool`/`.hic`, paramètres de reconstruction,
  ensembles, statistiques, export publication, sessions reproductibles.
- **Découverte** — parcours guidé, échelles comparées, glossaire, haltes commentées,
  accessibilité, FR/EN.

## Documentation

| Document | Contenu |
|----------|---------|
| [`docs/ROADMAP.md`](docs/ROADMAP.md) | Le plan sur 24 semaines, semaine par semaine, avec critères de fin |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Structure du dépôt, format `.g3d`, budget de rendu, stratégie deux-modes |
| [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md) | Catalogue des sources, accessions, licences, vérification |

## Principes

1. Une structure Hi-C unique est un artefact statistique — on affiche des ensembles.
2. Mesuré ≠ simulé ≠ interpolé, et ça se voit à l'écran.
3. Chaque bille remonte à une accession, une empreinte et une version de pipeline.
4. Validation contre données orthogonales, ou ce n'est qu'une illustration.
5. Coordonnées doubles hg38 / T2T-CHM13v2.0.

## Pile

Python (`cooler`, `cooltools`, `polychrom`/OpenMM, Snakemake) pour le pipeline.
TypeScript et Three.js/WebGL2 pour le moteur. `igv.js` et HiGlass pour les vues liées.
