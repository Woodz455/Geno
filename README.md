# Geno

Un modèle 3D navigable du génome humain, construit sur des données de conformation réelles,
où **chaque point cliqué répond de lui-même** : où suis-je sur le génome, dans quelle
structure nucléaire, quel gène, quelle régulation, quelle pathologie.

Deux mètres d'ADN diploïde repliés dans une sphère de dix micromètres — et il faut pouvoir
descendre jusqu'à la base.

## État

**Phase 0, semaine 2.** Lignée de référence figée (GM12878), machinerie de vérification des
données en place, socle 1D interrogeable.

```console
$ make build-store                              # construit depuis les fixtures
$ make query Q=chr7:5,527,000-5,530,600

chr7:5,527,000-5,530,600   3,601 pb   GRCh38   [fixtures]

  bande cytogénétique  (1)
    7p22.1  gneg                          chr7:5,000,001-7,200,000  2,200,000 pb

  gènes  (1)
    ACTB  -  protein_coding               chr7:5,527,151-5,530,601      3,451 pb
    ENSG00000075624.17

  exons  (6)
    …
```

Le magasin d'intervalles répond en **p99 0,82 ms** à l'échelle d'un locus sur un million
d'intervalles ([mesures](docs/ARCHITECTURE.md#8-socle-1d--index-dintervalles-et-mesures)).

> Les données affichées viennent de **fixtures**, pas de sources scientifiques — l'accès
> réseau aux dépôts est fermé dans l'environnement de développement. Toute fiche issue des
> fixtures le dit en en-tête. Voir [`docs/DATA_SOURCES.md` § 8](docs/DATA_SOURCES.md).

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

## Commandes

```
make help            toutes les cibles
make data            télécharge et vérifie les jeux de données du manifeste
make selftest        prouve la vérification d'empreintes, sans réseau
make build-store     construit le magasin d'intervalles
make query Q=…       interroge une région
make bench           latence de requête à l'échelle réelle
make test            suite de tests
```
