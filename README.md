# Geno

Un modèle 3D navigable du génome humain, construit sur des données de conformation réelles,
où **chaque point cliqué répond de lui-même** : où suis-je sur le génome, dans quelle
structure nucléaire, quel gène, quelle régulation, quelle pathologie.

Deux mètres d'ADN diploïde repliés dans une sphère de dix micromètres — et il faut pouvoir
descendre jusqu'à la base.

## État

**Phase 0–2, semaine 7.** Lignée de référence figée (GM12878), machinerie de vérification des
données en place, socle 1D interrogeable, callers de conformation validés contre une structure
plantée, reconstruction 3D calibrée contre une géométrie connue, moteur de rendu par
imposteurs dont la correction est prouvée dans un vrai navigateur, **un noyau diploïde
entier** — 46 chaînes, 8 082 billes d'échelle TAD, sans interpénétration — et un **ensemble de
200 repliements** du même génome.

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

Côté conformation, `make hic-validate` plante une structure Hi-C connue — compartiments,
TADs, boucles, biais de couverture — et vérifie que les callers la retrouvent : équilibrage
ICE r = +0,98 contre le biais planté, compartiments 100 %, frontières de TAD 100 / 100,
boucles 89 / 94. C'est la seule configuration où « le caller est correct » est vérifiable :
sur des données réelles, un désaccord avec les appels publiés ne dit pas lequel des deux a
tort. Détails et limites dans [`docs/VALIDATION.md`](docs/VALIDATION.md).

`make recon` va plus loin : il fabrique une conformation 3D, en dérive les contacts par un
modèle direct d'exposant connu, et mesure quel exposant de reconstruction la restitue. Résultat
consigné — **l'exposant optimal n'est pas 1/gamma** : il décroît vers lui avec la profondeur de
séquençage sans jamais l'atteindre, et son minimum s'aplatit justement là où les données sont
les plus creuses. Reprendre `alpha = 1/3` d'un article sans regarder sa profondeur n'est donc
pas une convention.

`make nucleus` monte d'un cran : les 46 chaînes du caryotype 46,XX de GM12878, découpées en
billes d'échelle TAD, placées dans une sphère de dix micromètres sous volume exclu, chaînes
tenues, confinement et rappel des LADs vers la périphérie. Le résultat de la semaine tient en
une distinction : **la part du génome qui *est* en LAD et la part qui *touche* la lamina ne
sont pas la même quantité**, et la seconde est bornée par la géométrie du modèle — 5,1 % à
densité uniforme, 38 % en empilement maximal à 750 kb par bille ; 1,2 % et 9 % à 10 kb. Les
deux bornes sont linéaires en rayon de bille, donc une « fraction de LADs à la lamina » ne se
compare pas d'un modèle à l'autre sans dire à quelle résolution. Le noyau produit tient le
critère de la semaine — chevauchement maximal sous 0,8 % du contact sur six graines, aucune
bille hors du noyau — et ce qu'il ne prouve pas est dit aussi clairement que ce qu'il prouve
([`docs/VALIDATION.md` § S6](docs/VALIDATION.md)).

`make ensemble` applique le principe n° 1 : une structure unique est un artefact statistique,
et voici de combien. Deux cents repliements du même génome, et une décomposition de variance sur
la profondeur nucléaire de chaque bille — **ICC 0,776**, donc **22 % de la position radiale
d'une bille tient au tirage** et non à la bille. Deux structures indépendantes ne s'accordent
qu'à r = +0,78 sur qui est profond et qui est superficiel. Afficher une bille à sa profondeur
sans afficher ses 276 nm de dispersion, ce serait présenter un tirage comme une mesure.

L'ensemble prédit aussi ce qu'une expérience Hi-C verrait — sans qu'aucune matrice de contacts
n'ait jamais été montrée au modèle. Fraction trans 22,4 % là où le hasard donnerait 98 %, pas
d'appariement des homologues, et une **P(s) qui s'aplatit au-delà de 15 Mb** quand le Hi-C réel
continue de décroître. C'est un désaccord, et c'est le résultat le plus utile de la semaine : le
modèle a le *fait* des territoires, pas leur organisation interne. La semaine 8 est là pour ça.

Côté moteur, `pnpm verify` lance 16 assertions sur les imposteurs de sphères dans un vrai
Chromium : la profondeur bombe, l'interpénétration se résout par fragment et non par quad, et
le picking GPU reste exact parmi 250 005 instances sur toute la plage d'identifiants 32 bits.
Ce verdict porte sur la correction seulement — l'environnement ne rend que via SwiftShader,
un rasteriseur logiciel, donc les images par seconde attendent du vrai matériel.

> Les données affichées viennent de **fixtures**, pas de sources scientifiques — l'accès
> réseau aux dépôts est fermé dans l'environnement de développement. Toute fiche issue des
> fixtures le dit en en-tête. Voir [`docs/DATA_SOURCES.md` § 8](docs/DATA_SOURCES.md).

## Ce qu'on construit

Un zoom continu sur quatre ordres de grandeur, du noyau (10 µm) à la paire de bases (2 nm),
en sept niveaux : noyau et repères nucléaires, territoires chromosomiques, compartiments A/B,
TADs, boucles CTCF, fibre de chromatine, double hélice. Le modèle de noyau et ses lois
d'échelle sont dans [`ARCHITECTURE.md` § 10](docs/ARCHITECTURE.md#10-modèle-de-noyau--billes-tad-territoires-périphérie).

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
| [`docs/VALIDATION.md`](docs/VALIDATION.md) | Ce qui est mesuré et prouvé, semaine par semaine |
| [`docs/SETUP.md`](docs/SETUP.md) | Installation, et pourquoi deux niveaux de dépendances |

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
make hic-validate    plante une structure Hi-C connue et valide les callers
make recon           balaie l'exposant contact → distance contre une géométrie connue
make nucleus         construit un noyau diploïde complet de billes TAD
make ensemble        200 repliements du même génome, et ce qui s'y reproduit
make test            suite de tests
```

Côté rendu (Node 22, pnpm) :

```
pnpm verify          16 assertions sur les imposteurs, dans un vrai Chromium
pnpm typecheck       TypeScript strict
pnpm build           produit dist/spike.html — la page de mesure à ouvrir sur du matériel
```
