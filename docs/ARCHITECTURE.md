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

Chiffres à confirmer par le spike de semaine 4. Le budget est ensuite **tenu par des tests
de non-régression** : un commit qui fait dépasser le budget casse la CI.

---

## 5. Rendu

**Imposteurs de sphères, pas de maillages.** Un quad orienté caméra, le fragment shader
calcule la normale de sphère et la profondeur. C'est ce que font Mol* et NGL, et c'est ce qui
permet le million de billes à 60 fps là où `InstancedMesh` avec une vraie géométrie de sphère
s'effondre.

**Picking GPU.** Rendu des identifiants dans une cible `RGBA32UI` hors écran, lecture d'un
pixel au clic. Coût O(1) quel que soit le nombre d'objets — le raycasting CPU de Three.js est
inutilisable à cette échelle.

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
