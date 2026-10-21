# Installation

Deux niveaux, volontairement séparés.

## Socle 1D — numpy seulement

Le magasin d'intervalles et `geno query` n'ont besoin que de numpy. Ils tournent avec
l'interpréteur système, sans environnement dédié.

```console
$ pip install numpy pytest
$ make build-store
$ make query Q=chr7:5,527,000-5,530,600
```

Cette séparation n'est pas cosmétique : le socle 1D est ce dont dépend le viewer, et il ne
doit pas traîner derrière lui scipy, h5py et une pile scientifique complète.

## Pile de conformation — cooler et cooltools

Nécessaire pour le Hi-C (`make hic-validate`, `make recon`), pour le modèle de noyau
(`make nucleus`) — qui n'a besoin que de numpy et de `scipy.spatial` —, pour les ensembles
(`make ensemble`), qui y ajoutent zarr, et pour l'échelle fine (`make fine`), qui y ajoute
OpenMM. Elle vit dans `pipeline/.venv`, que les cibles `make`
utilisent automatiquement si elle existe ; sinon elles retombent sur `python3` et signalent
proprement ce qui manque.

```console
$ cd pipeline
$ python3 -m venv .venv
$ .venv/bin/pip install --upgrade pip wheel
$ .venv/bin/pip install "setuptools==59.8.0"
$ .venv/bin/pip install --no-build-isolation asciitree
$ .venv/bin/pip install cooler cooltools "pandas<3" "zarr<3" openmm pytest
```

### Pourquoi ces épinglages

**`setuptools==59.8.0` puis `asciitree --no-build-isolation`.** `asciitree` est une
dépendance transitive de `cooler` qui utilise des API `setuptools` retirées depuis. Avec un
setuptools moderne, sa construction échoue sur `AttributeError: install_layout`, et
l'installation de cooler s'arrête là. On lui fournit donc un setuptools d'époque, une seule
fois, dans le venv.

**`zarr<3`.** zarr 3 a remplacé `create_dataset` par `create_array` et changé la sémantique
de `open_group`. Le magasin d'ensembles (`make ensemble`) écrit en API 2 ; l'épinglage évite
une rupture silencieuse au premier `pip install -U`.

**`openmm`, et pas `polychrom`.** `polychrom` est le wrapper habituel pour les polymères de
chromatine, mais il n'est pas publié sur PyPI et l'environnement de développement n'atteint
pas GitHub. Les cinq termes du champ de force de la semaine 8 sont donc écrits directement
sur l'API OpenMM — une quarantaine de lignes, avec l'avantage que chacun est lisible et
testable séparément. OpenMM 8.6 s'installe en revanche depuis PyPI sans conda, ce qui n'a pas
toujours été le cas.

Une contrainte de son API a façonné le format des instantanés d'extrusion :
`updateParametersInContext` modifie les longueurs et les raideurs d'une liaison, **jamais ses
particules**. Déplacer une liaison de cohésine au fil de l'extrusion lève `The set of
particles in a bond has changed`. D'où `fine/polymer.bond_catalogue`, qui déclare d'emblée
toutes les paires du parcours à raideur nulle et n'allume que celles tenues à l'instant
présent ([`ARCHITECTURE.md` § 12](ARCHITECTURE.md)).

**`pandas<3`.** `cooltools` 0.7.1 appelle `.idxmin()` sur une colonne entièrement NA dans
`api/dotfinder.py:977`. pandas 2 émet un `FutureWarning` ; **pandas 3 lève une
`ValueError`** et la détection de boucles échoue avec le message peu parlant
`Encountered all NA values`. Diagnostiqué, pas deviné : le warning de pandas 2 pointe la
ligne exacte.

## Format `.g3d` — aucune dépendance nouvelle

L'écriture n'utilise que numpy, déjà là, et `zlib` et `hashlib` de la bibliothèque standard. La lecture,
côté navigateur comme côté Node 22, n'utilise que des API natives : `DecompressionStream
('deflate-raw')` pour décompresser, `crypto.subtle` pour vérifier les empreintes, `fetch` pour
les requêtes Range. La mesure de temps pilote Chromium par le protocole DevTools sur le
`WebSocket` natif de Node 22, sans Playwright ni Puppeteer. Le dépôt garde donc ses deux seules
dépendances JavaScript, esbuild et TypeScript — jusqu'à la semaine 14, qui embarque igv.js et
HiGlass dans la seule page des vues liées (ci-dessous) ; le test de lecture croisée déclare les deux
symboles Node dont il a besoin plutôt que d'ajouter `@types/node`.

## Vues liées — igv.js et HiGlass (semaine 14)

Les premières dépendances JavaScript d'exécution du dépôt, toutes dans `apps/linked/` et toutes
épinglées : `igv` 3.8.10, `higlass` 2.3.5 et ses pairs `react` / `react-dom` 18.3.1 et `pixi.js`
6.5.10 (`pnpm install`). Elles ne touchent ni `first.html` ni `dive.html` : seule `vues.html` les
embarque (4,3 Mo). Quatre particularités, chacune vue avant d'être contournée :

- **igv** s'importe par `igv/dist/igv.esm.js` : le champ `browser` de son `package.json` désigne la
  version UMD, qu'esbuild prendrait sinon, sans export par défaut. Ses types ne couvrent pas ce
  chemin (`apps/linked/shims.d.ts`).
- **React** lit `process.env.NODE_ENV` : `tools/build.mjs` le définit à `"production"`.
- **HiGlass** a une feuille de style (`hglib.css`), incorporée à la page. Son data fetcher maison
  doit exposer un champ `dataConfig` — HiGlass le lit par son nom.
- **Le serveur local** (`tools/serve.mjs`) lit les fichiers de plus de 64 Mo par plages, sans les
  charger : igv lit le `.2bit` de 835 Mo par requêtes Range.

`make linked` (pipeline) n'ajoute aucune dépendance Python : numpy et scipy, déjà là.

## Vérifier

```console
$ make test            # 283 tests (… + 10 liftover, + 23 fiche, + 1 parseurs, + 7 vues liées)
$ make selftest        # 13 assertions sur la vérification d'empreintes, sans réseau
$ make hic-validate    # plante une structure Hi-C connue, valide les callers dessus
```

`make test` saute automatiquement les tests de conformation si cooler et cooltools sont
absents de l'interpréteur courant — le socle 1D reste testable sans la pile lourde.

## Réseau

`make data` a besoin d'un accès sortant vers `hgdownload.soe.ucsc.edu`, `ftp.ebi.ac.uk`,
`data.4dnucleome.org` et `ftp.ncbi.nlm.nih.gov`. Depuis le 2026-10-18, ils répondent ; voir
[`DATA_SOURCES.md` § 8](DATA_SOURCES.md) pour ce qui est récupéré et verrouillé.

La semaine 12 a besoin de `make data-core data-seq` (≈ 1,7 Go : `hg38.2bit` et `hs1.2bit`) et
de GENCODE (`gencode.annotation.gtf`, 80 Mo). Un téléchargement coupé par le réseau se reprend
par plage (`curl -C -`) ; c'est l'empreinte, confrontée au `md5sum.txt` du dépôt puis
verrouillée par `make data-lock`, qui dit que le fichier est entier. Sans `hg38.2bit`, `make g3d`
écrit un fichier sans séquence, et la plongée s'arrête aux nucléosomes.

La semaine 13 (la fiche) ajoute 24 sources, ≈ 750 Mo, dont ClinVar (200 Mo), le GWAS Catalog
(73 Mo), gnomAD (95 Mo), MONDO (53 Mo) et les bigWig de 4DN et d'ENCODE ; elles s'ajoutent à
`ensembl.org`, `encodeproject.org`, les S3 ouverts de 4DN, `storage.googleapis.com` (gnomAD,
GTEx), `purl.obolibrary.org` et `jaspar.genereg.net`. Le `.hic` de Rao 2014, plusieurs dizaines de
Go, est dans le même tier que des fichiers dont la fiche a besoin : `ONLY="id …"` dans
l'environnement de `make data-…` ne récupère que ces identifiants. `make fiche` a besoin de
`pyBigWig` (`pip install -e '.[hic,fiche]'`) ; sa première lecture des sources prend ≈ 4 min et
se met en cache dans `data/cache/fiche/`, indexé par les empreintes du manifeste.

Depuis le 2026-10-19, `hgdownload.soe.ucsc.edu` échoue au TLS (« unable to get local issuer
certificate ») ; rien de ce qui précède n'en dépend plus, sauf un premier `make data-core`.
