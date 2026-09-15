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

Nécessaire pour le Hi-C (`make hic-validate`, `make recon`) et pour le modèle de noyau
(`make nucleus`), qui n'a besoin que de numpy et de `scipy.spatial` mais vit dans le même
niveau. Elle est dans `pipeline/.venv`, que les cibles `make` utilisent automatiquement si
elle existe — sinon elles retombent sur `python3` et signalent proprement ce qui manque.

```console
$ cd pipeline
$ python3 -m venv .venv
$ .venv/bin/pip install --upgrade pip wheel
$ .venv/bin/pip install "setuptools==59.8.0"
$ .venv/bin/pip install --no-build-isolation asciitree
$ .venv/bin/pip install cooler cooltools "pandas<3" pytest
```

### Pourquoi ces deux lignes bizarres

**`setuptools==59.8.0` puis `asciitree --no-build-isolation`.** `asciitree` est une
dépendance transitive de `cooler` qui utilise des API `setuptools` retirées depuis. Avec un
setuptools moderne, sa construction échoue sur `AttributeError: install_layout`, et
l'installation de cooler s'arrête là. On lui fournit donc un setuptools d'époque, une seule
fois, dans le venv.

**`pandas<3`.** `cooltools` 0.7.1 appelle `.idxmin()` sur une colonne entièrement NA dans
`api/dotfinder.py:977`. pandas 2 émet un `FutureWarning` ; **pandas 3 lève une
`ValueError`** et la détection de boucles échoue avec le message peu parlant
`Encountered all NA values`. Diagnostiqué, pas deviné : le warning de pandas 2 pointe la
ligne exacte.

## Vérifier

```console
$ make test            # 55 tests (43 socle 1D + 12 conformation)
$ make selftest        # 13 assertions sur la vérification d'empreintes, sans réseau
$ make hic-validate    # plante une structure Hi-C connue, valide les callers dessus
```

`make test` saute automatiquement les tests de conformation si cooler et cooltools sont
absents de l'interpréteur courant — le socle 1D reste testable sans la pile lourde.

## Réseau

`make data` a besoin d'un accès sortant vers `hgdownload.soe.ucsc.edu`, `ftp.ebi.ac.uk`,
`data.4dnucleome.org` et `ftp.ncbi.nlm.nih.gov`. Voir
[`DATA_SOURCES.md` § 8](DATA_SOURCES.md) pour l'état actuel et les deux façons de débloquer.
PyPI reste joignable indépendamment de ces quatre hôtes.
