# Validation

Ce document grandit à chaque semaine qui produit un résultat vérifiable. Il devient le
rapport de validation de la **semaine 21**, qui confronte le modèle aux données
orthogonales — DNA-FISH, chromatin tracing, DamID, SPRITE, GAM.

Règle : on consigne les chiffres, **y compris les mauvais**. Un rapport où tout est vert est
un rapport qui n'a pas cherché.

---

## S7 — Ensembles : ce qui est une propriété du génome, ce qui est un tirage

### Le dispositif

Le principe n° 1 du projet — *une structure Hi-C unique est un artefact statistique* — est une
déclaration d'intention jusqu'à ce qu'on sache **de combien**. La semaine 7 produit deux cents
repliements du même génome et mesure, bille par bille, quelle part de « cette bille est à telle
profondeur » est un énoncé sur la bille et quelle part sur le tirage.

`make ensemble` — 200 structures, quatre cœurs, 36 minutes. Le résultat est un
`ensemble.zarr` dont la forme porte déjà l'essentiel : **un génome, N repliements**. Le
découpage en billes, les rayons et la piste LAD sont écrits une seule fois ; seules les
coordonnées portent l'indice de structure.

Ce n'est pas une économie de place. Si chaque structure portait son propre génome, la
variabilité mesurée mélangerait la variabilité de repliement et celle du génome, sans aucun
moyen de les séparer après coup — et c'est exactement ce que faisait le code de la semaine 6,
où `build(seed=k)` tirait la piste LAD et la conformation de la même graine. Deux cents appels
auraient donné deux cents génomes. `lad_seed` les sépare, un test le verrouille, et la
génération vérifie l'empreinte du génome à **chaque** structure plutôt que de s'y fier.

### Il n'y a pas de repère commun, et c'est mesuré

Deux noyaux recuits séparément ne partagent ni orientation — rien ne distingue un axe dans une
sphère — ni placement des territoires, puisque chr7:a atterrit ailleurs à chaque tirage. Une
« variance de la position de la bille i » en x, y, z serait donc un nombre sans objet.

Plutôt que de l'affirmer, on l'a mesuré. Trois écarts entre deux structures, tous après
rotation optimale (réflexion permise, cf. § S5) :

| | RMSD |
|---|---:|
| vraie correspondance, bille i contre bille i | **4,72 µm** |
| billes permutées à l'intérieur de leur couche radiale | 5,36 µm |
| billes permutées librement — référence sans information | 5,35 µm |

Dans un noyau de rayon 5 µm, **88 % de l'écart survit au meilleur alignement possible**.
Une bille est typiquement plus loin de sa contrepartie que du centre du noyau.

**Une hypothèse fausse, corrigée par son propre témoin.** Le témoin « couches radiales » avait
été ajouté pour vérifier que le peu rattrapé par l'alignement était la stratification radiale
partagée. Il dit le contraire : mêler les billes à profondeur constante ramène pratiquement à
la référence sans information — cette explication ne vaut que **−1 % du gain**, c'est-à-dire rien. Ce
que la correspondance porte est ailleurs, dans la **compacité des territoires** : deux
structures ont chacune 46 blobs, et une rotation bien choisie en superpose quelques-uns. Le
témoin est resté, parce qu'il dit ça.

### Ce qui se reproduit d'un tirage à l'autre

Reste la position radiale, invariante par rotation, donc comparable. Décomposition de variance
à un facteur — la bille est le sujet, la structure le juge :

| | |
|---|---:|
| ICC — part de variance attribuable à la bille | **0,776** |
| dispersion entre billes | 538 nm |
| dispersion d'une bille sur l'ensemble | 276 nm |
| corrélation entre les profils radiaux de deux structures | **+0,777** |

78 % de la variance de profondeur tient à la bille, 22 % au tirage. La dernière ligne
dit la même chose sans passer par l'ICC, et c'est celle qu'il faut retenir : deux structures
indépendantes s'accordent à r = +0,777 sur *qui est profond et qui est superficiel*, et
pas plus.

Conséquence directe sur ce que le viewer aura le droit d'afficher : une structure isolée porte
les deux composantes sans les distinguer. Montrer une bille à sa profondeur sans montrer
276 nm de dispersion serait afficher un tirage en le présentant comme une mesure.

### Distributions radiales

| Quartile de contenu LAD | Fraction LAD | Rayon moyen | Dispersion entre billes |
|---|---:|---:|---:|
| Q1 — le moins LAD | 0,00–0,00 | 0,671 | 0,071 |
| Q2 | 0,00–0,07 | 0,673 | 0,074 |
| Q3 | 0,07–0,73 | 0,799 | 0,061 |
| Q4 — le plus LAD | 0,73–1,00 | **0,870** | 0,056 |

Les bornes LAD comptent autant que les rayons : **plus de la moitié des billes ont une fraction
LAD nulle ou quasi nulle**, si bien que Q1 et Q2 ne se distinguent pas — 0,671 contre 0,673. Le
modèle ne les sépare pas, et il n'a rien pour les séparer. La stratification qu'on mesure vient
entièrement du tiers supérieur de la piste.

L'auto-cohérence avec la piste LAD d'entrée vaut r = +0,771. **Ce nombre ne valide
rien** : la piste LAD est ce qui a fixé les rayons visés du modèle. Il dit seulement que le
solveur a fait ce qu'on lui demandait, et il est ici pour ça.

### Le médoïde, et pour quelle distance

Structure 125, graine 1125 : distance cumulée 1419,4 contre 1503,6 pour la plus excentrée.

Il n'existe pas de médoïde « de l'ensemble » dans l'absolu, seulement un médoïde **pour une
distance donnée**, et celle-ci ne voit que la profondeur des billes. Deux structures aux
territoires disposés tout autrement peuvent avoir le même profil radial ; ce médoïde ne les
distinguera pas. C'est assumé — la profondeur est ce que cet ensemble mesure — mais ça se dit.

### Ce que l'ensemble prédit d'une expérience Hi-C, et où il se trompe

Aucune matrice de contacts n'a jamais été montrée au modèle. P(s), la fraction trans et les
contacts d'homologues sortent de la seule géométrie : chaîne de sphères tangentes, volume
exclu, confinement, territoires. Ce sont donc les **seules grandeurs de cet ensemble qu'une
expérience puisse contredire** — et l'une d'elles est contredite.

| Seuil de contact | Contacts / structure | Trans | Homologues | Plateau P(s) |
|---:|---:|---:|---:|---:|
| **1,50** | **41 874** | **22,4 %** | **2,0 %** | **0,0253** |
| 1,25 | 29 503 | 19,3 % | 2,0 % | 0,0164 |
| 2,00 | 99 167 | 28,4 % | 2,1 % | 0,0644 |

Sur l'ensemble entier au seuil 1,50 : 6 498 563 contacts cis, 1 876 214 trans, 38 442 entre
homologues.

Un seuil de contact sous 1,15 ne mesure rien : c'est l'allongement maximal d'une liaison, donc
en dessous même deux billes voisines de chaîne ne « se touchent » pas et il ne reste que les
chevauchements résiduels. Mesuré à 1,00 : 145 contacts par structure au lieu de 41 800. Le
seuil est un paramètre — une ligature Hi-C n'exige pas que deux nucléosomes se touchent — mais
il a un plancher, et ce plancher est une propriété du modèle de chaîne.

**P(s) n'a pas une pente, elle en a trois.**

| Régime | Domaine | Pente |
|---|---|---:|
| chaîne | 1,4–4 Mb | −1,31 |
| polymère | 2,2–9 Mb | **−0,85** |
| territoire | 15–150 Mb | **−0,10** |

- Le régime **chaîne** est une construction : deux billes liées sont toujours à portée dès que
  le seuil dépasse l'allongement maximal, donc P(1) = 1 exactement. Ce qu'on y lit n'est pas du
  repliement, c'est le modèle de liaison.
- Le régime **polymère** donne -0,85, à comparer au ~s^-1 du Hi-C réel. C'est le bon ordre
  de grandeur, obtenu sans qu'aucune donnée de contact n'ait été fournie.
- Le régime **territoire** s'aplatit : -0,10, avec un plateau à P = 0,0253 au-delà de
  15 Mb. **Le Hi-C réel ne fait pas ça** : il continue de décroître.

C'est le résultat de la semaine, et c'est un échec instructif. Le modèle reproduit le *fait*
des territoires — la fraction trans tombe à 22,4 % quand le hasard pur en donnerait 98 % —
mais pas leur **organisation interne**. Une fois deux loci séparés de plus d'une quinzaine de
mégabases, leur probabilité de contact ne dépend plus de leur distance génomique : le
territoire est devenu un sac bien mélangé. Il manque au modèle tout ce qui structure l'intérieur
d'un chromosome — boucles, TADs, ségrégation compartimentale à l'échelle sub-chromosomique.

C'est exactement le périmètre de la semaine 8 (extrusion de boucles, polymère fin), et ce
plateau en est la mesure de départ : il faudra le voir décroître.

Deux lectures secondaires :

- **Les homologues ne s'apparient pas** : 2,0 % des contacts trans, contre les
  2,2 % attendus si une copie contactait ses 45 voisines au hasard. Conforme à ce
  qu'on sait des cellules somatiques humaines, où l'appariement homologue est l'exception.
- **La fraction trans tombe dans la fourchette du Hi-C réel** (~25–40 % pour GM12878 selon le
  filtrage). Comme la territorialité est une entrée du modèle (§ S6), ce n'est pas une
  prédiction indépendante ; ce qui l'est, c'est la *valeur* — rien ne garantissait qu'elle
  tombe là plutôt qu'à 5 % ou à 70 %.

### Le livrable qui manque

La feuille de route demande la **corrélation entre position radiale modélisée et LADs DamID
publiés**. Elle n'est pas dans ce rapport, et il faut être précis sur pourquoi : aucune entrée
DamID du manifeste n'a pu être récupérée, l'egress de l'environnement restant fermé
([`DATA_SOURCES.md` § 8](DATA_SOURCES.md)).

Ce qui existe : `ensemble.damid` lit un bedGraph, attribue à chaque bille la moyenne des scores
pondérée par le recouvrement, exclut les billes non couvertes — une bille sans mesure est une
absence, pas un zéro, et la mettre à zéro fabriquerait de la corrélation là où il n'y a pas de
donnée — et rend la corrélation avec la profondeur moyenne. Le chemin complet, fichier compris,
tourne sous test. Il se lance par :

```console
$ make ensemble DAMID=chemin/vers/lads.bedGraph
```

Il lui faut un fichier, pas une ligne de code de plus.

### Défauts trouvés

**1. Un génome par structure.** Décrit plus haut : `build(seed=k)` tirait la piste LAD *et* la
conformation de la même graine. Un ensemble construit ainsi aurait mesuré la variabilité de
deux cents génomes différents en croyant mesurer celle d'un repliement, et rien dans les
coordonnées n'aurait permis de s'en apercevoir. Corrigé par `lad_seed`, verrouillé par un test,
et vérifié à l'exécution sur chaque structure.

**2. Une copie sans homologue était son propre homologue.** La table des partenaires rend
`partner[k] = k` pour une copie non appariée, si bien que le test `partner[ci] == cj` comptait
**tous ses contacts cis** comme des contacts d'homologues. Trouvé par un test sur une chaîne
unique, où le compte annonçait 234 contacts d'homologues pour 0 contact trans — une
contradiction dans les termes, puisqu'un contact d'homologues est trans par définition. Le
masque manquant est `& ~same`.

**3. Les quartiles calculés sur les valeurs, pas sur les rangs.** Une bille de 750 kb sans le
moindre LAD est fréquente — plus de la moitié du génome modélisé. Les quantiles 0 % et 25 % de
la piste valent donc tous deux zéro, le groupe Q1 sort vide, et le rapport affichait `nan`.
Découpés par rang, les quatre groupes sont toujours définis ; le tableau affiche en plus les
bornes LAD de chacun, pour que les ex æquo se voient au lieu de se deviner.

**4. Une fenêtre d'ajustement trop étroite pour ajuster quoi que ce soit.** Le régime « chaîne »
couvrait deux points de séparation, et la pente sortait `nan`. Ce n'est pas le garde-fou qui est
en cause — il refuse d'ajuster une droite sur moins de trois points, et il a bien fait — mais la
fenêtre, choisie en mégabases sans vérifier combien de billes ça faisait à cette résolution.

**5. Relancer la commande effaçait la série.** `zarr.open_group(mode="w")` recrée le magasin
sans rien demander : un second `make ensemble` aurait effacé une demi-heure de calcul. Reprendre
est devenu le défaut, recréer demande `--fresh`, et un magasin d'une autre taille est refusé
plutôt que silencieusement remplacé.

### Une erreur de méthode, commise et rattrapée

Le premier essai de P(s), sur un ensemble de rodage à 3 Mb par bille, donnait une pente unique
de **-0,82** sur une fenêtre de 1,5 à 50 Mb. J'ai failli l'écrire telle quelle : « P(s) ∝ s^-0,8,
proche du s^-1 du Hi-C réel, obtenu sans qu'aucune donnée de contact n'ait été fournie ». Ç'aurait
été un beau résultat, et il aurait été faux.

À 3 Mb par bille, cette fenêtre couvre les séparations de 1 à 17 billes, c'est-à-dire
essentiellement le régime **chaîne**, celui qui est une construction. En regardant enfin la
courbe au lieu de son ajustement, elle se casse en trois : raide jusqu'à 2 Mb, en s^-0,85
jusqu'à 9 Mb, **plate ensuite**. La conclusion honnête est l'inverse de celle que j'allais tirer :
le modèle ne reproduit pas P(s), il s'en écarte franchement au-delà de 15 Mb, et c'est ça qui est
intéressant.

La leçon est celle de la semaine 5 sous une autre forme : une fenêtre d'ajustement se choisit en
regardant la courbe, jamais avant. Les trois régimes sont maintenant nommés dans le code, et
le rapport donne les trois pentes plutôt qu'une moyenne qui n'a de sens nulle part.

### Ce que ça ne dit pas

**Rien sur un vrai noyau**, et pour les mêmes raisons qu'en semaine 6 : aucune donnée de
conformation ne contraint ces positions, la piste LAD est synthétique, les longueurs de
chromosomes viennent de la table interne. Deux cents structures fausses restent fausses — un
ensemble ne rachète pas ses entrées, il en chiffre la dispersion.

**L'ICC dépend de la force du rappel radial.** Avec un rappel plus fort, les billes seraient
plus reproductiblement placées et l'ICC monterait, sans que le modèle soit meilleur pour
autant. Le nombre caractérise ce modèle-ci à ce réglage-là, et il se cite avec.

**Deux cents n'est pas un nombre magique** — il vient de la feuille de route, pas d'un calcul
de puissance. Alors on a regardé :

| N | ICC | r entre deux | dispersion | RMSD aligné | trans | plateau | pente polymère | médoïde |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 0,773 | +0,773 | 276 nm | 4,80 µm | 22,5 % | 0,0253 | −0,84 | graine 1000 |
| 100 | 0,775 | +0,774 | 276 nm | 4,70 µm | 22,5 % | 0,0253 | −0,85 | graine 1000 |
| 150 | 0,776 | +0,777 | 276 nm | 4,77 µm | 22,4 % | 0,0253 | −0,85 | graine 1125 |
| 200 | 0,776 | +0,777 | 276 nm | 4,74 µm | 22,4 % | 0,0253 | −0,85 | graine 1125 |

**Tout est convergé dès cinquante structures** : même ICC à trois décimales, même plateau à
quatre, même pente. Les deux cents de la feuille de route sont confortables, pas justes.

Une exception, et elle est instructive : **le médoïde change**, de la graine 1000 à la graine
1125 entre cent et cent cinquante structures. C'est attendu — il suffit qu'une structure
arrivée plus tard soit un peu plus centrale — et ça rappelle qu'un médoïde est le membre le
plus typique d'un échantillon, pas une structure privilégiée. Les statistiques de l'ensemble
sont stables ; l'identité de son représentant ne l'est pas.

### Reproduire

```console
$ make ensemble                       # 200 structures puis le rapport, ~36 min
$ make ensemble N=50                  # plus court
$ make ensemble DAMID=lads.bedGraph   # avec la corrélation qui manque
$ cd pipeline && PYTHONPATH=. .venv/bin/python -m pytest tests/test_ensemble.py
```

Le magasin est repris s'il existe : relancer la commande complète la série au lieu de la
recommencer.

---

## S6 — Noyau diploïde complet : imposé, hérité, mesuré

### Le dispositif

La semaine 5 reconstruisait *une* chaîne à partir de contacts. La semaine 6 change d'échelle :
le caryotype entier de GM12878 — **46,XX**, donc deux exemplaires de chr1–22 plus un X actif et
un X inactif — découpé en billes d'échelle TAD et placé dans une sphère de dix micromètres sans
que rien ne se traverse.

Une bille est un TAD, encore faut-il dire lequel. Dixon 2012 en compte ~2 200 sur le génome
haploïde, taille moyenne ~880 kb ; Rao 2014 en compte 9 274, médiane 185 kb. Ce ne sont pas deux
mesures du même objet à quatre près, ce sont deux définitions. À **750 kb par bille** on est à
l'échelle Dixon, et le noyau diploïde tient en **8 082 billes** de 166,8 nm de rayon, pour une
fraction volumique `phi = 0,30` — dans la fourchette 12–52 % mesurée par ChromEMT (Ou 2017).

Trois conditions dures et une préférence molle :

| | Contrainte | Nature |
|---|---|---|
| 1 | volume exclu, `d ≥ r_i + r_j` | dure |
| 2 | longueur de liaison **maximale**, `d ≤ 1,15 · (r_i + r_j)` | dure, unilatérale |
| 3 | confinement, `\|x\| ≤ R − r` | dure |
| 4 | rappel radial ordonné par contenu LAD | molle, coupée au polissage |

Et trois étapes : placer 46 domaines sphériques, y faire pousser une chaîne auto-évitante,
recuire le tout en faisant croître les rayons (Lubachevsky–Stillinger).

### Résultat : le critère de la semaine est atteint

`make nucleus`, six graines, valeurs par défaut.

| Graine | Chevauchement max | p99 | Paires en contact | Liaison la plus tendue | Secousses | Territorialité av. → ap. | Temps |
|-------:|------------------:|----:|------------------:|----------------------:|----------:|-------------------------:|------:|
| 0 | **0,749 %** | 0,475 % | 202 | +0,88 % | 0 | 38,0 → 27,0 | 25 s |
| 1 | **0,696 %** | 0,552 % | 174 | +0,89 % | 0 | 37,7 → 26,8 | 24 s |
| 2 | **0,703 %** | 0,558 % | 164 | +0,98 % | 0 | 38,0 → 26,9 | 24 s |
| 3 | **0,661 %** | 0,411 % | 197 | +0,99 % | 1 | 38,0 → 26,4 | 38 s |
| 4 | **0,735 %** | 0,563 % | 169 | +0,95 % | 1 | 37,9 → 27,3 | 38 s |
| 5 | **0,776 %** | 0,685 % | 118 | +0,86 % | 1 | 38,0 → 27,0 | 39 s |

Aucune bille hors du noyau dans aucun cas, et **8 082 billes** dans la fourchette demandée.

Le critère de la feuille de route — « un noyau diploïde, ~6 000–10 000 billes TAD, sans
interpénétration » — demande un seuil. Une seule tolérance, **1 % en relatif**, appliquée aux
trois conditions dures : le chevauchement (1 % du contact, soit 3,3 nm sur des billes de 334 nm
de diamètre), la tension de liaison, et le confinement. Une chaîne tendue 20 % au-delà de sa
limite viole une condition autant qu'une interpénétration, et c'est précisément la signature des
cages octaédriques décrites plus bas — un critère qui ne regarderait que le chevauchement les
laisserait passer.

À lire correctement : le polissage **s'arrête dès que le critère est tenu**. Ces valeurs disent
donc que la barre est franchie, pas jusqu'où le solveur descendrait — resserrer la tolérance le
fait descendre plus bas et coûter plus cher, ce qui est le réglage attendu et non une propriété
du modèle.

### Ce qui est imposé, ce qui est hérité, ce qui est mesuré

Trois familles de nombres qu'il ne faut pas confondre, et c'est la leçon de méthode de la
semaine.

**Imposé.** L'absence d'interpénétration, le confinement, l'intégrité des chaînes sont des
conditions que le solveur applique. Les mesurer vérifie le solveur, pas la biologie.

**Hérité.** Une relaxation sous contraintes ne fait **jamais** se croiser deux chaînes : la
topologie du noyau est celle de son initialisation. La territorialité est donc une **entrée**
du modèle, pas un résultat. Ce qui la justifie est la mitose — les chromosomes se décondensent
là où la télophase les a laissés, ils n'ont pas à se trier. Le seul énoncé vérifiable est que
le recuit la conserve : indice de territorialité **38,0 à l'initialisation, 26,4 à 27,0 après
recuit**, soit environ **70 % de conservation**, pour un entremêlement de 31 %.

Un modèle qui afficherait « territoires chromosomiques reproduits » sur cette base mentirait.

**Mesuré.** Reste la stratification radiale, qui est demandée par le terme 4 et peut ne pas être
obtenue. Elle l'est : corrélation de Pearson entre fraction LAD et rayon **+0,66 à +0,71**,
quartile le plus LAD à **0,87** du rayon nucléaire contre **0,66** pour le moins LAD. Sans le
terme 4, le même noyau donne **+0,091** et 0,694 contre 0,654 — la stratification vient bien du
modèle, et le contrôle le montre.

### Le résultat de fond : « en LAD » n'est pas « à la lamina »

Le DamID dit qu'environ 35 % du génome est en LAD. Il est tentant d'en conclure qu'un noyau
correct doit mettre 35 % de sa chromatine au contact de l'enveloppe. La géométrie dit non, et
elle le dit avec des nombres.

Une bille **touche** l'enveloppe si l'écart entre sa surface et celle-ci est sous un demi-rayon.
Cette bande fait 0,5 r d'épaisseur en position de centre, soit moins que l'écart entre deux
couches empilées (~1,63 r) : elle ne contient donc qu'une monocouche. Deux bornes l'encadrent.

- **À densité uniforme**, la part des billes qui s'y trouve est la part de volume de la coquille
  dans la boule accessible aux centres : `[(R−r)³ − (R−1,5r)³] / (R−r)³`, soit environ
  `1,5 · r/(R−r)`. À 750 kb par bille : **5,1 %**.
- **En empilement maximal**, au plus `eta · 4(R−r)²/r²` billes touchent à la fois, avec
  `eta ≤ 0,9069` la densité hexagonale. À 750 kb : **38 %**.

Les deux sont **linéaires en r**, donc dépendent de la résolution du modèle :

| Résolution | Billes | Rayon | Densité uniforme | Empilement max |
|-----------:|-------:|------:|-----------------:|---------------:|
| 3 Mb | 2 021 | 264,7 nm | 8,2 % | 57 % |
| **750 kb** | **8 083** | **166,8 nm** | **5,1 %** | **38 %** |
| 250 kb | 24 248 | 115,6 nm | 3,5 % | 27 % |
| 100 kb | 60 621 | 85,2 nm | 2,6 % | 20 % |
| 10 kb | 606 208 | 39,5 nm | 1,2 % | 9 % |

Ce ne sont pas que des formules : le noyau a été construit pour de bon à 250 kb par bille
(`make nucleus N=250000`). **24 244 billes** — le tableau en annonce 24 248, la borne se
calculant sur le total génomique là où le découpage réel arrondit chromosome par chromosome —
chevauchement maximal 0,541 %, liaison la plus tendue +0,94 %, aucune bille hors du noyau. La
part uniforme mesurée y vaut 3,5 %, la valeur du tableau, et le contenu LAD au contact tombe à
**5 %** contre 7 % à 750 kb : le sens et l'ordre de grandeur annoncés. La corrélation LAD-rayon
y monte même à +0,749, le modèle plus fin ayant plus de latitude pour trier.

Trois conséquences.

1. **Les deux phrases ne sont pas la même.** « 35 % du génome est en LAD » décrit une
   *séquence* ; « 35 % du génome touche la lamina » décrit une *configuration*, et la seconde ne
   découle pas de la première. À 10 kb par bille, l'empilement maximal lui-même plafonne à 9 %.
2. **Une fraction de LADs à la lamina n'est pas comparable entre deux modèles** de granularité
   différente. Sans la résolution, le nombre ne veut rien dire.
3. **Atteindre l'empilement maximal demande un noyau à croûte dense et intérieur creux.** C'est
   exactement ce que produisait la première version de ce solveur, par accident, et elle se
   bloquait (§ défauts).

Le modèle, lui, donne deux nombres qu'il faut lire ensemble :

- **2,6 % des billes** touchent l'enveloppe, contre 5,1 % à densité uniforme. La coquille est
  donc **appauvrie d'un facteur 2**, et c'est attendu : le confinement est dur, et une chaîne a
  moins de latitude contre la paroi qu'au milieu.
- ces 2,6 % de billes portent **7 % du contenu LAD** du génome modélisé. À répartition
  indifférente elles en porteraient 2,6 % : le LAD y est donc **enrichi 2,7 fois**.

Autrement dit, le terme 4 **trie** sans réussir à **remplir** — et il ne le pourrait pas, la
coquille étant plus petite que le contenu LAD quelle que soit la force du rappel. Les 7 % sont
bien en dessous des 34 % de couverture LAD de la piste, et du même ordre que la lecture en
cellule unique de Kind 2013, où chaque cellule ne contacte qu'une fraction du répertoire LAD et
non l'ensemble.

Et ces nombres sont **plus sensibles à la densité d'empilement qu'au terme LAD lui-même**, ce
qui est la mise en garde à retenir. Toutes les lignes ci-dessous tiennent le critère de 1 % :

| Variante | Billes à la lamina | Contenu LAD à la lamina | Part à densité uniforme |
|---|---:|---:|---:|
| défaut (`phi` 0,30, plancher de bruit 0,02) | **2,6 %** | **7 %** | 5,1 % |
| plancher de bruit 0,05 | 6,1 % | 16 % | 5,1 % |
| sans rappel radial du tout | 3,6 % | 4 % | 5,1 % |
| `phi` 0,15 | 0,3 % | 1 % | 4,0 % |
| `phi` 0,45 | 14,9 % | 36 % | 5,8 % |

De 0,15 à 0,45 de fraction volumique, le contenu LAD au contact passe de 1 % à 36 % — un
facteur 36 — là où mettre ou retirer le terme LAD le fait passer de 4 % à 7 %. Le réglage qui
décide n'est donc pas celui qu'on croit.

Deux lectures secondaires, qui tiennent toutes deux à la comparaison avec la part uniforme :

- **La densité fait basculer le signe.** À `phi` 0,15 la coquille est très appauvrie (0,3 %
  contre 4,0 %), à 0,30 elle l'est encore d'un facteur 2, à 0,45 elle est **enrichie** 2,6 fois
  (14,9 % contre 5,8 %). Une chaîne confinée peu dense évite la paroi ; serrée, l'encombrement
  lui impose une couche de surface. Le modèle reproduit cette transition sans qu'on la lui ait
  demandée — c'est du volume exclu contre une paroi dure, rien de plus.
- **Le terme LAD réduit le nombre de billes à la lamina tout en y concentrant le LAD.** Sans
  lui, 3,6 % des billes touchent et portent 4 % du LAD, soit à peu près rien de particulier.
  Avec lui, 2,6 % touchent mais portent 7 %. Il tire vers le centre plus de billes qu'il n'en
  pousse vers la paroi, et ce qu'il gagne est un **tri**, pas un remplissage.

### Défauts trouvés, tous dans mon propre code

Chaque nombre cité est celui mesuré au moment où le défaut a été trouvé, avec le solveur dans
l'état où il était alors : ils disent l'ampleur de chaque défaut, ils ne se comparent pas entre
eux.

**1. Une poussée radiale au lieu d'un rappel.** Le biais LAD s'écrivait comme un déplacement
vers l'extérieur appliqué à chaque pas. Il s'accumule : son effet dépend du nombre de pas et
non du modèle. Mesuré, 500 pas d'une poussée même faible (gain 0,03) plaquaient toutes les
billes LAD contre l'enveloppe et y formaient une croûte bloquée — chevauchement résiduel
**21 %**, contre **5 %** sans poussée du tout. C'est exactement l'erreur corrigée en semaine 5
sur la ségrégation A/B, refaite dans un module neuf. Un test de régression la ferme maintenant :
la même cible doit sortir de 300 pas et de 3 000.

**2. Deux corrections d'initialisation qui n'ont rien corrigé.** La marche persistante de pas
`2r` reposait régulièrement une bille sur une précédente : chevauchement initial **99 %**, deux
billes confondues, 56 512 paires. Première correction, rendre la marche auto-évitante : les
paires tombent à 32 760, et le maximum reste à **99 %**. Deuxième correction, une grille de
hachage commune aux 46 copies au lieu d'une par copie — les territoires se recouvrent, donc les
collisions qui comptent sont celles entre chaînes différentes : 34 250 paires, maximum toujours
à **99 %**.

Les deux corrections sont justes et sont restées. Aucune ne résout le problème, parce que le
problème n'était pas là : une marche gloutonne à cette densité finit toujours par se piéger, et
une fois piégée, la « moins mauvaise » direction pose la bille sur sa voisine. Ce qui a résolu,
c'est d'**arrêter d'exiger une configuration initiale valide** — faire croître les rayons
pendant le recuit (Lubachevsky–Stillinger) plutôt que de partir à taille pleine. Une leçon de
diagnostic : deux hypothèses plausibles, chacune vérifiée, chacune fausse sur la cause.

**3. Un tube non étanche pendant l'inflation.** Le solveur fait croître les rayons pour éviter
de partir d'un empilement impossible. Mais à l'échelle `s`, deux billes liées s'excluent au
rayon `s·r` en restant séparées d'au plus `stretch · 2r` : une troisième bille passe entre elles
dès que **`s ≤ stretch / 2`**. À `s = 0,50` avec `stretch = 1,15` (seuil 0,575), une liaison
restait bloquée 20 % au-delà de sa limite avec 5 % de chevauchement, quand tout le reste était
déjà à 0,04 %. Le défaut par défaut est écarté en restant au-dessus du seuil.

**4. Une liaison bilatérale.** Une longueur de liaison *imposée* se bat contre le volume exclu
dans les replis serrés, et c'est la liaison qui gagne : écart aux liaisons 0,3 %, chevauchement
résiduel **6 %**. Or le volume exclu est la condition, la liaison ne fait qu'empêcher la chaîne
de casser. Rendue unilatérale — une longueur *maximale* — la frustration disparaît.

**5. Une projection de Jacobi à gain 1.** Elle stagne dès que les corrections demandées à une
bille se compensent : plafond à **21 %** de chevauchement. La sur-relaxation à 1,6 passe.

**6. Des cages octaédriques, et √2 comme signature.** Le polissage a des points fixes. Sur
certaines graines il tournait ses 4 000 itérations pour rester à 5,1 % de chevauchement avec
une liaison tendue à **+20,42 %** au-delà de sa limite, quand d'autres graines finissaient sous
0,5 %. Le même **+20,42 %** est ressorti à l'identique sur cinq configurations et trois graines
différentes — une constante exacte n'est jamais un accident local.

En regardant la géométrie : six billes de chr3:a, toutes à **0,98 × contact** les unes des
autres, et trois liaisons de la même copie mais éloignées le long de la chaîne (|i−j| = 30, 31,
45, 75, 89) tendues à **1,3849 × contact**. Or 1,3849 / 0,98 = **1,4132 ≈ √2**.

C'est un **octaèdre régulier**. Six billes en contact mutuel forment une cage rigide, et trois
liaisons tombées sur ses trois diagonales ne peuvent plus se raccourcir : raccourcir une
diagonale écarterait les quatre billes de l'équateur, ce que les deux autres liaisons
interdisent. La cage se verrouille elle-même, les liaisons compriment les arêtes de 2 %, et le
rapport diagonale/arête d'un octaèdre étant √2, la tension se fige à
`√2 · 0,98 / 1,15 − 1 = 20,4 %` — d'où la constante.

Ce qui ne marche pas : prolonger le polissage (c'est un point fixe), et réchauffer globalement.
Quatre réchauffages à 0,10 sur la graine 3 faisaient tomber la corrélation LAD-rayon de +0,665
à +0,596 **sans régler le problème** : on ne quitte pas une cage de contacts par agitation,
parce que ce n'est pas un puits peu profond, c'est une structure rigide.

Ce qui marche : une **secousse locale** d'amplitude comparable au rayon d'une bille, appliquée
aux seules billes fautives et à leur voisinage immédiat, suivie d'un polissage. Elle casse la
cage sans défaire le noyau.

**7. Une part de volume normalisée par `R³`.** La coquille de contact est une part du volume
accessible aux *centres*, c'est-à-dire de la boule de rayon `R − r`, pas de `R`. L'erreur
sous-estimait la référence de 10 % et sortait un enrichissement de 1,10 sur des points pourtant
tirés uniformément — trouvé par le test écrit exprès pour ça.

### Une erreur de méthode, commise et corrigée

Après avoir diagnostiqué le défaut n° 3, j'ai écrit dans le code que le défaut était
« topologique et définitif » et fait **échouer** `build()` en dessous du seuil. C'était une
généralisation à partir d'un seul réglage. Vérification faite, la même configuration à `s = 0,50`
polie avec un plancher de bruit de 0,05 au lieu de 0,02 et un gain de liaison de 1,0 au lieu de
0,5 retombe à **0,35 %**. Le passage crée un défaut *difficile* à défaire, pas indéfaisable.

L'exception a donc été retirée : la condition est calculée, rapportée dans la sortie, et les
valeurs par défaut restent au-dessus du seuil parce que c'est gratuit — mais rien n'interdit
d'aller en dessous, puisque la mesure dit que ça peut marcher.

### Ce que ça ne dit pas

**Rien sur un vrai noyau.** Aucune donnée de conformation ne contraint ces positions : le noyau
n'est pas ajusté sur du Hi-C, il est seulement admissible géométriquement. La piste LAD est
**synthétique** — chaîne de Markov à deux états calibrée sur les statistiques publiées de taille
et de couverture, pas un fichier DamID. Les longueurs de chromosomes viennent de la **table
interne**, pas du `hg38.chrom.sizes` officiel, faute de réseau (§ `DATA_SOURCES.md` § 8). Tout
ça est estampillé dans le `.json` frère de chaque structure écrite, et imprimé par la CLI.

**La territorialité est une entrée**, redite ici parce que c'est le piège principal.

**Ce qui n'est pas modélisé** : le nucléole et les bras courts acrocentriques (marqués dans les
données, pas traités), les NADs, les corps nucléaires, la différence de repliement entre Xa et
Xi — les deux X sont nommés séparément mais construits pareil. Et une seule structure n'est pas
un ensemble : c'est le principe n° 1 du projet, et c'est le livrable de la semaine 7.

### Reproduire

```console
$ make nucleus                       # noyau par défaut, ~30 s
$ make nucleus N=250000              # 24 244 billes, échelle des domaines de Rao
$ cd pipeline && PYTHONPATH=. .venv/bin/python -m pytest tests/test_nucleus.py
```

---

## S5 — Contact → distance : l'exposant n'est pas une constante

### Le dispositif

La semaine 3 plantait des *enrichissements de contacts*. Elle ne plantait aucune position dans
l'espace, et ne pouvait donc rien dire d'une reconstruction géométrique. Ici on inverse : on
fabrique d'abord une conformation 3D — marche persistante confinée, volume exclu, compartiments
A au centre et B en périphérie — puis on en **dérive** la matrice de contacts par un modèle
direct explicite :

```
f(i, j) ∝ d(i, j)^(-gamma)
```

`gamma` est planté. La reconstruction, elle, balaie `alpha` dans `d ∝ f^(-alpha)` sans le
connaître. Si la méthode était exacte, l'optimum tomberait sur **alpha = 1/gamma**. C'est une
prédiction chiffrée qui peut échouer.

C'est la mesure que les données réelles ne permettront jamais : dans du Hi-C réel, la structure
3D est précisément l'inconnue.

### Résultat : alpha* ne vaut pas 1/gamma, et dépend de la profondeur

`make recon`, gamma = 3 donc 1/gamma = 0,333. Cinq conformations de 300 billes par profondeur.

| Profondeur | Densité | alpha* médian | Étendue sur 5 tirages | nRMSD min | Plateau à +5 % |
|-----------:|--------:|--------------:|-----------------------|----------:|---------------:|
| 500 000 | 37 % | **0,525** | 0,47 – 0,75 | 0,271 | 0,100 |
| 2 000 000 | 65 % | 0,475 | 0,40 – 0,55 | 0,204 | 0,050 |
| 8 000 000 | 90 % | 0,425 | 0,42 – 0,47 | 0,135 | 0,050 |
| 40 000 000 | 100 % | 0,375 | 0,38 – 0,40 | 0,079 | 0,000 |
| 200 000 000 | 100 % | **0,350** | 0,35 – 0,35 | 0,037 | 0,000 |
| sans bruit | — | 0,325 | — | **0,006** | — |

Trois lectures :

1. **alpha\* décroît vers 1/gamma avec la profondeur, sans jamais l'atteindre à profondeur
   finie.** Il absorbe la compression de dynamique due au bruit de Poisson : prendre une
   puissance négative d'un petit comptage bruité gonfle la distance moyenne, et un alpha plus
   grand ré-étale la gamme.
2. **Sans bruit, l'inversion exacte gagne** : alpha* = 0,325 (pas de balayage à 0,025 près de
   0,333) et nRMSD 0,006. Le modèle direct est donc bien inversé — l'écart vient entièrement
   de l'échantillonnage.
3. **Le plateau s'élargit quand les données se creusent.** À 200 M de contacts, le minimum
   est un point unique reproductible sur les cinq tirages. À 500 000, il s'étale sur 0,10 et
   l'optimum varie de 0,47 à 0,75 selon la conformation. **Là où il faudrait le plus calibrer
   alpha, c'est là qu'il est le moins déterminé.**

Conséquence pratique pour le projet : reprendre `alpha = 1/3` d'un article sans regarder sa
profondeur de séquençage n'est pas une convention, c'est une approximation non chiffrée. Alpha
devra donc être calibré sur les données effectivement utilisées, pas sur la littérature — ce
qui n'a pas encore eu lieu : la semaine 6 construit un noyau **sans aucune contrainte Hi-C**,
et cette calibration attend les données réelles.

### La complétion géodésique n'est pas un raffinement

Les paires sans contact observé ne sont pas à distance infinie. ShRec3D (Lesne 2014) les
complète par le plus court chemin dans le graphe des contacts. Mesuré en remplaçant cette
complétion par une grande constante :

| Profondeur | nRMSD avec géodésiques | nRMSD sans |
|-----------:|-----------------------:|-----------:|
| 300 000 | 0,204 | 0,611 |
| 1 200 000 | 0,171 | 0,291 |
| 10 000 000 | 0,098 | 0,210 |
| 50 000 000 | 0,057 | 0,107 |

Un facteur 2 à 3 sur toute la plage. C'est aussi ce qui rend la matrice compatible avec
l'inégalité triangulaire, que le MDS classique suppose et qu'une matrice trouée ne respecte pas.

### La chiralité est irrécupérable

Une matrice de distances ne détermine la structure qu'à une isométrie près : rotation,
translation, et **réflexion**. Une reconstruction miroir est une reconstruction correcte.
L'alignement de Procruste autorise donc explicitement la réflexion — l'interdire compterait
la moitié des solutions valides comme des échecs. Un test le verrouille.

### Deux erreurs de méthode, commises et corrigées

**La ségrégation A/B tenait par chance.** La première version de `chain()` l'obtenait par une
dérive appliquée pendant la marche, qui luttait contre la diffusion du hasard. Écart radial
mesuré : +0,157 à 180 billes, **+0,009 à 350**. Un test à une seule taille l'aurait déclarée
acquise. C'est désormais une force de rappel dans la relaxation, donc une propriété convergée :
pire cas +0,159 sur six tailles et six graines.

**J'ai sur-généralisé le premier résultat.** Après un seul balayage, j'ai affirmé que l'écart
à 1/gamma est positif et décroît avec la profondeur. En changeant la longueur de chaîne,
l'optimum à faible profondeur est passé de 0,250 à 0,575 — au-dessous puis au-dessus de
1/gamma. Le minimum était plat et l'argmin instable ; la médiane sur cinq conformations a
rétabli la tendance, mais l'affirmation initiale reposait sur un seul tirage dans une zone
plate. La conclusion tient, la démarche qui y menait était insuffisante.

### Ce que ça ne dit pas

Le modèle direct est une loi de puissance pure. Le Hi-C réel n'en est pas une : il porte du
bruit de ligature aléatoire, des régions non mappables, et une décroissance P(s) qui ne suit
pas le même exposant à toutes les échelles. Ces résultats fixent la **méthode et ses pièges**,
pas la valeur d'alpha à utiliser sur GM12878.

### Reproduire

```console
$ make recon              # le tableau complet
$ make recon N=500        # autre longueur de chaîne
$ make test               # 19 assertions figeant les propriétés
```

---

## S4 — Imposteurs de sphères et picking GPU

### Ce qui est prouvé, et ce qui ne peut pas l'être ici

`pnpm verify` exécute 16 assertions dans un vrai Chromium, sur un vrai contexte WebGL2.
**Toutes passent.** Mais Chromium headless rend via **SwiftShader**, un rasteriseur
*logiciel* : ses images par seconde ne disent rien d'un GPU. Le verdict porte donc sur la
**correction**, pas sur la performance.

Le critère de fin de la semaine 4 — « go/no-go signé sur le budget de rendu » — demande des
chiffres sur trois cibles matérielles réelles. Ils ne peuvent pas sortir d'ici.
`dist/spike.html` est la page qui les produit ; elle s'ouvre dans un navigateur sur la
machine cible et rend le tableau à reporter en [`ARCHITECTURE.md` § 4](ARCHITECTURE.md).

### Pourquoi des imposteurs plutôt que des sphères

On ne dessine pas de sphères : un quad par bille, orienté face caméra, et le fragment shader
résout l'intersection rayon–sphère par pixel. Une vraie géométrie de sphère, même grossière à
80 triangles, ferait **80 millions de triangles** pour un million de billes. Un imposteur en
fait deux millions, et le fragment shader ne travaille que sur les pixels réellement couverts.

L'essentiel tient dans une ligne du fragment shader : `gl_FragDepth` est calculé depuis le
**point d'impact réel**, pas depuis le quad. Sans elle on obtient des vignettes plates — deux
billes qui s'interpénètrent se découpent selon l'arête du quad au lieu de la courbe
d'intersection, et un nuage dense devient un collage.

### Les 16 assertions

| Famille | Ce qui est vérifié |
|---------|--------------------|
| Contexte | WebGL2 disponible ; cible multiple couleur + identifiants R32UI + profondeur RGBA32F complète |
| Shaders | compilation et édition de liens en GLSL ES 3.0 |
| Forme | le centre de la bille est touché ; le **coin du quad est jeté** — un disque, pas un carré |
| Profondeur | elle **bombe** : au centre 0,98098, au bord 0,98244. L'imposteur a du relief |
| Interpénétration | la surface la plus proche gagne des deux côtés du plan d'intersection, **et le résultat ne dépend pas de l'ordre de dessin** |
| Picking | chaque bille rend son identifiant ; le fond et le hors-cadre rendent `NO_HIT` |
| Échelle | picking exact parmi **250 005 instances** ; un identifiant occupant les 32 bits (4 294 967 295) survit au transport |

Le test d'ordre de dessin est celui qui discrimine réellement. Deux billes de même profondeur
qui se chevauchent : si la profondeur était celle du quad, les deux quads seraient coplanaires
et l'ordre de dessin déciderait du vainqueur. Avec la profondeur du point d'impact, c'est la
surface la plus proche qui gagne — donc le même résultat dans les deux ordres.

### Un test mal posé, corrigé

La première version sondait une grille de 32 400 billes sur un canevas de 256 pixels. Chaque
bille y couvre **0,4 pixel** : elle n'est physiquement pas pointable, et les quatre sondes
échouaient. Ce n'était pas un défaut du picking mais la résolution de l'écran.

Le test sonde désormais des billes franches placées devant, le nuage de 250 000 servant à ce
qu'il doit servir : montrer que l'exactitude ne dépend pas du nombre d'instances. Au passage,
l'assertion « identifiant au-delà de 2¹⁶ » portait sur un identifiant de 32 400 — elle ne
testait rien. Elle porte maintenant sur 4 294 967 295.

### Pourquoi le picking GPU et pas le raycasting

Le raycasting CPU teste le rayon contre chaque objet : O(n). À six cent mille billes, un clic
coûte plus cher qu'une frame. Le picking GPU laisse le rasteriseur faire le travail qu'il fait
déjà — il a de toute façon déterminé quel fragment est devant — et relit **un pixel**. Coût
constant quel que soit le nombre d'objets.

C'est ce qui rend tenable la fiche de la semaine 13 : cliquer une bille parmi six cent mille
doit coûter la même chose que cliquer parmi dix.

Prix à payer, consigné : `readPixels` synchronise le CPU sur le GPU. Invisible sur un clic
isolé ; au survol continu il faudra passer par un `PIXEL_PACK_BUFFER` et une lecture
asynchrone.

### Reste à faire, sur du matériel

```console
$ pnpm build          # produit dist/spike.html
# ouvrir dist/spike.html sur GPU desktop, iGPU portable, téléphone
```

Trois tableaux à rapporter, puis le go/no-go. Tant qu'ils manquent, le budget LOD de
`ARCHITECTURE.md` § 4 reste une **arithmétique d'octets vérifiée, pas une mesure de rendu**.

---

## S3 — Callers de conformation contre structure plantée

### Pourquoi du synthétique, alors qu'on veut du réel

Le Hi-C réel n'a pas de vérité terrain. Les appels publiés de Rao 2014 sont un *point de
comparaison entre deux méthodes*, pas une vérité : quand notre caller diverge, rien ne dit
lequel des deux a tort. On ne peut donc jamais prouver qu'un caller est **correct** sur des
données réelles — seulement qu'il est **d'accord** avec un autre.

Sur une structure plantée, on connaît la réponse parce qu'on l'a écrite. C'est la seule étape
du projet où « correct » a un sens.

L'ordre est donc : prouver la correction ici, puis comparer à Rao sur les vraies données.
Le synthétique ne remplace pas les données réelles, **il les précède**.

### Le modèle génératif

Pour chaque paire de bins (i, j) d'un chromosome :

```
E[i,j] = profondeur
       × (|i-j| + 1)^(-1,1)      décroissance P(s) avec la distance
       × (1 ± 0,35)              compartiment : même type ou non
       × (1 + 0,9)               si i et j sont dans le même TAD
       × (1 + bosse)             aux coins des TADs, là où CTCF forme une ancre
       × biais_i × biais_j       biais de couverture log-normal, σ = 0,35

comptage[i,j] ~ Poisson(E[i,j])
```

Deux contraintes de réalisme que le premier jet n'avait pas, et qui changent la mesure :

1. **La fin du chromosome n'est pas une frontière.** La planter fabrique un vrai positif que
   personne ne peut appeler — il n'y a pas de bin derrière — et écrase le rappel mesuré sans
   que rien ne soit cassé dans le caller.
2. **Les transitions de compartiment tombent sur des frontières de TAD.** En tirant les deux
   indépendamment, un changement A/B au milieu d'un TAD crée une vraie insulation que le
   caller détecte correctement, mais qui n'est pas dans la vérité. On mesure alors comme une
   erreur du caller ce qui est une erreur du modèle. Avant correction : précision 68 %. Après,
   sans toucher au caller : **100 %**.

### Résultats

`make hic-validate` — 1 000 bins × 10 kb = 10 Mb, 3,5 M contacts, 22 TADs, 19 boucles.
Outils : `cooler` 0.10.4, `cooltools` 0.7.1.

| Étape | Mesure | Résultat |
|-------|--------|----------|
| Équilibrage ICE | CV des marginales | 0,422 → **0,039** |
| Équilibrage ICE | poids vs biais planté | **r = +0,982** |
| Compartiments A/B | accord sur 1 000 bins | **100 %** |
| Frontières de TAD | rappel / précision (fenêtre 100 kb, ±1 bin) | **100 % / 100 %** |
| Boucles | rappel / précision (±2 bins) | **89 % / 94 %** |

### Ce que la validation a appris

**La fenêtre d'insulation doit être nettement plus petite que le TAD cherché.** Sur des TADs
de 450 kb en médiane, le rappel s'effondre quand la fenêtre approche leur taille :

| Fenêtre | 100 kb | 150 kb | 200 kb | 300 kb | 400 kb | 600 kb |
|---------|--------|--------|--------|--------|--------|--------|
| Rappel  | 100 %  | 100 %  | 100 %  | 100 %  | 77 %   | 42 %   |

Une fenêtre trop large enjambe la frontière et lisse le minimum qu'on cherche. C'est le
réglage qu'on aurait copié d'un tutoriel sans regarder l'échelle de ses propres domaines.

**Les appels de frontière ne tombent pas au bin exact.** À tolérance nulle, le rappel chute à
38 % ; à ±1 bin il est de 100 %. L'insulation est un score lissé, pas un détecteur d'arête.
Toute comparaison de frontières — y compris celle à venir contre Rao 2014 — doit donc porter
une tolérance explicite, sinon elle mesure du bruit.

**Le signe du vecteur propre est arbitraire.** La décomposition rend une partition, pas une
étiquette : rien dans la matrice ne dit quel côté est actif. Il faut une information
extérieure — contenu GC ou densité génique — pour l'orienter. Deux propriétés distinctes,
donc deux tests distincts : que la partition soit retrouvée, et que la piste de phasage fixe
le bon sens.

### Reproduire

```console
$ make hic-validate          # le tableau ci-dessus
$ make test                  # 12 assertions figeant ces seuils
```

Les seuils des tests sont stricts à dessein. Un caller qui passe de 100 % à 85 % de rappel
sur une structure plantée a un bug, pas un mauvais jour.

### Ce que ça ne dit pas

Rien sur les données réelles. Une matrice plantée est propre : pas de régions non mappables,
pas d'aberration caryotypique, pas de variation de profondeur entre chromosomes, pas de
bruit de ligature aléatoire. Ces callers sont **corrects** ; reste à savoir s'ils sont
**robustes**. C'est la comparaison à Rao 2014, et elle attend l'accès réseau
([`DATA_SOURCES.md` § 8](DATA_SOURCES.md)).

---

## S2 — Socle 1D

Latence de requête sur 1 000 000 d'intervalles : p99 **0,819 ms** à l'échelle d'un locus,
attributs compris ; index seul sous 0,11 ms à toutes les échelles. Tableau complet et
conséquences pour le format `.g3d` dans
[`ARCHITECTURE.md` § 8](ARCHITECTURE.md#8-socle-1d--index-dintervalles-et-mesures).

## S1 — Vérification des données

13 assertions sans réseau (`make selftest`) : rejet d'une source altérée, refus d'une
accession non résolue, idempotence, non-conservation d'un fichier corrompu.
