# Validation

Ce document grandit à chaque semaine qui produit un résultat vérifiable. Il devient le
rapport de validation de la **semaine 21**, qui confronte le modèle aux données
orthogonales — DNA-FISH, chromatin tracing, DamID, SPRITE, GAM.

Règle : on consigne les chiffres, **y compris les mauvais**. Un rapport où tout est vert est
un rapport qui n'a pas cherché.

---

## S9 — Format et streaming : 1,41 s en Slow 4G, et une marge que TLS mangerait

### Le critère, et ce qu'il veut dire

La feuille de route demande : **niveau noyau < 500 ko, premier rendu < 1,5 s en 4G**. Le
premier chiffre est tenu de très loin. Le second dépend de ce qu'on appelle 4G, et la réponse
change le verdict. Les préréglages de DevTools en donnent deux :

| Profil | Latence par requête | Débit descendant |
|---|---:|---:|
| Fast 4G | 165 ms | 1 012,5 ko/s |
| Slow 4G | 562,5 ms | 180 ko/s |

`tools/loadtime.mjs` passe ces valeurs explicitement à `Network.emulateNetworkConditions` :
ce sont **ces nombres-là**, et non le nom du profil, qui définissent la mesure. En Slow 4G, les deux allers-retours incompressibles — la page, puis les données —
coûtent **1 125 ms avant le premier octet utile**. Il reste 375 ms pour transférer, décoder et
dessiner, soit ~34 ko sur le fil. Tout le format en découle : le premier rendu doit tenir en
**une** requête Range, et tout préfixe du fichier doit être utile.

### Ce qui est dans le fichier

`make g3d` écrit `data/g3d/gm12878.g3d` : **56 881 octets**, trois niveaux, toutes empreintes
vérifiées à la relecture.

| Niveau | Billes | Source | Erreur max | Rôle |
|---|---:|---|---:|---|
| aperçu, 3 Mb | 2 032 | médoïde S7, quatre billes fusionnées | 2,57 nm | premier rendu |
| noyau, 750 kb | 8 082 | médoïde S7 sur 200 structures | 2,80 nm | — |
| fin, 2 kb | 2 000 | médoïde S8 sur 1 200 conformations | 0,83 nm | chr7:4–8 Mb |

Le niveau noyau entier pèse **35 ko** contre 500 ko admis. Son premier rendu — copie et
positions — tient en 26,4 ko.

### La précision stockée suit l'incertitude, et c'est là que sont les octets

La première version quantifiait les positions sur toute la plage de 16 bits : pas de 0,14 nm,
**40,6 ko** de positions pour le noyau, à peine moins que les 48,5 ko bruts. Les positions sont
la colonne chère parce qu'à ce pas l'écart entre billes voisines porte de nombreux bits par axe
— des bits que le modèle n'a pas : la profondeur d'une bille varie de 276 nm (médiane) d'un
tirage à l'autre (S7). La règle retenue fixe l'erreur maximale à **1 % de
l'incertitude médiane du niveau** :

| Erreur max | Pas | Positions du noyau | Fraction de l'incertitude médiane |
|---:|---:|---:|---:|
| 0,07 nm | 0,14 nm | 40 642 o | 0,03 % |
| 0,50 nm | 1,00 nm | 32 528 o | 0,18 % |
| 1,00 nm | 2,00 nm | 28 960 o | 0,36 % |
| **2,76 nm** | **5,52 nm** | **26 387 o** | **0,99 %** |
| 5,00 nm | 10,0 nm | 23 872 o | 1,78 % |
| 10,0 nm | 20,0 nm | 20 938 o | 3,57 % |

La courbe est logarithmique : chaque moitié de précision rend environ un bit par axe. La règle
a été choisie **après** avoir vu la courbe, ce qu'il faut dire : 1 % n'est pas réglé sur un
budget d'octets, c'est le seuil usuel du négligeable devant une incertitude.

Le reste se compresse presque à rien, parce que le delta modulaire par plan transforme des
suites régulières en zéros : `start` et `end` passent de 32 ko à 0,7 ko chacun, l'index des
identifiants de 32 ko à 51 octets.

### Le premier rendu, mesuré

Chromium headless, page servie en gzip, `.g3d` servi brut avec Range, cache désactivé, dix
passages par case. Chaque temps est compté depuis l'origine de navigation, et **chaque
passage vérifie son rendu** — au moins 5 % des pixels couverts par une bille, et un picking au
centre qui tombe sur une bille — avant de publier son chiffre.

| Variante | Profil | Requêtes .g3d | Octets .g3d | Document | Données | **Premier rendu** | p90 |
|---|---|---:|---:|---:|---:|---:|---:|
| noyau, préfixe 16 Kio | Fast 4G | 2 | 42 792 | 179 ms | 579 ms | 819 ms | 837 ms |
| noyau, préfixe 16 Kio | Slow 4G | 2 | 42 792 | 607 ms | 1 999 ms | **2 253 ms** | 2 302 ms |
| noyau, préfixe 64 Kio | Fast 4G | 1 | 44 627 | 176 ms | 406 ms | 654 ms | 668 ms |
| noyau, préfixe 64 Kio | Slow 4G | 1 | 44 627 | 606 ms | 1 433 ms | **1 685 ms** | 1 750 ms |
| noyau, préfixe ajusté | Fast 4G | 1 | 29 153 | 177 ms | 393 ms | 637 ms | 674 ms |
| noyau, préfixe ajusté | Slow 4G | 1 | 29 153 | 607 ms | 1 353 ms | **1 602 ms** | 1 627 ms |
| aperçu, préfixe 16 Kio | Fast 4G | 1 | 16 384 | 178 ms | 380 ms | 521 ms | 558 ms |
| aperçu, préfixe 16 Kio | Slow 4G | 1 | 16 384 | 606 ms | 1 282 ms | **1 412 ms** | 1 428 ms |
| aperçu, préfixe 64 Kio | Fast 4G | 1 | 56 881 | 177 ms | 421 ms | 550 ms | 554 ms |
| aperçu, préfixe 64 Kio | Slow 4G | 1 | 56 881 | 606 ms | 1 508 ms | **1 646 ms** | 1 671 ms |
| aperçu, préfixe ajusté | Fast 4G | 1 | 11 149 | 177 ms | 374 ms | 512 ms | 574 ms |
| aperçu, préfixe ajusté | Slow 4G | 1 | 11 149 | 606 ms | 1 249 ms | **1 381 ms** | 1 389 ms |

Sans bridage, toutes les variantes rendent entre 164 et 294 ms. « Préfixe ajusté » : le
préfixe vaut exactement la fin du premier rendu, que donne le préambule ; une page déployée
avec son fichier connaît ce nombre, un client générique non.

**En Fast 4G, tout passe**, de 512 à 819 ms. **En Slow 4G, seul l'aperçu passe** : 1 381 ms
avec un préfixe ajusté, **1 412 ms avec un simple préfixe de 16 Kio** — c'est donc ce que fait
la page par défaut, sans rien savoir du fichier.

Trois lectures du tableau :

- **Une requête de trop coûte plus que tout le reste.** Le noyau lu avec un préfixe de 16 Kio,
  trop court, demande une seconde requête : +651 ms sur le préfixe ajusté. Une latence de
  562,5 ms ne se compense par aucun gain d'octets.
- **Un préfixe trop long se paie aussi.** 64 Kio sur un fichier de 57 ko, c'est télécharger
  tout le fichier pour en dessiner 11 ko : 1 646 ms au lieu de 1 381.
- **Le noyau n'échoue pas sur le réseau.** Ses données sont là à 1 353 ms, sous la barre ; ce
  sont les ~250 ms de dessin par SwiftShader qui le font passer au-dessus. Sur un vrai GPU, ce
  dessin coûterait probablement quelques dizaines de millisecondes — mais ce n'est pas mesuré
  ici, donc ce n'est pas revendiqué.

### La règle fixée d'avance, appliquée

Avant la mesure, le plan disait : l'aperçu n'entre dans le fichier par défaut que s'il fait
passer un profil nommé sous 1,5 s. Il fait passer Slow 4G et le noyau seul non : il est donc le
niveau du premier rendu par défaut, et `gm12878-noyau.g3d`, sans lui, reste écrit pour
comparaison. Il coûte 12 ko de fichier, qui ne sont lus que par qui les demande.

### La marge est de 88 ms, et voici ce qui la mangerait

Le bridage de DevTools agit **par requête** : il n'émule ni la poignée de main TCP, ni TLS,
ni le démarrage lent de TCP. Un vrai réseau Slow 4G en HTTPS ajoute au moins un aller-retour
de TLS 1.3 pour la connexion à l'origine — 562,5 ms à cette latence. **Avec lui, aucune
variante ne passerait sous 1,5 s.** Le chiffre mesuré est une **borne basse** côté réseau,
quand SwiftShader en fait une **borne haute** côté rendu ; les deux erreurs vont en sens
contraires et ne s'annulent pas forcément.

Le verdict honnête est donc : **critère tenu sous le profil Slow 4G de DevTools tel qu'il est
défini, avec 88 ms de marge ; tenu largement en Fast 4G ; non tenu sur un vrai Slow 4G en
HTTPS**, où il faudrait gagner un aller-retour — par exemple en servant le préfixe du premier
rendu dans la page elle-même, ce qui couple page et données et n'est pas fait ici.

### Le budget du § 4, en octets écrits

`make g3d-bench` remplace l'arithmétique sur float32 par des octets écrits, attributs compris :

| Bille | Billes | Structure | Chunks | Premier rendu | Niveau entier | o/bille |
|---:|---:|---|---:|---:|---:|---:|
| 1 Mb | 6 058 | noyau S6 réel | 1 | 19,6 Kio | 27,8 Kio | 4,69 |
| 250 kb | 24 244 | noyau S6 réel | 8 | 76,1 Kio | 109,6 Kio | 4,63 |
| 100 kb | 60 618 | noyau S6 réel | 8 | 173,5 Kio | 224,0 Kio | 3,78 |
| 25 kb | 242 472 | raffinement synthétique | 64 | 677,1 Kio | 929,6 Kio | 3,93 |
| 10 kb | 606 180 | raffinement synthétique | 190 | 1 610 Kio | 2 218 Kio | 3,75 |

Entre 3,8 et 4,7 octets par bille à toutes les échelles, contre 12 pour les seules positions en
float32. Sous 100 kb, les structures sont un **raffinement synthétique** du noyau à 100 kb —
chaque bille découpée en filles le long d'un pont brownien — qui suffit à compter des octets et
ne dit rien de la biologie.

### Le niveau fin, posé, et ce que la pose mesure

Le médoïde de la semaine 8 est posé dans le repère du noyau par Kabsch, réflexion permise,
**sans échelle**, sur six billes de chr7:a qu'il couvre à moitié au moins.

| | |
|---|---:|
| RMSD de la pose rigide | **394,7 nm** |
| échelle qui serait optimale | ×1,58 |
| RMSD si on l'appliquait | 382,5 nm |

Le chiffre important est le dernier : **une échelle ne rattrape presque rien**. L'écart entre
les deux modèles n'est pas une affaire de taille mais de **forme** — le modèle fin est une
pelote confinée dans 435 nm de rayon, le noyau de la semaine 6 une chaîne qui gonfle avec la
distance génomique. C'est la différence de pente de `R(s)` mesurée en semaine 8 (s^0,68 contre
s^−0,02), vue autrement. La pose ne corrige pas le raccord ; elle le **mesure**, et le fichier
transporte ce bilan avec la transformation. Poser sur chr7:a plutôt que chr7:b est arbitraire —
le modèle fin n'a pas d'haplotype — et le fichier le dit aussi.

### Défauts trouvés

**Une précision que le modèle n'a pas.** Le format initial quantifiait sur toute la plage de
16 bits ; c'était l'habitude du float32 sous une autre forme. Mesuré, ça coûtait 14 ko sur
les 34 que Slow 4G laisse au premier rendu.

**Une origine arrondie qui dépassait l'erreur annoncée.** L'en-tête arrondissait l'origine de
quantification au 1e-4 nm, assez pour que l'erreur relue dépasse l'erreur annoncée de
3·10⁻⁵ nm. Invisible à l'œil, vu par un test exact.

**Trois lectures au lieu de deux.** Avec un préfixe trop court pour l'en-tête, le lecteur
complétait l'en-tête, puis relisait les colonnes du premier rendu — alors que le préambule
donne déjà où elles finissent. Trouvé par la lecture croisée ; la lecture qui complète
l'en-tête va désormais jusqu'à la fin du premier rendu.

**Deux variantes de mesure confondues avec le banc.** La grille de temps a d'abord été lancée
pendant que le banc construisait un noyau à 100 kb ; les deux se disputaient le processeur, et
SwiftShader rend sur le processeur. La grille a été arrêtée et relancée seule.

### Ce que ça ne dit pas

**Rien sur un vrai GPU.** Les temps de dessin sont ceux de SwiftShader. **Rien sur TLS ni TCP.**
Le bridage est par requête. **Rien sur des données réelles** : les trois niveaux sont simulés
(S6–S8), et le fichier le dit dans son en-tête comme dans chaque niveau. **Rien au-delà du
premier rendu** : le chargement par octant visible des niveaux fins est écrit dans le lecteur
(`visibleChunks`) et testé, mais aucune navigation n'est encore mesurée — c'est la semaine 10.

### Reproduire

```console
$ cd pipeline && make g3d               # gm12878.g3d (avec aperçu) et gm12878-noyau.g3d
$ make g3d-bench                        # le budget en octets écrits (~10 min la première fois)
$ cd .. && node tools/test-g3d.mjs      # lecture croisée Python → TypeScript, HTTP compris
$ pnpm build && node tools/loadtime.mjs 10    # la grille de temps (~6 min)
```

---

## S8 — Échelle fine : les boucles sortent du mécanisme, le raccord ne tient pas

### Le dispositif

La semaine 7 s'est arrêtée sur un désaccord chiffré : au-delà de 15 Mb, la P(s) du noyau
entier s'aplatit à 0,0253 quand le Hi-C réel continue de décroître. Le modèle avait le *fait*
des territoires sans leur organisation interne. La semaine 8 descend là où cette organisation
naît — entre le kilobase et le mégabase — et y met un **mécanisme** plutôt qu'une
interpolation.

Le mécanisme est l'extrusion de boucles (Fudenberg 2016), en cinq règles dont aucune ne
mentionne le mot « domaine » :

1. un complexe se charge au hasard sur deux monomères voisins ;
2. ses deux jambes s'écartent, une vers la gauche, une vers la droite ;
3. elles ne se dépassent pas entre complexes ;
4. un motif CTCF **qui leur fait face** les arrête, avec sa probabilité d'occupation ;
5. le complexe se décroche au bout d'une durée de vie exponentielle.

Les domaines et les points d'angle sont censés *sortir* de la règle 4. C'est ce que la
validation mesure.

**La règle de convergence, écrite comme une asymétrie.** Un motif *forward* arrête une jambe
qui va vers la gauche ; un motif *reverse* arrête une jambe qui va vers la droite. Un cohésine
chargé entre un `+` en `i` et un `−` en `j > i` bloque donc ses deux jambes et la paire
accumule des contacts. Une paire **divergente** n'arrête rien. Inverser ces deux lignes
produirait des boucles exactement là où Rao 2014 n'en voit pas.

`make fine` — chr7:4 000 000–8 000 000 (la fenêtre qui contient ACTB, le gène que `make query`
affiche depuis la semaine 2), **2 000 monomères de 2 kb**, 20 cohésines de processivité 200 kb,
8 réplicats × 150 instantanés = **1 200 conformations**, 14,2 minutes sur quatre cœurs pour
l'exécution *et* son témoin. Moteur OpenMM, cinq termes, aucun ajusté après coup.

### Rien n'est réglé sur le modèle de noyau, et c'est ce qui rend le raccord falsifiable

Le rayon d'un monomère ne se choisit pas cette semaine : il sort de la loi de la semaine 6 —
volume de chromatine proportionnel aux paires de bases, constante fixée par la fraction
volumique `phi` de ChromEMT — et c'est littéralement la même fonction qui le calcule. À 2 kb,
**23,1 nm**, donc `sigma = 46,3 nm`. Le confinement suit : 2 000 monomères occupent la sphère
de **435 nm** de rayon, c'est-à-dire exactement le volume que le modèle de noyau alloue à 4 Mb
de chromatine.

Si la semaine 8 avait le droit de choisir la taille de ses monomères, on pourrait toujours
faire coïncider les deux modèles, et la coïncidence ne dirait rien.

### Le témoin fait partie du dispositif, pas du commentaire

Poser des barrières et constater qu'il en sort des domaines ne mesure rien : c'est vrai par
construction. Quatre témoins séparent ce qui est expliqué de ce qui est posé.

| Témoin | Ce qu'il retire | n |
|---|---|---:|
| domaines **divergents** | rien — mêmes barrières, orientation inversée | 3 |
| domaines **en tandem** | rien — deux `+` | 3 |
| paires **au hasard**, à distances appariées | l'emplacement | 400 |
| **CTCF inoccupé** | l'arrêt lui-même, positions conservées | tout |

Les trois classes de domaines partagent la même distribution de longueurs : un score plus élevé
aux convergents ne peut donc pas venir de la distance génomique. Le tirage plante 20 domaines —
14 convergents, 3 divergents, 3 tandem — et 19 frontières.

### Les barrières font des domaines

Score d'insulation, fenêtre 100 kb, tolérance ±1 casier de 10 kb :

| | appelées | rappel | précision |
|---|---:|---:|---:|
| frontières à deux sens bloqués (13/19) | 18 | **92 %** (12/13) | 67 % |
| frontières à un seul sens bloqué (6/19) | — | 5/6 | — |
| **témoin CTCF inoccupé** | 9 | **0 %** (0/13) | 0 % |

Douze frontières sur treize avec CTCF, **zéro sans**, alors que rien d'autre n'a changé entre
les deux exécutions — mêmes graines, mêmes positions de sites, même champ de force. Les
domaines viennent bien de l'arrêt de l'extrusion.

### Les points d'angle sont aux paires convergentes

Score de coin sur carte O/E, fond local en anneau, croix exclue :

| classe | n | médiane | moyenne |
|---|---:|---:|---:|
| **convergent** | 13 | **2,56** | 2,41 |
| divergent | 2 | 1,61 | 1,61 |
| tandem | 3 | 1,89 | 2,21 |
| hasard, distances appariées | 400 | **1,00** | 1,05 |
| **convergent, CTCF inoccupé** | 13 | **1,06** | — |

Convergents contre hasard : **2,57×**. Contre divergents : 1,59×. Et contre eux-mêmes sans
CTCF : 2,56 → **1,06**, c'est-à-dire le fond.

Le fond au hasard tombant exactement sur 1,00 est la vérification que le score ne fabrique
rien : sur une carte normalisée O/E avec fond local, une paire quelconque doit valoir 1.

### Un témoin qui n'en est qu'à moitié un, et pourquoi

Les scores domaine par domaine :

```
convergent   3,27  3,22  2,97  2,85  2,65  2,61  2,56  2,48  2,42  2,03  1,69  1,66  0,97
divergent    1,73  1,48
tandem       3,27  1,89  1,47
```

**Un domaine en tandem marque 3,27, autant que le meilleur convergent.** Ce n'est pas du bruit,
et ce n'est pas un défaut du modèle : c'est un défaut du témoin. Un domaine en tandem
(`+ +`) voit bien sa jambe gauche bloquée par son ancre gauche ; sa jambe droite, elle,
traverse son ancre droite — mais si le domaine **suivant** commence par un `−`, ce `−` se
trouve à un monomère de là et l'arrête. La boucle existe, ancrée un cran plus loin, et la
tolérance du score la ramasse.

Le seul témoin propre est donc **divergent**, où aucune des deux ancres n'arrête la jambe qui
l'atteint. Le rapport publié est celui-là, et le fond au hasard sert de référence à grand
effectif.

### Une erreur dans ma propre vérité, sortie par les résultats

La première version étiquetait une frontière comme « bloquante » si l'ancre à sa gauche était
`−` **ou** celle à sa droite `+`. Ça rate les deux autres cas, et le rapport annonçait
18 frontières bloquantes sur 19 — puis imprimait fièrement que la dix-neuvième, le « témoin
non bloquant », avait elle aussi été retrouvée. Un témoin qui contredit son propre énoncé n'en
est pas un.

En réalité **aucune frontière n'est non bloquante** dans ce dispositif : elle porte deux
ancres, et une ancre arrête toujours l'un des deux sens — un `+` les jambes qui vont à gauche,
un `−` celles qui vont à droite. Ce qui distingue deux frontières est leur **force** : 2 quand
les deux ancres sont de sens opposés, donc les deux sens arrêtés ; 1 sinon. 13 sur 19 ici.

### La fenêtre d'insulation : la semaine 3 avait raison, son optimum non

La semaine 3 avait mesuré qu'une fenêtre trop large tue le rappel (100 % à 100 kb contre 42 % à
600 kb sur des TADs de 450 kb) et j'en avais déduit un rapport de ~0,22, donc 50 kb pour des
domaines de 200 kb. Le balayage dit l'inverse à ce bout-là, et la colonne qui tranche est celle
du témoin :

| fenêtre | proéminence | appelées | rappel | précision | témoin : appelées | témoin : rappel |
|---:|---:|---:|---:|---:|---:|---:|
| 30 kb | 0,12 | 28 | 74 % | 50 % | 20 | 11 % |
| 50 kb | 0,05 | 32 | 95 % | 56 % | 33 | 42 % |
| 50 kb | 0,12 | 21 | 95 % | 86 % | 12 | 21 % |
| **100 kb** | **0,05** | **18** | **89 %** | **94 %** | **9** | **5 %** |
| 100 kb | 0,08 | 15 | 79 % | 100 % | 5 | 5 % |
| 100 kb | 0,18 | 2 | 11 % | 100 % | 0 | 0 % |
| 200 kb | 0,05 | 9 | 47 % | 100 % | 0 | 0 % |

Deux bornes, et l'optimum est **entre les deux**. Vers le bas, une fenêtre étroite moyenne
moins de pixels et appelle du bruit : à 50 kb, le témoin **sans aucun CTCF** retrouve encore
42 % des frontières, ce qui disqualifie le réglage quel que soit son rappel. Vers le haut, à
200 kb — l'échelle des domaines elle-même — le rappel s'effondre à 47 %, exactement ce que la
semaine 3 avait mesuré.

La contrainte de la semaine 3 tient donc, mais son optimum ne se transposait pas : ici le
rapport qui marche est 0,5 de la taille des domaines, pas 0,22.

Le seuil est donc choisi **en regardant ce tableau**, ce qui est une sélection ; ce qui la rend
légitime est qu'elle s'appuie sur une référence indépendante — le témoin sans CTCF — et pas sur
le rappel. Le tableau entier est publié, pas seulement sa ligne retenue.

### P(s), et pourquoi il ne faut pas se réjouir de −1,08

| régime | fenêtre | pente |
|---|---|---:|
| boucle | 20–150 kb | −0,69 |
| **domaine** | 150–800 kb | **−1,08** |
| confinement | 1–3 Mb | +0,26 |

Lieberman-Aiden 2009 rapporte `s^−1,08` — le même nombre. Il ne faut rien en conclure, pour
deux raisons qu'il vaut mieux écrire que taire.

D'abord la fenêtre : la pente publiée est ajustée sur **500 kb–7 Mb**, qu'une région de 4 Mb ne
peut pas couvrir. Les deux nombres se ressemblent sans être mesurés pareil, et la semaine 7 a
déjà payé une fois le prix d'une fenêtre d'ajustement prise pour acquise.

Ensuite la convergence : la pente **bouge encore** avec le nombre de réplicats (−0,95 à un
réplicat, −1,08 à huit, cf. ci-dessous). Ce qu'on peut dire honnêtement est qu'elle est entre
−0,9 et −1,1 dans le régime des domaines, ce qui est le bon ordre de grandeur. Pas qu'elle
vaut −1,08.

Le régime « confinement » à +0,26 n'a aucun sens physique — une probabilité de contact ne croît
pas avec la distance. Il mesure la paroi de la sphère et le petit nombre de paires disponibles
au-delà de 1 Mb sur 2 000 monomères. Il est nommé pour ça.

### Le raccord ne tient pas, et c'est le résultat de la semaine

| séparation | modèle fin | noyau | rapport | paires/conf. |
|---:|---:|---:|---:|---:|
| 0,75 Mb | 453 nm | 365 nm | **1,24** | 1 625 |
| 1,50 Mb | 448 nm | 612 nm | 0,73 | 1 250 |
| 2,25 Mb | 450 nm | 806 nm | 0,56 | 875 |
| 3,00 Mb | 430 nm | 958 nm | 0,45 | 500 |
| 3,75 Mb | 444 nm | 1 077 nm | 0,41 | 125 |

Écart le plus grand **2,42×**, contre un facteur 1,25 admis. **Le critère de la semaine 8 n'est
pas atteint.**

Ce que le tableau dit vraiment est une différence de pente, et elle est plus instructive qu'un
rapport :

| | `R(s) ∝ s^ν` |
|---|---:|
| modèle fin (2 kb, extrusion, confiné) | **−0,02** |
| noyau entier (750 kb, semaine 6) | **+0,68** |
| marche aléatoire idéale | 0,50 |
| marche auto-évitante gonflée | 0,59 |
| traçage de chromatine, au-dessus du Mb (Wang 2016, Bintu 2018, Su 2020) | **0,25 – 0,33** |

**Les deux modèles encadrent la mesure, et ils l'encadrent par les deux bouts.** Le noyau de la
semaine 6 gonfle comme une marche auto-évitante : ses billes ne sont liées que par une longueur
*maximale*, et rien ne les retient ensemble au-delà. Le modèle fin, lui, sature complètement,
parce qu'on lui a dit que 4 Mb tiennent dans 435 nm.

Ce dernier point est une réserve sérieuse et elle porte sur le haut du tableau : le confinement
du modèle fin vient de la loi de volume de la semaine 6, donc au-delà d'environ 1 Mb l'accord —
ou le désaccord — devient partiellement circulaire. **La ligne informative est la première** :
à 750 kb, l'échelle d'une seule bille grossière, les deux modèles s'accordent à 1,24×. Ils
divergent ensuite, et la divergence est dans le terme de chaîne, pas dans la taille des billes.

### Ce que ça dit du modèle de noyau, et c'est la même chose que la semaine 7

Une chaîne qui gonfle en `s^0,68` atteint la taille de son territoire beaucoup trop vite.
Passé ce point, tout est à portée de tout, et la probabilité de contact cesse de décroître —
exactement le plateau à 0,0253 au-delà de 15 Mb que la semaine 7 avait mesuré sans savoir d'où
il venait. Les deux observations sont la même : **le terme de chaîne de la semaine 6 n'a aucun
mécanisme de compaction à l'échelle du mégabase.**

C'est une prescription, pas seulement un constat. Il manque au modèle de noyau ce que le modèle
fin a et qu'il n'a pas : quelque chose qui tienne une chaîne compacte sans la confiner à la
main.

### Huit réplicats, est-ce assez ? Non pour tout

| réplicats | conformations | rappel | précision | convergents | divergents | conv./hasard | p. boucle | p. domaine |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 150 | 54 % | 28 % | 2,40 | 1,91 | 2,45× | −0,78 | −0,95 |
| 2 | 300 | 77 % | 40 % | 2,37 | 1,58 | 2,42× | −0,74 | −0,89 |
| 4 | 600 | 85 % | 50 % | 2,58 | 1,44 | 2,57× | −0,72 | −1,04 |
| 8 | 1 200 | 92 % | 67 % | 2,56 | 1,61 | 2,57× | −0,69 | −1,08 |

Le score de coin et le rapport au hasard sont posés dès quatre réplicats. **Le rappel des
frontières et les pentes de P(s), non** : le premier monte encore de 85 à 92 %, les secondes
dérivent régulièrement. Les conclusions qui portent sur les points d'angle sont donc solides ;
celles qui portent sur un chiffre de pente ne le sont pas, et c'est pour ça que la section P(s)
ne cite pas −1,08 comme un résultat.

### La trajectoire a-t-elle oublié son point de départ ?

Comparaison de `R(s)` entre la première et la seconde moitié des instantanés, **à l'intérieur de
chaque réplicat** — couper l'ensemble empilé en deux comparerait des réplicats entre eux, ce
qui est une autre question.

| séparation | 1re moitié | 2e moitié | dérive |
|---:|---:|---:|---:|
| 20 kb | 192 nm | 190 nm | −1,0 % |
| 200 kb | 364 nm | 347 nm | −4,6 % |
| 800 kb | 449 nm | 460 nm | +2,6 % |
| 2,00 Mb | 456 nm | 461 nm | +1,3 % |
| 3,60 Mb | 446 nm | 456 nm | +2,5 % |

Dérive maximale 4,6 % : stationnaire sur la fenêtre d'échantillonnage. C'est `R(s)` aux grandes
séparations qu'il faut regarder, jamais `Rg` — voir le défaut ci-dessous.

### `R(s)` du modèle fin, pour mémoire

| séparation | moyenne | écart-type |
|---:|---:|---:|
| 2 kb | 46 nm | 2 nm |
| 10 kb | 136 nm | 41 nm |
| 34 kb | 236 nm | 93 nm |
| 102 kb | 314 nm | 134 nm |
| 318 kb | 395 nm | 151 nm |
| 978 kb | 449 nm | 152 nm |
| 3,02 Mb | 429 nm | 151 nm |

Les 46 nm à 2 kb sont `sigma` : deux monomères liés sont au contact, par construction. L'ordre
de grandeur aux échelles intermédiaires — ~315 nm à 100 kb, ~395 nm à 320 kb — est celui que
le traçage de chromatine rapporte, mais rien ici ne le compare à une mesure : il n'y en a
aucune dans le dépôt, le réseau est fermé.

### Défauts et erreurs de méthode

**`Rg` est aveugle à la structure interne, et il m'a presque eu.** La première conformation de
départ était une marche aléatoire libre **comprimée** jusqu'à tenir dans la sphère. La courbe
de mise en place était alors parfaitement plate — `Rg` = 7,02 σ au pas 2 000, 7,05 σ au pas
40 000 — et la conclusion évidente était « l'équilibrage est immédiat ». Elle est fausse : une
compression divise **toutes** les distances par le même facteur, donc elle écrase la structure
interne exactement autant que la structure globale, et `Rg` ne peut pas voir la différence.
Mesuré sur `R(s)` à la place, le mode lent apparaît : `R(10)` et `R(100)` se stabilisent en
quelques dizaines de milliers de pas, `R(1000)` dérivait encore à 80 000. La conformation de
départ est maintenant une marche **réfléchie** sur la paroi, qui garde des pas de longueur
unité à toutes les échelles, et un test l'exige.

**Une étude d'équilibrage abandonnée au profit d'une mesure sur les données publiées.** La
première version mesurait la convergence dans une exécution séparée. C'est deux fois moins
bon : ça coûte le double et ça porte sur une trajectoire qu'on ne publie pas. La stationnarité
est maintenant mesurée sur les instantanés de production eux-mêmes (tableau ci-dessus).

**OpenMM ne déplace pas une liaison.** `updateParametersInContext` modifie longueurs et
raideurs, jamais les particules : une liaison par cohésine qu'on ferait suivre l'extrusion lève
`The set of particles in a bond has changed`. Découvert en le faisant. La parade — déclarer
d'emblée toutes les paires du parcours à raideur nulle et n'allumer que celles tenues — a une
conséquence de format : les instantanés d'extrusion gardent un emplacement par cohésine, `-1`
quand elle est décrochée, pour qu'une liaison ne saute jamais d'une paire à une paire sans
rapport.

**J'avais écrit l'inverse de la vraie raison d'exclure la croix.** Le commentaire de
`dot_score` disait qu'elle évite qu'une bande de la carte soit comptée comme un point. C'est
faux, et c'est un test qui l'a montré : une bande seule donne 1,00 à côté d'elle, exclusion ou
pas. La vraie raison est l'inverse — une ancre de boucle *émet* une bande, l'anneau de fond
tomberait dedans, le fond serait surestimé et le point réel **masqué**. C'est pour ça que
HiCCUPS retire la croix.

**Une troncature silencieuse dans le raccord.** La fenêtre de recouvrement était tronquée à la
plus courte des deux courbes. Sur 1 500 monomères, la séparation 1 500 n'existe pas — aucune
paire à cette distance — et le quatrième point disparaissait sans que rien ne le dise. La
fenêtre est maintenant calculée, et une région trop courte pour un seul point lève.

**Des témoins tirés à pile ou face ne sont pas des témoins.** Sur une vingtaine de domaines, un
tirage indépendant donnait couramment 1 divergent contre 5 tandem. Les deux classes alternent
désormais.

**Un fond tiré hors de la zone mesurable.** Les paires au hasard étaient tirées n'importe où, y
compris à moins de six casiers d'un bord où le score rend `NaN`. Sur une petite carte avec de
longues portées, le fond entier pouvait sortir vide. Les tirages sont maintenant contraints à
la zone scorable, et un test le vérifie.

**Une fausse alerte, et d'où elle venait.** Le chemin multiprocessus semblait cassé
(`FileNotFoundError: .../<stdin>`). C'est un artefact de `python - <<EOF` : `spawn` réimporte
`__main__` depuis son chemin, et `<stdin>` n'en a pas. Depuis un vrai fichier, il marche. Rien
à corriger dans le code, tout à corriger dans la façon de le tester.

### Ce que ça ne dit pas

**Rien de biologique.** Les sites CTCF sont **plantés**, pas lus : aucune piste publiée n'est
accessible, même mur réseau qu'aux semaines 1 à 7. Ce qui est validé est la méthode — le
mécanisme produit bien des domaines et des points d'angle, et seulement là où l'orientation le
prévoit. `make fine CTCF=motifs.bed` attend un fichier ; le lecteur est écrit, testé, et refuse
un site sans brin plutôt que d'en inventer un.

**Rien sur la topologie.** La répulsion est volontairement franchissable à 3 kT, parce que la
topoisomérase II laisse passer les brins à l'échelle de temps de l'extrusion. Le nombre
d'enlacements de la région n'est donc pas une prédiction du modèle.

**Rien au-delà du mégabase.** Le confinement sphérique impose la saturation de `R(s)`, et la
région de 4 Mb ne permet aucune mesure au-delà. Les deux dernières lignes du raccord reposent
sur 500 et 125 paires par conformation.

**Rien sur la vitesse réelle des cohésines.** Un pas d'extrusion est un monomère, pas une
seconde. Toutes les quantités rapportées sont des moyennes d'état stationnaire ; aucune ne
demande de convertir ces pas en temps.

### Reproduire

```console
$ cd pipeline && make fine                       # 8 réplicats + témoin, ~15 min sur 4 cœurs
$ make fine CTCF=motifs.bed                      # avec des motifs mesurés, quand il y en aura
$ .venv/bin/python -m geno_pipeline fine --report-only           # relit, sans recalculer
$ .venv/bin/python -m geno_pipeline fine --report-only --sweep   # le tableau fenêtre × seuil
```

Le balayage est dans la commande et pas dans un script à côté, précisément parce que c'est lui
qui justifie le réglage retenu : un tableau qu'on ne peut pas refaire ne justifie rien.

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
