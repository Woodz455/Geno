# Validation

Ce document grandit à chaque semaine qui produit un résultat vérifiable. Il devient le
rapport de validation de la **semaine 21**, qui confronte le modèle aux données
orthogonales — DNA-FISH, chromatin tracing, DamID, SPRITE, GAM.

Règle : on consigne les chiffres, **y compris les mauvais**. Un rapport où tout est vert est
un rapport qui n'a pas cherché.

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
profondeur de séquençage n'est pas une convention, c'est une approximation non chiffrée. La
semaine 6 calibrera alpha sur les données qu'elle utilise, pas sur la littérature.

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
