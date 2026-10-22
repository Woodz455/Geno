# Validation

Ce document grandit à chaque semaine qui produit un résultat vérifiable. Il devient le
rapport de validation de la **semaine 21**, qui confronte le modèle aux données
orthogonales — DNA-FISH, chromatin tracing, DamID, SPRITE, GAM.

Règle : on consigne les chiffres, **y compris les mauvais**. Un rapport où tout est vert est
un rapport qui n'a pas cherché.

---

## S15 — Le cas EPHA4 n'est pas rejoué : le modèle n'isole pas un TAD de 1,8 Mb

**Verdict : non fait.** La page `variant.html` existe — préréglages, délétions tracées à la main,
prédiction 1D en direct, 3D précalculée, surimpressions cliniques exactes — mais le critère de la
semaine, « la suppression d'une frontière rejoue en direct le cas EPHA4 », n'est pas atteint. Ni le
3D ni la prédiction 1D ne montrent l'effet de la frontière que Lupiáñez 2015 mesure, et le 3D
n'isole même pas le type sauvage dès qu'on lui donne assez de réplicats pour le juger. Ce qui suit
consigne les critères fixés d'avance, ce qu'ils donnent, les deux révisions faites après coup — et
pourquoi elles ne changent pas la conclusion.

### Le locus, et ce qui en est reconstruit

Région chr2:218–223 Mb (hg38) à 2 kb, 2 500 monomères. **Mesuré ou annoté** : les gènes (GENCODE
v50 : WNT6, WNT10A, IHH, NHEJ1, EPHA4, PAX3, SGPP2), les 105 pics CTCF d'ENCODE GM12878 dont 96
orientés par MA0139.1 (occupation = 0,95·signal / 90e centile du génome, bornée à [0,05 ; 0,95]),
la boucle de Rao 2014 qui borne le TAD d'EPHA4 (219,735 ↔ 221,575 Mb, levée hg19 → hg38), et la
frontière : la grappe CTCF orientée à ± 20 kb du TSS d'EPHA4 (221,572–221,574 Mb). Ce sont des
CTCF de **lignée lymphoblastoïde**, pas de bourgeon de membre.

**Reconstruit depuis le texte** — la Table S1 de Lupiáñez 2015 (points de cassure des patients)
n'a pas pu être lue : PMC derrière un reCAPTCHA, Europe PMC en 403/500, aucun fichier
supplémentaire exposé, aucune entrée ClinVar/dbVar. Le texte intégral (API BioC) contraint :
délétions de 1,75 à 1,9 Mb emportant EPHA4 et la frontière EPHA4/PAX3 ; grappe de trois enhancers
dans 30 kb, au contact du promoteur d'EPHA4 « across 1.66 Mb » (distance mesurée chez la souris,
appliquée telle quelle) ; DelBS, le témoin murin, garde la frontière. D'où non pas **un** variant
mais **la famille** que le texte autorise : 34 DEL-B (longueurs 1,75/1,80/1,85/1,90 Mb,
extrémité centromérique au pas de 50 kb au-delà des enhancers, extrémité télomérique entre la
frontière et la fin 3' de PAX3), et pour chacun deux témoins — **DelBS** (même extrémité
centromérique, s'arrête 2 kb avant la frontière ; 10 distincts) et **« frontière épargnée »**
(mêmes points de cassure, les sites CTCF de la frontière gardés au point de recollement : seule la
frontière diffère, pas la distance). Chaque coordonnée reconstruite est étiquetée comme telle dans
le fichier exporté (`evidence`).

### Les critères, fixés avant de regarder les mutants

1. **Le type sauvage isole** : en 3D, enhancers → promoteur d'EPHA4 ≥ 3 × enhancers → PAX3.
2. **La frontière décide** : sur le DEL-B central (220,20–222,00 Mb), enhancers → PAX3 ≥ 1,5 ×
   chacun des deux témoins, en 3D ; et en 1D pour **chacun** des 34 DEL-B.
3. **Les autres frontières tiennent** : DEL-B / sauvage, enhancers → IHH, dans [0,67 ; 1,5].
4. **La prédiction 1D suit le 3D** : Pearson des log du 4C virtuel depuis PAX3 (casiers de 20 kb)
   ≥ 0,9 sur les trois mutants, et le rapport du critère 2 du même côté de 1,5.

Plus, pour la page : geste → lecture affichée < 1 s au 95e centile, et les rapports 2 et 3 lus
dans la page.

### Plan pré-enregistré : un critère sur quatre

OpenMM, 8 réplicats × 150 instantanés par allèle, mêmes graines ; chaque allèle mutant est la
région moins la délétion (DEL-B 1 600 monomères, DelBS 1 815). Loi 1D `P = min(1, (d/d0)^−α + c)`
sur la distance effective (plus court chemin le long de la chaîne, boucles tenues comprises),
calée sur le 3D du type sauvage par déviance de Poisson : α = 1,253, d0 = 0,686, c = 0,00082,
r(log, 20 kb) = 0,651.

| | sauvage EPHA4/PAX3 | DEL-B / DelBS | DEL-B / épargnée | DEL-B / sauvage | IHH |
|---|---|---|---|---|---|
| 3D | **4,55** | 1,05 | 1,78 | 5,65 | 4,87 |
| 1D | 1,17 | 1,45 | 1,16 | 2,23 | 1,03 |
| 1D, famille (min sur 34) | — | 1,03 | 0,89 | — | — |

PASS 1 · FAIL 2 (3D) · FAIL 2 (famille 1D) · FAIL 3 · FAIL 4 (r 0,64–0,73).

### Deux défauts du modèle, trouvés en regardant pourquoi

**Le confinement dépend de la longueur.** La sphère est réglée sur la quantité de chromatine ; une
chaîne plus courte y a un plancher de contacts lointains plus haut — 3D, séparations de 500 à 800
monomères : 0,0012 pour le type sauvage, 0,0023 pour DEL-B. DEL-B et DelBS n'étaient donc pas
comparés à confinement égal, et l'IHH « monte » ×4,9 pour cette raison.

**Le plancher de la loi 1D est cet artefact.** c = 0,0008 reproduit le plateau que la sphère
impose au-delà d'un mégabase — déjà relevé à la semaine 7 (« la sphère décide : artefact ») — et il
dépasse les contacts lointains eux-mêmes : en 1D, le type sauvage n'isolait plus (1,17).

### Plan révisé (après coup, seuils inchangés)

- **Fenêtre constante** : chaque allèle garde 5 Mb (2 500 monomères), la délétion faisant entrer la
  séquence en aval (`Locus.sites(window=True)`). Le type sauvage ne change pas.
- **Loi sans plancher** : c = 0, calée sur les paires à moins de 1 Mb : α = 1,128, d0 = 0,535,
  r(log, 20 kb) = 0,801.

À 8 réplicats, il échoue aussi. Mais le détail par réplicat montre que **ni le plan pré-enregistré
ni la révision ne pouvaient trancher** : chaque contact enhancer → promoteur à ~1 Mb tient à un ou
deux réplicats (erreur type de 50 à 100 % de la moyenne). Un réplicat, c'est 150 000 pas de mise en
place puis 30 000 d'échantillonnage : ses 150 instantanés partagent à peu près **une** conformation
à grande échelle. Le tirage indépendant est le réplicat, pas l'instantané. IC à 95 % par
rééchantillonnage des réplicats, plan pré-enregistré : sauvage [0,94 ; 21,8], DEL-B/DelBS
[0,28 ; 6,1], DEL-B/épargnée [0,40 ; 99,9], IHH [0,47 ; 9,1]. Le « PASS » du critère 1 était un
tirage.

D'où `geno sv extend` : des réplicats ajoutés au plan révisé, 8 par tour, les quatre allèles
ensemble, graines à la suite. Visé : 48 par allèle. **Arrêté à 24** par la limite de durée des
tâches de fond de l'environnement (2 h) ; les tours étant enregistrés un à un, les quatre allèles
sont à 24, cohérents. Résultat, plan révisé, 24 réplicats :

| | sauvage EPHA4/PAX3 | DEL-B / DelBS | DEL-B / épargnée | DEL-B / sauvage | IHH |
|---|---|---|---|---|---|
| 3D | 0,90 | 1,02 | 2,66 | 1,03 | 1,46 |
| IC 95 % (réplicats) | [0,36 ; 2,17] | [0,37 ; 4,90] | [1,03 ; 9,56] | [0,44 ; 2,49] | [0,23 ; 5,89] |
| 1D | 1,67 | 1,86 | 1,18 | 5,10 | 0,93 |
| 1D, famille (min / médiane sur 34) | — | 1,05 / 1,59 | 0,87 / 1,03 | — | — |

FAIL 1 · FAIL 2 (3D) · FAIL 2 (famille 1D) · PASS 3 · FAIL 4 (r 0,67–0,81 ; le 3D contre lui-même,
moitiés de réplicats : 0,33–0,62).

### Ce que ça dit

- **Le 3D n'isole pas le TAD d'EPHA4.** À 24 réplicats, l'enhancer touche PAX3 autant qu'EPHA4
  (0,0020 contre 0,0018). À 1,7–2,4 Mb, les contacts du modèle sont le mélange du confinement, pas
  une structure de domaine : avec des cohésines de 200 kb de processivité, la boucle de 1,84 Mb que
  Rao 2014 mesure dans GM12878 ne se forme presque jamais. Sans type sauvage isolé, il n'y a pas
  d'isolement à supprimer.
- **Le seul rapport 3D dont l'IC exclut 1** est DEL-B / « frontière épargnée » (2,66 [1,03 ; 9,56]),
  à la limite.
- **La prédiction 1D ne voit pas la frontière.** Contre « frontière épargnée », le rapport médian
  sur la famille est 1,03. Contre DelBS il passe souvent 1,5, mais c'est la **distance** : DelBS
  enlève moins de séquence et laisse les enhancers plus loin de PAX3. Le témoin du papier mélange
  les deux effets ; le témoin in silico les sépare, et ne laisse rien.
- **Exploratoire, pas un résultat** : la processivité ne sauve pas le 1D. Balayée de 200 kb à 2 Mb
  (16 exécutions, loi inchangée), le rapport sauvage EPHA4/PAX3 reste entre 1,22 et 1,84, DEL-B /
  épargnée entre 1,01 et 1,51, sans tendance hors du bruit.

### Ce qui marche, et est vérifié

- **L'extrusion 1D en TypeScript** (`apps/variant/extrusion.ts`) est celle de Python : densité de
  pieds Pearson 0,961 (Python contre Python : 0,959), occupation 1,000 des deux côtés, lignes de
  contact Pearson des log 0,9991, et contacts enhancer → EPHA4/PAX3/IHH indiscernables sur 16 lots
  (Welch |t| ≤ 1,05). **Un biais trouvé en route** : avec le générateur mulberry32, enhancer → PAX3
  sortait 3,6 % sous Python (Welch p = 0,004 sur 20 lots) alors que les statistiques de boucles
  concordaient ; passé à sfc32, plus d'écart (p 0,07–0,49). Les premiers seuils de ce test (Pearson
  0,97, écart 10 %) avaient été posés avant de connaître le bruit d'un lot unique ; remplacés par
  « pas pire que Python contre Python » et Welch sur 16 lots.
- **La loi** : la première version, moindres carrés sur le log des seules paires non nulles (20 %),
  sélectionnait les paires chanceuses et dégénérait (α = 0,23, d0 → 0, r = 0,61) ; remplacée par la
  déviance de Poisson sur toutes les paires, zéros compris.
- **Surimpression clinique** dans la plongée (`dive.html`, menu « surimpression ») et la page du
  variant : P/LP ClinVar, associations GWAS, LOEUF gnomAD minimal, par bille. Relue dans le
  navigateur contre le lecteur Python du fichier de fiches : aperçu 2 032/2 032 billes égales pour
  chacune des trois couches, noyau 8 082/8 082 ; peinte en 21–29 ms. Sur le locus (809 monomères
  porteurs), recomptée depuis les fichiers bruts : ClinVar 809/809, GWAS 809/809 ; témoin décalé
  d'un monomère : 677 et 215.
- **La page** : 4 préréglages 3D (fondu sauvage → mutant, alignés par Kabsch sur les monomères
  communs ; les structures exportées sont celles du plan révisé, 24 réplicats, médoïde), délétions
  tracées sur la piste, deux 4C virtuels, et pour toute délétion qui emporte la frontière, le même
  allèle « frontière épargnée » calculé à côté.

### La page, mesurée (`node tools/variant.mjs 12`)

28 gestes CDP : 4 préréglages, 12 DEL-B tirés dans la famille et tracés à la souris, leurs 12 DelBS.

| | mesuré | critère |
|---|---|---|
| geste → lecture affichée | p50 1 828 ms · p95 2 800 ms · max 3 000 ms | p95 < 1 s — **échec** |
| DEL-B / frontière épargnée | min 0,80 · médiane 1,05 | ≥ 1,5 — **échec** |
| DEL-B / DelBS | min 1,12 · médiane 1,83 | ≥ 1,5 — **échec** (min) |
| DEL-B / sauvage, IHH | 0,93 – 1,10 | [0,67 ; 1,5] — tient |

La latence : une délétion qui emporte la frontière coûte deux prédictions (l'allèle et son témoin),
16 extrusions et ~20 lignes de plus court chemin chacune, sur un conteneur à 4 vCPU sans GPU.

### Ce que ça ne dit pas

- Rien sur la biologie d'EPHA4 : le modèle échoue, le cas de Lupiáñez reste ce qu'il est.
- Les points de cassure sont reconstruits ; les CTCF sont lymphoblastoïdes. Aucun des deux n'explique
  l'échec du type sauvage, qui ne dépend d'aucune délétion.
- SOX9 (Franke 2016) n'est pas commencé.

### Pistes, non tentées

Faire tenir un domaine de 1,8 Mb demande soit des cohésines de processivité locale bien plus
grande **avec** un 3D assez échantillonné pour le juger (≥ 48 réplicats par allèle : ~4 h de calcul
par allèle ici), soit un modèle contraint par le Hi-C mesuré du locus — le `.hic` de Rao 2014 reste
`unverified` (DATA_SOURCES § 8 ter). L'un ou l'autre est un choix de conception, pas un réglage.

### Reproduire

```console
$ cd pipeline
$ geno sv run                           # plan pré-enregistré : 4 allèles × 8 réplicats (~50 min)
$ geno sv run --window --only delb,delbs,delb-keep    # plan révisé (~45 min)
$ geno sv extend --replicates 24        # réplicats ajoutés, 8 par tour
$ geno sv check                         # critères, plan pré-enregistré (+ révisé si présent)
$ geno sv revised                       # plan révisé seul (loi sans plancher calée au besoin)
$ geno sv export                        # data/g3d/epha4.variant.json + recomptage clinique
$ cd .. && node tools/build.mjs && node tools/test-variant.mjs && node tools/overlay.mjs
$ node tools/variant.mjs 12
```

---

## S14 — Trois vues, une sélection : 90/90 vues ailleurs, zéro rebond

### Le critère, traduit

**Fait quand : une sélection faite dans l'une des trois vues apparaît dans les deux autres, sans
rebond.** Fixé avant de mesurer :

1. **« faite dans l'une des trois vues »** : par un vrai geste (événements souris CDP) dans la
   page `vues.html` — un pixel de la 3D, un élément d'une piste igv.js, un pixel de la carte
   HiGlass —, 30 fois par vue, à des endroits et des échelles variés ;
2. **« apparaît dans les deux autres »** : relu dans chacune, pas supposé — la 3D : des pixels de
   billes surlignées dans son tampon d'identifiants, de la teinte de la surbrillance ; igv : une
   région d'intérêt par intervalle, dessinée dans la vue ; HiGlass : les bandes de chaque intervalle
   dessinées sur la carte, dans la vue ;
3. **« sans rebond »** : **une** publication acceptée, venue de la vue du geste, **aucune**
   refusée (un écho ou une égalité refusés seraient un rebond rattrapé : on veut zéro tentative),
   chacune des deux autres vues appliquant la sélection **une** fois — jugé après le repos des
   trois vues et 400 ms de plus, pour laisser à un rebond tardif le temps de se produire. Avec un
   **témoin** : la conception que la feuille de route écarte, des vues synchronisées deux à deux,
   doit, elle, rebondir.

### Ce qui est montré, et ce que ça vaut

| Vue | Ce qu'elle montre | Évidence |
|---|---|---|
| 3D | la plongée des semaines 10 à 13, fiche comprise | modèle simulé, séquence hg38 |
| igv.js 3.8 | séquence hg38 (le `.2bit` verrouillé, lu par plages) ; gènes GENCODE v50 ; pics CTCF de GM12878, brin = motif ; boucles Rao 2014 levées en hg38 | annotation, mesuré, mesuré |
| | compartiments : l'E1 du modèle | simulé |
| HiGlass 2.3 | contacts du **modèle** : noyau (24 Mb → 750 kb), régions fines (512 kb → 2 kb) | simulé |

**La carte de HiGlass est celle du modèle, pas du Hi-C mesuré** — le choix de la 3D, qui reste
simulée tant que le rattrapage S1–S8 n'est pas fait ; la page le dit dans le cartouche de la carte.
Les deux cartes se superposent, la fine par-dessus (NaN hors de ses régions, donc transparente), et
leur frontière se voit : c'est le raccord S8, qui échoue, vu cette fois en contacts.

Aucune requête ne sort de la page : igv lit la référence et ses pistes sur le serveur local
(`/ref/hg38.2bit`, 835 Mo, par plages ; `/data/igv/`), HiGlass ses tuiles par un data fetcher
maison. La mesure relit la liste des ressources de la page : **148 requêtes, 0 hors du serveur
local** — après une correction : igv allait chercher sa liste de génomes sur `igv.org` à
l'ouverture (`loadDefaultGenomes: false`).

### La sélection : un store, deux défenses, une règle

Une sélection est **un ou deux intervalles haploïdes** (chrom, début, fin) : une bille, un
monomère, une base, un gène, un pic — ou un pixel de la carte, deux loci en contact. Haploïde,
parce qu'igv et HiGlass le sont ; la 3D surligne les deux homologues. Le store (`apps/linked/store.ts`)
la distribue aux deux autres vues, jamais à sa source, qui la montre d'elle-même. Deux défenses,
chacune journalisée :

- **égalité** : une sélection égale à la courante est refusée et n'est pas redistribuée ;
- **garde d'application** : pendant qu'une vue applique ce qu'elle reçoit, ce qu'elle publierait
  est refusé comme écho.

Le store seul est vérifié sous Node (`tools/test-store.mjs`, 7 assertions) : une sélection va aux
deux autres vues et pas à sa source ; la même, ou la même paire dans l'autre ordre, est refusée et
pas redistribuée ; une vue qui republie, élargie, ce qu'elle reçoit voit son écho refusé.

Et une règle, qui rend les deux défenses inutiles en pratique : **naviguer n'est pas
sélectionner**. Le déplacement d'une vue — igv `locuschange`, HiGlass `location`, l'orbite de la
3D — n'est jamais publié. Une vue qui reçoit une sélection hors de son champ s'y déplace : igv
ouvre une fenêtre par intervalle (côte à côte pour deux loci), HiGlass centre le pixel du contact,
la 3D vise les billes surlignées, en gardant son orientation.

**Dans la 3D**, une sélection : ses billes en jaune Okabe–Ito (absent des palettes de territoire et
de compartiment), les autres assombries à 40 %, le reste de leurs copies tramé à 30 %, les autres
copies à 15 % — sinon une bille de chr1 restait enfouie dans son propre territoire. Une copie qui
porte la sélection reste visible à toutes les échelles, contre l'effacement de la plongée. À
l'échelle de la fibre, ses monomères. Couleurs et groupes réécrits dans les tampons existants
(`ImpostorRenderer.recolor`), sans en créer.

### La carte de contacts : `gm12878.contacts`

`make linked` : les 200 repliements du noyau, comptés comme à la semaine 7 (centres à moins de
1,5 × (r_i + r_j)), sur la carte haploïde de 4 041 bins (3 875 404 pixels non nuls) ; les deux
régions fines, 2 000 monomères × 1 200 conformations chacune, comme à la semaine 8. **La grille de
HiGlass est uniforme en coordonnées absolues, les billes non** (747 à 753 kb, et un chromosome ne
commence pas sur un multiple de la résolution) : chaque bille est répartie sur les pixels qu'elle
recouvre, au prorata de sa longueur. Tuiles de 256 × 256 en float16, deflate, une empreinte chacune :
**21,1 Mo**, 325 tuiles (187 du noyau, 138 des régions fines ; 43,4 Mo en float32, 11 bits suffisent
à une couleur). Seules les tuiles x ≤ y sont écrites ; HiGlass produit les autres par miroir — et
l'orientation (colonnes = x, lignes = y) a été établie sur un motif asymétrique avant d'écrire un
octet.

| Contrôle (`geno linked check`) | Résultat | Témoin |
|---|---|---|
| conservation : somme des tuiles, miroir compris, contre la matrice de bins | écart relatif **≤ 1,5·10⁻⁵** sur les 15 résolutions | — |
| 200 pixels de 750 kb recomposés par distances directes (ni KD-tree ni matrice creuse) | **200/200** à 10⁻³ près | décalés d'un pixel : 39/200 (des zéros trans) |
| 200 pixels de 2 kb des régions fines, idem | **200/200** | décalés : 4/200 |
| lecture croisée Python → TypeScript (`tools/test-contacts.mjs`) | **80/80 tuiles au bit près**, 80/80 pixels, aller-retour chr17:7 673 802 → absolu | tuile corrompue d'un octet : refusée |

Pistes d'igv, écrites par la même commande depuis les sources verrouillées : gènes (1,04 Mo),
CTCF orientés (443 ko), boucles (128 ko ; 9 448 boucles, 22 perdues au levage), compartiments du
modèle (42 ko), cytobandes, tailles — chacune lue en entier, sans index.

### La mesure : 90/90

`tools/linked.mjs 30`, 1 440 × 860, SwiftShader :

| Geste dans | n | publiée une fois | sans rebond | vue dans les deux autres | geste → publication | → 3D | → igv | → HiGlass |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 3D (15 billes de 3 Mb, 15 monomères ou sites de la fibre) | 30 | 30 | 30 | **30** | 6 / 11 ms | — | 13 / 24 ms | 1 / 2 ms |
| igv (gènes et pics CTCF, 6 loci, 12 éléments distincts) | 30 | 30 | 30 | **30** | 501 / 502 ms | 13 / 14 ms | — | 3 / 5 ms |
| HiGlass (pixels de 6 Mb, 750 kb, 8 kb, 2 kb) | 30 | 30 | 30 | **30** | 0 / 1 ms | 31 / 353 ms | 23 / 46 ms | — |

(latences p50 / p95 ; « → vue » : de la publication à l'application — pour la 3D, l'image peinte
qui la montre ; une vue qui doit se déplacer y ajoute son déplacement.) **90/90** : chaque
sélection, publiée une fois, appliquée une fois par chacune des deux autres vues, vue dans chacune ;
**zéro** publication refusée. 29 des 90 sélections portent deux intervalles (pixels hors
diagonale), 61 un seul.

- **igv attend 500 ms** avant de traiter un clic simple : c'est son délai de double-clic
  (`doubleClickDelay`, le double-clic y zoome). Le geste → publication d'igv, c'est lui ; laissé
  par défaut.
- **HiGlass → 3D, 353 ms au p95** : la 3D peint sa première image après que HiGlass a décodé et
  dessiné ses tuiles sur le même fil, sous SwiftShader. Sur un vrai GPU, à mesurer.

**La navigation ne publie rien** : 15 gestes — glisser la piste d'igv, molette dans HiGlass,
tourner la 3D — : **0** publication.

**Le témoin rebondit.** La même page sans store (`?naive`) : igv `locuschange` déplace HiGlass et
la 3D, HiGlass `location` déplace igv et la 3D. Un seul geste (glisser igv) : **190 messages dans
la première demi-seconde**, 246 à 8 s, **52 loci distincts** — chaque vue arrondit ou ajuste la
fenêtre qu'on lui donne (le format d'igv, le rapport d'aspect de HiGlass), la renvoie, et la
fenêtre dérive. La boucle ne s'arrête pas ; on a coupé au bout de 8 s.

Quatre défauts trouvés en chemin, chacun avant que la mesure ne passe :

- **la carte du noyau ne s'affichait pas** sous celle des régions fines : le fond par défaut d'une
  heatmap HiGlass (gris clair) la couvrait — fond transparent ;
- **la surimpression de HiGlass restait sur l'ancienne vue** : `zoomTo` sans transition ne prévient
  pas ses écouteurs, et `getLocation()` lit des échelles qu'il ne remplace pas tout de suite —
  150 ms de transition, et la vue tenue à jour par l'événement `location` ;
- **deux clics igv sur trois ne sélectionnaient rien** : ils tombaient sur le cartouche du nom de
  piste, qu'igv dessine par-dessus la piste ;
- **une sélection venue d'ailleurs restait invisible dans la 3D** quand la caméra était à la fibre
  de chr17 : la plongée efface toutes les copies sauf chr17:a sous 12,5 µm de champ — une copie qui
  porte la sélection l'emporte désormais.

Et deux défauts de l'outil, corrigés : pixi.js réémet un `pointerup` synthétique 40 ms après le
vrai, qui écrasait l'instant du geste (seuls les événements de confiance comptent) ; et les pixels
de la 3D, tirés une fois avant le premier clic, tombaient ensuite sur des billes que la sélection
avait effacées.

### Non-régression

| Mesure | Semaine 13 | Semaine 14 |
|---|---:|---:|
| `make test` | 276 | **283** (+ 7 vues liées) |
| `pnpm verify` | 82 | **82** |
| lecture croisée `.g3d` · fiches | 46 · 6/6 | **46 · 6/6** ; contacts **5/5** ; store **7/7** |
| continuité, 600 images | 0 saut, transitions 0,94–1,00×, pire rapport 1,47 | **identique** (témoin : 3 sauts, 7,41×) |
| lisibilité S11 : ancres, anneau, arc, corps, flèches | 269 · 131 · 1 662 · 2 961 · 85 % | **identiques** ; raccords 45 % et 11 % |
| fidélité des bases | 90/90 (témoin 23 %) | **90/90** (témoin 23 %) |
| temps réel Slow 4G, images en attente | 0 (témoin 1) | **0** (témoin 1 à 9 s, 3 à 40 s) |
| clic → fiche peinte, chaud | p95 43,9 ms, 504/504 | p95 **26,9 ms**, 114/114 (120 demandés) |
| premier rendu S9, aperçu, préfixe 16 Kio, Slow 4G | 1 478 ms | **1 420 ms** (marge 80 ms) |
| page du premier rendu | 15,8 Kio | 16,1 Kio (+ `recolor`) |

La surbrillance ne change rien à la plongée sans sélection : mêmes couleurs, mêmes fondus. Le
premier rendu S9 **gagne** 58 ms en ajoutant 0,3 Kio à la page : c'est l'environnement qui varie
(± 70 ms d'un jour à l'autre, § S13), pas le code — la marge reste celle d'une mesure qui flotte.


### Ce que ça ne dit pas

- **Les contacts sont ceux du modèle.** Une carte de Hi-C mesuré de GM12878 (Rao 2014, `.hic` de
  dizaines de Go, ou un `.mcool` 4DN) n'est pas dans la page : le manifeste ne l'a pas verrouillée,
  et elle appartient au rattrapage S1–S8.
- **Le témoin n'est qu'une conception naïve parmi d'autres** : il montre que la mesure voit un
  rebond, pas que toute synchronisation deux-à-deux en produit.
- **SwiftShader, encore** : les latences sont des bornes hautes ; aucune vue liée n'a tourné sur un
  vrai GPU.
- **Une page lourde** : `vues.html` pèse 4,3 Mo (igv.js, HiGlass, React, pixi.js) — ce n'est pas la
  page du premier rendu, qui reste `first.html`, mais elle se charge en plusieurs secondes sous
  Slow 4G, non mesurées.
- **Un avertissement** reste dans la console : la dépréciation d'une API de pixi.js 6, émise une
  fois par HiGlass à la création de son rendu.

### Reproduire

```
cd pipeline && make linked       # gm12878.contacts + pistes igv, puis geno linked check
make test
cd .. && pnpm install && pnpm build
node tools/test-contacts.mjs     # lecture croisée des contacts, Python → TypeScript
node tools/test-store.mjs        # le store seul : distribution, égalité, écho
node tools/linked.mjs 30         # 3 × 30 sélections, navigation, témoin
node tools/serve.mjs             # puis http://127.0.0.1:8080/vues.html
```

---

## S13 — Un clic, une fiche : complète, exacte, peinte en 44 ms (p95)

### Le critère, traduit

**Fait quand : un clic sur n'importe quelle bille produit la fiche complète en moins de 100 ms.**
Trois mots, trois mesures, fixées avant de mesurer :

1. **« n'importe quelle bille »** : à chaque étape de la plongée — aperçu de 3 Mb, bille de
   750 kb, monomère de 2 kb, nucléosome, ADN, base de l'hélice — un pixel cliqué désigne le bon
   objet et la bonne position génomique, aussi fine que l'objet le permet ;
2. **« la fiche complète »** : pour **toutes** les billes de tous les niveaux, chaque champ de la
   feuille de route a une valeur ou une raison écrite de ne pas en avoir — aucun champ muet — et
   les valeurs sont exactes, recomposées par un autre chemin ;
3. **« en moins de 100 ms »** : du clic (événement souris réel) à l'image où la fiche est peinte,
   p95, dans la page interactive — avec un témoin sans préchargement qui, lui, attend.

### Les sources : réelles, sauf ce qui n'existe pas

Le réseau ouvert à la semaine 12 a permis de résoudre, télécharger et verrouiller 24 sources
(DATA_SOURCES § 8) ; leur md5 est celui que publie le dépôt partout où il en publie un (ENCODE ×9,
4DN ×4, ClinVar). La structure 3D, elle, reste celle du modèle : **la fiche annote une coordonnée
avec des données réelles, elle ne refait pas le modèle** — le rattrapage S1–S8 reste une étape à
part. Chaque ligne de la fiche porte son badge d'évidence, couleur permanente :

| Évidence | Champs |
|---|---|
| mesuré dans GM12878 | sous-compartiment, domaines et ancres (Rao 2014), cCRE, ChromHMM, CTCF, timing de réplication |
| mesuré dans une autre lignée | lamina LMNB1 et speckles SON (K562), nucléole (HFFc6), expression (LCL de GTEx) |
| simulé | compartiment E1, position radiale, LAD du modèle ; CTCF plantés et boucles aux niveaux fins |
| référence | position hg38 et T2T, cytobande, composition, bases |
| annotation | gènes GENCODE v50, contrainte gnomAD, ClinVar, GWAS |

Deux choix dits d'avance. **LAD, speckles et NAD n'ont aucune mesure en GM12878** : la fiche
montre la lignée la plus proche, étiquetée comme telle (décision posée en question). **Rao 2014 et
le Repli-seq sont en hg19** : levés par la chaîne Ensembl, ce qui ne se lève pas est compté et
perdu (47 sous-compartiments sur 3 849, 38 domaines sur 9 274, 22 boucles sur 9 448).

Le CTCF orienté vient d'un scan de MA0139.1 (JASPAR, seuil relatif 0,80) sur la séquence hg38 de
chaque pic IDR : 39 276 pics orientés sur 43 865 — 19 766 `+`, 19 510 `−`, 4 589 sans motif
reconnu, dits tels.

### Le fichier de fiches

`make fiche` : une fiche par intervalle distinct de chaque niveau — 1 016 pour les 2 032 billes de
l'aperçu, 4 041 pour les 8 082 du noyau, 2 000 + 2 000 aux régions fines —, plus une colonne par
copie pour ce que dit le modèle. **6,50 Mo** au total : tables partagées 2,11 Mo (gènes 1,70,
maladies 0,27, traits 0,14), aperçu 0,86, noyau 2,12, régions fines 0,51 et 0,87. C'est **au-dessus
du budget de 2,5 Mo** que je m'étais fixé ; la première version faisait 12,15 Mo, avec les gènes
recopiés dans chaque fiche (3,2 Mo à eux seuls). Ce que ça coûte se lit au préchargement, plus bas.

### « La fiche complète »

`geno fiche check`, sur le fichier écrit :

**Complétude — toutes les billes.** **0 champ muet** sur 2 032 + 8 082 + 2 000 + 2 000 = 14 114
fiches ; chaque absence porte sa raison. Les plus fréquentes : « les bases se lisent au zoom
maximal » (toute bille hors des régions fines), contrainte gnomAD sur des gènes non codants, X
inactif sans compartiment, trous d'assemblage sans ChromHMM ni Repli-seq. Témoin : une fiche dont
on retire `clinvar` est vue incomplète.

**Exactitude — un autre chemin.** 200 fiches tirées au hasard, recomposées **sans** le code qui
les a écrites : gènes, cCRE, CTCF et cytobandes par le magasin d'intervalles de la semaine 2
construit sur les fichiers verrouillés ; ClinVar et GWAS par un balayage ligne à ligne des
fichiers bruts.

| | gènes | cCRE | CTCF | cytobandes | ClinVar (comptes par classe) | GWAS |
|---|---:|---:|---:|---:|---:|---:|
| fiches | **200/200** | **200/200** | **200/200** | **200/200** | **200/200** | **200/200** |
| témoin, intervalles décalés d'une bille | 92 | 45 | 100 | 131 | 73 | 58 |

**Faits connus**, écrits avant de lire le fichier : TP53 c.818 en 17p13.1 ; son monomère porte
TP53, brin −, MANE Select ; pLI 0,998 et LOEUF 0,449 ; rs28934576 **C>T** (R273H) pathogène à
chr17:7 673 802 (5 variants ClinVar à cette base) ; la position T2T sur chr17 ; ACTB en 7p22.1,
8 062 TPM en LCL. Tous vérifiés.

**La position T2T, sans UCSC.** UCSC ne passait plus le TLS ; la chaîne hg38 → T2T est donc
contrôlée par les deux séquences : sur 10 000 positions tirées au hasard, 9 773 se lèvent ;
**9 558 ont le même 21-mer** dans `hs1.2bit` (97,8 % — deux assemblages qui diffèrent d'environ une
SNV par kb) et 9 764 la même base (99,91 %). Témoin décalé de 1 kb : 0 21-mer, 25,7 % de bases
— le hasard. TP53 c.818 tombe sur chr17:7 577 926 de T2T-CHM13.

**Lecture croisée** (`node tools/test-fiche.mjs`) : 160 fiches tirées dans les quatre niveaux, relues
par le lecteur TypeScript, sont celles de Python **champ pour champ** ; 1 000 positions des régions
fines levées par la chaîne embarquée dans le fichier sont celles de la chaîne complète ; un bloc
corrompu est refusé.

### « N'importe quelle bille » : le picking et le mapping inverse

La cible d'identifiants passe de R32UI à **RGBA32UI**, relue en `RGBA_INTEGER` — le seul format
entier que WebGL2 garantit pour `readPixels` — ; R porte l'identifiant, G le type de rendu.
Suite headless `pick` (`pnpm verify`), 12 assertions :

- chaque type de rendu écrit son identifiant et son type — sphère jusqu'à 4 294 967 294, tube
  (plage + identifiant local), ellipsoïde, hélice — ; le fond n'est rien ; une lecture coûte
  0,35 ms ;
- **un défaut de la semaine 11 trouvé et corrigé** : l'ADN de 1 nm entre les nucléosomes
  réutilisait la plage d'identifiants des cœurs. Sur une vue de deux monomères, 48 des 343
  identifiants d'ADN visibles se lisaient comme des nucléosomes ; avec sa propre plage : 0 ;
- **aller-retour base → pixel → base** sur 500 bases d'ADN visibles : 482 retrouvées à la
  résolution du pixel (≥ 1 pb ; vue de biais, un pixel couvre plusieurs bases), 464 à ± 1 pb,
  404 exactes. La base est celle dont **l'axe passe au plus près de la ligne de visée** ; deux
  témoins montrent pourquoi : l'identifiant seul (le segment de 6 pb) n'en donne que 32 exactes,
  et le point de surface touché — ce qu'utilisait la première version — glisse le long de l'axe
  de r / tan θ : 241 à ± 1 pb ;
- **un second défaut, plus ancien** : `Picker` arrondissait la coordonnée (`round`) et lisait le
  pixel voisin dès que le clic tombait dans la moitié droite d'un pixel. En `floor`, l'aller-retour
  passe de 252 à 404 bases exactes.

Ce que vaut la position, par objet, et que la fiche écrit : une bille, son intervalle ; un monomère
de la fibre, ses 2 kb — la plus petite unité du polymère, pas de position inventée dedans ; un
nucléosome, la base enroulée la plus proche de la ligne de visée ; l'ADN de liaison, la base de
l'axe ; un barreau de l'hélice, sa base, exacte. Au niveau de la base, la fiche du monomère se
complète de ce qui vaut pour elle : lettre et contexte (référence), position T2T exacte, variants
ClinVar à cette base, état ChromHMM et cCRE qui la couvrent.

### « En moins de 100 ms »

`node tools/pick.mjs 500` : la page interactive, de vrais événements souris (CDP
`Input.dispatchMouseEvent`), 1 280 × 720, SwiftShader. À chaque étape de la plongée, des pixels qui
portent un objet tirés à graine fixe, un clic à la fois ; la page chronomètre la lecture du pixel
et le mapping inverse, la fiche, le DOM, et l'image suivante — c'est ce temps-là qui est jugé.

| Étape | Clics | Objets | Pixel + mapping, p50 | Fiche | DOM | **Peinte, p50** | **p95** | max | Complètes |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| départ, aperçu (80 µm) | 84 | bille de 3 Mb | 22,7 | 0,5 | 1,3 | 28,2 | 36,8 | 182,7 | 84/84 |
| noyau | 84 | bille de 750 kb, fibre | 21,5 | 0,3 | 1,3 | 27,2 | 36,7 | 82,8 | 84/84 |
| fibre | 84 | monomère, site CTCF | 30,5 | 0,2 | 1,1 | 36,7 | 49,0 | 55,4 | 84/84 |
| nucléosomes | 84 | ADN, nucléosome | 22,1 | 0,2 | 1,6 | 28,3 | 38,3 | 69,4 | 84/84 |
| hélice | 84 | base | 23,4 | 1,2 | 1,6 | 31,2 | 37,3 | 46,3 | 84/84 |
| base | 84 | base | 23,5 | 1,2 | 1,4 | 30,6 | 39,8 | 67,7 | 84/84 |
| **tout** | **504** | | 23,2 | 0,4 | 1,4 | **29,5** | **43,9** | 182,7 | **504/504** |

(ms ; fiches préchargées, réseau non bridé.) **Critère tenu : p95 43,9 ms, et 504 fiches complètes
sur 504.** Un clic sur 504 dépasse 100 ms : le premier de la plongée, 182,7 ms, dont 176 dans la
lecture du pixel — la première lecture de la cible après l'ouverture. Composer la fiche coûte
0,2 à 1,2 ms (plus à la base, qui lit sa séquence et ses variants) ; l'écrire dans le DOM, 1 à 2 ms.

**Un défaut trouvé par cette mesure.** La première passe donnait 1,1 s par clic, dont 1,1 s dans
la lecture du pixel : `readPixels` attendait l'image en cours de rendu, et la page en rendait une à
chaque image — le repère du clic forçait un rendu WebGL, et l'échelle de l'hélice redessinait sans
cesse. SwiftShader met ~1 s à rendre 1 280 × 720. Corrigé : la page ne redessine que si des
nucléosomes se génèrent ou si de la séquence arrive (`DiveScene.busy`), et le repère se dessine sur
la surimpression 2D. Une lecture de pixel sur une vue immobile coûte alors 0,33 ms ; ce qui reste
dans la colonne « pixel » (~23 ms) est la synchronisation avec le compositeur du navigateur, qui
repeint la page après chaque clic — en logiciel ici. **Ce que ça ne couvre pas** : un clic pendant
une animation attend l'image en vol, soit une image de rendu (2,6 ms sur l'Iris Xe de la semaine
10, ~1 s sous SwiftShader).

**Ce que le préchargement achète, et ce qu'il coûte.** Il lit les tables, puis le noyau (le
chromosome visé d'abord), l'aperçu et les régions fines, une requête à la fois : 6,5 Mo en 185
requêtes, **3,4 s** sans bridage — mais **147 s sous Slow 4G** (tables à 44 s, noyau à 90 s). Pendant
ce temps, un clic montre aussitôt ce que le `.g3d` sait (objet, coordonnées, copie, précision :
28 ms) et complète la fiche quand son bloc arrive : sous Slow 4G, 13 clics sur 16 faits avant la
fin du préchargement ont attendu, **1,3 à 5,4 s**. Après lui, 60/60 prêts, peints en 42 ms (p95).
Témoin **sans préchargement**, Slow 4G : 15 clics sur 30 attendent leur bloc ; fiche complète en
1,7 s médian, 7,2 s au 95e centile. Un défaut trouvé en mesurant : les clics à froid attendaient
la fin du préchargement entier (150 s) au lieu de l'ouverture du fichier ; corrigé avant la
mesure ci-dessus.

### Sur du vrai matériel : le clic tient, les 60 images par seconde non

Le banc matériel republié (version 4) a été joué deux fois sur un **iPhone 14** (Apple GPU,
canevas de 716 × 403 pixels de rendu — dpr 3, plafonné à 2 par la page), le 2026-10-20 à 01:39 et 01:40 UTC — la plongée entière, du noyau à
la base c.818, puis 90 clics chronométrés chacun à trois étapes (bille, fibre, base) :

| Passe | rAF méd. / p95 | Images > 25 ms | Rendu seul méd. / p95 | Clic méd. / p95 / max |
|---|---|---:|---|---|
| 01:39:53 | 17 / **28 ms** | 7,5 % | 5 / 29 ms | **0 / 0 / 1 ms** |
| 01:40:07 | 17 / **31 ms** | 10,4 % | 7 / 31 ms | **0 / 1 / 1 ms** |
| rappel S10, même appareil, plongée jusqu'aux compartiments | 17 / 17 ms | 0 % | 3 / 7 ms | — |

**Le clic tient, largement** : lecture du pixel et mapping inverse coûtent 0 à 1 ms sur un vrai
GPU — Safari arrondit à la milliseconde —, contre ~23 ms sous SwiftShader, où l'essentiel était la
synchronisation avec le compositeur logiciel.

**Les 60 images par seconde ne tiennent pas** sur la plongée jusqu'à la base : p95 de 28 et 31 ms,
verdict « non tenu » de la page. Ce n'est pas un chargement — la passe suit un préchargement
complet —, c'est le rendu, et il est localisé : les images lentes tombent entre **59 % et 92 % du
trajet** (indices 141 à 220 sur 240), c'est-à-dire, d'après le tableau de continuité de 600
images, de l'arrivée des nucléosomes (170 nm de champ) jusqu'à l'apparition des lettres (6,3 nm).
Le rendu seul y vaut **22 à 27 ms** par décile du trajet, contre 2 à 3 ms à l'échelle du noyau ;
la vue finale de la base revient à 16,7 ms. C'est le premier chiffre matériel pour les échelles
des semaines 11 à 13 : les mesures de la semaine 10 s'arrêtaient aux compartiments. Pistes,
ni mesurées ni engagées — reportées à la semaine 22 par décision de l'utilisateur : la fenêtre de nucléosomes et d'ADN dessinée (jusqu'à 40 monomères de
contexte de part et d'autre de la base), et l'occlusion ambiante, calculée à la pleine résolution
du canevas — sa part n'est pas isolée par le banc.

### Non-régression

| Mesure | Semaine 12 | Semaine 13 |
|---|---:|---:|
| `make test` | 242 | **276** (+ 10 liftover, + 23 fiche, + 1 parseurs) |
| `pnpm verify` | 70 | **82** (+ 12 picking) |
| lecture croisée `.g3d` | 46 | **46** |
| lecture croisée des fiches | — | **6/6** (160 fiches, 1 000 positions T2T) |
| continuité, 600 images | 0 saut, transitions 0,94–1,00× | **identique** (témoin : 3 sauts, jusqu'à 7,41×) |
| lisibilité S11 : ancres, anneau, arc, corps, flèches | 269 · 131 · 1 662 · 2 961 · 85 % | **identiques** ; raccords 45 % et 11 % |
| fidélité des bases | 90/90 (témoin 23 %) | **90/90** (témoin 23 %) |
| temps réel Slow 4G, images en attente | 0 (témoin 1) | **0** (témoin 1), 11,8–31,8 Kio lus |
| premier rendu S9, aperçu, préfixe 16 Kio, Slow 4G | 1 418 ms | **1 478 ms** (marge 22 ms) |
| page du premier rendu | 15,1 Kio | 15,8 Kio |

**La marge de S9 a fondu, et ce n'est pas la semaine 13.** Le code de la semaine 12 (`ceca1ee`),
reconstruit et remesuré le même jour, donne **1 487 ms** : l'environnement a dérivé de ~70 ms
depuis le 18/10, et les phases réseau (document 606 ms, en-tête 1 278 ms, données 1 285 ms) sont
identiques à la milliseconde ; c'est le rendu logiciel qui est plus lent. Ce que la semaine 13
coûte, mesuré à part : la cible RGBA32UI, **9 ms** (1 470 ms en R32UI, expérience non conservée).
Deux corrections pour en arriver là : la page relisait **tout** le tampon d'identifiants dans la
fenêtre chronométrée pour se vérifier — quatre fois plus d'octets en RGBA32UI —, elle s'arrête
désormais sur la lecture d'un pixel, qui attend autant la fin du GPU, et vérifie ensuite ; et le
cache des intervalles de billes, ajouté à `G3D` que la page embarque, en est sorti (`BeadIndex`).
Avant elles : 1 511 ms, en échec.

### Ce que ça ne dit pas

- **La fiche est réelle ; la position 3D ne l'est pas.** Cliquer une bille donne les annotations
  réelles de son intervalle génomique ; l'endroit où la bille est dans le noyau reste celui du
  modèle simulé, et sa ligne « position radiale » le dit.
- **LAD, speckles et NAD viennent d'autres lignées.** Les LAD sont conservés à ~70 % entre types
  cellulaires ; un LAD de K562 n'est pas un LAD de GM12878, et le badge le dit.
- **Un seul appareil réel** a joué le banc de la semaine 13 (iPhone 14) : le clic y tient en 0–1 ms,
  les 60 fps de la plongée jusqu'à la base non. Ni GPU de bureau ni GPU intégré de PC n'ont encore
  rejoué cette plongée-là.
- **Le budget de taille n'est pas tenu** (6,5 Mo pour 2,5) : sous Slow 4G, le préchargement complet
  dure deux minutes et demie, et un clic qui le précède attend son bloc.
- **Pas de cytobande T2T** : UCSC injoignable.

### Reproduire

```
cd pipeline && ONLY="…" TIERS=all bash scripts/fetch.sh && make data-lock   # les 24 sources
pip install -e '.[hic,fiche]'
make fiche                       # écrit data/g3d/gm12878.fiche, puis geno fiche check
make test
cd .. && pnpm build && pnpm verify
node tools/test-fiche.mjs        # lecture croisée des fiches, Python → TypeScript
node tools/test-g3d.mjs
node tools/pick.mjs 500          # clic → fiche peinte : chaud, Slow 4G, témoin
```

---

## S12 — Du noyau à une base de TP53 : la vraie séquence, une hélice qui tient, aucun chargement visible

### Le critère, traduit

**Fait quand : un zoom continu enregistré, du noyau jusqu'à une base précise de TP53, sans
coupure ni chargement visible.** Quatre mots, quatre mesures, fixées avant de mesurer :

1. **« une base précise de TP53 »** : la base est désignée par l'annotation, pas recopiée — c.818
   du transcrit MANE de GENCODE, codon 273 — et ce que l'image dessine à cet endroit est la base
   du génome ;
2. **« la double hélice »** : la géométrie B est juste — pas, twist, sillons, sens — et le chemin
   de l'ADN qui la porte a une longueur, à la base près ;
3. **« sans coupure »** : la plongée, de 80 µm à 4 nm de champ, ne saute jamais — sur le rendu
   image par image, et sur la vidéo enregistrée elle-même ;
4. **« sans chargement visible »** : jouée en temps réel, réseau bridé, aucune image ne montre
   une hélice due sans sa séquence ou des nucléosomes dus pas encore générés — et un témoin sans
   préchargement, lui, en montre.

### Ce qui a changé sous le projet : le réseau

Le 2026-10-18, les dépôts primaires répondent ; ils renvoyaient 403 depuis la semaine 1. La
semaine 12 n'en prend que ce dont elle a besoin, verrouillé par le manifeste et confronté à
l'empreinte publiée par chaque dépôt ([`DATA_SOURCES.md` § 8](DATA_SOURCES.md)) : `hg38.2bit`
(md5 = celui d'UCSC), GENCODE v50 basic (md5 = celui d'EBI), les longueurs de chromosomes.

| Contrôle | Résultat |
|---|---|
| table interne des 24 longueurs contre `hg38.chrom.sizes` | **identiques**, à la base près |
| fenêtre chr17:7 673 797–7 673 808 lue par notre lecteur `.2bit`, contre l'API UCSC | `AAACACGCACCT` les deux fois |
| c.818 de TP53 recalculé depuis les CDS de ENST00000269305.9 (MANE Select) | chr17:7 673 802 = rs28934576 (R273H) |
| codon 273 lu dans le génome, brin codant | **CGT** (Arg) |
| c.1 de TP53 et d'ACTB | **ATG** les deux |

Un défaut trouvé en chemin : le caryotype cherchait `hg38.chrom.sizes` dans `data/core/`, où
`make data-core` ne dépose rien (il écrit dans `data/raw/core/`). Le fichier verrouillé n'aurait
jamais fait foi. Corrigé ; les longueurs étant identiques, aucun modèle n'en change.

### La région fine de TP53 : la méthode de S8, rejouée, et un piège de graine

`make fine-tp53` : chr17:5,67–9,67 Mb, 2 000 monomères de 2 kb, 1 200 conformations, sites CTCF
**plantés** comme en S8 (décision prise avec l'utilisateur : le rattrapage des données réelles
est une étape à part).

La première exécution reprenait les graines d'ACTB, et le bilan des boucles était **identique**
à celui d'ACTB, au dixième près. La cause : le 1D (extrusion, numpy à graines fixes) et le
placement des sites ne dépendent que des graines, pas de la région — les mêmes sites aux mêmes
rangs, les mêmes boucles. C'eût été afficher sous TP53 une copie de l'agencement d'ACTB.
Relancée avec ses propres graines (`--site-seed 17 --first-seed 17000`) :

| | ACTB (S8) | TP53 |
|---|---:|---:|
| frontières retrouvées (deux sens bloqués) | 12/13 | 11/13 |
| points d'angle convergents / hasard | 2,57× | 2,67× |
| … sans CTCF | 1,06 | 0,96 |
| arcs au seuil de 5 % : convergents / divergents / tandem | 11/14 · 0/3 · 1/3 | 13/15 · 0/5 · 1/4 |
| arcs sans CTCF | 0 (max 0,003) | 0 (max 0,011) |
| raccord, pire rapport R(s) | 2,42× | 2,54× |
| pose sur le noyau, RMSD | 426 nm | 375 nm |

Le mécanisme se reproduit sur une autre région, et le raccord échoue de la même façon.

### Le chemin de l'ADN : une longueur, à la base près

Mesuré sur 30 monomères (60 000 paires) d'une fibre synthétique à volume exclu, dans la suite
headless S12 :

| Propriété | Mesure | Critère fixé d'avance |
|---|---|---|
| pas sur les liens | 0,34 nm, écart max **0,37 %** | 1 % |
| pas sur l'ADN enroulé (imposé par 1,67 tour de 4,18 nm) | 0,3017 nm, écart max 0,02 % | 1 % |
| plus grand coude d'une base à la suivante | **16,7°**, 0 pas au-delà de 20° | 20° |
| twist sur l'ADN enroulé | 34,29°, écart max 0,45° (torsion de la superhélice) | 0,6° |
| twist moyen | 34,04° ; tous les pas à droite | 34,29 ± 0,5° |
| superhélice gauche | 300/300 cœurs | tous |
| cœurs à moins de 10 nm | 0 dans un monomère, 0 entre voisins | 0 ; ≤ 2 % entre monomères |
| déterministe, indépendant de l'ordre de calcul | oui ; une autre graine déplace tout | — |

Le chemin y est arrivé par quatre échecs, que le test a vus et que voici, parce que chacun
aurait produit une hélice fausse à l'écran :

1. **Liens allongés par leurs tangentes** : la première Hermite s'allongeait en grandissant ses
   tangentes, et rebroussait — 174° d'une base à la suivante, et des pas de 0,1 nm. Remplacée
   par une bosse en sin²(πu), qui allonge sans toucher aux tangentes d'extrémité.
2. **Corde visée trop courte** : les cœurs étaient placés au plus près du précédent, et les
   liens recevaient deux fois leur corde en longueur ; ils se repliaient en épingle (115°). La
   corde visée est maintenant la longueur moyenne d'un lien moins 15 % de mou, et les pas sont
   partagés en proportion des cordes.
3. **Tangentes à rebours** : une sortie de cœur tournée à l'envers de la corde annule la dérivée
   de l'Hermite près du bout — 97°. Avec des tangentes dans le sens de la corde et λ ≤ corde, la
   courbe avance toujours le long de celle-ci : c'est maintenant une quasi-contrainte du placement,
   avec une anticipation vers la cible suivante.
4. **Régler jusqu'à ce que ça passe** : chaque réglage déplaçait le défaut ailleurs (26°, puis
   40°, au gré du tirage). La construction juge désormais chaque réalisation de monomère et la
   retire si un lien se coude de plus de 18° par base. Le critère du test, 20°, n'a pas bougé.

Les cœurs de deux monomères voisins se touchaient dans 19 % des paires avant la contrainte des
demi-espaces (chaque cœur à au moins 5,5 nm des plans de jonction de son monomère) ; 0 après.

### La double hélice : géométrie, couleurs, picking

| Propriété | Mesure | Référence |
|---|---|---|
| rayon des phosphates | 8,91 Å | modèle de fibre B (Arnott) |
| pas de l'hélice | 3,5700 nm, 10,5 pb/tour | 3,57 nm |
| petit sillon / grand sillon (P–P le plus court − 5,8 Å) | **5,4 / 11,8 Å** | 5,7 / 11,7 Å publiés |
| témoin sans décalage axial ±2,85 Å | rapport des sillons 1,11 au lieu de 1,57 | doit s'égaliser |
| barreaux qui traversent un squelette | 0 | 0 |
| couleur de chaque pixel de barreau = sa base | **100 %** de 7 236 pixels | > 99 % |
| témoin : séquence décalée d'une base | 18,8 % | < 50 % |
| picking d'un barreau, contre un lancer de rayon analytique | **201/203** pixels non ambigus | ≥ 99 % |

Le picking a demandé deux corrections **de l'oracle**, pas du rendu : viser le milieu d'un
barreau comptait faux quand un voisin à 0,34 nm le cachait légitimement ; puis l'oracle donnait
tout un segment de squelette à sa première base, quand le tube le partage entre ses deux points
(« un pixel appartient au point le plus proche »). Ce qui est garanti est qu'un clic sur un
barreau rend sa base ; le squelette entre deux phosphates n'appartient à une base que par
convention, et il est compté à part.

**Fidélité dans la plongée elle-même** (`node tools/fidelity.mjs`) : la vue finale, 4 nm de
champ sur TP53 c.818, chaque pixel de barreau classé par teinte, une lettre votée par base
visible, contre la séquence relue par Python **directement dans `hg38.2bit`** — ni par le `.g3d`,
ni par le navigateur : **90/90 bases**, témoin décalé d'une base 23 %. La base visée est dessinée
C sur le brin + et G sur le brin −, comme dans le génome. Une première version éclaircissait la
base visée pour la distinguer : son C virait au cyan, lu A par la mesure — et par l'œil. Elle
garde sa couleur ; seule son épaisseur la marque.

### « Sans coupure »

`node tools/dive.mjs 600` : la plongée entière, de 80 µm à 4 nm de champ, 600 images à pas
constant sur le chemin de van Wijk et Nuij, écart moyen de chaque image à la précédente. Même
détecteur qu'en S10 et S11 — un saut est un écart trois fois plus grand que la médiane de ses 20
voisines. Le témoin fait les mêmes transitions en marches sèches.

| Transition | Champ | Fondu | Sèche |
|---|---:|---:|---:|
| niveau de détail 3 Mb → 750 kb | 19,8 µm | 0,94× | 1,37× |
| autres chromosomes effacés | 12,5 µm | 0,98× | 2,88× |
| territoires → compartiments | 5,26 µm | 0,96× | 3,00× |
| billes 750 kb → fibre 2 kb (raccord S8) | 2,10 µm | 1,00× | 7,41× |
| la boucle s'isole | 0,91 µm | 0,98× | 4,51× |
| fibre → nucléosomes | 170 nm | 0,96× | 2,39× |
| ADN de 1 nm → double hélice | 27,8 nm | **0,96×** | 2,76× |
| les lettres paraissent | 6,3 nm | 0,97× | 0,97× |

**0 saut** dans la plongée, 3 dans le témoin. Deux lectures honnêtes du tableau :

- le témoin de l'hélice (2,76×) passe **sous** le seuil de 3 fixé en S10, comme celui du niveau
  de détail (1,37×) : à 28 nm de champ, un ADN de 1 nm et une double hélice de 2 nm occupent
  presque les mêmes pixels, et une marche sèche ne change que ~3× l'écran. Le seuil, fixé
  d'avance, n'est pas abaissé pour l'attraper ; le rapport est publié ;
- les **lettres** sont une surimpression 2D, hors de l'image WebGL que ce banc relit : 0,97× des
  deux côtés ne prouve rien pour elles. Leur continuité est mesurée sur la vidéo, qui les
  contient (ci-dessous).

SwiftShader : 246 ms par image en médiane, p95 1 480 ms — une borne haute, sans portée pour les
fps.

### « Sans chargement visible »

`node tools/stream.mjs 9 40` : la plongée jouée **en temps réel** — la caméra suit l'horloge, une
image en retard a quand même avancé —, réseau bridé par CDP au profil Slow 4G de la semaine 9
(562,5 ms de latence, 180 ko/s), cache désactivé. Chaque image dit si elle montre un chargement :
hélice due sans sa séquence, ou nucléosomes dus pas encore générés.

| Durée | | Images | Chargement visible | dont séquence | dont nucléosomes | Séquence lue |
|---:|---|---:|---:|---:|---:|---:|
| 9 s (vitesse du bouton) | préchargement | 41 | **0** | 0 | 0 | 11,8 Kio |
| 9 s | témoin : à la demande | 41 | **1** (à 7,7 s, 18 nm de champ) | 1 | 0 | 7,8 Kio |
| 40 s | préchargement | 178 | **0** | 0 | 0 | 31,8 Kio |
| 40 s | témoin : à la demande | 183 | 0 | 0 | 0 | 7,8 Kio |

**Aucune image en attente avec préchargement**, aux deux vitesses ; et **seule la fenêtre est
lue** : 11,8 à 31,8 Kio de séquence sur les 2 Mo du fichier (la table d'empreintes, puis un à
sept blocs de 4 Ko). À dire avec le résultat :

- le témoin ne discrimine qu'à la vitesse rapide. À 40 s, la descente de 60 à 40 nm de champ
  laisse à une requête à la demande le temps d'arriver (≈ 0,7 s en Slow 4G) : il ne montre rien,
  et ce zéro-là ne prouve rien. À 9 s, il attend une fois — c'est cette ligne qui donne son sens
  au zéro de la plongée ;
- un premier essai, avant correction, comptait **1 image** de nucléosomes non générés à 248 nm
  de champ : la génération commençait à l'échelle de la boucle et coûtait 11 ms par monomère.
  Elle commence désormais dès la fibre, autour de la base visée dès l'ouverture, à 4,1 ms par
  monomère, sur un budget d'un quart d'intervalle entre images ;
- SwiftShader rend 4 à 5 images par seconde : la caméra avance par grands pas et la génération
  a moins d'images pour se faire. Sur du matériel à 60 images par seconde, la marge est plus
  grande, pas plus petite ; le réseau, lui, est simulé de façon réaliste, sans TCP ni TLS
  (même réserve qu'en S9).

### L'enregistrement

`node tools/record.mjs 45` : 1 350 images de la plongée, rendues une à une dans Chromium headless
en 1 280 × 720 avec leur surimpression (lettres, étiquette de la base visée, légende d'échelle
et de provenance), capturées par CDP, encodées par ffmpeg en H.264 à 30 i/s — 1,5 s tenue au
départ, 3 s à l'arrivée. **`dist/plongee-tp53.mp4`, 49,4 s, 21,1 Mo**, rendue en 17 min ; elle
n'entre pas dans git.

Sa continuité est mesurée **sur la vidéo encodée elle-même** — écart moyen de luminance d'une
image à la suivante (`tblend` + `signalstats` de ffmpeg), même détecteur qu'au banc : **0 saut**,
pire rapport 2,11, écart médian 2,49/255. Cette mesure-là voit la surimpression : les lettres
qui paraissent entre 8 et 5 nm de champ, que le banc WebGL ne voyait pas, ne font pas sauter
l'image.

### Le raccord S8, mesuré une seconde fois

Au milieu du fondu billes → fibre (2,1 µm de champ), sur l'image (`node tools/loop.mjs`) :

| Région | Part de la fibre dans la silhouette des billes qu'elle détaille | Écart des centres |
|---|---:|---:|
| chr7:a, ACTB | 45 % | 8 % de la hauteur |
| chr17:a, TP53 | **11 %** | 13 % |

Le fondu tramé rend la plongée continue à l'œil ; il ne rend pas la fibre juste. À TP53, ce qui
apparaît est presque entièrement **à côté** de ce qui disparaît. C'est le défaut S8, rendu
visible comme convenu, et la raison de la semaine de modélisation qui reste à planifier.

### Non-régression

| Mesure | Semaine 11 | Semaine 12 |
|---|---:|---:|
| `make test` | 233 | **242** (+ 5 séquence, + 4 format) |
| `pnpm verify` | 49 | **70** (+ 21 hélice) |
| lecture croisée `.g3d` | 37 | **46** (+ 9 séquence) |
| lisibilité S11 (`target=actb`) : ancres, anneau, arc, corps, flèches | 269 · 131 · 1 662 · 2 961 · 85 % · 0,95/0,93 | **identiques** |
| premier rendu S9, aperçu, préfixe 16 Kio, Slow 4G | 1 432 ms | **1 418 ms** (marge 82 ms) |
| page du premier rendu | 15,4 Kio | **15,1 Kio** |
| en-tête du `.g3d` / fin du premier rendu | 12 182 o | 13 854 o |

Deux défauts trouvés en chemin, l'un et l'autre corrigés avant la mesure :

- **La lisibilité S11 avait changé** — boucle 692–777 au lieu de 1527–1621, 0 pixel d'ancre +.
  ACTB porte désormais une cible (c.1), et la nouvelle règle « la boucle visée contient la base »
  remplaçait celle de S11. `target=actb` ne passe plus de base visée : la semaine 11 est rejouée
  à l'identique.
- **La page du premier rendu grossissait de 2 Kio** : les méthodes de séquence vivaient dans la
  classe `G3D`, qu'aucun élagage ne découpe. Elles sont dans `SequenceReader`, à part.

La variante « noyau seul » de S9, déjà en échec (1 602 ms), passe à 2 263 ms : son fichier porte
maintenant deux niveaux fins et la séquence dans l'en-tête, qui ne tient plus dans le préfixe de
16 Kio. C'est le fichier de comparaison, pas celui que sert la plongée ; la règle de S9 — l'aperçu
est le niveau du premier rendu — en sort renforcée.

### Ce que ça ne dit pas

- **La séquence est réelle ; sa position ne l'est pas.** La base c.818 est où le chemin d'ADN
  la pose : sur un nucléosome illustratif, sur une fibre simulée, aux CTCF plantés, posée sur un
  noyau simulé dont le raccord échoue. L'écran le dit à chaque échelle (ligne de provenance).
- **L'hélice est idéale.** 10,5 pb/tour partout, y compris sur le cœur (10,2 in vivo) ; aucun
  effet de séquence sur la géométrie locale (roulis, inclinaison, courbure intrinsèque) ; pas de
  sillons comprimés côté octamère.
- **Les nucléosomes sont placés, pas positionnés** : répétition de 200 pb en moyenne, liens de
  longueur variable, aucune MNase-seq lue.
- **L'enregistrement ne prouve pas l'absence de chargement** : rendu image par image, chaque
  image attend d'être complète. C'est la mesure en temps réel qui le prouve, sous SwiftShader —
  donc avec une caméra qui avance par grands pas, et un réseau qui, lui, est réaliste.
- **Les 60 fps de la plongée jusqu'à la base** ne pouvaient venir que du banc matériel. Ils en sont
  venus le 2026-10-20 (§ S13) : **non tenus** sur un iPhone 14 aux échelles nucléosome et hélice
  (rAF p95 28–31 ms, rendu seul 22–27 ms sur cette portion du trajet).

### Reproduire

```
cd pipeline && make data-core data-seq && make data-lock     # 9 entrées, empreintes contrôlées
make fine-tp53 && make g3d && make test
cd .. && pnpm build && pnpm verify
node tools/test-g3d.mjs          # lecture croisée, séquence comprise
node tools/dive.mjs 600          # continuité, 80 µm → 4 nm
node tools/stream.mjs 9 40       # temps réel, Slow 4G, avec et sans préchargement
node tools/fidelity.mjs          # les bases dessinées contre hg38.2bit
node tools/loop.mjs              # S11 inchangée, et le raccord S8 sur chr7 et chr17
node tools/record.mjs 45         # la vidéo, et sa continuité mesurée sur elle-même
```

---

## S11 — Une boucle CTCF se lit, ancres comprises ; le fondu cache le saut du raccord, pas l'écart

### Le critère, traduit

**Fait quand : une boucle CTCF se lit comme telle, ancres comprises, à l'échelle où elle
existe.** « Se lire » ne se mesure pas au goût. On l'a traduit, avant de mesurer, en propriétés
de l'image de la vue de boucle, relevées dans le **tampon d'identifiants** du picking — qui ne
garde que le fragment le plus proche, donc ce qui est *vu*, pas ce qui est dessiné :

1. les deux **ancres** — la flèche CTCF `+` au pied gauche, la `−` au pied droit — ont chacune
   au moins 50 pixels non occultés ; l'**anneau** de cohésine et l'**arc** de la paire aussi ;
2. au moins 80 % des monomères du **corps** de la boucle ont un pixel visible ;
3. chaque flèche suit **son brin vers l'intérieur** de la boucle, en 3D et à l'écran ;
4. les arcs dessinés sont les **bonnes paires** : jugés contre la vérité plantée, et contre le
   témoin sans CTCF où aucun ne doit passer le seuil ;
5. « à l'échelle où elle existe » : chaque transition de la plongée change l'écran autant que le
   mouvement de caméra seul, et le raccord noyau → fibre est **chiffré**, pas caché.

### Ce qui est dessiné

Sous 3 µm de champ, la bille de 750 kb cède à la **fibre de 2 kb** : un tube le long des 2 000
monomères du modèle de la semaine 8, posé dans le repère du noyau (semaine 9). Entre 1,4 µm et
600 nm, la boucle visée s'isole — le reste de la région s'efface par tramage — et ses annotations
apparaissent : **flèches CTCF** qui enveloppent la fibre dans le sens du motif, **anneaux** des
cohésines tenues dans la conformation affichée, **arcs** des paires d'ancres fréquentes dans la
population. Sous 250 nm, le tube cède aux **nucléosomes**.

Un arc n'est pas une liaison de la structure affichée : c'est une fréquence de population —
la fraction des 1 200 conformations où une cohésine tient la paire —, dessinée là où sont les
deux ancres dans cette structure-ci. Un anneau, lui, est l'état de cette conformation-là.

### Les jambes de cohésine, retrouvées sans refaire tourner OpenMM

Le fichier de la semaine 8 garde les positions, pas les jambes. L'étage 1D est du numpy pur à
graines fixes : le refaire tourner rend les jambes qui ont guidé le 3D — à condition de le
vérifier contre le 3D, puisqu'OpenMM, lui, ne se reproduit pas au bit (§ S10). Une cohésine tient
ses deux pieds à ~1 σ :

| | distance médiane des deux pieds |
|---|---:|
| jambes régénérées, paires tenues | **1,08 σ** |
| paires quelconques à même séparation génomique | 5,69 σ |
| témoin : jambes décalées d'une seule conformation | 2,18 σ |

Le garde-fou voit un décalage d'une seule conformation. Les exécutions futures sauvent les
jambes, avec `velocity` et `burnin_1d`, que le `.json` n'écrivait pas.

### Les arcs sont les bonnes paires

Les candidates ne doivent rien à la vérité : toutes les paires (site `+`, site `−` à sa droite)
à moins de 1 Mb, 66 après fusion des grappes de motifs voisins. Seuil fixé d'avance : une paire
tenue dans au moins **5 %** des conformations est une boucle de la population.

| | arcs dessinés |
|---|---:|
| domaines convergents plantés | **11 / 14** |
| domaines tandem plantés | 1 / 3 |
| domaines divergents plantés | 0 / 3 |
| paires hors de toute vérité | 1 |
| **témoin : sites CTCF inoccupés** | **0** (fréquence maximale 0,003) |

Les divergents (`− +`) ne sont pas candidats par construction : leur 0 ne prouve rien, et c'est
dit. Le tandem retrouvé est celui que la semaine 8 avait expliqué : le `−` du domaine suivant est
à un monomère de son ancre droite. Le seuil se lit dans les deux sens — à 10 %, 5 convergents
sur 14 seulement ; à 2 %, 12 sur 14 mais 4 paires hors vérité.

### La boucle visée, et une boucle de l'instant qui n'en était pas une

La première boucle visée était la plus proche d'ACTB parmi celles qu'une cohésine tient, ses
deux pieds sur des CTCF, dans la conformation affichée : 1341–1392, un domaine convergent planté.
Mais cette paire n'est tenue que dans **2,5 %** des conformations — une boucle de l'instant, pas
de la population, et elle n'avait pas d'arc. La règle est devenue : tenue ici **et** au-dessus du
seuil. La boucle montrée est **1527–1621** (188 kb, domaine convergent planté), tenue par
**deux** cohésines dans cette conformation et dans **18,7 %** des conformations de l'ensemble.
Elle est à 1,5 Mb d'ACTB, pas sur ACTB : c'est la plus proche qui remplisse les deux conditions.

### Elle se lit

`node tools/loop.mjs`, vue de boucle à 960 × 540 :

| | ancre `+` | ancre `−` | anneau | arc | corps visible | flèches vers l'intérieur | verdict |
|---|---:|---:|---:|---:|---:|---|---|
| vue de boucle | **269 px** | **131 px** | 1 662 px | 2 961 px | **85 %** | oui / oui (0,95 ; 0,93 à l'écran) | lisible |
| témoin : sans isolement | 0 | 0 | 0 | 0 | 7 % | — | non |

Sur les 26 ancres des 13 arcs dessinés, les 26 flèches suivent leur brin vers l'intérieur. Le
témoin dit ce que l'isolement achète : dans la pelote de 4 Mb, la boucle et ses ancres sont
entièrement cachées.

Il a fallu quatre corrections pour en arriver là, chacune vue sur l'image ou la mesure :
- les annotations de **toute la région** — 40 flèches, 20 anneaux, 13 arcs — flottaient autour de
  la boucle isolée, sans chromatine sous elles ; elles s'isolent maintenant avec elle ;
- la caméra regardait la base **à travers la pelote**, puis dans **l'axe** des deux flèches, fondues
  en une tache ; elle regarde maintenant de côté, perpendiculairement à l'axe base → boucle ;
- la boucle de l'instant ci-dessus ;
- la première mesure d'orientation comparait la flèche à la corde vers le monomère situé trois
  plus loin, et donnait −0,35 et −0,65 : à sa base, la chaîne fait demi-tour en 6 kb. Elle
  mesurait la pelote. La flèche se juge contre son brin, à un monomère.

### À l'échelle où elle existe : continue, et l'écart du raccord est sur l'image

`node tools/dive.mjs 360`, la plongée entière en 360 images à vitesse perçue constante :

| transition | champ | fondu | marche sèche |
|---|---:|---:|---:|
| aperçu 3 Mb → noyau 750 kb | 19,9 µm | 0,94× | 1,42× |
| autres chromosomes effacés | 12,5 µm | 1,06× | 3,06× |
| territoires → compartiments | 5,2 µm | 0,97× | 3,13× |
| **billes 750 kb → fibre 2 kb (raccord S8)** | **2,1 µm** | **0,98×** | **6,77×** |
| la boucle CTCF s'isole | 0,90 µm | 0,98× | 3,37× |
| fibre → nucléosomes | 0,17 µm | 0,97× | 2,91× |

**0 saut** dans la plongée, 4 dans le témoin. Mais le chiffre qui compte ici est ailleurs : le
raccord noyau → fibre ne « saute » pas, parce qu'un fondu étale n'importe quel changement sur
des dizaines d'images. Une mesure image à image ne distingue pas *la même chose qui change de
forme* de *une chose qui s'efface pendant qu'une autre apparaît ailleurs*. On a donc mesuré la
géométrie : au milieu du fondu, à la même caméra, les seules billes qui couvrent la région fine,
puis la seule fibre.

| au milieu du fondu (2,1 µm de champ) | |
|---|---:|
| part de la fibre dans la silhouette des 6 billes de 750 kb qu'elle détaille | **45 %** |
| écart des deux centres | **8 %** de la hauteur de l'image |

**Le fondu cache le saut, pas l'écart.** C'est le défaut de la semaine 8 rendu visible, et
mesuré, comme convenu : la fibre est une pelote de 435 nm de rayon, les six billes un chemin
étiré, et la pose rigide laisse 426 nm de RMSD. Ce chiffre décrit ; il n'a pas de témoin, faute
d'un niveau fin qui raccorderait.

### Le rendu, prouvé correct

Suite headless de la semaine 11 (`pnpm verify`, 13 assertions) :

| propriété | mesure | témoin |
|---|---|---|
| transport parallèle : jamais plus de rotation que la tangente | excès max 0,0001° | Frenet : 3 retournements de 180° |
| tube centré sur ses points, au bon rayon | 26,0 px pour 25,8 attendus | — |
| **aucune fissure** entre segments de tessellations différentes | **0 trou** sur 11 711 pixels sondés | côtés par segment : 257 trous |
| faces avant vers l'extérieur | 13 312 px sur 13 312 avec élimination de l'arrière | — |
| tessellation adaptative | 1 592 triangles à 60 unités, 15 328 à 4 | — |
| nucléosome de face, de profil | 104 × 104 px ; 104 × 52, rapport 0,500 | — |
| profondeur du point d'impact | au 10⁻⁶ près | — |
| ellipsoïdes qui s'interpénètrent | 0 pixel de différence selon l'ordre | — |

Deux de ces tests ont d'abord été faux. Le seuil de rotation du transport (12°) était arbitraire :
sur une courbe plane la normale *doit* tourner avec la tangente, et c'est l'excès de rotation qui
signe la faute. Le test de fissure ne voyait rien — témoin compris — tant que le tube était
dessiné des deux côtés : une fissure montre l'intérieur du tube, pas le fond. En éliminant les
faces arrière pour le voir, on a trouvé l'**enroulement des triangles inversé** — faces avant vers
l'intérieur —, corrigé. Le témoin ne fuyait alors que 2 pixels avec une sonde de 7 × 7 sur l'axe :
une fissure s'ouvre près de la silhouette. La sonde couvre maintenant 80 % du rayon projeté (pas
90 % : un tube à 8 côtés n'a que cos(π/8) = 0,92 r entre deux arêtes).

### Les nucléosomes, illustratifs

Environ dix par monomère de 2 kb — l'ordre de grandeur d'une répétition de 200 pb —, cœurs de
11 × 5,5 nm, liés par l'ADN (2 nm), en chaîne désordonnée confinée dans le volume du monomère,
distance minimale 10 nm. **Ce n'est pas un modèle** : aucune donnée de positionnement n'est
lisible, et les positions sont tirées, déterministes par monomère. La forme est juste — une chaîne
désordonnée, pas une fibre de 30 nm, que la ChromEMT n'observe pas *in situ* (Ou 2017) — et les
nombres sont des ordres de grandeur.

### Performance

SwiftShader, 960 × 540 : 271 ms par image en médiane, 2 277 ms au 95e centile — la génération des
19 990 nucléosomes à la première descente, et 20 000 ellipsoïdes sur un rasteriseur logiciel. Borne
haute sans portée. La page de banc matériel est republiée avec la plongée entière ; les mesures de
la semaine 10 y restent, étiquetées. **Aucune mesure de la semaine 11 sur vrai matériel encore.**

### La semaine 9 tient encore, de justesse

| aperçu 3 Mb, Slow 4G | semaine 9 | semaine 10 | semaine 11 |
|---|---:|---:|---:|
| préfixe 16 Kio | 1 412 ms | 1 416 ms | **1 432 ms** |
| préfixe ajusté | 1 381 ms | 1 381 ms | 1 409 ms |
| fin du premier rendu | 11 149 o | 11 619 o | 12 182 o |

La marge passe à 68 ms. Elle avait fondu à 62 ms : le shader du tube, gabarit à `${…}` au niveau
du module, entrait dans la page du premier rendu — 15,1 → 17,4 Kio, le piège de l'occlusion en
semaine 10. Corrigé de la même façon. Les descriptions longues des liens ont aussi quitté
l'en-tête, qui fait partie de la première requête : 12 389 → 12 182 octets.

### Ce que ça ne dit pas

**Rien d'une boucle mesurée** : les sites sont plantés, la boucle est celle du modèle. **Rien des
nucléosomes réels.** **Rien des fps de la vue de boucle sur vrai matériel**, tant que la page n'a
pas reçu de mesure. **Le raccord S8 n'est pas corrigé** : il est montré. La boucle visée n'est pas
sur ACTB.

### Reproduire

```console
$ cd pipeline && make g3d                       # liens, sites CTCF, bilan des boucles imprimé
$ cd .. && pnpm verify                          # 16 (S4) + 20 (S10) + 13 (S11) assertions
$ node tools/test-g3d.mjs                       # lecture croisée, liens compris (37)
$ node tools/loop.mjs                           # lisibilité de la boucle, témoin, raccord S8
$ node tools/dive.mjs 360                       # continuité de la plongée entière
$ node tools/loadtime.mjs 10                    # la semaine 9 tient-elle encore
```

---

## S10 — La plongée est continue, et tient 60 images par seconde sur un GPU intégré

### Le critère, et les deux moitiés qu'il contient

**Fait quand : on descend du noyau aux compartiments en plongée continue, à 60 fps sur la cible
matérielle médiane.** Deux affirmations de nature différente. « Continue » est une propriété de
l'algorithme : elle se mesure ici, image par image, contre un témoin construit pour échouer.
« 60 fps sur la cible médiane » est une propriété du matériel : l'environnement ne rend que par
SwiftShader, un rasteriseur logiciel, et aucun chiffre qui en sort ne la prouve ni ne l'infirme.
La première moitié est tenue ici. La seconde l'est sur deux vrais appareils, mesurés après coup
par la page de banc matériel (dernière section).

### Ce qui est montré

`dist/dive.html` lit le `.g3d` de la semaine 9 et descend en une plongée de la vue d'ensemble
(80 µm de champ) jusqu'à la bille d'ACTB sur chr7:a (2,6 µm). En chemin, trois transitions :

| transition | où | comment |
|---|---|---|
| aperçu 3 Mb → noyau 750 kb | bille du noyau entre 3 et 6 px de rayon | fondu tramé complémentaire |
| les 45 autres copies s'effacent | champ de 20 à 8 µm | extinction par groupe, tramée |
| territoires → compartiments | champ de 8 à 3,5 µm | seconde palette, mélangée dans le shader |

**L'état visuel est une fonction de la caméra, pas du temps** (`apps/dive/scene.ts`) : remonter
défait exactement ce que descendre a fait, et un arrêt en chemin montre un état stable. Seule
l'arrivée du noyau par le réseau est fondue sur 400 ms ; la mesure la neutralise.

### La caméra : vitesse perçue constante

Le chemin de van Wijk et Nuij (2003) minimise la longueur dans la métrique
`ds = √(ρ² du² + dw²/ρ²) / w`, où un déplacement compte en largeurs d'écran et un zoom en
facteurs. Parcouru à pas constants, il doit avoir une vitesse constante dans cette métrique — et
c'est cela qui est testé, pas seulement ses extrémités :

| | vitesse perçue, min – max |
|---|---:|
| chemin de van Wijk et Nuij, 20 → 0,8 µm de champ et 3 µm de côté | **1,0000 – 1,0000** |
| témoin : distance interpolée linéairement, même trajet | 0,712 – 16,747 |

Le témoin varie d'un facteur 23 : il passerait l'essentiel du trajet loin de la cible et
s'écraserait sur elle en quelques images.

### Le fondu tramé ne troue jamais l'image

Un fondu par alpha exigerait de trier 8 000 billes à chaque image et casserait la profondeur
écrite par les imposteurs. On trame : chaque pixel est gardé ou jeté selon un seuil fixe à
l'écran, et le niveau suivant, en mode inversé, garde exactement les pixels jetés. Quarante
billes qui se chevauchent, dessinées en double avec des identifiants distincts :

| t | trous | écart de profondeur | part du second niveau |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0,000 |
| 0,25 | 0 | 0 | 0,250 |
| 0,5 | 0 | 0 | 0,500 |
| 0,75 | 0 | 0 | 0,749 |
| 1 | 0 | 0 | 1,000 |
| témoin : deux trames identiques, t = 0,5 | **10 433 sur 20 886** | — | — |

Aucun trou, la profondeur visible identique au bit près à celle d'un niveau seul, et la part
du second niveau égale à `t`. Le témoin troue la moitié de l'image : la complémentarité est
la propriété qui compte, et le test la distingue. Un groupe éteint disparaît entièrement (5 619
pixels avant, 0 après) et ne coûte rien : son quad sort du volume de vue dans le vertex shader.

### L'occlusion ambiante : un pli, et la même image à toute échelle

| | |
|---|---:|
| bille seule, de face | 1,000 |
| fond | 1,000 |
| deux billes tangentes, flanc → face → pli (moyenne 5 × 5) | 1,000 ≥ 0,992 ≥ **0,802** |
| scène ×10 et rayon ×10 : pixels qui changent | **2 sur 7 687**, d'un échantillon chacun |
| témoin : scène ×10, rayon inchangé | 7 965 sur 13 986, de plus d'un échantillon |

L'invariance d'échelle est ce qui rend l'occlusion utilisable dans une plongée de deux ordres de
grandeur. Elle n'est pas exacte au bit : un échantillon tangent à une surface bascule d'un côté
ou de l'autre du test de profondeur selon l'arrondi flottant, et vaut 1/16. Le test compte donc
les pixels qui basculent et de combien, au lieu d'exiger l'égalité.

### Les compartiments du modèle, et un biais que le témoin a montré

Aucun compartiment n'est donné au modèle. On tire de l'ensemble de 200 structures la carte
qu'une expérience Hi-C aurait vue — homologues confondus, contacts au seuil de la semaine 7 — et
on y appelle E1 avec les appelants de la semaine 3, orienté par la piste LAD (A du côté pauvre
en LAD, faute de GC lisible). L'X inactif ne reçoit aucun compartiment : il n'en a pas in vivo.

**Le premier témoin a révélé un biais de mesure.** La force du damier, `(AA + BB) / (AB + BA)`
sur les quintiles extrêmes de E1, se calcule d'habitude sur la carte même dont E1 est tiré. Sur
un ensemble jouet **sans aucun compartiment**, elle vaut **2,80** au lieu de 1 : le vecteur
propre est par construction la partition qui fait le meilleur damier de cette carte, bruit
compris. On la croise donc : E1 d'une moitié des structures, carte de l'autre, dans les deux sens.

| ensemble jouet | force croisée | force sur l'échantillon même | moitié contre moitié |
|---|---:|---:|---:|
| damier planté | 8,37 | 26,59 | r = +0,879 |
| témoin sans damier | **0,97** | **2,80** | r = +0,003 |

Sur l'ensemble réel :

| | |
|---|---:|
| bins appelés | 3 976 sur 4 041 |
| part de A | 50 % |
| force du damier, croisée | **4,76** (5,09 sur l'échantillon même) |
| E1 d'une moitié des structures contre l'autre | **r = +0,980** |
| E1 contre la piste LAD | r = −0,49 |
| E1 contre la position radiale moyenne | **r = −0,60** |
| billes avec un compartiment | 7 747 sur 8 082 (Xi et bins non appelés exclus) |

Le damier est fort et reproductible. Mais il dit d'où il vient : E1 suit la **position
radiale** plus que la piste LAD elle-même. Le seul ingrédient du modèle capable de ségréger la
chromatine est le terme radial de la semaine 6 ; les compartiments du modèle sont sa
stratification radiale, vue à travers les contacts. C'est une prédiction, que le Hi-C de Rao
2014 confirmera ou non le jour où il sera lisible — `make g3d COMPARTMENTS=e1.bedGraph`
l'affiche alors à la place, et chiffre l'accord avec le modèle. La corrélation négative avec
LAD est en partie imposée par l'orientation ; son amplitude, la force et la reproductibilité ne
le sont pas.

### La plongée est continue

`node tools/dive.mjs 240` parcourt la plongée en 240 images à pas constants et mesure l'écart
moyen absolu entre images successives. Le mouvement seul produit un écart régulier (médiane
1,48/255) ; une transition qui saute fait un pic. Règle fixée avant la mesure : un **saut** est
une image dont l'écart dépasse 3 fois la médiane de ses 20 voisines.

| | sauts | pire rapport, toutes images |
|---|---:|---:|
| la plongée | **0** | **1,05** |
| témoin : les mêmes transitions en marches sèches | 2 | 3,84 |

Transition par transition, le rapport de l'écart à la médiane de ses voisines :

| transition | image | champ | fondu | sèche |
|---|---:|---:|---:|---:|
| aperçu 3 Mb → noyau 750 kb | 98 | 19,8 µm | **0,95×** | 1,69× |
| autres chromosomes effacés | 130 | 12,5 µm | **1,04×** | 3,24× |
| territoires → compartiments | 191 | 5,2 µm | **0,97×** | 3,84× |

Chaque transition fondue est indiscernable du mouvement de caméra seul. Le seuil de 3 fixé
d'avance attrape deux des trois marches du témoin et **manque la troisième** : le changement de
niveau sec ne fait « que » 1,69 fois le mouvement, parce qu'à 20 µm de champ une bille de 3 Mb
et ses quatre filles de 750 kb se ressemblent. Ce n'est pas une raison de baisser le seuil après
coup ; c'est la raison pour laquelle le tableau par transition est publié — il sépare les deux
séries partout, d'un facteur 1,8 au moins. Deux exécutions donnent les mêmes écarts au pixel près.

### Les 60 images par seconde : mesurées sur deux vrais appareils

| SwiftShader, 960 × 540 | médiane | p95 |
|---|---:|---:|
| rendu d'une image, GPU attendu | 120 ms | 248 ms |
| intervalle entre images (requestAnimationFrame) | 133 ms | — |

Environ 8 images par seconde, sur un rasteriseur qui exécute les fragment shaders sur quatre
cœurs de processeur. C'est une borne haute sans portée sur la cible : ne pas lire « 8 fps ».

Le chiffre qui manque vient de **la page de banc matériel** (`dist/materiel/index.html`, publiée
comme page privée sur claude.ai) : elle joue la même plongée, à la résolution de l'appareil qui
l'ouvre, en deux passes — synchronisée, pour la marge du GPU ; au rythme de
`requestAnimationFrame` sans aucune relecture, pour ce que voit l'utilisateur — et enregistre le
résultat là où Claude le relit. Verdict « tenu » si 95 % des intervalles sont sous 20 ms, seuil
écrit dans la page avant toute mesure. Deux appareils l'ont ouverte, 240 images par passe :

| appareil | GPU | canevas | intervalle rAF, médiane / p95 | images > 25 ms | rendu seul, médiane / p95 |
|---|---|---:|---:|---:|---:|
| iPhone 14 | Apple GPU | 716 × 403, dpr 3 | 17 / 17 ms | 0 % | 3 / 7 ms |
| PC de bureau | Intel Iris Xe, **intégré** | 960 × 540, dpr 1 | 16,7 / 16,8 ms | 0 % | 2,6 / 4,2 ms |

**Les 60 images par seconde tiennent sur les deux**, sans une image perdue : l'intervalle est celui
de l'écran à 60 Hz, et le rendu seul en occupe 15 à 20 %. L'Iris Xe est un GPU intégré, la classe
de la cible médiane (iGPU de portable) ; c'est lui qui tranche le critère. Le budget de la semaine
4 le laissait attendre — 10 114 billes au plus par image quand le spike en dessinait 620 000 —,
c'est maintenant mesuré.

Trois réserves, à garder avec ces chiffres. Le canevas de l'iPhone est celui de la mise en page
(716 × 403), pas l'écran entier (1 170 × 2 532) : la charge de fragments plein écran est environ
dix fois plus grande, et la marge mesurée (rendu seul 3 ms sur 16,7) dit qu'elle tiendrait
probablement, pas qu'elle tient. Safari arrondit `performance.now()` à la milliseconde, d'où les
temps entiers. Aucun GPU dédié n'a été mesuré ; il n'est pas la cible.

### La semaine 9 tient encore

Le viewer a grandi (tramage, groupes, palette double) et l'en-tête du `.g3d` porte les
métadonnées de l'attribut. Grille de la semaine 9 refaite, dix passages par case :

| aperçu 3 Mb, Slow 4G | semaine 9 | semaine 10 |
|---|---:|---:|
| préfixe 16 Kio (défaut) | 1 412 ms | **1 416 ms** |
| préfixe ajusté | 1 381 ms | 1 381 ms |
| fin du premier rendu | 11 149 o | 11 619 o |
| page | 14,1 Kio | 15,1 Kio |

La marge passe de 88 à 84 ms. Elle aurait fondu sans deux corrections ci-dessous : la page avait
d'abord grossi à 18,7 Kio.

### Défauts trouvés

**Le lockfile de la semaine 9 était périmé.** `apps/first/package.json` déclarait deux
dépendances absentes de `pnpm-lock.yaml` ; `pnpm install --frozen-lockfile`, ce que fait toute
CI, échouait. Vu en reconstruisant l'environnement après le recyclage du conteneur.

**L'occlusion ambiante voyageait dans la page du premier rendu.** Un gabarit de chaîne avec
`${…}` au niveau d'un module n'est pas éliminé par esbuild, qui ne peut pas prouver l'absence
d'effet de bord ; et les commentaires des shaders survivent à la minification. 18,7 Kio au lieu
de 13,7 ; 14,8 après les deux corrections (shaders construits à la demande, commentaires GLSL
retirés au bundle). Les suites headless sont construites pareil, donc elles prouvent que les
shaders dépouillés compilent.

**La force du damier, mesurée comme d'habitude, trouvait des compartiments dans le bruit.**
Ci-dessus. Le test du témoin la garde en vue : il exige qu'elle soit fausse sur l'échantillon
même et juste une fois croisée.

**cooltools refusait les vraies données.** Les billes ont la même longueur dans un chromosome
mais pas d'un chromosome à l'autre (747 à 753 kb) ; cooltools exige une taille de bin unique.
L'ensemble jouet du test avait des billes toutes égales et l'avait laissé passer. Le cooler
intermédiaire est maintenant indexé par bille ; le test a des longueurs par chromosome et
échoue sur l'ancien code.

**`make fine` ne reproduisait pas la semaine 8.** La commande de reproduction publiée tournait
en 4 réplicats × 120 instantanés, quand le rapport en utilisait 8 × 150. Défauts corrigés.

**La semaine 8 n'est reproductible qu'en statistique.** OpenMM sur processeur somme les forces
dans un ordre qui dépend des fils d'exécution : à graines égales, les trajectoires divergent.
Deux exécutions de 1 200 conformations :

| | semaine 8 | semaine 10 |
|---|---:|---:|
| frontières retrouvées, avec CTCF | 12/13, 18 appelées | 12/13, 19 appelées |
| frontières retrouvées, CTCF inoccupé | **0/13** | **2/13** |
| points d'angle aux convergents, médiane | 2,56 | 2,56 |
| idem, CTCF inoccupé | 1,06 | 1,00 |
| pente de P(s), régime domaine | −1,08 | −1,12 |
| pire écart du raccord | 2,42× | 2,35× |
| pente de R(s) du modèle fin | −0,02 | +0,02 |

Les conclusions tiennent ; le « 0 sur 13 » du témoin était le chiffre d'une exécution, et la
suivante en donne 2. La pose du niveau fin a bougé avec son médoïde : RMSD 426,2 nm, échelle
optimale ×1,05 (semaine 9 : 394,7 nm, ×1,58). La lecture ne change pas — une échelle ne
rattrape presque rien, l'écart est de forme. L'ensemble de la semaine 7, lui, se reproduit : même
médoïde, graine 1125. Le verdict de la pose citait « 2,42× » en dur ; il cite maintenant la plage.

**Deux de mes seuils de test étaient arbitraires.** Le test du pli exigeait 0,15 d'écart sur un
pixel seul, quantifié par pas de 1/16 ; l'invariance d'échelle exigeait l'égalité à 3/255 près,
quand un échantillon vaut 16/255. Les deux ont échoué au premier passage. Ils mesurent
maintenant ce qui est défendable — une décroissance monotone vers le contact sur des fenêtres
moyennées, et le nombre de pixels qui basculent d'un échantillon — avec des témoins qui, eux,
échouent.

### Ce que ça ne dit pas

**Rien des fps plein écran sur téléphone** ni sur un portable bas de gamme : deux appareils, pas
une flotte. **Rien
des compartiments réels** : ceux du modèle sont une prédiction. **L'Xi n'est pas replié
autrement que l'Xa** : il est montré sans compartiment parce qu'il n'en a pas in vivo, pas parce
que le modèle l'aurait trouvé. **Le tramage se voit** pendant une transition — un grain fixe à
l'écran, sur quelques dizaines d'images — ; c'est le prix d'un fondu sans tri. **La force de
l'occlusion (0,9) est réglée à l'œil.** La plongée s'arrête à la bille de 750 kb : en dessous,
c'est la semaine 11.

### Reproduire

```console
$ cd pipeline && make ensemble && make fine     # ~32 + ~15 min, si les données manquent
$ make g3d                                      # .g3d avec compartiments du modèle
$ cd .. && pnpm verify                          # 16 assertions S4 + 20 assertions S10
$ node tools/test-g3d.mjs                       # lecture croisée, attributs compris (34)
$ node tools/dive.mjs 240                       # continuité, et le témoin qui doit sauter
$ node tools/loadtime.mjs 10                    # la semaine 9 tient-elle encore
```

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
