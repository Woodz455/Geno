# Semaine 22 — budgets de performance, fixés avant toute mesure

Commité avant la première mesure de la semaine. Une modification ultérieure sera datée ici et
dite dans `VALIDATION.md` § S22, jamais faite en silence.

## Les trois cibles (décision de l'utilisateur, 2026-11-25)

| Cible | Appareil joué |
|---|---|
| téléphone | iPhone 14 (Apple GPU, Safari) |
| GPU intégré de portable | PC à Intel Iris Xe |
| GPU de bureau | une machine à GPU dédié |

Les temps viennent **uniquement** de la page de banc matériel (`apps/dive/materiel.ts`), jouée par
l'utilisateur sur chaque appareil. Ici, Chromium rend par SwiftShader, et aucun temps en ms mesuré
ici n'est jugé.

## Budgets sur matériel, par cible

La mesure : la plongée entière, du noyau à la base c.818 de TP53, en 240 images, dans la passe
rythmée par `requestAnimationFrame`, sans relecture.

| Budget | Seuil |
|---|---|
| image | p95 de l'intervalle rAF ≤ **20 ms** (le seuil déjà écrit dans la page de banc) |
| mémoire GPU allouée par la page | ≤ **256 Mo**, comptée par l'instrumentation (somme des `bufferData`, `texImage2D`, `texStorage2D`) |

## Budgets en CI, déterministes

Ils sont mesurés sous SwiftShader par un `WebGL2RenderingContext` instrumenté, sur les mêmes 240
images. Ils ne dépendent pas de la machine.

| Budget | Mesure |
|---|---|
| appels de dessin | maximum par image |
| instances | maximum par image |
| sommets soumis | maximum par image |
| pixels des passes plein écran | maximum par image (proxy des fragments de l'occlusion ambiante) |
| octets GPU alloués | en tout |
| pages de `dist/` | taille de chaque page |
| premier rendu | octets lus (préfixe de 16 Kio, S9) |

**Règle de cliquet** : chaque budget CI vaut la valeur mesurée après correction, plus 10 %. Un
commit qui la dépasse fait échouer `node tools/perf.mjs --check`, donc la CI. Un budget ne se
relève qu'à la main, par un commit qui le dit.

**Témoin** : une fenêtre de nucléosomes doublée (`CONTEXT = 80` dans `apps/dive/fiber.ts`) doit
faire échouer le contrôle, avec le budget dépassé nommé. Sans ce témoin, le contrôle n'est pas
réputé fonctionner.

## Le diagnostic, puis le choix de la correction

**Passage 1**, sur les trois cibles. Le segment lent, de 59 à 92 % du trajet, est joué sous quatre
variantes :
- qualité pleine ;
- occlusion ambiante coupée ;
- occlusion ambiante à demi-résolution ;
- fenêtre de nucléosomes et d'ADN réduite à ce qui est visible.

La plongée entière est ensuite jouée en qualité pleine.

**Règle de choix**, fixée ici :
- parmi les variantes (ou combinaisons de variantes) qui tiennent 20 ms au p95 sur les trois
  cibles, on retient celle qui s'écarte le moins de la qualité pleine ;
- l'écart est la moyenne absolue RGB sur les 240 images de la plongée, mesurée en headless ;
- à écart égal, on retient la plus simple ;
- couper l'occlusion ambiante n'est retenu que si aucune autre variante ne tient.

**Passage 2**, sur les trois cibles : la plongée entière avec la correction, jugée contre les deux
budgets matériel.

### Révision du 2026-11-25 — une cinquième variante, trouvée par l'instrumentation avant le passage 1

Datée et commitée **avant** tout passage sur matériel ; les seuils ne changent pas.

**Ce que les compteurs ont trouvé.** La mesure de référence, sous SwiftShader et sans aucun
temps, donne **27,3 millions de sommets** soumis par image au pire moment de la plongée, à
l'échelle de l'hélice. Un seul appel en soumet 19,3 millions : 1 536 indices × 12 544 segments du
tube de la double hélice.

**La cause.** `TubeRenderer.classify` (`packages/viewer/src/tube.ts`) n'écarte aucun segment hors
du champ. Pire, un segment derrière la caméra voit sa profondeur ramenée à 10⁻⁶ : il paraît
immense et part au niveau de détail maximal. Aucun de ces sommets ne fait un pixel.

**La variante ajoutée** : l'élimination des segments de tube hors du cône de vue (sphère
englobante du segment contre les quatre plans latéraux et le plan proche). Les compteurs, sur les
240 images de la plongée :

| Variante | Sommets par image (max) |
|---|---:|
| référence | 27,3 M |
| élimination hors champ | 2,49 M |
| fenêtre visible | 13,6 M |
| élimination + fenêtre visible | 2,05 M |

Le passage 1 joue donc cinq variantes et leurs combinaisons utiles :
- pleine ;
- élimination ;
- occlusion coupée ;
- occlusion à demi-résolution ;
- fenêtre visible.

La règle de choix ne change pas : la moins chère en écart d'image parmi celles qui tiennent.
L'élimination ne devrait rien changer à l'image ; c'est mesuré, pas supposé
(`tools/quality-diff.mjs`).

### Révision du 2026-11-25 (seconde) — l'écart de l'élimination, mesuré et isolé

Datée et commitée **avant** tout passage sur matériel. Les seuils et la règle de choix ne
changent pas.

**L'élimination change l'image, un peu.** Écart moyen 0,194 sur 255, max 1,34, avec 151 images
identiques sur 241. Elle n'était donc pas neutre, contrairement à ce qu'on supposait plus haut.
Écart de chaque variante contre la qualité pleine (moyenne absolue RGB sur 255, 241 images,
960 × 540, `tools/quality-diff.mjs`) :

| Variante | écart moyen | max | images identiques |
|---|---:|---:|---:|
| élimination, détail inchangé (ajoutée ici) | **0,0000** | 0,005 | 240 |
| élimination hors champ | 0,194 | 1,34 | 151 |
| occlusion à demi-résolution | 0,712 | 3,63 | 0 |
| fenêtre visible | 0,714 | 8,18 | 180 |
| élimination + occlusion moitié | 0,816 | 3,79 | 0 |
| élimination + fenêtre visible | 0,899 | 8,40 | 135 |
| occlusion coupée | 2,903 | 9,81 | 0 |

**La cause, isolée par une expérience.** Le nombre de côtés du tube suit le plus grand rayon
projeté parmi les segments classés. Sans élimination, ce rayon compte les segments **derrière** la
caméra, dont la profondeur est ramenée à 10⁻⁶. Le tube entier prend alors 16 côtés. Avec
l'élimination, ces segments ne comptent plus, et le tube visible peut en prendre moins.

Garder ce calcul tel quel, en éliminant quand même, rend l'image inchangée : 240 images
identiques sur 241, max 0,005. L'élimination elle-même est donc neutre. L'écart de 0,19 est
celui du maillage.

**D'où une variante de plus**, `cull: "same"` : élimination, détail calculé comme avant. Ses
compteurs, sur les 240 images :

| Réglage | sommets par image (max) | appels | instances |
|---|---:|---:|---:|
| référence | 27 295 758 | 28 | 36 886 |
| élimination, détail inchangé | 3 251 886 | 23 | 21 868 |
| élimination | 2 486 022 | 23 | 21 868 |

Le diagnostic du banc (v6) joue ces réglages, puis la plongée entière en qualité pleine :
- pleine ;
- élimination, détail inchangé ;
- élimination ;
- occlusion coupée ;
- occlusion moitié ;
- fenêtre visible ;
- élimination + fenêtre ;
- élimination + occlusion moitié.

Par la règle, si « détail inchangé » tient 20 ms sur les trois cibles, elle est retenue, car son
écart est nul.

Précision : le segment joué par le banc, depuis la v5, va de 55 à 95 % du trajet. C'est le
segment lent de 59 à 92 %, avec une marge de chaque côté.

### Cliquet provisoire (2026-11-25)

Pour que la CI tourne dès maintenant, `tools/perf-budgets.json` est inscrit à la qualité de
**référence**, plus 10 %. C'est un garde-fou contre la régression, pas le budget de la semaine.
Après le passage 1, le cliquet est réinscrit sur la variante retenue, donc abaissé, par un commit
qui le dit.

## Compatibilité

- Sans WebGL2, chaque page qui ouvre la scène montre un **repli** : un message qui nomme ce qui
  manque, et une vue 2D statique de l'aperçu. Jamais un écran noir, jamais une exception.
- Sans rendu en entiers 32 bits, le picking se désactive et la page le dit.
- Vérifié par Chromium lancé avec `--disable-webgl` puis `--disable-3d-apis` : 100 % des pages
  montrent le repli.
- WebGPU n'est pas porté. Sa présence est détectée et notée, sans être utilisée.

## Hors budget S22, dit

La fiche de 6,5 Mo dépasse son budget de 2,5 Mo (S13). Elle n'est pas réparée cette semaine. La CI
la garde sous cliquet : elle ne peut pas grossir.

## Ce qui fait la semaine

« Les budgets tiennent sur les trois cibles, et un commit qui les dépasse casse la CI. »
- **Faite** quand le passage 2 tient les deux budgets matériel sur les trois cibles, que la CI
  tourne sur GitHub avec ses budgets, et que le témoin la fait échouer.
- Si un passage matériel manque : « correction faite, verdict matériel en attente ».
- Si GitHub refuse le workflow : « contrôle local seulement ».
