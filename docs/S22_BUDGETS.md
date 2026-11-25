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
