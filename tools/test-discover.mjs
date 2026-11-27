/** Semaines 19–20 — le contenu du mode Découverte, vérifié hors navigateur, en français et en anglais.
 *
 *  Critères (S19, fixés avant la première mesure ; appliqués à chaque langue en S20) :
 *    1  structure   7 étapes dans l'ordre noyau → hélice, 4 haltes (HBB, FTO, HOXD, X inactif),
 *                   chaque étape avec un arrêt, un texte et une évidence connue ;
 *    2  glossaire   chaque [[terme]] a sa définition ; aucune définition orpheline ;
 *    3  jargon      sur chaque écran (introduction, chaque étape, chaque halte), la **première**
 *                   occurrence de chaque terme du lexique ci-dessous porte sa définition — 0 écart ;
 *    4  échelles    chaque valeur citée, recalculée ici par un autre chemin, à 1 % près ; aucune
 *                   marque {clé} inconnue ;
 *    5  lisibilité  indice de Kandel et Moles par écran (français), de Flesch (anglais) : rapporté,
 *                   sans seuil.
 *
 *  Ajoutés en S20, fixés avant la première exécution :
 *    6  parité      l'anglais a la même structure que le français : mêmes étapes, arrêts, évidences,
 *                   cartes de contacts, haltes et gènes, mêmes clés {…} dans chaque texte, même quiz
 *                   (nombre d'options, bonne réponse) ; même nombre de définitions ;
 *    7  interface   ui.fr.json et ui.en.json ont les mêmes clés, les mêmes {marques} dans chaque
 *                   chaîne, et les quatre phrases d'évidence de la scène (apps/dive/dive.ts) sont
 *                   toutes reconnues, à l'identique, par ui.fr.json ;
 *    8  nombres     chaque valeur d'échelle anglaise est égale à la française à 1 % près.
 *
 *    node tools/test-discover.mjs
 */

import { build } from "esbuild";
import { readFileSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/discover-test");
mkdirSync(OUT, { recursive: true });
await build({ entryPoints: [resolve(ROOT, "apps/discover/scale.ts")], bundle: true, format: "esm", platform: "node", outfile: resolve(OUT, "scale.mjs"), logLevel: "error" });
const { values, fill, FACTS } = await import(pathToFileURL(resolve(OUT, "scale.mjs")).href);
const json = (f) => JSON.parse(readFileSync(resolve(ROOT, "apps/discover/content", f), "utf8"));

/** Le lexique du jargon, pré-enregistré : terme du glossaire → formes à chercher. */
const LEXICON_FR = {
  "ADN": ["ADN"],
  "génome": ["génomes?"],
  "chromosome": ["chromosomes?"],
  "noyau": ["noyaux?"],
  "GM12878": ["GM12878"],
  "territoire chromosomique": ["territoires?"],
  "compartiment": ["compartiments?"],
  "gène": ["gènes?"],
  "TAD": ["TAD"],
  "contact": ["contacts?"],
  "cohésine": ["cohésines?"],
  "extrusion de boucle": ["extrusion"],
  "CTCF": ["CTCF"],
  "enhancer": ["enhancers?"],
  "promoteur": ["promoteurs?"],
  "nucléosome": ["nucléosomes?"],
  "chromatine": ["chromatine"],
  "double hélice": ["double hélice"],
  "base": ["bases?"],
  "mutation": ["mutations?"],
  "frontière": ["frontières?"],
  "X inactif": ["X inactif"],
  "superdomaine": ["superdomaines?"],
  "simulé": ["simulée?s?"],
  "prédit": ["prédite?s?"],
  "planté": ["plantée?s?"],
  "illustratif": ["illustrati(?:f|ve)s?"],
  "mesuré": ["mesurée?s?"],
};
/** Semaine 20 : le même lexique en anglais, pré-enregistré avant la première exécution. */
const LEXICON_EN = {
  "DNA": ["DNA"],
  "genome": ["genomes?"],
  "chromosome": ["chromosomes?"],
  "nucleus": ["nucleus", "nuclei"],
  "GM12878": ["GM12878"],
  "chromosome territory": ["territor(?:y|ies)"],
  "compartment": ["compartments?"],
  "gene": ["genes?"],
  "TAD": ["TADs?"],
  "contact": ["contacts?"],
  "cohesin": ["cohesins?"],
  "loop extrusion": ["extrusion"],
  "CTCF": ["CTCF"],
  "enhancer": ["enhancers?"],
  "promoter": ["promoters?"],
  "nucleosome": ["nucleosomes?"],
  "chromatin": ["chromatin"],
  "double helix": ["double helix"],
  "base": ["bases?"],
  "mutation": ["mutations?"],
  "border": ["borders?"],
  "inactive X": ["inactive X"],
  "superdomain": ["superdomains?"],
  "simulated": ["simulated"],
  "predicted": ["predicted"],
  "planted": ["planted"],
  "illustrative": ["illustrative"],
  "measured": ["measured"],
};

const LANGS = [
  { lang: "fr", C: json("fr.json"), G: json("glossary.fr.json"), UI: json("ui.fr.json"), LEX: LEXICON_FR, locale: "fr-FR" },
  { lang: "en", C: json("en.json"), G: json("glossary.en.json"), UI: json("ui.en.json"), LEX: LEXICON_EN, locale: "en-GB" },
];

const fails = [];
const ok = (cond, msg) => { if (!cond) fails.push(msg); return cond; };
const TAG = /\[\[([^\]|]+)(?:\|([^\]]+))?\]\]/g;
const KEYS = (t) => [...(t ?? "").matchAll(/\{([a-z_]+)\}/g)].map((m) => m[1]).sort().join(",");
const ORDER = ["noyau", "territoires", "compartiments", "tad", "boucle", "fibre", "helice"];
const report = {};

for (const { lang, C, G, LEX, locale } of LANGS) {
  const L = (m) => `[${lang}] ${m}`;
  // 1 — structure
  ok(JSON.stringify(C.steps.map((s) => s.id)) === JSON.stringify(ORDER), L(`étapes : ${C.steps.map((s) => s.id)}`));
  ok(JSON.stringify(C.halts.map((h) => h.id)) === JSON.stringify(["hbb", "fto", "hoxd", "xi"]), L("haltes"));
  for (const s of C.steps) ok(s.stop && s.text && s.show && ["mesuré", "simulé", "prédit", "illustratif"].includes(s.evidence), L(`étape ${s.id} incomplète`));
  for (const h of C.halts) ok(h.gene && h.text && h.limits, L(`halte ${h.id} incomplète`));
  ok(C.quiz.choices.length === 3 && C.quiz.choices.every((c) => c.answer >= 0 && c.answer < c.options.length), L("quiz"));

  const screens = [
    { id: "intro", texts: [C.intro] },
    ...C.steps.map((s) => ({ id: s.id, texts: [s.text, s.how, s.summary, s.show].filter(Boolean) })),
    ...C.halts.map((h) => ({ id: h.id, texts: [h.text, h.limits] })),
  ];

  // 2 — glossaire
  const used = new Set();
  for (const sc of screens) for (const t of sc.texts) for (const m of t.matchAll(TAG)) {
    used.add(m[1]);
    ok(m[1] in G, L(`${sc.id} : [[${m[1]}]] sans définition`));
  }
  const orphans = Object.keys(G).filter((k) => !used.has(k));
  ok(orphans.length === 0, L(`définitions jamais appelées : ${orphans.join(", ")}`));
  ok(Object.keys(LEX).every((k) => k in G), L(`lexique hors du glossaire : ${Object.keys(LEX).filter((k) => !(k in G)).join(", ")}`));

  // 3 — jargon : la première occurrence de chaque terme, sur chaque écran, est balisée.
  let checked = 0;
  const untagged = [];
  for (const sc of screens) {
    const raw = sc.texts.join("\n");
    const zones = [...raw.matchAll(TAG)].map((m) => [m.index, m.index + m[0].length]);
    for (const [term, forms] of Object.entries(LEX)) {
      const re = new RegExp(`(?<![\\p{L}\\p{N}])(?:${forms.join("|")})(?![\\p{L}\\p{N}])`, "u");
      const m = re.exec(raw);
      if (!m) continue;
      checked++;
      if (!zones.some(([a, b]) => m.index >= a && m.index < b)) untagged.push(`${sc.id} : « ${m[0]} »`);
    }
  }
  ok(untagged.length === 0, L(`jargon sans définition à sa première occurrence : ${untagged.join(" · ")}`));

  // 4 — échelles, recalculées par un autre chemin (nombres lus dans la mise en forme de la langue).
  const v = values(locale);
  const num = (s) => Number(lang === "fr" ? s.replace(/\s| | /g, "").replace(",", ".") : s.replace(/,/g, ""));
  const Ln = 6.2e9 * 0.34e-9;
  const EXPECT = {
    dna_m: Ln, nucleus_um: 10, fold: Ln / 10e-6, ball_dna_km: (Ln * 0.067) / 10e-6 / 1000,
    ball_helix_um: 2e-9 * (0.067 / 10e-6) * 1e6, bead_kb: 750, bead_novels: 750000 / 540000,
    genome_novels: 3.1e9 / 540000, tad_kb: 185, tad_um: 185000 * 0.34e-3, monomer_kb: 2, nuc_bp: 147,
    nuc_turns: 1.65, rise_nm: 0.34, helix_nm: 2,
  };
  const scale = [];
  for (const [k, want] of Object.entries(EXPECT)) {
    const got = num(v[k]);
    scale.push(`${k} ${v[k]}`);
    ok(Math.abs(got - want) / want <= 0.01, L(`${k} : ${v[k]} au lieu de ${want}`));
  }
  ok(Object.keys(v).length === Object.keys(EXPECT).length, L("valeurs d'échelle non recalculées"));
  for (const sc of screens) for (const t of sc.texts) { try { fill(t, v); } catch (e) { ok(false, L(`${sc.id} : ${e.message}`)); } }

  // 5 — lisibilité : Kandel et Moles (1958) en français, Flesch (1948) en anglais.
  const syll = (w) => {
    const x = w.toLowerCase();
    if (lang === "fr") return Math.max(1, (x.replace(/e$/, "").match(/[aeiouyàâäéèêëîïôöùûü]+/g) ?? []).length);
    return Math.max(1, (x.replace(/(?:[^laeiouy]es|ed|[^laeiouy]e)$/, "").replace(/^y/, "").match(/[aeiouy]{1,2}/g) ?? []).length);
  };
  const plain = (t) => fill(t, v).replace(TAG, (_, a, b) => b ?? a);
  const read = screens.map((sc) => {
    const s = plain(sc.texts.join(" ")), words = s.match(/[\p{L}\p{N}'’-]+/gu) ?? [], sentences = Math.max(1, (s.match(/[.!?:;]+(\s|$)/g) ?? []).length);
    const sy = words.reduce((a, w) => a + syll(w), 0);
    const score = lang === "fr" ? 207 - 1.015 * (words.length / sentences) - 73.6 * (sy / words.length)
                                : 206.835 - 1.015 * (words.length / sentences) - 84.6 * (sy / words.length);
    return { id: sc.id, words: words.length, score };
  });
  report[lang] = { C, G, screens, used, checked, untagged, scale, read, v };
}
ok(FACTS.risePerBp === 0.34e-9, "montée par paire");

// 6 — parité de structure FR / EN.
const [FR, EN] = [report.fr.C, report.en.C];
let parity = 0;
const same = (a, b, msg) => { parity++; ok(JSON.stringify(a) === JSON.stringify(b), `parité : ${msg}`); };
FR.steps.forEach((s, k) => {
  const e = EN.steps[k] ?? {};
  same([s.id, s.stop, s.evidence, !!s.contacts, !!s.how, !!s.summary], [e.id, e.stop, e.evidence, !!e.contacts, !!e.how, !!e.summary], `étape ${s.id}`);
  for (const f of ["text", "how", "summary", "show"]) same(KEYS(s[f]), KEYS(e[f]), `${s.id}.${f} : clés {…}`);
});
FR.halts.forEach((h, k) => {
  const e = EN.halts[k] ?? {};
  same([h.id, h.gene], [e.id, e.gene], `halte ${h.id}`);
  for (const f of ["text", "limits"]) same(KEYS(h[f]), KEYS(e[f]), `${h.id}.${f} : clés {…}`);
});
same(KEYS(FR.intro), KEYS(EN.intro), "intro : clés {…}");
same(FR.quiz.choices.map((c) => [c.options.length, c.answer]), EN.quiz.choices.map((c) => [c.options.length, c.answer]), "quiz");
same(Object.keys(report.fr.G).length, Object.keys(report.en.G).length, "nombre de définitions");

// 7 — interface : mêmes clés, mêmes marques ; les phrases d'évidence de la scène reconnues.
const flat = (o, p = "") => Object.entries(o).flatMap(([k, x]) => (x && typeof x === "object" ? flat(x, `${p}${k}.`) : [[`${p}${k}`, String(x)]]));
const [UF, UE] = [new Map(flat(LANGS[0].UI)), new Map(flat(LANGS[1].UI))];
ok(JSON.stringify([...UF.keys()].sort()) === JSON.stringify([...UE.keys()].sort()), `interface : clés différentes : ${[...UF.keys()].filter((k) => !UE.has(k)).concat([...UE.keys()].filter((k) => !UF.has(k))).join(", ")}`);
const marks = (t) => [...t.matchAll(/\{(\w+)\}/g)].map((m) => m[1]).sort().join(",");
const uiMarks = [...UF].filter(([k, t]) => UE.has(k) && marks(t) !== marks(UE.get(k))).map(([k]) => k);
ok(uiMarks.length === 0, `interface : marques {…} différentes : ${uiMarks.join(", ")}`);
const dive = readFileSync(resolve(ROOT, "apps/dive/dive.ts"), "utf8");
const e0 = dive.indexOf("evidence(): string {"), block = dive.slice(e0, dive.indexOf("highlight(regions", e0));
const scenePhrases = [...block.matchAll(/return "([^"]+)"/g)].map((m) => m[1]);
const known = Object.values(LANGS[0].UI.evidenceDetail);
const unknown = scenePhrases.filter((p) => !known.includes(p));
ok(scenePhrases.length === 4 && unknown.length === 0, `évidence : ${scenePhrases.length} phrases de la scène, inconnues : ${unknown.join(" | ")}`);

// 8 — nombres EN = FR.
const nfr = (s) => Number(s.replace(/\s| | /g, "").replace(",", ".")), nen = (s) => Number(s.replace(/,/g, ""));
const diffs = Object.keys(report.fr.v).filter((k) => Math.abs(nfr(report.fr.v[k]) - nen(report.en.v[k])) > 0.01 * Math.abs(nfr(report.fr.v[k])));
ok(diffs.length === 0, `nombres : ${diffs.map((k) => `${k} ${report.fr.v[k]} ≠ ${report.en.v[k]}`).join(", ")}`);

console.log("\n  Mode Découverte — contenu (semaines 19–20)");
for (const lang of ["fr", "en"]) {
  const r = report[lang];
  console.log(`  [${lang}]`);
  console.log(`  1  structure   ${r.C.steps.length} étapes, ${r.C.halts.length} haltes, quiz de ${r.C.quiz.choices.length} questions`);
  console.log(`  2  glossaire   ${Object.keys(r.G).length} définitions, ${r.used.size} appelées`);
  console.log(`  3  jargon      ${r.checked} premières occurrences contrôlées sur ${r.screens.length} écrans, ${r.untagged.length} sans définition`);
  console.log(`  4  échelles    ${r.scale.length} valeurs recalculées à 1 % : ${r.scale.join(" · ")}`);
  console.log(`  5  lisibilité  ${lang === "fr" ? "Kandel-Moles" : "Flesch"} (60–70 : standard ; plus haut : plus facile) : ${r.read.map((x) => `${x.id} ${x.score.toFixed(0)} (${x.words} mots)`).join(" · ")}`);
}
console.log(`  6  parité      ${parity} comparaisons FR/EN`);
console.log(`  7  interface   ${UF.size} chaînes par langue · ${scenePhrases.length} phrases d'évidence de la scène reconnues : ${scenePhrases.length - unknown.length}`);
console.log(`  8  nombres     ${Object.keys(report.fr.v).length} valeurs EN = FR à 1 %, écarts : ${diffs.length}`);
if (fails.length) console.log(`\n  \x1b[31méchec\x1b[0m\n    ${fails.join("\n    ")}`);
else console.log("\n  verdict : \x1b[32mfait\x1b[0m");
process.exit(fails.length ? 1 : 0);
