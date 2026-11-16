/** Semaine 19 — le contenu du mode Découverte, vérifié hors navigateur.
 *
 *  Critères (fixés avant la première mesure) :
 *    1  structure   7 étapes dans l'ordre noyau → hélice, 4 haltes (HBB, FTO, HOXD, X inactif),
 *                   chaque étape avec un arrêt, un texte et une évidence connue ;
 *    2  glossaire   chaque [[terme]] a sa définition ; aucune définition orpheline ;
 *    3  jargon      sur chaque écran (introduction, chaque étape, chaque halte), la **première**
 *                   occurrence de chaque terme du lexique ci-dessous porte sa définition — 0 écart ;
 *    4  échelles    chaque valeur citée, recalculée ici par un autre chemin, à 1 % près ; aucune
 *                   marque {clé} inconnue ;
 *    5  lisibilité  indice de Kandel et Moles par écran : rapporté, sans seuil.
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
const C = JSON.parse(readFileSync(resolve(ROOT, "apps/discover/content/fr.json"), "utf8"));
const G = JSON.parse(readFileSync(resolve(ROOT, "apps/discover/content/glossary.fr.json"), "utf8"));

/** Le lexique du jargon, pré-enregistré : terme du glossaire → formes à chercher. */
const LEXICON = {
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

const fails = [];
const ok = (cond, msg) => { if (!cond) fails.push(msg); return cond; };

// 1 — structure
const ORDER = ["noyau", "territoires", "compartiments", "tad", "boucle", "fibre", "helice"];
ok(JSON.stringify(C.steps.map((s) => s.id)) === JSON.stringify(ORDER), `étapes : ${C.steps.map((s) => s.id)}`);
ok(JSON.stringify(C.halts.map((h) => h.id)) === JSON.stringify(["hbb", "fto", "hoxd", "xi"]), "haltes");
for (const s of C.steps) ok(s.stop && s.text && s.show && ["mesuré", "simulé", "prédit", "illustratif"].includes(s.evidence), `étape ${s.id} incomplète`);
for (const h of C.halts) ok(h.gene && h.text && h.limits, `halte ${h.id} incomplète`);
ok(C.quiz.choices.length === 3 && C.quiz.choices.every((c) => c.answer >= 0 && c.answer < c.options.length), "quiz");

// Les écrans, et leur texte balisé.
const screens = [
  { id: "intro", texts: [C.intro] },
  ...C.steps.map((s) => ({ id: s.id, texts: [s.text, s.how, s.summary, s.show].filter(Boolean) })),
  ...C.halts.map((h) => ({ id: h.id, texts: [h.text, h.limits] })),
];

// 2 — glossaire
const TAG = /\[\[([^\]|]+)(?:\|([^\]]+))?\]\]/g;
const used = new Set();
for (const sc of screens) for (const t of sc.texts) for (const m of t.matchAll(TAG)) {
  used.add(m[1]);
  ok(m[1] in G, `${sc.id} : [[${m[1]}]] sans définition`);
}
const orphans = Object.keys(G).filter((k) => !used.has(k));
ok(orphans.length === 0, `définitions jamais appelées : ${orphans.join(", ")}`);

// 3 — jargon : la première occurrence de chaque terme, sur chaque écran, est balisée.
let checked = 0;
const untagged = [];
for (const sc of screens) {
  const raw = sc.texts.join("\n");
  // Les zones balisées, en positions du texte brut.
  const zones = [...raw.matchAll(TAG)].map((m) => [m.index, m.index + m[0].length]);
  for (const [term, forms] of Object.entries(LEXICON)) {
    const re = new RegExp(`(?<![\\p{L}\\p{N}])(?:${forms.join("|")})(?![\\p{L}\\p{N}])`, "u");
    const m = re.exec(raw);
    if (!m) continue;
    checked++;
    if (!zones.some(([a, b]) => m.index >= a && m.index < b)) untagged.push(`${sc.id} : « ${m[0]} »`);
  }
}
ok(untagged.length === 0, `jargon sans définition à sa première occurrence : ${untagged.join(" · ")}`);

// 4 — échelles, recalculées par un autre chemin.
const v = values();
const num = (s) => Number(s.replace(/\s| | /g, "").replace(",", "."));
const L = 6.2e9 * 0.34e-9;
const EXPECT = {
  dna_m: L, nucleus_um: 10, fold: L / 10e-6, ball_dna_km: (L * 0.067) / 10e-6 / 1000,
  ball_helix_um: 2e-9 * (0.067 / 10e-6) * 1e6, bead_kb: 750, bead_novels: 750000 / 540000,
  genome_novels: 3.1e9 / 540000, tad_kb: 185, tad_um: 185000 * 0.34e-3, monomer_kb: 2, nuc_bp: 147,
  nuc_turns: 1.65, rise_nm: 0.34, helix_nm: 2,
};
const scale = [];
for (const [k, want] of Object.entries(EXPECT)) {
  const got = num(v[k]);
  const rel = Math.abs(got - want) / want;
  scale.push(`${k} ${v[k]}`);
  ok(rel <= 0.01, `${k} : ${v[k]} au lieu de ${want}`);
}
ok(Object.keys(v).length === Object.keys(EXPECT).length, "valeurs d'échelle non recalculées");
for (const sc of screens) for (const t of sc.texts) { try { fill(t); } catch (e) { ok(false, `${sc.id} : ${e.message}`); } }
ok(FACTS.risePerBp === 0.34e-9, "montée par paire");

// 5 — lisibilité (Kandel et Moles 1958, adaptation de Flesch au français).
const syll = (w) => Math.max(1, (w.toLowerCase().replace(/e$/, "").match(/[aeiouyàâäéèêëîïôöùûü]+/g) ?? []).length);
const plain = (t) => fill(t).replace(TAG, (_, a, b) => b ?? a);
const km = (t) => {
  const s = plain(t), words = s.match(/[\p{L}\p{N}'’-]+/gu) ?? [], sentences = Math.max(1, (s.match(/[.!?:;]+(\s|$)/g) ?? []).length);
  const sy = words.reduce((a, w) => a + syll(w), 0);
  return { words: words.length, km: 207 - 1.015 * (words.length / sentences) - 73.6 * (sy / words.length) };
};
const read = screens.map((sc) => ({ id: sc.id, ...km(sc.texts.join(" ")) }));

console.log("\n  Mode Découverte — contenu (semaine 19)");
console.log(`  1  structure   ${C.steps.length} étapes, ${C.halts.length} haltes, quiz de ${C.quiz.choices.length} questions`);
console.log(`  2  glossaire   ${Object.keys(G).length} définitions, ${used.size} appelées`);
console.log(`  3  jargon      ${checked} premières occurrences contrôlées sur ${screens.length} écrans, ${untagged.length} sans définition`);
console.log(`  4  échelles    ${scale.length} valeurs recalculées à 1 % : ${scale.join(" · ")}`);
console.log(`  5  lisibilité  Kandel-Moles (60–70 : standard ; plus haut : plus facile) : ${read.map((r) => `${r.id} ${r.km.toFixed(0)} (${r.words} mots)`).join(" · ")}`);
if (fails.length) console.log(`\n  \x1b[31méchec\x1b[0m\n    ${fails.join("\n    ")}`);
else console.log("\n  verdict : \x1b[32mfait\x1b[0m");
process.exit(fails.length ? 1 : 0);
