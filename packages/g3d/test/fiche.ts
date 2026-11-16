/** Lecture croisée du fichier de fiches : le lecteur TypeScript contre le lecteur Python.
 *  Exécuté sous Node par `tools/test-fiche.mjs` ; les attentes viennent de Python. */

import { readFileSync } from "node:fs";
import { BufferSource, CorruptError, FicheReader } from "../src/index.js";

const [, , filePath, expectedPath] = process.argv;
const bytes = new Uint8Array(readFileSync(filePath));
const expected = JSON.parse(readFileSync(expectedPath, "utf8")) as {
  cards: { level: string; bead: number; card: unknown }[];
  lifts: { chrom: string; pos: number; hit: { chrom: string; pos: number; strand: string } | null }[];
};

let passed = 0;
let failed = 0;
const check = (name: string, ok: boolean, detail = "") => {
  if (ok) passed++;
  else failed++;
  console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${name}${detail ? `  \x1b[2m${detail}\x1b[0m` : ""}`);
};

/** JSON canonique : clés triées, pour comparer des objets construits dans deux langages. */
const canon = (v: unknown): string => JSON.stringify(v, (_, x) =>
  x && typeof x === "object" && !Array.isArray(x) ? Object.fromEntries(Object.entries(x).sort(([a], [b]) => (a < b ? -1 : 1))) : x);

const r = await FicheReader.open(new BufferSource(bytes));
check("en-tête lu, empreinte vérifiée", r.header.format === "geno-fiche", `${r.header.levels.length} niveaux, ${Object.keys(r.header.sources).length} sources`);

let same = 0;
let firstDiff = "";
for (const e of expected.cards) {
  const c = await r.cardAsync(e.level, e.bead);
  if (canon(c) === canon(e.card)) same++;
  else if (!firstDiff) {
    const a = c as Record<string, unknown>, b = e.card as Record<string, unknown>;
    const k = Object.keys(b).find((key) => canon(a[key]) !== canon(b[key])) ?? "?";
    firstDiff = `${e.level} #${e.bead} champ ${k} : ts ${canon(a[k]).slice(0, 160)} · py ${canon(b[k]).slice(0, 160)}`;
  }
}
check("chaque fiche relue par TypeScript est celle de Python, champ pour champ", same === expected.cards.length,
      `${same}/${expected.cards.length}${firstDiff ? ` — ${firstDiff}` : ""}`);

const sync = expected.cards.every((e) => r.has(e.level, e.bead) && canon(r.card(e.level, e.bead)) === canon(e.card));
check("une fois ses blocs lus, `card` rend la fiche sans attendre", sync);

let liftOk = 0;
for (const l of expected.lifts) {
  const t = r.liftHs1(l.chrom, l.pos);
  if (canon(t) === canon(l.hit)) liftOk++;
}
check("T2T à la base : la chaîne embarquée lève comme la chaîne complète de Python", liftOk === expected.lifts.length,
      `${liftOk}/${expected.lifts.length} positions des régions fines`);
check("hors des régions fines, pas de T2T à la base (pas d'interpolation)", r.liftHs1("chr1", 1_000_000) === null);

// Témoin : un octet changé dans un bloc doit être vu.
const lv = r.header.levels[1];
const base = 64 + new DataView(bytes.buffer).getUint32(12, true);
const bad = bytes.slice();
bad[base + lv.blocks[0].offset + 10] ^= 0xff;
const rb = await FicheReader.open(new BufferSource(bad));
let caught = false;
try {
  await rb.cardAsync(lv.name, 1 + Array.from({ length: lv.n_beads }, (_, i) => i).find(() => true)!);
} catch (e) {
  caught = e instanceof CorruptError;
}
// La bille 1 n'est pas forcément dans le bloc 0 : on cherche une bille qui y est.
if (!caught) {
  await rb.openLevel(lv.name);
  for (let b = 1; b <= lv.n_beads && !caught; b++) {
    try {
      await rb.cardAsync(lv.name, b);
    } catch (e) {
      caught = e instanceof CorruptError;
    }
  }
}
check("témoin : un bloc corrompu est refusé", caught);

console.log(`\n  ${passed} réussis, ${failed} échoués\n`);
process.exit(failed ? 1 : 0);
