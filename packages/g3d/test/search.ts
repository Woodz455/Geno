/** Lecture croisée de l'index de recherche : le lecteur TypeScript contre Python, sur les
 *  requêtes des critères de la semaine 16. Exécuté sous Node par `tools/test-search.mjs`. */

import { readFileSync } from "node:fs";
import { BufferSource, CorruptError, SearchIndex } from "../src/index.js";
import type { Resolution } from "../src/index.js";

const [, , filePath, expectedPath] = process.argv;
const bytes = new Uint8Array(readFileSync(filePath));
const expected = JSON.parse(readFileSync(expectedPath, "utf8")) as {
  queries: { group: string; query: string; expected: Resolution }[];
};

let passed = 0;
let failed = 0;
const check = (name: string, ok: boolean, detail = "") => {
  if (ok) passed++;
  else failed++;
  console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${name}${detail ? `  \x1b[2m${detail}\x1b[0m` : ""}`);
};

const key = (r: Resolution) =>
  JSON.stringify([r.kind, r.error, r.hits.map((h) => [h.label, h.chrom, h.start, h.end, h.in_model, h.sources ?? null])]);

const ix = await SearchIndex.open(new BufferSource(bytes));
check("en-tête lu, empreinte vérifiée", ix.header.format === "geno-search/1",
      `${ix.header.tables.genes.rows} clés de gènes, ${ix.header.tables.rs.distinct} rsID`);

const groups = new Map<string, { n: number; same: number; first: string }>();
for (const q of expected.queries) {
  const g = groups.get(q.group) ?? { n: 0, same: 0, first: "" };
  const got = await ix.resolve(q.query);
  g.n++;
  if (key(got) === key(q.expected)) g.same++;
  else if (!g.first) g.first = `${JSON.stringify(q.query)} → ${key(got)} ≠ ${key(q.expected)}`;
  groups.set(q.group, g);
}
for (const [name, g] of groups) {
  check(`${name} : la résolution TypeScript est celle de Python`, g.same === g.n, `${g.same}/${g.n}${g.first ? ` · d'abord ${g.first}` : ""}`);
}

const rs334 = await ix.resolve("rs334");
check("rs334 → chr11:5 227 002 (HBB)", rs334.hits.length === 1 && rs334.hits[0].chrom === "chr11" && rs334.hits[0].start === 5_227_001);

// Témoin : un octet corrompu dans un bloc de rsID doit être refusé.
const spec = ix.header.tables.rs.blocks[0];
const bad = bytes.slice();
const hlen = new DataView(bytes.buffer).getUint32(12, true);
bad[64 + hlen + spec.offset + 5] ^= 0xff;
const ib = await SearchIndex.open(new BufferSource(bad));
let refused = false;
try {
  await ib.rs(spec.first as number);
} catch (e) {
  refused = e instanceof CorruptError;
}
check("témoin : un bloc corrompu d'un octet est refusé", refused);

console.log(`\n${passed} passés, ${failed} échoués`);
process.exit(failed ? 1 : 0);
