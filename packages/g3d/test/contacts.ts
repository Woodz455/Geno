/** Lecture croisée du fichier de contacts : le lecteur TypeScript contre le lecteur Python.
 *  Exécuté sous Node par `tools/test-contacts.mjs` ; les attentes viennent de Python. */

import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { BufferSource, ContactsReader, CorruptError } from "../src/index.js";

const [, , filePath, expectedPath] = process.argv;
const bytes = new Uint8Array(readFileSync(filePath));
const expected = JSON.parse(readFileSync(expectedPath, "utf8")) as {
  tiles: { uid: string; z: number; x: number; y: number; sha256: string; nan: number }[];
  pixels: { uid: string; z: number; row: number; col: number; value: number | null }[];
};

let passed = 0;
let failed = 0;
const check = (name: string, ok: boolean, detail = "") => {
  if (ok) passed++;
  else failed++;
  console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${name}${detail ? `  \x1b[2m${detail}\x1b[0m` : ""}`);
};

const r = await ContactsReader.open(new BufferSource(bytes));
check("en-tête lu, empreinte vérifiée", r.header.format === "geno-contacts/1",
      `${r.header.tilesets.map((t) => `${t.uid} ${t.resolutions.length} résolutions`).join(", ")}`);

let same = 0;
let first = "";
for (const e of expected.tiles) {
  const t = await r.tile(e.uid, e.z, e.x, e.y);
  const h = t ? createHash("sha256").update(new Uint8Array(t.buffer, t.byteOffset, t.byteLength)).digest("hex") : "absente";
  if (h === e.sha256) same++;
  else if (!first) first = `${e.uid} ${e.z}.${e.x}.${e.y}`;
}
check("chaque tuile décodée par TypeScript est celle de Python, au bit près (float16 → float32)",
      same === expected.tiles.length, `${same}/${expected.tiles.length}${first ? ` · d'abord ${first}` : ""}`);

let px = 0;
for (const e of expected.pixels) {
  const v = await r.pixel(e.uid, e.z, e.row, e.col);
  if ((e.value === null && Number.isNaN(v)) || (e.value !== null && Object.is(Math.fround(e.value), v))) px++;
}
check("chaque pixel relu (miroir et tuiles absentes compris) est celui de Python", px === expected.pixels.length,
      `${px}/${expected.pixels.length}`);

const abs = r.absolute("chr17", 7_673_801);
const back = abs === null ? null : r.relative(abs);
check("chr17:7 673 802 → position absolue de HiGlass → retour", back?.chrom === "chr17" && back.pos === 7_673_801, `${abs}`);

// Témoin : un octet corrompu dans une tuile doit être refusé.
const ts = r.header.tilesets[0];
const spec = Object.values(ts.levels[ts.levels.length - 1].tiles)[0];
const bad = bytes.slice();
const hlen = new DataView(bytes.buffer).getUint32(12, true);
bad[64 + hlen + spec.offset + 10] ^= 0xff;
const rb = await ContactsReader.open(new BufferSource(bad));
const [z, x, y] = Object.keys(ts.levels[ts.levels.length - 1].tiles)[0].split(".").map(Number);
let refused = false;
try {
  await rb.tile(ts.uid, z, x, y);
} catch (e) {
  refused = e instanceof CorruptError;
}
check("témoin : une tuile corrompue d'un octet est refusée", refused);

console.log(`\n${passed} passés, ${failed} échoués`);
process.exit(failed ? 1 : 0);
