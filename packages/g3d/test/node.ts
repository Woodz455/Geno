/** Lecture croisée : un fichier écrit par Python, relu par le lecteur TypeScript.
 *
 *  Exécuté sous Node par `tools/test-g3d.mjs`. Les attentes viennent du lecteur Python
 *  sur le même fichier ; aucune n'est recopiée à la main.
 */

import { readFileSync } from "node:fs";
import { BufferSource, CorruptError, G3D, HttpSource, type Source, toNucleus } from "../src/index.js";

interface Expected {
  bytes: number;
  first_frame_end: number;
  levels: Record<
    string,
    {
      n: number;
      chunks: number;
      positions: number[];
      nucleus_frame: number[];
      copy: number[];
      start: number[];
      end: number[];
      variability: number[];
      variability_step: number;
      queries: { chrom: string; start: number; end: number; hits: { copy: string; ids: number[] }[] }[];
    }
  >;
}

const [, , filePath, expectedPath] = process.argv;
const bytes = new Uint8Array(readFileSync(filePath));
const expected: Expected = JSON.parse(readFileSync(expectedPath, "utf8"));

let failed = 0;
let passed = 0;
function check(name: string, ok: boolean, detail = ""): void {
  if (ok) passed++;
  else failed++;
  console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${name}${detail ? `  \x1b[2m${detail}\x1b[0m` : ""}`);
}

class Counting implements Source {
  reads: [number, number][] = [];
  constructor(private readonly inner: Source) {}
  read(offset: number, length: number) {
    this.reads.push([offset, length]);
    return this.inner.read(offset, length);
  }
}

const maxAbs = (a: ArrayLike<number>, b: ArrayLike<number>) => {
  let m = 0;
  for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
  return m;
};

async function main(): Promise<void> {
  console.log("\n  Lecture croisée .g3d — écrit par Python, relu par TypeScript\n");
  check("taille du fichier", bytes.length === expected.bytes, `${bytes.length} o`);

  // — Un préfixe qui couvre le premier rendu : une seule lecture en tout. —
  const one = new Counting(new BufferSource(bytes));
  const g = await G3D.open(one, { prefix: expected.first_frame_end });
  await g.firstFrame();
  check(
    "premier rendu en une seule lecture quand le préfixe le couvre",
    one.reads.length === 1 && g.requests === 1,
    `${one.reads.length} lecture(s)`,
  );

  // — Un préfixe trop court : une lecture de plus, pas davantage. —
  const two = new Counting(new BufferSource(bytes));
  const g2 = await G3D.open(two, { prefix: 4096 });
  await g2.firstFrame();
  check("préfixe trop court : exactement une lecture de plus", two.reads.length === 2, `${two.reads.length} lecture(s)`);

  for (const [name, exp] of Object.entries(expected.levels)) {
    const frame = await g.firstFrame(name);
    const lv = g.level(name);
    check(`${name} : ${exp.chunks} chunk(s), ${exp.n} billes`, lv.chunks.length === exp.chunks && frame.count === exp.n);
    check(`${name} : copies identiques`, frame.copy.every((c, i) => c === exp.copy[i]));
    // Positions : float32 côté TypeScript, float64 côté Python. Sur 10 µm, l'ulp d'un
    // float32 vaut ~0,001 nm ; tout écart au-delà serait une erreur de décodage.
    const dp = maxAbs(frame.positions, exp.positions);
    check(`${name} : positions au flottant près`, dp < 0.01, `écart max ${dp.toExponential(2)} nm`);
    const dn = maxAbs(toNucleus(lv, frame.positions), exp.nucleus_frame);
    check(`${name} : passage au repère du noyau`, dn < 0.01, `écart max ${dn.toExponential(2)} nm`);

    let idOk = true;
    let varErr = 0;
    for (const i of [0, 1, exp.n >> 1, exp.n - 1]) {
      const b = await g.identify(name, i + 1);
      idOk &&= b.start === exp.start[i] && b.end === exp.end[i];
      varErr = Math.max(varErr, Math.abs(b.variabilityNm - exp.variability[i]));
    }
    check(`${name} : identify rend début et fin exacts`, idOk);
    check(`${name} : variabilité à la précision du stockage`, varErr < 1e-6, `écart ${varErr.toExponential(1)} nm`);

    for (const q of exp.queries) {
      const hits = await g.locate(name, q.chrom, q.start, q.end);
      const same = JSON.stringify(hits) === JSON.stringify(q.hits);
      const n = hits.reduce((s, h) => s + h.ids.length, 0);
      check(`${name} : locate ${q.chrom}:${q.start}-${q.end}`, same, `${n} bille(s) sur ${hits.length} copie(s)`);
    }
  }

  // — Culling : un demi-espace qui contient tout garde tout, un autre en exclut. —
  const lv = g.level("noyau");
  const all = g.visibleChunks("noyau", [[1, 0, 0, 1e9]]);
  const xs = lv.chunks.map((c) => c.sphere[0]).sort((a, b) => a - b);
  const median = xs[xs.length >> 1];
  const half = g.visibleChunks("noyau", [[1, 0, 0, -median - 1e6]]);
  check("culling : tout dans un demi-espace qui contient tout", all.length === lv.chunks.length);
  check("culling : rien au-delà de toutes les sphères", half.length === 0);

  // — Corruption : un octet changé dans le premier rendu doit être vu. —
  const bad = bytes.slice();
  bad[expected.first_frame_end - 2] ^= 1;
  const gb = await G3D.open(new BufferSource(bad), { prefix: expected.first_frame_end });
  let caught = false;
  try {
    await gb.firstFrame();
  } catch (e) {
    caught = e instanceof CorruptError;
  }
  check("un octet altéré dans une colonne est détecté", caught);

  const badHeader = bytes.slice();
  badHeader[70] ^= 1;
  let headerCaught = false;
  try {
    await G3D.open(new BufferSource(badHeader));
  } catch (e) {
    headerCaught = e instanceof CorruptError;
  }
  check("un octet altéré dans l'en-tête est détecté", headerCaught);

  // — Par HTTP : le même fichier, servi avec Range par tools/serve.mjs. —
  const base = process.argv[4];
  if (base) {
    const http = new HttpSource(`${base}/g3d-fixture/fixture.g3d`);
    const gh = await G3D.open(http, { prefix: expected.first_frame_end });
    const fh = await gh.firstFrame();
    const ref = await g.firstFrame();
    check(
      "HTTP : premier rendu en une requête Range, réponse 206",
      http.log.length === 1 && http.log[0].status === 206 && http.log[0].received === expected.first_frame_end,
      `${http.log.length} requête(s), statut ${http.log[0]?.status}, ${http.log[0]?.received} o`,
    );
    check("HTTP : mêmes positions qu'en mémoire", maxAbs(fh.positions, ref.positions) === 0);
    await gh.identify("noyau", 5);
    check(
      "HTTP : identify ne relit qu'une plage, pas le fichier",
      http.log.length === 2 && http.log[1].received < bytes.length / 4,
      `${http.log[1]?.received} o sur ${bytes.length}`,
    );
    const past = await fetch(`${base}/g3d-fixture/fixture.g3d`, { headers: { Range: `bytes=${bytes.length + 10}-` } });
    check("HTTP : une plage hors du fichier rend 416", past.status === 416);
  }

  console.log(`\n  ${passed} réussis, ${failed} échoués\n`);
  process.exit(failed === 0 ? 0 : 1);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
