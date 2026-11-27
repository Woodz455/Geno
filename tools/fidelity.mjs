/** Fidélité à la séquence : ce que la vue finale dessine, contre ce que dit le génome.
 *
 *  La page (`dive.html?fidelity`) rend la vue de la base, classe la couleur de chaque pixel
 *  de barreau par teinte et vote une lettre par base visible. Ici, la séquence de référence
 *  est relue par Python **directement dans le .2bit** — pas dans le .g3d, pas par le
 *  navigateur — et chaque lettre dessinée y est comparée : brin + tel quel, brin −
 *  complémenté. Le témoin décale la référence d'une base : l'accord doit s'effondrer.
 *
 *    node tools/fidelity.mjs
 */

import { execFileSync } from "node:child_process";
import { existsSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const GENOME = resolve(ROOT, "pipeline/data/raw/seq/hg38.2bit");
if (!existsSync(GENOME)) {
  console.error("hg38.2bit absent — `make data-seq` dans pipeline/.");
  process.exit(2);
}
const COMP = { A: "T", C: "G", G: "C", T: "A", N: "N" };

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
let r = null;
try {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Page.navigate", { url: `${base}/dive.html?f=data/gm12878.g3d&fidelity` });
  for (let i = 0; i < 900 && !r; i++) {
    await new Promise((x) => setTimeout(x, 200));
    const { result } = await s("Runtime.evaluate", { expression: "window.__geno ? JSON.stringify(window.__geno) : null", returnByValue: true });
    if (result.value) r = JSON.parse(result.value);
  }
} finally {
  browser.close();
  server.close();
  await shutdown(proc, profile);
}
if (!r?.letters) {
  console.error(`pas de lettres : ${JSON.stringify(r).slice(0, 300)}`);
  process.exit(1);
}
const keys = Object.keys(r.letters).map((k) => k.split(":").map(Number));
const ranks = keys.map(([k]) => k);
const lo = Math.min(...ranks) - 2, hi = Math.max(...ranks) + 3;
const start = r.region_start + lo;
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const ref = execFileSync(existsSync(venv) ? venv : "python3", ["-m", "geno_pipeline.sequence", GENOME, r.target.chrom, String(start), String(r.region_start + hi)],
  { cwd: resolve(ROOT, "pipeline"), env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") } }).toString().trim();
const expect = (rank, strand, shift = 0) => {
  const b = ref[rank - lo + shift];
  return strand === 0 ? b : COMP[b];
};
let agree = 0, shifted = 0;
const wrong = [];
for (const [rank, strand] of keys) {
  const got = r.letters[`${rank}:${strand}`];
  if (got === expect(rank, strand)) agree++;
  else if (wrong.length < 8) wrong.push(`${r.target.chrom}:${(r.region_start + rank + 1).toLocaleString("fr")}/${strand ? "−" : "+"} dessiné ${got}, attendu ${expect(rank, strand)}`);
  if (got === expect(rank, strand, 1)) shifted++;
}
const n = keys.length;
const tr = r.target_rank;
const tPlus = r.letters[`${tr}:0`], tMinus = r.letters[`${tr}:1`];
console.log(`\n  Fidélité à la séquence — vue finale (${r.target.name}), ${n} bases visibles, référence relue par Python dans hg38.2bit\n`);
console.log(`  accord               ${agree}/${n}  (${((100 * agree) / n).toFixed(1)} %)`);
console.log(`  témoin décalé d'une base  ${shifted}/${n}  (${((100 * shifted) / n).toFixed(1)} %)`);
console.log(`  base visée           ${r.target.chrom}:${r.target.position.toLocaleString("fr")} — dessinée ${tPlus ?? "cachée"} (+) / ${tMinus ?? "cachée"} (−), génome ${ref[tr - lo]} (+)`);
for (const w of wrong) console.log(`    ✗ ${w}`);
const ok = n >= 10 && agree === n && shifted / n < 0.5;
console.log(`\n  Verdict : ${ok ? "\x1b[32mchaque base dessinée est la base du génome\x1b[0m" : "\x1b[31mnon\x1b[0m"}\n`);
process.stdout.write(`JSON ${JSON.stringify({ n, agree, shifted, target: r.target, drawn: { plus: tPlus, minus: tMinus }, genome: ref[tr - lo] })}\n`);
process.exit(ok ? 0 : 1);
