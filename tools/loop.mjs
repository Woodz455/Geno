/** Lisibilité de la boucle CTCF visée (semaine 11), et le témoin sans isolement.
 *
 *  Mesurée sur le tampon d'identifiants de la vue de boucle : pixels visibles des deux
 *  ancres, de l'anneau de cohésine, de l'arc ; part du corps de la boucle visible ;
 *  orientation des deux flèches à l'écran. Seuils dans `apps/dive/legibility.ts`.
 *
 *    node tools/loop.mjs
 */

import { existsSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
if (!existsSync(resolve(ROOT, "dist/dive.html"))) {
  console.error("dist/dive.html absent — lancer `pnpm build`.");
  process.exit(2);
}

async function run(browser, url) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Page.navigate", { url });
  let r = null;
  for (let i = 0; i < 600 && !r; i++) {
    await new Promise((x) => setTimeout(x, 200));
    const { result } = await s("Runtime.evaluate", { expression: "window.__geno ? JSON.stringify(window.__geno) : null", returnByValue: true });
    if (result.value) r = JSON.parse(result.value);
  }
  await browser.send("Target.closeTarget", { targetId });
  return r;
}

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
const rows = [];
let tp53 = null;
try {
  for (const [name, flag] of [["vue de boucle", ""], ["témoin : sans isolement", "&noisolate"]]) {
    const r = await run(browser, `${base}/dive.html?f=data/gm12878.g3d&target=actb&legibility${flag}`);
    rows.push({ name, ...r });
  }
  // Semaine 12 : le même raccord, mesuré sur la région de TP53 (chr17).
  tp53 = (await run(browser, `${base}/dive.html?f=data/gm12878.g3d&target=tp53&junction`))?.junction ?? null;
} finally {
  browser.close();
  server.close();
  await shutdown(proc, profile);
}
const L = rows[0].loop;
console.log(`\n  Boucle visée : monomères ${L?.i}–${L?.j} (${L?.kb} kb), paire tenue dans ${((rows[0].arcWeight ?? 0) * 100).toFixed(1)} % des conformations\n`);
console.log(`  ${"".padEnd(26)} ${"ancre +".padStart(8)} ${"ancre −".padStart(8)} ${"anneau".padStart(7)} ${"arc".padStart(6)} ${"corps".padStart(7)}   flèches vers l'intérieur   verdict`);
for (const r of rows) {
  console.log(
    `  ${r.name.padEnd(26)} ${String(r.anchorPx[0]).padStart(8)} ${String(r.anchorPx[1]).padStart(8)} ${String(r.ringPx).padStart(7)} ${String(r.arcPx).padStart(6)} ${(r.body * 100).toFixed(0).padStart(6)} %   ${r.inward.map((x, k) => `${x ? "oui" : "non"} (${r.inwardDot[k].toFixed(2)})`).join(" / ").padEnd(24)} ${r.ok ? "\x1b[32mlisible\x1b[0m" : "\x1b[31mnon\x1b[0m"}`,
  );
}
for (const J of [rows[0].junction, tp53]) {
  if (!J) continue;
  console.log(`\n  Raccord S8 sur ${J.region}, au milieu du fondu (${J.widthUm.toFixed(2)} µm de champ) : ${(J.insideShare * 100).toFixed(0)} % de la fibre`);
  console.log(`  dans la silhouette des ${J.beads} billes de 750 kb qu'elle détaille ; centres écartés de ${(J.centreOffset * 100).toFixed(0)} % de la hauteur de l'image.`);
}
console.log(`\n  Flèches qui suivent leur brin vers l'intérieur, sur les ${rows[0].arcAnchors} ancres des ${rows[0].arcAnchors / 2} arcs dessinés : ${rows[0].arcAnchorsInward}.`);
console.log("  Seuils fixés d'avance : 50 px par ancre, anneau et arc ; 80 % du corps ; deux flèches vers l'intérieur.\n");
process.stdout.write(`JSON ${JSON.stringify({ rows, tp53 })}\n`);
process.exit(rows[0].ok ? 0 : 1);
