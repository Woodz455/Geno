/** Semaine 22 — ce que chaque variante de qualité coûte à l'œil : l'écart d'image contre la
 *  qualité pleine, image par image, sur la plongée entière (`dive.html?compare`). C'est la mesure
 *  de la règle de choix de `docs/S22_BUDGETS.md`.
 *
 *    node tools/quality-diff.mjs [frames]
 */

import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const FRAMES = Number(process.argv[2] ?? 240);
const VARIANTS = [
  ["élimination hors champ", "cull=1"],
  ["fenêtre visible", "window=visible"],
  ["occlusion à demi-résolution", "ao=half"],
  ["occlusion coupée", "ao=off"],
  ["élimination + fenêtre visible", "cull=1&window=visible"],
  ["élimination + occlusion moitié", "cull=1&ao=half"],
];

async function run(browser, url) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Page.enable");
  await s("Page.navigate", { url });
  let out = null;
  const deadline = Date.now() + 40 * 60_000;
  while (!out && Date.now() < deadline) {
    await new Promise((r) => setTimeout(r, 500));
    const { result } = await s("Runtime.evaluate", { expression: "window.__geno ? JSON.stringify(window.__geno) : null", returnByValue: true });
    if (result.value) out = JSON.parse(result.value);
  }
  await browser.send("Target.closeTarget", { targetId });
  if (!out) throw new Error(`pas de résultat : ${url}`);
  return out;
}

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
const rows = [];
try {
  for (const [name, q] of VARIANTS) {
    const r = await run(browser, `${base}/dive.html?compare&frames=${FRAMES}&w=960&h=540&${q}`);
    const d = r.diffs;
    const mean = d.reduce((a, b) => a + b, 0) / d.length;
    rows.push({ name, query: q, mean, max: Math.max(...d), zero: d.filter((x) => x === 0).length, frames: d.length, diffs: d });
    console.log(`  ${name.padEnd(32)} écart moyen ${mean.toFixed(4)} · max ${Math.max(...d).toFixed(3)} · images identiques ${d.filter((x) => x === 0).length}/${d.length}`);
  }
} finally {
  browser.close?.();
  await shutdown(proc, profile).catch(() => {});
  server.close();
}
mkdirSync(resolve(ROOT, "dist/perf"), { recursive: true });
writeFileSync(resolve(ROOT, "dist/perf/quality-diff.json"), JSON.stringify(rows, null, 1));
