/** Semaine 22 — budgets de performance déterministes, ceux que la CI fait respecter.
 *
 *  La plongée entière (240 images, `dist/dive.html?bench=240`) est jouée sous SwiftShader avec un
 *  WebGL2 instrumenté (`tools/glcount.js`). On garde ce qui ne dépend pas de la machine :
 *  - par image (maximum sur la plongée) : appels de dessin, instances, sommets, pixels des passes
 *    plein écran ;
 *  - les octets GPU alloués (pic) ;
 *  - les octets du `.g3d` lus avant le premier rendu (`first.html`, préfixe de 16 Kio, S9) ;
 *  - la taille de chaque page de `dist/`.
 *  Les temps en ms ne sont jamais jugés ici : ils viennent du banc matériel. Budgets et règle de
 *  cliquet : `docs/S22_BUDGETS.md`.
 *
 *    node tools/perf.mjs                  mesure, écrit dist/perf/perf.json
 *    node tools/perf.mjs --check          échoue si un budget de tools/perf-budgets.json est dépassé
 *    node tools/perf.mjs --write          inscrit le cliquet : valeur mesurée + 10 %
 *    options : --ao=full|half|off  --window=context|visible  --cull=off|same|on  --frames=N
 *
 *  Données : `pipeline/data/g3d/` si présent, sinon `fixtures/` (empreinte vérifiée).
 */

import { createHash } from "node:crypto";
import { existsSync, mkdirSync, readdirSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { dirname, join, relative, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const arg = (k, d) => process.argv.find((a) => a.startsWith(`--${k}=`))?.split("=")[1] ?? d;
const CHECK = process.argv.includes("--check"), WRITE = process.argv.includes("--write");
const FRAMES = Number(arg("frames", 240));
const AO = arg("ao", ""), WINDOW = arg("window", ""), CULL = arg("cull", "");
const BUDGETS = resolve(ROOT, "tools/perf-budgets.json");

// Les données : le pipeline, sinon la fixture, dont l'empreinte est vérifiée.
if (!existsSync(resolve(ROOT, "pipeline/data/g3d/gm12878.g3d"))) {
  const fx = resolve(ROOT, "fixtures/gm12878.g3d");
  const want = readFileSync(resolve(ROOT, "fixtures/gm12878.g3d.sha256"), "utf8").trim().split(/\s+/)[0];
  const got = createHash("sha256").update(readFileSync(fx)).digest("hex");
  if (got !== want) { console.error(`fixtures/gm12878.g3d : empreinte ${got} au lieu de ${want}`); process.exit(2); }
  process.env.GENO_DATA = resolve(ROOT, "fixtures");
}
const { serve } = await import("./serve.mjs");
const GLCOUNT = readFileSync(resolve(ROOT, "tools/glcount.js"), "utf8");

async function run(browser, url) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Page.enable");
  await s("Page.addScriptToEvaluateOnNewDocument", { source: GLCOUNT });
  await s("Page.navigate", { url });
  const deadline = Date.now() + 30 * 60_000;
  let out = null;
  while (!out && Date.now() < deadline) {
    await new Promise((r) => setTimeout(r, 250));
    const { result } = await s("Runtime.evaluate", {
      expression: "window.__geno ? JSON.stringify({ geno: window.__geno, gl: window.__glcount }) : null", returnByValue: true });
    if (result.value) out = JSON.parse(result.value);
  }
  await browser.send("Target.closeTarget", { targetId });
  if (!out) throw new Error(`pas de résultat : ${url}`);
  return out;
}

function pages() {
  const out = {};
  const walk = (d) => {
    for (const f of readdirSync(d)) {
      const p = join(d, f);
      if (statSync(p).isDirectory()) { if (!/test|perf|discover-test/.test(f)) walk(p); }
      else if (p.endsWith(".html")) out[relative(resolve(ROOT, "dist"), p)] = statSync(p).size;
    }
  };
  walk(resolve(ROOT, "dist"));
  return out;
}

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
let result;
try {
  const q = `bench=${FRAMES}&w=960&h=540${AO ? `&ao=${AO}` : ""}${WINDOW ? `&window=${WINDOW}` : ""}${CULL ? `&cull=${CULL}` : ""}`;
  const dive = await run(browser, `${base}/dive.html?${q}`);
  const first = await run(browser, `${base}/first.html?f=data/gm12878.g3d&prefix=16384`);
  const fr = dive.gl.frames.filter(Boolean);
  const max = (k) => Math.max(...fr.map((f) => f[k]));
  const at = (k) => fr.findIndex((f) => f[k] === max(k));
  result = {
    quality: { ao: AO || "full", window: WINDOW || "context", cull: CULL || "off" },
    frames_counted: fr.length,
    unknown_formats: dive.gl.unknownFormats,
    budgets: {
      draws_per_frame: max("draws"),
      instances_per_frame: max("instances"),
      vertices_per_frame: max("vertices"),
      fullscreen_pixels_per_frame: max("fsPixels"),
      gpu_bytes_peak: dive.gl.peak,
      first_render_g3d_bytes: first.geno.g3d_bytes,
    },
    worst_frame: { draws: at("draws"), instances: at("instances"), vertices: at("vertices") },
    per_frame: fr,
    pages: pages(),
  };
} finally {
  browser.close?.();
  await shutdown(proc, profile).catch(() => {});
  server.close();
}

mkdirSync(resolve(ROOT, "dist/perf"), { recursive: true });
writeFileSync(resolve(ROOT, "dist/perf/perf.json"), JSON.stringify(result, null, 1));
const B = result.budgets;
console.log(`\n  Budgets déterministes — semaine 22 (qualité ${JSON.stringify(result.quality)}, ${result.frames_counted} images comptées)`);
for (const [k, v] of Object.entries(B)) console.log(`  ${k.padEnd(30)} ${v.toLocaleString("fr")}`);
console.log(`  pages : ${Object.entries(result.pages).map(([k, v]) => `${k} ${(v / 1024).toFixed(0)} Kio`).join(" · ")}`);

if (WRITE) {
  const up = (v) => Math.ceil(v * 1.1);
  const budgets = { rule: "valeur mesurée après correction + 10 % (docs/S22_BUDGETS.md)", measured_with: result.quality,
    ...Object.fromEntries(Object.entries(B).map(([k, v]) => [k, up(v)])),
    pages: Object.fromEntries(Object.entries(result.pages).map(([k, v]) => [k, up(v)])) };
  writeFileSync(BUDGETS, JSON.stringify(budgets, null, 1) + "\n");
  console.log(`\n  cliquet inscrit : ${relative(ROOT, BUDGETS)}`);
}
if (CHECK) {
  const budgets = JSON.parse(readFileSync(BUDGETS, "utf8"));
  const over = [];
  for (const [k, v] of Object.entries(B)) if (v > budgets[k]) over.push(`${k} ${v.toLocaleString("fr")} > ${budgets[k].toLocaleString("fr")}`);
  for (const [k, v] of Object.entries(result.pages)) {
    if (!(k in budgets.pages)) over.push(`page ${k} sans budget (${v} o)`);
    else if (v > budgets.pages[k]) over.push(`page ${k} ${v.toLocaleString("fr")} o > ${budgets.pages[k].toLocaleString("fr")} o`);
  }
  if (over.length) { console.log(`\n  \x1b[31mbudget dépassé\x1b[0m\n    ${over.join("\n    ")}`); process.exitCode = 1; }
  else console.log("\n  \x1b[32mbudgets tenus\x1b[0m");
}
