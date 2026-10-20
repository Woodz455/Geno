/** « Un clic sur n'importe quelle bille produit la fiche complète en moins de 100 ms. »
 *
 *  La page interactive (`dive.html?pick`), de vrais événements souris envoyés par CDP
 *  (`Input.dispatchMouseEvent`) : le chemin mesuré est celui d'un utilisateur — pointerup,
 *  lecture du pixel, mapping inverse, fiche, DOM. À chaque étape de la plongée (aperçu,
 *  noyau, fibre, nucléosomes, hélice, base), des pixels tirés au hasard parmi ceux qui
 *  portent un objet, à graine fixe. La page chronomètre chaque clic : pixel, fiche, DOM, et
 *  l'image suivante (requestAnimationFrame) — c'est ce temps-là qui est jugé.
 *
 *  Trois passes :
 *    chaud        les fiches préchargées (réseau non bridé), N clics : le critère ;
 *    Slow 4G      durée du préchargement, et un clic à froid juste après l'ouverture ;
 *    témoin       sans préchargement, sous Slow 4G : chaque clic attend son bloc — ce que
 *                 le préchargement achète.
 *
 *  SwiftShader est une borne haute : un vrai GPU lit un pixel plus vite.
 *
 *    node tools/pick.mjs [N=500]
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const N = Number(process.argv[2] ?? 500);
// Passes à jouer : warm, slow, control (défaut : les trois).
const PASSES = (process.argv[3] ?? "warm,slow,control").split(",");
const SLOW_4G = { latency: 562.5, downloadThroughput: (1.6e6 / 8) * 0.9, uploadThroughput: (750e3 / 8) * 0.9 };
const NONE = { latency: 0, downloadThroughput: -1, uploadThroughput: -1 };
const W = 1280, H = 720;
const STOP_NAMES = ["départ (aperçu, 80 µm)", "noyau (bille de 750 kb)", "fibre (2 kb)", "nucléosomes", "hélice", "base"];
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function open(browser, url, net) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Network.enable");
  await s("Network.setCacheDisabled", { cacheDisabled: true });
  await s("Network.emulateNetworkConditions", { offline: false, ...net });
  await s("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  const t0 = Date.now();
  await s("Page.navigate", { url });
  const ev = async (expr) => (await s("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true })).result.value;
  for (;;) {
    await sleep(100);
    if (await ev("!!window.__pick")) break;
    const g = await ev("window.__geno ? JSON.stringify(window.__geno) : null");
    if (g && !JSON.parse(g).ok) throw new Error(g);
    if (Date.now() - t0 > 10 * 60_000) throw new Error("page jamais prête");
  }
  const click = async (x, y) => {
    for (const type of ["mousePressed", "mouseReleased"]) {
      await s("Input.dispatchMouseEvent", { type, x, y, button: "left", buttons: type === "mousePressed" ? 1 : 0, clickCount: 1 });
    }
  };
  return { s, ev, click, close: () => browser.send("Target.closeTarget", { targetId }), openedMs: Date.now() - t0 };
}

async function clicks(page, perStop, seed, waitPrefetch) {
  if (waitPrefetch) {
    for (let k = 0; k < 6000 && !(await page.ev("window.__pick.prefetch.all > 0")); k++) await sleep(100);
  }
  const stops = await page.ev("window.__pick.stops");
  const rows = [];
  for (let i = 0; i < stops; i++) {
    await page.ev(`window.__pick.goto(${i}), true`);
    // Un clic sur une vue immobile : on attend que la page n'ait plus d'image à dessiner
    // (nucléosomes générés, séquence lue). Le cas « image en vol » est mesuré à part.
    for (let k = 0; k < 600 && !(await page.ev("window.__pick.idle()")); k++) await sleep(100);
    await sleep(100);
    const pts = await page.ev(`JSON.stringify(window.__pick.hits(${perStop}, ${seed + i}))`).then(JSON.parse);
    const before = await page.ev("window.__pick.log.length");
    for (const [x, y] of pts) {
      const n = await page.ev("window.__pick.log.length");
      await page.click(x, y);
      // Un clic à la fois : le suivant part quand l'image du précédent est peinte.
      for (let k = 0; k < 200 && !(await page.ev(`(window.__pick.log[${n}] || {}).frameMs !== undefined || window.__pick.log[${n}]?.hit === false`)); k++) await sleep(5);
    }
    // Laisse les images suivantes se peindre, et les fiches à froid arriver.
    for (let k = 0; k < 600; k++) {
      const pending = await page.ev(`window.__pick.log.slice(${before}).filter(e => e.hit && (e.frameMs === undefined || (!e.ready && e.waitMs === undefined))).length`);
      if (!pending) break;
      await sleep(100);
    }
    const log = await page.ev(`JSON.stringify(window.__pick.log.slice(${before}))`).then(JSON.parse);
    rows.push({ stop: STOP_NAMES[i] ?? `étape ${i}`, clicks: pts.length, log });
  }
  return rows;
}

const q = (a, p) => {
  if (!a.length) return NaN;
  const s = a.slice().sort((x, y) => x - y);
  return s[Math.min(s.length - 1, Math.floor(p * s.length))];
};
const f1 = (x) => (Number.isFinite(x) ? x.toFixed(1) : "—");

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
const out = {};
try {
  // 1. Chaud.
  if (PASSES.includes("warm")) {
    const page = await open(browser, `${base}/dive.html?f=data/gm12878.g3d&pick`, NONE);
    const rows = await clicks(page, Math.ceil(N / 6), 13, true);
    out.warm = { rows, prefetch: await page.ev("JSON.stringify(window.__pick.prefetch)").then(JSON.parse),
                 reader: await page.ev("JSON.stringify(window.__pick.reader())").then(JSON.parse) };
    await page.close();
  }
  // 2. Slow 4G : préchargement, et un clic à froid tout de suite.
  if (PASSES.includes("slow")) {
    const page = await open(browser, `${base}/dive.html?f=data/gm12878.g3d&pick`, SLOW_4G);
    const cold = await clicks(page, 3, 7, false);
    for (let k = 0; k < 6000 && !(await page.ev("window.__pick.prefetch.all > 0")); k++) await sleep(250);
    const warm = await clicks(page, 10, 99, false);
    out.slow = { opened_ms: page.openedMs, cold, warm, prefetch: await page.ev("JSON.stringify(window.__pick.prefetch)").then(JSON.parse),
                 reader: await page.ev("JSON.stringify(window.__pick.reader())").then(JSON.parse) };
    await page.close();
  }
  // 3. Témoin : sans préchargement, Slow 4G.
  if (PASSES.includes("control")) {
    const page = await open(browser, `${base}/dive.html?f=data/gm12878.g3d&pick&nofichepre`, SLOW_4G);
    await sleep(3000);       // le temps que l'en-tête des fiches arrive
    const rows = await clicks(page, 5, 21, false);
    out.control = { rows, reader: await page.ev("JSON.stringify(window.__pick.reader())").then(JSON.parse) };
    await page.close();
  }
} finally {
  browser.close();
  server.close();
  await shutdown(proc, profile);
}

if (!out.warm || !out.slow || !out.control) {
  process.stdout.write(`JSON ${JSON.stringify(out)}\n`);
  const cold = (out.slow?.cold ?? []).flatMap((r) => r.log.filter((e) => e.hit));
  console.log(`  passes partielles (${PASSES.join(", ")}) ; Slow 4G, clics avant le préchargement : ${cold.filter((e) => !e.ready).length}/${cold.length} attendent, attente ${cold.filter((e) => !e.ready).map((e) => f1(e.waitMs)).join(", ")} ms`);
  if (out.slow) console.log(`  préchargement Slow 4G : ${JSON.stringify(out.slow.prefetch)}`);
  process.exit(0);
}
const all = out.warm.rows.flatMap((r) => r.log.filter((e) => e.hit));
console.log(`\n  Clic → fiche, page interactive, événements souris CDP, SwiftShader (borne haute), ${W} × ${H}\n`);
console.log(`  ${"étape".padEnd(26)} ${"clics".padStart(5)} ${"objet".padEnd(26)} ${"pixel p50".padStart(9)} ${"fiche".padStart(6)} ${"DOM".padStart(6)} ${"peinte p50".padStart(10)} ${"p95".padStart(6)} ${"max".padStart(6)} ${"complètes".padStart(10)}`);
for (const r of out.warm.rows) {
  const L = r.log.filter((e) => e.hit);
  const kinds = [...new Set(L.map((e) => e.what))].join(", ");
  const fr = L.map((e) => e.frameMs);
  const full = L.filter((e) => e.ready && e.missing.length === 0).length;
  console.log(`  ${r.stop.padEnd(26)} ${String(L.length).padStart(5)} ${kinds.slice(0, 26).padEnd(26)} ${f1(q(L.map((e) => e.pickMs), 0.5)).padStart(9)} ${f1(q(L.map((e) => e.cardMs), 0.5)).padStart(6)} ${f1(q(L.map((e) => e.domMs), 0.5)).padStart(6)} ${f1(q(fr, 0.5)).padStart(10)} ${f1(q(fr, 0.95)).padStart(6)} ${f1(Math.max(...fr)).padStart(6)} ${`${full}/${L.length}`.padStart(10)}`);
}
const fr = all.map((e) => e.frameMs);
const complete = all.filter((e) => e.ready && e.missing.length === 0).length;
console.log(`  ${"tout".padEnd(26)} ${String(all.length).padStart(5)} ${"".padEnd(26)} ${f1(q(all.map((e) => e.pickMs), 0.5)).padStart(9)} ${f1(q(all.map((e) => e.cardMs), 0.5)).padStart(6)} ${f1(q(all.map((e) => e.domMs), 0.5)).padStart(6)} ${f1(q(fr, 0.5)).padStart(10)} ${f1(q(fr, 0.95)).padStart(6)} ${f1(Math.max(...fr)).padStart(6)} ${`${complete}/${all.length}`.padStart(10)}`);
const p = out.warm.prefetch;
console.log(`\n  Préchargement, réseau non bridé : tables ${f1(p.tables)} ms · noyau ${f1(p.noyau)} ms · tout ${f1(p.all)} ms · ${(out.warm.reader.bytes / 1e6).toFixed(2)} Mo en ${out.warm.reader.requests} requêtes`);
const sp = out.slow.prefetch;
console.log(`  Préchargement, Slow 4G : tables ${f1(sp.tables / 1000)} s · noyau ${f1(sp.noyau / 1000)} s · tout ${f1(sp.all / 1000)} s`);
const cold = out.slow.cold.flatMap((r) => r.log.filter((e) => e.hit));
const sw = out.slow.warm.flatMap((r) => r.log.filter((e) => e.hit));
console.log(`  Slow 4G, clics avant le préchargement : ${cold.filter((e) => !e.ready).length}/${cold.length} attendent leur bloc, attente p50 ${f1(q(cold.filter((e) => !e.ready).map((e) => e.waitMs), 0.5))} ms (coordonnées affichées tout de suite, ${f1(q(cold.map((e) => e.frameMs), 0.5))} ms)`);
console.log(`  Slow 4G, clics après : ${sw.filter((e) => e.ready).length}/${sw.length} prêts, peinte p95 ${f1(q(sw.map((e) => e.frameMs), 0.95))} ms`);
const ctl = out.control.rows.flatMap((r) => r.log.filter((e) => e.hit));
console.log(`  Témoin sans préchargement (Slow 4G) : ${ctl.filter((e) => !e.ready).length}/${ctl.length} attendent, fiche complète p50 ${f1(q(ctl.map((e) => e.waitMs ?? e.frameMs), 0.5))} ms, p95 ${f1(q(ctl.map((e) => e.waitMs ?? e.frameMs), 0.95))} ms`);
const ok = q(fr, 0.95) < 100 && complete === all.length && ctl.some((e) => !e.ready);
console.log(`\n  Verdict : ${ok ? "\x1b[32mfiche complète en moins de 100 ms (p95)\x1b[0m" : "\x1b[31mnon tenu\x1b[0m"} — et le témoin attend.\n`);
process.stdout.write(`JSON ${JSON.stringify(out)}\n`);
process.exit(ok ? 0 : 1);
