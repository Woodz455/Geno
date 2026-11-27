/** « La suppression d'une frontière rejoue en direct le cas EPHA4 des malformations des membres. »
 *
 *  La page `variant.html`, de vrais gestes envoyés par CDP :
 *
 *    préréglages   un clic sur chaque allèle simulé en 3D (DEL-B central, DelBS, frontière
 *                  épargnée) : le 3D change, les 4C se recalculent ;
 *    tracés        des délétions **glissées sur la piste**, tirées dans la famille DEL-B que le
 *                  texte autorise, et leurs témoins DelBS tracés de même ; la page calcule à
 *                  chaque fois l'allèle et, pour une délétion qui emporte la frontière, le même
 *                  allèle aux CTCF de la frontière épargnés.
 *
 *  Jugé, pour chaque geste : le temps du geste à la lecture affichée (critère : < 1 s au p95),
 *  et ce que la page affiche — enhancers → PAX3 contre le témoin (≥ 1,5 ×), enhancers → IHH contre
 *  le type sauvage (dans [0,67 ; 1,5]).
 *
 *    node tools/variant.mjs [N=12]
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";
import { writeFileSync } from "node:fs";

const N = Number(process.argv[2] ?? 12);
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const pct = (xs, q) => { const a = xs.filter(Number.isFinite).sort((x, y) => x - y); return a[Math.min(a.length - 1, Math.floor(q * a.length))]; };

const server = await serve(0);
const port = server.address().port;
const { proc, url, profile } = await launch();
let out;
try {
  const b = await cdp(url);
  const { targetId } = await b.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await b.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => b.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: 1440, height: 860, deviceScaleFactor: 1, mobile: false });
  await s("Page.navigate", { url: `http://127.0.0.1:${port}/variant.html` });
  const ev = async (e) => {
    const r = await s("Runtime.evaluate", { expression: e, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description ?? r.exceptionDetails.text);
    return r.result.value;
  };
  for (let k = 0; ; k++) {
    await sleep(200);
    const g = await ev("window.__geno ? JSON.stringify(window.__geno) : null");
    if (g) { if (!JSON.parse(g).ok) throw new Error(g); break; }
    if (k > 600) throw new Error("page jamais prête");
  }
  const move = (x, y, buttons = 0) => s("Input.dispatchMouseEvent", { type: "mouseMoved", x, y, button: buttons ? "left" : "none", buttons });
  const press = (x, y, type) => s("Input.dispatchMouseEvent", { type, x, y, button: "left", buttons: type === "mousePressed" ? 1 : 0, clickCount: 1 });
  const click = async (x, y) => { await move(x, y); await press(x, y, "mousePressed"); await press(x, y, "mouseReleased"); };
  const idle = async () => { for (let k = 0; k < 200 && !(await ev("window.__variant.idle()")); k++) await sleep(25); };
  const last = () => ev("window.__variant.timings.at(-1)");
  const count = () => ev("window.__variant.timings.length");

  const wt = (await ev("window.__variant.timings[0]")).contact;
  const results = [];

  // — Préréglages : un clic sur le bouton —
  for (const key of ["delb", "delbs", "delb-keep", "wt"]) {
    const box = await ev(`(() => { const b = [...document.querySelectorAll("header button")].find((x) => x.onclick && x.textContent && window.__variant.presets.includes(${JSON.stringify(key)}) && x === [...document.querySelectorAll("header button")][window.__variant.presets.indexOf(${JSON.stringify(key)})]); const r = b.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()`);
    const n0 = await count();
    await click(box[0], box[1]);
    while ((await count()) === n0) await sleep(10);
    await idle();
    results.push({ gesture: "préréglage", ...(await last()) });
  }

  // — Délétions tracées sur la piste —
  const fam = await ev("window.__variant.family()");
  const dels = fam.filter((v) => v.kind === "DEL-B");
  const delbs = fam.filter((v) => v.kind === "DelBS");
  const pick = [];
  for (let i = 0; i < N; i++) pick.push(dels[Math.floor((i * dels.length) / N)]);
  const track = await ev("window.__variant.track()");
  const y = track.top + track.height * 0.6;
  const drag = async (a, z) => {
    const xa = await ev(`window.__variant.xOf(${a})`), xz = await ev(`window.__variant.xOf(${z})`);
    const n0 = await count();
    await move(xa, y);
    await press(xa, y, "mousePressed");
    for (let k = 1; k <= 6; k++) { await move(xa + ((xz - xa) * k) / 6, y, 1); await sleep(8); }
    await press(xz, y, "mouseReleased");
    while ((await count()) === n0) await sleep(10);
    await idle();
    return last();
  };
  for (const v of pick) {
    const r = await drag(v.start, v.end);
    results.push({ gesture: "tracé", family: v.name, kind: "DEL-B", ...r });
    const c = delbs.find((x) => x.start === v.start);
    if (c) {
      const r2 = await drag(c.start, c.end);
      results.push({ gesture: "tracé", family: c.name, kind: "DelBS", ...r2, pairOf: v.name, delbPax3: r.contact.PAX3 });
    }
  }
  const errors = await ev("window.__variant.timings.length");

  // — Jugement —
  const tracedB = results.filter((r) => r.gesture === "tracé" && r.kind === "DEL-B");
  const tracedS = results.filter((r) => r.gesture === "tracé" && r.kind === "DelBS");
  const vsKeep = tracedB.map((r) => (r.control ? r.contact.PAX3 / r.control.PAX3 : NaN));
  const vsDelbs = tracedS.map((r) => r.delbPax3 / r.contact.PAX3);
  const ihh = tracedB.map((r) => r.contact.IHH / wt.IHH);
  const ms = results.map((r) => r.totalMs);
  console.log(`\nVariant structural au locus EPHA4 — ${results.length} gestes (${tracedB.length} DEL-B et ${tracedS.length} DelBS tracés, 4 préréglages)\n`);
  console.log(`  geste → lecture affichée     p50 ${pct(ms, 0.5).toFixed(0)} ms · p95 ${pct(ms, 0.95).toFixed(0)} ms · max ${Math.max(...ms).toFixed(0)} ms`);
  console.log(`  DEL-B / frontière épargnée   min ${Math.min(...vsKeep).toFixed(2)} · médiane ${pct(vsKeep, 0.5).toFixed(2)}  (enhancers → PAX3)`);
  console.log(`  DEL-B / DelBS                min ${Math.min(...vsDelbs).toFixed(2)} · médiane ${pct(vsDelbs, 0.5).toFixed(2)}`);
  console.log(`  DEL-B / sauvage, IHH         ${Math.min(...ihh).toFixed(2)} – ${Math.max(...ihh).toFixed(2)}`);
  for (const r of results.filter((x) => x.gesture === "préréglage")) {
    console.log(`  préréglage ${r.key.padEnd(10)} ${r.totalMs.toFixed(0)} ms · enhancers → PAX3 ×${(r.contact.PAX3 / wt.PAX3).toFixed(2)} du sauvage`);
  }
  const pass = pct(ms, 0.95) < 1000 && Math.min(...vsKeep) >= 1.5 && Math.min(...vsDelbs) >= 1.5 && ihh.every((x) => x >= 0.67 && x <= 1.5);
  console.log(`\n  verdict : ${pass ? "\x1b[32mrejoué en direct\x1b[0m" : "\x1b[31méchec\x1b[0m"}`);
  out = { n: N, results, wt, pass, timings: errors };
} finally {
  await shutdown(proc, profile).catch(() => {});
  server.close();
}
writeFileSync(new URL("../dist/variant-result.json", import.meta.url), JSON.stringify(out, null, 1));
process.exit(out?.pass ? 0 : 1);
