/** « Les trois formes de requête amènent la caméra au bon endroit, et l'URL rejoue la vue à
 *  l'identique. » (semaine 16)
 *
 *  La plongée (`dive.html?nav`), de vrais gestes CDP. Critères fixés avant de mesurer :
 *
 *    1  résolution   les 645 requêtes de `geno search dump` (200 gènes, 200 positions, 200 rsID,
 *                    20 rsID absents, 20 requêtes mal formées, 5 hors modèle ou hors bornes) :
 *                    la page rend la cible de Python, 645/645 ; les témoins ne volent pas et le
 *                    disent ;
 *    2  la caméra    30 requêtes **tapées** par forme, les trois de la feuille de route
 *                    comprises : à la fin du vol, la cible est dans le tampon d'identifiants
 *                    (≥ 1 pixel), son centre projeté dans le tiers central de l'image, les deux
 *                    homologues surlignés, la fiche ouverte sur la bille (ou le monomère) qui
 *                    porte la position ; 90/90 ;
 *    3  l'URL        30 vues (après une recherche, une orbite, un zoom, une surimpression, un
 *                    clic) relues dans un onglet neuf : état égal au bit près, fiche égale,
 *                    image égale au pixel près (WebGL et surimpression) ; 30/30. Témoin : le
 *                    lacet décalé de 10⁻³ rad change l'image. Signets : 5 enregistrés, page
 *                    rechargée, liste rendue, chacun rejoue sa vue.
 *
 *    node tools/search.mjs
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";
import { readFileSync, writeFileSync } from "node:fs";

const W = 1280, H = 760;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const pct = (xs, q) => { const a = xs.filter(Number.isFinite).sort((x, y) => x - y); return a[Math.min(a.length - 1, Math.floor(q * a.length))]; };
const expected = JSON.parse(readFileSync(new URL("../dist/search-test/expected.json", import.meta.url), "utf8")).queries;

const server = await serve(0);
const port = server.address().port;
const base = `http://127.0.0.1:${port}/dive.html?nav`;
const { proc, url, profile } = await launch();
let out = null;

async function open(browser, href) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  const ev = async (expr) => {
    const r = await s("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(`${expr.slice(0, 80)} : ${r.exceptionDetails.exception?.description ?? r.exceptionDetails.text}`);
    return r.result.value;
  };
  const ready = async () => {
    const t0 = Date.now();
    for (;;) {
      await sleep(200);
      const g = await ev("window.__geno && window.__nav ? JSON.stringify(window.__geno) : null").catch(() => null);
      if (g) {
        if (!JSON.parse(g).ok) throw new Error(g);
        break;
      }
      if (Date.now() - t0 > 180_000) throw new Error("page jamais prête");
    }
  };
  await s("Page.enable");
  await s("Page.navigate", { url: href });
  await ready();
  const idle = async (extra = 300) => {
    for (let k = 0; k < 2400 && !(await ev("window.__nav.idle()")); k++) await sleep(25);
    await sleep(extra);
    for (let k = 0; k < 2400 && !(await ev("window.__nav.idle()")); k++) await sleep(25);
  };
  const mouse = (type, x, y, buttons = 0) => s("Input.dispatchMouseEvent", { type, x, y, button: type === "mouseMoved" ? (buttons ? "left" : "none") : "left", buttons, clickCount: 1 });
  const click = async (x, y) => { await mouse("mouseMoved", x, y); await mouse("mousePressed", x, y, 1); await mouse("mouseReleased", x, y); };
  const drag = async (x, y, dx, dy) => {
    await mouse("mouseMoved", x, y);
    await mouse("mousePressed", x, y, 1);
    for (let k = 1; k <= 8; k++) { await mouse("mouseMoved", x + (dx * k) / 8, y + (dy * k) / 8, 1); await sleep(16); }
    await mouse("mouseReleased", x + dx, y + dy);
  };
  const wheel = (x, y, deltaY) => s("Input.dispatchMouseEvent", { type: "mouseWheel", x, y, deltaX: 0, deltaY });
  const box = (sel) => ev(`(() => { const r = document.querySelector(${JSON.stringify(sel)}).getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()`);
  /** Taper une requête comme un utilisateur : clic dans le champ, tout sélectionner, touches, Entrée. */
  const type = async (text) => {
    const [x, y] = await box("#search input");
    await click(x, y);
    await s("Input.dispatchKeyEvent", { type: "keyDown", key: "a", code: "KeyA", modifiers: 2, windowsVirtualKeyCode: 65 });
    await s("Input.dispatchKeyEvent", { type: "keyUp", key: "a", code: "KeyA", modifiers: 2, windowsVirtualKeyCode: 65 });
    await s("Input.dispatchKeyEvent", { type: "keyDown", key: "Backspace", code: "Backspace", windowsVirtualKeyCode: 8 });
    await s("Input.dispatchKeyEvent", { type: "keyUp", key: "Backspace", code: "Backspace", windowsVirtualKeyCode: 8 });
    for (const ch of text) {
      await s("Input.dispatchKeyEvent", { type: "keyDown", key: ch, text: ch, unmodifiedText: ch });
      await s("Input.dispatchKeyEvent", { type: "keyUp", key: ch });
    }
    const n0 = await ev("window.__nav.flights().length");
    const t0 = Date.now();
    await s("Input.dispatchKeyEvent", { type: "keyDown", key: "Enter", code: "Enter", windowsVirtualKeyCode: 13, text: "\r" });
    await s("Input.dispatchKeyEvent", { type: "keyUp", key: "Enter", code: "Enter", windowsVirtualKeyCode: 13 });
    return { n0, t0 };
  };
  return { s, ev, idle, click, drag, wheel, type, box, targetId };
}

try {
  const browser = await cdp(url);
  const A = await open(browser, base);

  // — 1. Résolution : la page contre Python —
  const key = (r) => JSON.stringify([r.kind, r.error, r.hits.map((h) => [h.label, h.chrom, h.start, h.end, h.in_model, h.sources ?? null])]);
  const groups = {};
  const firstBad = {};
  const tRes = [];
  for (const q of expected) {
    const t0 = Date.now();
    const got = await A.ev(`window.__nav.resolve(${JSON.stringify(q.query)})`);
    tRes.push(Date.now() - t0);
    const g = (groups[q.group] ??= { n: 0, same: 0 });
    g.n++;
    if (key(got) === key(q.expected)) g.same++;
    else firstBad[q.group] ??= { query: q.query, got, want: q.expected };
  }
  // Les témoins, lancés comme le formulaire les lance : aucun vol, un message.
  const controls = expected.filter((q) => ["rs-absent", "malformed", "other"].includes(q.group));
  let refused = 0;
  const said = new Set();
  for (const q of controls) {
    const n0 = await A.ev("window.__nav.flights().length");
    const r = await A.ev(`window.__nav.search(${JSON.stringify(q.query)})`);
    const n1 = await A.ev("window.__nav.flights().length");
    const msg = await A.ev("window.__nav.message()");
    const shouldFly = q.expected.hits.some((h) => h.in_model);
    if ((n1 === n0) === !shouldFly && (shouldFly || msg)) refused++;
    said.add(msg.split(" — ")[0]);
  }

  // — 2. La caméra au bon endroit : 30 requêtes tapées par forme —
  const pick = (group, first, n) => [first, ...expected.filter((q) => q.group === group && q.query !== first).slice(0, n - 1).map((q) => q.query)];
  const genes = pick("gene", "BRCA1", 30);
  const positions = pick("position", "chr17:43,044,295", 28).concat(["chr17:7,673,802", "chr17:9,000,000"]);
  const rss = pick("rs", "rs334", 30);
  const arrivals = [];
  const views = [];
  const viewAt = async (what) => {
    if (views.length >= 30) return;
    await A.idle(400);
    views.push({ what, url: await A.ev("window.__nav.url()"), snap: await A.ev("window.__nav.snapshot()"),
                 panel: await A.ev("window.__nav.panel()"), state: await A.ev("window.__nav.state()") });
  };
  let k = 0;
  for (const [form, list] of [["gene", genes], ["position", positions], ["rs", rss]]) {
    for (const text of list) {
      const { n0, t0 } = await A.type(text);
      for (let i = 0; i < 400 && (await A.ev("window.__nav.flights().length")) === n0; i++) await sleep(10);
      await A.idle(300);
      const ms = Date.now() - t0;
      const f = await A.ev("window.__nav.flight()");
      const a = await A.ev("JSON.stringify(window.__nav.arrived())").then(JSON.parse);
      const flew = (await A.ev("window.__nav.flights().length")) > n0 && f.query === text.trim();
      // Les deux homologues : a et b, ou Xa et Xi pour le X du modèle 46,XX.
      const both = a && a.copies.filter((c) => c.startsWith(`${f.hit.chrom}:`)).length === 2;
      const cardOk = a?.card && a.card.level === f.level && a.card.bead === f.bead && f.start <= f.pos && f.pos < f.end;
      const ok = flew && a && a.pixels >= 1 && a.central && both && cardOk;
      arrivals.push({ form, query: text, ok, ms, depth: a?.depth, pixels: a?.pixels, central: a?.central, copies: a?.copies,
                      card: a?.card, bead: f?.bead, pos: f?.pos, resolveMs: f?.resolveMs, flightMs: f?.arrivedMs });
      // Douze de ces vues pour le critère 3 : quatre par forme.
      if (k++ % 10 === 0 || text === "chr17:7,673,802") await viewAt(`recherche ${text}`);
    }
  }

  // — 3. D'autres vues : orbite, zoom, surimpression, clic —
  const free = [520, 470];
  for (let i = 0; i < 6; i++) {
    await A.drag(free[0], free[1], 60 + 15 * i, -25 + 9 * i);
    await viewAt(`orbite ${i + 1}`);
  }
  for (let i = 0; i < 4; i++) {
    await A.wheel(free[0], free[1], i % 2 ? 240 : -300);
    await viewAt(`molette ${i + 1}`);
  }
  for (const kind of ["clinvar", "gwas", "loeuf", ""]) {
    await A.type(["TP53", "rs334", "chr7:5,530,601", "HBB"][views.length % 4]);
    await A.idle(300);
    await A.ev(`(() => { const s = document.querySelector("#hud select"); s.value = ${JSON.stringify(kind)}; s.dispatchEvent(new Event("change")); })()`);
    await viewAt(`surimpression ${kind || "aucune"}`);
  }
  for (let i = 0; views.length < 30 && i < 20; i++) {
    const [[x, y] = []] = await A.ev(`window.__nav.hits(1, ${100 + i})`);
    if (x === undefined) continue;
    await A.click(x, y);
    await sleep(200);
    await viewAt(`clic ${i + 1}`);
  }

  // Chaque vue, relue dans un onglet neuf.
  const replays = [];
  for (const v of views) {
    const B = await open(browser, v.url);
    await B.idle(600);
    const snap = await B.ev("window.__nav.snapshot()");
    const panel = await B.ev("window.__nav.panel()");
    const state = await B.ev("window.__nav.state()");
    const control = replays.length < 5 ? await B.ev("window.__nav.perturbed(1e-3)") : null;
    replays.push({ what: v.what, state: state === v.state, gl: snap.gl === v.snap.gl, overlay: snap.overlay === v.snap.overlay,
                   panel: panel === v.panel, control });
    await browser.send("Target.closeTarget", { targetId: B.targetId });
  }

  // Signets : cinq, la page rechargée, chacun rejoué.
  await A.ev("localStorage.removeItem('geno.signets.v1')");
  const A2 = await open(browser, base);
  const saved = [];
  for (const text of ["BRCA1", "rs334", "chr17:43,044,295", "TP53", "chr17:7,673,802"]) {
    await A2.type(text);
    await A2.idle(400);
    const [bx, by] = await A2.ev(`(() => { const b = [...document.querySelectorAll("#marks button")].find((x) => x.textContent === "+ signet"); const r = b.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()`);
    await A2.click(bx, by);
    await sleep(100);
    saved.push({ text, state: await A2.ev("window.__nav.state()"), snap: await A2.ev("window.__nav.snapshot()") });
  }
  await A2.s("Page.reload");
  await sleep(500);
  for (let i = 0; i < 300 && !(await A2.ev("!!window.__nav").catch(() => false)); i++) await sleep(100);
  const listed = await A2.ev("window.__nav.bookmarks().length");
  const links = await A2.ev(`[...document.querySelectorAll("#marks a")].length`);
  const marksOk = [];
  for (let i = 0; i < saved.length; i++) {
    const [x, y] = await A2.ev(`(() => { const a = document.querySelectorAll("#marks a")[${i}]; const r = a.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()`);
    await A2.click(x, y);
    await sleep(300);
    await A2.idle(600);
    const snap = await A2.ev("window.__nav.snapshot()");
    marksOk.push({ text: saved[i].text, state: (await A2.ev("window.__nav.state()")) === saved[i].state, gl: snap.gl === saved[i].snap.gl });
  }

  // — Jugement —
  const total = Object.values(groups).reduce((a, g) => a + g.n, 0);
  const same = Object.values(groups).reduce((a, g) => a + g.same, 0);
  const okArr = arrivals.filter((a) => a.ok).length;
  const okRep = replays.filter((r) => r.state && r.gl && r.overlay && r.panel).length;
  const controlOk = replays.filter((r) => r.control).every((r) => r.control.differing > 0);
  const marksPass = listed === 5 && links === 5 && marksOk.every((m) => m.state && m.gl);
  console.log("\nRecherche et navigation — semaine 16\n");
  console.log(`  1  résolution   ${same}/${total} égales à Python  (${Object.entries(groups).map(([g, v]) => `${g} ${v.same}/${v.n}`).join(" · ")})`);
  console.log(`                  témoins : ${refused}/${controls.length} sans vol et dits · résolution p50 ${pct(tRes, 0.5)} ms, p95 ${pct(tRes, 0.95)} ms`);
  for (const [g, b] of Object.entries(firstBad)) console.log(`     d'abord en échec (${g}) : ${JSON.stringify(b).slice(0, 300)}`);
  const byForm = (f) => arrivals.filter((a) => a.form === f);
  console.log(`  2  la caméra    ${okArr}/${arrivals.length}  (${["gene", "position", "rs"].map((f) => `${f} ${byForm(f).filter((a) => a.ok).length}/${byForm(f).length}`).join(" · ")})`);
  console.log(`                  requête tapée → fin du vol : p50 ${pct(arrivals.map((a) => a.ms), 0.5)} ms, p95 ${pct(arrivals.map((a) => a.ms), 0.95)} ms · profondeurs ${JSON.stringify(arrivals.reduce((m, a) => ({ ...m, [a.depth]: (m[a.depth] ?? 0) + 1 }), {}))}`);
  for (const a of arrivals.filter((x) => !x.ok).slice(0, 5)) console.log(`     échec : ${JSON.stringify(a)}`);
  console.log(`  3  l'URL        ${okRep}/${replays.length} vues rejouées à l'identique (état, image WebGL, surimpression, fiche)`);
  for (const r of replays.filter((x) => !(x.state && x.gl && x.overlay && x.panel)).slice(0, 8)) console.log(`     écart : ${JSON.stringify(r)}`);
  console.log(`                  témoin, lacet + 10⁻³ rad : ${replays.filter((r) => r.control).map((r) => r.control.differing).join(", ")} pixels changés sur ${W * H}`);
  console.log(`                  signets : ${listed} rendus après rechargement (${links} liens), ${marksOk.filter((m) => m.state && m.gl).length}/5 rejoués à l'identique`);
  const pass = same === total && refused === controls.length && okArr === 90 && okRep === 30 && controlOk && marksPass;
  console.log(`\n  verdict : ${pass ? "\x1b[32mfait\x1b[0m" : "\x1b[31méchec\x1b[0m"}`);
  out = { groups, firstBad, refused, controls: controls.length, said: [...said], arrivals, replays, views: views.map((v) => ({ what: v.what, url: v.url })),
          bookmarks: { listed, links, marks: marksOk }, resolveMs: tRes, pass };
} finally {
  await shutdown(proc, profile).catch(() => {});
  server.close();
}
writeFileSync(new URL("../dist/search-result.json", import.meta.url), JSON.stringify(out, null, 1));
process.exit(out?.pass ? 0 : 1);
