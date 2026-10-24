/** « Une sélection faite dans l'une des trois vues apparaît dans les deux autres, sans rebond. »
 *
 *  La page des vues liées (`vues.html`), de vrais événements souris envoyés par CDP — le
 *  chemin mesuré est celui d'un utilisateur. Pour chaque vue, N sélections à des endroits et
 *  des échelles variés :
 *
 *    3D       un pixel qui porte un objet, au noyau entier puis à l'échelle de la fibre ;
 *    igv      un élément d'une piste (gène, pic CTCF), à plusieurs loci ;
 *    HiGlass  un pixel de la carte, du génome entier à la région fine.
 *
 *  Après chaque geste, on attend que les trois vues soient au repos, puis 400 ms de plus
 *  (un rebond tardif doit avoir le temps de se produire), et on juge :
 *
 *    - **une** publication acceptée, venue de la vue du geste, et **aucune** refusée (un
 *      écho ou une égalité refusés seraient un rebond rattrapé — on veut zéro tentative) ;
 *    - chacune des deux autres vues a appliqué la sélection **une** fois ;
 *    - elle s'y **voit** : 3D — des pixels de billes surlignées dans le tampon d'identifiants,
 *      de la couleur de la surbrillance ; igv — une ROI par intervalle, dans la vue ;
 *      HiGlass — les bandes de chaque intervalle dessinées sur la carte, dans la vue ;
 *    - la latence, du geste à l'application dans chaque vue.
 *
 *  Puis la navigation : glisser dans igv, molette dans HiGlass, tourner la 3D — rien ne doit
 *  être publié. Puis le **témoin** : la même page sans store (`?naive`), où les vues se
 *  synchronisent deux à deux sur leur déplacement ; un seul geste doit y déclencher une
 *  cascade — c'est ce que la mesure doit savoir voir.
 *
 *    node tools/linked.mjs [N=30]
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";
import { writeFileSync } from "node:fs";

const N = Number(process.argv[2] ?? 30);
const W = 1440, H = 860;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function open(browser, url) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  await s("Page.navigate", { url });
  const ev = async (expr) => {
    const r = await s("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(`${expr.slice(0, 80)} : ${r.exceptionDetails.exception?.description ?? r.exceptionDetails.text}`);
    return r.result.value;
  };
  const t0 = Date.now();
  for (;;) {
    await sleep(200);
    const g = await ev("window.__geno ? JSON.stringify(window.__geno) : null");
    if (g) {
      if (!JSON.parse(g).ok) throw new Error(g);
      break;
    }
    if (Date.now() - t0 > 120_000) throw new Error("page jamais prête");
  }
  const move = (x, y) => s("Input.dispatchMouseEvent", { type: "mouseMoved", x, y });
  const click = async (x, y) => {
    await move(x, y);
    for (const type of ["mousePressed", "mouseReleased"]) {
      await s("Input.dispatchMouseEvent", { type, x, y, button: "left", buttons: type === "mousePressed" ? 1 : 0, clickCount: 1 });
    }
  };
  const drag = async (x, y, dx, dy) => {
    await move(x, y);
    await s("Input.dispatchMouseEvent", { type: "mousePressed", x, y, button: "left", buttons: 1, clickCount: 1 });
    for (let k = 1; k <= 8; k++) {
      await s("Input.dispatchMouseEvent", { type: "mouseMoved", x: x + (dx * k) / 8, y: y + (dy * k) / 8, button: "left", buttons: 1 });
      await sleep(16);
    }
    await s("Input.dispatchMouseEvent", { type: "mouseReleased", x: x + dx, y: y + dy, button: "left", buttons: 0, clickCount: 1 });
  };
  const wheel = (x, y, dy) => s("Input.dispatchMouseEvent", { type: "mouseWheel", x, y, deltaX: 0, deltaY: dy });
  const idle = async (quiet = 400, max = 15_000) => {
    const t = Date.now();
    let since = null;
    while (Date.now() - t < max) {
      await sleep(50);
      if (await ev("window.__linked.idle()")) {
        since ??= Date.now();
        if (Date.now() - since >= quiet) return true;
      } else since = null;
    }
    return false;
  };
  return { s, ev, click, drag, wheel, move, idle, close: () => browser.send("Target.closeTarget", { targetId }) };
}

const pct = (xs, q) => {
  const a = xs.filter(Number.isFinite).sort((x, y) => x - y);
  return a.length ? a[Math.min(a.length - 1, Math.floor(q * a.length))] : NaN;
};
const keyOf = (r) => `${r.chrom}:${r.start}-${r.end}`;

/** Un geste et ce qu'il a produit. */
async function act(page, source, gesture) {
  const before = await page.ev(`(() => { const L = window.__linked; return { log: L.store.log.length, a: { "3d": L.v3d.shown.applied, igv: L.igv.shown.applied, higlass: L.hg.shown.applied } }; })()`);
  await gesture();
  // La publication, puis le repos des trois vues, puis 400 ms de calme.
  const t = Date.now();
  while (Date.now() - t < 5_000 && (await page.ev("window.__linked.store.log.length")) === before.log) await sleep(20);
  const settled = await page.idle(400);
  return page.ev(`(() => {
    const L = window.__linked;
    const log = L.store.log.slice(${before.log});
    const cur = L.store.current;
    const applied = { "3d": L.v3d.shown.applied - ${before.a["3d"]}, igv: L.igv.shown.applied - ${before.a.igv}, higlass: L.hg.shown.applied - ${before.a.higlass} };
    const keys = cur ? cur.regions.map((r) => r.chrom + ":" + r.start + "-" + r.end) : [];
    const vis = L.v3d.visible();
    const rois = L.igv.rois();
    const drawn = L.hg.drawn();
    const t0 = L.input.at;
    return {
      settled: ${settled},
      log, cur, applied, keys,
      shows: {
        "3d": { ok: vis.ids > 0 && vis.yellow >= 0.9 * vis.pixels, ...vis, seq: L.v3d.shown.seq, hl: L.v3d.shown.hl && { nucleus: L.v3d.shown.hl.nucleus.length, preview: L.v3d.shown.hl.preview.length, fiber: L.v3d.shown.hl.fiber, copies: L.v3d.shown.hl.copies } },
        igv: { ok: keys.length > 0 && keys.every((k) => rois.some((r) => r.visible && r.region === k.replace(":", "-"))), rois: rois.length, seq: L.igv.shown.seq, navigated: L.igv.shown.navigated },
        higlass: { ok: keys.length > 0 && keys.every((k) => drawn.some((d) => d.visible && d.region === k)), drawn: drawn.length, seq: L.hg.shown.seq, navigated: L.hg.shown.navigated },
      },
      // De la publication à l'application dans chaque vue (la 3D : l'image peinte qui la montre).
      latency: cur ? { "3d": L.v3d.shown.paintedAt - cur.at, igv: L.igv.shown.at - cur.at, higlass: L.hg.shown.at - cur.at } : {},
      publishMs: cur ? cur.at - t0 : null,
      errors: L.errors.length,
    };
  })()`).then((r) => ({ source, ...r }));
}

function judge(r) {
  const accepted = r.log.filter((e) => e.accepted);
  const refused = r.log.filter((e) => !e.accepted);
  const others = ["3d", "igv", "higlass"].filter((v) => v !== r.source);
  const once = others.every((v) => r.applied[v] === 1) && r.applied[r.source] === 0;
  const seen = others.every((v) => r.shows[v].ok && r.shows[v].seq === r.cur?.seq);
  return {
    published: accepted.length === 1 && accepted[0].source === r.source,
    noBounce: refused.length === 0 && accepted.length === 1 && once,
    seen,
    seenBy: Object.fromEntries(others.map((v) => [v, r.shows[v].ok])),
    ok: accepted.length === 1 && accepted[0].source === r.source && refused.length === 0 && once && seen && r.settled,
  };
}

async function selections(page) {
  const out = [];
  const L = (e) => page.ev(`window.__linked.${e}`);
  // — 3D : au noyau entier, puis à l'échelle de la fibre —
  const stops = await page.ev("window.__linked.v3d.scene.stops.length");
  for (const [stop, n] of [[0, Math.ceil(N / 2)], [2, Math.floor(N / 2)]]) {
    if (stop >= stops) continue;
    await page.ev(`(() => { const v = window.__linked.v3d; v.goto(v.scene.stops[${stop}]); return 1; })()`);
    await page.idle(300);
    for (let i = 0, seed = 11 + stop * 1000; i < n && seed < 11 + stop * 1000 + 40 * n; seed++) {
      // Tiré juste avant le clic (la sélection précédente a effacé des copies entières), et
      // ailleurs que la sélection courante : le même clic deux fois serait refusé comme égal,
      // ce qui n'est pas un rebond mais fausserait le compte.
      const [p] = await L(`v3d.hits(1, ${seed})`);
      if (!p) continue;
      const k = await L(`v3d.peek(${p[0]}, ${p[1]})`);
      const cur = await page.ev("(() => { const c = window.__linked.store.current; return c ? c.regions.map((r) => r.chrom + ':' + r.start + '-' + r.end).join(' ') : null; })()");
      if (!k || k === cur) continue;
      out.push(await act(page, "3d", () => page.click(p[0], p[1])));
      i++;
    }
  }
  // — igv : un gène ou un pic CTCF, à plusieurs loci (le déplacement n'est pas une sélection) —
  const loci = ["chr17:7,560,000-7,790,000", "chr7:5,400,000-5,650,000", "chr1:155,000,000-155,400,000",
                "chr11:5,200,000-5,320,000", "chrX:153,700,000-154,000,000", "chr2:60,400,000-60,800,000"];
  let k = 0;
  while (out.filter((r) => r.source === "igv").length < N && k < 4 * N) {
    const locus = loci[k % loci.length];
    const track = k % 3 === 2 ? "CTCF" : "Gènes";
    k++;
    await page.ev(`window.__linked.igv.goto(${JSON.stringify(locus)})`);
    await page.idle(300);
    const ts = await L(`igv.targets(${JSON.stringify(track)}, 6)`);
    const t = ts[(k * 7) % Math.max(1, ts.length)];
    if (!t) continue;
    out.push(await act(page, "igv", () => page.click(t.x, t.y)));
  }
  // — HiGlass : des pixels de la carte, du génome entier à la région fine —
  const windows = await page.ev(`(() => { const r = window.__linked.hg.reader; const a = (c, p) => r.absolute(c, p);
    return [[0, r.total], [a("chr7", 0), a("chr8", 0)], [a("chr17", 0), a("chr17", 83e6)], [a("chr17", 5.7e6), a("chr17", 9.6e6)],
            [a("chr17", 7.4e6), a("chr17", 7.9e6)], [a("chr7", 4.2e6), a("chr7", 7.8e6)]]; })()`);
  for (let i = 0; i < N; i++) {
    const [x0, x1] = windows[i % windows.length];
    await page.ev(`(() => { window.__linked.hg.goto(${x0}, ${x1}); return 1; })()`);
    await page.idle(300);
    const [[x, y]] = await L(`hg.targets(1, ${101 + i})`);
    out.push(await act(page, "higlass", () => page.click(x, y)));
  }
  return out;
}

async function navigation(page) {
  const before = await page.ev("window.__linked.store.log.length");
  const box = (sel) => page.ev(`(() => { const e = document.querySelector(${JSON.stringify(sel)}); const b = e.getBoundingClientRect(); return [b.left, b.top, b.width, b.height]; })()`);
  const [ix, iy, iw, ih] = await box("#pigv");
  const [hx, hy, hw, hh] = await box("#phg");
  const [tx, ty, tw, th] = await box("#p3d");
  for (let k = 0; k < 5; k++) {
    await page.drag(ix + iw * 0.5, iy + ih * 0.62, (k % 2 ? 1 : -1) * 120, 0);           // igv : glisser la piste
    await page.wheel(hx + hw * 0.5, hy + hh * 0.5, k % 2 ? 240 : -240);                 // HiGlass : molette
    await page.drag(tx + tw * 0.5, ty + th * 0.5, 90, 30);                               // 3D : tourner
    await page.idle(300);
  }
  return { gestures: 15, published: (await page.ev("window.__linked.store.log.length")) - before };
}

async function control(browser, port) {
  const page = await open(browser, `http://127.0.0.1:${port}/vues.html?naive`);
  await page.idle(600);
  const box = await page.ev(`(() => { const b = document.querySelector("#pigv").getBoundingClientRect(); return [b.left, b.top, b.width, b.height]; })()`);
  await page.drag(box[0] + box[2] * 0.5, box[1] + box[3] * 0.62, -120, 0);              // un seul geste
  // La cascade, relevée sur l'horloge de la page : messages cumulés après le geste.
  const t0 = await page.ev("window.__linked.input.at");
  const curve = [];
  for (const at of [500, 1000, 2000, 4000, 8000]) {
    while ((await page.ev("performance.now()")) - t0 < at) await sleep(50);
    curve.push([at, await page.ev("window.__linked.naive.hops.length")]);
  }
  const r = await page.ev(`(() => { const n = window.__linked.naive; return { hops: n.hops.length, cap: n.cap, loci: [...new Set(n.hops.map((h) => h.locus))].length,
    first: n.hops.slice(0, 6).map((h) => [Math.round(h.t), h.from, h.to, h.locus]) }; })()`);
  await page.close();
  return { ...r, curve, seconds: 8 };
}

const server = await serve(0);
const port = server.address().port;
const { proc, url, profile } = await launch();
let out;
try {
  const browser = await cdp(url);
  const page = await open(browser, `http://127.0.0.1:${port}/vues.html`);
  await page.idle(800);
  const runs = await selections(page);
  const nav = await navigation(page);
  const errors = await page.ev("window.__linked.errors");
  // Ce que la page a demandé (fetch, XHR, images…) hors du serveur local : on en attend zéro.
  const external = await page.ev("performance.getEntriesByType('resource').map((e) => e.name).filter((u) => !u.startsWith(location.origin))");
  const resources = await page.ev("performance.getEntriesByType('resource').length");
  await page.close();
  const ctl = await control(browser, port);

  const by = (src) => runs.filter((r) => r.source === src);
  const summary = {};
  console.log(`\nsélections (${runs.length}), jugées après repos + 400 ms\n`);
  console.log("source    n    publiée  sans rebond  vue ailleurs  geste→publ.  publ.→3D    publ.→igv   publ.→HiGlass  (p50 / p95, ms)");
  for (const src of ["3d", "igv", "higlass"]) {
    const rs = by(src);
    const j = rs.map(judge);
    const lat = (v) => rs.map((r) => r.latency[v]);
    const cell = (v) => (v === src ? "—".padEnd(10) : `${pct(lat(v), 0.5).toFixed(0)} / ${pct(lat(v), 0.95).toFixed(0)}`.padEnd(10));
    summary[src] = {
      n: rs.length, published: j.filter((x) => x.published).length, noBounce: j.filter((x) => x.noBounce).length,
      seen: j.filter((x) => x.seen).length, ok: j.filter((x) => x.ok).length,
      seenBy: Object.fromEntries(["3d", "igv", "higlass"].filter((v) => v !== src).map((v) => [v, j.filter((x) => x.seenBy[v]).length])),
      latency: Object.fromEntries(["3d", "igv", "higlass"].filter((v) => v !== src).map((v) => [v, { p50: pct(lat(v), 0.5), p95: pct(lat(v), 0.95), max: pct(lat(v), 1) }])),
      failures: rs.filter((_, i) => !j[i].ok).slice(0, 3).map((r) => ({ log: r.log, applied: r.applied, keys: r.keys, shows: r.shows, settled: r.settled })),
    };
    const S = summary[src];
    const pub = rs.map((r) => r.publishMs);
    S.publish = { p50: pct(pub, 0.5), p95: pct(pub, 0.95) };
    console.log(`${src.padEnd(8)} ${String(S.n).padStart(3)}   ${String(S.published).padStart(4)}     ${String(S.noBounce).padStart(5)}        ${String(S.seen).padStart(5)}       ` +
                `${`${S.publish.p50.toFixed(0)} / ${S.publish.p95.toFixed(0)}`.padEnd(11)}  ${cell("3d")}  ${cell("igv")}  ${cell("higlass")}`);
  }
  const total = runs.length;
  const okAll = runs.map(judge).filter((x) => x.ok).length;
  console.log(`\ncritère : ${okAll}/${total} sélections vues dans les deux autres vues, une publication, zéro rebond`);
  console.log(`navigation : ${nav.gestures} gestes (glisser igv, molette HiGlass, tourner la 3D) → ${nav.published} publication(s)`);
  console.log(`témoin (sans store, deux-à-deux sur le déplacement) : un geste → ${ctl.hops} messages en ${ctl.seconds} s` +
              `${ctl.hops >= ctl.cap ? " (plafond atteint)" : ""}, ${ctl.loci} loci distincts ; cumul ${ctl.curve.map(([t, n]) => `${t / 1000} s : ${n}`).join(", ")}`);
  console.log(`erreurs de la page : ${errors.length}${errors.length ? ` — ${errors.slice(0, 3).join(" | ").slice(0, 300)}` : ""}`);
  console.log(`requêtes : ${resources}, dont ${external.length} hors du serveur local${external.length ? ` — ${external.slice(0, 3).join(" ")}` : ""}`);
  out = { n: N, summary, ok: okAll, total, navigation: nav, control: ctl, errors, external, resources, runs };
} finally {
  await shutdown(proc, profile).catch(() => {});
  server.close();
}
writeFileSync(new URL("../dist/linked-result.json", import.meta.url), JSON.stringify(out, null, 1));
process.exit(out && out.ok === out.total && out.navigation.published === 0 && out.control.hops > 3 ? 0 : 1);
