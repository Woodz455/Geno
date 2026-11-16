/** Semaine 20 — l'accessibilité du mode Découverte, mesurée dans Chromium, en français et en anglais.
 *
 *  Critères (fixés avant la première exécution) — chaque écran = introduction, 7 étapes, 4 haltes,
 *  fin, soit 13 par langue et 26 en tout :
 *    1  axe       axe-core 4.10.3 (règles wcag2a, wcag2aa, wcag21a, wcag21aa, wcag22aa) sur chaque
 *                 écran à 1 280 × 720 : 0 violation ;
 *    2  contraste recalculé ici sur chaque élément de texte visible du panneau (et la bulle du
 *                 glossaire) : ≥ 4,5:1, ≥ 3:1 pour le grand texte ; bordures des boutons, champs et
 *                 groupes, et contours de focus : ≥ 3:1 contre leur fond ;
 *    3  reflow    à 320 × 640, et à 1 280 × 720 zoomé à 200 % (640 × 360 px CSS, densité 2) : aucun
 *                 élément du panneau ne sort de la fenêtre, aucun défilement horizontal, sur les 26
 *                 écrans ;
 *    4  clavier   le parcours entier à la flèche → : après chaque changement d'écran, le focus est sur
 *                 le titre (h2) du nouvel écran ; Échap ferme la bulle d'un terme sans déplacer le
 *                 focus ; chaque élément atteint par Tab montre un contour visible ;
 *    5  vue 3D    à chaque arrêt (7 étapes, 4 haltes), la description annoncée nomme l'objet, le
 *                 chromosome et la catégorie que la scène montre réellement (le picking au centre,
 *                 `categories()`), dans la langue de la page, et le détail d'évidence est reconnu :
 *                 11/11 par langue.
 *
 *  La non-régression des critères S19 dans chaque langue est `tools/discover.mjs fr|en`.
 *
 *    node tools/a11y.mjs        (après make g3d fiche linked search, et node tools/build.mjs)
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";
import { existsSync, mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/discover-test");
mkdirSync(OUT, { recursive: true });
for (const f of ["gm12878.g3d", "gm12878.fiche", "gm12878.contacts", "gm12878.search"]) {
  if (!existsSync(resolve(ROOT, "pipeline/data/g3d", f))) {
    console.error(`manque pipeline/data/g3d/${f} : make g3d fiche linked search`);
    process.exit(2);
  }
}
const AXE = readFileSync(resolve(ROOT, "node_modules/axe-core/axe.min.js"), "utf8");
const AXE_TAGS = ["wcag2a", "wcag2aa", "wcag21a", "wcag21aa", "wcag22aa"];
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const server = await serve(0);
const port = server.address().port;
const { proc, url, profile } = await launch();
const browser = await cdp(url);

async function open(lang, width, height, scale = 1) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width, height, deviceScaleFactor: scale, mobile: false });
  const ev = async (expr) => {
    const r = await s("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(`${expr.slice(0, 80)} : ${r.exceptionDetails.exception?.description ?? r.exceptionDetails.text}`);
    return r.result.value;
  };
  await s("Page.enable");
  await s("Page.navigate", { url: `http://127.0.0.1:${port}/decouverte.html?lang=${lang}&motion=0` });
  const t0 = Date.now();
  while (!(await ev("!!(window.__discover && window.__discover.ready)").catch(() => false))) {
    if (Date.now() - t0 > 240_000) throw new Error("page jamais prête");
    await sleep(200);
  }
  const key = async (k, code = k, vk = 0, mods = 0) => {
    await s("Input.dispatchKeyEvent", { type: "keyDown", key: k, code, windowsVirtualKeyCode: vk, modifiers: mods, ...(k === "Enter" ? { text: "\r" } : {}) });
    await s("Input.dispatchKeyEvent", { type: "keyUp", key: k, code, windowsVirtualKeyCode: vk, modifiers: mods });
  };
  const settle = async (id) => {
    for (let k = 0; k < 2400 && !(await ev("window.__discover.arrived()")); k++) await sleep(50);
    if (id === "tad") for (let k = 0; k < 100 && !(await ev("window.__discover.contacts()")); k++) await sleep(100);
    if (["hbb", "fto", "hoxd", "xi"].includes(id)) for (let k = 0; k < 100 && !(await ev(`document.querySelector("#halt-facts")?.textContent ?? ""`)); k++) await sleep(100);
    await sleep(300);
  };
  return { s, ev, key, settle, close: () => browser.send("Target.closeTarget", { targetId }) };
}

/** Les écrans d'une langue, parcourus au clavier (→) puis par les boutons de halte ; `visit` est
 *  appelé sur chacun, une fois arrivé. */
async function tour(P, visit) {
  const stops = await P.ev("window.__discover.stops()");
  await P.settle("intro");
  await visit("intro");
  for (let k = 0; k < 8; k++) {
    await P.key("ArrowRight", "ArrowRight", 39);
    const id = k < 7 ? stops[k] : "fin";
    await P.settle(id);
    await visit(id);
  }
  for (const id of stops.slice(7)) {
    await P.ev(`document.querySelector('[data-halt="${id}"]').click()`);
    await P.settle(id);
    await visit(id);
  }
}

// Le contraste, calculé dans la page : couleurs calculées, fond composé jusqu'au corps.
const CONTRAST = `(() => {
  const parse = (c) => { const m = c.match(/rgba?\\(([^)]+)\\)/); if (!m) return null; const v = m[1].split(/[ ,/]+/).filter(Boolean).map(Number); return { r: v[0], g: v[1], b: v[2], a: v.length > 3 ? v[3] : 1 }; };
  const over = (top, under) => ({ r: top.r * top.a + under.r * (1 - top.a), g: top.g * top.a + under.g * (1 - top.a), b: top.b * top.a + under.b * (1 - top.a), a: 1 });
  const lum = (c) => { const f = (x) => { x /= 255; return x <= 0.03928 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4; }; return 0.2126 * f(c.r) + 0.7152 * f(c.g) + 0.0722 * f(c.b); };
  const ratio = (a, b) => { const [x, y] = [lum(a), lum(b)].sort((p, q) => q - p); return (x + 0.05) / (y + 0.05); };
  const bgOf = (el) => {
    const layers = [];
    for (let e = el; e; e = e.parentElement) { const c = parse(getComputedStyle(e).backgroundColor); if (c && c.a > 0) layers.push(c); }
    let acc = { r: 11, g: 15, b: 20, a: 1 };          // --bg du corps
    for (const c of layers.reverse()) acc = over(c, acc);
    return acc;
  };
  const visible = (el) => { const r = el.getBoundingClientRect(); const cs = getComputedStyle(el); return r.width > 1 && r.height > 1 && cs.visibility !== "hidden" && cs.display !== "none" && !el.closest("[hidden], .sr, [aria-hidden=true]:not(#tip)"); };
  const out = { texts: 0, low: [], components: 0, lowComp: [], min: 99 };
  const roots = [document.getElementById("story"), document.getElementById("tip")].filter((r) => r && getComputedStyle(r).display !== "none");
  for (const root of roots) for (const el of [root, ...root.querySelectorAll("*")]) {
    if (!visible(el) || el.closest("button[disabled]")) continue;
    const own = [...el.childNodes].some((n) => n.nodeType === 3 && n.textContent.trim());
    if (own) {
      const cs = getComputedStyle(el), fg = parse(cs.color), bg = bgOf(el);
      const size = parseFloat(cs.fontSize), bold = Number(cs.fontWeight) >= 700;
      const need = size >= 24 || (size >= 18.66 && bold) ? 3 : 4.5;
      const r = ratio(over(fg, bg), bg);
      out.texts++; out.min = Math.min(out.min, r);
      if (r < need) out.low.push(el.tagName.toLowerCase() + (el.id ? "#" + el.id : "") + " " + r.toFixed(2) + " « " + el.textContent.trim().slice(0, 30) + " »");
    }
    if (el.matches("button:not(.term):not([disabled]), textarea, fieldset")) {
      const cs = getComputedStyle(el), bc = parse(cs.borderTopColor), bg = bgOf(el.parentElement);
      if (bc && parseFloat(cs.borderTopWidth) > 0) {
        const r = ratio(over(bc, bg), bg);
        out.components++;
        if (r < 3) out.lowComp.push(el.tagName.toLowerCase() + (el.id ? "#" + el.id : "") + " bordure " + r.toFixed(2));
      }
    }
  }
  return out;
})()`;

// Le contour de l'élément qui a le focus, contre son fond.
const FOCUS = `(() => {
  const el = document.activeElement; if (!el || el === document.body) return null;
  const cs = getComputedStyle(el);
  const parse = (c) => { const m = c.match(/rgba?\\(([^)]+)\\)/); if (!m) return null; const v = m[1].split(/[ ,/]+/).filter(Boolean).map(Number); return { r: v[0], g: v[1], b: v[2] }; };
  const lum = (c) => { const f = (x) => { x /= 255; return x <= 0.03928 ? x / 12.92 : ((x + 0.055) / 1.055) ** 2.4; }; return 0.2126 * f(c.r) + 0.7152 * f(c.g) + 0.0722 * f(c.b); };
  const oc = parse(cs.outlineColor), w = parseFloat(cs.outlineWidth), bg = { r: 14, g: 18, b: 24 };
  const r = oc ? ((Math.max(lum(oc), lum(bg)) + 0.05) / (Math.min(lum(oc), lum(bg)) + 0.05)) : 0;
  return { tag: el.tagName.toLowerCase(), cls: el.className, text: (el.textContent || el.getAttribute("aria-label") || "").trim().slice(0, 24), inStory: !!el.closest("#story"), style: cs.outlineStyle, width: w, ratio: r, ok: cs.outlineStyle !== "none" && w >= 1 && r >= 3 };
})()`;

const REFLOW = `(() => {
  const W = document.documentElement.clientWidth, out = [];
  for (const el of document.querySelectorAll("#story, #story *")) {
    const r = el.getBoundingClientRect(), cs = getComputedStyle(el);
    if (r.width < 1 || r.height < 1 || cs.display === "none" || el.closest(".sr, [hidden]")) continue;
    if (r.right > W + 0.5 || r.left < -0.5) out.push(el.tagName.toLowerCase() + (el.id ? "#" + el.id : "") + " " + Math.round(r.left) + "–" + Math.round(r.right));
  }
  const story = document.getElementById("story");
  return { width: W, out: out.slice(0, 8), n: out.length, hscroll: document.documentElement.scrollWidth > W + 0.5 || story.scrollWidth > story.clientWidth + 0.5 };
})()`;

const result = {};
try {
  for (const lang of ["fr", "en"]) {
    const R = (result[lang] = { axe: [], contrast: [], reflow: [], focusTitle: [], focusRing: [], escape: null, view: [] });

    // 1, 2, 4, 5 — à 1 280 × 720.
    const P = await open(lang, 1280, 720);
    await P.ev(AXE);
    await tour(P, async (id) => {
      const ax = await P.ev(`axe.run(document, { runOnly: { type: "tag", values: ${JSON.stringify(AXE_TAGS)} }, resultTypes: ["violations"] }).then((r) => r.violations.map((v) => ({ id: v.id, impact: v.impact, n: v.nodes.length, where: v.nodes.slice(0, 3).map((n) => n.target.join(" ")) })))`);
      R.axe.push({ id, violations: ax });
      R.contrast.push({ id, ...(await P.ev(CONTRAST)) });
      if (id !== "intro") {
        const f = await P.ev(`(() => { const a = document.activeElement; return { tag: a?.tagName.toLowerCase(), text: a?.textContent.trim().slice(0, 40), h2: document.querySelector("#story h2")?.textContent.trim().slice(0, 40) }; })()`);
        R.focusTitle.push({ id, ok: f.tag === "h2" && f.text === f.h2, got: `${f.tag} « ${f.text} »` });
      }
      // Tab depuis le titre : chaque arrêt du focus a son contour.
      const seen = new Set();
      for (let k = 0; k < 90; k++) {
        await P.key("Tab", "Tab", 9);
        const st = await P.ev(FOCUS);
        if (!st || !st.inStory) break;
        const sig = `${st.tag}.${st.cls}.${st.text}`;
        if (seen.has(sig)) break;
        seen.add(sig);
        if (!st.ok) R.focusRing.push(`${id} : ${st.tag}.${st.cls} « ${st.text} » ${st.style} ${st.width}px ${st.ratio.toFixed(2)}`);
      }
      (R.focusCount ??= 0), (R.focusCount += seen.size);
      // Échap : sur l'écran TAD, la bulle d'un terme se ferme sans déplacer le focus.
      if (id === "tad") {
        await P.ev(`document.querySelector("#story .term").focus()`);
        await sleep(50);
        const before = await P.ev(`({ shown: document.getElementById("tip").style.display === "block", term: document.activeElement.dataset.term })`);
        R.contrast.push({ id: "tad (bulle)", ...(await P.ev(CONTRAST)) });
        await P.key("Escape", "Escape", 27);
        await sleep(50);
        const after = await P.ev(`({ shown: document.getElementById("tip").style.display === "block", term: document.activeElement.dataset.term, current: window.__discover.current() })`);
        R.escape = { ok: before.shown && !after.shown && after.term === before.term && after.current === "tad", before, after };
        await P.ev(`document.querySelector("#story h2").focus()`);
      }
      // La vue en mots, contre ce que le picking voit au centre.
      if (!["intro", "fin"].includes(id)) {
        const d = await P.ev("window.__discover.described()");
        const c = await P.ev("window.__discover.centre()");
        const cats = await P.ev("window.__discover.categories()");
        const known = await P.ev("window.__discover.evidenceKnown()");
        const n = c ? c.chrom.replace(/^chr/, "") : "";
        const ok = !!d && !!c && d.what === c.what && d.chrom === c.chrom && cats.includes(d.cat) && d.known && known
          && new RegExp(`(?:chromosome) ${n}(?![0-9])`).test(d.text);
        R.view.push({ id, ok, text: d?.text ?? "", centre: c && `${c.what} ${c.chrom}`, known });
      }
    });
    await P.close();

    // 3 — reflow.
    for (const [w, h, sc, name] of [[320, 640, 1, "320 px"], [640, 360, 2, "zoom 200 %"]]) {
      const Q = await open(lang, w, h, sc);
      await tour(Q, async (id) => { R.reflow.push({ id, at: name, ...(await Q.ev(REFLOW)) }); });
      await Q.close();
    }
  }
} finally {
  browser.close?.();
  await shutdown(proc, profile).catch(() => {});
  server.close();
}

writeFileSync(resolve(OUT, "a11y.json"), JSON.stringify(result, null, 1));
const crit = {};
console.log("\n  Mode Découverte — accessibilité (semaine 20)");
for (const lang of ["fr", "en"]) {
  const R = result[lang];
  const axeN = R.axe.reduce((a, x) => a + x.violations.length, 0);
  const lowT = R.contrast.flatMap((x) => x.low.map((l) => `${x.id} ${l}`)), lowC = R.contrast.flatMap((x) => x.lowComp.map((l) => `${x.id} ${l}`));
  const minC = Math.min(...R.contrast.map((x) => x.min));
  const flow = R.reflow.filter((x) => x.n || x.hscroll);
  const titles = R.focusTitle.filter((x) => x.ok).length;
  const views = R.view.filter((x) => x.ok).length;
  crit[lang] = {
    axe: axeN === 0 && R.axe.length === 13,
    contraste: lowT.length === 0 && lowC.length === 0,
    reflow: flow.length === 0 && R.reflow.length === 26,
    clavier: titles === R.focusTitle.length && R.focusTitle.length === 12 && !!R.escape?.ok && R.focusRing.length === 0,
    vue: views === 11 && R.view.length === 11,
  };
  const c = crit[lang], m = (b) => (b ? "✓" : "✗");
  console.log(`  [${lang}]`);
  console.log(`  1  axe         ${R.axe.length} écrans, ${axeN} violations${axeN ? " : " + R.axe.filter((x) => x.violations.length).map((x) => `${x.id} ${x.violations.map((v) => `${v.id}×${v.n} (${v.where.join(", ")})`).join(" ")}`).join(" · ") : ""} ${m(c.axe)}`);
  console.log(`  2  contraste   ${R.contrast.reduce((a, x) => a + x.texts, 0)} textes, minimum ${minC.toFixed(2)}:1 · ${R.contrast.reduce((a, x) => a + x.components, 0)} bordures${lowT.length || lowC.length ? " · sous le seuil : " + [...lowT, ...lowC].slice(0, 8).join(" ; ") : ""} ${m(c.contraste)}`);
  console.log(`  3  reflow      ${R.reflow.length} écrans (320 px et zoom 200 %), ${flow.length} débordent${flow.length ? " : " + flow.slice(0, 6).map((x) => `${x.at} ${x.id} ${x.out.join(", ")}${x.hscroll ? " (défilement)" : ""}`).join(" · ") : ""} ${m(c.reflow)}`);
  console.log(`  4  clavier     titre focalisé ${titles}/${R.focusTitle.length} · Échap ${R.escape?.ok ? "ferme la bulle, focus gardé" : "✗ " + JSON.stringify(R.escape)} · ${R.focusCount} arrêts de Tab, ${R.focusRing.length} sans contour${R.focusRing.length ? " : " + R.focusRing.slice(0, 5).join(" ; ") : ""} ${m(c.clavier)}`);
  console.log(`  5  vue 3D      ${views}/${R.view.length} descriptions justes${R.view.filter((x) => !x.ok).map((x) => ` · ✗ ${x.id} « ${x.text} » (centre : ${x.centre})`).join("")} ${m(c.vue)}`);
  console.log(`                 ex. « ${R.view.find((x) => x.id === "helice")?.text ?? ""} »`);
}
const pass = Object.values(crit).every((c) => Object.values(c).every(Boolean));
console.log(`\n  verdict : ${pass ? "\x1b[32mfait\x1b[0m" : "\x1b[31méchec\x1b[0m"}`);
process.exitCode = pass ? 0 : 1;
