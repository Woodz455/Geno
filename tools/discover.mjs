/** Semaine 19 — le mode Découverte, mesuré dans la page par de vrais gestes (CDP).
 *
 *  Critères, fixés avant cette mesure (docs/VALIDATION.md § S19) :
 *    1  parcours    au clavier seul — → jusqu'à la fin, Tab + Entrée pour chaque halte —, chaque
 *                   arrêt (7 étapes, 4 haltes) met sa cible au centre de l'image : 11/11.
 *                   noyau : une bille ; territoires : une bille du chromosome 17 ; compartiments :
 *                   la bille de TP53 ; TAD : un monomère de la fenêtre de la carte ; boucle : un
 *                   monomère à ±10 de l'ancre ; fibre : à ±20 kb de la base visée ; hélice : la base
 *                   visée, exactement ; halte : la bille qui porte le TSS du gène.
 *                   Et la fin du parcours ramène à l'introduction.
 *    2  glossaire   sur chaque écran, chaque terme balisé montre sa définition au survol **et** au
 *                   focus clavier ; 100 %.
 *    3  évidence    à chaque arrivée, la catégorie annoncée par le texte est l'une de celles de la
 *                   couche réellement rendue (`scene.evidence()`) : 11/11.
 *    4  faits       la carte de contacts de l'étape TAD est dessinée sur la fenêtre du domaine ;
 *                   la fiche de chaque halte nomme son gène : 4/4.
 *    5  URL         chaque arrêt rejoué dans un onglet neuf par `#t=` : image identique au pixel
 *                   près : 11/11.
 *
 *    node tools/discover.mjs           (après make g3d fiche linked search, et node tools/build.mjs)
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";
import { existsSync, mkdirSync, writeFileSync } from "node:fs";
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
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const W = 1280, H = 720, PANEL = 400;
const server = await serve(0);
const port = server.address().port;
const base = `http://127.0.0.1:${port}/decouverte.html?w=${W - PANEL}&h=${H}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);

async function open(href) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  const ev = async (expr) => {
    const r = await s("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(`${expr.slice(0, 80)} : ${r.exceptionDetails.exception?.description ?? r.exceptionDetails.text}`);
    return r.result.value;
  };
  await s("Page.enable");
  await s("Page.navigate", { url: href });
  const t0 = Date.now();
  while (!(await ev("!!(window.__discover && window.__discover.ready)").catch(() => false))) {
    if (Date.now() - t0 > 240_000) throw new Error(`page jamais prête : ${await ev("document.body.innerText.slice(0, 300)").catch(() => "?")}`);
    await sleep(200);
  }
  const arrived = async () => {
    const t1 = Date.now();
    for (let k = 0; ; k++) {
      if (await ev("window.__discover.arrived()")) break;
      if (Date.now() - t1 > 120_000) throw new Error("jamais arrivé");
      await sleep(50);
    }
    await sleep(400);
    while (!(await ev("window.__discover.arrived()"))) await sleep(50);
  };
  const key = async (k, code = k, vk = 0) => {
    await s("Input.dispatchKeyEvent", { type: "keyDown", key: k, code, windowsVirtualKeyCode: vk });
    await s("Input.dispatchKeyEvent", { type: "keyUp", key: k, code, windowsVirtualKeyCode: vk });
  };
  const mouse = (x, y) => s("Input.dispatchMouseEvent", { type: "mouseMoved", x, y, button: "none", buttons: 0 });
  const shot = async () => (await s("Page.captureScreenshot", { format: "png", clip: { x: 0, y: 0, width: W - PANEL, height: H, scale: 1 } })).data;
  return { s, ev, arrived, key, mouse, shot, close: () => browser.send("Target.closeTarget", { targetId }) };
}

/** La cible attendue au centre de chaque arrêt. */
function judgeCentre(id, c, ctx) {
  if (!c) return { ok: false, why: "rien au centre" };
  const inside = (pos) => c.start <= pos && pos < c.end;
  switch (id) {
    case "noyau": return { ok: c.level === "noyau" || c.what.includes("bille"), why: c.what };
    case "territoires": return { ok: c.chrom === "chr17", why: `${c.chrom} ${c.what}` };
    case "compartiments": return { ok: c.chrom === "chr17" && inside(ctx.target), why: `${c.chrom}:${c.start}-${c.end} ${c.what}` };
    case "tad": {
      const d = ctx.domain, len = d.end - d.start;
      return { ok: c.level !== "noyau" && c.chrom === d.chrom && c.start >= d.start - len / 2 && c.end <= d.end + len / 2, why: `${c.what} ${c.start}` };
    }
    case "boucle": {
      const bp = ctx.bp, a = ctx.region.start + ctx.focus.i * bp;
      return { ok: c.level !== "noyau" && Math.abs(c.start - a) <= 10 * bp, why: `${c.what} ${c.start} (ancre ${a})` };
    }
    case "fibre": return { ok: c.level !== "noyau" && Math.abs((c.pos >= 0 ? c.pos : c.start) - ctx.target) <= 20_000, why: `${c.what} ${c.pos}` };
    case "helice": return { ok: c.pos === ctx.target, why: `${c.what} ${c.pos} (visée ${ctx.target})` };
    default: {
      const g = ctx.genes[id];
      return { ok: !!g && c.chrom === g.chrom && inside(g.tss), why: `${c.chrom}:${c.start}-${c.end} pour ${g?.name} ${g?.tss}` };
    }
  }
}

try {
  const P = await open(base);
  const stops = await P.ev("window.__discover.stops()");
  const steps = stops.slice(0, 7), halts = stops.slice(7);
  const ctx = {
    target: 7_673_801,                       // chr17:7,673,802, TP53 c.818 (semaine 12), 0-based
    domain: await P.ev("window.__discover.domain()"),
    focus: await P.ev("window.__discover.focus()"),
    region: await P.ev("window.__discover.region()"),
    genes: {},
  };
  ctx.bp = (ctx.domain.end - ctx.domain.start) / (ctx.focus.j - ctx.focus.i + 1);
  const results = [];
  const shots = {};
  const glossary = { terms: 0, hover: 0, focus: 0, misses: [] };

  /** Chaque terme balisé de l'écran courant : survol, puis focus clavier (Tab). */
  const checkTerms = async (screen) => {
    const terms = await P.ev(`[...document.querySelectorAll("#story .term")].map((e) => { const r = e.getBoundingClientRect(); return { term: e.dataset.term, x: r.left + Math.min(6, r.width / 2), y: r.top + r.height / 2 }; })`);
    for (const t of terms) {
      glossary.terms++;
      await P.mouse(t.x, t.y);
      await sleep(30);
      const hov = await P.ev(`(() => { const tip = document.getElementById("tip"); return tip.style.display === "block" ? tip.textContent : ""; })()`);
      if (hov.startsWith(t.term) && hov.length > t.term.length + 10) glossary.hover++;
      else glossary.misses.push(`${screen} survol ${t.term}`);
      await P.mouse(2, 2);
    }
    // Au clavier : Tab depuis le haut du panneau, jusqu'à avoir parcouru les termes.
    await P.ev(`document.activeElement && document.activeElement.blur(); document.querySelector("#story").setAttribute("tabindex", "-1"); document.querySelector("#story").focus()`);
    const seen = new Set();
    for (let k = 0; k < terms.length + 30 && seen.size < terms.length; k++) {
      await P.key("Tab", "Tab", 9);
      const st = await P.ev(`(() => { const a = document.activeElement; const tip = document.getElementById("tip"); return a && a.classList.contains("term") ? { term: a.dataset.term, shown: tip.style.display === "block" ? tip.textContent : "" } : null; })()`);
      if (!st) continue;
      const key = `${st.term}#${seen.size}`;
      seen.add(key);
      if (st.shown.startsWith(st.term)) glossary.focus++;
      else glossary.misses.push(`${screen} focus ${st.term}`);
    }
    if (seen.size < terms.length) glossary.misses.push(`${screen} : ${terms.length - seen.size} termes jamais atteints au clavier`);
    await P.ev(`document.activeElement && document.activeElement.blur()`);
  };

  await checkTerms("intro");
  // 1 — Les étapes, au clavier : →.
  for (const id of steps) {
    await P.key("ArrowRight", "ArrowRight", 39);
    await P.arrived();
    const cur = await P.ev("window.__discover.current()");
    const c = await P.ev("window.__discover.centre()");
    const j = judgeCentre(id, c, ctx);
    const cats = await P.ev("window.__discover.categories()"), declared = await P.ev("window.__discover.declared()");
    const contacts = id === "tad" ? await (async () => { for (let k = 0; k < 100; k++) { const x = await P.ev("window.__discover.contacts()"); if (x) return x; await sleep(100); } return null; })() : null;
    shots[id] = await P.shot();
    results.push({ id, cur, centre: j, evidence: { declared, cats, ok: cats.includes(declared) }, contacts, width: await P.ev("window.__discover.widthUm()") });
    await checkTerms(id);
  }
  // La fin, puis le retour à l'introduction.
  await P.key("ArrowRight", "ArrowRight", 39);
  await sleep(300);
  const atEnd = await P.ev("window.__discover.current()");
  const quizShown = await P.ev(`!!document.querySelector("#quiz textarea") && document.querySelectorAll("#quiz fieldset").length === 4`);
  await P.key("ArrowRight", "ArrowRight", 39);
  await sleep(300);
  const backToIntro = await P.ev("window.__discover.current()");

  // Les haltes : Tab jusqu'au bouton, Entrée.
  for (const id of halts) {
    let reached = false;
    for (let k = 0; k < 80; k++) {
      await P.key("Tab", "Tab", 9);
      if (await P.ev(`document.activeElement && document.activeElement.dataset.halt === ${JSON.stringify(id)}`)) { reached = true; break; }
    }
    if (reached) await P.key("Enter", "Enter", 13);
    await P.arrived();
    const c = await P.ev("window.__discover.centre()");
    const name = (await P.ev(`document.querySelector("#story h2").textContent`));
    // Le gène de la halte : son TSS, lu par l'index de la page (comme la recherche de la semaine 16).
    const gene = { hbb: "HBB", fto: "FTO", hoxd: "HOXD13", xi: "XIST" }[id];
    // Le TSS du gène, lu par l'index de la page (comme la recherche de la semaine 16).
    const tss = await P.ev(`window.__discover.resolve(${JSON.stringify(gene)})`);
    ctx.genes[id] = tss && { name: gene, chrom: tss.chrom, tss: tss.pos };
    const j = judgeCentre(id, c, ctx);
    const cats = await P.ev("window.__discover.categories()"), declared = await P.ev("window.__discover.declared()");
    let facts = "";
    for (let k = 0; k < 50 && !facts; k++) { facts = await P.ev(`document.querySelector("#halt-facts")?.textContent ?? ""`); if (!facts) await sleep(100); }
    shots[id] = await P.shot();
    results.push({ id, reached, title: name, centre: j, evidence: { declared, cats, ok: cats.includes(declared) }, facts, factsOk: facts.includes(gene) });
    await checkTerms(id);
  }
  await P.close();

  // 5 — Chaque arrêt rejoué par l'URL, dans un onglet neuf.
  const replay = [];
  for (const id of stops) {
    const Q = await open(`${base}#t=${id}`);
    await Q.arrived();
    if (id === "tad") for (let k = 0; k < 100 && !(await Q.ev("window.__discover.contacts()")); k++) await sleep(100);
    const img = await Q.shot();
    replay.push({ id, current: await Q.ev("window.__discover.current()"), same: img === shots[id] });
    if (img !== shots[id]) { writeFileSync(resolve(OUT, `${id}-a.png`), Buffer.from(shots[id], "base64")); writeFileSync(resolve(OUT, `${id}-b.png`), Buffer.from(img, "base64")); }
    await Q.close();
  }

  const centreOk = results.filter((r) => r.centre.ok).length;
  const evOk = results.filter((r) => r.evidence.ok).length;
  const tad = results.find((r) => r.id === "tad");
  const factsOk = results.filter((r) => r.factsOk).length;
  const replayOk = replay.filter((r) => r.same && r.current === r.id).length;
  const loopOk = atEnd === "fin" && quizShown && backToIntro === "intro" && results.filter((r) => r.reached !== false).length === 11;
  const crit = {
    parcours: centreOk === 11 && loopOk,
    glossaire: glossary.terms > 0 && glossary.hover === glossary.terms && glossary.focus === glossary.terms,
    evidence: evOk === 11,
    faits: !!tad?.contacts && tad.contacts.domain[0] === ctx.domain.start && factsOk === 4,
    url: replayOk === 11,
  };
  writeFileSync(resolve(OUT, "discover.json"), JSON.stringify({ ctx, results, replay, glossary, atEnd, backToIntro, crit }, null, 1));
  console.log("\n  Mode Découverte — semaine 19");
  console.log(`  1  parcours    centre juste ${centreOk}/11 · fin → quiz ${quizShown ? "oui" : "non"} → introduction ${backToIntro === "intro" ? "oui" : "non"} · haltes au clavier ${results.filter((r) => r.reached).length}/4 ${crit.parcours ? "✓" : "✗"}`);
  for (const r of results) console.log(`                 ${r.centre.ok ? "✓" : "✗"} ${r.id.padEnd(14)} ${r.centre.why}${r.width ? ` · champ ${r.width.toFixed(r.width < 0.1 ? 4 : 2)} µm` : ""}`);
  console.log(`  2  glossaire   ${glossary.terms} termes · survol ${glossary.hover} · focus ${glossary.focus}${glossary.misses.length ? ` · manqués : ${glossary.misses.slice(0, 6).join(" ; ")}` : ""} ${crit.glossaire ? "✓" : "✗"}`);
  console.log(`  3  évidence    ${evOk}/11 : ${results.map((r) => `${r.id} ${r.evidence.declared}${r.evidence.ok ? "" : ` ≠ ${r.evidence.cats.join("/")}`}`).join(" · ")} ${crit.evidence ? "✓" : "✗"}`);
  console.log(`  4  faits       carte du TAD ${tad?.contacts ? `${tad.contacts.bins} × ${tad.contacts.bins} bins de ${tad.contacts.resolution / 1000} kb` : "absente"} · haltes ${factsOk}/4 nomment leur gène ${crit.faits ? "✓" : "✗"}`);
  console.log(`  5  URL         ${replayOk}/11 arrêts rejoués au pixel près${replay.filter((r) => !r.same).length ? ` · différents : ${replay.filter((r) => !r.same).map((r) => r.id).join(", ")}` : ""} ${crit.url ? "✓" : "✗"}`);
  const pass = Object.values(crit).every(Boolean);
  console.log(`\n  verdict : ${pass ? "\x1b[32mfait\x1b[0m" : "\x1b[31méchec\x1b[0m"}`);
  process.exitCode = pass ? 0 : 1;
} finally {
  browser.close?.();
  await shutdown(proc, profile);
  server.close();
}
