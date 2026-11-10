/** « Un .mcool importé ressort en structure, et le manifeste de session la reproduit à
 *  l'identique. » (semaine 17)
 *
 *  La page `recherche.html`, servie par le service local (`geno research serve`), de vrais
 *  gestes CDP. Deux fichiers : un `.mcool` **planté** à structure connue (deux chromosomes, trois
 *  résolutions), et, s'il est là, un vrai `.mcool` public (Micro-C HFF, chr2 et chr17). Critères
 *  fixés avant de mesurer :
 *
 *    1  import      ce que la page affiche du fichier (résolutions, bins, chromosomes) est ce que
 *                   cooler y lit — fichier déposé (planté) et chemin local (réel) ;
 *    2  structure   pour chaque chromosome planté, à 50 kb : la page affiche la structure — toutes
 *                   les billes du manifeste dans le .g3d, au moins la moitié visibles dans le
 *                   tampon d'identifiants sous deux orientations, aucun identifiant étranger — et
 *                   elle est la structure plantée : nRMSD ≤ 0,20 après alignement (réflexion
 *                   permise) ; témoin : la vérité battue, nRMSD ≥ 0,8 ;
 *    3  manifeste   12 sessions (les deux fichiers, résolutions, régions, alpha, équilibrages
 *                   variés), manifestes gardés ; **service arrêté et relancé**, page rouverte,
 *                   fichiers réimportés, chaque manifeste rejoué par la page : mêmes octets,
 *                   12/12. Témoins : alpha + 0,01 dans le manifeste → « différent » ; le
 *                   manifeste rejoué sur l'autre fichier → refusé ; un fichier qui n'est pas du
 *                   JSON → refusé ;
 *    4  refus       une région trop grande pour la reconstruction dense est refusée avec sa raison.
 *
 *    node tools/research.mjs [vrai.mcool]
 */

import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/research-test");
mkdirSync(OUT, { recursive: true });
const REAL = process.argv[2] && existsSync(process.argv[2]) ? resolve(process.argv[2]) : null;
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
const env = { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") };
const geno = (...args) => execFileSync(python, ["-m", "geno_pipeline", "research", ...args], { cwd: resolve(ROOT, "pipeline"), env, encoding: "utf8", stdio: ["ignore", "pipe", "pipe"] });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const PLANTED = resolve(OUT, "plante.mcool");
if (!existsSync(PLANTED)) geno("plant", PLANTED);
const work = mkdtempSync(join(tmpdir(), "geno-recherche-"));

function startService() {
  const proc = spawn(python, ["-m", "geno_pipeline", "research", "serve", "--port", "0", "--work", work], { cwd: resolve(ROOT, "pipeline"), env, stdio: ["ignore", "pipe", "inherit"] });
  return new Promise((done, fail) => {
    let buf = "";
    proc.stdout.on("data", (d) => {
      buf += d;
      const m = /http:\/\/127\.0\.0\.1:(\d+)\/recherche\.html/.exec(buf);
      if (m) done({ proc, url: m[0], pid: proc.pid });
    });
    proc.on("exit", (c) => fail(new Error(`le service s'est arrêté (${c})`)));
  });
}

async function openPage(browser, url) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: 1280, height: 760, deviceScaleFactor: 1, mobile: false });
  await s("DOM.enable");
  await s("Page.navigate", { url });
  const ev = async (expr) => {
    const r = await s("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true });
    if (r.exceptionDetails) throw new Error(`${expr.slice(0, 80)} : ${r.exceptionDetails.exception?.description ?? r.exceptionDetails.text}`);
    return r.result.value;
  };
  for (let k = 0; k < 600 && !(await ev("!!window.__research").catch(() => false)); k++) await sleep(100);
  const idle = async () => {
    await sleep(150);
    for (let k = 0; k < 6000 && (await ev("window.__research.busy()")); k++) await sleep(50);
    await sleep(100);
  };
  const setFile = async (selector, path) => {
    const { root } = await s("DOM.getDocument", {});
    const { nodeId } = await s("DOM.querySelector", { nodeId: root.nodeId, selector });
    await s("DOM.setFileInputFiles", { nodeId, files: [path] });
    await idle();
  };
  const click = async (selector) => {
    const [x, y] = await ev(`(() => { const e = document.querySelector(${JSON.stringify(selector)}); e.scrollIntoView({ block: "center" }); const r = e.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; })()`);
    for (const type of ["mouseMoved", "mousePressed", "mouseReleased"]) await s("Input.dispatchMouseEvent", { type, x, y, button: type === "mouseMoved" ? "none" : "left", buttons: type === "mousePressed" ? 1 : 0, clickCount: 1 });
  };
  /** Un champ rempli comme par l'utilisateur : valeur, puis les événements que la page écoute. */
  const fill = (selector, value) => ev(`(() => { const e = document.querySelector(${JSON.stringify(selector)}); e.value = ${JSON.stringify(String(value))}; e.dispatchEvent(new Event("input", { bubbles: true })); e.dispatchEvent(new Event("change", { bubbles: true })); })()`);
  const importPath = async (path) => { await fill("#path", path); await click("#open"); await idle(); };
  const run = async (p) => {
    await fill("#resolution", p.resolution);
    await fill("#region", p.region);
    await fill("#balance", p.balance);
    await fill("#alpha", p.alpha);
    const t0 = Date.now();
    await click("#run");
    await idle();
    return { ms: Date.now() - t0, last: await ev("window.__research.last()"), message: await ev("window.__research.message()") };
  };
  return { s, ev, idle, setFile, click, fill, importPath, run, targetId };
}

let out = null;
let svc = await startService();
const { proc: chrome, url: wsUrl, profile } = await launch();
try {
  const browser = await cdp(wsUrl);
  let P = await openPage(browser, svc.url);

  // — 1. Import —
  const same = (page, py) => JSON.stringify([page.chroms, page.resolutions.map((r) => [r.resolution, r.bins, r.balanced, r.contacts])]) ===
    JSON.stringify([py.chroms, py.resolutions.map((r) => [r.resolution, r.bins, r.balanced, r.contacts])]);
  const imports = [];
  await P.setFile("#file", PLANTED);
  const filePlanted = await P.ev("window.__research.file()");
  imports.push({ file: "planté (déposé)", ok: !!filePlanted && same(filePlanted.inspect, JSON.parse(geno("inspect", PLANTED))), shown: await P.ev("document.getElementById('info').textContent") });

  // — 2. Structure, sur le planté —
  const structures = [];
  const sessions = [];             // { file, params, manifest }
  for (const region of ["chrP1", "chrP2"]) {
    const params = { resolution: 50000, region, alpha: 0.35, balance: "aucun" };
    const r = await P.run(params);
    const seen = await P.ev("window.__research.seen()");
    const beads = await P.ev("window.__research.beads().length");
    const judge = JSON.parse(geno("judge", PLANTED, "--out", join(work, "sessions", r.last.session)));
    const m = r.last.manifest;
    structures.push({ region, ms: r.ms, beads, manifestBeads: m.stats.beads, g3dBeads: judge.beads, ...seen, nrmsd: judge.nrmsd, shuffled: judge.nrmsd_shuffled,
                      ok: beads === m.stats.beads && judge.beads === beads && seen.distinct >= beads / 2 && seen.max <= beads && seen.covered > 0.02 * seen.W * seen.H && judge.nrmsd <= 0.2 && judge.nrmsd_shuffled >= 0.8 });
    sessions.push({ file: "planté", params, manifest: m });
  }

  // — 3. D'autres sessions, variées —
  const more = [
    { resolution: 100000, region: "chrP1", alpha: 0.35, balance: "fichier" },
    { resolution: 200000, region: "chrP2", alpha: 0.5, balance: "ice" },
    { resolution: 50000, region: "chrP1:5,000,001-20,000,000", alpha: 0.3, balance: "aucun" },
    { resolution: 100000, region: "génome", alpha: 0.4, balance: "aucun" },
    { resolution: 200000, region: "chrP1", alpha: 0.25, balance: "fichier" },
    { resolution: 50000, region: "chrP2", alpha: 0.45, balance: "ice" },
  ];
  for (const params of more) {
    const r = await P.run(params);
    if (!r.last || JSON.stringify(r.last.manifest.params.region) !== JSON.stringify(params.region)) throw new Error(`session non faite : ${JSON.stringify(params)} — ${r.message}`);
    sessions.push({ file: "planté", params, manifest: r.last.manifest });
  }
  let refusal = null;
  const timings = [];
  if (REAL) {
    await P.importPath(REAL);
    const fileReal = await P.ev("window.__research.file()");
    imports.push({ file: "réel (chemin local)", ok: !!fileReal && same(fileReal.inspect, JSON.parse(geno("inspect", REAL))), shown: await P.ev("document.getElementById('info').textContent") });
    for (const params of [
      { resolution: 100000, region: "chr17", alpha: 0.35, balance: "fichier" },
      { resolution: 1000000, region: "génome", alpha: 0.35, balance: "fichier" },
      { resolution: 10000, region: "chr17:7,000,001-9,000,000", alpha: 0.35, balance: "fichier" },
      { resolution: 100000, region: "chr2", alpha: 0.35, balance: "ice" },
    ]) {
      const before = await P.ev("window.__research.last()?.session ?? null");
      const r = await P.run(params);
      if (!r.last || r.last.session === before) throw new Error(`session réelle non faite : ${JSON.stringify(params)} — ${r.message}`);
      const seen = await P.ev("window.__research.seen()");
      sessions.push({ file: "réel", params, manifest: r.last.manifest });
      timings.push({ params, ms: r.ms, beads: r.last.manifest.stats.beads, masked: r.last.manifest.stats.masked, r: r.last.manifest.stats.r_log_distance, distinct: seen.distinct });
    }
    // — 4. Refus —
    const before = await P.ev("window.__research.last().session");
    const r = await P.run({ resolution: 10000, region: "chr2", alpha: 0.35, balance: "fichier" });
    refusal = { same: r.last.session === before, message: r.message[1] };
  } else {
    const before = await P.ev("window.__research.last().session");
    const r = await P.run({ resolution: 50000, region: "chrP1:1-200,000", alpha: 0.35, balance: "aucun" });
    refusal = { same: r.last.session === before, message: r.message[1] };
  }

  // — Service arrêté et relancé : rien ne survit que les fichiers et les manifestes —
  svc.proc.kill();
  await sleep(500);
  rmSync(join(work, "sessions"), { recursive: true, force: true });
  const oldPid = svc.pid;
  svc = await startService();
  P = await openPage(browser, svc.url);
  const replays = [];
  const controls = {};
  const tokens = {};
  for (const which of ["planté", "réel"]) {
    const list = sessions.filter((x) => x.file === which);
    if (!list.length) continue;
    if (which === "planté") await P.setFile("#file", PLANTED);
    else await P.importPath(REAL);
    tokens[which] = await P.ev("window.__research.file().token");
    for (const [k, x] of list.entries()) {
      const path = join(OUT, `manifeste-${which === "planté" ? "p" : "r"}${k}.json`);
      writeFileSync(path, JSON.stringify(x.manifest));
      const t0 = Date.now();
      await P.setFile("#replay", path);
      const verdict = await P.ev("window.__research.verdict()");
      const last = await P.ev("window.__research.last()");
      const bytes = readFileSync(join(work, "sessions", last.session, "structure.g3d"));
      const { createHash } = await import("node:crypto");
      const sha = createHash("sha256").update(bytes).digest("hex");
      replays.push({ file: which, params: x.params, ms: Date.now() - t0, verdict, identical: /identique/.test(verdict) && sha === x.manifest.output.sha256 && last.manifest.output.sha256 === sha });
    }
    if (which === "planté") {
      // Témoins, sur le planté.
      const m = JSON.parse(JSON.stringify(list[0].manifest));
      m.params.alpha = Math.round((m.params.alpha + 0.01) * 100) / 100;
      const pa = join(OUT, "manifeste-alpha.json");
      writeFileSync(pa, JSON.stringify(m));
      await P.setFile("#replay", pa);
      controls.alpha = await P.ev("window.__research.verdict()");
      const pj = join(OUT, "pas-du-json.json");
      writeFileSync(pj, "ceci n'est pas un manifeste");
      await P.setFile("#replay", pj);
      controls.notJson = await P.ev("window.__research.verdict()");
      const real = sessions.find((x) => x.file === "réel");
      if (real) {
        const pr = join(OUT, "manifeste-autre-fichier.json");
        writeFileSync(pr, JSON.stringify(real.manifest));
        await P.setFile("#replay", pr);
        controls.otherFile = await P.ev("window.__research.verdict()");
      }
    }
  }

  // — Jugement —
  const okImport = imports.every((i) => i.ok);
  const okStruct = structures.every((x) => x.ok);
  const okReplay = replays.length === sessions.length && replays.every((r) => r.identical);
  const okControls = /Différent/.test(controls.alpha ?? "") && /JSON/.test(controls.notJson ?? "") && (!REAL || /pas celui du manifeste/.test(controls.otherFile ?? ""));
  const okRefusal = refusal.same && /résolution plus grossière|trop peu/.test(refusal.message ?? "");
  console.log("\nMode Recherche — semaine 17\n");
  for (const i of imports) console.log(`  1  import      ${i.ok ? "égal à cooler" : "DIFFÉRENT"} — ${i.file}`);
  if (!REAL) console.log("                 (pas de vrai .mcool fourni : node tools/research.mjs chemin/vers/test.mcool)");
  for (const x of structures) console.log(`  2  structure   ${x.region} : ${x.beads} billes (manifeste ${x.manifestBeads}, .g3d ${x.g3dBeads}), ${x.distinct} vues, ${(100 * x.covered / (x.W * x.H)).toFixed(1)} % de l'image · nRMSD ${x.nrmsd.toFixed(3)} (vérité battue : ${x.shuffled.toFixed(2)}) · ${x.ms} ms  ${x.ok ? "✓" : "✗"}`);
  console.log(`  3  manifeste   ${replays.filter((r) => r.identical).length}/${sessions.length} rejoués à l'octet près, service relancé (pid ${oldPid} → ${svc.pid})`);
  for (const r of replays.filter((x) => !x.identical)) console.log(`     écart : ${JSON.stringify(r)}`);
  console.log(`                 témoins : alpha + 0,01 → « ${(controls.alpha ?? "").slice(0, 40)}… » · pas du JSON → « ${controls.notJson} »${REAL ? ` · autre fichier → « ${controls.otherFile} »` : ""}`);
  console.log(`  4  refus       ${refusal.same ? "aucune structure produite" : "UNE STRUCTURE A ÉTÉ PRODUITE"} — « ${refusal.message} »`);
  for (const t of timings) console.log(`     réel  ${t.params.region} à ${t.params.resolution} (${t.params.balance}) : ${t.beads} billes, ${t.masked} bins écartés, r ${t.r.toFixed(3)}, ${t.distinct} vues, ${t.ms} ms`);
  const pass = okImport && okStruct && okReplay && okControls && okRefusal;
  console.log(`\n  verdict : ${pass ? "\x1b[32mfait\x1b[0m" : "\x1b[31méchec\x1b[0m"}${REAL ? "" : " (sans fichier réel)"}`);
  out = { real: !!REAL, imports, structures, sessions: sessions.map((x) => ({ file: x.file, params: x.params, sha256: x.manifest.output.sha256, stats: x.manifest.stats, seconds: x.manifest.run.seconds })),
          replays, controls, refusal, timings, versions: sessions[0].manifest.versions, pass };
} finally {
  await shutdown(chrome, profile).catch(() => {});
  svc.proc.kill();
  rmSync(work, { recursive: true, force: true });
}
writeFileSync(resolve(OUT, "result.json"), JSON.stringify(out, null, 1));
process.exit(out?.pass ? 0 : 1);
