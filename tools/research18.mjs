/** Semaine 18 — l'explorateur d'ensemble, les statistiques, les exports, le .hic, dans la page.
 *
 *  `recherche.html` sur le service local, gestes CDP. Critères (fixés avant cette mesure, après
 *  les premiers essais en ligne de commande) :
 *
 *    1  ensemble     20 réplicats du chromosome planté : la page n'en lit un qu'à la demande (une
 *                    requête Range de billes × 12 octets par réplicat regardé), et ce qu'elle
 *                    dessine est le réplicat du fichier ; une carte dix fois moins profonde donne
 *                    une variabilité médiane au moins 1,5 fois plus grande ;
 *    2  statistiques deux courbes dessinées ; les valeurs affichées sont celles de statistics.json ;
 *    3  exports      PDB, 3DG et glTF téléchargés du service et relus **ici, en JavaScript** :
 *                    les positions sont celles que la page dessine, à la précision du format ;
 *                    SVG : une silhouette par bille devant la caméra, et au centre de chaque bille
 *                    testée le tampon d'identifiants rend cette bille ou une bille peinte après
 *                    (critère d'origine ; la révision « géométrie seule » de la première mesure,
 *                    77/80, est retirée après correction — VALIDATION § S18, reprise) ;
 *                    ajouté à la reprise : sur une grille de points du canevas, la bille du dessus
 *                    dans le SVG est celle de l'écran, hors bords (à 1 px d'une silhouette) ;
 *                    PNG : aux dimensions du canevas, signature PNG ;
 *    4  .hic         un .hic public importé par son chemin ressort en structure ;
 *    5  manifeste    chaque session (ensemble, .hic compris), service relancé et conversions
 *                    effacées : rejouée à l'octet près, structure **et** fichiers d'accompagnement.
 *
 *    node tools/research18.mjs vrai.mcool test.hic
 */

import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync, mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/research-test");
mkdirSync(OUT, { recursive: true });
const REAL = process.argv[2] && existsSync(process.argv[2]) ? resolve(process.argv[2]) : null;
const HIC = process.argv[3] && existsSync(process.argv[3]) ? resolve(process.argv[3]) : null;
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
const env = { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") };
const py = (code) => execFileSync(python, ["-W", "ignore", "-c", code], { cwd: resolve(ROOT, "pipeline"), env, encoding: "utf8", stdio: ["ignore", "pipe", "pipe"] });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const PLANTED = resolve(OUT, "plante.mcool"), SHALLOW = resolve(OUT, "plante-peu-profond.mcool");
const plant = (path, total) => py(`import logging; logging.disable(logging.INFO)\nfrom pathlib import Path\nfrom geno_pipeline.research.planted import plant\nplant(Path(${JSON.stringify(path)}), total=${total})`);
if (!existsSync(PLANTED)) plant(PLANTED, 40_000_000);
if (!existsSync(SHALLOW)) plant(SHALLOW, 4_000_000);
const work = mkdtempSync(join(tmpdir(), "geno-recherche18-"));

function startService() {
  const proc = spawn(python, ["-m", "geno_pipeline", "research", "serve", "--port", "0", "--work", work], { cwd: resolve(ROOT, "pipeline"), env, stdio: ["ignore", "pipe", "inherit"] });
  return new Promise((done, fail) => {
    let buf = "";
    proc.stdout.on("data", (d) => { buf += d; const m = /http:\/\/127\.0\.0\.1:(\d+)\/recherche\.html/.exec(buf); if (m) done({ proc, url: m[0], base: `http://127.0.0.1:${m[1]}` }); });
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
  const idle = async () => { await sleep(150); for (let k = 0; k < 12000 && (await ev("window.__research.busy()")); k++) await sleep(50); await sleep(100); };
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
  const fill = (selector, value) => ev(`(() => { const e = document.querySelector(${JSON.stringify(selector)}); e.value = ${JSON.stringify(String(value))}; e.dispatchEvent(new Event("input", { bubbles: true })); e.dispatchEvent(new Event("change", { bubbles: true })); })()`);
  const importPath = async (path) => { await fill("#path", path); await click("#open"); await idle(); };
  const run = async (p) => {
    const before = await ev("window.__research.last()?.session ?? null");
    for (const k of ["resolution", "region", "balance", "alpha", "ensemble", "seed"]) await fill(`#${k}`, p[k] ?? 0);
    const t0 = Date.now();
    await click("#run");
    await idle();
    const last = await ev("window.__research.last()");
    if (!last || last.session === before) throw new Error(`session non faite : ${JSON.stringify(p)} — ${await ev("window.__research.message()")}`);
    return { ms: Date.now() - t0, last };
  };
  return { s, ev, idle, setFile, click, fill, importPath, run };
}

const get = async (url, range) => { const r = await fetch(url, range ? { headers: { Range: range } } : {}); return { status: r.status, buf: Buffer.from(await r.arrayBuffer()) }; };
let out = null;
let svc = await startService();
const { proc: chrome, url: wsUrl, profile } = await launch();
try {
  const browser = await cdp(wsUrl);
  let P = await openPage(browser, svc.url);
  const sessions = [];

  const T0 = Date.now(), step = (m) => console.error(`  [${((Date.now() - T0) / 1000).toFixed(0)} s] ${m}`);
  // — 1. Ensemble —
  step("ensemble");
  await P.importPath(PLANTED);
  const deep = await P.run({ resolution: 50000, region: "chrP1", alpha: 0.35, balance: "aucun", ensemble: 20, seed: 1 });
  sessions.push({ file: PLANTED, name: "planté, ensemble 20", manifest: deep.last.manifest });
  const n = deep.last.manifest.stats.beads;
  const consensus = await P.ev("window.__research.shown()");
  const reads = [];
  for (const k of [7, 20, 7, 1]) {
    await P.fill("#replica", k);
    await P.idle();
    const shown = await P.ev("window.__research.shown()");
    const file = (await get(`${svc.base}/sessions/${deep.last.session}/ensemble.f32`, `bytes=${(k - 1) * n * 12}-${k * n * 12 - 1}`)).buf;
    let worst = 0;
    for (let i = 0; i < n * 3; i++) worst = Math.max(worst, Math.abs(shown[i] - Math.fround(file.readFloatLE(4 * i) / 1000)));
    reads.push({ k, worst, replica: await P.ev("window.__research.replica()"), bytes: await P.ev("window.__research.ensembleBytes()") });
  }
  const moved = reads.every((r) => r.worst === 0 && r.replica === r.k);
  const onDemand = reads[reads.length - 1].bytes === 3 * n * 12;      // 7, 20, 1 : trois lus, le second 7 vient du cache
  const total = deep.last.manifest.files.ensemble.bytes;
  await P.fill("#replica", 0);
  await P.fill("#colour", "var");
  await P.idle();
  const varDeep = (await P.ev("window.__research.statistics()")).variability_nm.median;
  await P.importPath(SHALLOW);
  const shallow = await P.run({ resolution: 50000, region: "chrP1", alpha: 0.35, balance: "aucun", ensemble: 20, seed: 1 });
  sessions.push({ file: SHALLOW, name: "planté peu profond, ensemble 20", manifest: shallow.last.manifest });
  const varShallow = (await P.ev("window.__research.statistics()")).variability_nm.median;

  step("statistiques");
  // — 2. Statistiques —
  const stFile = JSON.parse((await get(`${svc.base}/sessions/${shallow.last.session}/statistics.json`)).buf.toString());
  const stPage = await P.ev("window.__research.statistics()");
  const charts = await P.ev("window.__research.charts()");
  const okStats = charts === 2 && JSON.stringify(stFile) === JSON.stringify(stPage) && stFile.ps_slope < 0 && stFile.rs_nu > 0;

  step("exports");
  // — 3. Exports, sur la session profonde (rejouée par son manifeste : la page la réaffiche) —
  await P.importPath(PLANTED);
  const again = await P.run({ resolution: 50000, region: "génome", alpha: 0.35, balance: "aucun", ensemble: 0, seed: 0 });
  sessions.push({ file: PLANTED, name: "planté, génome", manifest: again.last.manifest });
  const shown = await P.ev("window.__research.shown()");        // µm, float32, ce que la page dessine
  const nb = again.last.manifest.stats.beads;
  const url = (f) => `${svc.base}/sessions/${again.last.session}/${f}`;
  const pdb = (await get(url("structure.pdb"))).buf.toString();
  const unit = Number(/REMARK {3}2 UNIT_NM (\S+)/.exec(pdb)[1]);
  const atoms = pdb.split("\n").filter((l) => l.startsWith("ATOM"));
  const worstOf = (coords) => { let w = 0; for (let i = 0; i < nb * 3; i++) w = Math.max(w, Math.abs(coords[i] - shown[i] * 1000)); return w; };
  const pdbXyz = atoms.flatMap((l) => [30, 38, 46].map((a) => Number(l.slice(a, a + 8)) * unit));
  const chains = [...new Set(atoms.map((l) => l[21]))];
  const tdg = (await get(url("structure.3dg"))).buf.toString().split("\n").filter((l) => l && !l.startsWith("#")).map((l) => l.split("\t"));
  const tdgXyz = tdg.flatMap((c) => c.slice(2).map(Number));
  const g = JSON.parse((await get(url("structure.gltf"))).buf.toString());
  const bin = Buffer.from(g.buffers[0].uri.split(",")[1], "base64");
  const gltfXyz = [];
  let gltfOk = g.asset.version === "2.0" && bin.length === g.buffers[0].byteLength;
  for (const m of g.meshes) {
    const a = g.accessors[m.primitives[0].attributes.POSITION];
    gltfOk &&= a.componentType === 5126 && a.type === "VEC3" && [3, 0].includes(m.primitives[0].mode);
    for (let i = 0; i < a.count * 3; i++) gltfXyz.push(bin.readFloatLE(a.byteOffset + 4 * i) * 1000);
  }
  const beads = await P.ev("window.__research.beads()");
  const exportsRes = {
    pdb: { atoms: atoms.length, chains: chains.length, unit, worst: worstOf(pdbXyz), ok: atoms.length === nb && chains.length === 2 && worstOf(pdbXyz) <= 0.0005 * unit + 1e-3 },
    tdg: { rows: tdg.length, worst: worstOf(tdgXyz), ok: tdg.length === nb && worstOf(tdgXyz) <= 0.0005 + 1e-3 && tdg.every((c, i) => c[0] === beads[i].chrom && Number(c[1]) === beads[i].start) },
    gltf: { points: gltfXyz.length / 3, worst: worstOf(gltfXyz), ok: !!gltfOk && gltfXyz.length === nb * 3 && worstOf(gltfXyz) <= 1e-3 },
  };
  // Silhouette d'une bille dans le SVG : ellipse (ex, ey, rx, ry) tournée de rot degrés ; tolérance en pixels.
  const inEllipse = (o, x, y, tol) => {
    const t = (o.rot * Math.PI) / 180, dx = x - o.ex, dy = y - o.ey;
    const a = dx * Math.cos(t) + dy * Math.sin(t), b = -dx * Math.sin(t) + dy * Math.cos(t);
    return (a / (o.rx + tol)) ** 2 + (b / (o.ry + tol)) ** 2 <= 1;
  };
  // Polygone (masque) : contient le point, et distance du point à ses arêtes.
  const inPoly = (pl, x, y) => {
    let inside = false;
    for (let i = 0, j = pl.length - 1; i < pl.length; j = i++) {
      const [xi, yi] = pl[i], [xj, yj] = pl[j];
      if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  };
  const edgeDist = (pl, x, y) => {
    let best = Infinity;
    for (let i = 0, j = pl.length - 1; i < pl.length; j = i++) {
      const [ax, ay] = pl[j], [bx, by] = pl[i], dx = bx - ax, dy = by - ay;
      const t = Math.max(0, Math.min(1, ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy || 1)));
      best = Math.min(best, Math.hypot(x - ax - t * dx, y - ay - t * dy));
    }
    return best;
  };
  step("svg");
  const svg = await P.ev("window.__research.svg()");
  step(`svg rendu, ${svg.svg.length} caractères`);
  const circlesInSvg = (svg.svg.match(/<ellipse /g) ?? []).length;
  const masked = (svg.svg.match(/ mask="url/g) ?? []).length;
  // L'ordre de dessin : du fond vers l'avant.
  const order = [...svg.svg.matchAll(/data-bead="(\d+)"/g)].map((m) => Number(m[1]));
  const rank = new Map(order.map((id, k) => [id, k]));
  const byId = new Map(svg.circles.map((o) => [o.id, o]));
  /** La bille que le SVG montre en (x, y) : la dernière peinte dont la silhouette contient le point et
   *  dont aucun masque ne l'ôte ; `edge` si le point est à moins d'1 px d'un bord qui décide. */
  const svgTop = (x, y) => {
    for (let k = order.length - 1; k >= 0; k--) {
      const o = byId.get(order[k]);
      if (!inEllipse(o, x, y, 1)) continue;
      if (!inEllipse(o, x, y, -1)) return { edge: true };
      if (o.hide.some((pl) => edgeDist(pl, x, y) < 1)) return { edge: true };
      if (o.hide.some((pl) => inPoly(pl, x, y))) continue;
      return { id: o.id };
    }
    return { id: 0 };
  };
  const judge = async (x, y, acc, extra) => {
    const t = svgTop(x, y);
    if (t.edge) { acc.edge++; return; }
    const id = await P.ev(`window.__research.idAt(${x}, ${y})`);
    const got = rank.has(id) ? id : 0;
    acc.tested++;
    if (got === t.id) acc.agree++;
    else acc.misses.push({ x, y, svg: t.id, screen: got, ...extra });
  };
  // Au centre de 80 billes réparties dans l'ordre, puis sur une grille de 20 × 15 points du canevas :
  // la bille visible du SVG est-elle celle de l'écran ?
  const centres = { tested: 0, agree: 0, edge: 0, misses: [] };
  let seen = 0;
  for (let k = 0; k < svg.circles.length && seen < 80; k += Math.max(1, Math.floor(svg.circles.length / 80))) {
    const c = svg.circles[k];
    if (c.x < 400 || c.x > 1270 || c.y < 5 || c.y > 755) continue;        // sous le panneau ou hors champ
    seen++;
    await judge(c.x, c.y, centres, { bead: c.id });
  }
  step(`centres ${centres.agree}/${centres.tested}`);
  const grid = { tested: 0, agree: 0, edge: 0, misses: [] };
  for (let gy = 0; gy < 15; gy++) for (let gx = 0; gx < 20; gx++) await judge(405 + (gx + 0.5) * (860 / 20), 10 + (gy + 0.5) * (740 / 15), grid, {});
  writeFileSync(resolve(OUT, "svg18.json"), JSON.stringify({ masked, centres, grid }, null, 1));
  await P.click("#ex-svg");
  await P.click("#ex-png");
  await sleep(600);
  const exported = await P.ev("window.__research.exported()");
  const png = await P.ev(`(async () => { const c = document.querySelector("canvas"); const b = await new Promise((r) => c.toBlob(r, "image/png")); return [c.width, c.height]; })()`);
  // Au centre de chaque bille testée et sur la grille, la bille que le SVG montre est celle de l'écran
  // (points à moins d'1 px d'un bord qui décide : comptés, non jugés).
  exportsRes.svg = { circles: circlesInSvg, masked, centres, grid, ok: circlesInSvg === svg.circles.length && circlesInSvg === nb && centres.tested >= 40 && centres.agree === centres.tested && grid.tested >= 100 && grid.agree === grid.tested };
  exportsRes.png = { exported, canvas: png, ok: exported.some((e) => e.name === "structure.png" && e.type === "image/png" && e.bytes > 1000) && exported.some((e) => e.name === "structure.svg" && e.bytes > 1000) };

  step("hic");
  // — 4. .hic, et une session réelle avec ensemble —
  let hic = null;
  if (HIC) {
    const t0 = Date.now();
    await P.importPath(HIC);
    const f = await P.ev("window.__research.file()");
    // Le fichier d'exemple est un sous-échantillon : à 250 kb aucune paire n'a de contact ; à 2,5 Mb, une sur cinq.
    const r = await P.run({ resolution: 2500000, region: f.inspect.chroms[0][0], alpha: 0.35, balance: "aucun", ensemble: 0, seed: 0 });
    const seen = await P.ev("window.__research.seen()");
    hic = { importMs: Date.now() - t0 - r.ms, runMs: r.ms, converted: !!f.inspect.converted_from_hic, beads: r.last.manifest.stats.beads, masked: r.last.manifest.stats.masked,
            r: r.last.manifest.stats.r_log_distance, distinct: seen.distinct, by: r.last.manifest.input.converted_by, ok: !!f.inspect.converted_from_hic && seen.distinct >= r.last.manifest.stats.beads / 2 };
    sessions.push({ file: HIC, name: ".hic", manifest: r.last.manifest });
  }
  let real = null;
  if (REAL) {
    await P.importPath(REAL);
    const r = await P.run({ resolution: 10000, region: "chr17:7,000,001-9,000,000", alpha: 0.35, balance: "fichier", ensemble: 12, seed: 5 });
    const st = await P.ev("window.__research.statistics()");
    real = { ms: r.ms, beads: r.last.manifest.stats.beads, variability: st.variability_nm, radius: r.last.manifest.stats.radius_nm, ps: st.ps_slope, nu: st.rs_nu };
    sessions.push({ file: REAL, name: "réel, ensemble 12", manifest: r.last.manifest });
  }

  step("manifestes");
  // — 5. Manifestes, service relancé, sessions et conversions effacées —
  svc.proc.kill();
  await sleep(500);
  rmSync(join(work, "sessions"), { recursive: true, force: true });
  rmSync(join(work, "converted"), { recursive: true, force: true });
  svc = await startService();
  P = await openPage(browser, svc.url);
  const replays = [];
  let imported = null;
  for (const [k, x] of sessions.entries()) {
    if (imported !== x.file) { await P.importPath(x.file); imported = x.file; }
    const path = join(OUT, `manifeste18-${k}.json`);
    writeFileSync(path, JSON.stringify(x.manifest));
    await P.setFile("#replay", path);
    const verdict = await P.ev("window.__research.verdict()");
    const last = await P.ev("window.__research.last()");
    const same = Object.entries(x.manifest.files).every(([key, f]) => last.manifest.files[key]?.sha256 === f.sha256) && last.manifest.output.sha256 === x.manifest.output.sha256;
    replays.push({ name: x.name, files: Object.keys(x.manifest.files).length, verdict, identical: /identique/.test(verdict) && same });
  }
  // Témoin : une autre graine dans le manifeste → l'ensemble diffère, et la page le dit.
  const m = JSON.parse(JSON.stringify(sessions[0].manifest));
  m.params.seed += 1;
  await P.importPath(sessions[0].file);
  const pc = join(OUT, "manifeste18-graine.json");
  writeFileSync(pc, JSON.stringify(m));
  await P.setFile("#replay", pc);
  const control = await P.ev("window.__research.verdict()");

  const okEnsemble = moved && onDemand && varShallow / varDeep >= 1.5;
  const okExports = Object.values(exportsRes).every((e) => e.ok);
  const okReplay = replays.every((r) => r.identical) && /Différent/.test(control) && /ensemble/.test(control);
  console.log("\nMode Recherche — semaine 18\n");
  console.log(`  1  ensemble     ${reads.length} lectures, réplicat affiché = réplicat du fichier : ${moved ? "oui" : "NON"} · à la demande : ${reads[reads.length - 1].bytes} octets lus sur ${total} (${onDemand ? "3 réplicats, 1 du cache" : "INATTENDU"})`);
  console.log(`                  variabilité médiane ${varDeep.toFixed(0)} nm (40 M de contacts) · ${varShallow.toFixed(0)} nm (4 M) · rapport ${(varShallow / varDeep).toFixed(2)} · ensemble de 20 : ${deep.ms} ms`);
  console.log(`  2  statistiques ${okStats ? "2 courbes, valeurs = statistics.json" : "ÉCART"} · P(s) pente ${stFile.ps_slope.toFixed(2)} · R(s) exposant ${stFile.rs_nu.toFixed(2)}`);
  console.log(`  3  exports      PDB ${exportsRes.pdb.atoms} atomes, ${exportsRes.pdb.chains} chaînes, unité ${unit} nm, écart max ${exportsRes.pdb.worst.toExponential(1)} nm ${exportsRes.pdb.ok ? "✓" : "✗"} · 3DG ${exportsRes.tdg.rows} lignes, ${exportsRes.tdg.worst.toExponential(1)} nm ${exportsRes.tdg.ok ? "✓" : "✗"} · glTF ${exportsRes.gltf.points} points, ${exportsRes.gltf.worst.toExponential(1)} nm ${exportsRes.gltf.ok ? "✓" : "✗"}`);
  console.log(`                  SVG ${exportsRes.svg.circles} silhouettes pour ${nb} billes, ${masked} masquées (sécantes) · bille visible = écran : centres ${centres.agree}/${centres.tested} (${centres.edge} sur un bord), grille ${grid.agree}/${grid.tested} (${grid.edge} sur un bord) ${exportsRes.svg.ok ? "✓" : "✗"} · PNG ${JSON.stringify(exported.find((e) => e.name.endsWith("png")) ?? null)} ${exportsRes.png.ok ? "✓" : "✗"}`);
  if (hic) console.log(`  4  .hic         converti (${hic.by}) en ${hic.importMs} ms, à 2,5 Mb : ${hic.beads} billes (${hic.masked} bins écartés), r ${hic.r.toFixed(3)}, ${hic.distinct} vues, ${hic.runMs} ms ${hic.ok ? "✓" : "✗"}`);
  else console.log("  4  .hic         non mesuré (pas de fichier fourni)");
  if (real) console.log(`     réel         chr17:7–9 Mb à 10 kb, ensemble 12 : ${real.beads} billes, variabilité médiane ${real.variability.median.toFixed(0)} nm (rayon ${real.radius.toFixed(0)}), P(s) ${real.ps.toFixed(2)}, nu ${real.nu.toFixed(2)}, ${real.ms} ms`);
  console.log(`  5  manifeste    ${replays.filter((r) => r.identical).length}/${replays.length} rejoués à l'octet près, fichiers d'accompagnement compris (${replays.map((r) => r.files).join(", ")}) · témoin, autre graine : « ${control.slice(0, 90)}… »`);
  for (const r of replays.filter((x) => !x.identical)) console.log(`     écart : ${JSON.stringify(r)}`);
  const pass = okEnsemble && okStats && okExports && (!HIC || hic.ok) && okReplay && !!HIC && !!REAL;
  console.log(`\n  verdict : ${pass ? "\x1b[32mfait\x1b[0m" : "\x1b[31méchec\x1b[0m"}${HIC && REAL ? "" : " (fichiers réels manquants)"}`);
  out = { reads, total, varDeep, varShallow, ensembleMs: deep.ms, statistics: stFile, okStats, exports: exportsRes, hic, real, replays, control, pass };
} finally {
  await shutdown(chrome, profile).catch(() => {});
  svc.proc.kill();
  rmSync(work, { recursive: true, force: true });
}
writeFileSync(resolve(OUT, "result18.json"), JSON.stringify(out, null, 1));
process.exit(out?.pass ? 0 : 1);
