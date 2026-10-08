/** Temps du premier rendu sous réseau bridé, dans un vrai Chromium.
 *
 *  Pilotage en CDP brut sur le `WebSocket` natif de Node 22 : pas de Playwright, pas de
 *  Puppeteer, aucune dépendance. Chaque passage ouvre un onglet neuf, désactive le
 *  cache, applique le profil réseau, navigue, et attend que la page publie un résultat
 *  **vérifié** (pixels couverts, picking au centre).
 *
 *  Jamais `--virtual-time-budget` ici : il accélère les horloges de la page, et c'est
 *  le temps qu'on mesure.
 *
 *  Deux bornes, dans deux directions opposées, à écrire avec chaque chiffre :
 *  - le bridage de DevTools agit **par requête** : il n'émule ni la poignée de main TCP,
 *    ni TLS, ni le démarrage lent de TCP — c'est une borne **basse** côté réseau ;
 *  - le rendu passe par SwiftShader, un rasteriseur logiciel — une borne **haute** côté
 *    GPU.
 *
 *    node tools/loadtime.mjs [passages]
 */

import { spawn } from "node:child_process";
import { existsSync, mkdtempSync, readFileSync, rmSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const RUNS = Number(process.argv[2] ?? 10);

// Valeurs des préréglages de DevTools, recopiées de Chromium
// (front_end/core/sdk/NetworkManager.ts). Elles n'ont pas pu être relues sur la source
// depuis cet environnement ; ce sont donc ces nombres-ci, publiés avec le résultat, qui
// définissent la mesure — pas le nom du profil.
const PROFILES = {
  "sans bridage": { latency: 0, downloadThroughput: -1, uploadThroughput: -1 },
  "Fast 4G": { latency: 165, downloadThroughput: (9e6 / 8) * 0.9, uploadThroughput: (1.5e6 / 8) * 0.9 },
  "Slow 4G": { latency: 562.5, downloadThroughput: (1.6e6 / 8) * 0.9, uploadThroughput: (750e3 / 8) * 0.9 },
};

// « ajusté » : le préfixe vaut exactement `fin du premier rendu`, lue dans le préambule.
// Une page déployée avec son fichier connaît ce nombre ; un client générique, non.
const VARIANTS = [
  { name: "noyau, préfixe 16 Kio", file: "gm12878.g3d", prefix: 16384 },
  { name: "noyau, préfixe 64 Kio", file: "gm12878.g3d", prefix: 65536 },
  { name: "noyau, préfixe ajusté", file: "gm12878.g3d", prefix: "ajusté" },
  { name: "aperçu 3 Mb, préfixe 64 Kio", file: "gm12878-apercu.g3d", prefix: 65536 },
  { name: "aperçu 3 Mb, préfixe ajusté", file: "gm12878-apercu.g3d", prefix: "ajusté" },
];

function firstFrameEnd(file) {
  const b = readFileSync(resolve(ROOT, "pipeline/data/g3d", file));
  return b.readUInt32LE(16);
}

const CHROME = [
  process.env.CHROME_PATH,
  "/opt/pw-browsers/chromium-1194/chrome-linux/chrome",
  "/usr/bin/chromium",
].filter(Boolean).find((p) => existsSync(p));
if (!CHROME) {
  console.error("Chromium introuvable. Définir CHROME_PATH.");
  process.exit(2);
}

// — CDP minimal —
function cdp(url) {
  const ws = new WebSocket(url);
  let next = 1;
  const pending = new Map();
  ws.onmessage = (ev) => {
    const msg = JSON.parse(ev.data);
    if (msg.id && pending.has(msg.id)) {
      const { ok, ko } = pending.get(msg.id);
      pending.delete(msg.id);
      msg.error ? ko(new Error(`${msg.error.message} (${msg.error.code})`)) : ok(msg.result);
    }
  };
  const send = (method, params = {}, sessionId) =>
    new Promise((ok, ko) => {
      const id = next++;
      pending.set(id, { ok, ko });
      ws.send(JSON.stringify({ id, method, params, ...(sessionId ? { sessionId } : {}) }));
    });
  return new Promise((ok, ko) => {
    ws.onopen = () => ok({ send, close: () => ws.close() });
    ws.onerror = () => ko(new Error("connexion CDP impossible"));
  });
}

function launch() {
  const profile = mkdtempSync(join(tmpdir(), "geno-chrome-"));
  const proc = spawn(CHROME, [
    "--headless",
    "--no-sandbox",
    "--use-gl=angle",
    "--use-angle=swiftshader",
    "--remote-debugging-port=0",
    `--user-data-dir=${profile}`,
    "--no-first-run",
    "about:blank",
  ]);
  return new Promise((ok, ko) => {
    let err = "";
    proc.stderr.on("data", (d) => {
      err += d;
      const m = /DevTools listening on (ws:\/\/\S+)/.exec(err);
      if (m) ok({ proc, url: m[1], profile });
    });
    proc.on("exit", (code) => ko(new Error(`Chromium s'est arrêté (${code}) :\n${err}`)));
  });
}

async function measure(browser, url, net) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Network.enable");
  await s("Network.setCacheDisabled", { cacheDisabled: true });
  await s("Network.emulateNetworkConditions", { offline: false, ...net });
  await s("Page.enable");
  await s("Page.navigate", { url });
  const deadline = Date.now() + 60_000;
  let result = null;
  while (!result && Date.now() < deadline) {
    await new Promise((r) => setTimeout(r, 25));
    const { result: r } = await s("Runtime.evaluate", {
      expression: "window.__geno ? JSON.stringify(window.__geno) : null",
      returnByValue: true,
    });
    if (r.value) result = JSON.parse(r.value);
  }
  await browser.send("Target.closeTarget", { targetId });
  if (!result) throw new Error(`pas de résultat dans le délai : ${url}`);
  return result;
}

const pct = (xs, q) => xs.slice().sort((a, b) => a - b)[Math.min(xs.length - 1, Math.floor(xs.length * q))];

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const page = resolve(ROOT, "dist/first.html");
if (!existsSync(page)) {
  console.error("dist/first.html absent — lancer `pnpm build`.");
  process.exit(2);
}
const { proc, url, profile } = await launch();
const browser = await cdp(url);

const pageBytes = statSync(page).size;
console.log(`\n  Premier rendu — Chromium headless, rendu SwiftShader, ${RUNS} passages par case`);
console.log(`  page ${(pageBytes / 1024).toFixed(1)} Kio non compressée (servie en gzip)\n`);
console.log(
  `  ${"variante".padEnd(30)} ${"profil".padEnd(13)} ${"requêtes".padStart(8)} ${"octets .g3d".padStart(11)}` +
    `  ${"document".padStart(9)} ${"en-tête".padStart(8)} ${"données".padStart(8)} ${"rendu méd.".padStart(11)} ${"p90".padStart(7)}  verdict`,
);

const rows = [];
let bad = 0;
try {
  for (const v of VARIANTS) {
    if (!existsSync(resolve(ROOT, "pipeline/data/g3d", v.file))) {
      console.log(`  ${v.name.padEnd(30)} — fichier absent (make g3d${v.file.includes("apercu") ? " PREVIEW=1 …" : ""})`);
      continue;
    }
    const prefix = v.prefix === "ajusté" ? firstFrameEnd(v.file) : v.prefix;
    for (const [pname, net] of Object.entries(PROFILES)) {
      const runs = [];
      for (let i = 0; i < RUNS; i++) {
        const r = await measure(browser, `${base}/first.html?f=data/${v.file}&prefix=${prefix}`, net);
        if (!r.ok) {
          bad++;
          console.log(`  ÉCHEC du rendu vérifié : ${JSON.stringify(r).slice(0, 300)}`);
          continue;
        }
        runs.push(r);
      }
      if (!runs.length) continue;
      const med = (k) => pct(runs.map((r) => r[k]), 0.5);
      const total = med("t_first_render");
      const p90 = pct(runs.map((r) => r.t_first_render), 0.9);
      const verdict = total < 1500 ? "\x1b[32m< 1,5 s\x1b[0m" : "\x1b[31m> 1,5 s\x1b[0m";
      rows.push({ variant: v.name, profile: pname, net, prefix, runs: runs.length, median: total, p90,
                  requests: runs[0].requests, bytes: runs[0].g3d_bytes,
                  document: med("t_document"), open: med("t_open"), data: med("t_data"),
                  beads: runs[0].beads, level: runs[0].level, coverage: runs[0].coverage });
      console.log(
        `  ${v.name.padEnd(30)} ${pname.padEnd(13)} ${String(runs[0].requests).padStart(8)} ` +
          `${runs[0].g3d_bytes.toLocaleString("fr").padStart(11)}  ${med("t_document").toFixed(0).padStart(7)} ms` +
          ` ${med("t_open").toFixed(0).padStart(6)} ms ${med("t_data").toFixed(0).padStart(6)} ms` +
          ` ${total.toFixed(0).padStart(8)} ms ${p90.toFixed(0).padStart(4)} ms  ${verdict}`,
      );
    }
  }
} finally {
  browser.close();
  proc.kill();
  server.close();
  rmSync(profile, { recursive: true, force: true });
}

console.log(`\n  Profils : ${Object.entries(PROFILES).map(([k, p]) =>
  `${k} = ${p.latency} ms, ${p.downloadThroughput < 0 ? "∞" : (p.downloadThroughput / 1000).toFixed(1) + " ko/s"}`).join(" · ")}`);
console.log("  Borne basse côté réseau (bridage par requête, ni TCP ni TLS) ; borne haute côté rendu (SwiftShader).\n");
process.stdout.write(`JSON ${JSON.stringify(rows)}\n`);
process.exit(bad ? 1 : 0);
