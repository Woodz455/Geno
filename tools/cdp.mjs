/** Chromium headless piloté en CDP brut, sur le `WebSocket` natif de Node 22.
 *
 *  Partagé par `loadtime.mjs` (semaine 9) et `dive.mjs` (semaine 10). Aucune dépendance :
 *  ni Playwright ni Puppeteer, quelques dizaines de lignes de protocole suffisent.
 */

import { spawn } from "node:child_process";
import { existsSync, mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

export const CHROME = [
  process.env.CHROME_PATH,
  "/opt/pw-browsers/chromium-1194/chrome-linux/chrome",
  "/usr/bin/chromium",
].filter(Boolean).find((p) => existsSync(p));
if (!CHROME) {
  console.error("Chromium introuvable. Définir CHROME_PATH.");
  process.exit(2);
}

// — CDP minimal —
export function cdp(url) {
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

export function launch(extra = []) {
  const profile = mkdtempSync(join(tmpdir(), "geno-chrome-"));
  const proc = spawn(CHROME, [
    "--headless",
    "--no-sandbox",
    "--use-gl=angle",
    "--use-angle=swiftshader",
    "--remote-debugging-port=0",
    `--user-data-dir=${profile}`,
    "--no-first-run",
    ...extra,
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

