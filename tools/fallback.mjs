/** Semaine 22 — le repli, vérifié : sans WebGL, chaque page qui ouvre une scène montre un message
 *  qui nomme ce qui manque, sans exception ni écran noir ; les pages pleines dessinent en plus
 *  l'aperçu en 2D. Chromium est lancé deux fois : `--disable-webgl`, puis `--disable-3d-apis`.
 *  Critère (`docs/S22_BUDGETS.md`) : 100 % des pages.
 *
 *    node tools/fallback.mjs
 */

import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
if (!existsSync(resolve(ROOT, "pipeline/data/g3d/gm12878.g3d"))) {
  const want = readFileSync(resolve(ROOT, "fixtures/gm12878.g3d.sha256"), "utf8").trim().split(/\s+/)[0];
  if (createHash("sha256").update(readFileSync(resolve(ROOT, "fixtures/gm12878.g3d"))).digest("hex") !== want) process.exit(2);
  process.env.GENO_DATA = resolve(ROOT, "fixtures");
}
const { serve } = await import("./serve.mjs");

const PAGES = [
  ["first.html", false], ["dive.html", true], ["decouverte.html", true], ["decouverte.html?lang=en", true],
  ["recherche.html", true], ["variant.html", true], ["vues.html", true],
];
const MODES = [["--disable-webgl"], ["--disable-3d-apis"]];

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const rows = [];
for (const flags of MODES) {
  const { proc, url, profile } = await launch(flags);
  const browser = await cdp(url);
  try {
    for (const [page, still] of PAGES) {
      const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
      const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
      const s = (m, p) => browser.send(m, p, sessionId);
      const errors = [];
      browser.on("Runtime.exceptionThrown", (p, sid) => { if (sid === sessionId) errors.push(p.exceptionDetails?.exception?.description ?? p.exceptionDetails?.text); });
      await s("Runtime.enable");
      await s("Page.enable");
      await s("Page.navigate", { url: `${base}/${page}` });
      const ev = async (e) => (await s("Runtime.evaluate", { expression: e, returnByValue: true })).result.value;
      let st = null;
      for (let k = 0; k < 120 && !st; k++) {
        await new Promise((r) => setTimeout(r, 250));
        st = await ev(`(() => { const b = document.getElementById("geno-fallback"); const n = window.__fallbackDrawn; return b && (${still} ? n !== undefined : true) ? { text: b.textContent, drawn: n ?? 0 } : null; })()`);
      }
      const ok = !!st && /WebGL2/.test(st.text) && (!still || st.drawn > 0) && errors.length === 0;
      rows.push({ mode: flags[0], page, ok, drawn: st?.drawn ?? 0, errors, text: st?.text.slice(0, 90) ?? "(pas de repli)" });
      await browser.send("Target.closeTarget", { targetId });
    }
  } finally {
    browser.close?.();
    await shutdown(proc, profile).catch(() => {});
  }
}
server.close();

console.log("\n  Repli sans WebGL — semaine 22");
for (const r of rows) console.log(`  ${r.ok ? "✓" : "✗"} ${r.mode.padEnd(18)} ${r.page.padEnd(26)} ${r.drawn ? `aperçu 2D ${r.drawn} billes · ` : ""}${r.errors.length ? `exceptions : ${r.errors.join(" | ").slice(0, 120)}` : r.text}`);
const pass = rows.every((r) => r.ok);
console.log(`\n  ${rows.filter((r) => r.ok).length}/${rows.length} pages · verdict : ${pass ? "\x1b[32mfait\x1b[0m" : "\x1b[31méchec\x1b[0m"}`);
process.exitCode = pass ? 0 : 1;
