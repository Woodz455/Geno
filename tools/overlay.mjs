/** La surimpression clinique de la plongée (semaine 15) : chaque bille peinte porte la valeur
 *  que le lecteur Python lit dans le fichier de fiches — P/LP ClinVar, associations GWAS, LOEUF
 *  minimal —, aperçu et noyau, toutes les billes. Puis le temps de peindre chaque couche.
 *
 *    node tools/overlay.mjs
 */

import { execFileSync } from "node:child_process";
import { existsSync, readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const expected = resolve(ROOT, "dist/overlay-test/expected.json");
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
execFileSync(existsSync(venv) ? venv : "python3", ["-m", "geno_pipeline", "sv", "clinical", "--dump", expected], {
  cwd: resolve(ROOT, "pipeline"), stdio: "inherit", env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") },
});
const E = JSON.parse(readFileSync(expected, "utf8"));
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const server = await serve(0);
const { proc, url, profile } = await launch();
let failed = 0;
try {
  const b = await cdp(url);
  const { targetId } = await b.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await b.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => b.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: 1280, height: 720, deviceScaleFactor: 1, mobile: false });
  await s("Page.navigate", { url: `http://127.0.0.1:${server.address().port}/dive.html?pick` });
  const ev = async (e) => (await s("Runtime.evaluate", { expression: e, returnByValue: true, awaitPromise: true })).result.value;
  for (let k = 0; k < 900 && !(await ev("!!(window.__pick && window.__pick.prefetch.all)")); k++) await sleep(200);
  for (const lv of Object.keys(E)) {
    for (const kind of ["clinvar", "gwas", "loeuf"]) {
      const got = await ev(`window.__pick.clinical(${JSON.stringify(lv)}, ${JSON.stringify(kind)})`);
      const want = E[lv][kind];
      let same = 0;
      for (let i = 0; i < want.length; i++) if (got[i] === want[i] || (got[i] !== null && want[i] !== null && Math.abs(got[i] - want[i]) < 1e-9)) same++;
      const ok = same === want.length;
      failed += ok ? 0 : 1;
      console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${lv.padEnd(7)} ${kind.padEnd(8)} ${same}/${want.length} billes égales à Python` +
                  `  \x1b[2m${want.filter((x) => x).length} non nulles\x1b[0m`);
    }
  }
  for (const kind of ["clinvar", "gwas", "loeuf", ""]) {
    const ms = await ev(`window.__pick.paint(${JSON.stringify(kind)})`);
    console.log(`  peinte : ${kind || "aucune"} en ${ms.toFixed(0)} ms`);
  }
} finally {
  await shutdown(proc, profile).catch(() => {});
  server.close();
}
process.exit(failed ? 1 : 0);
