/** Une image de la plongée à une étape donnée, en PNG — pour regarder, pas pour mesurer.
 *
 *    node tools/shot.mjs <dossier> [vue …]      vues : final fibre nuc helix base loop anchor
 *
 *  Chaque image passe par `dive.html?view=…&labels` : surimpression des lettres comprise.
 */

import { writeFileSync } from "node:fs";
import { resolve } from "node:path";
import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const [, , out = "dist/shots", ...views] = process.argv;
const list = views.length ? views : ["final", "fibre", "nuc", "helix", "base"];
const extra = process.env.SHOT_PARAMS ?? "";
const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
try {
  for (const v of list) {
    const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
    const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
    const s = (m, p) => browser.send(m, p, sessionId);
    await s("Emulation.setDeviceMetricsOverride", { width: 960, height: 540, deviceScaleFactor: 1, mobile: false });
    await s("Page.navigate", { url: `${base}/dive.html?f=data/gm12878.g3d&view=${v}&labels${extra}` });
    let r = null;
    for (let i = 0; i < 600 && !r; i++) {
      await new Promise((x) => setTimeout(x, 200));
      const { result } = await s("Runtime.evaluate", { expression: "window.__geno ? JSON.stringify(window.__geno) : null", returnByValue: true });
      if (result.value) r = JSON.parse(result.value);
    }
    const { data } = await s("Page.captureScreenshot", { format: "png" });
    const file = resolve(out, `${v}.png`);
    writeFileSync(file, Buffer.from(data, "base64"));
    const st = r?.state ?? {};
    console.log(`  ${v.padEnd(7)} ${file}  champ ${((st.widthUm ?? 0) * 1000).toFixed(1)} nm · fibre ${(st.fiber ?? 0).toFixed(2)} · nuc ${(st.nuc ?? 0).toFixed(2)} · hélice ${(st.helix ?? 0).toFixed(2)} · lettres ${(st.bases ?? 0).toFixed(2)}${r?.error ? ` · ERREUR ${r.error}` : ""}`);
    await browser.send("Target.closeTarget", { targetId });
  }
} finally {
  browser.close();
  server.close();
  await shutdown(proc, profile);
}
