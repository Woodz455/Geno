/** L'enregistrement : la plongée entière, du noyau à TP53 c.818, en vidéo.
 *
 *  Chaque image est rendue sur demande dans Chromium headless (`dive.html?record=N`), avec
 *  sa surimpression (lettres, étiquette, légende), capturée par CDP et passée à ffmpeg
 *  (H.264, 30 i/s). Puis la continuité est mesurée **sur la vidéo elle-même** : écart moyen
 *  de luminance d'une image à la suivante (`tblend` + `signalstats` de ffmpeg), et même
 *  détecteur de sauts que `tools/dive.mjs` — trois fois la médiane des 20 voisines.
 *
 *  Dit d'avance : une vidéo rendue image par image ne peut pas montrer un chargement, puisque
 *  chaque image attend d'être complète. Le chargement se mesure en temps réel, à part
 *  (`tools/stream.mjs`).
 *
 *    node tools/record.mjs [secondes] [sortie.mp4]      défaut : 45 s, dist/plongee-tp53.mp4
 */

import { execFile, spawn } from "node:child_process";
import { resolve } from "node:path";
import { promisify } from "node:util";
import { cdp, launch, shutdown } from "./cdp.mjs";
import { RATIO, ratios } from "./jumps.mjs";
import { serve } from "./serve.mjs";

const SECONDS = Number(process.argv[2] ?? 45);
const OUT = resolve(process.argv[3] ?? "dist/plongee-tp53.mp4");
const FPS = 30;
const HOLD_START = 1.5, HOLD_END = 3;
const N = Math.round(SECONDS * FPS);
const [W, H] = [1280, 720];

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
const ff = spawn("ffmpeg", ["-y", "-loglevel", "error", "-f", "image2pipe", "-framerate", String(FPS), "-i", "-",
  "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-preset", "medium", "-movflags", "+faststart", OUT], { stdio: ["pipe", "inherit", "inherit"] });
const t0 = Date.now();
try {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  await s("Page.navigate", { url: `${base}/dive.html?f=data/gm12878.g3d&record=${N}&w=${W}&h=${H}` });
  let ready = null;
  for (let i = 0; i < 900 && !ready; i++) {
    await new Promise((x) => setTimeout(x, 200));
    const { result } = await s("Runtime.evaluate", { expression: "window.__geno ? JSON.stringify(window.__geno) : null", returnByValue: true });
    if (result.value) ready = JSON.parse(result.value);
  }
  if (!ready?.ok) throw new Error(`page non prête : ${JSON.stringify(ready)}`);
  const write = (buf) => new Promise((ok) => (ff.stdin.write(buf) ? ok() : ff.stdin.once("drain", ok)));
  let shot = null;
  for (let k = 0; k < N; k++) {
    await s("Runtime.evaluate", { expression: `window.__record.render(${k})`, returnByValue: true });
    const { data } = await s("Page.captureScreenshot", { format: "png" });
    shot = Buffer.from(data, "base64");
    const reps = k === 0 ? Math.round(HOLD_START * FPS) : k === N - 1 ? Math.round(HOLD_END * FPS) : 1;
    for (let r = 0; r < reps; r++) await write(shot);
    if (k % 60 === 0) console.log(`  image ${k}/${N} · ${((Date.now() - t0) / 60000).toFixed(1)} min`);
  }
  await browser.send("Target.closeTarget", { targetId });
} finally {
  ff.stdin.end();
  browser.close();
  server.close();
  await shutdown(proc, profile);
}
await new Promise((ok) => ff.on("close", ok));

// Continuité mesurée sur la vidéo encodée : écart de luminance moyen d'une image à la suivante.
const { stderr } = await promisify(execFile)("ffmpeg", ["-hide_banner", "-i", OUT, "-vf",
  "tblend=all_mode=difference,signalstats,metadata=print:key=lavfi.signalstats.YAVG", "-f", "null", "-"], { maxBuffer: 256 * 1024 * 1024 });
const diffs = [...stderr.matchAll(/lavfi\.signalstats\.YAVG=([\d.]+)/g)].map((m) => Number(m[1]));
// Les images tenues au début et à la fin ont un écart nul : on ne juge que le mouvement.
const a = Math.round(HOLD_START * FPS), b = diffs.length - Math.round(HOLD_END * FPS);
const moving = diffs.slice(a, b);
const rs = ratios(moving);
const jumps = rs.flatMap((r, k) => (r > RATIO && moving[k] > 0.5 ? [{ frame: k + a, ratio: r, diff: moving[k] }] : []));
const { stdout: probe } = await promisify(execFile)("ffprobe", ["-v", "error", "-show_entries", "format=duration,size", "-of", "json", OUT]);
const info = JSON.parse(probe).format;
console.log(`\n  ${OUT}  ${(Number(info.size) / 1e6).toFixed(1)} Mo · ${Number(info.duration).toFixed(1)} s · ${W}×${H} · ${FPS} i/s · ${N} images rendues en ${((Date.now() - t0) / 60000).toFixed(1)} min`);
console.log(`  continuité sur la vidéo : ${jumps.length} saut(s) · pire rapport ${Math.max(...rs).toFixed(2)} · écart médian ${[...moving].sort((x, y) => x - y)[moving.length >> 1].toFixed(2)}/255`);
for (const j of jumps.slice(0, 10)) console.log(`    image ${j.frame} : ×${j.ratio.toFixed(1)} (écart ${j.diff.toFixed(2)})`);
process.stdout.write(`JSON ${JSON.stringify({ out: OUT, size: Number(info.size), duration: Number(info.duration), frames: N, jumps, worst: Math.max(...rs) })}\n`);
process.exit(jumps.length ? 1 : 0);
