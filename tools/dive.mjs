/** Mesure de la plongée : continuité de l'image, et temps par image.
 *
 *  « Plongée continue » se traduit en une propriété mesurable : d'une image à la suivante,
 *  l'écran change d'à peu près autant tout le long du trajet. Le chemin de van Wijk et Nuij
 *  est parcouru à vitesse perçue constante, donc le mouvement seul produit un écart régulier ;
 *  une transition de niveau, de sujet ou de palette qui « saute » se voit comme un pic.
 *
 *  Un **saut** est une image dont l'écart à la précédente dépasse 3 fois la médiane de ses
 *  20 voisines. Le témoin — les mêmes transitions, en marches sèches — doit en produire ;
 *  sinon le critère ne discrimine rien et la plongée « sans saut » ne prouve rien.
 *
 *  Temps par image : SwiftShader, rasteriseur logiciel, donc une borne **haute** sans valeur
 *  pour dire « 60 fps ». Ce chiffre-là ne vient que de vrai matériel (VALIDATION.md § S10).
 *
 *    node tools/dive.mjs [images]
 */

import { existsSync, rmSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { cdp, launch } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const FRAMES = Number(process.argv[2] ?? 240);
const RATIO = 3;
const WINDOW = 10;

if (!existsSync(resolve(ROOT, "dist/dive.html"))) {
  console.error("dist/dive.html absent — lancer `pnpm build`.");
  process.exit(2);
}
if (!existsSync(resolve(ROOT, "pipeline/data/g3d/gm12878.g3d"))) {
  console.error("pipeline/data/g3d/gm12878.g3d absent — lancer `make g3d`.");
  process.exit(2);
}

const median = (xs) => {
  const s = xs.slice().sort((a, b) => a - b);
  return s.length ? s[s.length >> 1] : NaN;
};
const pct = (xs, q) => xs.slice().sort((a, b) => a - b)[Math.min(xs.length - 1, Math.floor(xs.length * q))];

/** Rapport de chaque écart à la médiane de ses voisines (la valeur elle-même exclue). */
export function ratios(diffs) {
  return diffs.map((d, k) => {
    const around = [];
    for (let j = Math.max(0, k - WINDOW); j <= Math.min(diffs.length - 1, k + WINDOW); j++) if (j !== k) around.push(diffs[j]);
    return d / Math.max(median(around), 1e-3);
  });
}

/** Images dont l'écart à la précédente dépasse RATIO × la médiane de leurs voisines. */
export function jumps(diffs) {
  return ratios(diffs).flatMap((ratio, k) =>
    ratio > RATIO && diffs[k] > 0.5 ? [{ frame: k + 1, diff: diffs[k], ratio }] : []);
}

async function run(browser, url) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Page.enable");
  await s("Page.navigate", { url });
  const deadline = Date.now() + 30 * 60_000;
  let result = null;
  while (!result && Date.now() < deadline) {
    await new Promise((r) => setTimeout(r, 200));
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

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);

const rows = [];
let failed = false;
try {
  console.log(`\n  Plongée noyau → compartiments — ${FRAMES} images, Chromium headless (SwiftShader)\n`);
  for (const [name, flag] of [["fondus (la plongée)", ""], ["témoin : transitions sèches", "&hard"]]) {
    const r = await run(browser, `${base}/dive.html?f=data/gm12878.g3d&bench=${FRAMES}${flag}`);
    if (!r.ok) {
      console.log(`  ${name} : ÉCHEC ${JSON.stringify(r).slice(0, 300)}`);
      failed = true;
      continue;
    }
    const j = jumps(r.diffs);
    rows.push({ name, hard: r.hard, jumps: j, worst: Math.max(...ratios(r.diffs)), result: r });
    const worst = Math.max(...ratios(r.diffs));
    console.log(
      `  ${name.padEnd(30)} sauts ${String(j.length).padStart(2)}   pire rapport ${worst.toFixed(2).padStart(5)}   ` +
        `écart médian ${median(r.diffs).toFixed(2)}/255   image méd. ${median(r.frame_ms).toFixed(0)} ms` +
        `  p95 ${pct(r.frame_ms, 0.95).toFixed(0)} ms   rAF méd. ${median(r.raf_ms).toFixed(0)} ms`,
    );
    for (const x of j) {
      const st = r.states[x.frame];
      console.log(
        `      image ${String(x.frame).padStart(3)} : écart ${x.diff.toFixed(2)} (×${x.ratio.toFixed(1)})  ` +
          `champ ${st.w.toFixed(1)} µm, détail ${st.lod.toFixed(2)}, sujet ${st.focus.toFixed(2)}, palette ${st.mix.toFixed(2)}`,
      );
    }
  }
} finally {
  browser.close();
  proc.kill();
  server.close();
  rmSync(profile, { recursive: true, force: true });
}

const smooth = rows.find((r) => !r.hard);
const hard = rows.find((r) => r.hard);

// Transition par transition : l'image où chaque grandeur franchit 0,5, et le rapport de
// son écart à la médiane des 20 voisines. C'est la lecture la plus directe : une transition
// invisible a le rapport du mouvement seul, ~1.
const transitions = (r) =>
  ["lod", "focus", "mix"].map((key) => {
    const s = r.result.states;
    const k = s.findIndex((x, i) => i > 0 && s[i - 1][key] < 0.5 && x[key] >= 0.5);
    if (k < 1) return { key, frame: null };
    return { key, frame: k, width: s[k].w, ratio: ratios(r.result.diffs)[k - 1] };
  });
const NAMES = { lod: "niveau de détail 3 Mb → 750 kb", focus: "autres chromosomes effacés", mix: "territoires → compartiments" };
if (smooth && hard) {
  const ts = transitions(smooth);
  const th = transitions(hard);
  console.log(`\n  ${"transition".padEnd(34)} ${"image".padStart(5)} ${"champ".padStart(8)}   fondu   sèche`);
  ts.forEach((t, i) => {
    if (t.frame === null) return;
    console.log(`  ${NAMES[t.key].padEnd(34)} ${String(t.frame).padStart(5)} ${t.width.toFixed(1).padStart(5)} µm   ${t.ratio.toFixed(2)}×   ${th[i].ratio.toFixed(2)}×`);
  });
  smooth.transitions = ts;
  hard.transitions = th;
}
if (smooth && hard) {
  const ok = smooth.jumps.length === 0 && hard.jumps.length > 0;
  console.log(
    `\n  Verdict : ${ok ? "\x1b[32mcontinue\x1b[0m" : "\x1b[31mnon démontrée\x1b[0m"} — ` +
      `${smooth.jumps.length} saut(s) dans la plongée, ${hard.jumps.length} dans le témoin, qui doit en avoir.`,
  );
  if (!ok) failed = true;
}
console.log(`  Rendu ${rows[0]?.result.gpu ?? "?"} — les temps par image sont une borne haute, pas une mesure de fps.\n`);
process.stdout.write(`JSON ${JSON.stringify(rows.map((r) => ({ name: r.name, hard: r.hard, jumps: r.jumps, worst: r.worst, transitions: r.transitions,
  frames: r.result.frames, size: r.result.size, gpu: r.result.gpu, instances: r.result.instances,
  diff_median: median(r.result.diffs), frame_ms_median: median(r.result.frame_ms),
  frame_ms_p95: pct(r.result.frame_ms, 0.95), raf_ms_median: median(r.result.raf_ms),
  diffs: r.result.diffs, states: r.result.states })))}\n`);
process.exit(failed ? 1 : 0);
