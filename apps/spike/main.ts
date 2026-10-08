/** Spike de rendu — à exécuter sur du matériel réel.
 *
 *  La suite headless prouve que les imposteurs sont corrects. Elle ne peut rien
 *  dire de la performance : en headless, Chromium rend via SwiftShader, un
 *  rasteriseur logiciel. Cette page-ci produit les chiffres qui décident du
 *  budget de niveaux de détail de la semaine 4, et elle doit tourner sur les
 *  trois cibles : GPU desktop, iGPU portable, téléphone.
 *
 *  Ouvrir `dist/spike.html` dans un navigateur et attendre la fin.
 */

import { createTargets, ImpostorRenderer, NO_HIT, perspective, Picker } from "@geno/viewer";

const FOV = (50 * Math.PI) / 180;
const FRAMES = 90;
const WARMUP = 20;

/** Les paliers correspondent aux niveaux du budget LOD de ARCHITECTURE.md § 4. */
const STEPS = [
  { n: 6_200, label: "noyau", bin: "1 Mb" },
  { n: 24_800, label: "territoire", bin: "250 kb" },
  { n: 62_000, label: "compartiment", bin: "100 kb" },
  { n: 248_000, label: "TAD", bin: "25 kb" },
  { n: 620_000, label: "boucle", bin: "10 kb" },
  { n: 1_000_000, label: "au-delà du budget", bin: "—" },
];

const out = document.createElement("pre");
out.id = "results";
document.body.appendChild(out);
const lines: string[] = [];
const say = (s = "") => {
  lines.push(s);
  out.textContent = lines.join("\n");
};

/** Nuage sphérique : la forme d'un noyau, pas un cube. */
function nucleus(n: number, radius: number, depth: number) {
  const centers = new Float32Array(n * 3);
  const radii = new Float32Array(n);
  const ids = new Uint32Array(n);
  const colors = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) {
    const u = Math.random() * 2 - 1;
    const phi = Math.random() * Math.PI * 2;
    const r = radius * Math.cbrt(Math.random());
    const s = Math.sqrt(1 - u * u);
    centers[i * 3] = r * s * Math.cos(phi);
    centers[i * 3 + 1] = r * s * Math.sin(phi);
    centers[i * 3 + 2] = depth + r * u;
    // Rayon tel que le volume total occupé reste comparable d'un palier à l'autre.
    radii[i] = (radius * 0.9) / Math.cbrt(n);
    ids[i] = i + 1;
    colors[i * 3] = 0.7 - 0.4 * (i / n);
    colors[i * 3 + 1] = 0.2 + 0.3 * (i / n);
    colors[i * 3 + 2] = 0.25 + 0.5 * (i / n);
  }
  return { centers, radii, ids, colors, count: n };
}

const percentile = (sorted: number[], q: number) =>
  sorted[Math.min(sorted.length - 1, Math.floor(sorted.length * q))];

async function main(): Promise<void> {
  const canvas = document.createElement("canvas");
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const W = Math.round(1280 * dpr);
  const H = Math.round(720 * dpr);
  canvas.width = W;
  canvas.height = H;
  canvas.style.width = "640px";
  canvas.style.height = "360px";
  document.body.insertBefore(canvas, out);

  const gl = canvas.getContext("webgl2", { antialias: false, powerPreference: "high-performance" });
  if (!gl) {
    say("WebGL2 indisponible.");
    return;
  }

  const dbg = gl.getExtension("WEBGL_debug_renderer_info");
  const renderer = dbg ? (gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) as string) : "inconnu";
  say("  Spike de rendu Geno — imposteurs de sphères, WebGL2");
  say(`  ${renderer}`);
  say(`  cible ${W}×${H} (dpr ${dpr})   ${FRAMES} frames mesurées par palier`);
  say("");
  say("  instances      niveau            frame médiane      p95        fps médian");
  say("  ------------   ---------------   -------------   ---------   ----------");

  const targets = createTargets(gl, W, H);
  const picker = new Picker(gl, targets);
  const proj = perspective(FOV, W / H, 0.1, 100);
  let pickMs: number[] = [];

  for (const step of STEPS) {
    const r = new ImpostorRenderer(gl);
    r.upload(nucleus(step.n, 3.2, -9));

    const times: number[] = [];
    for (let f = 0; f < FRAMES + WARMUP; f++) {
      const t0 = performance.now();
      gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
      gl.viewport(0, 0, W, H);
      gl.clearBufferfv(gl.COLOR, 0, [0.055, 0.07, 0.095, 1]);
      gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
      gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
      gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
      r.draw(proj);
      // Une lecture d'un pixel force la fin du GPU : sans ça on mesure la file
      // d'attente, pas le rendu.
      picker.pick(W >> 1, H >> 1);
      const dt = performance.now() - t0;
      if (f >= WARMUP) times.push(dt);
      await new Promise((res) => requestAnimationFrame(() => res(null)));
    }
    times.sort((a, b) => a - b);
    const med = times[times.length >> 1];
    say(
      `  ${step.n.toLocaleString("fr").padStart(12)}   ${step.label.padEnd(15)}   ` +
        `${med.toFixed(2).padStart(10)} ms   ${percentile(times, 0.95).toFixed(2).padStart(6)} ms   ` +
        `${(1000 / med).toFixed(0).padStart(8)}`,
    );

    if (step.n === 620_000) {
      pickMs = [];
      for (let k = 0; k < 60; k++) {
        const t0 = performance.now();
        picker.pick(Math.random() * W, Math.random() * H);
        pickMs.push(performance.now() - t0);
      }
      pickMs.sort((a, b) => a - b);
    }
    r.dispose();
  }

  say("");
  if (pickMs.length) {
    say(
      `  picking GPU à 620 000 instances : médiane ${pickMs[pickMs.length >> 1].toFixed(3)} ms, ` +
        `p95 ${percentile(pickMs, 0.95).toFixed(3)} ms`,
    );
    say("  (coût constant : le rasteriseur a déjà fait le travail, on relit un pixel)");
  }
  say("");
  say("  Verdict : un palier tient à 60 fps sous 16,7 ms, à 30 fps sous 33,3 ms.");
  say("  Reporter ces chiffres dans docs/ARCHITECTURE.md § 4 et docs/VALIDATION.md.");
}

main().catch((e) => say(`\n  ÉCHEC : ${e}`));
