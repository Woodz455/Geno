/** La plongée scriptée, mesurée. Même code en headless et sur matériel réel.
 *
 *  Deux passes sur le même trajet, à pas de paramètre constants (donc à vitesse perçue
 *  constante, `camera.ts`) :
 *
 *  1. **image par image, synchronisée** : un pixel relu après chaque rendu force la fin du
 *     travail GPU, si bien que le temps mesuré est celui du rendu, pas de la file d'attente.
 *     Si `continuity`, l'image entière est aussi relue (hors chronomètre) et comparée à la
 *     précédente : c'est la mesure des sauts.
 *  2. **au rythme de requestAnimationFrame**, sans aucune relecture : ce que verrait un
 *     utilisateur. Sur du vrai matériel, c'est elle qui dit si 60 images par seconde tiennent.
 */

import { dive } from "@geno/viewer";
import { type DiveScene, FOV } from "./dive.js";

export interface BenchResult {
  frames: number;
  size: [number, number];
  gpu: string;
  pathLength: number;
  instances: { preview: number; noyau: number };
  /** Écart moyen absolu entre images successives, RGB 0–255 (vide sans `continuity`). */
  diffs: number[];
  frameMs: number[];
  rafMs: number[];
  states: { w: number; lod: number; focus: number; mix: number }[];
}

export interface BenchOptions {
  continuity?: boolean;
  onProgress?: (fraction: number, phase: "synchronisée" | "rAF") => void;
}

export async function runBench(scene: DiveScene, frames: number, opts: BenchOptions = {}): Promise<BenchResult> {
  const gl = scene.gl;
  const [W, H] = scene.size();
  const path = dive(scene.start, scene.final, FOV);
  const one = new Uint8Array(4);
  const px = opts.continuity ? new Uint8Array(W * H * 4) : null;
  const prev = opts.continuity ? new Uint8Array(W * H * 4) : null;
  const diffs: number[] = [];
  const frameMs: number[] = [];
  const states: BenchResult["states"] = [];
  scene.arrival = 1;            // le fondu d'arrivée réseau ne dépend pas de la caméra : neutralisé

  for (let k = 0; k < 5; k++) {
    scene.render(path.at(0));
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, one);
  }
  for (let k = 0; k <= frames; k++) {
    const t0 = performance.now();
    scene.render(path.at((k / frames) * path.length));
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, one);
    frameMs.push(performance.now() - t0);
    states.push({ w: scene.last.widthUm, lod: scene.last.lod, focus: scene.last.focus, mix: scene.last.mix });
    if (px && prev) {
      gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, px);
      if (k > 0) {
        let s = 0;
        for (let i = 0; i < px.length; i += 4) {
          s += Math.abs(px[i] - prev[i]) + Math.abs(px[i + 1] - prev[i + 1]) + Math.abs(px[i + 2] - prev[i + 2]);
        }
        diffs.push(s / ((px.length / 4) * 3));
      }
      prev.set(px);
    }
    if (k % 10 === 0) {
      opts.onProgress?.(k / frames / 2, "synchronisée");
      await new Promise((r) => setTimeout(r, 0));     // laisser la page respirer
    }
  }

  const rafMs: number[] = [];
  await new Promise<void>((done) => {
    let k = 0;
    let last = 0;
    const step = (now: number) => {
      if (last) rafMs.push(now - last);
      last = now;
      scene.render(path.at((k / frames) * path.length));
      if (k % 10 === 0) opts.onProgress?.(0.5 + k / frames / 2, "rAF");
      if (++k <= frames) requestAnimationFrame(step);
      else done();
    };
    requestAnimationFrame(step);
  });

  return {
    frames,
    size: [W, H],
    gpu: scene.gpu,
    pathLength: path.length,
    instances: { preview: scene.preview.frame.count, noyau: scene.fine?.frame.count ?? 0 },
    diffs,
    frameMs,
    rafMs,
    states,
  };
}

export const median = (xs: number[]): number => {
  const s = xs.slice().sort((a, b) => a - b);
  return s.length ? s[s.length >> 1] : NaN;
};

export const percentile = (xs: number[], q: number): number => {
  const s = xs.slice().sort((a, b) => a - b);
  return s.length ? s[Math.min(s.length - 1, Math.floor(s.length * q))] : NaN;
};
