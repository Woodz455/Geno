/** Page de mesure du premier rendu.
 *
 *  Elle fait ce que fera l'application, et rien d'autre : ouvrir un `.g3d` par requête
 *  Range, décoder le premier niveau, dessiner, et **vérifier** que le dessin a eu lieu
 *  avant de publier son temps. Un chronomètre arrêté sur un canevas noir mesurerait le
 *  chargement d'un script.
 *
 *  Tous les temps sont comptés depuis l'origine de navigation (`performance.timeOrigin`),
 *  c'est-à-dire **avant** la requête du document : ce que voit l'utilisateur qui vient
 *  de toucher le lien.
 */

import { G3D, HttpSource } from "@geno/g3d";
import { createTargets, ImpostorRenderer, NO_HIT, perspective, Picker } from "@geno/viewer";

const q = new URLSearchParams(location.search);
const file = q.get("f") ?? "data/gm12878.g3d";
const prefix = Number(q.get("prefix") ?? 64 * 1024);
const W = 960;
const H = 540;
const FOV = (45 * Math.PI) / 180;

declare global {
  interface Window {
    __geno?: Record<string, unknown>;
  }
}

function publish(r: Record<string, unknown>): void {
  window.__geno = r;
  const pre = document.createElement("pre");
  pre.id = "result";
  pre.textContent = JSON.stringify(r, null, 1);
  document.body.appendChild(pre);
}

async function main(): Promise<void> {
  const tScript = performance.now();
  const source = new HttpSource(file);
  const g3d = await G3D.open(source, { prefix });
  const tOpen = performance.now();
  const frame = await g3d.firstFrame();
  const tData = performance.now();

  // Caméra : le niveau entier dans le champ, centré, vu depuis l'axe z.
  const n = frame.count;
  const p = frame.positions;
  let cx = 0, cy = 0, cz = 0;
  for (let i = 0; i < n; i++) {
    cx += p[3 * i]; cy += p[3 * i + 1]; cz += p[3 * i + 2];
  }
  cx /= n; cy /= n; cz /= n;
  let reach = 0;
  for (let i = 0; i < n; i++) {
    reach = Math.max(reach, Math.hypot(p[3 * i] - cx, p[3 * i + 1] - cy, p[3 * i + 2] - cz));
  }
  const lv = frame.level;
  const rmax = Math.max(...lv.radius_nm);
  const dist = ((reach + rmax) / Math.sin(FOV / 2)) * 1.05;

  const centers = new Float32Array(n * 3);
  const radii = new Float32Array(n);
  const ids = new Uint32Array(n);
  const colors = new Float32Array(n * 3);
  const k = lv.copies.length;
  for (let i = 0; i < n; i++) {
    centers[3 * i] = (p[3 * i] - cx) / 1000;
    centers[3 * i + 1] = (p[3 * i + 1] - cy) / 1000;
    centers[3 * i + 2] = (p[3 * i + 2] - cz - dist) / 1000;
    const c = frame.copy[i];
    radii[i] = lv.radius_nm[c] / 1000;
    ids[i] = i + 1;
    const h = (c / k) * 6.283;
    colors[3 * i] = 0.55 + 0.4 * Math.cos(h);
    colors[3 * i + 1] = 0.55 + 0.4 * Math.cos(h - 2.094);
    colors[3 * i + 2] = 0.55 + 0.4 * Math.cos(h + 2.094);
  }

  const canvas = document.createElement("canvas");
  canvas.width = W;
  canvas.height = H;
  document.body.appendChild(canvas);
  const gl = canvas.getContext("webgl2", { antialias: false });
  if (!gl) throw new Error("WebGL2 indisponible");
  const targets = createTargets(gl, W, H);
  const renderer = new ImpostorRenderer(gl);
  renderer.upload({ centers, radii, ids, colors, count: n });

  gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
  gl.viewport(0, 0, W, H);
  gl.clearBufferfv(gl.COLOR, 0, [0.055, 0.07, 0.095, 1]);
  gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
  gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
  gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
  renderer.draw(perspective(FOV, W / H, dist / 1000 / 50, (dist + reach + rmax) / 1000 * 2));

  // Copie vers le canevas visible : c'est ce que l'utilisateur voit.
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
  gl.readBuffer(gl.COLOR_ATTACHMENT0);
  gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);
  gl.blitFramebuffer(0, 0, W, H, 0, 0, W, H, gl.COLOR_BUFFER_BIT, gl.NEAREST);

  // Vérification, qui force aussi la fin du travail GPU : on compte les pixels qui
  // portent une bille, et on relit l'identifiant au centre de l'image.
  const idbuf = new Uint32Array(W * H);
  gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
  gl.readBuffer(gl.COLOR_ATTACHMENT1);
  gl.readPixels(0, 0, W, H, gl.RED_INTEGER, gl.UNSIGNED_INT, idbuf);
  const tDraw = performance.now();

  let covered = 0;
  const seen = new Set<number>();
  for (const v of idbuf) {
    if (v !== NO_HIT) {
      covered++;
      seen.add(v);
    }
  }
  const centre = new Picker(gl, targets).pick(W / 2, H / 2);
  const nav = performance.getEntriesByType("navigation")[0] as PerformanceNavigationTiming | undefined;

  publish({
    ok: covered > 0.05 * W * H && centre !== NO_HIT && centre <= n,
    level: lv.name,
    beads: n,
    beads_visible: seen.size,
    coverage: covered / (W * H),
    centre_id: centre,
    prefix,
    requests: source.log.length,
    g3d_bytes: source.log.reduce((s, r) => s + r.received, 0),
    first_frame_end: g3d.firstFrameEnd,
    t_document: nav ? nav.responseEnd : null,
    t_script: tScript,
    t_open: tOpen,
    t_data: tData,
    t_first_render: tDraw,
    log: source.log,
  });
}

main().catch((e) => publish({ ok: false, error: String(e) }));
