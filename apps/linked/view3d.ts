/** La 3D des vues liées : la scène de la plongée (semaines 10 à 13), dans un panneau.
 *
 *  Un clic publie ce qu'il désigne — l'intervalle de la bille ou du monomère, la base pour
 *  l'hélice — et ouvre sa fiche ; une sélection reçue est surlignée (`DiveScene.highlight`).
 *  Le rendu est à la demande, comme dans `dive.html` : une image en vol ferait attendre la
 *  lecture du pixel du clic suivant.
 */

import { FicheReader, HttpSource } from "@geno/g3d";
import { distanceFor, dive, type Orbit, rotate, zoom } from "@geno/viewer";
import { divePath } from "../dive/bench.js";
import { type DiveScene, FOV, type Highlighted, HIGHLIGHT, openDive, type Selection as Pick } from "../dive/dive.js";
import { compose, FichePanel } from "../dive/fiche.js";
import { IDS } from "../dive/fiber.js";
import type { Interval, Selection } from "./store.js";

export interface View3D {
  scene: DiveScene;
  canvas: HTMLCanvasElement;
  /** Ce que la dernière sélection a surligné, quand l'image qui le montre a été peinte, et le
   *  nombre de sélections **reçues** du store (pas celles que la 3D a publiées). */
  shown: { seq: number; hl: Highlighted | null; paintedAt: number; applied: number; navigated: number };
  apply(sel: Selection | null): void;
  /** Pour la mesure : la vue va à une étape de la plongée, rendue tout de suite. */
  goto(o: Orbit): void;
  /** Billes surlignées visibles : identifiants distincts et pixels, relus sur le dernier rendu,
   *  avec la part de ces pixels dont la couleur est bien celle de la surbrillance. */
  visible(): { ids: number; pixels: number; yellow: number };
  /** Pixels CSS (dans la page) qui portent un objet, tirés à graine fixe. */
  hits(n: number, seed: number): [number, number][];
  idle(): boolean;
  /** Pour la mesure : l'intervalle qu'un clic en ce point de la page publierait, sans cliquer. */
  peek(x: number, y: number): string | null;
  /** Dernier clic : la sélection publiée, le temps du pixel. */
  last: { pick: Pick | null; at: number };
}

export async function openView3D(root: HTMLElement, file: string, ficheFile: string,
                                 publish: (regions: Interval[], label: string) => void): Promise<View3D> {
  const canvas = document.createElement("canvas");
  canvas.style.cssText = "display:block;width:100%;height:100%;touch-action:none";
  root.appendChild(canvas);
  const overlay = document.createElement("canvas");
  overlay.style.cssText = "position:absolute;left:0;top:0;width:100%;height:100%;pointer-events:none";
  root.appendChild(overlay);
  const fit = (): [number, number] => {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, Math.round(root.clientWidth * dpr));
    const h = Math.max(1, Math.round(root.clientHeight * dpr));
    canvas.width = overlay.width = w;
    canvas.height = overlay.height = h;
    return [w, h];
  };
  fit();
  const scene = await openDive(canvas, file, {});
  scene.arrival = 1;

  let orbit: Orbit = { ...scene.overview };
  let anim: { path: { length: number; at(s: number): Orbit }; t0: number; duration: number } | null = null;
  let dirty = true;
  let marker: number[] | null = null;
  const shown: View3D["shown"] = { seq: 0, hl: null, paintedAt: 0, applied: 0, navigated: 0 };
  let pendingSeq = -1;
  const last: View3D["last"] = { pick: null, at: 0 };

  const panel = new FichePanel(document);
  let reader: FicheReader | null = null;
  void FicheReader.open(new HttpSource(ficheFile)).then((r) => (reader = r)).catch(() => null);

  const drawMarker = () => {
    const ctx = overlay.getContext("2d")!;
    ctx.clearRect(0, 0, overlay.width, overlay.height);
    if (!marker) return;
    const [x, y] = scene.project(marker);
    const r = 9 * (overlay.width / overlay.clientWidth || 1);
    ctx.strokeStyle = `rgb(${HIGHLIGHT.map((c) => Math.round(c * 255)).join(",")})`;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(x, y, r, 0, 2 * Math.PI);
    ctx.stroke();
  };

  const onClick = (x: number, y: number) => {
    const sel = scene.select(x, y);
    last.pick = sel;
    last.at = performance.now();
    if (!sel) return;
    marker = sel.world;
    drawMarker();
    const regions = [sel.pos >= 0 ? { chrom: sel.chrom, start: sel.pos, end: sel.pos + 1 } : { chrom: sel.chrom, start: sel.start, end: sel.end }];
    // La 3D surligne sa propre sélection : elle la reçoit du store comme les autres vues ne
    // la reçoivent pas — le store ne la renvoie pas à sa source.
    applyLocal(regions, -1);
    publish(regions, sel.what);
    const r = reader;
    if (r) {
      const card = r.card(sel.level, sel.bead);
      panel.show(compose(sel, card, r.header, r, scene.sequence?.reader ?? null), r.header, !card);
      if (!card) void r.cardAsync(sel.level, sel.bead).then((full) => panel.show(compose(sel, full, r.header, r, scene.sequence?.reader ?? null), r.header));
    }
  };

  /** Si rien de la sélection n'est à l'écran, l'orbite qui la montre ; sinon null. À l'échelle
   *  de la fibre, une sélection qui ne la touche pas fait remonter au noyau ; ailleurs, la
   *  caméra garde son orientation et vise le centre des billes surlignées, avec un champ qui
   *  les contient toutes (les deux homologues peuvent être loin l'un de l'autre). */
  const frameSelection = (): Orbit | null => {
    const hl = shown.hl;
    const F = scene.fiber;
    if (!hl) return null;
    const [W, H] = scene.size();
    const onScreen = (p: number[]) => {
      const [x, y] = scene.project(p);
      return Number.isFinite(x) && x >= 0 && x <= W && y >= 0 && y <= H;
    };
    const fibre = scene.last.fiber > 0.5;
    if (fibre && hl.fiber && F) {
      const [lo, hi] = hl.fiber;
      const pts = Array.from({ length: hi - lo + 1 }, (_, k) => [F.points[3 * (lo + k)], F.points[3 * (lo + k) + 1], F.points[3 * (lo + k) + 2]]);
      if (pts.some(onScreen)) return null;
      return { ...orbit, target: pts[pts.length >> 1] as [number, number, number] };
    }
    const lay = scene.fine ?? scene.preview;
    const ids = scene.fine ? hl.nucleus : hl.preview;
    if (!ids.length) return null;
    const P = lay.frame.positions;
    const pts = ids.map((i) => [P[3 * (i - 1)] / 1000, P[3 * (i - 1) + 1] / 1000, P[3 * (i - 1) + 2] / 1000]);
    if (!fibre && pts.some(onScreen)) return null;
    const c = [0, 1, 2].map((k) => pts.reduce((s, p) => s + p[k], 0) / pts.length) as [number, number, number];
    const spread = Math.max(...pts.map((p) => Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2])));
    const field = Math.min(24, Math.max(3, 2.4 * spread + 2));
    return { ...orbit, target: field >= 24 ? [0, 0, 0] : c, distance: distanceFor(field, FOV) };
  };

  const applyLocal = (regions: Interval[] | null, seq: number) => {
    shown.hl = scene.highlight(regions);
    pendingSeq = seq;
    dirty = true;
  };

  let drag: { x: number; y: number; moved: number } | null = null;
  canvas.addEventListener("pointerdown", (e) => {
    drag = { x: e.clientX, y: e.clientY, moved: 0 };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x;
    const dy = e.clientY - drag.y;
    drag = { x: e.clientX, y: e.clientY, moved: drag.moved + Math.abs(dx) + Math.abs(dy) };
    orbit = rotate(orbit, -dx * 0.006, dy * 0.006);
    anim = null;
    dirty = true;
  });
  canvas.addEventListener("pointerup", (e) => {
    const wasClick = drag && drag.moved < 4;
    drag = null;
    if (!wasClick) return;
    const dpr = canvas.width / canvas.clientWidth;
    onClick(e.offsetX * dpr, e.offsetY * dpr);
  });
  canvas.addEventListener("wheel", (e) => {
    e.preventDefault();
    orbit = zoom(orbit, Math.exp(e.deltaY * 0.0015), 0.003, 150);
    anim = null;
    dirty = true;
  }, { passive: false });
  new ResizeObserver(() => {
    scene.resize(...fit());
    dirty = true;
  }).observe(root);

  // Deux boutons, comme la plongée : descendre vers la base visée, remonter.
  const bar = document.createElement("div");
  bar.className = "bar";
  const button = (label: string, to: () => { length: number; at(s: number): Orbit }) => {
    const b = document.createElement("button");
    b.textContent = label;
    b.onclick = () => {
      const path = to();
      anim = { path, t0: performance.now(), duration: 900 * Math.max(0.5, path.length) };
    };
    bar.appendChild(b);
  };
  button(`Plonger vers ${scene.target?.name ?? "la base"}`, () => divePath({ ...scene, stops: [orbit, ...scene.stops.slice(1)] }));
  button("Remonter", () => dive(orbit, scene.overview, FOV));
  root.appendChild(bar);

  const loop = (now: number) => {
    if (anim) {
      const k = Math.min(1, (now - anim.t0) / anim.duration);
      orbit = anim.path.at(k * anim.path.length);
      if (k >= 1) anim = null;
      dirty = true;
    }
    if (scene.busy) dirty = true;
    if (dirty) {
      scene.render(orbit);
      drawMarker();
      dirty = false;
      if (pendingSeq !== -2) {
        shown.seq = pendingSeq;
        shown.paintedAt = performance.now();
        pendingSeq = -2;
      }
    }
    requestAnimationFrame(loop);
  };
  pendingSeq = -2;
  requestAnimationFrame(loop);

  return {
    scene,
    canvas,
    shown,
    last,
    apply(sel) {
      // La fiche ouverte est celle du dernier clic dans la 3D : une sélection venue d'ailleurs
      // la rend caduque, et elle couvrirait ce qui est surligné.
      panel.hide();
      marker = null;
      shown.applied++;
      applyLocal(sel ? sel.regions : null, sel ? sel.seq : 0);
      const to = sel ? frameSelection() : null;
      if (to) {
        shown.navigated++;
        const path = dive(orbit, to, FOV);
        anim = { path, t0: performance.now(), duration: 600 * Math.max(0.5, path.length) };
      }
    },
    goto(o) {
      orbit = o;
      anim = null;
      scene.render(o);
      scene.render(o);
      dirty = true;
    },
    visible() {
      // Le rendu courant, relu dans la même tâche : l'image composée est encore là.
      scene.render(orbit);
      const ids = scene.ids();
      const gl = scene.gl;
      const [W, H] = scene.size();
      const rgb = new Uint8Array(W * H * 4);
      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
      gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, rgb);
      const hl = shown.hl;
      const want = new Set<number>();
      if (hl) {
        for (const i of hl.preview) want.add(1_000_000 + i);
        for (const i of hl.nucleus) want.add(i);
        if (hl.fiber) for (let k = hl.fiber[0]; k <= hl.fiber[1]; k++) want.add(IDS.fiber + k + 1);
      }
      const seen = new Set<number>();
      let pixels = 0;
      let yellow = 0;
      for (let i = 0; i < ids.length; i++) {
        if (!want.has(ids[i])) continue;
        seen.add(ids[i]);
        // Jaune : jugé sur la teinte, que l'ombrage et l'occlusion ne changent pas ; un pixel
        // trop sombre pour en avoir une n'est pas compté.
        const r = rgb[4 * i], g = rgb[4 * i + 1], b = rgb[4 * i + 2];
        const mx = Math.max(r, g, b), mn = Math.min(r, g, b);
        if (mx < 30) continue;
        pixels++;
        const hue = mx === mn ? -1 : mx === r ? (60 * (g - b)) / (mx - mn) : mx === g ? 120 + (60 * (b - r)) / (mx - mn) : 240 + (60 * (r - g)) / (mx - mn);
        if (hue >= 35 && hue <= 75 && (mx - mn) / mx > 0.35) yellow++;
      }
      return { ids: seen.size, pixels, yellow };
    },
    hits(n, seed) {
      const ids = scene.ids();
      const [W, H] = scene.size();
      const dpr = canvas.width / canvas.clientWidth;
      const box = canvas.getBoundingClientRect();
      let st = seed >>> 0;
      const rnd = () => ((st = (st * 1664525 + 1013904223) >>> 0) / 2 ** 32);
      const out: [number, number][] = [];
      // Assez d'essais pour une scène dont une sélection a effacé les copies à 85 %.
      for (let tries = 0; out.length < n && tries < Math.max(4000, 200 * n); tries++) {
        const x = Math.floor(rnd() * W), y = Math.floor(rnd() * H);
        const cx = box.left + (x + 0.5) / dpr, cy = box.top + (y + 0.5) / dpr;
        if (ids[(H - 1 - y) * W + x] && document.elementFromPoint(cx, cy) === canvas) out.push([cx, cy]);
      }
      return out;
    },
    idle: () => !dirty && !anim && !scene.busy,
    peek(x, y) {
      const box = canvas.getBoundingClientRect();
      const dpr = canvas.width / canvas.clientWidth;
      const sel = scene.select((x - box.left) * dpr, (y - box.top) * dpr);
      if (!sel) return null;
      return sel.pos >= 0 ? `${sel.chrom}:${sel.pos}-${sel.pos + 1}` : `${sel.chrom}:${sel.start}-${sel.end}`;
    },
  };
}
