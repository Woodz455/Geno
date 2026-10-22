/** Le mode variant structural (semaine 15) : supprimer une frontière, voir l'enhancer changer de cible.
 *
 *  Le cas EPHA4 de Lupiáñez 2015, sur le locus WNT6/IHH/EPHA4/PAX3 de GM12878. À gauche, le
 *  locus en 3D (OpenMM) pour les allèles précalculés ; à droite, les pistes et deux 4C virtuels
 *  — depuis le promoteur de PAX3, comme l'article, et depuis la grappe d'enhancers —, recalculés
 *  **en direct** par l'extrusion 1D pour n'importe quelle délétion tracée sur la piste.
 *
 *  Paramètres d'URL : f=… (défaut data/epha4.variant.json).
 */

import { createTargets, distanceFor, NO_HIT, type Orbit, perspective, rotate, Ssao, TubeRenderer, viewMatrix, zoom } from "@geno/viewer";
import { decodeCoords, type Deletion, predict, type Prediction, type VariantDoc } from "./model.js";

const q = new URLSearchParams(location.search);
const FILE = q.get("f") ?? "data/epha4.variant.json";
const FOV = (40 * Math.PI) / 180;

declare global {
  interface Window {
    __geno?: Record<string, unknown>;
    __variant?: unknown;
  }
}

// Okabe–Ito.
const C = {
  cen: [0.0, 0.62, 0.45],        // TAD de WNT6/IHH (vert bleuté)
  tad: [0.9, 0.62, 0.0],         // TAD d'EPHA4 (orange)
  tel: [0.0, 0.45, 0.7],         // TAD de PAX3 (bleu)
  enh: [0.94, 0.89, 0.26],       // enhancers (jaune)
  pax3: [0.84, 0.37, 0.0],       // promoteur de PAX3 (vermillon)
  epha4: [0.8, 0.47, 0.65],      // promoteur d'EPHA4 (pourpre)
} as const;
const css = (c: readonly number[], a = 1) => `rgba(${c.map((x) => Math.round(x * 255)).join(",")},${a})`;

const style = document.createElement("style");
style.textContent = `
  html, body { height: 100%; overflow: hidden; }
  body { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.15fr); grid-template-rows: auto 1fr; gap: 6px;
         padding: 6px; box-sizing: border-box; font: 13px/1.45 system-ui, sans-serif; }
  header { grid-column: 1 / 3; display: flex; flex-wrap: wrap; gap: 6px 10px; align-items: center; }
  header h1 { font-size: 14px; margin: 0 8px 0 0; font-weight: 600; }
  header button, header select { font: inherit; font-size: 12px; padding: 3px 10px; border-radius: 6px;
         border: 1px solid #3a4658; background: #1b2330; color: inherit; cursor: pointer; }
  header button[aria-pressed="true"] { border-color: rgb(240,228,66); color: rgb(240,228,66); }
  .k { color: #8b98aa; }
  #view { position: relative; border: 1px solid #2a3442; border-radius: 6px; overflow: hidden; min-height: 0; }
  #view canvas { position: absolute; inset: 0; width: 100%; height: 100%; }
  #view .cap { position: absolute; left: 8px; bottom: 6px; font-size: 11px; color: #c9d3df; background: rgba(14,18,24,.8);
               padding: 2px 8px; border-radius: 6px; max-width: calc(100% - 24px); }
  #side { display: grid; grid-template-rows: auto auto 1fr; gap: 6px; min-height: 0; }
  #tracks, #plots { border: 1px solid #2a3442; border-radius: 6px; position: relative; }
  #tracks canvas, #plots canvas { display: block; width: 100%; }
  #tracks canvas { cursor: crosshair; touch-action: none; }
  #readout { padding: 6px 8px; border: 1px solid #2a3442; border-radius: 6px; font-size: 12px; overflow: auto; }
  #readout b { color: rgb(240,228,66); font-weight: 600; }
  #readout .src { color: #8b98aa; font-size: 11px; margin-top: 4px; }
  @media (max-width: 760px) { body { grid-template-columns: 1fr; grid-template-rows: auto 45vh auto; overflow: auto; } header { grid-column: 1; } }
`;
document.head.appendChild(style);

function publish(r: Record<string, unknown>): void {
  window.__geno = r;
}

async function main(): Promise<void> {
  const doc = (await (await fetch(FILE)).json()) as VariantDoc;
  const R = doc.region;
  const header = document.createElement("header");
  header.innerHTML = `<h1>Geno — variant structural au locus EPHA4</h1>`;
  const view = document.createElement("section");
  view.id = "view";
  const side = document.createElement("section");
  side.id = "side";
  document.body.append(header, view, side);

  // — Les allèles précalculés et les préréglages —
  const dels = doc.variants;
  const central = doc.alleles.delb?.variant;
  const presets: { key: string; label: string; del: Deletion | null }[] = [
    { key: "wt", label: "Type sauvage", del: null },
  ];
  for (const [key, label] of [["delb", "DEL-B (brachydactylie)"], ["delbs", "DelBS (frontière gardée)"], ["delb-keep", "DEL-B, CTCF de la frontière épargnés"]] as const) {
    const a = doc.alleles[key];
    if (a?.variant) presets.push({ key, label, del: { ...a.variant, name: a.name, kind: a.kind ?? undefined } as Deletion });
  }
  const buttons = new Map<string, HTMLButtonElement>();
  for (const p of presets) {
    const b = document.createElement("button");
    b.textContent = p.label;
    b.onclick = () => choose(p.key, p.del);
    header.appendChild(b);
    buttons.set(p.key, b);
  }
  const layerSel = document.createElement("select");
  layerSel.innerHTML = `<option value="domains">couleurs : domaines</option><option value="clinvar">ClinVar P/LP</option>
    <option value="gwas">GWAS</option><option value="loeuf">contrainte gnomAD (LOEUF)</option>`;
  layerSel.setAttribute("aria-label", "Surimpression sur la structure");
  header.appendChild(layerSel);
  const hint = document.createElement("span");
  hint.className = "k";
  hint.textContent = "ou glisser sur la piste pour tracer une délétion";
  header.appendChild(hint);

  // — 3D —
  const canvas = document.createElement("canvas");
  const overlay = document.createElement("canvas");
  overlay.style.pointerEvents = "none";
  const cap = document.createElement("div");
  cap.className = "cap";
  view.append(canvas, overlay, cap);
  const gl = canvas.getContext("webgl2", { antialias: false })!;
  const fit = (): [number, number] => {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, Math.round(view.clientWidth * dpr)), h = Math.max(1, Math.round(view.clientHeight * dpr));
    canvas.width = overlay.width = w;
    canvas.height = overlay.height = h;
    return [w, h];
  };
  let [W, H] = fit();
  let targets = createTargets(gl, W, H);
  let ssao = new Ssao(gl, targets);
  const tube = new TubeRenderer(gl);

  const um = (doc.alleles.wt.sigma_nm ?? 46) / 1000;
  const wtXYZ = decodeCoords(doc.alleles.wt.coords);
  const n0 = R.n;
  const mid = (k: number) => R.start + (k + 0.5) * R.bp_per_bead;
  // Position, couleur et opacité courantes de chaque monomère du type sauvage.
  const pos = new Float32Array(n0 * 3);
  for (let i = 0; i < n0 * 3; i++) pos[i] = wtXYZ[i] * um;
  let from = pos.slice();
  let to = pos.slice();
  const alpha = new Float32Array(n0).fill(1);
  let alphaFrom = alpha.slice(), alphaTo = alpha.slice();
  let anim: { t0: number; ms: number } | null = null;
  let dirty = true;

  let cx = 0, cy = 0, cz = 0;
  for (let k = 0; k < n0; k++) { cx += pos[3 * k]; cy += pos[3 * k + 1]; cz += pos[3 * k + 2]; }
  cx /= n0; cy /= n0; cz /= n0;
  let ext = 0;
  for (let k = 0; k < n0; k++) ext = Math.max(ext, Math.hypot(pos[3 * k] - cx, pos[3 * k + 1] - cy, pos[3 * k + 2] - cz));
  let orbit: Orbit = { target: [cx, cy, cz], yaw: 0.4, pitch: 0.25, distance: distanceFor(2.4 * ext, FOV) };

  // Couleurs : domaines, ou une couche clinique.
  const inside = (m: number, a: number, b: number) => m >= a && m < b;
  const colors = new Float32Array(n0 * 3);
  const radii = new Float32Array(n0);
  const recolor = () => {
    const layer = layerSel.value;
    const cl = doc.clinical;
    const gloeuf = new Float32Array(n0).fill(NaN);
    if (layer === "loeuf") {
      for (const g of cl.genes) if (g.LOEUF !== null) {
        for (let k = 0; k < n0; k++) if (inside(mid(k), g.start, g.end)) gloeuf[k] = Math.min(Number.isNaN(gloeuf[k]) ? 9 : gloeuf[k], g.LOEUF);
      }
    }
    const maxOf = (a: number[]) => Math.max(1, ...a);
    const mPlp = maxOf(cl.clinvar_plp), mGw = maxOf(cl.gwas);
    for (let k = 0; k < n0; k++) {
      const m = mid(k);
      let c: readonly number[] =
        m < doc.tad[0] ? C.cen : m <= doc.boundary[1] ? C.tad : C.tel;
      if (layer !== "domains") {
        // Gris neutre pour l'absence, puis du clair au foncé (une seule teinte, magenta).
        const v = layer === "clinvar" ? Math.sqrt(cl.clinvar_plp[k] / mPlp) : layer === "gwas" ? Math.sqrt(cl.gwas[k] / mGw)
          : Number.isNaN(gloeuf[k]) ? 0 : Math.max(0, Math.min(1, 1 - gloeuf[k] / 1.0));
        c = v > 0 ? [0.95 - 0.35 * v, 0.75 - 0.6 * v, 0.95 - 0.35 * v] : [0.42, 0.45, 0.5];
      }
      radii[k] = 0.5 * um;
      if (inside(m, doc.enhancer[0], doc.enhancer[1])) { c = C.enh; radii[k] = 1.2 * um; }
      const pax3 = doc.genes.find((g) => g.name === "PAX3")!, epha4 = doc.genes.find((g) => g.name === "EPHA4")!;
      if (Math.abs(m - pax3.tss) < doc.promoter_half) { c = C.pax3; radii[k] = 1.2 * um; }
      if (Math.abs(m - epha4.tss) < doc.promoter_half) { c = C.epha4; radii[k] = 1.2 * um; }
      colors.set(c, 3 * k);
    }
    dirty = true;
  };
  recolor();
  layerSel.onchange = () => {
    recolor();
    drawTracks();
  };

  const upload = () => {
    // Le tube ne passe que par les monomères visibles : une délétion le coupe en deux, et les
    // deux bouts se rejoignent quand la conformation du mutant est là.
    const keep: number[] = [];
    for (let k = 0; k < n0; k++) if (alpha[k] > 0.02) keep.push(k);
    const pts = new Float32Array(keep.length * 3), rr = new Float32Array(keep.length), cc = new Float32Array(keep.length * 3);
    const ids = new Uint32Array(keep.length);
    const breaks: number[] = [];
    keep.forEach((k, i) => {
      pts.set(pos.subarray(3 * k, 3 * k + 3), 3 * i);
      rr[i] = radii[k];
      cc.set(colors.subarray(3 * k, 3 * k + 3), 3 * i);
      ids[i] = k + 1;
      if (i > 0 && k !== keep[i - 1] + 1 && !joined) breaks.push(i);
    });
    tube.upload({ points: pts, radii: rr, colors: cc, ids, breaks });
    const fade = new Float32Array(keep.length);
    keep.forEach((k, i) => (fade[i] = alpha[k]));
    tube.setFade(fade);
  };
  let joined = false;

  const render = () => {
    const view = viewMatrix(orbit);
    const proj = perspective(FOV, W / H, Math.max(1e-4, orbit.distance - 3 * ext), orbit.distance + 3 * ext);
    gl.bindFramebuffer(gl.FRAMEBUFFER, targets.fbo);
    gl.viewport(0, 0, W, H);
    gl.clearBufferfv(gl.COLOR, 0, [0.055, 0.07, 0.095, 1]);
    gl.clearBufferuiv(gl.COLOR, 1, [NO_HIT, 0, 0, 0]);
    gl.clearBufferfv(gl.COLOR, 2, [1, 0, 0, 1]);
    gl.clearBufferfi(gl.DEPTH_STENCIL, 0, 1.0, 0);
    tube.draw(proj, { view, height: H, idBase: 0 });
    ssao.occlusion(proj, { radius: 2 * um });
    ssao.composite(proj, { radius: 2 * um, strength: 0.8 });
    drawLabels(view, proj);
  };

  const project = (view: Float32Array, proj: Float32Array, p: number[]): [number, number, number] => {
    const e = [0, 1, 2].map((k) => view[k] * p[0] + view[4 + k] * p[1] + view[8 + k] * p[2] + view[12 + k]);
    const x = proj[0] * e[0] + proj[4] * e[1] + proj[8] * e[2] + proj[12];
    const y = proj[1] * e[0] + proj[5] * e[1] + proj[9] * e[2] + proj[13];
    const w = proj[3] * e[0] + proj[7] * e[1] + proj[11] * e[2] + proj[15];
    return [((x / w) * 0.5 + 0.5) * W, (1 - ((y / w) * 0.5 + 0.5)) * H, w];
  };
  const centerOf = (a: number, b: number): number[] | null => {
    let s = [0, 0, 0], c = 0;
    for (let k = 0; k < n0; k++) if (alpha[k] > 0.5 && mid(k) >= a && mid(k) < b) { s = [s[0] + pos[3 * k], s[1] + pos[3 * k + 1], s[2] + pos[3 * k + 2]]; c++; }
    return c ? s.map((v) => v / c) : null;
  };
  const drawLabels = (view: Float32Array, proj: Float32Array) => {
    const ctx = overlay.getContext("2d")!;
    ctx.clearRect(0, 0, W, H);
    const k = W / Math.max(1, overlay.clientWidth);
    const tag = (text: string, c: readonly number[], a: number, b: number) => {
      const p = centerOf(a, b);
      if (!p) return;
      const [x, y, w] = project(view, proj, p);
      if (w <= 0) return;
      ctx.font = `600 ${Math.round(12 * k)}px system-ui, sans-serif`;
      const tw = ctx.measureText(text).width;
      ctx.fillStyle = "rgba(14,18,24,0.82)";
      ctx.fillRect(x + 8 * k, y - 16 * k, tw + 10 * k, 18 * k);
      ctx.fillStyle = css(c);
      ctx.fillText(text, x + 13 * k, y - 3 * k);
      ctx.strokeStyle = css(c);
      ctx.lineWidth = 1.5 * k;
      ctx.beginPath();
      ctx.moveTo(x, y);
      ctx.lineTo(x + 8 * k, y - 7 * k);
      ctx.stroke();
    };
    const g = (name: string) => doc.genes.find((x) => x.name === name)!;
    tag("enhancers de membre", C.enh, doc.enhancer[0], doc.enhancer[1]);
    tag("PAX3", C.pax3, g("PAX3").tss - doc.promoter_half, g("PAX3").tss + doc.promoter_half);
    tag("EPHA4", C.epha4, g("EPHA4").tss - doc.promoter_half, g("EPHA4").tss + doc.promoter_half);
  };

  let drag: { x: number; y: number } | null = null;
  canvas.addEventListener("pointerdown", (e) => { drag = { x: e.clientX, y: e.clientY }; canvas.setPointerCapture(e.pointerId); });
  canvas.addEventListener("pointermove", (e) => {
    if (!drag) return;
    orbit = rotate(orbit, -(e.clientX - drag.x) * 0.006, (e.clientY - drag.y) * 0.006);
    drag = { x: e.clientX, y: e.clientY };
    dirty = true;
  });
  canvas.addEventListener("pointerup", () => (drag = null));
  canvas.addEventListener("wheel", (e) => { e.preventDefault(); orbit = zoom(orbit, Math.exp(e.deltaY * 0.0015), 0.05 * ext, 20 * ext); dirty = true; }, { passive: false });
  new ResizeObserver(() => {
    [W, H] = fit();
    ssao.dispose();
    targets = createTargets(gl, W, H);
    ssao = new Ssao(gl, targets);
    dirty = true;
    drawTracks();
    drawPlots();
  }).observe(view);

  // — Pistes, 4C, lecture —
  const tracksBox = document.createElement("div");
  tracksBox.id = "tracks";
  const tracks = document.createElement("canvas");
  tracks.setAttribute("role", "img");
  tracks.setAttribute("aria-label", "Pistes du locus : gènes, CTCF orientés, TAD, enhancers ; glisser pour tracer une délétion");
  tracksBox.appendChild(tracks);
  const plotsBox = document.createElement("div");
  plotsBox.id = "plots";
  const plots = document.createElement("canvas");
  plots.setAttribute("role", "img");
  plots.setAttribute("aria-label", "4C virtuels depuis PAX3 et depuis les enhancers, type sauvage et allèle courant");
  plotsBox.appendChild(plots);
  const readout = document.createElement("div");
  readout.id = "readout";
  readout.setAttribute("aria-live", "polite");
  side.append(tracksBox, plotsBox, readout);

  const TR_H = 150, PL_H = 250;
  const sizeCanvas = (c: HTMLCanvasElement, h: number) => {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    c.width = Math.round(c.parentElement!.clientWidth * dpr);
    c.height = Math.round(h * dpr);
    c.style.height = `${h}px`;
    return dpr;
  };
  const xOf = (g: number, w: number) => ((g - R.start) / (R.end - R.start)) * w;
  const gOf = (x: number, w: number) => R.start + (x / w) * (R.end - R.start);

  let wt: Prediction | null = null;
  let cur: Prediction | null = null;
  let curDel: Deletion | null = null;
  let draft: [number, number] | null = null;
  let current = "wt";

  const drawTracks = () => {
    const dpr = sizeCanvas(tracks, TR_H);
    const w = tracks.width, ctx = tracks.getContext("2d")!;
    ctx.clearRect(0, 0, w, tracks.height);
    ctx.save();
    ctx.scale(dpr, dpr);
    const W1 = w / dpr;
    const X = (g: number) => xOf(g, W1);
    // Domaines : la boucle de Rao 2014 qui borne le TAD d'EPHA4.
    ctx.fillStyle = css(C.cen, 0.18); ctx.fillRect(0, 4, X(doc.tad[0]), 14);
    ctx.fillStyle = css(C.tad, 0.25); ctx.fillRect(X(doc.tad[0]), 4, X(doc.boundary[1]) - X(doc.tad[0]), 14);
    ctx.fillStyle = css(C.tel, 0.2); ctx.fillRect(X(doc.boundary[1]), 4, W1 - X(doc.boundary[1]), 14);
    ctx.fillStyle = "#c9d3df";
    ctx.font = "11px system-ui, sans-serif";
    ctx.fillText("TAD d'EPHA4 (Rao 2014, GM12878)", X(doc.tad[0]) + 4, 15);
    // Échelle.
    ctx.strokeStyle = "#3a4658"; ctx.fillStyle = "#8b98aa";
    for (let m = Math.ceil(R.start / 5e5) * 5e5; m < R.end; m += 5e5) {
      ctx.beginPath(); ctx.moveTo(X(m), 20); ctx.lineTo(X(m), 26); ctx.stroke();
      ctx.fillText(`${(m / 1e6).toFixed(1)}`, X(m) + 2, 31);
    }
    ctx.fillText("Mb, chr2 (hg38)", W1 - 86, 31);
    // Gènes.
    for (const g of doc.genes) {
      const x0 = X(g.start), x1 = Math.max(X(g.end), x0 + 2);
      const c = g.name === "PAX3" ? C.pax3 : g.name === "EPHA4" ? C.epha4 : [0.6, 0.66, 0.74];
      ctx.fillStyle = css(c);
      ctx.fillRect(x0, 46, x1 - x0, 8);
      ctx.fillText(`${g.name} ${g.strand === "+" ? "›" : "‹"}`, x0, 68);
    }
    // CTCF orientés : flèche vers la droite (+) ou la gauche (−), hauteur ∝ occupation.
    for (const p of doc.ctcf) {
      if (!p.motif || p.occupancy === null) continue;
      const x = X((p.start + p.end) / 2), h = 4 + 14 * p.occupancy;
      ctx.fillStyle = p.motif === "+" ? "#56B4E9" : "#CC79A7";
      ctx.beginPath();
      if (p.motif === "+") { ctx.moveTo(x - 3, 92); ctx.lineTo(x + 3, 92 - h / 2); ctx.lineTo(x - 3, 92 - h); }
      else { ctx.moveTo(x + 3, 92); ctx.lineTo(x - 3, 92 - h / 2); ctx.lineTo(x + 3, 92 - h); }
      ctx.fill();
    }
    ctx.fillStyle = "#8b98aa";
    ctx.fillText("CTCF GM12878 (› +, ‹ −)", 4, 104);
    // Enhancers.
    ctx.fillStyle = css(C.enh);
    ctx.fillRect(X(doc.enhancer[0]) - 1, 110, Math.max(3, X(doc.enhancer[1]) - X(doc.enhancer[0])), 10);
    ctx.fillText("enhancers (reconstruits)", X(doc.enhancer[0]) + 6, 119);
    // Couche clinique.
    if (layerSel.value === "clinvar" || layerSel.value === "gwas") {
      const v = layerSel.value === "clinvar" ? doc.clinical.clinvar_plp : doc.clinical.gwas;
      const mx = Math.max(1, ...v);
      ctx.fillStyle = "#c06bc0";
      v.forEach((c, k) => { if (c) { const x = X(mid(k)); const h = 2 + 12 * Math.sqrt(c / mx); ctx.fillRect(x, 146 - h, 1.5, h); } });
      ctx.fillStyle = "#8b98aa";
      ctx.fillText(layerSel.value === "clinvar" ? "ClinVar P/LP par 2 kb" : "associations GWAS par 2 kb", 4, 134);
    }
    // Délétion courante, ou en cours de tracé.
    const d = draft ?? (curDel ? [curDel.start, curDel.end] : null);
    if (d) {
      ctx.fillStyle = "rgba(255,255,255,0.12)";
      ctx.fillRect(X(d[0]), 0, X(d[1]) - X(d[0]), TR_H);
      ctx.strokeStyle = "rgba(255,255,255,0.8)";
      ctx.strokeRect(X(d[0]) + 0.5, 0.5, X(d[1]) - X(d[0]) - 1, TR_H - 1);
      for (const [a, b] of curDel?.keep ?? []) { ctx.fillStyle = css(C.enh, 0.9); ctx.fillRect(X(a) - 1, 0, Math.max(3, X(b) - X(a)), TR_H); }
    }
    ctx.restore();
  };

  const drawPlots = () => {
    const dpr = sizeCanvas(plots, PL_H);
    const w = plots.width, ctx = plots.getContext("2d")!;
    ctx.clearRect(0, 0, w, plots.height);
    ctx.save();
    ctx.scale(dpr, dpr);
    const W1 = w / dpr;
    const X = (g: number) => xOf(g, W1);
    const panel = (y0: number, h: number, title: string, key: "PAX3" | "enhancer") => {
      ctx.fillStyle = "#c9d3df";
      ctx.font = "12px system-ui, sans-serif";
      ctx.fillText(title, 4, y0 + 12);
      const lo = -4.5, hi = 0; // log10 de la probabilité de contact
      const Y = (p: number) => y0 + 18 + (h - 22) * (1 - (Math.log10(Math.max(p, 10 ** lo)) - lo) / (hi - lo));
      ctx.strokeStyle = "#2a3442";
      for (let e = lo; e <= hi; e += 1.5) { ctx.beginPath(); ctx.moveTo(0, Y(10 ** e)); ctx.lineTo(W1, Y(10 ** e)); ctx.stroke(); }
      const curve = (p: Prediction, col: string, lw: number) => {
        ctx.strokeStyle = col;
        ctx.lineWidth = lw;
        ctx.beginPath();
        let started = false;
        const vp = p.lm[key === "PAX3" ? "PAX3" : "enhancer"];
        for (let i = 0; i < p.allele.n; i++) {
          const x = X(p.allele.mids[i]);
          if (vp.includes(i)) { started = false; continue; }
          const gap = i > 0 && p.allele.wild[i] !== p.allele.wild[i - 1] + 1;
          if (!started || gap) { ctx.moveTo(x, Y(p.from[key][i])); started = true; } else ctx.lineTo(x, Y(p.from[key][i]));
        }
        ctx.stroke();
      };
      if (wt) curve(wt, "rgba(200,208,220,0.55)", 1);
      if (cur && current !== "wt") curve(cur, "rgb(240,228,66)", 1.5);
      // Repères.
      const mark = (g: number, c: readonly number[]) => { ctx.fillStyle = css(c); ctx.fillRect(X(g) - 1, y0 + 18, 2, h - 22); };
      mark((doc.enhancer[0] + doc.enhancer[1]) / 2, C.enh);
      mark(doc.genes.find((g) => g.name === "PAX3")!.tss, C.pax3);
      ctx.fillStyle = "#8b98aa";
      ctx.font = "10px system-ui, sans-serif";
      ctx.fillText("10⁰", W1 - 22, Y(1) + 10);
      ctx.fillText("10⁻⁴·⁵", W1 - 34, Y(10 ** lo) - 2);
    };
    panel(0, PL_H / 2, "4C virtuel depuis le promoteur de PAX3 (vermillon) — gris : type sauvage, jaune : allèle", "PAX3");
    panel(PL_H / 2, PL_H / 2, "4C virtuel depuis les enhancers (jaune)", "enhancer");
    ctx.restore();
  };

  const fmt = (x: number) => (x >= 10 ? x.toFixed(0) : x.toFixed(x >= 1 ? 2 : 3));
  const writeReadout = () => {
    if (!wt || !cur) return;
    const r = cur.contact.PAX3 / wt.contact.PAX3;
    const ctl = cur.allele.deletion ? controls.get(`${cur.allele.deletion.start}`) : null;
    const d = cur.allele.deletion;
    readout.innerHTML =
      `<div>${d ? `<b>${d.name ?? "délétion tracée"}</b> · ${((d.end - d.start) / 1e6).toFixed(2)} Mb, chr2:${(d.start + 1).toLocaleString("fr")}–${d.end.toLocaleString("fr")}`
        : "<b>type sauvage</b>"}</div>` +
      `<div>enhancers → PAX3 : <b>×${fmt(r)}</b> par rapport au type sauvage` +
      (ctl ? ` · ×${fmt(cur.contact.PAX3 / ctl.contact.PAX3)} par rapport au même allèle, CTCF de la frontière épargnés` : "") + `</div>` +
      `<div>enhancers → EPHA4 : ${Number.isNaN(cur.contact.EPHA4) ? "EPHA4 supprimé" : `×${fmt(cur.contact.EPHA4 / wt.contact.EPHA4)}`} · ` +
      `enhancers → IHH : ×${fmt(cur.contact.IHH / wt.contact.IHH)}</div>` +
      `<div class="k">prédiction 1D : ${cur.ms.extrusion.toFixed(0)} + ${cur.ms.contacts.toFixed(0)} ms (extrusion, contacts) ·
        3D : ${doc.alleles[current] ? "OpenMM, conformation représentative de 1 200" : "non recalculée — type sauvage, segment supprimé effacé"}</div>` +
      `<div class="src">Coordonnées des variants et des enhancers reconstruites depuis le texte de Lupiáñez 2015 (Table S1 non lue) ;
        CTCF et TAD mesurés dans GM12878 ; structures et contacts simulés. Loi des contacts calée sur le 3D du type sauvage
        (α = ${doc.law.alpha}, d0 = ${doc.law.d0}).</div>`;
    cap.textContent = `${doc.alleles[current]?.name ?? (cur.allele.deletion?.name ?? "délétion tracée")} · ${doc.evidence.structures}`;
  };

  // Le témoin « frontière épargnée » d'une délétion qui emporte la frontière, calculé avec elle.
  const controls = new Map<string, Prediction>();

  const timings: Record<string, unknown>[] = [];
  const choose = (key: string, del: Deletion | null, t0 = performance.now()) => {
    current = key;
    for (const [k, b] of buttons) b.setAttribute("aria-pressed", String(k === key));
    curDel = del;
    draft = null;
    const t1 = performance.now();
    cur = del ? predict(doc, del) : wt;
    const removes = del && del.start < doc.boundary[0] && del.end > doc.boundary[1];
    if (del && removes && !del.keep?.length) controls.set(`${del.start}`, predict(doc, { ...del, keep: [doc.boundary] }));
    const t2 = performance.now();
    // 3D : la conformation simulée de l'allèle si elle existe ; sinon celle du type sauvage,
    // segment supprimé effacé.
    from = pos.slice();
    alphaFrom = alpha.slice();
    const a = doc.alleles[key];
    if (a) {
      const xyz = decodeCoords(a.coords);
      to = pos.slice();
      alphaTo = new Float32Array(n0);
      const first = Math.round((a.mids_first - R.start) / R.bp_per_bead - 0.5);
      let j = 0;
      for (let k = first; k < n0 && j < a.n; k++) {
        if (del && mid(k) >= del.start && mid(k) < del.end) continue;
        to.set([xyz[3 * j] * um, xyz[3 * j + 1] * um, xyz[3 * j + 2] * um], 3 * k);
        alphaTo[k] = 1;
        j++;
      }
      for (let k = 0; k < n0; k++) if (!alphaTo[k]) to.set(from.subarray(3 * k, 3 * k + 3), 3 * k);
      joined = true;
    } else {
      to = new Float32Array(n0 * 3);
      for (let i = 0; i < n0 * 3; i++) to[i] = wtXYZ[i] * um;
      alphaTo = new Float32Array(n0).fill(1);
      if (del) for (let k = 0; k < n0; k++) if (mid(k) >= del.start && mid(k) < del.end) alphaTo[k] = 0.0;
      joined = false;
    }
    anim = { t0: performance.now(), ms: 1200 };
    drawTracks();
    drawPlots();
    writeReadout();
    const t3 = performance.now();
    const entry = { key, del: del && { start: del.start, end: del.end }, predictMs: t2 - t1, totalMs: t3 - t0,
                    contact: cur!.contact, control: removes ? controls.get(`${del!.start}`)!.contact : null, ms: cur!.ms };
    timings.push(entry);
    return entry;
  };

  // Tracer une délétion sur la piste.
  let down: number | null = null;
  tracks.addEventListener("pointerdown", (e) => {
    down = gOf(e.offsetX, tracks.clientWidth);
    tracks.setPointerCapture(e.pointerId);
  });
  tracks.addEventListener("pointermove", (e) => {
    if (down === null) return;
    const g = gOf(e.offsetX, tracks.clientWidth);
    draft = [Math.min(down, g), Math.max(down, g)];
    drawTracks();
  });
  tracks.addEventListener("pointerup", (e) => {
    if (down === null) return;
    const t0 = performance.now();
    const g = gOf(e.offsetX, tracks.clientWidth);
    const [a, b] = [Math.min(down, g), Math.max(down, g)].map((v) => Math.round(v / 1000) * 1000);
    down = null;
    if (b - a < 4 * R.bp_per_bead) { draft = null; drawTracks(); return; }
    choose("drawn", { name: "délétion tracée", start: Math.max(R.start, a), end: Math.min(R.end, b) }, t0);
  });

  const t0 = performance.now();
  wt = predict(doc, null);
  const wtMs = performance.now() - t0;
  recolor();
  choose("wt", null);

  const loop = (now: number) => {
    if (anim) {
      const k = Math.min(1, (now - anim.t0) / anim.ms);
      const e = k * k * (3 - 2 * k);
      for (let i = 0; i < pos.length; i++) pos[i] = from[i] + (to[i] - from[i]) * e;
      for (let i = 0; i < n0; i++) alpha[i] = alphaFrom[i] + (alphaTo[i] - alphaFrom[i]) * e;
      upload();
      if (k >= 1) anim = null;
      dirty = true;
    }
    if (dirty) {
      render();
      dirty = false;
    }
    requestAnimationFrame(loop);
  };
  upload();
  requestAnimationFrame(loop);

  window.__variant = {
    doc,
    timings,
    wtMs,
    presets: presets.map((p) => p.key),
    choose: (key: string) => { const p = presets.find((x) => x.key === key)!; return choose(key, p.del); },
    family: () => dels,
    /** Pour la mesure : la boîte de la piste dans la page, et le x (CSS) d'une coordonnée. */
    track: () => { const b = tracks.getBoundingClientRect(); return { left: b.left, top: b.top, width: b.width, height: b.height }; },
    xOf: (g: number) => tracks.getBoundingClientRect().left + xOf(g, tracks.clientWidth),
    idle: () => !anim && !dirty,
  };
  void central;
  publish({ ok: true, variant: true, alleles: Object.keys(doc.alleles), wtMs });
}

main().catch((e) => publish({ ok: false, error: String((e as Error)?.stack ?? e) }));
