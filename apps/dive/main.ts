/** La plongée : du noyau entier aux compartiments d'un chromosome, sans coupure.
 *
 *  Deux niveaux du même `.g3d` — l'aperçu à 3 Mb, puis le noyau à 750 kb —, fondus l'un
 *  dans l'autre par tramage complémentaire ; les 45 autres copies qui s'effacent autour de
 *  chr7:a ; les couleurs de territoire qui glissent vers celles de compartiment ; une
 *  occlusion ambiante qui rend le volume lisible. Tout dépend de la caméra (`scene.ts`).
 *
 *  Paramètres d'URL :
 *    f=…        fichier .g3d (défaut data/gm12878.g3d)
 *    bench=N    plongée scriptée de N images, mesurée, puis résultat publié dans window.__geno
 *    hard       témoin : transitions sèches au lieu de fondus (la mesure doit y voir des sauts)
 *    w, h       taille du canevas en mode mesure (défaut 960 × 540)
 */

import { dive, type Orbit, rotate, zoom } from "@geno/viewer";
import { runBench } from "./bench.js";
import { FOV, openDive } from "./dive.js";

const q = new URLSearchParams(location.search);
const FILE = q.get("f") ?? "data/gm12878.g3d";
const BENCH = Number(q.get("bench") ?? 0);
const HARD = q.has("hard");

declare global {
  interface Window {
    __geno?: Record<string, unknown>;
  }
}

const style = document.createElement("style");
style.textContent = `
  html, body { height: 100%; overflow: hidden; }
  #hud { position: fixed; left: 12px; top: 10px; max-width: min(460px, calc(100vw - 24px));
         background: rgba(14,18,24,.82); padding: 10px 12px; border-radius: 8px; font-size: 12px; }
  #hud h1 { font-size: 13px; margin: 0 0 6px; font-weight: 600; }
  #hud button { font: inherit; margin: 6px 6px 0 0; padding: 4px 10px; border-radius: 6px;
                border: 1px solid #3a4658; background: #1b2330; color: inherit; cursor: pointer; }
  #hud .k { color: #8b98aa; }
  .sw { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin: 0 4px -1px 8px; }
  #info { margin-top: 6px; white-space: pre-wrap; }
  canvas { touch-action: none; }
`;
document.head.appendChild(style);

const canvas = document.createElement("canvas");
document.body.appendChild(canvas);

function fit(): [number, number] {
  if (BENCH) {
    const w = Number(q.get("w") ?? 960);
    const h = Number(q.get("h") ?? 540);
    canvas.width = w;
    canvas.height = h;
    canvas.style.width = `${w}px`;
    canvas.style.height = `${h}px`;
    return [w, h];
  }
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.round(window.innerWidth * dpr);
  canvas.height = Math.round(window.innerHeight * dpr);
  canvas.style.width = `${window.innerWidth}px`;
  canvas.style.height = `${window.innerHeight}px`;
  return [canvas.width, canvas.height];
}

function publish(r: Record<string, unknown>): void {
  window.__geno = r;
}

async function main(): Promise<void> {
  fit();
  const scene = await openDive(canvas, FILE, { hard: HARD });

  if (BENCH) {
    const r = await runBench(scene, BENCH, { continuity: true });
    publish({
      ok: r.diffs.length === BENCH && Number.isFinite(r.diffs[0]),
      hard: HARD,
      frames: r.frames,
      size: r.size,
      gpu: r.gpu,
      path_length: r.pathLength,
      first_render_ms: scene.firstRenderMs,
      instances: r.instances,
      diffs: r.diffs,
      frame_ms: r.frameMs,
      raf_ms: r.rafMs,
      states: r.states,
    });
    return;
  }

  // Mode interactif : glisser tourne, la molette zoome, un clic sans glisser identifie.
  let orbit: Orbit = { ...scene.overview };
  let anim: { path: ReturnType<typeof dive>; t0: number; duration: number } | null = null;
  let dirty = true;
  const arrivedAt = performance.now();
  const frameTimes: number[] = [];
  let lastFrame = performance.now();

  const hud = document.createElement("div");
  hud.id = "hud";
  hud.innerHTML = `<h1>Geno — plongée du noyau aux compartiments</h1>
    <div class="k">GM12878 · modèle simulé (médoïde de 200 repliements) · compartiments appelés sur la
    carte de contacts du modèle, pas sur du Hi-C mesuré</div>
    <div><span class="sw" style="background:rgb(224,107,51)"></span>A
         <span class="sw" style="background:rgb(56,107,219)"></span>B
         <span class="sw" style="background:rgb(128,133,140)"></span>sans compartiment (Xi)</div>`;
  document.body.appendChild(hud);
  const go = (label: string, to: () => Orbit) => {
    const b = document.createElement("button");
    b.textContent = label;
    b.onclick = () => {
      const path = dive(orbit, to(), FOV);
      anim = { path, t0: performance.now(), duration: 900 * Math.max(0.5, path.length) };
    };
    hud.appendChild(b);
  };
  go("Plonger vers ACTB", () => scene.final);
  go("Remonter", () => scene.overview);
  const status = document.createElement("div");
  const info = document.createElement("div");
  info.id = "info";
  hud.appendChild(status);
  hud.appendChild(info);

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
  canvas.addEventListener("pointerup", async (e) => {
    const wasClick = drag && drag.moved < 4;
    drag = null;
    if (!wasClick) return;
    const dpr = canvas.width / canvas.clientWidth;
    info.textContent = (await scene.describe(e.offsetX * dpr, e.offsetY * dpr)) ?? "";
  });
  canvas.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      orbit = zoom(orbit, Math.exp(e.deltaY * 0.0015), 0.4, 150);
      anim = null;
      dirty = true;
    },
    { passive: false },
  );
  window.addEventListener("resize", () => {
    scene.resize(...fit());
    dirty = true;
  });

  const loop = (now: number) => {
    if (scene.fine && scene.arrival < 1) {
      scene.arrival = Math.min(1, (now - arrivedAt) / 400);   // seul fondu réglé par le temps
      dirty = true;
    }
    if (anim) {
      const k = Math.min(1, (now - anim.t0) / anim.duration);
      orbit = anim.path.at(k * anim.path.length);
      if (k >= 1) anim = null;
      dirty = true;
    }
    if (dirty) {
      scene.render(orbit);
      dirty = false;
      frameTimes.push(now - lastFrame);
      if (frameTimes.length > 60) frameTimes.shift();
      const sorted = frameTimes.slice().sort((a, b) => a - b);
      const s = scene.last;
      status.innerHTML =
        `<span class="k">champ</span> ${s.widthUm.toFixed(1)} µm · ` +
        `<span class="k">750 kb</span> ${(s.lod * 100).toFixed(0)} % · ` +
        `<span class="k">autres copies</span> ${((1 - s.focus) * 100).toFixed(0)} % · ` +
        `<span class="k">compartiments</span> ${(s.mix * 100).toFixed(0)} % · ` +
        `<span class="k">image</span> ${sorted[sorted.length >> 1]?.toFixed(1) ?? "–"} ms<br>` +
        `<span class="k">${scene.gpu}</span>`;
    }
    lastFrame = now;
    requestAnimationFrame(loop);
  };
  requestAnimationFrame(loop);
}

main().catch((e) => publish({ ok: false, error: String(e) }));
