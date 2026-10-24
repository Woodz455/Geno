/** Les vues liées (semaine 14) : la 3D, igv.js et HiGlass, une seule sélection.
 *
 *  Chaque vue publie dans le store ce que l'utilisateur y sélectionne et applique ce que
 *  le store lui donne (`store.ts`) ; aucune ne parle à une autre.
 *
 *  Paramètres d'URL :
 *    f=…       fichier .g3d (défaut data/gm12878.g3d) ; les fiches, les contacts et les
 *              pistes igv sont à côté (gm12878.fiche, gm12878.contacts, igv/)
 *    naive     témoin : pas de store — les vues se synchronisent deux à deux sur leur
 *              déplacement (igv `locuschange` → HiGlass et 3D, HiGlass `location` → igv et
 *              3D), la conception que la semaine écarte. La mesure doit y voir des rebonds.
 */

import { openHiglass } from "./viewHiglass.js";
import { type IgvIndex, openIgv } from "./viewIgv.js";
import { openView3D } from "./view3d.js";
import { type Interval, key, type Selection, SelectionStore, type ViewName } from "./store.js";

// Les erreurs et avertissements des trois bibliothèques, gardés : la mesure les compte.
const errors: string[] = [];
for (const k of ["error", "warn"] as const) {
  const orig = console[k].bind(console);
  console[k] = (...a: unknown[]) => {
    errors.push(`${k}: ${a.map((x) => (x instanceof Error ? x.message : String(x))).join(" ").slice(0, 300)}`);
    orig(...a);
  };
}
window.addEventListener("error", (e) => errors.push(`exception: ${e.message}`));
window.addEventListener("unhandledrejection", (e) => errors.push(`rejet: ${String(e.reason?.stack ?? e.reason).slice(0, 600)}`));

// L'instant de chaque geste de l'utilisateur, sur l'horloge de performance.now() : d'où part
// la latence d'une sélection. Les événements de confiance seulement — pixi.js (sous HiGlass)
// réémet un `pointerup` synthétique 40 ms après le vrai, qui écrasait l'instant du geste.
const input = { at: 0, kind: "" };
for (const k of ["pointerup", "keydown", "wheel"]) {
  window.addEventListener(k, (e) => {
    if (!e.isTrusted) return;
    input.at = performance.now();
    input.kind = k;
  }, true);
}

const q = new URLSearchParams(location.search);
const FILE = q.get("f") ?? "data/gm12878.g3d";
const NAIVE = q.has("naive");

declare global {
  interface Window {
    __geno?: Record<string, unknown>;
    __linked?: unknown;
  }
}

const style = document.createElement("style");
style.textContent = `
  html, body { height: 100%; overflow: hidden; }
  body { display: grid; grid-template-columns: 1fr 1fr; grid-template-rows: 30px 1fr 1fr; gap: 6px; padding: 6px; box-sizing: border-box; }
  #status { grid-column: 1 / 3; display: flex; gap: 12px; align-items: center; font-size: 12px; white-space: nowrap; overflow: hidden; }
  #status .k { color: #8b98aa; }
  #status .sel { color: rgb(240,228,66); }
  .pane { position: relative; border: 1px solid #2a3442; border-radius: 6px; overflow: hidden; min-height: 0; }
  .pane > .cap { position: absolute; left: 8px; top: 6px; z-index: 6; font-size: 11px; color: #c9d3df;
                 background: rgba(14,18,24,.82); padding: 2px 8px; border-radius: 6px; pointer-events: none; }
  #p3d { grid-row: 2 / 4; }
  #p3d .bar { position: absolute; left: 8px; bottom: 8px; z-index: 6; }
  #p3d .bar button { font: inherit; font-size: 12px; margin-right: 6px; padding: 4px 10px; border-radius: 6px;
                     border: 1px solid #3a4658; background: #1b2330; color: inherit; cursor: pointer; }
  #pigv { background: #fff; color: #222; overflow: auto; }
  #pigv > .cap { position: sticky; display: inline-block; }
  #phg { background: #fff; }
  #phg > .cap { top: auto; bottom: 6px; }
  /* La fiche contre le bord gauche de la 3D, étroite : au centre, elle couvrait le noyau. */
  #fiche { right: auto !important; left: 12px !important; top: 46px !important; bottom: 56px !important; width: 340px !important; }
`;
document.head.appendChild(style);

const pane = (id: string, caption: string) => {
  const el = document.createElement("section");
  el.className = "pane";
  el.id = id;
  const cap = document.createElement("div");
  cap.className = "cap";
  cap.textContent = caption;
  el.appendChild(cap);
  const body = document.createElement("div");
  body.style.cssText = "position:absolute;inset:0";
  el.appendChild(body);
  return { el, body, cap };
};

const status = document.createElement("div");
status.id = "status";
document.body.appendChild(status);
const p3 = pane("p3d", "3D · noyau, fibre et bases — modèle simulé, séquence hg38");
const pi = pane("pigv", "igv.js · hg38 · GENCODE, CTCF et boucles mesurés, compartiments simulés");
const ph = pane("phg", "HiGlass · contacts du modèle — simulés, pas du Hi-C mesuré");
document.body.append(p3.el, pi.el, ph.el);
pi.body.style.cssText = "position:relative;padding-top:24px";

function publishResult(r: Record<string, unknown>) {
  window.__geno = r;
}

async function main() {
  const store = new SelectionStore();
  const base = FILE.replace(/[^/]*$/, "");
  const index = (await (await fetch(`${base}igv/index.json`)).json()) as IgvIndex;
  // Le témoin : une file de messages deux-à-deux, comptée. Plafonnée — sans plafond, la
  // boucle ne s'arrête pas, et c'est précisément ce qu'on veut montrer sans geler la page.
  const naive = { hops: [] as { t: number; from: string; to: string; locus: string }[], cap: 400 };
  const pub = (view: ViewName) => (regions: Interval[], label: string) => {
    if (NAIVE) return;
    store.publish(view, regions, label);
  };

  const [v3d, igvView, hg] = await Promise.all([
    openView3D(p3.body, FILE, FILE.replace(/\.g3d$/, ".fiche"), pub("3d")),
    openIgv(pi.body, index, pub("igv")),
    openHiglass(ph.body, FILE.replace(/\.g3d$/, ".contacts"), pub("higlass")),
  ]);

  const describe = (s: Selection | null) =>
    s ? `<span class="k">sélection</span> <span class="sel">${s.label || "—"}</span> · ${s.regions.map((r) => `${r.chrom}:${(r.start + 1).toLocaleString("fr")}–${r.end.toLocaleString("fr")}`).join(" × ")} <span class="k">· publiée par ${s.source} · n° ${s.seq}</span>`
      : `<span class="k">aucune sélection — cliquez une bille, un élément d'une piste igv ou un pixel de la carte</span>`;
  const show = (s: Selection | null) => (status.innerHTML = describe(s));
  show(null);

  if (!NAIVE) {
    store.subscribe("3d", (s) => v3d.apply(s));
    store.subscribe("igv", (s) => igvView.apply(s));
    store.subscribe("higlass", (s) => hg.apply(s));
    // La source ne reçoit pas sa propre sélection ; la barre d'état, elle, les voit toutes.
    store.observe(show);
  } else {
    // Deux-à-deux, sur le déplacement de chaque vue : ce que font beaucoup d'intégrations.
    const hop = (from: string, to: string, locus: string, fn: () => void) => {
      if (naive.hops.length >= naive.cap) return;
      naive.hops.push({ t: performance.now(), from, to, locus });
      fn();
    };
    const sel = (r: Interval): Selection => ({ regions: [r], source: "igv", label: "déplacement", seq: naive.hops.length, at: performance.now() });
    igvView.browser.on("locuschange", (loci: { chr: string; start: number; end: number }[]) => {
      const l = loci[0];
      if (!l) return;
      const r = { chrom: l.chr, start: Math.round(l.start), end: Math.round(l.end) };
      const a = hg.reader.absolute(r.chrom, r.start), b = hg.reader.absolute(r.chrom, r.end);
      hop("igv", "higlass", key([r]), () => a !== null && b !== null && hg.goto(a, b));
      hop("igv", "3d", key([r]), () => v3d.apply(sel(r)));
    });
    hg.api.on("location", (e: { xDomain: [number, number] }) => {
      const a = hg.reader.relative(e.xDomain[0]);
      const b = hg.reader.relative(e.xDomain[1]);
      if (!a || !b) return;
      const end = a.chrom === b.chrom ? b.pos : hg.reader.header.chromsizes.find(([c]) => c === a.chrom)![1];
      const r = { chrom: a.chrom, start: a.pos, end };
      hop("higlass", "igv", key([r]), () => void igvView.browser.search(`${r.chrom}:${r.start + 1}-${r.end}`));
      hop("higlass", "3d", key([r]), () => v3d.apply(sel(r)));
    }, "v");
  }

  window.__linked = {
    errors,
    input,
    naive: NAIVE ? naive : null,
    store,
    v3d,
    igv: igvView,
    hg,
    /** Toutes les vues au repos : rien à dessiner, rien à charger. */
    idle: () => v3d.idle() && igvView.idle() && hg.idle(),
  };
  publishResult({ ok: true, linked: true, naive: NAIVE, gpu: v3d.scene.gpu });
}

main().catch((e) => publishResult({ ok: false, error: String((e as Error)?.stack ?? e) }));
