/** La carte de contacts des vues liées : HiGlass, sur les tuiles du modèle (`make linked`).
 *
 *  Deux cartes superposées — le noyau (750 kb, tout le génome) dessous, les régions fines
 *  (2 kb, TP53 et ACTB) dessus, NaN ailleurs pour laisser voir le noyau. Toutes deux
 *  **simulées** : la carte est celle du modèle, pas du Hi-C mesuré, et la page le dit.
 *
 *  Les tuiles viennent d'un data fetcher maison (`GenoContacts`) qui lit `gm12878.contacts`
 *  par plages. Sélectionner, c'est cliquer un pixel de la carte : ses deux loci, à la
 *  résolution affichée. Naviguer n'est pas sélectionner ; une sélection reçue est dessinée
 *  par-dessus (bandes de ses intervalles, rectangle du pixel), et la vue s'y déplace si elle
 *  n'y est pas.
 */

import { ContactsReader, HttpSource, TILE, type Tileset } from "@geno/g3d";
import { viewer } from "higlass";
import type { Interval, Selection } from "./store.js";

let READER: Promise<ContactsReader> | null = null;

/** Les résolutions annoncées à HiGlass. La carte fine reçoit aussi celles, plus grossières, du
 *  noyau — servies vides, donc transparentes : sans elles, au génome entier, HiGlass lui
 *  demandait ses tuiles de 512 kb, 24 par côté (« too many visible tiles »), toutes vides
 *  sauf deux. */
function resolutionsOf(r: ContactsReader, ts: Tileset): number[] {
  if (ts.uid !== "fin") return ts.resolutions;
  return [...r.tileset("noyau").resolutions.filter((x) => x > Math.max(...ts.resolutions)), ...ts.resolutions];
}

interface HGC {
  utils: { DenseDataExtrema2D: new (a: Float32Array) => { minNonZeroInTile: number; maxNonZeroInTile: number } };
}

/** Le data fetcher que HiGlass instancie pour `data: { type: "geno-contacts", tileset }`. */
class GenoContacts {
  static pending = 0;
  /** Tuiles demandées par jeu, et ce qui a échoué : pour la mesure et le diagnostic. */
  static stats: Record<string, number> & { errors?: string[] } = {};

  // `dataConfig` : HiGlass lit ce champ du fetcher par son nom (`tilesetUid`, `server`).
  constructor(private readonly hgc: HGC, readonly dataConfig: { tileset: string }, _pubSub: unknown) {}

  async tilesetInfo(cb?: (info: unknown) => void): Promise<unknown> {
    const r = await READER!;
    const ts = r.tileset(this.dataConfig.tileset);
    // `chromsizes` : les pistes d'étiquettes de chromosomes s'en servent, sans serveur.
    const info = { min_pos: ts.min_pos, max_pos: ts.max_pos, resolutions: resolutionsOf(r, ts), bins_per_dimension: TILE,
                   name: ts.label, coordSystem: "hg38", chromsizes: r.header.chromsizes };
    cb?.(info);
    return info;
  }

  async fetchTilesDebounced(received: (tiles: Record<string, unknown>) => void, tileIds: string[]): Promise<Record<string, unknown>> {
    GenoContacts.pending++;
    try {
      const r = await READER!;
      const uid = this.dataConfig.tileset;
      const out: Record<string, unknown> = {};
      await Promise.all(tileIds.map(async (id) => {
        const [z, x, y] = id.split(".").slice(-3).map(Number);
        // Les niveaux prêtés à la carte fine (voir `resolutionsOf`) n'ont pas de tuile.
        const zf = z - (resolutionsOf(r, r.tileset(uid)).length - r.tileset(uid).resolutions.length);
        // Seules les tuiles x ≤ y sont écrites ; HiGlass n'en demande pas d'autres, mais une
        // autre se lirait transposée.
        const t = zf < 0 ? null : await r.tile(uid, zf, Math.min(x, y), Math.max(x, y));
        const dense = new Float32Array(TILE * TILE);
        if (!t) dense.fill(uid === "fin" ? NaN : 0);
        else if (x <= y) dense.set(t);
        else for (let i = 0; i < TILE; i++) for (let j = 0; j < TILE; j++) dense[i * TILE + j] = t[j * TILE + i];
        const dde = new this.hgc.utils.DenseDataExtrema2D(dense);
        out[id] = { dense, denseDataExtrema: dde, minNonZero: dde.minNonZeroInTile, maxNonZero: dde.maxNonZeroInTile,
                    tilePos: [x, y], zoomLevel: z, tilePositionId: id, tileId: id, tilesetUid: uid, server: "" };
        GenoContacts.stats[uid] = (GenoContacts.stats[uid] ?? 0) + 1;
      }));
      received(out);
      return out;
    } catch (e) {
      (GenoContacts.stats.errors ??= []).push(String(e));
      throw e;
    } finally {
      GenoContacts.pending--;
    }
  }
}

export interface HgApi {
  on(event: string, fn: (e: any) => void, viewId?: string): unknown;
  getLocation(viewId: string): { xDomain: [number, number]; yDomain: [number, number]; xRange: [number, number]; yRange: [number, number] };
  zoomTo(viewId: string, x0: number, x1: number, y0: number, y1: number, animate?: number): void;
}

export interface ViewHiglass {
  api: HgApi;
  reader: ContactsReader;
  root: HTMLElement;
  shown: { seq: number; at: number; applied: number; navigated: number };
  apply(sel: Selection | null): void;
  /** La boîte de la carte dans la page (coins en pixels CSS). */
  box(): { left: number; top: number; width: number; height: number } | null;
  /** Ce que la surimpression dessine : chaque intervalle, et s'il est dans la vue. */
  drawn(): { region: string; visible: boolean }[];
  /** Pour la mesure : des points de page sur la carte, à graine fixe. */
  targets(n: number, seed: number): [number, number][];
  idle(): boolean;
  goto(x0: number, x1: number): void;
  /** La résolution affichée, et la carte qui la sert. */
  resolution(at?: [number, number]): { res: number; tileset: string };
  stats: Record<string, unknown>;
}

const VIEW = "v";
/** Durée du déplacement programmé. Pas 0 : sans transition, HiGlass déplace la vue sans
 *  prévenir ses écouteurs (`location`), et la surimpression resterait sur l'ancienne. */
const ZOOM_MS = 150;
const fmtRes = (r: number) => (r >= 1e6 ? `${r / 1e6} Mb` : `${r / 1e3} kb`);

export async function openHiglass(root: HTMLElement, file: string, publish: (regions: Interval[], label: string) => void,
                                  start: Interval = { chrom: "chr17", start: 3_000_000, end: 13_000_000 }): Promise<ViewHiglass> {
  READER = ContactsReader.open(new HttpSource(file));
  const reader = await READER;
  const noyau = reader.tileset("noyau");
  const fin = reader.tileset("fin");
  const a0 = reader.absolute(start.chrom, start.start)!;
  const a1 = reader.absolute(start.chrom, start.end)!;
  const heat = (uid: string, tileset: string, name: string, colorbar: string) => ({
    uid, type: "heatmap", data: { type: "geno-contacts", tileset },
    // Fond transparent : la carte fine, NaN hors de ses régions, laisse voir le noyau dessous
    // (le fond par défaut d'une heatmap, gris clair, le cachait).
    options: { name, labelPosition: "hidden", colorbarPosition: colorbar, heatmapValueScaling: "log", backgroundColor: "transparent",
               colorRange: ["white", "rgba(245,166,35,1.0)", "rgba(208,2,27,1.0)", "black"], showTooltip: false },
  });
  const api = (await viewer(root, {
    editable: false,
    views: [{
      uid: VIEW,
      initialXDomain: [a0, a1], initialYDomain: [a0, a1],
      tracks: {
        top: [{ uid: "tl", type: "horizontal-chromosome-labels", data: { type: "geno-contacts", tileset: "noyau" }, height: 24 }],
        left: [{ uid: "ll", type: "vertical-chromosome-labels", data: { type: "geno-contacts", tileset: "noyau" }, width: 24 }],
        center: [{ uid: "c", type: "combined", contents: [
          heat("hn", "noyau", "contacts du noyau — simulés", "topRight"),
          heat("hf", "fin", "contacts des régions fines — simulés", "bottomRight"),
        ] }],
      },
      layout: { w: 12, h: 12, x: 0, y: 0 },
    }],
  } as never, { bounded: true, pluginDataFetchers: { "geno-contacts": { dataFetcher: GenoContacts } } } as never)) as unknown as HgApi;

  const shown: ViewHiglass["shown"] = { seq: 0, at: 0, applied: 0, navigated: 0 };
  let current: Selection | null = null;
  let settleUntil = 0;
  let trackOrigin: [number, number] | null = null;
  let cursor: { dataX: number; dataY: number; is2d: boolean } | null = null;
  let last: ReturnType<HgApi["getLocation"]> | null = null;
  const where = () => last ?? api.getLocation(VIEW);

  api.on("cursorLocation", (e) => {
    if (e.isFrom2dTrack) trackOrigin = [e.absX - e.relTrackX, e.absY - e.relTrackY];
    cursor = { dataX: e.dataX, dataY: e.dataY, is2d: !!e.isFrom2dTrack };
  });

  const overlay = document.createElement("div");
  overlay.style.cssText = "position:absolute;inset:0;pointer-events:none;overflow:hidden;z-index:5";
  root.appendChild(overlay);

  const box = () => {
    if (!trackOrigin) return null;
    const loc = where();
    const r = root.getBoundingClientRect();
    return { left: trackOrigin[0] - r.left, top: trackOrigin[1] - r.top, width: loc.xRange[1] - loc.xRange[0], height: loc.yRange[1] - loc.yRange[0] };
  };

  /** La résolution que HiGlass affiche en un point (absolu, x et y) : la plus fine qui a
   *  moins d'un bin par pixel, dans la carte fine si le point est dans le carré d'une région
   *  fine, dans celle du noyau sinon. Sans point : celle du centre de la vue. */
  const resolution = (at?: [number, number]): { res: number; tileset: string } => {
    const loc = where();
    const bpPerPx = (loc.xDomain[1] - loc.xDomain[0]) / Math.max(1, loc.xRange[1] - loc.xRange[0]);
    const pick = (ts: Tileset) => {
      const ok = ts.resolutions.filter((r) => r > bpPerPx);
      return ok.length ? Math.min(...ok) : Math.max(...ts.resolutions);
    };
    const [x, y] = at ?? [(loc.xDomain[0] + loc.xDomain[1]) / 2, (loc.yDomain[0] + loc.yDomain[1]) / 2];
    const inFine = (fin.regions ?? []).some((g) => x >= g.abs_start && x < g.abs_end && y >= g.abs_start && y < g.abs_end);
    const f = pick(fin);
    return inFine && f < Math.min(...noyau.resolutions) ? { res: f, tileset: "fin" } : { res: pick(noyau), tileset: "noyau" };
  };

  const toInterval = (abs: number, res: number): Interval | null => {
    const b0 = Math.floor(abs / res) * res;
    const p = reader.relative(b0);
    if (!p) return null;
    const size = reader.header.chromsizes.find(([c]) => c === p.chrom)![1];
    return { chrom: p.chrom, start: p.pos, end: Math.min(size, p.pos + res) };
  };

  const draw = () => {
    overlay.textContent = "";
    const b = box();
    if (!current || !b) return;
    const loc = where();
    const sx = (v: number) => b.left + ((v - loc.xDomain[0]) / (loc.xDomain[1] - loc.xDomain[0])) * b.width;
    const sy = (v: number) => b.top + ((v - loc.yDomain[0]) / (loc.yDomain[1] - loc.yDomain[0])) * b.height;
    const rect = (x0: number, x1: number, y0: number, y1: number, css: string, region: string) => {
      const l = Math.max(b.left, Math.min(sx(x0), sx(x1) - 2)), r = Math.min(b.left + b.width, Math.max(sx(x1), sx(x0) + 2));
      const t = Math.max(b.top, Math.min(sy(y0), sy(y1) - 2)), d = Math.min(b.top + b.height, Math.max(sy(y1), sy(y0) + 2));
      if (r <= l || d <= t) return;
      const el = document.createElement("div");
      el.className = "hg-sel";
      el.dataset.region = region;
      el.style.cssText = `position:absolute;left:${l}px;top:${t}px;width:${r - l}px;height:${d - t}px;${css}`;
      overlay.appendChild(el);
    };
    const abs = current.regions.map((g) => [reader.absolute(g.chrom, g.start), reader.absolute(g.chrom, g.end)] as const);
    current.regions.forEach((g, k) => {
      const [s, e] = abs[k];
      if (s === null || e === null) return;
      const key = `${g.chrom}:${g.start}-${g.end}`;
      rect(s, e, loc.yDomain[0], loc.yDomain[1], "background:rgba(240,228,66,0.22)", key);   // bande verticale (colonnes)
      rect(loc.xDomain[0], loc.xDomain[1], s, e, "background:rgba(240,228,66,0.22)", key);   // bande horizontale (lignes)
    });
    if (abs.length === 2 && abs.every(([s, e]) => s !== null && e !== null)) {
      const [[a0, a1], [b0, b1]] = abs as [number, number][];
      for (const [x0, x1, y0, y1] of [[a0, a1, b0, b1], [b0, b1, a0, a1]]) {
        rect(x0, x1, y0, y1, "outline:2px solid rgb(240,228,66);background:rgba(240,228,66,0.35)", "pixel");
      }
    }
  };
  // La vue telle que HiGlass la publie : `getLocation()` lit des échelles que `zoomTo` ne
  // remplace pas tout de suite — elle rendait l'ancienne vue après un déplacement.
  api.on("location", (e) => {
    last = e;
    settleUntil = performance.now() + 120;
    draw();
  }, VIEW);

  // Un clic sans glisser : le pixel sous le pointeur, à la résolution affichée.
  let down: { x: number; y: number } | null = null;
  root.addEventListener("pointerdown", (e) => (down = { x: e.clientX, y: e.clientY }), true);
  root.addEventListener("pointerup", (e) => {
    const d = down;
    down = null;
    if (!d || Math.abs(e.clientX - d.x) + Math.abs(e.clientY - d.y) > 4) return;
    const c = cursor;
    if (!c || !c.is2d || c.dataX < 0 || c.dataY < 0) return;
    const { res } = resolution([c.dataX, c.dataY]);
    const x = toInterval(c.dataX, res);
    const y = toInterval(c.dataY, res);
    if (!x || !y) return;
    const same = x.chrom === y.chrom && x.start === y.start;
    const regions = same ? [x] : [x, y];
    current = { regions, source: "higlass", label: "", seq: -1, at: performance.now() };
    draw();
    publish(regions, `pixel de la carte (${fmtRes(res)}${same ? ", diagonale" : ""})`);
  }, true);

  // Calage de la surimpression : un mouvement de souris synthétique au centre de la carte
  // donne la position de la piste centrale dans la page (HiGlass ne l'expose pas autrement).
  const calibrate = () => {
    const r = root.getBoundingClientRect();
    const x = r.left + r.width / 2, y = r.top + r.height / 2;
    document.elementFromPoint(x, y)?.dispatchEvent(new MouseEvent("mousemove", { clientX: x, clientY: y, bubbles: true }));
  };
  calibrate();
  new ResizeObserver(() => {
    trackOrigin = null;
    calibrate();
    draw();
  }).observe(root);

  const contained = (sel: Selection) => {
    const loc = where();
    const abs = sel.regions.map((g) => [reader.absolute(g.chrom, g.start), reader.absolute(g.chrom, g.end)]);
    if (abs.some(([s, e]) => s === null || e === null)) return true;          // hors du génome du modèle : rien à montrer
    const [x, y] = abs.length === 2 ? [abs[0], abs[1]] : [abs[0], abs[0]];
    return x[0]! >= loc.xDomain[0] && x[1]! <= loc.xDomain[1] && y[0]! >= loc.yDomain[0] && y[1]! <= loc.yDomain[1];
  };

  return {
    api,
    reader,
    root,
    shown,
    apply(sel) {
      shown.applied++;
      current = sel;
      if (sel && !contained(sel)) {
        // Une fenêtre autour de chaque intervalle : x sur le premier, y sur le second s'il y en
        // a deux — le pixel du contact au centre. Assez large pour voir son voisinage : dix
        // fois la sélection, au moins 400 kb dans une région fine, 20 Mb ailleurs.
        const abs = sel.regions.map((g) => [reader.absolute(g.chrom, g.start)!, reader.absolute(g.chrom, g.end)!]);
        const inFine = (s: number) => (fin.regions ?? []).some((g) => s >= g.abs_start && s < g.abs_end);
        const win = ([s, e]: number[]) => {
          const w = Math.max(10 * (e - s), inFine(s) ? 400_000 : 20_000_000);
          const c = (s + e) / 2;
          return [c - w / 2, c + w / 2];
        };
        const [x0, x1] = win(abs[0]);
        const [y0, y1] = win(abs[1] ?? abs[0]);
        shown.navigated++;
        settleUntil = performance.now() + ZOOM_MS + 250;
        api.zoomTo(VIEW, x0, x1, y0, y1, ZOOM_MS);
      }
      draw();
      shown.seq = sel?.seq ?? 0;
      shown.at = performance.now();
    },
    box,
    drawn() {
      return Array.from(overlay.querySelectorAll<HTMLElement>(".hg-sel")).map((el) => {
        const b = el.getBoundingClientRect();
        return { region: el.dataset.region ?? "", visible: b.width > 0 && b.height > 0 };
      });
    },
    targets(n, seed) {
      const b = box();
      if (!b) return [];
      const r = root.getBoundingClientRect();
      let st = seed >>> 0;
      const rnd = () => ((st = (st * 1664525 + 1013904223) >>> 0) / 2 ** 32);
      return Array.from({ length: n }, () => [r.left + b.left + 8 + rnd() * (b.width - 16), r.top + b.top + 8 + rnd() * (b.height - 16)] as [number, number]);
    },
    idle: () => GenoContacts.pending === 0 && performance.now() > settleUntil,
    goto(x0, x1) {
      settleUntil = performance.now() + ZOOM_MS + 250;
      api.zoomTo(VIEW, x0, x1, x0, x1, ZOOM_MS);
    },
    resolution,
    stats: GenoContacts.stats,
  };
}
