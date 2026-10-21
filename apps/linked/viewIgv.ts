/** Le navigateur linéaire des vues liées : igv.js, sur la séquence et les annotations
 *  verrouillées (`make linked`), sans réseau extérieur.
 *
 *  Sélectionner, c'est cliquer un élément d'une piste — un gène, un pic CTCF, un bin de
 *  compartiment, une boucle (ses deux ancres). Naviguer (glisser, zoomer, taper un locus)
 *  n'est **pas** sélectionner et n'est jamais publié. Une sélection reçue devient une
 *  région d'intérêt (ROI), et la vue s'y déplace si elle n'y est pas.
 */

import igv from "igv/dist/igv.esm.js";
import type { Interval, Selection } from "./store.js";

export interface IgvIndex {
  igv: Record<string, { path: string; bytes: number; sha256: string }>;
  evidence: Record<string, string>;
}

/** Ce que la page lit d'igv : un sous-ensemble de son API publique, typé ici. */
interface Frame {
  chr: string;
  start: number;
  end: number;
  bpPerPixel: number;
}
export interface IgvBrowser {
  referenceFrameList: Frame[];
  on(event: string, fn: (...a: any[]) => unknown): void;
  search(locus: string): Promise<unknown>;
  loadROI(roi: unknown): Promise<void>;
  clearROIs(): void;
  findTracks(fn: (t: any) => boolean): any[];
}

export interface ViewIgv {
  browser: IgvBrowser;
  root: HTMLElement;
  shown: { seq: number; at: number; applied: number; navigated: number };
  apply(sel: Selection | null): void;
  /** Les ROI dessinées : leur région (`chr:start-end`), et si elles tiennent dans la vue. */
  rois(): { region: string; visible: boolean }[];
  /** Pour la mesure : des points de page sur des éléments de la piste des gènes (ou d'une
   *  autre), avec l'élément qu'ils désignent. */
  targets(track: string, n: number): { x: number; y: number; feature: Interval & { name: string } }[];
  /** Plus rien en cours : pas de recherche ni de chargement de piste. */
  idle(): boolean;
  clicks: { track: string; type: string; features: number; data: number }[];
  goto(locus: string): Promise<void>;
}

const ROI_COLOR = "rgba(240, 228, 66, 0.45)";
const roiOf = (regions: Interval[]) => [{ name: "sélection", color: ROI_COLOR,
  features: regions.map((r) => ({ chr: r.chrom, start: r.start, end: r.end })) }];
const fmt = (n: number) => n.toLocaleString("en");

export async function openIgv(root: HTMLElement, index: IgvIndex, publish: (regions: Interval[], label: string) => void,
                              locus = "chr17:7,560,000-7,790,000"): Promise<ViewIgv> {
  const f = (k: string) => `data/igv/${index.igv[k].path}`;
  const browser = (await igv.createBrowser(root, {
    reference: {
      id: "geno-hg38", name: "GRCh38 (hg38.2bit verrouillé)",
      twoBitURL: "ref/hg38.2bit", chromSizesURL: f("chromsizes"), cytobandURL: f("cytoband"),
    },
    locus,
    showSampleNames: false,
    tracks: [
      { name: "Gènes · GENCODE v50 · annotation", type: "annotation", format: "bed", url: f("genes"), displayMode: "COLLAPSED", color: "#3a5e8c" },
      { name: "CTCF GM12878 · ENCODE · mesuré (brin = motif)", type: "annotation", format: "bed", url: f("ctcf"), displayMode: "COLLAPSED", color: "#7a4d9e" },
      { name: "Compartiments · E1 du modèle · simulé", type: "wig", format: "bedgraph", url: f("compartments"), color: "rgb(224,107,51)", altColor: "rgb(56,107,219)", height: 40 },
      { name: "Boucles · Rao 2014 · mesuré (hg19 → hg38)", type: "interact", format: "bedpe", url: f("loops"), height: 60, color: "rgb(0,114,178)" },
    ],
  })) as IgvBrowser;

  const shown: ViewIgv["shown"] = { seq: 0, at: 0, applied: 0, navigated: 0 };
  /** Chaque clic de piste qu'igv a signalé, sélection ou non : le diagnostic d'un clic perdu. */
  const clicks: { track: string; type: string; features: number; data: number }[] = [];
  let busy = 0;

  browser.on("trackclick", (track: any, data: any, _pos: number, features?: any[]) => {
    clicks.push({ track: String(track.name), type: String(track.type), features: features?.length ?? -1, data: Array.isArray(data) ? data.length : -1 });
    let regions: Interval[] = [];
    let label = "";
    const feat = features?.[0];
    if (feat && feat.chr !== undefined) {
      regions = [{ chrom: feat.chr, start: feat.start, end: feat.end }];
      label = `${String(track.name).split(" · ")[0].toLowerCase()} ${feat.name ?? ""}`.trim();
    } else if (track.type === "interact" && Array.isArray(data)) {
      // Une boucle : ses deux ancres. igv ne passe pas l'objet, on relit la bulle.
      const v = (k: string) => data.find((d: any) => d.name === k)?.value;
      const a = { chrom: v("chr1"), start: Number(v("start1")), end: Number(v("end1")) };
      const b = { chrom: v("chr2"), start: Number(v("start2")), end: Number(v("end2")) };
      if (a.chrom && b.chrom && Number.isFinite(a.start) && Number.isFinite(b.start)) {
        regions = [a, b];
        label = "boucle (deux ancres)";
      }
    } else if (track.type === "wig" && Array.isArray(data)) {
      const pos = data.find((d: any) => d.name === "Position" || d.name === "position")?.value;
      const m = typeof pos === "string" ? /^(chr\w+):([\d,]+)-([\d,]+)/.exec(pos) : null;
      if (m) {
        regions = [{ chrom: m[1], start: Number(m[2].replace(/,/g, "")) - 1, end: Number(m[3].replace(/,/g, "")) }];
        label = "bin de compartiment (E1 du modèle)";
      }
    }
    if (!regions.length) return undefined;           // pas un élément : la bulle par défaut
    // La source montre sa propre sélection : le store ne la lui renvoie pas.
    browser.clearROIs();
    void browser.loadROI(roiOf(regions));
    publish(regions, label);
    return false;                                    // sélectionné : pas de bulle par-dessus
  });

  const contains = (r: Interval) => browser.referenceFrameList.some((fr) => fr.chr === r.chrom && fr.start <= r.start && fr.end >= r.end);
  const view: ViewIgv = {
    browser,
    root,
    shown,
    apply(sel) {
      shown.applied++;
      browser.clearROIs();
      if (!sel) {
        shown.seq = 0;
        shown.at = performance.now();
        return;
      }
      busy++;
      const done = () => {
        busy--;
        shown.seq = sel.seq;
        shown.at = performance.now();
      };
      const roi = browser.loadROI(roiOf(sel.regions));
      if (sel.regions.every(contains)) {
        void roi.then(done, done);
        return;
      }
      // Se déplacer : une fenêtre par intervalle (igv affiche plusieurs loci côte à côte),
      // au moins 20 kb, trois fois la sélection.
      shown.navigated++;
      const loci = sel.regions.map((r) => {
        const w = Math.max(20_000, 3 * (r.end - r.start));
        const c = (r.start + r.end) / 2;
        return `${r.chrom}:${fmt(Math.max(1, Math.round(c - w / 2)))}-${fmt(Math.round(c + w / 2))}`;
      });
      void Promise.all([roi, browser.search(loci.join(" "))]).then(done, done);
    },
    rois() {
      const box = root.getBoundingClientRect();
      // igv.js 3 dessine dans un shadow DOM attaché au conteneur.
      return Array.from((root.shadowRoot ?? root).querySelectorAll<HTMLElement>(".igv-roi-region")).map((el) => {
        const b = el.getBoundingClientRect();
        return { region: el.dataset.region ?? "", visible: b.width > 0 && b.right > box.left && b.left < box.right };
      });
    },
    targets(name, n) {
      const t = browser.findTracks((tr) => String(tr.name).startsWith(name))[0];
      if (!t) return [];
      const fr = browser.referenceFrameList[0];
      // La piste dans la page : son élément de vue (une colonne par locus ; la première ici).
      const port = t.trackView?.viewports?.[0];
      const vp = port?.viewportElement as HTMLElement | undefined;
      if (!vp) return [];
      const box = vp.getBoundingClientRect();
      // Les éléments que la vue a dessinés (le cache de la fenêtre, interne à igv).
      const feats: any[] = (port.featureCache?.features ?? [])
        .filter((g: any) => g.chr === fr.chr && g.end > fr.start && g.start < fr.end);
      const out: ReturnType<ViewIgv["targets"]> = [];
      for (const g of feats) {
        const a = Math.max(g.start, fr.start), b = Math.min(g.end, fr.end);
        const mid = (a + b) / 2;
        const x = box.left + (mid - fr.start) / fr.bpPerPixel;
        if (x < box.left + 2 || x > box.right - 2) continue;
        // Un seul élément sous le pointeur (à trois pixels près) : sinon le clic en désigne un
        // autre, et la mesure ne saurait pas lequel attendre.
        const slack = 3 * fr.bpPerPixel;
        if (feats.some((o) => o !== g && o.start < mid + slack && o.end > mid - slack)) continue;
        // Et rien par-dessus : le cartouche du nom de la piste, en haut à gauche, prend le clic.
        const y = box.top + 12;
        const top = (root.shadowRoot ?? document).elementFromPoint(x, y);
        if (!top || !vp.contains(top) || top.closest(".igv-track-label")) continue;
        out.push({ x, y, feature: { chrom: g.chr, start: g.start, end: g.end, name: g.name ?? "" } });
        if (out.length >= n) break;
      }
      return out;
    },
    idle: () => busy === 0,
    clicks,
    async goto(l) {
      busy++;
      try {
        await browser.search(l);
      } finally {
        busy--;
      }
    },
  };
  return view;
}
