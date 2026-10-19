/** La double hélice dans la plongée : la séquence réelle, pour la seule fenêtre visible.
 *
 *  Trois choses, chacune une fonction de la caméra comme le reste de la plongée :
 *
 *  - **La fenêtre.** Les paires de bases dont l'axe tombe à moins de ~¾ de champ du point
 *    visé, dans les monomères déjà générés : quelques centaines à quelques milliers de
 *    paires, jamais la région. Leur géométrie est calculée quand la fenêtre change.
 *  - **Le préchargement.** Dès que la fibre apparaît (champ < 3 µm), les blocs de séquence
 *    des monomères proches du point visé sont demandés, les plus proches d'abord, un à la
 *    fois. À l'arrivée sur l'hélice, la séquence est là depuis longtemps.
 *  - **Le compte de ce qui manque.** Une image qui devrait montrer l'hélice sans avoir sa
 *    séquence, ou des nucléosomes non encore générés, est un **chargement visible** : elle
 *    est comptée, jamais cachée. C'est la mesure du critère de la semaine.
 */

import { type SequenceReader, type TargetSpec } from "@geno/g3d";
import { B_DNA, BASE_COLORS, baseOf, complement, HelixRenderer } from "@geno/viewer";
import { type Fiber, IDS } from "./fiber.js";
import type { State } from "./scene.js";

export interface Label {
  x: number;
  y: number;
  text: string;
  color: [number, number, number];
  size: number;
  alpha: number;
  /** Étiquette de la base visée, sur plusieurs lignes. */
  target?: boolean;
}

export interface LoadingStats {
  /** Images rendues depuis l'ouverture, et celles qui montraient un chargement. */
  frames: number;
  loading: number;
  /** Dont : hélice due sans sa séquence ; nucléosomes dus non générés. */
  sequence: number;
  nucleosomes: number;
}

type V = [number, number, number];
const sub = (a: ArrayLike<number>, b: ArrayLike<number>): V => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const len = (a: ArrayLike<number>) => Math.hypot(a[0], a[1], a[2]);

/** Au-delà, la fenêtre est rognée aux paires les plus proches. */
const MAX_BP = 4000;
const SAMPLE = 8;

export class SequenceLayer {
  readonly helix: HelixRenderer;
  readonly stats: LoadingStats = { frames: 0, loading: 0, sequence: 0, nucleosomes: 0 };
  /** Rang (dans la région) de la base visée, ou −1. */
  readonly targetRank: number;
  private key = "";
  private ready = false;
  private ranges: [number, number][] = [];
  private lastTarget: V = [Infinity, 0, 0];
  private lastWidth = 0;
  private inflight = false;
  private readonly wanted = new Set<number>();

  constructor(
    gl: WebGL2RenderingContext,
    readonly reader: SequenceReader,
    private readonly fiber: Fiber,
    readonly chrom: string,
    /** Coordonnée 0-based de la première base de la région fine. */
    readonly regionStart: number,
    readonly target: TargetSpec | null,
    private readonly prefetchOn = true,
  ) {
    this.helix = new HelixRenderer(gl, 0.001);
    this.targetRank = target ? target.position - 1 - regionStart : -1;
  }

  /** Octets de séquence lus depuis l'ouverture. */
  get bytes(): number {
    return this.reader.bytes;
  }

  private bpm(): number {
    return this.fiber.path.bpPerMonomer;
  }

  /** Préchargement : les blocs des monomères proches du point visé, un à la fois. */
  prefetch(p: number[], widthUm: number): void {
    if (!this.prefetchOn || widthUm > 3 || this.inflight) return;
    const r = Math.max(0.03, Math.min(0.15, widthUm / 8));
    const P = this.fiber.points;
    for (const k of this.fiber.windowNear(p)) {
      const d = Math.hypot(P[3 * k] - p[0], P[3 * k + 1] - p[1], P[3 * k + 2] - p[2]);
      if (d > r) break;
      const g0 = this.regionStart + k * this.bpm();
      if (this.reader.has(this.chrom, g0, g0 + this.bpm())) continue;
      this.inflight = true;
      this.reader.prefetch(this.chrom, g0, g0 + this.bpm()).finally(() => (this.inflight = false)).catch(() => {});
      return;
    }
  }

  /** Lit d'avance la séquence des monomères de la fenêtre à moins de `r` µm de `p`. Pour la
   *  mesure image par image seulement, qui ne rend jamais la main au réseau. */
  async preload(p: number[], r = 0.15): Promise<void> {
    const P = this.fiber.points;
    for (const k of this.fiber.windowNear(p)) {
      if (Math.hypot(P[3 * k] - p[0], P[3 * k + 1] - p[1], P[3 * k + 2] - p[2]) > r) break;
      const g0 = this.regionStart + k * this.bpm();
      await this.reader.prefetch(this.chrom, g0, g0 + this.bpm());
    }
  }

  /** Les plages de paires visibles autour de `p`, à `widthUm` de champ. */
  private windowAt(p: number[], widthUm: number): [number, number][] {
    const path = this.fiber.path;
    const R = 0.75 * widthUm + 0.006;
    const P = this.fiber.points;
    const hits: { b: number; d: number }[] = [];
    for (const k of this.fiber.windowNear(p)) {
      const d = Math.hypot(P[3 * k] - p[0], P[3 * k + 1] - p[1], P[3 * k + 2] - p[2]);
      if (d > R + 0.06) break;
      if (!this.fiber.generated(k)) continue;
      for (let b = k * this.bpm(); b < (k + 1) * this.bpm(); b += SAMPLE) {
        const q = len(sub(path.basePair(b).p, p));
        if (q < R) hits.push({ b, d: q });
      }
    }
    hits.sort((a, b) => a.d - b.d);
    const keep = hits.slice(0, Math.floor(MAX_BP / (2 * SAMPLE))).map((h) => h.b).sort((a, b) => a - b);
    const out: [number, number][] = [];
    for (const b of keep) {
      const lo = Math.max(0, b - SAMPLE), hi = Math.min(path.bpCount, b + SAMPLE);
      const last = out[out.length - 1];
      if (last && lo <= last[0] + last[1]) last[1] = Math.max(last[1], hi - last[0]);
      else out.push([lo, hi - lo]);
    }
    return out;
  }

  /** Une image : prépare la fenêtre si l'hélice est proche, compte ce qui manque. Rend
   *  la part d'hélice à dessiner — nulle si la séquence n'est pas là. */
  update(p: number[], widthUm: number, st: State): number {
    this.stats.frames++;
    let loading = false;
    if (st.nuc > 0 && this.fiber.missing(p, 0.5 * widthUm + 0.05) > 0) {
      this.stats.nucleosomes++;
      loading = true;
    }
    if (widthUm < 0.06) {
      const moved = len(sub(p, this.lastTarget)) > 0.1 * widthUm || Math.abs(Math.log(widthUm / this.lastWidth)) > 0.1;
      if (moved || !this.ready) {
        const ranges = this.windowAt(p, widthUm);
        const key = ranges.map((r) => r.join(":")).join(",");
        if (key !== this.key || !this.ready) {
          this.key = key;
          this.ranges = ranges;
          this.ready = this.upload();
        }
        this.lastTarget = [p[0], p[1], p[2]];
        this.lastWidth = widthUm;
      }
    }
    if (st.helix > 0 && !this.ready) {
      this.stats.sequence++;
      loading = true;
    }
    if (loading) this.stats.loading++;
    return this.ready ? st.helix : 0;
  }

  /** Charge la géométrie de la fenêtre si toute sa séquence est là ; sinon la demande. */
  private upload(): boolean {
    const g = (b: number) => this.regionStart + b;
    const seqs: string[] = [];
    for (const [b0, n] of this.ranges) {
      const s = this.reader.now(this.chrom, g(b0), g(b0 + n));
      if (s === null) {
        for (const [c0, m] of this.ranges) {
          if (!this.wanted.has(c0)) {
            this.wanted.add(c0);
            this.reader.prefetch(this.chrom, g(c0), g(c0 + m)).finally(() => this.wanted.delete(c0)).catch(() => {});
          }
        }
        return false;
      }
      seqs.push(s);
    }
    const total = this.ranges.reduce((s, r) => s + r[1], 0);
    const frames = new Float32Array(total * 9);
    let o = 0;
    for (const [b0, n] of this.ranges) {
      frames.set(this.fiber.path.basePairs(b0, b0 + n), o);
      o += n * 9;
    }
    this.helix.upload({ ranges: this.ranges, frames, bases: seqs.join(""), highlight: this.targetRank });
    return true;
  }

  draw(proj: Float32Array, view: Float32Array, height: number, helix: number): void {
    if (helix > 0 && this.ready) this.helix.draw(proj, { view, height, fade: helix, idBase: IDS.base });
  }

  /** Les lettres des bases assez grandes à l'écran et tournées vers la caméra, et
   *  l'étiquette de la base visée. */
  labels(st: State & { widthUm?: number }, project: (p: number[]) => [number, number], eye: number[], W: number, H: number): Label[] {
    const out: Label[] = [];
    const last = this.helix.last;
    if (!this.ready || !last || st.helix <= 0) return out;
    const { s1, s2, window } = last;
    const widthUm = st.widthUm ?? 0.004;
    // Une tranche de profondeur autour du plan de la base visée : sans elle, les bases de l'ADN
    // loin derrière recevaient leur lettre aussi, et le premier plan devenait illisible.
    let focusDepth = -1;
    const slab = 0.35 * widthUm;
    let q = 0;
    const found: { p: number[] } = { p: [] };
    for (const [b0, n] of window.ranges) {
      if (this.targetRank >= b0 && this.targetRank < b0 + n) {
        const t = q + (this.targetRank - b0);
        found.p = [(s1[3 * t] + s2[3 * t]) / 2, (s1[3 * t + 1] + s2[3 * t + 1]) / 2, (s1[3 * t + 2] + s2[3 * t + 2]) / 2];
        focusDepth = len(sub(eye, found.p));
      }
      q += n;
    }
    q = 0;
    for (const [b0, n] of window.ranges) {
      for (let i = 0; i < n; i++, q++) {
        const rank = b0 + i;
        const a = [s1[3 * q], s1[3 * q + 1], s1[3 * q + 2]];
        const b = [s2[3 * q], s2[3 * q + 1], s2[3 * q + 2]];
        const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
        const base = window.bases[q];
        if (rank === this.targetRank && this.target) {
          const h = [(a[0] + m[0]) / 2, (a[1] + m[1]) / 2, (a[2] + m[2]) / 2];
          const [x, y] = project(h);
          if (x >= 0 && y >= 0 && x <= W && y <= H) out.push({ x, y, text: this.targetText(), color: [1, 1, 1], size: 13, alpha: st.helix, target: true });
        }
        if (st.bases <= 0) continue;
        if (focusDepth >= 0 && Math.abs(len(sub(eye, m)) - focusDepth) > slab) continue;
        for (const [end, letter] of [[a, base], [b, complement(base)]] as const) {
          const h = [(end[0] + m[0]) / 2, (end[1] + m[1]) / 2, (end[2] + m[2]) / 2];
          // Seule la moitié tournée vers la caméra porte sa lettre : l'autre est derrière l'axe.
          if (len(sub(eye, h)) > len(sub(eye, m))) continue;
          const [x, y] = project(h);
          const [x0, y0] = project(end);
          const [x1, y1] = project(m);
          const px = Math.hypot(x1 - x0, y1 - y0);
          if (x < 0 || y < 0 || x > W || y > H || px <= 14) continue;
          out.push({ x, y, text: letter, color: BASE_COLORS[letter] ?? BASE_COLORS.N, size: Math.max(10, Math.min(30, 0.55 * px)), alpha: st.bases });
        }
      }
    }
    return out;
  }

  private targetText(): string {
    const t = this.target!;
    const plus = t.base_plus as string | undefined;
    const coding = t.base_coding as string | undefined;
    return [
      `${t.chrom}:${t.position.toLocaleString("fr")} · ${plus ?? "?"} sur le brin +`,
      `${t.gene} ${t.coding}, codon ${t.codon} · ${coding ?? "?"} sur le brin codant (${t.strand})`,
      `${t.transcript} · ${t.source}`,
    ].join("\n");
  }

  describe(id: number): string | null {
    if (id < IDS.base) return null;
    const { rank, strand } = baseOf(id - IDS.base);
    const pos = this.regionStart + rank;                  // 0-based
    const s = this.reader.now(this.chrom, pos, pos + 1);
    if (s === null) return null;
    const letter = strand === 0 ? s : complement(s);
    const where = `${this.chrom}:${(pos + 1).toLocaleString("fr")}`;
    const tag = rank === this.targetRank && this.target ? ` · ${this.target.gene} ${this.target.coding}, codon ${this.target.codon}` : "";
    return (
      `${where} · ${letter} sur le brin ${strand === 0 ? "+" : "−"} (${complement(letter)} en face)${tag}\n` +
      `séquence hg38 de référence · double hélice B idéale (${B_DNA.bpPerTurn} pb/tour) · position illustrative`
    );
  }
}
