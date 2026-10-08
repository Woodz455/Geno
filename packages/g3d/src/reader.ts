/** Lecteur `.g3d` : ouvrir avec un préfixe, dessiner, puis ne lire que ce qu'on regarde.
 *
 *  `open` lit un préfixe du fichier — 64 Kio par défaut — et en tire le préambule et
 *  l'en-tête. Si le fichier a été écrit pour ça, le même préfixe contient déjà les
 *  colonnes du premier rendu : `firstFrame` ne fait alors **aucune** requête de plus.
 *  Tout le reste — variabilité, coordonnées génomiques, index, autres niveaux — attend
 *  qu'on le demande.
 */

import {
  type Column,
  type ColumnSpec,
  type Header,
  type LevelSpec,
  PREAMBLE,
  decodeColumn,
  parseHeader,
  parsePreamble,
} from "./format.js";
import type { Source } from "./source.js";

export interface Frame {
  level: LevelSpec;
  count: number;
  /** xyz entrelacés, nm, dans le repère du niveau. L'identifiant de la bille i est i + 1. */
  positions: Float32Array;
  copy: Uint8Array;
}

export interface Bead {
  id: number;
  copy: string;
  start: number;
  end: number;
  variabilityNm: number;
}

export interface Hit {
  copy: string;
  ids: number[];
}

export interface OpenOptions {
  /** Octets lus d'emblée. Doit couvrir le premier rendu pour qu'il tienne en une requête. */
  prefix?: number;
  /** Vérifie l'empreinte de chaque bloc lu. */
  verify?: boolean;
}

export class G3D {
  /** Nombre de lectures émises vers la source depuis l'ouverture. */
  requests = 0;
  private readonly cache: Uint8Array;
  private readonly frames = new Map<string, Frame>();

  private constructor(
    private readonly source: Source,
    readonly header: Header,
    readonly dataStart: number,
    readonly firstFrameEnd: number,
    prefix: Uint8Array,
    private readonly verify: boolean,
  ) {
    this.cache = prefix;
  }

  static async open(source: Source, opts: OpenOptions = {}): Promise<G3D> {
    const prefixLength = opts.prefix ?? 64 * 1024;
    const verify = opts.verify ?? true;
    let prefix = await source.read(0, prefixLength);
    let requests = 1;
    const pre = parsePreamble(prefix);
    const dataStart = PREAMBLE + pre.headerLength;
    if (prefix.length < dataStart) {
      // En-tête plus long que le préfixe. La seconde lecture va **jusqu'à la fin du
      // premier rendu**, que le préambule donne déjà : compléter l'en-tête seul puis
      // relire les colonnes coûterait un aller-retour de plus pour quelques octets.
      // Le test de lecture croisée l'a montré — trois lectures au lieu de deux.
      const until = Math.max(dataStart, pre.firstFrameEnd);
      const rest = await source.read(prefix.length, until - prefix.length);
      requests++;
      const merged = new Uint8Array(prefix.length + rest.length);
      merged.set(prefix);
      merged.set(rest, prefix.length);
      prefix = merged;
    }
    const header = await parseHeader(prefix.subarray(PREAMBLE, dataStart), pre, verify);
    const g = new G3D(source, header, dataStart, pre.firstFrameEnd, prefix, verify);
    g.requests = requests;
    return g;
  }

  level(name: string = this.header.first): LevelSpec {
    const lv = this.header.levels.find((l) => l.name === name);
    if (!lv) throw new Error(`niveau inconnu : ${name}`);
    return lv;
  }

  /** Lit une plage absolue, depuis le préfixe si possible. */
  private async span(offset: number, length: number): Promise<Uint8Array> {
    if (offset + length <= this.cache.length) return this.cache.subarray(offset, offset + length);
    this.requests++;
    return this.source.read(offset, length);
  }

  /** Lit plusieurs colonnes en **une** lecture couvrant leur enveloppe. */
  private async columns(specs: ColumnSpec[]): Promise<Column[]> {
    const lo = Math.min(...specs.map((s) => s.offset));
    const hi = Math.max(...specs.map((s) => s.offset + s.length));
    const bytes = await this.span(this.dataStart + lo, hi - lo);
    return Promise.all(
      specs.map((s) => decodeColumn(bytes.subarray(s.offset - lo, s.offset - lo + s.length), s, this.verify)),
    );
  }

  /** Ce qu'il faut pour dessiner un niveau : positions et copie, rien d'autre. */
  async firstFrame(name: string = this.header.first): Promise<Frame> {
    const cached = this.frames.get(name);
    if (cached) return cached;
    const lv = this.level(name);
    const specs = lv.chunks.flatMap((c) => [c.columns.copy, c.columns.position]);
    const cols = await this.columns(specs);

    const positions = new Float32Array(lv.n_beads * 3);
    const copy = new Uint8Array(lv.n_beads);
    const [ox, oy, oz] = lv.quant.origin;
    const [sx, sy, sz] = lv.quant.step;
    lv.chunks.forEach((c, k) => {
      copy.set(cols[2 * k] as Uint8Array, c.first);
      const q = cols[2 * k + 1] as Uint16Array;            // trois plans : x…, y…, z…
      for (let i = 0; i < c.n; i++) {
        const o = (c.first + i) * 3;
        positions[o] = ox + q[i] * sx;
        positions[o + 1] = oy + q[c.n + i] * sy;
        positions[o + 2] = oz + q[2 * c.n + i] * sz;
      }
    });
    const frame = { level: lv, count: lv.n_beads, positions, copy };
    this.frames.set(name, frame);
    return frame;
  }

  /** Les billes d'un intervalle génomique, **sur toutes les copies** du chromosome.
   *
   *  Une coordonnée correspond à deux positions 3D chez un diploïde (ARCHITECTURE.md
   *  § 6) ; on rend les deux, étiquetées, plutôt que d'en choisir une.
   */
  async locate(name: string, chrom: string, start: number, end: number): Promise<Hit[]> {
    const lv = this.level(name);
    const ix = lv.index.columns;
    const [ids, starts, ends] = (await this.columns([ix.ids, ix.starts, ix.ends])) as Uint32Array[];
    const hits: Hit[] = [];
    lv.copies.forEach((label, c) => {
      if (label.split(":")[0] !== chrom) return;
      const a = lv.index.copy_offsets[c];
      const b = lv.index.copy_offsets[c + 1];
      // Les billes d'une copie ne se chevauchent pas : débuts et fins sont triés tous les deux.
      const lo = lowerBound(ends, start, a, b, (v, q) => v <= q);     // première fin > start
      const hi = lowerBound(starts, end, a, b, (v, q) => v < q);      // premier début ≥ end
      if (hi > lo) hits.push({ copy: label, ids: Array.from(ids.subarray(lo, hi), (i) => i + 1) });
    });
    return hits;
  }

  /** Ce qu'on sait d'une bille cliquée. L'identifiant vient du picking GPU (≥ 1). */
  async identify(name: string, id: number): Promise<Bead> {
    const lv = this.level(name);
    const i = id - 1;
    const c = lv.chunks.find((ch) => i >= ch.first && i < ch.first + ch.n);
    if (!c) throw new Error(`identifiant hors du niveau : ${id}`);
    const [start, end, vari] = await this.columns([c.columns.start, c.columns.end, c.columns.variability]);
    const frame = await this.firstFrame(name);
    const k = i - c.first;
    return {
      id,
      copy: lv.copies[frame.copy[i]],
      start: start[k],
      end: end[k],
      variabilityNm: vari[k] * lv.variability.step,
    };
  }

  /** Chunks dont la sphère englobante coupe le tronc de vue (plans ax+by+cz+d ≥ 0). */
  visibleChunks(name: string, planes: [number, number, number, number][]): number[] {
    return this.level(name).chunks.flatMap((c, k) => {
      const [x, y, z, r] = c.sphere;
      return planes.every(([a, b, cc, d]) => a * x + b * y + cc * z + d >= -r) ? [k] : [];
    });
  }
}

function lowerBound(
  a: Uint32Array,
  q: number,
  lo: number,
  hi: number,
  before: (v: number, q: number) => boolean,
): number {
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (before(a[mid], q)) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

/** Applique la transformation d'un niveau (4×4 en ligne) à des positions xyz entrelacées. */
export function toNucleus(level: LevelSpec, positions: Float32Array): Float32Array {
  const m = level.frame.transform;
  const out = new Float32Array(positions.length);
  for (let i = 0; i < positions.length; i += 3) {
    const x = positions[i], y = positions[i + 1], z = positions[i + 2];
    out[i] = m[0] * x + m[1] * y + m[2] * z + m[3];
    out[i + 1] = m[4] * x + m[5] * y + m[6] * z + m[7];
    out[i + 2] = m[8] * x + m[9] * y + m[10] * z + m[11];
  }
  return out;
}
