/** La sélection des vues liées (semaine 14) : une seule, publiée par un store.
 *
 *  Pas de synchronisation deux-à-deux. Chaque vue **publie** ce que l'utilisateur y a
 *  sélectionné et **applique** ce que le store lui donne ; aucune ne parle à une autre.
 *  Deux défenses contre le rebond, chacune comptée dans le journal :
 *
 *  - **égalité** : une sélection égale à la courante (mêmes intervalles) ne change rien et
 *    n'est pas redistribuée ;
 *  - **garde d'application** : pendant qu'une vue applique la sélection reçue, ce qu'elle
 *    publierait est refusé comme écho — une bibliothèque qui émet un événement quand on la
 *    pilote par programme ne peut pas relancer la boucle.
 *
 *  Et une règle de conception, qui rend les deux défenses inutiles en pratique : naviguer
 *  n'est pas sélectionner. Le déplacement d'une vue (igv `locuschange`, HiGlass `location`)
 *  n'est jamais publié — c'est la source connue des boucles infinies, et le témoin de la
 *  mesure (`naive`) la rejoue.
 */

export type ViewName = "3d" | "igv" | "higlass";

/** 0-based, demi-ouvert, haploïde. */
export interface Interval {
  chrom: string;
  start: number;
  end: number;
}

export interface Selection {
  /** Un intervalle (une bille, un gène, un pic) ou deux (un pixel de la carte : deux loci en contact). */
  regions: Interval[];
  source: ViewName;
  /** Ce que l'utilisateur a désigné, en clair : « bille de 750 kb », « gène TP53 »… */
  label: string;
  /** Numéro d'ordre : une vue sait si elle a déjà appliqué celle-ci. */
  seq: number;
  /** performance.now() de la publication. */
  at: number;
}

export interface LogEntry {
  t: number;
  source: ViewName;
  label: string;
  accepted: boolean;
  /** Pourquoi refusée : « égale », « écho » ; vide si acceptée. */
  reason: string;
  seq: number;
}

type Listener = (sel: Selection | null) => void;

export const key = (regions: Interval[]): string =>
  regions.map((r) => `${r.chrom}:${r.start}-${r.end}`).sort().join(" ");

export class SelectionStore {
  current: Selection | null = null;
  readonly log: LogEntry[] = [];
  private readonly listeners = new Map<ViewName, Listener>();
  private readonly observers: Listener[] = [];
  private applying: ViewName | null = null;
  private seq = 0;

  /** Une vue par nom : la remplacer est permis, en avoir deux non. */
  subscribe(view: ViewName, fn: Listener): void {
    this.listeners.set(view, fn);
  }

  /** Ce qui n'est pas une vue (la barre d'état) : appelé à chaque sélection acceptée, source
   *  comprise, et ne publie jamais. */
  observe(fn: Listener): void {
    this.observers.push(fn);
  }

  /** Ce que l'utilisateur a sélectionné dans `source`. Rend vrai si la sélection a changé. */
  publish(source: ViewName, regions: Interval[] | null, label = ""): boolean {
    const t = performance.now();
    const norm = (regions ?? []).map((r) => ({ chrom: r.chrom, start: Math.max(0, Math.floor(r.start)), end: Math.ceil(r.end) }))
      .filter((r) => r.end > r.start);
    const refuse = (reason: string) => {
      this.log.push({ t, source, label, accepted: false, reason, seq: this.seq });
      return false;
    };
    if (this.applying === source) return refuse("écho");
    const same = norm.length === 0 ? this.current === null : this.current !== null && key(this.current.regions) === key(norm);
    if (same) return refuse("égale");
    this.seq++;
    this.current = norm.length ? { regions: norm, source, label, seq: this.seq, at: t } : null;
    this.log.push({ t, source, label, accepted: true, reason: "", seq: this.seq });
    for (const [view, fn] of this.listeners) {
      if (view === source) continue;
      this.applying = view;
      try {
        fn(this.current);
      } finally {
        this.applying = null;
      }
    }
    for (const fn of this.observers) fn(this.current);
    return true;
  }
}
