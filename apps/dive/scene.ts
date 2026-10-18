/** Ce que montre la plongée, en fonction de la caméra seulement.
 *
 *  Règle de conception : **l'état visuel est une fonction de la caméra, pas du temps.**
 *  Niveau de détail, extinction des autres chromosomes, passage des territoires aux
 *  compartiments : tout se lit sur la largeur de vue et la taille des billes à l'écran.
 *  Une plongée est donc réversible — remonter défait exactement ce que descendre a fait —
 *  et un arrêt en chemin montre un état intermédiaire stable, jamais une animation figée.
 *  La seule exception est l'arrivée d'un niveau par le réseau, qui ne dépend pas de la
 *  caméra ; elle est fondue sur 400 ms, et la mesure de continuité la neutralise.
 */

export interface View {
  /** Largeur visible au point visé, µm. */
  widthUm: number;
  /** Rayon à l'écran d'une bille du niveau fin (750 kb), pixels. */
  beadPx: number;
}

export interface State {
  /** 0 : aperçu à 3 Mb ; 1 : noyau à 750 kb. Fondu tramé entre les deux. */
  lod: number;
  /** 0 : tous les chromosomes ; 1 : seule la copie visée. */
  focus: number;
  /** 0 : couleurs de territoire ; 1 : couleurs de compartiment. */
  mix: number;
  /** 0 → 1 : la fibre de 2 kb remplace les billes de 750 kb (semaine 11). Le saut S8 est ici. */
  fiber: number;
  /** 0 → 1 : la boucle visée s'isole ; sites CTCF, cohésines et arcs apparaissent. */
  loop: number;
  /** 0 → 1 : le tube cède aux nucléosomes et à l'ADN de liaison. */
  nuc: number;
  /** 0 → 1 : l'ADN de 1 nm cède à la double hélice, sur la séquence réelle (semaine 12). */
  helix: number;
  /** 0 → 1 : les lettres des bases paraissent. */
  bases: number;
}

/** Les seuils, en unités qui ont un sens : des pixels pour le détail, des µm pour le sujet. */
export const SCHEDULE = {
  // L'aperçu cède la place quand une bille du noyau dépasse 3 px de rayon, et a disparu à 6.
  lodPx: [3, 6] as const,
  // Les autres chromosomes s'effacent entre 20 et 8 µm de champ : avant que la caméra entre
  // dans le noyau (rayon 5 µm), sinon elle traverserait des billes opaques.
  focusUm: [20, 8] as const,
  // Puis les territoires laissent place aux compartiments, entre 8 et 3,5 µm.
  mixUm: [8, 3.5] as const,
  // Semaine 11. La fibre remplace les billes de 750 kb entre 3 et 1,5 µm : une bille y fait
  // 30 à 60 px, et c'est en dessous qu'elle cesse de dire quoi que ce soit.
  fiberUm: [3, 1.5] as const,
  // Une boucle de 100 à 200 kb occupe quelques centaines de nm : elle « existe » à l'écran
  // entre 1,4 µm et 600 nm de champ, et c'est là que ses annotations apparaissent.
  loopUm: [1.4, 0.6] as const,
  // Sous 250 nm, un monomère de 2 kb fait plus de 100 px : ses ~10 nucléosomes se voient.
  nucUm: [0.25, 0.12] as const,
  // Semaine 12. Sous 40 nm de champ, le diamètre de 2 nm de l'ADN dépasse 25 px sur un écran
  // de 540 px : les deux brins et leurs barreaux se séparent. À 20 nm, l'hélice est seule.
  helixUm: [0.04, 0.02] as const,
  // Sous 8 nm, un barreau fait plus de 60 px : une lettre y tient. À 5 nm, elles sont toutes là.
  basesUm: [0.008, 0.005] as const,
};

export function smoothstep(a: number, b: number, x: number): number {
  const t = Math.max(0, Math.min(1, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
}

/** Les transitions sont douces sur une plage. `hard` les remplace par des sauts à mi-plage :
 *  c'est le **témoin** de la mesure de continuité, qui doit, lui, produire des sauts. */
export function stateAt(v: View, hard = false): State {
  const step = (a: number, b: number, x: number) => (hard ? (x >= (a + b) / 2 ? 1 : 0) : smoothstep(a, b, x));
  // Sur la largeur, on raisonne en logarithme : 20 → 8 µm est un facteur, pas une différence.
  const lw = Math.log(v.widthUm);
  const [f0, f1] = SCHEDULE.focusUm.map(Math.log);
  const [m0, m1] = SCHEDULE.mixUm.map(Math.log);
  const zoomIn = (r: readonly [number, number]) => step(-Math.log(r[0]), -Math.log(r[1]), -lw);
  return {
    lod: step(SCHEDULE.lodPx[0], SCHEDULE.lodPx[1], v.beadPx),
    focus: step(-f0, -f1, -lw),
    mix: step(-m0, -m1, -lw),
    fiber: zoomIn(SCHEDULE.fiberUm),
    loop: zoomIn(SCHEDULE.loopUm),
    nuc: zoomIn(SCHEDULE.nucUm),
    helix: zoomIn(SCHEDULE.helixUm),
    bases: zoomIn(SCHEDULE.basesUm),
  };
}

// ---------------------------------------------------------------------------
// Palettes
// ---------------------------------------------------------------------------

function hsl(h: number, s: number, l: number): [number, number, number] {
  const a = s * Math.min(l, 1 - l);
  const f = (n: number) => {
    const k = (n + h / 30) % 12;
    return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1));
  };
  return [f(0), f(8), f(4)];
}

/** Couleur de territoire par copie : une teinte par chromosome (angle d'or, pour que deux
 *  numéros voisins ne se ressemblent pas), l'homologue « b » plus clair, l'Xi éteint — il
 *  est condensé en corps de Barr, et le dire par la couleur est une information. */
export function territoryPalette(copies: string[]): [number, number, number][] {
  return copies.map((label) => {
    const [chrom, hap] = label.split(":");
    const name = chrom.replace(/^chr/, "");
    const k = name === "X" ? 22 : name === "Y" ? 23 : Number(name) - 1;
    const hue = (k * 137.508) % 360;
    if (hap === "Xi") return hsl(hue, 0.12, 0.5);
    return hsl(hue, 0.55, hap === "b" ? 0.66 : 0.52);
  });
}

export const COMPARTMENT = {
  A: [0.88, 0.42, 0.2] as [number, number, number],
  B: [0.22, 0.42, 0.86] as [number, number, number],
  absent: [0.5, 0.52, 0.55] as [number, number, number],
};

/** Couleur de compartiment par bille : le signe de E1 choisit A ou B, sa valeur absolue
 *  la saturation — une bille à E1 ≈ 0 est entre les deux, et doit le paraître. */
export function compartmentColors(e1: Float32Array | null, n: number): Float32Array {
  const out = new Float32Array(n * 3);
  if (!e1) {
    for (let i = 0; i < n; i++) out.set(COMPARTMENT.absent, i * 3);
    return out;
  }
  const mags = Array.from(e1).filter(Number.isFinite).map(Math.abs).sort((a, b) => a - b);
  const q95 = mags.length ? mags[Math.floor(mags.length * 0.95)] || 1 : 1;
  for (let i = 0; i < n; i++) {
    const v = e1[i];
    if (!Number.isFinite(v)) {
      out.set(COMPARTMENT.absent, i * 3);
      continue;
    }
    const base = v > 0 ? COMPARTMENT.A : COMPARTMENT.B;
    const k = 0.35 + 0.65 * Math.min(1, Math.abs(v) / q95);
    for (let c = 0; c < 3; c++) out[i * 3 + c] = COMPARTMENT.absent[c] + (base[c] - COMPARTMENT.absent[c]) * k;
  }
  return out;
}
