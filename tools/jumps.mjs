/** Le détecteur de sauts des bancs de continuité (semaines 10 à 12).
 *
 *  Un **saut** est une image dont l'écart à la précédente dépasse 3 fois la médiane de ses
 *  20 voisines (et 0,5/255, pour ne pas compter le bruit d'une image presque immobile).
 */

export const RATIO = 3;
export const WINDOW = 10;

const median = (xs) => {
  const s = xs.slice().sort((a, b) => a - b);
  return s.length ? s[s.length >> 1] : NaN;
};

/** Rapport de chaque écart à la médiane de ses voisines (la valeur elle-même exclue). */
export function ratios(diffs) {
  return diffs.map((d, k) => {
    const around = [];
    for (let j = Math.max(0, k - WINDOW); j <= Math.min(diffs.length - 1, k + WINDOW); j++) if (j !== k) around.push(diffs[j]);
    return d / Math.max(median(around), 1e-3);
  });
}

/** Images dont l'écart à la précédente dépasse RATIO × la médiane de leurs voisines. */
export function jumps(diffs) {
  return ratios(diffs).flatMap((ratio, k) =>
    ratio > RATIO && diffs[k] > 0.5 ? [{ frame: k + 1, diff: diffs[k], ratio }] : []);
}
