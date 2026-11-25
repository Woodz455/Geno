/** Semaine 22 : la porte des pages qui ouvrent la scène. Elle mesure ce que le navigateur sait faire
 *  (`detect`) ; s'il manque quelque chose, elle montre le repli (`fallback.ts`) et la page ne se
 *  charge pas ; sinon, elle charge la page. Le bundle reste un seul fichier : esbuild n'initialise
 *  le module de la page qu'à l'appel de `load` (import dynamique sans découpage). */

import { detect } from "@geno/viewer";
import { showFallback } from "./fallback.js";

export function gate(load: () => Promise<unknown>): void {
  const support = detect();
  const w = window as unknown as Record<string, unknown>;
  w.__support = support;
  if (support.ok) {
    void load();
    return;
  }
  const q = new URLSearchParams(location.search);
  const f = q.get("f") ?? "";
  const lang = q.get("lang") === "en" || (!q.get("lang") && navigator.language.toLowerCase().startsWith("en")) ? "en" : "fr";
  // Pour les mesures : chaque page publie son état ; le repli publie le sien.
  const state = { ok: false, fallback: true, missing: support.missing };
  w.__geno = state;
  w.__discover = { ready: true, ...state };
  void showFallback(document.body, support, f.endsWith(".g3d") ? f : "data/gm12878.g3d", lang).then(() => {
    w.__fallbackDrawn = Number((document.querySelector("#geno-fallback canvas") as HTMLCanvasElement | null)?.dataset.beads ?? 0);
  });
}
