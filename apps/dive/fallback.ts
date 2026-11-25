/** Semaine 22 : le repli des pages qui ouvrent la scène, quand le navigateur n'a pas ce qu'il faut
 *  (`detect`, `packages/viewer/src/support.ts`). Jamais un écran noir ni une exception : un message
 *  qui nomme ce qui manque, et l'aperçu du noyau dessiné en Canvas2D, lu dans le même `.g3d`
 *  (positions de l'aperçu, couleur par chromosome, projection orthographique, billes de l'arrière
 *  vers l'avant). La vue 2D est une image fixe : ni rotation, ni plongée, ni clic.
 */

import { G3D, HttpSource, type Source } from "@geno/g3d";
import type { Support } from "@geno/viewer";
import { territoryPalette } from "./scene.js";

const TEXT = {
  fr: {
    title: "La vue 3D ne peut pas s'ouvrir dans ce navigateur",
    missing: "Il lui manque :",
    still: "Voici l'aperçu du noyau en image fixe : les 46 chromosomes, une couleur chacun. Un navigateur récent (Chrome, Edge, Firefox ou Safari à jour) affichera la vue 3D.",
    alt: "Aperçu fixe du noyau : les billes des 46 chromosomes, vues de face, une couleur par chromosome",
    noData: "L'aperçu n'a pas pu être lu non plus.",
  },
  en: {
    title: "The 3D view cannot open in this browser",
    missing: "It is missing:",
    still: "Here is a still image of the nucleus: the 46 chromosomes, one colour each. An up-to-date browser (Chrome, Edge, Firefox or Safari) will show the 3D view.",
    alt: "Still overview of the nucleus: the beads of the 46 chromosomes, seen from the front, one colour per chromosome",
    noData: "The overview could not be read either.",
  },
};

export async function showFallback(host: HTMLElement, support: Support, file: string | Source, lang: "fr" | "en" = "fr"): Promise<void> {
  const t = TEXT[lang];
  const box = document.createElement("section");
  box.id = "geno-fallback";
  box.setAttribute("role", "alert");
  box.style.cssText = "max-width:720px;margin:24px auto;padding:16px;font:15px/1.5 system-ui,sans-serif;color:#e8edf3;background:#0e1218;border:1px solid #6b7a90;border-radius:8px";
  box.innerHTML = `<h1 style="font-size:20px;margin:0 0 8px">${t.title}</h1><p style="margin:0 0 6px">${t.missing} <strong>${support.missing.join(", ")}</strong>.</p><p style="margin:0 0 12px">${t.still}</p>`;
  host.replaceChildren(box);
  const cv = document.createElement("canvas");
  cv.width = 640;
  cv.height = 640;
  cv.setAttribute("role", "img");
  cv.setAttribute("aria-label", t.alt);
  cv.style.cssText = "width:100%;max-width:640px;height:auto;display:block;background:#0b0f14;border-radius:6px";
  box.appendChild(cv);
  try {
    const g3d = await G3D.open(typeof file === "string" ? new HttpSource(file) : file);
    const frame = await g3d.firstFrame();
    const colors = territoryPalette(frame.level.copies);
    const P = frame.positions;
    let r = 0;
    for (let i = 0; i < frame.count; i++) r = Math.max(r, Math.hypot(P[3 * i], P[3 * i + 1], P[3 * i + 2]));
    const s = (cv.width / 2 - 8) / (r || 1);
    const order = Array.from({ length: frame.count }, (_, i) => i).sort((a, b) => P[3 * a + 2] - P[3 * b + 2]);
    const ctx = cv.getContext("2d")!;
    const rad = Math.max(1.5, (frame.level.radius_nm[0] ?? r / 60) * s);
    for (const i of order) {
      const [cr, cg, cb] = colors[frame.copy[i]].map((c) => Math.round(c * 255));
      const shade = 0.55 + 0.45 * ((P[3 * i + 2] / (r || 1) + 1) / 2);     // plus clair devant
      ctx.fillStyle = `rgb(${Math.round(cr * shade)},${Math.round(cg * shade)},${Math.round(cb * shade)})`;
      ctx.beginPath();
      ctx.arc(cv.width / 2 + P[3 * i] * s, cv.height / 2 - P[3 * i + 1] * s, rad, 0, 2 * Math.PI);
      ctx.fill();
    }
    cv.dataset.beads = String(frame.count);
  } catch {
    const p = document.createElement("p");
    p.textContent = t.noData;
    box.appendChild(p);
  }
}
