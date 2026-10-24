/** Picking GPU.
 *
 *  Le raycasting CPU teste le rayon contre chaque objet : O(n). À six cent mille
 *  billes, un clic coûte plus cher qu'une frame. Le picking GPU laisse le
 *  rasteriseur faire le travail qu'il fait déjà — il a de toute façon déterminé
 *  quel fragment est devant — et se contente de relire **un pixel**.
 *
 *  Le coût est donc constant quel que soit le nombre d'objets. C'est ce qui rend
 *  tenable la fiche d'information de la semaine 13 : cliquer une bille parmi
 *  six cent mille doit coûter la même chose que cliquer parmi dix.
 *
 *  Le prix à payer : `readPixels` synchronise le CPU sur le GPU. Sur un clic
 *  isolé c'est invisible — la semaine 13 le mesure (VALIDATION § S13) ; au survol
 *  continu il faudrait passer par un `PIXEL_PACK_BUFFER` et une lecture asynchrone.
 *
 *  La cible est RGBA32UI, lue en `RGBA_INTEGER` : le seul format entier que WebGL2
 *  garantit pour `readPixels`. R porte l'identifiant, G le type de rendu (`KIND`).
 */

import type { Targets } from "./gl.js";

/** Un pixel relu : identifiant, type de rendu, et profondeur du fragment dans [0, 1]. */
export interface Hit {
  id: number;
  kind: number;
  depth: number;
}

/** Identifiant rendu quand le pixel ne porte aucune bille. */
export const NO_HIT = 0;

export class Picker {
  private readonly scratch = new Uint32Array(4);

  constructor(
    private readonly gl: WebGL2RenderingContext,
    private readonly targets: Targets,
  ) {}

  /** Identifiant de la bille sous (x, y), en pixels, origine en haut à gauche.
   *  Rend `NO_HIT` sur le fond ou hors cadre. */
  pick(x: number, y: number): number {
    return this.read(x, y) ? this.scratch[0] : NO_HIT;
  }

  /** Identifiant, type de rendu et profondeur sous (x, y) — ce que le mapping inverse
   *  3D → génome demande. Deux lectures d'un pixel ; `null` sur le fond ou hors cadre. */
  hit(x: number, y: number): Hit | null {
    if (!this.read(x, y) || this.scratch[0] === NO_HIT) return null;
    return { id: this.scratch[0], kind: this.scratch[1], depth: this.depthAt(x, y) };
  }

  private read(x: number, y: number): boolean {
    const { gl, targets } = this;
    // Le pixel qui **contient** (x, y) : `floor`. `round` (jusqu'à la semaine 12) lisait le
    // voisin dès que x tombait dans la moitié droite d'un pixel — mesuré par l'aller-retour
    // de la semaine 13 : 404 bases exactes sur 500 en `floor`, 252 en `round`.
    const px = Math.floor(x);
    const py = targets.height - 1 - Math.floor(y); // WebGL compte depuis le bas
    if (px < 0 || py < 0 || px >= targets.width || py >= targets.height) return false;
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
    gl.readBuffer(gl.COLOR_ATTACHMENT1);
    gl.readPixels(px, py, 1, 1, gl.RGBA_INTEGER, gl.UNSIGNED_INT, this.scratch);
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    return true;
  }

  /** Profondeur écrite par le fragment shader en (x, y), dans [0, 1].
   *  Un depth buffer ne se relit pas en WebGL2 ; on la double donc dans un
   *  attachement flottant, ce qui rend la profondeur des imposteurs testable. */
  depthAt(x: number, y: number): number {
    const { gl, targets } = this;
    const px = Math.floor(x);
    const py = targets.height - 1 - Math.floor(y);
    const buf = new Float32Array(4);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
    gl.readBuffer(gl.COLOR_ATTACHMENT2);
    gl.readPixels(px, py, 1, 1, gl.RGBA, gl.FLOAT, buf);
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    return buf[0];
  }

  /** Couleur en (x, y), composantes dans [0, 255]. */
  colorAt(x: number, y: number): Uint8Array {
    const { gl, targets } = this;
    const px = Math.floor(x);
    const py = targets.height - 1 - Math.floor(y);
    const buf = new Uint8Array(4);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(px, py, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, buf);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    return buf;
  }
}
