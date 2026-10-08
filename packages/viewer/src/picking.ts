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
 *  isolé c'est invisible ; au survol continu il faudrait passer par un
 *  `PIXEL_PACK_BUFFER` et une lecture asynchrone.
 */

import type { Targets } from "./gl.js";

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
    const { gl, targets } = this;
    const px = Math.round(x);
    const py = Math.round(targets.height - 1 - y); // WebGL compte depuis le bas
    if (px < 0 || py < 0 || px >= targets.width || py >= targets.height) return NO_HIT;

    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
    gl.readBuffer(gl.COLOR_ATTACHMENT1);
    gl.readPixels(px, py, 1, 1, gl.RED_INTEGER, gl.UNSIGNED_INT, this.scratch);
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    return this.scratch[0];
  }

  /** Profondeur écrite par le fragment shader en (x, y), dans [0, 1].
   *  Un depth buffer ne se relit pas en WebGL2 ; on la double donc dans un
   *  attachement flottant, ce qui rend la profondeur des imposteurs testable. */
  depthAt(x: number, y: number): number {
    const { gl, targets } = this;
    const px = Math.round(x);
    const py = Math.round(targets.height - 1 - y);
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
    const px = Math.round(x);
    const py = Math.round(targets.height - 1 - y);
    const buf = new Uint8Array(4);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, targets.fbo);
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(px, py, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, buf);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    return buf;
  }
}
