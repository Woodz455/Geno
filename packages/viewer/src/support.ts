/** Semaine 22 : ce que le navigateur sait faire, mesuré avant d'ouvrir la scène.
 *
 *  Le moteur est WebGL2 seul. Il lui faut :
 *  - **WebGL2** ;
 *  - **EXT_color_buffer_float** : la cible de profondeur RGBA32F que lisent l'occlusion ambiante
 *    et le picking (`gl.ts`). C'est une extension, la seule vraiment optionnelle de la liste ;
 *  - une cible **RGBA32UI** rendable : les identifiants du picking. WebGL2 l'impose ; on la teste
 *    quand même, sur un vrai framebuffer.
 *  WebGPU est seulement noté (`navigator.gpu`) : il n'y a pas de moteur WebGPU. Sans le nécessaire,
 *  les pages montrent un repli (`apps/dive/fallback.ts`) au lieu d'un écran noir.
 */

export interface Support {
  webgl2: boolean;
  colorBufferFloat: boolean;
  uint32Target: boolean;
  maxTextureSize: number;
  webgpu: boolean;
  /** Ce qui manque, en mots, pour le message du repli. */
  missing: string[];
  ok: boolean;
}

export function detect(): Support {
  const out: Support = { webgl2: false, colorBufferFloat: false, uint32Target: false, maxTextureSize: 0,
    webgpu: typeof navigator !== "undefined" && "gpu" in navigator, missing: [], ok: false };
  let gl: WebGL2RenderingContext | null = null;
  try {
    gl = document.createElement("canvas").getContext("webgl2");
  } catch {
    gl = null;
  }
  if (!gl) {
    out.missing.push("WebGL2");
    return out;
  }
  out.webgl2 = true;
  out.maxTextureSize = Number(gl.getParameter(gl.MAX_TEXTURE_SIZE)) || 0;
  out.colorBufferFloat = !!gl.getExtension("EXT_color_buffer_float");
  if (!out.colorBufferFloat) out.missing.push("EXT_color_buffer_float");
  // Une cible RGBA32UI de 4 × 4, attachée : le framebuffer doit être complet.
  const t = gl.createTexture();
  const f = gl.createFramebuffer();
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA32UI, 4, 4);
  gl.bindFramebuffer(gl.FRAMEBUFFER, f);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, t, 0);
  out.uint32Target = gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE;
  if (!out.uint32Target) out.missing.push("rendu en entiers 32 bits (RGBA32UI)");
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.deleteFramebuffer(f);
  gl.deleteTexture(t);
  gl.getExtension("WEBGL_lose_context")?.loseContext();
  out.ok = out.webgl2 && out.colorBufferFloat && out.uint32Target;
  return out;
}
