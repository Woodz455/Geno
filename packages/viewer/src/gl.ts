/** Aides WebGL2 minimales. Pas un framework : juste ce qui évite de répéter. */

export function compile(gl: WebGL2RenderingContext, type: number, source: string): WebGLShader {
  const sh = gl.createShader(type);
  if (!sh) throw new Error("createShader a rendu null");
  gl.shaderSource(sh, source);
  gl.compileShader(sh);
  if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
    const log = gl.getShaderInfoLog(sh) ?? "(pas de log)";
    const kind = type === gl.VERTEX_SHADER ? "vertex" : "fragment";
    gl.deleteShader(sh);
    throw new Error(`shader ${kind} : ${log}`);
  }
  return sh;
}

export function link(gl: WebGL2RenderingContext, vs: string, fs: string): WebGLProgram {
  const p = gl.createProgram();
  if (!p) throw new Error("createProgram a rendu null");
  const v = compile(gl, gl.VERTEX_SHADER, vs);
  const f = compile(gl, gl.FRAGMENT_SHADER, fs);
  gl.attachShader(p, v);
  gl.attachShader(p, f);
  gl.linkProgram(p);
  gl.deleteShader(v);
  gl.deleteShader(f);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
    const log = gl.getProgramInfoLog(p) ?? "(pas de log)";
    gl.deleteProgram(p);
    throw new Error(`édition de liens : ${log}`);
  }
  return p;
}

/** Projection perspective, colonne-majeure comme l'attend WebGL. */
export function perspective(fovY: number, aspect: number, near: number, far: number): Float32Array {
  const f = 1 / Math.tan(fovY / 2);
  const nf = 1 / (near - far);
  // prettier-ignore
  return new Float32Array([
    f / aspect, 0, 0,                   0,
    0,          f, 0,                   0,
    0,          0, (far + near) * nf,  -1,
    0,          0, 2 * far * near * nf, 0,
  ]);
}

export interface Targets {
  fbo: WebGLFramebuffer;
  /** Couleur visible, RGBA8. */
  color: WebGLTexture;
  /** Identifiant d'instance, R32UI — c'est ce que lit le picking. */
  id: WebGLTexture;
  /** Profondeur écrite par le fragment shader, RGBA32F. Sert aux tests : on ne
   *  peut pas relire un depth buffer en WebGL2, mais on peut relire un float. */
  depth: WebGLTexture;
  depthBuffer: WebGLRenderbuffer;
  width: number;
  height: number;
}

/** Cible de rendu multiple : couleur + identifiants + profondeur lisible.
 *
 *  Les trois attachements sont écrits en une seule passe. Une seconde passe
 *  dédiée au picking doublerait le coût géométrique, ce qui à un million
 *  d'instances n'est pas une nuance. */
export function createTargets(gl: WebGL2RenderingContext, width: number, height: number): Targets {
  // RGBA32F n'est pas rendable en WebGL2 sans cette extension.
  if (!gl.getExtension("EXT_color_buffer_float")) {
    throw new Error("EXT_color_buffer_float absent : la cible de profondeur est impossible");
  }
  const tex = (internal: number, format: number, type: number) => {
    const t = gl.createTexture()!;
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texStorage2D(gl.TEXTURE_2D, 1, internal, width, height);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    void format;
    void type;
    return t;
  };

  const color = tex(gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE);
  const id = tex(gl.R32UI, gl.RED_INTEGER, gl.UNSIGNED_INT);
  const depth = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT);

  const depthBuffer = gl.createRenderbuffer()!;
  gl.bindRenderbuffer(gl.RENDERBUFFER, depthBuffer);
  gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT24, width, height);

  const fbo = gl.createFramebuffer()!;
  gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, color, 0);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, id, 0);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, depth, 0);
  gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, depthBuffer);
  gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2]);

  const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
  if (status !== gl.FRAMEBUFFER_COMPLETE) {
    throw new Error(`framebuffer incomplet : 0x${status.toString(16)}`);
  }
  return { fbo, color, id, depth, depthBuffer, width, height };
}
