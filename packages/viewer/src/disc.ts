/** Imposteurs d'ellipsoïdes de révolution : les nucléosomes.
 *
 *  Un cœur de nucléosome est un disque — 147 pb enroulés autour de l'octamère d'histones,
 *  environ 11 nm de diamètre pour 5,5 nm d'épaisseur. Une sphère le ferait passer pour une
 *  perle ; un ellipsoïde aplati garde l'information qui compte à cette échelle : chaque
 *  nucléosome a une **orientation**.
 *
 *  Même technique que les sphères de la semaine 4 : un quad orienté caméra par instance, et
 *  le fragment shader résout l'intersection rayon–ellipsoïde. On passe le rayon dans le
 *  repère propre de l'ellipsoïde, mis à l'échelle (1/a, 1/a, 1/c) : il y devient une sphère
 *  unité, et l'intersection est celle de la semaine 4. La profondeur écrite est celle du
 *  point d'impact réel, la normale celle de l'ellipsoïde (gradient de l'équation), pas
 *  celle de la sphère déformée.
 */

import { link } from "./gl.js";

const VERT = /* glsl */ `#version 300 es
precision highp float;

layout(location = 0) in vec2 aCorner;
layout(location = 1) in vec3 aCenter;    // repère monde
layout(location = 2) in vec3 aAxis;      // axe de révolution, unitaire, repère monde
layout(location = 3) in vec2 aRadii;     // a (équatorial), c (polaire)
layout(location = 4) in uint aId;
layout(location = 5) in vec3 aColor;
layout(location = 6) in float aAlpha;

uniform mat4 uProj;
uniform mat4 uView;

out vec3 vRay;
flat out vec3 vCenter;
flat out vec3 vAxis;
flat out vec2 vRadii;
flat out uint vId;
flat out vec3 vColor;
flat out float vAlpha;

void main() {
  vec3 c = (uView * vec4(aCenter, 1.0)).xyz;
  vCenter = c;
  vAxis = normalize(mat3(uView) * aAxis);
  vRadii = aRadii;
  vId = aId;
  vColor = aColor;
  vAlpha = aAlpha;
  float r = max(aRadii.x, aRadii.y);            // la sphère englobante
  float d = length(c);
  float k = clamp(r / max(d, 1e-6), 0.0, 0.99);
  float half_ = r / sqrt(1.0 - k * k);
  vec3 pos = c + vec3(aCorner * half_, 0.0);
  vRay = pos;
  gl_Position = aAlpha <= 0.0 ? vec4(2.0, 2.0, 2.0, 1.0) : uProj * vec4(pos, 1.0);
}
`;

const FRAG = /* glsl */ `#version 300 es
precision highp float;
precision highp int;

in vec3 vRay;
flat in vec3 vCenter;
flat in vec3 vAxis;
flat in vec2 vRadii;
flat in uint vId;
flat in vec3 vColor;
flat in float vAlpha;

uniform mat4 uProj;
uniform vec3 uLight;
uniform bool uInvert;
uniform float uFade;

layout(location = 0) out vec4 oColor;
layout(location = 1) out uvec4 oId;
layout(location = 2) out vec4 oDepth;

float threshold(vec2 p) {
  return fract(52.9829189 * fract(dot(p, vec2(0.06711056, 0.00583715))));
}

void main() {
  float a = vAlpha * uFade;
  float h = threshold(gl_FragCoord.xy);
  if (uInvert ? h < 1.0 - a : h >= a) discard;

  // Repère propre : w le long de l'axe, u et v dans le plan équatorial.
  vec3 w = vAxis;
  vec3 t = abs(w.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
  vec3 u = normalize(cross(w, t));
  vec3 v = cross(w, u);
  mat3 toLocal = transpose(mat3(u, v, w));
  vec3 inv = vec3(1.0 / vRadii.x, 1.0 / vRadii.x, 1.0 / vRadii.y);

  // Rayon depuis l'origine (la caméra), passé dans l'espace où l'ellipsoïde est la sphère unité.
  vec3 o = (toLocal * (-vCenter)) * inv;
  vec3 dir = (toLocal * normalize(vRay)) * inv;
  float A = dot(dir, dir);
  float B = dot(o, dir);
  float C = dot(o, o) - 1.0;
  float disc = B * B - A * C;
  if (disc < 0.0) discard;
  float s = (-B - sqrt(disc)) / A;
  if (s < 0.0) discard;

  vec3 hitView = normalize(vRay) * s;
  vec3 q = (toLocal * (hitView - vCenter));                 // point d'impact, repère propre
  vec3 nLocal = q * inv * inv;                              // gradient de x²/a² + y²/a² + z²/c²
  vec3 n = normalize(transpose(toLocal) * nLocal);

  vec4 clip = uProj * vec4(hitView, 1.0);
  gl_FragDepth = (clip.z / clip.w) * 0.5 + 0.5;
  float lambert = max(dot(n, normalize(uLight)), 0.0);
  oColor = vec4(vColor * (0.25 + 0.75 * lambert), 1.0);
  oId = uvec4(vId, 0u, 0u, 0u);
  oDepth = vec4(gl_FragDepth, n);
}
`;

export interface Ellipsoids {
  centers: Float32Array;
  /** Axe de révolution unitaire par instance, xyz. */
  axes: Float32Array;
  /** (a, c) par instance : rayon équatorial, demi-épaisseur. */
  radii: Float32Array;
  ids: Uint32Array;
  colors: Float32Array;
  /** Opacité par instance ; absente : 1. */
  alpha?: Float32Array;
  count: number;
}

export class EllipsoidRenderer {
  private readonly program: WebGLProgram;
  private readonly vao: WebGLVertexArrayObject;
  private readonly buffers: WebGLBuffer[] = [];
  private alphaBuffer: WebGLBuffer | null = null;
  private count = 0;

  constructor(private readonly gl: WebGL2RenderingContext) {
    this.program = link(gl, VERT, FRAG);
    this.vao = gl.createVertexArray()!;
    gl.bindVertexArray(this.vao);
    const cb = gl.createBuffer()!;
    gl.bindBuffer(gl.ARRAY_BUFFER, cb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    this.buffers.push(cb);
    gl.bindVertexArray(null);
  }

  upload(e: Ellipsoids): void {
    const gl = this.gl;
    this.count = e.count;
    gl.bindVertexArray(this.vao);
    const attr = (loc: number, data: ArrayBufferView, size: number, type: number, integer: boolean) => {
      const b = gl.createBuffer()!;
      gl.bindBuffer(gl.ARRAY_BUFFER, b);
      gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW);
      gl.enableVertexAttribArray(loc);
      if (integer) gl.vertexAttribIPointer(loc, size, type, 0, 0);
      else gl.vertexAttribPointer(loc, size, type, false, 0, 0);
      gl.vertexAttribDivisor(loc, 1);
      this.buffers.push(b);
      return b;
    };
    attr(1, e.centers, 3, gl.FLOAT, false);
    attr(2, e.axes, 3, gl.FLOAT, false);
    attr(3, e.radii, 2, gl.FLOAT, false);
    attr(4, e.ids, 1, gl.UNSIGNED_INT, true);
    attr(5, e.colors, 3, gl.FLOAT, false);
    this.alphaBuffer = attr(6, e.alpha ?? new Float32Array(e.count).fill(1), 1, gl.FLOAT, false);
    gl.bindVertexArray(null);
  }

  setAlpha(alpha: Float32Array): void {
    const gl = this.gl;
    if (!this.alphaBuffer) return;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.alphaBuffer);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, alpha);
  }

  draw(proj: Float32Array, o: { view: Float32Array; light?: [number, number, number]; fade?: number; invert?: boolean }): void {
    const gl = this.gl;
    if (!this.count) return;
    gl.useProgram(this.program);
    gl.bindVertexArray(this.vao);
    const loc = (n: string) => gl.getUniformLocation(this.program, n);
    gl.uniformMatrix4fv(loc("uProj"), false, proj);
    gl.uniformMatrix4fv(loc("uView"), false, o.view);
    gl.uniform3fv(loc("uLight"), o.light ?? [0.4, 0.6, 1]);
    gl.uniform1f(loc("uFade"), o.fade ?? 1);
    gl.uniform1i(loc("uInvert"), o.invert ? 1 : 0);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LESS);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, this.count);
    gl.bindVertexArray(null);
  }

  dispose(): void {
    for (const b of this.buffers) this.gl.deleteBuffer(b);
    this.gl.deleteVertexArray(this.vao);
    this.gl.deleteProgram(this.program);
  }
}
