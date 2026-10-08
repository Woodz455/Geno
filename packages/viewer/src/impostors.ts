/** Imposteurs de sphères.
 *
 *  On ne dessine pas des sphères. On dessine un quad par bille, orienté face à
 *  la caméra, et le fragment shader résout l'intersection rayon–sphère pour
 *  chaque pixel : il jette ce qui tombe hors du disque, calcule la normale, et
 *  **écrit la profondeur du point d'impact réel**.
 *
 *  Cette dernière ligne est tout l'intérêt de la technique. Sans elle, on a des
 *  vignettes plates : deux billes qui s'interpénètrent se découpent selon
 *  l'arête du quad au lieu de la courbe d'intersection, et un nuage dense
 *  devient un collage. Avec elle, la composition est celle de vraies sphères,
 *  pour le coût de deux triangles.
 *
 *  L'enjeu est le budget : une vraie géométrie de sphère, même grossière à 80
 *  triangles, fait 80 millions de triangles pour un million de billes. Un
 *  imposteur en fait deux millions, et le fragment shader ne travaille que sur
 *  les pixels réellement couverts.
 *
 *  Les shaders sont écrits en GLSL ES 3.0 brut et se transposent tels quels
 *  dans un `RawShaderMaterial` de Three.js le jour où le viewer l'utilisera.
 */

import { link, perspective } from "./gl.js";

const VERT = /* glsl */ `#version 300 es
precision highp float;

layout(location = 0) in vec2 aCorner;    // coin du quad, dans [-1, 1]
layout(location = 1) in vec3 aCenter;    // centre de la bille, espace vue
layout(location = 2) in float aRadius;
layout(location = 3) in uint aId;
layout(location = 4) in vec3 aColor;

uniform mat4 uProj;

out vec3 vCenter;
out float vRadius;
out vec3 vRay;        // point du quad, espace vue — origine du rayon
flat out uint vId;
flat out vec3 vColor;

void main() {
  vCenter = aCenter;
  vRadius = aRadius;
  vId = aId;
  vColor = aColor;

  // Sous perspective, la silhouette d'une sphère est une ellipse plus large que
  // son rayon : un quad de demi-taille r laisse dépasser le bord. Le facteur
  // exact à la profondeur du centre est r / sqrt(1 - (r/d)^2).
  float d = length(aCenter);
  float k = clamp(aRadius / max(d, 1e-4), 0.0, 0.99);
  float halfSize = aRadius / sqrt(1.0 - k * k);   // half est un mot réservé en GLSL ES

  vec3 pos = aCenter + vec3(aCorner * halfSize, 0.0);
  vRay = pos;
  gl_Position = uProj * vec4(pos, 1.0);
}
`;

const FRAG = /* glsl */ `#version 300 es
precision highp float;
precision highp int;

in vec3 vCenter;
in float vRadius;
in vec3 vRay;
flat in uint vId;
flat in vec3 vColor;

uniform mat4 uProj;
uniform vec3 uLight;

layout(location = 0) out vec4 oColor;
layout(location = 1) out uvec4 oId;
layout(location = 2) out vec4 oDepth;

void main() {
  // Caméra à l'origine en espace vue : le rayon part de 0 vers vRay.
  vec3 dir = normalize(vRay);
  vec3 oc = -vCenter;                       // origine - centre
  float b = dot(oc, dir);
  float c = dot(oc, oc) - vRadius * vRadius;
  float disc = b * b - c;
  if (disc < 0.0) discard;                  // le pixel est hors du disque

  float t = -b - sqrt(disc);                // racine la plus proche
  if (t < 0.0) discard;                     // sphère derrière la caméra

  vec3 hit = dir * t;
  vec3 normal = normalize(hit - vCenter);

  // LA ligne : la profondeur vient du point d'impact, pas du quad.
  vec4 clip = uProj * vec4(hit, 1.0);
  gl_FragDepth = (clip.z / clip.w) * 0.5 + 0.5;

  float lambert = max(dot(normal, normalize(uLight)), 0.0);
  float ambient = 0.25;
  oColor = vec4(vColor * (ambient + 0.75 * lambert), 1.0);
  oId = uvec4(vId, 0u, 0u, 0u);
  oDepth = vec4(gl_FragDepth, 0.0, 0.0, 1.0);
}
`;

export interface Beads {
  /** Centres en espace vue, xyz entrelacés. */
  centers: Float32Array;
  radii: Float32Array;
  /** 0 est réservé au fond : les identifiants utiles commencent à 1. */
  ids: Uint32Array;
  colors: Float32Array;
  count: number;
}

export class ImpostorRenderer {
  private readonly program: WebGLProgram;
  private readonly vao: WebGLVertexArrayObject;
  private readonly buffers: WebGLBuffer[] = [];
  private readonly uProj: WebGLUniformLocation;
  private readonly uLight: WebGLUniformLocation;
  private count = 0;

  constructor(private readonly gl: WebGL2RenderingContext) {
    this.program = link(gl, VERT, FRAG);
    this.uProj = gl.getUniformLocation(this.program, "uProj")!;
    this.uLight = gl.getUniformLocation(this.program, "uLight")!;

    this.vao = gl.createVertexArray()!;
    gl.bindVertexArray(this.vao);

    // Le quad : quatre coins, en triangle strip. Partagé par toutes les instances.
    const corners = new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]);
    const cb = gl.createBuffer()!;
    gl.bindBuffer(gl.ARRAY_BUFFER, cb);
    gl.bufferData(gl.ARRAY_BUFFER, corners, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    this.buffers.push(cb);

    gl.bindVertexArray(null);
  }

  /** Téléverse les billes. Les attributs sont par instance. */
  upload(beads: Beads): void {
    const gl = this.gl;
    this.count = beads.count;
    gl.bindVertexArray(this.vao);

    const attr = (
      loc: number,
      data: ArrayBufferView,
      size: number,
      type: number,
      integer: boolean,
    ) => {
      const b = gl.createBuffer()!;
      gl.bindBuffer(gl.ARRAY_BUFFER, b);
      gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(loc);
      if (integer) gl.vertexAttribIPointer(loc, size, type, 0, 0);
      else gl.vertexAttribPointer(loc, size, type, false, 0, 0);
      gl.vertexAttribDivisor(loc, 1);
      this.buffers.push(b);
    };

    attr(1, beads.centers, 3, gl.FLOAT, false);
    attr(2, beads.radii, 1, gl.FLOAT, false);
    attr(3, beads.ids, 1, gl.UNSIGNED_INT, true);
    attr(4, beads.colors, 3, gl.FLOAT, false);

    gl.bindVertexArray(null);
  }

  draw(proj: Float32Array, light: [number, number, number] = [0.4, 0.6, 1]): void {
    const gl = this.gl;
    gl.useProgram(this.program);
    gl.bindVertexArray(this.vao);
    gl.uniformMatrix4fv(this.uProj, false, proj);
    gl.uniform3fv(this.uLight, light);
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

export { perspective };
