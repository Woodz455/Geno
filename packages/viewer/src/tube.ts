/** Tubes le long d'une polyligne : la fibre de chromatine, les arcs de boucle, l'ADN de liaison.
 *
 *  La géométrie n'est pas construite sur le processeur. Chaque **segment** de la polyligne
 *  est une instance d'un gabarit (u, θ), et le vertex shader évalue lui-même la courbe :
 *  une spline d'Hermite cubique entre deux points, tangentes de Catmull-Rom, et un repère
 *  à **transport parallèle** précalculé une fois — un repère de Frenet tourne brutalement
 *  aux points d'inflexion, celui-ci jamais. Les points, normales et couleurs vivent dans
 *  des textures flottantes ; rien n'est retéléversé quand la caméra bouge.
 *
 *  **Tessellation adaptative.** À chaque image, chaque segment est rangé selon sa longueur
 *  à l'écran dans l'un de cinq gabarits (1 à 16 subdivisions le long du segment), ce qui
 *  coûte quatre petites listes d'instances. Le nombre de côtés du tube, lui, est choisi
 *  **une fois par image pour tout le tube** : si deux segments voisins avaient un nombre de
 *  côtés différent, l'anneau qu'ils partagent n'aurait pas les mêmes sommets des deux côtés
 *  et la surface se fendrait. Le long du segment, l'anneau de bord est le même quel que
 *  soit le gabarit : pas de fissure, et le test le vérifie pixel à pixel.
 *
 *  Plusieurs polylignes tiennent dans un même tube (`breaks`) : aucun segment ne relie
 *  la fin de l'une au début de la suivante. Les arcs de boucle et les liens d'ADN s'en
 *  servent.
 *
 *  Même cible multiple que les imposteurs : couleur, identifiant, profondeur + normale pour
 *  l'occlusion ; même opacité tramée, ici par point (fondu d'un intervalle génomique).
 */

import { link } from "./gl.js";

const TEX_W = 1024;
const SUBDIV = [1, 2, 4, 8, 16] as const;
const SIDES = [4, 6, 8, 12, 16] as const;

// Construit à la demande, comme celui de l'occlusion (ssao.ts) : un gabarit avec `${…}` au
// niveau du module n'est pas éliminé par esbuild, et le tube se retrouvait dans la page du
// premier rendu, qui n'en dessine pas — 2,3 Kio et 22 ms en Slow 4G.
const vertSource = () => /* glsl */ `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;

layout(location = 0) in vec2 aUV;          // u le long du segment, v = θ / 2π
layout(location = 1) in uint aSeg;         // index du premier point du segment

uniform sampler2D uPoint;    // xyz, rayon — négatif : premier point d'une polyligne
uniform sampler2D uNormal;   // normale transportée xyz, opacité
uniform sampler2D uColor;    // rgb première palette, identifiant (flottant exact < 2^24)
uniform sampler2D uColor2;   // rgb seconde palette, 1 : dernier point d'une polyligne
uniform mat4 uProj;
uniform mat4 uView;
uniform float uMix;
uniform uint uIdBase;
uniform int uCount;

out vec3 vNormal;
out vec3 vPos;
flat out uint vId;
flat out vec3 vColor;
flat out float vAlpha;

vec4 fetch(sampler2D t, int k) {
  k = clamp(k, 0, uCount - 1);
  return texelFetch(t, ivec2(k % ${TEX_W}, k / ${TEX_W}), 0);
}

void main() {
  int k = int(aSeg);
  vec4 p0 = fetch(uPoint, k), p1 = fetch(uPoint, k + 1);
  // Tangentes de Catmull-Rom, décentrées aux extrémités d'une polyligne : le premier point
  // n'a pas de voisin avant, le dernier pas de voisin après.
  vec4 pm = fetch(uPoint, k - 1), pp = fetch(uPoint, k + 2);
  bool startK = p0.w < 0.0;
  bool lastK1 = fetch(uColor2, k + 1).w > 0.5;
  vec3 m0 = startK ? (p1.xyz - p0.xyz) : 0.5 * (p1.xyz - pm.xyz);
  vec3 m1 = lastK1 ? (p1.xyz - p0.xyz) : 0.5 * (pp.xyz - p0.xyz);

  float u = aUV.x;
  float u2 = u * u, u3 = u2 * u;
  vec3 pos = (2.0 * u3 - 3.0 * u2 + 1.0) * p0.xyz + (u3 - 2.0 * u2 + u) * m0
           + (-2.0 * u3 + 3.0 * u2) * p1.xyz + (u3 - u2) * m1;
  vec3 tan = normalize((6.0 * u2 - 6.0 * u) * p0.xyz + (3.0 * u2 - 4.0 * u + 1.0) * m0
           + (-6.0 * u2 + 6.0 * u) * p1.xyz + (3.0 * u2 - 2.0 * u) * m1);

  vec4 n0 = fetch(uNormal, k), n1 = fetch(uNormal, k + 1);
  vec3 nrm = mix(n0.xyz, n1.xyz, u);
  nrm = normalize(nrm - dot(nrm, tan) * tan);
  vec3 bin = cross(tan, nrm);
  float th = 6.2831853 * aUV.y;
  vec3 dir = cos(th) * nrm + sin(th) * bin;
  float r = mix(abs(p0.w), abs(p1.w), u);

  // Un pixel appartient au point le plus proche : l'identifiant et l'opacité changent à
  // mi-segment, pas au milieu d'un point.
  bool first = u < 0.5;
  vec4 c = first ? fetch(uColor, k) : fetch(uColor, k + 1);
  vec4 c2 = first ? fetch(uColor2, k) : fetch(uColor2, k + 1);
  vId = uIdBase + uint(c.w);
  vColor = mix(c.rgb, c2.rgb, uMix);
  vAlpha = first ? n0.w : n1.w;

  vec4 world = vec4(pos + r * dir, 1.0);
  vec4 view = uView * world;
  vPos = view.xyz;
  vNormal = mat3(uView) * dir;
  gl_Position = vAlpha <= 0.0 ? vec4(2.0, 2.0, 2.0, 1.0) : uProj * view;
}
`;

const FRAG = /* glsl */ `#version 300 es
precision highp float;
precision highp int;

in vec3 vNormal;
in vec3 vPos;
flat in uint vId;
flat in vec3 vColor;
flat in float vAlpha;

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
  vec3 n = normalize(vNormal);
  if (!gl_FrontFacing) n = -n;
  float lambert = max(dot(n, normalize(uLight)), 0.0);
  oColor = vec4(vColor * (0.25 + 0.75 * lambert), 1.0);
  oId = uvec4(vId, 0u, 0u, 0u);
  oDepth = vec4(gl_FragCoord.z, n);
}
`;

export interface Polyline {
  /** Points xyz entrelacés. */
  points: Float32Array;
  /** Rayon par point. */
  radii: Float32Array;
  colors: Float32Array;
  colors2?: Float32Array;
  /** Identifiant local par point (< 2^24) ; le rendu ajoute `idBase`. */
  ids?: Uint32Array;
  /** Index des points qui commencent une nouvelle polyligne (le premier, 0, est implicite). */
  breaks?: number[];
}

export interface TubeDrawOptions {
  view: Float32Array;
  /** Hauteur du canevas en pixels, pour ranger les segments par taille à l'écran. */
  height: number;
  light?: [number, number, number];
  fade?: number;
  invert?: boolean;
  mix?: number;
  idBase?: number;
  /** **Témoin de test seulement** : un nombre de côtés par gabarit au lieu d'un par image —
   *  la conception qui fend la surface entre deux segments voisins. */
  sidesPerLevel?: boolean;
  /** Élimine les faces arrière. Un tube fermé n'en montre aucune ; une fissure, si. */
  cullBack?: boolean;
}

interface Template {
  vao: WebGLVertexArrayObject;
  instances: WebGLBuffer;
  indexCount: number;
}

/** Les gabarits et le rangement des segments sont exposés pour les tests. */
export function levelFor(px: number): number {
  return px < 3 ? 0 : px < 8 ? 1 : px < 20 ? 2 : px < 48 ? 3 : 4;
}

export function sidesFor(radiusPx: number): number {
  return radiusPx < 2 ? SIDES[0] : radiusPx < 5 ? SIDES[1] : radiusPx < 12 ? SIDES[2] : radiusPx < 30 ? SIDES[3] : SIDES[4];
}

/** Normales à transport parallèle le long d'une polyligne (une par point, unitaires, ⊥ tangente). */
export function transport(points: Float32Array, starts: number[], n: number): Float32Array {
  const out = new Float32Array(n * 3);
  const P = (k: number) => [points[3 * k], points[3 * k + 1], points[3 * k + 2]];
  const sub = (a: number[], b: number[]) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const dot = (a: number[], b: number[]) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm = (a: number[]) => {
    const l = Math.hypot(a[0], a[1], a[2]) || 1;
    return [a[0] / l, a[1] / l, a[2] / l];
  };
  const bounds = [...starts, n];
  for (let s = 0; s < bounds.length - 1; s++) {
    const a = bounds[s];
    const b = bounds[s + 1];
    const tangent = (k: number) => norm(sub(P(Math.min(k + 1, b - 1)), P(Math.max(k - 1, a))));
    const t = tangent(a);
    // Une normale initiale quelconque ⊥ t : l'axe le moins aligné. Puis transport : à chaque
    // point, la normale précédente projetée sur le plan ⊥ à la nouvelle tangente — la plus
    // petite rotation qui la garde perpendiculaire, donc aucune vrille ajoutée.
    const ax = Math.abs(t[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    let nrm = norm(sub(ax, t.map((x) => x * dot(ax, t))));
    for (let k = a; k < b; k++) {
      const tk = tangent(k);
      nrm = norm(sub(nrm, tk.map((x) => x * dot(nrm, tk))));
      out.set(nrm, 3 * k);
    }
  }
  return out;
}

export class TubeRenderer {
  private readonly program: WebGLProgram;
  private readonly templates = new Map<string, Template>();
  private readonly tex: Record<"point" | "normal" | "color" | "color2", WebGLTexture>;
  private readonly normalData: Float32Array[] = [];
  private segments: Uint32Array = new Uint32Array(0);
  private points: Float32Array = new Float32Array(0);
  private radii: Float32Array = new Float32Array(0);
  private count = 0;
  private normals: Float32Array = new Float32Array(0);
  /** Derniers chiffres de rangement, pour les tests et la mesure. */
  lastLevels = [0, 0, 0, 0, 0];
  lastSides = 0;

  constructor(private readonly gl: WebGL2RenderingContext) {
    this.program = link(gl, vertSource(), FRAG);
    const t = () => gl.createTexture()!;
    this.tex = { point: t(), normal: t(), color: t(), color2: t() };
  }

  private texture(tex: WebGLTexture, data: Float32Array): void {
    const gl = this.gl;
    const rows = Math.max(1, Math.ceil(this.count / TEX_W));
    const full = new Float32Array(TEX_W * rows * 4);
    full.set(data.subarray(0, Math.min(data.length, full.length)));
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, TEX_W, rows, 0, gl.RGBA, gl.FLOAT, full);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  }

  upload(line: Polyline): void {
    const n = line.radii.length;
    this.count = n;
    this.points = line.points;
    this.radii = line.radii;
    const starts = [0, ...(line.breaks ?? [])].filter((v, i, a) => a.indexOf(v) === i).sort((a, b) => a - b);
    const isStart = new Uint8Array(n + 1);
    for (const s of starts) isStart[s] = 1;

    const point = new Float32Array(n * 4);
    const color = new Float32Array(n * 4);
    const color2 = new Float32Array(n * 4);
    for (let k = 0; k < n; k++) {
      point.set([line.points[3 * k], line.points[3 * k + 1], line.points[3 * k + 2]], 4 * k);
      // Rayon négatif : premier point d'une polyligne ; color2.w = 1 : dernier point.
      point[4 * k + 3] = isStart[k] === 1 ? -Math.max(line.radii[k], 1e-9) : line.radii[k];
      color.set([line.colors[3 * k], line.colors[3 * k + 1], line.colors[3 * k + 2], line.ids ? line.ids[k] : k + 1], 4 * k);
      const c2 = line.colors2 ?? line.colors;
      const last = k === n - 1 || isStart[k + 1] === 1;
      color2.set([c2[3 * k], c2[3 * k + 1], c2[3 * k + 2], last ? 1 : 0], 4 * k);
    }
    this.normals = transport(line.points, starts, n);
    const normal = new Float32Array(n * 4);
    for (let k = 0; k < n; k++) {
      normal.set([this.normals[3 * k], this.normals[3 * k + 1], this.normals[3 * k + 2], 1], 4 * k);
    }
    this.normalData[0] = normal;
    this.texture(this.tex.point, point);
    this.texture(this.tex.normal, normal);
    this.texture(this.tex.color, color);
    this.texture(this.tex.color2, color2);

    const segs: number[] = [];
    for (let k = 0; k < n - 1; k++) if (!isStart[k + 1]) segs.push(k);
    this.segments = Uint32Array.from(segs);
  }

  /** Opacité par point (fondu d'un intervalle génomique). Une valeur ≤ 0 n'est pas dessinée. */
  setFade(alpha: ArrayLike<number>): void {
    const normal = this.normalData[0];
    if (!normal) return;
    for (let k = 0; k < this.count; k++) normal[4 * k + 3] = alpha[k];
    this.texture(this.tex.normal, normal);
  }

  private template(sub: number, sides: number): Template {
    const key = `${sub}x${sides}`;
    const hit = this.templates.get(key);
    if (hit) return hit;
    const gl = this.gl;
    const uv: number[] = [];
    for (let i = 0; i <= sub; i++) for (let j = 0; j <= sides; j++) uv.push(i / sub, j / sides);
    const idx: number[] = [];
    for (let i = 0; i < sub; i++) {
      for (let j = 0; j < sides; j++) {
        const a = i * (sides + 1) + j;
        const b = a + sides + 1;
        // Sens trigonométrique vu de l'extérieur : ∂θ × ∂u pointe hors du tube. La
        // première version tournait à l'envers, faces avant vers l'intérieur — invisible
        // tant que le tube est dessiné des deux côtés, fatal dès qu'on élimine l'arrière.
        idx.push(a, a + 1, b, a + 1, b + 1, b);
      }
    }
    const vao = gl.createVertexArray()!;
    gl.bindVertexArray(vao);
    const vb = gl.createBuffer()!;
    gl.bindBuffer(gl.ARRAY_BUFFER, vb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(uv), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    const ib = gl.createBuffer()!;
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(idx), gl.STATIC_DRAW);
    const instances = gl.createBuffer()!;
    gl.bindBuffer(gl.ARRAY_BUFFER, instances);
    gl.enableVertexAttribArray(1);
    gl.vertexAttribIPointer(1, 1, gl.UNSIGNED_INT, 0, 0);
    gl.vertexAttribDivisor(1, 1);
    gl.bindVertexArray(null);
    const t = { vao, instances, indexCount: idx.length };
    this.templates.set(key, t);
    return t;
  }

  /** Range les segments par taille à l'écran : longueur projetée et rayon projeté, en pixels. */
  classify(proj: Float32Array, view: Float32Array, height: number): { buckets: number[][]; sides: number } {
    const buckets: number[][] = SUBDIV.map(() => []);
    const f = proj[5] * 0.5 * height;      // pixels par unité à distance 1
    const P = this.points;
    let rMax = 0;
    const depth = (k: number) => -(view[2] * P[3 * k] + view[6] * P[3 * k + 1] + view[10] * P[3 * k + 2] + view[14]);
    const sx = (k: number, z: number) => (view[0] * P[3 * k] + view[4] * P[3 * k + 1] + view[8] * P[3 * k + 2] + view[12]) / z;
    const sy = (k: number, z: number) => (view[1] * P[3 * k] + view[5] * P[3 * k + 1] + view[9] * P[3 * k + 2] + view[13]) / z;
    for (const k of this.segments) {
      const z0 = Math.max(depth(k), 1e-6);
      const z1 = Math.max(depth(k + 1), 1e-6);
      const px = Math.hypot(sx(k, z0) - sx(k + 1, z1), sy(k, z0) - sy(k + 1, z1)) * f;
      buckets[levelFor(px)].push(k);
      rMax = Math.max(rMax, (this.radii[k] / z0) * f);
    }
    return { buckets, sides: sidesFor(rMax) };
  }

  draw(proj: Float32Array, o: TubeDrawOptions): void {
    const gl = this.gl;
    if (!this.segments.length) return;
    const { buckets, sides } = this.classify(proj, o.view, o.height);
    this.lastLevels = buckets.map((b) => b.length);
    this.lastSides = sides;

    gl.useProgram(this.program);
    const loc = (n: string) => gl.getUniformLocation(this.program, n);
    const bind = (unit: number, name: string, tex: WebGLTexture) => {
      gl.activeTexture(gl.TEXTURE0 + unit);
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.uniform1i(loc(name), unit);
    };
    bind(0, "uPoint", this.tex.point);
    bind(1, "uNormal", this.tex.normal);
    bind(2, "uColor", this.tex.color);
    bind(3, "uColor2", this.tex.color2);
    gl.uniformMatrix4fv(loc("uProj"), false, proj);
    gl.uniformMatrix4fv(loc("uView"), false, o.view);
    gl.uniform3fv(loc("uLight"), o.light ?? [0.4, 0.6, 1]);
    gl.uniform1f(loc("uMix"), o.mix ?? 0);
    gl.uniform1f(loc("uFade"), o.fade ?? 1);
    gl.uniform1i(loc("uInvert"), o.invert ? 1 : 0);
    gl.uniform1ui(loc("uIdBase"), o.idBase ?? 0);
    gl.uniform1i(loc("uCount"), this.count);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LESS);
    if (o.cullBack) {
      gl.enable(gl.CULL_FACE);
      gl.cullFace(gl.BACK);
    } else gl.disable(gl.CULL_FACE);

    buckets.forEach((segs, level) => {
      if (!segs.length) return;
      const t = this.template(SUBDIV[level], o.sidesPerLevel ? SIDES[level] : sides);
      gl.bindVertexArray(t.vao);
      gl.bindBuffer(gl.ARRAY_BUFFER, t.instances);
      gl.bufferData(gl.ARRAY_BUFFER, Uint32Array.from(segs), gl.STREAM_DRAW);
      gl.drawElementsInstanced(gl.TRIANGLES, t.indexCount, gl.UNSIGNED_SHORT, 0, segs.length);
    });
    gl.bindVertexArray(null);
    gl.disable(gl.CULL_FACE);
  }

  /** Triangles du dernier dessin : la mesure de la tessellation adaptative. */
  triangles(): number {
    return this.lastLevels.reduce((s, n, l) => s + n * SUBDIV[l] * this.lastSides * 2, 0);
  }

  dispose(): void {
    const gl = this.gl;
    for (const t of this.templates.values()) {
      gl.deleteVertexArray(t.vao);
      gl.deleteBuffer(t.instances);
    }
    for (const t of Object.values(this.tex)) gl.deleteTexture(t);
    gl.deleteProgram(this.program);
  }
}
