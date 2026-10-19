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
 *
 *  Semaine 10 : les centres sont en repère monde et passent par `uView` (identité
 *  par défaut, ce qui laisse la suite de la semaine 4 inchangée) ; chaque bille porte
 *  un **groupe** — sa copie de chromosome — dont l'opacité se règle d'un uniforme ; et
 *  l'opacité est rendue par **tramage** (voir `FRAG`), pas par mélange.
 */

import { KIND, link, perspective } from "./gl.js";

const VERT = /* glsl */ `#version 300 es
precision highp float;

layout(location = 0) in vec2 aCorner;    // coin du quad, dans [-1, 1]
layout(location = 1) in vec3 aCenter;    // centre de la bille, espace vue
layout(location = 2) in float aRadius;
layout(location = 3) in uint aId;
layout(location = 4) in vec3 aColor;
layout(location = 5) in uint aGroup;     // copie de chromosome ; 0 si non fourni
layout(location = 6) in vec3 aColor2;    // seconde palette, mélangée par uMix

uniform mat4 uProj;
uniform mat4 uView;
uniform float uFade;                      // opacité du niveau entier
uniform float uGroupFade[MAX_GROUPS];     // opacité par groupe
uniform float uMix;

out vec3 vCenter;
out float vRadius;
out vec3 vRay;        // point du quad, espace vue — origine du rayon
flat out uint vId;
flat out vec3 vColor;
flat out float vAlpha;

void main() {
  vec3 center = (uView * vec4(aCenter, 1.0)).xyz;
  vCenter = center;
  vRadius = aRadius;
  vId = aId;
  vColor = mix(aColor, aColor2, uMix);
  vAlpha = uFade * uGroupFade[min(aGroup, uint(MAX_GROUPS - 1))];
  if (vAlpha <= 0.0) {
    // Groupe éteint : le quad sort du volume de vue et n'est jamais rasterisé.
    // Éteindre 44 copies sur 46 doit coûter ce qu'elles coûtent, c'est-à-dire rien.
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
    return;
  }

  // Sous perspective, la silhouette d'une sphère est une ellipse plus large que
  // son rayon : un quad de demi-taille r laisse dépasser le bord. Le facteur
  // exact à la profondeur du centre est r / sqrt(1 - (r/d)^2).
  float d = length(center);
  float k = clamp(aRadius / max(d, 1e-4), 0.0, 0.99);
  float halfSize = aRadius / sqrt(1.0 - k * k);   // half est un mot réservé en GLSL ES

  vec3 pos = center + vec3(aCorner * halfSize, 0.0);
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
flat in float vAlpha;

uniform mat4 uProj;
uniform vec3 uLight;
uniform bool uInvert;

layout(location = 0) out vec4 oColor;
uniform uint uKind;
layout(location = 1) out uvec4 oId;
layout(location = 2) out vec4 oDepth;

// Seuil de tramage par pixel, dans [0, 1) : bruit à gradient entrelacé
// (Jimenez 2014). Fixe à l'écran, il ne scintille pas quand la caméra est immobile.
float threshold(vec2 p) {
  return fract(52.9829189 * fract(dot(p, vec2(0.06711056, 0.00583715))));
}

void main() {
  // Opacité par tramage : un pixel est gardé ou jeté, jamais mélangé. Il n'y a donc
  // ni tri ni ordre de dessin, la profondeur reste exacte, et le picking rend ce qui
  // est vu. Le mode inversé garde exactement les pixels que l'autre jette : un niveau
  // à 1 − t et le suivant à t, inversé, se partagent l'écran **sans trou ni doublon**.
  float h = threshold(gl_FragCoord.xy);
  if (uInvert ? h < 1.0 - vAlpha : h >= vAlpha) discard;

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
  oId = uvec4(vId, uKind, 0u, 0u);
  // La normale en espace vue voyage avec la profondeur : l'occlusion ambiante
  // (ssao.ts) a besoin des deux, et la cible flottante avait trois canaux libres.
  oDepth = vec4(gl_FragDepth, normal);
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
  /** Groupe de chaque bille (sa copie), < MAX_GROUPS. Absent : tout dans le groupe 0. */
  groups?: Uint8Array;
  /** Seconde palette, vers laquelle `uMix` fait glisser. Absente : noir, donc à ne pas mélanger. */
  colors2?: Float32Array;
}

/** Nombre de groupes adressables : 46 copies de chromosomes, avec de la marge. */
export const MAX_GROUPS = 64;

export interface DrawOptions {
  /** Repère monde → vue. Identité par défaut : centres déjà en espace vue. */
  view?: Float32Array;
  light?: [number, number, number];
  /** Opacité du niveau, rendue par tramage. */
  fade?: number;
  /** Garde les pixels complémentaires de ceux qu'un niveau à `1 − fade` garderait. */
  invert?: boolean;
  /** Opacité par groupe ; un groupe absent du tableau vaut 1. */
  groupFade?: ArrayLike<number>;
  /** 0 : première palette ; 1 : seconde. */
  mix?: number;
  /** Écrit dans le canal G de la cible d'identifiants. `KIND.sphere` par défaut. */
  kind?: number;
}

const IDENTITY = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);

export class ImpostorRenderer {
  private readonly program: WebGLProgram;
  private readonly vao: WebGLVertexArrayObject;
  private readonly buffers: WebGLBuffer[] = [];
  private readonly uProj: WebGLUniformLocation;
  private readonly uLight: WebGLUniformLocation;
  private readonly u: Record<"view" | "fade" | "groupFade" | "mix" | "invert" | "kind", WebGLUniformLocation>;
  private readonly fades = new Float32Array(MAX_GROUPS);
  private count = 0;

  constructor(private readonly gl: WebGL2RenderingContext) {
    this.program = link(gl, VERT.replace(/MAX_GROUPS/g, String(MAX_GROUPS)), FRAG);
    this.uProj = gl.getUniformLocation(this.program, "uProj")!;
    this.uLight = gl.getUniformLocation(this.program, "uLight")!;
    const loc = (n: string) => gl.getUniformLocation(this.program, n)!;
    this.u = {
      view: loc("uView"),
      fade: loc("uFade"),
      groupFade: loc("uGroupFade"),
      mix: loc("uMix"),
      invert: loc("uInvert"),
      kind: loc("uKind"),
    };

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
    // Attributs facultatifs. Désactivé, un attribut lit sa « valeur courante », qui
    // doit avoir le **type** déclaré par le shader : un `uint` lu sur la valeur
    // flottante par défaut est une INVALID_OPERATION au dessin. D'où les deux appels
    // explicites, refaits à chaque dessin parce que cette valeur est globale au
    // contexte et qu'un autre programme peut l'avoir changée.
    if (beads.groups) attr(5, beads.groups, 1, gl.UNSIGNED_BYTE, true);
    if (beads.colors2) attr(6, beads.colors2, 3, gl.FLOAT, false);
    this.hasGroups = !!beads.groups;
    this.hasColors2 = !!beads.colors2;

    gl.bindVertexArray(null);
  }

  private hasGroups = false;
  private hasColors2 = false;

  draw(proj: Float32Array, opts: DrawOptions | [number, number, number] = {}): void {
    const gl = this.gl;
    const o: DrawOptions = Array.isArray(opts) ? { light: opts } : opts;
    gl.useProgram(this.program);
    gl.bindVertexArray(this.vao);
    if (!this.hasGroups) gl.vertexAttribI4ui(5, 0, 0, 0, 0);
    if (!this.hasColors2) gl.vertexAttrib4f(6, 0, 0, 0, 1);
    gl.uniformMatrix4fv(this.uProj, false, proj);
    gl.uniform3fv(this.uLight, o.light ?? [0.4, 0.6, 1]);
    gl.uniformMatrix4fv(this.u.view, false, o.view ?? IDENTITY);
    gl.uniform1f(this.u.fade, o.fade ?? 1);
    gl.uniform1i(this.u.invert, o.invert ? 1 : 0);
    gl.uniform1f(this.u.mix, o.mix ?? 0);
    gl.uniform1ui(this.u.kind, o.kind ?? KIND.sphere);
    this.fades.fill(1);
    if (o.groupFade) for (let i = 0; i < Math.min(o.groupFade.length, MAX_GROUPS); i++) this.fades[i] = o.groupFade[i];
    gl.uniform1fv(this.u.groupFade, this.fades);
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
