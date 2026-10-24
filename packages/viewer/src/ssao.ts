/** Occlusion ambiante en espace écran.
 *
 *  Un nuage de billes éclairées d'une seule lumière se lit comme une surface : les
 *  billes du fond et celles de devant ont le même éclairage, et rien ne dit qu'entre
 *  deux territoires il y a un creux. L'occlusion ambiante assombrit ce qui est entouré —
 *  les plis, les contacts, l'intérieur d'une masse — et c'est elle qui rend la
 *  profondeur lisible sans faire tourner la vue.
 *
 *  Méthode : hémisphère orienté par la normale (Crytek 2007, dans la variante à
 *  normales de Chapman). Pour chaque pixel, on reconstruit le point vu depuis la
 *  profondeur écrite par les imposteurs, on lance `SAMPLES` points dans l'hémisphère
 *  de rayon `radius`, et on compte ceux qui tombent **derrière** la surface visible à
 *  leur position d'écran. Une passe de flou qui respecte la profondeur retire le bruit,
 *  et la composition multiplie la couleur par le résultat.
 *
 *  Le rayon est en unités de vue, et le biais lui est proportionnel : une scène
 *  agrandie dix fois avec un rayon dix fois plus grand donne la **même** image. Une
 *  plongée qui traverse deux ordres de grandeur en dépend — le test le vérifie.
 */

import { link, type Targets } from "./gl.js";

export const SAMPLES = 16;

const FULLSCREEN = /* glsl */ `#version 300 es
void main() {
  // Un triangle qui couvre l'écran, sans tampon de sommets.
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
`;

// Construites à la demande : un gabarit avec `${…}` au niveau du module n'est pas éliminé par
// esbuild (il ne peut pas prouver l'absence d'effet de bord), et l'occlusion se retrouvait
// dans des pages qui ne s'en servent pas, comme celle du premier rendu.
const ssaoSource = () => /* glsl */ `#version 300 es
precision highp float;
precision highp int;

uniform highp sampler2D uDepth;   // x : profondeur [0, 1] ; yzw : normale en espace vue
uniform mat4 uProj;
uniform mat4 uInvProj;
uniform float uRadius;
uniform vec3 uKernel[${SAMPLES}];

out vec4 oAo;

vec3 viewPos(ivec2 px, float depth, vec2 size) {
  vec4 ndc = vec4((vec2(px) + 0.5) / size * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
  vec4 v = uInvProj * ndc;
  return v.xyz / v.w;
}

// z en espace vue depuis une profondeur [0, 1], pour une perspective.
float viewZ(float depth) {
  return uProj[3][2] / (-(depth * 2.0 - 1.0) - uProj[2][2]);
}

float noise(vec2 p) {
  return fract(52.9829189 * fract(dot(p, vec2(0.06711056, 0.00583715))));
}

void main() {
  ivec2 px = ivec2(gl_FragCoord.xy);
  vec2 size = vec2(textureSize(uDepth, 0));
  vec4 g = texelFetch(uDepth, px, 0);
  if (g.x >= 1.0) { oAo = vec4(1.0); return; }    // fond : rien à occulter

  vec3 p = viewPos(px, g.x, size);
  vec3 n = normalize(g.yzw);

  // Repère tangent tourné d'un angle propre au pixel : le motif d'échantillonnage
  // change d'un pixel à l'autre, et le flou qui suit transforme ce bruit en douceur.
  float a = 6.2831853 * noise(gl_FragCoord.xy);
  vec3 r = vec3(cos(a), sin(a), 0.0);
  vec3 t = normalize(r - n * dot(r, n));
  if (any(isnan(t))) t = normalize(cross(n, vec3(0.0, 1.0, 0.0)));
  mat3 tbn = mat3(t, cross(n, t), n);

  float bias = 0.02 * uRadius;
  float occ = 0.0;
  for (int i = 0; i < ${SAMPLES}; i++) {
    vec3 s = p + tbn * uKernel[i] * uRadius;
    vec4 c = uProj * vec4(s, 1.0);
    vec2 uv = c.xy / c.w * 0.5 + 0.5;
    if (any(lessThan(uv, vec2(0.0))) || any(greaterThanEqual(uv, vec2(1.0)))) continue;
    float d = texelFetch(uDepth, ivec2(uv * size), 0).x;
    if (d >= 1.0) continue;
    float z = viewZ(d);
    // Le point est occulté si la surface vue à sa place est devant lui. La pondération
    // par la distance évite qu'une bille lointaine, loin derrière, assombrisse le bord
    // d'une bille proche.
    float range = smoothstep(0.0, 1.0, uRadius / abs(p.z - z));
    occ += (z >= s.z + bias ? 1.0 : 0.0) * range;
  }
  oAo = vec4(vec3(1.0 - occ / float(${SAMPLES})), 1.0);
}
`;

const COMPOSE = /* glsl */ `#version 300 es
precision highp float;
precision highp int;

uniform sampler2D uColor;
uniform sampler2D uAo;
uniform highp sampler2D uDepth;
uniform mat4 uProj;
uniform float uRadius;
uniform float uStrength;

out vec4 oColor;

float viewZ(float depth) {
  return uProj[3][2] / (-(depth * 2.0 - 1.0) - uProj[2][2]);
}

void main() {
  ivec2 px = ivec2(gl_FragCoord.xy);
  ivec2 size = textureSize(uAo, 0);
  vec4 col = texelFetch(uColor, px, 0);
  float d0 = texelFetch(uDepth, px, 0).x;
  if (d0 >= 1.0) { oColor = col; return; }
  float z0 = viewZ(d0);

  // Flou 4×4 qui ne mélange que des pixels de la même surface : une différence de
  // profondeur supérieure au rayon d'occlusion sépare deux objets, et leur occlusion
  // ne doit pas baver de l'un sur l'autre.
  float sum = 0.0;
  float wsum = 0.0;
  for (int y = -2; y < 2; y++) {
    for (int x = -2; x < 2; x++) {
      ivec2 q = clamp(px + ivec2(x, y), ivec2(0), size - 1);
      float d = texelFetch(uDepth, q, 0).x;
      if (d >= 1.0) continue;
      float w = abs(viewZ(d) - z0) < uRadius ? 1.0 : 0.0;
      sum += texelFetch(uAo, q, 0).r * w;
      wsum += w;
    }
  }
  float ao = wsum > 0.0 ? sum / wsum : 1.0;
  oColor = vec4(col.rgb * mix(1.0, ao, uStrength), col.a);
}
`;

/** Noyau d'échantillonnage déterministe : des points dans l'hémisphère z > 0, plus
 *  denses près du centre, où ils comptent le plus. Déterministe pour que deux rendus
 *  identiques donnent deux images identiques — les tests en dépendent. */
export function kernel(n = SAMPLES): Float32Array {
  const out = new Float32Array(n * 3);
  let seed = 0x9e3779b9;
  const rand = () => {
    seed ^= seed << 13;
    seed ^= seed >>> 17;
    seed ^= seed << 5;
    return (seed >>> 0) / 4294967296;
  };
  for (let i = 0; i < n; i++) {
    let x: number, y: number, z: number, l: number;
    do {
      x = rand() * 2 - 1;
      y = rand() * 2 - 1;
      z = rand();
      l = Math.hypot(x, y, z);
    } while (l > 1 || l < 1e-3 || z < 0.1);
    const t = (i + 1) / n;
    const k = (0.1 + 0.9 * t * t) / l;
    out.set([x * k, y * k, z * k], i * 3);
  }
  return out;
}

export interface SsaoOptions {
  /** Rayon de l'hémisphère, en unités de vue. */
  radius: number;
  /** 0 : pas d'occlusion ; 1 : occlusion pleine. */
  strength?: number;
}

export class Ssao {
  private readonly ssao: WebGLProgram;
  private readonly compose: WebGLProgram;
  private readonly vao: WebGLVertexArrayObject;
  readonly ao: WebGLTexture;
  readonly fbo: WebGLFramebuffer;
  private readonly kern = kernel();

  constructor(
    private readonly gl: WebGL2RenderingContext,
    private readonly targets: Targets,
  ) {
    this.ssao = link(gl, FULLSCREEN, ssaoSource());
    this.compose = link(gl, FULLSCREEN, COMPOSE);
    this.vao = gl.createVertexArray()!;
    this.ao = gl.createTexture()!;
    gl.bindTexture(gl.TEXTURE_2D, this.ao);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, targets.width, targets.height);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    this.fbo = gl.createFramebuffer()!;
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.ao, 0);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  }

  /** Calcule l'occlusion depuis la cible des imposteurs. */
  occlusion(proj: Float32Array, o: SsaoOptions): void {
    const { gl, targets } = this;
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo);
    gl.viewport(0, 0, targets.width, targets.height);
    gl.disable(gl.DEPTH_TEST);
    gl.useProgram(this.ssao);
    this.bind(this.ssao, "uDepth", 0, targets.depth);
    gl.uniformMatrix4fv(gl.getUniformLocation(this.ssao, "uProj"), false, proj);
    gl.uniformMatrix4fv(gl.getUniformLocation(this.ssao, "uInvProj"), false, invert(proj));
    gl.uniform1f(gl.getUniformLocation(this.ssao, "uRadius"), o.radius);
    gl.uniform3fv(gl.getUniformLocation(this.ssao, "uKernel"), this.kern);
    gl.bindVertexArray(this.vao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }

  /** Floute l'occlusion, la multiplie à la couleur, écrit dans `out` (null : le canevas). */
  composite(proj: Float32Array, o: SsaoOptions, out: WebGLFramebuffer | null = null): void {
    const { gl, targets } = this;
    gl.bindFramebuffer(gl.FRAMEBUFFER, out);
    gl.viewport(0, 0, targets.width, targets.height);
    gl.disable(gl.DEPTH_TEST);
    gl.useProgram(this.compose);
    this.bind(this.compose, "uColor", 0, targets.color);
    this.bind(this.compose, "uAo", 1, this.ao);
    this.bind(this.compose, "uDepth", 2, targets.depth);
    gl.uniformMatrix4fv(gl.getUniformLocation(this.compose, "uProj"), false, proj);
    gl.uniform1f(gl.getUniformLocation(this.compose, "uRadius"), o.radius);
    gl.uniform1f(gl.getUniformLocation(this.compose, "uStrength"), o.strength ?? 1);
    gl.bindVertexArray(this.vao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    gl.bindVertexArray(null);
  }

  /** Occlusion relue en (x, y), origine en haut à gauche, dans [0, 1]. */
  aoAt(x: number, y: number): number {
    const { gl, targets } = this;
    const buf = new Uint8Array(4);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, this.fbo);
    gl.readPixels(Math.round(x), Math.round(targets.height - 1 - y), 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, buf);
    gl.bindFramebuffer(gl.READ_FRAMEBUFFER, null);
    return buf[0] / 255;
  }

  private bind(p: WebGLProgram, name: string, unit: number, tex: WebGLTexture): void {
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.uniform1i(gl.getUniformLocation(p, name), unit);
  }

  dispose(): void {
    const gl = this.gl;
    gl.deleteProgram(this.ssao);
    gl.deleteProgram(this.compose);
    gl.deleteVertexArray(this.vao);
    gl.deleteTexture(this.ao);
    gl.deleteFramebuffer(this.fbo);
  }
}

/** Inverse d'une matrice 4×4 colonne-majeure (cofacteurs). */
export function invert(m: Float32Array): Float32Array {
  const inv = new Float32Array(16);
  inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
  inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
  inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
  inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
  inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
  inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
  inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
  inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
  inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
  inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
  inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
  inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
  inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
  inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
  inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
  inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];
  const det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
  if (det === 0) throw new Error("matrice singulière");
  for (let i = 0; i < 16; i++) inv[i] /= det;
  return inv;
}
