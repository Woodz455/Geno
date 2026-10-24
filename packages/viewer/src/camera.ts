/** Caméra orbitale et chemin de plongée.
 *
 *  Une plongée du noyau (10 µm) aux compartiments (quelques centaines de nm) traverse
 *  deux ordres de grandeur. Interpoler la distance linéairement passerait l'essentiel du
 *  trajet loin de la cible, puis s'écraserait en une image ; l'interpoler en logarithme
 *  corrige la vitesse du zoom mais pas le déplacement latéral, qui file à l'écran dès
 *  qu'on est près. Le chemin de van Wijk et Nuij (2003) résout les deux ensemble : il
 *  minimise la longueur dans l'espace (position, largeur de vue) muni de la métrique
 *  ds = √(ρ² du² + dw² / ρ²) / w, où le déplacement compte en **largeurs d'écran** et le
 *  zoom en **facteurs**. Parcouru à pas constants, il a une vitesse perçue constante :
 *  c'est ce que vérifie le test, pas seulement ses extrémités.
 */

export type Vec3 = [number, number, number];

export interface Orbit {
  target: Vec3;
  /** Rotation autour de l'axe y, radians. */
  yaw: number;
  /** Élévation, radians, bornée à ±89°. */
  pitch: number;
  distance: number;
}

export const PITCH_LIMIT = (89 * Math.PI) / 180;

export function eye(o: Orbit): Vec3 {
  const cp = Math.cos(o.pitch);
  return [
    o.target[0] + o.distance * cp * Math.sin(o.yaw),
    o.target[1] + o.distance * Math.sin(o.pitch),
    o.target[2] + o.distance * cp * Math.cos(o.yaw),
  ];
}

/** Matrice de vue, colonne-majeure comme l'attend WebGL. */
export function lookAt(from: Vec3, to: Vec3, up: Vec3 = [0, 1, 0]): Float32Array {
  const f = norm(sub(to, from));
  const s = norm(cross(f, up));
  const u = cross(s, f);
  // prettier-ignore
  return new Float32Array([
    s[0], u[0], -f[0], 0,
    s[1], u[1], -f[1], 0,
    s[2], u[2], -f[2], 0,
    -dot(s, from), -dot(u, from), dot(f, from), 1,
  ]);
}

export function viewMatrix(o: Orbit): Float32Array {
  return lookAt(eye(o), o.target);
}

/** Tourne l'orbite d'un glissement écran, en radians. */
export function rotate(o: Orbit, dYaw: number, dPitch: number): Orbit {
  return {
    ...o,
    yaw: o.yaw + dYaw,
    pitch: Math.max(-PITCH_LIMIT, Math.min(PITCH_LIMIT, o.pitch + dPitch)),
  };
}

/** Zoom multiplicatif : un cran de molette vaut un facteur, pas une longueur. */
export function zoom(o: Orbit, factor: number, min = 0, max = Infinity): Orbit {
  return { ...o, distance: Math.max(min, Math.min(max, o.distance * factor)) };
}

// ---------------------------------------------------------------------------
// Chemin de van Wijk et Nuij
// ---------------------------------------------------------------------------

/** ρ = √2 : le compromis que van Wijk et Nuij retiennent après essais utilisateurs. */
export const RHO = Math.SQRT2;

export interface ZoomPath {
  /** Longueur du chemin dans la métrique ci-dessus. Un pas de s constant = vitesse perçue constante. */
  length: number;
  /** Vue au paramètre s ∈ [0, length]. */
  at(s: number): { center: Vec3; width: number };
}

/** Chemin optimal entre deux vues, chacune donnée par son centre et la largeur visible au centre. */
export function zoomPath(c0: Vec3, w0: number, c1: Vec3, w1: number, rho = RHO): ZoomPath {
  const d = sub(c1, c0);
  const u1 = Math.hypot(d[0], d[1], d[2]);
  const along = (u: number): Vec3 => (u1 > 0 ? add(c0, scale(d, u / u1)) : [...c0]);

  // Sans déplacement latéral, le chemin dégénère en zoom pur, exponentiel en s.
  if (u1 < 1e-9 * Math.max(w0, w1)) {
    const k = Math.log(w1 / w0) / rho;
    const length = Math.abs(k);
    return {
      length,
      at: (s) => ({ center: [...c0], width: w0 * Math.exp(Math.sign(k) * rho * s) }),
    };
  }

  const r2 = rho * rho;
  const r4 = r2 * r2;
  const b0 = (w1 * w1 - w0 * w0 + r4 * u1 * u1) / (2 * w0 * r2 * u1);
  const b1 = (w1 * w1 - w0 * w0 - r4 * u1 * u1) / (2 * w1 * r2 * u1);
  const rr0 = Math.log(-b0 + Math.sqrt(b0 * b0 + 1));
  const rr1 = Math.log(-b1 + Math.sqrt(b1 * b1 + 1));
  const length = (rr1 - rr0) / rho;
  const a = (w0 / r2) * Math.cosh(rr0);
  return {
    length,
    at: (s) => {
      const t = rho * s + rr0;
      return {
        center: along(a * Math.tanh(t) - (w0 / r2) * Math.sinh(rr0)),
        width: (w0 * Math.cosh(rr0)) / Math.cosh(t),
      };
    },
  };
}

/** Largeur visible au centre d'une orbite, et l'inverse. */
export function widthAt(distance: number, fovY: number): number {
  return 2 * distance * Math.tan(fovY / 2);
}

export function distanceFor(width: number, fovY: number): number {
  return width / (2 * Math.tan(fovY / 2));
}

/** Plongée d'une orbite à une autre. L'orientation suit la progression, adoucie aux bouts. */
export function dive(from: Orbit, to: Orbit, fovY: number): { length: number; at(s: number): Orbit } {
  const p = zoomPath(from.target, widthAt(from.distance, fovY), to.target, widthAt(to.distance, fovY));
  // L'angle de lacet le plus court : on ne fait pas un tour complet pour 10°.
  let dy = to.yaw - from.yaw;
  dy -= 2 * Math.PI * Math.round(dy / (2 * Math.PI));
  return {
    length: p.length,
    at: (s) => {
      const v = p.at(Math.max(0, Math.min(p.length, s)));
      const x = p.length > 0 ? Math.max(0, Math.min(1, s / p.length)) : 1;
      const k = x * x * (3 - 2 * x);
      return {
        target: v.center,
        distance: distanceFor(v.width, fovY),
        yaw: from.yaw + dy * k,
        pitch: from.pitch + (to.pitch - from.pitch) * k,
      };
    },
  };
}

// ---------------------------------------------------------------------------

function sub(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}
function add(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}
function scale(a: Vec3, k: number): Vec3 {
  return [a[0] * k, a[1] * k, a[2] * k];
}
function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
function cross(a: Vec3, b: Vec3): Vec3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
function norm(a: Vec3): Vec3 {
  const l = Math.hypot(a[0], a[1], a[2]) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
}
