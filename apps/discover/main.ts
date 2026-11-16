/** Geno — mode Découverte (semaines 19–20) : « Du noyau à une lettre ».
 *
 *  La plongée de la semaine 12 (`openDive`, cible TP53), avec un préréglage et une couche
 *  narrative : sept étapes du noyau à une base, quatre haltes commentées, un glossaire au survol et
 *  au focus, des échelles comparées calculées (`scale.ts`), une fin de parcours avec une question
 *  ouverte et trois questions à choix. Aucun moteur nouveau : les arrêts sont ceux de la scène,
 *  les haltes passent par `destination()` comme la recherche de la semaine 16, et le badge
 *  d'évidence lit `scene.evidence()` — le texte ne peut pas dire « mesuré » quand la couche
 *  affichée est simulée sans que la mesure (`tools/discover.mjs`) le voie.
 *
 *  Semaine 20 : français et anglais (le contenu et l'interface sont des fichiers, `content/`),
 *  et l'accessibilité WCAG 2.2 AA — lien d'évitement, focus sur le titre de chaque écran, une
 *  description en mots de chaque vue 3D (lue par la scène, pas écrite à la main), le glossaire
 *  en boutons décrits, une description chiffrée de la carte de contacts, des animations qu'on
 *  peut couper. Mesuré par `tools/a11y.mjs`.
 *
 *  URL : `#t=<étape ou halte>` reprend le parcours où il était.
 *  Paramètres : lang=fr|en, f=… (.g3d), w=, h= (taille fixe, pour la mesure), motion=0 (sans vol).
 */

import { ContactsReader, FicheReader, HttpSource, SearchIndex } from "@geno/g3d";
import { distanceFor, dive, type Orbit, widthAt } from "@geno/viewer";
import { type DiveScene, FOV, openDive } from "../dive/dive.js";
import { destination } from "../dive/navigate.js";
import contentEn from "./content/en.json";
import contentFr from "./content/fr.json";
import glossaryEn from "./content/glossary.en.json";
import glossaryFr from "./content/glossary.fr.json";
import uiEn from "./content/ui.en.json";
import uiFr from "./content/ui.fr.json";
import { fieldCompare, fieldText, fill, type Locale, values } from "./scale.js";

const q = new URLSearchParams(location.search);
const store = {
  get: (k: string): string | null => { try { return localStorage.getItem(k); } catch { return null; } },
  set: (k: string, v: string): void => { try { localStorage.setItem(k, v); } catch { /* navigation privée */ } },
};

// — Langue : ?lang=, sinon le choix gardé, sinon celle du navigateur —
const asked = q.get("lang");
const LANG: "fr" | "en" = asked === "en" || asked === "fr" ? asked
  : (store.get("geno.decouverte.lang") as "fr" | "en" | null) ?? (navigator.language.toLowerCase().startsWith("en") ? "en" : "fr");
if (asked === "en" || asked === "fr") store.set("geno.decouverte.lang", asked);
const content: typeof contentFr = LANG === "en" ? (contentEn as typeof contentFr) : contentFr;
const GLOSS: Record<string, string> = LANG === "en" ? glossaryEn : glossaryFr;
const UI: typeof uiFr = LANG === "en" ? (uiEn as typeof uiFr) : uiFr;
const LOC = UI.locale as Locale;
const VALS = values(LOC);
/** Une chaîne de l'interface, ses `{clés}` remplies. */
const T = (s: string, v: Record<string, string | number> = {}): string => s.replace(/\{(\w+)\}/g, (m, k: string) => (k in v ? String(v[k]) : m));
const num = (x: number, d = 0): string => x.toLocaleString(LOC, { maximumFractionDigits: d, minimumFractionDigits: d });

type Step = (typeof content.steps)[number];
type Halt = (typeof content.halts)[number];
type Stop = { kind: "step"; index: number; step: Step } | { kind: "halt"; index: number; halt: Halt };

const FILE = q.get("f") ?? "data/gm12878.g3d";
const CONTACTS = FILE.replace(/\.g3d$/, ".contacts");
const SEARCH = FILE.replace(/\.g3d$/, ".search");
const FICHE = FILE.replace(/\.g3d$/, ".fiche");
const FIXED = q.has("w");
// Mouvement : ?motion=0 et la préférence du système d'abord ; puis le bouton de la page, gardé.
const savedMotion = store.get("geno.decouverte.motion");
let reduced = q.get("motion") === "0" || (savedMotion === null ? matchMedia("(prefers-reduced-motion: reduce)").matches : savedMotion === "0");

// — Style —
// Jetons de couleur choisis pour les seuils AA (semaine 20) : texte ≥ 4,5:1 sur le panneau,
// bordures de composants et contours de focus ≥ 3:1. `tools/a11y.mjs` les recalcule.
const style = document.createElement("style");
style.textContent = `
  :root { --bg: #0b0f14; --panel: #0e1218; --ink: #e8edf3; --mute: #a9b4c4; --line: #2b3646;
          --edge: #6b7a90; --accent: #f0e342; --link: #9cc3ff; }
  html, body { height: 100%; margin: 0; overflow: hidden; background: var(--bg); color: var(--ink);
               font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }
  canvas { display: block; touch-action: none; }
  #overlay { position: fixed; left: 0; top: 0; pointer-events: none; }
  .skip { position: fixed; left: 8px; top: -60px; z-index: 20; background: var(--panel); color: var(--ink);
          padding: 8px 12px; border: 2px solid var(--link); border-radius: 8px; }
  .skip:focus { top: 8px; }
  #story { position: fixed; right: 0; top: 0; bottom: 0; width: min(400px, 100vw); box-sizing: border-box;
           background: var(--panel); border-left: 1px solid var(--line); padding: 18px 20px 14px;
           display: flex; flex-direction: column; gap: 10px; overflow-y: auto; overflow-x: hidden; overflow-wrap: anywhere; }
  #story:focus { outline: none; }
  #story header { display: flex; flex-wrap: wrap; justify-content: space-between; gap: 6px 12px; align-items: baseline; }
  #story header .kicker { color: var(--mute); font-size: 12px; letter-spacing: .04em; text-transform: uppercase; }
  #story header .tools { display: flex; gap: 8px; align-items: center; font-size: 13px; }
  #story header .tools a, #story header .tools button { min-height: 24px; }
  #story h1 { font-size: 15px; margin: 2px 0 0; font-weight: 600; color: var(--mute); flex-basis: 100%; }
  #story h2 { font-size: 22px; line-height: 1.25; margin: 4px 0 0; }
  #story h2:focus { outline: none; }
  #story h2:focus-visible { outline: 2px solid var(--link); outline-offset: 2px; }
  #story p { margin: 0; }
  #story .show { color: var(--mute); font-size: 13.5px; }
  #story .summary { border-left: 3px solid var(--accent); padding: 4px 0 4px 10px; font-weight: 500; }
  a { color: var(--link); }
  a:focus-visible { outline: 2px solid var(--link); outline-offset: 2px; }
  button.term { font: inherit; color: inherit; background: none; border: 0; padding: 0; margin: 0; border-radius: 2px;
                text-decoration: underline dotted 1.5px; text-underline-offset: 3px; cursor: help; text-align: inherit; }
  button.term:hover { background: none; }
  #tip { position: fixed; z-index: 10; max-width: min(300px, calc(100vw - 16px)); box-sizing: border-box; background: #1b2330; color: var(--ink);
         border: 1px solid var(--edge); border-radius: 8px; padding: 8px 10px; font-size: 13.5px; line-height: 1.45;
         display: none; box-shadow: 0 6px 20px rgba(0,0,0,.4); }
  #tip b { display: block; margin-bottom: 2px; }
  .badge { display: inline-flex; gap: 6px; align-items: baseline; font-size: 12.5px; padding: 3px 9px;
           border-radius: 999px; border: 1px solid var(--edge); background: #121822; }
  .badge .dot { width: 8px; height: 8px; border-radius: 50%; display: inline-block; align-self: center; }
  #evidence-detail { color: var(--mute); font-size: 12px; }
  #ruler { display: flex; align-items: center; gap: 10px; font-size: 13px; color: var(--mute); }
  #ruler .bar { flex: 0 0 auto; height: 3px; background: var(--ink); width: 80px; position: relative; }
  #ruler .bar::before, #ruler .bar::after { content: ""; position: absolute; top: -4px; width: 2px; height: 11px; background: var(--ink); }
  #ruler .bar::before { left: 0; } #ruler .bar::after { right: 0; }
  #ruler b { color: var(--ink); font-weight: 600; }
  #contacts { width: 100%; height: auto; border-radius: 6px; background: #0e131b; }
  #contacts-caption, #contacts-summary { color: var(--mute); font-size: 12.5px; }
  nav.steps { display: flex; flex-wrap: wrap; gap: 8px; margin-top: auto; padding-top: 8px; }
  button { font: inherit; padding: 7px 14px; border-radius: 8px; border: 1px solid var(--edge);
           background: #1b2330; color: inherit; cursor: pointer; min-height: 24px; }
  button:hover { background: #243042; }
  button:focus-visible, input:focus-visible, textarea:focus-visible { outline: 2px solid var(--link); outline-offset: 2px; }
  button.primary { background: #2c4a73; border-color: #7f9cc4; }
  button[disabled] { opacity: .45; cursor: default; }
  button.small { font-size: 13px; padding: 2px 8px; }
  ol.dots { display: flex; gap: 5px; align-items: center; list-style: none; margin: 0; padding: 0; }
  ol.dots li { width: 8px; height: 8px; border-radius: 50%; background: #56647a; }
  ol.dots li.on { background: var(--accent); }
  #halts { display: flex; flex-wrap: wrap; gap: 6px; }
  #halts button { font-size: 13px; padding: 4px 10px; }
  #halts-label { color: var(--mute); font-size: 12px; }
  .quiz fieldset { border: 1px solid var(--edge); border-radius: 8px; margin: 0 0 8px; padding: 8px 10px; min-width: 0; }
  .quiz legend { padding: 0 4px; font-weight: 500; }
  .quiz label { display: flex; gap: 8px; align-items: center; padding: 3px 0; cursor: pointer; }
  .quiz input[type=radio] { width: 18px; height: 18px; margin: 0; flex: 0 0 auto; }
  .quiz .verdict { font-size: 13px; margin-top: 4px; }
  .quiz textarea { width: 100%; box-sizing: border-box; min-height: 90px; font: inherit; color: inherit;
                   background: #0e131b; border: 1px solid var(--edge); border-radius: 6px; padding: 6px 8px; }
  .sr { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); clip-path: inset(50%); white-space: nowrap; margin: -1px; padding: 0; border: 0; }
  @media (max-width: 760px) { #story { top: auto; height: 55vh; width: 100vw; border-left: 0; border-top: 1px solid var(--line); } }
`;
document.head.appendChild(style);
document.documentElement.lang = LANG;
document.title = `Geno — ${content.title}`;

// — Squelette —
const skip = document.createElement("a");
skip.className = "skip";
skip.href = "#story";
skip.textContent = UI.skip;
document.body.appendChild(skip);
const canvas = document.createElement("canvas");
canvas.setAttribute("aria-hidden", "true");
document.body.appendChild(canvas);
const overlay = document.createElement("canvas");
overlay.id = "overlay";
overlay.setAttribute("aria-hidden", "true");
document.body.appendChild(overlay);
overlay.getContext("2d", { willReadFrequently: true });
const story = document.createElement("main");
story.id = "story";
story.tabIndex = -1;
story.setAttribute("aria-label", content.title);
document.body.appendChild(story);
// La vue 3D en mots : réécrite à chaque arrivée depuis ce que la scène montre au centre.
const view = document.createElement("div");
view.id = "view";
view.className = "sr";
view.setAttribute("role", "status");
view.setAttribute("aria-live", "polite");
view.setAttribute("aria-label", UI.viewLabel);
document.body.appendChild(view);
// La bulle du glossaire est visuelle : chaque terme porte déjà sa définition (aria-describedby).
const tip = document.createElement("div");
tip.id = "tip";
tip.setAttribute("aria-hidden", "true");
document.body.appendChild(tip);

const panelWidth = (): number => (window.innerWidth <= 760 ? 0 : Math.min(400, window.innerWidth));
function fit(): [number, number] {
  if (FIXED) {
    const w = Number(q.get("w") ?? 1280), h = Number(q.get("h") ?? 720);
    for (const c of [canvas, overlay]) {
      c.width = w;
      c.height = h;
      c.style.width = `${w}px`;
      c.style.height = `${h}px`;
    }
    return [w, h];
  }
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const cssW = Math.max(1, window.innerWidth - panelWidth());
  const cssH = window.innerWidth <= 760 ? Math.round(window.innerHeight * 0.45) : window.innerHeight;
  for (const c of [canvas, overlay]) {
    c.width = Math.round(cssW * dpr);
    c.height = Math.round(cssH * dpr);
    c.style.width = `${cssW}px`;
    c.style.height = `${cssH}px`;
  }
  return [canvas.width, canvas.height];
}

// — Glossaire : [[terme]] ou [[terme|texte affiché]] —
// Chaque terme est un bouton ; sa définition est dans le DOM (cachée), liée par aria-describedby :
// un lecteur d'écran la lit au focus, sans dépendre de la bulle.
const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]!);
let tipFor: HTMLElement | null = null;
let termSeq = 0;
function rich(text: string): string {
  return esc(fill(text, VALS)).replace(/\[\[([^\]|]+)(?:\|([^\]]+))?\]\]/g, (_, term: string, shown?: string) => {
    if (!(term in GLOSS)) throw new Error(`terme absent du glossaire : ${term}`);
    const id = `def-${++termSeq}`;
    return `<button type="button" class="term" aria-describedby="${id}" data-term="${esc(term)}">${shown ?? term}</button><span id="${id}" hidden>${esc(GLOSS[term])}</span>`;
  });
}
function showTip(el: HTMLElement): void {
  const term = el.dataset.term!;
  tip.innerHTML = `<b>${esc(term)}</b>${esc(GLOSS[term])}`;
  tip.style.display = "block";
  const r = el.getBoundingClientRect();
  const w = tip.offsetWidth;
  tip.style.left = `${Math.max(8, Math.min(window.innerWidth - w - 8, r.left))}px`;
  const below = r.bottom + 8, h = tip.offsetHeight;
  tip.style.top = `${below + h < window.innerHeight ? below : Math.max(8, r.top - h - 8)}px`;
  tipFor = el;
}
function hideTip(): void {
  tip.style.display = "none";
  tipFor = null;
}
// La bulle reste tant que le pointeur est sur le terme **ou sur elle** (WCAG 1.4.13), et Échap
// la ferme sans déplacer le focus.
const leaving = (e: MouseEvent) => { const to = e.relatedTarget as Node | null; if (!to || !(tip.contains(to) || (tipFor && tipFor.contains(to)))) hideTip(); };
story.addEventListener("mouseover", (e) => { const t = (e.target as HTMLElement).closest<HTMLElement>(".term"); if (t) showTip(t); });
story.addEventListener("mouseout", (e) => { if ((e.target as HTMLElement).closest(".term")) leaving(e); });
tip.addEventListener("mouseout", leaving);
story.addEventListener("focusin", (e) => { const t = (e.target as HTMLElement).closest<HTMLElement>(".term"); if (t) showTip(t); else hideTip(); });
story.addEventListener("focusout", (e) => { if (!(e.relatedTarget as HTMLElement | null)?.closest?.(".term")) hideTip(); });
story.addEventListener("click", (e) => { const t = (e.target as HTMLElement).closest<HTMLElement>(".term"); if (t) (tipFor === t && tip.style.display === "block" ? hideTip() : showTip(t)); });

// — Évidence : ce que la couche affichée est vraiment —
const EVIDENCE_COLORS: Record<string, string> = { "mesuré": "#5ec27a", "simulé": "#e0a84a", "prédit": "#c58be0", "illustratif": "#7aa8e6" };
/** Les catégories que couvre le texte d'évidence de la scène (`DiveScene.evidence`). */
export function categories(evidence: string): string[] {
  if (evidence.startsWith("séquence")) return ["mesuré"];
  if (evidence.startsWith("nucléosomes")) return ["illustratif"];
  if (evidence.startsWith("fibre")) return ["simulé"];
  return ["simulé", "prédit"];
}
/** Le détail d'évidence dans la langue de la page : la scène rend l'une de quatre phrases
 *  françaises fixes ; on la reconnaît **entière** et on rend sa traduction. Une phrase inconnue
 *  passe telle quelle — et la mesure la compte comme une faute. */
const EV_FR = uiFr.evidenceDetail as Record<string, string>;
function evidenceDetail(ev: string): { text: string; known: boolean } {
  const key = Object.keys(EV_FR).find((k) => EV_FR[k] === ev);
  return key ? { text: (UI.evidenceDetail as Record<string, string>)[key], known: true } : { text: ev, known: false };
}
const evLabel = (cat: string): string => (UI.evidenceLabels as Record<string, string>)[cat] ?? cat;

/** L'en-tête du panneau : la langue (le lien garde l'étape) et le bouton des animations. */
function headerHtml(): string {
  const other = UI.switchLang;
  const p = new URLSearchParams(location.search);
  p.set("lang", other);
  const href = `${location.pathname}?${p.toString()}${location.hash}`;
  return `<header><div class="kicker">${esc(UI.kicker)}</div><div class="tools"><a id="lang" href="${esc(href)}" lang="${other}" hreflang="${other}" title="${esc(UI.switchLabel)}">${esc(UI.switchTo)}</a>` +
    `<button type="button" class="small" id="motion">${esc(reduced ? UI.motionOff : UI.motionOn)}</button></div><h1>${esc(content.title)}</h1></header>`;
}
story.addEventListener("click", (e) => {
  const b = (e.target as HTMLElement).closest<HTMLButtonElement>("#motion");
  if (!b) return;
  reduced = !reduced;
  store.set("geno.decouverte.motion", reduced ? "0" : "1");
  b.textContent = reduced ? UI.motionOff : UI.motionOn;
});
// Le lien de langue suit l'étape courante (le #t= change sans recharger).
window.addEventListener("hashchange", () => { const a = story.querySelector<HTMLAnchorElement>("#lang"); if (a) a.href = a.href.replace(/#.*$/, "") + location.hash; });

async function main(): Promise<void> {
  const [w0, h0] = fit();
  story.innerHTML = `${headerHtml()}<p role="status">${esc(UI.loading)}</p>`;
  const scene: DiveScene = await openDive(canvas, FILE, { target: "tp53" });
  scene.resize(w0, h0);
  // Toute la fenêtre de la fibre, générée une fois : la fenêtre de séquence de l'hélice ne retient
  // que les monomères déjà générés, et un vol en génère d'autres qu'un saut direct — l'image d'un
  // arrêt dépendait du chemin (5 pixels d'un niveau à l'hélice, semaine 20).
  scene.fiber?.grow(scene.final.target as number[], Number.POSITIVE_INFINITY);
  const index = SearchIndex.open(new HttpSource(SEARCH)).catch(() => null);
  const fiche = FicheReader.open(new HttpSource(FICHE)).catch(() => null);
  const contacts = ContactsReader.open(new HttpSource(CONTACTS)).catch(() => null);

  // Les arrêts des étapes : ceux de la plongée de la semaine 12, et trois vues réglées sur le
  // calendrier des fondus (`SCHEDULE`, scene.ts) pour que chaque étape montre sa couche.
  const S = scene.stops;
  const at = (o: Orbit, widthUm: number): Orbit => ({ ...o, distance: distanceFor(widthUm, FOV) });
  const chr17 = scene.region?.chrom ?? "chr17";
  const level = scene.region ? scene.g3d.level(scene.region.level) : null;
  // Le TAD de l'étape : le domaine CTCF de la **population** qui contient la base visée — la
  // paire d'ancres convergentes la plus étroite autour d'elle, au-dessus du seuil des arcs
  // (fiber.arcs). Dans la conformation affichée, aucune cohésine ne tient de boucle autour de
  // la base de TP53 (scene.fiber.focus est nul, comme à la semaine 12) : le domaine vient donc
  // de l'ensemble des conformations, pas de l'instant montré.
  const F = scene.fiber;
  const tm = scene.target && scene.region && level ? Math.floor((scene.target.position - 1 - scene.region.start) / level.bp_per_bead) : -1;
  const focus: { i: number; j: number } | null = F?.focus
    ?? (F && tm >= 0 ? F.arcs.filter((a) => a.i <= tm && tm <= a.j).sort((a, b) => (a.j - a.i) - (b.j - b.i))[0] ?? null : null);
  const pt = (k: number): [number, number, number] => [F!.points[3 * k], F!.points[3 * k + 1], F!.points[3 * k + 2]];
  let domainView: Orbit | null = null, anchorView: Orbit | null = null;
  if (F && focus) {
    const c: [number, number, number] = [0, 0, 0];
    for (let k = focus.i; k <= focus.j; k++) for (let q = 0; q < 3; q++) c[q] += F.points[3 * k + q] / (focus.j - focus.i + 1);
    let extent = 0;
    for (let k = focus.i; k <= focus.j; k++) extent = Math.max(extent, Math.hypot(...pt(k).map((v, q) => v - c[q])));
    // Même cadrage que la vue de boucle de la plongée (dive.ts) : 2,4 × l'étendue, entre 0,3 et 0,58 µm.
    domainView = { ...scene.final, target: c, distance: distanceFor(Math.max(0.3, Math.min(0.58, 2.4 * extent)), FOV) };
    anchorView = { ...scene.final, target: pt(focus.i), distance: distanceFor(0.3, FOV) };
  }
  const orbits: Record<string, Orbit> = {
    noyau: scene.overview,
    territoires: at(scene.final, 16),             // les autres chromosomes encore là (20 → 8 µm), assombris
    compartiments: at(scene.final, 3.2),          // compartiments pleins (< 3,5 µm), pas encore de fibre (> 3 µm)
    tad: scene.loopView ?? domainView ?? S[2],
    boucle: scene.anchorView ? at(scene.anchorView, 0.3) : anchorView ?? S[2],   // l'ancre, avant les nucléosomes (< 0,25 µm)
    fibre: S[3],
    helice: S[S.length - 1],
  };
  // Une vue dont le centre montre bien l'objet de l'étape. À 16 et 3,2 µm, la bille de TP53 est
  // dans le noyau, derrière d'autres : vue de la plongée (scene.final), le centre de l'image
  // tombait sur chrX ou sur une autre région de chr17. On tourne autour de la même cible — même
  // champ — jusqu'à ce que le picking du centre désigne ce qu'il faut : l'orbite d'origine
  // d'abord, puis 12 azimuts × 5 hauteurs, dans un ordre fixe (la reprise par l'URL refait le
  // même calcul et trouve la même vue).
  const clearView = (o: Orbit, ok: (s: NonNullable<ReturnType<DiveScene["select"]>>) => boolean): Orbit => {
    const tries = [o];
    for (const dp of [0, -0.5, 0.5, -1, 1]) for (let k = dp ? 0 : 1; k < 12; k++) tries.push({ ...o, yaw: o.yaw + (k * Math.PI) / 6, pitch: Math.max(-1.2, Math.min(1.2, o.pitch + dp)) });
    for (const q of tries) {
      scene.render(q);
      const s = scene.select(canvas.width / 2, canvas.height / 2);
      if (s && ok(s)) return q;
    }
    return o;
  };
  const tpos = scene.target ? scene.target.position - 1 : -1;
  orbits.territoires = clearView(orbits.territoires, (s) => s.chrom === chr17);
  orbits.compartiments = clearView(orbits.compartiments, (s) => s.chrom === chr17 && s.start <= tpos && tpos < s.end);
  const domain = focus && scene.region && level
    ? { chrom: scene.region.chrom, start: scene.region.start + focus.i * level.bp_per_bead, end: scene.region.start + (focus.j + 1) * level.bp_per_bead }
    : null;

  const stops: Stop[] = [
    ...content.steps.map((step, index) => ({ kind: "step" as const, index, step })),
    ...content.halts.map((halt, index) => ({ kind: "halt" as const, index, halt })),
  ];
  const idOf = (s: Stop) => (s.kind === "step" ? s.step.id : s.halt.id);
  let cur = -1;          // −1 : l'introduction ; content.steps.length : la fin
  let orbit: Orbit = scene.start;
  let anim: { path: { length: number; at(s: number): Orbit }; to: Orbit; t0: number; duration: number } | null = null;
  let dirty = true;
  const arrivals: { id: string; orbit: Orbit; arrivedAt: number; dest?: unknown }[] = [];
  let arrived = true;
  let settled = true;

  const flyTo = (to: Orbit, jump = false): void => {
    arrived = false;
    settled = false;
    if (reduced || jump) {
      orbit = to;
      anim = null;
      dirty = true;
      return;
    }
    const path = dive(orbit, to, FOV);
    anim = { path, to, t0: performance.now(), duration: 900 * Math.max(0.6, path.length) };
    dirty = true;
  };

  // — Rendu du panneau —
  const dots = () => `<ol class="dots" aria-label="${esc(UI.stepsLabel)}">${content.steps.map((st, k) => `<li class="${k <= cur ? "on" : ""}"${k === cur ? ' aria-current="step"' : ""}><span class="sr">${esc(st.title)}</span></li>`).join("")}</ol>`;
  const haltsHtml = () => `<div id="halts-label">${esc(UI.haltsLabel)}</div><div id="halts" role="group" aria-labelledby="halts-label">${content.halts.map((h) => `<button type="button" data-halt="${h.id}">${esc(h.title.split(",")[0])}</button>`).join("")}</div>`;
  const evidenceHtml = () => `<div><span class="badge" id="badge"><span class="dot" aria-hidden="true"></span><span id="badge-text"></span></span></div><div id="evidence-detail"></div>`;
  const rulerHtml = () => `<div id="ruler" role="group" aria-label="${esc(UI.rulerLabel)}"><span class="bar" aria-hidden="true"></span><span id="ruler-text"></span></div>`;
  const nav = (prevLabel: string, nextLabel: string, prevOff = false, nextOff = false) =>
    `<nav class="steps" aria-label="${esc(UI.navLabel)}"><button type="button" id="prev" ${prevOff ? "disabled" : ""}><span aria-hidden="true">← </span>${esc(prevLabel)}</button><button type="button" class="primary" id="next" ${nextOff ? "disabled" : ""}>${esc(nextLabel)}<span aria-hidden="true"> →</span></button></nav>`;

  let declared = "";
  const render = (): void => {
    hideTip();
    let html = headerHtml();
    if (cur < 0) {
      declared = "";
      html += `<h2 tabindex="-1">${esc(content.subtitle)}</h2><p>${rich(content.intro)}</p>${rulerHtml()}${haltsHtml()}${nav(UI.back, UI.start, true)}`;
    } else if (cur < content.steps.length) {
      const st = content.steps[cur];
      declared = st.evidence;
      html += `<div class="kicker">${esc(T(UI.stepOf, { n: cur + 1, total: content.steps.length }))}</div>${dots()}<h2 tabindex="-1">${esc(st.title)}</h2><p>${rich(st.text)}</p>`;
      if ("how" in st && st.how) html += `<p>${rich(st.how)}</p>`;
      if ("contacts" in st && st.contacts) html += `<canvas id="contacts" width="360" height="190" role="img" aria-label="${esc(UI.contactsLabel)}" aria-describedby="contacts-caption contacts-summary"></canvas><div id="contacts-caption">${esc(UI.contactsCaption)}</div><div id="contacts-summary"></div>`;
      if ("summary" in st && st.summary) html += `<p class="summary">${rich(st.summary)}</p>`;
      html += `<p class="show">${rich(st.show)}</p>${evidenceHtml()}${rulerHtml()}${nav(UI.prev, cur === content.steps.length - 1 ? UI.finish : UI.next, false)}${haltsHtml()}`;
    } else if (cur >= content.steps.length && cur < stops.length) {
      const h = (stops[cur] as Extract<Stop, { kind: "halt" }>).halt;
      declared = "simulé";
      html += `<div class="kicker">${esc(UI.halt)}</div><h2 tabindex="-1">${esc(h.title)}</h2><p>${rich(h.text)}</p><p class="show"><b>${esc(UI.notShown)}</b> ${rich(h.limits)}</p><div id="halt-facts" class="show"></div>${evidenceHtml()}${rulerHtml()}${nav(UI.backToTour, UI.endOfTour)}${haltsHtml()}`;
    } else {
      declared = "";
      const Q = content.quiz;
      html += `<div class="kicker">${esc(UI.endOfTour)}</div><h2 tabindex="-1">${esc(Q.intro)}</h2><form class="quiz" id="quiz"><fieldset><legend>${esc(Q.open)}</legend><label class="sr" for="open">${esc(UI.yourAnswer)}</label><textarea id="open" aria-describedby="open-note"></textarea><div class="show" id="open-note">${esc(Q.openNote)}</div></fieldset>` +
        Q.choices.map((c, k) => `<fieldset aria-describedby="v${k}"><legend>${esc(c.q)}</legend>${c.options.map((o, j) => `<label><input type="radio" name="q${k}" value="${j}"> ${esc(o)}</label>`).join("")}<div class="verdict" id="v${k}" aria-live="polite"></div></fieldset>`).join("") +
        `<button type="submit" class="primary">${esc(UI.check)}</button></form>${haltsHtml()}${nav(UI.reviewLast, UI.restart)}`;
    }
    story.innerHTML = html;
    wire();
    updateLive();
    if (cur >= 0 && cur < content.steps.length && "contacts" in content.steps[cur]) void drawContacts();
    if (cur >= content.steps.length && cur < stops.length) void haltFacts((stops[cur] as Extract<Stop, { kind: "halt" }>).halt);
  };

  const wire = (): void => {
    story.querySelector<HTMLButtonElement>("#prev")?.addEventListener("click", () => back());
    story.querySelector<HTMLButtonElement>("#next")?.addEventListener("click", () => forward());
    story.querySelectorAll<HTMLButtonElement>("[data-halt]").forEach((b) => b.addEventListener("click", () => void goHalt(b.dataset.halt!)));
    const quiz = story.querySelector<HTMLFormElement>("#quiz");
    if (quiz) {
      const saved = store.get("geno.decouverte.reponse") ?? "";
      const ta = quiz.querySelector<HTMLTextAreaElement>("#open")!;
      ta.value = saved;
      ta.addEventListener("input", () => store.set("geno.decouverte.reponse", ta.value));
      quiz.addEventListener("submit", (e) => {
        e.preventDefault();
        content.quiz.choices.forEach((c, k) => {
          const v = quiz.querySelector<HTMLInputElement>(`input[name=q${k}]:checked`);
          const out = quiz.querySelector<HTMLDivElement>(`#v${k}`)!;
          out.textContent = !v ? UI.noAnswer : Number(v.value) === c.answer ? `${UI.right} ${c.explain}` : `${UI.wrong} ${c.explain}`;
        });
      });
    }
  };

  /** Le badge et la règle suivent la scène à chaque image ; le reste du panneau ne bouge pas. */
  const updateLive = (): void => {
    const ev = scene.evidence();
    const badge = story.querySelector<HTMLSpanElement>("#badge-text");
    if (badge) {
      const cats = categories(ev);
      const shown = declared && cats.includes(declared) ? declared : cats[0];
      badge.textContent = T(UI.badge, { cat: evLabel(shown) });
      story.querySelector<HTMLSpanElement>("#badge .dot")!.style.background = EVIDENCE_COLORS[shown] ?? "#888";
      story.querySelector<HTMLDivElement>("#evidence-detail")!.textContent = evidenceDetail(ev).text;
    }
    const ruler = story.querySelector<HTMLSpanElement>("#ruler-text");
    if (ruler) {
      const width = widthAt(orbit.distance, FOV);
      const bar = width * (80 / Math.max(1, canvas.clientWidth));       // la barre fait 80 px de l'image
      ruler.innerHTML = `<b>${fieldText(bar, LOC)}</b> · ${esc(T(UI.ruler, { width: fieldText(width, LOC), compare: fieldCompare(width, LOC) }))}`;
    }
  };

  // — La carte de contacts de l'étape TAD : triangle, diagonale en bas —
  let contactsDrawn: { bins: number; resolution: number; domain: [number, number] } | null = null;
  let contactRatio = NaN;
  const drawContacts = async (): Promise<void> => {
    const cv = story.querySelector<HTMLCanvasElement>("#contacts");
    const r = await contacts;
    if (!cv || !r || !domain) return;
    const ts = r.tileset("fin");
    const len = domain.end - domain.start;
    const lo = Math.max(0, domain.start - len / 2), hi = domain.end + len / 2;
    let z = ts.levels.length - 1;
    while (z > 0 && (hi - lo) / ts.levels[z].resolution > 120) z--;
    const res = ts.levels[z].resolution;
    const a0 = r.absolute(domain.chrom, lo)!, a1 = r.absolute(domain.chrom, hi)!;
    const b0 = Math.floor(a0 / res), n = Math.max(2, Math.ceil(a1 / res) - b0);
    const M = new Float32Array(n * n);
    let max = 0;
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      const v = await r.pixel("fin", z, b0 + j, b0 + i);
      M[i * n + j] = v;
      if (j > i + 1) max = Math.max(max, v);
    }
    // En mots, pour qui ne voit pas le triangle : à séparation égale le long du fil, le contact
    // moyen entre deux bins du domaine, rapporté au contact moyen d'un bin du domaine avec un bin
    // hors du domaine. Les séparations comparées sont celles qui existent des deux côtés.
    const i0 = Math.floor(r.absolute(domain.chrom, domain.start)! / res) - b0, i1 = Math.ceil(r.absolute(domain.chrom, domain.end)! / res) - b0;
    const inD = (k: number) => k >= i0 && k < i1;
    let sumIn = 0, sumOut = 0;
    for (let d = 1; d < i1 - i0; d++) {
      let a = 0, na = 0, c = 0, nc = 0;
      for (let i = 0; i + d < n; i++) {
        const j = i + d, v = M[i * n + j];
        if (inD(i) && inD(j)) { a += v; na++; } else if (inD(i) !== inD(j)) { c += v; nc++; }
      }
      if (na && nc) { sumIn += a / na; sumOut += c / nc; }
    }
    contactRatio = sumOut > 0 ? sumIn / sumOut : NaN;
    const summary = story.querySelector<HTMLDivElement>("#contacts-summary");
    if (summary && Number.isFinite(contactRatio)) summary.textContent = T(UI.contactsSummary, { bins: n, res: num(res / 1000), ratio: num(contactRatio, 1) });
    const ctx = cv.getContext("2d")!;
    const W = cv.width, H = cv.height, cell = W / n;
    ctx.clearRect(0, 0, W, H);
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      const v = M[i * n + j];
      if (!(v > 0)) continue;
      const t = Math.min(1, Math.log1p(v) / Math.log1p(max || 1));
      ctx.fillStyle = `rgb(${Math.round(250 - 205 * t)},${Math.round(245 - 215 * t)},${Math.round(240 - 190 * t)})`;
      // (i, j) en losange : centre à x = (i + j + 1) / 2, hauteur (j − i) / 2.
      const cx = ((i + j + 1) / 2) * cell, cy = H - ((j - i) / 2) * cell;
      ctx.beginPath();
      ctx.moveTo(cx, cy - cell / 2);
      ctx.lineTo(cx + cell / 2, cy);
      ctx.lineTo(cx, cy + cell / 2);
      ctx.lineTo(cx - cell / 2, cy);
      ctx.fill();
    }
    const d0 = (r.absolute(domain.chrom, domain.start)! / res - b0) * cell, d1 = (r.absolute(domain.chrom, domain.end)! / res - b0) * cell;
    ctx.strokeStyle = "#f0e342";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(d0, H);
    ctx.lineTo((d0 + d1) / 2, H - (d1 - d0) / 2);
    ctx.lineTo(d1, H);
    ctx.stroke();
    contactsDrawn = { bins: n, resolution: res, domain: [domain.start, domain.end] };
  };

  // — Les faits d'une halte, lus dans la fiche de sa bille —
  const haltDest: Record<string, Awaited<ReturnType<typeof destination>>> = {};
  const haltFacts = async (h: Halt): Promise<void> => {
    const out = story.querySelector<HTMLDivElement>("#halt-facts");
    const d = haltDest[h.id];
    const f = await fiche;
    if (!out || !d || !f) return;
    const card = f.card(d.level, d.bead) ?? (await f.cardAsync(d.level, d.bead));
    const genes: string[] = ((card?.genes?.v?.genes ?? []) as { name: string }[]).map((g) => g.name).filter(Boolean);
    // Les gènes nommés d'abord ceux dont parle la halte, s'ils sont dans la bille ; puis les
    // autres symboles, sans les identifiants Ensembl sans nom.
    const named = [...new Set(genes.filter((g) => !/^ENSG\d/.test(g)))];
    const told = named.filter((g) => new RegExp(`(?<![\\w-])${g.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}(?![\\w-])`).test(`${h.gene} ${h.text}`));
    const shown = [...told, ...named.filter((g) => !told.includes(g))].slice(0, 6);
    const list = `${shown.join(", ")}${named.length > shown.length ? "…" : ""}`;
    out.textContent = T(UI.haltFacts, { kb: num(Math.round((d.end - d.start) / 1000)), copy: d.copy, start: num(d.start + 1), end: num(d.end) }) +
      (genes.length ? ` · ${T(genes.length > 1 ? UI.genesMany : UI.genesOne, { n: num(genes.length), list })}` : "");
  };

  // — Navigation —
  const enter = async (k: number, jump = false, focusTitle = !jump): Promise<void> => {
    cur = k;
    scene.fiber?.highlight(-1, -1);
    scene.highlight(null);
    if (k >= 0 && k < content.steps.length) {
      const st = content.steps[k];
      if (st.id === "territoires") scene.highlight([{ chrom: chr17, start: 0, end: 1e9 }]);
      if ((st.id === "tad" || st.id === "boucle") && focus) scene.fiber?.highlight(focus.i, focus.j);
      flyTo(orbits[st.stop], jump);
      arrivals.push({ id: st.id, orbit: orbits[st.stop], arrivedAt: -1 });
    } else if (k >= content.steps.length && k < stops.length) {
      const h = (stops[k] as Extract<Stop, { kind: "halt" }>).halt;
      const d = await haltDestination(h);
      if (d) {
        flyTo(d.orbit, jump);
        arrivals.push({ id: h.id, orbit: d.orbit, arrivedAt: -1, dest: d });
      }
    }
    render();
    // Le focus va au titre du nouvel écran : un lecteur d'écran l'annonce, et Tab repart de là.
    if (focusTitle) story.querySelector<HTMLElement>("h2")?.focus();
    view.textContent = "";
    try { history.replaceState(null, "", `${location.pathname}${location.search}#t=${k < 0 ? "intro" : k < stops.length ? idOf(stops[k]) : "fin"}`); } catch { /* fichier local */ }
    const a = story.querySelector<HTMLAnchorElement>("#lang");
    if (a) a.href = a.href.replace(/#.*$/, "") + location.hash;
  };
  const haltDestination = async (h: Halt) => {
    if (haltDest[h.id]) {
      scene.highlight(haltRegions[h.id]);
      return haltDest[h.id];
    }
    const ix = await index;
    if (!ix) return null;
    const res = await ix.resolve(h.gene);
    const hit = res.hits.find((x) => x.in_model);
    if (!hit) return null;
    // Le X inactif : les deux X surlignés en entier, la caméra sur la bille de XIST.
    const regions = h.id === "xi" ? [{ chrom: hit.chrom, start: 0, end: 1e9 }] : [{ chrom: hit.chrom, start: hit.start, end: hit.end }];
    haltRegions[h.id] = regions;
    scene.highlight(regions);
    const pos = hit.strand === "-" ? hit.end - 1 : hit.start;
    // L'orientation part de la vue d'ensemble, pas de la vue courante : une halte arrive à la même
    // vue d'où qu'on vienne, et la reprise par l'URL la retrouve.
    const d = await destination(scene, { chrom: hit.chrom, start: hit.start, end: hit.end, pos }, scene.overview);
    // Même règle que les étapes : le centre de l'image sur la bille qui porte le gène.
    if (d) d.orbit = clearView(d.orbit, (s) => s.chrom === hit.chrom && s.start <= pos && pos < s.end);
    if (d) haltDest[h.id] = d;
    return d;
  };
  const haltRegions: Record<string, { chrom: string; start: number; end: number }[]> = {};
  const goHalt = async (id: string): Promise<void> => {
    const k = stops.findIndex((s) => idOf(s) === id);
    if (k >= 0) await enter(k);
  };
  let lastStep = 0;
  const forward = (): void => {
    if (cur < 0) void enter(0);
    else if (cur < content.steps.length - 1) void enter(cur + 1);
    else if (cur === content.steps.length - 1) { lastStep = cur; void enter(stops.length); }
    else if (cur >= content.steps.length && cur < stops.length) void enter(stops.length);
    else void enter(-1);
  };
  const back = (): void => {
    if (cur <= 0) void enter(-1);
    else if (cur < content.steps.length) void enter(cur - 1);
    else void enter(Math.max(0, Math.min(lastStep, content.steps.length - 1)));
  };
  document.addEventListener("keydown", (e) => {
    const t = e.target as HTMLElement;
    if (t.closest("textarea, input")) return;
    if (e.key === "ArrowRight") { e.preventDefault(); forward(); }
    else if (e.key === "ArrowLeft") { e.preventDefault(); back(); }
    else if (e.key === "Escape") {
      // Échap ferme d'abord la bulle, sans toucher au focus (WCAG 1.4.13) ; sinon, quitte la halte.
      if (tip.style.display === "block") { e.preventDefault(); hideTip(); }
      else if (cur >= content.steps.length && cur < stops.length) back();
    }
  });
  window.addEventListener("resize", () => {
    if (FIXED) return;
    scene.resize(...fit());
    dirty = true;
  });

  // — La vue en mots, à chaque arrivée : ce que le picking voit au centre de l'image —
  let described: { what: string; chrom: string; cat: string; known: boolean } | null = null;
  const describe = (): void => {
    const field = fieldText(widthAt(orbit.distance, FOV), LOC);
    const s = cur >= 0 && cur < stops.length ? scene.select(canvas.width / 2, canvas.height / 2) : null;
    if (!s) { described = null; view.textContent = cur >= 0 && cur < stops.length ? T(UI.describeNone, { field }) : ""; return; }
    const cats = categories(scene.evidence());
    const cat = declared && cats.includes(declared) ? declared : cats[0];
    const whatMap = UI.what as Record<string, string>;
    const what = whatMap[s.what] ?? s.what;
    const chrom = T(UI.chromosome, { n: s.chrom.replace(/^chr/, "") });
    const range = s.pos >= 0 ? num(s.pos + 1) : `${num(s.start + 1)}–${num(s.end)}`;
    view.textContent = T(UI.describe, { what, chrom, range, field, cat: evLabel(cat) });
    described = { what: s.what, chrom: s.chrom, cat, known: s.what in whatMap };
  };

  // — Boucle d'images —
  const labels = (): void => {
    const ctx = overlay.getContext("2d")!;
    ctx.clearRect(0, 0, overlay.width, overlay.height);
    for (const l of scene.labels()) {
      if (l.target) continue;
      ctx.globalAlpha = l.alpha;
      const [r, g, b] = l.color.map((c) => Math.round(c * 255));
      ctx.font = `700 ${Math.round(l.size)}px ui-monospace, monospace`;
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.lineWidth = Math.max(2, l.size / 6);
      ctx.strokeStyle = "rgba(10,12,16,0.85)";
      ctx.strokeText(l.text, l.x, l.y);
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fillText(l.text, l.x, l.y);
    }
    ctx.globalAlpha = 1;
  };
  const loop = (now: number): void => {
    if (anim) {
      const k = Math.min(1, (now - anim.t0) / anim.duration);
      // L'arrivée est l'orbite visée, exactement : la reprise par l'URL rejoue la même image.
      orbit = k >= 1 ? anim.to : anim.path.at(k * anim.path.length);
      if (k >= 1) anim = null;
      dirty = true;
    }
    if (scene.busy) dirty = true;
    if (dirty) {
      scene.render(orbit);
      labels();
      updateLive();
      dirty = false;
    }
    if (!anim && !arrived && !scene.busy && !settled) {
      // À l'arrêt, la fenêtre de séquence est recalculée là où est la caméra, puis l'image
      // refaite : l'image d'un arrêt ne dépend plus du chemin qui y a mené.
      settled = true;
      scene.sequence?.settle();
      dirty = true;
    } else if (!anim && !arrived && !scene.busy) {
      arrived = true;
      describe();
      const a = arrivals[arrivals.length - 1];
      if (a && a.arrivedAt < 0) a.arrivedAt = now;
    }
    requestAnimationFrame(loop);
  };

  // L'URL reprend le parcours où il était.
  const t = new URLSearchParams(location.hash.slice(1)).get("t");
  const start = t === "fin" ? stops.length : t && t !== "intro" ? stops.findIndex((s) => idOf(s) === t) : -1;
  // Reprise : on y est d'emblée, sans vol ; sinon l'introduction, sur le noyau entier.
  if (start >= 0) await enter(start, true);
  else {
    await enter(t === "fin" ? stops.length : -1, false, false);
    orbit = scene.overview;
  }
  dirty = true;
  requestAnimationFrame(loop);

  // — Pour la mesure (`tools/discover.mjs`) —
  (window as unknown as Record<string, unknown>).__discover = {
    ready: true,
    stops: () => stops.map(idOf),
    current: () => (cur < 0 ? "intro" : cur < stops.length ? idOf(stops[cur]) : "fin"),
    declared: () => declared,
    evidence: () => scene.evidence(),
    categories: () => categories(scene.evidence()),
    arrived: () => arrived && !anim && !scene.busy,
    arrivals: () => arrivals,
    orbit: () => orbit,
    widthUm: () => widthAt(orbit.distance, FOV),
    /** La sélection sous le centre de l'image (pixels du canevas) : ce que la caméra vise. */
    centre: () => { scene.render(orbit); const s = scene.select(canvas.width / 2, canvas.height / 2); return s && { ...s, world: null }; },
    /** La sélection sous un pixel quelconque du canevas (mesure). */
    pick: (x: number, y: number) => { scene.render(orbit); const s = scene.select(x, y); return s && { ...s, world: null }; },
    domain: () => domain,
    focus: () => focus && { i: focus.i, j: focus.j },
    region: () => scene.region,
    haltDest: () => haltDest,
    contacts: () => contactsDrawn,
    contactRatio: () => contactRatio,
    lang: LANG,
    described: () => described && { ...described, text: view.textContent },
    evidenceKnown: () => evidenceDetail(scene.evidence()).known,
    ids: () => Array.from(scene.ids()).length,
    /** Le TSS d'un gène, par l'index de la page : chrom et position 0-based. */
    resolve: async (gene: string) => {
      const ix = await index;
      const h = ix && (await ix.resolve(gene)).hits.find((x) => x.in_model);
      return h ? { chrom: h.chrom, pos: h.strand === "-" ? h.end - 1 : h.start } : null;
    },
  };
}

main().catch((e) => {
  story.innerHTML = `<p role="alert">${esc(T(UI.failed, { error: String(e) }))}</p>`;
  console.error(e);
});
