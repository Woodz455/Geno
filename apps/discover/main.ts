/** Geno — mode Découverte (semaine 19) : « Du noyau à une lettre ».
 *
 *  La plongée de la semaine 12 (`openDive`, cible TP53), avec un préréglage et une couche
 *  narrative : sept étapes du noyau à une base, quatre haltes commentées, un glossaire au survol et
 *  au focus, des échelles comparées calculées (`scale.ts`), une fin de parcours avec une question
 *  ouverte et trois questions à choix. Aucun moteur nouveau : les arrêts sont ceux de la scène,
 *  les haltes passent par `destination()` comme la recherche de la semaine 16, et le badge
 *  d'évidence lit `scene.evidence()` — le texte ne peut pas dire « mesuré » quand la couche
 *  affichée est simulée sans que la mesure (`tools/discover.mjs`) le voie.
 *
 *  Le contenu est hors du code (`content/fr.json`, `content/glossary.fr.json`) : l'anglais de la
 *  semaine 20 est une traduction de fichiers, pas une réécriture.
 *
 *  URL : `#t=<étape ou halte>` reprend le parcours où il était.
 *  Paramètres : f=… (.g3d), w=, h= (taille fixe, pour la mesure), motion=0 (sans vol).
 */

import { ContactsReader, FicheReader, HttpSource, SearchIndex } from "@geno/g3d";
import { distanceFor, dive, type Orbit, widthAt } from "@geno/viewer";
import { type DiveScene, FOV, openDive } from "../dive/dive.js";
import { destination } from "../dive/navigate.js";
import content from "./content/fr.json";
import glossary from "./content/glossary.fr.json";
import { fieldCompare, fieldText, fill } from "./scale.js";

type Step = (typeof content.steps)[number];
type Halt = (typeof content.halts)[number];
type Stop = { kind: "step"; index: number; step: Step } | { kind: "halt"; index: number; halt: Halt };

const q = new URLSearchParams(location.search);
const FILE = q.get("f") ?? "data/gm12878.g3d";
const CONTACTS = FILE.replace(/\.g3d$/, ".contacts");
const SEARCH = FILE.replace(/\.g3d$/, ".search");
const FICHE = FILE.replace(/\.g3d$/, ".fiche");
const FIXED = q.has("w");
const REDUCED = q.get("motion") === "0" || matchMedia("(prefers-reduced-motion: reduce)").matches;

// — Style —
const style = document.createElement("style");
style.textContent = `
  :root { --bg: #0b0f14; --panel: rgba(14,18,24,.92); --ink: #e8edf3; --mute: #9aa7b8; --line: #2b3646;
          --accent: #f0e342; --link: #9cc3ff; }
  html, body { height: 100%; margin: 0; overflow: hidden; background: var(--bg); color: var(--ink);
               font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }
  canvas { display: block; touch-action: none; }
  #overlay { position: fixed; left: 0; top: 0; pointer-events: none; }
  #story { position: fixed; right: 0; top: 0; bottom: 0; width: min(400px, 100vw); box-sizing: border-box;
           background: var(--panel); border-left: 1px solid var(--line); padding: 18px 20px 14px;
           display: flex; flex-direction: column; gap: 10px; overflow-y: auto; }
  #story header .kicker { color: var(--mute); font-size: 12px; letter-spacing: .04em; text-transform: uppercase; }
  #story h1 { font-size: 15px; margin: 2px 0 0; font-weight: 600; color: var(--mute); }
  #story h2 { font-size: 22px; line-height: 1.25; margin: 4px 0 0; }
  #story p { margin: 0; }
  #story .show { color: var(--mute); font-size: 13.5px; }
  #story .summary { border-left: 3px solid var(--accent); padding: 4px 0 4px 10px; font-weight: 500; }
  .term { text-decoration: underline dotted 1.5px; text-underline-offset: 3px; cursor: help; color: inherit; }
  .term:focus-visible { outline: 2px solid var(--link); outline-offset: 2px; border-radius: 2px; }
  #tip { position: fixed; z-index: 10; max-width: 300px; background: #1b2330; color: var(--ink);
         border: 1px solid #3a4658; border-radius: 8px; padding: 8px 10px; font-size: 13.5px; line-height: 1.45;
         pointer-events: none; display: none; box-shadow: 0 6px 20px rgba(0,0,0,.4); }
  #tip b { display: block; margin-bottom: 2px; }
  .badge { display: inline-flex; gap: 6px; align-items: baseline; font-size: 12.5px; padding: 3px 9px;
           border-radius: 999px; border: 1px solid var(--line); background: #121822; }
  .badge .dot { width: 8px; height: 8px; border-radius: 50%; display: inline-block; align-self: center; }
  #evidence-detail { color: var(--mute); font-size: 12px; }
  #ruler { display: flex; align-items: center; gap: 10px; font-size: 13px; color: var(--mute); }
  #ruler .bar { flex: 0 0 auto; height: 3px; background: var(--ink); width: 80px; position: relative; }
  #ruler .bar::before, #ruler .bar::after { content: ""; position: absolute; top: -4px; width: 2px; height: 11px; background: var(--ink); }
  #ruler .bar::before { left: 0; } #ruler .bar::after { right: 0; }
  #ruler b { color: var(--ink); font-weight: 600; }
  #contacts { width: 100%; height: auto; border-radius: 6px; background: #0e131b; }
  #contacts-caption { color: var(--mute); font-size: 12px; }
  nav.steps { display: flex; gap: 8px; margin-top: auto; padding-top: 8px; }
  button { font: inherit; padding: 7px 14px; border-radius: 8px; border: 1px solid #3a4658;
           background: #1b2330; color: inherit; cursor: pointer; }
  button:hover { background: #243042; }
  button:focus-visible, input:focus-visible, textarea:focus-visible { outline: 2px solid var(--link); outline-offset: 2px; }
  button.primary { background: #2c4a73; border-color: #3d6396; }
  button[disabled] { opacity: .45; cursor: default; }
  .dots { display: flex; gap: 5px; align-items: center; }
  .dots span { width: 8px; height: 8px; border-radius: 50%; background: #3a4658; }
  .dots span.on { background: var(--accent); }
  #halts { display: flex; flex-wrap: wrap; gap: 6px; }
  #halts button { font-size: 13px; padding: 4px 10px; }
  #halts-label { color: var(--mute); font-size: 12px; }
  .quiz fieldset { border: 1px solid var(--line); border-radius: 8px; margin: 0 0 8px; padding: 8px 10px; }
  .quiz legend { padding: 0 4px; font-weight: 500; }
  .quiz label { display: block; padding: 2px 0; cursor: pointer; }
  .quiz .verdict { font-size: 13px; margin-top: 4px; }
  .quiz textarea { width: 100%; box-sizing: border-box; min-height: 90px; font: inherit; color: inherit;
                   background: #0e131b; border: 1px solid #3a4658; border-radius: 6px; padding: 6px 8px; }
  .sr { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap; }
  @media (max-width: 760px) { #story { top: auto; height: 55vh; width: 100vw; border-left: 0; border-top: 1px solid var(--line); } }
`;
document.head.appendChild(style);
document.documentElement.lang = "fr";
document.title = `Geno — ${content.title}`;

// — Squelette —
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
story.setAttribute("aria-label", content.title);
document.body.appendChild(story);
const tip = document.createElement("div");
tip.id = "tip";
tip.setAttribute("role", "tooltip");
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
const GLOSS = glossary as Record<string, string>;
const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]!);
let tipFor: HTMLElement | null = null;
function rich(text: string): string {
  return esc(fill(text)).replace(/\[\[([^\]|]+)(?:\|([^\]]+))?\]\]/g, (_, term: string, shown?: string) => {
    if (!(term in GLOSS)) throw new Error(`terme absent du glossaire : ${term}`);
    return `<span class="term" tabindex="0" role="button" aria-describedby="tip" data-term="${esc(term)}">${shown ?? term}</span>`;
  });
}
function showTip(el: HTMLElement): void {
  const term = el.dataset.term!;
  tip.innerHTML = `<b>${esc(term)}</b>${esc(GLOSS[term])}`;
  tip.style.display = "block";
  const r = el.getBoundingClientRect();
  const w = Math.min(300, window.innerWidth - 16);
  tip.style.left = `${Math.max(8, Math.min(window.innerWidth - w - 8, r.left))}px`;
  const below = r.bottom + 8, h = tip.offsetHeight;
  tip.style.top = `${below + h < window.innerHeight ? below : Math.max(8, r.top - h - 8)}px`;
  tipFor = el;
}
function hideTip(): void {
  tip.style.display = "none";
  tipFor = null;
}
story.addEventListener("mouseover", (e) => { const t = (e.target as HTMLElement).closest<HTMLElement>(".term"); if (t) showTip(t); });
story.addEventListener("mouseout", (e) => { if ((e.target as HTMLElement).closest(".term")) hideTip(); });
story.addEventListener("focusin", (e) => { const t = (e.target as HTMLElement).closest<HTMLElement>(".term"); if (t) showTip(t); else hideTip(); });
story.addEventListener("focusout", () => hideTip());

// — Évidence : ce que la couche affichée est vraiment —
const EVIDENCE_COLORS: Record<string, string> = { "mesuré": "#5ec27a", "simulé": "#e0a84a", "prédit": "#c58be0", "illustratif": "#7aa8e6" };
/** Les catégories que couvre le texte d'évidence de la scène (`DiveScene.evidence`). */
export function categories(evidence: string): string[] {
  if (evidence.startsWith("séquence")) return ["mesuré"];
  if (evidence.startsWith("nucléosomes")) return ["illustratif"];
  if (evidence.startsWith("fibre")) return ["simulé"];
  return ["simulé", "prédit"];
}

async function main(): Promise<void> {
  const [w0, h0] = fit();
  story.innerHTML = `<header><div class="kicker">Geno · mode Découverte</div><h1>${esc(content.title)}</h1></header><p>Chargement du noyau…</p>`;
  const scene: DiveScene = await openDive(canvas, FILE, { target: "tp53" });
  scene.resize(w0, h0);
  const index = SearchIndex.open(new HttpSource(SEARCH)).catch(() => null);
  const fiche = FicheReader.open(new HttpSource(FICHE)).catch(() => null);
  const contacts = ContactsReader.open(new HttpSource(CONTACTS)).catch(() => null);

  // Les arrêts des étapes : ceux de la plongée de la semaine 12, et trois vues réglées sur le
  // calendrier des fondus (`SCHEDULE`, scene.ts) pour que chaque étape montre sa couche.
  const S = scene.stops;
  const at = (o: Orbit, widthUm: number): Orbit => ({ ...o, distance: distanceFor(widthUm, FOV) });
  const focus = scene.fiber?.focus ?? null;
  const orbits: Record<string, Orbit> = {
    noyau: scene.overview,
    territoires: at(scene.final, 16),             // les autres chromosomes encore là (20 → 8 µm), assombris
    compartiments: at(scene.final, 3.2),          // compartiments pleins (< 3,5 µm), pas encore de fibre (> 3 µm)
    tad: scene.loopView ?? S[2],
    boucle: scene.anchorView ? at(scene.anchorView, 0.3) : S[2],   // l'ancre, avant les nucléosomes (< 0,25 µm)
    fibre: S[3],
    helice: S[S.length - 1],
  };
  const chr17 = scene.region?.chrom ?? "chr17";
  const level = scene.region ? scene.g3d.level(scene.region.level) : null;
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
  let anim: { path: { length: number; at(s: number): Orbit }; t0: number; duration: number } | null = null;
  let dirty = true;
  const arrivals: { id: string; orbit: Orbit; arrivedAt: number; dest?: unknown }[] = [];
  let arrived = true;

  const flyTo = (to: Orbit, jump = false): void => {
    arrived = false;
    if (REDUCED || jump) {
      orbit = to;
      anim = null;
      dirty = true;
      return;
    }
    const path = dive(orbit, to, FOV);
    anim = { path, t0: performance.now(), duration: 900 * Math.max(0.6, path.length) };
    dirty = true;
  };

  // — Rendu du panneau —
  const dots = () => `<div class="dots" aria-hidden="true">${content.steps.map((_, k) => `<span class="${k <= cur ? "on" : ""}"></span>`).join("")}</div>`;
  const haltsHtml = () => `<div id="halts-label">Haltes commentées</div><div id="halts">${content.halts.map((h) => `<button type="button" data-halt="${h.id}">${esc(h.title.split(",")[0])}</button>`).join("")}</div>`;
  const evidenceHtml = () => `<div><span class="badge" id="badge"><span class="dot"></span><span id="badge-text"></span></span></div><div id="evidence-detail"></div>`;
  const rulerHtml = () => `<div id="ruler" aria-live="off"><span class="bar"></span><span id="ruler-text"></span></div>`;
  const nav = (prevLabel: string, nextLabel: string, prevOff = false, nextOff = false) =>
    `<nav class="steps" aria-label="Parcours"><button type="button" id="prev" ${prevOff ? "disabled" : ""}>← ${prevLabel}</button><button type="button" class="primary" id="next" ${nextOff ? "disabled" : ""}>${nextLabel} →</button></nav>`;

  let declared = "";
  const render = (): void => {
    hideTip();
    let html = `<header><div class="kicker">Geno · mode Découverte</div><h1>${esc(content.title)}</h1></header>`;
    if (cur < 0) {
      declared = "";
      html += `<h2>${esc(content.subtitle)}</h2><p>${rich(content.intro)}</p>${rulerHtml()}${haltsHtml()}${nav("Retour", "Commencer", true)}`;
    } else if (cur < content.steps.length) {
      const st = content.steps[cur];
      declared = st.evidence;
      html += `<div class="kicker">Étape ${cur + 1} sur ${content.steps.length}</div>${dots()}<h2>${esc(st.title)}</h2><p>${rich(st.text)}</p>`;
      if ("how" in st && st.how) html += `<p>${rich(st.how)}</p>`;
      if ("contacts" in st && st.contacts) html += `<canvas id="contacts" width="360" height="190" role="img" aria-label="Carte de contacts de la région : un triangle sombre marque le TAD"></canvas><div id="contacts-caption">Carte de contacts du modèle autour du TAD : plus c'est sombre, plus les deux morceaux se touchent. Les traits marquent ses frontières.</div>`;
      if ("summary" in st && st.summary) html += `<p class="summary">${rich(st.summary)}</p>`;
      html += `<p class="show">${rich(st.show)}</p>${evidenceHtml()}${rulerHtml()}${nav("Précédent", cur === content.steps.length - 1 ? "Terminer" : "Suivant", false)}${haltsHtml()}`;
    } else if (cur >= content.steps.length && cur < stops.length) {
      const h = (stops[cur] as Extract<Stop, { kind: "halt" }>).halt;
      declared = "simulé";
      html += `<div class="kicker">Halte</div><h2>${esc(h.title)}</h2><p>${rich(h.text)}</p><p class="show"><b>Ce que le modèle ne montre pas.</b> ${rich(h.limits)}</p><div id="halt-facts" class="show"></div>${evidenceHtml()}${rulerHtml()}${nav("Retour au parcours", "Fin du parcours")}${haltsHtml()}`;
    } else {
      declared = "";
      const Q = content.quiz;
      html += `<div class="kicker">Fin du parcours</div><h2>${esc(Q.intro)}</h2><form class="quiz" id="quiz"><fieldset><legend>${esc(Q.open)}</legend><label class="sr" for="open">Votre réponse</label><textarea id="open"></textarea><div class="show">${esc(Q.openNote)}</div></fieldset>` +
        Q.choices.map((c, k) => `<fieldset><legend>${esc(c.q)}</legend>${c.options.map((o, j) => `<label><input type="radio" name="q${k}" value="${j}"> ${esc(o)}</label>`).join("")}<div class="verdict" id="v${k}" aria-live="polite"></div></fieldset>`).join("") +
        `<button type="submit" class="primary">Vérifier</button></form>${haltsHtml()}${nav("Revoir la dernière étape", "Recommencer")}`;
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
      const saved = (() => { try { return localStorage.getItem("geno.decouverte.reponse") ?? ""; } catch { return ""; } })();
      const ta = quiz.querySelector<HTMLTextAreaElement>("#open")!;
      ta.value = saved;
      ta.addEventListener("input", () => { try { localStorage.setItem("geno.decouverte.reponse", ta.value); } catch { /* navigation privée */ } });
      quiz.addEventListener("submit", (e) => {
        e.preventDefault();
        content.quiz.choices.forEach((c, k) => {
          const v = quiz.querySelector<HTMLInputElement>(`input[name=q${k}]:checked`);
          const out = quiz.querySelector<HTMLDivElement>(`#v${k}`)!;
          out.textContent = !v ? "Pas de réponse." : Number(v.value) === c.answer ? `Juste. ${c.explain}` : `Pas tout à fait. ${c.explain}`;
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
      badge.textContent = `ce que vous voyez : ${shown}`;
      story.querySelector<HTMLSpanElement>("#badge .dot")!.style.background = EVIDENCE_COLORS[shown] ?? "#888";
      story.querySelector<HTMLDivElement>("#evidence-detail")!.textContent = ev;
    }
    const ruler = story.querySelector<HTMLSpanElement>("#ruler-text");
    if (ruler) {
      const width = widthAt(orbit.distance, FOV);
      const bar = width * (80 / Math.max(1, canvas.clientWidth));       // la barre fait 80 px de l'image
      ruler.innerHTML = `<b>${fieldText(bar)}</b> · image : ${fieldText(width)}, ${fieldCompare(width)}`;
    }
  };

  // — La carte de contacts de l'étape TAD : triangle, diagonale en bas —
  let contactsDrawn: { bins: number; resolution: number; domain: [number, number] } | null = null;
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
    const fmt = (x: number) => x.toLocaleString("fr-FR");
    out.innerHTML = `Dans le modèle : bille de ${fmt((d.end - d.start) / 1000)} kb, ${esc(d.copy)}, ${fmt(d.start + 1)}–${fmt(d.end)} (hg38)` +
      (genes.length ? ` · ${genes.length} gène${genes.length > 1 ? "s" : ""} dont ${esc(genes.slice(0, 6).join(", "))}${genes.length > 6 ? "…" : ""}` : "");
  };

  // — Navigation —
  const enter = async (k: number, jump = false): Promise<void> => {
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
    try { history.replaceState(null, "", `${location.pathname}${location.search}#t=${k < 0 ? "intro" : k < stops.length ? idOf(stops[k]) : "fin"}`); } catch { /* fichier local */ }
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
    const d = await destination(scene, { chrom: hit.chrom, start: hit.start, end: hit.end, pos }, orbit);
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
    else if (e.key === "Escape") { hideTip(); if (cur >= content.steps.length && cur < stops.length) back(); }
    else if (e.key === "Enter" && t.classList.contains("term")) { tipFor === t ? hideTip() : showTip(t); }
  });
  window.addEventListener("resize", () => {
    if (FIXED) return;
    scene.resize(...fit());
    dirty = true;
  });

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
      orbit = anim.path.at(k * anim.path.length);
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
    if (!anim && !arrived && !scene.busy) {
      arrived = true;
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
    await enter(t === "fin" ? stops.length : -1);
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
    domain: () => domain,
    focus: () => focus && { i: focus.i, j: focus.j },
    region: () => scene.region,
    haltDest: () => haltDest,
    contacts: () => contactsDrawn,
    ids: () => Array.from(scene.ids()).length,
  };
}

main().catch((e) => {
  story.innerHTML = `<p>La page n'a pas pu s'ouvrir : ${esc(String(e))}</p>`;
  console.error(e);
});
