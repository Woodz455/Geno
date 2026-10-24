/** La fiche d'information (semaine 13) : ce que l'on sait de ce qu'on vient de cliquer.
 *
 *  Deux temps, séparés pour être mesurés :
 *  - `compose` : sélection + fiche du fichier → sections et lignes de texte. Pur, synchrone,
 *    sans DOM ; c'est lui qui dit si la fiche est complète (chaque champ attendu a une
 *    ligne : une valeur, ou la raison de son absence).
 *  - `FichePanel.show` : ces lignes dans un panneau DOM, avec pour chacune le badge de son
 *    évidence — le code couleur permanent du principe n° 2 : mesuré ≠ simulé ≠ référence.
 *
 *  À l'échelle de la base, la fiche du monomère se complète de ce qui vaut pour **cette**
 *  base : sa lettre, sa position T2T exacte, les variants ClinVar qui la touchent, l'état
 *  ChromHMM et le cCRE qui la couvrent, son contexte de séquence.
 */

import type { Card, FicheField, FicheHeader, FicheReader, SequenceReader } from "@geno/g3d";
import type { Selection } from "./dive.js";

/** Palette Okabe–Ito : lisible par les trois daltonismes courants. */
export const EVIDENCE_COLORS: Record<string, string> = {
  mesure: "#009E73",
  autre: "#E69F00",
  simule: "#CC79A7",
  reference: "#56B4E9",
  annotation: "#0072B2",
};
const EVIDENCE_SHORT: Record<string, string> = {
  mesure: "mesuré GM12878",
  autre: "autre lignée",
  simule: "simulé",
  reference: "référence",
  annotation: "annotation",
};

export interface Line {
  key: string;
  label: string;
  text: string;
  evidence: string;
  /** Lignée ou tissu, pour une mesure faite ailleurs que dans GM12878. */
  cell?: string | undefined;
  source: string;
  absent: boolean;
}

export interface Composed {
  title: string;
  subtitle: string;
  sections: { name: string; lines: Line[] }[];
  /** Champs attendus sans ligne : vide = fiche complète. */
  missing: string[];
  /** Champs attendus, et combien ont une valeur (le reste : une absence motivée). */
  expected: number;
  valued: number;
}

const SECTION_NAMES: Record<string, string> = {
  coordonnees: "Coordonnées",
  structure: "Contexte structural",
  genes: "Gènes",
  regulation: "Régulation",
  clinique: "Clinique",
  sequence: "Séquence",
};

const fmt = (n: number) => n.toLocaleString("fr-FR");
const pos1 = (p: number) => fmt(p + 1);
const kb = (bp: number) => (bp >= 1e6 ? `${(bp / 1e6).toFixed(2)} Mb` : bp >= 1e3 ? `${(bp / 1e3).toFixed(1)} kb` : `${bp} pb`);
const pct = (x: number) => `${Math.round(x * 100)} %`;

function text(key: string, f: FicheField, h: FicheHeader, sel: Selection, card: Card): string {
  const v = f.v;
  // Les tables d'expression et de contrainte sont indexées par identifiant Ensembl ; on
  // affiche le nom du gène, que porte la ligne Gènes.
  const gv = card.genes?.v;
  const named = new Map<string, string>([...(gv?.genes ?? []), ...(gv?.nearest ?? [])].map((g: any) => [g.id, g.name]));
  const nameOf = (id: string) => named.get(id) ?? id;
  switch (key) {
    case "hg38":
      return `${v.chrom}:${pos1(v.start)}–${fmt(v.end)} (${kb(v.end - v.start)})`;
    case "t2t":
      return `${v.chrom}:${pos1(v.start)}–${fmt(v.end)}${v.strand === "-" ? " (retourné)" : ""} · ${pct(v.mapped)} aligné${v.split ? " · coupé entre deux chaînes" : ""}`;
    case "cytobande":
      return v.length > 3 ? `${v[0]} … ${v[v.length - 1]} (${v.length} bandes)` : v.join(", ");
    case "compartiment":
      return `${v.ab} (E1 ${v.e1.toFixed(2)})`;
    case "sous_compartiment":
      return v.slice(0, 3).map((x: any) => `${x.name} ${pct(x.frac)}`).join(" · ");
    case "tad": {
      const c = v.containing ? `dans le domaine ${pos1(v.containing[0])}–${fmt(v.containing[1])} (${kb(v.containing[1] - v.containing[0])})` : "hors domaine au centre";
      return `${v.n} domaine${v.n > 1 ? "s" : ""} · ${c} · frontière la plus proche à ${kb(v.boundary_distance)}`;
    }
    case "boucles":
      if (!v.n) return "aucune ancre HiCCUPS dans l'intervalle";
      return `${v.n} ancre${v.n > 1 ? "s" : ""} · ` + v.anchors.slice(0, 3).map((a: any) =>
        `${kb(a.anchor[1] - a.anchor[0])} ↔ ${a.partner[0]}:${pos1(a.partner[1])} (motif ${a.motif ?? "—"}/${a.partner_motif ?? "—"})`).join(" · ");
    case "radial":
      return `affichée ${v.shown.toFixed(2)}${v.mean != null ? ` · ensemble ${v.mean.toFixed(2)} ± ${v.sd.toFixed(2)} (${v.over})` : ""} (0 centre, 1 lamina)`;
    case "lad_modele":
      return `fraction LAD ${pct(v.frac)} (piste synthétique)`;
    case "lad":
      return `LMNB1 log2 ${v.log2.toFixed(2)} · ${pct(v.frac)} des bins > 0${v.frac > 0.5 ? " → LAD" : ""}`;
    case "speckles":
      return `SON ${v.score.toFixed(2)}${v.decile ? ` · décile ${v.decile}/10 du niveau` : ""}${v.decile >= 9 ? " → proche des speckles" : ""}`;
    case "nad":
      return `nucléole / contrôle log2 ${v.log2.toFixed(2)}${v.log2 > 0.5 ? " → associé au nucléole" : ""}`;
    case "replication":
      return `${v.class} (score ${v.score}, terciles ${h.repli_terciles.join(" / ")})`;
    case "ctcf_plante": {
      const s = Object.entries(v.sites);
      return s.length ? s.map(([k, o]) => `site ${k} occupation ${(o as number).toFixed(2)}`).join(" · ") : "aucun site planté sur ce monomère";
    }
    case "boucles_modele":
      return v.loops.length ? v.loops.slice(0, 4).map((l: any) => `→ ${l.partner} (${l.kind}${l.weight != null ? `, ${pct(l.weight)}` : ""})`).join(" · ") : "aucune boucle du modèle sur ce monomère";
    case "genes": {
      const show = (g: any) => `${g.name}${g.mane ? "*" : ""} (${g.strand}, ${g.type.replace(/_/g, " ")})`;
      if (!v.genes.length) {
        return "aucun gène GENCODE · plus proches codants : " + (v.nearest.map((g: any) => `${show(g)} à ${kb(g.distance)} en ${g.side}`).join(" · ") || "aucun");
      }
      const at = sel.pos >= 0 ? v.genes.filter((g: any) => g.start <= sel.pos && sel.pos < g.end) : [];
      const lead = at.length ? `à cette base : ${at.map(show).join(", ")} · ` : "";
      const pc = v.genes.filter((g: any) => g.type === "protein_coding");
      const list = (pc.length ? pc : v.genes).slice(0, 8).map(show).join(", ");
      return `${lead}${v.genes.length} gène${v.genes.length > 1 ? "s" : ""} (${pc.length} codant${pc.length > 1 ? "s" : ""}) : ${list}${(pc.length || v.genes.length) > 8 ? "…" : ""} · * MANE Select`;
    }
    case "expression": {
      const rows = Object.entries(v as Record<string, any>).sort((a, b) => b[1].lcl - a[1].lcl).slice(0, 3);
      return rows.map(([id, r]) => `${nameOf(id)} : ${r.lcl} TPM en LCL, max ${r.max} (${h.tissues[r.tissue].replace(/_/g, " ")})`).join(" · ");
    }
    case "contrainte": {
      const rows = Object.entries(v as Record<string, any>).filter(([, c]) => c.loeuf != null).sort((a, b) => a[1].loeuf - b[1].loeuf).slice(0, 3);
      return rows.map(([id, c]) => `${nameOf(id)} : LOEUF ${c.loeuf}${c.loeuf < 0.35 ? " (intolérant)" : ""} · pLI ${c.pli}`).join(" · ") || "aucune valeur publiée";
    }
    case "ccre": {
      const cl = Object.entries(v.classes as Record<string, number>).sort((a, b) => b[1] - a[1]).slice(0, 4);
      return `${v.n} actif${v.n > 1 ? "s" : ""} dans GM12878${cl.length ? " : " + cl.map(([k, n]) => `${k} ${n}`).join(", ") : ""} · ${v.inactive} inactif${v.inactive > 1 ? "s" : ""} (Low-DNase)`;
    }
    case "chromhmm":
      return v.states.slice(0, 3).map((s: any) => `${s.state} ${pct(s.frac)}`).join(" · ");
    case "ctcf":
      return `${v.n} pic${v.n > 1 ? "s" : ""} · motif + ${v.plus} · − ${v.minus} · sans motif ${v.none}`;
    case "clinvar": {
      const c = v.counts;
      const lead = v.plp.slice(0, 3).map((x: any) => {
        const d = x.diseases[0];
        return `${x.genes.join("/") || "—"} ${x.hgvs?.split(":")[1] ?? ""} ${"★".repeat(x.stars)}${d ? ` (${d[1]})` : ""}`;
      }).join(" · ");
      return `P/LP ${c.plp} · VUS ${c.vus} · B/LB ${c.blb} · contradictoires ${c.conf}${lead ? ` — ${lead}` : ""}${v.plp_total > 3 ? ` … +${v.plp_total - 3}` : ""}`;
    }
    case "gwas":
      if (!v.n) return "aucune association publiée dans l'intervalle";
      return `${v.n} association${v.n > 1 ? "s" : ""} · ` + v.top.slice(0, 3).map((x: any) => `${x.rs} ${x.trait} (p ${x.p.toExponential(0)})`).join(" · ");
    case "composition":
      return `GC ${pct(v.gc)}${v.n ? ` · N ${pct(v.n)}` : ""}`;
    case "bases":
      return "lisibles au zoom maximal (référence hg38)";
    default:
      return JSON.stringify(v).slice(0, 120);
  }
}

/** Ce qui ne vaut que pour la base cliquée, ajouté à la fiche de son monomère. */
function baseLines(sel: Selection, card: Card, reader: FicheReader | null, seq: SequenceReader | null): Line[] {
  if (sel.pos < 0) return [];
  const out: Line[] = [];
  const p = sel.pos;
  const ctx = seq?.now(sel.chrom, p - 10, p + 11);
  if (ctx) {
    const letter = ctx[10];
    const comp: Record<string, string> = { A: "T", C: "G", G: "C", T: "A", N: "N" };
    const strand = sel.strand === 1 ? ` — barreau du brin − : ${comp[letter]}` : "";
    out.push({ key: "base", label: "Base", text: `${letter} sur le brin +${strand} · ${ctx.slice(0, 10).toLowerCase()}[${letter}]${ctx.slice(11).toLowerCase()}`, evidence: "reference", source: "hg38.2bit", absent: false });
  } else {
    out.push({ key: "base", label: "Base", text: "", evidence: "reference", source: "hg38.2bit", absent: true, cell: undefined });
    out[out.length - 1].text = "séquence pas encore lue à cette échelle (elle l'est au zoom de l'hélice)";
  }
  const t = reader?.liftHs1(sel.chrom, p);
  out.push(t
    ? { key: "t2t_base", label: "Base en T2T", text: `${t.chrom}:${pos1(t.pos)}${t.strand === "-" ? " (brin retourné)" : ""}`, evidence: "reference", source: "chain.hg38.to.hs1", absent: false }
    : { key: "t2t_base", label: "Base en T2T", text: "aucun bloc de chaîne ne couvre cette base", evidence: "reference", source: "chain.hg38.to.hs1", absent: true });
  const vs = (card.clinvar?.v?.variants ?? []).filter((x: any) => x.pos === p);
  out.push({
    key: "clinvar_base", label: "ClinVar à cette base",
    text: vs.length ? vs.map((x: any) => `${x.ref}>${x.alt} ${x.clnsig} ${"★".repeat(x.stars)} (VCV ${x.id}${x.rs ? `, ${x.rs}` : ""})`).join(" · ") : "aucun variant publié à cette base",
    evidence: "annotation", source: "clinvar.vcf", absent: false,
  });
  const seg = (card.chromhmm?.v?.segments ?? []).find((s: number[]) => s[0] <= p && p < s[1]);
  out.push({ key: "chromhmm_base", label: "ChromHMM à cette base", text: seg ? `${seg[2]} (${pos1(seg[0])}–${fmt(seg[1])})` : "non segmentée", evidence: "mesure", source: "chromhmm.gm12878", absent: !seg });
  const cc = (card.ccre?.v?.list ?? []).find((c: any) => c.start <= p && p < c.end);
  out.push({ key: "ccre_base", label: "cCRE à cette base", text: cc ? `${cc.id} ${cc.class}` : "aucun élément actif ne la couvre", evidence: "mesure", source: "encode.ccres", absent: false });
  return out;
}

export function compose(sel: Selection, card: Card | null, header: FicheHeader, reader: FicheReader | null = null, seq: SequenceReader | null = null): Composed {
  const title = sel.pos >= 0
    ? `${sel.what} · ${sel.chrom}:${pos1(sel.pos)}`
    : `${sel.what} · ${sel.chrom}:${pos1(sel.start)}–${fmt(sel.end)}`;
  const subtitle = `${sel.copy} · niveau ${sel.level} · ${sel.precision}${sel.detail ? ` · ${sel.detail}` : ""}`;
  const sections: Composed["sections"] = [];
  const missing: string[] = [];
  let expected = 0;
  let valued = 0;
  if (!card) return { title, subtitle, sections, missing: ["fiche"], expected: 0, valued: 0 };
  const fine = header.levels.find((l) => l.name === sel.level)?.fine ?? false;
  for (const [sec, keys] of Object.entries(header.sections)) {
    const all = [...keys, ...(fine ? header.fine_only[sec] ?? [] : [])];
    const lines: Line[] = [];
    if (sec === "coordonnees" || sec === "sequence") lines.push(...(sec === "sequence" ? baseLines(sel, card, reader, seq).filter((l) => l.key === "base") : []));
    for (const key of all) {
      expected++;
      const f = card[key];
      if (!f || (f.v === undefined && !f.absent)) {
        missing.push(key);
        continue;
      }
      const absent = f.v === undefined;
      if (!absent) valued++;
      lines.push({
        key, label: header.labels[key] ?? key, evidence: f.e, cell: f.cell, source: f.s, absent,
        text: absent ? f.absent! : text(key, f, header, sel, card),
      });
    }
    if (sec === "coordonnees") lines.push(...baseLines(sel, card, reader, seq).filter((l) => l.key === "t2t_base"));
    if (sec === "regulation") lines.push(...baseLines(sel, card, reader, seq).filter((l) => l.key === "chromhmm_base" || l.key === "ccre_base"));
    if (sec === "clinique") lines.push(...baseLines(sel, card, reader, seq).filter((l) => l.key === "clinvar_base"));
    sections.push({ name: SECTION_NAMES[sec] ?? sec, lines });
  }
  return { title, subtitle, sections, missing, expected, valued };
}

const CSS = `
  #fiche { position: fixed; right: 12px; top: 10px; bottom: 10px; width: min(420px, calc(100vw - 24px));
           overflow-y: auto; background: rgba(14,18,24,.9); border-radius: 8px; padding: 10px 12px;
           font-size: 12px; line-height: 1.35; display: none; }
  #fiche.open { display: block; }
  #fiche h2 { font-size: 13px; margin: 0; font-weight: 600; }
  #fiche .sub { color: #8b98aa; margin: 2px 0 6px; }
  #fiche details { margin: 4px 0; }
  #fiche summary { cursor: pointer; font-weight: 600; color: #dbe3ee; }
  #fiche .row { display: grid; grid-template-columns: 112px 1fr; gap: 6px; margin: 3px 0; }
  #fiche .lab { color: #8b98aa; }
  #fiche .abs { color: #8b98aa; font-style: italic; }
  #fiche .ev { display: inline-block; font-size: 10px; padding: 0 5px; border-radius: 8px; margin-right: 4px;
               color: #0e1218; font-weight: 600; white-space: nowrap; }
  #fiche .legend { margin-top: 8px; color: #8b98aa; }
  #fiche button.close { float: right; font: inherit; background: none; border: 1px solid #3a4658; color: inherit;
                        border-radius: 6px; cursor: pointer; padding: 0 6px; }
  #fiche .wait { color: #E69F00; }
`;

export class FichePanel {
  readonly el: HTMLElement;
  /** Appelé quand l'utilisateur ferme la fiche. */
  onHide: (() => void) | null = null;

  constructor(doc: Document) {
    const st = doc.createElement("style");
    st.textContent = CSS;
    doc.head.appendChild(st);
    this.el = doc.createElement("aside");
    this.el.id = "fiche";
    this.el.setAttribute("aria-live", "polite");
    this.el.setAttribute("aria-label", "Fiche d'information du point cliqué");
    doc.body.appendChild(this.el);
  }

  hide(): void {
    this.el.classList.remove("open");
  }

  /** Écrit la fiche. `waiting` : le bloc de la fiche n'est pas encore là ; les coordonnées
   *  s'affichent tout de suite, le reste quand il arrive. */
  show(c: Composed, header: FicheHeader, waiting = false): void {
    const esc = (s: string) => s.replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[ch]!);
    const badge = (l: Line) => {
      const src = header.sources[l.source];
      const tip = `${header.evidence[l.evidence] ?? l.evidence}${l.cell ? ` — ${l.cell}` : ""} · source : ${l.source}${src ? ` (${src.license}) ${src.url}` : ""}`;
      return `<span class="ev" style="background:${EVIDENCE_COLORS[l.evidence] ?? "#999"}" title="${esc(tip)}">${esc(EVIDENCE_SHORT[l.evidence] ?? l.evidence)}${l.cell ? ` · ${esc(l.cell.split(" ")[0])}` : ""}</span>`;
    };
    const parts = [
      `<button class="close" aria-label="Fermer la fiche">×</button>`,
      `<h2>${esc(c.title)}</h2><div class="sub">${esc(c.subtitle)}</div>`,
    ];
    if (waiting) parts.push(`<div class="wait">fiche en cours de lecture…</div>`);
    for (const s of c.sections) {
      parts.push(`<details open><summary>${esc(s.name)}</summary>`);
      for (const l of s.lines) {
        parts.push(`<div class="row"><div class="lab">${esc(l.label)}</div><div>${badge(l)}<span class="${l.absent ? "abs" : ""}">${esc(l.text)}</span></div></div>`);
      }
      parts.push("</details>");
    }
    if (c.missing.length && !waiting) parts.push(`<div class="wait">champs sans valeur ni raison : ${esc(c.missing.join(", "))}</div>`);
    parts.push(`<div class="legend">${Object.entries(EVIDENCE_SHORT).map(([k, v]) => `<span class="ev" style="background:${EVIDENCE_COLORS[k]}">${v}</span>`).join(" ")}</div>`);
    this.el.innerHTML = parts.join("");
    this.el.querySelector("button.close")!.addEventListener("click", () => {
      this.hide();
      this.onHide?.();
    });
    this.el.classList.add("open");
  }
}
