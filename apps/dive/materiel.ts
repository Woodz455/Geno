/** Banc matériel : la plongée mesurée sur l'appareil de la personne qui ouvre la page.
 *
 *  Publiée comme Artifact (privée par défaut). Le `.g3d` voyage **dans** la page, en base64 :
 *  les Artifacts ne servent pas de fichier binaire à côté d'une page. Chaque mesure est
 *  enregistrée dans la base de l'Artifact (capacité `db`), où Claude la relit : c'est le
 *  seul chemin par lequel un temps par image mesuré sur un vrai GPU revient au projet.
 *  Sans base (page ouverte hors de claude.ai), la mesure s'affiche et n'est pas enregistrée.
 */

import { BufferSource } from "@geno/g3d";
import { distanceFor, type Orbit } from "@geno/viewer";
import { median, percentile, runBench } from "./bench.js";
import { FOV, openDive } from "./dive.js";

interface Doc {
  id: string;
  data(): Record<string, unknown> | undefined;
}
interface Db {
  collection(path: string): {
    add(data: Record<string, unknown>): Promise<unknown>;
    orderBy(field: string, dir?: "asc" | "desc"): {
      limit(n: number): { onSnapshot(next: (s: { docs: Doc[] }) => void, error?: (e: { code: string }) => void): () => void };
    };
  };
}

declare global {
  interface Window {
    claude?: { use(name: string): Promise<unknown> };
  }
}

const FRAMES = 240;
const $ = <T extends HTMLElement = HTMLElement>(id: string) => document.getElementById(id) as T;
const canvas = $<HTMLCanvasElement>("scene");
const status = $("status");
const button = $<HTMLButtonElement>("run");

function fit(): [number, number] {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(320, Math.round(canvas.clientWidth * dpr));
  const h = Math.round((w * 9) / 16);
  canvas.width = w;
  canvas.height = h;
  return [w, h];
}

const fmt = (ms: number) => `${ms.toFixed(1)} ms`;
const verdictOf = (p95: number) => p95 <= 20;

function renderRuns(docs: Doc[]): void {
  const body = $("runs");
  if (!docs.length) return;
  body.replaceChildren(
    ...docs.map((d) => {
      const r = d.data() ?? {};
      const tr = document.createElement("tr");
      const cell = (text: string, cls = "") => {
        const td = document.createElement("td");
        td.textContent = text;
        if (cls) td.className = cls;
        tr.appendChild(td);
      };
      const num = (k: string) => (typeof r[k] === "number" ? fmt(r[k] as number) : "–");
      // Les mesures d'avant la semaine 11 ne portent pas d'étiquette : elles s'arrêtaient aux
      // compartiments.
      cell(String(r.scene || "semaine 10"));
      cell(String(r.device || "—"));
      cell(String(r.gpu ?? "?"));
      cell(`${r.width}×${r.height}`);
      cell(num("raf_ms_median"), "n");
      cell(num("raf_ms_p95"), "n");
      cell(num("frame_ms_median"), "n");
      const ok = typeof r.raf_ms_p95 === "number" && verdictOf(r.raf_ms_p95 as number);
      const td = document.createElement("td");
      const pill = document.createElement("span");
      pill.className = `pill ${ok ? "ok" : "ko"}`;
      pill.textContent = ok ? "tenu" : "non tenu";
      td.appendChild(pill);
      tr.appendChild(td);
      return tr;
    }),
  );
}

async function main(): Promise<void> {
  fit();
  const embedded = document.getElementById("g3d")?.textContent?.trim();
  if (!embedded) throw new Error("le fichier du modèle manque à la page");
  const bytes = Uint8Array.from(atob(embedded), (c) => c.charCodeAt(0));
  const scene = await openDive(canvas, new BufferSource(bytes));
  $("cap").textContent = `${scene.preview.frame.count.toLocaleString("fr")} puis ${(scene.fine?.frame.count ?? 0).toLocaleString("fr")} billes · ${canvas.width}×${canvas.height}`;

  // Arrivée du noyau en fondu, puis une rotation lente : la scène est vivante sans être mesurée.
  const t0 = performance.now();
  let running = false;
  const idle = (now: number) => {
    if (!running) {
      scene.arrival = Math.min(1, (now - t0) / 400);
      const o: Orbit = { ...scene.overview, yaw: scene.overview.yaw + (now - t0) / 9000 };
      scene.render(o);
    }
    requestAnimationFrame(idle);
  };
  requestAnimationFrame(idle);
  button.disabled = false;

  // La base, si la page est servie par claude.ai. Absente : on mesure quand même.
  let db: Db | null = null;
  try {
    db = ((await window.claude?.use("db")) as Db | null) ?? null;
  } catch {
    db = null;
  }
  if (db) {
    db.collection("runs").orderBy("at", "desc").limit(50).onSnapshot(
      (s) => renderRuns(s.docs),
      () => (status.textContent = "Le tableau des mesures est indisponible pour cette visite."),
    );
  } else {
    status.textContent = "Mesure possible ; elle ne sera pas enregistrée (pas de base de données dans cette vue).";
  }

  button.addEventListener("click", async () => {
    button.disabled = true;
    running = true;
    try {
      const [w, h] = fit();
      scene.resize(w, h);
      const r = await runBench(scene, FRAMES, {
        onProgress: (f, phase) => (status.textContent = `Mesure en cours, passe ${phase} : ${Math.round(f * 100)} %`),
      });
      const raf = median(r.rafMs);
      const rafP95 = percentile(r.rafMs, 0.95);
      const dropped = r.rafMs.filter((x) => x > 25).length / r.rafMs.length;
      const sync = median(r.frameMs);
      const ok = verdictOf(rafP95);
      $("raf").textContent = fmt(raf);
      $("rafp95").textContent = fmt(rafP95);
      $("dropped").textContent = `${(dropped * 100).toFixed(1)} %`;
      $("sync").textContent = fmt(sync);
      $("gpu").textContent = `${r.gpu} · canevas ${w}×${h} · ${FRAMES} images par passe`;
      const v = $("verdict");
      v.className = `pill ${ok ? "ok" : "ko"}`;
      v.textContent = ok ? "60 fps tenus" : "60 fps non tenus";
      $("result").hidden = false;

      const record = {
        at: new Date().toISOString(),
        scene: "semaine 11",
        device: ($<HTMLInputElement>("device").value || "").slice(0, 80),
        gpu: r.gpu.slice(0, 200),
        ua: navigator.userAgent.slice(0, 300),
        width: w,
        height: h,
        dpr: window.devicePixelRatio || 1,
        frames: FRAMES,
        instances: r.instances,
        raf_ms_median: raf,
        raf_ms_p95: rafP95,
        raf_dropped_share: dropped,
        frame_ms_median: sync,
        frame_ms_p95: percentile(r.frameMs, 0.95),
        // Les séries complètes, arrondies : de quoi refaire l'analyse sans refaire la mesure.
        raf_ms: r.rafMs.map((x) => Math.round(x * 10) / 10),
        frame_ms: r.frameMs.map((x) => Math.round(x * 10) / 10),
      };
      if (db) {
        await db.collection("runs").add(record);
        status.textContent = "Mesure enregistrée. Vous pouvez la refaire sur un autre appareil en rouvrant ce lien.";
      } else {
        status.textContent = "Mesure terminée (non enregistrée).";
      }
    } catch (e) {
      status.textContent = `La mesure a échoué : ${e instanceof Error ? e.message : String(e)}`;
    } finally {
      running = false;
      button.disabled = false;
      scene.render({ ...scene.overview, distance: distanceFor(24, FOV) });
    }
  });
}

main().catch((e) => {
  status.textContent = `Impossible de démarrer : ${e instanceof Error ? e.message : String(e)}`;
});
