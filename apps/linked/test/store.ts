/** Le store des vues liées, seul : ses deux défenses et sa règle de distribution.
 *  Exécuté sous Node par `tools/test-store.mjs`. */

import { type Interval, SelectionStore, type ViewName } from "../store.js";

let passed = 0;
let failed = 0;
const check = (name: string, ok: boolean, detail = "") => {
  if (ok) passed++;
  else failed++;
  console.log(`  ${ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m"}  ${name}${detail ? `  \x1b[2m${detail}\x1b[0m` : ""}`);
};
const A: Interval = { chrom: "chr17", start: 7_668_401, end: 7_687_550 };
const B: Interval = { chrom: "chr7", start: 5_527_151, end: 5_530_601 };

{
  const s = new SelectionStore();
  const got: Record<ViewName, number> = { "3d": 0, igv: 0, higlass: 0 };
  for (const v of ["3d", "igv", "higlass"] as const) s.subscribe(v, () => got[v]++);
  s.publish("igv", [A], "TP53");
  check("une sélection va aux deux autres vues, pas à sa source", got["3d"] === 1 && got.higlass === 1 && got.igv === 0, JSON.stringify(got));
  const again = s.publish("3d", [{ ...A }], "TP53 encore");
  check("la même sélection, d'une autre vue, est refusée et pas redistribuée", !again && got["3d"] === 1 && got.higlass === 1 && got.igv === 0,
        s.log.at(-1)?.reason);
  const swapped = s.publish("higlass", [B, A], "paire");
  const same = s.publish("igv", [A, B], "paire, autre ordre");
  check("deux intervalles : l'ordre ne fait pas une autre sélection", swapped && !same);
  check("numéro d'ordre : un par sélection acceptée", s.current?.seq === 2, `seq ${s.current?.seq}`);
  s.publish("3d", null);
  // igv a reçu la paire (venue de HiGlass) et l'effacement — jamais sa propre sélection.
  check("null efface, et se distribue", s.current === null && got.igv === 2 && got.higlass === 2, `igv ${got.igv}, HiGlass ${got.higlass}`);
}

{
  // Une vue naïve : elle republie ce qu'elle reçoit, un peu élargi — ce que fait une
  // bibliothèque qui émet un événement quand on la pilote. Sans la garde, la boucle tourne.
  const s = new SelectionStore();
  let applied = 0;
  s.subscribe("igv", (sel) => {
    applied++;
    if (sel) s.publish("igv", sel.regions.map((r) => ({ ...r, start: r.start - 1000, end: r.end + 1000 })), "écho élargi");
  });
  s.subscribe("higlass", () => applied++);
  s.publish("3d", [A], "bille");
  const echoes = s.log.filter((e) => e.reason === "écho").length;
  check("garde d'application : l'écho d'une vue pendant qu'elle applique est refusé", echoes === 1 && applied === 2 && s.current?.source === "3d",
        `${echoes} écho refusé, ${applied} applications`);
}

{
  const s = new SelectionStore();
  const seen: string[] = [];
  s.observe((sel) => seen.push(sel?.source ?? "—"));
  s.publish("higlass", [A]);
  s.publish("higlass", [A]);
  s.publish("3d", [B]);
  check("un observateur voit chaque sélection acceptée, source comprise, une fois", seen.join(",") === "higlass,3d", seen.join(","));
}

console.log(`\n${passed} passés, ${failed} échoués`);
process.exit(failed ? 1 : 0);
