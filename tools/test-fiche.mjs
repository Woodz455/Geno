/** Lecture croisée Python → TypeScript du fichier de fiches (semaine 13).
 *
 *  Le pipeline relit le vrai fichier (`data/g3d/gm12878.fiche`) et écrit, pour des billes
 *  tirées au hasard dans chaque niveau, la fiche que **son** lecteur compose, et des
 *  positions levées en T2T par **sa** chaîne. Le lecteur TypeScript doit rendre les mêmes
 *  fiches, champ pour champ, et les mêmes positions. Puis un octet corrompu doit être vu.
 */

import { build } from "esbuild";
import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/fiche-test");
mkdirSync(OUT, { recursive: true });
const FICHE = resolve(ROOT, "pipeline/data/g3d/gm12878.fiche");
if (!existsSync(FICHE)) {
  console.error("pas de fichier de fiches : cd pipeline && make fiche");
  process.exit(2);
}
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
const expected = resolve(OUT, "expected.json");
execFileSync(python, ["-m", "geno_pipeline", "fiche", "dump", "--dump", expected], {
  cwd: resolve(ROOT, "pipeline"),
  stdio: "inherit",
  env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") },
});
const bundle = resolve(OUT, "fiche-test.mjs");
await build({ entryPoints: [resolve(ROOT, "packages/g3d/test/fiche.ts")], bundle: true, platform: "node", format: "esm",
              target: "node22", outfile: bundle, logLevel: "warning" });
const code = await new Promise((done) => {
  spawn(process.execPath, [bundle, FICHE, expected], { stdio: "inherit" }).on("exit", done);
});
process.exit(code ?? 1);
