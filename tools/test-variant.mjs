/** Lecture croisée Python → TypeScript de l'extrusion 1D et des contacts prédits (semaine 15).
 *
 *  Python fait tourner son extrusion sur l'allèle sauvage du locus EPHA4 (`geno sv dump`) ; le
 *  port TypeScript, celui que la page fait tourner en direct, doit en rendre les mêmes
 *  statistiques — densité de pieds, occupation, contacts — à l'erreur d'échantillonnage près.
 */

import { build } from "esbuild";
import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/variant-test");
mkdirSync(OUT, { recursive: true });
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
const expected = resolve(OUT, "expected.json");
if (!process.argv.includes("--reuse") || !existsSync(expected)) {
  execFileSync(python, ["-m", "geno_pipeline", "sv", "dump", "--dump", expected], {
    cwd: resolve(ROOT, "pipeline"), stdio: "inherit", env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") },
  });
}
const bundle = resolve(OUT, "variant-test.mjs");
await build({ entryPoints: [resolve(ROOT, "apps/variant/test/extrusion.ts")], bundle: true, platform: "node", format: "esm",
              target: "node22", outfile: bundle, logLevel: "warning" });
const code = await new Promise((done) => spawn(process.execPath, [bundle, expected], { stdio: "inherit" }).on("exit", done));
process.exit(code ?? 1);
