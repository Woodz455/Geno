/** Lecture croisée Python → TypeScript de l'index de recherche (semaine 16).
 *
 *  Le pipeline résout les requêtes des critères — gènes, positions, rsID, et leurs témoins —
 *  par son lecteur ; le lecteur TypeScript, celui de la page, doit rendre les mêmes cibles.
 *  Puis un octet corrompu doit être vu.
 */

import { build } from "esbuild";
import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/search-test");
mkdirSync(OUT, { recursive: true });
const FILE = resolve(ROOT, "pipeline/data/g3d/gm12878.search");
if (!existsSync(FILE)) {
  console.error("pas d'index de recherche : cd pipeline && make search");
  process.exit(2);
}
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
const expected = resolve(OUT, "expected.json");
execFileSync(python, ["-m", "geno_pipeline", "search", "dump", "--dump", expected], {
  cwd: resolve(ROOT, "pipeline"),
  stdio: "inherit",
  env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") },
});
const bundle = resolve(OUT, "search-test.mjs");
await build({ entryPoints: [resolve(ROOT, "packages/g3d/test/search.ts")], bundle: true, platform: "node", format: "esm",
              target: "node22", outfile: bundle, logLevel: "warning" });
const code = await new Promise((done) => {
  spawn(process.execPath, [bundle, FILE, expected], { stdio: "inherit" }).on("exit", done);
});
process.exit(code ?? 1);
