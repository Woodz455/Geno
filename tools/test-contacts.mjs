/** Lecture croisée Python → TypeScript du fichier de contacts (semaine 14).
 *
 *  Le pipeline relit le vrai fichier (`data/g3d/gm12878.contacts`) et écrit l'empreinte de
 *  tuiles tirées au hasard, décodées par **son** lecteur, et des pixels. Le lecteur
 *  TypeScript — celui que HiGlass interroge dans la page — doit rendre les mêmes, au bit
 *  près. Puis un octet corrompu doit être vu.
 */

import { build } from "esbuild";
import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/contacts-test");
mkdirSync(OUT, { recursive: true });
const FILE = resolve(ROOT, "pipeline/data/g3d/gm12878.contacts");
if (!existsSync(FILE)) {
  console.error("pas de fichier de contacts : cd pipeline && make linked");
  process.exit(2);
}
const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
const expected = resolve(OUT, "expected.json");
execFileSync(python, ["-m", "geno_pipeline", "linked", "dump", "--dump", expected], {
  cwd: resolve(ROOT, "pipeline"),
  stdio: "inherit",
  env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") },
});
const bundle = resolve(OUT, "contacts-test.mjs");
await build({ entryPoints: [resolve(ROOT, "packages/g3d/test/contacts.ts")], bundle: true, platform: "node", format: "esm",
              target: "node22", outfile: bundle, logLevel: "warning" });
const code = await new Promise((done) => {
  spawn(process.execPath, [bundle, FILE, expected], { stdio: "inherit" }).on("exit", done);
});
process.exit(code ?? 1);
