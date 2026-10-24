/** Lecture croisée Python → TypeScript du format `.g3d`.
 *
 *  1. le pipeline écrit un fichier témoin et ce que son propre lecteur en relit ;
 *  2. esbuild empaquette le test TypeScript (les imports `.js` du dépôt ne passent pas
 *     par le simple retrait de types de Node 22) ;
 *  3. Node l'exécute. `DecompressionStream` et `crypto.subtle` y sont natifs.
 */

import { build } from "esbuild";
import { execFileSync, spawn } from "node:child_process";
import { existsSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { serve } from "./serve.mjs";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/g3d-fixture");
mkdirSync(OUT, { recursive: true });

const venv = resolve(ROOT, "pipeline/.venv/bin/python");
const python = existsSync(venv) ? venv : "python3";
execFileSync(python, ["-m", "geno_pipeline.export.fixture", OUT], {
  cwd: resolve(ROOT, "pipeline"),
  stdio: "inherit",
  env: { ...process.env, PYTHONPATH: resolve(ROOT, "pipeline") },
});

const bundle = resolve(OUT, "node-test.mjs");
await build({
  entryPoints: [resolve(ROOT, "packages/g3d/test/node.ts")],
  bundle: true,
  platform: "node",
  format: "esm",
  target: "node22",
  outfile: bundle,
  logLevel: "warning",
});

// Le test HTTP passe par le vrai serveur : `/` y sert `dist/`, donc le témoin.
const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const code = await new Promise((done) => {
  const child = spawn(process.execPath, [bundle, resolve(OUT, "fixture.g3d"), resolve(OUT, "expected.json"), base], {
    stdio: "inherit",
  });
  child.on("exit", done);
});
server.close();
process.exit(code ?? 1);
