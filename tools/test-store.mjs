/** Le store des vues liées (semaine 14), sous Node : `apps/linked/test/store.ts`. */

import { build } from "esbuild";
import { spawn } from "node:child_process";
import { mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const OUT = resolve(ROOT, "dist/store-test");
mkdirSync(OUT, { recursive: true });
const bundle = resolve(OUT, "store-test.mjs");
await build({ entryPoints: [resolve(ROOT, "apps/linked/test/store.ts")], bundle: true, platform: "node", format: "esm",
              target: "node22", outfile: bundle, logLevel: "warning" });
const code = await new Promise((done) => spawn(process.execPath, [bundle], { stdio: "inherit" }).on("exit", done));
process.exit(code ?? 1);
