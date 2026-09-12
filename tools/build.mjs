/** Bundle esbuild → un seul HTML autonome.
 *
 *  Autonome parce que les modules ES ne se chargent pas depuis `file://` (CORS),
 *  et qu'un serveur de dev juste pour lancer un test headless est une pièce
 *  mobile de plus qui peut tomber en panne.
 */

import { build } from "esbuild";
import { mkdir, writeFile } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");

const page = (title, script, body = "") => `<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><title>${title}</title>
<style>
  :root { color-scheme: dark }
  body { margin:0; background:#0e1218; color:#e8ecf1;
         font:14px/1.6 ui-monospace,SFMono-Regular,Menlo,monospace }
  canvas { display:block }
  #results { white-space:pre-wrap; word-break:break-all; padding:12px }
</style></head><body>${body}<script>${script}</script></body></html>`;

async function bundle(entry) {
  const out = await build({
    entryPoints: [resolve(ROOT, entry)],
    bundle: true,
    format: "iife",
    target: "es2022",
    write: false,
    logLevel: "warning",
  });
  return out.outputFiles[0].text;
}

const targets = [
  { entry: "packages/viewer/test/headless.ts", out: "dist/headless.html", title: "Geno — correction des imposteurs" },
  { entry: "apps/spike/main.ts", out: "dist/spike.html", title: "Geno — spike de rendu" },
];

await mkdir(resolve(ROOT, "dist"), { recursive: true });
for (const t of targets) {
  const js = await bundle(t.entry);
  await writeFile(resolve(ROOT, t.out), page(t.title, js));
  console.log(`  ${t.out}  ${(js.length / 1024).toFixed(1)} ko`);
}
