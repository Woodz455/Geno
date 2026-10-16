/** Bundle esbuild → un seul HTML autonome.
 *
 *  Autonome parce que les modules ES ne se chargent pas depuis `file://` (CORS),
 *  et qu'un serveur de dev juste pour lancer un test headless est une pièce
 *  mobile de plus qui peut tomber en panne.
 */

import { build } from "esbuild";
import { copyFile, mkdir, readFile, writeFile } from "node:fs/promises";
import { existsSync } from "node:fs";
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

/** Retire les commentaires des shaders GLSL. La minification ne touche pas au contenu des
 *  chaînes, et les shaders commentés pesaient 2 Kio dans la page du premier rendu — dont
 *  chaque kilo-octet coûte ~6 ms en Slow 4G (§ S9). Les suites headless sont construites
 *  de la même façon : elles prouvent donc que les shaders dépouillés compilent encore. */
const stripGlsl = {
  name: "strip-glsl-comments",
  setup(b) {
    b.onLoad({ filter: /packages[\\/]viewer[\\/]src[\\/].*\.ts$/ }, async (args) => {
      const src = await readFile(args.path, "utf8");
      const contents = src.replace(/`#version 300 es[\s\S]*?`/g, (glsl) =>
        glsl.replace(/[ \t]*\/\/[^\n]*/g, "").replace(/\n\s*\n/g, "\n"));
      return { contents, loader: "ts" };
    });
  },
};

async function bundle(entry, minify = false) {
  const out = await build({
    entryPoints: [resolve(ROOT, entry)],
    bundle: true,
    format: "iife",
    target: "es2022",
    write: false,
    minify,
    plugins: [stripGlsl],
    logLevel: "warning",
  });
  return out.outputFiles[0].text;
}

const targets = [
  { entry: "packages/viewer/test/headless.ts", out: "dist/headless.html", title: "Geno — correction des imposteurs" },
  { entry: "packages/viewer/test/macro.ts", out: "dist/macro.html", title: "Geno — correction de la plongée" },
  { entry: "apps/spike/main.ts", out: "dist/spike.html", title: "Geno — spike de rendu" },
  // Minifiée : sous Slow 4G, chaque kilo-octet de la page coûte ~6 ms avant la première
  // requête de données. Les deux autres pages sont des outils, celle-ci est mesurée.
  { entry: "apps/first/main.ts", out: "dist/first.html", title: "Geno", minify: true },
  { entry: "apps/dive/main.ts", out: "dist/dive.html", title: "Geno — plongée", minify: true },
];

await mkdir(resolve(ROOT, "dist"), { recursive: true });
for (const t of targets) {
  const js = await bundle(t.entry, t.minify);
  await writeFile(resolve(ROOT, t.out), page(t.title, js));
  console.log(`  ${t.out}  ${(js.length / 1024).toFixed(1)} ko`);
}

// Banc matériel : une page au format Artifact (pas de squelette html/head/body, la
// publication l'ajoute), et le .g3d à côté. Publiée pour être ouverte sur de vrais appareils.
{
  const js = await bundle("apps/dive/materiel.ts", true);
  const tpl = await readFile(resolve(ROOT, "apps/dive/materiel.html"), "utf8");
  await mkdir(resolve(ROOT, "dist/materiel"), { recursive: true });
  await writeFile(resolve(ROOT, "dist/materiel/index.html"), tpl.replace("/*SCRIPT*/", () => js.replace(/<\/script/gi, "<\\/script")));
  const g3d = resolve(ROOT, "pipeline/data/g3d/gm12878.g3d");
  if (existsSync(g3d)) await copyFile(g3d, resolve(ROOT, "dist/materiel/gm12878.g3d"));
  console.log(`  dist/materiel/index.html  ${(js.length / 1024).toFixed(1)} ko${existsSync(g3d) ? "  (+ gm12878.g3d)" : "  (sans .g3d : make g3d)"}`);
}
