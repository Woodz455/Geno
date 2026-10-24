/** Bundle esbuild → un seul HTML autonome.
 *
 *  Autonome parce que les modules ES ne se chargent pas depuis `file://` (CORS),
 *  et qu'un serveur de dev juste pour lancer un test headless est une pièce
 *  mobile de plus qui peut tomber en panne.
 */

import { build } from "esbuild";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { existsSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");

const page = (title, script, body = "", css = "") => `<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><title>${title}</title>${css ? `\n<style>${css}</style>` : ""}
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
    // React (sous HiGlass, semaine 14) lit `process.env.NODE_ENV` : sans valeur, la page
    // tomberait sur le build de développement, ou sur une ReferenceError.
    define: { "process.env.NODE_ENV": '"production"', global: "globalThis" },
    logLevel: "warning",
  });
  return out.outputFiles[0].text.replace(/<\/script/gi, "<\\/script");
}

const targets = [
  { entry: "packages/viewer/test/headless.ts", out: "dist/headless.html", title: "Geno — correction des imposteurs" },
  { entry: "packages/viewer/test/macro.ts", out: "dist/macro.html", title: "Geno — correction de la plongée" },
  { entry: "packages/viewer/test/fiber.ts", out: "dist/fiber.html", title: "Geno — correction de la fibre" },
  { entry: "packages/viewer/test/helix.ts", out: "dist/helix.html", title: "Geno — correction de l'hélice" },
  { entry: "packages/viewer/test/pick.ts", out: "dist/pick.html", title: "Geno — correction du picking" },
  { entry: "apps/spike/main.ts", out: "dist/spike.html", title: "Geno — spike de rendu" },
  // Minifiée : sous Slow 4G, chaque kilo-octet de la page coûte ~6 ms avant la première
  // requête de données. Les deux autres pages sont des outils, celle-ci est mesurée.
  { entry: "apps/first/main.ts", out: "dist/first.html", title: "Geno", minify: true },
  { entry: "apps/dive/main.ts", out: "dist/dive.html", title: "Geno — plongée", minify: true },
  // Semaine 15 : le mode variant structural (extrusion 1D en direct, allèles OpenMM).
  { entry: "apps/variant/main.ts", out: "dist/variant.html", title: "Geno — variant structural", minify: true },
  // Semaine 14 : igv.js et HiGlass embarqués (React, pixi) — 4 Mo, pas une page de premier rendu.
  { entry: "apps/linked/main.ts", out: "dist/vues.html", title: "Geno — vues liées", minify: true,
    css: "apps/linked/node_modules/higlass/dist/hglib.css" },
];

await mkdir(resolve(ROOT, "dist"), { recursive: true });
for (const t of targets) {
  const js = await bundle(t.entry, t.minify);
  const css = t.css ? await readFile(resolve(ROOT, t.css), "utf8") : "";
  await writeFile(resolve(ROOT, t.out), page(t.title, js, "", css));
  console.log(`  ${t.out}  ${(js.length / 1024).toFixed(1)} ko`);
}

// Banc matériel : une page au format Artifact (pas de squelette html/head/body, la
// publication l'ajoute), avec le .g3d **dans** la page en base64 — un Artifact ne sert
// pas de fichier binaire à côté. Publiée pour être ouverte sur de vrais appareils.
{
  const js = await bundle("apps/dive/materiel.ts", true);
  const tpl = await readFile(resolve(ROOT, "apps/dive/materiel.html"), "utf8");
  const g3d = resolve(ROOT, "pipeline/data/g3d/gm12878.g3d");
  const data = existsSync(g3d) ? (await readFile(g3d)).toString("base64") : "";
  await mkdir(resolve(ROOT, "dist/materiel"), { recursive: true });
  const html = tpl.replace("<script>/*SCRIPT*/</script>", () =>
    `<script id="g3d" type="text/plain">${data}</script>\n<script>${js.replace(/<\/script/gi, "<\\/script")}</script>`);
  await writeFile(resolve(ROOT, "dist/materiel/index.html"), html);
  console.log(`  dist/materiel/index.html  ${(html.length / 1024).toFixed(1)} ko${data ? "  (.g3d intégré)" : "  (sans .g3d : make g3d)"}`);
}
