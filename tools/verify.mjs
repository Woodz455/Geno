/** Exécute la suite de correction dans Chromium headless et rend son verdict.
 *
 *  Rappel permanent : headless rend via SwiftShader, un rasteriseur logiciel.
 *  Ce script prouve que la technique est CORRECTE. Il ne dit rien de sa
 *  performance — pour ça, ouvrir `dist/spike.html` sur du vrai matériel.
 */

import { execFile } from "node:child_process";
import { existsSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const exec = promisify(execFile);
const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");

const CANDIDATES = [
  process.env.CHROME_PATH,
  "/opt/pw-browsers/chromium-1194/chrome-linux/chrome",
  "/usr/bin/chromium",
  "/usr/bin/google-chrome",
].filter(Boolean);

const chrome = CANDIDATES.find((p) => existsSync(p));
if (!chrome) {
  console.error("Chromium introuvable. Définir CHROME_PATH.");
  process.exit(2);
}

// Une suite par page : celle de la semaine 4 (imposteurs, picking) et celle de la
// semaine 10 (caméra, fondu tramé, occlusion), celle de la semaine 11 (tubes,
// nucléosomes) et celle de la semaine 12 (chemin de l'ADN, double hélice). Chacune se
// compte à part.
const SUITES = [
  { page: "dist/headless.html", title: "Correction des imposteurs (semaine 4)" },
  { page: "dist/macro.html", title: "Plongée : caméra, fondu tramé, occlusion (semaine 10)" },
  { page: "dist/fiber.html", title: "Fibre : tubes, nucléosomes (semaine 11)" },
  { page: "dist/helix.html", title: "Hélice : chemin de l'ADN, B-ADN, séquence (semaine 12)" },
];

const decode = (s) =>
  s.replace(/&quot;/g, '"').replace(/&#39;/g, "'").replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&amp;/g, "&");

let failed = 0;
for (const suite of SUITES) {
  const pagePath = resolve(ROOT, suite.page);
  if (!existsSync(pagePath)) {
    console.error(`${suite.page} absent — lancer \`pnpm build\` d'abord.`);
    process.exit(2);
  }
  const { stdout } = await exec(
    chrome,
    [
      "--headless",
      "--no-sandbox",
      "--use-gl=angle",
      "--use-angle=swiftshader",
      "--virtual-time-budget=20000",
      "--dump-dom",
      `file://${pagePath}`,
    ],
    { maxBuffer: 64 * 1024 * 1024 },
  );
  const match = stdout.match(/<pre id="results">([\s\S]*?)<\/pre>/);
  if (!match) {
    console.error(`${suite.page} n'a pas rendu de résultats — le script a probablement levé.`);
    process.exit(1);
  }
  const payload = JSON.parse(decode(match[1]));
  console.log(`\n  ${suite.title} — Chromium headless (SwiftShader, rendu logiciel)\n`);
  for (const r of payload.results) {
    const mark = r.ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m";
    console.log(`  ${mark}  ${r.name}${r.detail ? `  \x1b[2m${r.detail}\x1b[0m` : ""}`);
  }
  console.log(`\n  ${payload.passed} réussis, ${payload.failed} échoués`);
  failed += payload.failed;
}
console.log(
  "\n  \x1b[2mCe verdict porte sur la correction, pas sur la performance :\n" +
    "  SwiftShader est un rasteriseur logiciel. Pour les fps, ouvrir dist/spike.html\n" +
    "  sur du matériel réel.\x1b[0m\n",
);

process.exit(failed === 0 ? 0 : 1);
