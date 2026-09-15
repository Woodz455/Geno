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

const pagePath = resolve(ROOT, "dist/headless.html");
if (!existsSync(pagePath)) {
  console.error("dist/headless.html absent — lancer `pnpm build` d'abord.");
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
  console.error("la page n'a pas rendu de résultats — le script a probablement levé.");
  process.exit(1);
}

const decode = (s) =>
  s.replace(/&quot;/g, '"').replace(/&#39;/g, "'").replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&amp;/g, "&");

const payload = JSON.parse(decode(match[1]));

console.log("\n  Correction des imposteurs — Chromium headless (SwiftShader, rendu logiciel)\n");
for (const r of payload.results) {
  const mark = r.ok ? "\x1b[32mPASS\x1b[0m" : "\x1b[31mFAIL\x1b[0m";
  console.log(`  ${mark}  ${r.name}${r.detail ? `  \x1b[2m${r.detail}\x1b[0m` : ""}`);
}
console.log(`\n  ${payload.passed} réussis, ${payload.failed} échoués`);
console.log(
  "\n  \x1b[2mCe verdict porte sur la correction, pas sur la performance :\n" +
    "  SwiftShader est un rasteriseur logiciel. Pour les fps, ouvrir dist/spike.html\n" +
    "  sur du matériel réel.\x1b[0m\n",
);

process.exit(payload.failed === 0 ? 0 : 1);
