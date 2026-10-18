/** « Sans chargement visible » : la plongée jouée en temps réel, réseau bridé.
 *
 *  La caméra suit l'horloge (`dive.html?stream=S`) : si une image arrive tard, la caméra
 *  a avancé quand même, comme pour un utilisateur. Le réseau est bridé par CDP au profil
 *  Slow 4G de la semaine 9. Chaque image dit si elle montre un chargement — hélice due sans
 *  sa séquence, ou nucléosomes dus pas encore générés — et ces images sont comptées.
 *
 *  Le témoin est la même plongée **sans préchargement** : la séquence n'est demandée que
 *  quand l'hélice en a besoin. Il doit compter des images en attente ; sinon la mesure ne
 *  voit rien et le zéro de la vraie plongée ne prouve rien.
 *
 *  SwiftShader rend 4 à 8 images par seconde : la caméra avance donc par grands pas, et le
 *  temps de calcul y est une borne haute. Le réseau, lui, est le vrai sujet de cette mesure.
 *
 *    node tools/stream.mjs [secondes …]      défaut : 9 (la vitesse du bouton) et 40
 */

import { cdp, launch, shutdown } from "./cdp.mjs";
import { serve } from "./serve.mjs";

const SECONDS = process.argv.slice(2).map(Number).filter((x) => x > 0);
const DURATIONS = SECONDS.length ? SECONDS : [9, 40];
const SLOW_4G = { latency: 562.5, downloadThroughput: (1.6e6 / 8) * 0.9, uploadThroughput: (750e3 / 8) * 0.9 };

async function run(browser, url, net) {
  const { targetId } = await browser.send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await browser.send("Target.attachToTarget", { targetId, flatten: true });
  const s = (m, p) => browser.send(m, p, sessionId);
  await s("Network.enable");
  await s("Network.setCacheDisabled", { cacheDisabled: true });
  await s("Network.emulateNetworkConditions", { offline: false, ...net });
  await s("Emulation.setDeviceMetricsOverride", { width: 960, height: 540, deviceScaleFactor: 1, mobile: false });
  await s("Page.navigate", { url });
  let r = null;
  const deadline = Date.now() + 20 * 60_000;
  while (!r && Date.now() < deadline) {
    await new Promise((x) => setTimeout(x, 250));
    const { result } = await s("Runtime.evaluate", { expression: "window.__geno ? JSON.stringify(window.__geno) : null", returnByValue: true });
    if (result.value) r = JSON.parse(result.value);
  }
  await browser.send("Target.closeTarget", { targetId });
  if (!r) throw new Error(`pas de résultat : ${url}`);
  return r;
}

const server = await serve(0);
const base = `http://127.0.0.1:${server.address().port}`;
const { proc, url, profile } = await launch();
const browser = await cdp(url);
const rows = [];
try {
  for (const sec of DURATIONS) {
    for (const [name, flag] of [["préchargement", ""], ["témoin : à la demande", "&noprefetch"]]) {
      const r = await run(browser, `${base}/dive.html?f=data/gm12878.g3d&stream=${sec}${flag}`, SLOW_4G);
      const f = r.frames;
      const firstHelix = f.find((x) => x.helix > 0);
      const waits = f.filter((x) => x.loading);
      rows.push({ name, sec, frames: f.length, loading: waits.length, stats: r.stats, bytes: r.sequence_bytes, requests: r.requests,
                  first: waits[0] ?? null, helixAt: firstHelix?.t ?? null });
    }
  }
} finally {
  browser.close();
  server.close();
  await shutdown(proc, profile);
}
console.log(`\n  Plongée jusqu'à TP53 c.818 en temps réel, réseau bridé Slow 4G (562,5 ms, 180 ko/s), SwiftShader\n`);
console.log(`  ${"durée".padStart(6)}  ${"".padEnd(24)} ${"images".padStart(7)} ${"chargement visible".padStart(19)} ${"dont séquence".padStart(14)} ${"dont nucléosomes".padStart(17)} ${"séquence lue".padStart(13)} ${"requêtes".padStart(9)}`);
for (const r of rows) {
  console.log(
    `  ${String(r.sec).padStart(4)} s  ${r.name.padEnd(24)} ${String(r.frames).padStart(7)} ${String(r.loading).padStart(19)} ` +
      `${String(r.stats?.sequence ?? "?").padStart(14)} ${String(r.stats?.nucleosomes ?? "?").padStart(17)} ${(r.bytes / 1024).toFixed(1).padStart(10)} Kio ${String(r.requests).padStart(9)}` +
      (r.first ? `   1re attente à ${r.first.t.toFixed(1)} s, champ ${(r.first.w * 1000).toFixed(0)} nm` : ""),
  );
}
const ok = rows.filter((r) => r.name === "préchargement").every((r) => r.loading === 0) && rows.filter((r) => r.name !== "préchargement").some((r) => r.loading > 0);
console.log(`\n  Verdict : ${ok ? "\x1b[32maucun chargement visible\x1b[0m" : "\x1b[31mnon démontré\x1b[0m"} — le témoin doit en montrer, la plongée aucun.\n`);
process.stdout.write(`JSON ${JSON.stringify(rows)}\n`);
process.exit(ok ? 0 : 1);
