/** Serveur statique minimal, qui parle Range.
 *
 *  `/` sert `dist/`, `/data/` sert `pipeline/data/g3d/`. Les `.g3d` partent bruts avec
 *  `Accept-Ranges` et répondent 206 à une plage — jamais compressés à la volée, sans
 *  quoi la plage désignerait des octets d'un flux gzip et ne servirait à rien. Le HTML
 *  part en gzip, comme sur n'importe quel hébergement statique : mesurer une page non
 *  compressée gonflerait le temps du document et ferait accuser le format.
 *
 *    node tools/serve.mjs [port]
 */

import { createServer } from "node:http";
import { readFile, stat } from "node:fs/promises";
import { extname, join, normalize, resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { gzipSync } from "node:zlib";

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const MOUNTS = [
  ["/data/", resolve(ROOT, "pipeline/data/g3d")],
  ["/", resolve(ROOT, "dist")],
];
const TYPES = { ".html": "text/html; charset=utf-8", ".g3d": "application/octet-stream", ".js": "text/javascript" };

function locate(url) {
  const path = decodeURIComponent(new URL(url, "http://x").pathname);
  for (const [prefix, dir] of MOUNTS) {
    if (path.startsWith(prefix)) {
      const full = normalize(join(dir, path.slice(prefix.length)));
      return full.startsWith(dir) ? full : null;       // pas de sortie du dossier monté
    }
  }
  return null;
}

export function serve(port = 0) {
  const server = createServer(async (req, res) => {
    const file = locate(req.url);
    let info;
    try {
      info = file && (await stat(file));
    } catch {
      info = null;
    }
    if (!info || !info.isFile()) {
      res.writeHead(404).end();
      return;
    }
    const type = TYPES[extname(file)] ?? "application/octet-stream";
    const body = await readFile(file);
    const head = { "Content-Type": type, "Cache-Control": "no-store", "Accept-Ranges": "bytes" };

    const range = req.headers.range && /^bytes=(\d+)-(\d*)$/.exec(req.headers.range);
    if (range) {
      const start = Number(range[1]);
      const end = Math.min(range[2] ? Number(range[2]) : body.length - 1, body.length - 1);
      if (start >= body.length || start > end) {
        res.writeHead(416, { ...head, "Content-Range": `bytes */${body.length}` }).end();
        return;
      }
      res.writeHead(206, {
        ...head,
        "Content-Range": `bytes ${start}-${end}/${body.length}`,
        "Content-Length": end - start + 1,
      });
      res.end(body.subarray(start, end + 1));
      return;
    }
    if (type.startsWith("text/") && /gzip/.test(req.headers["accept-encoding"] ?? "")) {
      const gz = gzipSync(body, { level: 9 });
      res.writeHead(200, { ...head, "Content-Encoding": "gzip", "Content-Length": gz.length });
      res.end(gz);
      return;
    }
    res.writeHead(200, { ...head, "Content-Length": body.length });
    res.end(body);
  });
  return new Promise((ok) => server.listen(port, "127.0.0.1", () => ok(server)));
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const server = await serve(Number(process.argv[2] ?? 8080));
  console.log(`http://127.0.0.1:${server.address().port}/first.html`);
}
