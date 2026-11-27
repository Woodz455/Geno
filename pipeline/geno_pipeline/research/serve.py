"""Le service local du mode Recherche : `geno research serve`.

Il n'écoute que la boucle locale (127.0.0.1) et refuse toute requête d'écriture dont l'origine
n'est pas lui-même : une page ouverte ailleurs dans le navigateur ne peut pas lui faire lire un
fichier de la machine. Il sert la page (`dist/`), reçoit le fichier de l'utilisateur ou lit un
chemin local, lance les reconstructions une à une dans un fil, et sert le `.g3d` par plages.

    POST /api/upload?name=…        le corps est le fichier ; rend {file, sha256, inspect}
    POST /api/open     {path}      un fichier déjà sur la machine ; rend la même chose
    POST /api/run      {file, params}             → {job}
    POST /api/reproduce {file, manifest}          → {job}
    GET  /api/job/<id>             {state, step, progress, manifest | error | reproduce}
    GET  /sessions/<id>/structure.g3d | .pdb | .3dg | .gltf | session.json | statistics.json | ensemble.f32
"""

from __future__ import annotations

import json
import re
import shutil
import threading
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import session as S

_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,120}$")


class Service:
    def __init__(self, work: Path, static: Path | None):
        self.work = Path(work)
        self.static = Path(static) if static else None
        (self.work / "files").mkdir(parents=True, exist_ok=True)
        (self.work / "sessions").mkdir(parents=True, exist_ok=True)
        self.files: dict[str, dict] = {}        # jeton → {path, sha256, inspect}
        self.jobs: dict[str, dict] = {}
        self.lock = threading.Lock()            # une reconstruction à la fois : elles sont denses

    def register(self, path: Path) -> dict:
        cool = S.as_cooler(path, self.work / "converted")
        info = S.inspect(cool)
        info["name"], info["bytes"] = Path(path).name, Path(path).stat().st_size
        info["converted_from_hic"] = cool != Path(path)
        token = uuid.uuid4().hex[:16]
        self.files[token] = {"path": Path(path), "cool": cool, "sha256": S.sha256_file(path), "inspect": info}
        return {"file": token, "sha256": self.files[token]["sha256"], "inspect": info}

    def start(self, kind: str, token: str, payload: dict) -> str:
        if token not in self.files:
            raise S.SessionError("fichier inconnu de ce service : l'importer d'abord")
        job = uuid.uuid4().hex[:16]
        self.jobs[job] = {"state": "en attente", "step": "", "progress": 0.0, "kind": kind}
        threading.Thread(target=self._work, args=(job, kind, token, payload), daemon=True).start()
        return job

    def _work(self, job: str, kind: str, token: str, payload: dict) -> None:
        st = self.jobs[job]
        f = self.files[token]
        out = self.work / "sessions" / job

        def progress(step: str, frac: float) -> None:
            st.update(step=step, progress=frac)

        with self.lock:
            st["state"] = "en cours"
            try:
                if kind == "run":
                    try:
                        p = S.Params(**payload)
                    except TypeError as e:
                        raise S.SessionError(f"paramètres illisibles : {e}") from None
                    st["manifest"] = S.run(f["cool"], p, out, progress, input_sha256=f["sha256"], source=f["path"])
                else:
                    r = S.reproduce(payload, f["path"], out, progress, cache=self.work / "converted")
                    st["manifest"] = r.pop("manifest")
                    st["reproduce"] = r
                st["session"] = job
                st["state"] = "fait"
            except S.SessionError as e:
                st.update(state="refusé", error=str(e))
            except Exception as e:                           # une erreur du service, dite telle quelle
                st.update(state="échec", error=f"{type(e).__name__} : {e}")


def handler(svc: Service):
    class H(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *a):                           # silencieux : la page montre l'état
            pass

        def _json(self, code: int, obj) -> None:
            body = json.dumps(obj, ensure_ascii=False).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _same_origin(self) -> bool:
            host = self.headers.get("Host", "")
            origin = self.headers.get("Origin")
            if not re.fullmatch(r"(127\.0\.0\.1|localhost)(:\d+)?", host):
                return False                                 # contre le DNS rebinding
            return origin is None or origin in (f"http://{host}",)

        def _file(self, path: Path, ctype: str) -> None:
            if not path.is_file():
                self._json(404, {"error": "introuvable"})
                return
            size = path.stat().st_size
            rng = re.fullmatch(r"bytes=(\d+)-(\d*)", self.headers.get("Range", ""))
            start, end = 0, size - 1
            if rng:
                start = int(rng.group(1))
                end = min(int(rng.group(2)) if rng.group(2) else size - 1, size - 1)
                if start >= size or start > end:
                    self.send_response(416)
                    self.send_header("Content-Range", f"bytes */{size}")
                    self.send_header("Content-Length", "0")
                    self.end_headers()
                    return
            self.send_response(206 if rng else 200)
            self.send_header("Content-Type", ctype)
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("Cache-Control", "no-store")
            if rng:
                self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            self.send_header("Content-Length", str(end - start + 1))
            self.end_headers()
            with open(path, "rb") as f:
                f.seek(start)
                left = end - start + 1
                while left > 0:
                    b = f.read(min(1 << 20, left))
                    if not b:
                        break
                    self.wfile.write(b)
                    left -= len(b)

        def do_GET(self):
            if not self._same_origin():
                self._json(403, {"error": "origine refusée"})
                return
            u = urlparse(self.path)
            m = re.fullmatch(r"/api/job/([0-9a-f]{16})", u.path)
            if m:
                st = svc.jobs.get(m.group(1))
                self._json(200 if st else 404, st or {"error": "tâche inconnue"})
                return
            m = re.fullmatch(r"/sessions/([0-9a-f]{16})/(structure\.(?:g3d|pdb|3dg|gltf)|session\.json|statistics\.json|ensemble\.f32)", u.path)
            if m:
                self._file(svc.work / "sessions" / m.group(1) / m.group(2),
                           "application/json" if m.group(2).endswith(".json") else "application/octet-stream")
                return
            if svc.static and re.fullmatch(r"/[A-Za-z0-9._-]*", u.path):
                name = u.path[1:] or "recherche.html"
                self._file(svc.static / name, "text/html; charset=utf-8" if name.endswith(".html") else "application/octet-stream")
                return
            self._json(404, {"error": "introuvable"})

        def do_POST(self):
            if not self._same_origin():
                self._json(403, {"error": "origine refusée"})
                return
            u = urlparse(self.path)
            n = int(self.headers.get("Content-Length", "0"))
            try:
                if u.path == "/api/upload":
                    name = parse_qs(u.query).get("name", ["import.mcool"])[0]
                    if not _NAME.match(name):
                        raise S.SessionError("nom de fichier refusé")
                    d = svc.work / "files" / uuid.uuid4().hex[:16]
                    d.mkdir(parents=True)
                    dest = d / name
                    with open(dest, "wb") as f:
                        left = n
                        while left > 0:
                            b = self.rfile.read(min(1 << 20, left))
                            if not b:
                                break
                            f.write(b)
                            left -= len(b)
                    try:
                        self._json(200, svc.register(dest))
                    except S.SessionError:
                        shutil.rmtree(d, ignore_errors=True)
                        raise
                    return
                body = json.loads(self.rfile.read(n) or b"{}")
                if u.path == "/api/open":
                    p = Path(str(body.get("path", ""))).expanduser()
                    if not p.is_file():
                        raise S.SessionError("ce chemin n'est pas un fichier de cette machine")
                    self._json(200, svc.register(p))
                elif u.path == "/api/run":
                    self._json(200, {"job": svc.start("run", body.get("file", ""), body.get("params", {}))})
                elif u.path == "/api/reproduce":
                    self._json(200, {"job": svc.start("reproduce", body.get("file", ""), body.get("manifest", {}))})
                else:
                    self._json(404, {"error": "introuvable"})
            except S.SessionError as e:
                self._json(400, {"error": str(e)})
            except json.JSONDecodeError:
                self._json(400, {"error": "corps JSON illisible"})

    return H


def serve(work: Path, static: Path | None, port: int = 0) -> ThreadingHTTPServer:
    """Le serveur, prêt ; l'appelant choisit `serve_forever()` ou un fil."""
    srv = ThreadingHTTPServer(("127.0.0.1", port), handler(Service(work, static)))
    srv.daemon_threads = True
    return srv
