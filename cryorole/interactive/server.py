"""Loopback-only standard-library HTTP server for interactive exploration."""

from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import threading
from typing import Any

from cryorole.interactive.session import ExploreSession


_ASSET_ROOT = Path(__file__).with_name("assets")
_STATIC = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/assets/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/assets/styles.css": ("styles.css", "text/css; charset=utf-8"),
}


class ExploreServer:
    def __init__(self, session: ExploreSession, *, port: int = 0) -> None:
        self.session = session
        self.host = "127.0.0.1"
        handler = _handler_for(session)
        self._server = ThreadingHTTPServer((self.host, int(port)), handler)
        self.port = int(self._server.server_address[1])
        self.url = f"http://{self.host}:{self.port}/"
        self._thread: threading.Thread | None = None

    def serve_forever(self) -> None:
        self._server.serve_forever(poll_interval=0.1)

    def start_in_thread(self) -> threading.Thread:
        if self._thread and self._thread.is_alive():
            return self._thread
        self._thread = threading.Thread(target=self.serve_forever, name="cryorole-explore", daemon=True)
        self._thread.start()
        return self._thread

    def shutdown(self) -> None:
        self._server.shutdown()
        self._server.server_close()


def create_explore_server(session: ExploreSession, *, port: int = 0) -> ExploreServer:
    return ExploreServer(session, port=port)


def _handler_for(session: ExploreSession):
    class Handler(BaseHTTPRequestHandler):
        server_version = "cryoROLEExplore/1.0"

        def do_GET(self) -> None:  # noqa: N802
            path = self.path.split("?", 1)[0]
            if path == "/api/health":
                self._json(200, {"status": "ok", "run_id": session.run_id})
                return
            if path == "/api/session":
                self._json(200, session.session_payload())
                return
            asset = _STATIC.get(path)
            if asset is None:
                self._json(404, {"error": "not_found"})
                return
            filename, content_type = asset
            data = (_ASSET_ROOT / filename).read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header(
                "Content-Security-Policy",
                "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'",
            )
            self.end_headers()
            self.wfile.write(data)

        def do_POST(self) -> None:  # noqa: N802
            path = self.path.split("?", 1)[0]
            if path not in {"/api/evaluate", "/api/confirm", "/api/shutdown"}:
                self._json(404, {"error": "not_found"})
                return
            try:
                payload = self._payload()
                if payload.get("session_token") != session.session_token:
                    raise ValueError("invalid session token")
                if path == "/api/evaluate":
                    result = session.evaluate(
                        center=payload["center"],
                        representation=payload.get("representation", "rotvec"),
                        radius_deg=float(payload["radius_deg"]),
                        expected_run_id=payload.get("run_id"),
                        expected_landscape_sha256=payload.get("landscape_sha256"),
                    )
                    self._json(200, result)
                    return
                if path == "/api/confirm":
                    result = session.confirm(
                        str(payload.get("selection_id", "")),
                        expected_run_id=payload.get("run_id"),
                        expected_landscape_sha256=payload.get("landscape_sha256"),
                    )
                    self._json(201, result)
                    return
                self._json(200, {"status": "shutting_down"})
                threading.Thread(target=self.server.shutdown, daemon=True).start()
            except (KeyError, TypeError, ValueError, FileExistsError) as exc:
                self._json(409 if isinstance(exc, FileExistsError) else 400, {"error": str(exc)})

        def _payload(self) -> dict[str, Any]:
            length = int(self.headers.get("Content-Length", "0"))
            if length < 1 or length > 1_000_000:
                raise ValueError("request body must be between 1 byte and 1 MB")
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("JSON request must be an object")
            return payload

        def _json(self, status: int, payload: dict[str, Any]) -> None:
            data = (json.dumps(payload, sort_keys=True) + "\n").encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def log_message(self, _format: str, *_args) -> None:
            return

    return Handler
