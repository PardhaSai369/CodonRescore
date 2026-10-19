#!/usr/bin/env python3
"""One-command launcher: runs the full analysis pipeline, builds the dashboard and opens it.

    python run.py                # run pipeline (~1 min), then serve dashboard at http://localhost:8000
    python run.py --skip-pipeline  # reuse existing results, just serve the dashboard
    python run.py --port 9000 --no-browser
"""
import argparse, functools, http.server, socketserver, subprocess, sys, threading, time, webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [("Benchmarking gene callers against RefSeq annotation", "code/evaluate.py"),
         ("Sensitivity analysis, bootstrap CIs, C. diphtheriae run, figures", "code/analysis.py"),
         ("Exporting dashboard data", "code/export_dashboard.py")]


def run_step(msg, script):
    print(f"[run] {msg} ...", flush=True)
    t = time.time()
    r = subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, capture_output=True, text=True)
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-3000:], sep="\n")
        sys.exit(f"[run] step failed: {script}")
    print(f"[run]   done in {time.time() - t:.0f}s", flush=True)


def build_dashboard():
    d = ROOT / "dashboard"
    html = (d / "template.html").read_text(encoding="utf-8").replace("/*DATA*/", (d / "data.json").read_text(encoding="utf-8"))
    (d / "index.html").write_text(html, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-pipeline", action="store_true", help="reuse existing results/")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    try:
        import numpy, pyrodigal, Bio, matplotlib  # noqa
    except ImportError as e:
        sys.exit(f"Missing dependency ({e.name}). Run: pip install -r requirements.txt")
    if not a.skip_pipeline or not (ROOT / "dashboard/data.json").exists():
        for msg, s in STEPS:
            run_step(msg, s)
    build_dashboard()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT / "dashboard"))
    handler.log_message = lambda *args, **kw: None
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("127.0.0.1", a.port), handler) as srv:
        url = f"http://localhost:{a.port}/index.html"
        print(f"[run] Dashboard running at {url}   (Ctrl+C to stop)", flush=True)
        if not a.no_browser:
            threading.Timer(0.6, lambda: webbrowser.open(url)).start()
        try:
            srv.serve_forever()
        except KeyboardInterrupt:
            print("\n[run] stopped")


if __name__ == "__main__":
    main()
