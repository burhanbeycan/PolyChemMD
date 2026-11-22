"""Install a local Studio environment and launch the loopback web application."""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
from urllib.request import urlopen
import venv
import webbrowser


REPO = Path(__file__).resolve().parents[1]


def setup(path, requirements):
    path = path.expanduser().resolve()
    bindir = path / ("Scripts" if os.name == "nt" else "bin")
    python = bindir / ("python.exe" if os.name == "nt" else "python")
    if not python.exists():
        venv.EnvBuilder(with_pip=True).create(path)
    environment = {
        **os.environ,
        "PATH": str(bindir) + os.pathsep + os.environ.get("PATH", ""),
    }
    subprocess.run(
        [str(python), "-m", "pip", "install", *requirements],
        env=environment,
        check=True,
    )
    return python, environment


def open_when_ready(port):
    url = f"http://127.0.0.1:{port}/"
    for _ in range(60):
        try:
            with urlopen(url, timeout=0.5) as r:
                if r.status == 200:
                    webbrowser.open(url)
                    return
        except OSError:
            time.sleep(0.25)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--full",
        action="store_true",
        help="Install both packers, LAMMPS, MDAnalysis and isolated RadonPy; prepare public IFF + OPLS.",
    )
    p.add_argument("--env-dir", type=Path, default=REPO / ".studio-env")
    p.add_argument(
        "--atlas",
        type=Path,
        help="Import the existing PVA/Atlas archive for optional MATCH/CGenFF parameters (full setup).",
    )
    p.add_argument("--port", type=int, default=8765)
    p.add_argument(
        "--open",
        action="store_true",
        help="Open the browser when the application is ready.",
    )
    p.add_argument("--setup-only", action="store_true")
    args = p.parse_args()
    if sys.version_info[:2] < (3, 11) or (
        args.full and sys.version_info[:2] != (3, 12)
    ):
        p.error("Use Python 3.12 for the tested full stack, or 3.11+ for the core GUI.")
    if args.atlas and (not args.full or not args.atlas.is_file()):
        p.error("--atlas needs --full and a local Atlas ZIP file")
    if not 1 <= args.port <= 65535:
        p.error("Invalid port")
    extras = (
        "studio,studio-packers,md-smoke,studio-analysis,studio-models"
        if args.full
        else "studio"
    )
    python, env = setup(args.env_dir, ["-e", str(REPO) + f"[{extras}]"])
    if args.full:
        rpython, _ = setup(
            args.env_dir.with_name(args.env_dir.name + "-radon"),
            ["radonpy-pypi==0.2.11", "rdkit==2026.9.1"],
        )
        env["POLYCHEMMD_RADON_PYTHON"] = str(rpython)
        prepare = [
            str(python),
            "-m",
            "polychemmd.studio.iff_sources",
            "--fetch",
            "--opls",
        ]
        if args.atlas:
            prepare += ["--atlas", str(args.atlas.resolve())]
        subprocess.run(prepare, env=env, check=True)
    else:
        other = (
            args.env_dir.with_name(args.env_dir.name + "-radon")
            / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
        ).resolve()
        if other.is_file():
            env.setdefault("POLYCHEMMD_RADON_PYTHON", str(other))
    if args.setup_only:
        return
    if args.open:
        threading.Thread(target=open_when_ready, args=(args.port,), daemon=True).start()
    raise SystemExit(
        subprocess.call(
            [str(python), "-m", "polychemmd.studio", "serve", "--port", str(args.port)],
            env=env,
        )
    )


if __name__ == "__main__":
    main()
