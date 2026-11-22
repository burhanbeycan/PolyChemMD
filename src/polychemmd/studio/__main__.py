"""Independent local application entry point."""

import argparse
import json
import os
from pathlib import Path
import sys


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="polychemmd-studio", description="Local MXene–polymer simulation studio"
    )
    p.add_argument(
        "--data-dir",
        type=Path,
        default=Path(
            os.environ.get(
                "POLYCHEMMD_STUDIO_DATA", Path.home() / ".polychemmd" / "studio"
            )
        ),
    )
    sub = p.add_subparsers(dest="command", required=True)
    server = sub.add_parser("serve")
    server.add_argument("--port", type=int, default=8765)
    server.add_argument(
        "--external-worker",
        action="store_true",
        help="Run a worker in a separate terminal.",
    )
    worker = sub.add_parser("worker")
    m = worker.add_mutually_exclusive_group()
    m.add_argument("--once", action="store_true")
    m.add_argument("--drain", action="store_true")
    sub.add_parser("doctor")
    run = sub.add_parser("run")
    from .jobs import KINDS

    run.add_argument("kind", choices=KINDS)
    run.add_argument("recipe", type=Path)
    sub.add_parser("jobs")
    plan = sub.add_parser("plan")
    plan.add_argument(
        "--objective",
        choices=("termination", "hydrolysis", "model"),
        default="termination",
    )
    plan.add_argument("--model", help="Optional installed Ollama tag on localhost")
    campaign = sub.add_parser(
        "campaign", help="Queue paired, independently seeded study runs"
    )
    campaign.add_argument("recipe", type=Path)
    inspect = sub.add_parser(
        "campaign-status", help="Inspect completed pairs and uncertainty"
    )
    inspect.add_argument("id")
    args = p.parse_args(argv)
    try:
        if args.command == "serve":
            import uvicorn
            from .server import app_factory

            uvicorn.run(
                app_factory(args.data_dir, start_worker=not args.external_worker),
                host="127.0.0.1",
                port=args.port,
            )
            return 0
        if args.command == "doctor":
            from .operations import status

            value = status()
        elif args.command == "plan":
            from .planner import suggest

            value = suggest(
                "ollama" if args.model else "rules", args.model, args.objective
            )
        else:
            from .jobs import JobStore, execute_next, work

            store = JobStore(args.data_dir)
            if args.command == "campaign":
                from .campaign import submit_campaign

                value = submit_campaign(
                    store, json.loads(args.recipe.read_text(encoding="utf-8"))
                )
            elif args.command == "campaign-status":
                from .campaign import campaign_summary

                value = campaign_summary(store, args.id)
            elif args.command == "jobs":
                value = store.list()
            elif args.command == "worker":
                if args.once:
                    execute_next(store)
                elif args.drain:
                    while execute_next(store):
                        pass
                else:
                    work(store)
                return 0
            else:
                item = store.submit(
                    args.kind, json.loads(args.recipe.read_text(encoding="utf-8"))
                )
                execute_next(store, item["id"])
                value = store.get(item["id"])
                # A separate process can claim a job between submission and this call.
                import time

                while value["status"] in ("queued", "running"):
                    time.sleep(0.25)
                    execute_next(store, item["id"])
                    value = store.get(item["id"])
                if value["status"] == "failed":
                    print(
                        json.dumps(value, ensure_ascii=False, indent=2), file=sys.stderr
                    )
                    return 1
        print(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False))
        return 0
    except (ValueError, KeyError, ImportError, RuntimeError, OSError) as e:
        print(f"polychemmd-studio: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
