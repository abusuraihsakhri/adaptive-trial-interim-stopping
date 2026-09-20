"""Command-line interface for adaptive trial interim monitoring."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
import sys
from typing import Iterable, Sequence

from adaptive_boundaries import (
    FutilityAssessor,
    LanDeMetsSpending,
    OBFlemingBoundary,
    PocockBoundary,
    SampleSizeReestimator,
)


def _fractions(value: str | None) -> Sequence[float] | None:
    if value is None:
        return None
    items = [x.strip() for x in value.split(",") if x.strip()]
    if not items:
        raise argparse.ArgumentTypeError("fractions must contain at least one value")
    try:
        return [float(x) for x in items]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("fractions must be comma-separated numbers") from exc


def _boundary_rows(method: str, alpha: float, looks: int, fractions=None):
    if method == "obf":
        cls = OBFlemingBoundary
    elif method == "pocock":
        cls = PocockBoundary
    else:
        raise ValueError("method must be 'obf' or 'pocock'")
    return cls(alpha, looks, information_fractions=fractions).compute_boundaries()


def _print_boundary_table(rows: Iterable) -> None:
    rows = list(rows)
    print(f"{'Look':>4} {'Info':>8} {'Z boundary':>12} {'Alpha incr.':>12} {'Alpha cum.':>12}")
    for row in rows:
        print(
            f"{row.look_number:>4} {row.information_fraction:>8.3f} "
            f"{row.z_bound:>12.4f} {row.alpha_spent:>12.7f} {row.cumulative_alpha:>12.7f}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="adaptive-trial",
        description="Group-sequential efficacy boundaries and adaptive-trial monitoring utilities.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    def add_boundary_args(p):
        p.add_argument("--method", choices=("obf", "pocock"), default="obf")
        p.add_argument("--alpha", type=float, default=0.05)
        p.add_argument("--looks", type=int, default=5)
        p.add_argument("--fractions", type=_fractions, help="Comma-separated information fractions")
        p.add_argument("--json", action="store_true")

    add_boundary_args(sub.add_parser("boundaries", help="Calculate group-sequential efficacy boundaries"))
    add_boundary_args(sub.add_parser("audit", help="Compatibility alias for boundaries"))

    p_spend = sub.add_parser("spending", help="Calculate Lan-DeMets alpha spending")
    p_spend.add_argument("--method", choices=("obf", "pocock"), default="obf")
    p_spend.add_argument("--alpha", type=float, default=0.05)
    p_spend.add_argument("--fractions", type=_fractions, required=True)
    p_spend.add_argument("--json", action="store_true")

    p_fut = sub.add_parser("futility", help="Calculate conditional power from an interim Z statistic")
    p_fut.add_argument("--interim-z", type=float, required=True)
    p_fut.add_argument("--information-fraction", type=float, required=True)
    p_fut.add_argument("--critical-z", type=float, default=1.959963984540054)
    p_fut.add_argument("--assumed-final-mean-z", type=float)
    p_fut.add_argument("--threshold", type=float, default=0.10)
    p_fut.add_argument("--json", action="store_true")

    p_ssr = sub.add_parser("ssr", help="Simple sample-size re-estimation")
    p_ssr.add_argument("--mode", choices=("blinded", "unblinded"), required=True)
    p_ssr.add_argument("--n", type=int, required=True)
    p_ssr.add_argument("--current-variance", type=float)
    p_ssr.add_argument("--expected-variance", type=float)
    p_ssr.add_argument("--observed-effect", type=float)
    p_ssr.add_argument("--expected-effect", type=float)
    p_ssr.add_argument("--json", action="store_true")

    p_batch = sub.add_parser("batch", help="Calculate boundaries for rows in a CSV file")
    p_batch.add_argument("-i", "--input", required=True)
    p_batch.add_argument("-o", "--output", default="results.csv")

    p_serve = sub.add_parser("serve", help="Launch the optional FastAPI server")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8000)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.command in {"boundaries", "audit"}:
        rows = _boundary_rows(args.method, args.alpha, args.looks, args.fractions)
        if args.json:
            print(json.dumps([asdict(x) for x in rows], indent=2))
        else:
            print(rows[0].method)
            _print_boundary_table(rows)
        return 0

    if args.command == "spending":
        result = LanDeMetsSpending().compute_looks(args.fractions, args.alpha, args.method)
        payload = asdict(result)
        if args.json:
            print(json.dumps(payload, indent=2))
        else:
            print(result.function_name)
            for i, (inc, cum) in enumerate(zip(result.alpha_at_looks, result.cumulative_alpha), 1):
                print(f"Look {i}: incremental={inc:.8f} cumulative={cum:.8f}")
        return 0

    if args.command == "futility":
        result = FutilityAssessor.assess_from_z(
            args.interim_z,
            args.information_fraction,
            args.critical_z,
            args.assumed_final_mean_z,
            args.threshold,
        )
        payload = asdict(result)
        if args.json:
            print(json.dumps(payload, indent=2))
        else:
            print(f"Conditional power: {result.conditional_power:.2%}")
            print(f"Futility flag: {'yes' if result.is_futile else 'no'}")
            print(result.recommendation)
        return 0

    if args.command == "ssr":
        ssr = SampleSizeReestimator()
        if args.mode == "blinded":
            if args.current_variance is None or args.expected_variance is None:
                raise SystemExit("blinded mode requires --current-variance and --expected-variance")
            result = ssr.blinded_reestimate(args.n, args.current_variance, args.expected_variance)
        else:
            if args.observed_effect is None or args.expected_effect is None:
                raise SystemExit("unblinded mode requires --observed-effect and --expected-effect")
            result = ssr.unblinded_reestimate(args.n, args.observed_effect, args.expected_effect)
        if args.json:
            print(json.dumps(asdict(result), indent=2))
        else:
            print(f"Original N: {result.original_n}")
            print(f"Revised N: {result.revised_n}")
            print(f"Inflation factor: {result.inflation_factor:.4f}")
        return 0

    if args.command == "batch":
        output_rows = []
        with open(args.input, newline="", encoding="utf-8-sig") as handle:
            for record in csv.DictReader(handle):
                method = (record.get("method") or "obf").strip().lower()
                alpha = float(record.get("alpha") or 0.05)
                looks = int(record.get("looks") or 5)
                for row in _boundary_rows(method, alpha, looks):
                    output_rows.append({
                        "method": row.method,
                        "alpha": alpha,
                        "look": row.look_number,
                        "information_fraction": row.information_fraction,
                        "z_boundary": row.z_bound,
                        "alpha_increment": row.alpha_spent,
                        "cumulative_alpha": row.cumulative_alpha,
                    })
        with open(args.output, "w", newline="", encoding="utf-8") as handle:
            fields = ["method", "alpha", "look", "information_fraction", "z_boundary", "alpha_increment", "cumulative_alpha"]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(output_rows)
        print(f"Wrote {len(output_rows)} boundary rows to {args.output}")
        return 0

    if args.command == "serve":
        try:
            import uvicorn
        except ImportError:
            print("Server dependencies are not installed. Install with: pip install -e '.[server]'", file=sys.stderr)
            return 1
        from .server import create_app
        uvicorn.run(create_app(), host=args.host, port=args.port)
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
