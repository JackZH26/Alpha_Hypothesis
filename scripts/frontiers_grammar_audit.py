#!/usr/bin/env python3
"""Finite grammar audit for the Bost-Connes alpha ansatz.

The audit is intentionally conservative and reproducible. It enumerates a
small expression space around the manuscript's three-term grammar and also
reports a fixed-seed Monte Carlo sample from the same space. The resulting
fractions are diagnostics, not formal p-values.
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

from frontiers_quadratic_field_screen import (
    L_minus_3_for_quadratic_discriminant,
    fundamental_discriminant,
    fundamental_unit_data,
    is_squarefree,
    trace_epsilon_squared,
)


getcontext().prec = 50

TARGET = Decimal("137.035999177")
HYPOTHESIS_VALUE = Decimal("137.03599916476563934505723564140907572836137437744")
HYPOTHESIS_ABS_RESIDUAL = abs(HYPOTHESIS_VALUE - TARGET)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit a finite arithmetic expression grammar.")
    parser.add_argument("--max-d", type=int, default=120, help="largest squarefree d to include")
    parser.add_argument("--max-exponent", type=int, default=8, help="largest regulator exponent")
    parser.add_argument("--max-coeff", type=int, default=3, help="largest absolute torsion coefficient")
    parser.add_argument("--top-n", type=int, default=20, help="number of closest expressions to record")
    parser.add_argument("--unit-search-limit", type=int, default=500000)
    parser.add_argument("--monte-carlo-samples", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260429)
    parser.add_argument("--summary-out", type=Path)
    parser.add_argument("--top-out", type=Path)
    return parser.parse_args()


def decimal_fraction(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def field_records(max_d: int, unit_search_limit: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for d in [n for n in range(2, max_d + 1) if is_squarefree(n)]:
        D = fundamental_discriminant(d)
        try:
            eps, a, b, denom = fundamental_unit_data(d, limit=unit_search_limit)
        except RuntimeError as exc:
            print(f"skipping {exc}", file=sys.stderr)
            continue
        Lm3 = L_minus_3_for_quadratic_discriminant(D)
        if Lm3 == 0:
            continue
        scale = decimal_fraction(Fraction(720, 1) / Lm3)
        trace_eps2 = trace_epsilon_squared(d, a, b, denom)
        trace_dec = decimal_fraction(trace_eps2)
        records.append({"d": d, "D": D, "eps": eps, "scale": scale, "trace": trace_dec})
    return records


def evaluate(record: dict[str, object], p: int, coeff: int, q: int, sign: int, r: int) -> Decimal:
    eps = record["eps"]
    scale = record["scale"]
    trace = record["trace"]
    assert isinstance(eps, Decimal)
    assert isinstance(scale, Decimal)
    assert isinstance(trace, Decimal)
    return (
        scale / (eps**p)
        + Decimal(coeff) * Decimal(2) / (eps**q)
        + Decimal(sign) / ((trace * eps) ** r)
    )


def expression_label(record: dict[str, object], p: int, coeff: int, q: int, sign: int, r: int) -> str:
    return (
        f"(3!/zeta_K(-3))*eps^-{p} "
        f"{coeff:+d}*w_K*eps^-{q} "
        f"{sign:+d}*(Tr(eps^2)*eps)^-{r}"
        f"; d={record['d']}, Delta={record['D']}"
    )


def maybe_add_top(top_rows: list[dict[str, object]], row: dict[str, object], top_n: int) -> None:
    top_rows.append(row)
    top_rows.sort(key=lambda item: item["abs_residual_decimal"])
    del top_rows[top_n:]


def write_summary(summary: dict[str, object], path: Path | None) -> None:
    fieldnames = list(summary)
    if path is None:
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(summary)
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(summary)


def write_top(top_rows: list[dict[str, object]], path: Path | None) -> None:
    if path is None:
        return
    fieldnames = [
        "rank",
        "d",
        "Delta",
        "p",
        "torsion_coeff",
        "q",
        "local_sign",
        "r",
        "value",
        "signed_residual",
        "abs_residual",
        "expression",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for rank, row in enumerate(top_rows, start=1):
            writer.writerow(
                {
                    "rank": rank,
                    "d": row["d"],
                    "Delta": row["Delta"],
                    "p": row["p"],
                    "torsion_coeff": row["torsion_coeff"],
                    "q": row["q"],
                    "local_sign": row["local_sign"],
                    "r": row["r"],
                    "value": row["value"],
                    "signed_residual": row["signed_residual"],
                    "abs_residual": row["abs_residual"],
                    "expression": row["expression"],
                }
            )


def sci_decimal(value: Decimal) -> str:
    if value == 0:
        return "0.000000000000000000E+0"
    return f"{value:.18E}"


def main() -> None:
    args = parse_args()
    records = field_records(args.max_d, args.unit_search_limit)
    exponents = range(1, args.max_exponent + 1)
    coeffs = range(-args.max_coeff, args.max_coeff + 1)
    signs = (-1, 0, 1)

    total = 0
    better_or_equal = 0
    strictly_better = 0
    top_rows: list[dict[str, object]] = []

    for record in records:
        for p in exponents:
            for coeff in coeffs:
                for q in exponents:
                    for sign in signs:
                        for r in exponents:
                            value = evaluate(record, p, coeff, q, sign, r)
                            signed_residual = value - TARGET
                            abs_residual = abs(signed_residual)
                            total += 1
                            if abs_residual < HYPOTHESIS_ABS_RESIDUAL:
                                strictly_better += 1
                            if abs_residual <= HYPOTHESIS_ABS_RESIDUAL:
                                better_or_equal += 1
                            if len(top_rows) < args.top_n or abs_residual < top_rows[-1]["abs_residual_decimal"]:
                                maybe_add_top(
                                    top_rows,
                                    {
                                        "abs_residual_decimal": abs_residual,
                                        "d": record["d"],
                                        "Delta": record["D"],
                                        "p": p,
                                        "torsion_coeff": coeff,
                                        "q": q,
                                        "local_sign": sign,
                                        "r": r,
                                        "value": f"{value:.18f}",
                                        "signed_residual": f"{signed_residual:.18E}",
                                        "abs_residual": f"{abs_residual:.18E}",
                                        "expression": expression_label(record, p, coeff, q, sign, r),
                                    },
                                    args.top_n,
                                )

    rng = random.Random(args.seed)
    mc_better_or_equal = 0
    if records and args.monte_carlo_samples:
        for _ in range(args.monte_carlo_samples):
            record = rng.choice(records)
            p = rng.choice(list(exponents))
            coeff = rng.choice(list(coeffs))
            q = rng.choice(list(exponents))
            sign = rng.choice(signs)
            r = rng.choice(list(exponents))
            if abs(evaluate(record, p, coeff, q, sign, r) - TARGET) <= HYPOTHESIS_ABS_RESIDUAL:
                mc_better_or_equal += 1

    summary = {
        "target": str(TARGET),
        "hypothesis_value": str(HYPOTHESIS_VALUE),
        "hypothesis_abs_residual": str(HYPOTHESIS_ABS_RESIDUAL),
        "max_d": args.max_d,
        "fields_evaluated": len(records),
        "max_exponent": args.max_exponent,
        "max_abs_torsion_coeff": args.max_coeff,
        "total_expressions": total,
        "strictly_better_than_hypothesis": strictly_better,
        "better_or_equal_to_hypothesis": better_or_equal,
        "hypothesis_rank": strictly_better + 1,
        "exhaustive_fraction_better_or_equal": sci_decimal(Decimal(better_or_equal) / Decimal(total)),
        "monte_carlo_seed": args.seed,
        "monte_carlo_samples": args.monte_carlo_samples,
        "monte_carlo_better_or_equal": mc_better_or_equal,
        "monte_carlo_fraction_better_or_equal": (
            sci_decimal(Decimal(mc_better_or_equal) / Decimal(args.monte_carlo_samples))
            if args.monte_carlo_samples
            else "0"
        ),
    }
    write_summary(summary, args.summary_out)
    write_top(top_rows, args.top_out)


if __name__ == "__main__":
    main()
