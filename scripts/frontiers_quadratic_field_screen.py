#!/usr/bin/env python3
"""Pilot screen for the quadratic-field arithmetic scales used in the draft.

The script is intentionally small and dependency-free.  It is not a proof of
statistical significance; it records the reproducible field-level screen used
in the manuscript's look-elsewhere discussion.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from math import isqrt


getcontext().prec = 50


def is_squarefree(n: int) -> bool:
    for p in range(2, isqrt(n) + 1):
        if n % (p * p) == 0:
            return False
    return True


def fundamental_discriminant(d: int) -> int:
    return d if d % 4 == 1 else 4 * d


def kronecker(a: int, n: int) -> int:
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        return (-1 if a < 0 else 1) * kronecker(a, -n)

    result = 1
    twos = 0
    while n % 2 == 0:
        twos += 1
        n //= 2

    if twos:
        if a % 2 == 0:
            return 0
        result *= (1 if a % 8 in (1, 7) else -1) ** twos

    if n == 1:
        return result

    a %= n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def bernoulli4(x: Fraction) -> Fraction:
    return x**4 - 2 * x**3 + x**2 - Fraction(1, 30)


def L_minus_3_for_quadratic_discriminant(D: int) -> Fraction:
    f = abs(D)
    total = Fraction(0, 1)
    for a in range(1, f + 1):
        total += kronecker(D, a) * bernoulli4(Fraction(a, f))
    generalized_bernoulli = f**3 * total
    return -generalized_bernoulli / 4


def fundamental_unit_data(d: int, limit: int = 20000) -> tuple[Decimal, int, int, int]:
    """Return epsilon and its exact representation (a + b sqrt(d)) / denom."""
    best = None
    if d % 4 == 1:
        for b in range(1, limit + 1):
            for sign in (-4, 4):
                a2 = d * b * b + sign
                a = isqrt(a2)
                if a * a == a2 and (a - b) % 2 == 0:
                    unit = (Decimal(a) + Decimal(b) * Decimal(d).sqrt()) / 2
                    if unit > 1 and (best is None or unit < best):
                        best = (unit, a, b, 2)
                        break
            if best is not None:
                return best
    else:
        for b in range(1, limit + 1):
            for sign in (-1, 1):
                a2 = d * b * b + sign
                a = isqrt(a2)
                if a * a == a2:
                    unit = Decimal(a) + Decimal(b) * Decimal(d).sqrt()
                    if unit > 1 and (best is None or unit < best):
                        best = (unit, a, b, 1)
                        break
            if best is not None:
                return best
    raise RuntimeError(f"no unit found for d={d}")


def trace_epsilon_squared(d: int, a: int, b: int, denom: int) -> Fraction:
    if denom == 2:
        return Fraction(a * a + d * b * b, 2)
    return Fraction(2 * (a * a + d * b * b), 1)


def main() -> None:
    # CODATA 2022 adjusted central value for alpha(0)^(-1), used as the
    # empirical target in the field-level negative-control screen.
    target = Decimal("137.035999177")
    print(
        "d,D,epsilon,trace_epsilon_squared,L(-3,chi_D),"
        "3!/zeta_K(-3),leading_scale,leading_residual,"
        "uniform_three_term,uniform_residual"
    )
    for d in [n for n in range(2, 31) if is_squarefree(n)]:
        D = fundamental_discriminant(d)
        eps, a, b, denom = fundamental_unit_data(d)
        trace_eps2 = trace_epsilon_squared(d, a, b, denom)
        Lm3 = L_minus_3_for_quadratic_discriminant(D)
        if Lm3 == 0:
            continue
        scale = Fraction(720, 1) / Lm3
        leading = Decimal(scale.numerator) / Decimal(scale.denominator)
        leading = leading / (eps * eps)
        leading_residual = leading - target
        trace_dec = Decimal(trace_eps2.numerator) / Decimal(trace_eps2.denominator)
        uniform = leading - Decimal(2) / (eps**3) + Decimal(1) / ((trace_dec * eps) ** D)
        uniform_residual = uniform - target
        print(
            f"{d},{D},{eps:.12f},{trace_eps2},{Lm3},"
            f"{float(scale):.12g},{leading:.12f},{leading_residual:.12f},"
            f"{uniform:.12f},{uniform_residual:.12f}"
        )


if __name__ == "__main__":
    main()
