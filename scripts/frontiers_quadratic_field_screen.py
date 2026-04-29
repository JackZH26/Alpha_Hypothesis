#!/usr/bin/env python3
"""Pilot screen for the quadratic-field leading scale used in the Frontiers draft.

The script is intentionally small and dependency-free.  It is not a proof of
statistical significance; it records the reproducible field-level screen used
in the manuscript's look-elsewhere discussion.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from math import gcd, isqrt, sqrt


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


def fundamental_unit_float(d: int, limit: int = 20000) -> float:
    best = None
    if d % 4 == 1:
        for b in range(1, limit + 1):
            for sign in (-4, 4):
                a2 = d * b * b + sign
                a = isqrt(a2)
                if a * a == a2 and (a - b) % 2 == 0:
                    unit = (a + b * sqrt(d)) / 2
                    if unit > 1 and (best is None or unit < best):
                        best = unit
                        break
            if best is not None:
                return best
    else:
        for b in range(1, limit + 1):
            for sign in (-1, 1):
                a2 = d * b * b + sign
                a = isqrt(a2)
                if a * a == a2:
                    unit = a + b * sqrt(d)
                    if unit > 1 and (best is None or unit < best):
                        best = unit
                        break
            if best is not None:
                return best
    raise RuntimeError(f"no unit found for d={d}")


def main() -> None:
    target = Decimal("137.03599916476563934505723564140")
    print("d,D,epsilon,L(-3,chi_D),3!/zeta_K(-3),leading_scale,residual")
    for d in [n for n in range(2, 31) if is_squarefree(n)]:
        D = fundamental_discriminant(d)
        eps = fundamental_unit_float(d)
        Lm3 = L_minus_3_for_quadratic_discriminant(D)
        if Lm3 == 0:
            continue
        scale = Fraction(720, 1) / Lm3
        leading = Decimal(scale.numerator) / Decimal(scale.denominator)
        eps_dec = Decimal(str(eps))
        leading = leading / (eps_dec * eps_dec)
        residual = leading - target
        print(
            f"{d},{D},{eps:.12g},{Lm3},{float(scale):.12g},"
            f"{leading:.12f},{residual:.12f}"
        )


if __name__ == "__main__":
    main()
