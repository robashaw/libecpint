"""Closed forms for the type 2 radial integrals, unrolled from the recurrence relations of

    [Shaw2017] R. A. Shaw and J. G. Hill, J. Chem. Phys. 147, 074108 (2017),
               doi: 10.1063/1.4986887

Not used by the library at present: see README.md in this directory for why, and for what
would be needed to use the closed forms again.

For a primitive triple (ECP exponent n, basis exponents a and b at distances A and B from the
ECP), with p = n + a + b, x = aA and y = bB,

    Q_ijk = int_0^inf r^k exp(-p r^2) M_i(2xr) M_j(2yr) dr,

where M_l is the modified spherical Bessel function of the first kind. The recurrences
(eqs. 28, 29 and 33 of [Shaw2017]) reduce Q_ijk (i <= j) to the base integrals Q_00N, and those
with N < 1 are eliminated in favour of G1A, G1B and H2 (eqs. 37-39). For each (i, j, k) this
prints the coefficients of the base integrals, either as C++ switch cases (--format code, case
label i*10000 + j*100 + k, with values[N - 2] holding Q_00N for N >= 2), or as a readable
listing (--format simple).

usage: python unrol_radial.py [--max-am 4] [--format code|simple] > output
Requires sympy.
"""

import argparse
import re
import sys

from sympy import Symbol, simplify

x = Symbol('x')  # k_A / 2 = aA
y = Symbol('y')  # k_B / 2 = bB
p = Symbol('p')  # n + a + b


class Qijk:
    """A term Q_ijk and its reduction to base integrals.

    Members:
      i, j, k     - the integral indices
      size        - number of base integrals in the reduction, 2i + j + 1
      start, end  - the lowest and highest N of the base integrals Q_00N
      subq, terms - the Qijk in one step of the recurrence, and the recurrence coefficient of each
      bases       - for each N in [start, end], the coefficient of Q_00N
      f, ga, gb, h - coefficients of the base integrals of each type after sort() and eliminate()
    """

    def __init__(self, i=0, j=0, k=0):
        self.i, self.j, self.k = i, j, k
        self.size = 2 * i + j + 1
        self.start = k - i - j
        self.end = k + i
        self.subq = []
        self.terms = []
        self.bases = [[] for _ in range(self.size)]
        self.f, self.ga, self.gb, self.h = [], [], [], []

    def format_simple(self):
        """Readable listing of the simplified coefficient of each base integral."""
        lines = ["Q %d %d %d" % (self.i, self.j, self.k)]
        for name, coefs, offset in (("F", self.f, 0), ("GB", self.gb, 1), ("GA", self.ga, 1),
                                    ("H", self.h, 2)):
            for n, c in enumerate(coefs):
                N = 2 * n + self.start + offset
                if N > 0:
                    lines.append("%s %d : %s" % (name, N, simplify(c)))
        return "\n".join(lines) + "\n"

    def format_code(self):
        """C++ switch case accumulating result from the base integrals."""
        out = ["case %d: {" % (self.i * 10000 + self.j * 100 + self.k)]

        def emit(c, target, first=False):
            c = simplify(c)
            if c != 0:
                out.append("  result %s (%s) * %s;" % ("=" if first else "+=", powers(str(c)),
                                                       target))

        for n, c in enumerate(self.f):  # F_N = Q_00N, N even
            N = 2 * n + self.start
            if N >= 2:
                emit(c, "values[%d]" % (N - 2), first=(N == 2))
        for n, c in enumerate(self.gb):  # Q_00N, N odd
            N = 2 * n + self.start + 1
            if N == 1:
                emit(c, "G1B")
            elif N > 1:
                emit(c, "values[%d]" % (N - 2))
        for n, c in enumerate(self.ga):
            if 2 * n + self.start + 1 == 1:
                emit(c, "G1A")
        for n, c in enumerate(self.h):
            if 2 * n + self.start + 2 == 2:
                emit(c, "H2")
        out += ["  break;", "}"]
        return "\n".join(out) + "\n"

    def combine(self):
        """Sums the coefficient products collected for each base integral."""
        self.bases = [sum((coefficient(t) for t in terms), 0) for terms in self.bases]

    def sort(self):
        """Sorts the base integrals by type: F for even N, Gb for odd N."""
        for n, c in enumerate(self.bases):
            (self.f if n % 2 == 0 else self.gb).append(c)

    def eliminate(self):
        """Eliminates base integrals with N < 1, as described by eqs. 37-39 of [Shaw2017]."""
        if self.start >= 1:
            return
        # append a zero term if max N < 1, and another if it is exactly 1
        if self.end < 1:
            self.bases.append(0)
        if self.end == 1:
            self.bases.append(0)

        N = self.start
        ix = 0  # start at the lowest N and recur up
        gaix = -1  # where we are up to in the Ga expansion (eq. 38)
        hix = -1  # and the H expansion (eq. 39)
        while N < 1:
            # all base integrals share the first term, X_N: (2p / (N - 1)) X_{N+2}
            z = self.bases[ix]
            self.bases[ix + 2] += (2 * p / (N - 1)) * z
            # and X_{N+1} gains a -k_B / (N - 1)
            self.bases[ix + 1] -= (2 * y / (N - 1)) * z

            if ix % 2 == 0:  # even: a Ga then an H
                if gaix > -1:
                    self.ga[gaix] -= (2 * x / (N - 1)) * z
                else:
                    self.ga.append(-(2 * x / (N - 1)) * z)
                    gaix += 1
                if hix > -1:
                    z = self.h[hix]
                    self.h.append((2 * p / (N - 1)) * z)
                    hix += 1
                    self.ga[gaix] -= (2 * y / (N - 1)) * z
                    self.bases[ix + 1] -= (2 * x / (N - 1)) * z
            else:  # odd: an H then a Ga
                if hix > -1:
                    self.h[hix] -= (2 * x / (N - 1)) * z
                else:
                    self.h.append(-(2 * x / (N - 1)) * z)
                    hix += 1
                if gaix > -1:
                    z = self.ga[gaix]
                    self.ga.append((2 * p / (N - 1)) * z)
                    gaix += 1
                    self.h[hix] -= (2 * y / (N - 1)) * z
                    self.bases[ix + 1] -= (2 * x / (N - 1)) * z
            N += 1
            ix += 1


def coefficient(terms):
    """Product of the recurrence coefficients along one path of the unrolling (eqs. 28-33)."""
    c = 1
    for t in terms:
        if t[0] == "mu":
            _, i, j, k = t
            c *= (2 + j - i - k) / (2 * x)
        elif t[0] == "nu":
            c *= -y / x
        elif t[0] == "xi":
            c *= p / x
        elif t[0] == "rho":
            c *= (1 - 2 * t[1]) / (2 * y)
        elif t[0] == "om":
            c *= -1 / (2 * y)
        # "sigma" and "ups" have coefficient 1
    return c


def powers(s):
    """Writes x**n, y**n and p**n with the squares x2, y2, p2 (declared by the caller)."""
    def rep(m):
        v, n = m.group(1), int(m.group(2))
        f = [v + "2"] * (n // 2) + ([v] if n % 2 else [])
        return f[0] if len(f) == 1 else "(" + "*".join(f) + ")"
    return re.sub(r"\b([xyp])\*\*(\d+)", rep, s)


def unrol(q):
    """Recursively unrolls the recurrences of [Shaw2017], eqs. 28-33, for q."""
    if q.i == 0 and q.j == 0:  # already a base integral
        return
    if q.i > 0:  # reduce i, eq. 28
        q.subq += [Qijk(q.i - 1, q.j, q.k - 1), Qijk(q.i - 1, q.j - 1, q.k),
                   Qijk(q.i - 1, q.j, q.k + 1)]
        q.terms += [("mu", q.i, q.j, q.k), ("nu",), ("xi",)]
    elif q.j > 1:  # reduce j, eq. 29
        q.subq += [Qijk(0, q.j - 2, q.k), Qijk(0, q.j - 1, q.k - 1)]
        q.terms += [("sigma",), ("rho", q.j)]
    else:  # eq. 33
        q.subq += [Qijk(0, 0, q.k), Qijk(0, 0, q.k - 1)]
        q.terms += [("ups",), ("om",)]
    for s in q.subq:
        unrol(s)


def collect(q, Q, path=()):
    """Appends to Q.bases the coefficient path from Q down to each base integral in q."""
    if q.i == 0 and q.j == 0:
        Q.bases[q.k - Q.start].append(path)
    else:
        for s, t in zip(q.subq, q.terms):
            collect(s, Q, path + (t,))


def algebraic_unrol(i, j, k):
    """Q_ijk reduced to base integrals, with simplified coefficients."""
    q = Qijk(i, j, k)
    unrol(q)
    collect(q, q)
    q.combine()
    return q


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--max-am", type=int, default=4,
                        help="largest angular momentum j to unroll (default 4)")
    parser.add_argument("--format", choices=("code", "simple"), default="code")
    args = parser.parse_args()
    L = args.max_am
    for j in range(L + 1):
        for i in range(j + 1):
            for k in range(1, 3 * L + 1 - i - j):
                if (i + j + k) % 2 == 0:
                    q = algebraic_unrol(i, j, k)
                    q.eliminate()
                    q.sort()
                    sys.stdout.write(q.format_code() if args.format == "code" else
                                     q.format_simple())
                    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
