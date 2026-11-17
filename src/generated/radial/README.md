# Closed forms for the type 2 radial integrals (not used)

`unrol_radial.py` unrolls the recurrence relations of Shaw and Hill, J. Chem. Phys. 147, 074108 (2017) into closed forms for the type 2 radial integrals of a primitive triple, in terms of a few base integrals. Up to f6f4db6 the library used them, as generated C++, for angular momenta up to 4. They are kept here so that the idea can be revisited. The library now calculates every radial integral by quadrature (`src/lib/radial_gen.cpp`).

## Why they are not used

**Accuracy.** The closed forms are exact, but in double precision they lose accuracy to cancellation when a basis primitive is diffuse or close to the ECP, i.e. when x = aA or y = bB is small. The base integrals are sums of alternating terms in X1 ∝ exp(p P1²) and X2 ∝ exp(p P2²), which nearly cancel when x is small, and the recurrence coefficients grow like powers of 1/x and 1/y.

- `check_closed_forms.py --digits 16` with aA = 0.2 gives relative errors up to 2e-3, against 2e-13 at aA = 1.7.
- On Ag basis sets, errors reached 1e5 relative for l1 + l2 ≥ 5, and 2% for a tight primitive paired with a moderate one.
- Radials are multiplied downstream by powers of the distances (the binomial factors of the shell displacements and derivatives), so small absolute errors in the radials are not small in the integrals.

**Speed.** An a posteriori error estimate with a quadrature fallback (branch `exp/closed-forms-estimate`, 0bf8f29) made them safe, but then quadrature was still needed for 93% of primitive triples on Ag, so there was no speed-up. Since then, the shared radial grid has made quadrature cheap for all but tight primitives (exponent > 30). It evaluates each shell's contracted radial factor once and reuses it for every partner shell and ECP channel, whereas a closed form costs the same for every primitive triple. The remaining per-triple quadrature for tight primitives, where x is large and the closed forms are accurate, took under 1% of the run time on a 25-atom Ag benchmark (2026-11-13).

## A bug in the generated code used up to f6f4db6

The original script encoded each recurrence coefficient as a string and read the indices back one character at a time, so k ≥ 10 was misread. That changed one generated case, Q(1, 1, 10) (switch label 10110), which was wrong by O(1). This case is reached for l1 = l2 = 1 with N = 8 and an r⁰ ECP term (n = 2 in the ECP file, as in every ecp28mdf term), i.e. for LA + LB ≥ 8: g–g shell pairs, or f–f pairs in first-derivative integrals. The script now keeps coefficients as tuples. With that fixed, all 63 cases up to angular momentum 4 agree with direct quadrature to the working precision (`check_closed_forms.py`, 40 digits). The old generated code differed only in that case.

## Files

- **`unrol_radial.py`**: writes, for each (i, j, k), the coefficients of the base integrals.
  - `--format code` gives C++ switch cases with label i·10000 + j·100 + k, where `values[N - 2]` holds Q_00N, alongside `G1A`, `G1B` and `H2`.
  - `--format simple` gives a readable listing.
  - `--max-am` sets the largest angular momentum (default 4).
- **`check_closed_forms.py`**: evaluates every generated case as the library did, against direct quadrature of the radial integral.
  - The base integrals are computed as in the removed C++ (`radial_gen.part1`, in git history before this commit): F_N and G_N by recurrences in X1 and X2; G1A and G1B by Dawson functions; H2 directly.
  - `--digits` sets the working precision.

Both need sympy; the checker also needs mpmath.

```bash
python unrol_radial.py > cases.txt
python check_closed_forms.py cases.txt                                # exact arithmetic: ~1e-37
python check_closed_forms.py cases.txt --digits 16 1.3 0.1 0.15 2 2.5  # diffuse: ~1e-3
```

## Using them again

- `RadialIntegral::type2` still takes `nbase`, the number of base integrals the closed forms need.
- `RadialIntegral::estimate_type2` is the screening estimate the old generated code used. It is not a bound.
- The natural place for closed forms is the per-triple path for tight primitives in `RadialIntegral::type2`, where x and y are large. They would need a rigorous accuracy criterion; see the error estimate on `exp/closed-forms-estimate`.
- Validate any change against the independent integral reference (the `Suite.Live` tests, `ctest -L reference`; see `tests/README.md`) as well as against the radial-level references (`tests/radial`). Errors that are harmless in a radial integral can be amplified in the final integrals.
