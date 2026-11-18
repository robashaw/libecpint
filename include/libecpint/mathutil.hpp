/*
 *      Copyright (c) 2020 Robert Shaw
 *		This file is a part of Libecpint.
 *
 *      Permission is hereby granted, free of charge, to any person obtaining
 *      a copy of this software and associated documentation files (the
 *      "Software"), to deal in the Software without restriction, including
 *      without limitation the rights to use, copy, modify, merge, publish,
 *      distribute, sublicense, and/or sell copies of the Software, and to
 *      permit persons to whom the Software is furnished to do so, subject to
 *      the following conditions:
 *
 *      The above copyright notice and this permission notice shall be
 *      included in all copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *      NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *      LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *      OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *      WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef MATHUTILHEADERDEF
#define MATHUTILHEADERDEF

/**
 * \file  mathutil.hpp
 * \brief Mathematical constants, special function tables, and utility functions
 */

#include <cmath>
#include <numeric>
#include <vector>

#include "multiarr.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_EULER
#define M_EULER 2.71828182845904523536
#endif

#define MAX_FAC 100   ///< the maximum factorial needed
#define MAX_DFAC 200  ///< the maximum double factorial needed

namespace libecpint {

constexpr double ROOT_PI = 1.772453850905516;  ///< square root of PI
constexpr double SINH_1 = 1.1752011936;
constexpr int MAX_POW = 24;
/// Largest positive exponent n for which FAST_POW[n] is a valid power routine (entries above are
/// 1/z, 1/z^2)
constexpr int FAST_POW_MAX_EXPONENT = 24;

extern double FAC[MAX_FAC];    ///< Array of factorials
extern double DFAC[MAX_DFAC];  ///< Array of double factorials

double pow_m2(double);
double pow_m1(double);
double pow_0(double);
double pow_1(double);
double pow_2(double);
double pow_3(double);
double pow_4(double);
double pow_5(double);
double pow_6(double);
double pow_7(double);
double pow_8(double);
double pow_9(double);
double pow_10(double);
double pow_11(double);
double pow_12(double);
double pow_13(double);
double pow_14(double);
double pow_15(double);
double pow_16(double);
double pow_17(double);
double pow_18(double);
double pow_19(double);
double pow_20(double);
double pow_21(double);
double pow_22(double);
double pow_23(double);
double pow_24(double);

/// Array of function pointers to hand-coded x**n routines
static double (*FAST_POW[MAX_POW + 3])(double){
    pow_0,  pow_1,  pow_2,  pow_3,  pow_4,  pow_5,  pow_6,  pow_7,  pow_8,
    pow_9,  pow_10, pow_11, pow_12, pow_13, pow_14, pow_15, pow_16, pow_17,
    pow_18, pow_19, pow_20, pow_21, pow_22, pow_23, pow_24, pow_m1, pow_m2};

/**
 * Gamma function tabulation, where GAMMA[i] = Gamma((i+1)/2)
 * e.g. GAMMA[0] = Gamma(1/2) = sqrt(Pi), GAMMA[1] = 0! = 1, GAMMA[2] = Gamma(3/2) = sqrt(Pi) / 2,
 * etc.
 */
const double GAMMA[30] = {1.772453850905516,  1.0,          0.886226925452758,  1.0,
                          1.329340388179137,  2.0,          3.3233509704478426, 6.0,
                          11.631728396567448, 24.0,         52.34277778455352,  120.0,
                          287.88527781504433, 720.0,        1871.2543057977884, 5040.0,
                          14034.407293483413, 40320.0,      119292.46199460901, 362880.0,
                          1133278.3889487856, 3628800.0,    11899423.083962249, 39916800.0,
                          136843365.46556586, 479001600.0,  1710542068.3195732, 6227020800.0,
                          23092317922.31424,  87178291200.0};

/**
 * Calculates real spherical harmonics S_lm(theta, phi) for all l, m up to lmax, writing them
 * into a caller-provided buffer (resized as needed) to avoid per-call allocation in hot loops.
 * @param lmax - the maximum angular momentum needed
 * @param x - cos(theta), where theta is the polar angle in spherical coordinates
 * @param phi - the azimuth angle in spherical coordinates
 * @param values - matrix S(l, l+m) of the spherical harmonic values, filled on return
 */
void realSphericalHarmonics(int lmax, double x, double phi, TwoIndex<double>& values);

/**
 * As above, with sin(theta) given as well. Computing it from the Cartesian components, as
 * sqrt(x^2 + y^2) / r, avoids the cancellation in sqrt(1 - cos^2(theta)) near the poles, whose
 * error (up to ~1e-8 for a direction on the axis whose cosine is rounded) would be amplified
 * in the integrals.
 * @param sx - sin(theta) >= 0
 */
void realSphericalHarmonics(int lmax, double x, double sx, double phi, TwoIndex<double>& values);

/**
 * Calculates real spherical harmonics S_lm(theta, phi) for all l, m up to lmax
 * @param lmax - the maximum angular momentum needed
 * @param x - cos(theta), where theta is the polar angle in spherical coordinates
 * @param phi - the azimuth angle in spherical coordinates
 * @return a matrix S(l, l+m) of the spherical harmonic values
 */
TwoIndex<double> realSphericalHarmonics(int lmax, double x, double phi);

/**
 * @param mat - a reference to a TwoIndex array
 * @return the Frobenius norm of mat
 */
double frobenius_norm(const TwoIndex<double>& mat);

/**
 * Initialises the global factorial and double factorial arrays
 */
void initFactorials();
}  // namespace libecpint

#endif
