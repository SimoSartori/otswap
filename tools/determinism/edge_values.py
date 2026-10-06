# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Write tests/detmath_reference.h: the reference values of the edge-case
tests of src/detmath.h and of the pixel lookup in src/arc.h.

Every function value is computed with mpmath at a precision far above
double, then rounded to the nearest double, ties to even, subnormals
included; it is the correctly rounded result. Special cases (signed zeros,
infinities, NaN) are written from C99 Annex F rather than computed.

The pixel references replay Healpix's vec2pix, ang2pix and loc2pix in
Python, whose floats are IEEE doubles evaluated one operation at a time,
like the C++ built with -ffp-contract=off. Only atan2, cos and sin are taken
from mpmath, correctly rounded. A test failure on a boundary point therefore
means that the deterministic function did not return the correctly rounded
value there, or that loc2pix was not compiled as written.

Every double is written as an exact hexadecimal literal.

Run by hand when the cases change; it needs mpmath:

    python -B tools/determinism/edge_values.py tests/detmath_reference.h
"""

import math
import sys

import mpmath
from mpmath import mp, mpf

# Enough bits to reduce trigonometric arguments up to 2^1024 exactly and
# still hold more than a thousand bits of the result.
mp.prec = 2400

INF = math.inf
NAN = math.nan
TINY = 5e-324                       # 2^-1074, the smallest subnormal
DBL_MAX = sys.float_info.max


def up(x):
    return math.nextafter(x, INF)


def down(x):
    return math.nextafter(x, -INF)


def rounded(v):
    """The double nearest to the mpf v, ties to even; subnormals and
    overflow handled without double rounding."""
    if mpmath.isnan(v):
        return NAN
    if mpmath.isinf(v):
        return INF if v > 0 else -INF
    if v == 0:
        return 0.0
    sign, man, exp, bc = v._mpf_
    msb = exp + bc - 1
    lsb = max(msb - 52, -1074)
    if exp >= lsb:
        q = man << (exp - lsb)
    else:
        shift = lsb - exp
        q, rem = man >> shift, man & ((1 << shift) - 1)
        half = 1 << (shift - 1)
        if rem > half or (rem == half and q & 1):
            q += 1
    if q and lsb + q.bit_length() > 1024:
        r = INF
    else:
        r = math.ldexp(q, lsb)
    return -r if sign else r


def lit(x):
    """An exact C++ literal for the double x."""
    if math.isnan(x):
        return "kNaN"
    if math.isinf(x):
        return "kInf" if x > 0 else "-kInf"
    return x.hex()


# ---------------------------------------------------------------- functions

def F(name, *args, special=None):
    return (name, args, special)


PI = rounded(mp.pi)
HALF_PI = rounded(mp.pi / 2)
QUARTER_PI = rounded(mp.pi / 4)


def multiples_of_half_pi():
    xs = []
    for k in (1, 2, 3, 4, 5, 6, 7, 8, 9, 100, 1000003):
        x = rounded(k * mp.pi / 2)
        xs += [x, -x]
    for k in (1, 2):
        x = rounded(k * mp.pi / 2)
        xs += [up(x), down(x)]
    for k in (3, 5, 7, 9):        # the argument reduction's own thresholds
        xs.append(rounded(k * mp.pi / 4))
    return xs


TRIG_ARGS = (
    [TINY, -TINY, 2.0 ** -30, 2.0 ** -27, up(2.0 ** -27), 0.5, 1.0, -1.0, 2.0]
    + multiples_of_half_pi()
    + [1647099.0, 1647100.0, 2.0 ** 20 * PI, 1e22, -1e22, 2.0 ** 60, 1e300,
       6381956970095103.0 * 2.0 ** 797, 2.0 ** 1023, DBL_MAX, -DBL_MAX]
)

SIN = ([F("sin", 0.0, special=0.0), F("sin", -0.0, special=-0.0),
        F("sin", INF, special=NAN), F("sin", -INF, special=NAN), F("sin", NAN, special=NAN)]
       + [F("sin", x) for x in TRIG_ARGS])

COS = ([F("cos", 0.0, special=1.0), F("cos", -0.0, special=1.0),
        F("cos", INF, special=NAN), F("cos", -INF, special=NAN), F("cos", NAN, special=NAN)]
       + [F("cos", x) for x in TRIG_ARGS])

O_THRESHOLD = float.fromhex("0x1.62e42fefa39efp+9")    # 709.78..., largest finite exp
DBL_MIN_ARG = float.fromhex("-0x1.6232bdd7abcd2p+9")   # -708.39..., exp near DBL_MIN
U_THRESHOLD = float.fromhex("-0x1.74910d52d3051p+9")   # -745.13..., smallest non-zero exp

EXP = ([F("exp", 0.0, special=1.0), F("exp", -0.0, special=1.0),
        F("exp", INF, special=INF), F("exp", -INF, special=0.0), F("exp", NAN, special=NAN),
        F("exp", up(O_THRESHOLD), special=INF), F("exp", 1000.0, special=INF),
        F("exp", -1000.0, special=0.0)]
       + [F("exp", x) for x in (
           1.0, -1.0, 0.5, TINY, -TINY, 2.0 ** -54, -(2.0 ** -54), 1e-300,
           rounded(mp.log(2) / 2), rounded(3 * mp.log(2) / 2), 700.0,
           O_THRESHOLD, down(O_THRESHOLD),
           DBL_MIN_ARG, up(DBL_MIN_ARG), down(DBL_MIN_ARG),
           -720.0, -740.0, -744.0, U_THRESHOLD, up(U_THRESHOLD), down(U_THRESHOLD), -745.2)])

LOG1P = ([F("log1p", 0.0, special=0.0), F("log1p", -0.0, special=-0.0),
          F("log1p", -1.0, special=-INF), F("log1p", down(-1.0), special=NAN),
          F("log1p", -2.0, special=NAN), F("log1p", -INF, special=NAN),
          F("log1p", INF, special=INF), F("log1p", NAN, special=NAN)]
         + [F("log1p", x) for x in (
             TINY, -TINY, 2.0 ** -60, 2.0 ** -54, 2.0 ** -29, up(2.0 ** -29),
             1e-10, -1e-10, 1e-5, -1e-5,
             rounded(mp.sqrt(2) - 1), rounded(mp.sqrt(2) / 2 - 1), 1.0, -0.5,
             up(-1.0), -0.999999, -0.99999999999, 1e300, DBL_MAX)])

ASIN = ([F("asin", 0.0, special=0.0), F("asin", -0.0, special=-0.0),
         F("asin", up(1.0), special=NAN), F("asin", 2.0, special=NAN),
         F("asin", INF, special=NAN), F("asin", NAN, special=NAN)]
        + [F("asin", x) for x in (
            TINY, 2.0 ** -30, 2.0 ** -27, 0.5, -0.5, 0.975, up(0.975), down(0.975),
            down(1.0), 1.0, -1.0)])

THREE_QUARTER_PI = rounded(3 * mp.pi / 4)

ATAN2 = ([F("atan2", y, x, special=s) for (y, x, s) in (
            (0.0, 0.0, 0.0), (-0.0, 0.0, -0.0), (0.0, -0.0, PI), (-0.0, -0.0, -PI),
            (0.0, -1.0, PI), (-0.0, -1.0, -PI), (0.0, 1.0, 0.0), (-0.0, 1.0, -0.0),
            (0.0, -INF, PI), (-0.0, -INF, -PI), (0.0, INF, 0.0), (-0.0, INF, -0.0),
            (1.0, 0.0, HALF_PI), (1.0, -0.0, HALF_PI),
            (-1.0, 0.0, -HALF_PI), (-1.0, -0.0, -HALF_PI),
            (1.0, -INF, PI), (-1.0, -INF, -PI), (1.0, INF, 0.0), (-1.0, INF, -0.0),
            (INF, 1.0, HALF_PI), (-INF, 1.0, -HALF_PI), (INF, -1.0, HALF_PI),
            (INF, INF, QUARTER_PI), (-INF, INF, -QUARTER_PI),
            (INF, -INF, THREE_QUARTER_PI), (-INF, -INF, -THREE_QUARTER_PI),
            (NAN, 1.0, NAN), (1.0, NAN, NAN), (NAN, NAN, NAN))]
         + [F("atan2", y, x) for (y, x) in (
             (1.0, 1.0), (1.0, -1.0), (-1.0, -1.0), (-1.0, 1.0),
             (3.0, 4.0), (4.0, -3.0), (-3.0, -4.0), (-4.0, 3.0),
             (1.0, up(1.0)), (up(1.0), 1.0),
             (1e-300, 1e300), (-1e-300, 1e300), (1e-300, -1e300), (-1e-300, -1e300),
             (1e300, 1e-300), (TINY, TINY), (1.0, TINY), (TINY, -1.0),
             (-1e-17, 1.0), (1e-17, -1.0), (0.7071067811865476, -0.7071067811865476),
             (DBL_MAX, DBL_MAX), (-DBL_MAX, TINY))])

POW = ([F("pow", x, y, special=s) for (x, y, s) in (
          (2.0, 0.0, 1.0), (2.0, -0.0, 1.0), (NAN, 0.0, 1.0), (NAN, -0.0, 1.0),
          (1.0, NAN, 1.0), (1.0, INF, 1.0), (1.0, -INF, 1.0), (1.0, 3.5, 1.0),
          (-1.0, INF, 1.0), (-1.0, -INF, 1.0),
          (0.0, -3.0, INF), (-0.0, -3.0, -INF), (0.0, -2.0, INF), (-0.0, -2.0, INF),
          (0.0, -0.5, INF), (-0.0, -0.5, INF), (0.0, -INF, INF), (-0.0, -INF, INF),
          (0.0, 3.0, 0.0), (-0.0, 3.0, -0.0), (0.0, 2.0, 0.0), (-0.0, 2.0, 0.0),
          (-0.0, 0.5, 0.0), (0.0, INF, 0.0), (-0.0, INF, 0.0),
          (-2.0, 0.5, NAN), (-2.0, 1.0 / 3.0, NAN), (-INF, 0.5, INF),
          (0.5, -INF, INF), (-0.5, -INF, INF), (2.0, -INF, 0.0), (-2.0, -INF, 0.0),
          (0.5, INF, 0.0), (-0.5, INF, 0.0), (2.0, INF, INF), (-2.0, INF, INF),
          (-INF, -3.0, -0.0), (-INF, -2.0, 0.0), (-INF, 3.0, -INF), (-INF, 2.0, INF),
          (INF, -1.0, 0.0), (INF, -0.5, 0.0), (INF, 2.0, INF), (INF, 0.5, INF),
          (NAN, 1.0, NAN), (2.0, NAN, NAN), (NAN, NAN, NAN),
          (2.0, 1024.0, INF), (-2.0, 1025.0, -INF), (2.0, -1080.0, 0.0), (-2.0, -1081.0, -0.0))]
       + [F("pow", x, y) for (x, y) in (
           (2.0, 0.5), (4.0, 0.5), (2.0, 10.0), (-2.0, 3.0), (-2.0, -3.0), (-2.0, 2.0),
           (10.0, -5.0), (10.0, 22.0), (2.0, -1074.0), (2.0, -1075.0), (0.5, 1074.0),
           (2.0, 1023.0), (2.0, 1023.5), (up(1.0), 2.0 ** 53), (down(1.0), 2.0 ** 60),
           (1.0000001, 1e10), (1e-4, -1.0 / 3.0), (3.7e-5, -1.0 / 3.0), (2.3e-7, -1.0 / 3.0),
           (0.5, 0.6), (0.25, -0.3), (0.9, -2.7), (DBL_MAX, 0.5), (TINY, 0.5), (TINY, -0.5))])

FUNCTIONS = [("sin", SIN), ("cos", COS), ("exp", EXP), ("log1p", LOG1P),
             ("asin", ASIN), ("atan2", ATAN2), ("pow", POW)]

MP = {
    "sin": lambda x: mp.sin(mpf(x)),
    "cos": lambda x: mp.cos(mpf(x)),
    "exp": lambda x: mp.exp(mpf(x)),
    "log1p": lambda x: mp.log1p(mpf(x)),
    "asin": lambda x: mp.asin(mpf(x)),
    "atan2": lambda y, x: mp.atan2(mpf(y), mpf(x)),
    "pow": lambda x, y: mp.power(mpf(x), mpf(y)),
}


def evaluate(case):
    name, args, special = case
    if special is not None:
        return special, True
    v = MP[name](*args)
    if isinstance(v, mpmath.mpc):
        raise ValueError(f"{name}{args} is complex; give it as a special case")
    r = rounded(v)
    if r == 0.0 and v < 0:
        r = -0.0
    return r, False


# ------------------------------------------------------------------- pixels

INV_HALFPI = 0.6366197723675813430755350534900574
TWOTHIRD = 2.0 / 3.0
NSIDES = (1, 8, 16, 1024, 8192)
# Checked on the 64-bit base only, which holds their pixel counts.
LARGE_NSIDES = (1 << 20, 1 << 29)


def fmodulo(v1, v2):
    if v1 >= 0:
        return v1 if v1 < v2 else math.fmod(v1, v2)
    tmp = math.fmod(v1, v2) + v2
    return 0.0 if tmp == v2 else tmp


def spread_bits(v):
    r = 0
    for b in range(32):
        r |= ((v >> b) & 1) << (2 * b)
    return r


def xyf2nest(ix, iy, face, order):
    return (face << (2 * order)) + spread_bits(ix) + (spread_bits(iy) << 1)


def trunc(x):
    return int(x)


def loc2pix(nside, nest, z, phi, sth, have_sth):
    """healpix_base.cc, T_Healpix_Base<I>::loc2pix, operation for
    operation. Python's integers are exact, so this is the 64-bit base's
    lookup, and the 32-bit one's wherever the pixel count fits an int."""
    order = nside.bit_length() - 1
    ns = float(nside)
    za = abs(z)
    tt = fmodulo(phi * INV_HALFPI, 4.0)
    npix = 12 * nside * nside
    ncap = 2 * (nside * nside - nside)
    if not nest:
        if za <= TWOTHIRD:
            nl4 = 4 * nside
            temp1 = ns * (0.5 + tt)
            temp2 = ns * z * 0.75
            jp = trunc(temp1 - temp2)
            jm = trunc(temp1 + temp2)
            ir = nside + 1 + jp - jm
            kshift = 1 - (ir & 1)
            t1 = jp + jm - nside + kshift + 1 + nl4 + nl4
            ip = ((t1 >> 1) & (nl4 - 1)) if order > 0 else ((t1 >> 1) % nl4)
            return ncap + (ir - 1) * nl4 + ip
        tp = tt - trunc(tt)
        tmp = ns * math.sqrt(3 * (1 - za)) if (za < 0.99 or not have_sth) \
            else ns * sth / math.sqrt((1. + za) / 3.)
        jp = trunc(tp * tmp)
        jm = trunc((1.0 - tp) * tmp)
        ir = jp + jm + 1
        ip = trunc(tt * ir)
        assert 0 <= ip < 4 * ir
        return 2 * ir * (ir - 1) + ip if z > 0 else npix - 2 * ir * (ir + 1) + ip
    if za <= TWOTHIRD:
        temp1 = ns * (0.5 + tt)
        temp2 = ns * (z * 0.75)
        jp = trunc(temp1 - temp2)
        jm = trunc(temp1 + temp2)
        ifp = jp >> order
        ifm = jm >> order
        face = (ifp | 4) if ifp == ifm else (ifp if ifp < ifm else ifm + 8)
        ix = jm & (nside - 1)
        iy = nside - (jp & (nside - 1)) - 1
        return xyf2nest(ix, iy, face, order)
    ntt = min(3, trunc(tt))
    tp = tt - ntt
    tmp = ns * math.sqrt(3 * (1 - za)) if (za < 0.99 or not have_sth) \
        else ns * sth / math.sqrt((1. + za) / 3.)
    jp = min(trunc(tp * tmp), nside - 1)
    jm = min(trunc((1.0 - tp) * tmp), nside - 1)
    return xyf2nest(nside - jm - 1, nside - jp - 1, ntt, order) if z >= 0 \
        else xyf2nest(jp, jm, ntt + 8, order)


def vec2pix(nside, nest, x, y, z):
    xl = 1. / math.sqrt(x * x + y * y + z * z)
    phi = 0.0 if (x == 0. and y == 0.) else atan2_reference(y, x)
    nz = z * xl
    if abs(nz) > 0.99:
        return loc2pix(nside, nest, nz, phi, math.sqrt(x * x + y * y) * xl, True)
    return loc2pix(nside, nest, nz, phi, 0, False)


def ang2pix(nside, nest, theta, phi):
    c = rounded(mp.cos(mpf(theta)))
    s = rounded(mp.sin(mpf(theta)))
    if theta < 0.01 or theta > 3.14159 - 0.01:
        return loc2pix(nside, nest, c, phi, s, True)
    return loc2pix(nside, nest, c, phi, 0., False)


def atan2_reference(y, x):
    """atan2 correctly rounded, x and y not both zero. mpmath has no signed
    zero, so a zero y is resolved here: the result is 0 or pi, with the sign
    of y."""
    if y == 0.0:
        return math.copysign(0.0 if x > 0 else PI, y)
    r = rounded(mp.atan2(mpf(y), mpf(x)))
    return math.copysign(r, y) if r == 0.0 else r


S75 = math.sqrt(0.75)
S375 = math.sqrt(0.375)
S5_3 = rounded(mp.sqrt(5) / 3)
S099 = math.sqrt(1 - 0.99 * 0.99)

VECTORS = [
    ("north pole", 0.0, 0.0, 1.0), ("south pole", 0.0, 0.0, -1.0),
    ("north pole, not normalised", 0.0, 0.0, 7.5), ("north pole, negative zeros", -0.0, -0.0, 1.0),
    ("next to the north pole", 1e-17, 0.0, 1.0), ("next to the north pole, third quadrant", -1e-17, -1e-17, 1.0),
    ("next to the south pole", 1e-9, 1e-9, -1.0),
    ("RA 0, equator", 1.0, 0.0, 0.0), ("RA -0, equator", 1.0, -0.0, 0.0),
    ("RA 0, z 0.3", 1.0, 0.0, 0.3), ("RA just below 2 pi, by a subnormal", 1.0, -TINY, 0.0),
    ("RA just below 2 pi, z 0.3", 1.0, -1e-17, 0.3), ("RA just below 2 pi, polar", 1.0, -1e-17, -0.9),
    ("RA just above 0", 1.0, 1e-17, 0.5),
    ("RA pi", -1.0, 0.0, 0.0), ("RA -pi", -1.0, -0.0, 0.0), ("RA just past pi", -1.0, -1e-300, 0.2),
    ("RA pi/2", 0.0, 1.0, 0.0), ("RA -pi/2", 0.0, -1.0, 0.0), ("RA pi/2, x -0", -0.0, 1.0, 0.1),
    ("RA pi/4", 1.0, 1.0, 0.0), ("RA 3pi/4", -1.0, 1.0, 0.0), ("RA -3pi/4", -1.0, -1.0, 0.0),
    ("RA -pi/4", 1.0, -1.0, 0.0), ("RA pi/4, z 0.5", 1.0, 1.0, 0.5), ("RA pi/4, z -0.5", 1.0, 1.0, -0.5),
    ("3-4-5", 3.0, 4.0, 0.0), ("diagonal", 1.0, 1.0, 1.0),
    ("z 2/3", S5_3, 0.0, TWOTHIRD), ("z just above 2/3", S5_3, 0.0, up(TWOTHIRD)),
    ("z just below 2/3", S5_3, 0.0, down(TWOTHIRD)), ("z -2/3", S5_3, 0.0, -TWOTHIRD),
    ("z 0.99", S099, 0.0, 0.99), ("z just above 0.99", S099, 0.0, up(0.99)),
    ("ring edge at z 0.5 for NSIDE 8", S75, 0.0, 0.5), ("ring edge at z 0.5, RA pi/4", S375, S375, 0.5),
    ("ring edge at z -0.5", S75, 0.0, -0.5),
]

THETA_2_3 = rounded(mp.acos(mpf(2) / 3))
THIRD_PI = rounded(mp.pi / 3)
TWO_PI = rounded(2 * mp.pi)

ANGLES = [
    ("north pole", 0.0, 0.0), ("north pole, RA 1", 0.0, 1.0),
    ("south pole", PI, 0.0), ("south pole, RA 2", PI, 2.0),
    ("theta 0.01, where the sine is first used", 0.01, 0.5),
    ("theta just below 0.01", down(0.01), 0.5),
    ("theta 3.14159 - 0.01", 3.14159 - 0.01, 0.5), ("theta just above 3.14159 - 0.01", up(3.14159 - 0.01), 0.5),
    ("equator, RA 0", HALF_PI, 0.0), ("equator, RA -0", HALF_PI, -0.0),
    ("equator, RA subnormal", HALF_PI, TINY), ("equator, RA pi/4", HALF_PI, QUARTER_PI),
    ("equator, RA pi/2", HALF_PI, HALF_PI), ("equator, RA pi", HALF_PI, PI),
    ("equator, RA 3pi/2", HALF_PI, rounded(3 * mp.pi / 2)),
    ("equator, RA 2pi", HALF_PI, TWO_PI), ("equator, RA just below 2pi", HALF_PI, down(TWO_PI)),
    ("equator, just south", up(HALF_PI), 0.0), ("equator, just north", down(HALF_PI), 0.0),
    ("z 0.5, RA 0", THIRD_PI, 0.0), ("z 0.5, RA pi/4", THIRD_PI, QUARTER_PI),
    ("z 0.5, one ulp north", down(THIRD_PI), 0.0), ("z 0.5, one ulp south", up(THIRD_PI), 0.0),
    ("z 2/3", THETA_2_3, 0.3), ("z 2/3, one ulp north", down(THETA_2_3), 0.3),
    ("z 2/3, one ulp south", up(THETA_2_3), 0.3),
    ("z -0.5", rounded(2 * mp.pi / 3), 1.0),
]


# ------------------------------------------------------------------- output

def main(path):
    out = []
    w = out.append
    w("// SPDX-License-Identifier: GPL-2.0-or-later")
    w("// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it")
    w("//")
    w("// Generated by tools/determinism/edge_values.py (mpmath " + mpmath.__version__ +
      "); do not edit.")
    w("// Reference values of tests/test_detmath.cpp: correctly rounded function")
    w("// values, and the pixels Healpix's lookup gives with correctly rounded")
    w("// atan2, cos and sin.")
    w("")
    w("#ifndef OTSWAP_TESTS_DETMATH_REFERENCE_H")
    w("#define OTSWAP_TESTS_DETMATH_REFERENCE_H")
    w("")
    w("#include <limits>")
    w("")
    w("namespace detref {")
    w("")
    w("  const double kInf = std::numeric_limits<double>::infinity();")
    w("  const double kNaN = std::numeric_limits<double>::quiet_NaN();")
    w("")
    w("  // special: a value C99 Annex F prescribes, to be matched exactly (sign")
    w("  // of zero included; any NaN for NaN). Otherwise the correctly rounded")
    w("  // value, which the function must match within one ulp, exactly when")
    w("  // exact is set.")
    w("  struct Case { double x, y, expected; bool special, exact; };")
    w("")
    for name, cases in FUNCTIONS:
        w(f"  const Case k_{name}[] = {{")
        for case in cases:
            args = case[1]
            value, special = evaluate(case)
            exact = special or (not math.isnan(value) and not math.isinf(value)
                                and mpf(value) == MP[name](*args))
            y = args[1] if len(args) > 1 else 0.0
            w(f"    {{{lit(args[0])}, {lit(y)}, {lit(value)}, "
              f"{'true' if special else 'false'}, {'true' if exact else 'false'}}},")
        w("  };")
        w("")
    w("  const int kNsides[] = {" + ", ".join(str(n) for n in NSIDES) + "};")
    w(f"  constexpr int kNNsides = {len(NSIDES)};")
    w("")
    w("  struct VecCase { const char* label; double x, y, z; int ring[kNNsides], nest[kNNsides]; };")
    w("  struct AngCase { const char* label; double theta, phi; int ring[kNNsides], nest[kNNsides]; };")
    w("")
    w("  const VecCase kVectors[] = {")
    for label, x, y, z in VECTORS:
        ring = [vec2pix(n, False, x, y, z) for n in NSIDES]
        nest = [vec2pix(n, True, x, y, z) for n in NSIDES]
        w(f'    {{"{label}", {lit(x)}, {lit(y)}, {lit(z)},')
        w(f"     {{{', '.join(map(str, ring))}}}, {{{', '.join(map(str, nest))}}}}},")
    w("  };")
    w("")
    w("  const AngCase kAngles[] = {")
    for label, theta, phi in ANGLES:
        ring = [ang2pix(n, False, theta, phi) for n in NSIDES]
        nest = [ang2pix(n, True, theta, phi) for n in NSIDES]
        w(f'    {{"{label}", {lit(theta)}, {lit(phi)},')
        w(f"     {{{', '.join(map(str, ring))}}}, {{{', '.join(map(str, nest))}}}}},")
    w("  };")
    w("")
    w("  // The same points at NSIDE 2^20 and 2^29, on the 64-bit base.")
    w("  const long long kLargeNsides[] = {" + ", ".join(str(n) for n in LARGE_NSIDES) + "};")
    w(f"  constexpr int kNLargeNsides = {len(LARGE_NSIDES)};")
    w("")
    w("  struct LargeVecCase { const char* label; double x, y, z; "
      "long long ring[kNLargeNsides], nest[kNLargeNsides]; };")
    w("  struct LargeAngCase { const char* label; double theta, phi; "
      "long long ring[kNLargeNsides], nest[kNLargeNsides]; };")
    w("")
    w("  const LargeVecCase kLargeVectors[] = {")
    for label, x, y, z in VECTORS:
        ring = [vec2pix(n, False, x, y, z) for n in LARGE_NSIDES]
        nest = [vec2pix(n, True, x, y, z) for n in LARGE_NSIDES]
        w(f'    {{"{label}", {lit(x)}, {lit(y)}, {lit(z)},')
        w(f"     {{{', '.join(map(str, ring))}}}, {{{', '.join(map(str, nest))}}}}},")
    w("  };")
    w("")
    w("  const LargeAngCase kLargeAngles[] = {")
    for label, theta, phi in ANGLES:
        ring = [ang2pix(n, False, theta, phi) for n in LARGE_NSIDES]
        nest = [ang2pix(n, True, theta, phi) for n in LARGE_NSIDES]
        w(f'    {{"{label}", {lit(theta)}, {lit(phi)},')
        w(f"     {{{', '.join(map(str, ring))}}}, {{{', '.join(map(str, nest))}}}}},")
    w("  };")
    w("")
    w("}")
    w("")
    w("#endif")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
