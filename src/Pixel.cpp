// SPDX-License-Identifier: GPL-2.0-or-later
/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * This program is free software; you can redistribute it and/or    *
 * modify it under the terms of the GNU General Public License as   *
 * published by the Free Software Foundation; either version 2 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU General Public License for more details.                     *
 *                                                                  *
 * You should have received a copy of the GNU General Public        *
 * License along with this program; if not, write to the Free       *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,      *
 * Boston, MA 02110-1301 USA.                                       *
 ********************************************************************/

/**
 *  @file src/Pixel.cpp
 *
 *  @brief The Healpix pixelisation on a 64-bit pixel index.
 *
 *  Derived from HEALPix 3.82 (Healpix_cxx/healpix_base.cc,
 *  healpix_tables.cc and cxxsupport/math_utils.h), Copyright (C) 2003-2016
 *  Max-Planck-Society, author Martin Reinecke, distributed under the GNU
 *  General Public License, version 2 or later.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>

#include "detmath.h"
#include "internal.h"
#include "pixel.h"

namespace {

  using I = std::int64_t;

  /// Healpix's inv_halfpi and twothird, from the same literals.
  constexpr double kInvHalfPi = 0.6366197723675813430755350534900574;
  constexpr double kTwoThird = 2.0/3.0;

  /// Ring and position offsets of the 12 base faces.
  constexpr int kJrll[12] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
  constexpr int kJpll[12] = { 1,3,5,7,0,2,4,6,1,3,5,7 };

  /// Offsets of the eight neighbours in the face's x and y, and, for a
  /// neighbour across a face edge, the face it lies on and the flips and
  /// swap of x and y that face needs. The rows of the last two are the
  /// directions S, SE, E, SW, centre, NE, W, NW and N, in that order.
  constexpr int kNbXOffset[8] = { -1,-1, 0, 1, 1, 1, 0,-1 };
  constexpr int kNbYOffset[8] = {  0, 1, 1, 1, 0,-1,-1,-1 };
  constexpr int kNbFaceArray[9][12] =
    { {  8, 9,10,11,-1,-1,-1,-1,10,11, 8, 9 },
      {  5, 6, 7, 4, 8, 9,10,11, 9,10,11, 8 },
      { -1,-1,-1,-1, 5, 6, 7, 4,-1,-1,-1,-1 },
      {  4, 5, 6, 7,11, 8, 9,10,11, 8, 9,10 },
      {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11 },
      {  1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4 },
      { -1,-1,-1,-1, 7, 4, 5, 6,-1,-1,-1,-1 },
      {  3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7 },
      {  2, 3, 0, 1,-1,-1,-1,-1, 0, 1, 2, 3 } };
  constexpr int kNbSwapArray[9][3] =
    { { 0,0,3 },
      { 0,0,6 },
      { 0,0,0 },
      { 0,0,5 },
      { 0,0,0 },
      { 5,0,0 },
      { 0,0,0 },
      { 6,0,0 },
      { 3,0,0 } };

  /// v1 modulo v2, in [0, v2), for v2 > 0.
  double fmodulo (const double v1, const double v2)
  {
    if (v1 >= 0) return (v1 < v2) ? v1 : std::fmod(v1, v2);
    const double tmp = std::fmod(v1, v2) + v2;
    return (tmp == v2) ? 0. : tmp;
  }

  /// The integer n with n*n <= arg < (n+1)*(n+1).
  std::uint32_t isqrt (const I arg)
  {
    I res = (I)std::sqrt(double(arg) + 0.5);
    if (arg < (I(1) << 50)) return std::uint32_t(res);
    if (res*res > arg)
      --res;
    else if ((res+1)*(res+1) <= arg)
      ++res;
    return std::uint32_t(res);
  }

  /// a / b for a result in [0, 3].
  I special_div (I a, const I b)
  {
    const I t = (a >= (b<<1));
    a -= t*(b<<1);
    return (t<<1) + (a >= b);
  }

  /// The bits of v moved to the even positions of the result.
  I spread_bits (const int v)
  {
    std::uint64_t x = (std::uint32_t)v;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
    x = (x | (x <<  8)) & 0x00FF00FF00FF00FFull;
    x = (x | (x <<  4)) & 0x0F0F0F0F0F0F0F0Full;
    x = (x | (x <<  2)) & 0x3333333333333333ull;
    x = (x | (x <<  1)) & 0x5555555555555555ull;
    return (I)x;
  }

  /// The even bits of v, packed: the inverse of spread_bits.
  int compress_bits (const I v)
  {
    std::uint64_t x = (std::uint64_t)v & 0x5555555555555555ull;
    x = (x | (x >>  1)) & 0x3333333333333333ull;
    x = (x | (x >>  2)) & 0x0F0F0F0F0F0F0F0Full;
    x = (x | (x >>  4)) & 0x00FF00FF00FF00FFull;
    x = (x | (x >>  8)) & 0x0000FFFF0000FFFFull;
    x = (x | (x >> 16)) & 0x00000000FFFFFFFFull;
    return (int)x;
  }

}


// ============================================================================


void otswap::internal::HealpixBase::set (const I nside, const bool nested)
{
  m_order = -1;
  if ((nside & (nside - 1)) == 0) {
    m_order = 0;
    while ((I(1) << m_order) < nside) ++m_order;
  }
  m_nside = nside;
  m_npface = m_nside*m_nside;
  m_ncap = (m_npface - m_nside) << 1;
  m_npix = 12*m_npface;
  m_nested = nested;
}


// ============================================================================


I otswap::internal::HealpixBase::loc2pix (const double z, const double phi, const double sth,
                                          const bool haveSth) const
{
  const double za = std::fabs(z);
  const double tt = fmodulo(phi*kInvHalfPi, 4.0);

  if (!m_nested) {
    if (za <= kTwoThird) {
      const I nl4 = 4*m_nside;
      const double temp1 = m_nside*(0.5+tt);
      const double temp2 = m_nside*z*0.75;
      const I jp = I(temp1-temp2);
      const I jm = I(temp1+temp2);

      const I ir = m_nside + 1 + jp - jm;
      const I kshift = 1 - (ir&1);

      const I t1 = jp + jm - m_nside + kshift + 1 + nl4 + nl4;
      const I ip = (m_order > 0) ? (t1>>1)&(nl4-1) : ((t1>>1)%nl4);

      return m_ncap + (ir-1)*nl4 + ip;
    }

    const double tp = tt - I(tt);
    const double tmp = ((za < 0.99) || (!haveSth)) ?
                       m_nside*std::sqrt(3*(1-za)) :
                       m_nside*sth/std::sqrt((1.+za)/3.);

    const I jp = I(tp*tmp);
    const I jm = I((1.0-tp)*tmp);

    const I ir = jp + jm + 1;
    const I ip = I(tt*ir);
    if (!((ip >= 0) && (ip < 4*ir)))
      throw Error("internal error in the pixel lookup: in-ring index " + std::to_string(ip) +
                  " outside [0, " + std::to_string(4*ir) + ") at z = " + std::to_string(z) +
                  ", phi = " + std::to_string(phi) + ", NSIDE " + std::to_string(m_nside));

    return (z > 0) ? 2*ir*(ir-1) + ip : m_npix - 2*ir*(ir+1) + ip;
  }

  if (za <= kTwoThird) {
    const double temp1 = m_nside*(0.5+tt);
    const double temp2 = m_nside*(z*0.75);
    const I jp = I(temp1-temp2);
    const I jm = I(temp1+temp2);
    const I ifp = jp >> m_order;
    const I ifm = jm >> m_order;
    const int face = (int)((ifp == ifm) ? (ifp|4) : ((ifp < ifm) ? ifp : (ifm+8)));

    const int ix = (int)(jm & (m_nside-1));
    const int iy = (int)(m_nside - (jp & (m_nside-1)) - 1);
    return xyf2nest(ix, iy, face);
  }

  const int ntt = std::min(3, int(tt));
  const double tp = tt - ntt;
  const double tmp = ((za < 0.99) || (!haveSth)) ?
                     m_nside*std::sqrt(3*(1-za)) :
                     m_nside*sth/std::sqrt((1.+za)/3.);

  I jp = I(tp*tmp);
  I jm = I((1.0-tp)*tmp);
  jp = std::min(jp, m_nside-1);
  jm = std::min(jm, m_nside-1);
  return (z >= 0) ?
    xyf2nest((int)(m_nside-jm-1), (int)(m_nside-jp-1), ntt) : xyf2nest((int)jp, (int)jm, ntt+8);
}


// ============================================================================


void otswap::internal::HealpixBase::neighbors (const I pix, std::array<I, 8>& result) const
{
  int ix, iy, face;
  if (m_nested)
    nest2xyf(pix, ix, iy, face);
  else
    ring2xyf(pix, ix, iy, face);

  const I nsm1 = m_nside - 1;
  if ((ix > 0) && (ix < nsm1) && (iy > 0) && (iy < nsm1)) {
    if (!m_nested)
      for (int m = 0; m < 8; ++m)
        result[(std::size_t)m] = xyf2ring(ix + kNbXOffset[m], iy + kNbYOffset[m], face);
    else {
      const I fpix = I(face) << (2*m_order),
        px0 = spread_bits(ix  ), py0 = spread_bits(iy  ) << 1,
        pxp = spread_bits(ix+1), pyp = spread_bits(iy+1) << 1,
        pxm = spread_bits(ix-1), pym = spread_bits(iy-1) << 1;

      result[0] = fpix+pxm+py0; result[1] = fpix+pxm+pyp;
      result[2] = fpix+px0+pyp; result[3] = fpix+pxp+pyp;
      result[4] = fpix+pxp+py0; result[5] = fpix+pxp+pym;
      result[6] = fpix+px0+pym; result[7] = fpix+pxm+pym;
    }
    return;
  }

  for (int i = 0; i < 8; ++i) {
    int x = ix + kNbXOffset[i], y = iy + kNbYOffset[i];
    int nbnum = 4;
    if (x < 0)
      { x += (int)m_nside; nbnum -= 1; }
    else if (x >= m_nside)
      { x -= (int)m_nside; nbnum += 1; }
    if (y < 0)
      { y += (int)m_nside; nbnum -= 3; }
    else if (y >= m_nside)
      { y -= (int)m_nside; nbnum += 3; }

    const int f = kNbFaceArray[nbnum][face];
    if (f >= 0) {
      const int bits = kNbSwapArray[nbnum][face>>2];
      if (bits&1) x = (int)(m_nside - x - 1);
      if (bits&2) y = (int)(m_nside - y - 1);
      if (bits&4) std::swap(x, y);
      result[(std::size_t)i] = m_nested ? xyf2nest(x, y, f) : xyf2ring(x, y, f);
    }
    else
      result[(std::size_t)i] = -1;
  }
}


// ============================================================================


I otswap::internal::HealpixBase::xyf2nest (const int ix, const int iy, const int face) const
{
  return (I(face) << (2*m_order)) + spread_bits(ix) + (spread_bits(iy) << 1);
}


// ============================================================================


void otswap::internal::HealpixBase::nest2xyf (I pix, int& ix, int& iy, int& face) const
{
  face = (int)(pix >> (2*m_order));
  pix &= (m_npface-1);
  ix = compress_bits(pix);
  iy = compress_bits(pix >> 1);
}


// ============================================================================


I otswap::internal::HealpixBase::xyf2ring (const int ix, const int iy, const int face) const
{
  const I nl4 = 4*m_nside;
  const I jr = (kJrll[face]*m_nside) - ix - iy - 1;

  I nr, n_before;
  bool shifted;
  ring_info_small(jr, n_before, nr, shifted);
  nr >>= 2;
  const I kshift = 1 - shifted;
  I jp = (kJpll[face]*nr + ix - iy + 1 + kshift) / 2;
  if (!(jp <= 4*nr))
    throw Error("internal error in the neighbour lookup: position " + std::to_string(jp) +
                " beyond the " + std::to_string(4*nr) + " pixels of ring " + std::to_string(jr) +
                ", NSIDE " + std::to_string(m_nside));
  if (jp < 1) jp += nl4;

  return n_before + jp - 1;
}


// ============================================================================


void otswap::internal::HealpixBase::ring2xyf (const I pix, int& ix, int& iy, int& face) const
{
  I iring, iphi, kshift, nr;
  const I nl2 = 2*m_nside;

  if (pix < m_ncap) {
    iring = (1 + isqrt(1 + 2*pix)) >> 1;
    iphi = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    face = (int)special_div(iphi-1, nr);
  }
  else if (pix < (m_npix-m_ncap)) {
    const I ip = pix - m_ncap;
    const I tmp = (m_order >= 0) ? ip >> (m_order+2) : ip/(4*m_nside);
    iring = tmp + m_nside;
    iphi = ip - tmp*4*m_nside + 1;
    kshift = (iring + m_nside)&1;
    nr = m_nside;
    const I ire = tmp + 1,
            irm = nl2 + 1 - tmp;
    I ifm = iphi - (ire>>1) + m_nside - 1,
      ifp = iphi - (irm>>1) + m_nside - 1;
    if (m_order >= 0)
      { ifm >>= m_order; ifp >>= m_order; }
    else
      { ifm /= m_nside; ifp /= m_nside; }
    face = (int)((ifp == ifm) ? (ifp|4) : ((ifp < ifm) ? ifp : (ifm+8)));
  }
  else {
    const I ip = m_npix - pix;
    iring = (1 + isqrt(2*ip - 1)) >> 1;
    iphi = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2 - iring;
    face = (int)special_div(iphi-1, nr) + 8;
  }

  const I irt = iring - ((2 + (face>>2))*m_nside) + 1;
  I ipt = 2*iphi - kJpll[face]*nr - kshift - 1;
  if (ipt >= nl2) ipt -= 8*m_nside;

  ix = (int)(( ipt-irt) >> 1);
  iy = (int)((-ipt-irt) >> 1);
}


// ============================================================================


void otswap::internal::HealpixBase::ring_info_small (const I ring, I& startpix, I& ringpix,
                                                     bool& shifted) const
{
  if (ring < m_nside) {
    shifted = true;
    ringpix = 4*ring;
    startpix = 2*ring*(ring-1);
  }
  else if (ring < 3*m_nside) {
    shifted = ((ring-m_nside) & 1) == 0;
    ringpix = 4*m_nside;
    startpix = m_ncap + (ring-m_nside)*ringpix;
  }
  else {
    shifted = true;
    const I nr = 4*m_nside - ring;
    ringpix = 4*nr;
    startpix = m_npix - 2*nr*(nr+1);
  }
}


// ============================================================================


I otswap::internal::vec2pix (const HealpixBase& base, const Vec3& v)
{
  const double xl = 1./v.length();
  const double phi = (v.x == 0. && v.y == 0.) ? 0.0 : det_atan2(v.y, v.x);
  const double nz = v.z*xl;
  if (std::abs(nz) > 0.99)
    return base.loc2pix(nz, phi, std::sqrt(v.x*v.x + v.y*v.y)*xl, true);
  else
    return base.loc2pix(nz, phi, 0, false);
}


// ============================================================================


I otswap::internal::ang2pix (const HealpixBase& base, const double theta, const double phi)
{
  return ((theta < 0.01) || (theta > 3.14159-0.01)) ?
    base.loc2pix(det_cos(theta), phi, det_sin(theta), true) :
    base.loc2pix(det_cos(theta), phi, 0., false);
}
