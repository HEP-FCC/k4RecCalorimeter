/*
  Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file RecCaloCommon/include/RecCaloCommon/phihelper.h
 *
 * Adapted from CxxUtils/phihelper.h from ATLAS:
 * @author Frank Winklmeier
 * @date May 2019
 * @brief Helper for azimuthal angle calculations
 */
#ifndef RECCALOCOMMON_PHIHELPER_H
#define RECCALOCOMMON_PHIHELPER_H

#include <cmath>
#include <numbers>
#include <type_traits>

namespace k4::recCalo {

/**
 * Wrap angle in radians to [-pi, pi]
 *
 * Odd positive (negative) multiples of pi map to (-)pi.
 */
template <std::floating_point T>
inline T wrapToPi(T phi) {
  constexpr T PI = std::numbers::pi_v<T>;
  // For large values this is faster:
  if (phi < -100 || phi > 100) {
    return std::remainder(phi, 2 * PI);
  }
  while (phi > PI)
    phi -= 2 * PI;
  while (phi < -PI)
    phi += 2 * PI;
  return phi;
}

/**
 * Return difference phiA - phiB in range [-pi, pi]
 */
template <std::floating_point T>
inline T deltaPhi(T phiA, T phiB) {
  return wrapToPi(phiA - phiB);
}

/**
 * Calculate average of two angles
 *
 * The average is calculated as the angle of the sum of the two unit vectors
 * with angles phiA and phiB. As such it always returns the angle in the
 * narrower of the two complimentary regions spanned by phiA and phiB. If the
 * vector sum is zero, return the angle within [-pi/2, pi/2]. This method is
 * symmetric in its arguments except if |mean| equals pi.
 *
 * The returned value is within the range [-pi, pi].
 */
template <std::floating_point T>
inline T phiMean(T phiA, T phiB) {
  const T diff = wrapToPi(phiA - phiB);
  return wrapToPi(phiB + 0.5 * diff);
}

/**
 * Bisect (average) the angle spanned by phiA and phiB
 *
 * The average is calculated by rotating phiA half-way counter-clockwise onto
 * phiB. Note that his method is not symmetric in its arguments. Typical
 * use-cases are to determine the centre of a region in the detector.
 *
 * The returned value is within the range [-pi, pi].
 */
template <std::floating_point T>
inline T phiBisect(T phiA, T phiB) {
  T phi = 0.5 * (phiA + phiB);
  if (phiA > phiB)
    phi += std::numbers::pi;
  return wrapToPi(phi);
}

} // namespace k4::recCalo

#endif // not RECCALOCOMMON_PHIHELPER_H
