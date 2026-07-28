/*
  Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
*/
/**
 * @file RecCaloCommon/tests/phihelper_test.cxx
 *
 * Adapted from CxxUtils/tests/phihelper_test.h from ATLAS:
 * @file CxxUtils/test/phihelper_test.cxx
 * @author Frank Winklmeier
 * @date May 2019
 * @brief Regression tests for phihelper
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE phihelper_test
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include "RecCaloCommon/phihelper.h"
#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

namespace utf = boost::unit_test;
using k4::recCalo::deltaPhi;
using k4::recCalo::phiBisect;
using k4::recCalo::phiMean;
using k4::recCalo::wrapToPi;

// Types we want to test
typedef boost::mpl::list<float, double, long double> test_types;

// Tolerance for floating point comparisons (float, double, long double)
#define TOLERANCE *utf::tolerance(1e-4f) * utf::tolerance(1e-8) * utf::tolerance(1e-8L)

// cppcheck-suppress unknownMacro
BOOST_TEST_DECORATOR(TOLERANCE)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_wrap, T, test_types) {
  // Needed for boost to apply the correct tolerance in the comparison.
  // Also note the use of floats for all literals (0.0f) when used
  // in a comparison to ensure the lowest common type is float.
  constexpr T PI = std::numbers::pi_v<T>;

  // No wrapping for values within |pi|
  BOOST_TEST(wrapToPi<T>(3.0) == 3.0f);
  BOOST_TEST(wrapToPi<T>(-3.0) == -3.0f);
  // Corner cases
  BOOST_TEST(wrapToPi<T>(0.0) == 0.0f);
  BOOST_TEST(wrapToPi<T>(PI) == PI);
  BOOST_TEST(wrapToPi<T>(-PI) == -PI);
  BOOST_TEST(wrapToPi<T>(2 * PI) == 0.0f);
  BOOST_TEST(wrapToPi<T>(-2 * PI) == 0.0f);
  if constexpr (sizeof(T) <= sizeof(double)) {
    // With long double rounding, we have PI+PI+PI-PI-PI > PI, so this
    // test doesn't work out in that case.
    BOOST_TEST(wrapToPi<T>(3 * PI) == PI);
    BOOST_TEST(wrapToPi<T>(-3 * PI) == -PI);
  }
  // Wrap once
  BOOST_TEST(wrapToPi<T>(4.0) == 4.0f - 2 * PI);
  BOOST_TEST(wrapToPi<T>(-4.0) == -4.0f + 2 * PI);
  // Wrap many times
  BOOST_TEST(wrapToPi<T>(10 * PI) == 0.0f);
  BOOST_TEST(wrapToPi<T>(100 * PI) == 0.0f);
  // Check consistency between loop and remainder solution
  BOOST_TEST(wrapToPi<T>(40.1 * PI) == wrapToPi<T>(10.1 * PI));
  BOOST_TEST(wrapToPi<T>(-40.1 * PI) == wrapToPi<T>(-10.1 * PI));
}

// cppcheck-suppress unknownMacro
BOOST_TEST_DECORATOR(TOLERANCE)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_delta, T, test_types) {
  constexpr T PI = std::numbers::pi_v<T>;
  BOOST_TEST(deltaPhi<T>(3.0, 2.0) == 1.0f);
  BOOST_TEST(deltaPhi<T>(3 * M_PI, 2 * PI) == PI);
  BOOST_TEST(deltaPhi<T>(30 * M_PI, 25 * PI) == PI);
}

BOOST_TEST_DECORATOR(TOLERANCE)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_mean, T, test_types) {
  constexpr T PI = std::numbers::pi_v<T>;
  // Check values against unit vector addition
  std::vector<std::pair<T, T>> v = {{0, 0},
                                    {-1, 1},
                                    {-1, -1},
                                    {1, 2},
                                    {5, 8},
                                    {0.25 * PI, 0.75 * PI},
                                    {-0.25 * PI, -0.75 * PI},
                                    {-0.25 * PI, -0.75 * PI},
                                    {9 * PI / 4, 11 * PI / 4}};
  for (const auto& [a, b] : v) {
    // Check symmetry
    const T r1 = phiMean<T>(a, b);
    const T r2 = phiMean<T>(b, a);
    BOOST_TEST(r1 == r2, "phiMean(" << a << "," << b << ") not symmetric: " << r1 << " != " << r2);

    // Compare to trigonometric result
    const T mean = atan2(sin(a) + sin(b), cos(a) + cos(b));
    BOOST_TEST(phiMean<T>(a, b) == mean);
  }

  // If mean is |pi| we can only check the absolute value
  BOOST_TEST(std::abs(phiMean<T>(0.75 * PI, -0.75 * PI)) == PI);

  // Check values that result in zero unit vector, i.e. atan2(0,0)
  // Our algorithm always returns values within [-pi/2, pi/2]
  std::vector<std::tuple<T, T, T>> w = {{0.0, 1.0 * PI, 0.5 * PI},
                                        {0.25 * PI, -0.75 * PI, -0.25 * PI},
                                        {0.5 * PI, -0.5 * PI, 0.0},
                                        {0.75 * PI, -0.25 * PI, 0.25 * PI},
                                        {0.0, -1.0 * PI, -0.5 * PI}};
  for (const auto& [a, b, r] : w) {
    // Check symmetry and result
    const auto result = phiMean<T>(a, b);
    BOOST_TEST(result == phiMean<T>(b, a));
    BOOST_TEST(result == r, "phiMean(" << a << "," << b << ") = " << result << " != " << r);
  }
}

BOOST_TEST_DECORATOR(TOLERANCE)
BOOST_AUTO_TEST_CASE_TEMPLATE(test_bisect, T, test_types) {
  constexpr T PI = std::numbers::pi_v<T>;
  BOOST_TEST(phiBisect<T>(-1.0, 1.0) == 0.0f);
  BOOST_TEST(phiBisect<T>(1.0, -1.0) == PI);
  BOOST_TEST(phiBisect<T>(-1.0, -1.0) == -1.0f);
  BOOST_TEST(phiBisect<T>(0.0, PI) == 0.5f * PI);
  BOOST_TEST(phiBisect<T>(PI, 0.0) == -0.5f * PI);
  BOOST_TEST(phiBisect<T>(0.25 * PI, 0.75 * PI) == 0.5f * PI);
  BOOST_TEST(phiBisect<T>(-0.25 * PI, -0.75 * PI) == 0.5f * PI);
  BOOST_TEST(phiBisect<T>(0.75 * PI, -0.75 * PI) == PI);
  BOOST_TEST(phiBisect<T>(-0.75 * PI, 0.75 * PI) == 0.0f);
  BOOST_TEST(phiBisect<T>(9 * PI / 4, 11 * PI / 4) == 0.5f * PI);
}
