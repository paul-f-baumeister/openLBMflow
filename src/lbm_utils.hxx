#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf
#include <cmath> // std::atan
#include <vector> // std::vector<T>

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

  template <typename real_t> inline
  real_t pow2(real_t const x) { return x*x; }

  template <typename real_t> inline
  real_t pow3(real_t const x) { return x*x*x; }

  template <typename real_t> inline
  real_t pow4(real_t const x) { return pow2(pow2(x)); }
  
namespace lbm_utils {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += STATUS_TEST_NOT_INCLUDED;
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_utils
