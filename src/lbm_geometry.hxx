#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // uint32_t, size_t

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
// (3*21) bit interleaved coordinates
#include "global_coordinates.hxx" // ::get

class Geometry {
  // global cell description
private:
  // members
  uint32_t periodic_[3];

public:
  
  Geometry(size_t const period[3]=nullptr) {

      uint32_t constexpr MaxPeriod = (1ul << 21); // wrap round for global_coordinates
      size_t product{1};
      for(int d = 0; d < 3; ++d) {
          size_t const peri = (period) ? period[d] : MaxPeriod;
          assert(0 < peri); assert(peri <= MaxPeriod);
          periodic_[d] = peri;
          product *= peri;
      } // d
      std::printf("# %s: max %.3f k x %.3f k x %.3f k = %.6f M lattice sites\n", __func__
                , periodic_[0]*.001, periodic_[1]*.001, periodic_[2]*.001, product*1e-6);

  } // [default] constructor

}; // class Geometry


namespace lbm_geometry {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_geometry(int const echo=0) {
      status_t stat(0);
      size_t const p111[] = {1, 1, 1};
      Geometry gmax, gmin(p111); // envoke constructor for largest and smallest possible cell
      return stat;
  } // test_geometry

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_geometry(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_geometry
