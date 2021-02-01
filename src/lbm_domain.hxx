#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstring> // std::strncpy

#include "status.hxx" // status_t

class Domain {
  // rectangular 3D domain of lattice points
private:
  // static members
  int static constexpr id_Length = 7;
  // members
  int32_t n_[3]; // number of grid points
  int32_t offset_[3];
  char id_[id_Length + 1];

public:
  
  Domain(int32_t const extent=0, char const *id="?^3") {
      int32_t const extent_cubic[] = {extent, extent, extent};
      constructor(extent_cubic, id);
  } // [default] constructor

  Domain(int32_t const extent[3], char const *id="???") {
      constructor(extent, id);
  } // constructor

  // number of all grid points
  size_t n() const { return size_t(n_[0])*size_t(n_[1])*size_t(n_[2]); }

  int32_t n(int d) const {
      assert(0 <= d); assert(d < 3);
      return n_[d];
  } // n

  int32_t operator [](char direction) const {
      return n((direction | 32) - 'x');
  } // []
  
  int32_t operator [](int d) const {
      assert(0 <= d); assert(d < 3);
      return n_[d];
  } // []
  
  int32_t m(int d) const {
      assert(0 <= d); assert(d < 3);
      return n_[d] - 1;
  } // m == n - 1

  int32_t const * global_coordinate_offset() const { return offset_; }
  int32_t const * extent() const { return n_; }

  // local coordinates of the neighboring +/-1 (or more) lattice points
  template <char Direction>
  inline int32_t neighbor(int32_t coordinate, int shift) {
      int const dir = (Direction | 32) - 'x'; // ignore case with |32
      int32_t const n = n_[dir];
      assert(0 <= coordinate); assert(coordinate < n);
      int32_t const sc = coordinate + shift;
      assert(shift*shift <= n*n); // the following return fomula is incorrect for |shift| > n
      return sc + (sc < 0)*n - (sc >= n)*n; // modulo operation
  } // neighbor

  inline int32_t neighbor(int32_t coordinate, int shift, int dir) {
      switch (dir) {
          case  0: return neighbor<'x'>(coordinate, shift);
          case  1: return neighbor<'y'>(coordinate, shift);
          case  2: return neighbor<'z'>(coordinate, shift);
          default: assert(0 <= dir); assert(dir < 3); return 0;
      } // dir
  } // neighbor

private:

  void constructor(int32_t const extent[3], char const *id) {

      int32_t constexpr MaxExtent = (1l << 21); // wrap round for global_coordinates
      size_t product{1};
      for(int d = 0; d < 3; ++d) {
          auto const ext = (extent) ? extent[d] : MaxExtent;
          assert(0 <= ext); assert(ext <= MaxExtent);
          n_[d] = ext;
          product *= ext;
          offset_[d] = 0;
      } // d

      if (product > 0) {
          std::strncpy(id_, id, id_Length); // copy the identifyer
          id_[id_Length] = '\0'; // ensure null-termination

          std::printf("# Domain: %d x %d x %d = %.3f k lattice sites, id=%s\n",
                          n_[0], n_[1], n_[2],  product*.001,         id_);
      } else {
          id_[0] = '_'; id_[1] = '\0'; // mark that there were no lattice sites
      } // default constructor is silent

  } // constructor

}; // class Domain


namespace lbm_domain {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_constructors(int const echo=0) {
      int32_t const n111[] = {1, 1, 1};
      Domain d000default, d000(0), dmin(n111, "1x1x1_long_name"); // envoke different constructors
      return 0;
  } // test_constructors

  inline status_t test_domain(int const echo=0) {
      status_t stat(0);
      int32_t const nnn[] = {3, 5, 7};
      Domain domain(nnn, "3x5x7");
      for(int d = 0; d < 3; ++d) { // spatial direction
          auto const n = domain[d];
          stat += (n != domain.n(d));
          stat += (n != domain[char('x' + d)]);
          stat += (n != domain[char('X' + d)]);
          for(int i = 0; i < n; ++i) { // local coordinate
              stat += (domain.neighbor(i,  0, d) !=  i            );
              stat += (domain.neighbor(i, +1, d) != (i +1    ) % n);
              stat += (domain.neighbor(i, -1, d) != (i -1 + n) % n);
          } // i
      } // d direction
      return stat;
  } // test_domain

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_constructors(echo);
      stat += test_domain(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_domain
