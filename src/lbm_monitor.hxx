#pragma once

#include <cstdio> // std::printf, std::sprintf, std::snprintf
#include <cassert> // assert
#include <cstring> // std::strncpy
#include <cmath> // std::round
#include <algorithm> // std::min, std::max

#include "status.hxx" // status_t

#include "lbm_domain.hxx" // Domain

namespace lbm_monitor {

  char constexpr esc = 27; // escape char
  char const def[4] = {esc,'[','m','\0'};

// ! In 256 color mode (ESC[38;5;<fgcode>m and ESC[48;5;<bgcode>m), the color-codes are the following:[citation needed]
// ! 
// !  0x00-0x07:  standard colors (as in ESC [ 30..37 m)
// !  0x08-0x0f:  high intensity colors (as in ESC [ 90..97 m)
// !  0x10-0xe7:  6*6*6=216 colors: 16 + 36*r + 6*g + b (0<=r,g,b<=5)
// !  0xe8-0xff:  grayscale from black to white in 24 steps
//   
// 
// ! \e[31m --> red font color
// ! \e[31;1m --> red bold font color
// ! \e[41m --> red background color
// 
// ! Intensity       0       1       2       3       4       5       6       7
// ! Normal          Black   Red     Green   Yellow  Blue    Magenta Cyan    White
// ! Bright          Black   Red     Green   Yellow  Blue    Magenta Cyan    White
// 

  template<unsigned Nchars=8>
  class string_t {
  private:
      char c[Nchars];
  public:
      string_t(char const *string=nullptr) { 
          if (Nchars) c[0] = '\0'; 
          if (nullptr != string) {
              std::snprintf(c, Nchars, "%s", string);
          }
      } // default constructor
      char const* c_str() const { return (Nchars) ? c : nullptr; }
      char*        data()       { return (Nchars) ? c : nullptr; }
      
      bool operator==(string_t const & rhs) const {
          bool same{true};
          for(int k = 0; k < Nchars; ++k) {
              same &= (c[k] == rhs.c[k]);
          } // k
          return same;
      } // ==

  }; // class string_t

  
  typedef string_t<8>  string7_t;
  typedef string_t<16> string15_t;

  template <int N=5>
  inline int scala(float const value) { return std::min(std::max(0, int(std::round(value*N))), N); }

  inline string15_t colorchar(int const rd6, int const gr6, int const bl6) { // rgb < 6
      int const i = 16 + (bl6*6 + gr6)*6 + rd6;
      string15_t s; std::sprintf(s.data(), "%c[48;5;%im", esc, i);
      return s;
  } // colorchar

  inline string15_t colorchar(uint8_t const rgb[3]) { return colorchar(rgb[0], rgb[1], rgb[2]); }
  
//   inline double ssqrt(double const x) { return x; }
  // signed squareroot
  inline double ssqrt(double const x) { return ((x < 0) ? -1: 1)*std::sqrt(std::abs(x)); }

  union RGBA32 {
      int32_t i32; // makes it possible to easily compare if twpo colors are the same
      uint8_t rgb[4]; // red, green, blue and alpha (opacity) in [0, 255]
  }; // RGBA32

  template <typename real_t>
  status_t show_vector_field(
        real_t const ux[]
      , real_t const uy[]
      , real_t const uz[]
      , Domain const & domain
      , int const z=0
      , double const *const maxval=nullptr
      , double const *const minval=nullptr
  ) {

      int const Nx = domain['x'], Ny = domain['y'], Nz = domain['z'];
      assert(0 <= z); assert(z < Nz);
      // determine minimum and maximum in the data set
      double mini{9e33}, maxi{-9e33}; // init min and max
      for(int xy = 0; xy < Ny*Nx; ++xy) {
          int const xyz = z*Ny*Nx + xy;
          double const value = pow2(ux[xyz]) + pow2(uy[xyz]) + pow2(uz[xyz]);
          mini = std::min(mini, value);
          maxi = std::max(maxi, value);
      } // xy

      auto const f = 1/std::sqrt(std::max(1e-12, maxi));
      std::printf("# %s values in [%g %g] --> sqrt([%g %g]) on %d x %d\n", 
                    __func__, mini, maxi, mini*f*f, maxi*f*f,   Ny, Nx);

      for(int y = 0; y < Ny; ++y) {
          RGBA32 prev{1 << 24}; // initialize with a value outside of 24bit
          for(int x = 0; x < Nx; ++x) {
              int const xyz = (z*Ny + y)*Nx + x;
              RGBA32 c{0};
              c.rgb[0] = scala(.5 + .5*ssqrt(ux[xyz]*f));
              c.rgb[1] = scala(.5 + .5*ssqrt(uy[xyz]*f));
              c.rgb[2] = scala(.5 + .5*ssqrt(uz[xyz]*f));
              // print "  ", i.e. two blank per cell
              if (c.i32 == prev.i32) {
              } else {
                  auto const color = colorchar(c.rgb);
                  std::printf("%s", color.c_str());
                  prev = c;
              }
              std::printf("  ");
          } // x
          std::printf("%s\n", def); // reset and newline 
      } // y
      std::printf("%s\n", def); // reset and newline 
      return 0;
  } // show_vector_field
  
  template <typename real_t>
  status_t show_scalar(
        real_t const rho[]
      , int const Nx
      , int const Ny
      , int const Nz=1
      , int const z=0
      , double const *const maxval=nullptr
      , double const *const minval=nullptr
  ) {

      // determine minimum and maximum in the data set
      double mini{9e33}, maxi{-9e33}; // init min and max
      for(int xy = 0; xy < Ny*Nx; ++xy) {
          int const xyz = z*Ny*Nx + xy;
          double const value = rho[xyz];
          mini = std::min(mini, value);
          maxi = std::max(maxi, value);
      } // xy

      auto const f = 1/std::max(1e-12, maxi - mini);
      std::printf("# %s values in [%g %g] --> [%g %g] on %d x %d\n", 
                  __func__, mini, maxi, f, (maxi - mini)*f, Ny, Nx);

      for(int y = 0; y < Ny; ++y) {
          int prev{-1};
          for(int x = 0; x < Nx; ++x) {
              int const xyz = (z*Ny + y)*Nx + x;
              float const grey = (rho[xyz] - mini)*f;
              auto const igrey = scala(grey);
              // print "  ", i.e. two blank per cell
              if (igrey == prev) {
                  // we do not have to change the color
                  std::printf("  ");
              } else {
                  auto const color = colorchar(igrey, igrey, igrey);
                  std::printf("%s  ", color.c_str());
                  prev = igrey;
              }
          } // x
          std::printf("%s\n", def); // reset and newline 
      } // y
      std::printf("%s\n", def); // reset and newline 
      return 0;
  } // show_scalar
  


#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_simple_colors(int const echo=0) {

      // the 6x6x6 color cube
      for (int ib = 5; ib >= 0; --ib) {
          for (int ig = 5; ig >= 0; --ig) {
              for (int ir = 5; ir >= 0; --ir) {
                  // print a single ' ' with the color formatter in front
                  std::printf("%s ", colorchar(ir/5., ig/5., ib/5.).c_str());
              } // ir
          } // ig
          std::printf("%s\n", def); // reset and newline 
      } // ib
      std::printf("%s\n", def); // reset and newline 

      return 0;
  } // test_simple_colors

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_simple_colors(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_monitor
