#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf
// #include <cmath> // std::tanh

//   template <typename real_t>
//   inline real_t pow2(real_t const v) { return v*v; } // ToDo: move to utils

  class Droplet {
  public:

      Droplet(
            int const x=0, int const y=0, int const z=0
          , int const radius=-999
          , int const drop=0
          , double const rho_low=0
          , double const rho_high=1
          , int const echo=0 // increase echo if also the default constructor should print to stdout
      )
          : rho_low_(rho_low), rho_high_(rho_high), drop_(drop)
      {
          xyzr_[0] = x; xyzr_[1] = y; xyzr_[2] = z; xyzr_[3] = radius;
          if (echo > 0) {
              std::printf("# %s at %d %d %d with radius %d, drop=%d\n", __func__,
                    xyzr_[0], xyzr_[1], xyzr_[2], xyzr_[3], drop_);
          } // drop==0 (default constructor should not print)
          assert(drop*drop <= 1 && "meaning of |drop| > 1 unknown!");
      } // constructor and default constructor

      int const* position() const { return xyzr_; }
      int position(int d) const { assert(0 <= d); assert(d < 4); return xyzr_[d]; }
      int radius() const { return xyzr_[3]; }
      int drop() const { return drop_; }

      template <typename real_t>
      double density(
            real_t const x
          , real_t const y
          , real_t const z
          , double const inverse_interface_width=1 // 1/ifaceW
      ) {
          auto const p = position();
          double const dist2 = pow2(x - p[0]) + pow2(y - p[1]) + pow2(z - p[2]);
          auto const dist = std::sqrt(dist2);
          auto const arg = 2*(dist - radius())*inverse_interface_width;
          return 0.5*((rho_high_ + rho_low_) - drop_*(rho_high_ - rho_low_)*std::tanh(arg));
      } // density
      
  private:
      // members
      double rho_low_, rho_high_;
      int xyzr_[4];
      int drop_; // 1:droplet, -1:bubble, 0:not_initialized

  }; // Droplet
