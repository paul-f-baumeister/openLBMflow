#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf
#include <cmath> // std::atan
#include <vector> // std::vector<T>

#include "data_view.hxx" // view4D<T>
#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#include "lbm_utils.hxx" // pow2
#include "lbm_monitor.hxx" // ::show_scalar
  
  double factorial(int const n) { return (n < 2) ? (n >= 0) : n*factorial(n - 1); }

  template <typename real_t>
  inline void Hermite_polynomials(
        real_t H[]
      , double const x
      , int const numax
      , double const sigma2=1 // sigma^2
  ) {
      double Hnum1{0}, Hnu{1}, Hnup1;
      for(int nu = 0; nu <= numax; ++nu) {
          H[nu] = Hnu; // store
          // use the two step recurrence relation to create the 1D Hermite polynomials
          // H[nu+1] = x * H[nu] - nu * H[nu-1];
          Hnup1 = x * Hnu - nu * Hnum1 * sigma2;
          // rotate registers for the next iteration
          Hnum1 = Hnu; Hnu = Hnup1; // ordering is important here
      } // nu
  } // Hermite_polynomials

  template <typename real_t>
  inline void Hermite_polynomials_saving( // 4*numax flop
        real_t H[]
      , double const x
      , int const numax
      , double const sigma2=1 // sigma^2
  ) {
      // same as Hermite_polynomials with some less Flops
      H[0] = 1;
      double Hnum1, Hnu{0}, Hnup1{1}; // init for nu=-1
      for(int nu = 0; nu < numax; ++nu) {
          // rotate registers
          Hnum1 = Hnu; Hnu = Hnup1; // ordering is important here
          // use the two step recurrence relation to create the 1D Hermite polynomials
          // H[nu+1] = x * H[nu] - nu * H[nu-1];
          Hnup1 = x * Hnu - nu * Hnum1 * sigma2;
          H[nu + 1] = Hnup1;
      } // nu
  } // Hermite_polynomials_saving

  inline int test_Hermite_orthogonality(int const echo=5, int const numax=7) {
      double const dx = 0.25;
      double const sigma = 2;
      int const mx = 19.0*sigma/dx;
      int const M = 1 + numax;
      double const pi = 4*std::atan(1.0);
      double const sqrt2pi = std::sqrt(2*pi);
      std::vector<double> ovl(M*M, 0);
      double H[96];
      for(int ix = -mx; ix <= mx; ++ix) {
          double const xi = ix*dx/sigma;
          Hermite_polynomials(H, xi, numax);
          double const Gauss = std::exp(-.5*xi*xi)*dx;
          for(int i = 0; i <= numax; ++i) {
              for(int j = 0; j <= numax; ++j) {
                  ovl[i*M + j] += H[i]*H[j]*Gauss;
              } // j
          } // i
      } // ix

      double dev[] = {0, 0};
      for(int i = 0; i <= numax; ++i) {
          if (echo > 5) std::printf("# %s: row%3i ", __func__, i);
          double const f = 1./(sigma*factorial(i)*sqrt2pi);
          for(int j = 0; j <= numax; ++j) {
              double const ovl_ij = ovl[i*M + j]*f;
              if (echo > 5) std::printf("%16.9f", ovl_ij);
              int const diag = (i == j);
              dev[diag] += std::abs(ovl_ij - diag);
          } // j
          if (echo > 5) std::printf("\n");
      } // i

      if (echo > 2) std::printf("# %s: deviations %.1e and %.1e for diagonal and"
          " off-diagonal elements\n", __func__,  dev[1], dev[0]);
      return (dev[0] > 1e-14) + (dev[1] > 1e-14);
  } // test_Hermite_orthogonality


  inline int number_of_moments(int const maxM, int const D, int const Q, int const echo=0) {
      int nom{1}, denom{1};
      for(int d = 0; d < D; ++d) {
          nom   *= (1 + d + maxM);
          denom *= (1 + d);
      } // d
      assert(denom > 0);
      assert(0 == nom % denom);
      int const n_momenta = nom/denom;
      if (echo > 7) std::printf("# D%dQ%d   %d moments for order %d \n", D, Q, n_momenta, maxM);
      return n_momenta;
  } // number_of_moments
           
           
  inline int find_maximum_moments(int const D, int const Q, int const echo=0) {
      // find the maximum number of moments
      int maxM{0};
      while (number_of_moments(maxM + 1, D, Q, echo) <= Q) { ++maxM; }
      
      if (echo > 2) std::printf("# D%dQ%d  %d moments up to order %d \n",
                      D, Q, number_of_moments(maxM, D, Q, echo), maxM);
      return maxM;
  } // find_maximum_moments


  
    template <typename real_t>
    inline double check_overlap(
          real_t const projector[]  // [Q][M + 1] // Hermite polynomials + kinetic energy in the last column
        , real_t const expansion[]  // [M][Q]     // weighted transposed Hermite polynomials
        , int const Q               // number of populations
        , int const M               // number of moments
        , int const echo=0          // log-level
    ) {
        double dev{0}; // deviation from a unit matrix

        for(int n = 0; n < M; ++n) {
            for(int m = 0; m < M; ++m) {
                double t{0};
                for(int q = 0; q < Q; ++q) {
                    t += projector[q*(M + 1) + n] * expansion[m*Q + q];
                } // m
                auto const dt = std::abs(t - int(n == m));
                if (echo > 7 || dt > 1e-12) std::printf("# %s: m=%d m\'=%d %g\n", __func__, n, m, t);
                dev = std::max(dev, dt);
            } // m
        } // n

        return dev;
    } // check_overlap
  
  
#ifdef  STANDALONE_TEST
    #define restrict
#endif

    template <typename real_t>
    inline double collide_general(
          real_t       *restrict const f_out    // result populations[Q] (optional)
        , real_t       *restrict const moment   // result moments, optional input
        , real_t const *restrict const f_in     // input  populations[Q] (optional, take moments otherwise)
        , double const projector[]  // [Q][M + 1] // Hermite polynomials + kinetic energy in the last column
        , double const expansion[]  // [M][Q]     // weighted transposed Hermite polynomials
        , int const D               // dimensions
        , int const Q               // number of populations
        , int const M               // number of moments
        , int const numax           // maximum order (should match M == sum_nu^numax number_of_moments(nu))
        , double const omega[]      // relaxation parameters (omega=1/tau)
//      , double const body_force_xyz[3]=nullptr
        , int const echo=0          // log-level
        , double const vel[][4]=nullptr
    ) {
        if (nullptr != f_in) {
            // transform input populations into momenta
            for(int m = 0; m <= M; ++m) moment[m] = 0;
            for(int q = 0; q < Q; ++q) {
                for(int m = 0; m <= M; ++m) {
                    moment[m] += f_in[q] * projector[q*(M + 1) + m]; // 1184 Flop for D2Q37
                } // m
            } // q
        } // f_in

        if (echo > 11) {
            std::printf("\n# %s: %s %d moments\n", __func__, f_in?"projectored":"input", M);
            for(int m = 0; m < M; ++m) {
                std::printf("# found moment[%i]= %g\n", m, moment[m]);
            } // m
            std::printf("# rho*kinetic_energy = %g\n\n", moment[M]);
        } // echo
        
        real_t u2{0};
        auto const rho = moment[0];
        if (rho > 0) {
            auto const inv_rho = 1/rho; // 1 inversion
            auto const kinetic_energy = moment[M]*inv_rho; // 1 flop
            if (echo > 7) std::printf("# kinetic_energy = %g\n", kinetic_energy/D);
            real_t u[3];
            for(int d = 0; d < D; ++d) {
                u[d] = moment[1 + d]*inv_rho; // velocity, D flop
                u2 += pow2(u[d]); // 2*D flop
            } // d
            

    //         if (nullptr != body_force_xyz) {
    //             add the body force
    //             for(int d = 0; d < D; ++d) {
    //                 mom[1 + d] += tau[1 + d]*body_force_xyz[d]*inv_rho;
    //             } // d
    //         } // body_force?

            // the temperature T is defined as 1/D*( <v^2> - <v>^2 ), the variance of the distribution
            auto const temperature = kinetic_energy - u2*(1./D); // 2 flop
            if (echo > 7) std::printf("# temperature = %g\n", temperature);

            

            // prepare equilibrium Hermite polynomials
            real_t Hermite_eq[3][8];
            Hermite_eq[2][0] = 1; // for the 2D case
            for(int d = 0; d < D; ++d) {
                Hermite_polynomials(Hermite_eq[d], u[d], 7, 1 - temperature);
            } // d

#define CHECK_ME
#ifdef  CHECK_ME
            double meq_ref[35];
#endif // CHECK_ME

            { // scope: relax towards equilibrium momenta
                int m{0};
                for(int nu = 0; nu <= numax; ++nu) { // order
                    auto const omega_nu = omega[nu];
    //                 if (nu < 2) assert(0 == omega_nu); // we must not relax the momenta of 0th and 1st order for conservation laws
                    auto const one_minus_omega_nu = 1 - omega_nu;
                    for(int nz = 0; nz <= nu*(D > 2); ++nz) { // nz runs from 0 to 0 if D==2
                        for(int ny = 0; ny <= nu - nz; ++ny) {
                            int const nx = nu - ny - nz;
                            auto const equilibrium_moment = rho*Hermite_eq[0][nx]
                                                               *Hermite_eq[1][ny]
                                                               *Hermite_eq[2][nz];
#ifdef  CHECK_ME
                            meq_ref[m] = equilibrium_moment; // store
#endif // CHECK_ME
                            auto const old_moment = moment[m];
                            moment[m] = one_minus_omega_nu*old_moment + omega_nu*equilibrium_moment;
                            if (echo > 9) {
                                std::printf("# %d (%1x%1x%1x) moments:  old %.6f\tequ %.6f\tnew %.6f\n", 
                                               nu, nz,ny,nx, old_moment, equilibrium_moment, moment[m]);
                            } // echo
                            assert(m < M);
                            ++m;
                        } // ny
                    } // nz
                } // nu
                assert(M == m);
            } // scope
            
#ifdef  CHECK_ME
            if (D == 2 && Q == 37) {
                // compute the equilibrium populations the classical way
                // according to https://www.imacm.uni-wuppertal.de/fileadmin/imacm/preprints/2020/imacm_20_60.pdf
                if (nullptr != vel) {
                    double const Tm1 = temperature - 1;
//                  double const inv_cs2 = 1;
                    double const inv_cs2 = pow2(1.19697977039307435897239);
                    double const u4 = pow2(u2);
                    double feq[37];
                    for(int q = 0; q < Q; ++q) {
                        double const uc = u[0]*vel[q][0]  + u[1]*vel[q][1];
                        double const c2 = pow2(vel[q][0]) + pow2(vel[q][1]);
                        double const uc2 = pow2(uc);
                        
                        double const f0 = 1; // zeroth order
                        double const f1 = uc; // first order
                        double const f2 = inv_cs2/2.*( uc2 - u2 + Tm1*(c2 - D) ); // seconds order
                        double const f3 = pow2(inv_cs2)/6.*uc*( uc2 - 3*u2 + 3*Tm1*(c2 - D - 2) ); // third order
                        double const f4 = pow3(inv_cs2)/24.*( pow2(uc2) - 6*uc2*u2 + 3*u4
                                                            + 6*Tm1*( uc2*(c2 - D - 2) + u2*(2 + D - c2) )
                                                            + 3*pow2(Tm1)*( pow2(c2) - 2*(D + 2)*c2 + D*(D + 2) )
                                                            ); // fourth order
                        double const wq = expansion[q]; // pure weights
                        feq[q] = wq*rho*(f0 + f1 + f2 + f3 + f4);
                    } // q

                    // project them
                    double meq[15]; assert(15 == M);
                    if (echo > 1) std::printf("\n# compare %d moments: Gabbana vs. here\n", M);
                    for(int m = 0; m < M; ++m) {
                        meq[m] = 0;
                        for(int q = 0; q < Q; ++q) {
                            meq[m] += feq[q] * projector[q*(M + 1) + m];
                        } // q
                        if (echo > 1) {
                            double const dif = meq[m] - meq_ref[m];
                            double const rel = (std::abs(meq_ref[m]) > 1e-16) ? dif/meq_ref[m] : 0;
                            std::printf("# %dth moment: Gabbana %g here %g diff %g relative %g\n", 
                                        m, meq[m], meq_ref[m], dif, rel);
                        } // echo
                    } // m
                } // nullptr != vel
            } // D2Q37
#endif // CHECK_ME
            
            
            
        } else {
            for(int m = 0; m <= M; ++m) {
                moment[m] = 0; // clear all moments
            } // m
        } // rho > 0

        if (nullptr != f_out) {
            // transform momenta into output populations
            for(int q = 0; q < Q; ++q) f_out[q] = 0; // init results
            for(int m = 0; m < M; ++m) {
                for(int q = 0; q < Q; ++q) {
                    f_out[q] += moment[m] * expansion[m*Q + q]; // 1110 Flop for D2Q37
                } // q
            } // m
        } // f_out

        return u2;
    } // collide_general

    template <typename real_t>
    inline void propagate_general(
          view4D<real_t>       & f_out    // result populations[z][y][x][Q]
        , real_t const *restrict const f_in // input         populations[Q]
        , int const D               // dimensions
        , int const Q               // number of populations
        , int8_t const vel[][4]     // lattice-matched velocities
        , int const ix
        , int const iy
        , int const iz
        , int const nx
        , int const ny
        , int const nz=1
        , int const echo=0          // log-level
    ) {
        for(int q = 0; q < Q; ++q) {
            int const jx = (ix + vel[q][0] + 9*nx) % nx; 
            int const jy = (iy + vel[q][1] + 9*ny) % ny; 
            int const jz = (iz + vel[q][2] + 9*nz) % nz;
            f_out(jz,jy,jx,q) = f_in[q];
        } // q
    } // propagate_general

    
    template <typename real_t> 
    void initalize_Taylor_Green_2D(
          view4D<real_t> & moments
        , int const n[3]
        , int const iz=0
        , double const rho=1
    ) {
        double const pi = std::acos(-1.);
        double const TG2D[2][2] = {{1, -1}, {2*pi/n[0], 2*pi/n[1]}}; // Taylor-Green vortex flow
        assert( std::abs(TG2D[0][0]*TG2D[1][0] + TG2D[0][1]*TG2D[1][1]) < 1e-9 );
        for (int iy = 0; iy < n[1]; ++iy) {
            double const cos_y = std::cos(TG2D[1][1]*iy);
            double const sin_y = std::sin(TG2D[1][1]*iy);
            for (int ix = 0; ix < n[0]; ++ix) {
                moments(iz,iy,ix,0) = rho;
                moments(iz,iy,ix,1) = TG2D[0][0] * std::sin(TG2D[1][0]*ix) * cos_y;
                moments(iz,iy,ix,2) = TG2D[0][1] * std::cos(TG2D[1][0]*ix) * sin_y;
            } // ix
        } // iy
    } // initalize_Taylor_Green_2D

    
    template <typename real_t> 
    void initalize_droplet(
          view4D<real_t> & moments
        , int const n[3]
        , double const xyzr[4]
        , double const rho=1
    ) {
        auto const radius2 = pow2(xyzr[3]);
        for         (int iz = 0; iz < n[2]; ++iz) {   auto const z2 = pow2(iz - xyzr[2]);
            for     (int iy = 0; iy < n[1]; ++iy) {   auto const y2 = pow2(iy - xyzr[1]);
                for (int ix = 0; ix < n[0]; ++ix) {   auto const x2 = pow2(ix - xyzr[0]);
                    auto const dist2 = x2 + y2 + z2;
                    moments(iz,iy,ix,0) += (dist2 < radius2)*rho;
                } // ix
            }     // iy
        }         // iz
    } // initalize_droplet
    
namespace lbm_general {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t=double>
  inline status_t test_collide(int const D, int const echo=0, int const dist2=10) {
      status_t stat(0);
      if (echo > 0) std::printf("\n# %s start D=%d!\n", __func__, D);

      int constexpr v2 = 3;
      double  vel[999][4];
      int8_t ivel[999][4]; // integer velocities
      int     iv2[999];
      std::vector<double> weight_d2(99, 0.0);
      std::vector<int> multiplicity_d2(99, 0);

      int constexpr X=0, Y=1, Z=2;
      
      double const weights_D2Q37[11] = { // from Biferale et al., similar to Shan & Chen (2007)
        .23315066913235250228660,
        .10730609154221900241246,
        .05766785988879488203006,
        .0                      ,
        .01420821615845075026469,
        .00535304900051377523273,
        .0                      ,
        .0                      ,
        .00101193759267357547541,
        .00024530102775771734547,
        .00028341425299419821740 };
      double const r_D2Q37 = 1.19697977039307435897239, r = r_D2Q37;
      double const *w = weights_D2Q37;
      
      
      // create a stencil
      int qi{0};
      double weight_sum{0};
      for(int d2 = 0; d2 <= dist2; ++d2) {
          int const mx = std::sqrt(1.001*d2);
          int const my = mx*(D > 1);
          int const mz = mx*(D > 2);
          for(int z = -mz; z <= mz; ++z) {
          for(int y = -my; y <= my; ++y) {
          for(int x = -mx; x <= mx; ++x) {
              if (pow2(x) + pow2(y) + pow2(z) == d2) {
                  vel[qi][0] = x*r;
                  vel[qi][1] = y*r;
                  vel[qi][2] = z*r;
                  vel[qi][v2] = d2*r*r;
                  iv2[qi]     = d2; // integer distance^2
                  ivel[qi][0] = x;
                  ivel[qi][1] = y;
                  ivel[qi][2] = z;
                  ivel[qi][3] = d2;
                  ++qi;
                  ++multiplicity_d2[d2];
              } // v^2 matches d2
          } // x
          } // y
          } // z
          weight_d2[d2] = (nullptr != w) ? w[d2] : std::exp(-0.5*d2);
          weight_sum += weight_d2[d2] * multiplicity_d2[d2];
          std::printf("# D%d distance^2 %d multiplicity %d weight %g\n", D, d2, multiplicity_d2[d2], weight_d2[d2]);
      } // d2
      int const Q = qi;

      if (nullptr != w) {
          std::printf("# weights sum up to %.15f\n", weight_sum);
      } else {
          weight_sum = 1./weight_sum;
          for(int d2 = 0; d2 <= dist2; ++d2) weight_d2[d2] *= weight_sum; // normalize weights
      }
      

      double const weights_D3Q45[45] = {0.20740740740740618, 
        0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 0.05787037037037047, 
        0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 0.00462962962962958, 
        0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939, 0.0004629629629629939};
      // weights(|v|^2/3) = {.20740740740740618, .05787037037037047, 0, .00462962962962958, 0, .0004629629629629939}; // 0,3,9,15 appearing as v^2 with multiplicity 1,12,20,12
      double const velocity_D3Q45[3][45] = {
      {0.0, 0.06386083877343968, -0.06386083877343968, 1.2239121278243665, -1.2239121278243665, 1.2239121278243665, -1.2239121278243665, 1.5766994272507744, -1.5766994272507744, 0.5069610024977665, -0.5069610024977665, -0.5069610024977665, 0.5069610024977665, 2.403092127540177, -2.403092127540177, -2.403092127540177, 2.403092127540177, -0.8892242114059369, 0.8892242114059369, -0.8892242114059369, 0.8892242114059369, -1.5602655313772367, -1.5602655313772367, 1.5602655313772367, 1.5602655313772367, 0.4744978678080795, 0.4744978678080795, -0.4744978678080795, -0.4744978678080795, 2.9239876105912574, -2.9239876105912574, 1.7320508075688787, -1.7320508075688787, -2.7367507163016924, 2.7367507163016924, 2.7367507163016924, -2.7367507163016924, 0.14279717659756475, -0.14279717659756475, -3.5256070994177073, 3.5256070994177073, 1.1335992635264445, -1.1335992635264445, 1.1335992635264445, -1.1335992635264445},
      {0.0, -1.2239121278243665, 1.2239121278243665, -0.06386083877343968, 0.06386083877343968, 1.2239121278243665, -1.2239121278243665, -0.5069610024977665, 0.5069610024977665, 0.5069610024977665, -0.5069610024977665, 1.5766994272507744, -1.5766994272507744, 0.8892242114059369, -0.8892242114059369, 1.5602655313772367, -1.5602655313772367, -2.403092127540177, 2.403092127540177, 1.5602655313772367, -1.5602655313772367, 0.8892242114059369, 2.403092127540177, -0.8892242114059369, -2.403092127540177, 0.4744978678080795, 2.9239876105912574, -0.4744978678080795, -2.9239876105912574, 0.4744978678080795, -0.4744978678080795, 1.7320508075688787, -1.7320508075688787, 0.14279717659756475, 2.7367507163016924, -0.14279717659756475, -2.7367507163016924, -2.7367507163016924, 2.7367507163016924, 1.1335992635264445, -1.1335992635264445, -3.5256070994177073, 3.5256070994177073, 1.1335992635264445, -1.1335992635264445},
      {0.0, -1.2239121278243665, 1.2239121278243665, 1.2239121278243665, -1.2239121278243665, -0.06386083877343968, 0.06386083877343968, -0.5069610024977665, 0.5069610024977665, -1.5766994272507744, 1.5766994272507744, -0.5069610024977665, 0.5069610024977665, -1.5602655313772367, 1.5602655313772367, -0.8892242114059369, 0.8892242114059369, 1.5602655313772367, -1.5602655313772367, -2.403092127540177, 2.403092127540177, 2.403092127540177, 0.8892242114059369, -2.403092127540177, -0.8892242114059369, 2.9239876105912574, 0.4744978678080795, -2.9239876105912574, -0.4744978678080795, 0.4744978678080795, -0.4744978678080795, 1.7320508075688787, -1.7320508075688787, -2.7367507163016924, -0.14279717659756475, 2.7367507163016924, 0.14279717659756475, -2.7367507163016924, 2.7367507163016924, 1.1335992635264445, -1.1335992635264445, 1.1335992635264445, -1.1335992635264445, -3.5256070994177073, 3.5256070994177073}
      };
      // D3Q45
      // double const phi = (std::sqrt(5) + 1)/2; // golden ratio
      // double const phi_inv = 2/(sqrt(5) + 1) = (sqrt(5) - 1)/2
      //  --> phi + 1/phi = sqrt(5)
      //  --> 1/phi + 1 = phi
      // we can construct a rhombic dodecahedron with radius sqrt(2) as
      //    (+/-phi, +/-1/phi, 0), cyclic
      //    plus a unit cube
      //    (+/-1, +/-1, +/-1)
      // then, the centers of the triangular faces form an icosahedron
      // and the face centers are at the sum of five, e.g.
      //  (+phi, +/-1/phi, 0) (2x)
      //  (1, +/-1, 1)        (2x)
      //  (1/phi, 0, phi)     (1x)
      //  added gives (2*phi + 1 + 1/phi, 0, phi + 1)
      //              (+/-3*phi, 0, +/-(phi + 1)) --> (12x)


      
      for(int q45 = 0; q45 < 45; ++q45) {
          auto const v2 = pow2(velocity_D3Q45[0][q45]) + pow2(velocity_D3Q45[1][q45]) + pow2(velocity_D3Q45[2][q45]);
          std::printf("# D3Q45 q=%2i v %16.9f%16.9f%16.9f v^2 %g weight %g\n", 
                                  q45, velocity_D3Q45[0][q45], velocity_D3Q45[1][q45], velocity_D3Q45[2][q45], v2, weights_D3Q45[q45]);
      } // q45
      
      find_maximum_moments(D, Q, echo);
      int const numax = 4; // for D2Q37
      assert(numax < 8);
      int const M = number_of_moments(numax, D, Q, echo);
      
      
      uint32_t binom_coeff[8][8];
      { // scope: setup the binomial coefficients
          for(int k = 0; k < 8; ++k) binom_coeff[0][k] = 0;
          binom_coeff[0][0] = 1;
          for(int n = 1; n < 8; ++n) {
              binom_coeff[n][0] = 1;
              for(int k = 1; k < 8; ++k) {
                  binom_coeff[n][k] = binom_coeff[n-1][k] + binom_coeff[n-1][k-1];
                  if (k <= n) {
                      assert(binom_coeff[n][k]*factorial(k)*factorial(n - k) == factorial(n));
//                    std::printf("# n=%d k=%d %d\n", n, k, binom_coeff[n][k]);
                  } // k <= n
              } // k
          } // n
      } // scope

      std::vector<double> omega(1 + numax, 0.0); // omega = 1/tau (how much of the equilibrium_moment is added?)
//       omega[0] = 0.0; // conserve mass
//       omega[1] = 0.0; // conserve linear momentum

      omega[2] = 0.125; // relax the higher moments
      omega[3] = 0.25; // relax order 3 momenta
      omega[4] = 0.5; // fully relax order 4 momenta

      std::vector<double> projector(Q*(M + 1), 0.0);
      std::vector<double> expansion(M*Q      , 0.0);
      for(int q = 0; q < Q; ++q) {

            double H[3][8]; H[Z][0] = 1;
            for(int d = 0; d < D; ++d) {
                Hermite_polynomials(H[d], vel[q][d], numax); // if we need to expand into a rotated velocity set, apply a rotation to vel[q] first
            } // d

            auto const d2 = iv2[q]; // integer distance^2
            int m{0};
            for(int nu = 0; nu <= numax; ++nu) { // order
                int k{0};
                for(int nz = 0; nz <= nu*(D > 2); ++nz) {
                    for(int ny = 0; ny <= nu - nz; ++ny) {
                        int const nx = nu - ny - nz;
                        double const Hxyz = H[X][nx]*H[Y][ny]*H[Z][nz];
                        // if (0 == q) std::printf("# D%d m=%i nu=%d xyz= %i %i %i\n", D, m, nu, nx, ny, nz);
                        assert(m < M);
                        double const geometry_factors = binom_coeff[nu][k]/factorial(nu);
                        projector[q*(M + 1) + m] = Hxyz;
                        expansion[m*Q + q] = w[d2]*Hxyz*geometry_factors;
                        ++m;
                        ++k;
                    } // nx
                } // nz
            } // nu
            assert(M == m);

            // to compute the kinetic energy:
            projector[q*(M + 1) + M] = vel[q][v2]/D;

      } // q
 
 
      // Factorization:
      //      projector and expansion matrix can be factorized. For eaxmple D2Q37
      //      the projection operations take 37*2*15 = 1110 Flop if performed 
      //      as matrix multiplication like implemented here.
      //      The projector matrix has 15*37 = 555 elements, the populations have 37 and the moments 15 elements
      //      Using factorization, we have a 2 step procedure:
      //        Step 1) transforming from 37 populations in 37*2*5 = 370 Flop to 7*5 f_nx(y) 
      //        Step 2) transforming from 35 intermediates in 7*2*15 = 210 Flop to 15 moments
      //        ... so we would save about half of the Flops
      //      The factorized matrices have 2*7*5 = 70 elements or could even be computed on-the-fly
      
      // for the reconstruction/expansion, we have to weight 
      //      the 15 moments with geometry_factors[m] and 
      //      the 37 populations with population_weights[q]
      //      so an additional 15+37=52 numbers are needed.
      //      
      //      Geometry factors are
      //              1
      //            1   1
      //          1/2   2/2   1/2
      //        1/6   3/6   3/6   1/6
      //      1/24   4/24   6/24   4/24   1/24
      //      could be scaled up by a common denominator 24:
      //              24
      //            24  24
      //          12  24  12
      //        4   12  12  4
      //      1   4   6   4   1
      //      so they fit into 15 int8_t
      //      and the factor 1/24 is merged into the population weights
      //      
      //      population_weights can be compressed to 8 doubles using an indirection
      //      indirection indices could be compressed into 37*3 bit < 16 Byte
      //      or 37*4 bit < 19 Byte
      
      //      if we do not use omega_nu for nu == 0 and nu == 1
      //      but branch for that case, we do not need to store the (trivial) geometry_factors (3x 24)
      //      then, 19 Byte and 12 Byte fit into one cache line
      
      //      for the projection and expansion with Hermite polynomials, we also need the 3bit info, about the y-velocity
      //        
      //      double tmp[5][8];
      //      {
      //          for (int q = 0; q < 37; ++q) {
      //              int const cy =  xy[q]       & 0x7;
      //              int const cx = (xy[q] >> 4) & 0x7;
      //              for (int nx = 0; nx <= 4; ++nx) {
      //                  tmp[nx][cy] += f_in[q] * H[cx][nx]; // 37*5*2 = 370 Flop
      //              } // nx
      //          } // q
      //      }
      //      double moment[15];
      //      {
      //          int m{0};
      //          for (int nu = 0; nu <= 4; ++nu) {
      //              for (int ny = 0; ny < nu; ++ny) {
      //                  int const nx = nu - ny;
      //                  double mom{0}
      //                  for (int cy = 0; cy < 7; ++cy) {
      //                      mom += tmp[nx][cy] * H[cy][ny]; // 15*7*2 = 210 Flop
      //                  } // cy
      //                  moment[m] = mom;
      //                  ++m;
      //              } // ny
      //          } // nu
      //          assert(15 == m);
      //      }
      //      //
      //      // 1 inversion + 7 Flop
      //      // 
      //      // How to compute the kinetic energy? (moment[3] + moment[5] + 2*moment[0])/(rho*D);
      //      // // moment[3], moment[5] are cx^2 - 1 and cy^2 - 1
      //      //
//             double Hermite_eq[2][5];
//             double const mTm1 = 1 - temperature; // 1 Flop
//             Hermite_polynomials_saving(Hermite_eq[0], u[0], 4, mTm1); // 16 Flop
//             Hermite_polynomials_saving(Hermite_eq[1], u[1], 4, mTm1); // 16 Flop
//             { // scope: relax towards equilibrium momenta
//                 int m{0};
//                 for(int nu = 0; nu <= numax; ++nu) { // order
//                     auto const omega_nu = (nu > 1) ? omega[nu] : 0;
//                     auto const one_minus_omega_nu = 1 - omega_nu; // 5 Flop
//                     for(int ny = 0; ny <= nu - nz; ++ny) {
//                         int const nx = nu - ny - nz;
//                         auto const equilibrium_moment = rho*Hermite_eq[0][nx]
//                                                            *Hermite_eq[1][ny]; // 15*2 Flop
//                         auto const old_moment = moment[m];
//                         moment[m] = one_minus_omega_nu*old_moment + omega_nu*equilibrium_moment; // 15*3 Flop
//                         ++m;
//                     } // ny
//                 } // nu
//                 assert(M == m);
//             } // scope
      //
      //      // and expansion (double loops nu,ny can be fused)
      //      {
      //          int m{0};
      //          for (int nu = 0; nu <= 4; ++nu) {
      //              for (int ny = 0; ny < nu; ++ny) {
      //                  int const nx = nu - ny;
      //                  moment[m] *= geometry_factors[m]; // 15 Flop
      //                  for (int cy = 0; cy < 7; ++cy) {
      //                      tmp[nx][cy] += H[cy][ny] * moment[m]; // 15*7*2 = 210 Flop
      //                  } // cy
      //                  ++m;
      //              } // ny
      //          } // nu
      //          assert(15 == m);
      //      }
      //      {
      //          for (int q = 0; q < 37; ++q) {
      //              int const cy =  xy[q]       & 0x7;
      //              int const cx = (xy[q] >> 4) & 0x7;
      //              double f{0};
      //              for (int nx = 0; nx <= 4; ++nx) {
      //                  f += H[cx][nx] * tmp[nx][cy]; // 37*5*2 = 370 Flop
      //              } // nx
      //              double const w = weights_d2[cx*cx + cy*cy];
      //              f_out[q] = w*f; // 37 Flop
      //          } // q
      //      }
      //      
      
      // Flop count D2Q37 is 370 + 210 + 7 + 4 + 1 + 16 + 16 + 5 + 30 + 45 + 15 + 210 + 370 + 37 = 1336 Flop
      // plus 1 inversion
      // most but not all flops are FMA
      // Arithmetic intensity:
      //    How much loads and stores? (only those quantities that are site-intensive)
      //    37*2*8 Byte = 592 Byte
      //    AI = 2.25 Flop/Byte
      //    A100: 9.5 TFlop/s over 1 TByte/s --> memory bound
      //  ==> complexity of the custom version for D2Q37 with Hermite polynomials is not necessary
      //      yes, we do save some bandwidth since H is much smaller than projector and expansion matrix
      //      but that is not so much if we use long vectors (e.g. 64 CUDA threads per block)
      //      since the matrices will be cached

      // Flop count D2Q37 with matrices 1184 + 7 + 4 + 1 + 16 + 16 + 5 + 30 + 45 + 1110 = 2418 Flop
      //  37*2*8 + 37*(16+15)*8/vector = 735.375 for vector=64
      //  ==> AI = 3.288
      //  if we consider the matrices to be cached
      //  --> AI = 4.08
      //  otherwise, we could consider to constuct the matrices 
      //      on-the-fly from ivelocity[][2],weight_d2[],geometry_factors[],r=1.19697977

      if (1) {
          auto const dev = check_overlap(projector.data(), expansion.data(), Q, M, echo/2);
          if (echo > 1) std::printf("# %s: projector*expansion deviates max %.1e from unity\n", __func__, dev);
//        std::exit(__LINE__);
      } // DEBUG


      
      if (0) {
          std::vector<real_t> moments(M + 1, 0.0);

          // inital moments
          moments[0] = 1.0; // density
          for(int d = 0; d < D; ++d) {
              moments[1 + d] = 0.0 + d*0.1; // linear momentum
          } // d
          moments[M] = moments[0]; // kinetic energy

          int constexpr Ntime=5;
          int input{0};
          int const reference=2;

          
          std::vector<real_t> f[3]; for(int i = 0; i < 3; ++i) f[i].resize(Q, 0.0);

          int output{-1};
          for(int time = 0; time < Ntime; ++time) {
              output = 1 - input;
              collide_general(
                      f[output].data()
                    , moments.data()
                    , time ? f[input].data() : nullptr
                    , projector.data()
                    , expansion.data()
                    , D
                    , Q
                    , M
                    , numax
                    , omega.data()
                    , echo
                    , vel // for CHECK_ME
                    );

              if (0 == time) {
                  for(int q = 0; q < Q; ++q) {
                      f[reference][q] = f[output][q];
                  } // q
              } // after 1st call

              input = 1 - input; // switch between 0 and 1 back and forth
          } // time
          assert(-1 != output && "Run at least one time step!"); 

          double const threshold = (sizeof(real_t) < 8) ? 2e-9 : 2e-15;
          // check that the populations have kept their content exactly
          double maxdev{0};
          for(int q = 0; q < Q; ++q) {
              auto const v_ref = f[reference][q];
              auto const v_out = f[output   ][q]; // get last output
              if (echo > 6) {
                  std::printf("# D%dQ%d collide: init f[%i]=\t%16.9f\treference=%16.9f\n",
                                D, Q,                   q, v_out,         v_ref);
              } // echo
              double const dev = std::abs(v_out - v_ref);
              stat += (dev > threshold);
              maxdev = std::max(maxdev, dev);
          } // q
          std::printf("# D%dQ%d collide: largest deviation is %.1e\n", D, Q, maxdev);
          
      } // test stability
      
      

      
      
      if (1) {
          int const n1 = 1 << 8;
          int const n[3] = {n1, n1, 1};
          view4D<real_t> moments(n[2], n[1], n[0], M + 1, 0.0);

          // initalize moments

          if (0) {
              double const droplet0_position_and_radius[] = {n[0]*.5 + n[0]*.125, n[1]*.5, 0, n[0]*.125};
              double const droplet1_position_and_radius[] = {n[0]*.5 - n[0]*.125, n[1]*.5, 0, n[0]*.125};
              initalize_droplet(moments, n, droplet0_position_and_radius);
              initalize_droplet(moments, n, droplet1_position_and_radius);
          } else {
              initalize_Taylor_Green_2D(moments, n); assert(1 == n[2]);
          }

          for (int xyz = 0; xyz < n[2]*n[1]*n[0]; ++xyz) {
              moments(0,0,xyz,M) = moments(0,0,xyz,0); // kinetic energy = rho
          } // xyz

          view4D<real_t> f[2]; // population arrays
          f[0] = view4D<real_t>(n[2], n[1], n[0], Q + 3, 0.0);
          f[1] = view4D<real_t>(n[2], n[1], n[0], Q + 3, 0.0);
          int input{0}, output{-1};

          int constexpr Ntime = 333;
          for(int time = 0; time < Ntime; ++time) {
              output = 1 - input;
              
              if (0 == (time & 0x0)) { // monitor

                  view3D<float> rho_u(4, n[0], n[1], 0.f); // {rho,ux,uy,uz}[Ny][Nx]
                  int const iz = 0;
                  for (int iy = 0; iy < n[1]; ++iy) {
                      for (int ix = 0; ix < n[0]; ++ix) {
                          auto const rho = moments(iz,iy,ix,0);
                          rho_u(0,iy,ix) = rho;
                          if (rho > 0) {
                              auto const rho_inv = 1/rho;
                              for(int d = 0; d < D; ++d) {
                                  rho_u(1 + d,iy,ix) = moments(iz,iy,ix,1 + d)*rho_inv;
                              } // d
                          } // rho > 0
                      } // ix
                  } // iy
//                   std::printf("\n# time step %i density\n", time);
//                   lbm_monitor::show_scalar(rho_u.data(), n[0], n[1]); // density

                  
                  if (1) {
                      view3D<float> vorticity(n[2], n[1], n[0], 0.f);
                      for (int iy = 0; iy < n[1]; ++iy) {
                          int const iym = (iy - 1 + n[1]) % n[1];
                          int const iyp = (iy + 1 + n[1]) % n[1];
                          for (int ix = 0; ix < n[0]; ++ix) {
                              // compute the vorticity:
                              //          dv_z/dy - dv_y/dz
                              //    v =   dv_x/dz - dv_z/dx
                              //          dv_y/dx - dv_x/dy    (only this component is non-zero in 2D flow)
                              int const ixm = (ix - 1 + n[0]) % n[0];
                              int const ixp = (ix + 1 + n[0]) % n[0];
                              vorticity(iz,iy,ix) = moments(iz,iy,ixp,2) - moments(iz,iy,ixm,2)
                                                  - moments(iz,iyp,ix,1) + moments(iz,iym,ix,1);
                          } // ix
                      } // iy
//                       std::printf("\n# time step %i vorticity\n", time);
//                       lbm_monitor::show_scalar(vorticity.data(), n[0], n[1]);

                      std::printf("\n# time step %i velocity field and vorticity\n", time);
                      lbm_monitor::show_vector_field(
                                            rho_u[1].data(), rho_u[2].data(), rho_u[3].data(),
                                            n[0], n[1], 1, 0, nullptr, nullptr, vorticity.data());
                  } else { 

                      std::printf("\n# time step %i velocity field and density\n", time);
                      double const mass = lbm_monitor::show_vector_field(
                                            rho_u[1].data(), rho_u[2].data(), rho_u[3].data(), 
                                            n[0], n[1], 1, 0, nullptr, nullptr, rho_u.data());
                      std::printf("# mass at time step %i is %.15g\n\n", time, mass);
                  } // compute vorticity?

              } // monitor

              
              for (int iz = 0; iz < n[2]; ++iz) {
              for (int iy = 0; iy < n[1]; ++iy) {
              for (int ix = 0; ix < n[0]; ++ix) {
                  real_t tmp_f[128];
                  real_t const *const f_in = time ? &f[input](iz,iy,ix,0) : nullptr;
                  collide_general(
                          tmp_f
                        , &moments(iz,iy,ix,0)
                        , f_in
                        , projector.data(), expansion.data()
                        , D, Q
                        , M, numax, omega.data()
                        , echo/2
                        , (time > 999)? vel : nullptr // for CHECK_ME
                        );
                  propagate_general(f[output], tmp_f, D, Q, ivel, ix, iy, iz, n[0], n[1], n[2], echo);
              } // ix
              } // iy
              } // iz

              input = 1 - input; // switch between 0 and 1 back and forth
          } // time
          assert(-1 != output && "Run at least one time step!"); 
          
      } // test run
      

      if (0 != int(stat) && echo > 0) std::printf("# %s<D%dQ%d>: found %d errors!\n",
                                                    __func__, D, Q, int(stat));
      return stat;
  } // test_collide

  template <typename real_t>
  inline status_t test_stencils(int const echo=0) {
      status_t stat(0);
      stat += test_collide<real_t>(2, echo);
//       stat += test_collide<real_t>(3, echo);
      return 0;
  } // test_stencils

  inline status_t test_Hermite_version(int const echo=0) {
      int constexpr numax = 15;
      double H[1 + numax], Hs[1 + numax], x=.6;
      Hermite_polynomials(H, x, numax);
      Hermite_polynomials_saving(Hs, x, numax);
      double dev{0};
      for(int nu = 0; nu <= numax; ++nu) {
          dev = std::max(dev, std::abs(H[nu] - Hs[nu]));
      } // nu
      return (dev > 1e-12);
  } // test_Hermite_version
  
  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_Hermite_version(echo);
      stat += test_Hermite_orthogonality(echo);
      stat += test_stencils<double>(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_general
