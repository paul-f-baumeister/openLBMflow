#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "lbm_stencil.hxx" // q_ooo, ..., q_nop, ...


#ifdef  STANDALONE_TEST
        #define restrict
#endif

    template <typename real_t, class stencil_t>
    inline int collide_D3Q19(
          real_t       *restrict const f_out // result vector[Q]
        , real_t const *restrict const f_in  // input  vector[Q]
        , stencil_t const & stencil
        , double const tau=0.5
        , double const body_force_xyz[3]=nullptr
    ) {
#ifdef  STANDALONE_TEST
        assert(f_out != f_in); // never pass the same pointer (forbidden by restrict keywords)
#endif
        int constexpr Q = stencil_t::Q;
        assert(19 == Q);

        // relaxation time constant tau
        double const inv_tau = 1.0/tau;
        double const min_tau = 1.0 - inv_tau;
        double const half_wi0 = 0.5*stencil.weight(0), 
                     half_wi1 = 0.5*stencil.weight(1),
                     half_wi2 = 0.5*stencil.weight(2);

                    double const f_ooo = f_in[q_ooo];
                    double const f_poo = f_in[q_poo];
                    double const f_ppo = f_in[q_ppo];
                    double const f_opo = f_in[q_opo];
                    double const f_npo = f_in[q_npo];
                    double const f_noo = f_in[q_noo];
                    double const f_nno = f_in[q_nno];
                    double const f_ono = f_in[q_ono];
                    double const f_pno = f_in[q_pno];
                    double const f_opp = f_in[q_opp];
                    double const f_oop = f_in[q_oop];
                    double const f_onp = f_in[q_onp];
                    double const f_onn = f_in[q_onn];
                    double const f_oon = f_in[q_oon];
                    double const f_opn = f_in[q_opn];
                    double const f_pop = f_in[q_pop];
                    double const f_nop = f_in[q_nop];
                    double const f_non = f_in[q_non];
                    double const f_pon = f_in[q_pon];

                    // calculate rho and ux, uy, uz
                    double const tmp_rho = f_ooo 
                                    + (f_poo + f_noo) 
                                    + (f_opo + f_ono) 
                                    + (f_oop + f_oon)
                                    
                                    + (f_opp + f_onn)
                                    + (f_onp + f_opn)
                                    + (f_pop + f_non)
                                    + (f_nop + f_pon)
                                    + (f_ppo + f_nno)
                                    + (f_npo + f_pno);
                                
                    double const inv_rho = 1.0/tmp_rho;
                    
                    // x-current        p      n
                    double tmp_ux = ( (f_poo - f_noo) 
                                    + (f_ppo - f_nno) 
                                    + (f_pno - f_npo)
                                    + (f_pop - f_non)
                                    + (f_pon - f_nop) )*inv_rho;

                    // y-current         p      n
                    double tmp_uy = ( (f_opo - f_ono) 
                                    + (f_ppo - f_nno) 
                                    + (f_npo - f_pno)
                                    + (f_opp - f_onn)
                                    + (f_opn - f_onp) )*inv_rho;

                    // z-current          p      n
                    double tmp_uz = ( (f_oop - f_oon) 
                                    + (f_pop - f_non) 
                                    + (f_nop - f_pon)
                                    + (f_opp - f_onn)
                                    + (f_onp - f_opn) )*inv_rho;
 
                if (nullptr != body_force_xyz) {
                    // add the body force
                    tmp_ux += tau*body_force_xyz[0];
                    tmp_uy += tau*body_force_xyz[1];
                    tmp_uz += tau*body_force_xyz[2];
                }
//                     rho[xyz] = tmp_rho; // store density
//                     ux[xyz] = tmp_ux;
//                     uy[xyz] = tmp_uy; // store current directions
//                     uz[xyz] = tmp_uz;
                    auto const ux2 = pow2(tmp_ux);
                    auto const uy2 = pow2(tmp_uy);
                    auto const uz2 = pow2(tmp_uz);
                    auto const uxyz2 = ux2 + uy2 + uz2;

                    // equilibrium(hw, uv, u2) = hw*(2 + 6*uv - 3*u2 + 9*uv*uv);
                    
                    auto const tmp_rho_inv_tau = tmp_rho*inv_tau;
                    auto const w0 = half_wi0*tmp_rho_inv_tau;
                    f_out[q_ooo] = f_ooo*min_tau + w0*(2 - 3*uxyz2);

                    auto const w1 = half_wi1*tmp_rho_inv_tau;
                    f_out[q_poo] = f_poo*min_tau + w1*(2 + 6*tmp_ux + 9*ux2 - 3*uxyz2);
                    f_out[q_noo] = f_noo*min_tau + w1*(2 - 6*tmp_ux + 9*ux2 - 3*uxyz2);
                    f_out[q_opo] = f_opo*min_tau + w1*(2 + 6*tmp_uy + 9*uy2 - 3*uxyz2);
                    f_out[q_ono] = f_ono*min_tau + w1*(2 - 6*tmp_uy + 9*uy2 - 3*uxyz2);
                    f_out[q_oop] = f_oop*min_tau + w1*(2 + 6*tmp_uz + 9*uz2 - 3*uxyz2);
                    f_out[q_oon] = f_oon*min_tau + w1*(2 - 6*tmp_uz + 9*uz2 - 3*uxyz2);

                    auto const uxy2 = ux2 + uy2;
                    auto const uyz2 = uy2 + uz2;
                    auto const uzx2 = uz2 + ux2;
                    auto const uxy = 2*tmp_ux*tmp_uy;
                    auto const uyz = 2*tmp_uy*tmp_uz;
                    auto const uzx = 2*tmp_uz*tmp_ux;

                    auto const w2 = half_wi2*tmp_rho_inv_tau;
                    f_out[q_ppo] = f_ppo*min_tau + w2*(2 + 6*(+tmp_ux + tmp_uy) + 9*(uxy2 + uxy) - 3*uxyz2);
                    f_out[q_npo] = f_npo*min_tau + w2*(2 + 6*(-tmp_ux + tmp_uy) + 9*(uxy2 - uxy) - 3*uxyz2);
                    f_out[q_nno] = f_nno*min_tau + w2*(2 + 6*(-tmp_ux - tmp_uy) + 9*(uxy2 + uxy) - 3*uxyz2);
                    f_out[q_pno] = f_pno*min_tau + w2*(2 + 6*(+tmp_ux - tmp_uy) + 9*(uxy2 - uxy) - 3*uxyz2);

                    f_out[q_opp] = f_opp*min_tau + w2*(2 + 6*(+tmp_uy + tmp_uz) + 9*(uyz2 + uyz) - 3*uxyz2);
                    f_out[q_onp] = f_onp*min_tau + w2*(2 + 6*(-tmp_uy + tmp_uz) + 9*(uyz2 - uyz) - 3*uxyz2);
                    f_out[q_onn] = f_onn*min_tau + w2*(2 + 6*(-tmp_uy - tmp_uz) + 9*(uyz2 + uyz) - 3*uxyz2);
                    f_out[q_opn] = f_opn*min_tau + w2*(2 + 6*(+tmp_uy - tmp_uz) + 9*(uyz2 - uyz) - 3*uxyz2);
                    
                    f_out[q_pop] = f_pop*min_tau + w2*(2 + 6*(+tmp_ux + tmp_uz) + 9*(uzx2 + uzx) - 3*uxyz2);
                    f_out[q_nop] = f_nop*min_tau + w2*(2 + 6*(-tmp_ux + tmp_uz) + 9*(uzx2 - uzx) - 3*uxyz2);
                    f_out[q_non] = f_non*min_tau + w2*(2 + 6*(-tmp_ux - tmp_uz) + 9*(uzx2 + uzx) - 3*uxyz2);
                    f_out[q_pon] = f_pon*min_tau + w2*(2 + 6*(+tmp_ux - tmp_uz) + 9*(uzx2 - uzx) - 3*uxyz2);
        
        
        std::printf("# D3Q%d collide(rho= %g, u= %g %g %g)\n", 
                          Q, 0., 0., 0., 0.);

        return 0;
    } // collide_D3Q19

namespace lbm_collide {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <int D, int Q, typename real_t=double>
  inline status_t test_collide(int const echo=0) {
      status_t stat(0);

      BKG_stencil<3,19> stencil(9);

      int constexpr Ntime=9;
      int input{0};
      int const reference=2;
      
      real_t f[3][Q];
      double const rho = 1;
      double const u[3] = {0.5, 0.2, -.15};
      auto const u2 = pow2(u[0]) + pow2(u[1]) + pow2(u[2]);
      for(int q = 0; q < Q; ++q) {
          auto const v = stencil.velocity(q);
          double const uv = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
          double const weight = stencil.weight(v[3]);
          double const value = stencil.equilibrium(0.5*rho*weight, uv, u2);
          f[reference][q] = value;
          f[input    ][q] = value;
          if (echo > 21) {
              std::printf("# D%dQ%d collide: init f[%i]= %g\n",
                             D, Q,                   q, value);
          } // echo
      } // q

      int output{-1};
      for(int time = 0; time < Ntime; ++time) { // run an even number of time steps
          output = 1 - input;
          int const n0 = collide_D3Q19<real_t> (f[output],
                                      f[input], stencil);
          assert(0 == n0);
          input = 1 - input; // switch between 0 and 1 back and forth
      } // time
      assert(-1 != output && "Run at least one time step!"); 

      double const threshold = (sizeof(real_t) < 8) ? 2e-9 : 2e-15;
      // check that the populations have kept their content exactly
      for(int q = 0; q < Q; ++q) {
          auto const v_ref = f[reference][q];
          auto const v_out = f[output   ][q]; // get last output
          if (echo > 11) {
              std::printf("# D%dQ%d propagate: init f[%i]= %g  reference= %g\n",
                             D, Q,                     q, v_out,         v_ref);
          } // echo
          stat += (std::abs(v_out - v_ref) > threshold);
      } // q

      if (0 != int(stat) && echo > 0) std::printf("# %s<D%dQ%d>: found %d errors!\n",
                                                    __func__, D, Q, int(stat));
      return stat;
  } // test_collide

  template <typename real_t>
  inline status_t test_stencils(int const echo=0) {
      status_t stat(0);
      stat += test_collide<3,19,real_t>(99);
      return 0;
  } // test_stencils

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_stencils<double>(echo);
//       stat += test_stencils<float> (echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_collide
