#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "lbm_stencil.hxx" // q_ooo, ..., q_nop, ...

    template <typename real_t, class Stencil>
    inline int propagate(
#ifdef  STANDALONE_TEST
        #define restrict
          real_t       *restrict const fnew  // output
        #define f_out(x, y, z, q)      fnew[((x*Ny + y)*Nz + z)*Stencil::Q + q]
#else  // STANDALONE_TEST
          view4D<real_t> & f_out // output --> view4D is an object that can transpose at runtime
#endif // STANDALONE_TEST
        , real_t const *restrict const f_in // input vector[Q]
        , int const x, int const y, int const z
        , int const Nx, int const Ny, int const Nz // local domain limits
    ) {
        int number_of_populations{0}; // result (should match Q on exit in STANDALONE_TEST mode)
#ifdef  STANDALONE_TEST
        assert(fnew != f_in); // never pass the same pointer (forbidden by restrict keywords)
        #define count_pop(number) number_of_populations += number
#else
        #define count_pop(number) // do nothing (and allow the compiler to reorder)
#endif

        // exploit that we can assume the local domain to be periodic
        int const oz = z;
        int const oy = y;
        int const ox = x;
        int const px = (x >= Nx - 1) ? 0 : x + 1; // periodic
        int const nx = (x <= 0) ? Nx - 1 : x - 1; // periodic
        
        // depending of their velocities, the polulations in f_in[Q] at (x,y,z)
        // are streamed to the other lattice sites

        // all stencils: D1Q3, D2Q5, D3Q7, D2Q9, D3Q15, D3Q19, D3Q27
        f_out(ox, oy, oz, q_ooo) = f_in[q_ooo]; //
        f_out(px, oy, oz, q_poo) = f_in[q_poo]; // +x
        f_out(nx, oy, oz, q_noo) = f_in[q_noo]; // -x
        count_pop(3);

        if (Stencil::D >= 2) {
            int const py = (y >= Ny - 1) ? 0 : y + 1; // periodic
            int const ny = (y <= 0) ? Ny - 1 : y - 1; // periodic

            // stencils: D2Q5, D3Q7, D2Q9, D3Q15, D3Q19, D3Q27
            f_out(ox, py, oz, q_opo) = f_in[q_opo]; //   +y
            f_out(ox, ny, oz, q_ono) = f_in[q_ono]; //   -y
            count_pop(2);

            if (Stencil::Q > 5) {

                int constexpr k12 = (Stencil::is_D3Q15() || Stencil::is_D3Q7()) ? 12 : 0;
                if (0 == k12) {
                    // stencils: D2Q9, D3Q19, D3Q27
                    f_out(px, py, oz, q_ppo) = f_in[q_ppo]; // +x+y
                    f_out(nx, py, oz, q_npo) = f_in[q_npo]; // -x+y
                    f_out(nx, ny, oz, q_nno) = f_in[q_nno]; // -x-y
                    f_out(px, ny, oz, q_pno) = f_in[q_pno]; // +x-y
                    count_pop(4);
                } // not for D3Q7 nor for D3Q15

                if (3 == Stencil::D) {
                    int const pz = (z >= Nz - 1) ? 0 : z + 1; // periodic
                    int const nz = (z <= 0) ? Nz - 1 : z - 1; // periodic

                    if (0 == k12) {
                        // stencils: D3Q19, D3Q27
                        f_out(ox, py, pz, q_opp) = f_in[q_opp]; //   +y+z
                        f_out(ox, ny, pz, q_onp) = f_in[q_onp]; //   -y+z
                        f_out(ox, ny, nz, q_onn) = f_in[q_onn]; //   -y-z
                        f_out(ox, py, nz, q_opn) = f_in[q_opn]; //   +y-z
                        
                        f_out(px, oy, pz, q_pop) = f_in[q_pop]; // +x  +z
                        f_out(nx, oy, pz, q_nop) = f_in[q_nop]; // -x  +z
                        f_out(nx, oy, nz, q_non) = f_in[q_non]; // -x  -z
                        f_out(px, oy, nz, q_pon) = f_in[q_pon]; // +x  -z
                        count_pop(8);
                    } // not for D3Q7 nor for D3Q15

                    // stencils: D3Q7, D3Q15, D3Q19, D3Q27 (all 3D stencils)  
                    f_out(ox, oy, pz, q_oop -k12) = f_in[q_oop -k12]; //     +z
                    f_out(ox, oy, nz, q_oon -k12) = f_in[q_oon -k12]; //     -z
                    count_pop(2);

                    if (19 != Stencil::Q) {
                        // stencils: D3Q15, D3Q27  
                        f_out(px, py, pz, q_ppp -k12) = f_in[q_ppp -k12]; // +x+y+z
                        f_out(nx, py, pz, q_npp -k12) = f_in[q_npp -k12]; // -x+y+z
                        f_out(px, ny, pz, q_pnp -k12) = f_in[q_pnp -k12]; // +x-y+z
                        f_out(nx, ny, pz, q_nnp -k12) = f_in[q_nnp -k12]; // -x-y+z
                        f_out(px, py, nz, q_ppn -k12) = f_in[q_ppn -k12]; // +x+y-z
                        f_out(nx, py, nz, q_npn -k12) = f_in[q_npn -k12]; // -x+y-z
                        f_out(px, ny, nz, q_pnn -k12) = f_in[q_pnn -k12]; // +x-y-z
                        f_out(nx, ny, nz, q_nnn -k12) = f_in[q_nnn -k12]; // -x-y-z
                        count_pop(8);
                    } // not D3Q19
                } // D == 3
            } // Q > 5
        } // D >= 2

//     q_ooo =  0, //         wi0
//     q_poo =  1, // +x      wi1
//     q_noo =  2, // -x      wi1
//     up to here for D1Q3
//     q_opo =  3, //   +y    wi1
//     q_ono =  4, //   -y    wi1
//     up to here for D2Q5
//     q_ppo =  5, // +x+y    wi2
//     q_npo =  6, // -x+y    wi2
//     q_nno =  7, // -x-y    wi2
//     q_pno =  8, // +x-y    wi2
//     up to here for D2Q9
//     q_opp =  9, //   +y+z  wi2
//     q_onp = 10, //   -y+z  wi2
//     q_onn = 11, //   -y-z  wi2
//     q_opn = 12, //   +y-z  wi2
//     q_pop = 13, // +x  +z  wi2
//     q_nop = 14, // -x  +z  wi2
//     q_non = 15, // -x  -z  wi2
//     q_pon = 16, // +x  -z  wi2
//     q_oop = 17, //     +z  wi1
//     q_oon = 18, //     -z  wi1
//     up to here for D3Q19
//     q_ppp = 19, // +x+y+z  wi3
//     q_npp = 20, // -x+y+z  wi3
//     q_pnp = 21, // +x-y+z  wi3
//     q_nnp = 22, // -x-y+z  wi3
//     q_ppn = 23, // +x+y-z  wi3
//     q_npn = 24, // -x+y-z  wi3
//     q_pnn = 25, // +x-y-z  wi3
//     q_nnn = 26; // -x-y-z  wi3
//     up to here for D3Q27

//         std::printf("# D%dQ%d propagate(%i, %i, %i)\n", 
//            Stencil::D, number_of_populations, x, y, z);

        return number_of_populations;
    } // propagate

typedef int status_t;

namespace lbm_propagate {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <int D, int Q, typename real_t=double>
  inline status_t test_propagation(int const echo=0) {
      status_t stat(0);

      int constexpr Nx=3, Ny=5, Nz=7;
      int constexpr Ntime=Nx*Ny*Nz; // make sure the number of time steps is a (or the least) common multiple
      int input{0};
      int const reference=2;
      real_t f[3][Nx*Ny*Nz*Q];
      for(int xyz = 0; xyz < Nx*Ny*Nz; ++xyz) {
          for(int q = 0; q < Q; ++q) {
              int const index = xyz*Q + q;
              int const value = xyz*(1 + Q) + (1 + q); // initialize with 1+q (to avoid zero values)
              f[reference][index] = value;
              f[input    ][index] = value;
              if (echo > 21) {
                  std::printf("# D%dQ%d propagate: init(xyz=%i, q=%i)= f[%i]= %i\n",
                                 D, Q, xyz, q, index, value);
              } // echo
          } // q
      } // xyz

      int output{-1};
      for(int time = 0; time < Ntime; ++time) { // run an even number of time steps
          output = 1 - input;
          for(int x = 0; x < Nx; ++x) {
              for(int y = 0; y < Ny; ++y) {
                  for(int z = 0; z < Nz; ++z) {
                      int const xyz = (x*Ny + y)*Nz + z;
                      int const nq = propagate<real_t,BKG_stencil<D,Q>> (f[output],
                                            f[input] + xyz*Q, x, y, z, Nx, Ny, Nz);
                      assert(Q == nq);
                  } // z
              } // y
          } // x
          input = 1 - input; // switch between 0 and 1 back and forth
      } // time
      assert(-1 != output && "Run at least one time step!"); 

      // check that the populations have kept their content exactly
      for(int xyz = 0; xyz < Nx*Ny*Nz; ++xyz) {
          for(int q = 0; q < Q; ++q) {
              int const index = xyz*Q + q;
              auto const v_out = f[output][index]; // get last output
              auto const v_ref = f[reference][index]; 
              if (echo > 11) {
                  std::printf("# D%dQ%d propagate: final(xyz=%i, q=%i)= f[%i]= %g reference= %i\n", 
                                 D, Q,                   xyz,    q,     index, double(v_out), double(v_ref));
              } // echo
              stat += (v_out != v_ref);
          } // q
      } // xyz
      
      if (0 != int(stat) && echo > 0) std::printf("# %s<D%dQ%d>: found %d errors!\n",
                                                    __func__, D, Q, int(stat));
      return stat;
  } // test_propagation

  template <typename real_t>
  inline status_t test_stencils(int const echo=0) {
      status_t stat(0);
      stat += test_propagation<1, 3,real_t>(echo);
      stat += test_propagation<2, 5,real_t>(echo);
      stat += test_propagation<2, 9,real_t>(echo);
      stat += test_propagation<3,15,real_t>(echo);
      stat += test_propagation<3,19,real_t>(echo);
      stat += test_propagation<3,27,real_t>(echo);
      return 0;
  } // test_stencils

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_stencils<double>(echo);
      stat += test_stencils<float> (echo);
      stat += test_stencils<int>   (echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_propagate
