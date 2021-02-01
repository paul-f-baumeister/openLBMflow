#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf
#include <cstdint> // uint8_t, int8_t
#include <cmath> // std::abs

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

template <typename real_t>
inline real_t pow2(real_t const v) { return v*v; } // ToDo: move to utils

  // order should, in principle, not matter, but this ordering 
  // can be kept constant for D1Q3, D2Q5, D2Q9, D3Q19, D3Q27
  uint8_t constexpr
    q_ooo =  0, //         wi0
    q_poo =  1, // +x      wi1
    q_noo =  2, // -x      wi1
    // up to here for D1Q3
    q_opo =  3, //   +y    wi1
    q_ono =  4, //   -y    wi1
    // up to here for D2Q5
    q_ppo =  5, // +x+y    wi2
    q_npo =  6, // -x+y    wi2
    q_nno =  7, // -x-y    wi2
    q_pno =  8, // +x-y    wi2
    // up to here for D2Q9
    q_opp =  9, //   +y+z  wi2
    q_onp = 10, //   -y+z  wi2
    q_onn = 11, //   -y-z  wi2
    q_opn = 12, //   +y-z  wi2
    q_pop = 13, // +x  +z  wi2
    q_nop = 14, // -x  +z  wi2
    q_non = 15, // -x  -z  wi2
    q_pon = 16, // +x  -z  wi2
    q_oop = 17, //     +z  wi1
    q_oon = 18, //     -z  wi1
    // up to here for D3Q19
    q_ppp = 19, // +x+y+z  wi3
    q_npp = 20, // -x+y+z  wi3
    q_pnp = 21, // +x-y+z  wi3
    q_nnp = 22, // -x-y+z  wi3
    q_ppn = 23, // +x+y-z  wi3
    q_npn = 24, // -x+y-z  wi3
    q_pnn = 25, // +x-y-z  wi3
    q_nnn = 26; // -x-y-z  wi3
    // up to here for D3Q27

//   uint8_t constexpr // special for D3Q7 and D3Q15
//     q_oop_D3Q15 =  5, //     +z  wi1
//     q_oon_D3Q15 =  6, //     -z  wi1
//     q_ppp_D3Q15 =  7, // +x+y+z  wi3
//     q_npp_D3Q15 =  8, // -x+y+z  wi3
//     q_pnp_D3Q15 =  9, // +x-y+z  wi3
//     q_nnp_D3Q15 = 10, // -x-y+z  wi3
//     q_ppn_D3Q15 = 11, // +x+y-z  wi3
//     q_npn_D3Q15 = 12, // -x+y-z  wi3
//     q_pnn_D3Q15 = 13, // +x-y-z  wi3
//     q_nnn_D3Q15 = 14, // -x-y-z  wi3
//     q_oop_D3Q7 = q_oop_D3Q15,
//     q_oon_D3Q7 = q_oon_D3Q15;

  // put D and Q into one number to make a faster switching
  inline constexpr int D_Q_code(int const D, int const Q) { return D*4 + Q; }

template <int Dimensions, int Populations>
class BKG_stencil {
    // a Lattice Boltzmann stencil for the BKG equation
public:
    static int constexpr D = Dimensions; // number of dimensions
    static int constexpr Q = Populations; // number of populations

public:

    BKG_stencil(int const echo=3) { // default constructor

        // setup of quadrature weights
        for(int v2 = 0; v2 < 4; ++v2) {
            weight_v2_[v2] = 0; // init clear
        } // v2

        auto constexpr D_Q_ = D_Q_code(D, Q);

        std::printf("\n# D%dQ%d\n", D, Q);
        
        int not_implemented{0};
        // Gauss-Lobatto quadrature weights in 1D are {1/6., 4/6., 1/6.}

        switch ( D_Q_ ) {
            case D_Q_code(3,27):
                weight_v2_[0] = 64/216.; // weight for the center (1x), resting populations
                weight_v2_[1] = 16/216.; // weights for +/-1 in 1 Cartesian direction (6x)
                weight_v2_[2] =  4/216.; // weights for +/-1 in 2 Cartesian directions (12x)
                weight_v2_[3] =  1/216.; // weights for +/-1 in 3 Cartesian directions (8x)
                not_implemented = 1; // ToDo: adjust equilibrium function
            break;
            case D_Q_code(3,19):
                weight_v2_[0] = 72/216.; // weight for the center (1x), resting populations
                weight_v2_[1] = 12/216.; // weights for +/-1 in 1 Cartesian direction (6x)
                weight_v2_[2] =  6/216.; // weights for +/-1 in 2 Cartesian directions (12x)
                // ok, equilibrium function is taken from openLBMflow
            break;
            case D_Q_code(3,15):
                weight_v2_[0] = 48/216.; // weight for the center (1x), resting populations
                weight_v2_[1] = 24/216.; // weights for +/-1 in 1 Cartesian direction (6x)
                weight_v2_[3] =  3/216.; // weights for +/-1 in 3 Cartesian directions (8x)
                not_implemented = 1; // ToDo: adjust equilibrium function
            break;
            case D_Q_code(3, 7): 
                weight_v2_[0] = 120/216.; // weight for the center (1x), resting populations
                weight_v2_[1] =  16/216.; // weights for +/-1 in 1 Cartesian direction (6x)
                not_implemented = 2; // ToDo: check weights (diffusion only), adjust equilibrium function
            break;
            case D_Q_code(2, 9): 
                weight_v2_[0] = 16/36.; // weight for the center (1x), resting populations
                weight_v2_[1] =  4/36.; // weights for +/-1 in 1 Cartesian direction (4x)
                weight_v2_[2] =  1/36.; // weights for +/-1 in 2 Cartesian directions (4x)
                not_implemented = 1; // ToDo: adjust equilibrium function
            break;
            case D_Q_code(2, 5): 
                weight_v2_[0] = 20/36.; // weight for the center (1x), resting populations
                weight_v2_[1] =  4/36.; // weights for +/-1 in 1 Cartesian direction (4x)
                not_implemented = 2; // ToDo: check weights (diffusion only), adjust equilibrium function
            break;
            case D_Q_code(1, 3): 
                weight_v2_[0] = 4/6.; // weight for the center (1x), resting populations
                weight_v2_[1] = 1/6.; // weights for +/-1 in 1 Cartesian direction (2x)
                not_implemented = 1; // ToDo: adjust equilibrium function
            break;
            default:
                std::printf("\n# ERROR: D%dQ%d unknown!\n\n", D, Q);
                assert(0 && "Bad choices of D and Q for the LBM stencil!");
            break;
        } // switch D_Q_

        if (not_implemented > 0) {
            std::printf("# Warning: D%dQ%d needs an adjusted equilibrium function%s!\n",
                                    D, Q, (not_implemented > 1) ? " and weights" : "");
        } // not_implemented

        // check sanity of the quadrature weights
        int const w_class[4] = { 1,   // 1x center weight (resting population q000)
                               2*D,   // pure Cartesian axis, +/-1 per dimension
                       2*D*(D - 1),   // two Cartesian coordinates nonzero
                  (1 << D)*(D > 2) }; // three Cartesian axis --> 2^D combinations
        auto const class_sum = w_class[0] + w_class[1] + w_class[2] + w_class[3];
        auto const three_pow_D = 3*(1 + 2*(D > 1))*(1 + 2*(D > 2));
        assert(three_pow_D == class_sum);
        double const weight_sum = w_class[0]*weight_v2_[0]
                                + w_class[1]*weight_v2_[1]
                                + w_class[2]*weight_v2_[2]
                                + w_class[3]*weight_v2_[3];

        std::printf("# D%dQ%d weight_sum= %.15f, dev= %.1e, class_sum= %d\n",
                       D, Q,  weight_sum,   weight_sum - 1, class_sum);
        assert(std::abs(weight_sum - 1) < 1e-15);

        for(int q = 0; q < 32; ++q) {
//          ex_[q] = 0; ey_[q] = 0; ez_[q] = 0; // init
            for(int d = 0; d < 4; ++d) { // d=3 stores v^2
                velocity_[q][d] = 0; // init
            } // d
            opposite_[q] = 0;
        } // q

        { // scope: set velocity vectors
            int constexpr k12 = (D_Q_code(3,15) == D_Q_
                              || D_Q_code(3, 7) == D_Q_) ? 12 : 0; // special case orderings

    // decode the velocity vector components from the characters of the qxyz constants
    #define set_velocity_vector(         q_xyz) \
            _set_velocity_vector(#q_xyz, q_xyz)

            set_velocity_vector(q_ooo); //         wi0
            set_velocity_vector(q_poo); // +x      wi1
            set_velocity_vector(q_noo); // -x      wi1
            set_velocity_vector(q_opo); //   +y    wi1
            set_velocity_vector(q_ono); //   -y    wi1
          if (0 == k12) {
            set_velocity_vector(q_ppo); // +x+y    wi2
            set_velocity_vector(q_npo); // -x+y    wi2
            set_velocity_vector(q_nno); // -x-y    wi2
            set_velocity_vector(q_pno); // +x-y    wi2
            set_velocity_vector(q_opp); //   +y+z  wi2
            set_velocity_vector(q_onp); //   -y+z  wi2
            set_velocity_vector(q_onn); //   -y-z  wi2
            set_velocity_vector(q_opn); //   +y-z  wi2
            set_velocity_vector(q_pop); // +x  +z  wi2
            set_velocity_vector(q_nop); // -x  +z  wi2
            set_velocity_vector(q_non); // -x  -z  wi2
            set_velocity_vector(q_pon); // +x  -z  wi2
          } // 0 == k12
            set_velocity_vector(q_oop -k12); //     +z  wi1
            set_velocity_vector(q_oon -k12); //     -z  wi1

            set_velocity_vector(q_ppp -k12); // +x+y+z  wi3
            set_velocity_vector(q_npp -k12); // -x+y+z  wi3
            set_velocity_vector(q_pnp -k12); // +x-y+z  wi3
            set_velocity_vector(q_nnp -k12); // -x-y+z  wi3
            set_velocity_vector(q_ppn -k12); // +x+y-z  wi3
            set_velocity_vector(q_npn -k12); // -x+y-z  wi3
            set_velocity_vector(q_pnn -k12); // +x-y-z  wi3
            set_velocity_vector(q_nnn -k12); // -x-y-z  wi3

    #undef  set_velocity_vector
        } // scope

        int8_t cube[3][3][3];
        for(int xyz = 0; xyz < 27; ++xyz) {
            cube[0][0][xyz] = -1; // initialize
        } // xyz in all 27 combinations of +/-1

        int count_velocity2[4] = {0, 0, 0, 0};
        for(int q = 0; q < Q; ++q) {
            // set SoA arrays
//          ex_[q] = velocity_[q][0];
//          ey_[q] = velocity_[q][1];
//          ez_[q] = velocity_[q][2];

            // length^2 of the velocity vector
            auto const *const vel = velocity_[q];
            auto const velocity2 = pow2(vel[0]) + pow2(vel[1]) + pow2(vel[2]);
//          auto const velocity2 = pow2(ex_[q]) + pow2(ey_[q]) + pow2(ez_[q]);
            wi[q] = weight_v2_[velocity2]; // store in private member (error when we do not)
            velocity_[q][3] = velocity2; // store here
            ++count_velocity2[velocity2];

            // register in the velocity cube
//          cube[1+ex_[q]][1+ey_[q]][1+ez_[q]] = q;
            cube[1+vel[0]][1+vel[1]][1+vel[2]] = q;
        } // q
 
        { // scope: additional weight sum check
            double w8_sum{0};
            for(int v2 = 0; v2 < 4; ++v2) {
                w8_sum += count_velocity2[v2] * weight_v2_[v2];
                assert(count_velocity2[v2] <= w_class[v2]);
//              std::printf("# D%dQ%d: count= %d w8= %g\n", D, Q, count_velocity2[v2], weight_v2_[v2]);
            } // v2
            std::printf("# D%dQ%d: w8_sum= %g\n", D, Q, w8_sum);
            assert(std::abs(w8_sum - 1) < 1e-15);
        } // scope

        // check opposite index array (self-check)
        for (int q = 0; q < Q; ++q) {
            auto const vel = velocity(q);
            auto & c = cube[1-vel[0]][1-vel[1]][1-vel[2]];
//          if (c < 0) std::printf("# D%dQ%d: q=%i e= %i %i %i  c= %i\n", D, Q, q, ex_[q], ey_[q], ez_[q], c);
            if (c < 0) std::printf("# D%dQ%d: q=%i e= %i %i %i  c= %i\n", D, Q, q, vel[0], vel[1], vel[2], c);
            assert(c >= 0); // make sure we do not hit uninitialized values
            opposite_[q] = c;
            c = -9; // mark used (c is a reference to a cube entry)

            { // scope: additional checks
//                 assert(ex_[q] == velocity_[q][0]);
//                 assert(ey_[q] == velocity_[q][1]);
//                 assert(ez_[q] == velocity_[q][2]);

                auto const q_anti = opposite_[q];
        //         printf("# q = %2i  velocity vector%3i%3i%3i opposite%3i%3i%3i\n",
        //                q, ex_[q], ey_[q], ez_[q],   ex_[q_anti], ey_[q_anti], ez_[q_anti] );
                for(int d = 0; d < 3; ++d) {
                    assert(velocity_[q][d] == -velocity_[q_anti][d]);
                } // d
//                 assert(ex_[q] == -ex_[q_anti]);
//                 assert(ey_[q] == -ey_[q_anti]);
//                 assert(ez_[q] == -ez_[q_anti]);
            } // scope
        } // q
        
        // count if all opposites have been hit
        {
            int nhit{0};
            for(int xyz = 0; xyz < 27; ++xyz) {
                nhit += (-9 == cube[0][0][xyz]);
            } // xyz
            assert(Q == nhit);
        }

    } // default constructor
    
    template <typename real_t>
    inline real_t equilibrium(
          real_t half_weight // half weight
        , real_t uv // vec u dot vector v
        , real_t u2 // vec u^2
    ) const {
        return half_weight*(2 + 6*uv - 3*u2 + 9*uv*uv);
    } // equilibrium

private:
  
    void _set_velocity_vector(char const *q_xyz, int const q) {
        assert(0 <= q); assert(q < 27);
        assert(nullptr != q_xyz);
        assert('q' == q_xyz[0] && '_' == q_xyz[1]); // must start from "q_"
        for(int d = 0; d < 3; ++d) {
            char const c = q_xyz[d + 2]; // read chars after leading "q_"
            assert('\0' != c); // the name must have at least 4 non-null chars
            auto const c2v = int(c | 32) - 'o'; // ignore case and subtract center
            // now modify AoS member variable
            velocity_[q][d] = c2v;
        } // d
        char const t = q_xyz[5]; // tail char
        assert('\0' == t || '_' == t || '-' == t || ' ' == t); // the name should end or contnue with underscore
    } // _set_velocity_vector

public:
  
    static bool constexpr is_D3Q15() { return (D == 3) && (Q == 15); }
    static bool constexpr is_D3Q7()  { return (D == 3) && (Q ==  7); }

private:
    // velocity directions
//     int8_t ex_[32], ey_[32], ez_[32]; // SoA layout
    int8_t velocity_[Q][4]; // same velocities, but transposed to AoS, hold v2 in addition

    // Half-Way bounce back for all non-zero velocities, flip velocities with their opposites
    uint8_t opposite_[32]; // needs to be initialized
    double wi[27]; // weights as ordered in the populations, not needed, but errors when eliminated
    double weight_v2_[4]; // = {wi0, wi1, wi2, wi3}; // weights as a function of velocity^2

public:
    int opposite(int const q) const { assert(0 <= q); assert(q < Q); return opposite_[q]; }
    double weight(int const v2) const { assert(0 <= v2); assert(v2 < 4); return weight_v2_[v2]; }
    int8_t const* velocity(int const q) const { assert(0 <= q); assert(q < Q); return velocity_[q]; }
    uint16_t velocity_squared(int const q) const { assert(0 <= q); assert(q < Q); return velocity_[q][3]; }

//     deprecated
//     int8_t ex(int const q) const { assert(0 <= q); assert(q < Q); return ex_[q]; }
//     int8_t ey(int const q) const { assert(0 <= q); assert(q < Q); return ey_[q]; }
//     int8_t ez(int const q) const { assert(0 <= q); assert(q < Q); return ez_[q]; }

}; // class BKG_stencil


typedef int status_t;

namespace lbm_stencil {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_construction(int const echo=0) {
      BKG_stencil<1,3> d1q3;
      BKG_stencil<2,5> d2q5;
      BKG_stencil<2,9> d2q9;
      BKG_stencil<3,15> d3q15;
      BKG_stencil<3,19> d3q19;
      BKG_stencil<3,27> d3q27;
      return 0;
  } // test_construction

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_construction(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_stencil
