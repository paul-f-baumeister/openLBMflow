#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf
#include <cstdint> // uint8_t, int8_t
#include <cmath> // std::abs

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace lbm_initialize {

  void initialize_body_force(
        double body_force_xyz[3] // result: body force
      , double const body_force
      , double const body_force_dir=0 // 0: negative y-direction
  ) {
      // set body_force vector
      double const arg = body_force_dir*(M_PI/180.);
      body_force_xyz[0] =  body_force*std::sin(arg);
      body_force_xyz[1] = -body_force*std::cos(arg);
      body_force_xyz[2] =  0;

  } // initialize_body_force

  template <typename real_t, class Stencil>
  void initialize_distrFunc(
        view2D<real_t> & populations // result: mover populations(xyz,q)
      , Stencil const & stencil
      , int const NzNyNx
      , CellInfo const cell[]
      , double const *restrict const rho
      , double const *restrict const ux
      , double const *restrict const uy
      , double const *restrict const uz
  ) {

      for (index_t xyz = 0; xyz < NzNyNx; ++xyz) {
          if (cell[xyz].is_liquid()) {
              double const rho_tmp = rho[xyz];
              double const ux_tmp = ux[xyz];
              double const uy_tmp = uy[xyz]; // here, a more complex flow field could be initialized
              double const uz_tmp = uz[xyz];

              double const u_squared = pow2(ux_tmp) + pow2(uy_tmp) + pow2(uz_tmp);
              for (int q = 0; q < stencil.Q; q++) {
                  auto const e = stencil.velocity(q);
                  double const vel = e[0]*ux_tmp + e[1]*uy_tmp + e[2]*uz_tmp;
                  auto const v2 = stencil.velocity_squared(q);
                  assert(v2 == pow2(e[0]) + pow2(e[1]) + pow2(e[2]));
                  double const weight = stencil.weight(v2);
  //              double const feq = 0.5*weight*rho_tmp*(2 + 6*vel + 9*pow2(vel) - 3*u_squared);
                  double const feq = stencil.equilibrium(0.5*weight*rho_tmp, vel, u_squared);
                  populations(xyz, q) = feq;
              } // q

          } // is_liquid
      } // xyz

  } // initialize_distrFunc
  
  
  template <typename real_t>
  void initialize_boundary(
        CellInfo cell[]
      , real_t *restrict const rho
      , real_t *restrict const ux
      , real_t *restrict const uy
      , real_t *restrict const uz
      , int    const boundary[3][2] // {{lef,rig},{bot,top},{fro,bac}}
      , double const rho_solid
      , int const Nx, int const Ny, int const Nz
      , double const wall_speed[3][2] // {{lef,rig},{bot,top},{fro,bac}}
  ) {

      // initialize type of cells
      for (index_t xyz = 0; xyz < Nz*Ny*Nx; ++xyz) {
          assert(cell[xyz].is_liquid());
      } // xyz

      // define Half way Bounce Back Boundary condition
      for (int lu = 0; lu <= 1; ++lu) { // {0:lower, 1:upper}
        
          {   int constexpr dir = 0; // x-direction (left/right)
              if (boundary[dir][lu] > 0) {
                  int const x = lu*(Nx - 1);
                  for (int z = 0; z < Nz; ++z) {
                      for (int y = 0; y < Ny; ++y) {
                          index_t const xyz = indexyz(x, y, z, Nx, Ny); // node on left/right boundary, x==min/max
                          cell[xyz].set_solid();
                          rho[xyz] = rho_solid;
                      } // y
                  } // z
              } // boundary
          } // x-direction

          {   int constexpr dir = 1; // y-direction (down/up)
              if (boundary[dir][lu] > 0) {
                  int const y = lu*(Ny - 1);
                  for (int z = 0; z < Nz; ++z) {
                      for (int x = 0; x < Nx; ++x) {
                          index_t const xyz = indexyz(x, y, z, Nx, Ny); // node on bottom/top boundary, y==min/max
                          cell[xyz].set_solid();
                          rho[xyz] = rho_solid;
                          // specialty of the y-direction: for the lid-driven cavity flow
                          if (0 != wall_speed[dir][lu]) {
                              ux[xyz] = wall_speed[dir][lu];
                              uy[xyz] = 0;
                              uz[xyz] = 0;
//                               std::printf("# lbm: set wall speed [%g %g %g] on xyz=(%i %i %i)\n",
//                                                             ux[xyz], uy[xyz], uz[xyz], x, y, z);
                          } // wall_speed flowing in x-direction
                      } // x
                  } // z
              } // boundary
          } // y-direction

          {   int constexpr dir = 2; // z-direction (front/back)
              if (boundary[dir][lu] > 0) {
                  int const z = lu*(Nz - 1);
                  for (int y = 0; y < Ny; ++y) {
                      for (int x = 0; x < Nx; ++x) {
                          index_t const xyz = indexyz(x, y, z, Nx, Ny); // node on front/back boundary, z==min/max
                          cell[xyz].set_solid();
                          rho[xyz] = rho_solid;
                      } // x
                  } // y
              } // boundary
          } // z-direction

      } // lu

  } // initialize_boundary

  void initialize_density(
        double *restrict const rho // result
      , CellInfo const cell[]
      , double const value
      , int const NzNyNx
  ) {
      for (index_t xyz = 0; xyz < NzNyNx; ++xyz) {
          if (cell[xyz].is_liquid()) rho[xyz] = value;
      } // xyz
  } // initialize_density

  void initialize_droplet(
        double     *restrict const rho // in/output density
      , CellInfo const cell[]
      , int const Nx, int const Ny, int const Nz
      , double const dx, double const dy, double const dz // droplet center
      , double const dr // droplet radius
      , int    const drop
      , double const ifaceW
      , double const rho_high
      , double const rho_low
  ) {
      auto const rho_diff = (rho_high - rho_low)*drop;
      auto const rho_sum  =  rho_high + rho_low;
      for (int z = 0; z < Nz; ++z) {
          for (int y = 0; y < Ny; ++y) {
              for (int x = 0; x < Nx; ++x) {
                  index_t const xyz = indexyz(x, y, z, Nx, Ny);
                  if (cell[xyz].is_liquid()) {
                      auto const radius = std::sqrt(pow2(x - dx) + pow2(y - dy) + pow2(z - dz));
                      auto const tmp_rho = 0.5*(rho_sum - rho_diff*std::tanh(2.0*(radius - dr)/ifaceW));
                      if (tmp_rho > rho[xyz]) rho[xyz] = tmp_rho;
                  } // is_liquid
              } // x
          } // y
      } // z

  } // initialize_droplet

  
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_construction(int const echo=0) {
      return 0;
  } // test_construction

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_construction(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_initialize
