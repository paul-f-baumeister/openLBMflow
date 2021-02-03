// openLBMflow v1.0.1 Copyright (C) 2013 LBMflow
// Open Source Lattice Boltzmann Solver
// www.lbmflow.com
// open@lbmflow.com

// LICENSE
// The openLBMflow code is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// DISCLAIMER OF WARRANTY
// The code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// This software may contain errors that could cause failures or loss of data,
// and may be incomplete or contain inaccuracies.  You expressly acknowledge and agree
// that use of the openLBMflow software is at your sole risk.
// The openLBMflow software is provided 'AS IS' and without warranty of any kind.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cassert> // assert
#include <numeric> // std::accumulate

typedef size_t index_t; // defines the integer data type for direct indexing

#include "get_memory.hxx" // get_memory<T>
#include "lbm_visualization.hxx" // write_collection_pvd, writeVTK
#include "warnings.hxx" // warn, show_warnings

#include "lbm_stencil.hxx" // BKG_stencil<D,Q>, qooo, q..., pow2
#include "lbm_cell.hxx" // CellInfo

// === begin global variables

// include file with initial parameters
#include "openLBMFlow_conf.c"

int constexpr D = 3, Q = 19;
BKG_stencil<D,Q> const stencil;

int Nx, Ny, Nz, NxNyNz; // local lattice sizes, NxNyNz >= Nx*Ny*Nz

// === end global variables

#ifndef Lattice3D
    #error "only 3D version available"
#endif  

#define restrict __restrict__

#define wall_clock(noarg) ((double) clock()/((double) CLOCKS_PER_SEC))


inline index_t indexyz(int const x, int const y, int const z) { return (x*Ny + y)*Nz + z; }
inline index_t phindex(int const x, int const y, int const z) { return ((x+1)*(Ny+2) + (y+1))*(Nz+2) + (z+1); } // enlarged with a halo of thickness 1
int constexpr Q_aligned = Q + 1; // Q numbers are always odd on regular lattices, so we get better memory alignment by +1
inline index_t ipop(index_t const xyz, int const q) {
    return xyz*Q_aligned + q; // activate for AoS (array of structs) data layout
//     return q*NxNyNz + xyz; // activate for SoA (structure of arrays) data layout
} // ipop

#ifdef MultiPhase
template <typename real_t>
void transfer_phi_halos(
      real_t *restrict const phi
//  , Nx, Ny, Nz                  from global variables
) {

    int const Np[] = {Nx, Ny, Nz};
    // ============================================================================================
    // these operations may not be reordered or performed in parallel!
    // ============================================================================================
    for (int dir = 0; dir < 3; ++dir) { // ordering dir=0, dir=1, dir=2 may not be changed!
        int const dis = (2 - dir) >> 1; // slow index of the one of the two perpendicular directions
        int const dif = 2 - (dir >> 1); // fast index of the one of the two perpendicular directions
        int const os = (dir + 1) >> 1;  // slow offset that ensure to copy and send only meaningful data
        int const of = dir >> 1;        // fast offset that ensure to copy and send only meaningful data
        int const Nr = Np[dir], Ns = Np[dis], Nf = Np[dif];
        int iup[3], idn[3], upi[3], dni[3];

        { // serial and periodic
            // fill halo-cells with periodic values
                    dni[dir] = Nr; upi[dir] = -1; idn[dir] =  0; iup[dir] = Nr - 1;
            for (int si = -os; si < Ns + os; ++si) {      
                    dni[dis] = si; upi[dis] = si; idn[dis] = si; iup[dis] = si; // slow index
                for (int fi = -of; fi < Nf + of; ++fi) { 
                    dni[dif] = fi; upi[dif] = fi; idn[dif] = fi; iup[dif] = fi; // fast index

                    phi[phindex(dni[0], dni[1], dni[2])] = phi[phindex(idn[0], idn[1], idn[2])];
                    phi[phindex(upi[0], upi[1], upi[2])] = phi[phindex(iup[0], iup[1], iup[2])];
                } // f
            } // s

        } // serial and periodic

    } // dir

} // transfer_phi_halos
#endif // MultiPhase

template <typename real_t>
inline void solid_cell_treatment(
      index_t const xyz
    , double const *restrict const rho // rho[nxyz]
    , real_t const *restrict const fn  // fn[nxyz*Q] or fn[Q*nxyz]
    , real_t       *restrict const tmp_fn // inout: tmp_fn[Q]
//  , Q, stencil, top_wall_speed, bot_wall_speed           from global variables
) {

    for (int q = 0; q < Q; ++q) {
        tmp_fn[q] = fn[ipop(xyz, stencil.opposite(q))]; // reflection
    } // q


    if (0 != top_wall_speed) { // moving top wall (y=Ny-1, speed in x-direction)
//      singlephase_couette_flow:      top_wall_speed == 0.5
//      singlephase_lid_driven_cavity: top_wall_speed == 0.5
        tmp_fn[q_nno] -= top_wall_speed*(rho[xyz]/6.);
        tmp_fn[q_pno] += top_wall_speed/(rho[xyz]*6.);
    } // top wall speed
    if (0 != bot_wall_speed) { // moving bottom wall (y=0, speed in x-direction)
        tmp_fn[q_ppo] += bot_wall_speed*(rho[xyz]/6.);
        tmp_fn[q_npo] -= bot_wall_speed/(rho[xyz]*6.);
    } // bottom wall speed

} // solid_cell_treatment


// propagate-kernel (performance critical code section)
template <typename real_t>
inline void propagate(
      int const x, int const y, int const z
    , real_t const *restrict const tmp_fn // input vector[Q]
    , real_t       *restrict const fn     // output
//  , Nx, Ny, Nz                            from global variables
) {
    // exploit that we can assume the local domain to be periodic
    int const px = (x >= Nx - 1) ? 0 : x + 1;
    int const ox = x;
    int const nx = (x <= 0) ? Nx - 1 : x - 1;
    
    int const py = (y >= Ny - 1) ? 0 : y + 1;
    int const oy = y;
    int const ny = (y <= 0) ? Ny - 1 : y - 1;

    int const pz = (z >= Nz - 1) ? 0 : z + 1;
    int const oz = z;
    int const nz = (z <= 0) ? Nz - 1 : z - 1;

    // depending of their velocities, the polulations in tmp_fn[Q] at (x,y,z)
    // are streamed to the other lattice sites
    fn[ipop(indexyz(ox, oy, oz), q_ooo)] = tmp_fn[q_ooo]; //
    fn[ipop(indexyz(px, oy, oz), q_poo)] = tmp_fn[q_poo]; // +x
    fn[ipop(indexyz(px, py, oz), q_ppo)] = tmp_fn[q_ppo]; // +x+y
    fn[ipop(indexyz(ox, py, oz), q_opo)] = tmp_fn[q_opo]; //   +y
    fn[ipop(indexyz(nx, py, oz), q_npo)] = tmp_fn[q_npo]; // -x+y
    fn[ipop(indexyz(nx, oy, oz), q_noo)] = tmp_fn[q_noo]; // -x
    fn[ipop(indexyz(nx, ny, oz), q_nno)] = tmp_fn[q_nno]; // -x-y
    fn[ipop(indexyz(ox, ny, oz), q_ono)] = tmp_fn[q_ono]; //   -y
    fn[ipop(indexyz(px, ny, oz), q_pno)] = tmp_fn[q_pno]; // +x-y
    fn[ipop(indexyz(ox, py, pz), q_opp)] = tmp_fn[q_opp]; //   +y+z
    fn[ipop(indexyz(ox, oy, pz), q_oop)] = tmp_fn[q_oop]; //     +z
    fn[ipop(indexyz(ox, ny, pz), q_onp)] = tmp_fn[q_onp]; //   -y+z
    fn[ipop(indexyz(ox, ny, nz), q_onn)] = tmp_fn[q_onn]; //   -y-z
    fn[ipop(indexyz(ox, oy, nz), q_oon)] = tmp_fn[q_oon]; //     -z
    fn[ipop(indexyz(ox, py, nz), q_opn)] = tmp_fn[q_opn]; //   +y-z
    fn[ipop(indexyz(px, oy, pz), q_pop)] = tmp_fn[q_pop]; // +x  +z
    fn[ipop(indexyz(nx, oy, pz), q_nop)] = tmp_fn[q_nop]; // -x  +z
    fn[ipop(indexyz(nx, oy, nz), q_non)] = tmp_fn[q_non]; // -x  -z
    fn[ipop(indexyz(px, oy, nz), q_pon)] = tmp_fn[q_pon]; // +x  -z

} // propagate

// template <typename real_t>
// inline real_t equilibrium(
//       real_t half_weight // half weight
//     , real_t uv // vec u dot vector v
//     , real_t u2 // vec u^2
// ) {
//     return half_weight*(2 + 6*uv - 3*u2 + 9*uv*uv);
// } // equilibrium

// kernel (performance critical code section)
template <typename real_t>
void update(
      char   const *restrict const solid
    , double const *restrict const body_force_xyz
    , real_t const *restrict const fp //  input populations
    , real_t       *restrict const fn // output populations
    , double       *restrict const rho // output density
    , double       *restrict const ux
    , double       *restrict const uy // output macroscopic velocities
    , double       *restrict const uz
    , double       *restrict const phi // output phase field in the case of multiphase flow
//  , Nx, Ny, Nz                            from global variables
//  , stencil                               from global variables

) {
    // relaxation time constant tau
    double const inv_tau = 1.0/tau;
    double const min_tau = 1.0 - inv_tau;
    double const half_wi0 = 0.5*stencil.weight(0), 
                 half_wi1 = 0.5*stencil.weight(1),
                 half_wi2 = 0.5*stencil.weight(2);

    // calculate rho, ux, uy, uz and phi
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                index_t const xyz = indexyz(x, y, z);
                if (!solid[xyz]) {
                    // load
                    double const f_ooo = fp[ipop(xyz, q_ooo)];
                    double const f_poo = fp[ipop(xyz, q_poo)];
                    double const f_ppo = fp[ipop(xyz, q_ppo)];
                    double const f_opo = fp[ipop(xyz, q_opo)];
                    double const f_npo = fp[ipop(xyz, q_npo)];
                    double const f_noo = fp[ipop(xyz, q_noo)];
                    double const f_nno = fp[ipop(xyz, q_nno)];
                    double const f_ono = fp[ipop(xyz, q_ono)];
                    double const f_pno = fp[ipop(xyz, q_pno)];
                    double const f_opp = fp[ipop(xyz, q_opp)];
                    double const f_oop = fp[ipop(xyz, q_oop)];
                    double const f_onp = fp[ipop(xyz, q_onp)];
                    double const f_onn = fp[ipop(xyz, q_onn)];
                    double const f_oon = fp[ipop(xyz, q_oon)];
                    double const f_opn = fp[ipop(xyz, q_opn)];
                    double const f_pop = fp[ipop(xyz, q_pop)];
                    double const f_nop = fp[ipop(xyz, q_nop)];
                    double const f_non = fp[ipop(xyz, q_non)];
                    double const f_pon = fp[ipop(xyz, q_pon)];

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
 
                    // add the body force
                    tmp_ux += tau*body_force_xyz[0];
                    tmp_uy += tau*body_force_xyz[1];
                    tmp_uz += tau*body_force_xyz[2];
                    rho[xyz] = tmp_rho; // store density
#ifdef MultiPhase
// ----------------------------------------------------// multi-phase //---------------------------------------begin
                    ux[xyz] = tmp_ux; //
                    uy[xyz] = tmp_uy; // store macroscopic velocities
                    uz[xyz] = tmp_uz; //
                } // solid
                phi[phindex(x, y, z)] = 1 - std::exp(-rho[xyz]); // calculate interparticular force in multiphase Shan-Chen model
            } // z
        } // y
    } // x

    transfer_phi_halos(phi); // make the halo-enlarged array periodic

    double const inv_w2 = 1/36., inv_w1 = 2/36.; // weights for D3Q19

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                index_t const xyz = indexyz(x, y, z);
                if (!solid[xyz]) {

                    double const tmp_rho = rho[xyz]; // load density
                    double const inv_rho = 1/tmp_rho;
          #define ph(X,Y,Z) phi[phindex((X), (Y), (Z))]
                    double const tmp_phi = ph(x, y, z);
                    // calculate phi-gradients
                    double grad_phi_x = (ph(x+1, y, z) - ph(x-1, y, z))*inv_w1;
                    double grad_phi_y = (ph(x, y+1, z) - ph(x, y-1, z))*inv_w1;
                    double grad_phi_z = (ph(x, y, z+1) - ph(x, y, z-1))*inv_w1;
                    // every phi value is used twice, maybe buffer them
                    grad_phi_x += (ph(x+1, y+1, z) - ph(x-1, y+1, z) + ph(x+1, y-1, z) - ph(x-1, y-1, z))*inv_w2;
                    grad_phi_y += (ph(x+1, y+1, z) + ph(x-1, y+1, z) - ph(x-1, y-1, z) - ph(x+1, y-1, z))*inv_w2;
                    grad_phi_z += (ph(x+1, y, z+1) + ph(x-1, y, z+1) - ph(x-1, y, z-1) - ph(x+1, y, z-1))*inv_w2;
                    grad_phi_x += (ph(x+1, y, z+1) - ph(x-1, y, z+1) + ph(x+1, y, z-1) - ph(x-1, y, z-1))*inv_w2;
                    grad_phi_y += (ph(x, y+1, z+1) + ph(x, y+1, z-1) - ph(x, y-1, z+1) - ph(x, y-1, z-1))*inv_w2;
                    grad_phi_z += (ph(x, y+1, z+1) + ph(x, y-1, z+1) - ph(x, y-1, z-1) - ph(x, y+1, z-1))*inv_w2;
          #undef ph // abbreviation

                    // interparticular potential in equilibrium velocity
                    double const tmp_ux = ux[xyz] - tau*(G*tmp_phi*grad_phi_x)*inv_rho;
                    double const tmp_uy = uy[xyz] - tau*(G*tmp_phi*grad_phi_y)*inv_rho;
                    double const tmp_uz = uz[xyz] - tau*(G*tmp_phi*grad_phi_z)*inv_rho;

                    // load again
                    double const f_ooo = fp[ipop(xyz, q_ooo)];
                    double const f_poo = fp[ipop(xyz, q_poo)];
                    double const f_ppo = fp[ipop(xyz, q_ppo)];
                    double const f_opo = fp[ipop(xyz, q_opo)];
                    double const f_npo = fp[ipop(xyz, q_npo)];
                    double const f_noo = fp[ipop(xyz, q_noo)];
                    double const f_nno = fp[ipop(xyz, q_nno)];
                    double const f_ono = fp[ipop(xyz, q_ono)];
                    double const f_pno = fp[ipop(xyz, q_pno)];
                    double const f_opp = fp[ipop(xyz, q_opp)];
                    double const f_oop = fp[ipop(xyz, q_oop)];
                    double const f_onp = fp[ipop(xyz, q_onp)];
                    double const f_onn = fp[ipop(xyz, q_onn)];
                    double const f_oon = fp[ipop(xyz, q_oon)];
                    double const f_opn = fp[ipop(xyz, q_opn)];
                    double const f_pop = fp[ipop(xyz, q_pop)];
                    double const f_nop = fp[ipop(xyz, q_nop)];
                    double const f_non = fp[ipop(xyz, q_non)];
                    double const f_pon = fp[ipop(xyz, q_pon)];

// ----------------------------------------------------// multi-phase //---------------------------------------end
#endif // MultiPhase

                    ux[xyz] = tmp_ux;
                    uy[xyz] = tmp_uy; // store current directions
                    uz[xyz] = tmp_uz;
                    auto const ux2 = pow2(tmp_ux);
                    auto const uy2 = pow2(tmp_uy);
                    auto const uz2 = pow2(tmp_uz);
                    auto const uxyz2 = ux2 + uy2 + uz2;

                    real_t tmp_fn[Q];
                    
                    // equilibrium(hw, uv, u2) = hw*(2 + 6*uv - 3*u2 + 9*uv*uv);
                    
                    auto const tmp_rho_inv_tau = tmp_rho*inv_tau;
                    auto const w0 = half_wi0*tmp_rho_inv_tau;
                    tmp_fn[q_ooo] = f_ooo*min_tau + w0*(2 - 3*uxyz2);

                    auto const w1 = half_wi1*tmp_rho_inv_tau;
                    tmp_fn[q_poo] = f_poo*min_tau + w1*(2 + 6*tmp_ux + 9*ux2 - 3*uxyz2);
                    tmp_fn[q_noo] = f_noo*min_tau + w1*(2 - 6*tmp_ux + 9*ux2 - 3*uxyz2);
                    tmp_fn[q_opo] = f_opo*min_tau + w1*(2 + 6*tmp_uy + 9*uy2 - 3*uxyz2);
                    tmp_fn[q_ono] = f_ono*min_tau + w1*(2 - 6*tmp_uy + 9*uy2 - 3*uxyz2);
                    tmp_fn[q_oop] = f_oop*min_tau + w1*(2 + 6*tmp_uz + 9*uz2 - 3*uxyz2);
                    tmp_fn[q_oon] = f_oon*min_tau + w1*(2 - 6*tmp_uz + 9*uz2 - 3*uxyz2);

                    auto const uxy2 = ux2 + uy2;
                    auto const uyz2 = uy2 + uz2;
                    auto const uzx2 = uz2 + ux2;
                    auto const uxy = 2*tmp_ux*tmp_uy;
                    auto const uyz = 2*tmp_uy*tmp_uz;
                    auto const uzx = 2*tmp_uz*tmp_ux;
                            
                    auto const w2 = half_wi2*tmp_rho_inv_tau;
                    tmp_fn[q_ppo] = f_ppo*min_tau + w2*(2 + 6*(+tmp_ux + tmp_uy) + 9*(uxy2 + uxy) - 3*uxyz2);
                    tmp_fn[q_npo] = f_npo*min_tau + w2*(2 + 6*(-tmp_ux + tmp_uy) + 9*(uxy2 - uxy) - 3*uxyz2);
                    tmp_fn[q_nno] = f_nno*min_tau + w2*(2 + 6*(-tmp_ux - tmp_uy) + 9*(uxy2 + uxy) - 3*uxyz2);
                    tmp_fn[q_pno] = f_pno*min_tau + w2*(2 + 6*(+tmp_ux - tmp_uy) + 9*(uxy2 - uxy) - 3*uxyz2);

                    tmp_fn[q_opp] = f_opp*min_tau + w2*(2 + 6*(+tmp_uy + tmp_uz) + 9*(uyz2 + uyz) - 3*uxyz2);
                    tmp_fn[q_onp] = f_onp*min_tau + w2*(2 + 6*(-tmp_uy + tmp_uz) + 9*(uyz2 - uyz) - 3*uxyz2);
                    tmp_fn[q_onn] = f_onn*min_tau + w2*(2 + 6*(-tmp_uy - tmp_uz) + 9*(uyz2 + uyz) - 3*uxyz2);
                    tmp_fn[q_opn] = f_opn*min_tau + w2*(2 + 6*(+tmp_uy - tmp_uz) + 9*(uyz2 - uyz) - 3*uxyz2);
                    
                    tmp_fn[q_pop] = f_pop*min_tau + w2*(2 + 6*(+tmp_ux + tmp_uz) + 9*(uzx2 + uzx) - 3*uxyz2);
                    tmp_fn[q_nop] = f_nop*min_tau + w2*(2 + 6*(-tmp_ux + tmp_uz) + 9*(uzx2 - uzx) - 3*uxyz2);
                    tmp_fn[q_non] = f_non*min_tau + w2*(2 + 6*(-tmp_ux - tmp_uz) + 9*(uzx2 + uzx) - 3*uxyz2);
                    tmp_fn[q_pon] = f_pon*min_tau + w2*(2 + 6*(+tmp_ux - tmp_uz) + 9*(uzx2 - uzx) - 3*uxyz2);

                    // the loops for collide and propate are merged, otherwise we would need to store tmp_fn back into fp
                    propagate(x, y, z, tmp_fn, fn); // writes into fn, also propagates into solid boundary cells
                } // solid
            } // z
        } // y
    } // x

    // we have to treat the solid cells since some velocities may have penetrated them
    for (int x = 0 ; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                index_t const xyz = indexyz(x, y, z);
                if (!solid[xyz]) {
                    // invoke propagate here, if propagate is separated from collide: copy fp into tmp_fp and call propagate(x, y, z, tmp_fp, fn)
                } else {
                    real_t tmp_fn[Q];
                    solid_cell_treatment(xyz, rho, fn, tmp_fn); // reflect velocities by swapping the corresponding populations
                    propagate(x, y, z, tmp_fn, fn);
                } // solid
            } // z
        } // y
    } // x

} // update


 //
 //
 // initialization functions rho, ux, uy, uz (not critical for performance)
 //
 //

template <typename real_t>
void initialize_boundary(
      int    const boundary[3][2]
    , double const rho_boundary
    , char   *restrict const solid
    , real_t *restrict const rho
    , real_t *restrict const ux
    , real_t *restrict const uy
    , real_t *restrict const uz
//  , Nx, Ny, Nz                            from global variables
) {

    // initialize type of cells
    for (index_t xyz = 0; xyz < Nx*Ny*Nz; ++xyz) {
        assert(0 == solid[xyz]);
    } // xyz

    double const rho_solid = rho_boundary*(rhoh - rhol) + rhol;
    
    // define Bounce Back Boundary condition
    for (int x = 0; x < Nx; x++) {
        for (int z = 0; z < Nz; z++) {
            if (0 < boundary[1][0]) {
                index_t const xyz = indexyz(x, 0, z); // node on bottom boundary, y==min
                solid[xyz] = 1;
                rho[xyz] = rho_solid;
                if (0 != bot_wall_speed) { ux[xyz] = bot_wall_speed; uy[xyz] = 0; uz[xyz] = 0; }
            }
            if (0 < boundary[1][1]) {
                index_t const xyz = indexyz(x, Ny - 1, z); // node on top boundary, y==max
                solid[xyz] = 2;
                rho[xyz] = rho_solid;
                if (0 != top_wall_speed) { ux[xyz] = top_wall_speed; uy[xyz] = 0; uz[xyz] = 0; }
            }
        } // z
    } // x

    for (int y = 0; y < Ny; y++) {
        for (int z = 0; z < Nz; z++) {
            if (0 < boundary[0][0]) {
                index_t const xyz = indexyz(0, y, z); // node on left boundary, x==min
                solid[xyz] = 3;
                rho[xyz] = rho_solid;
            } 
            if (0 < boundary[0][1]) {
                index_t const xyz = indexyz(Nx - 1, y, z); // node on right boundary, x==max
                solid[xyz] = 4;
                rho[xyz] = rho_solid;
            }
        } // z
    } // y

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            if (0 < boundary[2][0]) {
                index_t const xyz = indexyz(x, y, 0); // node on front boundary, z==min
                solid[xyz] = 5;
                rho[xyz] = rho_solid;
            }
            if (0 < boundary[2][1]) {
                index_t const xyz = indexyz(x, y, Nz - 1); // node on back boundary, z==max
                solid[xyz] = 6;
                rho[xyz] = rho_solid;
            }
        } // y
    } // x

} // initialize_boundary

void initialize_density(
      double *restrict const rho // result
    , char const *restrict const solid
    , double const value
//  , Nx*Ny*Nz             from global variables
) {
    for (index_t xyz = 0; xyz < Nx*Ny*Nz; ++xyz) {
        if (!solid[xyz]) rho[xyz] = value;
    } // xyz
} // initialize_density

template <typename real_t>
void initialize_distrFunc(
      real_t *restrict const f  // f[nx*ny*nz*Q_aligned] mover populations
    , char const *restrict const solid
    , double const *restrict const rho
    , double *restrict const ux
    , double *restrict const uy
    , double *restrict const uz
    , double body_force_xyz[3]
//  , Nx*Ny*Nz             from global variables
) {

    // set body_force vector (global variable body_force_xyz)
    double const arg = body_force_dir*(M_PI/180.);
    body_force_xyz[0] =  body_force*std::sin(arg);
    body_force_xyz[1] = -body_force*std::cos(arg);
    body_force_xyz[2] =  0;

    double const uvec[3] = {0, 0, 0}; // input current

    for (index_t xyz = 0; xyz < Nx*Ny*Nz; ++xyz) {
        if (!solid[xyz]) {
            double const rho_tmp = rho[xyz];
            double const ux_tmp = uvec[0];
            double const uy_tmp = uvec[1]; // here, a more complex flow field could be initialized
            double const uz_tmp = uvec[2];

            double const u_squared = pow2(ux_tmp) + pow2(uy_tmp) + pow2(uz_tmp);
            for (int q = 0; q < Q; q++) {
                auto const e = stencil.velocity(q);
                double const vel = e[0]*ux_tmp + e[1]*uy_tmp + e[2]*uz_tmp;
//                auto const v2 = e[3];
                auto const v2 = stencil.velocity_squared(q);
                assert(v2 == pow2(e[0]) + pow2(e[1]) + pow2(e[2]));
                double const weight = stencil.weight(v2);
//              double const feq = 0.5*weight*rho_tmp*(2 + 6*vel + 9*pow2(vel) - 3*u_squared);
                double const feq = stencil.equilibrium(0.5*weight*rho_tmp, vel, u_squared);
                f[ipop(xyz, q)] = feq;
            } // q

            ux[xyz] = ux_tmp;
            uy[xyz] = uy_tmp;
            uz[xyz] = uz_tmp;
        } // solid
    } // xyz

} // initialize_distrFunc


void initialize_droplet(
      double const dx, double const dy, double const dz // droplet center
    , double const dr // droplet radius
    , int const drop
    , char const *restrict const solid
    , double     *restrict const rho // in/output density
    , double const rho_high
    , double const rho_low=0
//  , Nx, Ny, Nz                            from global variables
) {
    auto const rho_diff = (rho_high - rho_low)*drop;
    auto const rho_sum  =  rho_high + rho_low;
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                index_t const xyz = indexyz(x, y, z);
                if (!solid[xyz]) {
                    auto const radius = std::sqrt(pow2(x - dx) + pow2(y - dy) + pow2(z - dz));
                    auto const tmp = 0.5*(rho_sum - rho_diff*std::tanh(2.0*(radius - dr)/ifaceW));
                    if (tmp > rho[xyz]) rho[xyz] = tmp;
                } // solid
            } // z
        } // y
    } // x
} // initialize_droplet


template <typename real_t=double>
double outputSave(
      int const time
    , double const rho[]
    , double const ux[]
    , double const uy[]
    , double const uz[]
    , double const phi[]
    , int const ranks[3]
    , int const myrank
//  , Nx, Ny, Nz                            from global variables
//  , nx, ny, nz                            from global variables
) {
    static double timer_start, step_start{0};

    double time_stop = wall_clock(); // stop internal timer

    // calculate performance in units of Mega Lattice Site Updates per second: MLUP/s
    double const Speed = ((nx*ny)*(nz*1e-6)*(time - step_start)/(time_stop - timer_start)); // use global lattice sizes nx,ny,nz
    step_start = time;

    double const mass = std::accumulate(rho, rho + Nx*Ny*Nz, 0.0);

    auto const is_master = (0 == myrank);
    if (is_master) printf("t=%d\tSpeed=%f MLUP/s mass=%f\n", time, Speed, mass);

    int const nall = nx*ny*nz;
    int const offs[3] = {ranks[0]*Nx, ranks[1]*Ny, ranks[2]*Nz};
    auto const rho_all = (save_rho) ? get_memory<real_t>(nall) : nullptr;
    auto const pre_all = (save_pre) ? get_memory<real_t>(nall) : nullptr;
    real_t* vel_all[3] = {nullptr, nullptr, nullptr};
    for(int d = 0; d < D*save_vel; ++d) {
        vel_all[d] = get_memory<real_t>(nall);
    } // d

    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int z = 0; z < Nz; z++) {
                index_t const xyz = indexyz(x, y, z); // local index into rho, ux, uy, uz
                size_t const gxyz = ((x + offs[0])*ny + (y + offs[1]))*nz + (z + offs[2]); // global index
                if (save_rho) rho_all[gxyz] = rho[xyz];
                if (save_pre) pre_all[gxyz] = rho[xyz]/3.0 + ((phi) ? G*pow2(phi[phindex(x, y, z)])/6.0 : 0);
                if (save_vel) {
                    vel_all[0][gxyz] = ux[xyz];
                    vel_all[1][gxyz] = uy[xyz];
                    vel_all[2][gxyz] = uz[xyz];
                } // velocities
            } // z
        } // y
    } // x

    if (is_master) {
#ifndef SuppressIO
        lbm_visualization::writeVTK(time, nx, ny, nz, "output", "openLBMflow", 
                        rho_all, pre_all, vel_all[0], vel_all[1], vel_all[2]);
#else
        std::printf("# SuppressIO for writeVTK\n");
#endif
    } // is_master

    if (rho_all) delete[] rho_all;
    if (pre_all) delete[] pre_all;
    for(int d = 0; d < D*save_vel; ++d) {
        if (vel_all[d]) delete[] vel_all[d];
    } // d

    timer_start = wall_clock(); // start internal timer again
    return Speed;
} // outputSave

template <typename real_t> // floating point type of populations
double run(
    int const myrank=0
//  , Nx, Ny, Nz                            from global variables (write)
//  , NxNyNz                                from global variables (write)
//  , nx, ny, nz                            from global variables (read-only)
//  , boundary_*                            from global variables (read-only)
//  , D, Q, Q_aligned                       from global variables (read-only)
//  , rhoh, rhol, rho_boundary, drop...     from global variables (read-only)
//  , total_time, time_save                 from global variables (read-only)
) {
    assert(2*(1/real_t(2)) == 1); // real_t must be a floating point type

    Nx = nx; Ny = ny; Nz = nz; // local lattice sizes equal to global
    NxNyNz = Nx*Ny*Nz; // here is some potential for memory alignment, but watch out for loop limits
    auto const NxNyNz_aligned = (((NxNyNz - 1) >> 2) + 1) << 2;
    assert(NxNyNz_aligned >= NxNyNz);

    // allocate memory
    size_t const t_mem = (  sizeof(char)  *1 
                          + sizeof(double)*4
                          + sizeof(real_t)*2*Q_aligned
#ifdef MultiPhase
                          + sizeof(double)*1 // does not account for phi-halos
#endif
                           ) * nx*ny*nz;
    printf("LBM  needs %.6f GiByte total for D%dQ%d with (%d x %d x %d) cells.\n", 
                t_mem/double(1ull << 30),    D, Q,        nx,  ny,  nz);

    // cell info
    auto const solid = get_memory<char>  (NxNyNz); // one Byte per cell
    // populations (two copies)
    auto const f0    = get_memory<real_t>(NxNyNz*Q_aligned);
    auto const f1    = get_memory<real_t>(NxNyNz*Q_aligned);

    // observables
    // ToDo: group together into a view2D<double> that can be switched between SoA[4][NxNyNz_aligned] and AoS[NxNyNz][4];
    auto const rho   = get_memory<double>(NxNyNz_aligned, 0.0);
    auto const ux    = get_memory<double>(NxNyNz_aligned);
    auto const uy    = get_memory<double>(NxNyNz_aligned);
    auto const uz    = get_memory<double>(NxNyNz_aligned);
#ifdef MultiPhase
    auto const phi = get_memory<double>((Nx+2l)*(Ny+2l)*(Nz+2l)); // this array has halo-borders and needs to be indexed using phindex(x,y,z)
#else
    double* const phi = nullptr;
#endif

    int const boundary[3][2] = { {boundary_lef, boundary_rig},
                                 {boundary_bot, boundary_top},
                                 {boundary_fro, boundary_bac} };
    initialize_boundary(boundary, rho_boundary, solid, rho, ux, uy, uz);

    initialize_density(rho, solid, rhol); // low density value

    if (d1r > 0) initialize_droplet(d1x, d1y, d1z, d1r, drop1, solid, rho, rhoh, rhol);  // droplet 1
    if (d2r > 0) initialize_droplet(d2x, d2y, d2z, d2r, drop2, solid, rho, rhoh, rhol);  // droplet 2
    // extend:   initialize_droplet(d3x, d3y, d3z, d3r, drop3, solid, rho, rhoh, rhol);  // droplet 3

    double body_force_xyz[3];
    initialize_distrFunc(f0, solid, rho, ux, uy, uz, body_force_xyz);

    double speed_stats[] = {0, 0, 0};

#ifdef TIME_TOTAL
    #define time_total TIME_TOTAL  // control the maximum number of time steps at compile time using -DTIME_TOTAL
#endif

#ifdef TIME_SAVE
    #define time_save  TIME_SAVE   // control the .vti-file output frequency at compile time using -DTIME_SAVE
#endif

    assert(0 == time_save % 2); // if the check point interval is an odd number, we cannot combine two time steps

    int const ranks[] = {0, 0, 0}; // 3D ranks for parallelization

    // main iteration loop
    for (int t = 0; t <= time_total; t+=2) {
        if (0 == t % time_save) {
            auto const speed = outputSave(t, rho, ux, uy, uz, phi, ranks, myrank); // save output to VTK image file
            if (speed > 0) {
                speed_stats[0] += 1;
                speed_stats[1] += speed;
                speed_stats[2] += speed*speed;
            } // speed
        } // measure and dump data as .vti files

        // calculate the distribution function for the next two time steps
        update(solid, body_force_xyz, f0, f1, rho, ux, uy, uz, phi);
        update(solid, body_force_xyz, f1, f0, rho, ux, uy, uz, phi);
    } // t

#ifndef SuppressIO
    lbm_visualization::write_collection_pvd(nx, ny, nz, "openLBMflow", "output", time_total, time_save);
#endif    

    // compute mean and deviation    
    auto const denom = ((speed_stats[0] > 0) ? 1/speed_stats[0] : 1.); // inverse denominator
    speed_stats[1] *= denom;
    speed_stats[2] *= denom; speed_stats[2] -= speed_stats[1]*speed_stats[1]; speed_stats[2] = std::abs(speed_stats[2]);
    speed_stats[2] *= ((speed_stats[0] > 1)?speed_stats[0]/(speed_stats[0] - 1.0):1.0); // Gaussian variance
    printf("LBM  avgSpeed %.3f +/- %.3f MLUP/s\n", speed_stats[1], std::sqrt(speed_stats[2]));

    // free allocated memory
    delete[] solid;
    delete[] f0;
    delete[] f1;
    delete[] rho;
    delete[] ux;
    delete[] uy;
    delete[] uz;
    if (phi) delete[] phi;

    return speed_stats[1];
} // run

int main(int argc, char *argv[]) {

    printf("openLBMflow v2.0.0 (c) 2010 www.lbmflow.com\n");

#ifdef _GIT_KEY
    // stringify the value of a macro, two expansion levels needed
    #define macro2string(a) stringify(a)
    #define stringify(b) #b
    std::printf("# %s git checkout " macro2string(_GIT_KEY) "\n\n", argc ? argv[0] : "");
    #undef  stringify
    #undef  macro2string
#endif // _GIT_KEY

    run<double>();
//     run<float>();

    return 0;
} // main
