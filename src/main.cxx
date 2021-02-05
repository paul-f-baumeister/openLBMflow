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
#include <vector> // std::vector<T>

typedef size_t index_t; // defines the integer data type for direct indexing

#include "get_memory.hxx" // get_memory<T>

#include "warnings.hxx" // warn, show_warnings

#include "lbm_stencil.hxx" // BKG_stencil<D,Q>, qooo, q..., pow2
#include "lbm_cell.hxx" // CellInfo

#include "data_view.hxx" // view2D<T>

inline index_t indexyz(int const x, int const y, int const z, int const Nx, int const Ny) { return (z*Ny + y)*Nx + x; }
#define phindex(x,y,z) ((x+1)*(Ny+2) + (y+1))*(Nz+2) + (z+1) // enlarged with a halo of thickness 1

#define restrict __restrict__

#include "lbm_initialize.hxx" // ::initialize_body_force, initialize_distrFunc
#include "lbm_visualize.hxx" // ::write_collection_pvd, ::writeVTK, ::outputSave


template <typename real_t>
void transfer_phi_halos(
      real_t *restrict const phi // modified
    , int const Nx
    , int const Ny
    , int const Nz
) 
    // only used in multiphase case
{
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

template <typename real_t>
inline void solid_cell_treatment( // ToDo: reorder argument list
      real_t        *restrict const tmp_fn // inout: tmp_fn[Q]
    , view2D<real_t> const & fn  // fn[nxyz*Q] or fn[Q*nxyz]
    , index_t const xyz
    , double  const *restrict const rho // rho[nxyz]
    , uint8_t const *restrict const opposite
    , int     const Q
    , double  const top_wall_speed
    , double  const bot_wall_speed=0
) {

    for (int q = 0; q < Q; ++q) {
        tmp_fn[q] = fn(xyz, opposite[q]); // reflection
    } // q

    double constexpr sixth = 1/6.;
    if (0 != top_wall_speed) { // moving top wall (y=Ny-1, speed in x-direction)
//      singlephase_couette_flow:      top_wall_speed == 0.5
//      singlephase_lid_driven_cavity: top_wall_speed == 0.5
        tmp_fn[q_nno] -= top_wall_speed*sixth*rho[xyz];
        tmp_fn[q_pno] += top_wall_speed*sixth/rho[xyz];
    } // top wall speed
    if (0 != bot_wall_speed) { // moving bottom wall (y=0, speed in x-direction)
        tmp_fn[q_ppo] += bot_wall_speed*sixth*rho[xyz];
        tmp_fn[q_npo] -= bot_wall_speed*sixth/rho[xyz];
    } // bottom wall speed

} // solid_cell_treatment


// propagate-kernel (performance critical code section)
template <typename real_t>
inline void propagate(
      view2D<real_t> & fn // output f(xyz, q)
    , real_t const *restrict const tmp_fn // input vector[Q]
    , int const  x, int const  y, int const  z
    , int const Nx, int const Ny, int const Nz
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
    fn(indexyz(ox, oy, oz, Nx, Ny), q_ooo) = tmp_fn[q_ooo]; //
    fn(indexyz(px, oy, oz, Nx, Ny), q_poo) = tmp_fn[q_poo]; // +x
    fn(indexyz(px, py, oz, Nx, Ny), q_ppo) = tmp_fn[q_ppo]; // +x+y
    fn(indexyz(ox, py, oz, Nx, Ny), q_opo) = tmp_fn[q_opo]; //   +y
    fn(indexyz(nx, py, oz, Nx, Ny), q_npo) = tmp_fn[q_npo]; // -x+y
    fn(indexyz(nx, oy, oz, Nx, Ny), q_noo) = tmp_fn[q_noo]; // -x
    fn(indexyz(nx, ny, oz, Nx, Ny), q_nno) = tmp_fn[q_nno]; // -x-y
    fn(indexyz(ox, ny, oz, Nx, Ny), q_ono) = tmp_fn[q_ono]; //   -y
    fn(indexyz(px, ny, oz, Nx, Ny), q_pno) = tmp_fn[q_pno]; // +x-y
    fn(indexyz(ox, py, pz, Nx, Ny), q_opp) = tmp_fn[q_opp]; //   +y+z
    fn(indexyz(ox, oy, pz, Nx, Ny), q_oop) = tmp_fn[q_oop]; //     +z
    fn(indexyz(ox, ny, pz, Nx, Ny), q_onp) = tmp_fn[q_onp]; //   -y+z
    fn(indexyz(ox, ny, nz, Nx, Ny), q_onn) = tmp_fn[q_onn]; //   -y-z
    fn(indexyz(ox, oy, nz, Nx, Ny), q_oon) = tmp_fn[q_oon]; //     -z
    fn(indexyz(ox, py, nz, Nx, Ny), q_opn) = tmp_fn[q_opn]; //   +y-z
    fn(indexyz(px, oy, pz, Nx, Ny), q_pop) = tmp_fn[q_pop]; // +x  +z
    fn(indexyz(nx, oy, pz, Nx, Ny), q_nop) = tmp_fn[q_nop]; // -x  +z
    fn(indexyz(nx, oy, nz, Nx, Ny), q_non) = tmp_fn[q_non]; // -x  -z
    fn(indexyz(px, oy, nz, Nx, Ny), q_pon) = tmp_fn[q_pon]; // +x  -z

} // propagate


// kernel (performance critical code section)
template <typename real_t, bool multiphase>
void update(
      view2D<real_t> & fn // output populations fn(xyz, q)
    , double       *restrict const rho // output density
    , double       *restrict const ux
    , double       *restrict const uy // output macroscopic velocities
    , double       *restrict const uz
    , view2D<real_t> const & f_previous // input previous populations fp(xyz, q)
    , BKG_stencil<3,19> const & stencil
    , double const tau
    , int const Nx, int const Ny, int const Nz
    , CellInfo const cell[]
    , double const body_force[3]
    , double const top_wall_speed
    , double const bot_wall_speed
    , double       *restrict const phi=nullptr // output phase field in the case of multiphase flow
    , double const G=0 // interparticular interaction potential
) {
    int constexpr Q = 19; assert(Q == stencil.Q);
    // relaxation time constant tau
    double const inv_tau = 1.0/tau;
    double const min_tau = 1.0 - inv_tau;
    double const half_wi0 = 0.5*stencil.weight(0), 
                 half_wi1 = 0.5*stencil.weight(1),
                 half_wi2 = 0.5*stencil.weight(2);


    if (multiphase) {
        // update rho to calculate phi
        assert(nullptr != phi);

        for (int z = 0; z < Nz; ++z) {
            for (int y = 0; y < Ny; ++y) {
                for (int x = 0; x < Nx; ++x) {
                    index_t const xyz = indexyz(x, y, z, Nx, Ny);
                    if (cell[xyz].is_liquid()) {

                        auto const fp = f_previous[xyz]; // get a 1D subview
                        // calculate only rho
// #define ALLOW_DEVIATIONS
#ifndef ALLOW_DEVIATIONS
                        double const tmp_rho = double(fp[q_ooo])
                                        + (double(fp[q_poo]) + double(fp[q_noo])) 
                                        + (double(fp[q_opo]) + double(fp[q_ono]))
                                        + (double(fp[q_oop]) + double(fp[q_oon]))
                                        
                                        + (double(fp[q_opp]) + double(fp[q_onn]))
                                        + (double(fp[q_onp]) + double(fp[q_opn]))
                                        + (double(fp[q_pop]) + double(fp[q_non]))
                                        + (double(fp[q_nop]) + double(fp[q_pon]))
                                        + (double(fp[q_ppo]) + double(fp[q_nno]))
                                        + (double(fp[q_npo]) + double(fp[q_pno]));
#else  // ALLOW_DEVIATIONS
                        // if we allow some deviations, it can be written much shorter 
                        //          and much more flexible w.r.t. different stencils:
                        double tmp_rho{0};
                        for(int q = 0; q < Q; ++q) {
                            tmp_rho += double(fp[q]);
                        } // q
#endif // ALLOW_DEVIATIONS
                        rho[xyz] = tmp_rho; // store density
                    } // is_liquid
                    // calculate interparticular force in multiphase Shan-Chen model
                    phi[phindex(x, y, z)] = 1 - std::exp(-rho[xyz]); 
                } // x
            } // y
        } // z

        transfer_phi_halos(phi, Nx, Ny, Nz); // make the halo-enlarged array periodic

    } // multiphase

    for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0; x < Nx; ++x) {
                index_t const xyz = indexyz(x, y, z, Nx, Ny);
                if (cell[xyz].is_liquid()) {
  
                    double G_phi_grad_phi[3];
                    if (multiphase) {
                        double constexpr inv_w2 = 1/36., inv_w1 = 2/36.; // weights for D3Q19

                        double grad_phi[3];
                        // calculate phi-gradients
                        #define ph(X,Y,Z) phi[phindex((X), (Y), (Z))]
                        grad_phi[0] = (ph(x+1, y, z) - ph(x-1, y, z))*inv_w1;
                        grad_phi[1] = (ph(x, y+1, z) - ph(x, y-1, z))*inv_w1;
                        grad_phi[2] = (ph(x, y, z+1) - ph(x, y, z-1))*inv_w1;

                        // in the next three sections, every phi value is used twice
                        double const ph_ppo = ph(x+1, y+1, z);
                        double const ph_npo = ph(x-1, y+1, z);
                        double const ph_pno = ph(x+1, y-1, z);
                        double const ph_nno = ph(x-1, y-1, z);
                        grad_phi[0] += (ph_ppo - ph_npo + ph_pno - ph_nno)*inv_w2;
                        grad_phi[1] += (ph_ppo + ph_npo - ph_nno - ph_pno)*inv_w2;

                        double const ph_pop = ph(x+1, y, z+1);
                        double const ph_nop = ph(x-1, y, z+1);
                        double const ph_pon = ph(x+1, y, z-1);
                        double const ph_non = ph(x-1, y, z-1);
                        grad_phi[2] += (ph_pop + ph_nop - ph_non - ph_pon)*inv_w2;
                        grad_phi[0] += (ph_pop - ph_nop + ph_pon - ph_non)*inv_w2;
                        
                        double const ph_opp = ph(x, y+1, z+1);
                        double const ph_onp = ph(x, y-1, z+1);
                        double const ph_opn = ph(x, y+1, z-1);
                        double const ph_onn = ph(x, y-1, z-1);
                        grad_phi[1] += (ph_opp + ph_opn - ph_onp - ph_onn)*inv_w2;
                        grad_phi[2] += (ph_opp + ph_onp - ph_onn - ph_opn)*inv_w2;

                        // central point
                        double const G_phi = G*ph(x, y, z);
                        #undef ph // abbreviation

                        for(int d = 0; d < 3; ++d) {
                            G_phi_grad_phi[d] = G_phi*grad_phi[d];
                        } // d

                    } // multiphase

                    auto const fp = f_previous[xyz]; // get a 1D subview
                    // load
                    double const f_ooo = fp[q_ooo];
                    double const f_poo = fp[q_poo];
                    double const f_ppo = fp[q_ppo];
                    double const f_opo = fp[q_opo];
                    double const f_npo = fp[q_npo];
                    double const f_noo = fp[q_noo];
                    double const f_nno = fp[q_nno];
                    double const f_ono = fp[q_ono];
                    double const f_pno = fp[q_pno];
                    double const f_opp = fp[q_opp];
                    double const f_oop = fp[q_oop];
                    double const f_onp = fp[q_onp];
                    double const f_onn = fp[q_onn];
                    double const f_oon = fp[q_oon];
                    double const f_opn = fp[q_opn];
                    double const f_pop = fp[q_pop];
                    double const f_nop = fp[q_nop];
                    double const f_non = fp[q_non];
                    double const f_pon = fp[q_pon];
                    
                    // calculate rho and ux, uy, uz
                    auto const tmp_rho = f_ooo 
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

                    //                p:positive, n:negative, o:origin
                    // x-current         p       n      
                    double tmp_ux = ( (f_poo - f_noo)
                                    + (f_ppo - f_nno) 
                                    + (f_pno - f_npo)
                                    + (f_pop - f_non)
                                    + (f_pon - f_nop) )*inv_rho;

                    // y-current          p       n
                    double tmp_uy = ( (f_opo - f_ono) 
                                    + (f_ppo - f_nno) 
                                    + (f_npo - f_pno)
                                    + (f_opp - f_onn)
                                    + (f_opn - f_onp) )*inv_rho;

                    // z-current           p       n
                    double tmp_uz = ( (f_oop - f_oon) 
                                    + (f_pop - f_non) 
                                    + (f_nop - f_pon)
                                    + (f_opp - f_onn)
                                    + (f_onp - f_opn) )*inv_rho;
                                
                    double force[3] = {body_force[0], body_force[1], body_force[2]};
                    if (multiphase) {

        #ifndef ALLOW_DEVIATIONS
                        // make sure the previous triple loop has produced the same density
                        assert(tmp_rho == rho[xyz]);
        #endif // ALLOW_DEVIATIONS

                        // interparticular potential in equilibrium velocity
                        for(int d = 0; d < 3; ++d) {
                            force[d] -= G_phi_grad_phi[d]*inv_rho;
                        } // d
                    } else {
                        rho[xyz] = tmp_rho; // store the density
                    } // multiphase

                    // add the body force and phi-gradients in multiphase
                    tmp_ux += tau*force[0];
                    tmp_uy += tau*force[1];
                    tmp_uz += tau*force[2];

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
                    propagate(fn, tmp_fn, x, y, z, Nx, Ny, Nz); // writes into fn, also propagates into solid boundary cells
                } // is_liquid
            } // x
        } // y
    } // z

    // we have to treat the solid cells since some velocities may have penetrated them
    for (int z = 0; z < Nz; ++z) {
        for (int y = 0; y < Ny; ++y) {
            for (int x = 0 ; x < Nx; ++x) {
                index_t const xyz = indexyz(x, y, z, Nx, Ny);
                if (cell[xyz].is_liquid()) {
                    // invoke propagate here, if propagate is separated from collide: copy fp into tmp_fp and call propagate(x, y, z, tmp_fp, fn)
                } else {
                    real_t tmp_fn[Q];
                    // reflect velocities by swapping the corresponding populations
                    solid_cell_treatment(tmp_fn, fn, xyz, rho, stencil.opposite(), stencil.Q, top_wall_speed, bot_wall_speed);
                    propagate(fn, tmp_fn, x, y, z, Nx, Ny, Nz);
                } // is_liquid
            } // x
        } // y
    } // z

} // update

template <typename real_t> // floating point type of populations
double run(
    int const myrank=0
) {
    assert(2*(1/real_t(2)) == 1); // real_t must be a floating point type

#include "openLBMFlow_conf.c"
//  , nx, ny, nz
//  , boundary_*
//  , rhoh, rhol, rho_boundary, drop...
//  , bot_wall_speed, top_wall_speed
//  , total_time, time_save  
//  , ...
#ifndef Lattice3D
    #error "only 3D version available"
#endif  
    
#ifdef MultiPhase
    bool constexpr multiphase = true;
#else
    bool constexpr multiphase = false;
#endif

    int const Nx = nx;
    int const Ny = ny; 
    int const Nz = nz; // local lattice sizes equal to global
    int const NzNyNx = Nz*Ny*Nx; // here is some potential for memory alignment, but watch out for loop limits
    auto const NzNyNx_aligned = (((NzNyNx - 1) >> 2) + 1) << 2; // align to 4 real_t
    assert(NzNyNx_aligned >= NzNyNx);

    
    BKG_stencil<3,19> const stencil;
    
    int constexpr Q_aligned = stencil.Q + 1; // Q numbers are always odd on regular lattices, so we get better memory alignment by +1

    // allocate memory
    size_t const t_mem = (  sizeof(char)  *1 
                          + sizeof(double)*4
                          + sizeof(real_t)*2*Q_aligned
                          + sizeof(double)*int(multiphase) // does not account for phi-halos
                           ) * nz*ny*nx;
    printf("LBM  needs %.6f GiByte total for D%dQ%d with (%d x %d x %d) cells.\n", 
                t_mem/double(1ull << 30), stencil.D, stencil.Q, nx, ny, nz);

    // cell info
    auto const solid = get_memory<CellInfo>(NzNyNx, 0); // one Byte per cell

    // observables
    // ToDo: group together into a view2D<double> that can be switched between SoA[4][NzNyNx_aligned] and AoS[NzNyNx][4];
    auto const rho   = get_memory<double>(NzNyNx_aligned, 0.0);
    auto const ux    = get_memory<double>(NzNyNx_aligned, 0.0);
    auto const uy    = get_memory<double>(NzNyNx_aligned, 0.0);
    auto const uz    = get_memory<double>(NzNyNx_aligned, 0.0);
    // this array has halo-borders and needs to be indexed using phindex(x,y,z);
    auto const phi   = multiphase ? get_memory<double>((Nx+2l)*(Ny+2l)*(Nz+2l)) : nullptr; 

    int const boundary[3][2] = { {boundary_lef, boundary_rig},
                                 {boundary_bot, boundary_top},
                                 {boundary_fro, boundary_bac} };
    double const wall_speed[3][2] = { {0, 0}, {bot_wall_speed, top_wall_speed}, {0, 0} };

    double const rho_solid = rho_boundary*(rhoh - rhol) + rhol;
    lbm_initialize::initialize_boundary(solid, rho, ux, uy, uz, boundary, rho_solid, Nx, Ny, Nz, wall_speed);

    lbm_initialize::initialize_density(rho, solid, rhol, NzNyNx); // low density value
    
    if (d1r > 0) lbm_initialize::initialize_droplet(rho, solid, Nx, Ny, Nz, d1x, d1y, d1z, d1r, drop1, ifaceW, rhoh, rhol);  // droplet 1
    if (d2r > 0) lbm_initialize::initialize_droplet(rho, solid, Nx, Ny, Nz, d2x, d2y, d2z, d2r, drop2, ifaceW, rhoh, rhol);  // droplet 2
    // extend:   lbm_initialize::initialize_droplet(rho, solid, Nx, Ny, Nz, d3x, d3y, d3z, d3r, drop3, ifaceW, rhoh, rhol);  // droplet 3

    double body_force_xyz[3];
    lbm_initialize::initialize_body_force(body_force_xyz, body_force, body_force_dir);

    
    // populations of the velocity distribution functions (two copies)
//     auto const populations0 = get_memory<real_t>(NzNyNx*Q_aligned);
//     auto const populations1 = get_memory<real_t>(NzNyNx*Q_aligned);
//     ipop(xyz, q): xyz*Q_aligned + q; // activate for AoS (array of structs) data layout
//     view2D<real_t> f0(populations0, Q_aligned); // warp
//     view2D<real_t> f1(populations1, Q_aligned); // warp
    // use memory owning constructors
#define AoS
#ifdef AoS
    view2D<real_t> f0(NzNyNx, Q_aligned); // warp
    view2D<real_t> f1(NzNyNx, Q_aligned); // warp
#else
    view2D<real_t> f0_memory(Q_aligned, NzNyNx); // warp
    view2D<real_t> f1_memory(Q_aligned, NzNyNx); // warp
    auto f0 = f0_memory.transpose();
    auto f1 = f1_memory.transpose();
#endif

    
    
    lbm_initialize::initialize_distrFunc(f0, stencil, NzNyNx, solid, rho, ux, uy, uz);

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
            // save output to VTK image file
            auto const speed = lbm_visualize::outputSave(t, rho, ux, uy, uz, phi, ranks, myrank, Nx, Ny, Nz, nx, ny, nz, save_rho, save_pre, save_vel, G);
            if (speed > 0) {
                speed_stats[0] += 1;
                speed_stats[1] += speed;
                speed_stats[2] += speed*speed;
            } // speed
        } // measure and dump data as .vti files

        // calculate the distribution function for the next two time steps
        update<real_t, multiphase>(f1, rho, ux, uy, uz, f0, stencil, tau, Nx, Ny, Nz, solid, body_force_xyz, top_wall_speed, bot_wall_speed, phi, G);
        update<real_t, multiphase>(f0, rho, ux, uy, uz, f1, stencil, tau, Nx, Ny, Nz, solid, body_force_xyz, top_wall_speed, bot_wall_speed, phi, G);
    } // t

#ifndef SuppressIO
    lbm_visualize::write_collection_pvd(nx, ny, nz, "openLBMflow", "output", time_total, time_save);
#endif    

    // compute mean and deviation    
    auto const denom = ((speed_stats[0] > 0) ? 1/speed_stats[0] : 1.); // inverse denominator
    speed_stats[1] *= denom;
    speed_stats[2] *= denom; speed_stats[2] -= speed_stats[1]*speed_stats[1]; speed_stats[2] = std::abs(speed_stats[2]);
    speed_stats[2] *= ((speed_stats[0] > 1)?speed_stats[0]/(speed_stats[0] - 1.0):1.0); // Gaussian variance
    printf("LBM  avgSpeed %.3f +/- %.3f MLUP/s\n", speed_stats[1], std::sqrt(speed_stats[2]));

    // free allocated memory
    delete[] solid;
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
