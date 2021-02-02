#pragma once

#include <cassert> // assert
#include <cstdio> // std::printf
#include <vector> // std::vector<T>
#include <algorithm> // std::max
#include <cmath> // M_PI

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#include "lbm_droplet.hxx" // Droplet

namespace lbm_examples {

  // assemble all the examples that come with openLBMflow v1.0.1

  struct default_example {
      #include "../source/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = true;
  };

  struct multiphase_coalescence {
      #include "../examples/multiphase coalescence/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = true;
  };

  struct multiphase_coalescence_impact {
      #include "../examples/multiphase coalescence impact/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = true;
  };

  struct multiphase_drop_impact {
      #include "../examples/multiphase drop impact/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = true;
  };

  struct multiphase_drop_on_drop_impact {
      #include "../examples/multiphase drop on drop impact/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = true;
  };

  struct multiphase_rising_bubble {
      #include "../examples/multiphase rising bubble/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = true;
  };

  struct singlephase_couette_flow {
      #include "../examples/singlephase couette flow/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = false;
  };

  struct singlephase_lid_driven_cavity {
      #include "../examples/singlephase lid driven cavity/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = false;
  };

  struct singlephase_poiseuille_flow_3D_channel {
      #include "../examples/singlephase poiseuille flow 3D channel/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = false;
  };

  struct singlephase_poiseuille_flow_plate {
      #include "../examples/singlephase poiseuille flow plate/openLBMFlow_conf.c"
      static bool constexpr is_multiphase = false;
  };

#ifdef  Lattice3D
 #undef Lattice3D
#endif

#ifdef  Float
 #undef Float
#endif

#ifdef  MultiPhase
 #undef MultiPhase
#endif

//   class Droplet {
//   public:
// 
//       Droplet(
//             int const x=0, int const y=0, int const z=0
//           , int const radius=-999
//           , int const drop=0
//           , int const echo=0 // increase echo if also the default constructor should print to stdout
//       )
//           : drop_(drop)
//       {
//           xyzr_[0] = x; xyzr_[1] = y; xyzr_[2] = z; xyzr_[3] = radius;
//           if (echo > 0) {
//               std::printf("# %s at %d %d %d with radius %d, drop=%d\n", __func__,
//                     xyzr_[0], xyzr_[1], xyzr_[2], xyzr_[3], drop_);
//           } // drop==0 (default constructor should not print)
//       } // constructor and default constructor
// 
//       int const* position() const { return xyzr_; }
//       int position(int d) const { assert(0 <= d); assert(d < 4); return xyzr_[d]; }
//       int radius() const { return xyzr_[3]; }
//       int drop() const { return drop_; }
// 
//   private:
//       // members
//       int xyzr_[4];
//       int drop_; // 1:droplet, -1:bubble, 0:not_initialized
// 
//   }; // Droplet
  
  
  class Problem {
  private:
    static int constexpr DENSITY_LOW = 0,
                         DENSITY_BOUNDARY = 1,
                         DENSITY_HIGH = 2;
  public:
    


#define WARNING(CRITICAL_CONDITION) \
            if (CRITICAL_CONDITION) { \
                std::printf("# Warning:  " #CRITICAL_CONDITION "  critical!\n"); \
            }
    
    template <class Example>
    Problem(
          Example const & example
        , char const example_name[]=nullptr
        , int const echo=9
    ) {
        if (echo > 0) std::printf("\n#\n# Example %s\n#\n", example_name);

// // //  examples must define the following quantities:
//     int nx = 30;                //lattice size x
//     int ny = 30;                //lattice size y
//     int nz = 30;                //lattice size z (only for 3D code)
//     double tau = 1.0;           //relaxation time
//     double rhoh = 2.6429;       //high density fluid (rhoh=1 for singlephase)
//     double rhol = 0.0734;       //low  density fluid (rhoh=1 for singlephase)
//     double rho_boundary = 0.2;  //define density of solid walls
//     double ifaceW = 3.0;        //interface width
//     double G = -6.0;            //interparticular interaction potential (G=0 for singlephase)
//     double body_force = 0.0;    //gravity force
//     double body_force_dir = 0;  //gravity direction (0=down 90=right 180=top 270=left)
//     int time_total = 5000;      //total time step
//     int time_save = 200;        //save result to VTK image file (*.vti can be open in Paraview)
//     int save_rho = 0;           //save density  in output file (0=don't save, 1=save)
//     int save_pre = 0;           //save pressure in output file (0=don't save, 1=save)
//     int save_vel = 1;           //save velocity in output file (0=don't save, 1=save)
//     int boundary_bot = 1;       //0=periodic, 1=HBB, set half way bounce back on bottom wall
//     int boundary_top = 1;       //0=periodic, 1=HBB, set half way bounce back on top wall
//     int boundary_lef = 1;       //0=periodic, 1=HBB, set half way bounce back on left wall
//     int boundary_rig = 1;       //0=periodic, 1=HBB, set half way bounce back on right wall
//     int boundary_fro = 0;       //0=periodic, 1=HBB, set half way bounce back on front wall (only for 3D)
//     int boundary_bac = 0;       //0=periodic, 1=HBB, set half way bounce back on back wall (only for 3D)
//     double top_wall_speed = 0.5;//speed of top wall for lid driven cavity
//     double bot_wall_speed = 0.0;//speed of bottom wall
//     int d1r = 0;                //droplet1 radius
//     int d1x = 0;                //droplet1 position x (for multiphase model)
//     int d1y = 0;                //droplet1 position y (for multiphase model)
//     int d1z = 0;                //droplet1 position z (for multiphase model)
//     int drop1 = 1.0;            //1=droplet, -1=buble
//     int d2r = 0;                //droplet2 radius
//     int d2x = 0;                //droplet2 position x (for multiphase model)
//     int d2y = 0;                //droplet2 position y (for multiphase model)
//     int d2z = 0;                //droplet2 position z (for multiphase model)
//     int drop2 = 1.0;            //1=droplet, -1=buble
        
        
        // copy the data into the problem description
        nxyz_[0] = example.nx; assert(example.nx > 0); // or better nz here? ToDo
        nxyz_[1] = example.ny; assert(example.ny > 0); //
        nxyz_[2] = example.nz; assert(example.nz > 0); // or better nx here? ToDo
        if (echo > 0) std::printf("# nx= %d ny= %d nz= %d\n", n(0), n(1), n(2));

        WARNING(example.tau < 1 && "over-relaxation");
        tau_ = example.tau; assert(tau_ > 0);

        WARNING(example.rhoh < example.rhol);
        WARNING(example.rhol < 0);
        WARNING(example.rho_boundary < 0);
        rho_low_boundary_high_[DENSITY_HIGH]     = example.rhoh;
        rho_low_boundary_high_[DENSITY_LOW]      = example.rhol;
        rho_low_boundary_high_[DENSITY_BOUNDARY] = example.rho_boundary;
        
        WARNING(example.ifaceW <= 0); // negatives will turn around the tanh argument
        assert (example.ifaceW != 0);
        ifaceW_ = example.ifaceW;
        
        // interparticular_interaction_potential
        G_ = 0;
        if (example.is_multiphase) {
            G_ = example.G;
        } else {
            WARNING(example.G != 0 && "in singlephase");
        } // multiphase

        
        // init
        for(int d = 0; d < 3; ++d) {
            wall_speed_[d][0] = 0;
            wall_speed_[d][1] = 0;
            body_force_xyz_[d] = 0;
        } // d

        { // scope: set body force / gravity
            double body_force{example.body_force}, body_force_dir{example.body_force_dir}; // non-const temporaries
            if (body_force < 0) {
                WARNING("negative body_force, better use positive and opposite direction, ");
                body_force = -body_force; body_force_dir -= 180;
            } // body_force defined negative
            double const arg = body_force_dir*(M_PI/180.);
            body_force_xyz_[0] =  body_force*std::sin(arg);
            body_force_xyz_[1] = -body_force*std::cos(arg);
        } // scope

        WARNING(example.time_total < 1);
        total_time_steps_ = example.time_total;
        WARNING(example.time_save < 1);
        save_every_time_steps_ = std::max(1, example.time_save);

        // which observables do we want in the VTK files?
        save_ruuup_[0] = (example.save_rho != 0) ? 1 : 0;
        save_ruuup_[4] = (example.save_pre != 0) ? 1 : 0;
        for(int ud = 1; ud < 4; ++ud) { 
            save_ruuup_[ud] = (example.save_vel != 0) ? 1 : 0;
        } // ud
        
        // boundary conditions (0: periodic, 1: half way bounce back)
        boundaries_[0][0] = example.boundary_lef;
        boundaries_[0][1] = example.boundary_rig;
        boundaries_[1][0] = example.boundary_bot;
        boundaries_[1][1] = example.boundary_top;
        boundaries_[2][0] = example.boundary_fro;
        boundaries_[2][1] = example.boundary_bac;

        wall_speed_[1][0] = example.bot_wall_speed;
        wall_speed_[1][1] = example.top_wall_speed;
        
        if (echo > 0) {
            for(int d = 0; d < 3; ++d) {
                std::printf("# boundaries in %c-direction %d %d  wall speeds %g %g\n", 'x' + d,
                    boundaries_[d][0], boundaries_[d][1], wall_speed_[d][0], wall_speed_[d][1]);
            } // d
        } // echo
        
        droplets_.clear();
        if (example.is_multiphase) {
            droplets_.reserve(2);
            droplets_.push_back(Droplet(example.d1x, example.d1y, example.d1z, example.d1r, example.drop1, echo));
            droplets_.push_back(Droplet(example.d2x, example.d2y, example.d2z, example.d2r, example.drop2, echo));
            is_multiphase_ = true;
        } else {
            is_multiphase_ = false;
        } // multiphase

        if (echo > 0) std::printf("# %sphase\n#\n", is_multiphase_ ? "multi" : "single");
    } // constructor
    
    
    // accessors

    int n(int const d) const { assert(0 <= d); assert(d < 3); return nxyz_[d]; }
    int operator [](int const d) const { return n(d); } // [int]
    int operator [](char const dir) const { return n((dir | 32) - 'x'); } // [char]
    template <int d> int n() const { return nxyz_[d]; }
    template <char dir> int n() const { return nxyz_[(dir | 32) - 'x']; }
    uint32_t const* n() const { return nxyz_; }

    double rho_boundary() const { return rho_low_boundary_high_[DENSITY_BOUNDARY]; } // density of solid walls
    double rhol()         const { return rho_low_boundary_high_[DENSITY_LOW ]; } // low
    double rhoh()         const { return rho_low_boundary_high_[DENSITY_HIGH]; } // high

    size_t n_droplets() const { return droplets_.size(); }
    Droplet const & droplet(int const i) const { assert(0 <= i); assert(i < droplets_.size()); return droplets_[i]; }
    double const* body_force_xyz() const { return body_force_xyz_; }
    double interparticular_interaction_potential() const { return G_; }
    double interface_width() const { return ifaceW_; }
    double relaxation_time_parameter() const { return tau_; }

    int  time_total() const { return total_time_steps_; }
    int  time_save()  const { return save_every_time_steps_; }
    int  save_rho()   const { return save_ruuup_[0]; } // density
    int  save_vel()   const { return save_ruuup_[1] + save_ruuup_[2] + save_ruuup_[3]; }
    int  save_pre()   const { return save_ruuup_[4]; } // pressure

    bool is_multiphase() const { return is_multiphase_; }

  private:
    // members
    double    rho_low_boundary_high_[3];
    double    tau_;
    double    ifaceW_;
    double    G_;
    double    body_force_xyz_[3];
    uint32_t  nxyz_[3]; // lattice size in x,y,z (nxyz[2] should be 1 for 2D code) 
    int8_t    boundaries_[3][2]; // 0=periodic, 1=HBB, set half way bounce back
    double    wall_speed_[3][2]; // //speed of e.g. the top wall for lid driven cavity
    int       total_time_steps_;
    int       save_every_time_steps_;
    uint8_t   save_ruuup_[5];
    std::vector<Droplet> droplets_;
    bool      is_multiphase_;
  }; // class Problem



#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_example_setups(int const echo=0) {

      #define make_example(EXAMPLE_STRUCT) \
      {                                    \
          EXAMPLE_STRUCT example;          \
          Problem problem(example, #EXAMPLE_STRUCT, echo); \
      }

      make_example(multiphase_coalescence);
      make_example(multiphase_coalescence_impact);
      make_example(multiphase_drop_impact);
      make_example(multiphase_drop_on_drop_impact);
      make_example(multiphase_rising_bubble);
      make_example(singlephase_couette_flow);
      make_example(singlephase_lid_driven_cavity);
      make_example(singlephase_poiseuille_flow_3D_channel);
      make_example(singlephase_poiseuille_flow_plate);
      make_example(default_example);
   
      #undef make_example
      return 0;
  } // test_example_setups

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_example_setups(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  
  
} // namespace lbm_examples
