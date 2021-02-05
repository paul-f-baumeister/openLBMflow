#pragma once
 //
 //
 // visualization functions (not critical for performance)
 //
 //

#include <cstdio> // std::fprintf, std::snprintf, std::fopen, std::fclose
#include <sys/stat.h> // mkdir

namespace lbm_visualize {

  void write_collection_pvd(
        int const nx, int const ny, int const nz
      , char const *filename
      , char const *directory="."
      , int const time_total=1
      , int const time_save=1
  ) {
      int const dir = mkdir(directory, 0777); // delete the second argument in WIN32
      if (false && (0 == dir)) printf("Error: Cannot create output directory!\n");
      char fName[256];
      std::snprintf(fName, 255, "%s/%s.pvd", directory, filename);
      auto const f = std::fopen(fName,"w");
      std::fprintf(f, "<?xml version=\"1.0\"?>\n"
                        "<!-- openLBMflow v1.0.1, www.lbmflow.com -->\n"
                        "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                        "  <Collection>\n");
      for (int tt = 0; tt <= time_total; tt += time_save) {
          std::fprintf(f, "    <DataSet  timestep=\"%d\" group=\"\" part=\"%d\" file=\"%s_%07d.vti\"/>\n",
                                                    tt,                     0,    filename, tt);
      } // tt
      std::fprintf(f, "  </Collection>\n"
                        "</VTKFile>\n");
      std::fclose(f);
  } // write_collection_pvd

  template <typename real_t>
  void writeVTK(
        int const time // time
      , int const nx, int const ny, int const nz // sizes
      , char const *directory="."
      , char const *filename="out"
      , real_t const rho[]=nullptr // density
      , real_t const pre[]=nullptr // pressure 
      , real_t const ux[]=nullptr
      , real_t const uy[]=nullptr // velocities
      , real_t const uz[]=nullptr
  ) { // file
      int const dir = mkdir(directory, 0777); // delete the second argument in WIN32
      if (false && (0 == dir)) printf("Error: Cannot create output directory!\n");
      char datafilename[256];
      std::snprintf(datafilename, 255, "%s/%s_%07d.vti", directory, filename, time);
      auto const f = std::fopen(datafilename, "w");
      std::fprintf(f, "<?xml version=\"1.0\"?>\n"
                        "<!-- openLBMflow v1.0.1, www.lbmflow.com -->\n"
                        "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
      std::fprintf(f, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", nx-1, ny-1, nz-1);
      std::fprintf(f, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n", nx-1, ny-1, nz-1);
      std::fprintf(f, "    <PointData Scalars=\"scalars\">\n");

      auto const DataArrayFormat = "%s<DataArray type=\"Float%d\" Name=\"%s\""
                                   " NumberOfComponents=\"%d\" format=\"ascii\">\n";
      if (rho) { // write density
          std::fprintf(f, DataArrayFormat, "      ", sizeof(real_t)*8, "Density", 1);
          for (int z = 0; z < nz; z++) {
              for (int y = 0; y < ny; y++) {
                  for (int x = 0; x < nx; x++) {
                      std::fprintf(f,"%.4e ", rho[(x*ny + y)*nz + z]);
                  } // x
                  std::fprintf(f, "\n");
              } // y
          } // z
          std::fprintf(f, "      </DataArray>\n");
      } // rho
      
      if (pre) { // write pressure
          std::fprintf(f, DataArrayFormat, "      ", sizeof(real_t)*8, "Pressure", 1);
          for (int z = 0; z < nz; z++) {
              for (int y = 0; y < ny; y++) {
                  for (int x = 0; x < nx; x++) {
                      std::fprintf(f,"%.4e ", pre[(x*ny + y)*nz + z]);
                  } // x
                  std::fprintf(f, "\n");
              } // y
          } // z
          std::fprintf(f, "      </DataArray>\n");
      } // pre
      
      int const nuc = (nullptr != ux) + (nullptr != uy) + (nullptr != uz);
      if (ux || uy || uz) { // write velocity
          std::fprintf(f, DataArrayFormat, "      ", sizeof(real_t)*8, "Velocity", nuc);
          for (int z = 0; z < nz; z++) {
              for (int y = 0; y < ny; y++) {
                  for (int x = 0; x < nx; x++) {
                      const size_t xyz = (x*ny + y)*nz + z; // global index
                      if (ux) std::fprintf(f,"%.4e ", ux[xyz]);
                      if (uy) std::fprintf(f,"%.4e ", uy[xyz]);
                      if (uz) std::fprintf(f,"%.4e ", uz[xyz]);
                  } // x
                  std::fprintf(f, "\n");
              } // y
          } // z
          std::fprintf(f, "      </DataArray>\n");
      } // any velocity

      std::fprintf(f, "    </PointData>\n"
                        "    <CellData>\n"
                        "    </CellData>\n"
                        "  </Piece>\n"
                        "  </ImageData>\n"
                        "</VTKFile>\n");
      std::fclose(f);
  } // writeVTK

  
inline double wall_clock() { return clock()/double(CLOCKS_PER_SEC); }

template <typename real_t=double>
double outputSave(
      int const time
    , double const rho[]
    , double const ux[]
    , double const uy[]
    , double const uz[]
    , double const phi[] // index using phi[phindex(x, y, z)]
    , int const ranks[3]
    , int const myrank
    , int const Nx, int const Ny, int const Nz //  local lattice sites
    , int const nx, int const ny, int const nz // global lattice sites
    , int const save_rho
    , int const save_pre
    , int const save_vel
    , double const G
) {
    static double timer_start, step_start{0};

    double time_stop = wall_clock(); // stop internal timer

    // calculate performance in units of Mega Lattice Site Updates per second: MLUP/s
    double const Speed = ((nz*ny)*(nx*1e-6)*(time - step_start)/(time_stop - timer_start));
    step_start = time;

    double const mass = std::accumulate(rho, rho + Nz*Ny*Nx, 0.0);
    // ToDo: MPI_Allreduce(mass)

    auto const is_master = (0 == myrank);
    if (is_master) {
        printf("t= %d\tSpeed= %.1f MLUP/s mass= %.9f\n", time, Speed, mass);
#ifndef SuppressIO
    } // is_master

    int const nall = nx*ny*nz;
    int const offs[3] = {ranks[0]*Nx, ranks[1]*Ny, ranks[2]*Nz};
    std::vector<real_t> rho_all(save_rho*nall);
    std::vector<real_t> pre_all(save_pre*nall);
    std::vector<real_t> vel_all[3];
    for(int d = 0; d < 3; ++d) {
        vel_all[d] = std::vector<real_t>(save_vel*nall);
    } // d

    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            for (int z = 0; z < Nz; ++z) {
                index_t const xyz = indexyz(x, y, z, Nx, Ny); // local index into rho, ux, uy, uz
                size_t const gxyz = ((x + offs[0])*ny + (y + offs[1]))*nz + (z + offs[2]); // global index
                if (save_rho) rho_all[gxyz] = rho[xyz];
                if (save_pre) {
                    pre_all[gxyz] = rho[xyz]/3.0;
                    if (nullptr != phi) {
                        pre_all[gxyz] += G*pow2(phi[phindex(x, y, z)])/6.0;
                    }
                } // save_pre
                if (save_vel) {
                    vel_all[0][gxyz] = ux[xyz];
                    vel_all[1][gxyz] = uy[xyz];
                    vel_all[2][gxyz] = uz[xyz];
                } // velocities
            } // z
        } // y
    } // x

    if (is_master) {
        lbm_visualize::writeVTK(time, nx, ny, nz, "output", "openLBMflow",
                        rho_all.data(), pre_all.data(), vel_all[0].data(), vel_all[1].data(), vel_all[2].data());
#else
        std::printf("# SuppressIO for writeVTK\n");
#endif
    } // is_master

    timer_start = wall_clock(); // start internal timer again
    return Speed;
} // outputSave
  
} // namespace lbm_visualize
