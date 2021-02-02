#pragma once
 //
 //
 // vizualization functions (not critical for performance)
 //
 //

#include <cstdio> // std::fprintf, std::snprintf, std::fopen, std::fclose
#include <sys/stat.h> // mkdir

void write_collection_pvd(
      int const nx, int const ny, int const nz
    , char const *filename
    , char const *directory="."
    , int const time_total=1
    , int const time_save=1
) {
#ifdef SuppressIO
    return;
#endif                
    int const dir = mkdir(directory, 0777); // delete the second argument in WIN32
    // if (0 == dir) printf("Error: Cannot create output directory!\n");
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
#ifdef SuppressIO
    return;
#endif
    int const dir = mkdir(directory, 0777); // delete the second argument in WIN32
    if (0 == dir) printf("Error: Can't create output directory!\n");
    char datafilename[256];
    std::snprintf(datafilename, 255, "%s/%s_%07d.vti", directory, filename, time);
    auto const f = std::fopen(datafilename, "w");
    std::fprintf(f, "<?xml version=\"1.0\"?>\n"
                      "<!-- openLBMflow v1.0.1, www.lbmflow.com -->\n"
                      "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    std::fprintf(f, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", nx-1, ny-1, nz-1);
    std::fprintf(f, "  <Piece Extent=\"0 %d 0 %d 0 %d\">\n", nx-1, ny-1, nz-1);
    std::fprintf(f, "    <PointData Scalars=\"scalars\">\n");

    if (rho) { // write density
        std::fprintf(f, "      <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
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
        std::fprintf(f, "      <DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n");
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
        std::fprintf(f, "      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"%d\" format=\"ascii\">\n", nuc);
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

