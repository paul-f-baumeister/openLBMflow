// c++ -std=c++11 diff_vti_files_numerically.cxx -o numdiff
// use rapidxml

// from https://sourceforge.net/projects/rapidxml/files/latest/download
#include "rapidxml/rapidxml.hpp" // ::xml_document<>
#include "rapidxml/rapidxml_utils.hpp" // ::file<>

#include <cassert> // assert
#include <cstdio> // std::fprintf, stderr
#include <cmath> // std::abs
#include <cstring> // std::strcmp
#include <cstdlib> // std::strtod
#include <vector> // std::vector<T>
#include <algorithm> // std::min, std::max


    char const empty_string[] = "";
  
    inline char const * find_attribute(
          rapidxml::xml_node<> const *node
        , char const *const name
        , char const *const default_value=""
    ) { 
        if (nullptr == node) return empty_string;
        for (auto attr = node->first_attribute(); attr; attr = attr->next_attribute()) {
            if (0 == std::strcmp(name, attr->name())) {
                return attr->value();
            } // found
        } // attr
        return default_value;
    } // find_attribute 
  
    inline rapidxml::xml_node<> const * find_child(
          rapidxml::xml_node<> const *node
        , char const *const name
    ) { 
        if (nullptr == node) return nullptr;
        for (auto child = node->first_node(); child; child = child->next_sibling()) {
            if (0 == std::strcmp(name, child->name())) {
                return child;
            } // found
        } // attr
        return nullptr;
    } // find_child 

    template <typename real_t>
    std::vector<real_t> read_sequence(
          char const *sequence
        , size_t const reserve=0
    ) {
        char *end;
        char const *seq{sequence};
        std::vector<real_t> v;
        v.reserve(reserve);
        for (double f = std::strtod(seq, &end); seq != end; f = std::strtod(seq, &end)) {
            seq = end;
            if (errno == ERANGE){
                std::fprintf(stderr, "range error, got %g", f);
                errno = 0;
            } else {
                v.push_back(real_t(f));
            }
        } // f
        return v;
    } // read_sequence

#define print(...) if (log) { std::fprintf(log, __VA_ARGS__); }

    template <typename real_t=double>
    struct Observables {
        std::vector<real_t> density;
        std::vector<real_t> pressure;
        std::vector<real_t> velocity;
        size_t n;
    }; // Observables
    

    template <typename real_t=double>
    Observables<real_t> read_VTK_file(char const *filename, FILE* log=nullptr) {
      
        assert(filename != nullptr && "Filename invalid");
        assert(*filename && "Filename empty");

        rapidxml::file<> infile(filename);
        rapidxml::xml_document<> doc;
        doc.parse<0>(infile.data());

    // <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
    //   <ImageData WholeExtent="0 31 0 23 0 23" Origin="0 0 0" Spacing="1 1 1">
    //     <Piece Extent="0 31 0 23 0 23">
    //       <PointData Scalars="scalars">
    //         <DataArray type="Float32" Name="Density" NumberOfComponents="1" format="ascii">
        auto const VTKFile = doc.first_node("VTKFile");
        assert(VTKFile);
        auto const ImageData = find_child(VTKFile, "ImageData");
        assert(ImageData);
        auto const whole_extend = find_attribute(ImageData, "WholeExtent");
        print("# whole extent is %s\n", whole_extend);
        auto const Piece = find_child(ImageData, "Piece");
        assert(Piece);
        auto const piece_extend = find_attribute(Piece, "Extent");
        print("# Piece extent is %s\n", piece_extend);
        auto const ext = read_sequence<int>(piece_extend);
        assert(ext.size() == 6);
        int const lim[][2] = {{ext[0], ext[1]}, {ext[2], ext[3]}, {ext[4], ext[5]}};
        auto const total = size_t(lim[0][1] + 1 - lim[0][0])
                                *(lim[1][1] + 1 - lim[1][0])
                                *(lim[2][1] + 1 - lim[2][0]);
        auto const PointData = find_child(Piece, "PointData");
        assert(PointData);
        Observables<real_t> obs;
        obs.n = total;
        for(auto DataArray = PointData->first_node(); DataArray; 
                 DataArray = DataArray->next_sibling()) {
            assert(DataArray);
            auto const array_name = find_attribute(DataArray, "Name");
            auto const n_components = find_attribute(DataArray, "NumberOfComponents");
            auto const ncomp = std::atoi(n_components);
            print("# %d component%c in array '%s'\n", ncomp, (ncomp != 1)*'s', array_name);
            auto const values = DataArray->value();
            // print("# values: %s\n", values); // this can be very long
            auto const array = read_sequence<real_t>(values, 1.125*total*ncomp);
            print("# found %ld array entries, expect %ld\n", array.size(), total*ncomp);
            assert(total*ncomp == array.size());
                 if (0 == std::strcmp(array_name, "Density" )) { obs.density  = array; }
            else if (0 == std::strcmp(array_name, "Pressure")) { obs.pressure = array; }
            else if (0 == std::strcmp(array_name, "Velocity")) { obs.velocity = array; }
        } // Data
        return obs;
    } // read_VTK_file

    void add_stats(double value, double stats[5]) {
        stats[0] += 1;
        stats[1] += value;
        stats[2] += value*value;
        stats[3] = std::min(stats[3], value);
        stats[4] = std::max(stats[4], value);
    } // add_stats

    void show_stats(double stats[5], FILE* log, char const *title="") {
        size_t const n = std::max(1., stats[0]);
        auto const mean = stats[1]/n;
        print("# stats %s n=%ld%c", title, n, n?' ':'\n');
        if (n < 1) return;
        print("average %g +/- %g in [%g, %g]\n", mean, 0.0, stats[3], stats[4]);
    } // show_stats

    template <typename real_t=double>
    void analyse_deviation(
          real_t const array_now[]
        , real_t const array_ref[]
        , size_t const n
        , char const *name
        , FILE* log
        , int const ncomp=1
    ) {
        for(int icomp = 0; icomp < ncomp; ++icomp) {
            double ddev[] = {0, 0, 0, 9e33, -9e33};
            double adev[] = {0, 0, 0, 9e33, -9e33};
            double rdev[] = {0, 0, 0, 9e33, -9e33};
            for(size_t i = 0; i < n; ++i) {
                int const j = i*ncomp + icomp;
                auto const ref = array_ref[j];
                auto const now = array_now[j];
        //      print("%g ", ref);
                auto const dif = now - ref;
                auto const abs = std::abs(dif);
                if (std::abs(ref) > 1e-37) {
                    auto const rel = std::abs(dif/ref);
                    add_stats(rel, rdev);
                }
                add_stats(dif, ddev);
                add_stats(abs, adev);
            } // i
            char dir[4]; dir[0] = '-'*(ncomp > 1); dir[1] = 'x' + icomp; dir[2] = 0;
            print("# for array '%s%s':\n", name, dir);
            show_stats(ddev, log, "   plain difference");
            show_stats(adev, log, "absolute difference");
            show_stats(rdev, log, "relative difference");
        } // icomp
    } // analyse_deviation

int main(int argc, char *argv[]) {
    assert(argc > 2 && "2 filenames required as command line argument");
    auto const filename_ref = argv[1];
    auto const filename_now = argv[2];

    auto const log = nullptr; // stderr;
    auto const ref = read_VTK_file(filename_ref, log);
    auto const now = read_VTK_file(filename_now, log);
    auto const n = ref.n;
    assert(   n == now.n);
    
    auto const out = stdout;
    analyse_deviation(now.density.data(),  ref.density.data(),  n, "Density",  out);
    analyse_deviation(now.pressure.data(), ref.pressure.data(), n, "Pressure", out);
    analyse_deviation(now.velocity.data(), ref.velocity.data(), n, "Velocity", out, 3);

} // main
