// g++ -std=c++11 -O0 -g -pedantic -Wall -Wno-format-security -Wno-format test_MODULE.cxx && ./a.out
#include <cstdio> // std::printf
#include "MODULE.hxx" // ::all_tests
int main() {
    auto const stat = MODULE::all_tests(9);
    std::printf("\n# %s: all_tests = %i\n", __FILE__, int(stat));
    return int(stat);
} // main
