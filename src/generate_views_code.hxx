#pragma once

#include <cstdio> // std::printf, std::fprintf, stderr
#include <cassert> // assert

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

namespace generate_views_code {

  
  // configure here
  int constexpr MaxRank = 4;
  bool constexpr LowerBoundsCheck = true;
  bool constexpr UpperBoundsCheck = true;
  bool constexpr debug = true;
  bool constexpr warn_nullptr = true;
  bool constexpr write_access_on_transpose = false;

  status_t generate_default_constructor(int const rank, int const echo=0) {
      status_t stat(0);

      std::printf(R"STRING(
        
    // default constructor
    view%dD()
    : _data(nullptr))STRING", rank);
      for (int r = 0; r < rank; ++r) {
          std::printf(", _s%d(0)", r);
      } // r
      std::printf(", _mem(0)\n    {", rank);

      if (UpperBoundsCheck) {
          std::printf(R"STRING(
        for (int r = 0; r < Rank; ++r) {
            _n[r] = views::RangeUnknown;
        } // r
)STRING");
      } // UpperBoundsCheck

      if (debug) {
          std::printf(R"STRING(        std::printf("#  view%dD() default constructor\n");
)STRING", rank);
      } // debug
      std::printf("    } // default constructor\n\n");

      return stat;
  } // generate_default_constructor

  status_t generate_view_constructor(int const rank, int const echo=0) {
      status_t stat(0);

      std::printf(R"STRING(
        
    // data view constructor
    view%dD(T *const ptr)STRING", rank);
      for (int r = rank - 1; r >= 0; --r) {
          std::printf(", int_t s%d", r);
      } // r
      std::printf("=1)\n"
      "    : _data(ptr)");
      for (int r = 0; r < rank; ++r) {
          std::printf(", _s%d(s%d)", r, r);
      } // r
      std::printf(", _mem(0)\n    {");

      if (UpperBoundsCheck) {
          std::printf(R"STRING(
        for (int r = 0; r < Rank; ++r) {
            _n[r] = views::RangeUnknown;
        } // r
)STRING");
      } // UpperBoundsCheck

      if (warn_nullptr) {
          std::printf(R"STRING(        if (nullptr == ptr) std::printf("#  view%dD() received a nullptr!\n");
)STRING", rank);         
      } // warn_nullptr

      if (debug) {
          std::printf("        std::printf(\"#  view%dD(ptr", rank);
          for (int r = rank - 1; r >= 0; --r) {
              std::printf(", s%i=%%lld", r);
          } // r
          std::printf(") data view constructor\\n\"\n          ");
          for (int r = rank - 1; r >= 0; --r) {
              std::printf(", int64_t(s%i)", r);
          } // r
          std::printf(");\n");
      } // debug
      std::printf("    } // data view constructor\n\n");

      return stat;
  } // generate_view_constructor
  
  status_t generate_memory_constructor(int const rank, int const echo=0) {
      status_t stat(0);
      
      std::printf(R"STRING(
        
    // memory owner constructor
    view%dD()STRING", rank);
      for (int r = rank - 1; r >= 0; --r) {
          std::printf("size_t n%i, ", r);
      } // r
      std::printf("T init={0})\n"
      "    : _data(new T[");
      for (int r = rank - 1; r >= 0; --r) {
          std::printf("n%i%s", r, r?"*":"])");
      } // r
      for (int r = 0; r < rank; ++r) {
          std::printf(", _s%i(%s", r, r?"":"1)");
          for (int s = r - 1; s >= 0; --s) {
              std::printf("n%i%c", s, s?'*':')');
          } // s
      } // r
      std::printf(", _mem(");
      for (int r = rank - 1; r >= 0; --r) {
          std::printf("n%i*", r);
      } // r
      std::printf("sizeof(T))\n    {\n");

      if (UpperBoundsCheck) {
          for (int r = 0; r < rank; ++r) {
              std::printf("        _n[%i] = n%i;\n", r, r);
          } // r
      } // UpperBoundsCheck

      if (warn_nullptr) {
          std::printf(R"STRING(         if (nullptr == _data) std::printf("#  view%dD() found a nullptr after allocation!\n");
)STRING", rank);         
      } // warn_nullptr

      if (debug) {
          std::printf("        std::printf(\"#  view%dD(", rank);
          for (int r = rank - 1; r >= 0; --r) {
              std::printf("n%i=%%llu, ", r);
          } // r
          std::printf("init) memory owner constructor\\n\"\n          ");
          for (int r = rank - 1; r >= 0; --r) {
              std::printf(", n%i", r);
          } // r
          std::printf(");\n");
      } // debug
      std::printf("    } // memory owner constructor\n\n");

      return stat;
  } // generate_memory_constructor

  
  
  
  
  status_t generate_view_rank(int const rank, int const echo=0) {
      status_t stat(0);

      std::printf("\n\n\n\n\n// =============================================================\n");
      if (debug) {
          std::printf("#define debug_printf(...) std::printf(__VA_ARGS__)\n");
      } else {
          std::printf("#define debug_printf(...)\n");
      }
      std::printf(R"STRING(

template <typename T, typename int_t=int32_t>
class view%dD
{
public:
  typedef T value_type;
  static int constexpr Rank = %d;
)STRING", rank, rank);

      // add member functions
      stat += generate_default_constructor(rank, echo);
      stat += generate_view_constructor(rank, echo);
      stat += generate_memory_constructor(rank, echo);

      std::printf(R"STRING(
      
    // destructor
    ~view%dD() {
        if (_data && (_mem > 0)) {
            debug_printf("# ~view%dD() destructor tries to free %%g kByte\n", _mem*.001);
            delete[] _data;
        } else {
            debug_printf("# ~view%dD() destructor\n");
        } // is memory owner
    } // destructor

)STRING", rank, rank, rank);


      std::printf("    view%dD(view%dD<T> && rhs) {\n", rank, rank);
      if (debug) {
          std::printf("        std::printf(\"#  view%dD(view%dD<T> && rhs);\\n\");\n", rank, rank);
      } // debug
      std::printf(R"STRING(
        *this = std::move(rhs);
    } // move constructor

    view%dD(view%dD<T> const & rhs) = delete;

    view%dD & operator= (view%dD<T> const & rhs) = delete;

    view%dD & operator= (view%dD<T> && rhs) {
        _data = rhs._data;
)STRING", rank, rank, rank, rank, rank, rank);
      for (int r = 0; r < rank; ++r) {
          std::printf("        _s%i = rhs._s%i;\n", r, r);
      } // r
      if (UpperBoundsCheck) {
          std::printf("        for (int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }\n");
      } // UpperBoundsCheck
      if (debug) {
          std::printf("        std::printf(\"#  view%dD & operator= (view%dD<T> && rhs);\\n\");\n", rank, rank);
      } // debug
      std::printf(R"STRING(
        _mem  = rhs._mem; rhs._mem = 0; // steal ownership
        return *this;
    } // move assignment
    
private:
  
  inline size_t _index(
)STRING");
      for (int r = rank - 1; r >= 0; --r) {
          std::printf("      %c int_t i%i\n", (r < rank - 1)?',':' ', r);
      } // r
      std::printf("  ) const {\n");
      for (int r = 0; r < rank; ++r) {
          std::printf("      views::check_index(i%i, %i, ", r, r);
          std::printf(UpperBoundsCheck ? "_n[%i]" : "0", r);
          std::printf(", Rank, (void*)this);\n");
      } // r
      std::printf("      return ");
      for (int r = rank - 1; r > 0; --r) {
          std::printf("i%i*_s%i + ", r, r);
      } // r
      std::printf(R"STRING(i0*_s0;
  } // _index

public:

  T* data() const { return _data; }

  // access operators
)STRING");      

      for (int cst = 0; cst <= 1; ++cst) {
          auto const is_const = cst ? "const" : "     ";
          
          std::printf("  T %s & operator () (", is_const);
          for (int r = rank - 1; r >= 0; --r) {
              std::printf("int_t i%i%s", r, r?", ":")");
          } // r
          std::printf(" %s {\n            return _data[_index(", is_const);
          for (int r = rank - 1; r >= 0; --r) {
              std::printf("i%i%s", r, r?",":")]; }\n\n");
          } // r
          
      } // cst

 
      for (int cst = 0; cst <= 1; ++cst) {
          auto const is_const = cst ? "const" : "     ";
          if (rank > 1) {
              std::printf("  view%dD<T %s, int_t> operator [] (int_t i%i) %s {\n"
                          "            return view%dD<T %s, int_t>(_data + _index(i%i"
                          , rank - 1, is_const, rank - 1, is_const, rank - 1, is_const, rank - 1);
              for (int r = 0; r < rank - 1; ++r) {
                  std::printf(",0");
              } // r
              std::printf(")");
              for (int r = rank - 2; r >= 0; --r) {
                  std::printf(", _s%i", r);
              } // r
              std::printf("); }\n\n");
          } else {
              std::printf("  T %s & operator [] (int_t i0) %s { return _data[_index(i0)]; }\n"
                              , is_const                  , is_const);
          } // rank > 1
      } // cst
      

      std::printf("  // more advanced slicing operators\n");
      for (int cst = 0; cst <= 1; ++cst) {
          auto const is_const = cst ? "const" : "     ";
          for (int subrank = 1; subrank < rank; ++subrank) {
              int const fix = rank - subrank; // this may indices are fixed
              std::printf("  view%dD<T %s, int_t> operator () (", subrank, is_const);
              for (int r = fix - 1; r >= 0; --r) {
                  std::printf("int_t i%i, ", r);
              } // r
              std::printf("char const *slice=\"");
              for (int r = 0; r < rank; ++r) {
                  std::printf("%c", (r >= fix)?':':'_');
              } // r
              std::printf("\")%c      %s\n  {\n", cst?'\n':' ', is_const);
              std::printf("      int_t const strides[Rank] = {");
              for (int r = 0; r < rank; ++r) {
                  std::printf("_s%i%s", r, (r < rank - 1)?", ":"};\n");
              } // r
              std::printf("      int_t const indices[%d] = {", fix);
              for (int r = 0; r < fix; ++r) {
                  std::printf("i%i%s", r, (r < fix - 1)?", ":"};\n");
              } // r
              std::printf("      int_t sp[%d], ip[Rank];\n", subrank);
              std::printf("      views::subview_from_string(sp, ip, slice, Rank, %d, strides, indices);\n", subrank);
              std::printf("      return view%dD<T %s, int_t>(_data + _index(", subrank, is_const);
              for(int r = rank - 1; r >= 0; --r) {
                  std::printf("ip[%i]%c", r, r?',':')');
              } // r
              for(int r = subrank - 1; r >= 0; --r) {
                  std::printf(", sp[%i]", r);
              } // r
              std::printf(");\n");
              std::printf("  } // slicing operator\n");
          } // subrank
      } // cst

      // transpose
      if (rank > 1) {
          std::printf(R"STRING(

  view%dD transpose(uint64_t permutation=views::DefaultTranspose)
      const
  {
      int_t const strides[Rank] = {)STRING", rank);
          for (int r = 0; r < rank; ++r) {
              std::printf("_s%i%s", r, (r < rank - 1)? ", ": "};\n");
          } // r
          std::printf(R"STRING(
      uint8_t rank_count[Rank];
      int_t sp[Rank]; // permuted strides
      for (int r = 0; r < Rank; ++r) {
          rank_count[r] = 0;
          sp[r] = 0;
      } // r
      bool non_transposed{true};
      uint64_t const mask = (Rank < 16) ? ~((~uint64_t(0)) << (4*Rank)) : ~uint64_t(0);
      auto perm{permutation & mask}; // non-const copy
      for (int r = 0; r < Rank; ++r) {
          uint8_t const permuted_rank = perm & 0xf;
          assert(permuted_rank < Rank);
          perm >>= 4; // shift 4 bits out (one hexadecimal digit)
          non_transposed &= (permuted_rank == r);
          ++rank_count[permuted_rank];
          sp[r] = strides[permuted_rank];
      } // r
      if (non_transposed) {
          debug_printf("#  view%dD transpose with non-transposed permutation %%%d.%dllx\n", permutation & mask);
      } // warn
      int non_standard{0};
      for (int r = 0; r < Rank; ++r) {
          non_standard += (1 != rank_count[r]);
      } // r
      if (non_standard > 0) {
          debug_printf("#  view%dD transpose with non-standard permutation %%%d.%dllx\n", permutation & mask);
      } // warn

      return view%dD(_data)STRING", rank, rank, rank, rank, rank, rank, rank);
          for (int r = rank - 1; r >= 0; --r) {
              std::printf(", sp[%i]", r);
          } // r
          std::printf(");\n  } // transpose\n\n");
      } // rank > 1

      
      std::printf(R"STRING(
private:

  // private data members
  T * _data;
  int_t)STRING");
      for (int r = 0; r < rank; ++r) {
          std::printf(" _s%i%c", r, (r == rank - 1)?';':',');
      } // r
      std::printf(" // strides\n  size_t _mem; // only > 0 if memory owner\n");
      if (UpperBoundsCheck) {
          std::printf("  int_t _n[Rank];\n");
      } // UpperBoundsCheck
      
      std::printf("\n}; // class view%dD\n", rank);

      std::printf("\n#undef debug_printf\n\n");
      
      return stat;
  } // generate_view_rank
  
  
  
  status_t generate_test_rank(int const rank, int const echo=0) {
      status_t stat(0);

      std::printf("      { // open scope for Rank=%d tests\n", rank);
      std::printf("          std::printf(\"\\n# Tests for view%dD\\n\");\n\n", rank);
      
      std::printf("          { view%dD<int8_t> D; } // test default constructor\n\n", rank);

      std::printf("          // test memory owner constructor\n");
      std::printf("          view%dD<int8_t> A(", rank);
      for (int r = rank - 1; r >= 0; --r) {
          std::printf("n[%i]%c", r, r?',':')');
      } // r
      std::printf("; assert(nullptr != A.data());\n\n");
      
      std::printf("          // test read-write access\n          auto & a = A(");
      for (int r = rank - 1; r >= 0; --r) {
          std::printf("n[%i]-1%c", r, r?',':')'); // just an example to pass n[] as strides, not really meaningful
      } // r
      std::printf(";\n"
                  "          a = -1; assert(-1 == a);\n");

      
      std::printf("\n          // test [] operator\n");
      std::printf("          { auto const S = A[n[%d]-1]; ", rank - 1);
      if (rank > 1) {
          std::printf("assert(%d == S.Rank); } \n\n", rank - 1);
      } else {
          std::printf("assert(S == S); } // scalar value\n\n");
      } // rank > 1
      
      
      if (rank > 1) {
          std::printf("\n          // test () subview operators\n");
          for (int subrank = 1; subrank < rank; ++subrank) {
              std::printf("          { auto const SV = A(");
              for (int r = rank - 1; r >= subrank; --r) {
                  std::printf("n[%i]-1%c", r, (r > subrank)?',':')');
              } // r
              std::printf("; assert(%d == SV.Rank); }\n", subrank);
          } // subrank
          
          std::printf("\n          { auto const At = A.transpose(); } // test default transpose\n");
          std::printf("\n          { auto const A1 = A.transpose(views::NonTransposed); } // test identical view\n");
          
      } // rank > 1
      
      
      std::printf("\n      } // close scope for Rank=%d tests\n\n", rank);
      
      return stat;
  } // generate_test_rank
  
  
  
  
  
  
  
  
  
  status_t generate(int const echo=0) {
      status_t stat(0);

      std::printf(R"STRING(#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // uint_t
#include <utility> // std::move
#include <algorithm> // std::fill

namespace views {

  uint64_t constexpr NonTransposed    = 0xfedcba9876543210;
  uint64_t constexpr DefaultTranspose = 0xfedcba9876543201;
)STRING");

    if (UpperBoundsCheck) std::printf("  int64_t  constexpr RangeUnknown     = -1;\n\n");
    
    std::printf(R"STRING(
  // auto-generated code, please modify generator code %s

  template <typename int_t>
  inline void check_index(
        int_t const index
      , int const rank
      , int64_t const number
      , int const Rank
      , void* const at=nullptr
  ) {
)STRING", __FILE__);
    
      if (LowerBoundsCheck) {
          std::printf(R"STRING(
          if (index < 0) {
              std::fprintf(stderr, "\nERROR: view%%dD at %%p(i%%i= %%lli < 0)\n\n",
                            Rank, at, rank, int64_t(index));
          }
          assert(index >= 0);
)STRING");
      } // LowerBoundsCheck
      if (UpperBoundsCheck) {
          std::printf(R"STRING(
          if (number > RangeUnknown) {
              if (index >= number) {
                  std::fprintf(stderr, "\nERROR: view%%dD at %%p(i%%i= %%lli >= %%lld)\n\n", 
                                Rank, at, rank, int64_t(index), int64_t(number));
              }
              assert(index < number);
          } // RangeUnknown
)STRING");
      } // UpperBoundsCheck
  std::printf(R"STRING(
  } // check_index

  template <typename int_t>
  void subview_from_string(
        int_t sp[]
      , int_t ip[]
      , char const *const slice
      , int const nRank
      , int const nColon
      , int_t const strides[]
      , int_t const indices[]
  ) {
      assert(nullptr != slice);
      assert(nColon < nRank);
      int colon{nColon}, index{nRank - nColon};
      if (%d) std::printf("#  view%%dD --> sub view%%dD", nRank, nColon);
      for (int c = 0; c < nRank; ++c) {
          char const s = slice[c];
          assert('\0' != s && "slicing string too short");
          int const r = nRank - 1 - c;
          if (':' == s) {
              --colon;
              assert(colon >= 0 && "slicing string has too many colons");
              sp[colon] = strides[r];
              ip[r] = 0;
              if (%d) std::printf("%%c%%c", c?',':'(', s);
          } else {
              --index;
              assert(index >= 0 && "slicing string has not enough colons");
              ip[r] = indices[index];
              if (%d) std::printf("%%c%%i", c?',':'(', ip[r]);
          }
      } // r
      if (%d) std::printf(")\n");
      assert(0 == colon);
      assert(0 == index);
      assert('\0' == slice[nRank] || ' ' == slice[nRank]); // string stops or a space separator (consistency)
  } // subview_from_string

} // namespace views
)STRING", debug, debug, debug, debug);

  
  

      for (int rank = 1; rank <= MaxRank; ++rank) {
          stat += generate_view_rank(rank, echo);
      } // rank

      

      
      
      
      // tests
      std::printf(R"STRING(

namespace views {
  
  inline int all_tests(int const echo=0) {
      int stat{0};

      size_t const n[] = {8,7,6,5,4,3,2,1,4,3,2,1,4,3,2,1}; // test dimensions

)STRING");
    
      for (int rank = 1; rank <= MaxRank; ++rank) {
          stat += generate_test_rank(rank, echo);
      } // rank

      std::printf(R"STRING(
      return stat;
  } // all_tests

} // namespace views
)STRING");
      
      return stat;
  } // generate

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += generate(echo);
      return stat;
  } // all_tests

} // namespace generate_views_code
