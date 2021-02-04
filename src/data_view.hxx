#pragma once

#include <cstdio> // printf, std::fflush, stdout
#include <cassert> // assert
#include <cstdint> // uint_t --> replace size_t with this?
#include <algorithm> // std::fill
#include <utility> // std::move

typedef int status_t;

uint64_t constexpr NonTransposed = 0xfedcba9876543210;

#define UpperBoundsCheck
#ifdef  UpperBoundsCheck
  int64_t constexpr RangeUnknown = -1;
  #define UpperBound(UB) (UB)
#else  // UpperBoundsCheck
  #define UpperBound(UB) 0
#endif // UpperBoundsCheck

#ifdef  LowerBoundsCheck
  bool constexpr LowerBoundCheck = true;
#else
  bool constexpr LowerBoundCheck = false;
#endif

// for view1D
#define debug_printf(...) printf(__VA_ARGS__)
// #define debug_printf(...)

namespace data_view {
  
  template <typename int_t>
  inline void check_index(
        int_t const index
      , int const rank
      , int64_t const number
      , int const Rank
      , void* const at=nullptr
  ) {
      if (LowerBoundCheck) {
          if (index < 0) {
              std::fprintf(stderr, "\nERROR: view%dD at %p(i%i= %li < 0)\n\n",
                                                Rank, at, rank, index);
          }
          assert(0 <= index);
      } // LowerBoundCheck
      #ifdef UpperBoundsCheck
          if (number > RangeUnknown) {
              if (index >= number) {
                  std::fprintf(stderr, "\nERROR: view%dD at %p(i%i= %li >= %ld)\n\n", 
                                                    Rank, at, rank, index, number);
              }
              assert(index < number);
          } // RangeUnknown
      #endif // UpperBoundsCheck
  } // check_index

} // namespace data_view  
  
template <typename T, typename int_t=int64_t>
class view1D {
  // the data view can be transposed at no cost as we include only the two strides
  // we have to check if the smaller stride is 1 before we use it in dgemm
public:
  static int constexpr Rank = 1;

  view1D() : _data(nullptr), _s0(0), _mem(0) {
#ifdef  UpperBoundsCheck
      _n0 = RangeUnknown;
#endif // UpperBoundsCheck
      debug_printf("# view1D() default constructor\n");
  } // default constructor

  view1D(T* const ptr, int_t const stride=1) 
    : _data(ptr), _s0(stride), _mem(0) { 
#ifdef  UpperBoundsCheck
      _n0 = RangeUnknown;
#endif // UpperBoundsCheck
      debug_printf("# view1D(ptr, stride=%i) constructor\n", stride);
  } // data view constructor

  view1D(size_t const n0, T const init_value={0})
    : _data(new T[n0]), _s0(1), _mem(n0*sizeof(T)) {
#ifdef  UpperBoundsCheck
      _n0 = n0;
#endif // UpperBoundsCheck
      debug_printf("# view1D(%d [,init_value]) constructor allocates %g kByte\n", n0, _mem*.001);
      std::fill(_data, _data + n0, init_value); // warning! first touch here!
  } // memory owning constructor, ToDo: move this to the derived type

  ~view1D() {
      if (_data && (_mem > 0)) {
          debug_printf("# ~view1D() destructor tries to free %g kByte\n", _mem*.001);
          delete[] _data;
      } else {
          debug_printf("# ~view1D() destructor\n");
      } // is memory owner
  } // destructor

  // view1D(view1D<T> && rhs) = delete;
  view1D(view1D<T> && rhs) {
      debug_printf("# view1D(view1D<T> && rhs);\n");
      *this = std::move(rhs);
  } // move constructor

  view1D(view1D<T> const & rhs) = delete;
  // view1D(view1D<T> const & rhs) {
  //     // copy constructor (may lead to problems, who owns the memory afterwards?)
  //     debug_printf("# view1D(view1D<T> const & rhs);\n");
  //     *this = rhs;
  // } // copy constructor

  view1D& operator= (view1D<T> && rhs) {
      debug_printf("# view1D& operator= (view1D<T> && rhs);\n");
      _data = rhs._data;
      _s0   = rhs._s0;
      _mem  = rhs._mem; rhs._mem = 0; // steal ownership
#ifdef  UpperBoundsCheck
      _n0 = rhs._n0;
#endif // UpperBoundsCheck
      return *this;
  } // move assignment

  view1D& operator= (view1D<T> const & rhs) = delete;
  // view1D& operator= (view1D<T> const & rhs) {
  //     debug_printf("# view1D& operator= (view1D<T> const & rhs);\n");
  //     _data = rhs._data;
  //     _s0   = rhs._s0;
  //     _mem  = 0; // we are just a shallow copy
#ifdef  UpperBoundsCheck
//       _n0 = rhs._n0;
#endif // UpperBoundsCheck
  //     return *this;
  // } // move assignment

private:
  
  inline size_t _index(
        int_t const i0
  ) const {
      data_view::check_index(i0, 0, UpperBound(_n0), Rank, (void*)this);
      return i0*_s0;
  } // _index

public:
  
  T const & operator () (int_t const i0) const { return _data[_index(i0)]; }
  T       & operator () (int_t const i0)       { return _data[_index(i0)]; }
  
  // no slicing, [] also mean data access (like vectors)
  T const & operator [] (int_t const i0) const { return _data[_index(i0)]; }
  T       & operator [] (int_t const i0)       { return _data[_index(i0)]; }

  T*      data()   const { return _data; }
  
//size_t  stride() const { return 1; } // view1D did not exist in old versions

private:
  // private data members
  T * _data;
  int_t _s0;
  size_t _mem; // only > 0 if memory owner
#ifdef  UpperBoundsCheck
  int64_t _n0;
#endif // UpperBoundsCheck

}; // view1D

#undef  debug_printf













// for view2D
#define debug_printf(...) printf(__VA_ARGS__)
// #define debug_printf(...)

template <typename T, typename int_t=int64_t>
class view2D {
  // the data view can be transposed at no cost as we include only the two strides
  // we have to check if the smaller stride is 1 before we use it in dgemm
public:
  static int constexpr Rank = 2;

  view2D() : _data(nullptr), _s0(0), _s1(0), _mem(0) {
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = RangeUnknown; }
#endif // UpperBoundsCheck
      debug_printf("# view2D() default constructor\n");
  } // default constructor

  view2D(T* const ptr, int_t const stride1, int_t const stride0=1) 
    : _data(ptr), _s0(stride0), _s1(stride1), _mem(0) { 
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = RangeUnknown; }
#endif // UpperBoundsCheck
      debug_printf("# view2D(ptr, s1=%i, s0=%i) constructor\n", stride1, stride0);
  } // data view constructor

  view2D(size_t const n1, size_t const n0, T const init_value={0})
    : _data(new T[n1*n0]), _s0(1), _s1(n0), _mem(n1*n0*sizeof(T)) {
#ifdef  UpperBoundsCheck
      _n[0] = n0; _n[1] = n1;
#endif // UpperBoundsCheck
      debug_printf("# view2D(%d, %d [,init_value]) constructor allocates %g kByte\n", n1, n0, _mem*.001);
      std::fill(_data, _data + n1*n0, init_value); // warning! first touch here!
  } // memory owning constructor, ToDo: move this to the derived type

  ~view2D() {
      if (_data && (_mem > 0)) {
          debug_printf("# ~view2D() destructor tries to free %g kByte\n", _mem*.001);
          delete[] _data;
      } else {
          debug_printf("# ~view2D() destructor\n");
      } // is memory owner
  } // destructor

  // view2D(view2D<T> && rhs) = delete;
  view2D(view2D<T> && rhs) {
      debug_printf("# view2D(view2D<T> && rhs);\n");
      *this = std::move(rhs);
  } // move constructor

  view2D(view2D<T> const & rhs) = delete;
  // view2D(view2D<T> const & rhs) {
  //     // copy constructor (may lead to problems, who owns the memory afterwards?)
  //     debug_printf("# view2D(view2D<T> const & rhs);\n");
  //     *this = rhs;
  // } // copy constructor

  view2D& operator= (view2D<T> && rhs) {
      debug_printf("# view2D& operator= (view2D<T> && rhs);\n");
      _data = rhs._data;
      _s0   = rhs._s0;
      _s1   = rhs._s1;
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }
#endif // UpperBoundsCheck
      _mem  = rhs._mem; rhs._mem = 0; // steal ownership
      return *this;
  } // move assignment

  view2D& operator= (view2D<T> const & rhs) = delete;
  // view2D& operator= (view2D<T> const & rhs) {
  //     debug_printf("# view2D& operator= (view2D<T> const & rhs);\n");
  //     _data = rhs._data;
  //     _s0   = rhs._s0; 
  //     _s1   = rhs._s1;
  //     _mem  = 0; // we are just a shallow copy
#ifdef  UpperBoundsCheck
//       for(int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }
#endif // UpperBoundsCheck
  //     return *this;
  // } // move assignment

private:
  
  inline size_t _index(
        int_t const i1
      , int_t const i0
  ) const {
      data_view::check_index(i0, 0, UpperBound(_n[0]), Rank, (void*)this);
      data_view::check_index(i1, 1, UpperBound(_n[1]), Rank, (void*)this);
      return i1*_s1 + i0*_s0;
  } // _index

public:

  T const & operator () (int_t const i1, int_t const i0) const { return _data[_index(i1,i0)]; }
  T       & operator () (int_t const i1, int_t const i0)       { return _data[_index(i1,i0)]; }

  // slicing
  view1D<T const,int_t> operator [] (int_t const i1) const { return view1D<T const,int_t>(_data + _index(i1,0), _s0); }
  view1D<T      ,int_t> operator [] (int_t const i1)       { return view1D<T      ,int_t>(_data + _index(i1,0), _s0); }

  T*      data()   const { return _data; }

  // for backwards compatability with older versions
  size_t  stride() const { assert(_s0 >= _s1 && _s1 >= 0); return _s1; }

// #define WRITE_ACCESS_ON_TRANSPOSE // does not work!

  view2D transpose(uint64_t const permutation=0x01) // 0x10 is return an exact view like *this
#ifndef WRITE_ACCESS_ON_TRANSPOSE
      const
#endif // WRITE_ACCESS_ON_TRANSPOSE
  {
      int_t const new_s1 = ((permutation >> 4) & 0xf) ? _s1 : _s0;
      int_t const new_s0 = ( permutation       & 0xf) ? _s1 : _s0;
      if (new_s0 == _s0 || new_s1 == _s1) {
          debug_printf("# view2D non-transposing permutation %2.2x\n", permutation);
      } else
      if (new_s0 != _s1 || new_s1 != _s0) {
          debug_printf("# view2D transpose with non-standard permutation %2.2x\n", permutation);
      } // warn
      return view2D(_data, new_s1, new_s0); // UpperBoundsCheck not possible on transpose
  } // transpose

private:
  // private data members
  T * _data;
  int_t _s0, _s1;
  size_t _mem; // only > 0 if memory owner
#ifdef  UpperBoundsCheck
  int64_t _n[Rank];
#endif // UpperBoundsCheck

}; // view2D




  template<typename Ta, typename Tb, typename Tc>
  void gemm(
        view2D<Tc> & c //            result matrix, shape(N,M)
      , int const N
      , view2D<Tb> const & b //  left input matrix, shape(N,K)
      , int const K
      , view2D<Ta> const & a // right input matrix, shape(K,M)
      , int const aM=-1 // user specify M different from min(a.stride(), c.stride())
      , char const beta='0' // '0':overwrite elements of c, else add to c
  ) {
      // define a generic matrix-matrix multiplication: c(N,M) = b(N,K) * a(K,M)

      int const M = (-1 == aM) ? std::min(c.stride(), a.stride()) : aM;
      if (M > a.stride()) error("M= %d > %ld =a.stride", M, a.stride());
      if (K > b.stride()) error("K= %d > %ld =b.stride", M, b.stride());
      if (M > c.stride()) error("M= %d > %ld =c.stride", M, c.stride());
      assert( M <= a.stride() );
      assert( K <= b.stride() );
      assert( M <= c.stride() );
      for(int n = 0; n < N; ++n) {
          for(int m = 0; m < M; ++m) {
              Tc t(0);
              for(int k = 0; k < K; ++k) { // contraction index
                  t += b(n,k) * a(k,m);
              } // k
              if ('0' == beta) { c(n,m) = t; } else { c(n,m) += t; } // store
          } // m
      } // n
  } // gemm
         
         
#undef  debug_printf












// for view3D
#define debug_printf(...) printf(__VA_ARGS__)
// #define debug_printf(...)

template <typename T, typename int_t=int64_t>
class view3D {
public:
  static int constexpr Rank = 3;

  view3D() : _data(nullptr), _s0(0), _s1(0), _s2(0), _mem(0) { 
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = RangeUnknown; }
#endif // UpperBoundsCheck
      debug_printf("# view3D() default constructor\n");
  } // default constructor

  view3D(T* const ptr, int_t const stride2, int_t const stride1, int_t const stride0=1) 
    : _data(ptr), _s0(stride0), _s1(stride1), _s2(stride2), _mem(0) { 
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = RangeUnknown; }
#endif // UpperBoundsCheck
      debug_printf("# view3D(ptr, s2=%i, s1=%i, s0=%i) constructor\n", stride2, stride1, stride0);
  } // data view constructor

  view3D(size_t const n2, size_t const n1, size_t const n0, T const init_value={0})
    : _data(new T[n2*n1*n0]), _s0(1), _s1(n0), _s2(n1*n0), _mem(n2*n1*n0*sizeof(T)) {
#ifdef  UpperBoundsCheck
      _n[0] = n0; _n[1] = n1; _n[2] = n2;
#endif // UpperBoundsCheck
      debug_printf("# view3D(%d, %d, %d [,init_value]) constructor allocates %g kByte\n", n2, n1, n0, _mem*.001);
      std::fill(_data, _data + n2*n1*n0, init_value); // warning! first touch here!
  } // memory owning constructor, ToDo: move this to the derived type

  ~view3D() {
      if (_data && (_mem > 0)) {
          debug_printf("# ~view3D() destructor tries to free %g kByte\n", _mem*.001);
          delete[] _data;
      } else {
          debug_printf("# ~view3D() destructor\n");
      } // is memory owner
  } // destructor

  // view3D(view3D<T> && rhs) = delete;
  view3D(view3D<T> && rhs) {
      debug_printf("# view3D(view3D<T> && rhs);\n");
      *this = std::move(rhs);
  } // move constructor

  view3D(view3D<T> const & rhs) = delete;
  // view3D(view3D<T> const & rhs) {
  //     // copy constructor (may lead to problems, who owns the memory afterwards?)
  //     debug_printf("# view3D(view3D<T> const & rhs);\n");
  //     *this = rhs;
  // } // copy constructor

  view3D& operator= (view3D<T> && rhs) {
      debug_printf("# view3D& operator= (view3D<T> && rhs);\n");
      _data = rhs._data;
      _s0   = rhs._s0;
      _s1   = rhs._s1;
      _s2   = rhs._s2;
      _mem  = rhs._mem; rhs._mem = 0; // steal ownership
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }
#endif // UpperBoundsCheck
      return *this;
  } // move assignment

  view3D& operator= (view3D<T> const & rhs) = delete;
  // view3D& operator= (view3D<T> const & rhs) {
  //     debug_printf("# view3D& operator= (view3D<T> const & rhs);\n");
  //     _data = rhs._data;
  //     _s0   = rhs._s0; 
  //     _s1   = rhs._s1;
  //     _s2   = rhs._s2;
  //     _mem  = 0; // we are just a shallow copy
#ifdef  UpperBoundsCheck
//       for(int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }
#endif // UpperBoundsCheck
  //     return *this;
  // } // move assignment

private:
  
  inline size_t _index(
        int_t const i2
      , int_t const i1
      , int_t const i0
  ) const {
      data_view::check_index(i0, 0, UpperBound(_n[0]), Rank, (void*)this);
      data_view::check_index(i1, 1, UpperBound(_n[1]), Rank, (void*)this);
      data_view::check_index(i2, 2, UpperBound(_n[2]), Rank, (void*)this);
      return i2*_s2 + i1*_s1 + i0*_s0;
  } // _index

public:
  
  T const & operator () (int_t const i2, int_t const i1, int_t const i0) const { return _data[_index(i2,i1,i0)]; }
  T       & operator () (int_t const i2, int_t const i1, int_t const i0)       { return _data[_index(i2,i1,i0)]; }

  // slicing
  view2D<T const,int_t> operator [] (int_t const i2) const {
            return view2D<T const,int_t>(_data + _index(i2,0,0), _s1, _s0); }
  view2D<T      ,int_t> operator [] (int_t const i2)       {
            return view2D<T      ,int_t>(_data + _index(i2,0,0), _s1, _s0); }

  T* data() const { return _data; }
  
  // for backwards compatability with older versions
  // where view3D v(n2,n1,n0); v(i2,i1,i0) --> (i2*n1 + i1)*n0 + i0;
  //        then _s2 == n1*n0, _s1 == n0, _s0 == 1
  size_t  stride() const {
      assert(is_ordered()); assert(_s1 >= 0);
      return _s1; // == n0
  } // stride

  size_t  dim1()   const {
      assert(is_ordered()); assert(_s2 >= 0); assert(_s1 > 0);
      return _s2/_s1; // == n1
  } // dim1

  view3D transpose(uint64_t const permutation=0x201) // e.g. 0x210 returns an exact view like *this
#ifndef WRITE_ACCESS_ON_TRANSPOSE
      const
#endif // WRITE_ACCESS_ON_TRANSPOSE
  {
      int_t const strides[Rank] = {_s0, _s1, _s2};  // extend here for higher Rank
      uint8_t rank_count[Rank] = {0, 0, 0};         // extend here for higher Rank
      int_t sp[Rank]; // permuted strides
//    uint8_t pr[Rank]; // permuted rank
      bool non_transposed{true};
      auto perm{permutation}; // non-const copy
      for(int r = 0; r < Rank; ++r) {
          uint8_t const permuted_rank = perm & 0xf;
          perm >>= 4; // shift 4 bits out (one hexadecimal digit)
          assert(permuted_rank < Rank);
          non_transposed &= (permuted_rank == r);
          ++rank_count[permuted_rank];
          auto const permuted_stride = strides[permuted_rank];
//        pr[r] = permuted_rank;
          sp[r] = permuted_stride;
      } // r
      if (non_transposed) {
          debug_printf("# view3D transpose with non-transposed permutation %3.3x\n", permutation);
      } // warn
      int non_standard{0};
      for(int r = 0; r < Rank; ++r) {
          non_standard += (1 != rank_count[r]);
      } // r
      if (non_standard) {
          debug_printf("# view3D transpose with non-standard permutation %3.3x\n", permutation);
      } // warn

      return view3D(_data, sp[2], sp[1], sp[0]);    // extend here for higher Rank
  } // transpose

private:

  bool is_ordered() const { return (_s0 <= _s1 && _s1 <= _s2); }
  
  // private data members
  T * _data;
  int_t _s0, _s1, _s2;
  size_t _mem; // only > 0 if memory owner
#ifdef  UpperBoundsCheck
  int64_t _n[Rank];
#endif // UpperBoundsCheck

}; // view3D

#undef  debug_printf









// for view4D
#define debug_printf(...) printf(__VA_ARGS__)
// #define debug_printf(...)

template <typename T, typename int_t=int64_t>
class view4D {
public:
  static int constexpr Rank = 4;

  view4D() : _data(nullptr), _s0(0), _s1(0), _s2(0), _s3(0), _mem(0) { 
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = RangeUnknown; }
#endif // UpperBoundsCheck
      debug_printf("# view4D() default constructor\n");
  } // default constructor

  view4D(T* const ptr, int_t const stride3, int_t const stride2, int_t const stride1, int_t const stride0=1) 
    : _data(ptr), _s0(stride0), _s1(stride1), _s2(stride2), _s3(stride3), _mem(0) { 
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = RangeUnknown; }
#endif // UpperBoundsCheck
      debug_printf("# view4D(ptr, s3=%i, s2=%i, s1=%i, s0=%i) constructor\n", stride3, stride2, stride1, stride0);
  } // data view constructor

  view4D(size_t const n3, size_t const n2, size_t const n1, size_t const n0, T const init_value={0})
    : _data(new T[n3*n2*n1*n0]), _s0(1), _s1(n0), _s2(n1*n0), _s3(n2*n1*n0), _mem(n3*n2*n1*n0*sizeof(T)) {
#ifdef  UpperBoundsCheck
      _n[0] = n0; _n[1] = n1; _n[2] = n2; _n[3] = n3;
#endif // UpperBoundsCheck
      debug_printf("# view4D(%d, %d, %d, %d [,init_value]) constructor allocates %g kByte\n", n3, n2, n1, n0, _mem*.001);
      std::fill(_data, _data + n3*n2*n1*n0, init_value); // warning! first touch here!
  } // memory owning constructor, ToDo: move this to the derived type

  ~view4D() {
      if (_data && (_mem > 0)) {
          debug_printf("# ~view4D() destructor tries to free %g kByte\n", _mem*.001);
          delete[] _data;
      } else {
          debug_printf("# ~view4D() destructor\n");
      } // is memory owner
  } // destructor

  // view4D(view4D<T> && rhs) = delete;
  view4D(view4D<T> && rhs) {
      debug_printf("# view4D(view4D<T> && rhs);\n");
      *this = std::move(rhs);
  } // move constructor

  view4D(view4D<T> const & rhs) = delete;
  // view4D(view4D<T> const & rhs) {
  //     // copy constructor (may lead to problems, who owns the memory afterwards?)
  //     debug_printf("# view4D(view4D<T> const & rhs);\n");
  //     *this = rhs;
  // } // copy constructor

  view4D& operator= (view4D<T> && rhs) {
      debug_printf("# view4D& operator= (view4D<T> && rhs);\n");
      _data = rhs._data;
      _s0   = rhs._s0;
      _s1   = rhs._s1;
      _s2   = rhs._s2;
      _s3   = rhs._s3;
      _mem  = rhs._mem; rhs._mem = 0; // steal ownership
#ifdef  UpperBoundsCheck
      for(int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }
#endif // UpperBoundsCheck
      return *this;
  } // move assignment

  view4D& operator= (view4D<T> const & rhs) = delete;
  // view4D& operator= (view4D<T> const & rhs) {
  //     debug_printf("# view4D& operator= (view4D<T> const & rhs);\n");
  //     _data = rhs._data;
  //     _s0   = rhs._s0; 
  //     _s1   = rhs._s1;
  //     _s2   = rhs._s2;
  //     _s3   = rhs._s3;
  //     _mem  = 0; // we are just a shallow copy
#ifdef  UpperBoundsCheck
//       for(int r = 0; r < Rank; ++r) { _n[r] = rhs._n[r]; }
#endif // UpperBoundsCheck
  //     return *this;
  // } // move assignment

private:
  
  inline size_t _index(
        int_t const i3
      , int_t const i2
      , int_t const i1
      , int_t const i0
  ) const {
      data_view::check_index(i0, 0, UpperBound(_n[0]), Rank, (void*)this);
      data_view::check_index(i1, 1, UpperBound(_n[1]), Rank, (void*)this);
      data_view::check_index(i2, 2, UpperBound(_n[2]), Rank, (void*)this);
      data_view::check_index(i3, 3, UpperBound(_n[3]), Rank, (void*)this);
      return i3*_s3 + i2*_s2 + i1*_s1 + i0*_s0;
  } // _index

public:

  T const & operator () (int_t const i3, int_t const i2, int_t const i1, int_t const i0) const {
            return _data[_index(i3,i2,i1,i0)]; }
  T       & operator () (int_t const i3, int_t const i2, int_t const i1, int_t const i0)       {
            return _data[_index(i3,i2,i1,i0)]; }

  // slicing
  view3D<T const,int_t> operator [] (int_t const i3) const {
            return view3D<T const,int_t>(_data + _index(i3,0,0,0), _s2, _s1, _s0); }
  view3D<T      ,int_t> operator [] (int_t const i3)       {
            return view3D<T      ,int_t>(_data + _index(i3,0,0,0), _s2, _s1, _s0); }

  T* data() const { return _data; }
  
  // for backwards compatability with older versions
  // where view4D v(n3,n2,n1,n0); v(i3,i2,i1,i0) --> ((i3*n2 + i2)*n1 + i1)*n0 + i0;
  //        then _s3 == n2*n1*n0, _s2 == n1*n0, _s1 == n0, _s0 == 1
  size_t  stride() const {
      assert(is_ordered()); assert(_s1 >= 0);
      return _s1; // == n0
  } // stride

  size_t  dim1()   const {
      assert(is_ordered()); assert(_s2 >= 0); assert(_s1 > 0);
      return _s2/_s1; // == n1
  } // dim1

  size_t  dim2()   const {
      assert(is_ordered()); assert(_s3 >= 0); assert(_s2 > 0);
      return _s3/_s2; // == n2
  } // dim1

  view4D transpose(uint64_t const permutation=0x201) // e.g. 0x210 returns an exact view like *this
#ifndef WRITE_ACCESS_ON_TRANSPOSE
      const
#endif // WRITE_ACCESS_ON_TRANSPOSE
  {
      int_t const strides[Rank] = {_s0, _s1, _s2, _s3};  // extend here for higher Rank
      uint8_t rank_count[Rank] = {0, 0, 0, 0};           // extend here for higher Rank
      int_t sp[Rank]; // permuted strides
//    uint8_t pr[Rank]; // permuted rank
      bool non_transposed{true};
      auto perm{permutation}; // non-const copy
      for(int r = 0; r < Rank; ++r) {
          uint8_t const permuted_rank = perm & 0xf;
          perm >>= 4; // shift 4 bits out (one hexadecimal digit)
          assert(permuted_rank < Rank);
          non_transposed &= (permuted_rank == r);
          ++rank_count[permuted_rank];
          auto const permuted_stride = strides[permuted_rank];
//        pr[r] = permuted_rank;
          sp[r] = permuted_stride;
      } // r
      if (non_transposed) {
          debug_printf("# view4D transpose with non-transposed permutation %4.4x\n", permutation);
      } // warn
      int non_standard{0};
      for(int r = 0; r < Rank; ++r) {
          non_standard += (1 != rank_count[r]);
      } // r
      if (non_standard) {
          debug_printf("# view4D transpose with non-standard permutation %4.4x\n", permutation);
      } // warn

      return view4D(_data, sp[3], sp[2], sp[1], sp[0]);    // extend here for higher Rank
  } // transpose

private:
  
  bool is_ordered() const { return (_s0 <= _s1 && _s1 <= _s2 && _s2 <= _s3); }
  
  // private data members
  T * _data;
  int_t _s0, _s1, _s2, _s3;
  size_t _mem; // only > 0 if memory owner
#ifdef  UpperBoundsCheck
  int64_t _n[Rank];
#endif // UpperBoundsCheck

}; // view4D

#undef  debug_printf















namespace data_view {

#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <typename real_t=double>
  inline status_t test_view1D(int const echo=9) {
      // this test is to see if everything compiles
      status_t stat(0);
      int const n0 = 5;
      if (echo > 0) printf("\n# %s(%d)\n", __func__, n0);

      auto a = view1D<real_t>(n0, 3.14); // test memory allocating constructor
      for(int i0 = 0; i0 < n0; ++i0) {
          if (echo > 7) printf("# a1D(%i) = %g\n", i0, a(i0)); // test const & ()
          a(i0) = i0 + .1; // test & ()
          if (echo > 7) printf("# a1D(%i) = %g\n", i0, a[i0]); // test const & []
          a[i0] *= 1.1; // test & []
          if (echo > 7) printf("# a1D(%i) = %g\n", i0, a(i0));
      } // test i0

      view1D<real_t> b; // test default constructor
      b = view1D<real_t>(a.data(), 2); // test data view constructor & move assignment
      for(int i0 = 0; 2*i0 < n0; ++i0) {
          stat += (a(2*i0) != b[i0]); // test const & () and []
          b(i0) = i0 + .7;
          if (echo > 7) printf("# a1D(%i) = %g\n", 2*i0, a(2*i0));
          stat += (a[2*i0] != b(i0));
      } // test i0

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view1D


  template <typename real_t=double>
  inline status_t test_view2D(int const echo=9) {
      status_t stat(0);
      int constexpr n1 = 5, n0 = 3;
      if (echo > 0) printf("\n# %s(%d,%d)\n", __func__, n0, n1);

      view2D<real_t> a(n1, n0, 3.14); // test memory allocating constructor
      for(int i1 = 0; i1 < n1; ++i1) {
          auto a_i0 = a[i1]; // create a view1D
          for(int i0 = 0; i0 < n0; ++i0) {
              real_t const v = i1 + .1*i0;
              a(i1,i0) = v;
              if (echo > 0) printf("# a2D(%i,%i) = %g\n", i1, i0, a(i1,i0));
              stat += (v != a[i1][i0]);
              stat += (v != a[i1](i0));
              stat += (v != a(i1,i0));
              stat += (v != a_i0[i0]);
              stat += (v != a_i0(i0));
              a_i0[i0] = -3.14;   // can modify via the subview?
              a_i0(i0) = -3.1415; // can modify via the subview?
          } // i0
      } // i1

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view2D


  template <typename real_t=double>
  inline status_t test_view3D(int const echo=9) {
      status_t stat(0);
      int constexpr n2 = 5, n1 = 3, n0 = 2;
      if (echo > 0) printf("\n# %s(%d,%d,%d)\n", __func__, n2, n1, n0);

      view3D<real_t> a(n2, n1, n0, 3.14); // test memory allocating constructor
      auto const a_t021 = a.transpose(); // test the default transposition a(i0,j0,i2) --> a_transposed(i0,i2,j0)
      for(int i2 = 0; i2 < n2; ++i2) {
          auto a_i2 = a[i2]; // create a view2D
          for(int i1 = 0; i1 < n1; ++i1) {
              auto a_i2_i1 = a_i2[i1]; // create a view1D
              for(int i0 = 0; i0 < n0; ++i0) {
                  real_t const v = i2 + .1*i1 + .01*i0;
                  a(i2,i1,i0) = v;
                  if (echo > 0) printf("# a3D(%i,%i,%i) = %g\n", i2, i1, i0, a(i2,i1,i0));
                  stat += (v != a(i2,i1,i0)); // value arrived
                  stat += (v != a_i2(i1,i0));
                  stat += (v != a_i2_i1(i0));
                  stat += (v != a_i2_i1[i0]);
                  a_i2(i1,i0) = -3.14;   // can modify via the subview?
                  a_i2_i1(i0) = -3.1415; // can modify via the subview?
                  a_i2_i1[i0] = -3.14;   // can modify via the subview?
              } // i0
          } // i1
      } // i2

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view3D

  
  template <typename real_t=double>
  inline status_t test_view4D(int const echo=9) {
      status_t stat(0);
      int constexpr n3 = 4, n2 = 5, n1 = 3, n0 = 2;
      if (echo > 0) printf("\n# %s(%d,%d,%d,%d)\n", __func__, n3, n2, n1, n0);

      view4D<real_t> a(n3, n2, n1, n0, 3.14); // test memory allocating constructor
      for(int i3 = 0; i3 < n3; ++i3) {
          auto a_i3 = a[i3]; // create a view3D subview
          for(int i2 = 0; i2 < n2; ++i2) {
              for(int i1 = 0; i1 < n1; ++i1) {
                  for(int i0 = 0; i0 < n0; ++i0) {
                      real_t const v = i3 + .1*i2 + .01*i1 + .001*i0;
                      a(i3,i2,i1,i0) = v;
                      if (echo > 0) printf("# a4D(%i,%i,%i,%i) = %g\n", i3, i2, i1, i0, a(i3,i2,i1,i0));
                      stat += (v != a(i3,i2,i1,i0)); // value arrived
                      stat += (v != a_i3(i2,i1,i0));
                  } // i0
              } // i1
          } // i2
      } // i3

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view4D
 
  
  
  template <typename real_t=double>
  inline status_t test_view(int const echo=0) {
      status_t stat = 0;
      if (echo > 0) printf("\n\n# %s<%s>\n\n", __func__,
            (4 == sizeof(real_t)) ? "float" : "double");

      stat += test_view1D<real_t>(echo);
      stat += test_view2D<real_t>(echo);
      stat += test_view3D<real_t>(echo);
      stat += test_view4D<real_t>(echo);

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view

  
  template <typename integer_t=int>
  inline status_t test_view2D_transpose(int const echo=9) {
      status_t stat(0);
      int constexpr n1 = 3, n0 = 5;
      if (echo > 0) printf("\n# %s(%d,%d)\n", __func__, n1, n0);

      view2D<integer_t> a(n1, n0, -3); // test memory allocating constructor
      auto const a_t01 = a.transpose(); // only the view is transposed, not the data, default_transpose
#define TEST_SPECIAL_TRANSPOSE
#ifdef  TEST_SPECIAL_TRANSPOSE
      auto const a_t10 = a.transpose(0x10); // non_transposed
      auto const a_t00 = a.transpose(0x00); // special transpose, check if warnings are generated
      auto const a_t11 = a.transpose(0x11); // special transpose, check if warnings are generated
#endif // TEST_SPECIAL_TRANSPOSE
      for(int i1 = 0; i1 < n1; ++i1) {
          for(int i0 = 0; i0 < n0; ++i0) {
              integer_t const v = i1*10 + i0;
              a(i1,i0) = v;
              if (echo > 0) printf("# a2D(%i,%i) =%3d\n", i1, i0, a(i1,i0));
              stat += (v != a(i1,i0));
              stat += (v != a_t01(i0,i1));

#ifdef  TEST_SPECIAL_TRANSPOSE
              stat += (v != a_t10(i1,i0));
              stat += (v != a_t10(i1,i0));
#endif // TEST_SPECIAL_TRANSPOSE
              
#ifdef  WRITE_ACCESS_ON_TRANSPOSE
              a_t01(i0,i1) = -v;
              stat += (v != -a(i0,i1));
#endif // WRITE_ACCESS_ON_TRANSPOSE

          } // i0
      } // i1

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view2D_transpose

  
  template <typename integer_t=int>
  inline status_t test_view3D_transpose(int const echo=9) {
      status_t stat(0);
      int constexpr n2 = 2, n1 = 3, n0 = 5;
      if (echo > 0) printf("\n# %s(%d,%d,%d)\n", __func__, n2, n1, n0);
      
//    int64_t const perm[6] = {0x012, 0x021, 0x102, 0x120, 0x201, 0x210};

      view3D<integer_t> a(n2, n1, n0, -3); // test memory allocating constructor
      auto const a_t012 = a.transpose(0x012);
      auto const a_t021 = a.transpose(0x021);
      auto const a_t102 = a.transpose(0x102);
      auto const a_t120 = a.transpose(0x120);
      auto       a_t201 = a.transpose(); // default_transpose
#ifdef  TEST_SPECIAL_TRANSPOSE
      auto const a_t210 = a.transpose(0x210); // non_transposed
      auto const a_t000 = a.transpose(0x000); // special transpose, check if warnings are generated
      auto const a_t111 = a.transpose(0x111); // special transpose, check if warnings are generated
#endif // TEST_SPECIAL_TRANSPOSE
      for(int i2 = 0; i2 < n2; ++i2) {
          for(int i1 = 0; i1 < n1; ++i1) {
              for(int i0 = 0; i0 < n0; ++i0) {
                  integer_t const v = i2*100 + i1*10 + i0;
                  a(i2,i1,i0) = v;
                  if (echo > 0) printf("# a3D(%i,%i,%i) =%4d\n", i2, i1, i0, a(i2,i1,i0));
                  stat += (v != a(i2,i1,i0)); // verify that the writing worked
                  stat += (v != a_t012(i0,i1,i2));
                  stat += (v != a_t021(i0,i2,i1));
                  stat += (v != a_t102(i1,i0,i2));
                  stat += (v != a_t120(i1,i2,i0));
                  stat += (v != a_t201(i2,i0,i1));
#ifdef  TEST_SPECIAL_TRANSPOSE
                  stat += (v != a_t210(i2,i1,i0));
#endif // TEST_SPECIAL_TRANSPOSE

#ifdef  WRITE_ACCESS_ON_TRANSPOSE
                  a_t201(i2,i0,i1) = -v;
                  stat += (v != -a_t201(i2,i0,i1));
#endif // WRITE_ACCESS_ON_TRANSPOSE

              } // i0
          } // i1
      } // i2

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view3D_transpose
  

  template <typename integer_t=int>
  inline status_t test_view4D_transpose(int const echo=9) {
      status_t stat(0);
      int constexpr n3 = 4, n2 = 2, n1 = 3, n0 = 5;
      if (echo > 0) printf("\n# %s(%d,%d,%d,%d)\n", __func__, n3, n2, n1, n0);

      view4D<integer_t> a(n3, n2, n1, n0, -3); // test memory allocating constructor

      view4D<integer_t> a_trans[24]; // test the default constructor
      int8_t perm[24][4]; // index permutations
      uint16_t perm_hex[24]; // permutation key
      { // scope: generate all 24 permutations
          int i24{0};
          bool b[4] = {true, true, true, true};
          for(int i3 = 0; i3 < 4; ++i3) {              b[i3] = false;
          for(int i2 = 0; i2 < 4; ++i2) { if (b[i2]) { b[i2] = false;
          for(int i1 = 0; i1 < 4; ++i1) { if (b[i1]) { b[i1] = false;
          for(int i0 = 0; i0 < 4; ++i0) { if (b[i0])
          {
              perm[i24][0] = i0;
              perm[i24][1] = i1;
              perm[i24][2] = i2;
              perm[i24][3] = i3;
              assert(i0 != i3); assert(i0 != i2); assert(i0 != i1);
              assert(i1 != i3); assert(i1 != i2);
              assert(i2 != i3);
              uint64_t const permutation = ((i3*16 + i2)*16 + i1)*16 + i0; // convert to hexadecimal digits
              a_trans[i24] = a.transpose(permutation);
              perm_hex[i24] = permutation;
              ++i24;
          }              } // i0 
          b[i1] = true; }} // i1
          b[i2] = true; }} // i2 
          b[i3] = true;  } // i3
          assert(24 == i24);
      } // scope

      for(int i3 = 0; i3 < n3; ++i3) {
          for(int i2 = 0; i2 < n2; ++i2) {
              for(int i1 = 0; i1 < n1; ++i1) {
                  for(int i0 = 0; i0 < n0; ++i0) {
                      integer_t const v = i3*1000 + i2*100 + i1*10 + i0;
                      a(i3,i2,i1,i0) = v;
                      if (echo > 0) printf("# a4D(%i,%i,%i,%i) =%5d\n", i3, i2, i1, i0, a(i3,i2,i1,i0));
                      stat += (v != a(i3,i2,i1,i0)); // verify that the writing worked

                      int const iv[4] = {i0, i1, i2, i3}; // as vector
                      for(int i24 = 0; i24 < 24; ++i24) {
                          auto const p = perm[i24]; // abbreviate
                          int const jv[4] = {iv[p[0]], iv[p[1]], iv[p[2]], iv[p[3]]};
                          auto const vp = a_trans[i24](jv[3],jv[2],jv[1],jv[0]);
                          int const error = (v != vp);
                          stat += error;
                          if (error) std::printf("# %s: a4D(%i,%i,%i,%i)= %d but a4D.transpose(0x%4.4x)(%i,%i,%i,%i)= %d\n",
                                __func__, i3,i2,i1,i0, v, perm_hex[i24], jv[3],jv[2],jv[1],jv[0], vp);
                      } // i24
                      
                  } // i0
              } // i1
          } // i2
      } // i3

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view3D_transpose
  

  template <typename integer_t=int>
  inline status_t test_view_transpose(int const echo=0) {
      status_t stat(0);
      if (echo > 0) printf("\n# %s\n", __func__);

      stat += test_view2D_transpose<integer_t>(echo);
      stat += test_view3D_transpose<integer_t>(echo);
      stat += test_view4D_transpose<integer_t>(echo);
      
      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // test_view_transpose


  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_view<double>(echo);
      stat += test_view<float> (echo);
      stat += test_view_transpose(echo);

      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace data_view
