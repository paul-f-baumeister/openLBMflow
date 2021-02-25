#pragma once

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // int64_t, int32_t, uint32_t
#include <cstdlib> // std::rand, RAND_MAX

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED

#define debug_printf(...) std::printf(__VA_ARGS__)
// #define debug_printf(...)


namespace lbm_index {

  template <unsigned Bz, unsigned By, unsigned Bx, typename base_int_t=uint32_t>
  class MultiIndex {
      // group 4 integers of fixed bit width into one standard integer type
      // the width parameters need to be chose at compile time

  private:
      static unsigned constexpr BitsPerByte = 8;
  public:
      // a flexible entitiy that can hold up to 4 integer numbers (z,y,x,w)
      // with free choice of how many bits may be taken by each component

      static unsigned constexpr xBits = Bx, yBits = By, zBits = Bz;
      // the remaining number of bits
      static int constexpr wBits = BitsPerByte*sizeof(base_int_t) - (Bx + By + Bz);
      
  private:
      static base_int_t constexpr wMax_ = base_int_t(1) << wBits, wMask_ = wMax_ - 1,
                                  xMax_ = base_int_t(1) << xBits, xMask_ = xMax_ - 1,
                                  yMax_ = base_int_t(1) << yBits, yMask_ = yMax_ - 1,
                                  zMax_ = base_int_t(1) << zBits, zMask_ = zMax_ - 1;
  public:
    
      // default constructor
      MultiIndex(base_int_t const bits=0) : bits_(bits) { assert(wBits >= 0); }

      MultiIndex(int const z, int const y, int const x, int const w=0) {
          assert(wBits >= 0);
          base_int_t const wm = w & wMask_,
                           xm = x & xMask_,
                           ym = y & yMask_,
                           zm = z & zMask_;
          assert(zm == z); assert(ym == y); assert(xm == x); assert(wm == w);
          // synthesis
          bits_ = (((((zm << yBits) | ym) << xBits) | xm) << wBits) | wm;
      } // constructor

      inline uint32_t z() const { return (bits_ >> (yBits + xBits + wBits)) & zMask_; }
      inline uint32_t y() const { return (bits_ >> (        xBits + wBits)) & yMask_; }
      inline uint32_t x() const { return (bits_ >> (                wBits)) & xMask_; }
      inline uint32_t w() const { return (bits_ >> (                    0)) & wMask_; }

      base_int_t retrieve(int & iz, int & iy, int & ix) const {
          iz = z(); iy = y(); ix = x();
          return bits_;
      } // retrieve indices back

      base_int_t retrieve(int & iz, int & iy, int & ix, int & iw) const {
          iw = w();
          return retrieve(iz, iy, ix);
      } // retrieve indices back
      
      // these numbers could, in principle, be constants, e.g. making xMax public,
      // but then, we cannot exchange this container by a dynamic bit container
      static base_int_t constexpr w_max() { return wMax_; }
      static base_int_t constexpr x_max() { return xMax_; }
      static base_int_t constexpr y_max() { return yMax_; }
      static base_int_t constexpr z_max() { return zMax_; }

      base_int_t data() const { return bits_; }

  private:
      // member
      base_int_t bits_;

  }; // MultiIndex
  
  
  
  
  
  
  template <unsigned Rank, typename base_int_t=uint32_t, typename index_t=uint32_t>
  class BitCompressor {
    // this class is a vector of unsigned integers
    // of specified bit width that are compressed
    // into a standard integer type
  private:
      static unsigned constexpr BitsPerByte = 8;
      static base_int_t constexpr one = 1;
  public:
      static unsigned constexpr Rank_ = Rank;

      // [default] constructor
      BitCompressor(base_int_t const max_range[Rank]=nullptr) {
          assert(Rank > 0);
          // how many bits can we use?
          int const n_all = sizeof(base_int_t)*BitsPerByte;
          debug_printf("# %s<%d>: distribute %d bits as", __func__, Rank, n_all);
          int const n_each = n_all/Rank; // integer divide for default distribution
          shift_[0] = 0;
          for(int r = 0; r < Rank - 1; ++r) {
              if (nullptr != max_range) {
                  assert(max_range[r] >= 0);
                  nbits_[r] = ceil_log2(max_range[r]);
              } else {
                  nbits_[r] = n_each;
              } // nullptr != max_range
              debug_printf(" %d +", nbits_[r]);
              shift_[r + 1] = shift_[r] + nbits_[r]; // prefix sum
          } // r
          int const n_remaining_bits = n_all - shift_[Rank - 1];
          debug_printf(" %d\n", n_remaining_bits);
          assert(n_remaining_bits >= 0 && "Not enough bits!");
          nbits_[Rank - 1] = n_remaining_bits;

          // prepare masks for decoding
          for(int r = 0; r < Rank; ++r) {
              mask_[r] = (one << nbits_[r]) - one;
          } // r
          if (nullptr != max_range) {
              assert(mask_[Rank - 1] >= max_range[Rank - 1] && "Not enogh bits for highest rank!");
          }
          
      } // [default] constructor

      base_int_t encode(index_t const indices[Rank]) const {
          base_int_t bits{0};
          for(int r = 0; r < Rank; ++r) {
              bits |= ((indices[r] & mask_[r]) << shift_[r]);
          } // r
          return bits;
      } // encode

      void get(index_t indices[Rank], base_int_t const bits) const {
          for(int r = 0; r < Rank; ++r) {
              indices[r] = (bits >> shift_[r]) & mask_[r];
          } // r
      } // get all

      template <unsigned r>
      index_t get(base_int_t const bits) const {
          assert(r < Rank);
          return (bits >> shift_[r]) & mask_[r];
      } // get one

      index_t max(unsigned const r) const { assert(r < Rank); return mask_[r]; }

  private:

      unsigned ceil_log2(index_t const max_range) { //  == int(std::ceil(std::log2(max_range)))
          unsigned log2{0};
          index_t mrange{max_range}; // a non-const copy
          while (mrange > 0) {
              mrange >>= 1; // divide by 2
              ++log2;
          } // while mrange > 0
          return log2;
      } // ceil_log2

  private:
      // members
      index_t  mask_[Rank];
      uint8_t nbits_[Rank];
      uint8_t shift_[Rank];
      
  }; // BitCompressor
  
  
#undef debug_printf
  
  
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  template <class Index_t>
  inline status_t test_range(int const echo=9) {
      // this test is to see if everything compiles
      status_t stat(0);
  
      int const nB[4] = {Index_t::wBits, Index_t::xBits, Index_t::yBits, Index_t::zBits};
      int const total = nB[0] + nB[1] + nB[2] + nB[3];
      if (echo > 0) {
          std::printf("# %s:  MultiIndex<%d,%d,%d>\tholds %d + %d + %d + %d \t= %d bit = %ld Byte\n", 
              __func__, nB[3], nB[2], nB[1],  nB[3], nB[2], nB[1], nB[0],   total, sizeof(Index_t));
      } // echo
      assert(sizeof(Index_t)*8 == total);

      size_t number_of_checks{0};
      for(int64_t z = 0; z < Index_t::z_max(); z = z*2+1) {
      for(int64_t y = 0; y < Index_t::y_max(); y = y*2+1) {
      for(int64_t x = 0; x < Index_t::x_max(); x = x*2+1) {
      for(int64_t w = 0; w < Index_t::w_max(); w = w*2+1) {
          Index_t const idx(z, y, x, w);
          int iv[4];
          auto const idx_data = idx.retrieve(iv[3], iv[2], iv[1], iv[0]);
          stat += (iv[3] != z);
          stat += (iv[2] != y);
          stat += (iv[1] != x);
          stat += (iv[0] != w);
          Index_t const jdx(iv[3], iv[2], iv[1], iv[0]);
          stat += (idx.data() != idx_data);
          stat += (idx.data() != jdx.data());
          stat += (iv[3] != idx.z());
          stat += (iv[2] != idx.y());
          stat += (iv[1] != idx.x());
          stat += (iv[0] != idx.w());
//           stat += (jdx.get<'z'>() != z);
//           stat += (jdx.get<'y'>() != y);
//           stat += (jdx.get<'x'>() != x);
//           stat += (jdx.get<'w'>() != w);
          ++number_of_checks;
      }}}} // w x y z

      if (number_of_checks < 1) {
          if (echo > 0) std::printf("# %s:  MultiIndex<%d,%d,%d> no checks performed\n", __func__, nB[3], nB[2], nB[1]);
          return stat;
      } // echo
      if (int(stat) > 0 && echo > 0) {
          std::printf("# %s:  MultiIndex<%d,%d,%d> %d errors found!\n", __func__, nB[3], nB[2], nB[1], int(stat));
      } // echo
      if (echo > 7) {
          std::printf("# %s:  MultiIndex<%d,%d,%d> %ld checks performed, %d errors found\n",
                        __func__, nB[3], nB[2], nB[1], number_of_checks, int(stat));
      } // echo

      return stat;
  } // test_range

  
  inline status_t test_basic(int const echo=9) {
      // this test is to see if everything compiles
      status_t stat(0);

      stat += test_range<MultiIndex<9,9,9>>(echo); // a typical example
//    stat += test_range<MultiIndex<32,0,0>>(echo); // does not compile, needs at least one bit in a different rank
//    stat += test_range<MultiIndex<31,1,1>>(echo); // does not compile, too many bits required
//    stat += test_range<MultiIndex<33,1,1,int64_t>>(echo); // fails because the return type is uint32_t which cannot host 33 bits
      stat += test_range<MultiIndex<31,1,0>>(echo);
      stat += test_range<MultiIndex<31,0,1>>(echo);
      stat += test_range<MultiIndex<31,0,0>>(echo);

      stat += test_range<MultiIndex<1,21,21,int64_t>>(echo); // equivalent to global coordinates but not bit-interleaved
      stat += test_range<MultiIndex<16,16,16,int64_t>>(echo); // equivalent to union { int64_t; uint16_t[4]; };
      stat += test_range<MultiIndex<4,2,1,int8_t>>(echo); //
      stat += test_range<MultiIndex<2,2,2,int8_t>>(echo); //

      return stat;
  } // test_basic

  
  
  inline status_t test_BitCompressor(int const echo=9) {
      // this test is to see if everything compiles
      status_t stat(0);

      // test default constructors
//    BitCompressor<0> bc0; // ? --> compiler needs zero-length-array
      BitCompressor<1> bc1; // 32
      BitCompressor<2> bc2; // 16 + 16
      BitCompressor<3> bc3; // 10 + 10 + 12
      BitCompressor<4> bc4; // 8 + 8 + 8 + 8
      BitCompressor<5> bc5; // 6 + 6 + 6 + 6 + 8
      BitCompressor<6> bc6; // 5 + 5 + 5 + 5 + 5 + 7
      BitCompressor<7> bc7; // 4 + 4 + 4 + 4 + 4 + 4 + 8
      BitCompressor<8> bc8; // 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4
//    BitCompressor<33> bc33; // 0 + ... + 0 + 32 (not what we want)
      BitCompressor<3,int64_t> gc3; // 21 + 21 + 22 // equivalent to global coordinates but without interleaved bits

      // test user constructors:
      uint32_t const b9995[] = {511, 511, 511, 31}; // typical 3D + q-speeds
      BitCompressor<4> bc9995(b9995);
      for(int r = 0; r < bc9995.Rank_; ++r) assert(bc9995.max(r) == b9995[r]);

      uint32_t const b13136[] = {8191, 8191, 63}; // typical 2D + q-speeds
      BitCompressor<3> bc13136(b13136);
      for(int r = 0; r < bc13136.Rank_; ++r) assert(bc13136.max(r) == b13136[r]);
      
      uint32_t const b111110[] = {2047, 2047, 1023}; // typical 3D with potentially shorter z-direction
      BitCompressor<3> bc111110(b111110);
      for(int r = 0; r < bc111110.Rank_; ++r) assert(bc111110.max(r) == b111110[r]);

//       uint32_t const b1617[] = {65535, 131071}; // this 2D example must fail!
//       BitCompressor<2> bc1617(b1617);
//       for(int r = 0; r < bc1617.Rank_; ++r) assert(bc1617.max(r) == b1617[r]);
      
      int constexpr MaxRank = 8; // increase if we test a compressor with a higher rank

      // test encoding and decoding (get)
      for(int it = 0; it < (1 << 12); ++it) {
        
          uint32_t i[MaxRank], j[MaxRank], rn[MaxRank]; 
          for(int r = 0; r < MaxRank; ++r) rn[r] = std::rand(); // eight random numbers

#define   TestBitCompressorRank(BC)                                        \
          {   for(int r = 0; r < BC.Rank_; ++r) i[r] = rn[r] & BC.max(r);  \
              BC.get(j, BC.encode(i));                                     \
              for(int r = 0; r < BC.Rank_; ++r) stat += (i[r] != j[r]);    \
          }

          TestBitCompressorRank(bc1);
          TestBitCompressorRank(bc2);
          TestBitCompressorRank(bc3);
          TestBitCompressorRank(bc4);
          TestBitCompressorRank(bc5);
          TestBitCompressorRank(bc6);
          TestBitCompressorRank(bc7);
          TestBitCompressorRank(bc8);

          TestBitCompressorRank(bc9995);
          TestBitCompressorRank(bc13136);
          TestBitCompressorRank(bc111110);

      } // i

      return stat;
  } // test_BitCompressor
  
  
  
  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_basic(echo);
      stat += test_BitCompressor(echo);
      if (0 != stat && echo > 0) std::printf("# %s: failed with %d errors!\n", __func__, int(stat));
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace lbm_index
