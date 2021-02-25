#pragma once

#include <cassert> // assert
#include <cmath> // std::abs
#include <cstdlib> // std::rand, RAND_MAX
#include <cstdio> // std::printf

namespace heads_and_tails {

  union h_t { double d; int32_t i[2]; };
  
  inline float head(double const x) { 
      h_t u;
      u.d = x;
      u.i[0] &= 0xe0000000; // switch off the least significant 29 bits of the tail
      // this is different from a simple conversion double --> float
      auto const f = float(u.d);
      // ensures that the exponent ranges are ok
      assert( std::abs(x - f) <= 1e-6*std::abs(x) ); 
      return f;
  } // head

  inline int32_t tail(double const x) {
      h_t u;
      u.d = x;
      return u.i[0];
  } // tail

  static constexpr bool CheckMatchingBits = true;

  inline double join(float const head, int32_t const tail) {

      if (CheckMatchingBits) {
          // the IEEE 754 standard forsees for sign, exponent, mantissa
          //      1 + 11 + 52 bits in fp64 (aka double) and
          //      1 +  8 + 23 bits in fp32 (aka float)
          // therefore, converting a double into head and tail
          // requires to truncate 29 bits of the 52-bit mantissa.
          // 
          // When we join head+tail into a double again,
          // the most significant 3 bits of tail should
          // be equal to the least significant 3 bits of head
          // if head and tail were generated from the same double.
          union { float f; uint32_t i; } u4;
          u4.f = head;
          auto const check_head = u4.i & 0x7; // least significant 3 bits
          // check that the most significant 3 bits of tail match with check_head
          auto const check_tail = uint32_t(tail) >> 29;
//               std::printf("# %s: head= %.8e = 0%o tail= 0x%8.8x check_tail= 0%o check_head= 0%o\n",
//                             __func__, head, u.i, tail, check_tail, check_head);
          assert( check_tail == check_head );
      } // check

      h_t u;
      u.d = double(head);
      u.i[0] = tail;
      return u.d;
  } // join

  
  
  
  
  // same functionality collected in a class
  union head_and_tail_union {
  private:
      double x_;
      int32_t ht_[2];
      // since this is a union, the 8 Byte of x_
      // are physically the same memory with the 2x4 Byte of ht_
  public:

      head_and_tail_union(double const x=0) : x_(x) {}
      head_and_tail_union(float const head, int32_t const tail) 
        : x_(double(head))
      {
          if (CheckMatchingBits) {
              union { float f; int32_t i; } u; u.f = head;
              auto const check_head = u.i & 0x7; // least significant 3 bits
              // check that the most significant 3 bits of tail match with check_head
              auto const check_tail = uint32_t(tail) >> 29;
//               std::printf("# %s: head= %.8e = 0%o tail= 0x%8.8x check_tail= 0%o check_head= 0%o\n",
//                             __func__, head, u.i, tail, check_tail, check_head);
              assert( check_tail == check_head );
          } // check
          ht_[0] = tail;
      } // joining constructor

      int32_t tail() const { return ht_[0]; }

      float head() const { 
          union { double f; uint64_t i; } u;
          u.f = x_;
          u.i &= 0xffffffffe0000000; // switch off the least significant 29 bits 
          double const x = u.f;
          // ensures that the exponent ranges are ok
          assert( std::abs(x_ - float(x)) <= 1e-6*std::abs(x_) ); 
          return float(x);
      } // head

      double operator()() const { return x_; }

  }; // head_and_tail_union

  
  inline int all_tests(int const echo=0) {
      int stat(0);
      for(int i = 0; i < (1ul << 14); ++i) {
          auto const x = -0.5 + rand()/double(RAND_MAX);
          
//           head_and_tail_union const ht(x);
//           auto const h = ht.head();
//           auto const t = ht.tail();
//           head_and_tail_union const jn(h, t);
//           auto const j = jn();

          auto const h = head(x);
          auto const t = tail(x);
          auto const j = join(h, t);
          
          if (x != j) {
              std::printf("# %s: x= %.16e\t\thead= %.8e\ttail= %d\n#         join= %.16e\n", 
                            __func__, x, h, t, j);
          } // deviation
          stat += (x != j);
      } // i
      return stat;
  } // test

} // namespace heads_and_tails
