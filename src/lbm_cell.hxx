#pragma once

#include <cstdint> // uint8_t

class CellInfo {
  // Cell information for a single cubic lattice cell:
  //    in this version: is_solid and is_liquid == !is_solid are offered
  //    but if we want to allow thin walls, we need 6 bits to
  //    encode if a cell has walls on its 6 faces
  
  // furthermore, for parallelization, we have to mark liquid cells with a priority tag
  
private:
    static int constexpr SOLID = 0x1;
    static int constexpr PRIORITY = 0x80;
public:

    // default constructor
    CellInfo(uint8_t const bits=0) : bits_(bits) {}
    
    bool const is_solid()  const { return bits_ & SOLID; }
    bool const is_liquid() const { return !(is_solid()); }
    void set_solid() { bits_ |= SOLID; }

    bool const priority()  const { return bits_ & PRIORITY; }
    void set_priority() { bits_ |= PRIORITY; }

    uint8_t data() const { return bits_; }
  
private:
    uint8_t bits_;
}; // CellInfo

