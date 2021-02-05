#pragma once

#include <cstdint> // uint8_t

class CellInfo {
private:
    static int constexpr SOLID = 0x1;
public:

    // default constructor
    CellInfo(uint8_t const bits=0) : bits_(bits) {}
    
    bool const is_solid()  const { return bits_ & SOLID; }
    bool const is_liquid() const { return !(is_solid()); }
    void set_solid() { bits_ |= SOLID; }

    uint8_t data() const { return bits_; }
  
private:
    uint8_t bits_;
}; // CellInfo

