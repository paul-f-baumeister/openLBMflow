#pragma once

template <typename T>
T* get_memory(size_t const n, T const init_value=T(0)) {
    auto p = new T[n];
    assert(p);
    for(int i = 0; i < n; ++i) {
        p[i] = init_value;
    } // i
    return p;
} // get_memory

template <typename T>
int free_memory(T* const & p) {
    if (nullptr == p) return 1;
    delete[] p;
    return 0;
} // free_memory
