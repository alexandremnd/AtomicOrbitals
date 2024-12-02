#include "include/Atom/atom.hpp"
#include <iostream>

#ifdef __AVX__
    #include <immintrin.h>
#else
    #warning AVX is not available. Code will not compile
#endif

#ifdef __SSE2__
    #include <emmintrin.h>
#else
    #warning SSE is not available
#endif