#ifndef VECMATH_DEFINES_INCLUDED
#define VECMATH_DEFINES_INCLUDED

/// @brief used macro definitions
#define VECMATH_EPSILON_FZERO 1e-6f
#define VECMATH_EPSILON_DZERO 1e-12
#define VECMATH_EPSILON_PI 3.14159265358979323846
#define VECMATH_EPSILON_FLT 1.192092896e-07F
#define VECMATH_EPSILON_INT_MIN (-2147483647 - 1)

/// @brief compilation options
#if defined(VECMATH_REQUESTING_HEADER_ONLY)
    #define VECMATH_API static inline // Header-only version - use static inline
    #undef VECMATH_BUILD_SHARED // header-only overrides shared library build
#elif defined(VECMATH_BUILD_SHARED) // Shared library build
    #if defined(_WIN32) || defined(_WIN64)
        #if defined(VECMATH_EXPORTS)
            #define VECMATH_API __declspec(dllexport)
        #else
            #define VECMATH_API __declspec(dllimport)
        #endif
    #elif defined(__linux__) && !defined(__ANDROID__)
        #define VECMATH_API __attribute__((visibility("default")))
    #else
        #define VECMATH_API
    #endif
#else
    #define VECMATH_API // static library
#endif

#endif // VECMATH_DEFINES_INCLUDED
