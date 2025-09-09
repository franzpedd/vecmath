#ifndef VECMATH_TYPES_INCLUDED
#define VECMATH_TYPES_INCLUDED

/// @brief boolean type for checking true/false
typedef signed char vecbool;
#define vec_true 1
#define vec_false 0

/// @brief single precision unions
typedef union float2_t {
    struct { float x, y; } xy;
    struct { float u, v; } uv;
    float data[2];
} float2;

typedef union float3_t {
    struct { float x, y, z; } xyz;
    struct { float r, g, b; } rgb;
    float data[3];
} float3;

typedef union float4_t {
    struct { float x, y, z, w; } xyzw;
    struct { float r, g, b, a; } rgba;
    float data[4];
} float4;

typedef union fmat2_t {
    float data[2][2];
    struct { float2 c0; float2 c1; } column;
    struct { 
        float m00, m01; 
        float m10, m11; 
    } matrix;
} fmat2;

typedef union fmat3_t {
    float data[3][3];
    struct { float3 c0; float3 c1; float3 c2 } column;
    struct { 
        float m00, m01, m02; 
        float m10, m11, m12; 
        float m20, m21, m22; 
    } matrix;
} fmat3;

typedef union fmat4_t {
    float data[4][4];
    struct { float4 c0; float4 c1; float4 c2; float4 c3 } column;
    struct { 
        float m00, m01, m02, m03; 
        float m10, m11, m12, m13; 
        float m20, m21, m22, m23; 
        float m30, m31, m32, m33; 
    } matrix;
} fmat4;

typedef union fquat_t {
	float data[4];
	struct { float x, y, z, w; } vector;
} fquat;

/// @brief double precision unions
typedef union double2_t {
    struct { double x, y; } xy;
    struct { double u, v; } uv;
    double data[2];
} double2;

typedef union double3_t {
    struct { double x, y, z; } xyz;
    struct { double r, g, b; } rgb;
    double data[3];
} double3;

typedef union double4_t {
    struct { double x, y, z, w; } xyzw;
    struct { double r, g, b, a; } rgba;
    double data[4];
} double4;

typedef union dmat2_t {
    double data[2][2];
    struct { double2 c0; double2 c1; } column;
    struct { 
        double m00, m01; 
        double m10, m11; 
    } matrix;
} dmat2;

typedef union dmat3_t {
    double data[3][3];
    struct { double3 c0; double3 c1; double3 c2 } column;
    struct { 
        double m00, m01, m02; 
        double m10, m11, m12; 
        double m20, m21, m22; 
    } matrix;
} dmat3;

typedef union dmat4_t {
    double data[4][4];
    struct { double4 c0; double4 c1; double4 c2; double4 c3 } column;
    struct { 
        double m00, m01, m02, m03; 
        double m10, m11, m12, m13; 
        double m20, m21, m22, m23; 
        double m30, m31, m32, m33; 
    } matrix;
} dmat4;

typedef union dquat_t {
	double data[4];
	struct { double x, y, z, w; } vector;
} dquat;

#endif // VECMATH_INCLUDED