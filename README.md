# vecmath
vecmath is a C math library with vector and matrix operations with no SIMD operations yet, made to be used on graphics projects (the reason for yet). 

All functions are available both in single and double precision variants. 

## types
* vecbool: used for boolean checking, use vec_true/vec_false to check.
* float2, float3, float4, fmat2, fmat3, fmat4, fquat: single-precision structs;
* double2, double3, double4, dmat2, dmat3, dmat4, dquat : double-precision structs;
## functions
* **basic**: add, sub, mul, div, equals, aprox_equals;
* **vector**: mul_mat, scalar, length, length_sqrt, normalize, dot, cross, scale, lerp, distance, distance_sqrt, reflect, project;
* **matrix**: identity, mul_vec, transpose, determinant, inverse;
    * **row/col major**: get_translation, get_rotation, get_scale, decompose, translate, rotate, scale;
    * **vulkan/opengl/directx**: lookat, perspective, orthographic;
* **quaternion**: identity, mul, length, conjugate, normalize, dot, lerp, slerp, from_euler, to_euler;
    * **row/col major**: to_mat;
* **ray**: 
    * **vulkan**: from_screen_point, to_world_point;
* **utils**: to_radians, to_degrees, cos, sin, tan, max, min, power, log2, floor;

## todo
* **general**: SIMD instruction;
* **funcs**: ray_sphere_intersect, ray_triangle_intersect, ray_square_intersect, aabb_intersects_aabb, rgb_to_hsv, hsv_to_rgb, srgb_to_linear, linear_to_srgb, apply_gamma, perlin_noise, simplex_noise, extract_frustum, sphere_in_frustum, bezier_curve, catmull_rom_spline;

## build
* static: vecmath doesn't need any library, therefore just include it's header and source files in your project and start using it.
* dynamic: same as static, but define VECMATH_BUILD_SHARED on your project's setting, the visibility will be handled by the api.
* header-only: use the files on 'headeronly' folder, set "#define VECMATH_IMPLEMENTATION" on exactly ONE source file before including vecmath.h, don't compile vecmath.c but it must be present alongside vecmath.h