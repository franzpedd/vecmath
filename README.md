# vecmath
vecmath is a C math library with vector and matrix operations with no SIMD operations yet to be used on graphics projects (the reason for yet).

## types
* vecbool: used for boolean checking, use vec_true/vec_false to check.
* float2, float3, float4, fmat2, fmat3, fmat4, fquat: single-precision unions;
* double2, double3, double4, dmat2, dmat3, dmat4, dquat_t : double-precision unions;
## functions
* basic-operations: add, sub, mul, div, equal, aprox_equals

## todo
* vector funcs: len, normalize, dot, cross, scale, lerp, distance, reflect, project;
* matrix funcs: identity, vec_mul, transpose, determinant, inverse, translation, rotation, scale, lookat, perspective, ortographic, decompose, get_translation, get_rotation, get_scale
* quat funcs: mul, conjugate, normalize, dot, slerp, from_axis_angle, to_mat4
* util funcs: clamp, lerp, deg_to_rad, rad_to_deg, min, max