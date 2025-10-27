#include "vecmath_ray_op.h"
#include "vecmath_vec_op.h"
#include <math.h>

VECMATH_API fray fray_from_screen_point_vulkan(const float2* screenPos, const float2* viewportSize, const fmat4* inverseProj, const fmat4* inverseView)
{
    fray ray = { 0 };
    
    // convert to Vulkan NDC
    float ndcX = (2.0f * screenPos->xy.x) / viewportSize->xy.x - 1.0f;
    float ndcY = 1.0f - (2.0f * screenPos->xy.y) / viewportSize->xy.y; 
    // float ndcY = (2.0f * screenPos->xy.y) / viewportSize->xy.y - 1.0f; // TEST THIS IF DOESNT WORK
    
    // near and far points in NDC (Vulkan: Z range [0,1])
    float4 near = { ndcX, ndcY, 0.0f, 1.0f };
    float4 far = { ndcX, ndcY, 1.0f, 1.0f };
    
    // convert to view space
    float4 nearPointView = float4_mul_fmat4(&near, inverseProj);
    float4 farPointView = float4_mul_fmat4(&far, inverseProj);
    
    // perspective divide
    nearPointView = float4_scalar(&nearPointView, 1.0f / nearPointView.xyzw.w);
    farPointView = float4_scalar(&farPointView, 1.0f / farPointView.xyzw.w);
    
    // convert to world space
    float4 nearPointWorld = float4_mul_fmat4(&nearPointView, inverseView);
    float4 farPointWorld = float4_mul_fmat4(&farPointView, inverseView);
    
    // set ray origin and direction
    ray.origin = (float3){ nearPointWorld.xyzw.x, nearPointWorld.xyzw.y, nearPointWorld.xyzw.z };
    float3 farPos = (float3){ farPointWorld.xyzw.x, farPointWorld.xyzw.y, farPointWorld.xyzw.z };
    float3 subtraction = (float3){farPos.xyz.x - ray.origin.xyz.x, farPos.xyz.y - ray.origin.xyz.y, farPos.xyz.z - ray.origin.xyz.z };
    ray.direction = float3_normalize(&subtraction);
    
    return ray;
}

VECMATH_API float3 fray_get_point(const fray* ray, float distance)
{
    float3 mul = float3_scale(&ray->direction, distance);
    float3 point = (float3){ ray->origin.xyz.x + mul.xyz.x, ray->origin.xyz.y + mul.xyz.y, ray->origin.xyz.z + mul.xyz.z };
    return point;
}

VECMATH_API vecbool fray_sphere_intersection(const fray* ray, const float3* center, float radius, float* outDistance)
{
    float3 oc = (float3){ ray->origin.xyz.x - center->xyz.x, ray->origin.xyz.y - center->xyz.y, ray->origin.xyz.z - center->xyz.z };
    float a = float3_dot(&ray->direction, &ray->direction);
    float b = 2.0f * float3_dot(&oc, &ray->direction);
    float c = float3_dot(&oc, &oc) - radius * radius;
    
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return vec_false;
    
    float sqrtDisc = sqrtf(discriminant);
    float t0 = (-b - sqrtDisc) / (2.0f * a);
    float t1 = (-b + sqrtDisc) / (2.0f * a);

    if (t0 >= 0.0f) {
        *outDistance = t0;
        return vec_true;
    }

    else if (t1 >= 0.0f) {
        *outDistance = t1;
        return vec_true;
    }

    return vec_false;
}

VECMATH_API vecbool fray_quad_intersection(const fray *ray, const float3 *v0, const float3 *v1, const float3 *v2, const float3 *v3, float *outDistance)
{
    // compute normal of the quad (using two edges)
    float3 edge1 = (float3){ v1->xyz.x - v0->xyz.x, v1->xyz.y - v0->xyz.y, v1->xyz.z - v0->xyz.z };
    float3 edge2 = (float3){ v3->xyz.x - v0->xyz.x, v3->xyz.y - v0->xyz.y, v3->xyz.z - v0->xyz.z };
    float3 cross = float3_cross(&edge1, &edge2);
    float3 normal = float3_normalize(&cross);

    // ray-plane intersection
    float denom = float3_dot(&normal, &ray->direction);
    if (fabsf(denom) < VECMATH_EPSILON_FZERO) return vec_false; // ray is parallel to the quad plane

    float3 diff = (float3){ v0->xyz.x - ray->origin.xyz.x, v0->xyz.y - ray->origin.xyz.y, v0->xyz.z - ray->origin.xyz.z };
    float t = float3_dot(&diff, &normal) / denom;
    if (t < 0.0f) return vec_false; // intersection behind the ray

    // compute intersection point
    float3 hitPoint = fray_get_point(ray, t);

    // use barycentric technique for point-in-quad test, split the quad into two triangles (v0,v1,v2) and (v0,v2,v3)
    float3 u, v, w, sub0, sub1;
    float denom1, denom2;

    sub0 = (float3){ v1->xyz.x - v0->xyz.x, v1->xyz.y - v0->xyz.y, v1->xyz.z - v0->xyz.z };
    sub1 = (float3){ hitPoint.xyz.x - v0->xyz.x, hitPoint.xyz.y - v0->xyz.y, hitPoint.xyz.z - v0->xyz.z };
    float3 c0 = float3_cross(&sub0, &sub1);

    sub0 = (float3){ v2->xyz.x - v1->xyz.x, v2->xyz.y - v1->xyz.y, v2->xyz.z - v1->xyz.z };
    sub1 = (float3){ hitPoint.xyz.x - v1->xyz.x, hitPoint.xyz.y - v1->xyz.y, hitPoint.xyz.z - v1->xyz.z };
    float3 c1 = float3_cross(&sub0, &sub1);

    sub0 = (float3){ v3->xyz.x - v2->xyz.x, v3->xyz.y - v2->xyz.y, v3->xyz.z - v2->xyz.z };
    sub1 = (float3){ hitPoint.xyz.x - v2->xyz.x, hitPoint.xyz.y - v2->xyz.y, hitPoint.xyz.z - v2->xyz.z };
    float3 c2 = float3_cross(&sub0, &sub1);
    
    sub0 = (float3){ v0->xyz.x - v3->xyz.x, v0->xyz.y - v3->xyz.y, v0->xyz.z - v3->xyz.z };
    sub1 = (float3){ hitPoint.xyz.x - v3->xyz.x, hitPoint.xyz.y - v3->xyz.y, hitPoint.xyz.z - v3->xyz.z };
    float3 c3 = float3_cross(&sub0, &sub1);

    // all cross products should point roughly in the same direction as the normal
    if (float3_dot(&c0, &normal) < 0.0f) return vec_false;
    if (float3_dot(&c1, &normal) < 0.0f) return vec_false;
    if (float3_dot(&c2, &normal) < 0.0f) return vec_false;
    if (float3_dot(&c3, &normal) < 0.0f) return vec_false;

    *outDistance = t;
    return vec_true;
}

VECMATH_API vecbool fray_aabb_intersection(const fray *ray, const float3 *minBound, const float3 *maxBound, float *outDistance)
{
    float3 invDir = (float3){ 1.0f / ray->direction.xyz.x,  1.0f / ray->direction.xyz.y, 1.0f / ray->direction.xyz.z };

    float tmin = (minBound->xyz.x - ray->origin.xyz.x) * invDir.xyz.x;
    float tmax = (maxBound->xyz.x - ray->origin.xyz.x) * invDir.xyz.x;
    if (tmin > tmax) { 
        float tmp = tmin;
        tmin = tmax; 
        tmax = tmp; 
    }

    float tymin = (minBound->xyz.y - ray->origin.xyz.y) * invDir.xyz.y;
    float tymax = (maxBound->xyz.y - ray->origin.xyz.y) * invDir.xyz.y;
    if (tymin > tymax) { 
        float tmp = tymin;
        tymin = tymax;
        tymax = tmp;
    }

    if ((tmin > tymax) || (tymin > tmax)) return vec_false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    float tzmin = (minBound->xyz.z - ray->origin.xyz.z) * invDir.xyz.z;
    float tzmax = (maxBound->xyz.z - ray->origin.xyz.z) * invDir.xyz.z;
    if (tzmin > tzmax) {
        float tmp = tzmin;
        tzmin = tzmax;
        tzmax = tmp;
    }

    if ((tmin > tzmax) || (tzmin > tmax)) return vec_false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;

    *outDistance = (tmin >= 0.0f) ? tmin : tmax;
    return (*outDistance >= 0.0f) ? vec_true : vec_false;
}

VECMATH_API vecbool fray_obb_intersection(const fray *ray, const float3 *center, const fmat3 *rotation, const float3 *halfExtents, float *outDistance)
{
    // compute R^T (transpose)
    fmat3 rotT;
    rotT.data[0][0] = rotation->data[0][0]; rotT.data[0][1] = rotation->data[1][0]; rotT.data[0][2] = rotation->data[2][0];
    rotT.data[1][0] = rotation->data[0][1]; rotT.data[1][1] = rotation->data[1][1]; rotT.data[1][2] = rotation->data[2][1];
    rotT.data[2][0] = rotation->data[0][2]; rotT.data[2][1] = rotation->data[1][2]; rotT.data[2][2] = rotation->data[2][2];

    // transform ray into OBB local space
    float3 oc = (float3){ ray->origin.xyz.x - center->xyz.x, ray->origin.xyz.y - center->xyz.y, ray->origin.xyz.z - center->xyz.z };

    // origin' = R^T * (origin - center)
    float3 localOrigin = { 
        rotT.data[0][0] * oc.xyz.x + rotT.data[0][1] * oc.xyz.y + rotT.data[0][2] * oc.xyz.z, 
        rotT.data[1][0] * oc.xyz.x + rotT.data[1][1] * oc.xyz.y + rotT.data[1][2] * oc.xyz.z,
        rotT.data[2][0] * oc.xyz.x + rotT.data[2][1] * oc.xyz.y + rotT.data[2][2] * oc.xyz.z };

    // dir' = R^T * dir
    float3 localDir = { 
        rotT.data[0][0] * ray->direction.xyz.x + rotT.data[0][1] * ray->direction.xyz.y + rotT.data[0][2] * ray->direction.xyz.z,
        rotT.data[1][0] * ray->direction.xyz.x + rotT.data[1][1] * ray->direction.xyz.y + rotT.data[1][2] * ray->direction.xyz.z,
        rotT.data[2][0] * ray->direction.xyz.x + rotT.data[2][1] * ray->direction.xyz.y + rotT.data[2][2] * ray->direction.xyz.z
    };

    // build a temporary ray in local space
    fray localRay = { 0 };
    localRay.origin = localOrigin;
    localRay.direction = localDir;

    // intersect with local AABB (-E to +E)
    float3 minB = { -halfExtents->xyz.x, -halfExtents->xyz.y, -halfExtents->xyz.z };
    float3 maxB = {  halfExtents->xyz.x,  halfExtents->xyz.y,  halfExtents->xyz.z };

    return fray_aabb_intersection(&localRay, &minB, &maxB, outDistance);
}
