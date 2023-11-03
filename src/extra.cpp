#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include "iostream"

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
//int comparePrimitivesX(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
//{
//    return computePrimitiveCentroid(p1)[0] < computePrimitiveCentroid(p2)[0];
//}
//int comparePrimitivesY(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
//{
//    return computePrimitiveCentroid(p1)[1] < computePrimitiveCentroid(p2)[1];
//}
//int comparePrimitivesZ(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
//{
//    return computePrimitiveCentroid(p1)[2] < computePrimitiveCentroid(p2)[2];
//}
//glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
//{
//    float x_center = (primitive.v0.position[0] + primitive.v1.position[0] + primitive.v2.position[0]) / 3.f;
//    float y_center = (primitive.v0.position[1] + primitive.v1.position[1] + primitive.v2.position[1]) / 3.f;
//    float z_center = (primitive.v0.position[2] + primitive.v1.position[2] + primitive.v2.position[2]) / 3.f;
//
//    return glm::vec3 { x_center, y_center, z_center };
//}

float calculateVolumeAABB(AxisAlignedBox aabb) {
    return (aabb.upper[0] - aabb.lower[0]) * (aabb.upper[1] - aabb.lower[1]) * (aabb.upper[2] - aabb.lower[2]);
}


bool comparePrimitivesX(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
{
    return computePrimitiveCentroid(p1)[0] < computePrimitiveCentroid(p2)[0];
}
bool comparePrimitivesY(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
{
    return computePrimitiveCentroid(p1)[1] < computePrimitiveCentroid(p2)[1];
}
bool comparePrimitivesZ(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
{
    return computePrimitiveCentroid(p1)[2] < computePrimitiveCentroid(p2)[2];
}

size_t findIndexOfSplit(int bin, int bins, uint32_t axis, const AxisAlignedBox& aabb, const std::span<BVH::Primitive> primitives) {
    /*AxisAlignedBox leftAABB = { aabb.lower,
        { aabb.lower[0] + (aabb.upper[0] - aabb.lower[0]) * 1.f * bin / bins, aabb.upper[1], aabb.upper[2] } }*/
    float delimiter = aabb.lower[axis] + (aabb.upper[axis] - aabb.lower[axis]) * 1.f * bin / bins;
    int count = 1;
    for (int i = 1; i < primitives.size() - 1; i++) {
        if (computePrimitiveCentroid(primitives[i])[axis] <= delimiter)
            count++;
    }

    return count;
}

size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    size_t triangle_count = primitives.size();

    //size_t f_axis;
    size_t triangle_i= - 1;
    float min_cost = FLT_MAX;
    //for (int axis = 0; axis < 3; axis++) {
    /* std::sort(primitives.begin(), primitives.end(), [axis](const Primitive& p1, const Primitive& p2){
        return computePrimitiveCentroid(p1)[axis] < computePrimitiveCentroid(p2)[axis]; 
    });*/
    if (axis == 0)
        std::sort(primitives.begin(), primitives.end(), comparePrimitivesX);
    else if (axis == 1)
        std::sort(primitives.begin(), primitives.end(), comparePrimitivesY);
    else
        std::sort(primitives.begin(), primitives.end(), comparePrimitivesZ);

    //SAH no binning
    //for (int i = 1; i < triangle_count - 1; i++) {
    //    auto aabb_left = computeSpanAABB(primitives.subspan(0, i));
    //    auto aabb_right = computeSpanAABB(primitives.subspan(i, std::dynamic_extent));
    //    float volume_left = calculateVolumeAABB(aabb_left) / calculateVolumeAABB(aabb);
    //    float volume_right = calculateVolumeAABB(aabb_right) / calculateVolumeAABB(aabb);
    //    float cost = ( volume_left * i + volume_right * (triangle_count - i));
    //   
    //    //std::cout << " cost: " << cost << " left: " << volume_left << " right: " << volume_right << std::endl;

    //    if (cost < min_cost) {
    //        //f_axis = axis;
    //        triangle_i = i;
    //        min_cost = cost;
    //    }
    //}
    //}
    int totalBins = 20;
    for (int i = 0; i < totalBins - 1; i++) {
        size_t split = findIndexOfSplit(i + 1, totalBins, axis, aabb, primitives);

        auto span_left = primitives.subspan(0, split);
        auto span_right = primitives.subspan(split, std::dynamic_extent);
        auto aabb_left = computeSpanAABB(span_left);
        auto aabb_right = computeSpanAABB(span_right);
        float volume_left = calculateVolumeAABB(aabb_left) / calculateVolumeAABB(aabb);
        float volume_right = calculateVolumeAABB(aabb_right) / calculateVolumeAABB(aabb);
        float cost = ( volume_left * span_left.size() + volume_right * span_right.size());

        if (cost <= min_cost) {
        //f_axis = axis;
        triangle_i = span_left.size();
        min_cost = cost;
        }
    }


    return triangle_i; // This is clearly not the solution
}

