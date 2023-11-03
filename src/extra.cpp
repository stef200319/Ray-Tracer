#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <texture.h>
#include <algorithm>

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

    float focalLength = features.extra.focalLength;
    float aperture = features.extra.aperture;
    int n = features.extra.raysDoF;

    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };
            auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());

            Ray initialRay = rays.at(0);

            glm::vec3 intersection = initialRay.origin + initialRay.direction * focalLength;

            for (int i = 1; i < n; i++) {
                float sampleX = state.sampler.next_1d();
                float sampleY = state.sampler.next_1d();

                sampleX = sampleX * 2 * aperture - aperture;
                sampleY = sampleY * 2 * aperture - aperture;
                glm::vec3 offset = sampleX * camera.left() + sampleY * camera.up();

                Ray newRay(initialRay.origin + offset, glm::normalize(intersection - initialRay.origin - offset));
                rays.push_back(newRay);
            }
            
            auto L = renderRays(state, rays);
            screen.setPixel(x, y, L);
        }
    }
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

//this is so weird, it becomes 10000 times faster by not using this boxFilter method
// i suppose cache misses, idk
// it was a box filter because I was using it in early testing, in the main code a binomial gradient is used
glm::vec3 boxFilter(const Screen& image, int x, int y) {
    glm::vec3 edgeColor { 0.f, 0.f, 0.f };
    int filterSize = 1;
    auto pixels = image.pixels();

    //std::cout << x << " - " << image.resolution().x << " " << y << " - " << image.resolution().y 
    //    << " " << x * y << " out of: " << image.resolution().x * image.resolution().y << std::endl
    //          << std::endl
    //          << std::endl;
    glm::vec3 sum { 0.f, 0.f, 0.f };
    /*for (int i = -filterSize; i < filterSize + 1; ++i) {
        for (int j = -filterSize; j < filterSize + 1; ++j) {
            int index = image.indexAt(x + j, y + i);
            if (! (index < 0 || index >= pixels.size()) )
                sum += pixels[index];
        }*/
    for (int i = -filterSize; i < filterSize + 1; ++i){
        int index = image.indexAt(x + i, y);
        if (!(index < 0 || index >= pixels.size()))
            sum += pixels[index];
    }

    //sum *= 1.f / ((filterSize * 2 + 1) * (filterSize * 2 + 1));
    sum /= (filterSize * 2 + 1);
    return sum;
}

float choose(int n, int k) {
    if (k == 0)
        return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

std::vector<float> bloomFilter(int size) {
    std::vector<float> res;
    res.reserve(size);

    for (int i = 1; i <= size; i++)
        res.push_back(choose(size, i));

    return res;
}

void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }
    
    Screen brightMask(image.resolution(), true);
    brightMask.clear(glm::vec3 { 0.f, 0.f, 0.f });
    //glm::vec3 tresholdColor { 0.7f, 0.7f, 0.7f };
    auto pixels = image.pixels();
    //initialize mask
    for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {
            auto currPixel = pixels[image.indexAt(x, y)];
            //if (currPixel.x >= tresholdColor.x && currPixel.y >= tresholdColor.y && currPixel.z >= tresholdColor.z)
                 if (currPixel.x + currPixel.y + currPixel.z >= features.extra.bloomTreshold)
                brightMask.setPixel(x, y, currPixel);
            else
                brightMask.setPixel(x, y, glm::vec3 { 0.f, 0.f, 0.f });
        }
    }

    int filterSize = features.extra.bloomFilterSize;
    auto bloom = bloomFilter(filterSize * 2 + 1);
    float bloomDivider = 0;
    for (auto filterV : bloom)
        bloomDivider += filterV;

    Screen horizontalMask(image.resolution(), true);
    horizontalMask.clear(glm::vec3 { 0.f, 0.f, 0.f });
    auto brightMaskP = brightMask.pixels();
    //box filter and scale on mask
    for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {

            //calling another function seems to slow it down by a lot, even tho only references are passed, maybe optimization problems?
            /*auto currPixel = pixels[brightMask.indexAt(x, y)];
            if (currPixel.x != 0.f && currPixel.y != 0.f && currPixel.z != 0.f)*/
            // mask2.setPixel(x, y, boxFilter(mask, x, y));

            glm::vec3 sum { 0.f, 0.f, 0.f };
            for (int i = -filterSize, j = 0; i < filterSize + 1; i++, j++) {
                int index = brightMask.indexAt(x + i, y);
                if (!(index < 0 || index >= pixels.size()))
                     sum += brightMaskP[index] * bloom[j];
                    //sum += brightMaskP[index];
            }

             //sum /= (filterSize * 2 + 1);
            sum /= bloomDivider;
            horizontalMask.setPixel(x, y, sum);
        }
    }

     Screen verticalMask(image.resolution(), true);
     verticalMask.clear(glm::vec3 { 0.f, 0.f, 0.f });
     auto horizontalMaskP = horizontalMask.pixels();
    // box filter and scale on mask
    /* for (int x = 0; x < image.resolution().x; x++) {
         for (int y = 0; y < image.resolution().y; y++) {*/
    for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {

            glm::vec3 sum { 0.f, 0.f, 0.f };
            for (int i = -filterSize, j = 0; i < filterSize + 1; i++, j++) {
                int index = horizontalMask.indexAt(x, y + i);
                if (!(index < 0 || index >= pixels.size()))
                    //sum += horizontalMaskP[index];
                     sum += brightMaskP[index] * bloom[j];
            }

            // sum *= 1.f / ((filterSize * 2 + 1) * (filterSize * 2 + 1));
            //sum /= (filterSize * 2 + 1);
            sum /= bloomDivider;
            verticalMask.setPixel(x, y, sum);
        }
    }


    //add to original
    auto maskPixels = verticalMask.pixels();
     for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {
            auto currPixel = pixels[image.indexAt(x, y)];
            auto maskPixel =maskPixels[image.indexAt(x, y)];

            image.setPixel(x, y, currPixel + maskPixel);
            //image.setPixel(x, y, horizontalMaskP[image.indexAt(x, y)]);
        }
    }
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
         //Part of your implementation should go here

        float x = ray.direction.x;
        float y = ray.direction.y;
        float z = ray.direction.z;
        float u;
        float v;
        float a = 0.3333332f;

        if (x < 0 && abs(x) >= abs(y) && abs(x) >= abs(z)) {
            u = z / abs(x) / 2 + 0.5f;
            v = y / abs(x) / 2 + 0.5f;

            //these are coords within respective face, now we map to correct position of face within the texture

            u = u / 4;
            v = v / 3 + 1.0f/3.0f;

        } else if (x > 0 && abs(x) >= abs(y) && abs(x) >= abs(z)) {
            u = -z / abs(x) / 2 + 0.5f;
            v = y / abs(x) / 2 + 0.5f;

            // these are coords within respective face, now we map to correct position of face within the texture

            u = u / 4 + 0.5f;
            v = v / 3 + 1.0f / 3.0f;

        } else if (y < 0 && abs(y) >= abs(x) && abs(y) >= abs(z)) {
            u = x / abs(y) / 2 + 0.5f;
            v = z / abs(y) / 2 + 0.5f;

            // these are coords within respective face, now we map to correct position of face within the texture

            u = u / 4 + 0.25f;
            v = v / 3;

        } else if (y > 0 && abs(y) >= abs(x) && abs(y) >= abs(z)) {
            u = x / abs(y) / 2 + 0.5f;
            v = -z / abs(y) / 2 + 0.5f;

            // these are coords within respective face, now we map to correct position of face within the texture

            u = u / 4 + 0.25f;
            v = v / 3 + 2.0f/3.0f;

        } else if (z < 0 && abs(z) >= abs(x) && abs(z) >= abs(y)) {
            u = -x / abs(z) / 2 + 0.5f;
            v = y / abs(z) / 2 + 0.5f;

            // these are coords within respective face, now we map to correct position of face within the texture

            u = u / 4 + 0.75f;
            v = v / 3 + 1.0f / 3.0f;
        } else { //z>0, abs(z) is max
            u = x / abs(z) / 2 + 0.5f;
            v = y / abs(z) / 2 + 0.5f;

            // these are coords within respective face, now we map to correct position of face within the texture

            u = u / 4 + 0.25f;
            v = v / 3 + 1.0f / 3.0f;
        } 
      

        glm::vec2 texCoord(u, v);

        return sampleTextureNearest(state.scene.environment, texCoord);
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

