#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <iostream>

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

//this is so weird, it becomes 10000 times faster by not using this boxFilter method
// i suppose cache misses, idk
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

void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }
    
    Screen mask(image.resolution(), true);
    mask.clear(glm::vec3 { 0.f, 0.f, 0.f });
    glm::vec3 tresholdColor { 0.9f, 0.9f, 0.9f };

    auto pixels = image.pixels();
    //initialize mask
    for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {
            auto currPixel = pixels[image.indexAt(x, y)];
            //if ( currPixel.x >= tresholdColor.x && currPixel.y >= tresholdColor.y  && currPixel.z >= tresholdColor.z)
            if (currPixel.x + currPixel.y + currPixel.z >= 0.6f)
                mask.setPixel(x, y, currPixel);
        }
    }
    int filterSize = 10;
    Screen maskHorizontal(image.resolution(), true);
    maskHorizontal.clear(glm::vec3 { 0.f, 0.f, 0.f });
    //box filter and scale on mask
   /* for (int x = 0; x < image.resolution().x; x++) {
        for (int y = 0; y < image.resolution().y; y++) {*/
    for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {
            //auto currPixel = pixels[mask.indexAt(x, y)];
            //if (currPixel.x >= tresholdColor.x && currPixel.y >= tresholdColor.y && currPixel.z >= tresholdColor.z)
                //mask2.setPixel(x, y, boxFilter(mask, x, y));
            glm::vec3 sum { 0.f, 0.f, 0.f };
            for (int i = -filterSize; i < filterSize + 1; ++i) {
                int index = image.indexAt(x + i, y);
                if (!(index < 0 || index >= pixels.size()))
                    sum += pixels[index];
            }

            // sum *= 1.f / ((filterSize * 2 + 1) * (filterSize * 2 + 1));
            sum /= (filterSize * 2 + 1);
            mask2.setPixel(x, y, sum);
        }
    }

    //add to original
    auto maskP = mask2.pixels();
    auto testMask = mask.pixels();
     for (int y = 0; y < image.resolution().y; y++) {
        for (int x = 0; x < image.resolution().x; x++) {
            auto currPixel = pixels[image.indexAt(x, y)];
            auto maskPixel =maskP[image.indexAt(x, y)];
                /* if (currPixel.x >= tresholdColor.x && currPixel.y >= tresholdColor.y && currPixel.z >= tresholdColor.z)
                */
            image.setPixel(x, y, testMask[image.indexAt(x, y)]);
            //image.setPixel(x, y, currPixel + maskPixel);
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
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    return 0; // This is clearly not the solution
}