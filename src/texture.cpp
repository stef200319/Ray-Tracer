#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    float i = (1 - texCoord.y) * image.height;
    float j = texCoord.x * image.width;

    i = floor(i);
    j = floor(j);

    int pos = i * image.width + j;

    return image.pixels[pos];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    float i = (1 - texCoord.y) * image.height;
    float j = texCoord.x * image.width;

    int a = floor(i);
    int b = floor(j);
    int c;
    int d;

    float k = i - 0.5;
    float l = j - 0.5;

    if (floor(k) == floor(i))
        c = floor(i) + 1;
    else
        c = floor(i) - 1;

    if (floor(l) == floor(j))
        d = floor(j) + 1;
    else
        d = floor(j) - 1;


    if ((a == 0 || a == image.height-1) && (b == 0 || b == image.height-1)) {
        return image.pixels[a * image.width + b];
    }

    if (a == 0 || a == image.height-1) {
        int pos1 = a * image.width + b;
        int pos2 = a * image.width + d;

        float w1 = glm::abs(d + 0.5 - j);
        float w2 = glm::abs(b + 0.5 - j);

        return image.pixels[pos1] * w1 + image.pixels[pos2] * w2;
    }
        

    if (b == 0 || b == image.height-1) {
        int pos1 = a * image.width + b;
        int pos2 = c * image.width + b;

        float w1 = glm::abs(c + 0.5 - i);
        float w2 = glm::abs(a + 0.5 - i);

        return image.pixels[pos1] * w1 + image.pixels[pos2] * w2;
    }

    int pos1 = a * image.width + b;
    int pos2 = c * image.width + b;
    int pos3 = a * image.width + d;
    int pos4 = c * image.width + d;

    //pos2 = glm::clamp(1, 0, image.width * image.height);
    //pos3 = glm::clamp(1, 0, image.width * image.height);
    //pos4 = glm::clamp(1, 0, image.width * image.height);

    float w1 = glm::abs(c + 0.5 - i) * glm::abs(d + 0.5 - j);
    float w2 = glm::abs(a + 0.5 - i) * glm::abs(d + 0.5 - j);
    float w3 = glm::abs(c + 0.5 - i) * glm::abs(b + 0.5 - j);
    float w4 = glm::abs(a + 0.5 - i) * glm::abs(b + 0.5 - j);

    return image.pixels[pos1] * w1 + image.pixels[pos2] * w2 + image.pixels[pos3] * w3 + image.pixels[pos4] * w4;
}