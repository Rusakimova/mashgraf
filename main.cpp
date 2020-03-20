#include <iostream>
#include <cstdint>
#include <limits>
#include <cmath>
#include <fstream>
#include "geometry.h"
#include <string>
#include <vector>
#include <unordered_map>

const uint32_t RED   = 0x000000FF;
const uint32_t GREEN = 0x0000FF00;
const uint32_t BLUE  = 0x00FF0000;

class sphere {
public:
    float radius;
    Vec3f center;
    sphere (const Vec3f &c, const float &r): radius(r), center(c){}
    bool intersect(Vec3f center, Vec3f point, Vec3f dir, float r) {
        float dist;
        Vec3f cp(center[0]- point[0], center[1]- point[1], center[2]- point[2]);
        Vec3f dotvec = cross(dir, cp);
        dist = dotvec.norm()/dir.norm();
        if (dist > r) {
            return false;
        } else {
            return true;
        }
    }
};

void render(void) {
    const int width = 512;
    const int height = 512;
    const float fov = M_PI/2.;
    std::vector<Vec3f> framebuffer(width*height);
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            framebuffer[i+j*width] = Vec3f(j/float(height),i/float(width), 0);
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, sphere);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main(int argc, const char** argv)
{
    Sphere sphere(Vec3f(-3, 0, -16), 2);
    render(sphere);
    return 0;
  }
