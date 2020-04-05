#include "geometry.h"

class image_texture {
    public:
	unsigned char *data;
        int nx, ny;
        image_texture() {}
        image_texture(unsigned char *pixels, int A, int B) : data(pixels), nx(A), ny(B) {}
	Vec3f value(double u, double v, const Vec3f& p) const {            
            if (data == nullptr)
                return Vec3f(0,1,1);

            auto i = static_cast<int>((  u)*nx);
            auto j = static_cast<int>((1-v)*ny-0.001);
            if (i < 0) i = 0;
            if (j < 0) j = 0;
            if (i > nx-1) i = nx-1;
            if (j > ny-1) j = ny-1;	    

            auto r = static_cast<int>(data[3*i + 3*nx*j+0]) / 255.0;
            auto g = static_cast<int>(data[3*i + 3*nx*j+1]) / 255.0;
            auto b = static_cast<int>(data[3*i + 3*nx*j+2]) / 255.0;

            return Vec3f(r, g, b);
        }
        ~image_texture() {
            delete data;
        }
};
