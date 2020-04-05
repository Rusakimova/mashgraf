#include <iostream>
#include <cstdint>
#include <limits>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include "geometry.h"
#include "Bitmap.h"
#include "texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace std;
const int width = 1024;
const int height = 1024;

struct Material {
    Material(const float &r, const Vec4f &a, const Vec3f &color) : ref(r), ambient(a), diffuse_color(color) {}
    Material() : ref(), ambient(1, 0, 0, 0), diffuse_color() {}
    float ref;
    Vec4f ambient;
    Vec3f diffuse_color;    
};

struct ray {
    Vec3f orig;
    Vec3f dir;
    ray(const Vec3f &o, const Vec3f &d):orig(o), dir(d) {}
};

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

class Sphere {
public:
    float radius;
    Vec3f center;
    Material material;
    int specular;
    Sphere (const Vec3f &c, const float &r, const Material &m, const int &s): radius(r), center(c), material(m), specular(s) {}
    bool intersect(const ray &r, float &t0) const {
        Vec3f L = center - r.orig;
        float tca = L*r.dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0 = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

class Triangle {
public:
    Vec3f V0, V1, V2;
    Material material;
    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Material &m): V0(v0), V1(v1), V2(v2), material(m){}   
    bool intersect(const ray &r, float &t) const {
	float eps = 0.0001;
	Vec3f A = V1-V0;
	Vec3f B = V2-V1;
	Vec3f C = V0-V2;
	Vec3f N = cross(A, B);
	float s = N*r.dir;
	if (fabs(s) < eps) return false;
	t = (N*V0 + N*r.orig)/s;
	if(t < 0) return false;
	Vec3f point = r.orig + r.dir*t;
	Vec3f D; 
	Vec3f vp0 = point - V0; 
	D = cross(A, vp0); 
	if (N*D < 0) return false;  
	Vec3f vp1 = point - V1; 
	D = cross(B, vp1); 
	if (N*D < 0)  return false;  
	Vec3f vp2 = point - V2; 
	D = cross(C, vp2); 
	if (N*D < 0) return false;  
	return true;
    }    
};    
     
Vec3f refract(const Vec3f &dir, const Vec3f &N, const float &ref) { 
    float cosi = - max(-1.f, std::min(1.f, dir*N));
    float etai = 1, etat = ref;
    Vec3f n = N;
    if (cosi < 0) { 
	cosi = -cosi;
        swap(etai, etat); n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k < 0 ? Vec3f(0,0,0) : dir*eta + n*(eta * cosi - sqrtf(k));
}

void get_sphere_uv(const Vec3f &p, float &u, float &v) {
    float phi = atan2(p.z, p.x);
    float theta = asin(p.y);
    u = 1-(phi+M_PI)/(2 * M_PI);
    v = (theta + M_PI/2)/M_PI;
}

bool scene_intersect(const ray &r, const vector<Sphere> &spheres, const vector <Triangle> &trians, int &specular, Vec3f &hit, Vec3f &N, Material &material, int &flag) {
    float dist = numeric_limits<float>::max();
    for (size_t i = 0; i < trians.size(); i++) {
	float dist_i;
	if (trians[i].intersect(r, dist_i) && dist_i < dist) {
	    flag = 0;
	    N = Vec3f(0, 0, 0);
	    dist = dist_i;
	    hit = r.orig + r.dir*dist;
	    material = trians[i].material;	
	} 
    }   
    for (size_t i=0; i < spheres.size(); i++) {	
        float dist_i;
        if (spheres[i].intersect(r, dist_i) && dist_i < dist) {
            dist = dist_i;
            hit = r.orig + r.dir*dist_i;
            N = (hit - spheres[i].center).normalize();
	    if (i == 0) {
		flag = 1;
	    } else {
		material = spheres[i].material;
		flag = 0;
	    }	        
        }
    }    
    float checkerboard_dist = numeric_limits<float>::max();
    if (fabs(r.dir.y)>1e-3)  {	
	float d = -(r.orig.y+4)/r.dir.y; 
        Vec3f pt = r.orig + r.dir*d;
        if (d>0 && d<dist) {
	    flag = 0;
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
	    specular = 125;
            material.diffuse_color = (int(.3*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(0.3,0.55,0.65) : Vec3f(0.4, 0.7, 0.95);
	    material.ambient = Vec4f(0.6,  0.3, 0.1, 0.0);
	    material.ref = 1.0;	    
        }
    }
    return min(dist, checkerboard_dist) < 1000;
}


Vec3f calc_ray(const ray &r, const vector<Sphere> &spheres, const vector<Light> &lights, const vector <Triangle> &trians, const image_texture &image_map, const image_texture &image_sphere, int depth = 0) {
    Vec3f point, N;
    float t0, u, v;
    int specular = 0, flag;
    Material material;
    if (depth>3 || !scene_intersect(r, spheres, trians, specular, point, N, material, flag)) {
	Material blue(1.0, Vec4f(0.6,  0.3, 0.0, 0.0), Vec3f(0, 0, 0.8));
	Sphere sphere(Vec3f(0,    0,   0), 100, blue, 0);
	sphere.intersect(r, t0);		
	get_sphere_uv((r.dir*t0-Vec3f(0, 0, 0))/100, u, v);
	return image_map.value(u, v, point);	
    }
    if (flag == 1) {	
	get_sphere_uv((point-spheres[0].center)/spheres[0].radius, u, v);
	material.diffuse_color = image_sphere.value(u, v, point);
    }
    float light_intens = 0, specular_intens = 0;    
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();
        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; 
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(ray(shadow_orig, light_dir), spheres, trians, specular, shadow_pt, shadow_N, tmpmaterial, flag) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;
        light_intens  += lights[i].intensity * max(0.f, light_dir*N);
	int s;  
	if (specular == 0) {      
	    s = spheres[i].specular;
	} else {
	    s = specular;
	}
	if (s != -1) {
	    Vec3f R = light_dir - N*2.f*(light_dir*N);
	    specular_intens += powf(max(0.f, R*r.dir), s)*lights[i].intensity;	    	    
	}	
    } 
    if (N == Vec3f(0, 0, 0)) return material.diffuse_color; 
    Vec3f reflect_dir = r.dir - N*2.f*(r.dir*N);
    Vec3f refract_dir = refract(r.dir, N, material.ref).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; 
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;    
    Vec3f reflect_color = calc_ray(ray(reflect_orig, reflect_dir), spheres, lights, trians, image_map, image_sphere, depth + 1);
    Vec3f refract_color = calc_ray(ray(refract_orig, refract_dir), spheres, lights, trians, image_map, image_sphere, depth + 1);      
    return material.diffuse_color * light_intens * material.ambient[0] + Vec3f(1., 1., 1.)*specular_intens * material.ambient[1] + reflect_color*material.ambient[2] + refract_color*material.ambient[3];
}

void render(const vector<Sphere> &spheres, const vector<Light> &lights, const vector<Triangle> &trians, const image_texture &image_map, const image_texture &image_sphere, const string file) {
    const int width = 512;
    const int height = 512;
    const float fov = M_PI/2.;    
    vector<Vec3f> framebuffer(width*height);
    for (int j = 0; j<height; j++) {
        for (int i = 0; i<width; i++) {
	    Vec3f res(0, 0, 0);
	    float x1, y1, x2, y2;
	    x1 =  (2*(i+0.3)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
	    y1 = -(2*(j+0.3)/(float)height - 1)*tan(fov/2.);
	    Vec3f dir = Vec3f(x1, y1, -1).normalize();
	    res = res + calc_ray(ray(Vec3f(0,0,0), dir), spheres, lights, trians, image_map, image_sphere);
	    y2 = -(2*(j+0.8)/(float)height - 1)*tan(fov/2.);
	    dir = Vec3f(x1, y2, -1).normalize();
	    res = res + calc_ray(ray(Vec3f(0,0,0), dir), spheres, lights, trians, image_map, image_sphere);
	    x2 = (2*(i+0.8)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
	    dir = Vec3f(x2, y2, -1).normalize();
	    res = res + calc_ray(ray(Vec3f(0,0,0), dir), spheres, lights, trians, image_map, image_sphere);
	    dir = Vec3f(x2, y1, -1).normalize();
	    res = res + calc_ray(ray(Vec3f(0,0,0), dir), spheres, lights, trians, image_map, image_sphere);
	    framebuffer[i+j*width] = res/4;
	}
    }
    SaveBMP(file.c_str(), framebuffer, width, height);
}


int main(int argc, const char** argv) {
    const int width = 1024;
    const int height = 1024;
    string scene;
    string outFilePath;
    for (int i = 0; i<argc; i++) {
	string option(argv[i]);
	if (option == "-scene") {
	    scene = argv[i+1];
	} else if (option == "-out") {
	    outFilePath = argv[i + 1];
	}
    }
    if (scene == "1") {
	int w1, h1, n, w2, h2;	
	unsigned char *texture_map = stbi_load("envmap.jpg", &w1, &h1, &n, 0);
	image_texture image_map(texture_map, w1, h1);
	unsigned char *texture_sphere = stbi_load("asphalt.JPG", &w2, &h2, &n, 0);
	image_texture image_sphere(texture_sphere, w2, h2);	
	Material orange(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(1, 0.3, 0));
	Material green(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.1, 0.4, 0.3));
	Material mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0));
	Material glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8));
	vector<Sphere> spheres;
	vector<Light> lights;
	vector<Triangle> triangles;
	lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
	lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));
	Material blue(1.0, Vec4f(0.6,  0.3, 0.0, 0.0), Vec3f(0, 0, 0.8));
	spheres.push_back(Sphere(Vec3f(-5,    -4,   -10), 1.5, blue, 100));
	spheres.push_back(Sphere(Vec3f( 3,    -1.5,   -7), 2, mirror, 10));
	spheres.push_back(Sphere(Vec3f( -1.0,    -3,   -8), 0.5, orange, 50));
	spheres.push_back(Sphere(Vec3f( -2.0,    -3,   -8), 0.5, orange, 50));
	spheres.push_back(Sphere(Vec3f( 0,    -3,   -7.5), 0.5, orange, 50));
	spheres.push_back(Sphere(Vec3f( -0.5,    -3,   -7), 0.5, orange, 50));
	spheres.push_back(Sphere(Vec3f( -1.5,    -3,   -7), 0.5, orange, 50));
	spheres.push_back(Sphere(Vec3f( -1,    -3,   -6), 0.5, orange, 50));
	spheres.push_back(Sphere(Vec3f(6, 4, -14), 1, glass, 125));	
	spheres.push_back(Sphere(Vec3f(-5, 2, -14), 1, glass, 125));
	triangles.push_back(Triangle(Vec3f(-3, -4, -10), Vec3f(-2.5, 0, -10), Vec3f(-2, -4, -10), green));
	triangles.push_back(Triangle(Vec3f(-2, -4, -10), Vec3f(-1.5, 0, -10), Vec3f(-0.5, -4, -10), green));
	triangles.push_back(Triangle(Vec3f(-0.5, -4, -10), Vec3f(0, 0, -10), Vec3f(1, -4, -9), green));
	render(spheres, lights, triangles, image_map, image_sphere, outFilePath);	
    } else {
	printf("there is only 1 scene\n");
    } 
    return 0;
}
