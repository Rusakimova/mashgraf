#include <vector>
#include <fstream>
#include <cstring>
#include "geometry.h"
using namespace std;
struct Pixel { unsigned char r, g, b; };

void WriteBMP(const char* fname, Pixel *framebuffer, int width, int height)
{
  int paddedsize = (width*height) * sizeof(Pixel);

  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};

  bmpfileheader[ 2] = (unsigned char)(paddedsize    );
  bmpfileheader[ 3] = (unsigned char)(paddedsize>> 8);
  bmpfileheader[ 4] = (unsigned char)(paddedsize>>16);
  bmpfileheader[ 5] = (unsigned char)(paddedsize>>24);

  bmpinfoheader[ 4] = (unsigned char)(width    );
  bmpinfoheader[ 5] = (unsigned char)(width>> 8);
  bmpinfoheader[ 6] = (unsigned char)(width>>16);
  bmpinfoheader[ 7] = (unsigned char)(width>>24);
  bmpinfoheader[ 8] = (unsigned char)(height    );
  bmpinfoheader[ 9] = (unsigned char)(height>> 8);
  bmpinfoheader[10] = (unsigned char)(height>>16);
  bmpinfoheader[11] = (unsigned char)(height>>24);

  std::ofstream out(fname, std::ios::out | std::ios::binary);
  out.write((const char*)bmpfileheader, 14);
  out.write((const char*)bmpinfoheader, 40);
  out.write((const char*)framebuffer, paddedsize);
  out.flush();
  out.close();
}

void SaveBMP(const char* fname, vector<Vec3f> framebuffer, int w, int h)
{
	
  vector<Pixel> pixels(w*h);
  int left = w * (h-1);
  int right = w*h-1;
  int i = 0, k = 0;
  Pixel px;
  while (left >= 0) {      
    px.b       = (unsigned char)(255 * max(0.f, min(1.f, framebuffer[left+i].x)));
    px.g       = (unsigned char)(255 * max(0.f, min(1.f, framebuffer[left+i].y)));
    px.r       = (unsigned char)(255 * max(0.f, min(1.f, framebuffer[left+i].z)));
    pixels[k] = px;
    if((i+1)%w == 0) {
	left -= w;
	i = -1;
    }
    i++;
    k++;
  }
  WriteBMP(fname, &pixels[0], w, h);  
}
