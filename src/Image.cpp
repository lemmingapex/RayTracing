#include <math.h>
#include <fstream>
#include <iostream>

#include "Image.h"

using namespace std;

Image::Image(int m, int n): xsize(m), ysize(n) {
	rgb = new RGB[m*n];
}

RGB &Image::pixel(int i, int j) {
	return rgb[i+xsize*j];
}

static unsigned char clampnround(double x) {
	return (unsigned char)floor(max(min(x, 255.0), 0.0)+.5);
}

bool Image::save_to_ppm_file(const char *filename) {
	ofstream ofs(filename,ios::binary);
	if(ofs) {
		ofs << "P6" << endl;
		ofs << xsize << " " << ysize << endl << 255 << endl;
		for(int i=0; i<xsize*ysize; i++)	{
			unsigned char r = clampnround(256*rgb[i].r);
			unsigned char g = clampnround(256*rgb[i].g);
			unsigned char b = clampnround(256*rgb[i].b);
			ofs.write((char*)&r, sizeof(char));
			ofs.write((char*)&g, sizeof(char));
			ofs.write((char*)&b, sizeof(char));
		}
		return true;
	}
	return false;
}
