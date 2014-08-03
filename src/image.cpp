// Scott Wiedemann
// 08/28/2009
// image.cpp
// Ray Tracing

#include "image.h"

#include <math.h>
#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

image::image ( int m, int n ) : xsize(m), ysize(n) {
	rgb = new RGB[m*n];
}

RGB &image::pixel ( int i, int j ) {
	return rgb[i+xsize*j];
}


static unsigned char clampnround ( double x ) {
	if (x>255)
		x = 255;
	if (x<0) 
		x = 0;
	return (unsigned char)floor(x+.5);
}

void image::save_to_ppm_file ( const char *filename ) {
	ofstream ofs(filename,ios::binary);
	assert(ofs);
	ofs << "P6" << endl;
	ofs << xsize << " " << ysize << endl << 255 << endl;
	for ( int i=0; i<xsize*ysize; i++ )	{
		unsigned char r = clampnround(256*rgb[i].r);
		unsigned char g = clampnround(256*rgb[i].g);
		unsigned char b = clampnround(256*rgb[i].b);
		ofs.write((char*)&r,sizeof(char));
		ofs.write((char*)&g,sizeof(char));
		ofs.write((char*)&b,sizeof(char));
	}
}
