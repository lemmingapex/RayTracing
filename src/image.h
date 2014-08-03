// Scott Wiedemann
// 08/28/2009
// image.h
// Ray Tracing

#include <math.h>
#include <fstream>
#include <iostream>
#include <assert.h>

typedef struct {
	double r,g,b;
}
RGB;

class image {
	public:
		image ( int m, int n );       // allocates image of specified size
		RGB &pixel ( int i, int j );  // access to a specific pixel
		void save_to_ppm_file ( const char *filename );
	private:
		int xsize,ysize; // resolution
		RGB *rgb;        // pixel intensities
};