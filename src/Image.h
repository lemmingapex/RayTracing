typedef struct {
	double r,g,b;
}
RGB;

class Image {
	public:
		Image ( int m, int n );       // allocates image of specified size
		RGB &pixel ( int i, int j );  // access to a specific pixel
		bool save_to_ppm_file ( const char *filename );
	private:
		int xsize,ysize; // resolution
		RGB *rgb;        // pixel intensities
};
