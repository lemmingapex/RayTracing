typedef struct {
	double r,g,b;
}
RGB;

class Image {
	public:
		// allocates image of specified size
		Image(int m, int n);

		// access to a specific pixel
		RGB &pixel(int i, int j);

		bool save_to_ppm_file(const char *filename);
	private:
		int xsize, ysize;

		// pixel intensities
		RGB *rgb;
};
