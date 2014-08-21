// Scott Wiedemann
// 08/28/2009
// raytrace.cpp
// Ray Tracing

#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Image.h"
#include "Point.h"
#include "Material.h"
#include "Primitive.h"
#include "Sphere.h"
#include "Triangle.h"

using namespace std;

struct intersectionInfo {
	unsigned int index;
	double t;
};

// global variables
int resolution_x, resolution_y;
vector <Primitive*> Primitives;

// scene stuff
Point view_point, eye_ray, light_source, lower_left_corner, horizontal_point, vertical_point;
double light_intensity, ambient_light_intensity;

Point eyeRay(const int& X, const int& Y) {
	return (lower_left_corner+horizontal_point*(((double)X+0.5)/resolution_x)+vertical_point*(((double)Y+0.5)/resolution_y))-view_point;
}

Point viewPointVector(const Point& P) {
	return (view_point-P).normalize();
}

Point lightSourceVector(const Point& P) {
	return (light_source-P).normalize();
}

intersectionInfo nearestIntersection(Point E) {
	intersectionInfo I;
	I.t=99999;

	for(unsigned int i=0; i<Primitives.size(); i++) {
		double intersection = Primitives[i]->Intersection(view_point, E);
		if(intersection<I.t && intersection != -1) {
			I.t=intersection;
			I.index=i;
		}
	}

	if(I.t==99999) {
		I.t=-1;
	}
	return I;
}

bool inShadow(intersectionInfo II) {
	Point p = view_point + eye_ray*II.t;
	Point shadow = lightSourceVector(p);

	for(unsigned int i=0; i<Primitives.size(); i++) {
		if(i!=II.index) {
			double t = Primitives[i]->Intersection(p, shadow);
			if(t>=0) {
				return true;
			}
		}
	}

	return false;
}

/**
* Uses Phong shading
*/
RGB illumination(intersectionInfo II) {
	// p = o + dt
	Point p = view_point + eye_ray*II.t;

	// output colors
	RGB colors;

	Primitive* P = Primitives[II.index];
	Point normal=P->Normal(view_point, p);

	if( (((normal).dot(light_source-p))<0) || inShadow(II) ) {
		// p in shadow of its primitive, return the value of the Ambient term
		colors.r=P->m.k_ambient_R*ambient_light_intensity;
		colors.g=P->m.k_ambient_G*ambient_light_intensity;
		colors.b=P->m.k_ambient_B*ambient_light_intensity;
	} else {
		Point N = normal.normalize();
		Point L = lightSourceVector(p);
		Point V = viewPointVector(p);
		Point H = (L+V).normalize();

		// some values not actually points, but using point class for methods!!!

		// Ambient term
		// Ia=I*ka
		Point Ia = Point(P->m.k_ambient_R, P->m.k_ambient_G, P->m.k_ambient_B)*ambient_light_intensity;

		// Diffuse term
		// Id = I*kd*(N dot L);
		Point Id = Point(P->m.k_diff_R, P->m.k_diff_G, P->m.k_diff_B)*(N).dot(L)*light_intensity;
		
		// Specular term
		// Is = I*ks*(H dot N)^n;
		double Is = pow(((H).dot(N)),(P->m.n_spec))*P->m.k_spec*light_intensity;

		// Iout = ambient term + diffuse term + specular term 
		// Iout = Ia + Id + Is
		Point Iout = Ia + Id + Is;
		colors.r=Iout.X();
		colors.g=Iout.Y();
		colors.b=Iout.Z();
	}
	return colors;
}

/**
 * Populates Primitives and the View
*/
bool read_input_file(string Filename) {
	ifstream ifs(Filename.c_str());
	if(ifs) {
		double temp_point[3];
		int number_of_primitives;

		ifs >> resolution_x >> resolution_y;

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		view_point=Point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		lower_left_corner=Point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		horizontal_point=Point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		vertical_point=Point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		light_source=Point(temp_point);

		ifs >> light_intensity;
		ifs >> ambient_light_intensity;
		ifs >> number_of_primitives;

		Material temp_material;
		double k_diffuse[3], k_ambient[3];
		double k_specular, n_specular;
		char primitive_type;

		for(int i=0; i<number_of_primitives; i++) {
			ifs >> primitive_type;
			primitive_type = toupper(primitive_type);
			double center[3];
			double radius;
			double a1[3], a2[3], a3[3];

			switch(primitive_type) {
				case 'S':
					ifs >> center[0] >> center[1] >> center[2];
					ifs >> radius;
					break;
				case 'T':
					ifs >> a1[0] >> a1[1] >> a1[2];
					ifs >> a2[0] >> a2[1] >> a2[2];
					ifs >> a3[0] >> a3[1] >> a3[2];
					break;
				default:
					cerr << "Unrecognized primitive: " << toupper(primitive_type) << endl;
					continue;
					break;
			}

			ifs >> k_diffuse[0] >> k_diffuse[1] >> k_diffuse[2];
			ifs >> k_ambient[0] >> k_ambient[1] >> k_ambient[2];
			ifs >> k_specular >> n_specular;
			temp_material=Material(k_diffuse, k_ambient, k_specular, n_specular);

			switch(primitive_type) {
				case 'S':
					Primitives.push_back(new Sphere(center, radius, temp_material));
					break;
				case 'T':
					Primitives.push_back(new Triangle(a1, a2, a3, temp_material));
					break;
				default:
					continue;
					break;
			}
		}
		return true;
	}
	return false;
}

int main(int argc, const char *argv[]) {
	string inputFileName = "input.txt";
	string outputFileName = "output.ppm";

	if (argc > 2) {
		cerr << "Usage: " << argv[0] << " <" << inputFileName << ">" << endl;
		return 1;
	} else if (argc == 2) {
		inputFileName = argv[1];
		size_t dot = inputFileName.find_last_of(".");
		if (dot == string::npos) {
			dot = inputFileName.length();
		}

		outputFileName = inputFileName.substr(0, dot)+".ppm";
	}
	if(!read_input_file(inputFileName)) {
		cerr << "File " << inputFileName << " does not exist, or can not be read." << endl;
		return 2;
	}

	Image img(resolution_x, resolution_y);
	for(int x=0; x<resolution_x; x++) {
		for(int y=0; y<resolution_y; y++) {
			eye_ray = eyeRay(x, y);
			intersectionInfo ii = nearestIntersection(eye_ray);
			
			if(ii.t>-1) {
				RGB &pix = img.pixel(x, resolution_y-y-1);
				pix = illumination(ii);
			}
		}
	}

	if(!img.save_to_ppm_file(outputFileName.c_str())) {
		cerr << "Problem writing " << outputFileName << "." << endl;
		return 3;
	}
	return 0;
}
