// Scott Wiedemann
// 08/28/2009
// RayTrace.cpp
// Ray Tracing

#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <vecmath.h>

#include "Image.h"
#include "Material.h"
#include "Primitive.h"
#include "Triangle.h"
#include "Sphere.h"
#include "Ellipsoid.h"

using namespace std;

struct IntersectionInformation {
	unsigned int index;
	double t;
};

// global variables
vector <Primitive*> Primitives;

// scene stuff
int resolution_x, resolution_y;
Vector3f view_point, eye_ray, light_source, lower_left_corner, horizontal_point, vertical_point;
double light_intensity, ambient_light_intensity;

Vector3f eyeRay(const int& X, const int& Y) {
	return (lower_left_corner+horizontal_point*(((double)X+0.5)/resolution_x)+vertical_point*(((double)Y+0.5)/resolution_y))-view_point;
}

Vector3f viewPointVector(const Vector3f& P) {
	return (view_point-P).normalized();
}

Vector3f lightSourceVector(const Vector3f& P) {
	return (light_source-P).normalized();
}

IntersectionInformation nearestIntersection(Vector3f E) {
	IntersectionInformation intersectionInformation;
	intersectionInformation.t = numeric_limits<double>::max();

	for(unsigned int i = 0; i < Primitives.size(); i++) {
		double intersection = Primitives[i]->Intersection(view_point, E);
		if(intersection < intersectionInformation.t && intersection != -1) {
			intersectionInformation.t = intersection;
			intersectionInformation.index = i;
		}
	}

	if(intersectionInformation.t == numeric_limits<double>::max()) {
		intersectionInformation.t = -1.0;
	}
	return intersectionInformation;
}

bool isInShadow(IntersectionInformation intersectionInformation) {
	Vector3f p = view_point + eye_ray*intersectionInformation.t;
	Vector3f shadow = lightSourceVector(p);

	for(unsigned int i=0; i<Primitives.size(); i++) {
		if(i != intersectionInformation.index) {
			double t = Primitives[i]->Intersection(p, shadow);
			if(t > 0) {
				return true;
			}
		}
	}
	return false;
}

/**
* Uses Phong shading
*/
RGB illumination(IntersectionInformation intersectionInformation) {
	// p = o + dt
	Vector3f p = view_point + eye_ray*intersectionInformation.t;

	// output colors
	RGB colors;

	Primitive* P = Primitives[intersectionInformation.index];
	Vector3f N = P->Normal(view_point, p).normalized();
	Vector3f L = lightSourceVector(p);
	Vector3f Iout;

	// check if the intersection is in shadow
	// in shadow can occur in 2 cases
	// 1) the primitive can block itself
	// 2) a different primitive can block the intersection
	if((N).dot(L) < 0 || isInShadow(intersectionInformation)) {
		// return the value of the Ambient term
		Iout = P->m.kAmbient*ambient_light_intensity;
	} else {
		Vector3f V = viewPointVector(p);
		Vector3f H = (L+V).normalized();

		// some values not actually points, but using point class for methods!

		// Ambient term
		// Ia=I*ka
		Vector3f Ia = P->m.kAmbient*ambient_light_intensity;

		// Diffuse term
		// Id = I*kd*(N dot L);
		Vector3f Id = P->m.kDiff*(N).dot(L)*light_intensity;
		
		// Specular term
		// Is = I*ks*(H dot N)^n;
		double Is = pow(((H).dot(N)),(P->m.nSpec))*P->m.kSpec*light_intensity;

		// Iout = ambient term + diffuse term + specular term 
		// Iout = Ia + Id + Is
		Iout = Ia + Id + Is;
	}
	colors.r = Iout.x();
	colors.g = Iout.y();
	colors.b = Iout.z();
	return colors;
}

/**
 * Populates Primitives and the View
*/
bool read_input_file(string Filename) {
	ifstream ifs(Filename.c_str());
	if(ifs) {
		float temp_point[3];
		int number_of_primitives;

		ifs >> resolution_x >> resolution_y;

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		view_point=Vector3f(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		lower_left_corner=Vector3f(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		horizontal_point=Vector3f(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		vertical_point=Vector3f(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		light_source=Vector3f(temp_point);

		ifs >> light_intensity;
		ifs >> ambient_light_intensity;
		ifs >> number_of_primitives;

		Material temp_material;
		float k_diffuse[3], k_ambient[3];
		double k_specular, n_specular;
		char primitive_type;

		for(int i=0; i<number_of_primitives; i++) {
			ifs >> primitive_type;
			primitive_type = toupper(primitive_type);
			float center[3];
			float axisLengths[3];
			double radius;
			float a1[3], a2[3], a3[3];

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
				case 'E':
					ifs >> center[0] >> center[1] >> center[2];
					ifs >> axisLengths[0] >> axisLengths[1] >> axisLengths[2];
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
				case 'E':
					Primitives.push_back(new Ellipsoid(center, axisLengths, temp_material));
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
			IntersectionInformation intersectionInformation = nearestIntersection(eye_ray);
			
			if(intersectionInformation.t > 0) {
				RGB &pix = img.pixel(x, resolution_y-y-1);
				pix = illumination(intersectionInformation);
			}
		}
	}

	if(!img.save_to_ppm_file(outputFileName.c_str())) {
		cerr << "Problem writing " << outputFileName << "." << endl;
		return 3;
	}
	return 0;
}
