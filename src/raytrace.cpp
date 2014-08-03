// Scott Wiedemann
// 08/28/2009
// raytrace.cpp
// Ray Tracing

#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <string>
#include <vector>

#include "image.h"

using namespace std;

class point
{
	private:
		double x,y,z;
	public:	
		double X() { return x; }
		double Y() { return y; }
		double Z() { return z; }
	
		point()
		{ }

		point(double A[3])
		{
			x=A[0];
			y=A[1];
			z=A[2];
		}

		point(const double& A, const double& B, const double& C)
		{
			x=A;
			y=B;
			z=C;
		}

		point operator+ (const point& a) const
		{
			point result;
			result.x = x+a.x;
			result.y = y+a.y;
			result.z = z+a.z;
			return result;
		}

		point operator+ (const double& a) const
		{
			point result;
			result.x = x+a;
			result.y = y+a;
			result.z = z+a;
			return result;
		}

		point operator- (const point& a) const
		{
			point result;
			result.x = x-a.x;
			result.y = y-a.y;
			result.z = z-a.z;
			return result;
		}

		point operator- (const double& a) const
		{
			point result;
			result.x = x-a;
			result.y = y-a;
			result.z = z-a;
			return result;
		}

		point operator* (const point& a) const
		{
			point result;
			result.x = x*a.x;
			result.y = y*a.y;
			result.z = z*a.z;
			return result;
		}

		point operator* (const double& a) const
		{
			point result;
			result.x = x*a;
			result.y = y*a;
			result.z = z*a;
			return result;
		}

		point operator/ (const point& a) const
		{
			point result;
			result.x = x/a.x;
			result.y = y/a.y;
			result.z = z/a.z;
			return result;
		}

		point operator/ (const double& a) const
		{
			point result;
			result.x = x/a;
			result.y = y/a;
			result.z = z/a;
			return result;
		}

		point normalize()
		{
			point result = point(x,y,z);
			result = result/(result.magnitude());
			return result;
		}

		double magnitude()
		{
			return sqrt((x*x)+(y*y)+(z*z));
		}

		double dot(const point& P)
		{
			return x*P.x+y*P.y+z*P.z;
		}

		point cross(const point& P)
		{
			return point((y*P.z-z*P.y),(z*P.x-x*P.z),(x*P.y-y*P.x));
		}
};

class material
{
	public:
		double k_diff_R, k_diff_G, k_diff_B;
		double k_ambient_R, k_ambient_G, k_ambient_B;
		double k_spec;
		double n_spec;

		material()
		{ }

		material(double K_Diff[3], double K_Ambient[3], double K_Spec, double N_Spec)
		{
			k_diff_R=K_Diff[0];
			k_diff_G=K_Diff[1];
			k_diff_B=K_Diff[2];
			k_ambient_R=K_Ambient[0];
			k_ambient_G=K_Ambient[1];
			k_ambient_B=K_Ambient[2];
			k_spec=K_Spec;
			n_spec=N_Spec;
		}
};

class triangle
{
	public:
		point a1, a2, a3;
		material m;

		triangle()
		{ }

		triangle(double A1[3], double A2[3], double A3[3], material M)
		{
			a1=point(A1);
			a2=point(A2);
			a3=point(A3);
			m=M;
		}
};

class sphere
{
	public:
		point center;
		double radius;
		material m;

		sphere()
		{ }

		sphere(double Center[3], double Radius, material M)
		{
			center=point(Center);
			radius=Radius;
			m=M;
		}
};

struct intersectionInfo
{
	char type;
	unsigned int index;
	double t;
};

// global variables
int resolution_x, resolution_y;
vector <triangle> Triangles;
vector <sphere> Spheres;
point view_point, eye_ray, light_source, lower_left_corner, horizontal_point, vertical_point;
double light_intensity, ambient_light_intensity;

point eyeRay(const int& X, const int& Y)
{
	return (lower_left_corner+horizontal_point*(((double)X+0.5)/resolution_x)+vertical_point*(((double)Y+0.5)/resolution_y))-view_point;
}

double sphereIntersection(point initialPoint, point Point, sphere Sphere)
{
	point co = initialPoint-Sphere.center;
	double t =-1.0;
	// t=-B+/-sqrt(B^2-4AC)/2A
	double A = (Point.magnitude())*(Point.magnitude());
	double B = 2*(Point.dot(co));
	double C = ((co.magnitude()*co.magnitude())-(Sphere.radius*Sphere.radius));
	double D = B*B-4*A*C;

	if(D>0)
	{
		double t1=(-B-sqrt(D))/(2*A);
		double t2=(-B+sqrt(D))/(2*A);

		if(t2<t1 && t2>0)
		{
			t1=t2;
		}
		// positive?
		if(t1>=0)
		{
			t=t1;
		}
	}
	return t;
}

double triangleIntersection(point initialPoint, point Point, triangle Triangle)
{
	double t =-1;
	// normal;
	point n = (Triangle.a2-Triangle.a1).cross(Triangle.a3-Triangle.a1);
	t=-1*(((initialPoint-Triangle.a1).dot(n))/(Point.dot(n)));
	if(t<0)
	{
		return -1;
	}
	
	point p=initialPoint + (Point*t);

	point v1 = (Triangle.a1-p).cross(Triangle.a2-p);
	point v2 = (Triangle.a2-p).cross(Triangle.a3-p);
	point v3 = (Triangle.a3-p).cross(Triangle.a1-p);

	if((v1.dot(v2)<0)||(v2.dot(v3)<0)||(v3.dot(v1)<0))
	{
		return -1;
	}
	return t;
}

intersectionInfo nearestIntersection(point E)
{
	intersectionInfo I;
	I.t=99999;

	// spheres
	for(unsigned int i=0; i<Spheres.size(); i++)
	{
		double intersection = sphereIntersection(view_point, E, Spheres[i]);
		if(intersection<I.t && intersection !=-1)
		{
			I.t=intersection;
			I.index=i;
			I.type='S';
		}
	}

	// triangles
	for(unsigned int i=0; i<Triangles.size(); i++)
	{
		double intersection = triangleIntersection(view_point, E, Triangles[i]);
		if(intersection<I.t && intersection !=-1)
		{
			I.t=intersection;
			I.index=i;
			I.type='T';
		}
	}

	if(I.t==99999)
	{
		I.t=-1;
	}
	return I;
}

point viewPointVector(const point& P)
{
	return (view_point-P).normalize();
}

point lightSourceVector(const point& P)
{
	return (light_source-P).normalize();
}

point normalTriangle(const triangle& T)
{
	point normal = (T.a3-T.a1).cross(T.a2-T.a1);
	if((normal).dot(view_point-T.a1)<0)
	{
		normal=normal*-1.0;
	}
	return normal;
}

point normalSphere(const sphere& S, point P)
{
	return P-S.center;
}

bool inShadow(intersectionInfo II) {
	point p = view_point + eye_ray*II.t;
	point shadow = lightSourceVector(p);

	for(unsigned int i=0; i<Triangles.size(); i++) {
		if(i!=II.index || II.type!='T') {
			double t = triangleIntersection(p, shadow, Triangles[i]);
			if(t>=0) {
				return true;
			}
		}
	}

	for(unsigned int i=0; i<Spheres.size(); i++) {
		if(i!=II.index || II.type!='S') {
			double t = sphereIntersection(p, shadow, Spheres[i]);
			if(t>=0) {
				return true;
			}
		}
	}

	return false;
}

RGB illumination(intersectionInfo II)
{
	// p = o + dt
	point p = view_point + eye_ray*II.t;

	// output colors	
	RGB colors;
	point normal;

	if(II.type=='T')
	{
		triangle T=Triangles[II.index];
		normal=normalTriangle(T);

		if( (((normal).dot(light_source-p))<0) || inShadow(II) )
		{
			// p in shadow of its primitive, return the value of the Ambient term
			colors.r=T.m.k_ambient_R*ambient_light_intensity;
			colors.g=T.m.k_ambient_G*ambient_light_intensity;
			colors.b=T.m.k_ambient_B*ambient_light_intensity;
		}
		else
		{
			point N = normal.normalize();
			point L = lightSourceVector(p);
			point V = viewPointVector(p);
			point H = (L+V).normalize();

			// some values not actually points, but using point class for methods!!!

			// Ambient term
			// Ia=I*ka
			point Ia = point(T.m.k_ambient_R, T.m.k_ambient_G, T.m.k_ambient_B)*ambient_light_intensity;

			// Diffuse term
			// Id = I*kd*(N dot L);
			point Id = point(T.m.k_diff_R, T.m.k_diff_G, T.m.k_diff_B)*(N).dot(L)*light_intensity;
			
			// Specular term
			// Is = I*ks*(H dot N)^n;
			double Is = pow(((H).dot(N)),(T.m.n_spec))*T.m.k_spec*light_intensity;

			// Iout = ambient term + diffuse term + specular term 
			// Iout = Ia + Id + Is
			point Iout = Ia + Id + Is;
			colors.r=Iout.X();
			colors.g=Iout.Y();
			colors.b=Iout.Z();
		}
	}
	else
	{
		sphere S=Spheres[II.index];
		normal=normalSphere(S, p);

		if( (((normal).dot(light_source-p))<0) || inShadow(II) )
		{
			// p in shadow of its primitive, return the value of the Ambient term
			colors.r=S.m.k_ambient_R*ambient_light_intensity;
			colors.g=S.m.k_ambient_G*ambient_light_intensity;
			colors.b=S.m.k_ambient_B*ambient_light_intensity;
		}
		else
		{
			point N = normal.normalize();
			point L = lightSourceVector(p);
			point V = viewPointVector(p);
			point H = (L+V).normalize();

			point Ia = point(S.m.k_ambient_R, S.m.k_ambient_G, S.m.k_ambient_B)*ambient_light_intensity;
			point Id = point(S.m.k_diff_R, S.m.k_diff_G, S.m.k_diff_B)*(N).dot(L)*light_intensity;
			double Is = pow(((H).dot(N)),(S.m.n_spec))*S.m.k_spec*light_intensity;

			point Iout = Ia + Id + Is;
			colors.r=Iout.X();
			colors.g=Iout.Y();
			colors.b=Iout.Z();
		}
	}
	return colors;
}

bool read_input_file(string Filename)
{
	ifstream ifs(Filename.c_str());
	if(ifs) {
		double temp_point[3];
		int number_of_primitives;

		ifs >> resolution_x >> resolution_y;

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		view_point=point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		lower_left_corner=point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		horizontal_point=point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		vertical_point=point(temp_point);

		ifs >> temp_point[0] >> temp_point[1] >> temp_point[2];
		light_source=point(temp_point);

		ifs >> light_intensity;
		ifs >> ambient_light_intensity;
		ifs >> number_of_primitives;

		material temp_material;
		double k_diffuse[3], k_ambient[3];
		double k_specular, n_specular;
		char primitive_type;

		for(int i=0; i<number_of_primitives; i++) {
			ifs >> primitive_type;
			double center[3];
			double radius;
			double a1[3], a2[3], a3[3];

			switch(toupper(primitive_type))
			{
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
					assert(0);
					break;
			}

			ifs >> k_diffuse[0] >> k_diffuse[1] >> k_diffuse[2];
			ifs >> k_ambient[0] >> k_ambient[1] >> k_ambient[2];
			ifs >> k_specular >> n_specular;
			temp_material=material(k_diffuse, k_ambient, k_specular, n_specular);

			switch(toupper(primitive_type))
			{
				case 'S':
					Spheres.push_back(sphere(center, radius, temp_material));
					break;
				case 'T':
					Triangles.push_back(triangle(a1, a2, a3, temp_material));
					break;
				default:
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
		outputFileName = inputFileName.substr(0,inputFileName.length()-4)+".ppm";
	}
	if(!read_input_file(inputFileName)) {
		cerr << "File " << inputFileName << " does not exist, or can't be read." << endl;
		return 2;
	}

	image img(resolution_x, resolution_y);
	for(int x=0; x<resolution_x; x++) {
		for(int y=0; y<resolution_y; y++) {
			eye_ray = eyeRay(x, y);
			intersectionInfo ii = nearestIntersection(eye_ray);
			
			RGB &pix = img.pixel(x, resolution_y-y-1);
			
			if(ii.t>-1) {
				pix = illumination(ii);
			}
		}
	}

	img.save_to_ppm_file(outputFileName.c_str());
	return 0;
}
