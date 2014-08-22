#include "Material.h"
#include "Point.h"
#include "Primitive.h"

class Sphere: public Primitive {
public:
	Point center;
	double radius;

	Sphere()
	{ }

	Sphere(double Center[3], double Radius, Material M) {
		center = Point(Center);
		radius = Radius;
		m = M;
	}

	virtual double Intersection(Point viewPoint, Point eyeRay) {
		double t = -1.0;
		Point co = viewPoint-center;

		// ray sphere intersection
		// sphere:
		// mag(x - c)^2 = r^2

		// line:
		// x = t*l + l0

		// substitute and expand

		// t=-B+/-sqrt(B^2-4AC)/2A
		double A = (eyeRay.magnitude())*(eyeRay.magnitude());
		double B = 2*(eyeRay.dot(co));
		double C = ((co.magnitude()*co.magnitude())-(radius*radius));
		double D = B*B-4*A*C;

		// ignore complex solutions
		if(D>0) {
			// get both roots
			double t1=(-B-sqrt(D))/(2*A);
			double t2=(-B+sqrt(D))/(2*A);

			// find closest, positive, interestion
			if(t2<t1 && t2>0) {
				t1=t2;
			}
			// positive?
			if(t1>=0) {
				t=t1;
			}
		}
		return t;
	}

	virtual Point Normal(Point viewPoint, Point intersectionPoint) {
		return intersectionPoint-center;
	}
};