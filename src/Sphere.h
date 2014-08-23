#include <vecmath.h>

#include "Material.h"
#include "Primitive.h"

class Sphere: public Primitive {
public:
	Vector3f center;
	double radius;

	Sphere()
	{ }

	Sphere(float Center[3], double Radius, Material M) {
		center = Vector3f(Center);
		radius = Radius;
		m = M;
	}

	virtual double Intersection(Vector3f viewPoint, Vector3f l) {
		double t = -1.0;
		Vector3f diff = viewPoint-center;

		// ray sphere intersection
		// sphere:
		// mag(x - c)^2 = r^2

		// line:
		// x = t*l + l0

		// substitute and expand

		// A*t^2 + B*t + C = 0
		// t = -B+/-sqrt(B^2-4AC)/2A
		double A = (l.abs())*(l.abs());
		double B = 2*(l.dot(diff));
		double C = ((diff.abs()*diff.abs())-(radius*radius));
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

	virtual Vector3f Normal(Vector3f viewPoint, Vector3f intersectionPoint) {
		Vector3f normal = intersectionPoint-center;

		if((intersectionPoint-center).abs() > (viewPoint-center).abs()) {
			normal = -1.0*normal;
		}

		return normal;
	}
};