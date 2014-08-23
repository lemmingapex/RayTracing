#include <vecmath.h>

#include "Material.h"
#include "Primitive.h"

class Ellipsoid: public Sphere {
public:
	Matrix3f A;

	Ellipsoid()
	{ }

	Ellipsoid(float Center[3], double Radius, float AxisLengths[3], Material M): Sphere(Center, Radius, M) {
		A = Matrix3f(1.0f/(AxisLengths[0]*AxisLengths[0]), 0.0f, 0.0f, 0.0f, 1.0f/(AxisLengths[1]*AxisLengths[1]), 0.0f, 0.0f, 0.0f, 1.0f/(AxisLengths[2]*AxisLengths[2]));
	}

	virtual double Intersection(Vector3f viewPoint, Vector3f l) {
		double t = -1.0;
		Vector3f diff = viewPoint - center;

		// ray ellipsoid intersection
		// ellipsoid:
		// (x - c)^T * A * (x - c) = 1

		// line:
		// x = t*l + l0

		// substitute and expand

		// A*t^2 + B*t + C = 0  
		// t = -B+/-sqrt(B^2-4AC)/2A
		Vector3f matDir = A*l;
		Vector3f matDiff = A*diff;
		double A = l.dot(matDir);
		double B = 2*(l.dot(matDiff));
		double C = (diff.dot(matDiff)) - 1.0;
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
};