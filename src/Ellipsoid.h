#ifndef ELLIPSOID_H
#define	ELLIPSOID_H

#include <vecmath.h>

#include "Material.h"
#include "Primitive.h"

class Ellipsoid: public Primitive {
protected:
	Vector3f center;
	// AxisLengths encoded
	Matrix3f A;

public:
	Ellipsoid()
	{ }

	Ellipsoid(float Center[3], float AxisLengths[3], Material M) {
		center = Vector3f(Center);
		m = M;
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
		double a = l.dot(matDir);
		double b = 2*(l.dot(matDiff));
		double c = (diff.dot(matDiff)) - 1.0;
		double d = b*b-4*a*c;

		// ignore complex solutions
		if(d>0) {
			// get both roots
			double t1=(-b-sqrt(d))/(2*a);
			double t2=(-b+sqrt(d))/(2*a);

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

#endif	/* ELLIPSOID_H */