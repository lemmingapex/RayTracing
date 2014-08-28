#ifndef ELLIPSOID_H
#define	ELLIPSOID_H

#include <vecmath.h>
#include <math.h>

#include "Material.h"
#include "Primitive.h"

#define PI 3.14159265

class Ellipsoid: public Primitive {
protected:
	Vector3f center;
	// AxisLengths encoded
	Matrix3f A;
	Matrix3f rotation;

public:
	Ellipsoid()
	{ }

	Ellipsoid(float Center[3], float AxisLengths[3], float theta, float phi, Material M) {
		center = Vector3f(Center);
		m = M;
		A = Matrix3f(1.0f/(AxisLengths[0]*AxisLengths[0]), 0.0f, 0.0f, 0.0f, 1.0f/(AxisLengths[1]*AxisLengths[1]), 0.0f, 0.0f, 0.0f, 1.0f/(AxisLengths[2]*AxisLengths[2]));
		// convert to radians
		theta = (theta * PI) / 180.0;
		phi = (phi * PI) / 180.0;
		rotation = Matrix3f(1.0f, 0.0f, 0.0f, 0.0f, cos(theta), -1.0*sin(theta), 0.0f, sin(theta), cos(theta)) * Matrix3f(cos(phi), 0.0f, sin(phi), 0.0f, 1.0f, 0.0f, -1.0*sin(phi), 0.0f, cos(phi));
	}

	virtual double Intersection(Vector3f viewPoint, Vector3f l) {
		double t = -1.0;

		// apply the rotaion matrix to A
		Matrix3f RTAT = rotation.transposed()*A*rotation;
		Vector3f diff = viewPoint - center;

		// ray ellipsoid intersection
		// ellipsoid:
		// (x - c)^T * A * (x - c) = 1

		// line:
		// x = t*l + l0

		// substitute and expand

		// a*t^2 + b*t + c = 0
		// t = -b+/-sqrt(b^2-4ac)/2a
		Vector3f matDir = RTAT*l;
		Vector3f matDiff = RTAT*diff;
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