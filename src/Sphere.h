#ifndef SPHERE_H
#define	SPHERE_H

#include <vecmath.h>

#include "Material.h"
#include "Ellipsoid.h"

class Sphere: public Ellipsoid {
protected:
	double radius;

public:
	Sphere()
	{ }

	Sphere(float Center[3], double Radius, Material M) {
		center = Vector3f(Center);
		m = M;
		A = Matrix3f(1.0f/(Radius*Radius), 0.0f, 0.0f, 0.0f, 1.0f/(Radius*Radius), 0.0f, 0.0f, 0.0f, 1.0f/(Radius*Radius));
		radius = Radius;
	}

	/**
	* Although the ellipsoid intersection algorithm could easily be used, the degeneracies that 
	* occur in the sphere case save a considerable ammount of computation time.
	*/
 	virtual double Intersection(Vector3f viewPoint, Vector3f l) {
		double t = -1.0;
		Vector3f diff = viewPoint - center;

		// ray sphere intersection
		// sphere:
		// mag(x - c)^2 = r^2

		// line:
		// x = t*l + l0

		// substitute and expand

		// a*t^2 + b*t + c = 0
		// t = -b+/-sqrt(b^2-4ac)/2a
		double a = (l.abs())*(l.abs());
		double b = 2*(l.dot(diff));
		double c = ((diff.abs()*diff.abs())-(radius*radius));
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
};

#endif	/* SPHERE_H */