#include "Material.h"
#include "Point.h"
#include "Primitive.h"

class Triangle: public Primitive {
public:
	Point a1, a2, a3;

	Triangle()
	{ }

	Triangle(double A1[3], double A2[3], double A3[3], Material M) {
		a1=Point(A1);
		a2=Point(A2);
		a3=Point(A3);
		m=M;
	}

	virtual double Intersection(Point viewPoint, Point eyeRay) {
		double t = -1;
		// normal
		Point n = (a2-a1).cross(a3-a1);
		
		// ray plane intersection
		// plane:
		// (p - p0) dot n = 0

		// line:
		// p = t*l + l0

		// substitute:
		// (t*l + l0 - p0) dot n = 0

		t=-1*(((viewPoint-a1).dot(n))/(eyeRay.dot(n)));
		if(t<0)
		{
			return -1;
		}
		
		// find point on plane
		Point p = viewPoint + (eyeRay*t);

		// check if point is inside triangle
		Point v1 = (a1-p).cross(a2-p);
		Point v2 = (a2-p).cross(a3-p);
		Point v3 = (a3-p).cross(a1-p);

		// if the dot product of all of the vectors is the same (should be 1 or -1 depending on triangle orientation)
		// then all the vectors point in the same direction, and the point is inside the triangle.
		// don't check equivalence becasue of numerical precision issues, use < or >
		if((v1.dot(v2)<0)||(v2.dot(v3)<0)||(v3.dot(v1)<0)) {
			return -1;
		}
		return t;
	}

	virtual Point Normal(Point viewPoint, Point intersectionPoint) {
		Point normal = (a3-a1).cross(a2-a1);
		if((normal).dot(viewPoint-a1)<0) {
			normal=normal*-1.0;
		}
		return normal;
	}
};