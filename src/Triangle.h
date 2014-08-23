#include <vecmath.h>

#include "Material.h"
#include "Primitive.h"

class Triangle: public Primitive {
public:
	Vector3f a1, a2, a3;

	Triangle()
	{ }

	Triangle(float A1[3], float A2[3], float A3[3], Material M) {
		a1=Vector3f(A1);
		a2=Vector3f(A2);
		a3=Vector3f(A3);
		m=M;
	}

	virtual double Intersection(Vector3f viewPoint, Vector3f l) {
		double t = -1;
		// normal
		Vector3f n = (a2-a1).cross(a3-a1);
		
		// ray plane intersection
		// plane:
		// (p - p0) dot n = 0

		// line:
		// p = t*l + l0

		// substitute:
		// (t*l + l0 - p0) dot n = 0

		t=-1.0*(((viewPoint-a1).dot(n))/(l.dot(n)));
		if(t < 0) {
			return -1;
		}
		
		// find point on plane
		Vector3f p = viewPoint + (l*t);

		// check if point is inside triangle
		Vector3f v1 = (a1-p).cross(a2-p);
		Vector3f v2 = (a2-p).cross(a3-p);
		Vector3f v3 = (a3-p).cross(a1-p);

		// if the dot product of all of the vectors is the same (should be 1 or -1 depending on triangle orientation)
		// then all the vectors point in the same direction, and the point is inside the triangle.
		// don't check equivalence becasue of numerical precision issues, use < or >
		if((v1.dot(v2)<0)||(v2.dot(v3)<0)||(v3.dot(v1)<0)) {
			return -1;
		}
		return t;
	}

	virtual Vector3f Normal(Vector3f viewPoint, Vector3f intersectionPoint) {
		Vector3f normal = (a3-a1).cross(a2-a1);
		if((normal).dot(viewPoint-a1)<0) {
			normal=normal*-1.0;
		}
		return normal;
	}
};