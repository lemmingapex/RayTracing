#ifndef PRIMITIVE_H
#define	PRIMITIVE_H

#include <vecmath.h>

#include "Material.h"

class Primitive {
	public:
		Material m;
		virtual double Intersection(Vector3f viewPoint, Vector3f eyeRay) = 0;
		virtual Vector3f Normal(Vector3f viewPoint, Vector3f intersectionPoint) = 0;
};

#endif	/* PRIMITIVE_H */