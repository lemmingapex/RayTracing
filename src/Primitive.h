#ifndef PRIMITIVE_H
#define	PRIMITIVE_H

#include "Material.h"
#include "Point.h"

class Primitive {
	public:
		Material m;
		virtual double Intersection(Point viewPoint, Point eyeRay) = 0;
		virtual Point Normal(Point viewPoint, Point intersectionPoint) = 0;
};

#endif	/* PRIMITIVE_H */