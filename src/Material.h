#ifndef MATERIAL_H
#define	MATERIAL_H

#include "Point.h"

using namespace std;

class Material {
public:
	Material();
	Material(double KDiff[3], double KAmbient[3], double KSpec, double NSpec);
	Point kDiff;
	Point kAmbient;
	double kSpec;
	double nSpec;
};

#endif	/* MATERIAL_H */