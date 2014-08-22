#ifndef MATERIAL_H
#define	MATERIAL_H

#include <vecmath.h>

using namespace std;

class Material {
public:
	Material();
	Material(float KDiff[3], float KAmbient[3], double KSpec, double NSpec);
	Vector3f kDiff;
	Vector3f kAmbient;
	double kSpec;
	double nSpec;
};

#endif	/* MATERIAL_H */