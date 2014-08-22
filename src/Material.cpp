#include "Material.h"

using namespace std;

Material::Material() {
	
}

Material::Material(float KDiff[3], float KAmbient[3], double KSpec, double NSpec) {
	kDiff = Vector3f(KDiff);
	kAmbient = Vector3f(KAmbient);
	kSpec = KSpec;
	nSpec = NSpec;
}