#include "Material.h"

using namespace std;

Material::Material() {
	
}

Material::Material(double KDiff[3], double KAmbient[3], double KSpec, double NSpec) {
	kDiff = Point(KDiff);
	kAmbient = Point(KAmbient);
	kSpec = KSpec;
	nSpec = NSpec;
}