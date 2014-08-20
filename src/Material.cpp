#include "Material.h"

using namespace std;

Material::Material() {
	
}

Material::Material(double K_Diff[3], double K_Ambient[3], double K_Spec, double N_Spec) {
	k_diff_R=K_Diff[0];
	k_diff_G=K_Diff[1];
	k_diff_B=K_Diff[2];
	k_ambient_R=K_Ambient[0];
	k_ambient_G=K_Ambient[1];
	k_ambient_B=K_Ambient[2];
	k_spec=K_Spec;
	n_spec=N_Spec;
}