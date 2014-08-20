#ifndef MATERIAL_H
#define	MATERIAL_H

using namespace std;

class Material {
public:
	Material();
	Material(double K_Diff[3], double K_Ambient[3], double K_Spec, double N_Spec);
	double k_diff_R, k_diff_G, k_diff_B;
	double k_ambient_R, k_ambient_G, k_ambient_B;
	double k_spec;
	double n_spec;
};

#endif	/* MATERIAL_H */