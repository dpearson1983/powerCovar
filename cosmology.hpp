#ifndef _COSMOLOGY_HPP_
#define _COSMOLOGY_HPP_

#include <cmath>
#include <gsl/gsl_integration.h>

double c = 299792.458; // Speed of light in km/s

class cosmology{
    double H_0, Omega_M, Omega_L;
    gsl_integration_workspace *w;

public:
    cosmology(double H_0, double Omega_M, double Omega_L);

    double redshiftToDist(double z);

    ~cosmology();
};

double f(double x, void *params) {
    double *p = (double *)params;
    return 1/(std::sqrt(p[0]*(1.0 + x)*(1.0 + x)*(1.0 + x) + p[1]));
}

cosmology::cosmology(double H_0, double Omega_M, double Omega_L) {
    this->H_0 = H_0;
    this->Omega_M = Omega_M;
    this->Omega_L = Omega_L;
    this->w = gsl_integration_workspace_alloc(1000000);
}

double cosmology::redshiftToDist(double z) {
    gsl_function F;
    double p[] = {this->Omega_M, this->Omega_L};
    F.function = &f;
    F.params = &p;
    double result, error;
    gsl_integration_qags(&F, 0.0, z, 1E-6, 1E-6, 1000000, this->w, &result, &error);
    return c*result/this->H_0;
}

cosmology::~cosmology() {
    gsl_integration_workspace_free(this->w);
}
#endif
