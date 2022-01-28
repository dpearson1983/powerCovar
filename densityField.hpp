#ifndef _DENSITY_FIELD_HPP_
#define _DENSITY_FIELD_HPP_

#include <vector>
#include <cmath>
#include "galaxy.hpp"

class densityField{
    std::vector<double> delta, L, dr, nbw;
    std::vector<int> N;

public:
    densityField(std::vector<int> N, std::vector<double> L);

    void binCIC(std::vector<galaxy> &gals, std::vector<double> r_min);

    double getDeltaR(int i, int j, int k);

    std::vector<double> getNbw();
};

densityField::densityField(std::vector<int> N, std::vector<double> L) {
    this->N = N;
    this->L = L;
    this->dr = {L[0]/double(N[0]), L[1]/double(N[1]), L[2]/double(N[2])};
    this->nbw = std::vector<double>(3);

    this->delta = std::vector<double>(N[0]*N[1]*N[2]);
}

void densityField::binCIC(std::vector<galaxy> &gals, std::vector<double> r_min) {
    double V = this->dr[0]*this->dr[1]*this->dr[2];
    for (int gal = 0; gal < gals.size(); ++gal) {
        std::vector<double> r = gals[gal].getCartPos();
        r[0] -= r_min[0];
        r[1] -= r_min[1];
        r[2] -= r_min[2];

        double w = gals[gal].getWeight();

        this->nbw[0] += w;
        this->nbw[1] += w*w;
        this->nbw[2] += gals[gal].getNbar()*w*w;

        double x_c = std::round(r[0]/this->dr[0])*this->dr[0];
        double y_c = std::round(r[1]/this->dr[1])*this->dr[1];
        double z_c = std::round(r[2]/this->dr[2])*this->dr[2];

        for (int i = -1; i <= 1; i += 2) {
            double x = r[0] + i*this->dr[0]/2.0;
            int i_p = int(x/this->dr[0]);
            for (int j = -1; j <= 1; j += 2) {
                double y = r[1] + j*this->dr[1]/2.0;
                int j_p = int(y/this->dr[1]);
                for (int k = -1; k <= 1; k += 2) {
                    double z = r[2] + k*this->dr[2]/2.0;
                    int k_p = int(z/this->dr[2]);
                    int index = k_p + this->N[2]*(j_p + this->N[1]*i_p);
                    double dV = std::abs(x - x_c)*std::abs(y - y_c)*std::abs(z - z_c);
                    this->delta[index] += w*dV/V;
                }
            }
        }
    }
}

double densityField::getDeltaR(int i, int j, int k) {
    int index = k + this->N[2]*(j + this->N[1]*i);
    return this->delta[index];
}

std::vector<double> densityField::getNbw() {
    return this->nbw;
}

#endif
