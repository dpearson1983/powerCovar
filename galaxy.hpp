#ifndef _GALAXY_HPP_
#define _GALAXY_HPP_

#include <vector>
#include <cmath>
#include "cosmology.hpp"

class galaxy{
    std::vector<double> astro, cart;
    double w, nbar;

public:
    galaxy(double ra, double dec, double red, double w, double nbar);

    void astroToCart(cosmology &cosmo);

    std::vector<double> getCartPos();

    std::vector<double> getAstroPos();

    double getWeight();

    double getNbar();
};

galaxy::galaxy(double ra, double dec, double red, double w, double nbar) {
    this->w = w;
    this->nbar = nbar;
    this->astro = std::vector<double>(3);
    this->cart = std::vector<double>(3);

    this->astro[0] = ra;
    this->astro[1] = dec;
    this->astro[2] = red;
}

void galaxy::astroToCart(cosmology &cosmo) {
    double r = cosmo.redshiftToDist(this->astro[2]);
    this->cart[0] = r*std::cos(this->astro[1]*M_PI/180.0)*std::cos(this->astro[0]*M_PI/180.0);
    this->cart[1] = r*std::cos(this->astro[1]*M_PI/180.0)*std::sin(this->astro[0]*M_PI/180.0);
    this->cart[2] = r*std::sin(this->astro[1]*M_PI/180.0);
}

std::vector<double> galaxy::getCartPos() {
    return this->cart;
}

std::vector<double> galaxy::getAstroPos() {
    return this->astro;
}

double galaxy::getWeight() {
    return this->w;
}

double galaxy::getNbar() {
    return this->nbar;
}

#endif
