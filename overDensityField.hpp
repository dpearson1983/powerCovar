#ifndef _OVER_DENSITY_FIELD_HPP_
#define _OVER_DENSITY_FIELD_HPP_

#include <vector>
#include <fftw3.h>
#include "densityField.hpp"

class overDensityField{
    std::vector<int> N;
    fftw_plan dr2dk, dk2dr;
    bool realSpace;

    void fftFreq();

    double sinc(double x);

public:
    std::vector<double> deltaR, L, dr, kx, ky, kz;
    std::vector<fftw_complex> deltaK;

    double shotNoise;

    overDensityField(std::vector<int> N, std::vector<double> L);

    void calculate(densityField &n_gals, densityField &n_rans);

    void fourierTransform();

    double getDeltaR(int i, int j, int k);

    std::vector<double> getDeltaK(int i, int j, int k);

    double getFreq(int i, int j, int k);

    std::vector<double> getWaveVec(int i, int j, int k);

    double gridCorrection(int i, int j, int k);

    ~overDensityField();
};

void overDensityField::fftFreq() {
    double dkx = (2.0*M_PI)/this->L[0];
    double dky = (2.0*M_PI)/this->L[1];
    double dkz = (2.0*M_PI)/this->L[2];

    for (int i = 0; i <= this->N[0]/2; ++i)
        this->kx[i] = i*dkx;
    for (int i = this->N[0]/2 + 1; i < this->N[0]; ++i)
        this->kx[i] = (i - this->N[0])*dkx;

    for (int i = 0; i <= this->N[1]/2; ++i)
        this->ky[i] = i*dky;
    for (int i = this->N[1]/2 + 1; i < this->N[1]; ++i)
        this->ky[i] = (i - this->N[1])*dky;

    for (int i = 0; i <= this->N[2]/2; ++i)
        this->kz[i] = i*dkz;
    for (int i = this->N[2]/2 + 1; i < this->N[2]; ++i)
        this->kz[i] = (i - this->N[2])*dkz;
}

double overDensityField::sinc(double x) {
    double val = 1.0;
    if (x != 0) {
        val = std::sin(x)/x;
    }
    return val;
}

overDensityField::overDensityField(std::vector<int> N, std::vector<double> L) {
    this->N = N;
    this->L = L;
    this->dr = {L[0]/double(N[0]), L[1]/double(N[1]), L[2]/double(N[2])};
    this->realSpace = true;
    this->kx = std::vector<double>(N[0]);
    this->ky = std::vector<double>(N[1]);
    this->kz = std::vector<double>(N[2]);

    std::cout << "    Getting wave numbers..." << std::endl;
    overDensityField::fftFreq();

    std::cout << "    Initializing arrays..." << std::endl;
    this->deltaR = std::vector<double>(N[0]*N[1]*N[2]);
    this->deltaK = std::vector<fftw_complex>(N[0]*N[1]*(N[2]/2 + 1));

    std::cout << "    Planning FFTs..." << std::endl;
    fftw_init_threads();
    fftw_plan_with_nthreads(4);
    fftw_import_wisdom_from_filename("fftw_wisdom.dat");
    this->dr2dk = fftw_plan_dft_r2c_3d(N[0], N[1], N[2], this->deltaR.data(), this->deltaK.data(),
                                       FFTW_MEASURE);
    this->dk2dr = fftw_plan_dft_c2r_3d(N[0], N[1], N[2], this->deltaK.data(), this->deltaR.data(),
                                       FFTW_MEASURE);
    fftw_export_wisdom_to_filename("fftw_wisdom.dat");

    std::cout << "    Rezeroing..." << std::endl;
    for (int i = 0; i < this->deltaR.size(); ++i) {
        this->deltaR[i] = 0.0;
        if (i < this->deltaK.size()) {
            this->deltaK[i][0] = 0.0;
            this->deltaK[i][1] = 0.0;
        }
    }
}

void overDensityField::calculate(densityField &n_gals, densityField &n_rans) {
    std::vector<double> gal_nbw = n_gals.getNbw();
    std::vector<double> ran_nbw = n_rans.getNbw();

    double alpha = gal_nbw[0]/ran_nbw[0];
    std::cout << "    alpha = " << alpha << std::endl;
    this->shotNoise = gal_nbw[1] + alpha*alpha*ran_nbw[1];
    std::cout << "    shotNoise = " << this->shotNoise << std::endl;

    for (int i = 0; i < this->N[0]; ++i) {
        for (int j = 0; j < this->N[1]; ++j) {
            for (int k = 0; k < this->N[2]; ++k) {
                int index = k + N[2]*(j + N[1]*i);
                this->deltaR[index] = n_gals.getDeltaR(i,j,k) - alpha*n_rans.getDeltaR(i,j,k);
            }
        }
    }
}

void overDensityField::fourierTransform() {
    if (this->realSpace) {
        fftw_execute(this->dr2dk);
        this->realSpace = false;
    } else {
        fftw_execute(this->dk2dr);
        for (int i = 0; i < this->deltaR.size(); ++i) {
            this->deltaR[i] /= (N[0]*N[1]*N[2]);
        }
        this->realSpace = true;
    }
}

double overDensityField::getDeltaR(int i, int j, int k) {
    int index = k + N[2]*(j + N[1]*i);
    return this->deltaR[index];
}

std::vector<double> overDensityField::getDeltaK(int i, int j, int k) {
    int index = k + (N[2]/2 + 1)*(j + N[1]*i);
    std::vector<double> val = {this->deltaK[index][0], this->deltaK[index][1]};
    return val;
}

double overDensityField::getFreq(int i, int j, int k) {
    return std::sqrt(this->kx[i]*this->kx[i] + this->ky[j]*this->ky[j] + this->kz[k]*this->kz[k]);
}

std::vector<double> overDensityField::getWaveVec(int i, int j, int k) {
    return {this->kx[i], this->ky[j], this->kz[k]};;
}

double overDensityField::gridCorrection(int i, int j, int k) {
    double sincx = overDensityField::sinc(this->kx[i]*this->L[0]/double(this->N[0]));
    double sincy = overDensityField::sinc(this->ky[j]*this->L[1]/double(this->N[1]));
    double sincz = overDensityField::sinc(this->kz[k]*this->L[2]/double(this->N[2]));
    double prodSinc = sincx*sincy*sincz;
    return 1.0/(prodSinc*prodSinc);
}

overDensityField::~overDensityField() {
    fftw_destroy_plan(this->dr2dk);
    fftw_destroy_plan(this->dk2dr);
    fftw_cleanup_threads();
}

#endif
