#ifndef _POWER_SPECTRUM_HPP_
#define _POWER_SPECTRUM_HPP_

#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <omp.h>
#include "overDensityField.hpp"

class powerSpectrum{
    std::vector<double> P, k;
    std::vector<int> Nk;
    double k_min, k_max, delta_k;
    int N;

public:
    powerSpectrum(double k_min, double k_max, double delta_k);

    void calculate(overDensityField &delta, std::vector<int> N, double norm);

    void writeFile(std::string file);

    std::vector<std::vector<double>> getBootstrapCov(overDensityField &delta, std::vector<int> N,
                                                     int N_bs, double norm);

    std::vector<std::vector<double>> getRawCov(overDensityField &delta, std::vector<int> N,
                                               int N_s, double Norm);

    std::vector<std::vector<double>> getParExCov(overDensityField &delta, std::vector<int> N,
                                                 int N_s, double norm);
};

powerSpectrum::powerSpectrum(double k_min, double k_max, double delta_k) {
    int N = (k_max - k_min)/delta_k;
    this->N = N;
    this->k_min = k_min;
    this->k_max = k_max;
    this->delta_k = delta_k;
    this->P = std::vector<double>(N);
    this->k = std::vector<double>(N);
    this->Nk = std::vector<int>(N);
}

void powerSpectrum::calculate(overDensityField &delta, std::vector<int> N, double norm) {
    std::cout << "    Binning frequencies..." << std::endl;
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            for (int k = 0; k <= N[2]/2; ++k) {
                double k_mag = delta.getFreq(i,j,k);

                if (k_mag >= this->k_min && k_mag < this->k_max) {
                    int bin = (k_mag - this->k_min)/this->delta_k;
                    std::vector<double> dk = delta.getDeltaK(i,j,k);
                    this->P[bin] += (dk[0]*dk[0] + dk[1]*dk[1])*delta.gridCorrection(i,j,k) - delta.shotNoise;
                    this->k[bin] += k_mag;
                    this->Nk[bin] += 1;
                }
            }
        }
    }

    std::cout << "    Normalizing..." << std::endl;
    for (int i = 0; i < this->P.size(); ++i) {
        if (this->Nk[i] > 0) {
            P[i] /= (this->Nk[i]*norm);
            k[i] /= this->Nk[i];
        }
    }
}

void powerSpectrum::writeFile(std::string file) {
    std::ofstream fout(file);
    fout.precision(15);
    for (int i = 0; i < this->P.size(); ++i) {
        fout << this->k[i] << " " << this->P[i] << "\n";
    }
    fout.close();
}

std::vector<std::vector<double>> powerSpectrum::getBootstrapCov(overDensityField &delta,
                                                                std::vector<int> N, int N_bs,
                                                                double norm) {
    std::vector<std::vector<double>> cov(this->N, std::vector<double>(this->N));
    std::vector<std::vector<double>> Pks(this->N);
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            for (int k = 0; k <= N[2]/2; ++k) {
                double k_mag = delta.getFreq(i,j,k);

                if (k_mag >= this->k_min && k_mag < this->k_max) {
                    int bin = (k_mag - this->k_min)/this->delta_k;
                    std::vector<double> dk = delta.getDeltaK(i,j,k);
                    double Pk= (dk[0]*dk[0] + dk[1]*dk[1])*delta.gridCorrection(i,j,k) - delta.shotNoise;
                    Pks[bin].push_back(Pk);
                }
            }
        }
    }

    std::random_device seeder;
    std::mt19937_64 gen(seeder());

    for (int i = 0; i < N_bs; ++i) {
        std::vector<double> P_bs(this->N);
        for (int j = 0; j < this->N; ++j) {
            std::uniform_int_distribution<int> dist(0, this->Nk[j]);
            for (int k = 0; k < this->Nk[j]; ++k) {
                P_bs[j] += Pks[j][dist(gen)];
            }
            P_bs[j] /= (Nk[j]*norm);
        }
        for (int j = 0; j < this->N; ++j) {
            for (int k = 0; k < this->N; ++k) {
                cov[j][k] += (P_bs[j] - this->P[j])*(P_bs[k] - this->P[k])/double(N_bs - 1);
            }
        }
    }

    return cov;
}

std::vector<std::vector<double>> powerSpectrum::getRawCov(overDensityField &delta, std::vector<int> N,
                                                          int N_s, double norm) {
    std::vector<std::vector<double>> cov(this->N, std::vector<double>(this->N));
    std::vector<std::vector<double>> Pks(this->N);
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            for (int k = 0; k <= N[2]/2; ++k) {
                double k_mag = delta.getFreq(i,j,k);

                if (k_mag >= this->k_min && k_mag < this->k_max) {
                    int bin = (k_mag - this->k_min)/this->delta_k;
                    std::vector<double> dk = delta.getDeltaK(i,j,k);
                    double Pk= (dk[0]*dk[0] + dk[1]*dk[1])*delta.gridCorrection(i,j,k) - delta.shotNoise;
                    Pks[bin].push_back(Pk/norm);
                }
            }
        }
    }

    std::random_device seeder;
    std::mt19937_64 gen(seeder());

    for (int i = 0; i < this->N; ++i) {
        std::uniform_int_distribution<int> i_dist(0, this->Nk[i]-1);
        for (int j = 0; j < this->N; ++j) {
            std::uniform_int_distribution<int> j_dist(0, this->Nk[j]-1);
            for (int s = 0; s < N_s; ++s) {
                if (i != j) {
                    cov[i][j] += (Pks[i][i_dist(gen)] - this->P[i])*(Pks[j][j_dist(gen)] - this->P[j])/(N_s-1);
                } else {
                    int val = i_dist(gen);
                    cov[i][i] += (Pks[i][val] - this->P[i])*(Pks[i][val] - this->P[i])/(N_s-1);
                }
            }
        }
    }
    return cov;
}

std::vector<std::vector<double>> powerSpectrum::getParExCov(overDensityField &delta, std::vector<int> N,
                                                          int N_s, double norm) {
    std::vector<std::vector<double>> cov(this->N, std::vector<double>(this->N));
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_int_distribution<int> xdist(0,N[0]-1);
    std::uniform_int_distribution<int> ydist(0,N[1]-1);
    std::uniform_int_distribution<int> zdist(0,N[2]-1);
    double theta_min = std::cos(1.0*M_PI/180.0);
    double theta_max = std::cos(179*M_PI/180.0);
    for (int s = 0; s < N_s; ++s) {
        std::cout << "    Sample #" << s+1 << std::endl;
        std::vector<double> Pi(this->N);
        std::vector<double> Ni(this->N);
        int i_ex = xdist(gen);
        int j_ex = ydist(gen);
        int l_ex = zdist(gen);
        std::vector<double> k_ex = delta.getWaveVec(i_ex, j_ex, l_ex);
        double kex_mag = delta.getFreq(i_ex, j_ex, l_ex);
        for (int i = 0; i < N[0]; ++i) {
            double kx = delta.kx[i];
            for (int j = 0; j < N[1]; ++j) {
                double ky = delta.ky[j];
                for (int k = 0; k <= N[2]/2; ++k) {
                    double k_mag = std::sqrt(kx*kx + ky*ky + delta.kz[k]*delta.kz[k]);
                    double kdot = kx*k_ex[0] + ky*k_ex[1] + delta.kz[k]*k_ex[2];
                    double theta = kdot/(k_mag*kex_mag);

                    if (k_mag >= this->k_min && k_mag < this->k_max && theta < theta_min && theta > theta_max) {
                        int bin = (k_mag - this->k_min)/this->delta_k;
                        int index = k + (N[2]/2 + 1)*(j + N[1]*i);
                        Pi[bin] += ((delta.deltaK[index][0]*delta.deltaK[index][0] + delta.deltaK[index][1]*delta.deltaK[index][1])*delta.gridCorrection(i,j,k) - delta.shotNoise);
                        Ni[bin] += 1;
                    }
                }
            }
        }

        for (int i = 0; i < this->N; ++i) {
            if (Ni[i] > 0) {
                Pi[i] /= (norm*Ni[i]);
            }
        }
        double end = omp_get_wtime();

        for (int i = 0; i < this->N; ++i) {
            for (int j = 0; j < this->N; ++j) {
                cov[i][j] += (Pi[i] - this->P[i])*(Pi[j] - this->P[j])/(N_s - 1);
            }
        }
    }
    return cov;
}

#endif
