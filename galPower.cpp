#include <iostream>
#include "cosmology.hpp"
#include "galaxy.hpp"
#include "densityField.hpp"
#include "overDensityField.hpp"
#include "fitsRead.hpp"
#include "powerSpectrum.hpp"
#include <omp.h>

int main() {
    double H_0 = 100.0; // Hubble constant in km/s/Mpc
    double Omega_M = 0.3111;
    double Omega_L = 0.6889;

    std::vector<int> N = {512,1024,512};
    std::vector<double> L = {1792, 3584, 1792};
    std::vector<double> r_min(3), r_max(3);

    cosmology cosmo(H_0, Omega_M, Omega_L);

    std::cout.precision(15);
    std::cout << "Reading galaxy file..." << std::endl;
    std::vector<galaxy> gals = readGalFitsFile("galaxy_DR12v5_CMASS_North.fits", cosmo, r_min, r_max);
    std::cout << "    Minimum position: (" << r_min[0] << ", " << r_min[1] << ", " << r_min[2] << ")\n";
    std::cout << "    Maximum position: (" << r_max[0] << ", " << r_max[1] << ", " << r_max[2] << ")\n";
    std::cout << "    Extent of box:    (" << r_max[0] - r_min[0] << ", " << r_max[1] - r_min[1] << ", "
              << r_max[2] - r_min[2] << ")\n";

    std::cout << "Reading randoms file..." << std::endl;
    std::vector<galaxy> rans = readRanFitsFile("random0_DR12v5_CMASS_North.fits", cosmo, r_min, r_max);

    std::cout << "    Minimum position: (" << r_min[0] << ", " << r_min[1] << ", " << r_min[2] << ")\n";
    std::cout << "    Maximum position: (" << r_max[0] << ", " << r_max[1] << ", " << r_max[2] << ")\n";
    std::cout << "    Extent of box:    (" << r_max[0] - r_min[0] << ", " << r_max[1] - r_min[1] << ", "
              << r_max[2] - r_min[2] << ")\n";

    r_min[0] -= (L[0] - r_max[0] + r_min[0])/2;
    r_min[1] -= (L[1] - r_max[1] + r_min[1])/2;
    r_min[2] -= (L[2] - r_max[2] + r_min[2])/2;
    std::cout << "Minimum position: (" << r_min[0] << ", " << r_min[1] << ", " << r_min[2] << ")\n";

    std::cout << "Binning galaxies..." << std::endl;
    densityField n_gal(N, L);
    n_gal.binCIC(gals, r_min);


    std::cout << "Binning randoms..." << std::endl;
    densityField n_ran(N, L);
    n_ran.binCIC(rans, r_min);

    std::cout << "Setting up overdensity field..." << std::endl;
    overDensityField delta(N, L);
    std::cout << "Calculating overdensity..." << std::endl;
    delta.calculate(n_gal, n_ran);

    std::cout << "Fourier transforming..." << std::endl;
    delta.fourierTransform();

    std::cout << "Setting up power spectrum..." << std::endl;
    powerSpectrum Pk(0.012, 0.3, 0.008);

    double start = omp_get_wtime();
    std::cout << "Calculating power spectrum..." << std::endl;
    std::cout << "    normalization = " << n_gal.getNbw()[2] << std::endl;
    Pk.calculate(delta, N, n_gal.getNbw()[2]);
    double end = omp_get_wtime();
    std::cout << "Time to compute power spectrum: " << end - start << std::endl;

    std::cout << "Writing power spectrum file..." << std::endl;
    Pk.writeFile("bossPk.dat");

    std::cout << "Getting bootstrapped covariance..." << std::endl;
    std::vector<std::vector<double>> C = Pk.getParExCov(delta, N, 2048, n_gal.getNbw()[2]);

    std::cout << "Writing correlation matrix to file..." << std::endl;
    std::ofstream fout("ParExCor.dat");
    std::ofstream dout("ParExCov.dat");
    fout.precision(15);
    dout.precision(15);
    for (int i = 0; i < C.size(); ++i) {
        for (int j = 0; j < C[i].size(); ++j) {
            fout.width(25);
            dout.width(25);
            fout << C[i][j]/sqrt(C[i][i]*C[j][j]);
            dout << C[i][j];
        }
        fout << "\n";
        dout << "\n";
    }
    fout.close();
    dout.close();

    return 0;
}
