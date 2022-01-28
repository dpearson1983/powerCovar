#ifndef _FITSREAD_HPP_
#define _FITSREAD_HPP_

#include <vector>
#include <string>
#include <CCfits/CCfits>
#include "galaxy.hpp"
#include "cosmology.hpp"

std::vector<galaxy> readGalFitsFile(std::string file, cosmology &cosmo, std::vector<double> &r_min,
                                    std::vector<double> &r_max) {
    //std::vector<std::string> primaryKeys;
    //primaryKeys.push_back(std::string("PRIMARY"));
    std::unique_ptr<CCfits::FITS> pInFile(new CCfits::FITS(file, CCfits::Read, false));

    CCfits::ExtHDU &table = pInFile->extension(1);

    std::vector<double> RA, DEC, RED, w_sys, w_noz, w_cp, w_fkp, nbar;

    long start = 1L;
    long end = table.rows();

    table.column("RA").read(RA, start, end);
    table.column("DEC").read(DEC, start, end);
    table.column("Z").read(RED, start, end);
    table.column("WEIGHT_FKP").read(w_fkp, start, end);
    table.column("WEIGHT_NOZ").read(w_noz, start, end);
    table.column("WEIGHT_CP").read(w_cp, start, end);
    table.column("WEIGHT_SYSTOT").read(w_sys, start, end);
    table.column("NZ").read(nbar, start, end);

    r_min[0] = 1E6;
    r_min[1] = 1E6;
    r_min[2] = 1E6;
    r_max[0] = -1E6;
    r_max[1] = -1E6;
    r_max[2] = -1E6;
    std::vector<galaxy> gals;
    for (long i = 0; i < RA.size(); ++i) {
        if (RED[i] >= 0.43 && RED[i] <= 0.7) {
            double w = w_fkp[i]*w_sys[i]*(w_noz[i] + w_cp[i] - 1);
            galaxy gal(RA[i], DEC[i], RED[i], w, nbar[i]);
            gal.astroToCart(cosmo);
            gals.push_back(gal);

            std::vector<double> pos = gal.getCartPos();
            for (int j = 0; j < 3; ++j) {
                if (pos[j] < r_min[j]) r_min[j] = pos[j];
                if (pos[j] > r_max[j]) r_max[j] = pos[j];
            }
        }
    }

    return gals;
}

std::vector<galaxy> readRanFitsFile(std::string file, cosmology &cosmo, std::vector<double> &r_min,
                                    std::vector<double> &r_max) {
    //std::vector<std::string> primaryKeys;
    //primaryKeys.push_back(std::string("PRIMARY"));
    std::unique_ptr<CCfits::FITS> pInFile(new CCfits::FITS(file, CCfits::Read, false));

    CCfits::ExtHDU &table = pInFile->extension(1);

    std::vector<double> RA, DEC, RED, w_fkp, nbar;

    long start = 1L;
    long end = table.rows();

    table.column("RA").read(RA, start, end);
    table.column("DEC").read(DEC, start, end);
    table.column("Z").read(RED, start, end);
    table.column("WEIGHT_FKP").read(w_fkp, start, end);
    table.column("NZ").read(nbar, start, end);

    r_min[0] = 1E6;
    r_min[1] = 1E6;
    r_min[2] = 1E6;
    r_max[0] = -1E6;
    r_max[1] = -1E6;
    r_max[2] = -1E6;
    std::vector<galaxy> rans;
    for (long i = 0; i < RA.size(); ++i) {
        if (RED[i] >= 0.43 && RED[i] <= 0.7) {
            galaxy ran(RA[i], DEC[i], RED[i], w_fkp[i], nbar[i]);
            ran.astroToCart(cosmo);
            rans.push_back(ran);

            std::vector<double> pos = ran.getCartPos();
            for (int j = 0; j < 3; ++j) {
                if (pos[j] < r_min[j]) r_min[j] = pos[j];
                if (pos[j] > r_max[j]) r_max[j] = pos[j];
            }
        }
    }

    return rans;
}

#endif
