#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
// #include "functions_amu.hpp"
// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"
#include "resampling_new.hpp"
#include "global.hpp"
#include "fit_all.hpp"
// #include "do_analysis_charm.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <map>
constexpr double hbarc = 197.326963;



double read_single_Nobs(FILE* stream, int header_size, int Njack) {
    int Nobs;
    long int tmp;
    int s = header_size;

    // size_t i = fread(&Njack, sizeof(int), 1, stream);


    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);

    return Nobs;

}

data_single read_single_dataj(FILE* stream) {

    int Njack;
    int Nobs;


    data_single dj;
    dj.header.read_header_jack(stream);
    dj.header.print_header();
    dj.Nobs = read_single_Nobs(stream, dj.header.struct_size, dj.header.Njack);
    dj.Njack = dj.header.Njack;
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;

}

data_all read_all_the_files(std::vector<std::string> files, const char* resampling) {
    data_all jackall;
    jackall.resampling = resampling;
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    for (int i = 0; i < files.size();i++) {
        printf("reading: %s\n", files[i].c_str());
        FILE* f = open_file(files[i].c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[i] = read_single_dataj(f);
        jackall.en[i].resampling = resampling;

        fclose(f);

        // printing
        printf("file: %s\n", files[i].c_str());
        jackall.en[i].header.print_header();

    }
    return jackall;

}

int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];

    std::vector<std::string> files;

    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th1_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th2_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th3_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th4_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th6_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th7_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th8_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th9_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th9.5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_B64 = files.size();

    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th1_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th2_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th3_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th4_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th6_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th7_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th8_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th9_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th9.5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_B96 = files.size() - Ntheta_B64;

    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th1_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th2_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th3_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th4_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th5_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th6_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th7_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th8_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th9_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th9.5_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_C80 = files.size() - Ntheta_B64 - Ntheta_B96;

    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th1_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th2_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th3_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th4_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th5_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th6_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th7_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th8_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th9_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th9.5_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_D96 = files.size() - Ntheta_B64 - Ntheta_B96 - Ntheta_C80;

    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th1_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th2_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th3_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th4_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th5_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th6_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th7_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th8_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th9_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cE211.044.112_th9.5_t91_72.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_E112 = files.size() - Ntheta_B64 - Ntheta_B96 - Ntheta_C80 - Ntheta_D96;

     mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th1_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th2_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th3_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th4_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th5_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th6_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th7_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th8_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th9_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.48_th9.5_t48_36.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_B48 = files.size() - Ntheta_B64 - Ntheta_B96 - Ntheta_C80 - Ntheta_D96 - Ntheta_E112 ;

    std::vector<int> all_en(files.size());
    for (int e = 0; e < files.size(); e++) {
        all_en[e] = e;
    }

    ////////////////////////////////////////////////////////////////////////////
    //// read
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    //// generalised resampling
    ////////////////////////////////////////////////////////////////////////////
    data_all jackall = read_all_the_files(files, argv[1]);
    jackall.create_generalised_resampling();

    int Njack = jackall.en[0].Njack;
    if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        myres = new resampling_boot(Njack - 1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // fits sigma4
    /////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<std::string>th_map = { "1","2","3","4","5","6","7","8","9","9.5" };
    error(th_map.size() != Ntheta_B64, 1, "check the thetas", "");
    double** fit_result_dG = malloc_2<double>(Ntheta_B64, Njack);
    for (int th = 0;th < Ntheta_B64;th++) {
        fit_type fit_info;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Nvar = 3;
        fit_info.Njack = Njack;

        fit_info.myen.resize(6);
        fit_info.myen[0] = th;
        fit_info.myen[1] = fit_info.myen[0] + Ntheta_B64;
        fit_info.myen[2] = fit_info.myen[1] + Ntheta_B96;
        fit_info.myen[3] = fit_info.myen[2] + Ntheta_C80;
        fit_info.myen[4] = fit_info.myen[3] + Ntheta_D96;
        fit_info.myen[5] = fit_info.myen[4] + Ntheta_E112;

        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();// Nvar, entot, Njack
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e = 0;e < fit_info.myen.size();e++) {
            for (int e : fit_info.myen) {
                for (int j = 0;j < Njack;j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[58][j], 2);
                    fit_info.x[1][count][j] = pow(jackall.en[e].header.thetas[0] * M_PI / (jackall.en[e].header.L * jackall.en[e].jack[58][j] * 1000.0 / hbarc), 2);// GeV
                    fit_info.x[2][count][j] = jackall.en[e].header.L;
                }
                count++;
            }
        }
        fit_info.function = rhs_const;
        fit_info.corr_id = { 69 };//
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "dGamma_th%s_const", th_map[th].c_str());
        fit_result dGamma_th = fit_all_data(argv, jackall, lhs_identity, fit_info, namefit);
        fit_info.band_range = { 0,0.0081 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", dGamma_th, dGamma_th, 0, fit_info.myen.size() - 1, 0.001);

        fit_info.restore_default();
        for (int j = 0;j < Njack;j++)
            fit_result_dG[th][j] = dGamma_th.P[0][j];

        dGamma_th.clear();
    }
    double** Cov = myres->comp_cov(Ntheta_B64, fit_result_dG);
    char namecov[NAMESIZE];
    mysprintf(namecov, NAMESIZE, "%s/dGamma_const_corr.dat", argv[3]);
    FILE* fcov = open_file(namecov, "w+");
    for (int i = 0;i < Ntheta_B64;i++) {
        for (int j = 0;j < Ntheta_B64;j++) {
            fprintf(fcov, "%-22.4g", Cov[i][j] / sqrt(Cov[i][i] * Cov[j][j]));
        }
        fprintf(fcov, "\n");
    }
    fclose(fcov);
    free_2(Ntheta_B64, Cov);
    free_2(Ntheta_B64, fit_result_dG);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // fits sigma4 Z0 part
    /////////////////////////////////////////////////////////////////////////////////////////////////
    error(th_map.size() != Ntheta_B64, 1, "check the thetas", "");
    fit_result_dG = malloc_2<double>(Ntheta_B64, Njack);
    for (int th = 0;th < Ntheta_B64;th++) {
        fit_type fit_info;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Nvar = 3;
        fit_info.Njack = Njack;

        fit_info.myen.resize(6);
        fit_info.myen[0] = th;
        fit_info.myen[1] = fit_info.myen[0] + Ntheta_B64;
        fit_info.myen[2] = fit_info.myen[1] + Ntheta_B96;
        fit_info.myen[3] = fit_info.myen[2] + Ntheta_C80;
        fit_info.myen[4] = fit_info.myen[3] + Ntheta_D96;
        fit_info.myen[5] = fit_info.myen[4] + Ntheta_E112;


        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();// Nvar, entot, Njack
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e = 0;e < fit_info.myen.size();e++) {
            for (int e : fit_info.myen) {
                for (int j = 0;j < Njack;j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[58][j], 2);
                    fit_info.x[1][count][j] = pow(jackall.en[e].header.thetas[0] * M_PI / (jackall.en[e].header.L * jackall.en[e].jack[58][j] * 1000.0 / hbarc), 2);// GeV
                    fit_info.x[2][count][j] = jackall.en[e].header.L;
                }
                count++;
            }
        }
        fit_info.function = rhs_const;
        fit_info.corr_id = { 71 };//
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "dGamma_th%s_part_const", th_map[th].c_str());
        fit_result dGamma_th = fit_all_data(argv, jackall, lhs_identity, fit_info, namefit);
        fit_info.band_range = { 0,0.0081 };
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", dGamma_th, dGamma_th, 0, fit_info.myen.size() - 1, 0.001);

        fit_info.restore_default();
        for (int j = 0;j < Njack;j++)
            fit_result_dG[th][j] = dGamma_th.P[0][j];

        dGamma_th.clear();
    }
    Cov = myres->comp_cov(Ntheta_B64, fit_result_dG);
    mysprintf(namecov, NAMESIZE, "%s/dGamma_const_corr.dat", argv[3]);
    
    free_2(Ntheta_B64, Cov);
    free_2(Ntheta_B64, fit_result_dG);

}