#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"
#include "fit_all.hpp"
#include "functions_heavy_quarks.hpp"
#include "tower.hpp"
#include "HLT.hpp"


struct kinematic kinematic_2pt;
constexpr double hbarc = 197.326963;
char** argv_to_options(char** argv) {
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[2]);         // path
    mysprintf(option[4], NAMESIZE, argv[6]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
    return option;
}

void init_global_head(generic_header head) {
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[1];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = head.mus[1];
    file_head.k[3] = head.mus[1];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

void read_twopt(FILE* stream, double*** to_write, generic_header head) {
    // write your function to read the data
    int fi = 0;
    int id;
    int i = fread(&id, sizeof(int), 1, stream);
    for (int k = 0; k < head.ncorr; k++) {
        for (int t = 0; t < head.T; t++) {
            fi += fread(to_write[k][t], sizeof(double), 2, stream);
        }
    }
}

int main(int argc, char** argv) {
    error(argc != 7, 1, "main ",
        "usage:./heavy_quarks -p path file -bin $bin  jack/boot \n separate "
        "path and file please");

    char** option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    FILE* infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head;
    head.read_header(infile);
    head.print_header();
    init_global_head(head);
    //////////////////////////////////// setup jackboot and binning
    int confs = head.Njack;
    int bin = atoi(argv[5]);
    int Neff = confs / bin;
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = 500;//(Neff * 2 + 1);
        myres = new resampling_boot(Njack - 1);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    // now Njack need to be the number of jacks
    head.Njack = Njack;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4],
        option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    head.write_header(jack_file);
    fflush(jack_file);

    // FILE* rjack_file = open_file(namefile, "r");
    // head.read_header(rjack_file);
    // head.print_header();

    //////////////////////////////////// confs
    int ncorr_max = head.ncorr + 56;
    double**** data = calloc_corr(confs, head.ncorr, head.T);
    double**** data_4pt = calloc_corr(confs, ncorr_max, head.T);

    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < confs; iconf++) {
        read_twopt(infile, data[iconf], head);
    }
    //////////////////////////////////////////////////////////////
    // build 4pt
    //////////////////////////////////////////////////////////////
    for (int iconf = 0; iconf < confs; iconf++) {
        for (int t = 0; t < head.T;t++) {
            for (int n = 0; n < head.ncorr;n++) {
                data_4pt[iconf][n][t][0] = data[iconf][n][t][0];
                data_4pt[iconf][n][t][1] = data[iconf][n][t][1];
            }
        }

    }
    int TJW, TDs, myerr, myseed;
    line_read_param(option, "TJW", TJW, myerr, myseed, namefile_plateaux);
    line_read_param(option, "TDs", TDs, myerr, myseed, namefile_plateaux);

    int t0 = TJW;
    for (int j = 0; j < confs;j++) {
        for (int t = 0; t < head.T;t++) {
            for (int g1 = 0;g1 < 8;g1++) {
                for (int g2 = 0;g2 < 8;g2++) {
                    int idg1 = 2 * head.gammas.size() + g1;
                    int idg2 = 3 * head.gammas.size() + g2;
                    int id_final = 2 * head.gammas.size() + g1 + g2 * 8;
                    std::complex<double> cp(data[j][idg1][t][0], data[j][idg1][t][1]);
                    std::complex<double> cm(data[j][idg2][t0][0], data[j][idg2][t0][1]);
                    // we put a minus sign because respect to the connected 
                    // part there is a sign change due to the fermion loops
                    std::complex<double>  c4 = - cp * cm; 
                    data_4pt[j][id_final][t][0] = c4.real();
                    data_4pt[j][id_final][t][1] = c4.imag();
                }
            }

        }
    }

    //////////////////////////////////////////////////////////////
    // binning
    //////////////////////////////////////////////////////////////
    double**** data_bin = binning(confs, ncorr_max, head.T, data_4pt, bin);
    double**** conf_jack = myres->create(Neff, ncorr_max, head.T, data_bin);
    free_corr(Neff, ncorr_max, head.T, data_bin);
    free_corr(confs, head.ncorr, head.T, data);
    free_corr(confs, ncorr_max, head.T, data_4pt);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3],
        option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3],
        option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3],
        option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3],
        option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_kernel", option[3],
        option[6]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_AoverB", option[3],
        option[6]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");

    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
            M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);
    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;

    ///////////////////// symmetrize
    // for (int i = 0; i < head.ncorr;i++) {
    //     symmetrise_jackboot(Njack, i, head.T, conf_jack);
    // }

    int id_Ds_ss_2pt = 0;
    int id_Ds_ls_2pt = 0 + head.gammas.size();
    int id_K_ss_2pt = 0 + 2 * head.gammas.size();
    int id_K_ls_2pt = 0 + 3 * head.gammas.size();

    symmetrise_jackboot(Njack, id_Ds_ss_2pt, head.T, conf_jack, 1);
    symmetrise_jackboot(Njack, id_Ds_ls_2pt, head.T, conf_jack, 1);
    // symmetrise_jackboot(Njack, id_K_ss_2pt, head.T, conf_jack, 1);
    // symmetrise_jackboot(Njack, id_K_ls_2pt, head.T, conf_jack, 1);


    ////////////////////


    //////////////////////////////
    // start fitting
    /////////////////////////////
    corr_counter = -1;
    double* M_PS, * M_Ds;
    M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_Ds_ls_2pt, "M_{Ds}_sl", M_eff_T, jack_file);
    check_correlatro_counter(0);
    free(M_PS);

    M_Ds = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_Ds_ss_2pt, "M_{Ds}_ss", M_eff_T, jack_file);
    check_correlatro_counter(1);


    M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_K_ss_2pt, "M_{K}_sl", M_eff_T, jack_file);
    check_correlatro_counter(2);
    free(M_PS);

    M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_K_ls_2pt, "M_{K}_ss", M_eff_T, jack_file);
    check_correlatro_counter(3);
    free(M_PS);


    ///////////////////////////////////////////////
    /// compute ZDs smeared smeared
    struct fit_type fit_info;
    struct fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
    fit_info.ext_P[0] = M_Ds;
    // for (int j = 0; j < fit_info.Njack; j++) {
    //     fit_info.ext_P[0][j] = M_Ds[j];
    // }
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = head.T;
    fit_info.corr_id = { id_Ds_ss_2pt };

    // c++ 0 || r 1
    struct fit_result ZDs = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_ZPS, "Z_{Ds}_ss", fit_info, jack_file);
    // free_fit_result(fit_info, fit_out);
    check_correlatro_counter(4);
    fit_info.restore_default();

    // getting the Z factor with a 2 parameters fit
    fit_info.Nvar = 1;
    fit_info.Npar = 2;
    fit_info.N = 1;
    fit_info.Njack = Njack;

    fit_info.function = rhs_2fit_par;
    fit_info.linear_fit = false;
    fit_info.T = head.T;
    fit_info.corr_id = { id_Ds_ss_2pt };

    // c++ 0 || r 1
    struct fit_result f2parZ = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_2fit_par, "2parZ_{Ds}_ss", fit_info, jack_file);
    free_fit_result(fit_info, f2parZ);
    check_correlatro_counter(5);
    fit_info.restore_default();


    //// struct to ignore the plateaux range
    fit_type tmp_info;
    tmp_info.codeplateaux = true;
    tmp_info.tmin = 20;
    tmp_info.tmax = 20;
    //////////////////////
    std::vector<std::vector<int>> id_Mmunu(4, std::vector<int>(4));
    int ncorr_new = ncorr_max;

    double ZA, ZV, dZA, dZV;

    line_read_param(option, "ZA", ZA, dZA, myseed, namefile_plateaux);
    double* ZA_jack = myres->create_fake(ZA, dZA, myseed);
    line_read_param(option, "ZV", ZV, dZV, myseed, namefile_plateaux);
    double* ZV_jack = myres->create_fake(ZV, dZV, myseed);
    double a_val, da_val;
    line_read_param(option, "a", a_val, da_val, myseed, namefile_plateaux);
    double* a_fm_jack = myres->create_fake(a_val, da_val, myseed);
    double* a_MeV_jack = (double*)malloc(sizeof(double) * Njack);
    double* a_GeV_jack = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0; j < Njack;j++) {
        a_MeV_jack[j] = a_fm_jack[j] / hbarc;
        a_GeV_jack[j] = a_MeV_jack[j] * 1000.0;
    }

    for (int mu = 0; mu < 4;mu++) {  // C_munu   mu runs faster // mu is given by the gamma // nu by the insertion
        for (int nu = 0; nu < 4;nu++) {


            fit_info.N = 1;
            fit_info.Njack = Njack;

            fit_info.T = head.T;
            fit_info.myen = { TDs, TJW ,mu, nu };
            fit_info.n_ext_P = 3;
            fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
            fit_info.ext_P[0] = M_Ds;
            fit_info.ext_P[1] = ZA_jack;
            fit_info.ext_P[2] = ZV_jack;
            // fit_info.ext_P[1] = (double*)malloc(sizeof(double) * Njack);
            // for (size_t j = 0; j < Njack; j++) {
            //     fit_info.ext_P[1][j] = 
            // }


            //                  2pt_files      +   A_4pt              +  gamma=Vmu + insertions Vmu     
            int id_VV = 2 * head.gammas.size() + (mu + 4) + (nu + 4) * 8;
            //                   2pt_files       +  gamma=Amu + insertions Amu     
            int id_AA = 2 * head.gammas.size() + (mu + 0) + (nu + 0) * 8;
            //                  2pt_files      +A_4pt+  gamma=Vmu + insertions Vmu     
            int id_VA = 2 * head.gammas.size() + (mu + 4) + (nu + 0) * 8;
            //                  2pt_files      +   A_4pt                +  gamma=Amu + insertions Vmu     
            int id_AV = 2 * head.gammas.size() + (mu + 0) + (nu + 4) * 8;
            fit_info.corr_id = { id_VV, id_AA, id_VA, id_AV, id_Ds_ss_2pt };//diag{ ll, ss}



            fit_info.verbosity = 0;
            id_Mmunu[mu][nu] = ncorr_new;
            printf("C_%d,%d\n", mu, nu);
            add_correlators(option, ncorr_new, conf_jack, compute_Mmunu, fit_info);
            printf(" ncorr after C_munu %d\n", ncorr_new);
            char name[NAMESIZE];
            mysprintf(name, NAMESIZE, "M_{%d,%d}", mu, nu);

            double* Mmunu = plateau_correlator_function(
                option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
                namefile_plateaux, outfile, ncorr_new - 1, name, identity, jack_file, tmp_info);
            check_correlatro_counter(6 + (nu + mu * 4) * 2);
            free(Mmunu);
            mysprintf(name, NAMESIZE, "ImM_{%d,%d}", mu, nu);
            double* Mmunu_im = plateau_correlator_function(
                option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
                namefile_plateaux, outfile, ncorr_new - 1, name, identity_im, jack_file, tmp_info);
            check_correlatro_counter(7 + (nu + mu * 4) * 2);
            free(Mmunu_im);

            fit_info.restore_default();
        }
    }

    std::vector<int> id_Y(6);
    id_Y[0] = -1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.T = head.T;
    // fit_info.n_ext_P = 1;
    // fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, Njack);
    // for (size_t j = 0; j < Njack; j++) {
    //     fit_info.ext_P[0] = head.thetas[0];
    // }


    // Y1
    id_Y[1] = ncorr_new;
    fit_info.corr_id = { id_Mmunu[1][1],id_Mmunu[2][2] };
    add_correlators(option, ncorr_new, conf_jack, compute_Y1, fit_info);

    // Y2
    id_Y[2] = ncorr_new;
    fit_info.corr_id = { id_Mmunu[0][0] };
    add_correlators(option, ncorr_new, conf_jack, compute_Y2, fit_info);

    // Y3
    id_Y[3] = ncorr_new;
    fit_info.corr_id = { id_Mmunu[3][3] };
    add_correlators(option, ncorr_new, conf_jack, compute_Y3, fit_info);

    // Y4
    id_Y[4] = ncorr_new;
    fit_info.corr_id = { id_Mmunu[0][3],id_Mmunu[3][0] };
    add_correlators(option, ncorr_new, conf_jack, compute_Y4, fit_info);

    // Y5
    id_Y[5] = ncorr_new;
    fit_info.corr_id = { id_Mmunu[1][2],id_Mmunu[2][1] };
    add_correlators(option, ncorr_new, conf_jack, compute_Y5, fit_info);

    for (int i = 1;i < 6;i++) {
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "Y_{%d}", i);
        double* Y = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Y[i], name, identity, jack_file, tmp_info);
        check_correlatro_counter(38 + (i - 1) * 2);
        free(Y);
        mysprintf(name, NAMESIZE, "ImY_{%d}", i);
        Y = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Y[i], name, identity_im, jack_file, tmp_info);
        check_correlatro_counter(39 + (i - 1) * 2);
        free(Y);
    }

    // Z0 Z1 Z2
    std::vector<int> id_Z = { ncorr_new, ncorr_new + 1, ncorr_new + 2 };
    fit_info.N = 3;
    fit_info.Njack = Njack;
    fit_info.T = head.T;
    fit_info.corr_id = id_Y;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P); // pass only the pointer
    fit_info.ext_P[0] = M_Ds;

    add_correlators(option, ncorr_new, conf_jack, compute_Z_factors, fit_info);

    fit_info.restore_default();

    for (int i = 0;i < 3;i++) {
        char name[NAMESIZE];
        mysprintf(name, NAMESIZE, "Z_{%d}", i);
        double* Z = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Z[i], name, identity, jack_file, tmp_info);
        check_correlatro_counter(48 + i * 2);
        free(Z);
        mysprintf(name, NAMESIZE, "ImZ_{%d}", i);
        Z = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Z[i], name, identity_im, jack_file, tmp_info);
        check_correlatro_counter(49 + i * 2);
        free(Z);
    }

    ///////////// structure for fits alpha
    std::vector<double> alphas = { 2.0,0,-1.99 };
    // std::vector<double> alphas = { 2.0 };
    data_all jackall;
    jackall.resampling = argv[6];
    jackall.ens = alphas.size();
    jackall.en = new data_single[jackall.ens];
    for (int i = 0;i < jackall.ens;i++) {
        jackall.en[i].header = head;
        jackall.en[i].Nobs = 1;
        jackall.en[i].Njack = head.Njack;
        jackall.en[i].jack = (double**)malloc(sizeof(double*) * jackall.en[i].Nobs);
    }
    ///////////// structure for fits alpha
    // std::vector<double> sigmas = { 0.12 ,   0.06 , 0.02, 0.05, 0.02 };
    // std::vector<double> sigmas = { 0.120, 0.100, 0.075, 0.050, 0.037, 0.025, 0.017, 0.010, 0.005, 0.0035, 0.002, 0.0005 };
    std::vector<double> sigmas = { 0.005,  0.010,  0.020,  0.030,  0.050,  0.070,  0.090,  0.110,  0.130,  0.150,  0.170,
      0.190,  0.210,  0.230,  0.250,  0.270,  0.290,  0.310,  0.330,  0.350,  0.400,  0.500 };
    // std::vector<double> sigmas = { 0.0005 };

    data_all jackall_sigma;
    jackall_sigma.resampling = argv[6];
    jackall_sigma.ens = sigmas.size();
    jackall_sigma.en = new data_single[jackall_sigma.ens];
    for (int i = 0;i < jackall_sigma.ens;i++) {
        jackall_sigma.en[i].header = head;
        jackall_sigma.en[i].Nobs = 1;
        jackall_sigma.en[i].Njack = head.Njack;
        jackall_sigma.en[i].jack = (double**)malloc(sizeof(double*) * jackall_sigma.en[i].Nobs);

    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// HLT
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Z0
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    HLT_type_input HLT_info;
    HLT_info.tmin = 1;
    HLT_info.T = head.T;
    HLT_info.type_b = HLT_EXP_b;
    HLT_info.prec = 50 * 3.33;
    HLT_info.integration_deg_limit = 1e+6;//1e+3;
    HLT_info.integration_eval_limit = 1e+9;//1e+6;
    HLT_info.integration_depth_limit = 1e+9;//1e+6;
    double omega = head.thetas[0] * M_PI / (head.L * M_Ds[Njack - 1]);
    double sigma1, dsigma1, E0_HLT;
    line_read_param(option, "sigma1", sigma1, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "E0_HLT", E0_HLT, dsigma1, myseed, namefile_plateaux);
    int dtmax;
    line_read_param(option, "HLT_tmax", HLT_info.tmax, dtmax, dtmax, namefile_plateaux);

    std::vector<double>  theta_p(3);
    // theta_p[0] = 5.739387e-01;//(1 - omega) * M_Ds[Njack - 1];
    // theta_p[1] = 9.598966e-02;//sigma1 * M_Ds[Njack - 1];
    // theta_p[2] = 2.254511e-02;//pow(omega, 3);
    // HLT_info.E0 = 3.9995692470e-02;//E0_HLT * M_Ds[Njack - 1];

    theta_p[0] = (1 - omega) * M_Ds[Njack - 1];// omega0max
    theta_p[2] = pow(omega, 3); // const
    HLT_info.E0 = E0_HLT * M_Ds[Njack - 1];

    printf("T=%d\n", HLT_info.T);
    printf("omega0=%g\n", theta_p[0]);
    printf("sigma=%g\n", theta_p[1]);
    printf("c_0=%g\n", theta_p[2]);
    printf("E_0=%g\n", HLT_info.E0);
    // HLT_type HLT_space(HLT_info);


    fit_type_HLT fit_info_HLT;
    fit_info_HLT.Njack = Njack;
    fit_info_HLT.corr_id = { id_Z[0] };
    fit_info_HLT.outfile_kernel = outfile_HLT_kernel;
    fit_info_HLT.outfile_AoverB = outfile_HLT_AoverB;
    fit_info_HLT.outfile = outfile;
    fit_info_HLT.maxE_check_reconstuct = 10;
    fit_info_HLT.stepsE_check_reconstuct = 100;
    fit_info_HLT.lambda_start = pow(2, 25);//33554432;
    fit_info_HLT.nsame = 4;
    fit_info_HLT.plateau_strategy = HLT_plateau_strategy_left;
    fit_info_HLT.strict_checks = false;
    // fit_info_HLT.nlambda_max = 10;
    line_read_param(option, "HLT_nlambda_max", fit_info_HLT.nlambda_max, dtmax, dtmax, namefile_plateaux);

    fit_info_HLT.reduce_lambda = 0.5;
    fit_info_HLT.diag_cov = false;

    double Max = 0;
    for (int t = HLT_info.tmin;t < HLT_info.tmax;t++) {
        theta_p[1] = sigmas[0] * M_Ds[Njack - 1];
        double a = alphas[0] + t + 1 / theta_p[1];
        // printf("a=%g\n", a);
        double precision = pow(10, -50);
        double b = log(precision * a / theta_p[2]);
        double max = (theta_p[0] / theta_p[1] - b) / a;
        if (max > Max) Max = max;
        // printf("maxE %d =%g\n", t, max);

    }
    printf("maxE all t =%g\n", Max);
    HLT_info.integration_maxE = 100;//Max*10;

    std::vector<HLT_type*> HLT_space;
    for (int ai = 0; ai < alphas.size();ai++) {
        HLT_info.alpha = alphas[ai];

        HLT_space.emplace_back(new HLT_type(HLT_info));
    }

    for (int si = 0;si < sigmas.size(); si++) {
        theta_p[1] = sigmas[si] * M_Ds[Njack - 1];
        for (int ai = 0; ai < alphas.size();ai++) {
            // HLT_info.alpha = alphas[ai];
            // HLT_type HLT_space(HLT_info);
            wrapper_smearing Delta(c_theta_s_HLT, theta_p, HLT_space[ai]);

            char nameHLT[NAMESIZE];
            mysprintf(nameHLT, NAMESIZE, "%s/out/%s_HLT_f_Z0-sig%f-alpha%2.2f", option[3], option[6], sigmas[si], alphas[ai]);
            if (std::filesystem::exists(nameHLT)) {
                HLT_space[ai]->load_f_EXP_b(nameHLT);
            }
            else {
                HLT_space[ai]->compute_f_EXP_b(Delta);
                HLT_space[ai]->store_f_EXP_b(nameHLT);
            }

            char namefit[NAMESIZE];
            mysprintf(namefit, NAMESIZE, "HLT_Z0-sig%f-alpha%2.2f", sigmas[si], alphas[ai]);
            fit_result tmp = HLT_space[ai]->HLT_of_corr(option, conf_jack, namefile_plateaux, namefit, Delta, jack_file, fit_info_HLT);
            jackall.en[ai].jack[0] = tmp.P[0];
            // check_correlatro_counter(54 + ai + (alphas.size() + 1) * si);
        }
        ///// fit in alpha
        fit_info.restore_default();
        fit_info.Npar = 1;
        fit_info.Nvar = 1;
        fit_info.Njack = head.Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(alphas.size());
        for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e : fit_info.myen) {
            for (int e = 0; e < fit_info.myen.size(); e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][count][j] = alphas[e];
                }
                count++;
            }
        }
        fit_info.corr_id = { 0 };
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.covariancey = false;
        fit_info.verbosity = 0;
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z0-sig%f-alpha", option[6], sigmas[si]);
        char** temp_argv = malloc_2<char>(5, NAMESIZE);
        mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
        mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

        fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);

        fit_info.band_range = { -2,2.1 };
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "alpha", fit_inter_alpha, fit_inter_alpha, 0, fit_info.myen.size() - 1, 0.5);
        // write_jack(fit_inter_alpha.P[0], Njack, jack_file);
        // check_correlatro_counter(57 + (alphas.size() + 1) * si);
        jackall_sigma.en[si].jack[0] = fit_inter_alpha.P[0];
    }
    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    int offset = 6;
    fit_info.myen = std::vector<int>(sigmas.size() - offset);
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n + offset;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    int count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e + offset]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma2_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    char namefit[NAMESIZE];
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z0-sigma", option[6]);
    char** temp_argv = malloc_2<char>(5, NAMESIZE);
    mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

    fit_result fit_Z0_sigma = fit_all_data(temp_argv, jackall_sigma, lhs_identity, fit_info, namefit);

    fit_info.band_range = { 0,0.5 };
    print_fit_band(temp_argv, jackall_sigma, fit_info, fit_info, namefit, "sigma", fit_Z0_sigma, fit_Z0_sigma, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z0_sigma.P[0], Njack, jack_file);
    check_correlatro_counter(54);




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Z1
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_all jackall_sigma_Z1;
    jackall_sigma_Z1.resampling = argv[6];
    jackall_sigma_Z1.ens = sigmas.size();
    jackall_sigma_Z1.en = new data_single[jackall_sigma_Z1.ens];
    for (int i = 0;i < jackall_sigma_Z1.ens;i++) {
        jackall_sigma_Z1.en[i].header = head;
        jackall_sigma_Z1.en[i].Nobs = 1;
        jackall_sigma_Z1.en[i].Njack = head.Njack;
        jackall_sigma_Z1.en[i].jack = (double**)malloc(sizeof(double*) * jackall_sigma_Z1.en[i].Nobs);

    }
    if (strcmp(option[6], "cC211.06.80_th9.5_t65_51.dat") == 0 || strcmp(option[6], "cC211.06.80_th9_t65_51.dat") == 0 ||
        strcmp(option[6], "cE211.044.112_th9.5_t91_72.dat") == 0 || strcmp(option[6], "cE211.044.112_th9_t91_72.dat") == 0) {
        fit_info_HLT.lambda_start = pow(2, 23);
    }
    theta_p.resize(4);
    theta_p[0] = (1 - omega) * M_Ds[Njack - 1];// omega_0^max
    // Poly
    theta_p[2] = pow(omega, 2) * (1 - omega);//  const
    theta_p[3] = -pow(omega, 2) / M_Ds[Njack - 1]; // \omega_0
    HLT_info.E0 = E0_HLT * M_Ds[Njack - 1];// 3.9995692470e-02;//E0_HLT * M_Ds[Njack - 1];
    printf("T=%d\n", HLT_info.T);
    printf("omega0=%g\n", theta_p[0]);
    printf("c_0=%g\n", theta_p[2]);
    printf("c_1=%g\n", theta_p[3]);
    printf("E_0=%g\n", HLT_info.E0);

    fit_info_HLT.corr_id = { id_Z[1] };
    for (int si = 0;si < sigmas.size(); si++) {
        theta_p[1] = sigmas[si] * M_Ds[Njack - 1];
        for (int ai = 0; ai < alphas.size();ai++) {
            // HLT_info.alpha = alphas[ai];
            // HLT_type HLT_space(HLT_info);
            wrapper_smearing Delta(c1_theta_s_HLT, theta_p, HLT_space[ai]);
            char nameHLT[NAMESIZE];
            mysprintf(nameHLT, NAMESIZE, "%s/out/%s_HLT_f_Z1-sig%f-alpha%2.2f", option[3], option[6], sigmas[si], alphas[ai]);
            if (std::filesystem::exists(nameHLT)) {
                HLT_space[ai]->load_f_EXP_b(nameHLT);
            }
            else {
                HLT_space[ai]->compute_f_EXP_b(Delta);
                HLT_space[ai]->store_f_EXP_b(nameHLT);
            }


            char namefit[NAMESIZE];
            mysprintf(namefit, NAMESIZE, "HLT_Z1-sig%f-alpha%2.2f", sigmas[si], alphas[ai]);
            fit_result tmp = HLT_space[ai]->HLT_of_corr(option, conf_jack, namefile_plateaux, namefit, Delta, jack_file, fit_info_HLT);
            jackall.en[ai].jack[0] = tmp.P[0];
            // check_correlatro_counter(75 + ai + (alphas.size() + 1) * si);
            // if (si==0 && ai==1) exit(1);
        }

        fit_info.restore_default();
        fit_info.Npar = 1;
        fit_info.Nvar = 1;
        fit_info.Njack = head.Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(alphas.size());
        for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e : fit_info.myen) {
            for (int e = 0; e < fit_info.myen.size(); e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][count][j] = alphas[e]; // mus1
                }
                count++;
            }
        }
        fit_info.corr_id = { 0 };
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.covariancey = false;
        fit_info.verbosity = 0;
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z1-sig%f-alpha", option[6], sigmas[si]);
        char** temp_argv = malloc_2<char>(5, NAMESIZE);
        mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
        mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

        fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);

        fit_info.band_range = { -2,2.1 };
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "alpha", fit_inter_alpha, fit_inter_alpha, 0, fit_info.myen.size() - 1, 0.5);
        // write_jack(fit_inter_alpha.P[0], Njack, jack_file);
        // check_correlatro_counter(78 + (alphas.size() + 1) * si);

        jackall_sigma_Z1.en[si].jack[0] = fit_inter_alpha.P[0];

    }

    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    fit_info.myen = std::vector<int>(sigmas.size());
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma2_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z1-sigma", option[6]);

    fit_result fit_Z1_sigma = fit_all_data(temp_argv, jackall_sigma_Z1, lhs_identity, fit_info, namefit);

    fit_info.band_range = { 0,0.5 };
    print_fit_band(temp_argv, jackall_sigma_Z1, fit_info, fit_info, namefit, "sigma", fit_Z1_sigma, fit_Z1_sigma, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z1_sigma.P[0], Njack, jack_file);
    check_correlatro_counter(55);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Z2
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_all jackall_sigma_Z2;
    jackall_sigma_Z2.resampling = argv[6];
    jackall_sigma_Z2.ens = sigmas.size();
    jackall_sigma_Z2.en = new data_single[jackall_sigma_Z2.ens];
    for (int i = 0;i < jackall_sigma_Z2.ens;i++) {
        jackall_sigma_Z2.en[i].header = head;
        jackall_sigma_Z2.en[i].Nobs = 1;
        jackall_sigma_Z2.en[i].Njack = head.Njack;
        jackall_sigma_Z2.en[i].jack = (double**)malloc(sizeof(double*) * jackall_sigma_Z2.en[i].Nobs);

    }

    fit_info_HLT.lambda_start = pow(2, 23);

    theta_p.resize(5);
    theta_p[0] = (1 - omega) * M_Ds[Njack - 1];// omega_0^max
    // Poly
    theta_p[2] = omega * (1 - omega) * (1 - omega);//  const= |omega| * omega0max^2
    theta_p[3] = -2.0 * omega * (1 - omega) / M_Ds[Njack - 1]; // -2 |omega| * omega0max    omega0
    theta_p[4] = omega / (M_Ds[Njack - 1] * M_Ds[Njack - 1]); //  |omega|  omega0^2
    printf("T=%d\n", HLT_info.T);
    printf("omega0=%g\n", theta_p[0]);
    printf("c_0=%g\n", theta_p[2]);
    printf("c_1=%g\n", theta_p[3]);
    printf("c_2=%g\n", theta_p[4]);
    printf("E_0=%g\n", HLT_info.E0);

    fit_info_HLT.corr_id = { id_Z[2] };
    for (int si = 0;si < sigmas.size(); si++) {
        theta_p[1] = sigmas[si] * M_Ds[Njack - 1];
        for (int ai = 0; ai < alphas.size();ai++) {

            wrapper_smearing Delta(c2_theta_s_HLT, theta_p, HLT_space[ai]);
            char nameHLT[NAMESIZE];
            mysprintf(nameHLT, NAMESIZE, "%s/out/%s_HLT_f_Z2-sig%f-alpha%2.2f", option[3], option[6], sigmas[si], alphas[ai]);
            if (std::filesystem::exists(nameHLT)) {
                HLT_space[ai]->load_f_EXP_b(nameHLT);
            }
            else {
                HLT_space[ai]->compute_f_EXP_b(Delta);
                HLT_space[ai]->store_f_EXP_b(nameHLT);
            }

            char namefit[NAMESIZE];
            mysprintf(namefit, NAMESIZE, "HLT_Z2-sig%f-alpha%2.2f", sigmas[si], alphas[ai]);
            fit_result tmp = HLT_space[ai]->HLT_of_corr(option, conf_jack, namefile_plateaux, namefit, Delta, jack_file, fit_info_HLT);
            jackall.en[ai].jack[0] = tmp.P[0];
            // check_correlatro_counter(75 + ai + (alphas.size() + 1) * si);
            // if (si==0 && ai==1) exit(1);
        }

        fit_info.restore_default();
        fit_info.Npar = 1;
        fit_info.Nvar = 1;
        fit_info.Njack = head.Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(alphas.size());
        for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e : fit_info.myen) {
            for (int e = 0; e < fit_info.myen.size(); e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][count][j] = alphas[e]; // mus1
                }
                count++;
            }
        }
        fit_info.corr_id = { 0 };
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.covariancey = false;
        fit_info.verbosity = 0;
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z2-sig%f-alpha", option[6], sigmas[si]);
        char** temp_argv = malloc_2<char>(5, NAMESIZE);
        mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
        mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

        fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);

        fit_info.band_range = { -2,2.1 };
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "alpha", fit_inter_alpha, fit_inter_alpha, 0, fit_info.myen.size() - 1, 0.5);
        // write_jack(fit_inter_alpha.P[0], Njack, jack_file);
        // check_correlatro_counter(78 + (alphas.size() + 1) * si);

        jackall_sigma_Z2.en[si].jack[0] = fit_inter_alpha.P[0];

    }

    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    fit_info.myen = std::vector<int>(sigmas.size());
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma2_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z2-sigma", option[6]);

    fit_result fit_Z2_sigma = fit_all_data(temp_argv, jackall_sigma_Z2, lhs_identity, fit_info, namefit);

    fit_info.band_range = { 0,0.5 };
    print_fit_band(temp_argv, jackall_sigma_Z2, fit_info, fit_info, namefit, "sigma", fit_Z2_sigma, fit_Z2_sigma, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z2_sigma.P[0], Njack, jack_file);
    check_correlatro_counter(56);

    ////////////////////////////////////////////////////////////////////////////////////////////
    // sum the Z
    /// 

    double* dGammadq = (double*)malloc(sizeof(double) * Njack);
    double* q2 = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        dGammadq[j] = fit_Z0_sigma.P[0][j] + fit_Z1_sigma.P[0][j] + fit_Z2_sigma.P[0][j];
        if (j == Njack - 1) printf(" sumZ=%g\n", dGammadq[j]);
        dGammadq[j] *= (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) / (2 * M_PI);
        if (j == Njack - 1) printf(" dgammaZ=%g\n", dGammadq[j]);
        if (j == Njack - 1) printf(" a=%g\n", a_fm_jack[j]);
        if (j == Njack - 1) printf(" a GeV=%g\n", a_GeV_jack[j]);
        if (j == Njack - 1) printf(" MDS GeV=%g\n", M_Ds[j] / a_GeV_jack[j]);
        q2[j] = head.thetas[0] * M_PI / (head.L * a_GeV_jack[j]);
        q2[j] *= q2[j];
    }
    fprintf(outfile, "\n\n#dGamma/dq_GEV3 summign Z\n");
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z0_sigma.P[0][Njack - 1], myres->comp_error(fit_Z0_sigma.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z1_sigma.P[0][Njack - 1], myres->comp_error(fit_Z1_sigma.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z2_sigma.P[0][Njack - 1], myres->comp_error(fit_Z2_sigma.P[0]));
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", "dGammadq_GEV3", 0, 2, 0.0, 0.0);
    fprintf(outfile, "%.15g    %.15g     %.15g    %.15g\n", dGammadq[Njack - 1], myres->comp_error(dGammadq), q2[Njack - 1], myres->comp_error(q2));
    printf("dGammadq_GEV3 = %.15g    %.15g    \n", dGammadq[Njack - 1], myres->comp_error(dGammadq));


    write_jack(dGammadq, Njack, jack_file);
    check_correlatro_counter(57);
    write_jack(a_fm_jack, Njack, jack_file);
    check_correlatro_counter(58);

    {
        double* M_Z = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Z[0], "M_{Z_0}", M_eff_log, jack_file);
        check_correlatro_counter(59);
        free(M_Z);

        M_Z = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Z[1], "M_{Z_1}", M_eff_log, jack_file);
        check_correlatro_counter(60);
        free(M_Z);

        M_Z = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Z[2], "M_{Z_2}", M_eff_log, jack_file);
        check_correlatro_counter(61);
        free(M_Z);

    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Z0 by part  (omega-omega0) theta(omega-omega0)
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_all jackall_sigma_Z0_part;
    jackall_sigma_Z0_part.resampling = argv[6];
    jackall_sigma_Z0_part.ens = sigmas.size();
    jackall_sigma_Z0_part.en = new data_single[jackall_sigma_Z0_part.ens];
    for (int i = 0;i < jackall_sigma_Z0_part.ens;i++) {
        jackall_sigma_Z0_part.en[i].header = head;
        jackall_sigma_Z0_part.en[i].Nobs = 1;
        jackall_sigma_Z0_part.en[i].Njack = head.Njack;
        jackall_sigma_Z0_part.en[i].jack = (double**)malloc(sizeof(double*) * jackall_sigma_Z0_part.en[i].Nobs);

    }
    HLT_info;
    HLT_info.tmin = 1;
    HLT_info.T = head.T;
    HLT_info.prec = 50 * 3.33;
    HLT_info.integration_deg_limit = 1e+6;//1e+3;
    HLT_info.integration_eval_limit = 1e+9;//1e+6;
    HLT_info.integration_depth_limit = 1e+9;//1e+6;
    line_read_param(option, "sigma1", sigma1, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "E0_HLT", E0_HLT, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "HLT_tmax", HLT_info.tmax, dtmax, dtmax, namefile_plateaux);


    theta_p.resize(4);
    theta_p[0] = (1 - omega) * M_Ds[Njack - 1];// omega_0^max
    // Poly
    theta_p[2] = +pow(omega, 3) * theta_p[0];   // constant
    theta_p[3] = -pow(omega, 3); // coeff of omega term
    HLT_info.type_b = HLT_EXP_bovert;
    HLT_info.E0 = E0_HLT * M_Ds[Njack - 1];



    fit_info_HLT.Njack = Njack;
    fit_info_HLT.corr_id = { id_Z[0] };
    fit_info_HLT.outfile_kernel = outfile_HLT_kernel;
    fit_info_HLT.outfile_AoverB = outfile_HLT_AoverB;
    fit_info_HLT.outfile = outfile;
    fit_info_HLT.maxE_check_reconstuct = 10;
    fit_info_HLT.stepsE_check_reconstuct = 100;
    fit_info_HLT.lambda_start = pow(2, 25);//3.09454500212;
    fit_info_HLT.nsame = 4;
    fit_info_HLT.strict_checks = false;
    // fit_info_HLT.nlambda_max = 10;
    line_read_param(option, "HLT_nlambda_max", fit_info_HLT.nlambda_max, dtmax, dtmax, namefile_plateaux);

    fit_info_HLT.reduce_lambda = 0.5;
    fit_info_HLT.diag_cov = false;

    HLT_info.integration_maxE = 100;//Max*10;
    std::vector<HLT_type*> HLT_space_part;
    for (int ai = 0; ai < alphas.size();ai++) {
        HLT_info.alpha = alphas[ai];

        HLT_space_part.emplace_back(new HLT_type(HLT_info));
    }

    for (int si = 0;si < sigmas.size(); si++) {
        theta_p[1] = sigmas[si] * M_Ds[Njack - 1];
        for (int ai = 0; ai < alphas.size();ai++) {
            wrapper_smearing Delta(c1_theta_s_HLT, theta_p, HLT_space_part[ai]);

            char nameHLT[NAMESIZE];
            mysprintf(nameHLT, NAMESIZE, "%s/out/%s_HLT_f_part_Z0-sig%f-alpha%2.2f", option[3], option[6], sigmas[si], alphas[ai]);
            if (std::filesystem::exists(nameHLT)) {
                HLT_space_part[ai]->load_f_EXP_b(nameHLT);
            }
            else {
                HLT_space_part[ai]->compute_f_EXP_b(Delta);
                HLT_space_part[ai]->store_f_EXP_b(nameHLT);
            }

            char namefit[NAMESIZE];
            mysprintf(namefit, NAMESIZE, "HLT_part_Z0-sig%f-alpha%2.2f", sigmas[si], alphas[ai]);
            fit_result tmp = HLT_space_part[ai]->HLT_of_corr(option, conf_jack, namefile_plateaux, namefit, Delta, jack_file, fit_info_HLT);
            jackall.en[ai].jack[0] = tmp.P[0];
            // check_correlatro_counter(54 + ai + (alphas.size() + 1) * si);
        }
        ///// fit in alpha
        fit_info.restore_default();
        fit_info.Npar = 1;
        fit_info.Nvar = 1;
        fit_info.Njack = head.Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(alphas.size());
        for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e : fit_info.myen) {
            for (int e = 0; e < fit_info.myen.size(); e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][count][j] = alphas[e];
                }
                count++;
            }
        }
        fit_info.corr_id = { 0 };
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.covariancey = false;
        fit_info.verbosity = 0;
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "%s_fit_HLT_part_Z0-sig%f-alpha", option[6], sigmas[si]);
        char** temp_argv = malloc_2<char>(5, NAMESIZE);
        mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
        mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

        fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);

        fit_info.band_range = { -2,2.1 };
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "alpha", fit_inter_alpha, fit_inter_alpha, 0, fit_info.myen.size() - 1, 0.5);
        // write_jack(fit_inter_alpha.P[0], Njack, jack_file);
        // check_correlatro_counter(57 + (alphas.size() + 1) * si);
        jackall_sigma_Z0_part.en[si].jack[0] = fit_inter_alpha.P[0];
    }
    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    fit_info.myen = std::vector<int>(sigmas.size() - offset);
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n + offset;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e + offset]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma2_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    namefit[NAMESIZE];
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_part_Z0-sigma", option[6]);
    mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// outfile

    fit_result fit_part_Z0_sigma = fit_all_data(temp_argv, jackall_sigma_Z0_part, lhs_identity, fit_info, namefit);

    fit_info.band_range = { 0,0.5 };
    print_fit_band(temp_argv, jackall_sigma_Z0_part, fit_info, fit_info, namefit, "sigma", fit_part_Z0_sigma, fit_part_Z0_sigma, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_part_Z0_sigma.P[0], Njack, jack_file);
    check_correlatro_counter(62);


    ////////////////////////////////////////////////////////////////////////////////////////////
    // sum the Z
    ////////////////////////////////////////////////////////////////////////////////////////////


    for (int j = 0;j < Njack;j++) {
        dGammadq[j] = fit_part_Z0_sigma.P[0][j] + fit_Z1_sigma.P[0][j] + fit_Z2_sigma.P[0][j];
        if (j == Njack - 1) printf(" sumZ=%g\n", dGammadq[j]);
        dGammadq[j] *= (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) / (2 * M_PI);
        if (j == Njack - 1) printf(" dgammaZ=%g\n", dGammadq[j]);
        if (j == Njack - 1) printf(" a=%g\n", a_fm_jack[j]);
        if (j == Njack - 1) printf(" a GeV=%g\n", a_GeV_jack[j]);
        if (j == Njack - 1) printf(" MDS GeV=%g\n", M_Ds[j] / a_GeV_jack[j]);
        q2[j] = head.thetas[0] * M_PI / (head.L * a_GeV_jack[j]);
        q2[j] *= q2[j];
    }
    fprintf(outfile, "\n\n#dGamma/dq_GEV3 summign Z, Z0 by part\n");
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_part_Z0_sigma.P[0][Njack - 1], myres->comp_error(fit_part_Z0_sigma.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z1_sigma.P[0][Njack - 1], myres->comp_error(fit_Z1_sigma.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z2_sigma.P[0][Njack - 1], myres->comp_error(fit_Z2_sigma.P[0]));
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", "dGammadq_part_GEV3", 0, 2, 0.0, 0.0);
    fprintf(outfile, "%.15g    %.15g     %.15g    %.15g\n", dGammadq[Njack - 1], myres->comp_error(dGammadq), q2[Njack - 1], myres->comp_error(q2));
    printf("dGammadq_part_GEV3 = %.15g    %.15g    \n", dGammadq[Njack - 1], myres->comp_error(dGammadq));


    write_jack(dGammadq, Njack, jack_file);
    check_correlatro_counter(63);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Z0 bounding plus
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_all jackall_sigma_Z0_plus;
    jackall_sigma_Z0_plus.resampling = argv[6];
    jackall_sigma_Z0_plus.ens = sigmas.size();
    jackall_sigma_Z0_plus.en = new data_single[jackall_sigma_Z0_plus.ens];
    for (int i = 0;i < jackall_sigma_Z0_plus.ens;i++) {
        jackall_sigma_Z0_plus.en[i].header = head;
        jackall_sigma_Z0_plus.en[i].Nobs = 1;
        jackall_sigma_Z0_plus.en[i].Njack = head.Njack;
        jackall_sigma_Z0_plus.en[i].jack = (double**)malloc(sizeof(double*) * jackall_sigma_Z0_plus.en[i].Nobs);

    }
    HLT_info.tmin = 1;
    HLT_info.T = head.T;
    HLT_info.prec = 50 * 3.33;
    HLT_info.integration_deg_limit = 1e+6;//1e+3;
    HLT_info.integration_eval_limit = 1e+9;//1e+6;
    HLT_info.integration_depth_limit = 1e+9;//1e+6;
    line_read_param(option, "sigma1", sigma1, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "E0_HLT", E0_HLT, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "HLT_tmax", HLT_info.tmax, dtmax, dtmax, namefile_plateaux);


    theta_p.resize(3);
    theta_p[0] = (1 - omega) * M_Ds[Njack - 1];// omega0max
    theta_p[2] = pow(omega, 3); // const

    HLT_info.type_b = HLT_EXP_bovert;
    HLT_info.E0 = E0_HLT * M_Ds[Njack - 1];



    fit_info_HLT.Njack = Njack;
    fit_info_HLT.corr_id = { id_Z[0] };
    fit_info_HLT.outfile_kernel = outfile_HLT_kernel;
    fit_info_HLT.outfile_AoverB = outfile_HLT_AoverB;
    fit_info_HLT.outfile = outfile;
    fit_info_HLT.maxE_check_reconstuct = 10;
    fit_info_HLT.stepsE_check_reconstuct = 100;
    fit_info_HLT.lambda_start = pow(2, 25);//3.09454500212;
    fit_info_HLT.nsame = 4;
    fit_info_HLT.strict_checks = false;
    // fit_info_HLT.nlambda_max = 10;
    line_read_param(option, "HLT_nlambda_max", fit_info_HLT.nlambda_max, dtmax, dtmax, namefile_plateaux);

    // fit_info_HLT.reduce_lambda = 0.5;
    // fit_info_HLT.diag_cov = false;

    HLT_info.integration_maxE = 200;//Max*10;
    std::vector<HLT_type*> HLT_space_Z0_plus;

    alphas = { 0 };
    for (int ai = 0; ai < alphas.size();ai++) {
        HLT_info.alpha = alphas[ai];

        HLT_space_Z0_plus.emplace_back(new HLT_type(HLT_info));
    }

    for (int si = 0;si < sigmas.size(); si++) {
        theta_p[1] = sigmas[si] * M_Ds[Njack - 1];
        for (int ai = 0; ai < alphas.size();ai++) {
            wrapper_smearing Delta(c_thetap_s_HLT, theta_p, HLT_space_Z0_plus[ai]);

            char nameHLT[NAMESIZE];
            mysprintf(nameHLT, NAMESIZE, "%s/out/%s_HLT_f_plus_Z0-sig%f-alpha%2.2f", option[3], option[6], sigmas[si], alphas[ai]);
            if (std::filesystem::exists(nameHLT)) {
                double t1 = timestamp();
                HLT_space_Z0_plus[ai]->load_f_EXP_b(nameHLT);
                printf("time to read f:  %g s\n", timestamp() - t1);
            }
            else {
                double t1 = timestamp();
                HLT_space_Z0_plus[ai]->compute_f_EXP_b(Delta);
                HLT_space_Z0_plus[ai]->store_f_EXP_b(nameHLT);
                printf("time to compute f:  %g s\n", timestamp() - t1);
            }

            char namefit[NAMESIZE];
            mysprintf(namefit, NAMESIZE, "HLT_plus_Z0-sig%f-alpha%2.2f", sigmas[si], alphas[ai]);
            fit_result tmp = HLT_space_Z0_plus[ai]->HLT_of_corr(option, conf_jack, namefile_plateaux, namefit, Delta, jack_file, fit_info_HLT);
            jackall.en[ai].jack[0] = tmp.P[0];
            // check_correlatro_counter(54 + ai + (alphas.size() + 1) * si);
        }
        ///// fit in alpha
        fit_info.restore_default();
        fit_info.Npar = 1;
        fit_info.Nvar = 1;
        fit_info.Njack = head.Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(alphas.size());
        for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e : fit_info.myen) {
            for (int e = 0; e < fit_info.myen.size(); e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][count][j] = alphas[e];
                }
                count++;
            }
        }
        fit_info.corr_id = { 0 };
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.covariancey = false;
        fit_info.verbosity = 0;
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "%s_fit_HLT_plus_Z0-sig%f-alpha", option[6], sigmas[si]);
        char** temp_argv = malloc_2<char>(5, NAMESIZE);
        mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
        mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

        fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);

        fit_info.band_range = { -2,2.1 };
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "alpha", fit_inter_alpha, fit_inter_alpha, 0, fit_info.myen.size() - 1, 0.5);
        // write_jack(fit_inter_alpha.P[0], Njack, jack_file);
        // check_correlatro_counter(57 + (alphas.size() + 1) * si);
        jackall_sigma_Z0_plus.en[si].jack[0] = fit_inter_alpha.P[0];
    }
    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    fit_info.myen = std::vector<int>(sigmas.size() - offset);
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n + offset;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e + offset]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma2_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    namefit[NAMESIZE];
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_plus_Z0-sigma", option[6]);
    mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// outfile

    fit_result fit_plus_Z0_sigma = fit_all_data(temp_argv, jackall_sigma_Z0_plus, lhs_identity, fit_info, namefit);

    fit_info.band_range = { 0,0.5 };
    print_fit_band(temp_argv, jackall_sigma_Z0_plus, fit_info, fit_info, namefit, "sigma", fit_plus_Z0_sigma, fit_plus_Z0_sigma, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_plus_Z0_sigma.P[0], Njack, jack_file);
    check_correlatro_counter(64);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Z0 bounding minus
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    data_all jackall_sigma_Z0_minus;
    jackall_sigma_Z0_minus.resampling = argv[6];
    jackall_sigma_Z0_minus.ens = sigmas.size();
    jackall_sigma_Z0_minus.en = new data_single[jackall_sigma_Z0_minus.ens];
    for (int i = 0;i < jackall_sigma_Z0_minus.ens;i++) {
        jackall_sigma_Z0_minus.en[i].header = head;
        jackall_sigma_Z0_minus.en[i].Nobs = 1;
        jackall_sigma_Z0_minus.en[i].Njack = head.Njack;
        jackall_sigma_Z0_minus.en[i].jack = (double**)malloc(sizeof(double*) * jackall_sigma_Z0_minus.en[i].Nobs);

    }
    HLT_info.tmin = 1;
    HLT_info.T = head.T;
    HLT_info.prec = 50 * 3.33;
    HLT_info.integration_deg_limit = 1e+6;//1e+3;
    HLT_info.integration_eval_limit = 1e+9;//1e+6;
    HLT_info.integration_depth_limit = 1e+9;//1e+6;
    line_read_param(option, "sigma1", sigma1, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "E0_HLT", E0_HLT, dsigma1, myseed, namefile_plateaux);
    line_read_param(option, "HLT_tmax", HLT_info.tmax, dtmax, dtmax, namefile_plateaux);


    theta_p.resize(3);
    theta_p[0] = (1 - omega) * M_Ds[Njack - 1];// omega0max
    theta_p[2] = pow(omega, 3); // const

    HLT_info.type_b = HLT_EXP_bovert;
    HLT_info.E0 = E0_HLT * M_Ds[Njack - 1];



    fit_info_HLT.Njack = Njack;
    fit_info_HLT.corr_id = { id_Z[0] };
    fit_info_HLT.outfile_kernel = outfile_HLT_kernel;
    fit_info_HLT.outfile_AoverB = outfile_HLT_AoverB;
    fit_info_HLT.outfile = outfile;
    fit_info_HLT.maxE_check_reconstuct = 10;
    fit_info_HLT.stepsE_check_reconstuct = 100;
    fit_info_HLT.lambda_start = pow(2, 25);//3.09454500212;
    fit_info_HLT.nsame = 4;
    fit_info_HLT.strict_checks = false;
    // fit_info_HLT.nlambda_max = 10;
    line_read_param(option, "HLT_nlambda_max", fit_info_HLT.nlambda_max, dtmax, dtmax, namefile_plateaux);

    // fit_info_HLT.reduce_lambda = 0.5;
    // fit_info_HLT.diag_cov = false;

    HLT_info.integration_maxE = 200;//Max*10;
    std::vector<HLT_type*> HLT_space_Z0_minus;

    alphas = { 0 };
    for (int ai = 0; ai < alphas.size();ai++) {
        HLT_info.alpha = alphas[ai];

        HLT_space_Z0_minus.emplace_back(new HLT_type(HLT_info));
    }

    for (int si = 0;si < sigmas.size(); si++) {
        theta_p[1] = sigmas[si] * M_Ds[Njack - 1];
        for (int ai = 0; ai < alphas.size();ai++) {
            wrapper_smearing Delta(c_thetam_s_HLT, theta_p, HLT_space_Z0_minus[ai]);

            char nameHLT[NAMESIZE];
            mysprintf(nameHLT, NAMESIZE, "%s/out/%s_HLT_f_minus_Z0-sig%f-alpha%2.2f", option[3], option[6], sigmas[si], alphas[ai]);
            if (std::filesystem::exists(nameHLT)) {
                double t1 = timestamp();
                HLT_space_Z0_minus[ai]->load_f_EXP_b(nameHLT);
                printf("time to read f:  %g s\n", timestamp() - t1);

            }
            else {
                double t1 = timestamp();
                HLT_space_Z0_minus[ai]->compute_f_EXP_b(Delta);
                HLT_space_Z0_minus[ai]->store_f_EXP_b(nameHLT);
                printf("time to compute f:  %g s\n", timestamp() - t1);
            }

            char namefit[NAMESIZE];
            mysprintf(namefit, NAMESIZE, "HLT_minus_Z0-sig%f-alpha%2.2f", sigmas[si], alphas[ai]);
            fit_result tmp = HLT_space_Z0_minus[ai]->HLT_of_corr(option, conf_jack, namefile_plateaux, namefit, Delta, jack_file, fit_info_HLT);
            jackall.en[ai].jack[0] = tmp.P[0];
            // check_correlatro_counter(54 + ai + (alphas.size() + 1) * si);
        }
        ///// fit in alpha
        fit_info.restore_default();
        fit_info.Npar = 1;
        fit_info.Nvar = 1;
        fit_info.Njack = head.Njack;
        fit_info.N = 1;
        fit_info.myen = std::vector<int>(alphas.size());
        for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
        fit_info.entot = fit_info.myen.size() * fit_info.N;
        fit_info.malloc_x();
        int count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            // for (int e : fit_info.myen) {
            for (int e = 0; e < fit_info.myen.size(); e++) {
                for (int j = 0;j < fit_info.Njack;j++) {
                    fit_info.x[0][count][j] = alphas[e];
                }
                count++;
            }
        }
        fit_info.corr_id = { 0 };
        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.covariancey = false;
        fit_info.verbosity = 0;
        char namefit[NAMESIZE];
        mysprintf(namefit, NAMESIZE, "%s_fit_HLT_minus_Z0-sig%f-alpha", option[6], sigmas[si]);
        char** temp_argv = malloc_2<char>(5, NAMESIZE);
        mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
        mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

        fit_result fit_inter_alpha = fit_all_data(temp_argv, jackall, lhs_identity, fit_info, namefit);

        fit_info.band_range = { -2,2.1 };
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "alpha", fit_inter_alpha, fit_inter_alpha, 0, fit_info.myen.size() - 1, 0.5);
        // write_jack(fit_inter_alpha.P[0], Njack, jack_file);
        // check_correlatro_counter(57 + (alphas.size() + 1) * si);
        jackall_sigma_Z0_minus.en[si].jack[0] = fit_inter_alpha.P[0];
    }
    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 2;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    fit_info.myen = std::vector<int>(sigmas.size() - offset);
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n + offset;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e + offset]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma2_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    namefit[NAMESIZE];
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_minus_Z0-sigma", option[6]);
    mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// outfile

    fit_result fit_minus_Z0_sigma = fit_all_data(temp_argv, jackall_sigma_Z0_minus, lhs_identity, fit_info, namefit);

    fit_info.band_range = { 0,0.5 };
    print_fit_band(temp_argv, jackall_sigma_Z0_minus, fit_info, fit_info, namefit, "sigma", fit_minus_Z0_sigma, fit_minus_Z0_sigma, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_minus_Z0_sigma.P[0], Njack, jack_file);
    check_correlatro_counter(65);

    //////////////////////////////////////////////////////////////
    // fit in sigma^2 + sigma^4
    //////////////////////////////////////////////////////////////
    // fit in sigma
    fit_info.restore_default();
    fit_info.Npar = 3;
    fit_info.Nvar = 1;
    fit_info.Njack = head.Njack;
    fit_info.N = 1;
    offset = 6;
    int esclude_last = 1;
    fit_info.myen = std::vector<int>(sigmas.size() - offset - esclude_last);
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n + offset;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e + offset]; // mus1
            }
            count++;
        }
    }
    fit_info.corr_id = { 0 };
    fit_info.function = sigma4_fit;//constant_fit;
    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling
    fit_info.band_range = { 0,0.5 };

    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z0-sigma4", option[6]);
    fit_result fit_Z0_sigma4 = fit_all_data(temp_argv, jackall_sigma, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_sigma, fit_info, fit_info, namefit, "sigma", fit_Z0_sigma4, fit_Z0_sigma4, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z0_sigma4.P[0], Njack, jack_file);
    // fit_Z0_sigma4.clear();
    check_correlatro_counter(66);

    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z1-sigma4", option[6]);
    fit_result fit_Z1_sigma4 = fit_all_data(temp_argv, jackall_sigma_Z1, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_sigma_Z1, fit_info, fit_info, namefit, "sigma", fit_Z1_sigma4, fit_Z1_sigma4, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z1_sigma4.P[0], Njack, jack_file);
    // fit_Z1_sigma4.clear();
    check_correlatro_counter(67);

    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_Z2-sigma4", option[6]);
    fit_result fit_Z2_sigma4 = fit_all_data(temp_argv, jackall_sigma_Z2, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_sigma_Z2, fit_info, fit_info, namefit, "sigma", fit_Z2_sigma4, fit_Z2_sigma4, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z2_sigma4.P[0], Njack, jack_file);
    // fit_Z2_sigma4.clear();
    check_correlatro_counter(68);


    ////////////////////////////////////////////////////////////////////////////////////////////
    // sum the Z
    ////////////////////////////////////////////////////////////////////////////////////////////


    for (int j = 0;j < Njack;j++) {
        dGammadq[j] = fit_Z0_sigma4.P[0][j] + fit_Z1_sigma4.P[0][j] + fit_Z2_sigma4.P[0][j];
        if (j == Njack - 1) printf(" sumZ=%g\n", dGammadq[j]);
        dGammadq[j] *= (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) / (2 * M_PI);
        if (j == Njack - 1) printf(" dgammaZ=%g\n", dGammadq[j]);
        if (j == Njack - 1) printf(" a=%g\n", a_fm_jack[j]);
        if (j == Njack - 1) printf(" a GeV=%g\n", a_GeV_jack[j]);
        if (j == Njack - 1) printf(" MDS GeV=%g\n", M_Ds[j] / a_GeV_jack[j]);
        q2[j] = head.thetas[0] * M_PI / (head.L * a_GeV_jack[j]);
        q2[j] *= q2[j];
    }
    fprintf(outfile, "\n\n#dGamma/dq_GEV3 summign Z, sigma4 extrapolation\n");
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z0_sigma4.P[0][Njack - 1], myres->comp_error(fit_Z0_sigma4.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z1_sigma4.P[0][Njack - 1], myres->comp_error(fit_Z1_sigma4.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z2_sigma4.P[0][Njack - 1], myres->comp_error(fit_Z2_sigma4.P[0]));
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", "dGammadq_GEV3_sigma4", 0, 2, 0.0, 0.0);
    fprintf(outfile, "%.15g    %.15g     %.15g    %.15g\n", dGammadq[Njack - 1], myres->comp_error(dGammadq), q2[Njack - 1], myres->comp_error(q2));
    printf("dGammadq_part_GEV3 = %.15g    %.15g    \n", dGammadq[Njack - 1], myres->comp_error(dGammadq));


    write_jack(dGammadq, Njack, jack_file);
    check_correlatro_counter(69);

    offset = 3;
    esclude_last = 2;
    fit_info.myen = std::vector<int>(sigmas.size() - offset - esclude_last);
    for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n + offset;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        // for (int e : fit_info.myen) {
        for (int e = 0; e < fit_info.myen.size(); e++) {
            for (int j = 0;j < fit_info.Njack;j++) {
                fit_info.x[0][count][j] = sigmas[e + offset]; // mus1
            }
            count++;
        }
    }
    mysprintf(namefit, NAMESIZE, "%s_fit_HLT_part_Z0-sigma4", option[6]);
    fit_result fit_Z0_part_sigma4 = fit_all_data(temp_argv, jackall_sigma_Z0_part, lhs_identity, fit_info, namefit);
    print_fit_band(temp_argv, jackall_sigma_Z0_part, fit_info, fit_info, namefit, "sigma", fit_Z0_part_sigma4, fit_Z0_part_sigma4, 0, fit_info.myen.size() - 1, 0.01);
    write_jack(fit_Z0_part_sigma4.P[0], Njack, jack_file);
    // fit_Z0_part_sigma4.clear();
    check_correlatro_counter(70);

    for (int j = 0;j < Njack;j++) {
        dGammadq[j] = fit_Z0_part_sigma4.P[0][j] + fit_Z1_sigma4.P[0][j] + fit_Z2_sigma4.P[0][j];
        if (j == Njack - 1) printf(" sumZ=%g\n", dGammadq[j]);
        dGammadq[j] *= (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) * (M_Ds[j] / a_GeV_jack[j]) / (2 * M_PI);
        if (j == Njack - 1) printf(" dgammaZ=%g\n", dGammadq[j]);
        if (j == Njack - 1) printf(" a=%g\n", a_fm_jack[j]);
        if (j == Njack - 1) printf(" a GeV=%g\n", a_GeV_jack[j]);
        if (j == Njack - 1) printf(" MDS GeV=%g\n", M_Ds[j] / a_GeV_jack[j]);
        q2[j] = head.thetas[0] * M_PI / (head.L * a_GeV_jack[j]);
        q2[j] *= q2[j];
    }
    fprintf(outfile, "\n\n#dGamma/dq_GEV3 summign Z, Z0 by part\n");
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z0_part_sigma4.P[0][Njack - 1], myres->comp_error(fit_Z0_part_sigma4.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z1_sigma4.P[0][Njack - 1], myres->comp_error(fit_Z1_sigma4.P[0]));
    fprintf(outfile, "%d   %.15e    %.15e\t", 0, fit_Z2_sigma4.P[0][Njack - 1], myres->comp_error(fit_Z2_sigma4.P[0]));
    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", "dGammadq_GEV3_part_sigma4", 0, 2, 0.0, 0.0);
    fprintf(outfile, "%.15g    %.15g     %.15g    %.15g\n", dGammadq[Njack - 1], myres->comp_error(dGammadq), q2[Njack - 1], myres->comp_error(q2));
    printf("dGammadq_part_GEV3 = %.15g    %.15g    \n", dGammadq[Njack - 1], myres->comp_error(dGammadq));

    write_jack(dGammadq, Njack, jack_file);
    check_correlatro_counter(71);


}