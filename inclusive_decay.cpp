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
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
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

    //////////////////////////////////// confs
    double**** data = calloc_corr(confs, head.ncorr, head.T);

    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < confs; iconf++) {
        read_twopt(infile, data[iconf], head);
    }

    double**** data_bin = binning(confs, head.ncorr, head.T, data, bin);
    double**** conf_jack = myres->create(Neff, head.ncorr, head.T, data_bin);
    free_corr(Neff, head.ncorr, head.T, data_bin);
    free_corr(confs, head.ncorr, head.T, data);

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
    symmetrise_jackboot(Njack, id_K_ss_2pt, head.T, conf_jack, 1);
    symmetrise_jackboot(Njack, id_K_ls_2pt, head.T, conf_jack, 1);


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


    //////////////////////
    std::vector<std::vector<int>> id_Mmunu(4, std::vector<int>(4));
    int ncorr_new = head.ncorr;
    int TJW, TDs, myerr, myseed;
    line_read_param(option, "TJW", TJW, myerr, myseed, namefile_plateaux);
    line_read_param(option, "TDs", TDs, myerr, myseed, namefile_plateaux);
    double ZA, ZV, dZA, dZV;

    line_read_param(option, "ZA", ZA, dZA, myseed, namefile_plateaux);
    double* ZA_jack = myres->create_fake(ZA, dZA, myseed);
    line_read_param(option, "ZV", ZV, dZV, myseed, namefile_plateaux);
    double* ZV_jack = myres->create_fake(ZV, dZV, myseed);

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
            int id_VV = 4 * head.gammas.size() + 4 * head.gammas.size() + (mu + 5) + nu * head.gammas.size();
            //                   2pt_files       +  gamma=Amu + insertions Amu     
            int id_AA = 4 * head.gammas.size() + 0 + (mu + 1) + nu * head.gammas.size();
            //                  2pt_files      +A_4pt+  gamma=Vmu + insertions Vmu     
            int id_VA = 4 * head.gammas.size() + 0 + (mu + 5) + nu * head.gammas.size();
            //                  2pt_files      +   A_4pt                +  gamma=Amu + insertions Vmu     
            int id_AV = 4 * head.gammas.size() + 4 * head.gammas.size() + (mu + 1) + nu * head.gammas.size();
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
                namefile_plateaux, outfile, ncorr_new - 1, name, identity, jack_file);
            check_correlatro_counter(6 + (nu + mu * 4) * 2);
            free(Mmunu);
            mysprintf(name, NAMESIZE, "ImM_{%d,%d}", mu, nu);
            double* Mmunu_im = plateau_correlator_function(
                option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
                namefile_plateaux, outfile, ncorr_new - 1, name, identity_im, jack_file);
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
            namefile_plateaux, outfile, id_Y[i], name, identity, jack_file);
        check_correlatro_counter(38 + (i - 1) * 2);
        free(Y);
        mysprintf(name, NAMESIZE, "ImY_{%d}", i);
        Y = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Y[i], name, identity_im, jack_file);
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
            namefile_plateaux, outfile, id_Z[i], name, identity, jack_file);
        check_correlatro_counter(48 + i * 2);
        free(Z);
        mysprintf(name, NAMESIZE, "ImZ_{%d}", i);
        Z = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_Z[i], name, identity_im, jack_file);
        check_correlatro_counter(49 + i * 2);
        free(Z);
    }

    HLT_type_input HLT_info;
    HLT_info.tmax = 32;
    HLT_info.T = head.T;
    HLT_info.E0 = 0.1;
    HLT_info.type_b = HLT_EXP_b;
    HLT_info.prec = 50 * 3.33;
    HLT_type HLT_space(HLT_info);
    std::vector<double>  theta_p = { 0.5, 0.1 };
    // wrapper_smearing Delta(theta_s_HLT, theta_p, &HLT_space);
    wrapper_smearing Delta(gaussian_for_HLT, theta_p, &HLT_space);
    HLT_space.compute_f_EXP_b(Delta);

    fit_type_HLT fit_info_HLT;
    fit_info_HLT.Njack = Njack;
    fit_info_HLT.corr_id = { id_Z[0] };
    fit_info_HLT.lambdas = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9  };
    fit_info_HLT.outfile_kernel = outfile_HLT_kernel;
    fit_info_HLT.outfile_AoverB = outfile_HLT_AoverB;
    fit_info_HLT.maxE_check_reconstuct =2;
    fit_info_HLT.stepsE_check_reconstuct=100;
    double** tmp = HLT_space.HLT_of_corr(option, conf_jack, namefile_plateaux, "HLT_Z0", Delta, jack_file, fit_info_HLT);


    // eg of fit to correlator
  //   struct fit_type fit_info;
  //   struct fit_result fit_out;

  //   fit_info.Nvar = 1;
  //   fit_info.Npar = 1;
  //   fit_info.N = 1;
  //   fit_info.Njack = Njack;
  //   fit_info.n_ext_P = 3;
  //   fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
  //   for (int j = 0; j < fit_info.Njack; j++) {
  //     fit_info.ext_P[0][j] = M_PS[j];
  //     fit_info.ext_P[1][j] = head.mus[0];
  //     fit_info.ext_P[2][j] = head.mus[0];
  //   }
  //   fit_info.function = constant_fit;
  //   fit_info.linear_fit = true;
  //   fit_info.T = head.T;
  //   fit_info.corr_id = { 0 };

  //   // c++ 0 || r 1
  //   struct fit_result f_PS = fit_fun_to_fun_of_corr(
  //     option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
  //     outfile, lhs_function_f_PS, "f_{Ds}", fit_info, jack_file);
  //   // free_fit_result(fit_info, fit_out);
  //   check_correlatro_counter(1);
  //   fit_info.restore_default();



}