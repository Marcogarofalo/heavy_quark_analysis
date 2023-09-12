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
  for (int i = 0; i < head.ncorr;i++) {
    symmetrise_jackboot(Njack, i, head.T, conf_jack);
  }

  // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

  ////////////////////

  id_contraction  ids(head.mus.size(), head.gammas.size());

  //////////////////////////////
  // start fitting
  /////////////////////////////
  corr_counter = -1;
  double* M_PS = plateau_correlator_function(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    namefile_plateaux, outfile, 0, "M_{PS}", M_eff_T, jack_file);
  check_correlatro_counter(0);

  // eg of fit to correlator
  struct fit_type fit_info;
  struct fit_result fit_out;

  fit_info.Nvar = 1;
  fit_info.Npar = 1;
  fit_info.N = 1;
  fit_info.Njack = Njack;
  fit_info.n_ext_P = 3;
  fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
  for (int j = 0; j < fit_info.Njack; j++) {
    fit_info.ext_P[0][j] = M_PS[j];
    fit_info.ext_P[1][j] = head.mus[0];
    fit_info.ext_P[2][j] = head.mus[0];
  }
  fit_info.function = constant_fit;
  fit_info.linear_fit = true;
  fit_info.T = head.T;
  fit_info.corr_id = { 0 };

  // c++ 0 || r 1
  struct fit_result f_PS = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS, "f_{PS}", fit_info, jack_file);
  // free_fit_result(fit_info, fit_out);
  check_correlatro_counter(1);
  free(fit_info.ext_P[0]);free(fit_info.ext_P[1]);free(fit_info.ext_P[2]);
  fit_info.restore_default();



  ////////////////////////////////////////////
  // GEVP
  ////////////////////////////////////////////
  int ncorr_new = head.ncorr;
  fit_info.N = 2;
  fit_info.Njack = Njack;

  fit_info.corr_id = { ids.return_id(0,0,0,0,0,0), ids.return_id(0,0,0,0,0,1),
                       ids.return_id(0,0,0,0,0,2), ids.return_id(0,0,0,0,0,3) };//diag{ ll, ss}
  fit_info.value_or_vector = 0; // 0= values
  fit_info.t0_GEVP = 3;
  fit_info.GEVP_ignore_warning_after_t = 1;
  fit_info.verbosity = 0;
  printf("GEVP_00\n");
  add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
  printf(" ncorr after GEVP %d\n", ncorr_new);
  int id_GEVP_00 = ncorr_new - 2;

  fit_info.value_or_vector = 2; // 2= Amplitudes
  fit_info.N = 2 * 2;
  add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
  int id_GEVP_00_a0_0 = ncorr_new - 4; // gives <0|P5|P5>
  int id_GEVP_00_a0_1 = ncorr_new - 3; // gives <0|P5smeared|P5>
  int id_GEVP_00_a1_0 = ncorr_new - 2; // gives <0|P5|P5next-state>
  int id_GEVP_00_a1_1 = ncorr_new - 1; // gives <0|P5smeared|P5next-state>

  fit_info.restore_default();

  double* lambda0 = plateau_correlator_function(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    namefile_plateaux, outfile, id_GEVP_00, "lambda0_00", identity, jack_file);
  check_correlatro_counter(2);

  double* M_PS_GEVP = plateau_correlator_function(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    namefile_plateaux, outfile, id_GEVP_00, "M_{PS-00-GEVP}", M_eff_T, jack_file);
  check_correlatro_counter(3);


  fit_info.Nvar = 1;
  fit_info.Npar = 1;
  fit_info.N = 1;
  fit_info.Njack = Njack;
  fit_info.n_ext_P = 3;
  fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
  for (int j = 0; j < fit_info.Njack; j++) {
    fit_info.ext_P[0][j] = M_PS_GEVP[j];
    fit_info.ext_P[1][j] = head.mus[0];
    fit_info.ext_P[2][j] = head.mus[0];
  }
  fit_info.function = constant_fit;
  fit_info.linear_fit = true;
  fit_info.T = head.T;
  fit_info.corr_id = { id_GEVP_00_a0_0 };


  struct fit_result f_PS_GEVP = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS_GEVP, "f_{PS-00-GEVP}", fit_info, jack_file);
  // free_fit_result(fit_info, fit_out);
  check_correlatro_counter(4);

  // fit_info.corr_id[0] = id_GEVP_00_a0_1;
  // struct fit_result f_PS_GEVP1 = fit_fun_to_fun_of_corr(
  //   option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
  //   outfile, lhs_function_f_PS_GEVP, "f_{PS-001-GEVP}", fit_info, jack_file);
  // // free_fit_result(fit_info, fit_out);
  // check_correlatro_counter(4);

  fit_info.restore_default();

  ///////////// structure for fits
  data_all jackall;
  jackall.resampling = argv[6];
  jackall.ens = 3;
  jackall.en = new data_single[jackall.ens];
  for (int i = 0;i < jackall.ens;i++) {
    jackall.en[i].header = head;
    jackall.en[i].Nobs = 1;
    jackall.en[i].Njack = head.Njack;
    jackall.en[i].jack = (double**)malloc(sizeof(double*) * jackall.en[i].Nobs);
  }

  ///// KAON

  for (int is = 0;is < 3;is++) {
    fit_info.Njack = Njack;
    fit_info.N = 2;
    fit_info.corr_id = { ids.return_id(is + 1,0,0,0,0,0), ids.return_id(is + 1,0,0,0,0,1),
                         ids.return_id(is + 1,0,0,0,0,2), ids.return_id(is + 1,0,0,0,0,3) };//diag{ ll, ss}
    fit_info.value_or_vector = 0; // 0= values
    fit_info.t0_GEVP = 3;
    fit_info.GEVP_ignore_warning_after_t = 1;
    fit_info.verbosity = 0;

    printf("GEVP_00\n");
    add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
    printf(" ncorr after GEVP %d\n", ncorr_new);
    int id_GEVP_l0 = ncorr_new - 2;
    char name_obs[NAMESIZE];
    mysprintf(name_obs, NAMESIZE, "M_{K}(s%d)", is);
    jackall.en[is].jack[0] = plateau_correlator_function(
      option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
      namefile_plateaux, outfile, id_GEVP_l0, name_obs, M_eff_T, jack_file);
    check_correlatro_counter(5 + is);
  }
  // fit the kaon
  fit_info.restore_default();
  fit_info.Npar = 2;
  fit_info.Nvar = 1;
  fit_info.Njack = head.Njack;
  fit_info.N = 1;
  fit_info.myen = { 0,1,2 };
  fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
  int count = 0;
  for (int n = 0;n < fit_info.N;n++) {
    for (int e : fit_info.myen) {
      for (int j = 0;j < fit_info.Njack;j++) {
        fit_info.x[0][count][j] = head.mus[e + 1]; // mus1
      }
      count++;
    }
  }
  fit_info.corr_id = { 0 };
  fit_info.function = rhs_M_K_linear;
  fit_info.linear_fit = true;
  fit_info.covariancey = false;
  fit_info.verbosity = 0;
  char namefit[NAMESIZE];
  mysprintf(namefit, NAMESIZE, "fit_MK_vs_mu");
  char** temp_argv = malloc_2<char>(5, NAMESIZE);
  mysprintf(temp_argv[1], NAMESIZE, "%s", argv[6]);// resampling
  mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

  fit_result fit_inter_MK = fit_all_data(temp_argv, jackall, lhs_M_K_inter, fit_info, namefit);

  fit_info.band_range = { head.mus[1] * 0.9,head.mus[3] * 1.1 };
  print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "amu", fit_inter_MK, fit_inter_MK, 0, fit_info.myen.size() - 1, 0.0005);

  fit_info.band_range = { .01692 ,0.01692 * 1.1 };
  print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "MK_interpolated", fit_inter_MK, fit_inter_MK, 0, 0, 100);

  ///////////// structure for fits
  data_all jackall_c;
  jackall_c.resampling = argv[6];
  jackall_c.ens = 3 * 4;
  jackall_c.en = new data_single[jackall_c.ens];
  for (int i = 0;i < jackall_c.ens;i++) {
    jackall_c.en[i].header = head;
    jackall_c.en[i].Nobs = 1;
    jackall_c.en[i].Njack = head.Njack;
    jackall_c.en[i].jack = (double**)malloc(sizeof(double*) * jackall_c.en[i].Nobs);
  }

  ///// Ds
  int Ns = 3, Nc = 3;
  for (int is = 0;is < Ns;is++) {
    for (int ic = 0;ic < Nc;ic++) {
      fit_info.Njack = Njack;
      fit_info.N = 2;
      fit_info.corr_id = { ids.return_id(ic + 5,is + 1,0,0,0,0), ids.return_id(ic + 5,is + 1,0,0,0,1),
                           ids.return_id(ic + 5,is + 1,0,0,0,2), ids.return_id(ic + 5,is + 1,0,0,0,3) };//diag{ ll, ss}
      fit_info.value_or_vector = 0; // 0= values
      fit_info.t0_GEVP = 3;
      fit_info.GEVP_ignore_warning_after_t = 1;
      fit_info.verbosity = 0;

      printf("GEVP_00\n");
      add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
      printf(" ncorr after GEVP %d\n", ncorr_new);
      int id_GEVP_l0 = ncorr_new - 2;
      char name_obs[NAMESIZE];
      mysprintf(name_obs, NAMESIZE, "M_{Ds}^{GEVP}(c%d,s%d)", ic, is);
      jackall_c.en[ic + is * 3].jack[0] = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_GEVP_l0, name_obs, M_eff_T, jack_file);
      check_correlatro_counter(8 + ic + is * Nc);
    }
  }
  // fit the Ds
  fit_info.restore_default();
  fit_info.Npar = 3;
  fit_info.Nvar = 2;
  fit_info.Njack = head.Njack;
  fit_info.N = 1;
  fit_info.myen = std::vector<int>(Ns * Nc);
  for (int n = 0;n < fit_info.myen.size();n++) fit_info.myen[n] = n;
  fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
  count = 0;
  for (int n = 0;n < fit_info.N;n++) {
    for (int e : fit_info.myen) {
      for (int j = 0;j < fit_info.Njack;j++) {
        int ic = e % Nc;
        int is = e / Nc;
        fit_info.x[0][count][j] = head.mus[ic + 5]; // mus1
        fit_info.x[1][count][j] = head.mus[is + 1]; // mus1
      }
      count++;
    }
  }
  fit_info.corr_id = { 0 };
  fit_info.function = rhs_M_Ds_linear;
  fit_info.linear_fit = true;
  fit_info.covariancey = false;
  fit_info.verbosity = 0;
  mysprintf(namefit, NAMESIZE, "fit_MDs_vs_mu");

  fit_result fit_inter_MDs = fit_all_data(temp_argv, jackall_c, lhs_M_K_inter, fit_info, namefit);

  fit_info.band_range = { head.mus[5] * 0.9,head.mus[7] * 1.1 };
  for (int is = 0;is < Ns;is++) {
    char name_band[NAMESIZE];
    mysprintf(name_band, NAMESIZE, "amuc_%d", is);
    print_fit_band(temp_argv, jackall_c, fit_info, fit_info, namefit, name_band, fit_inter_MDs, fit_inter_MDs, 0, 0 + is * Nc, 0.0005);
  }
  fit_info.band_range = { 0.2368 ,0.2368 * 1.1 };
  std::vector<double> xcont = { 0.2368, 0.01692 };
  print_fit_band(temp_argv, jackall_c, fit_info, fit_info, namefit, "MDs_interpolated", fit_inter_MDs, fit_inter_MDs, 0, 0, 100, xcont);


  double* M_Ds_ll_c0s0 = plateau_correlator_function(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    namefile_plateaux, outfile, ids.return_id(5, 1, 0, 0, 0, 0), "M_{Ds}^{ll}(c0,s0)", M_eff_T, jack_file);
  check_correlatro_counter(17);

  double* M_Ds_sl_c0s0 = plateau_correlator_function(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    namefile_plateaux, outfile, ids.return_id(5, 1, 0, 0, 0, 1), "M_{Ds}^{sl}(c0,s0)", M_eff_T, jack_file);
  check_correlatro_counter(18);

  double* M_Ds_ss_c0s0 = plateau_correlator_function(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
    namefile_plateaux, outfile, ids.return_id(5, 1, 0, 0, 0, 3), "M_{Ds}^{ss}(c0,s0)", M_eff_T, jack_file);
  check_correlatro_counter(19);
  /////////////////////////////////////


  fit_info.Nvar = 1;
  fit_info.Npar = 1;
  fit_info.N = 1;
  fit_info.Njack = Njack;
  fit_info.n_ext_P = 3;
  fit_info.malloc_ext_P();
  for (int j = 0; j < fit_info.Njack; j++) {
    fit_info.ext_P[0][j] = M_Ds_ll_c0s0[j];
    fit_info.ext_P[1][j] = head.mus[5];
    fit_info.ext_P[2][j] = head.mus[1];
  }
  fit_info.function = constant_fit;
  fit_info.linear_fit = true;
  fit_info.T = head.T;

  fit_info.corr_id = { ids.return_id(5, 1, 0, 0, 0, 0) };
  for (int j = 0; j < fit_info.Njack; j++)     fit_info.ext_P[0][j] = M_Ds_ll_c0s0[j];
  struct fit_result f_PS_Ds = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS, "f_{Ds}_ll(c0,s0)", fit_info, jack_file);
  free_fit_result(fit_info, f_PS_Ds);
  check_correlatro_counter(20);

  fit_info.corr_id = { ids.return_id(5, 1, 0, 0, 0, 1) };
  for (int j = 0; j < fit_info.Njack; j++)     fit_info.ext_P[0][j] = M_Ds_sl_c0s0[j];
  f_PS_Ds = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS, "f_{Ds}_sl(c0,s0)", fit_info, jack_file);
  free_fit_result(fit_info, f_PS_Ds);
  check_correlatro_counter(21);

  fit_info.corr_id = { ids.return_id(5, 1, 0, 0, 0, 3) };
  for (int j = 0; j < fit_info.Njack; j++)     fit_info.ext_P[0][j] = M_Ds_ss_c0s0[j];
  struct fit_result f_PS_Ds_ss = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS, "f_{Ds}_ss(c0,s0)", fit_info, jack_file);
  free_fit_result(fit_info, f_PS_Ds_ss);
  check_correlatro_counter(22);


  struct fit_result me_Ds_ss = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_me, "me_{Ds}_ss(c0,s0)", fit_info, jack_file);
  check_correlatro_counter(23);
  fit_info.restore_default();


  ///////////////////////////////////////////////////// GEVP Ds operator
  fit_info.value_or_vector = 2; // 2= Amplitudes
  fit_info.N = 2 * 2;
  int ic = 0, is = 0;
  fit_info.corr_id = { ids.return_id(ic + 5,is + 1,0,0,0,0), ids.return_id(ic + 5,is + 1,0,0,0,1),
                           ids.return_id(ic + 5,is + 1,0,0,0,2), ids.return_id(ic + 5,is + 1,0,0,0,3) };//diag{ ll, ss}

  fit_info.t0_GEVP = 3;
  fit_info.GEVP_ignore_warning_after_t = 1;
  fit_info.verbosity = 0;
  add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
  int id_GEVP_Ds_a0_0 = ncorr_new - 4;

  fit_info.Nvar = 1;
  fit_info.Npar = 1;
  fit_info.N = 1;
  fit_info.Njack = Njack;
  fit_info.n_ext_P = 3;
  fit_info.ext_P = malloc_2<double>(fit_info.n_ext_P, fit_info.Njack);
  for (int j = 0; j < fit_info.Njack; j++) {
    fit_info.ext_P[0][j] = jackall_c.en[ic + is * 3].jack[0][j];
    fit_info.ext_P[1][j] = head.mus[5];
    fit_info.ext_P[2][j] = head.mus[1];
  }
  fit_info.function = constant_fit;
  fit_info.linear_fit = true;
  fit_info.T = head.T;
  fit_info.corr_id = { id_GEVP_Ds_a0_0 };


  struct fit_result f_Ds_00_GEVP = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS_GEVP, "f_{Ds-00-GEVP}", fit_info, jack_file);
  check_correlatro_counter(24);
  fit_info.restore_default();
  for (int i = 0; i < fit_info.n_ext_P; i++) free(fit_info.ext_P[i]);

  ///////////////////////////////////////////////////// ss and sl for matrix element
  fit_info.Nvar = 1;
  fit_info.Npar = 1;
  fit_info.N = 1;
  fit_info.Njack = Njack;
  fit_info.n_ext_P = 4;
  fit_info.malloc_ext_P();
  for (int j = 0; j < fit_info.Njack; j++) {
    fit_info.ext_P[0][j] = M_Ds_ss_c0s0[j];
    fit_info.ext_P[1][j] = me_Ds_ss.P[0][j];
    fit_info.ext_P[2][j] = head.mus[5];
    fit_info.ext_P[3][j] = head.mus[1];
  }
  fit_info.function = constant_fit;
  fit_info.linear_fit = true;
  fit_info.T = head.T;
  fit_info.corr_id = { ids.return_id(ic + 5,is + 1,0,0,0,1), ids.return_id(ic + 5,is + 1,0,0,0,3) };


  struct fit_result f_Ds_00_ss_ls = fit_fun_to_fun_of_corr(
    option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    outfile, lhs_function_f_PS_ss_ls, "f_{Ds-00-ss-ls}", fit_info, jack_file);
  check_correlatro_counter(25);

  fit_info.restore_default();


}