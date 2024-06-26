#ifndef functions_heavy_quarks_H
#define functions_heavy_quarks_H
#include "non_linear_fit.hpp"
#include "arb_calc.h"

double lhs_function_heavy_quarks_eg(int j, double**** in, int t, struct fit_type fit_info);

class mass_index {
public:
    int nk;
    int**** table;
    int tot_mu_r;
    mass_index(int nmus);
    ~mass_index();
};

class id_contraction: mass_index {
public:
    id_contraction(int nmus, int ngamma);
    int return_id(int k1, int k2, int r1, int r2, int ig, int is);
    int Ngamma;
};

double lhs_M_K_inter(int n, int e, int j, data_all gjack, struct fit_type fit_info);
double rhs_M_K_linear(int n, int Nvar, double* x, int Npar, double* P);
double rhs_M_Ds_linear(int n, int Nvar, double* x, int Npar, double* P);

double** compute_Mmunu(int j, double**** in, int t, struct fit_type fit_info);
double** compute_Y1(int j, double**** in, int t, struct fit_type fit_info);
double** compute_Y2(int j, double**** in, int t, struct fit_type fit_info);
double** compute_Y3(int j, double**** in, int t, struct fit_type fit_info);
double** compute_Y4(int j, double**** in, int t, struct fit_type fit_info);
double** compute_Y5(int j, double**** in, int t, struct fit_type fit_info);
double** compute_Z_factors(int j, double**** in, int t, struct fit_type fit_info);

double lhs_function_ZPS(int j, double**** in, int t, struct fit_type fit_info);
double lhs_function_f_PS_ss_ls(int j, double**** in, int t, struct fit_type fit_info);
double lhs_function_me(int j, double**** in, int t, struct fit_type fit_info);

double lhs_2fit_par(int j, double**** in, int t, struct fit_type fit_info) ;
double rhs_2fit_par(int n, int Nvar, double* x, int Npar, double* P);

double sigma2_fit(int n, int Nvar, double* x, int Npar, double* P);
double sigma4_fit(int n, int Nvar, double* x, int Npar, double* P);

int c_thetap_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
int c_thetam_s_HLT(acb_ptr res, const acb_t z, void* param, slong order, slong prec);
#endif
