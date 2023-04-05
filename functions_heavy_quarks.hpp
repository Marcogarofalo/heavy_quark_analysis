#ifndef functions_heavy_quarks_H
#define functions_heavy_quarks_H
#include "non_linear_fit.hpp"

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
#endif
