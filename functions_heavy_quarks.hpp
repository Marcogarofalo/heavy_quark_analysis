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

#endif
