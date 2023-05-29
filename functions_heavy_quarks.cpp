#define functions_heavy_quarks_C
#include "functions_heavy_quarks.hpp"
#include "tower.hpp"

double lhs_function_heavy_quarks_eg(int j, double**** in, int t, struct fit_type fit_info) {
    double r = in[j][0][t][0];
    return r;
}


mass_index::mass_index(int nmus) {
    int k1, k2, r1, r2, i;
    nk = nmus;

    table = (int****)malloc(sizeof(int***) * nk);
    for (k1 = 0;k1 < nk;k1++) {
        table[k1] = (int***)malloc(sizeof(int**) * 2);
        for (r1 = 0;r1 < 2;r1++) {
            table[k1][r1] = (int**)malloc(sizeof(int*) * (k1 + 1));
            for (k2 = 0;k2 <= k1;k2++) {
                table[k1][r1][k2] = (int*)malloc(sizeof(int) * 2);
            }
        }
    }

    tot_mu_r = 0;
    for (k1 = 0;k1 < nk;k1++)
        for (r1 = 0;r1 < 2;r1++)
            for (k2 = 0;k2 <= k1;k2++)
                for (r2 = 0;r2 < 2;r2++) {
                    table[k1][r1][k2][r2] = tot_mu_r;
                    tot_mu_r++;
                }


}


mass_index::~mass_index() {

    for (int k1 = 0;k1 < nk;k1++) {
        for (int r1 = 0;r1 < 2;r1++) {
            for (int k2 = 0;k2 <= k1;k2++) {
                free(table[k1][r1][k2]);
            }
            free(table[k1][r1]);
        }
        free(table[k1]);
    }
    free(table);

}

id_contraction::id_contraction(int nmus, int ngamma) : mass_index(nmus) {
    Ngamma = ngamma;
};
int id_contraction::return_id(int k1, int k2, int r1, int r2, int ig, int is) {
    return table[k1][r1][k2][r2] + tot_mu_r * (ig + Ngamma * (is));

}


double lhs_M_K_inter(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[0]][j];

    return r;
}

double rhs_M_K_linear(int n, int Nvar, double* x, int Npar, double* P) {
    return P[0] + x[0] * P[1];
}

double rhs_M_Ds_linear(int n, int Nvar, double* x, int Npar, double* P) {
    return P[0] + x[0] * P[1] + x[1] * P[2];
}

double** compute_Mmunu(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = malloc_2<double>(1, 2);
    int VV = fit_info.corr_id[0];
    int AA = fit_info.corr_id[1];
    int VA = fit_info.corr_id[2];
    int AV = fit_info.corr_id[3];
    int C2 = fit_info.corr_id[4];

    double ZDs = fit_info.ext_P[0][j];
    r[0][0] = in[j][VV][t][0] + in[j][AA][t][0] - in[j][VA][t][0] - in[j][AV][t][0];
    r[0][0] /= in[j][C2][t][0] * in[j][C2][fit_info.t0_GEVP][0];
    r[0][0] *= ZDs;
    r[0][1] = in[j][VV][t][1] + in[j][AA][t][1] - in[j][VA][t][1] - in[j][AV][t][1];
    r[0][1] /= in[j][C2][t][0] * in[j][C2][fit_info.t0_GEVP][0];
    r[0][1] *= ZDs;
    

    return r;
}


double lhs_function_ZPS(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double corr = in[j][id][t][0];
    double M = fit_info.ext_P[0][j];
    double mu1 = fit_info.ext_P[1][j];
    double mu2 = fit_info.ext_P[2][j];

    double me = corr / (exp(-t * M) + exp(-(fit_info.T - t) * M));
    return me;
}