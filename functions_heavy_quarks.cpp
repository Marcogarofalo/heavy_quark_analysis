#define functions_heavy_quarks_C
#include "functions_heavy_quarks.hpp"
#include "tower.hpp"
#include "mutils.hpp"

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
    int dt = fit_info.myen[0] - fit_info.myen[1];
    int mu = fit_info.myen[2];
    int nu = fit_info.myen[3];

    double M = fit_info.ext_P[0][j];
    int t1 = (fit_info.myen[1] - t);
    if (t1 < 0) { r[0][0] = 0; r[0][1] = 0; return r; }
    r[0][0] = in[j][VV][t1][0] + in[j][AA][t1][0];//- in[j][VA][t][0] - in[j][AV][t][0];
    r[0][0] /= in[j][C2][t1][0] * exp(-M * dt);
    r[0][0] *= M;

    r[0][1] = in[j][VV][t1][1] + in[j][AA][t1][1];//- in[j][VA][t][0] - in[j][AV][t][0];
    r[0][1] /= in[j][C2][t1][0] * exp(-M * dt);
    r[0][1] *= M;
    
    if ((mu == 1 && nu == 2) || (mu == 2 && nu == 1)) {
        
        r[0][0] = - in[j][VA][t][0] - in[j][AV][t][0];
        r[0][0] /= in[j][C2][t1][0] * exp(-M * dt);
        r[0][0] *= M;

        r[0][1] = - in[j][VA][t][0] - in[j][AV][t][0];
        r[0][1] /= in[j][C2][t1][0] * exp(-M * dt);
        r[0][1] *= M;
    }

    return r;
}


double** compute_Y1(int j, double**** in, int t, struct fit_type fit_info) {
    double** Y1 = malloc_2<double>(1, 2);
    int W11 = fit_info.corr_id[0];
    int W22 = fit_info.corr_id[1];
    int id_input = 2;
    error(fit_info.corr_id.size() != id_input, 1, "compute Y1", "fit_info.corr_id.size() must be %d, intead it is %d",
        id_input, fit_info.corr_id.size());
    Y1[0][0] = (in[j][W11][t][0] + in[j][W22][t][0])/2.0 ;
    Y1[0][1] = (in[j][W11][t][1] + in[j][W22][t][1])/2.0 ;
    return Y1;
}

double** compute_Y2(int j, double**** in, int t, struct fit_type fit_info) {
    double** Y = malloc_2<double>(1, 2);
    int W00 = fit_info.corr_id[0];
    int id_input = 1;
    error(fit_info.corr_id.size() != id_input, 1, "compute Y2", "fit_info.corr_id.size() must be %d, intead it is %d",
        id_input, fit_info.corr_id.size());
    Y[0][0] = (in[j][W00][t][0]);
    Y[0][1] = (in[j][W00][t][1]);
    return Y;
}
double** compute_Y3(int j, double**** in, int t, struct fit_type fit_info) {
    double** Y = malloc_2<double>(1, 2);
    int W33 = fit_info.corr_id[0];
    int id_input = 1;
    error(fit_info.corr_id.size() != id_input, 1, "compute Y3", "fit_info.corr_id.size() must be %d, intead it is %d",
        id_input, fit_info.corr_id.size());
    Y[0][0] = -(in[j][W33][t][0]);
    Y[0][1] = -(in[j][W33][t][1]);
    return Y;
}
double** compute_Y4(int j, double**** in, int t, struct fit_type fit_info) {
    double** Y = malloc_2<double>(1, 2);
    int W03 = fit_info.corr_id[0];
    int W30 = fit_info.corr_id[1];
    int id_input = 2;
    error(fit_info.corr_id.size() != id_input, 1, "compute Y4", "fit_info.corr_id.size() must be %d, intead it is %d",
        id_input, fit_info.corr_id.size());
    Y[0][0] = -(in[j][W03][t][1] + in[j][W30][t][1])/2.0;
    Y[0][1] = (in[j][W03][t][0] + in[j][W30][t][0])/2.0;
    return Y;
}
double** compute_Y5(int j, double**** in, int t, struct fit_type fit_info) {
    double** Y = malloc_2<double>(1, 2);
    int W12 = fit_info.corr_id[0];
    int W21 = fit_info.corr_id[1];
    int id_input = 2;
    error(fit_info.corr_id.size() != id_input, 1, "compute Y5", "fit_info.corr_id.size() must be %d, intead it is %d",
        id_input, fit_info.corr_id.size());
    // there is an i in front
    Y[0][0] = (in[j][W12][t][1] - in[j][W21][t][1]) / 2.0;
    Y[0][1] = -(in[j][W12][t][0] - in[j][W21][t][0]) / 2.0;
    return Y;
}

double** compute_Z_factors(int j, double**** in, int t, struct fit_type fit_info) {
    double** Z = malloc_2<double>(3, 2);
    int Y1 = fit_info.corr_id[1];
    int Y2 = fit_info.corr_id[2];
    int Y3 = fit_info.corr_id[3];
    int Y4 = fit_info.corr_id[4];
    int id_input = 6;
    error(fit_info.corr_id.size() != id_input, 1, "compute Y1", "fit_info.corr_id.size() must be %d, intead it is %d",
        id_input, fit_info.corr_id.size());
    Z[0][0] = (in[j][Y2][t][0] + in[j][Y3][t][0] - 2 * in[j][Y4][t][0]);
    Z[0][1] = (in[j][Y2][t][1] + in[j][Y3][t][1] - 2 * in[j][Y4][t][1]);

    Z[1][0] = 2.0 * (in[j][Y3][t][0] - 2.0 * in[j][Y1][t][0] - in[j][Y4][t][0]);
    Z[1][1] = 2.0 * (in[j][Y3][t][1] - 2.0 * in[j][Y1][t][1] - in[j][Y4][t][1]);

    Z[2][0] = (in[j][Y3][t][0] - 2.0 * in[j][Y1][t][0]);
    Z[2][1] = (in[j][Y3][t][1] - 2.0 * in[j][Y1][t][1]);
    return Z;
}

double lhs_function_ZPS(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double corr = in[j][id][t][0];
    double M = fit_info.ext_P[0][j];


    double me = corr / (exp(-t * M) + exp(-(fit_info.T - t) * M));
    return me;
}

double lhs_2fit_par(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    return in[j][id][t][0];
}
double rhs_2fit_par(int n, int Nvar, double* x, int Npar, double* P) {
    int t = x[0];
    return 0.5 * P[1] * P[1] * (exp(-(P[0] * t)) + exp(-(P[0] * (file_head.l0 - t))));
}


double lhs_function_f_PS_ss_ls(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double amp = 0;
    double M = fit_info.ext_P[0][j];
    double me_ss = fit_info.ext_P[1][j];
    double mu1 = fit_info.ext_P[2][j];
    double mu2 = fit_info.ext_P[3][j];


    double me = me_ss * in[j][fit_info.corr_id[0]][t][0] / in[j][fit_info.corr_id[1]][t][0];
    return (mu1 + mu2) * me / (M * sinh(M));
}


double lhs_function_me(int j, double**** in, int t, struct fit_type fit_info) {
    int id = fit_info.corr_id[0];
    double amp = 0;
    double M = fit_info.ext_P[0][j];


    double me = in[j][fit_info.corr_id[0]][t][0] * 2 * M / (exp(-t * M) + exp(-(fit_info.T - t) * M));
    return sqrt(me);
}