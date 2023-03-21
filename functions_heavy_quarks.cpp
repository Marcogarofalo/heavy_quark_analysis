#define functions_heavy_quarks_C
#include "functions_heavy_quarks.hpp"

double lhs_function_heavy_quarks_eg(int j, double**** in, int t, struct fit_type fit_info){
    double r=in[j][0][t][0];
    return r;
}

double lhs_function_f_PS(int j, double**** in, int t, struct fit_type fit_info){
    int id=fit_info.corr_id[0];
    double corr=in[j][id][t][0];
    double M=fit_info.ext_P[0][j];
    double mu1=fit_info.ext_P[1][j];
    double mu2=fit_info.ext_P[2][j];

    double me=sqrt(corr*2*M/(exp(-t*M)+exp(-(fit_info.T-t)*M)));
    return (mu1+mu2)*me/(M*sinh(M));
}
