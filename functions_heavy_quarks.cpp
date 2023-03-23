#define functions_heavy_quarks_C
#include "functions_heavy_quarks.hpp"

double lhs_function_heavy_quarks_eg(int j, double**** in, int t, struct fit_type fit_info){
    double r=in[j][0][t][0];
    return r;
}


mass_index::mass_index(int nmus)
{
     int k1, k2,r1,r2,i;
     nk=nmus;

     table=(int****) malloc(sizeof(int***)*nk);
     for (k1=0;k1<nk;k1++){
     	table[k1]=(int***) malloc(sizeof(int**)*2);
     		for (r1=0;r1<2;r1++){
     			table[k1][r1]=(int**) malloc(sizeof(int*)*(k1+1));
     			for (k2=0;k2<=k1;k2++){
	     			table[k1][r1][k2]=(int*) malloc(sizeof(int)*2);		
			}
		}
     }

     tot_mu_r=0;
     for (k1=0;k1<nk;k1++)
     for (r1=0;r1<2;r1++)
     for (k2=0;k2<=k1;k2++)
     for (r2=0;r2<2;r2++){
    	table[k1][r1][k2][r2]=i;
    	tot_mu_r++;
    }


}


mass_index::~mass_index()
{
    
     for (int k1=0;k1<nk;k1++){
     		for (int r1=0;r1<2;r1++){
     			for (int k2=0;k2<=k1;k2++){
                    free(table[k1][r1][k2]);
			}
            free(table[k1][r1]);
		}
        free(table[k1]);
     }
     free(table);

}

id_contraction::id_contraction(int nmus, int ngamma): mass_index(nmus){
    Ngamma=ngamma;
};
int id_contraction::return_id(int k1, int k2, int r1, int r2, int ig, int is){
    return table[k1][r1][k2][r2]+tot_mu_r*(ig+Ngamma*(is) );

}

