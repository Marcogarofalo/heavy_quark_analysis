#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
// #include "functions_amu.hpp"
// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"
#include "resampling_new.hpp"
#include "global.hpp"
#include "fit_all.hpp"
// #include "do_analysis_charm.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <map>



double read_single_Nobs(FILE* stream, int header_size, int Njack) {
    int Nobs;
    long int tmp;
    int s = header_size;

    // size_t i = fread(&Njack, sizeof(int), 1, stream);


    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);

    return Nobs;

}

data_single read_single_dataj(FILE* stream) {

    int Njack;
    int Nobs;

    
    data_single dj;
    dj.header.read_header_jack(stream);
    dj.header.print_header();
    dj.Nobs = read_single_Nobs(stream, dj.header.struct_size, dj.header.Njack);
    dj.Njack = dj.header.Njack;
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;

}

data_all read_all_the_files(std::vector<std::string> files, const char* resampling) {
    data_all jackall;
    jackall.resampling = resampling;
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    for (int i=0; i<files.size();i++) {
        printf("reading: %s\n",files[i].c_str());
        FILE* f = open_file(files[i].c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[i] = read_single_dataj(f);
        jackall.en[i].resampling = resampling;
        
        fclose(f);

        // printing
        printf("file: %s\n",files[i].c_str());
        jackall.en[i].header.print_header();

    }
    return jackall;

}

int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];
    char namefit[NAMESIZE];

    std::vector<std::string> files;

    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th1_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th2_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th3_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th4_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th6_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th7_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th8_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th9_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.64_th9.5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_B64 = files.size();

    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th1_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th2_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th3_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th4_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th6_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th7_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th8_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th9_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB211.072.96_th9.5_t56_44.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_B96 = files.size() - Ntheta_B64;

    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th1_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th2_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th3_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th4_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th5_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th6_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th7_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th8_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th9_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC211.06.80_th9.5_t65_51.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_C80 = files.size() - Ntheta_B64 - Ntheta_B96;

    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th1_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th2_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th3_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th4_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th5_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th6_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th7_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th8_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th9_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD211.054.96_th9.5_t78_62.dat", argv[2], argv[1]);
    files.emplace_back(namefile);

    int Ntheta_D96 = files.size() - Ntheta_B64 - Ntheta_B96 - Ntheta_C80;


    std::vector<int> all_en(files.size());
    for (int e = 0; e < files.size(); e++) {
        all_en[e] = e;
    }

    data_all jackall = read_all_the_files(files, argv[1]);
    jackall.create_generalised_resampling();






}