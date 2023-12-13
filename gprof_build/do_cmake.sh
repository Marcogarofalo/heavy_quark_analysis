rm -rv CMakeCache.txt
rm -rv CMakeFiles

cmake -DCMAKE_PREFIX_PATH=/home/garofalo/analysis/analysis_program/gproof_build/install_dir \
      -DCMAKE_CXX_FLAGS="-g -O2 -pg -no-pie"  \
      -DWITH_ARB=ON \
      -DCMAKE_VERBOSE_MAKEFILE=ON \
      \
      ..

