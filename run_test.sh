#!/bin/bash

ROOT=$(cd $(dirname $BASH_SOURCE) && pwd)

cd ${ROOT}/build
make

echo "========================================================================"
echo "                       Rossby-Haurwitz test                             "
if [[ ! -d ${ROOT}/test/rh_01 ]]; then
  mkdir -p ${ROOT}/test/rh_01
fi
cd ${ROOT}/test/rh_01
cp ${ROOT}/run/namelist.rh_test .
rh_start_time=$(date '+%s')
${ROOT}/build/dycore_test.exe namelist.rh_test
rh_end_time=$(date '+%s')
ncl -Q ${ROOT}/src/test_cases/barotropic/plot_rossby_haurwitz_wave_test.ncl \
  file_prefix=\"rh_test.360x181.dt240\"

echo "========================================================================"
echo "                      Mountain zonal flow test                          "
if [[ ! -d ${ROOT}/test/mz_01 ]]; then
  mkdir -p ${ROOT}/test/mz_01
fi
cd ${ROOT}/test/mz_01
cp ${ROOT}/run/namelist.mz_test .
mz_start_time=$(date '+%s')
${ROOT}/build/dycore_test.exe namelist.mz_test
mz_end_time=$(date '+%s')
ncl -Q ${ROOT}/src/test_cases/barotropic/plot_mountain_zonal_flow_test.ncl \
  file_prefix=\"mz_test.360x181.dt240\"

echo "========================================================================"
echo "                          Jet zonal flow test                           "
if [[ ! -d ${ROOT}/test/jz_01 ]]; then
  mkdir -p ${ROOT}/test/jz_01
fi
cd ${ROOT}/test/jz_01
cp ${ROOT}/run/namelist.jz_test .
jz_start_time=$(date '+%s')
${ROOT}/build/dycore_test.exe namelist.jz_test
jz_end_time=$(date '+%s')
ncl -Q ${ROOT}/src/test_cases/barotropic/plot_jet_zonal_flow_test.ncl \
  file_prefix=\"jz_test.360x181.diffused.dt240\"

echo "========================================================================"
echo "                          Run time report                               "
echo
echo "Running duration (minutes):"
echo "Rossby-Haurwitz test:     $(((rh_end_time - rh_start_time) / 60))"
echo "Mountain zonal flow test: $(((mz_end_time - mz_start_time) / 60))"
echo "Jet zonal flow test:      $(((jz_end_time - jz_start_time) / 60))"
