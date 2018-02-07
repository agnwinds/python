#!/bin/bash
nlim=10
fraclim=0.5
Setup_Py_Dir
$1 -z py_hydro > output
mv py_heatcool.dat py_heatcool_init.dat
./test_heatcool.py py_heatcool_init.dat $nlim $fraclim > init_output.dat 2> tmp
if [ $(tail -1 tmp) -eq 0 ]
then
  echo "Initial test ran OK "
else
  echo "Initial test had problems, check init_output.dat"
fi
rm tmp
cp py_hydro.wind_save py_hydro_restart.wind_save
$1 -z -r py_hydro_restart > output2
mv py_heatcool.dat py_heatcool_restart.dat
./test_heatcool.py py_heatcool_restart.dat $nlim $fraclim > restart_output.dat 2> tmp
if [ $(tail -1 tmp) -eq 0 ]
then
  echo "Restart test ran OK "
else
  echo "Restart test had problems, check restart_output.dat"
fi
