#!/bin/bash
rm -fr mod
rm cells3d
mkdir mod
dir='./src'
if [ $1 == '-g' ]; then
	ifort -g -traceback -check bounds -check all -stand f03 -r8 -module mod $dir/global_m.F90 $dir/misc_m.F90 $dir/derivatives_m.F90 $dir/sim_init_m.F90 $dir/run_cells_m.F90  $dir/main.F90 -o cells3d_dbg
else
	ifort -g -traceback -stand f03 -O3 -fast -r8 -module mod $dir/global_m.F90 $dir/misc_m.F90 $dir/derivatives_m.F90 $dir/sim_init_m.F90 $dir/run_cells_m.F90  $dir/main.F90 -o cells3d
fi
# debug
#ifort -g -traceback -check bounds -check all -stand f03 -r8 -module mod $dir/global_m.F90 $dir/misc_m.F90 $dir/derivatives_m.F90 $dir/sim_init_m.F90 $dir/run_cells_m.F90  $dir/main.F90 -o cells3d
