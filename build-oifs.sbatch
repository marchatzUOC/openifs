#!/bin/bash
#SBATCH -A bsc32
###SBATCH -q gp_bsces
#SBATCH -q gp_debug
#SBATCH -c 112
#SBATCH --time 01:30:00

source ./oifs-config.bsc-mn5-gpp.sh
$OIFS_HOME/scripts/build_test/openifs-test.sh -b -j 112 #--clean 
