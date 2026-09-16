#
#
#   oifs-config.edit_me.sh
#
#
#   This script sets the environment for OpenIFS 48r1
#
#
#   Read this script using the command:
#
#   source ./oifs-config.edit_me.sh
#
#
#--- set machine specific settings -----------------------------

export OIFS_HOST="bsc"
export OIFS_PLATFORM="mn5-gpp"

#--- if OIFS_DEBUGGER is set, oifs-run can use it to activate debugger ---------
export OIFS_DEBUGGER="DDT"

#--- set principal OIFS variables ----------------------------------------------

export OIFS_CYCLE=48r1

#---Base code assumes openifs-48r1 and openifs-expt are installed
#---in $HOME. Either these can be changed by the user------------
export OIFS_HOME="/gpfs/scratch/bsc32/bsc032198/gitlab_bsc/openifs-48r1"

#---It is recommended that the openifs-expt and oifs_data dir
#---exist in a location designed for permanent storage-----------
export OIFS_EXPT="/gpfs/scratch/bsc32/bsc032198/gitlab_bsc/exp"
export OIFS_DATA_DIR="/gpfs/projects/bsc32/models/48r1"

#---Set the path for the arch directory. Depending on system,i.e.,
#---all libs are installed on the sytem, this is not required,
#---so set to an empty string OIFS_ARCH=""
export OIFS_ARCH="./arch/${OIFS_HOST}/${OIFS_PLATFORM}"

#---Path to the executable for 3d global model. This is the
#---default path for the exe, produced by openifs-test.sh.
#---DP means double precision. To run single precision change
#---DP to SP
export OIFS_EXEC="${OIFS_HOME}/build/bin/ifsMASTER.DP"

#---Default assumed paths, only change if you know what you are doing
export OIFS_TEST="${OIFS_HOME}/scripts/build_test"
export OIFS_RUN_SCRIPT="${OIFS_HOME}/scripts/exp_3d"
export OIFS_LOGFILE="${OIFS_HOME}/oifs_test_log.txt"

alias oenv="env -0 | sort -z | tr '\0' '\n' | grep -a OIFS_"

echo -e "\nOpenIFS environment variables are:"
echo "------------------------------------------------------"
env -0 | sort -z | tr '\0' '\n' | grep -a OIFS_
echo

#---Path to the executable for the SCM. This is the
#---default path for the exe, produced by openifs-test.sh.
#---DP means double precision. To run single precision change
#---DP to SP
export SCM_EXEC="${OIFS_HOME}/build/bin/MASTER_scm.DP"

#---Default assumed paths, only change if you know what you are doing
export SCM_TEST="${OIFS_HOME}/scripts/scm"
export SCM_VERSIONDIR="${OIFS_EXPT}/scm_openifs/48r1"
export SCM_PROJDIR="${SCM_VERSIONDIR}/scm-projects"
export SCM_RUNDIR="${SCM_PROJDIR}/ref48r1"
export SCM_LOGFILE="${SCM_RUNDIR}/scm_run_log.txt"

alias scm_env="env -0 | sort -z | tr '\0' '\n' | grep -a SCM_"

echo -e "\nSCM environment variables are:"
echo "------------------------------------------------------"
env -0 | sort -z | tr '\0' '\n' | grep -a SCM_
echo
