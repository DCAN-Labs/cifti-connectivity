#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MATLAB Runtime environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`dirname "$0"`

# Added 2019-11-26 to ensure that $TMPDIR is nonempty and points to a real path
if [ "x$TMPDIR" = "x" ]; then
    TMPDIR=${exe_dir}/temp
    if [ ! -d $TMPDIR ]; then
        mkdir $TMPDIR;
    fi
fi 

if [ ! -d $TMPDIR/$USER ]; then
    mkdir $TMPDIR/$USER
fi
export MCR_CACHE_ROOT=$TMPDIR/$USER
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=$1
      args="${args} \"${token}\"" 
      shift
  done
  eval "\"${exe_dir}/cifti_conn_matrix_for_wrapper\"" $args
fi
exit

