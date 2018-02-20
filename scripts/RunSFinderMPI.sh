#!/bin/bash


#######################################
##         CHECK ENVIRONMENT
#######################################
## Check if CAESAR environment variable exists and is properly set
CURRENTDIR="$PWD"
if [ "$CAESAR_DIR" = "" ]; then
	echo "ERROR: Missing CAESAR_DIR environment variable, please set it to your CAESAR installation path."
	exit 1
fi


#######################################
##         CHECK ARGS
#######################################
NARGS="$#"
echo "INFO: NARGS= $NARGS"

if [ "$NARGS" -lt 2 ]; then
	
	echo "ERROR: Invalid number of arguments...see script usage!"
  echo ""
	echo "**************************"
  echo "***     USAGE          ***"
	echo "**************************"
 	echo "$0 [ARGS]"
	echo ""
	echo "=========================="
	echo "==    ARGUMENT LIST     =="
	echo "=========================="
	echo "*** MANDATORY ARGS ***"
	echo "--config=[CONFIGFILE]"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"	
	echo "--nproc=[NPROC] - Number of MPI processor to run on (default=1)"
	echo "--hostfile=[HOSTFILE] - File listing hostnames to run MPI (default=all available hosts)"
	echo "--containerrun - Run inside Caesar container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with CAESAR installed software"
  echo "****************"
  exit 1
fi

#######################################
##         PARSE ARGS
#######################################
CONFIGFILE=""
CONFIGFILE_GIVEN=false
NPROC=1
HOSTFILE=""
HOSTFILE_GIVEN=false

for item in $*
do
	case $item in 	
    --config=*)
    	CONFIGFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			CONFIGFILE_GIVEN=true
    ;;
		--nproc=*)
      NPROC=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--hostfile=*)
    	HOSTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			HOSTFILE_GIVEN=true
    ;;
		--containerimg=*)
    	CONTAINER_IMG=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--containerrun*)
    	RUN_IN_CONTAINER=true
    ;;

    *)
    # Unknown option
    echo "ERROR: Unknown option...exit!"
    exit 1
    ;;
	esac
done

echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "CONFIGFILE: $CONFIGFILE"
echo "NPROC: $NPROC"
echo "HOSTFILE_GIVEN? $HOSTFILE_GIVEN, HOSTFILE: $HOSTFILE"
echo "RUN_IN_CONTAINER? $RUN_IN_CONTAINER, CONTAINER_IMG=$CONTAINER_IMG"
echo "****************************"
echo ""

if [ "$CONFIGFILE" = "" ] || [ "$CONFIGFILE_GIVEN" = false ]; then
  echo "ERROR: Empty CONFIGFILE argument (hint: you must specify a path to a CAESAR configuration file)!"
  exit 1
fi

#######################################
##         RUN COMMAND
#######################################
#CMD="mpirun -np $NPROC -f $HOSTFILE  $EXE --config=$CONFIGFILE"
CMD="mpirun -np $NPROC "
if [ "$HOSTFILE_GIVEN" = true ] ; then
	CMD="$CMD -f $HOSTFILE "
fi

if [ "$RUN_IN_CONTAINER" = true ] ; then
	EXE="singularity run --app sfinder $CONTAINER_IMG"		
else
	EXE="$CAESAR_DIR/bin/FindSourceMPI"
fi
EXE_ARGS="--config=$CONFIGFILE"

CMD="$CMD $EXE $EXE_ARGS"


echo "INFO: Running cmd: $CMD"
#eval $CMD
exec $CMD &

echo "INFO: End source finder script run."


