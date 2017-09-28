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

EXE="$CAESAR_DIR/bin/FindSourceMPI"


#######################################
##         CHECK ARGS
#######################################
NARGS="$#"
echo "INFO: NARGS= $NARGS"

if [ "$NARGS" -lt 2 ]; then
	echo "ERROR: Invalid number of arguments...see script usage!"
  echo ""
  echo "**** USAGE ****"
 	echo "$0 --nproc=[NPROCESSOR] [--hostfile=[HOSTFILE]] --config=[CONFIGFILE]"
  echo "****************"
  exit 1
fi

#######################################
##         PARSE ARGS
#######################################
CONFIGFILE=""
NPROC=1
HOSTFILE=""
HOSTFILE_GIVEN=false

for item in $*
do
	case $item in 	
    --config=*)
    	CONFIGFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nproc=*)
      NPROC=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--hostfile=*)
    	HOSTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			HOSTFILE_GIVEN=true
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
echo "****************************"
echo ""

#######################################
##         RUN COMMAND
#######################################
#CMD="mpirun -np $NPROC -f $HOSTFILE  $EXE --config=$CONFIGFILE"
CMD="mpirun -np $NPROC "
if [ "$HOSTFILE_GIVEN" = true ] ; then
	CMD="$CMD -f $HOSTFILE "
fi
CMD="$CMD $EXE --config=$CONFIGFILE"

echo "INFO: Running cmd: $CMD"
eval $CMD


echo "INFO: End source finder run."


