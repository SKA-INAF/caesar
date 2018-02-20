#!/bin/bash

NARGS="$#"

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
	echo "--start=[STARTJOBID] - Killed started job id"
	echo "--end=[ENDJOBID] - Killed end job id"
	echo "=========================="
	exit 1
fi

STARTJOBID=""
STARTJOBID_GIVEN=false
ENDJOBID=""
ENDJOBID_GIVEN=false
STEP=1

for item in "$@"
do
	case $item in 
		## MANDATORY ##	
		--start=*)
    	STARTJOBID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$STARTJOBID" != "" ]; then
				STARTJOBID_GIVEN=true
			fi
    ;;
		--end=*)
    	ENDJOBID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
			if [ "$ENDJOBID" != "" ]; then
				ENDJOBID_GIVEN=true
			fi
    ;;	
		
    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done

if [ "$STARTJOBID_GIVEN" = false ] || [ "$ENDJOBID_GIVEN" = false ]; then
  echo "ERROR: Missing or empty STARTJOBID/ENDJOBID args (hint: you should specify them)!"
  exit 1
fi

echo "INFO: Killing jobs from id $STARTJOBID to $ENDJOBID ..."

for ((jobid=$STARTJOBID; jobid<=$ENDJOBID; jobid++))
do
	echo "INFO: Killing job $jobid ..."
  qdel $jobid

done

echo "### END RUN ###"
