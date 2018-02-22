#!/bin/bash

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
	echo "--filelist=[FILELIST] - Ascii file with list of true source file in root format"
	echo "--filelist-rec=[FILELIST_REC] - Ascii file with list of FITS rec maps" 
	echo "--inputfile=[FILENAME] - Skymodel root filename. If the --filelist option is given this option is skipped."
	echo "--inputfile-rec=[FILENAME_REC] - FITS rec map filename. If the --filelist-rec option is given this option is skipped."
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"
	echo "=== CORRELATOR OPTIONS ==="
	echo "--overlapthr=[OVERLAP_THRESHOLD] - Fraction of matching pixels to consider source matching (default=0.1) "
	echo "--posthr=[POS_THRESHOLD] - Pixel distance below which two sources are matched (default=2.5) "
	echo "--applyfluxoverlapthr - Apply flux overlap threshold (default: not applied)"
	echo "--fluxoverlapthr=[FLUX_OVERLAP_THRESHOLD] - Flux overlap threshold (default: 0.3)"
	echo "--filtersourcetype - Consider only true sources with given type when searching the match (default=no)"
	echo "--filteredsourcetype - True source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED) (default=-1)"
	echo "--filtersourcesimtype - Consider only true sources with given type when searching the match (default=no)"
	echo "--filteredsourcesimtype -  True source sim types to be crossmatched (eRingLike=1,eBubbleLike=2,eEllipseLike=3,eDiskLike=4,eBlobLike=5) (default=-1)"
	echo "--no-compactsourcecorr - Disable correlation search for compact sources (default=enabled)"
	echo "--no-extendedsourcecorr - Disable correlation search for extended sources (default=enabled)"
	echo "--correctflux - Correct rec integrated flux by beam area (default=no correction)"
	echo ""

	echo "=== RUN OPTIONS ==="
	echo "--startid=[START_ID] - Run start id (default: 1)"		
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"		
	echo "--containerrun - Run inside Caesar container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with CAESAR installed software"
	echo ""

	echo "=== SUBMISSION OPTIONS ==="
	echo "--submit - Submit the script to the batch system using queue specified"
	echo "--queue=[BATCH_QUEUE] - Name of queue in batch system"
	echo "--jobwalltime=[JOB_WALLTIME] - Job wall time in batch system (default=96:00:00)"	
	echo "--jobmemory=[JOB_MEMORY] - Memory in GB required for the job (default=4)"
	echo "--jobusergroup=[JOB_USER_GROUP] - Name of job user group batch system (default=empty)" 
	echo "=========================="
	exit 1
fi

#######################################
##         PARSE ARGS
#######################################
NMAX_PROCESSED_FILES=-1
START_ID=1
SUBMIT=false
BATCH_QUEUE=""
JOB_WALLTIME="96:00:00"
JOB_MEMORY="4"
JOB_USER_GROUP=""
JOB_USER_GROUP_OPTION=""
ENV_FILE=""
FILELIST=""
FILELIST_GIVEN=false
FILELIST_REC=""
FILELIST_REC_GIVEN=false
INPUTFILE=""
INPUTFILE_GIVEN=false
INPUTFILE_REC=""
INPUTFILE_REC_GIVEN=false
APPLY_FLUX_OVERLAP_THR_OPTION=""
FILTER_SOURCE_TYPE_OPTION=""
FILTER_SOURCE_SIM_TYPE_OPTION=""
ENABLE_COMPACT_SOURCE_CORR_OPTION=""
ENABLE_EXTENDED_SOURCE_CORR_OPTION=""
CORRECT_FLUX_OPTION=""
OVERLAP_THRESHOLD="0.1"
POS_THRESHOLD="2.5"
FLUX_OVERLAP_THRESHOLD="0.3"
FILTERED_SOURCE_TYPE="-1"
FILTERED_SOURCE_SIM_TYPE="-1"


for item in "$@"
do
	case $item in 
		## MANDATORY ##	
		--filelist=*)
    	FILELIST=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$FILELIST" != "" ]; then
				FILELIST_GIVEN=true
			fi
    ;;
		--filelist-rec=*)
    	FILELIST_REC=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$FILELIST_REC" != "" ]; then
				FILELIST_REC_GIVEN=true
			fi
    ;;
		--inputfile=*)
    	INPUTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
			if [ "$INPUTFILE" != "" ]; then
				INPUTFILE_GIVEN=true
			fi
    ;;	
		--inputfile-rec=*)
    	INPUTFILE_REC=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
			if [ "$INPUTFILE_REC" != "" ]; then
				INPUTFILE_REC_GIVEN=true
			fi
    ;;	
		
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		

		## OPTIONAL ##
		--applyfluxoverlapthr*)
    	APPLY_FLUX_OVERLAP_THR_OPTION="--applyFluxOverlapThr"
    ;;
		--filtersourcetype*)
    	FILTER_SOURCE_TYPE_OPTION="--filterByType"
    ;;
		--filtersourcesimtype*)
    	FILTER_SOURCE_SIM_TYPE_OPTION="--filterBySimType"
    ;;
		--no-compactsourcecorr*)
    	ENABLE_COMPACT_SOURCE_CORR_OPTION="--no-compactSourceCorrelation"
    ;;
		--no-extendedsourcecorr*)
    	ENABLE_EXTENDED_SOURCE_CORR_OPTION="--no-extendedSourceCorrelation"
    ;;
		--correctflux*)
    	CORRECT_FLUX_OPTION="--correctFlux"
    ;;

		--overlapthr=*)
    	OVERLAP_THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--posthr=*)
    	POS_THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fluxoverlapthr=*)
    	FLUX_OVERLAP_THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--filteredsourcetype=*)
    	FILTERED_SOURCE_TYPE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--filteredsourcesimtype=*)
    	FILTERED_SOURCE_SIM_TYPE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		--startid=*)
    	START_ID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--maxfiles=*)
    	NMAX_PROCESSED_FILES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;	
		
		--queue=*)
    	BATCH_QUEUE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--jobwalltime=*)
			JOB_WALLTIME=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;
		--jobmemory=*)
			JOB_MEMORY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;
		--jobusergroup=*)
			JOB_USER_GROUP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
			JOB_USER_GROUP_OPTION="#PBS -A $JOB_USER_GROUP"
		;;
		--submit*)
    	SUBMIT=true
    ;;
		--containerimg=*)
    	CONTAINER_IMG=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--containerrun*)
    	RUN_IN_CONTAINER=true
    ;;

    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done


## Check arguments parsed
if [ "$FILELIST_GIVEN" = false ] && [ "$INPUTFILE_GIVEN" = false ]; then
  echo "ERROR: Missing or empty FILELIST and INPUTFILE args (hint: you should specify at least one)!"
  exit 1
fi

if [ "$FILELIST_REC_GIVEN" = false ] && [ "$INPUTFILE_REC_GIVEN" = false ]; then
  echo "ERROR: Missing or empty FILELIST_REC and INPUTFILE_REC args (hint: you should specify at least one)!"
  exit 1
fi

echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "INPUTFILE: $INPUTFILE"
echo "FILELIST: $FILELIST, NMAX_PROCESSED_FILES: $NMAX_PROCESSED_FILES (START_ID=$START_ID)"
echo "INPUTFILE_REC: $INPUTFILE_REC, FILELIST_REC: $FILELIST_REC"
echo "APPLY_FLUX_OVERLAP_THR_OPTION: $APPLY_FLUX_OVERLAP_THR_OPTION"
echo "FILTER_SOURCE_TYPE_OPTION: $FILTER_SOURCE_TYPE_OPTION"
echo "FILTER_SOURCE_SIM_TYPE_OPTION: $FILTER_SOURCE_SIM_TYPE_OPTION"
echo "ENABLE_COMPACT_SOURCE_CORR_OPTION: $ENABLE_COMPACT_SOURCE_CORR_OPTION"
echo "ENABLE_EXTENDED_SOURCE_CORR_OPTION: $ENABLE_EXTENDED_SOURCE_CORR_OPTION"
echo "CORRECT_FLUX_OPTION: $CORRECT_FLUX_OPTION"
echo "OVERLAP_THRESHOLD: $OVERLAP_THRESHOLD"
echo "POS_THRESHOLD: $POS_THRESHOLD"
echo "FLUX_OVERLAP_THRESHOLD: $FLUX_OVERLAP_THRESHOLD"
echo "FILTERED_SOURCE_TYPE: $FILTERED_SOURCE_TYPE"
echo "FILTERED_SOURCE_SIM_TYPE: $FILTERED_SOURCE_SIM_TYPE"
echo "ENV_FILE: $ENV_FILE"
echo "SUBMIT? $SUBMIT, QUEUE=$BATCH_QUEUE, JOB_WALLTIME: $JOB_WALLTIME, JOB_MEMORY: $JOB_MEMORY, JOB_USER_GROUP: $JOB_USER_GROUP"
echo "RUN_IN_CONTAINER? $RUN_IN_CONTAINER, CONTAINER_IMG=$CONTAINER_IMG"
echo "****************************"
echo ""


#######################################
##     DEFINE & LOAD ENV VARS
#######################################
export BASEDIR="$PWD"
export OUTPUT_DATADIR="$PWD"
export DATADIR=""

## Load env file
echo "INFO: Loading environment variables defined in file $ENV_FILE ..."
source "$ENV_FILE"

#######################################
##   DEFINE GENERATE EXE SCRIPT FCN
#######################################
generate_exec_script(){

	local shfile=$1
	local jobindex=$2
	local exe=$3
	local exe_args=$4
	local jobdir=$5

	echo "INFO: Creating sh file $shfile (jobindex=$jobindex, exe=$exe, exe_args=$exe_args)..."
	( 
			echo "#!/bin/bash"
			echo "#PBS -N CorrJob$jobindex"			
			echo "#PBS -j oe"
  		echo "#PBS -o $BASEDIR"
			echo "#PBS -l select=1:ncpus=1:mem=$JOB_MEMORY"'GB'
			echo "#PBS -l walltime=$JOB_WALLTIME"
    	echo '#PBS -r n'
      echo '#PBS -S /bin/bash'
      echo '#PBS -p 1'
			echo "$JOB_USER_GROUP_OPTION"

      echo " "
      echo " "
			echo 'echo "INFO: Running on host $HOSTNAME ..."'

      echo 'echo "*************************************************"'
      echo 'echo "****         PREPARE JOB                     ****"'
      echo 'echo "*************************************************"'

      echo 'echo ""'
            
      echo " "

      echo 'echo ""'
      echo 'echo "INFO: Source the software environment ..."'
      echo "source $ENV_FILE"

      echo 'echo ""'
      
			echo "JOBDIR=$jobdir"
			echo 'echo "INFO: Entering job directory $JOBDIR ..."'
			echo 'cd $JOBDIR'
	
     	echo ""
      echo " "

           
      echo 'echo ""'

      echo " "

      echo " "
      echo 'echo "*************************************************"'
      echo 'echo "****         RUN SIMULATION                  ****"'
      echo 'echo "*************************************************"'
      echo 'echo ""'
      echo '  cd $JOBDIR'

      echo "  $exe $exe_args"
      
      echo '  echo ""'

      echo " "
      echo " "
      
      echo 'echo "*** END RUN ***"'

 	) > $shfile

	chmod +x $shfile
}
## close function generate_exec_script()


#######################################
##   GENERATE AND SUBMIT SCRIPT JOBS
#######################################


if [ "$FILELIST_GIVEN" = true ]; then

	#######################################
	##     LOOP OVER FILES IN LIST
	#######################################
	echo "INFO: Looping over input files listed in file $FILELIST ..."
	cd $BASEDIR
	file_counter=0
	index=$START_ID

	while IFS=$'\t' read -r filename filename_rec 
	do

		## Extract base filename from file given in list 
		filename_base=$(basename "$filename")
		file_extension="${filename_base##*.}"
		filename_base_noext="${filename_base%.*}"
		filename_dir=$(dirname "${filename}")
		echo "INFO: Processing item $filename_base_noext in list ..."

		## Create job top directory
		JOB_DIR="$BASEDIR"

		## Define args
		INPUTFILE_REC_OPTION="--input-rec=$filename_rec"
		##outputfile="CorrOut_$filename_base_noext"'-'"RUN$index"'.root'
		outputfile="CorrOut_$filename_base_noext"'.root'
			
		## Define executable & args variables and generate script
		#shfile="RunSourceCorrelator_$filename_base_noext"'-RUN'"$index.sh"
		shfile="RunSourceCorrelator_$filename_base_noext.sh"
		if [ "$RUN_IN_CONTAINER" = true ] ; then
			EXE="singularity run --app catalogcorr $CONTAINER_IMG"
		else
			EXE="$CAESAR_DIR/bin/CorrelateSourceCatalogs"
		fi

		EXE_ARGS="--input=$filename $INPUTFILE_REC_OPTION --output=$outputfile $ENABLE_COMPACT_SOURCE_CORR_OPTION $ENABLE_EXTENDED_SOURCE_CORR_OPTION $FILTER_SOURCE_TYPE_OPTION $FILTER_SOURCE_SIM_TYPE_OPTION --selectedType=$FILTERED_SOURCE_TYPE --selectedSimType=$FILTERED_SOURCE_SIM_TYPE $CORRECT_FLUX_OPTION $APPLY_FLUX_OVERLAP_THR_OPTION --overlapThr=$OVERLAP_THRESHOLD --posThr=$POS_THRESHOLD --fluxOverlapThr=$FLUX_OVERLAP_THRESHOLD "

		echo "INFO: Creating script file $shfile for input file: $filename_base ..."
		generate_exec_script "$shfile" "$index" "$EXE" "$EXE_ARGS" "$JOB_DIR"

		# Submits the job to batch system
		if [ "$SUBMIT" = true ] ; then
			echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
			qsub -q $BATCH_QUEUE $BASEDIR/$shfile
		fi

		(( file_counter= $file_counter + 1 ))
		(( index= $index + 1 ))

		## If total number of jobs exceeds the maximum stop everything!
		if [ "$NMAX_PROCESSED_FILES" != "-1" ] && [ $file_counter -ge $NMAX_PROCESSED_FILES ]; then
  		echo "INFO: Maximum number of processed files ($NMAX_PROCESSED_FILES) reached, exit loop..."
  		break;
		fi

	done < <(paste $FILELIST $FILELIST_REC)


else
	
	################################################
	##     USE INPUT FILE PROVIDED IN SCRIPT ARGS
	################################################
	## Extract base filename from file given in input 
	filename_base=$(basename "$INPUTFILE")
	file_extension="${filename_base##*.}"
	filename_base_noext="${filename_base%.*}"
	filename_dir=$(dirname "${INPUTFILE}")
	echo "INFO: Processing item $filename_base_noext in list ..."
	
	JOB_DIR="$BASEDIR"

	## Define args
	INPUTFILE_REC_OPTION=""
	if [ "$INPUTFILE_REC_GIVEN" = true ]; then
  	INPUTFILE_REC_OPTION="--input-rec=$INPUTFILE_REC"
	fi

	outputfile="CorrOut_$filename_base_noext"'.root'
	
	## Define executable & args variables and generate script
	shfile="RunSourceCorrelator_$filename_base_noext"'.sh'

	if [ "$RUN_IN_CONTAINER" = true ] ; then
		EXE="singularity run --app catalogcorr $CONTAINER_IMG"
	else
		EXE="$CAESAR_DIR/bin/CorrelateSourceCatalogs"
	fi

	EXE_ARGS="--input=$INPUTFILE $INPUTFILE_REC_OPTION --output=$outputfile $ENABLE_COMPACT_SOURCE_CORR_OPTION $ENABLE_EXTENDED_SOURCE_CORR_OPTION $FILTER_SOURCE_TYPE_OPTION $FILTER_SOURCE_SIM_TYPE_OPTION --selectedType=$FILTERED_SOURCE_TYPE --selectedSimType=$FILTERED_SOURCE_SIM_TYPE $CORRECT_FLUX_OPTION $APPLY_FLUX_OVERLAP_THR_OPTION --overlapThr=$OVERLAP_THRESHOLD --posThr=$POS_THRESHOLD --fluxOverlapThr=$FLUX_OVERLAP_THRESHOLD"

	echo "INFO: Creating script file $shfile for input file: $filename_base ..."
	jobId=" "
	generate_exec_script "$shfile" "$jobId" "$EXE" "$EXE_ARGS" "$JOB_DIR"

	# Submits the job to batch system
	if [ "$SUBMIT" = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
		qsub -q $BATCH_QUEUE $BASEDIR/$shfile
	fi

fi

echo "*** END SUBMISSION ***"

