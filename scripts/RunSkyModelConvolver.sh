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
	echo "--filelist=[FILELIST] - Ascii file with list of ROOT skymodels"
	echo "--filelist-rec=[FILELIST_REC] - Ascii file with list of FITS rec maps" 
	echo "--inputfile=[FILENAME] - Skymodel root filename. If the --filelist option is given this option is skipped."
	echo "--inputfile-rec=[FILENAME_REC] - FITS rec map filename. If the --filelist-rec option is given this option is skipped."
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"
	echo "--bmaj=[BMAJ] - Convolution beam bmaj in arcsec. NB: Ineffective if filelist-rec or inputfile-rec are given."
	echo "--bmin=[BMIN] - Convolution beam bmin in arcsec. NB: Ineffective if filelist-rec or inputfile-rec are given."
	echo "--bpa=[BMIN] - Convolution beam position angle in degrees. NB: Ineffective if filelist-rec or inputfile-rec are given."
	echo "--truncthreshold=[TRUNC_THRESHOLD] - Flux loss threshold for source truncation (default=0.001) "
	echo "--userthreshold - Use fixed thresholds for all sources provided by options --threshold/--threshold-ext (default=no)"
	echo "--threshold=[THRESHOLD] - Flux threshold in Jy below which pixels are removed from compact sources (default=0)"
	echo "--threshold-ext=[THRESHOLD_EXT] - Flux threshold in Jy below which pixels are removed from extended sources (default=0)"
	echo "--mergesources - Merge overlapping sources (default=no)"
	echo "--mergecompactsources - Merge compact to compact sources (default=no)"
	echo "--mergeextsources - Merge extended to extended or compact sources (default=no)"
	echo "--nsigmas=[NSIGMAS] - Number of gaussian sgmas used in convolution kernel (default=10)"
	echo "--startid=[START_ID] - Run start id (default: 1)"	
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"		
	echo "--submit - Submit the script to the batch system using queue specified"
	echo "--containerrun - Run inside Caesar container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with CAESAR installed software"
	echo "--queue=[BATCH_QUEUE] - Name of queue in batch system"
	echo "--jobwalltime=[JOB_WALLTIME] - Job wall time in batch system (default=96:00:00)"
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
ENV_FILE=""
FILELIST=""
FILELIST_GIVEN=false
FILELIST_REC=""
FILELIST_REC_GIVEN=false
INPUTFILE=""
INPUTFILE_GIVEN=false
INPUTFILE_REC=""
INPUTFILE_REC_GIVEN=false
RECMAP_GIVEN=false
TRUNC_THRESHOLD=0.001
THRESHOLD=0
THRESHOLD_EXT=0
USER_THRESHOLD_OPTION=""
BMAJ=""
BMIN=""
BPA=""
BEAM_GIVEN=false
MERGE_SOURCES=false
MERGE_COMPACT_SOURCES=false
MERGE_EXTENDED_SOURCES=false
MERGE_SOURCES_OPTION=""
MERGE_COMPACT_SOURCES_OPTION=""
MERGE_EXTENDED_SOURCES_OPTION=""
NSIGMAS=10

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
		
		--bmaj=*)
    	BMAJ=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bmin=*)
    	BMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bpa=*)
    	BPA=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;

		## OPTIONAL ##
		--truncthreshold=*)
    	TRUNC_THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--userthreshold*)
    	USER_THRESHOLD_OPTION="--userthreshold"
    ;;
		--threshold=*)
    	THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--threshold-ext=*)
    	THRESHOLD_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--mergesources*)
    	MERGE_SOURCES=true
			MERGE_SOURCES_OPTION="--mergesources"
    ;;
		--mergecompactsources*)
    	MERGE_COMPACT_SOURCES=true
			MERGE_COMPACT_SOURCES_OPTION="--mergecompactsources"
    ;;
		--mergeextsources*)
    	MERGE_EXTENDED_SOURCES=true
			MERGE_EXTENDED_SOURCES_OPTION="--mergeextsources"
    ;;
		--nsigmas=*)
    	NSIGMAS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
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

if [ "$BMAJ" != "" ] && [ "$BMIN" != "" ] && [ "$BPA" != "" ]; then
  BEAM_GIVEN=true
fi
if [ "$INPUTFILE_REC_GIVEN" = true ] || [ "$FILELIST_REC_GIVEN" = true ]; then
  RECMAP_GIVEN=true
fi

if [ "$RECMAP_GIVEN" = false ] && [ "$BEAM_GIVEN" = false ]; then
  echo "ERROR: Missing or empty beam args (hint: you should specify them when no recmap is given)!"
  exit 1
fi

echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "INPUTFILE: $INPUTFILE"
echo "FILELIST: $FILELIST, NMAX_PROCESSED_FILES: $NMAX_PROCESSED_FILES (START_ID=$START_ID)"
echo "INPUTFILE_REC: $INPUTFILE_REC, FILELIST_REC: $FILELIST_REC"
echo "BEAM ($BMAJ,$BMIN,$BPA), BEAM_GIVEN? $BEAM_GIVEN"
echo "USER_THRESHOLD_OPTION: $USER_THRESHOLD_OPTION, THRESHOLD: $THRESHOLD, THRESHOLD_EXT=$THRESHOLD_EXT"
echo "TRUNC_THRESHOLD: $TRUNC_THRESHOLD"
echo "MERGE_SOURCES? $MERGE_SOURCES, MERGE_COMPACT_SOURCES? $MERGE_COMPACT_SOURCES, MERGE_EXTENDED_SOURCES? $MERGE_EXTENDED_SOURCES"
echo "ENV_FILE: $ENV_FILE"
echo "SUBMIT? $SUBMIT, QUEUE=$BATCH_QUEUE, JOB_WALLTIME: $JOB_WALLTIME"
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
			echo "#PBS -N ConvJob$jobindex"			
			echo "#PBS -j oe"
  		echo "#PBS -o $BASEDIR"
			echo "#PBS -l select=1:ncpus=1"
			echo "#PBS -l walltime=$JOB_WALLTIME"
    	echo '#PBS -r n'
      echo '#PBS -S /bin/bash'
      echo '#PBS -p 1'

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


	if [ "$FILELIST_REC_GIVEN" = true ]; then

		while IFS=$'\t' read -r filename filename_rec 
		do

			## Extract base filename from file given in list 
			filename_base=$(basename "$filename")
			file_extension="${filename_base##*.}"
			filename_base_noext="${filename_base%.*}"
			filename_dir=$(dirname "${filename}")
			sim_dir=$(dirname "${filename_dir}")
			echo "INFO: Processing item $filename_base_noext in list ..."

			## Create job top directory
			JOB_DIR="$BASEDIR"

			## Define args
			INPUTFILE_REC_OPTION="--input-rec=$filename_rec"
			BEAM_OPTIONS=""
			outputfile="$filename_base_noext"'_conv-'"RUN$index"'.root'
			outputfile_ds9="ds9regions_$filename_base_noext"'_conv-'"RUN$index"'.reg'
			outputfile_map="$filename_base_noext"'_conv-'"RUN$index"'.fits'
		
		
			## Define executable & args variables and generate script
			shfile="RunConvolver_$filename_base_noext"'-RUN'"$index.sh"
			if [ "$RUN_IN_CONTAINER" = true ] ; then
				EXE="singularity run --app convolver $CONTAINER_IMG"
			else
				EXE="$CAESAR_DIR/bin/SkyModelConvolver"
			fi

			EXE_ARGS="--input=$filename $INPUTFILE_REC_OPTION --output=$outputfile --output-ds9=$outputfile_ds9 --output-map=$outputfile_map $BEAM_OPTIONS --fluxtruncthr=$TRUNC_THRESHOLD $USER_THRESHOLD_OPTION --threshold=$THRESHOLD --threshold-ext=$THRESHOLD_EXT $MERGE_SOURCES_OPTION $MERGE_COMPACT_SOURCES_OPTION $MERGE_EXTENDED_SOURCES_OPTION"

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

	else ## No rec filelist given

		while read filename 
		do

			## Extract base filename from file given in list 
			filename_base=$(basename "$filename")
			file_extension="${filename_base##*.}"
			filename_base_noext="${filename_base%.*}"
			filename_dir=$(dirname "${filename}")
			sim_dir=$(dirname "${filename_dir}")
			echo "INFO: Processing item $filename_base_noext in list ..."

			## Create job top directory
			JOB_DIR="$BASEDIR"

			## Define args
			INPUTFILE_REC_OPTION=""
			if [ "$INPUTFILE_REC_GIVEN" = true ]; then
  			INPUTFILE_REC_OPTION="--input-rec=$INPUTFILE_REC"
			fi

			BEAM_OPTIONS=""
			if [ "$BEAM_GIVEN" = true ]; then
  			BEAM_OPTIONS="--bmaj=$BMAJ --bmin=$BMIN --bpa=$BPA"
			fi

			outputfile="$filename_base_noext"'_conv-'"RUN$index"'.root'
			outputfile_ds9="ds9regions_$filename_base_noext"'_conv-'"RUN$index"'.reg'
			outputfile_map="$filename_base_noext"'_conv-'"RUN$index"'.fits'
		
		
			## Define executable & args variables and generate script
			shfile="RunConvolver_$filename_base_noext"'-RUN'"$index.sh"
			if [ "$RUN_IN_CONTAINER" = true ] ; then
				EXE="singularity run --app convolver $CONTAINER_IMG"
			else
				EXE="$CAESAR_DIR/bin/SkyModelConvolver"
			fi

			EXE_ARGS="--input=$filename $INPUTFILE_REC_OPTION --output=$outputfile --output-ds9=$outputfile_ds9 --output-map=$outputfile_map $BEAM_OPTIONS --fluxtruncthr=$TRUNC_THRESHOLD $USER_THRESHOLD_OPTION --threshold=$THRESHOLD --threshold-ext=$THRESHOLD_EXT $MERGE_SOURCES_OPTION $MERGE_COMPACT_SOURCES_OPTION $MERGE_EXTENDED_SOURCES_OPTION"

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

		done < "$FILELIST"
	
	fi

else
	
	################################################
	##     USE INPUT FILE PROVIDED IN SCRIPT ARGS
	################################################
	## Extract base filename from file given in input 
	filename_base=$(basename "$INPUTFILE")
	file_extension="${filename_base##*.}"
	filename_base_noext="${filename_base%.*}"
	filename_dir=$(dirname "${INPUTFILE}")
	sim_dir=$(dirname "${filename_dir}")
	echo "INFO: Processing item $filename_base_noext in list ..."
	
	JOB_DIR="$BASEDIR"

	## Define args
	INPUTFILE_REC_OPTION=""
	if [ "$INPUTFILE_REC_GIVEN" = true ]; then
  	INPUTFILE_REC_OPTION="--input-rec=$INPUTFILE_REC"
	fi

	BEAM_OPTIONS=""
	if [ "$BEAM_GIVEN" = true ]; then
  	BEAM_OPTIONS="--bmaj=$BMAJ --bmin=$BMIN --bpa=$BPA"
	fi

	outputfile="$filename_base_noext"'_conv.root'
	outputfile_ds9="ds9regions_$filename_base_noext"'_conv.reg'
	outputfile_map="$filename_base_noext"'_conv.fits'

	## Define executable & args variables and generate script
	shfile="RunConvolver_$filename_base_noext"'.sh'

	if [ "$RUN_IN_CONTAINER" = true ] ; then
		EXE="singularity run --app convolver $CONTAINER_IMG"
	else
		EXE="$CAESAR_DIR/bin/SkyModelConvolver"
	fi

	EXE_ARGS="--input=$INPUTFILE $INPUTFILE_REC_OPTION --output=$outputfile --output-ds9=$outputfile_ds9 --output-map=$outputfile_map $BEAM_OPTIONS --fluxtruncthr=$TRUNC_THRESHOLD $USER_THRESHOLD_OPTION --threshold=$THRESHOLD --threshold-ext=$THRESHOLD_EXT $MERGE_SOURCES_OPTION $MERGE_COMPACT_SOURCES_OPTION $MERGE_EXTENDED_SOURCES_OPTION"

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

