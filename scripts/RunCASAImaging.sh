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
	echo "--filelist=[FILELIST] - Ascii file with list of CASA visibility MS data to be imaged" 
	echo "--inputfile=[FILENAME] - Input CASA visibility MS name to be imaged (NB: it's a full path directory not a file). If the --filelist option is given this option is skipped."
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo "--mapsize=[MAP_SIZE] - Map size (in pixels)"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"
	echo "--startid=[START_ID] - Run start id (default: 1)"	
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"
	echo "--outproject=[OUTPUT_PROJECT] - Name of CASA output project (default=rec)"
	echo "--no-mosaic - Disable mosaic generation (default=enabled)"
	echo "--mosaicimgname=[MOSAIC_IMAGE_NAME] - Name of CASA output mosaic image (default=mosaic)"
	echo "--no-fitsout - Disable export of mosaic image to fits (default=enabled)"
	echo "--fitsout=[FITSOUT] - Name of exported mosaic fits (default=output_[FILENAME].fits)"
	echo "--pixsize=[PIX_SIZE] - Pixel size of cleaned map (in arcsec) (default=1arcsec)"
	echo "--phasecenter=[PHASE_CENTER] - Phase center of cleaned fields (default=J2000 254.85067091580214deg -41.47631111052697deg)"
	echo "--deconvolver=[DECONVOLVER] - Deconvolver to be used in clean (hogbom | clark | multiscale | mtmfs | mem | clarkstokes) (default=clark)"
	echo "--gridder=[GRIDDER] - Gridder to be used in clean (standard | wproject | widefield | mosaic | awproject) (default=standard)"
	echo "--weighting=[WEIGHTING] - Weighting model to be used in clean (natural | uniform | briggs | superuniform | radial) (default=briggs)"
	echo "--projection=[PROJECTION] - Astro projection to be used for output cleaned maps (default=SIN)"
	echo "--niter=[NITER] - Max number of iterations to stop clean (default=100000)"
	echo "--cycleniter=[CYCLENITER] - Max number of minor clean cycle iterations (default=10000)"
	echo "--threshold=[THRESHOLD] - Clean residual threshold to stop clean (default=0)"
	echo "--usemask - Use a mask in cleaning (default=no)"
	echo "--mask=[MASK] - Name of mask file (default=search for mask file in simulation dir))"
	echo "--scales=[SCALES] - List of scales (in pixels) for multiscale cleaning (default: [0])"
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
ENV_FILE=""
MAP_SIZE=""
MAP_SIZE_GIVEN=false
PIX_SIZE="1arcsec"
FILELIST_GIVEN=false
INPUTFILE=""
INPUTFILE_GIVEN=false
OUTPUT_PROJECT="rec"
MOSAIC_FLAG=""
MOSAIC_IMAGE_NAME="mosaic"
FITSOUT_FLAG=""
FITSOUT_GIVEN=false
FITSOUT="output.fits"
PHASE_CENTER="J2000 254.85067091580214deg -41.47631111052697deg"
DECONVOLVER="clark"
GRIDDER="standard"
WEIGHTING="briggs"
PROJECTION="SIN"
NITER=100000
CYCLENITER=10000
THRESHOLD=0
MASK=""
USE_MASK=false
MASK_GIVEN=false
##SCALES="['0']"
SCALES="'0'"
JOB_WALLTIME="96:00:00"

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
		--inputfile=*)
    	INPUTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
			if [ "$INPUTFILE" != "" ]; then
				INPUTFILE_GIVEN=true
			fi
    ;;	
		--mapsize=*)
    	MAP_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$MAP_SIZE" != "" ]; then
				MAP_SIZE_GIVEN=true
			fi
    ;;
		
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		

		## OPTIONAL ##
		--startid=*)
    	START_ID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--maxfiles=*)
    	NMAX_PROCESSED_FILES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;	
		--pixsize=*)
    	PIX_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--outproject=*)
    	OUTPUT_PROJECT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--no-mosaic*)
    	MOSAIC_FLAG="--no-mosaic"
    ;;
		--mosaicimgname=*)
    	MOSAIC_IMAGE_NAME=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--no-fitsout*)
    	FITSOUT_FLAG="--no-fitsout"
    ;;
		--fitsout=*)
    	FITSOUT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			FITSOUT_GIVEN=true
    ;;
		--phasecenter=*)
    	PHASE_CENTER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--deconvolver=*)
    	DECONVOLVER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--gridder=*)
    	GRIDDER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--weighting=*)
    	WEIGHTING=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--projection=*)
    	PROJECTION=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--niter=*)
    	NITER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cycleniter=*)
    	CYCLENITER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--threshold=*)
    	THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--usemask*)
    	USE_MASK=true
    ;;
		--mask=*)
    	MASK=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			MASK_GIVEN=true
    ;;
		--scales=*)
    	SCALES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
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



echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "INPUTFILE: $INPUTFILE"
echo "FILELIST: $FILELIST, NMAX_PROCESSED_FILES: $NMAX_PROCESSED_FILES (START_ID=$START_ID)"
echo "PHASE_CENTER: $PHASE_CENTER"
echo "WEIGHTING: $WEIGHTING"
echo "DECONVOLVER: $DECONVOLVER"
echo "GRIDDER: $GRIDDER"
echo "SCALES: $SCALES"
echo "SUBMIT? $SUBMIT, QUEUE=$BATCH_QUEUE, JOB_WALLTIME: $JOB_WALLTIME"
echo "ENV_FILE: $ENV_FILE"
echo "RUN_IN_CONTAINER? $RUN_IN_CONTAINER, CONTAINER_IMG=$CONTAINER_IMG"
echo "MAP_SIZE: $MAP_SIZE (pixels), PIX_SIZE=$PIX_SIZE (arcsec)"
echo "****************************"
echo ""

## Check arguments parsed
if [ "$FILELIST_GIVEN" = false ] && [ "$INPUTFILE_GIVEN" = false ]; then
  echo "ERROR: Missing or empty FILELIST and INPUTFILE args (hint: you should specify at least one)!"
  exit 1
fi

if [ "$MAP_SIZE_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty MAP_SIZE args (hint: you should specify how big is your skymodel & sim map in pixels)!"
  exit 1
fi



if [ "$BATCH_QUEUE" = "" ] && [ "$SUBMIT" = true ]; then
  echo "ERROR: Empty BATCH_QUEUE argument (hint: you must specify a queue if submit option is activated)!"
  exit 1
fi

if [ "$ENV_FILE" = "" ]; then
  echo "ERROR: Empty ENV_FILE arg!"
  exit 1
fi

if [ "$CONTAINER_IMG" = "" ] && [ "$RUN_IN_CONTAINER" = true ]; then
  echo "ERROR: Empty CONTAINER_IMG argument (hint: you must specify a container image if run in container option is activated)!"
  exit 1
fi

## Define flags
MASK_FLAG="$MASK"
if [ "$USE_MASK" = false ] ; then
	MASK_FLAG=""
fi

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
	local recdir=$6

	echo "INFO: Creating sh file $shfile (jobindex=$jobindex, exe=$exe, exe_args=$exe_args)..."
	( 
			echo "#!/bin/bash"
			echo "#PBS -N RecJob$jobindex"			
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
			echo 'echo "INFO: Creating job top directory $JOBDIR ..."'
			echo 'mkdir -p "$JOBDIR"'
			echo 'echo ""'

			echo "RECDIR=$recdir"
			echo 'echo "INFO: Creating job rec directory $RECDIR ..."'
			echo 'mkdir -p "$RECDIR"'
			echo 'echo ""'

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
	#index=1
	index=$START_ID

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
		JOB_DIR="$BASEDIR/Rec_$filename_base_noext-RUN$index"
		CASA_REC_DIR="$JOB_DIR/$OUTPUT_PROJECT"
		#echo "INFO: Creating job top directory $JOB_DIR ..."
		#mkdir -p "$JOB_DIR"
		#echo "INFO: Creating rec directory $CASA_REC_DIR ..."
		#mkdir -p "$CASA_REC_DIR"

		## Define rec output files	
  	recmapfile='recmap_'"$filename_base_noext-RUN$index"'.fits'
		echo "INFO: Set rec output fits file to $recmapfile ..."

		## Define mask file
		maskimg=""
		if [ "$USE_MASK" = true ] ; then
			## If mask is given in arguments use that, otherwise search in sim directory
			if [ "$MASK_GIVEN" = true ] ; then
				maskimg=$MASK
			else
				echo "INFO: Searching for mask file in dirs [$filename_dir,$sim_dir] "
				maskimg=`find $filename_dir $sim_dir -name "casamask*.dat" -type f -print -quit`
			fi
			
			## Check if empty
			if [ "$maskimg" = "" ] ; then
				echo "WARN: No mask found in sim dir or empty map name passed in arguments, no maps will be applied in clean!"
			fi
		fi
		
		## Define executable & args variables and generate script
		shfile="Rec_$filename_base_noext"'-RUN'"$index.sh"
		if [ "$RUN_IN_CONTAINER" = true ] ; then
			EXE="singularity run --app imaging $CONTAINER_IMG"
		else
			EXE="$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/imaging_observation.py"
		fi

		EXE_ARGS="--vis=$filename --mask=$maskimg --mapsize=$MAP_SIZE $MOSAIC_FLAG --outimage=$MOSAIC_IMAGE_NAME --outproject=$OUTPUT_PROJECT $FITSOUT_FLAG --fitsout=$recmapfile --pixsize=$PIX_SIZE --phasecenter='""$PHASE_CENTER""'"" --deconvolver=$DECONVOLVER --gridder=$GRIDDER --weighting=$WEIGHTING --scales=$SCALES --projection=$PROJECTION --niter=$NITER --cycleniter=$CYCLENITER --threshold=$THRESHOLD "
		
		echo "INFO: Creating script file $shfile for input file: $filename_base ..."
		generate_exec_script "$shfile" "$index" "$EXE" "$EXE_ARGS" "$JOB_DIR" "$CASA_REC_DIR"


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

	## Create job top directory
	JOB_DIR="$BASEDIR/Rec_$filename_base_noext"
	CASA_REC_DIR="$JOB_DIR/$OUTPUT_PROJECT"
	#echo "INFO: Creating job top directory $JOB_DIR ..."
	#mkdir -p "$JOB_DIR"
	#echo "INFO: Creating rec directory $CASA_REC_DIR ..."
	#mkdir -p "$CASA_REC_DIR"

	## Define rec output files	
  recmapfile='recmap_'"$filename_base_noext"'.fits'
	echo "INFO: Set rec output fits file to $recmapfile ..."

	## Define mask file
	maskimg=""
	if [ "$USE_MASK" = true ] ; then
		## If mask is given in arguments use that, otherwise search in sim directory
		if [ "$MASK_GIVEN" = true ] ; then
			maskimg=$MASK
		else
			echo "INFO: Searching for mask file in dirs [$filename_dir,$sim_dir] "
			maskimg=`find $filename_dir $sim_dir -name "casamask*.dat" -type f -print -quit`
		fi
			
		## Check if empty
		if [ "$maskimg" = "" ] ; then
			echo "WARN: No mask found in sim dir or empty map name passed in arguments, no maps will be applied in clean!"
		fi
	fi
		
	## Define executable & args variables and generate script
	shfile="Rec_$filename_base_noext"'.sh'

	if [ "$RUN_IN_CONTAINER" = true ] ; then
		EXE="singularity run --app imaging $CONTAINER_IMG"
	else
		EXE="$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/imaging_observation.py"
	fi

	EXE_ARGS="--vis=$INPUTFILE --mask=$maskimg --mapsize=$MAP_SIZE $MOSAIC_FLAG --outimage=$MOSAIC_IMAGE_NAME --outproject=$OUTPUT_PROJECT $FITSOUT_FLAG --fitsout=$recmapfile --pixsize=$PIX_SIZE --phasecenter='""$PHASE_CENTER""'"" --deconvolver=$DECONVOLVER --gridder=$GRIDDER --weighting=$WEIGHTING --scales=$SCALES --projection=$PROJECTION --niter=$NITER --cycleniter=$CYCLENITER --threshold=$THRESHOLD "

	echo "INFO: Creating script file $shfile for input file: $filename_base ..."
	jobId=" "
	generate_exec_script "$shfile" "$jobId" "$EXE" "$EXE_ARGS" "$JOB_DIR" "$CASA_REC_DIR"

	# Submits the job to batch system
	if [ "$SUBMIT" = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
		qsub -q $BATCH_QUEUE $BASEDIR/$shfile
	fi

fi

echo "*** END SUBMISSION ***"



