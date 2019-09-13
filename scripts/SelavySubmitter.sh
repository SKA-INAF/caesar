#!/bin/bash

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
	echo "--filelist=[FILELIST] - Ascii file with list of image files (.fits/.root, full path) to be processed" 
	echo "--inputfile=[FILENAME] - Input file name to be searched (.fits/.root). If the --filelist option is given this option is skipped."
	echo ""

	echo "*** OPTIONAL ARGS ***"
	echo "=== SFINDER OUTPUT OPTIONS ==="
	echo "--save-regions - Save DS9 regions (default=no)"
	echo "--save-bkgmap - Save bkg map in output ROOT file (default=no)"
	echo "--save-rmsmap - Save rms map in output ROOT file (default=no)"	
	echo "--save-significancemap - Save significance map in output ROOT file (default=no)"
	echo "--save-residualmap - Save residual map in output ROOT file (default=no)"
	echo ""

	echo "=== INPUT OPTIONS ==="
	echo "--split-processing - Distributed processing in multiple tiles (default=no)"
	echo "--nsubx=[NSUBX] - The number of subdivisions in the x-direction when making the subimages (default=1=no tile split)"
	echo "--nsuby=[NSUBY] - The number of subdivisions in the y-direction when making the subimages (default=1=no tile split)"	
	echo "--overlapx=[OVERLAPX] - The number of pixels of overlap between neighbouring subimages in the x-direction (default=no=no overlap)"	
	echo "--overlapy=[OVERLAPY] - The number of pixels of overlap between neighbouring subimages in the y-direction (default=no=no overlap)"
	echo "--readsubimage - Read subimage instead of full image (default=no)"
	echo "--xmin=[XMIN] - Read sub-image of input image starting from pixel x=xmin in 1-based notation (default=0=read full image)"
	echo "--xmax=[XMAX] - Read sub-image of input image up to pixel x=xmax in 1-based notation (default=0=read full image)"
	echo "--ymin=[YMIN] - Read sub-image of input image starting from pixel y=xmin in 1-based notation (default=0=read full image)"
	echo "--ymax=[YMAX] - Read sub-image of input image up to pixel y=ymax in 1-based notation (default=0=read full image)"
	echo ""

	echo "=== SFINDER BKG OPTIONS ==="
	echo "--bmaj=[BMAJ] - User-supplied beam Bmaj in arcsec (NB: used only when beam info is not available in input map) (default: 10 arcsec)"
	echo "--bmin=[BMIN] - User-supplied beam Bmin in arcsec (NB: used only when beam info is not available in input map) (default: 5 arcsec)"
	echo "--bpa=[BMIN] - User-supplied beam position angle in degrees (NB: used only when beam info is not available in input map) (default: 0 deg)"
	echo "--globalbkg - Use global bkg (default=use local bkg)"
	echo "--bkgbox=[BKG_BOXSIZE] - Box size in pixels used to compute local bkg (default=101)"
	echo ""

	echo "=== SFINDER COMPACT SOURCE OPTIONS ==="
	echo "--npixmin=[NPIX_MIN] - Minimum number of pixel to form a compact source (default=5 pixels)"
	echo "--seedthr=[SEED_THR] - Seed threshold (in nsigmas) used in flood-fill (default=5 sigmas)"
	echo "--mergethr=[MERGE_THR] - Merge threshold (in nsigmas) used in flood-fill (default=2.6 sigmas)"
	echo ""

	
	echo "=== SFINDER SOURCE FITTING OPTIONS ==="
	echo "--fitsources - Fit compact point-like sources found (default=false)"
	echo "--fit-maxcomponents=[FIT_MAX_COMPONENTS] - Maximum number of components fitted in a blob (default=3)"	
	
	echo ""

	
	echo "=== SFINDER RUN OPTIONS ==="	
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"
	echo "--addrunindex - Append a run index to submission script (in case of list execution) (default=no)"
	echo "--outdir=[OUTPUT_DIR] - Output directory where to put run output file (default=pwd)"
	echo "--no-logredir - Do not redirect logs to output file in script "
	echo "--no-mpi - Disable MPI run (even with 1 proc) (default=enabled)"
	echo "--mpioptions - Options to be passed to MPI (e.g. --bind-to {none,hwthread, core, l1cache, l2cache, l3cache, socket, numa, board}) (default=)"
	echo "--nproc=[NPROC] - Number of MPI processors per node used (NB: mpi tot nproc=nproc x nnodes) (default=1)"
	echo "--hostfile=[HOSTFILE] - Ascii file with list of hosts used by MPI (default=no hostfile used)"
	echo "--containerrun - Run inside Caesar container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with CAESAR installed software"
	echo "--containeroptions=[CONTAINER_OPTIONS] - Options to be passed to container run (e.g. -B /home/user:/home/user) (default=none)"
	echo ""

	echo "=== SFINDER SUBMISSION OPTIONS ==="
	echo "--submit - Submit the script to the batch system using queue specified"
	echo "--batchsystem - Name of batch system. Valid choices are {PBS,SLURM} (default=PBS)"
	echo "--queue=[BATCH_QUEUE] - Name of queue in batch system" 
	echo "--jobwalltime=[JOB_WALLTIME] - Job wall time in batch system (default=96:00:00)"
	echo "--jobcpus=[JOB_NCPUS] - Number of cpu per node requested for the job (default=1)"
	echo "--jobnodes=[JOB_NNODES] - Number of nodes requested for the job (default=1)"
	echo "--jobmemory=[JOB_MEMORY] - Memory in GB required for the job (default=4)"
	echo "--jobusergroup=[JOB_USER_GROUP] - Name of job user group batch system (default=empty)" 	
	echo "=========================="
  exit 1
fi

#######################################
##         PARSE ARGS
#######################################
## INPUT OPTIONS
FILELIST=""
FILELIST_GIVEN=false
INPUTFILE=""
INPUTFILE_GIVEN=false

SPLIT_PROCESSING=false
NSUBX="1"
NSUBY="1"
OVERLAPX="0"
OVERLAPY="0"
READ_SUBIMAGE="false"
XMIN=0
XMAX=0
YMIN=0
YMAX=0

## OUTPUT OPTIONS
SAVE_DS9REGIONS="false"
SAVE_BKG_MAP="false"
SAVE_RMS_MAP="false"
SAVE_SIGNIFICANCE_MAP="false"
SAVE_RESIDUAL_MAP="false"

## RUN OPTIONS
NMAX_PROCESSED_FILES=-1
ENV_FILE=""
ENV_FILE_GIVEN=false
CONTAINER_IMG=""
CONTAINER_OPTIONS=""
RUN_IN_CONTAINER=false
APPEND_RUN_INDEX=false
MPI_ENABLED=true
NPROC=1
MPI_OPTIONS=""
HOSTFILE=""
HOSTFILE_GIVEN=false
REDIRECT_LOGS=true

## SUBMIT OPTIONS
SUBMIT=false
BATCH_SYSTEM="PBS"
BATCH_QUEUE=""
JOB_WALLTIME="96:00:00"
JOB_MEMORY="4"
JOB_USER_GROUP=""
JOB_USER_GROUP_OPTION=""
JOB_NNODES="1"
JOB_NCPUS="1"
OUTPUT_DIR=$PWD

## BACKGROUND OPTIONS
USE_LOCAL_BKG="true"
BKG_BOXSIZE="101"

## SFINDER OPTIONS
NPIX_MIN="3"
SEED_THR="5"
MERGE_THR="3"

## SOURCE FIT OPTIONS
FIT_SOURCES="false"
FIT_MAX_COMPONENTS="3"



for item in "$@"
do
	case $item in 
		## INPUT OPTIONS ##
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
		
		--split-processing*)
			SPLIT_PROCESSING=true
		;;
    --nsubx=*)
    	NSUBX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		 --nsuby=*)
    	NSUBY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--overlapx=*)
    	OVERLAPX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--overlapy=*)
    	OVERLAPY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;

		--readsubimage*)
			READ_SUBIMAGE="true"
		;;
		--xmin=*)
    	XMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--xmax=*)
    	XMAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ymin=*)
    	YMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ymax=*)
    	YMAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		

		## OUTPUT OPTIONS
		--outdir=*)
    	OUTPUT_DIR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--save-regions*)
    	SAVE_DS9REGIONS="true"
    ;;
		--save-bkgmap*)
    	SAVE_BKG_MAP="true"
    ;;
		--save-rmsmap*)
    	SAVE_RMS_MAP="true"
    ;;
		--save-significancemap*)
    	SAVE_SIGNIFICANCE_MAP="true"
    ;;
		--save-residualmap*)
    	SAVE_RESIDUAL_MAP="true"
    ;;

		## COMPACT SOURCE OPTIONS
		--npixmin=*)
    	NPIX_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		--seedthr=*)
    	SEED_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		--mergethr=*)
    	MERGE_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;


		## RUN OPTIONS
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$ENV_FILE" != "" ]; then
				ENV_FILE_GIVEN=true
			fi
    ;;
		--containerimg=*)
    	CONTAINER_IMG=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--containerrun*)
    	RUN_IN_CONTAINER=true
    ;;
		--containeroptions=*)
    	CONTAINER_OPTIONS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--loglevel=*)
    	LOG_LEVEL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--maxfiles=*)
    	NMAX_PROCESSED_FILES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--addrunindex*)
			APPEND_RUN_INDEX=true
		;;

		
		--no-mpi*)
    	MPI_ENABLED=false
    ;;
		--no-logredir*)
			REDIRECT_LOGS=false
		;;
		--nproc=*)
      NPROC=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--mpioptions=*)
      MPI_OPTIONS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--hostfile=*)
    	HOSTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			HOSTFILE_GIVEN=true
    ;;
		
		## BKG OPTIONS
		--globalbkg*)
    	USE_LOCAL_BKG="false"
    ;;
		--bkgbox=*)
    	BKG_BOXSIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			USE_LOCAL_BKG="true"
    ;;
		
		## SOURCE FITTING OPTIONS
		--fitsources*)
    	FIT_SOURCES="true"
    ;;	
		--fit-maxcomponents=*)
    	FIT_MAX_COMPONENTS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		
		## SUBMISSION OPTIONS
		--submit*)
    	SUBMIT="true"
    ;;
		--batchsystem=*)
    	BATCH_SYSTEM=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--queue=*)
    	BATCH_QUEUE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--jobwalltime=*)
			JOB_WALLTIME=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;	
		--jobcpus=*)
      JOB_NCPUS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--jobnodes=*)
      JOB_NNODES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--jobmemory=*)
			JOB_MEMORY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;
		--jobusergroup=*)
			JOB_USER_GROUP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
			JOB_USER_GROUP_OPTION="#PBS -A $JOB_USER_GROUP"
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

if [ "$BATCH_QUEUE" = "" ] && [ "$SUBMIT" = true ]; then
  echo "ERROR: Empty BATCH_QUEUE argument (hint: you must specify a queue if submit option is activated)!"
  exit 1
fi

if [ "$BATCH_SYSTEM" = "" ] && [ "$SUBMIT" = true ]; then
  echo "ERROR: Empty BATCH_SYSTEM argument (hint: you must specify a batch systen if submit option is activated)!"
  exit 1
fi

if [ "$BATCH_SYSTEM" != "PBS" ] && [ "$BATCH_SYSTEM" != "SLURM" ]; then
  echo "ERROR: Unknown/not supported BATCH_SYSTEM argument $BATCH_SYSTEM (hint: PBS/SLURM are supported)!"
  exit 1
fi


if [ "$CONTAINER_IMG" = "" ] && [ "$RUN_IN_CONTAINER" = true ]; then
  echo "ERROR: Empty CONTAINER_IMG argument (hint: you must specify a container image if run in container option is activated)!"
  exit 1
fi
RUN_IN_CONTAINER_FLAG=""
CONTAINER_IMG_FLAG=""
if [ "$RUN_IN_CONTAINER" = true ]; then
	RUN_IN_CONTAINER_FLAG="--containerrun "
	CONTAINER_IMG_FLAG="--containerimg=$CONTAINER_IMG "	
fi



## Compute total number of MPI processor to be given to mpirun
NPROC_TOT=$(($NPROC * $JOB_NNODES))
MPI_ENABLED_OPTION=""
if [ "$MPI_ENABLED" = false ]; then
  MPI_ENABLED_OPTION="--no-mpi"
fi


## Compute bkg half box size
BOXSIZE=$(((BKG_BOXSIZE-1)/2))

#######################################
##     DEFINE & LOAD ENV VARS
#######################################
export BASEDIR="$PWD"
export OUTPUT_DATADIR="$PWD"
export DATADIR=""

## Load env file
echo "INFO: Loading environment variables defined in file $ENV_FILE ..."
source $ENV_FILE

## Define batch run options
if [ "$BATCH_SYSTEM" = "PBS" ]; then
  BATCH_SUB_CMD="qsub"
	BATCH_QUEUE_NAME_OPTION="-q"
	BATCH_JOB_NAME_DIRECTIVE="#PBS -N"
	BATCH_JOB_OUTFILE_DIRECTIVE="#PBS -o $BASEDIR"
	BATCH_JOB_ERRFILE_DIRECTIVE="#PBS -e $BASEDIR"
	BATCH_JOB_JOINOUTERR_DIRECTIVE="#PBS -j oe"
	BATCH_JOB_WALLTIME_DIRECTIVE="#PBS -l walltime=$JOB_WALLTIME"
	BATCH_JOB_SHELL_DIRECTIVE="#PBS -S /bin/bash"
	BATCH_JOB_USERGRP_DIRECTIVE="#PBS -A $JOB_USER_GROUP"
	BATCH_JOB_PRIORITY="#PBS -p 1"
	BATCH_JOB_NOREQUEUE_DIRECTIVE="#PBS -r n"
	BATCH_JOB_SCATTER_DIRECTIVE="#PBS -l place=scatter"
	BATCH_JOB_NNODES_DIRECTIVE="#PBS -l select=$JOB_NNODES"':'"ncpus=$JOB_NCPUS"':'"mpiprocs=$NPROC"':'"mem=$JOB_MEMORY"'gb'
	BATCH_JOB_NPROC_DIRECTIVE=""
	BATCH_JOB_MEM_DIRECTIVE=""
	BATCH_JOB_NCORE_DIRECTIVE=""

elif [ "$BATCH_SYSTEM" = "SLURM" ]; then
  BATCH_SUB_CMD="sbatch"
	BATCH_QUEUE_NAME_OPTION="-p"
	BATCH_JOB_NAME_DIRECTIVE="#SBATCH -J"
	BATCH_JOB_OUTFILE_DIRECTIVE="#SBATCH -o $BASEDIR"
	BATCH_JOB_ERRFILE_DIRECTIVE="#SBATCH -e $BASEDIR"
	BATCH_JOB_JOINOUTERR_DIRECTIVE="" # There is no such option in SLURM
	BATCH_JOB_WALLTIME_DIRECTIVE="#SBATCH --time=$JOB_WALLTIME"
	BATCH_JOB_SHELL_DIRECTIVE="" # Equivalent SLURM directive not found
	BATCH_JOB_USERGRP_DIRECTIVE="#SBATCH -A $JOB_USER_GROUP"
	BATCH_JOB_PRIORITY="" # Equivalent SLURM directive not found
	BATCH_JOB_NOREQUEUE_DIRECTIVE="#SBATCH --no-requeue"
	BATCH_JOB_SCATTER_DIRECTIVE="#SBATCH --spread-job"
	BATCH_JOB_NNODES_DIRECTIVE="#SBATCH --nodes=$JOB_NNODES"
	BATCH_JOB_NPROC_DIRECTIVE="#SBATCH --ntasks-per-node=$NPROC"
	BATCH_JOB_MEM_DIRECTIVE="#SBATCH --mem=$JOB_MEMORY"'gb'
	BATCH_JOB_NCORE_DIRECTIVE="#SBATCH --ntasks-per-node=$JOB_NCPUS"
else 
	echo "ERROR: Unknown/not supported BATCH_SYSTEM argument $BATCH_SYSTEM (hint: PBS/SLURM are supported)!"
  exit 1
fi





#######################################
##   DEFINE GENERATE CONFIG FCN
#######################################
generate_config(){

	local configfile=$1

	echo "INFO: Creating config file $configfile ..."
	(
  
		echo '############################################'
		echo '##     INPUT OPTIONS'
		echo '############################################'
		echo '# The image to be searched'
		echo "Selavy.image = $inputfile"

		echo '# This is how we divide it up for distributed processing, with the number of subdivisions in each direction, and the size of the overlap region in pixels'
		echo "Selavy.nsubx = $NSUBX"
		echo "Selavy.nsuby = $NSUBY"
		echo "Selavy.overlapx = $OVERLAPX"
		echo "Selavy.overlapy = $OVERLAPY"
		
		echo ""

		echo '# We will search just a subsection of the image'
		echo "Selavy.flagSubsection = $READ_SUBIMAGE"
		echo '# The subsection shows the pixel range for each axis, giving the first & last pixel to be used (1-based)'
		echo "Selavy.subsection = [$XMIN:$XMAX,$YMIN:$YMAX,*,*]"

		echo ""

		echo '############################################'
		echo '##     OUTPUT OPTIONS'
		echo '############################################'
		echo "Selavy.flagDS9 = $SAVE_DS9REGIONS"
		echo "Selavy.ds9File = $ds9region_file"
		echo "Selavy.resultsFile = $catalog_file"
		echo "Selavy.votFile = $vot_file"
		echo "Selavy.writeFitResults = true"

		echo ""

		echo "############################################"
		echo "##     RUN OPTIONS"
		echo "############################################"
		echo "## Detection"
		echo "Selavy.findSpectralTerms = [false, false]"
		echo ""
		echo "# The search threshold, in the flux units of the image. Takes precedence over snrCut"
		echo "#Selavy.threshold = 0.1"
		echo "Selavy.snrCut = $SEED_THR"

		echo ""

		echo "# Aggregation/merge threshold. growthThreshold takes precedence over growthCut"
		echo "Selavy.flagGrowth = true"
		echo "#Selavy.growthThreshold = 0.05"
		echo "Selavy.growthCut = $MERGE_THR"
		
		echo ""
		
		echo "# Size criteria for the final list of detected islands"
		echo "Selavy.minPix = $NPIX_MIN"
		echo "Selavy.minVoxels = 3"
		echo "Selavy.minChannels = 1"

		echo ""
		
		echo "## Save threshold maps"
		echo "Selavy.VariableThreshold = $USE_LOCAL_BKG"
		echo "Selavy.VariableThreshold.boxSize = $BOXSIZE"
		echo "###Selavy.VariableThreshold.ThresholdImageName=detThresh"
		
		if [ "$SAVE_BKG_MAP" = "true" ]; then
			echo "Selavy.VariableThreshold.AverageImageName=$bkgmap_file"
		fi
		if [ "$SAVE_RMS_MAP" = "true" ]; then
			echo "Selavy.VariableThreshold.NoiseImageName=$rmsmap_file"
		fi
		if [ "$SAVE_SIGNIFICANCE_MAP" = "true" ]; then
			echo "Selavy.VariableThreshold.SNRimageName=$snrmap_file"
		fi
		
		echo '###Selavy.Weights.weightsImage = ""'
		echo "###Selavy.Weights.weightsCutoff = 0.15"
		echo "Selavy.Weights.weightsCutoff = -1"

		echo ""

		echo "# Parameters to switch on and control the Gaussian fitting"
		echo "Selavy.Fitter.doFit = $FIT_SOURCES"
		echo "Selavy.Fitter.fitTypes = [full]"
		echo "Selavy.Fitter.numGaussFromGuess = true"
		echo "Selavy.Fitter.maxNumGauss = $FIT_MAX_COMPONENTS"
		echo "Selavy.Fitter.maxReducedChisq = 10."

		echo ""

		echo "# Force the component maps to be casa images for now (casa/fits)"
		echo "Selavy.Fitter.imagetype = fits"

		echo ""
	
		echo "Selavy.threshSpatial = 5"
		echo "Selavy.flagAdjacent = false"

		

		echo "# How the islands are sorted in the final catalogue - by integrated flux in this case"
		echo "#Selavy.sortingParam = iflux"
		echo "Selavy.sortingParam = -pflux"
	
		echo ""

		echo "# Not performing RM Synthesis for this case"
		echo "Selavy.RMSynthesis = false"

		
		
 ) > $configfile

}
## close function

#######################################
##   DEFINE GENERATE EXE SCRIPT FCN
#######################################
generate_run_script(){

	local shfile=$1
	local parsetfile=$2
	local exe="selavy"
	local exe_args="-c $parsetfile"

	echo "INFO: Creating sh file $shfile (exe=$exe, exe_args=$exe_args)..."
	( 
			echo "#!/bin/bash"

			echo ""
			echo "$exe $exe_args"

	) > $shfile

	chmod +x $shfile

}
## close function generate_run_script()


generate_exec_script(){

	local shfile=$1
	local jobindex=$2
	local exe=$3
	local exe_args=$4
	local logfile=$5
	
	echo "INFO: Creating sh file $shfile (jobindex=$jobindex, exe=$exe, exe_args=$exe_args)..."
	( 
			echo "#!/bin/bash"

			echo "$BATCH_JOB_NAME_DIRECTIVE SelavyJob$jobindex"
			echo "$BATCH_JOB_OUTFILE_DIRECTIVE"
			echo "$BATCH_JOB_ERRFILE_DIRECTIVE"
			echo "$BATCH_JOB_JOINOUTERR_DIRECTIVE"
			echo "$BATCH_JOB_WALLTIME_DIRECTIVE"
			echo "$BATCH_JOB_SHELL_DIRECTIVE"
			echo "$BATCH_JOB_USERGRP_DIRECTIVE"
			echo "$BATCH_JOB_PRIORITY"
			echo "$BATCH_JOB_NOREQUEUE_DIRECTIVE"
			echo "$BATCH_JOB_SCATTER_DIRECTIVE"
			echo "$BATCH_JOB_NNODES_DIRECTIVE"
			echo "$BATCH_JOB_NPROC_DIRECTIVE"
			echo "$BATCH_JOB_MEM_DIRECTIVE"
			echo "$BATCH_JOB_NCORE_DIRECTIVE"

      echo " "
      echo " "

      echo 'echo "*************************************************"'
      echo 'echo "****         PREPARE JOB                     ****"'
      echo 'echo "*************************************************"'

      echo 'echo ""'
            
      echo " "

      echo 'echo ""'

			if [ "$ENV_FILE_GIVEN" = true ]; then
				echo 'echo "INFO: Source the software environment ..."'
      	echo "source $ENV_FILE"
			fi
      

      echo 'echo ""'
      
      echo "JOBDIR=$BASEDIR"
     
      echo " "

           
      echo 'echo ""'

      echo " "

      echo " "
      echo 'echo "*************************************************"'
      echo 'echo "****         RUN                             ****"'
      echo 'echo "*************************************************"'
      echo 'echo ""'
      echo '  cd $JOBDIR'
			if [ $REDIRECT_LOGS = true ]; then			
      	echo "  $exe $exe_args >& $logfile"
			else
				echo "  $exe $exe_args"
      fi

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


if [ $FILELIST_GIVEN = true ]; then

	#######################################
	##     LOOP OVER FILES IN LIST
	#######################################
	echo "INFO: Looping over input files listed in file $FILELIST ..."
	cd $BASEDIR
	file_counter=0
	index=1

	while read filename 
	do

		## Extract base filename from file given in list 
		filename_base=$(basename "$filename")
		file_extension="${filename_base##*.}"
		filename_base_noext="${filename_base%.*}"

  	export CURRENTJOBDIR=$BASEDIR  

		## Define input/output filenames
  	inputfile=$filename
		outputfile="out-$filename_base_noext"'.root'
 		ds9region_file="ds9-$filename_base_noext"'.reg'
		ds9fitregion_file="ds9_fitcomp-$filename_base_noext"'.reg'
		catalog_file="catalog-$filename_base_noext"'.txt'
		catalog_fitcomp_file="catalog_fitcomp-$filename_base_noext"'.txt'
		vot_file="catalog-$filename_base_noext"'.xml'
		bkgmap_file="bkg-$filename_base_noext"'.fits'
		rmsmap_file="rms-$filename_base_noext"'.fits'
		snrmap_file="snr-$filename_base_noext"'.fits'

		## Define output log filename
		logfile="output_$filename_base_noext"'.log'

		## Define config & run script file names 
		if [ $APPEND_RUN_INDEX = true ]; then
			configfile="selavy_$filename_base_noext"'_'"$index.parset"
			runfile="Run_$filename_base_noext"'_'"$index.sh"
			shfile="Submit_$filename_base_noext"'_'"$index.sh"
		else
			configfile="selavy_$filename_base_noext"'.parset'
			runfile="Run_$filename_base_noext"'.sh'
			shfile="Submit_$filename_base_noext"'.sh'
		fi

		
		## Generate config file
  	echo "INFO: Creating config file $configfile for input file: $inputfile ..."
		generate_config $configfile

		## Generate run script file (if running on container)
		if [ $RUN_IN_CONTAINER = true ] ; then
			echo "INFO: Generating run script file $runfile ..."
			generate_run_script "$runfile" "$CURRENTJOBDIR/$configfile"
		fi

		## Generate script
		echo "INFO: Creating submit script file $shfile for input file: $inputfile ..."

		CMD=""
		if [ $MPI_ENABLED = true ]; then	
			CMD="mpirun -np $NPROC_TOT $MPI_OPTIONS "
			if [ "$HOSTFILE_GIVEN" = true ] ; then
				CMD="$CMD -f $HOSTFILE "
			fi
		fi
		if [ $RUN_IN_CONTAINER = true ] ; then
			##EXE="$CMD singularity exec $CONTAINER_OPTIONS $CONTAINER_IMG selavy "
			EXE="$CMD singularity exec $CONTAINER_OPTIONS $CONTAINER_IMG $CURRENTJOBDIR/$runfile "
			EXE_ARGS=" "	
		else
			EXE="$CMD selavy "
			EXE_ARGS="-c $BASEDIR/$configfile "
		fi
		
		generate_exec_script "$shfile" "$index" "$EXE" "$EXE_ARGS" "$logfile"


		# Submits the job to batch system
		if [ $SUBMIT = true ] ; then
			echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE using $BATCH_SYSTEM batch system ..."
			$BATCH_SUB_CMD $BATCH_QUEUE_NAME_OPTION $BATCH_QUEUE $CURRENTJOBDIR/$shfile
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

	## Extract base filename from file given in list 
	filename_base=$(basename "$INPUTFILE")
	file_extension="${filename_base##*.}"
	filename_base_noext="${filename_base%.*}"

  export CURRENTJOBDIR=$BASEDIR  

	## Define input/output filenames
  inputfile=$INPUTFILE
	outputfile="out-$filename_base_noext"'.root'
 	ds9region_file="ds9-$filename_base_noext"'.reg'
	ds9fitregion_file="ds9_fitcomp-$filename_base_noext"'.reg'
	catalog_file="catalog-$filename_base_noext"'.txt'
	catalog_fitcomp_file="catalog_fitcomp-$filename_base_noext"'.txt'
	vot_file="catalog-$filename_base_noext"'.xml'	
	bkgmap_file="bkg-$filename_base_noext"'.fits'
	rmsmap_file="rms-$filename_base_noext"'.fits'
	snrmap_file="snr-$filename_base_noext"'.fits'

	## Define output log filename
	logfile="output_$filename_base_noext"'.log'

	## Define and generate config file
	configfile="selavy_$filename_base_noext"'.parset'
  echo "INFO: Creating config file $configfile for input file: $inputfile ..."
	generate_config $configfile

	## Generate run script file (if running on container)
	runfile="Run_$filename_base_noext"'.sh'	
	if [ $RUN_IN_CONTAINER = true ] ; then
		echo "INFO: Generating run script file $runfile ..."
		generate_run_script "$runfile" "$CURRENTJOBDIR/$configfile"
	fi

	## Define executable & args variables and generate script
	shfile="Submit_$filename_base_noext"'.sh'
	
	CMD=""
	if [ $MPI_ENABLED = true ]; then	
		CMD="mpirun -np $NPROC_TOT $MPI_OPTIONS "
		if [ "$HOSTFILE_GIVEN" = true ] ; then
			CMD="$CMD -f $HOSTFILE "
		fi
	fi

	if [ $RUN_IN_CONTAINER = true ] ; then
		##EXE="$CMD singularity exec $CONTAINER_OPTIONS $CONTAINER_IMG selavy "
		EXE="$CMD singularity exec $CONTAINER_OPTIONS $CONTAINER_IMG $CURRENTJOBDIR/$runfile "		
		EXE_ARGS=" "
	else
		EXE="$CMD selavy "
		EXE_ARGS="-c $BASEDIR/$configfile "
	fi
	
	echo "INFO: Creating script file $shfile for input file: $inputfile ..."
	jobId=" "
	generate_exec_script "$shfile" "$jobId" "$EXE" "$EXE_ARGS" "$logfile"

	# Submits the job to batch system
	if [ $SUBMIT = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE using $BATCH_SYSTEM batch system ..."
		$BATCH_SUB_CMD $BATCH_QUEUE_NAME_OPTION $BATCH_QUEUE $CURRENTJOBDIR/$shfile
	fi

fi

echo "*** END SUBMISSION ***"

