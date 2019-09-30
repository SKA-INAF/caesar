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
	echo ""

	echo "=== SFINDER BKG OPTIONS ==="
	echo "--bmaj=[BMAJ] - User-supplied beam Bmaj in degrees (NB: used only when beam info is not available in input map) (default: 10 arcsec)"
	echo "--bmin=[BMIN] - User-supplied beam Bmin in degrees (NB: used only when beam info is not available in input map) (default: 5 arcsec)"
	echo "--bpa=[BMIN] - User-supplied beam position angle in degrees (NB: used only when beam info is not available in input map) (default: 0 deg)"
	#echo "--bkgbox=[BKG_BOXSIZE] - The [x,y] size of the grid to use. Default = ~4* beam size square (default=101)"
	#echo "--bkggridsize=[BKG_BOX_GRID_SIZE] - The [x,y] size of the box over which the rms/bkg is calculated. Default = 5*grid (default=101)"
	echo ""

	echo "=== SFINDER COMPACT SOURCE OPTIONS ==="
	echo "--seedthr=[SEED_THR] - Seed threshold (in nsigmas) used in flood-fill (default=5 sigmas)"
	echo "--mergethr=[MERGE_THR] - Merge threshold (in nsigmas) used in flood-fill (default=3 sigmas)"
	echo ""
	
	echo "=== SFINDER SOURCE FITTING OPTIONS ==="
	echo "--fit-maxcomponents=[FIT_MAX_COMPONENTS] - Maximum number of components fitted in a blob (default=5)"	
	echo ""
	
	echo "=== SFINDER RUN OPTIONS ==="	
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"
	echo "--addrunindex - Append a run index to submission script (in case of list execution) (default=no)"
	echo "--outdir=[OUTPUT_DIR] - Output directory where to put run output file (default=pwd)"
	echo "--no-logredir - Do not redirect logs to output file in script "
	echo "--ncores=[NCORES] - Number of processing cores used (default=1)"
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


## OUTPUT OPTIONS
SAVE_DS9REGIONS="false"
SAVE_BKG_MAP="false"

## RUN OPTIONS
NMAX_PROCESSED_FILES=-1
ENV_FILE=""
ENV_FILE_GIVEN=false
APPEND_RUN_INDEX=false
NCORES=1
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
BKG_BOX_SIZE="101"
BKG_GRID_SIZE="101"

## SFINDER OPTIONS
SEED_THR="5"
MERGE_THR="3"

## SOURCE FIT OPTIONS
FIT_MAX_COMPONENTS="5"



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

		## COMPACT SOURCE OPTIONS
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
		
		--loglevel=*)
    	LOG_LEVEL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--maxfiles=*)
    	NMAX_PROCESSED_FILES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--addrunindex*)
			APPEND_RUN_INDEX=true
		;;

		
		--no-logredir*)
			REDIRECT_LOGS=false
		;;
		--ncores=*)
      NCORES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		## BKG OPTIONS
		--bkgbox=*)
    	BKG_BOX_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--bkggrid=*)
    	BKG_GRID_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		## SOURCE FITTING OPTIONS
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









generate_run_script(){

	local shfile=$1
	local jobindex=$2
	local filename=$3
	
	## Extract base filename
	filename_base=$(basename "$filename")
	file_extension="${filename_base##*.}"
	filename_base_noext="${filename_base%.*}"

	## Set RMS & background map filenames
	rms_file="$filename_base_noext"'_rms.fits'
	bkg_file="$filename_base_noext"'_bkg.fits'
	
	## Set catalog filename
	catalog_file="catalog-$filename_base_noext"'.dat'
	catalog_tab_file="catalog-$filename_base_noext"'.tab'

	## Set DS9 region filename
	ds9_file="ds9-$filename_base_noext"'.reg'

	## Set logfile
	logfile="output_$filename_base_noext"'.log'
	
	echo "INFO: Creating sh file $shfile (jobindex=$jobindex)..."
	( 
			echo "#!/bin/bash"

			echo "$BATCH_JOB_NAME_DIRECTIVE AegeanJob$jobindex"
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
			
			echo 'echo "INFO: Copying data file to job dir ..."'
			echo "cp $filename $BASEDIR"
      echo " "
			
			echo "python << 'EOF'"
			
			echo " "			
			echo 'import sys'
			echo 'import bdsf'
			echo " "

			echo "input_image='""$filename_base""'"
      echo "img = bdsf.process_image(input_image,adaptive_rms_box=True,output_all=True,thresh_pix=$SEED_THR,thresh_isl=$MERGE_THR,ncores=$NCORES)"
			echo " "

			echo 'print("INFO: Write the gaussian source component list... ")'
			echo "img.write_catalog(format='ascii', catalog_type='gaul')"
			echo "img.write_catalog(format='ds9', catalog_type='gaul')"
			echo " "

			echo 'print("INFO: Write the source list... ")'
			echo "img.write_catalog(format='ascii', catalog_type='srl')"
			echo "img.write_catalog(format='ds9', catalog_type='srl')"
			echo " "

			echo 'print("INFO: Export image... ")'
			echo "img.export_image(img_type='gaus_resid')"
			echo " "

			echo 'print("INFO: Export image... ")'
			echo "img.export_image(img_type='gaus_model', outfile=input_image+'.model')"
			echo " "
	
			echo 'EOF'

      echo " "
      echo " "
				
      #if [ $REDIRECT_LOGS = true ]; then			
      #	echo "  $exe $exe_args >& $logfile"
			#else
			#	echo "  $exe $exe_args"
      #fi

      echo '  echo ""'

      echo " "
      echo " "
      
      echo 'echo "*** END RUN ***"'

 	) > $shfile

	chmod +x $shfile
}
## close function generate_run_script()








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


		## Define run script file names 
		if [ $APPEND_RUN_INDEX = true ]; then
			shfile="Run_$filename_base_noext"'_'"$index.sh"
		else
			shfile="Run_$filename_base_noext"'.sh'
		fi

		generate_run_script "$shfile" "$index" "$filename"

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

	## Generate run script
	echo "INFO: Creating script file $shfile for input file: $inputfile ..."
	jobId=" "
	shfile="Run_$filename_base_noext"'.sh'
	generate_run_script "$shfile" "$jobId" "$INPUTFILE"

	# Submits the job to batch system
	if [ $SUBMIT = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE using $BATCH_SYSTEM batch system ..."
		$BATCH_SUB_CMD $BATCH_QUEUE_NAME_OPTION $BATCH_QUEUE $CURRENTJOBDIR/$shfile
	fi

fi

echo "*** END SUBMISSION ***"

