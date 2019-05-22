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
	echo "--mosaic=[MOSAIC_FILE] - Input mosaic map filename (.fits/casa)"
	echo ""	

	echo "*** OPTIONAL ARGS ***"
	echo "=== SKYMODEL SIM OPTIONS ==="	
	echo "--no-gensources - Turn off point source generation (default=no)"
	echo "--smin=[SMIN] - Min generated compact source flux density in Jy (default: 1.e-4)"
	echo "--smax=[SMAX] - Max generated compact source flux density in Jy (default: 1)"
	echo "--nsources=[NSOURCES] - Number of compact sources generated in mosaic (if >0 override source density options) (default: 0)"
	echo "--sourcedensity=[SOURCE_DENSITY] - Compact source density in sources/deg^2 (default: 1000)"
	echo "--genextsources - Generate extended sources (default=no)"
	echo "--smin-ext=[SMIN_EXT] - Min generated extended source flux density in Jy (default: 1.e-4)"
	echo "--smax-ext=[SMAX_EXT] - Max generated extended source flux density in Jy (default: 1)"
	echo "--smodel=[SMODEL] - Source flux random model (uniform, exp) (default: uniform)"
	echo "--sslope=[SSLOPE] - Source flux random exp model slope (default: 1.6)"
	echo "--nsources-ext=[NSOURCES_EXT] - Number of extended sources generated in mosaic (if >0 override source density options) (default: 0)"
	echo "--sourcedensity-ext=[SOURCE_DENSITY_EXT] - Extended source density in sources/deg^2 (default: 50)"
	echo "--bmaj-min=[BMAJ_MIN] - Min bmaj for extended source generation (default: 4)"
	echo "--bmaj-max=[BMAJ_MAX] - Max bmaj for extended source generation (default: 10)"
	echo "--bmin-min=[BMIN_MIN] - Min bmin for extended source generation (default: 4)"
	echo "--bmin-max=[BMIN_MAX] - Max bmin for extended source generation (default: 10)"
	echo "--pa-min=[BMIN_MIN] - Min pos angle for extended source generation (default: -90)"
	echo "--pa-max=[BMIN_MAX] - Max pos angle extended source generation (default: 90)"
	echo ""

	echo "=== RUN OPTIONS ==="	
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo "--nruns=[RUN_ID] - Number of runs (default: 1)"	
	echo "--startid=[START_ID] - Run start id (default: 1)"	
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"
	echo "--addrunindex - Append a run index to submission script (in case of list execution) (default=no)"	
	echo "--containerrun - Run inside ASKAPSoft container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with ASKAPSoft installed software"
	echo "--containeroptions=[CONTAINER_OPTIONS] - Options to be passed to container run (e.g. -B /home/user:/home/user) (default=none)"	
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
## MANDATORY OPTIONS
MOSAIC_FILE=""

## SKYMODEL GEN DEFAULT OPTIONS
SMIN=1.e-4
SMAX=1
SOURCE_DENSITY=1000
NSOURCES=0
SMODEL="uniform"
SSLOPE="1.6"
GEN_PS_SOURCES=" "
GEN_EXT_SOURCES=" "
SMIN_EXT=1.e-4
SMAX_EXT=1
SOURCE_DENSITY_EXT=1000
NSOURCES_EXT=0
BMAJ_MIN="4"
BMAJ_MAX="10"
BMIN_MIN="4"
BMIN_MAX="10"
PA_MIN="-90"
PA_MAX="90"

## RUN DEFAULT OPTIONS
ENV_FILE=""
NRUNS=1
NRUNS_GIVEN=false
START_ID=1
CONTAINER_IMG=""
RUN_IN_CONTAINER=false
CONTAINER_OPTIONS=""

## SUBMIT DEFAULT OPTIONS
SUBMIT=false
BATCH_QUEUE=""
JOB_WALLTIME="96:00:00"
JOB_MEMORY="4"
JOB_USER_GROUP=""
JOB_USER_GROUP_OPTION=""

##for item in $*
for item in "$@"
do
	case $item in 
		## MANDATORY ##	
		--mosaic=*)
			MOSAIC_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;


		## SKY MODEL SIM ##	
		--no-gensources*)
    	GEN_PS_SOURCES="--no-sources"
    ;;	
		--nsources=*)
			NSOURCES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;
		--sourcedensity=*)
    	SOURCE_DENSITY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--smin=*)
    	SMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--smax=*)
    	SMAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--smodel=*)
    	SMODEL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--sslope=*)
    	SSLOPE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--genextsources*)
    	GEN_EXT_SOURCES="--extsources"
    ;;
		--nsources-ext=*)
			NSOURCES_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
		;;
		--sourcedensity-ext=*)
    	SOURCE_DENSITY_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--smin-ext=*)
    	SMIN_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--smax-ext=*)
    	SMAX_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bmaj-min=*)
    	BMAJ_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bmaj-max=*)
    	BMAJ_MAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bmin-min=*)
    	BMIN_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bmin-max=*)
    	BMIN_MAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--pa-min=*)
    	PA_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--pa-max=*)
    	PA_MAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
			
		## RUN OPTIONS	
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nruns=*)
    	NRUNS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$NRUNS" != "" ]; then
				NRUNS_GIVEN=true
			fi
    ;;
		--startid=*)
    	START_ID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
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
		
		
		## SUBMISSION OPTIONS	
		--submit*)
    	SUBMIT=true
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
		
		
    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done



#######################################
##         CHECK ARGS
#######################################
if [ "$MOSAIC_FILE" = "" ]; then
	echo "ERROR: Missing or empty input mosaic filename!"
	exit 1				
fi



#######################################
##     DEFINE & LOAD ENV VARS
#######################################
export BASEDIR="$PWD"


#######################################
##     RUN
#######################################
RUN_ID=$START_ID
echo "Generate $NRUNS runs starting from RUN_ID=$RUN_ID ..."

for ((index=1; index<=$NRUNS; index=$index+1))
do
	## Create job top directory
	JOB_DIR="$BASEDIR/RUN$RUN_ID"

	skymodel_file="skymodel-RUN$RUN_ID.fits"
	ds9_file="ds9-RUN$RUN_ID.reg"
	source_file="sources-RUN$RUN_ID.dat"
	img_fits="simmap-RUN$RUN_ID.fits"
	sourcepar_file="sourcepars-RUN$RUN_ID.dat"
	extsourcemodel_file="extsourcemap-RUN$RUN_ID.fits"
	
	## Define run script
	shfile="Sim-RUN$RUN_ID.sh"
	echo "INFO: Creating sh file $shfile ..."
	(
		echo "#!/bin/bash"
		echo "#PBS -N SimJob$RUN_ID"
		echo "#PBS -j oe"
		echo "#PBS -o $BASEDIR"
    echo "#PBS -l select=1:ncpus=1:mem=$JOB_MEMORY"'GB'
		echo "#PBS -l walltime=$JOB_WALLTIME"
    echo '#PBS -r n'
    echo '#PBS -S /bin/bash'    
    echo '#PBS -p 1'
		echo "$JOB_USER_GROUP_OPTION"

    echo " "
		echo 'echo "INFO: Running on host $HOSTNAME ..."'
		echo " "

    echo 'echo "*************************************************"'
    echo 'echo "****         PREPARE JOB                     ****"'
    echo 'echo "*************************************************"'
		echo "JOBDIR=$JOB_DIR" 
		echo 'echo "INFO: Creating job top directory $JOBDIR ..."'
		echo 'mkdir -p "$JOBDIR"'
	
		echo " "
	
		echo 'echo "INFO: Entering job directory $JOBDIR ..."'
		echo 'cd $JOBDIR'
		
		echo " "	
		   
		if [ "$ENV_FILE" != "" ]; then
			echo 'echo "INFO: Source the software environment ..."'
			echo "source $ENV_FILE"		
		fi
    
		echo " "
	
		echo 'echo "*************************************************"'
    echo 'echo "****         RUN SKYMODEL SIMULATION         ****"'
    echo 'echo "*************************************************"'
		
		echo " "
		echo 'EXE="'"$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/generate_ps_map.py"'"'
		echo 'EXE_ARGS="'"--no-convolve --mosaic=$MOSAIC_FILE --outfile_img=$JOB_DIR/$img_fits --outfile_model=$JOB_DIR/$skymodel_file --outfile_model_ext=$JOB_DIR/$extsourcemodel_file --outfile_sources=$JOB_DIR/$source_file --outfile_sourcepars=$JOB_DIR/$sourcepar_file --outfile_ds9region=$JOB_DIR/$ds9_file $GEN_PS_SOURCES --source_density=$SOURCE_DENSITY --nsources=$NSOURCES --Smin=$SMIN --Smax=$SMAX $GEN_EXT_SOURCES --Smodel=$SMODEL --Sslope=$SSLOPE --source_density_ext=$SOURCE_DENSITY_EXT --nsources_ext=$NSOURCES_EXT --Smin_ext=$SMIN_EXT --Smax_ext=$SMAX_EXT --bmaj_min=$BMAJ_MIN --bmaj_max=$BMAJ_MAX --bmin_min=$BMIN_MIN --bmin_max=$BMIN_MAX --pa_min=$PA_MIN --pa_max=$PA_MAX "'"'
		
		echo " "

		echo 'echo "INFO: Running command $EXE $EXE_ARGS ..."'
		echo '$EXE $EXE_ARGS'
			
		echo " "
		
    echo 'echo "*** END RUN ***"'

	) > $shfile
	chmod +x $shfile

	####mv $shfile $CURRENTJOBDIR
	(( RUN_ID= $RUN_ID + 1 ))
	
	
	# Submits the job to batch system
	if [ "$SUBMIT" = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
		qsub -q $BATCH_QUEUE $BASEDIR/$shfile
	fi

done

echo "*** END SUBMISSION ***"


