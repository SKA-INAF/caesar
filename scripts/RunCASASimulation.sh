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
	echo "--nruns=[NRUNS] - Number of simulation runs"
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo "--mapsize=[MAP_SIZE] - Map size (in pixels)"
	echo "--pixsize=[PIX_SIZE] - Pixel size (in arcsec)"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"	
	echo "--startid=[START_ID] - Run start id (default: 1)"
	echo "--sourcegenmargin=[SOURCE_GEN_MARGIN_SIZE] - Left/right margin in skymodel map for source generation (default: 0)"
	echo "--bmaj=[BMAJ] - Beam Bmaj of sky model map in arcsec (default: 9.8 arcsec)"
	echo "--bmin=[BMIN] - Beam Bmin of sky model map in arcsec (default: 5.8 arcsec)"
	echo "--bpa=[BMIN] - Beam position angle of sky model map in deg (default: -3 deg)"
	echo "--bkglevel=[BKG_LEVEL] - Bkg level of sky model toy sim map in Jy (default: 10.e-6)"
	echo "--bkgrms=[BKG_RMS] - RMS level of sky model toy sim map in Jy (default: 100.e-6)"
	echo "--sources - Generate compact sources in skymodel (default: true)"
	echo "--no-sources - Do not generate compact sources in skymodel. NB: If not given compact sources are generated."
	echo "--zmin=[ZMIN] - Min generated compact source significance level wrt to bkg level & rms (default: 1)"
	echo "--zmax=[ZMAX] - Max generated compact source significance level wrt to bkg level & rms (default: 10000)"
	##echo "--zmin-model=[ZMIN_MODEL] - Minimum source significance level in sigmas above the bkg below which source data are set to 0 (default: 1)"
	echo "--truncthr=[TRUNC_THRESHOLD] - Flux loss with respect to total at which model is truncated (default: 0.001)"
	echo "--extsources - Generate extended sources in skymodel (default=no)"
	echo "--no-extsources - Do not generate extended sources in skymodel. NB: By default extended sources are not generated."
	echo "--sourcedensity=[SOURCE_DENSITY] - Compact source density in sources/deg^2 (default: 1000)"
	echo "--ext-zmin=[ZMIN_EXT] - Min generated extended source significance level wrt to bkg level & rms (default: 1)"
	echo "--ext-zmax=[ZMAX_EXT] - Max generated extended source significance level wrt to bkg level & rms (default: 5)"
	echo "--ext-sourcedensity=[EXT_SOURCE_DENSITY] - Extended source density in sources/deg^2 (default: 100)"
	echo "--ext-scalemin=[EXT_SCALE_MIN] - Minimum extended source size in arcsec (default: 10)"
	echo "--ext-scalemax=[EXT_SCALE_MAX] - Maximum extended source size in arcsec (default: 100)"
	echo "--ring-wmin=[RING_WIDTH_MIN] - Minimum source ring size in arcsec (default: 5)"
	echo "--ring-wmax=[RING_WIDTH_MAX] - Maximum source ring size in arcsec (default: 20)"
	echo "--simproject=[SIM_PROJECT] - Name of CASA simulation project (default: sim)"
	echo "--visimagename=[VIS_IMAGE_NAME] - Name of CASA sim visibility image name (default: vis.ms)"
	echo "--simtottime=[SIM_TOT_TIME] - Simulation total time in seconds (default: 43200)"
	echo "--telconfigs=[TELCONFIGS] - Antenna configurations (default: [atca_all.cfg])"
	echo "--addnoise - Add noise to CASA simulation (default: no)"
	echo "--maptype=[MAP_TYPE] - Simulated map type (square|hexagonal) (default=square)"
	echo "--frequency=[FREQUENCY_CENTER] - Frequency centroid of simulated data with units (default=2.1GHz)"
	echo "--frequencybw=[FREQUENCY_BANDWIDTH] - Frequency bandwidth of simulated data with units (default=10MHz)"
	echo "--submit - Submit the script to the batch system using queue specified"
	echo "--containerrun - Run inside Caesar container"
	echo "--containerimg=[CONTAINER_IMG] - Singularity container image file (.simg) with CAESAR installed software"
	echo "--queue=[BATCH_QUEUE] - Name of queue in batch system" 
	echo "--jobwalltime=[JOB_WALLTIME] - Job wall time in batch system (default=96:00:00)" 
	echo "--with-graphics - Enable graphics in simobserve. NB: Container run crashes with CASA graphics enabled (default=disabled)" 	
	echo "=========================="
	exit 1
fi



#######################################
##         PARSE ARGS
#######################################
NRUNS=1
NRUNS_GIVEN=false
START_ID=1
SUBMIT=false
BATCH_QUEUE=""
ENV_FILE=""
MAP_SIZE=""
MAP_SIZE_GIVEN=false
PIX_SIZE=""
PIX_SIZE_GIVEN=false
SOURCE_GEN_MARGIN_SIZE=0
GEN_SOURCES=true
GEN_EXT_SOURCES=false
CONTAINER_IMG=""
RUN_IN_CONTAINER=false
BMAJ=9.8
BMIN=5.8
BPA=-3
BKG_LEVEL=10e-6 # Jy
BKG_RMS=100e-6 # Jy
ZMIN=1
ZMAX=10000
ZMIN_EXT=1
ZMAX_EXT=5
##ZMIN_MODEL=1
TRUNC_THRESHOLD=0.001
SOURCE_DENSITY=1000
EXT_SOURCE_DENSITY=100
EXT_SCALE_MIN=10
EXT_SCALE_MAX=100
RING_WIDTH_MIN=5
RING_WIDTH_MAX=20
SIM_PROJECT="" # "sim"
SIM_PROJECT_GIVEN=false
VIS_IMAGE_NAME="" # "vis.ms"
VIS_IMAGE_NAME_GIVEN=false
SIM_TOT_TIME=43200
#TELCONFIGS="['atca_all.cfg']"
TELCONFIGS="'atca_all.cfg'"
ADD_NOISE=false
MAP_TYPE="square"
FREQUENCY="2.1GHz"
FREQUENCYBW="10MHz"
GRAPHICS_FLAG="--no-graphics"
JOB_WALLTIME="96:00:00"

##for item in $*
for item in "$@"
do
	case $item in 
		## MANDATORY ##	
		--nruns=*)
    	NRUNS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$NRUNS" != "" ]; then
				NRUNS_GIVEN=true
			fi
    ;;
		--mapsize=*)
    	MAP_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$MAP_SIZE" != "" ]; then
				MAP_SIZE_GIVEN=true
			fi
    ;;
		--pixsize=*)
    	PIX_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			if [ "$PIX_SIZE" != "" ]; then
				PIX_SIZE_GIVEN=true
			fi
    ;;
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		

		## OPTIONAL ##	
		--startid=*)
    	START_ID=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;	
		--sourcegenmargin=*)
    	SOURCE_GEN_MARGIN_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
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
		--bkglevel=*)
    	BKG_LEVEL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--bkgrms=*)
    	BKG_RMS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;	

		## COMPACT SOURCES OPTIONS
		--sources*)
    	GEN_SOURCES=true
    ;;
		--no-sources*)
    	GEN_SOURCES=false
    ;;
		--zmin=*)
    	ZMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--zmax=*)
    	ZMAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
			
		--truncthr=*)
    	TRUNC_THRESHOLD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--sourcedensity=*)
    	SOURCE_DENSITY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;

		## EXTENDED SOURCES OPTIONS
		--extsources*)
    	GEN_EXT_SOURCES=true
    ;;
		--no-extsources*)
    	GEN_EXT_SOURCES=false
    ;;
		--ext-zmin=*)
    	ZMIN_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--ext-zmax=*)
    	ZMAX_EXT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--ext-sourcedensity=*)
    	EXT_SOURCE_DENSITY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--ext-scalemin=*)
    	EXT_SCALE_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--ext-scalemax=*)
    	EXT_SCALE_MAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--ring-wmin=*)
    	RING_WIDTH_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		--ring-wmax=*)
    	RING_WIDTH_MAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
    ;;
		
		## SUBMISSION OPTIONS
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
		--simproject=*)
    	SIM_PROJECT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			SIM_PROJECT_GIVEN=true
    ;;
		--visimagename=*)
    	VIS_IMAGE_NAME=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			VIS_IMAGE_NAME_GIVEN=true
    ;;
		--simtottime=*)
    	SIM_TOT_TIME=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--telconfigs=*)
    	TELCONFIGS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--addnoise*)
    	ADD_NOISE=true
    ;;
		--maptype=*)
    	MAP_TYPE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--frequency=*)
    	FREQUENCY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--frequencybw=*)
    	FREQUENCYBW=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--with-graphics*)
    	GRAPHICS_FLAG="--graphics"
    ;;

    *)
    # Unknown option
    echo "ERROR: Unknown option ($item)...exit!"
    exit 1
    ;;
	esac
done

# Compute other parameters
CRPIX=`expr $MAP_SIZE / 2`

echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "SUBMIT? $SUBMIT, QUEUE=$BATCH_QUEUE, JOB_WALLTIME: $JOB_WALLTIME"
echo "NRUNS: $NRUNS (START_ID=$START_ID)"
echo "ENV_FILE: $ENV_FILE"
echo "RUN_IN_CONTAINER? $RUN_IN_CONTAINER, CONTAINER_IMG=$CONTAINER_IMG"
echo "MAP_SIZE: $MAP_SIZE (pixels), PIX_SIZE=$PIX_SIZE (arcsec)"
echo "BEAM ($BMAJ arcsec, $BMIN arcsec, $BPA deg)"
echo "CRPIX: $CRPIX"
echo "BKG_LEVEL: $BKG_LEVEL, BKG_RMS: $BKG_RMS"
echo "ZMIN/ZMAX: $ZMIN/$ZMAX"
echo "ZMIN_EXT/ZMAX_EXT: $ZMIN_EXT/$ZMAX_EXT"
echo "EXT_SCALE MIN/MAX: $EXT_SCALE_MIN/$EXT_SCALE_MAX"
echo "RING WIDTH MIN/MAX: $RING_WIDTH_MIN/$RING_WIDTH_MAX"
echo "SOURCE_DENSITY: $SOURCE_DENSITY, EXT_SOURCE_DENSITY: $EXT_SOURCE_DENSITY"
echo "TRUNC_THRESHOLD: $TRUNC_THRESHOLD"
echo "SOURCE_GEN_MARGIN_SIZE: $SOURCE_GEN_MARGIN_SIZE (pixels)"
echo "GEN_SOURCES? $GEN_EXT_SOURCES, GEN_EXT_SOURCES? $GEN_EXT_SOURCES"
echo "SIM_PROJECT: $SIM_PROJECT, VIS: $VIS_IMAGE_NAME"
echo "SIM_TOT_TIME: $SIM_TOT_TIME"
echo "TELCONFIGS: $TELCONFIGS"
echo "ADD_NOISE: $ADD_NOISE"
echo "MAP_TYPE: $MAP_TYPE"
echo "FREQUENCY CENTER/BW: $FREQUENCY/$FREQUENCYBW"
echo "GRAPHICS_FLAG: $GRAPHICS_FLAG"
echo "****************************"
echo ""

## Check arguments parsed
if [ "$NRUNS_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty NRUNS args (hint: you should specify how many simulation to be performed)!"
  exit 1
fi

if [ "$MAP_SIZE_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty MAP_SIZE args (hint: you should specify how big is your skymodel & sim map in pixels)!"
  exit 1
fi

if [ "$PIX_SIZE_GIVEN" = false ] ; then
  echo "ERROR: Missing or empty PIX_SIZE args (hint: you should specify how big is your skymodel & sim pixel size in arcsec)!"
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
GEN_SOURCE_FLAG="--compactsources"
GEN_EXT_SOURCE_FLAG="--extsources"
if [ "$GEN_SOURCES" = false ] ; then
	GEN_SOURCE_FLAG="--no-compactsources"
fi
if [ "$GEN_EXT_SOURCES" = false ] ; then
	GEN_EXT_SOURCE_FLAG="--no-extsources"
fi
ADD_NOISE_FLAG=""
if [ "$ADD_NOISE" = true ] ; then
	ADD_NOISE_FLAG="--addnoise"
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

###################################################
###     GENERATE SUBMISSION SCRIPT
###################################################
RUN_ID=$START_ID
echo "Generate submission scripts (start run numbering from RUN_ID=$RUN_ID)"

for ((index=1; index<=$NRUNS; index=$index+1))
	do

	## Create job top directory
	JOB_DIR="$BASEDIR/RUN$RUN_ID"
	##CASA_SIM_DIR="$JOB_DIR/sim"
	##echo "INFO: Creating job top directory $JOB_DIR ..."
	##mkdir -p "$JOB_DIR"

	
	## Define skymodel simulation files
  simmapfile='simmap-RUN'"$RUN_ID"'.fits'
	skymodelfile='skymodel-RUN'"$RUN_ID"'.fits'
	sourcefile='sources-RUN'"$RUN_ID"'.root'
	ds9regionfile='ds9regions-RUN'"$RUN_ID"'.reg' 
	casaregionfile='casamask-RUN'"$RUN_ID"'.dat'

	## Define CASA simulation image
	##simproject='sim-RUN'"$RUN_ID"
	simproject='sim'
	if [ "$SIM_PROJECT_GIVEN" = true ] ; then
		simproject=$SIM_PROJECT
	fi
	

	## Define CASA visibility image
	##visimg='vis-RUN'"$RUN_ID"'.ms'
	visimg='vis.ms'
	if [ "$VIS_IMAGE_NAME_GIVEN" = true ] ; then
		visimg=$VIS_IMAGE_NAME	
	fi

  echo ""

	## Generate script
	shfile="Sim-RUN$RUN_ID.sh"
	echo "INFO: Creating sh file $shfile ..."
	(
		echo "#!/bin/bash"
		echo "#PBS -N SimJob$RUN_ID"
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
		echo "JOBDIR=$JOB_DIR" 
		echo 'echo "INFO: Creating job top directory $JOBDIR ..."'
		echo 'mkdir -p "$JOBDIR"'
		echo 'echo ""'
	
		echo 'echo "INFO: Entering job directory $JOBDIR ..."'
		echo 'cd $JOBDIR'
    echo 'echo "INFO: Source the software environment ..."'
    echo "source $ENV_FILE"
		echo 'echo ""'
	
		echo " "
    echo " "

    echo 'echo "*************************************************"'
    echo 'echo "****         RUN SKYMODEL SIMULATION         ****"'
    echo 'echo "*************************************************"'
		echo 'echo ""'
    echo 'cd $JOBDIR'

		if [ "$RUN_IN_CONTAINER" = true ] ; then
			echo 'EXE="'"singularity run --app skymodel $CONTAINER_IMG"'"'
		else
			echo 'EXE="'"$CAESAR_SCRIPTS_DIR/skymodel_generator.py"'"'
		fi

		echo 'EXE_ARGS="'"--nx=$MAP_SIZE --ny=$MAP_SIZE --pixsize=$PIX_SIZE --marginx=$SOURCE_GEN_MARGIN_SIZE --marginy=$SOURCE_GEN_MARGIN_SIZE $GEN_SOURCE_FLAG $GEN_EXT_SOURCE_FLAG --bmaj=$BMAJ --bmin=$BMIN --bpa=$BPA --crpix1=$CRPIX --crpix2=$CRPIX --bkg --bkg_level=$BKG_LEVEL --bkg_rms=$BKG_RMS --zmin=$ZMIN --zmax=$ZMAX --zmin_ext=$ZMIN_EXT --zmax_ext=$ZMAX_EXT --source_density=$SOURCE_DENSITY --trunc_thr=$TRUNC_THRESHOLD --ext_source_density=$EXT_SOURCE_DENSITY --ext_scale_min=$EXT_SCALE_MIN --ext_scale_max=$EXT_SCALE_MAX --ring_wmin=$RING_WIDTH_MIN --ring_wmax=$RING_WIDTH_MAX --outputfile=$simmapfile --outputfile_model=$skymodelfile --outputfile_sources=$sourcefile --outputfile_ds9region=$ds9regionfile --outputfile_casaregion=$casaregionfile "'"'

		echo 'echo "Running command $EXE $EXE_ARGS"'
		echo '$EXE $EXE_ARGS'

		echo 'echo ""'

	 	echo " "
    echo " "
    
		echo 'echo "*************************************************"'
    echo 'echo "****         RUN CASA SIMULATION             ****"'
    echo 'echo "*************************************************"'
		echo 'echo ""'
    echo 'cd $JOBDIR'

		if [ "$RUN_IN_CONTAINER" = true ] ; then
			echo 'EXE="'"singularity run --app simulation $CONTAINER_IMG"'"'
		else
			##echo 'EXE="$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/simulate_observation.py"'
			echo 'EXE="'"$CASAPATH/bin/casa --nologger --log2term --nogui -c $CAESAR_SCRIPTS_DIR/simulate_observation.py"'"'
		fi

		echo 'EXE_ARGS="'"--outproject=$simproject --vis=$visimg --skymodel=$skymodelfile --total_time=$SIM_TOT_TIME --telconfigs=$TELCONFIGS $ADD_NOISE_FLAG --maptype=$MAP_TYPE --frequency_center=$FREQUENCY --frequency_bandwidth=$FREQUENCYBW $GRAPHICS_FLAG "'"'
		
		echo 'echo "Running command $EXE $EXE_ARGS"'
		echo '$EXE $EXE_ARGS'
    
    echo 'echo ""'
    
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

