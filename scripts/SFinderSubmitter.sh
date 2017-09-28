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
	echo "--envfile=[ENV_FILE] - File (.sh) with list of environment variables to be loaded by each processing node"
	echo ""
	echo ""
	echo "*** OPTIONAL ARGS ***"	
	echo "--submit - Submit the script to the batch system using queue specified"
	echo "--queue=[BATCH_QUEUE] - Name of queue in batch system" 
	echo "--loglevel=[LOG_LEVEL] - Logging level string {INFO, DEBUG, WARN, ERROR, OFF} (default=INFO)"
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"
	echo "--outdir=[OUTPUT_DIR] - Output directory where to put run output file (default=pwd)"
	echo "--nproc=[NPROC] - Number of processors to be used in MPI run (default=1)"
	echo "--hostfile=[HOSTFILE] - Ascii file with list of hosts used by MPI (default=no hostfile used)"
	echo "--nthreads=[NTHREADS] - Number of threads to be used in OpenMP (default=-1=all available)"
	echo "--tilesize=[TILE_SIZE] - Size (in pixels) of tile used to partition input image in distributed processing (default=0=no tile split)"		
	echo "--tilestep=[TILE_STEP] - Tile step size (range 0-1) expressed as tile fraction used in tile overlap (default=1=no overlap)"
	echo "--xmin=[XMIN] - Read sub-image of input image starting from pixel x=xmin (default=0=read full image)"
	echo "--xmax=[XMAX] - Read sub-image of input image up to pixel x=xmax (default=0=read full image)"
	echo "--ymin=[YMIN] - Read sub-image of input image starting from pixel y=xmin (default=0=read full image)"
	echo "--ymax=[YMAX] - Read sub-image of input image up to pixel y=ymax (default=0=read full image)"	
	echo "--globalbkg - Use global bkg (default=use local bkg)"
	echo "--bkgbox=[BKG_BOXSIZE] - Box size (muliple of beam size) used to compute local bkg (default=20 x beam)"
	echo "--bkggrid=[BKG_GRIDSIZE] - Grid size (fraction of bkg box) used to compute local bkg (default=0.2 x box)"
	echo "--no-compactsearch - Do not search compact sources"
	echo "--npixmin=[NPIX_MIN] - Minimum number of pixel to consider a source (default=5 pixels)"
	echo "--seedthr=[SEED_THR] - Seed threshold (in nsigmas) used in flood-fill (default=5 sigmas)"
	echo "--mergethr=[MERGE_THR] - Merge threshold (in nsigmas) used in flood-fill (default=2.6 sigmas)"
	echo "--brightseedthr=[BRIGHT_SEED_THR] - Seed threshold (in nsigmas) used in flood-fill for bright source removal (default=10 sigmas)"
	echo "--no-extendedsearch - Do not search extended sources"
	echo "--extsfinder=[EXT_SFINDER_METHOD] - Extended source search method {1=WT-thresholding,2=SPSegmentation,3=ActiveContour,4=Saliency thresholding} (default=3)"	
	echo "--activecontour=[AC_METHOD] - Active contour method {1=Chanvese, 2=LRAC} (default=2)"
	echo "--no-nestedsearch - Do not search nested sources (default=search)"		
	echo "--selectsources - Apply selection to compact sources found (default=false)"
	echo "--fitsources - Fit compact point-like sources found (default=false)"	
  echo "--spsize - Superpixel size (in pixels) used in hierarchical clustering (default=20)"	
	echo "--spbeta - Superpixel regularization par (beta) used in hierarchical clustering (default=1)"	
	echo "--spminarea - Superpixel min area (in pixels) used in hierarchical clustering (default=10)"	
	echo "--saliencythr - Saliency map threshold factor wrt optimal/median threshold (default=2.8)"
	echo "--saliencyminreso - Superpixel size (in pixels) used in multi-reso saliency map smallest scale (default=20 pixels)"
	echo "--saliencymaxreso - Superpixel size (in pixels) used in multi-reso saliency map highest scale (default=60 pixels)"
	echo "--saliencyresostep - Superpixel size step (in pixels) used in multi-reso saliency map computation (default=10 pixels)"
	echo "--saliencynn - Fraction of most similar region neighbors used in saliency map computation (default=1)"
	echo "=========================="
  exit 1
fi

#######################################
##         PARSE ARGS
#######################################
SUBMIT=false
FILELIST_GIVEN=false
INPUTFILE=""
INPUTFILE_GIVEN=false
NPROC=1
HOSTFILE=""
HOSTFILE_GIVEN=false
TILE_SIZE=0
SPLIT_IN_TILES="false"
TILE_OVERLAP="false"
TILE_STEP=1
FILELIST=""
NMAX_PROCESSED_FILES=-1
BATCH_QUEUE=""
LOG_LEVEL="INFO"
OUTPUT_DIR=$PWD
NTHREADS=1
READ_TILE="false"
XMIN=0
XMAX=0
YMIN=0
YMAX=0
USE_LOCAL_BKG="true"
BKG_BOXSIZE="20"
BKG_GRIDSIZE="0.2"
NPIX_MIN="5"
SEED_THR="5"
BRIGHT_SEED_THR="10"
MERGE_THR="2.6"
EXT_SFINDER_METHOD="3"
AC_METHOD="2"
SEARCH_NESTED_SOURCES="true"
SELECT_SOURCES="false"
SP_SIZE="20"
SP_BETA="1"
SP_MINAREA="10"
SALIENCY_THR="2.8"
SALIENCY_MIN_RESO="20"
SALIENCY_MAX_RESO="60"
SALIENCY_RESO_STEP="10"
SALIENCY_NN_PAR="1"
SEARCH_COMPACT_SOURCES="true"
SEARCH_EXTENDED_SOURCES="true"
FIT_SOURCES="false"

for item in $*
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
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
	
		## OPTIONAL ##	
		--submit*)
    	SUBMIT="true"
    ;;
		--queue=*)
    	BATCH_QUEUE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--loglevel=*)
    	LOG_LEVEL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--maxfiles=*)
    	NMAX_PROCESSED_FILES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--outdir=*)
    	OUTPUT_DIR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		--nproc=*)
      NPROC=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--hostfile=*)
    	HOSTFILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			HOSTFILE_GIVEN=true
    ;;
		--nthreads=*)
    	NTHREADS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
			
    --tilesize=*)
    	TILE_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			SPLIT_IN_TILES="true"
    ;;
		--tilestep=*)
    	TILE_STEP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			TILE_OVERLAP="true"
    ;;
		--xmin=*)
    	XMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			READ_TILE="true"
    ;;
		--xmax=*)
    	XMAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			READ_TILE="true"
    ;;
		--ymin=*)
    	YMIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			READ_TILE="true"
    ;;
		--ymax=*)
    	YMAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			READ_TILE="true"
    ;;
		
		
		--globalbkg*)
    	USE_LOCAL_BKG="false"
    ;;
		--bkgbox=*)
    	BKG_BOXSIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			USE_LOCAL_BKG="true"
    ;;
		--bkggrid=*)
    	BKG_GRIDSIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
			USE_LOCAL_BKG="true"
    ;;


		--no-compactsearch*)
    	SEARCH_COMPACT_SOURCES="false"
    ;;
		--no-nestedsearch*)
    	SEARCH_NESTED_SOURCES="false"
    ;;

		--npixmin=*)
    	NPIX_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		--seedthr=*)
    	SEED_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--brightseedthr=*)
    	BRIGHT_SEED_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--mergethr=*)
    	MERGE_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
	
		--no-extendedsearch*)
    	SEARCH_EXTENDED_SOURCES="false"
    ;;

		--extsfinder=*)
    	EXT_SFINDER_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--activecontour=*)
    	AC_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--selectsources*)
    	SELECT_SOURCES="true"
    ;;
		--fitsources*)
    	FIT_SOURCES="true"
    ;;

		--saliencythr*)
    	SALIENCY_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliencyminreso*)
    	SALIENCY_MIN_RESO=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliencymaxreso*)
    	SALIENCY_MAX_RESO=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliencyresostep*)
    	SALIENCY_RESO_STEP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliencynn*)
    	SALIENCY_NN_PAR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--spsize*)
    	SP_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--spbeta*)
    	SP_BETA=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--spminarea*)
    	SP_MINAREA=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
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
echo "SUBMIT? $SUBMIT, QUEUE=$BATCH_QUEUE"
echo "ENV_FILE: $ENV_FILE"
echo "INPUTFILE: $INPUTFILE"
echo "FILELIST: $FILELIST, NMAX_PROCESSED_FILES: $NMAX_PROCESSED_FILES"
echo "SPLIT_IN_TILES? $SPLIT_IN_TILES, TILE_SIZE: $TILE_SIZE, TILE_STEP: $TILE_STEP"
echo "NPROC: $NPROC, NTHREADS: $NTHREADS"
echo "HOSTFILE_GIVEN? $HOSTFILE_GIVEN, HOSTFILE: $HOSTFILE"
echo "LOG_LEVEL: $LOG_LEVEL"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "BKG BOX: $BKG_BOXSIZE, GRID: $BKG_GRID_SIZE"
echo "NPIX_MIN: $NPIX_MIN, SEED_THR: $SEED_THR, BRIGHT_SEED_THR: $BRIGHT_SEED_THR, MERGE_THR: $MERGE_THR"
echo "EXT_SFINDER_METHOD: $EXT_SFINDER_METHOD"
echo "AC_METHOD: $AC_METHOD"
echo "SEARCH_NESTED_SOURCES: $SEARCH_NESTED_SOURCES"
echo "SELECT_SOURCES: $SELECT_SOURCES"
echo "FIT_SOURCES: $FIT_SOURCES"
echo "SP PARS: ($SP_SIZE, $SP_BETA, $SP_MINAREA)"
echo "SALIENCY_THR: $SALIENCY_THR"
echo "SALIENCY_RESO: ($SALIENCY_MIN_RESO, $SALIENCY_MAX_RESO, $SALIENCY_RESO_STEP)"
echo "SALIENCY_NN_PAR: $SALIENCY_NN_PAR"
echo "****************************"
echo ""


## Check arguments parsed
if [ "$FILELIST_GIVEN" = false ] && [ "$INPUTFILE_GIVEN" = false ]; then
  echo "ERROR: Missing or empty FILELIST and INPUTFILE args (hint: you should specify at least one)!"
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



#######################################
##     DEFINE & LOAD ENV VARS
#######################################
export BASEDIR="$PWD"
export OUTPUT_DATADIR="$PWD"
export DATADIR=""

## Load env file
echo "INFO: Loading environment variables defined in file $ENV_FILE ..."
source $ENV_FILE


#######################################
##         CHECK ENVIRONMENT
#######################################
## Check if CAESAR environment variable exists and is properly set
if [ "$CAESAR_DIR" = "" ]; then
	echo "ERROR: Missing CAESAR_DIR environment variable, please set it to your CAESAR installation path."
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
		echo '###    CAESAR CONFIG OPTIONS'
		echo '############################################'
		echo '###'
		echo '//============================'
		echo '//==      INPUT             =='
		echo '//============================'
		echo "inputFile = $inputfile 					        | Input image filename (.root/.fits)"
		echo 'inputImage = img											  | Input image name in ROOT file'
		echo "readTileImage = $READ_TILE              | Read sub-image (T/F)"
		echo "tileMinX = $XMIN                        | Min x coords to be read in image"
		echo "tileMaxX = $XMAX												| Max x coords to be read in image"
		echo "tileMinY = $YMIN												| Min y coords to be read in image"
		echo "tileMaxY = $YMAX												| Max y coords to be read in image"
		echo '###'
		echo '###'
		echo '//============================'
		echo '//==      BEAM INFO         =='
		echo '//============================'
		echo 'beamFWHM = 6.5                          | User-supplied circular beam FWHM in arcsec (beamFWHM=BMAJ=BMIN, default=6.5 arcsec)'
		echo 'pixSize = 1.0   												| User-supplied map pixel area in arcsec (pixSize=CDELT, default=1 arcsec)'
		echo 'beamTheta = 0.0                         | User-supplied beam theta in deg (default=0)'
		echo '###'
		echo '###'
		echo '//============================================='
		echo '//==      DISTRIBUTED PROCESSING             =='
		echo '//============================================='
		echo "nThreads = $NTHREADS                               | Number of threads used if OPENMP is enabled (-1=all available threads)"
		echo "splitInTiles = $SPLIT_IN_TILES                     | Split input image in tiles (default=false)"
		echo "tileSizeX = $TILE_SIZE                             | Size of tile X (in pixels) to partition the input image"
    echo "tileSizeY = $TILE_SIZE									           | Size of tile Y (in pixels) to partition the input image"
		echo "useTileOverlap = $TILE_OVERLAP                     | Allow for tile overlap"
    echo "tileStepSizeX = $TILE_STEP                         | Tile step size fraction X to partition the input image (1=no overlap,0.5=half overlap, ...)"
    echo "tileStepSizeY = $TILE_STEP												 | Tile step size fraction Y to partition the input image (1=no overlap,0.5=half overlap, ...)"
		echo "mergeSourcesAtEdge = true 							           | Merge sources found at tile edge by each workers (default=true)"
    echo '###'
    echo '###'
    echo '//=============================='
    echo '//==      LOGGING             =='
    echo '//=============================='
    echo 'loggerTarget = 1                          | Logger target (1=CONSOLE, 2=FILE, 3=SYSLOG)'
    echo 'loggerTag = logger                        | Tag given to the log messages'
    echo "logLevel = $LOG_LEVEL                     | Log level threshold (DEBUG>INFO>WARN>ERROR>FATAL)"
    echo 'logFile = out.log                         | Log file name (for FILE target only)'
    echo 'appendToLogFile = false                   | If false a new log file is created, otherwise logs are appended (T/F)'
    echo 'maxLogFileSize = 10MB                     | Max size of log file before rotation (e.g. 10MB, 1KB)'
    echo 'maxBackupLogFiles = 2                     | Max number of backup files created after file threshold is reached'
    echo 'consoleTarget = System.out                | Console target (System.out/System.err)'
    echo 'syslogFacility = local6                   | Syslog facility used with syslog target'
    echo '###'
    echo '###'
    echo '//============================='
    echo '//==      OUTPUT             =='
    echo '//============================='
    echo 'isInteractiveRun = false                          | Is interactive run (graph plots enabled) (T/F)'
    echo "outputFile = $outputfile                         | Output filename (.root)"
    echo "ds9RegionFile = $ds9region_file					          | DS9 region file (.reg) where to store source catalog"
    echo 'DS9RegionFormat = 2                               | DS9 region format (1=ellipse, 2=polygon)'
    echo 'inputMapFITSFile = 	input_map.fits				        | Output filename where to store input map in FITS format (.fits)'
    echo 'residualMapFITSFile = residual_map.fits		        | Output filename where to store residual map in FITS format (.fits)'
    echo 'saliencyMapFITSFile = saliency_map.fits		        | Output filename where to store saliency map in FITS format (.fits)'
    echo 'bkgMapFITSFile = bkg_map.fits							        | Output filename where to store bkg map in FITS format (.fits)'
		echo 'noiseMapFITSFile = noise_map.fits					        | Output filename where to store noise map in FITS format (.fits)'
    echo 'significanceMapFITSFile = significance_map.fits	  | Output filename where to store significance map in FITS format (.fits)' 
    echo 'saveToFile = true																	| Save results & maps to output ROOT file (T/F)'
    echo 'saveToFITSFile = false														| Save results to output FITS file(s) (T/F)'
    echo 'saveInputMap = false															| Save input map to ROOT file (T/F)'
    echo 'saveConfig = true																	| Save config options to ROOT file (T/F)'
    echo 'saveResidualMap = false														| Save residual map to ROOT file (T/F)'
    echo 'saveBkgMap = false																| Save bkg map to ROOT file (T/F)'
    echo 'saveNoiseMap = false															| Save noise map to ROOT file (T/F)'
    echo 'saveSignificanceMap = false												| Save significance map to ROOT file (T/F)'
    echo 'saveSaliencyMap = false														| Save saliency map to ROOT file (T/F)'
    echo 'saveSources = true																| Save sources to ROOT file (T/F)'
    echo 'saveEdgenessMap = false                           | Save edgeness map computed in extended source search to ROOT file (T/F)'
    echo 'saveCurvatureMap = false                          | Save curvature map to ROOT file (T/F)'
    echo 'saveSegmentedMap = false                          | Save segmented map computed in extended source search to ROOT file (T/F)'
    echo '###'
    echo '###'
		echo '//==========================='
		echo '//==   BKG OPTIONS         =='
		echo '//==========================='
    echo "useLocalBkg = $USE_LOCAL_BKG							| Use local background calculation instead of global bkg (T/F)"
    echo 'localBkgMethod = 1												| Local background method (1=Grid, 2=Superpixel)'
    echo 'use2ndPassInLocalBkg = true								| Use 2nd pass to refine noise calculation in local bkg (T/F)'
		echo 'skipOutliersInLocalBkg = false						| Skip outliers (e.g. bright point sources) in local bkg computation (T/F)'
		echo 'bkgEstimator = 2													| Background estimator (1=Mean,2=Median,3=BiWeight,4=ClippedMedian)'
    echo 'useBeamInfoInBkg = true                   | Use beam information in bkg box definition (if available) (T/F)'
		echo "boxSizeX = $BKG_BOXSIZE										| X Size of local background box in #pixels"
		echo "boxSizeY = $BKG_BOXSIZE										| Y Size of local background box in #pixels"
		echo "gridSizeX = $BKG_GRIDSIZE									| X Size of local background grid used for bkg interpolation"
		echo "gridSizeY = $BKG_GRIDSIZE									| Y Size of local background grid used for bkg interpolation"
		echo '###'
		echo '###'
		echo '//==============================='
		echo '//==  FILTERING OPTIONS        =='
		echo '//==============================='
		echo 'usePreSmoothing = true										          | Use a pre-smoothing stage to filter input image for extended source search (T/F)'
		echo 'smoothFilter = 2													          | Smoothing filter (1=gaus,2=guided)'
		echo 'gausFilterKernSize = 5										          | Gaussian filter kernel size'
		echo 'gausFilterSigma = 1																	| Gaussian filter sigma'
		echo 'guidedFilterRadius = 12															| Guided filter radius par'
		echo 'guidedFilterColorEps = 0.04													| Guided filter color epsilon parameter'
		echo '###'
		echo '###'
		echo '//===================================='
		echo '//==  SOURCE FINDING OPTIONS        =='
		echo '//===================================='
		echo "searchCompactSources = $SEARCH_COMPACT_SOURCES			| Search compact sources (T/F)"
		echo "minNPix = $NPIX_MIN																	| Minimum number of pixel to consider a source"
		echo "seedBrightThr = $BRIGHT_SEED_THR										| Seed threshold in flood-filling algo for bright sources"
		echo "seedThr = $SEED_THR 																| Seed threshold in flood filling algo for faint sources"
		echo "mergeThr = $MERGE_THR																| Merge/aggregation threshold in flood filling algo"
		echo 'mergeBelowSeed = false                              | Aggregate to seed only pixels above merge threshold but below seed threshold (T/F)'
		echo 'searchNegativeExcess = false												| Search negative excess together with positive in compact source search'
		echo "compactSourceSearchNIters = 10                      | Number of iterations to be performed in compact source search (default=10)"
		echo '###'
		echo '###'
		echo '//==========================================='
		echo '//==  NESTED SOURCE FINDING OPTIONS        =='
		echo '//==========================================='
		echo "searchNestedSources = $SEARCH_NESTED_SOURCES 			  | Search for nested sources inside candidate sources (T/F)"
		echo 'nestedBlobThrFactor = 1                             | Threshold (multiple of curvature rms) used for nested blob finding'
		echo 'minNestedMotherDist = 2                             | Minimum distance in pixels (in x or y) between nested and parent blob below which nested is skipped'
		echo 'maxMatchingPixFraction = 0.5                        | Maximum fraction of matching pixels between nested and parent blob above which nested is skipped'
		echo '###'
		echo '###'
		echo '//=================================='
		echo '//==  Source fitting options   =='
		echo '//=================================='
		echo "fitSources = $FIT_SOURCES                      | Deblend point-like sources with multi-component gaus fit (T/F)"
		echo 'fitMaxNComponents = 3                          | Maximum number of components fitted in a blob (T/F)'
		echo '###'
		echo '###'
		echo '//============================================='
		echo '//==  EXTENDED SOURCE FINDING OPTIONS        =='
		echo '//============================================='
		echo "searchExtendedSources = $SEARCH_EXTENDED_SOURCES		| Search extended sources after bright source removal (T/F)"
		echo "extendedSearchMethod = $EXT_SFINDER_METHOD					| Extended source search method (1=WT-thresholding,2=SPSegmentation,3=ActiveContour,4=Saliency thresholding)"
		echo 'useResidualInExtendedSearch = true									| Use residual image (with selected sources dilated) as input for extended source search'
		echo "activeContourMethod = $AC_METHOD										| Active contour method (1=Chanvese, 2=LRAC)"
		echo 'wtScaleExtended = 6																	| Wavelet scale to be used for extended source search'
		echo '###'
		echo '###'
		echo '//================================'
		echo '//==  SOURCE SELECTION OPTIONS  =='
		echo '//================================'
		echo "applySourceSelection = $SELECT_SOURCES  						| Apply selection cuts to sources (T/F)"
		echo 'useMinBoundingBoxCut = true													| Use bounding box cut (T(F)'
		echo 'sourceMinBoundingBox = 2														| Minimum bounding box cut (source tagged as bad if below this threshold)'
		echo 'useCircRatioCut = false															| Use circularity ratio cut (T/F)'
		echo 'psCircRatioThr = 0.4																| Circular ratio threshold (source passes point-like cut if above this threshold)'
		echo 'useElongCut = false																	| Use elongation cut (T/F)'
		echo 'psElongThr = 0.7																		| Elongation threshold (source passes point-like cut if below this threshold'
		echo 'useEllipseAreaRatioCut = false											| Use Ellipse area ratio cut (T/F)'
		echo 'psEllipseAreaRatioMinThr = 0.6											| Ellipse area ratio min threshold'
		echo 'psEllipseAreaRatioMaxThr = 1.4											| Ellipse area ratio max threshold'
		echo 'useMaxNPixCut = true																| Use max npixels cut (T/F)'
		echo 'psMaxNPix = 1000																		| Max number of pixels for point-like sources (source passes point-like cut if below this threshold)'
		echo '###'
		echo '###'
		echo '//================================'
		echo '//==  SOURCE RESIDUAL OPTIONS   =='
		echo '//================================'
		echo 'dilateNestedSources = true 													| Dilate sources nested inside bright sources (T/F)'
		echo 'dilateKernelSize = 9																| Size of kernel (odd) to be used in dilation operation'
		echo 'dilatedSourceType = 2																| Type of bright sources to be dilated from the input image (-1=ALL,1=COMPACT,2=POINT-LIKE,3=EXTENDED) '
		echo 'dilateSourceModel = 1																| Dilate source model  (1=bkg,2=median)'
		echo 'dilateRandomize = false															| Randomize dilated values (T/F)'
		echo '###'
		echo '###'
		echo '//==================================='
		echo '//==  CHAN-VESE ALGORITHM OPTIONS  =='
		echo '//==================================='
		echo 'cvNIters = 1000 																		| Number of iterations'
		echo 'cvTimeStepPar = 0.007																| Chan-Vese time step par'
		echo 'cvWindowSizePar = 1																	| Chan-Vese window size par'
		echo 'cvLambda1Par = 1																		| Chan-Vese lambda1 par'
		echo 'cvLambda2Par = 2																		| Chan-Vese lambda2 par'
		echo 'cvMuPar = 0.5																				| Chan-Vese mu par'
		echo 'cvNuPar = 0																					|	Chan-Vese nu par'
		echo 'cvPPar = 1																					| Chan-Vese p par'
		echo '###'
		echo '###'
		echo '//==================================='
		echo '//==  LRAC ALGORITHM OPTIONS       =='
		echo '//==================================='
		echo 'lracNIters = 1000 																	| Number of iterations'
		echo 'lracLambdaPar = 0.1																  | Regularization par'
		echo 'lracRadiusPar = 1																	  | Radius of locatization ball par'
		echo 'lracEpsPar = 0.1																	  | Convergence par'
		echo '###'
		echo '###'
		echo '//==============================='
		echo '//==  SALIENCY FILTER OPTIONS  =='
		echo '//==============================='
		echo "saliencyThrFactor = $SALIENCY_THR										| Saliency threshold factor for tagging signal regions (thr=<saliency>*factor)"
		echo 'saliencyBkgThrFactor = 1														| Saliency threshold factor for tagging bkg regions (thr=<saliency>*factor)'
		echo 'saliencyImgThrFactor = 1														| Threshold factor to consider a region as significant (thr=<img>*factor)'
		echo "saliencyResoMin = $SALIENCY_MIN_RESO								| Saliency min reso par"
		echo "saliencyResoMax = $SALIENCY_MAX_RESO								| Saliency max reso par"
		echo "saliencyResoStep = $SALIENCY_RESO_STEP							| Saliency reso step par"
		echo 'saliencyUseCurvInDiss = false 											| Use curvature parameter in dissimilarity estimation (T/F)'
		echo 'saliencyUseRobustPars = false												| Use robust pars in saliency map computation (T/F)'
		echo 'saliencyUseBkgMap = true														| Use bkg map in saliency map computation (T/F)'
		echo 'saliencyUseNoiseMap = true					 								| Use noise map in saliency map computation (T/F)'
		echo "saliencyNNFactor = $SALIENCY_NN_PAR   							| Fraction of most similar neighbors used in saliency map computation"
		echo 'saliencySpatialRegFactor = 6												| Spatial regularization factor (ruling exp decay in saliency spatial weighting)'
		echo 'saliencyMultiResoCombThrFactor = 0.7								| Fraction of combined salient multi-resolution maps to consider global saliency'
		echo 'saliencyDissExpFalloffPar = 100                     | Dissimilarity exponential cutoff parameter (value)'
		echo 'saliencySpatialDistRegPar = 1                       | Spatial-color distance regularization par (value, 1=equal weights)'
		echo '###'
		echo '###'
		echo '//=================================='
		echo '//==    SP SEGMENTATION OPTIONS   =='
		echo '//=================================='
		echo "spSize = $SP_SIZE													  | Initial superpixel size"
		echo "spBeta = $SP_BETA														| Initial superpixel regularization parameter"
		echo "spMinArea = $SP_MINAREA											| Initial superpixel min area"
		echo 'spUseLogContrast = false										| Use intensity log scale to generate initial superpixel partition (T/F)'
		echo '###'
		echo '###'
		echo '//=================================='
		echo '//==    SP MERGING OPTIONS        =='
		echo '//=================================='
		echo 'spMergingNSegmentsToStop = 1                 | Number of segments below which the segmentation is stopped (default=1)'
		echo 'spMergingRatio = 0.3                         | Fraction of similar segments merged per each hierarchical level (default=0.3)'
		echo 'spMergingRegPar = 0.5                        | Regularization parameters balancing region edge and similarity (default=0.5)'
		echo 'spMergingMaxDissRatio = 1000                 | Mutual segment dissimilarity ratio (R=Diss(ji)/Diss(ij)) above which 1st-order neighbors are not merged (default=1000)'
		echo 'spMergingMaxDissRatio2ndNeighbours = 1.05    | Mutual segment dissimilarity ratio (R=Diss(ji)/Diss(ij)) above which 2nd-order neighbors are not merged (default=1.05)'
		echo 'spMergingDissThreshold = 3                   | Absolute dissimilarity ratio (wrt to initial average diss) threshold above which segments are not merged (default=3)'
		echo 'spMergingEdgeModel = 2                       | Superpixel edge model (1=Kirsch,2=Chan-Vese) '
		echo 'spMergingIncludeSpatialPars = false          | Include spatial pars in region dissimilarity measure (default=false)'
		echo '###'
		echo '###'
		echo '//spHierMergingPars = 1  0.3  0.5  3         | Hierarch algo pars: MIN_SEGMENTS, MERGING_RATIO, DIST_REGULARIZATION, DIST_THRESHOLD es 0.25 (value,value,value,value)'
		echo '//spHierMergingMaxDissRatio = 1000  1.05     | Maximum mutual dissimilarity among regions for merging 1st and 2nd neighbors es. 1.15 (value)'
		echo '//spMergingEdgeModel = 2                     | Edge model (1=Kirsch,2=Chan-Vese) (value)'
		echo '//use2ndNeighborsInSPMerging  T              | Use 2nd-order neighbors in superpixel merging (T/F)'
		echo '//useCurvatureInSPMerging  T                 | Use curvature params in superpixel merging (T/F)'
		echo '//useLogContrastInSPGeneration  F            | Use logarithmic contrast to generate initial partition (T/F)'
		echo '//usePixelRatioCut  T  0.4                   | Use pixel ratio to tag background regions (if npix/npixtot>cut)  (T/F,value)'
		echo '//tagSignificantSP  T  1  0.5                | Tag significant superpixels if fraction of significant subregions is above cut (1=saliency,2=Mahalanobis) (T/F,method,cut)'
		echo '###'
		echo '###'

 ) > $configfile

}
## close function

#######################################
##   DEFINE GENERATE EXE SCRIPT FCN
#######################################
generate_exec_script(){

	local shfile=$1
	local jobindex=$2
	local exe=$3
	local exe_args=$4

	echo "INFO: Creating sh file $shfile (jobindex=$jobindex, exe=$exe, exe_args=$exe_args)..."
	( 
  		echo "#PBS -o $BASEDIR"
    	echo "#PBS -o $BASEDIR"
    	echo '#PBS -r n'
      echo '#PBS -S /bin/sh'
      echo "#PBS -N SFinderJob$jobindex"
      echo '#PBS -p 1'

      echo " "
      echo " "

      echo 'echo "*************************************************"'
      echo 'echo "****         PREPARE JOB                     ****"'
      echo 'echo "*************************************************"'

      echo 'echo ""'
            
      echo " "

      echo 'echo ""'
      echo 'echo "INFO: Source the software environment ..."'
      echo "source $ENV_FILE"

      echo 'echo ""'
      
      echo "export JOBDIR=$BASEDIR"
     
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
		outputfile="Out_$filename_base_noext"'.root'
 		ds9region_file="DS9_$filename_base_noext"'.reg'

		## Define and generate config file
		configfile="config_$filename_base_noext"'_'"$index.cfg"
  	echo "INFO: Creating config file $configfile for input file: $inputfile ..."
		generate_config $configfile


		## Define executable & args variables and generate script
		shfile="Run_$filename_base_noext"'_'"$index.sh"
		EXE="$CAESAR_DIR/scripts/RunSFinderMPI.sh"
		EXE_ARGS="--nproc=$NPROC --config=$configfile"
		if [ "$HOSTFILE_GIVEN" = true ] ; then
			EXE_ARGS="$EXE_ARGS --hostfile=$HOSTFILE"
		fi
		
		echo "INFO: Creating script file $shfile for input file: $inputfile ..."
		generate_exec_script "$shfile" "$index" "$EXE" "$EXE_ARGS"


		# Submits the job to batch system
		if [ "$SUBMIT" = true ] ; then
			echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
			qsub -q $BATCH_QUEUE $CURRENTJOBDIR/$shfile
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
	outputfile="Out_$filename_base_noext"'.root'
 	ds9region_file="DS9_$filename_base_noext"'.reg'

	## Define and generate config file
	configfile="config_$filename_base_noext"'.cfg'
  echo "INFO: Creating config file $configfile for input file: $inputfile ..."
	generate_config $configfile


	## Define executable & args variables and generate script
	shfile="Run_$filename_base_noext"'.sh'
	EXE="$CAESAR_DIR/scripts/RunSFinderMPI.sh"
	EXE_ARGS="--nproc=$NPROC --config=$configfile"
	if [ "$HOSTFILE_GIVEN" = true ] ; then
		EXE_ARGS="$EXE_ARGS --hostfile=$HOSTFILE"
	fi

	echo "INFO: Creating script file $shfile for input file: $inputfile ..."
	jobId=" "
	generate_exec_script "$shfile" "$jobId" "$EXE" "$EXE_ARGS"

	# Submits the job to batch system
	if [ "$SUBMIT" = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE ..."
		qsub -q $BATCH_QUEUE $CURRENTJOBDIR/$shfile
	fi

fi

echo "*** END SUBMISSION ***"

