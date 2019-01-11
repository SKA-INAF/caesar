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

	echo "*** OPTIONAL ARGS ***"
	echo "=== SFINDER OUTPUT OPTIONS ==="
	echo "--save-inputmap - Save input map in output ROOT file (default=no)"
	echo "--save-bkgmap - Save bkg map in output ROOT file (default=no)"
	echo "--save-rmsmap - Save rms map in output ROOT file (default=no)"	
	echo "--save-significancemap - Save significance map in output ROOT file (default=no)"
	echo "--save-residualmap - Save residual map in output ROOT file (default=no)"
	echo "--save-saliencymap - Save saliency map in output ROOT file (default=no)"
	echo "--save-segmentedmap - Save segmented map in output ROOT file (default=no)"
	echo "--save-regions - Save DS9 regions (default=no)"
	echo "--convertregionstowcs - Save DS9 regions in WCS format (default=no)"
	echo "--regionwcs=[DS9REGION_WCSTYPE] - DS9 region WCS output format (0=J2000,1=B1950,2=GALACTIC) (default=0)"
	echo ""

	echo "=== SFINDER IMG READ OPTIONS ==="
	echo "--tilesize=[TILE_SIZE] - Size (in pixels) of tile used to partition input image in distributed processing (default=0=no tile split)"		
	echo "--tilestep=[TILE_STEP] - Tile step size (range 0-1) expressed as tile fraction used in tile overlap (default=1=no overlap)"
	echo "--xmin=[XMIN] - Read sub-image of input image starting from pixel x=xmin (default=0=read full image)"
	echo "--xmax=[XMAX] - Read sub-image of input image up to pixel x=xmax (default=0=read full image)"
	echo "--ymin=[YMIN] - Read sub-image of input image starting from pixel y=xmin (default=0=read full image)"
	echo "--ymax=[YMAX] - Read sub-image of input image up to pixel y=ymax (default=0=read full image)"
	echo ""

	echo "=== SFINDER STATS OPTIONS ==="
	echo "--no-parallelmedian - Switch off parallel median algorithm (based on parallel nth-element)"
	echo ""

	echo "=== SFINDER BKG OPTIONS ==="
	echo "--bmaj=[BMAJ] - User-supplied beam Bmaj in arcsec (NB: used only when beam info is not available in input map) (default: 10 arcsec)"
	echo "--bmin=[BMIN] - User-supplied beam Bmin in arcsec (NB: used only when beam info is not available in input map) (default: 5 arcsec)"
	echo "--bpa=[BMIN] - User-supplied beam position angle in degrees (NB: used only when beam info is not available in input map) (default: 0 deg)"
	echo "--mappixsize=[MAP_PIXSIZE] - Map pixel size in arcsec (NB: used only when info is not available in input map) (default=1 arcsec)"
	echo "--globalbkg - Use global bkg (default=use local bkg)"
	echo "--bkgestimator=[BKG_ESTIMATOR] - Stat estimator used for bkg (1=Mean,2=Median,3=BiWeight,4=ClippedMedian) (default=2)"
	echo "--bkgbox=[BKG_BOXSIZE] - Box size (muliple of beam size) used to compute local bkg (default=20 x beam)"
	echo "--bkggrid=[BKG_GRIDSIZE] - Grid size (fraction of bkg box) used to compute local bkg (default=0.2 x box)"
	echo "--no-bkg2ndpass - Do not perform a 2nd pass in bkg estimation (default=true)"
	echo "--bkgskipoutliers - Remove bkg outliers (blobs above seed thr) when estimating bkg (default=no)"
	echo "--sourcebkgboxborder=[SOURCE_BKGBOX_BORDER] - Border size (in pixels) of box around source bounding box used to estimate bkg for fitting (default=20)"
	echo ""

	echo "=== SFINDER SOURCE FINDING OPTIONS ==="
	echo "--mergeedgesources - Merge sources at tile edges. NB: Used for multitile processing. (default=no)"
	echo "--no-mergesources - Disable source merging in each tile (default=enabled)."
	echo "--no-mergecompactsources - Disable compact-compact source merging in each tile (default=enabled)."
	echo "--no-mergeextsources - Disable extended-extended and extended-compact source merging in each tile (default=enabled)."
	echo ""

	echo "=== SFINDER COMPACT SOURCE OPTIONS ==="
	echo "--no-compactsearch - Do not search compact sources"
	echo "--npixmin=[NPIX_MIN] - Minimum number of pixel to form a compact source (default=5 pixels)"
	echo "--seedthr=[SEED_THR] - Seed threshold (in nsigmas) used in flood-fill (default=5 sigmas)"
	echo "--mergethr=[MERGE_THR] - Merge threshold (in nsigmas) used in flood-fill (default=2.6 sigmas)"
	echo "--compactsearchiters=[COMPACT_SOURCE_SEARCH_NITERS] - Maximum number of compact source search iterations (default=5)"	
	echo "--seedthrstep=[SEED_THR_STEP] - Seed thr decrease step across iterations (default=1)"
	echo ""

	echo "=== SFINDER COMPACT SOURCE SELECTION OPTIONS ==="
	echo "--selectsources - Apply selection to compact sources found (default=false)"	
	echo "--no-boundingboxcut - Do not apply bounding box cut (default=apply)"
	echo "--minboundingbox=[MIN_BOUNDING_BOX] - Minimum bounding box cut in pixels (NB: source tagged as bad if below this threshold) (default=2)"
	echo "--no-circratiocut - Do not apply circular ratio parameter cut (default=apply)"
	echo "--circratiothr=[CIRC_RATIO_THR] - Circular ratio threshold (0=line, 1=circle) (source passes point-like cut if above this threshold) (default=0.4)"
	echo "--no-elongationcut - Do not apply elongation parameter cut (default=apply)"
	echo "--elongationthr=[ELONGATION_THR] - Elongation threshold (source passes point-like cut if below this threshold (default=0.7)"
	echo "--ellipsearearatiocut - Apply ellipse area ratio parameter cut (default=not applied)"	
	echo "--ellipsearearatiominthr=[ELLIPSE_AREA_RATIO_MIN_THR] - Ellipse area ratio min threshold (default=0.6)"	
	echo "--ellipsearearatiomaxthr=[ELLIPSE_AREA_RATIO_MAX_THR] - Ellipse area ratio max threshold (default=1.4)"	
	echo "--maxnpixcut - Apply max pixels cut (NB: source below this thr passes the point-like cut) (default=not applied)"	
	echo "--maxnpix=[MAX_NPIX] - Max number of pixels for point-like sources (source passes point-like cut if below this threshold) (default=1000)"
	echo "--no-nbeamscut - Use number of beams in source cut (default=applied)"
	echo "--nbeamsthr=[NBEAMS_THR] - nBeams threshold (sources passes point-like cut if nBeams<thr) (default=3)"
	echo ""

	echo "=== SFINDER EXTENDED SOURCE OPTIONS ==="
	echo "--no-extendedsearch - Do not search extended sources"
	echo "--extsfinder=[EXT_SFINDER_METHOD] - Extended source search method {1=WT-thresholding,2=SPSegmentation,3=ActiveContour,4=Saliency thresholding} (default=3)"	
	echo "--activecontour=[AC_METHOD] - Active contour method {1=Chanvese, 2=LRAC} (default=2)"
	echo ""

	echo "=== SFINDER COMPACT NESTED SOURCE OPTIONS ==="
	echo "--no-nestedsearch - Do not search nested sources (default=search)"		
	echo "--blobmaskmethod=[BLOB_MASK_METHOD] - Blob mask method (1=gaus smooth+Laplacian,2=multi-scale LoG) (default=2)"
	echo "--nested-sourcetobeamthr=[NESTED_SOURCE_TO_BEAM_THR] - Source area/beam thr to add nested sources (e.g. npix>thr*beamArea). NB: thr=0 means always if searchNestedSources is enabled (default=5)"
	echo "--nested-blobthr=[NESTED_BLOB_THR] - Threshold (multiple of curvature median) used for nested blob finding (default=0)"
	echo "--nested-minmotherdist=[NESTED_MIN_MOTHER_DIST] - Minimum distance in pixels (in x or y) between nested and parent blob below which nested is skipped (default=2)"
	echo "--nested-maxmotherpixmatch=[NESTED_MAX_MOTHER_PIX_MATCH] - Maximum fraction of matching pixels between nested and parent blob above which nested is skipped (default=0.5)"
	echo "--nested-blobpeakzthr=[NESTED_BLOB_PEAK_ZTHR] - Nested blob peak significance threshold (in scale curv map) (default=5)"
	echo "--nested-blobpeakzthrmerge=[NESTED_BLOB_PEAK_ZTHR_MERGE] - Nested blob significance merge threshold (in scale curv map) (default=2.5)"
	echo "--nested-blobminscale=[NESTED_BLOB_MIN_SCALE] - Nested blob min scale search factor f (blob sigma_min=f x beam width) (default=1)"
	echo "--nested-blobmaxscale=[NESTED_BLOB_MAX_SCALE] - Nested blob max scale search factor f (blob sigma_max=f x beam width) (default=3)"
	echo "--nested-blobscalestep=[NESTED_BLOB_SCALE_STEP] - Nested blob scale step (sigma=sigma_min + step) (default=1)"
	echo "--nested-blobkernfactor=[NESTED_BLOB_KERN_FACTOR] - Nested blob curvature/LoG kernel size factor f (kern size=f x sigma) (default=1)"
	
	echo ""
	
	echo "=== SFINDER SOURCE RESIDUAL OPTIONS ==="
	echo "--res-removenested - When a source has nested sources remove only on nested (default=false)"
	echo "--res-zthr=[RESIDUAL_ZTHR] - Seed threshold (in nsigmas) used to dilate sources (default=5 sigmas)"	
	echo "--res-zhighthr=[RESIDUAL_ZHIGHTHR] - Seed threshold (in nsigmas) used to dilate sources (even if they have nested components or different dilation type) (default=10 sigmas)"		
	echo "--dilatekernsize=[DILATE_KERNEL_SIZE] - Size of dilating kernel in pixels (default=9)"
	echo "--res-removedsourcetype=[RESIDUAL_REMOVED_SOURCE_TYPE] - Type of source dilated from the input image (-1=ALL,1=COMPACT,2=POINT-LIKE,3=EXTENDED) (default=2)"
	echo "--res-pssubtractionmethod=[PS_SUBTRACTION_METHOD] - Method used to subtract point-sources in residual map (1=DILATION, 2=FIT MODEL REMOVAL)"
	echo ""

	echo "=== SFINDER SOURCE FITTING OPTIONS ==="
	echo "--fitsources - Fit compact point-like sources found (default=false)"
	echo "--fit-usethreads - Enable multithread in source fitting (NB: use Minuit2 minimizer if enabled) (default=disabled)"
	echo "--fit-minimizer=[FIT_MINIMIZER] - Fit minimizer {Minuit,Minuit2} (default=Minuit2)"
	echo "--fit-minimizeralgo=[FIT_MINIMIZER_ALGO] - Fit minimizer algo {migrad,simplex,minimize,scan,fumili (Minuit2)} (default=minimize)"
	echo "--fit-printlevel=[FIT_PRINTLEVEL] - Fit print level (default=1)"
	echo "--fit-strategy=[FIT_STRATEGY] - Fit strategy (default=2)"
	echo "--fit-maxnbeams=[FIT_MAX_NBEAMS] - Maximum number of beams for fitting if compact source (default=20)"
	echo "--fit-maxcomponents=[FIT_MAX_COMPONENTS] - Maximum number of components fitted in a blob (default=3)"	
	echo "--fit-usenestedascomponents - Initialize fit components to nested sources found in source (default=no)"
	echo "--fit-freebkg - Fit with bkg offset parameter free to vary (default=fixed)"
	echo "--fit-estimatedbkg - Set bkg par starting value to estimated bkg (average over source pixels by default, box around source if --fit-estimatedboxbkg is given) (default=use fixed bkg start value)"
	echo "--fit-usebkgboxestimate - Set bkg par starting value to estimated bkg (from box around source)"
	echo "--fit-bkg=[FIT_BKG] - Bkg par starting value (NB: ineffective when -fit-estimatedbkg is enabled) (default=0)"
	echo "--fit-ampllimit=[FIT_AMPL_LIMIT] - Limit amplitude range par (Speak*(1+-FIT_AMPL_LIMIT)) (default=0.3)"
	echo "--prefit-freeampl - Set free amplitude par in pre-fit (default=fixed)"
	echo "--fit-sigmalimit=[FIT_SIGMA_LIMIT] - Gaussian sigma limit around psf or beam (Bmaj*(1+-FIT_SIGMA_LIMIT)) (default=0.3)"
	echo "--fit-thetalimit=[FIT_THETA_LIMIT] - Gaussian theta limit around psf or beam in degrees (e.g. Bpa +- FIT_THETA_LIMIT) (default=90)"
	echo "--fit-nobkglimits - Do not apply limits in bkg offset parameter in fit (default=fit with limits when par is free)"
	echo "--fit-noampllimits - Do not apply limits in Gaussian amplitude parameters in fit (default=fit with limits)"
	echo "--fit-nosigmalimits - Do not apply limits in Gaussian sigma parameters in fit (default=fit with limits)"	
	echo "--fit-noposlimits - Do not apply limits in Gaussian mean parameters in fit (default=fit with limits)"
	echo "--fit-poslimit=[FIT_POS_LIMIT] - Source centroid limits in pixel (default=3)"
	echo "--prefit-freepos - Set free centroid pars in pre-fit (default=fixed)"
	echo "--fit-nothetalimits - Do not apply limits in Gaussian ellipse pos angle parameters in fit (default=fit with limits)"
	echo "--fit-fixsigma - Fit with sigma parameters fixed to start value (beam bmaj/bmin) (default=fit with sigma free and constrained)"
	echo "--prefit-fixsigma - Fix sigma parameters in pre-fit (default=free)"	
	echo "--fit-fixtheta - Fit with theta parameters fixed to start value (beam bpa) (default=fit with theta free and constrained)"
	echo "--prefit-fixtheta - Fix theta parameter in pre-fit (default=free)"	
	echo "--fit-peakminkern=[PEAK_MIN_KERNEL_SIZE] - Minimum dilation kernel size (in pixels) used to detect peaks (default=3)"
	echo "--fit-peakmaxkern=[PEAK_MAX_KERNEL_SIZE] - Maximum dilation kernel size (in pixels) used to detect peaks (default=7)"
	echo "--fit-peakmultiplicitythr=[PEAK_KERNEL_MULTIPLICITY_THR] - Requested peak multiplicity across different dilation kernels (-1=peak found in all given kernels,1=only in one kernel, etc) (default=1)"
	echo "--fit-peakshifttol=[PEAK_SHIFT_TOLERANCE] - Shift tolerance (in pixels) used to compare peaks in different dilation kernels (default=2 pixels)"
	echo "--fit-peakzthrmin=[PEAK_ZTHR_MIN] - Minimum peak flux significance (in nsigmas above avg source bkg & noise) below which peak is skipped (default=1)"
	echo "--fit-fcntol=[FIT_FCNTOL] - Fit function tolerance for convergence (default 1.e-2)"
	echo "--fit-maxniters=[FIT_MAXNITERS] - Maximum number of fit iterations or function calls performed (default 10000)"
	echo "--fit-noimproveconvergence - Do not use iterative fitting to try to achieve fit convergence (default=use)"
	echo "--fit-noretry - Do not iteratively retry fit with less components in case of failed convergence (default=retry)"
	echo "--fit-nretries=[FIT_NRETRIES] - Maximum number of fit retries if fit failed or has parameters at bound (default 10)"
	echo "--fit-parboundincreasestep - Fit par bound increase step size (e.g. parmax= parmax_old+(1+nretry)*fitParBoundIncreaseStepSize*0.5*|max-min|). Used in iterative fitting. (default=0.1)"
	echo "--fit-improveerrors - Run final minimizer step (e.g. HESS) to improve fit error estimates (default=no)"
	echo "--fit-scaledatatomax - Scale source data to max pixel flux for fitting. Otherwise scale to mJy (default=no)"
	echo "--fit-nochi2cut - Do not apply reduced chi2 cut to fitted sources (default=apply)"
	echo "--fit-chi2cut=[FIT_CHI2_CUT] - Chi2 cut value (default=5)"
	echo "--fit-useellipsecuts - Apply ellipse cuts to fitted sources (default=not applied)"
	echo ""

	echo "=== SFINDER SMOOTHING FILTER OPTIONS ==="
	echo "--no-presmoothing - Do not smooth input/residual map before extended source search (default=yes)"
	echo "--smoothfilter=[SMOOTH_FILTER] - Smoothing filter to be used (1=gaussian, 2=guided filter) (default=2)"
	echo "--guidedfilter-radius=[GUIDED_FILTER_RADIUS] - Guided filter radius par (default=12)"
	echo "--guidedfilter-eps=[GUIDED_FILTER_EPS] - Guided filter eps par (default=0.04)"
	echo ""
	
	echo "=== SFINDER SALIENCY FILTER OPTIONS ==="
  echo "--sp-size - Superpixel size (in pixels) used in hierarchical clustering (default=20)"	
	echo "--sp-beta - Superpixel regularization par (beta) used in hierarchical clustering (default=1)"	
	echo "--sp-minarea - Superpixel min area (in pixels) used in hierarchical clustering (default=10)"
	echo "--saliency-nooptimalthr - Do not use optimal threshold in multiscale saliency estimation (e.g. use median thr) (default=use optimal)"	
	echo "--saliency-thr - Saliency map threshold factor wrt optimal/median threshold (default=2.8)"
	echo "--saliency-minreso - Superpixel size (in pixels) used in multi-reso saliency map smallest scale (default=20 pixels)"
	echo "--saliency-maxreso - Superpixel size (in pixels) used in multi-reso saliency map highest scale (default=60 pixels)"
	echo "--saliency-resostep - Superpixel size step (in pixels) used in multi-reso saliency map computation (default=10 pixels)"
	echo "--saliency-nn - Fraction of most similar region neighbors used in saliency map computation (default=1)"
	echo "--saliency-usebkgmap - Use bkg map in saliency computation (default=not used)"
	echo "--saliency-usermsmap - Use noise map in saliency computation (default=not used)"
	echo "--saliency-userobustpars - Use robust pars in saliency computation (default=no)"
	echo ""

	echo "=== SFINDER ACTIVE-CONTOUR MAIN OPTIONS ==="
	echo "--ac-niters=[AC_NITERS] - Maximum number of iterations in active-contour algorithms (default=1000)"
	echo "--ac-levelset=[AC_LEVELSET_METHOD] - Init level set method in active-contour algorithms (1=circle,2=checkerboard,3=saliency) (default=1)"
	echo "--ac-levelsetsize=[AC_LEVELSET_SIZE] - Init level set size par in active-contour algorithms (default=0.1ximage size)"
	echo "--ac-tolerance=[AC_TOLERANCE] - Tolerance par in active-contour algorithms (default=0.1)"
	echo ""

	echo "=== SFINDER CHAN-VESE OPTIONS ==="
	#echo "--cv-niters=[CV_NITERS] - Maximum number of iterations in ChanVese algorithms (default=1000)"
	echo "--cv-nitersinner=[CV_NITERS_INNER] - Maximum number of inner iterations in ChanVese algorithm (default=1000)"
	echo "--cv-nitersreinit=[CV_NITERS_REINIT] - Maximum number of re-init iterations in ChanVese algorithm (default=1000)"
	echo "--cv-timestep=[CV_TIMESTEP] - Chan-Vese time step parameter (default=0.007)"
	echo "--cv-wsize=[CV_WINDOWSIZE] - Chan-Vese window size parameter (default=1)"
	echo "--cv-lambda1=[CV_LAMBDA1] - Chan-Vese lambda1 parameter (default=1)"
	echo "--cv-lambda2=[CV_LAMBDA2] - Chan-Vese lambda2 parameter (default=2)"	
	echo "--cv-mu=[CV_MU] - Chan-Vese mu parameter (default=0.5)"
	echo "--cv-nu=[CV_NU] - Chan-Vese nu parameter (default=0)"	
	echo "--cv-p=[CV_P] - Chan-Vese p parameter (default=1)"		
	echo ""	

	echo "=== SFINDER WAVELET TRANSFORM FILTER OPTIONS ==="
	echo "--wtscalemin - Minimum Wavelet Transform scale for extended source search (default=3)"
	echo "--wtscalemax - Maximum Wavelet Transform scale for extended source search (default=6)"
	echo ""
	
	echo "=== SFINDER RUN OPTIONS ==="	
	echo "--loglevel=[LOG_LEVEL] - Logging level string {INFO, DEBUG, WARN, ERROR, OFF} (default=INFO)"
	echo "--maxfiles=[NMAX_PROCESSED_FILES] - Maximum number of input files processed in filelist (default=-1=all files)"
	echo "--addrunindex - Append a run index to submission script (in case of list execution) (default=no)"
	echo "--outdir=[OUTPUT_DIR] - Output directory where to put run output file (default=pwd)"
	echo "--no-mpi - Disable MPI run (even with 1 proc) (default=enabled)"
	echo "--mpioptions - Options to be passed to MPI (e.g. --bind-to {none,hwthread, core, l1cache, l2cache, l3cache, socket, numa, board}) (default=)"
	echo "--nproc=[NPROC] - Number of MPI processors per node used (NB: mpi tot nproc=nproc x nnodes) (default=1)"
	echo "--nthreads=[NTHREADS] - Number of threads to be used in OpenMP (default=-1=all available in node)"
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
ENV_FILE=""
SUBMIT=false
BATCH_SYSTEM="PBS"
CONTAINER_IMG=""
CONTAINER_OPTIONS=""
RUN_IN_CONTAINER=false
FILELIST_GIVEN=false
INPUTFILE=""
INPUTFILE_GIVEN=false
APPEND_RUN_INDEX=false
MPI_ENABLED=true
NPROC=1
MPI_OPTIONS=""
HOSTFILE=""
HOSTFILE_GIVEN=false
TILE_SIZE=0
SPLIT_IN_TILES="false"
TILE_OVERLAP="false"
TILE_STEP=1
FILELIST=""
NMAX_PROCESSED_FILES=-1
BATCH_QUEUE=""
JOB_WALLTIME="96:00:00"
JOB_MEMORY="4"
JOB_USER_GROUP=""
JOB_USER_GROUP_OPTION=""
JOB_NNODES="1"
JOB_NCPUS="1"
LOG_LEVEL="INFO"
OUTPUT_DIR=$PWD
NTHREADS=1
READ_TILE="false"
XMIN=0
XMAX=0
YMIN=0
YMAX=0
USE_PARALLEL_MEDIAN_ALGO="true"
USE_LOCAL_BKG="true"
BKG_ESTIMATOR="2"
BKG_BOXSIZE="20"
BKG_GRIDSIZE="0.2"
BKG_USE_2ND_PASS="true"
BKG_SKIP_OUTLIERS="false"
SOURCE_BKGBOX_BORDER="20"

NPIX_MIN="5"
SEED_THR="5"
MERGE_THR="2.6"
SEARCH_COMPACT_SOURCES="true"
COMPACT_SOURCE_SEARCH_NITERS="5"
SEED_THR_STEP="1"

##DILATE_NESTED="false"
RESIDUAL_REMOVE_NESTED="false"
##DILATE_BRIGHT_THR="10"
RESIDUAL_ZHIGHTHR="10"
##DILATE_THR="5"
RESIDUAL_ZTHR="5"
##DILATED_SOURCE="2"
RESIDUAL_REMOVED_SOURCE_TYPE="2"
PS_SUBTRACTION_METHOD="1"
DILATE_KERNEL_SIZE="9"
USE_PRESMOOTHING="true"
SMOOTH_FILTER="2"
GUIDED_FILTER_RADIUS="12"
GUIDED_FILTER_EPS="0.04"
EXT_SFINDER_METHOD="3"
AC_METHOD="2"


SELECT_SOURCES="false"
USE_BOUNDING_BOX_CUT="true"
USE_ELONGATION_CUT="true"
USE_CIRC_RATIO_CUT="true"
USE_ELLIPSE_AREA_RATIO_CUT="false"
USE_NMAXPIX_CUT="false"
USE_NBEAMS_CUT="true"
MIN_BOUNDING_BOX="2"
CIRC_RATIO_THR="0.4"
ELONGATION_THR="0.7"
ELLIPSE_AREA_RATIO_MIN_THR="0.6"
ELLIPSE_AREA_RATIO_MAX_THR="1.4"
MAX_NPIX="1000"
NBEAMS_THR="3"

SEARCH_NESTED_SOURCES="true"
NESTED_SOURCE_TO_BEAM_THR="5"
NESTED_BLOB_THR="0"
NESTED_MIN_MOTHER_DIST="2"
NESTED_MAX_MOTHER_PIX_MATCH="0.5"
NESTED_BLOB_PEAK_ZTHR="5"
NESTED_BLOB_PEAK_ZTHR_MERGE="2.5"
NESTED_BLOB_MIN_SCALE="1"
NESTED_BLOB_MAX_SCALE="3"
NESTED_BLOB_SCALE_STEP="1"
NESTED_BLOB_KERN_FACTOR="1"
BLOB_MASK_METHOD="2"
	
SP_SIZE="20"
SP_BETA="1"
SP_MINAREA="10"
SALIENCY_THR="2.8"
SALIENCY_MIN_RESO="20"
SALIENCY_MAX_RESO="60"
SALIENCY_RESO_STEP="10"
SALIENCY_NN_PAR="1"
SALIENCY_USE_ROBUST_PARS="false"
USE_BKG_MAP_IN_SALIENCY="false"
USE_RMS_MAP_IN_SALIENCY="false"
USE_OPTIMAL_THR_IN_SALIENCY="true"

SEARCH_EXTENDED_SOURCES="true"

FIT_SOURCES="false"
FIT_USETHREADS="false"
FIT_MINIMIZER="Minuit2"
FIT_MINIMIZER_ALGO="minimize"
FIT_PRINTLEVEL="1"
FIT_STRATEGY="2"
FIT_MAX_NBEAMS="20"
FIT_MAX_COMPONENTS="3"
FIT_WITH_FIXED_BKG="true"
FIT_WITH_ESTIMATED_BKG="false"
FIT_USE_BKGBOX_ESTIMATE="false"
FIT_BKG="0"
FIT_WITH_BKG_LIMITS="true"
FIT_WITH_AMPL_LIMITS="true"
FIT_AMPL_LIMIT="0.3"
FIT_WITH_POS_LIMITS="true"
FIT_POS_LIMIT="3"
PREFIT_FIX_POS="true"
FIT_WITH_SIGMA_LIMITS="true"
FIT_SIGMA_LIMIT="0.3"
FIT_WITH_THETA_LIMITS="true"
FIT_THETA_LIMIT="90"
FIT_WITH_SIGMA_FIXED="false"
FIT_WITH_THETA_FIXED="false"
PEAK_MIN_KERN_SIZE="3"
PEAK_MAX_KERN_SIZE="7"
PEAK_KERNEL_MULTIPLICITY_THR="1"
PEAK_SHIFT_TOLERANCE="2"
PEAK_ZTHR_MIN="1"
FIT_USE_NESTED_AS_COMPONENTS="false"
FIT_IMPROVE_CONVERGENCE="true"
FIT_PARBOUNDINCREASE_STEPSIZE="0.1"
FIT_MAXNITERS="10000"
FIT_NRETRIES="10"
FIT_RETRY_WITH_LESS_COMPONENTS="true"
FIT_FCNTOL="1.e-2"
FIT_IMPROVE_ERRORS="false"
FIT_SCALE_DATA_TO_MAX="false"
PREFIT_FIX_AMPL="true"
PREFIT_FIX_SIGMA="false"
PREFIT_FIX_THETA="false"
FIT_USE_CHI2_CUT="true"
FIT_USE_ELLIPSE_CUTS="false"
FIT_CHI2_CUT="5"

WTSCALE_MIN="3"
WTSCALE_MAX="6"
MERGE_EDGE_SOURCES="false"
MERGE_SOURCES="true"
MERGE_COMPACT_SOURCES="true"
MERGE_EXTENDED_SOURCES="true"
SAVE_INPUT_MAP="false"
SAVE_BKG_MAP="false"
SAVE_RMS_MAP="false"
SAVE_SIGNIFICANCE_MAP="false"
SAVE_RESIDUAL_MAP="false"
SAVE_SALIENCY_MAP="false"
SAVE_SEGMENTED_MAP="false"
BMAJ=10
BMIN=5
BPA=0
MAP_PIXSIZE="1.0"
SAVE_DS9REGIONS="false"
CONVERT_DS9REGIONS_TO_WCS="false"
DS9REGION_WCSTYPE="0"

AC_NITERS="1000"
AC_LEVELSET_METHOD="1"
AC_LEVELSET_SIZE="0.1"
AC_TOLERANCE="0.1"

##CV_NITERS="1000"
CV_NITERS_INNER="5"
CV_NITERS_REINIT="5"
CV_TIMESTEP="0.007"
CV_WINDOWSIZE="1"
CV_LAMBDA1="1"
CV_LAMBDA2="2"
CV_MU="0.5"
CV_NU="0"
CV_P="1"


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
		--envfile=*)
    	ENV_FILE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
	
		## OPTIONAL ##	
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

		--outdir=*)
    	OUTPUT_DIR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		--no-mpi*)
    	MPI_ENABLED=false
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
		--nthreads=*)
    	NTHREADS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
			
		--save-inputmap*)
    	SAVE_INPUT_MAP="true"
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
		--save-saliencymap*)
    	SAVE_SALIENCY_MAP="true"
    ;;
		--save-segmentedmap*)
    	SAVE_SEGMENTED_MAP="true"
    ;;
		--save-regions*)
    	SAVE_DS9REGIONS="true"
    ;;
		

		--convertregionstowcs*)
			CONVERT_DS9REGIONS_TO_WCS="true"
		;;
		--regionwcs=*)
			DS9REGION_WCSTYPE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`	
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
		--mappixsize=*)
    	MAP_PIXSIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`		
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
		
		## STATS OPTIONS
		--no-parallelmedian*)
    	USE_PARALLEL_MEDIAN_ALGO="false"
    ;;

		## BKG OPTIONS
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
		--bkgestimator=*)
			BKG_ESTIMATOR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--no-bkg2ndpass*)
			BKG_USE_2ND_PASS="false"
		;;
		--bkgskipoutliers*)
			BKG_SKIP_OUTLIERS="true"
		;;
		--sourcebkgboxborder=*)
			SOURCE_BKGBOX_BORDER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;

		## SOURCE MERGING
		--mergeedgesources*)
    	MERGE_EDGE_SOURCES="true"
    ;;
		--no-mergesources*)
    	MERGE_SOURCES="false"
    ;;
		--no-mergecompactsources*)
    	MERGE_COMPACT_SOURCES="false"
    ;;
		--no-mergeextsources*)
    	MERGE_EXTENDED_SOURCES="false"
    ;;

		## COMPACT SOURCE OPTIONS
		--no-compactsearch*)
    	SEARCH_COMPACT_SOURCES="false"
    ;;
		
		--compactsearchiters=*)
			COMPACT_SOURCE_SEARCH_NITERS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--seedthrstep=*)		
			SEED_THR_STEP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;

		--npixmin=*)
    	NPIX_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		
		--seedthr=*)
    	SEED_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		--mergethr=*)
    	MERGE_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		## SOURCE SELECTION
		--selectsources*)
    	SELECT_SOURCES="true"
    ;;
		--no-boundingboxcut*)
    	USE_BOUNDING_BOX_CUT="false"
    ;;
		--no-circratiocut*)
    	USE_CIRC_RATIO_CUT="false"
    ;;
		--no-elongationcut*)
    	USE_ELONGATION_CUT="false"
    ;;
		--ellipsearearatiocut*)
    	USE_ELLIPSE_AREA_RATIO_CUT="true"
    ;;
		--maxnpixcut*)
    	USE_NMAXPIX_CUT="true"
    ;;
		--no-nbeamscut*)
    	USE_NBEAMS_CUT="false"
    ;;
		
		--minboundingbox=*)
    	MIN_BOUNDING_BOX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--circratiothr=*)
    	CIRC_RATIO_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--elongationthr=*)
    	ELONGATION_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ellipsearearatiominthr=*)
    	ELLIPSE_AREA_RATIO_MIN_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ellipsearearatiomaxthr=*)
    	ELLIPSE_AREA_RATIO_MAX_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--maxnpix=*)
    	MAX_NPIX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;	
		--nbeamsthr=*)
    	NBEAMS_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
	
		## NESTED OPTIONS
		--no-nestedsearch*)
    	SEARCH_NESTED_SOURCES="false"
    ;;
		--blobmaskmethod=*)
			BLOB_MASK_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--nested-sourcetobeamthr=*)
    	NESTED_SOURCE_TO_BEAM_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nested-blobthr=*)
    	NESTED_BLOB_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nested-minmotherdist=*)
    	NESTED_MIN_MOTHER_DIST=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nested-maxmotherpixmatch=*)
    	NESTED_MAX_MOTHER_PIX_MATCH=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--nested-blobpeakzthr=*)
			NESTED_BLOB_PEAK_ZTHR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--nested-blobpeakzthrmerge=*)
			NESTED_BLOB_PEAK_ZTHR_MERGE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--nested-blobminscale=*)
			NESTED_BLOB_MIN_SCALE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--nested-blobmaxscale=*)
			NESTED_BLOB_MAX_SCALE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--nested-blobscalestep=*)
			NESTED_BLOB_SCALE_STEP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--nested-blobkernfactor=*)
			NESTED_BLOB_KERN_FACTOR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;

		## RESIDUAL OPTIONS
		--res-zthr=*)
    	RESIDUAL_ZTHR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--res-zhighthr=*)
    	RESIDUAL_ZHIGHTHR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--res-removenested*)
			RESIDUAL_REMOVE_NESTED="true"		
		;;
		--res-removedsourcetype=*)
			RESIDUAL_REMOVED_SOURCE_TYPE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--dilatekernsize=*)
			DILATE_KERNEL_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--res-pssubtractionmethod=*)
			PS_SUBTRACTION_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;

		## SMOOTHING FILTER OPTIONS
		--no-presmoothing*)
    	USE_PRESMOOTHING="false"
    ;;
		--smoothfilter=*)
			SMOOTH_FILTER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--guidedfilter-radius=*)
			GUIDED_FILTER_RADIUS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--guidedfilter-eps=*)
			GUIDED_FILTER_EPS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		

		## EXTENDED SOURCE OPTIONS
		--no-extendedsearch*)
    	SEARCH_EXTENDED_SOURCES="false"
    ;;

		--extsfinder=*)
    	EXT_SFINDER_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--activecontour=*)
    	AC_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		

		## SOURCE FITTING OPTIONS
		--fitsources*)
    	FIT_SOURCES="true"
    ;;	
		--fit-usethreads*)
    	FIT_USETHREADS="true"
    ;;
		--fit-minimizer=*)
			FIT_MINIMIZER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-minimizeralgo=*)
			FIT_MINIMIZER_ALGO=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-printlevel=*)
			FIT_PRINTLEVEL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-strategy=*)
			FIT_STRATEGY=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-maxnbeams=*)
			FIT_MAX_NBEAMS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-maxcomponents=*)
    	FIT_MAX_COMPONENTS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fit-usenestedascomponents*)
			FIT_USE_NESTED_AS_COMPONENTS="true"
		;;
		--fit-freebkg*)
    	FIT_WITH_FIXED_BKG="false"
    ;;
		--fit-estimatedbkg*)
    	FIT_WITH_ESTIMATED_BKG="true"
    ;;
		--fit-usebkgboxestimate*)
			FIT_USE_BKGBOX_ESTIMATE="true"
		;;
		--fit-bkg=*)
    	FIT_BKG=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fit-nobkglimits*)
    	FIT_WITH_BKG_LIMITS="false"
    ;;
		--fit-ampllimit=*)
    	FIT_AMPL_LIMIT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;	
		--fit-noampllimits*)
    	FIT_WITH_AMPL_LIMITS="false"
    ;;
		--prefit-freeampl*)
    	PREFIT_FIX_AMPL="false"
    ;;	
		--fit-sigmalimit=*)
    	FIT_SIGMA_LIMIT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--prefit-fixsigma*)
    	PREFIT_FIX_SIGMA="true"
    ;;
		--fit-nosigmalimits*)
    	FIT_WITH_SIGMA_LIMITS="false"
    ;;
		--fit-noposlimits*)
    	FIT_WITH_POS_LIMITS="false"
    ;;
		--fit-poslimit=*)
    	FIT_POS_LIMIT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--prefit-freepos*)
			PREFIT_FIX_POS="false"
		;;
		--fit-thetalimit=*)
    	FIT_THETA_LIMIT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--prefit-fixtheta*)
    	PREFIT_FIX_THETA="true"
    ;;
		--fit-nothetalimits*)
    	FIT_WITH_THETA_LIMITS="false"
    ;;	
		--fit-fixsigma*)
    	FIT_WITH_SIGMA_FIXED="true"
    ;;
		--fit-fixtheta*)
    	FIT_WITH_THETA_FIXED="true"
    ;;
		--fit-peakminkern=*)
    	PEAK_MIN_KERN_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fit-peakmaxkern=*)
    	PEAK_MAX_KERN_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fit-peakmultiplicitythr=*)
    	PEAK_KERNEL_MULTIPLICITY_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fit-peakshifttol=*)
    	PEAK_SHIFT_TOLERANCE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--fit-peakzthrmin=*)
    	PEAK_ZTHR_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		--fit-noimproveconvergence*)
			FIT_IMPROVE_CONVERGENCE="false"
		;;
		--fit-noretry*)
			FIT_RETRY_WITH_LESS_COMPONENTS="false"	
		;;
		--fit-improveerrors*)
			FIT_IMPROVE_ERRORS="true"
		;;
		--fit-parboundincreasestep=*)
			FIT_PARBOUNDINCREASE_STEPSIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-maxniters=*)
			FIT_MAXNITERS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-nretries=*)
			FIT_NRETRIES=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-fcntol=*)
			FIT_FCNTOL=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-scaledatatomax*)
			FIT_SCALE_DATA_TO_MAX="true"
		;;
		--fit-nochi2cut*)
			FIT_USE_CHI2_CUT="false"
		;;
		--fit-chi2cut=*)
			FIT_CHI2_CUT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
		;;
		--fit-useellipsecuts*)
			FIT_USE_ELLIPSE_CUTS="true"
		;;


		## SALIENCY FILTER OPTIONS
		--saliency-nooptimalthr*)
			USE_OPTIMAL_THR_IN_SALIENCY="false"
		;;

		--saliency-thr=*)
    	SALIENCY_THR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliency-minreso=*)
    	SALIENCY_MIN_RESO=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliency-maxreso=*)
    	SALIENCY_MAX_RESO=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliency-resostep=*)
    	SALIENCY_RESO_STEP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliency-nn=*)
    	SALIENCY_NN_PAR=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--saliency-usebkgmap*)
    	USE_BKG_MAP_IN_SALIENCY="true"
    ;;
		--saliency-usermsmap*)
    	USE_RMS_MAP_IN_SALIENCY="true"
    ;;
		--saliency-userobustpars*)
			SALIENCY_USE_ROBUST_PARS="true"
		;;

		--sp-size=*)
    	SP_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--sp-beta=*)
    	SP_BETA=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--sp-minarea=*)
    	SP_MINAREA=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		## WAVELET TRANSFORM FILTER OPTIONS
		--wtscalemin=*)
    	WTSCALE_MIN=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--wtscalemax=*)
    	WTSCALE_MAX=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		## ACTIVE CONTOURS OPTIONS
		--ac-niters=*)
    	AC_NITERS=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ac-levelset=*)
    	AC_LEVELSET_METHOD=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ac-levelsetsize=*)
    	AC_LEVELSET_SIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--ac-tolerance=*)
    	AC_TOLERANCE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;

		## CHAN-VESE OPTIONS
		--cv-nitersinner=*)
    	CV_NITERS_INNER=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-nitersreinit=*)
    	CV_NITERS_REINIT=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-timestep=*)
    	CV_TIMESTEP=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-wsize=*)
    	CV_WINDOWSIZE=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-lambda1=*)
    	CV_LAMBDA1=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-lambda2=*)
    	CV_LAMBDA2=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-mu=*)
    	CV_MU=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-nu=*)
    	CV_NU=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
    ;;
		--cv-p=*)
    	CV_P=`echo $item | sed 's/[-a-zA-Z0-9]*=//'`
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



echo ""
echo "*****  PARSED ARGUMENTS ****"
echo "SUBMIT? $SUBMIT, BATCH_SYSTEM=$BATCH_SYSTEM, QUEUE=$BATCH_QUEUE, JOB_CPU: $JOB_NCPUS, JOB_NODES=$JOB_NNODES, JOB_WALLTIME: $JOB_WALLTIME, JOB_MEMORY: $JOB_MEMORY, JOB_USER_GROUP: $JOB_USER_GROUP"
echo "RUN_IN_CONTAINER? $RUN_IN_CONTAINER, CONTAINER_IMG=$CONTAINER_IMG, CONTAINER_OPTIONS=$CONTAINER_OPTIONS"
echo "ENV_FILE: $ENV_FILE"
echo "INPUTFILE: $INPUTFILE"
echo "FILELIST: $FILELIST, NMAX_PROCESSED_FILES: $NMAX_PROCESSED_FILES"
echo "SAVE_INPUT_MAP? $SAVE_INPUT_MAP, SAVE_BKG_MAP: $SAVE_BKG_MAP, SAVE_RMS_MAP? $SAVE_RMS_MAP, SAVE_SIGNIFICANCE_MAP? $SAVE_SIGNIFICANCE_MAP, SAVE_RESIDUAL_MAP: $SAVE_RESIDUAL_MAP"
echo "SAVE_SALIENCY_MAP? $SAVE_SALIENCY_MAP, SAVE_SEGMENTED_MAP? $SAVE_SEGMENTED_MAP"
echo "SPLIT_IN_TILES? $SPLIT_IN_TILES, TILE_SIZE: $TILE_SIZE, TILE_STEP: $TILE_STEP"
echo "NPROC: $NPROC, NTHREADS: $NTHREADS, MPI_ENABLED? $MPI_ENABLED, MPI_OPTIONS: $MPI_OPTIONS"
echo "HOSTFILE_GIVEN? $HOSTFILE_GIVEN, HOSTFILE: $HOSTFILE"
echo "LOG_LEVEL: $LOG_LEVEL"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "USE_PARALLEL_MEDIAN_ALGO: $USE_PARALLEL_MEDIAN_ALGO"
echo "BKG BOX: $BKG_BOXSIZE, GRID: $BKG_GRID_SIZE, BKG_ESTIMATOR: $BKG_ESTIMATOR, BKG_USE_2ND_PASS: $BKG_USE_2ND_PASS, BKG_SKIP_OUTLIERS: $BKG_SKIP_OUTLIERS"
echo "NPIX_MIN: $NPIX_MIN, SEED_THR: $SEED_THR, MERGE_THR: $MERGE_THR, NITERS: $COMPACT_SOURCE_SEARCH_NITERS. SEED_THR_STEP=$SEED_THR_STEP"
echo "DILATE_NESTED? $DILATE_NESTED, DILATE_BRIGHT_THR: $DILATE_BRIGHT_THR, DILATE_THR: $DILATE_THR, DILATED_SOURCE: $DILATED_SOURCE, DILATE_KERNEL_SIZE: $DILATE_KERNEL_SIZE"
echo "USE_PRESMOOTHING? $USE_PRESMOOTHING, SMOOTH_FILTER: $SMOOTH_FILTER"
echo "GUIDED_FILTER PARS: ($GUIDED_FILTER_RADIUS, $GUIDED_FILTER_EPS)"
echo "EXT_SFINDER_METHOD: $EXT_SFINDER_METHOD"
echo "AC_METHOD: $AC_METHOD"
echo "SEARCH_NESTED_SOURCES: $SEARCH_NESTED_SOURCES"
echo "BLOB_MASK_METHOD: $BLOB_MASK_METHOD, NESTED_SOURCE_TO_BEAM_THR: $NESTED_SOURCE_TO_BEAM_THR, NESTED_BLOB_THR: $NESTED_BLOB_THR, NESTED_MIN_MOTHER_DIST: $NESTED_MIN_MOTHER_DIST, NESTED_MAX_MOTHER_PIX_MATCH: $NESTED_MAX_MOTHER_PIX_MATCH, NESTED_BLOB_PEAK_ZTHR: $NESTED_BLOB_PEAK_ZTHR, NESTED_BLOB_PEAK_ZTHR_MERGE: $NESTED_BLOB_PEAK_ZTHR_MERGE"
echo "SELECT_SOURCES: $SELECT_SOURCES"
echo "MERGE_EDGE_SOURCES: $MERGE_EDGE_SOURCES"
echo "FIT_SOURCES: $FIT_SOURCES"
echo "SP PARS: ($SP_SIZE, $SP_BETA, $SP_MINAREA)"
echo "SALIENCY_THR: $SALIENCY_THR"
echo "SALIENCY_RESO: ($SALIENCY_MIN_RESO, $SALIENCY_MAX_RESO, $SALIENCY_RESO_STEP)"
echo "SALIENCY_NN_PAR: $SALIENCY_NN_PAR"
echo "WTSCALE_MIN/WTSCALE_MAX: $WTSCALE_MIN/$WTSCALE_MAX"
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

if [ "$BATCH_SYSTEM" = "" ] && [ "$SUBMIT" = true ]; then
  echo "ERROR: Empty BATCH_SYSTEM argument (hint: you must specify a batch systen if submit option is activated)!"
  exit 1
fi

if [ "$BATCH_SYSTEM" != "PBS" ] && [ "$BATCH_SYSTEM" != "SLURM" ]; then
  echo "ERROR: Unknown/not supported BATCH_SYSTEM argument $BATCH_SYSTEM (hint: PBS/SLURM are supported)!"
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
RUN_IN_CONTAINER_FLAG=""
CONTAINER_IMG_FLAG=""
if [ "$RUN_IN_CONTAINER" = true ]; then
	RUN_IN_CONTAINER_FLAG="--containerrun "
	CONTAINER_IMG_FLAG="--containerimg=$CONTAINER_IMG "	
fi

## Check minimizer options
if [ "$FIT_MINIMIZER" != "Minuit" ] && [ "$FIT_MINIMIZER" != "Minuit2" ]; then
  echo "ERROR: Empty or invalid fit minimizer option given (hint: valid choices are {Minuit, Minuit2})!"
  exit 1
fi
if [ "$FIT_MINIMIZER_ALGO" != "migrad" ] && [ "$FIT_MINIMIZER_ALGO" != "minimize" ] && [ "$FIT_MINIMIZER_ALGO" != "simplex" ] && [ "$FIT_MINIMIZER_ALGO" != "fumili" ] && [ "$FIT_MINIMIZER_ALGO" != "scan" ]; then
  echo "ERROR: Empty or invalid fit minimizer algo option given (hint: valid choices are {migrad, minimize, simplex, scan, fumili})!"
  exit 1
fi


## Compute total number of MPI processor to be given to mpirun
NPROC_TOT=$(($NPROC * $JOB_NNODES))
MPI_ENABLED_OPTION=""
if [ "$MPI_ENABLED" = false ]; then
  MPI_ENABLED_OPTION="--no-mpi"
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
	#BATCH_JOB_NPROC_DIRECTIVE="#PBS -l mpiprocs="
	#BATCH_JOB_MEM_DIRECTIVE="#PBS -l mem="
	#BATCH_JOB_NCORE_DIRECTIVE="#PBS -l ncpus="
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
##         CHECK ENVIRONMENT
#######################################
## Check if CAESAR environment variable exists and is properly set
if [ "$CAESAR_DIR" = "" ] && [ "$RUN_IN_CONTAINER" = false ]; then
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
		echo "pixSize = $MAP_PIXSIZE   								| User-supplied map pixel area in arcsec (pixSize=CDELT, default=1 arcsec)"
		echo "beamFWHM = $BMIN                          | User-supplied circular beam FWHM in arcsec (beamFWHM=BMAJ=BMIN, default=6.5 arcsec)"
		echo "beamBmaj = $BMAJ                        | User-supplied elliptical beam bmaj FWHM in arcsec (default=10 arcsec)"
		echo "beamBmin = $BMIN                        | User-supplied elliptical beam bmin FWHM in arcsec (default=5 arcsec)"
		echo "beamTheta = $BPA                        | User-supplied beam theta in deg (default=0)"
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
		echo "mergeSourcesAtEdge = $MERGE_EDGE_SOURCES           | Merge sources found at tile edge by each workers (default=true)"
		echo "mergeSources = $MERGE_SOURCES								       | Merge overlapping sources found by each workers (default=false)"
		echo "mergeCompactSources = $MERGE_COMPACT_SOURCES       | Merge overlapping compact sources found by each workers (default=false)"
		echo "mergeExtendedSources = $MERGE_EXTENDED_SOURCES     | Merge overlapping extended/compact sources found by each workers (default=false)"
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
    echo 'isInteractiveRun = false                            | Is interactive run (graph plots enabled) (T/F)'
    echo "outputFile = $outputfile                            | Output filename (.root)"
		echo "outputCatalogFile = $catalog_file										| Output catalog ascii filename"
		echo "outputComponentCatalogFile = $catalog_fitcomp_file	| Output fitted component catalog ascii filename"
    echo "ds9RegionFile = $ds9region_file					            | DS9 region file (.reg) where to store source catalog"	
		echo "ds9FitRegionFile = $ds9fitregion_file					      | DS9 region file (.reg) where to store fitted source catalog"
    echo 'ds9RegionFormat = 2                                 | DS9 region format (1=ellipse, 2=polygon)'
		echo "convertDS9RegionsToWCS = $CONVERT_DS9REGIONS_TO_WCS | Convert DS9 regions (contours & ellipses) to WCS (default=false)"
		echo "ds9WCSType = $DS9REGION_WCSTYPE                     | DS9 region WCS output format (0=J2000,1=B1950,2=GAL) (default=0)"
    echo 'inputMapFITSFile = 	input_map.fits				          | Output filename where to store input map in FITS format (.fits)'
    echo 'residualMapFITSFile = residual_map.fits		          | Output filename where to store residual map in FITS format (.fits)'
    echo 'saliencyMapFITSFile = saliency_map.fits		          | Output filename where to store saliency map in FITS format (.fits)'
    echo 'bkgMapFITSFile = bkg_map.fits							          | Output filename where to store bkg map in FITS format (.fits)'
		echo 'noiseMapFITSFile = noise_map.fits					          | Output filename where to store noise map in FITS format (.fits)'
    echo 'significanceMapFITSFile = significance_map.fits	    | Output filename where to store significance map in FITS format (.fits)' 
    echo 'saveToFile = true																	  | Save results & maps to output ROOT file (T/F)'
		echo 'saveToCatalogFile = true														! Save sources to catalog files (island, fitted components) (T/F)'
    echo 'saveToFITSFile = false														  | Save results to output FITS file(s) (T/F)'
		echo "saveDS9Region = $SAVE_DS9REGIONS									  | Save DS9 region files (T/F) (default=T)"
    echo 'saveConfig = true																	  | Save config options to ROOT file (T/F)'
		echo 'saveSources = true																  | Save sources to ROOT file (T/F)'
		echo "saveInputMap = $SAVE_INPUT_MAP										  | Save input map to ROOT file (T/F)"
		echo "saveBkgMap = $SAVE_BKG_MAP												  | Save bkg map to ROOT file (T/F)"
    echo "saveNoiseMap = $SAVE_RMS_MAP											  | Save noise map to ROOT file (T/F)"
    echo "saveResidualMap = $SAVE_RESIDUAL_MAP  						  | Save residual map to ROOT file (T/F)"
    echo "saveSignificanceMap = $SAVE_SIGNIFICANCE_MAP  		  | Save significance map to ROOT file (T/F)"
    echo "saveSaliencyMap = $SAVE_SALIENCY_MAP							  | Save saliency map to ROOT file (T/F)"
    echo "saveSegmentedMap = $SAVE_SEGMENTED_MAP              | Save segmented map computed in extended source search to ROOT file (T/F)"
    echo 'saveEdgenessMap = false                             | Save edgeness map computed in extended source search to ROOT file (T/F)'
    echo 'saveCurvatureMap = false                            | Save curvature map to ROOT file (T/F)'
    echo '###'
    echo '###'
		echo '//==========================='
		echo '//==   STATS OPTIONS         =='
		echo '//==========================='
		echo "useParallelMedianAlgo = $USE_PARALLEL_MEDIAN_ALGO   | Use parallel median algo (based on nth_parallel) (default=true) (T/F)"
		echo '###'
    echo '###'
		echo '//==========================='
		echo '//==   BKG OPTIONS         =='
		echo '//==========================='
    echo "useLocalBkg = $USE_LOCAL_BKG							        | Use local background calculation instead of global bkg (T/F)"
    echo 'localBkgMethod = 1												        | Local background method (1=Grid, 2=Superpixel)'
    echo "use2ndPassInLocalBkg = $BKG_USE_2ND_PASS	        | Use 2nd pass to refine noise calculation in local bkg (T/F)"
		echo "skipOutliersInLocalBkg = $BKG_SKIP_OUTLIERS				| Skip outliers (e.g. bright point sources) in local bkg computation (T/F)"
		echo "bkgEstimator = $BKG_ESTIMATOR                     | Background estimator (1=Mean,2=Median,3=BiWeight,4=ClippedMedian)"
    echo 'useBeamInfoInBkg = true                           | Use beam information in bkg box definition (if available) (T/F)'
		echo "boxSizeX = $BKG_BOXSIZE										        | X Size of local background box in #pixels"
		echo "boxSizeY = $BKG_BOXSIZE										        | Y Size of local background box in #pixels"
		echo "gridSizeX = $BKG_GRIDSIZE									        | X Size of local background grid used for bkg interpolation"
		echo "gridSizeY = $BKG_GRIDSIZE									        | Y Size of local background grid used for bkg interpolation"
		echo "sourceBkgBoxBorderSize = $SOURCE_BKGBOX_BORDER    | Border size in pixels of box around source bounding box used to estimate bkg for fitting"
		echo '###'
		echo '###'
		echo '//==============================='
		echo '//==  FILTERING OPTIONS        =='
		echo '//==============================='
		echo "usePreSmoothing = $USE_PRESMOOTHING								  | Use a pre-smoothing stage to filter input image for extended source search (T/F)"
		echo "smoothFilter = $SMOOTH_FILTER												| Smoothing filter (1=gaus,2=guided)"
		echo 'gausFilterKernSize = 5										          | Gaussian filter kernel size'
		echo 'gausFilterSigma = 1																	| Gaussian filter sigma'
		echo "guidedFilterRadius = $GUIDED_FILTER_RADIUS					| Guided filter radius par"
		echo "guidedFilterColorEps = $GUIDED_FILTER_EPS						| Guided filter color epsilon parameter"
		echo '###'
		echo '###'
		echo '//===================================='
		echo '//==  SOURCE FINDING OPTIONS        =='
		echo '//===================================='
		echo "searchCompactSources = $SEARCH_COMPACT_SOURCES			            | Search compact sources (T/F)"
		echo "minNPix = $NPIX_MIN																	            | Minimum number of pixel to consider a source"
		echo "seedThr = $SEED_THR 																            | Seed threshold in flood filling algo for faint sources"
		echo "mergeThr = $MERGE_THR																            | Merge/aggregation threshold in flood filling algo"
		echo 'mergeBelowSeed = false                                          | Aggregate to seed only pixels above merge threshold but below seed threshold (T/F)'
		echo 'searchNegativeExcess = false												            | Search negative excess together with positive in compact source search'
		echo "compactSourceSearchNIters = $COMPACT_SOURCE_SEARCH_NITERS       | Number of iterations to be performed in compact source search (default=10)"
		echo "seedThrStep = $SEED_THR_STEP                                    | Seed threshold decrease step size between iteration (default=1)"
		echo '###'
		echo '###'
		echo '//==========================================='
		echo '//==  NESTED SOURCE FINDING OPTIONS        =='
		echo '//==========================================='
		echo "searchNestedSources = $SEARCH_NESTED_SOURCES 			               | Search for nested sources inside candidate sources (T/F)"
		echo "blobMaskMethod = $BLOB_MASK_METHOD									             | Blob mask method (1=gaus smooth + laplacian,2=multi-scale LoG)"
		echo "sourceToBeamAreaThrToSearchNested = $NESTED_SOURCE_TO_BEAM_THR   | Source area/beam thr to add nested sources (e.g. npix>thr*beamArea). NB: thr=0 means always if searchNestedSources is enabled (default=0)"
		echo "nestedBlobThrFactor = $NESTED_BLOB_THR                           | Threshold (multiple of curvature rms) used for nested blob finding"
		echo "minNestedMotherDist = $NESTED_MIN_MOTHER_DIST                    | Minimum distance in pixels (in x or y) between nested and parent blob below which nested is skipped"
		echo "maxMatchingPixFraction = $NESTED_MAX_MOTHER_PIX_MATCH            | Maximum fraction of matching pixels between nested and parent blob above which nested is skipped"
		echo "nestedBlobPeakZThr = $NESTED_BLOB_PEAK_ZTHR											 | Nested blob peak significance thr (in curv map) (below thr nested blob is skipped) (default=5 sigmas) "
		echo "nestedBlobPeakZMergeThr = $NESTED_BLOB_PEAK_ZTHR_MERGE        	 | Nested blob significance merge thr (in curv map) (default=2.5 sigmas)"
		echo "nestedBlobMinScale = $NESTED_BLOB_MIN_SCALE                      | Nested blob min search scale (sigma=minscale x beam width) (default=1)"
		echo "nestedBlobMaxScale = $NESTED_BLOB_MAX_SCALE                      | Nested blob max search scale (sigma=maxscale x beam width) (default=3)"
		echo "nestedBlobScaleStep = $NESTED_BLOB_SCALE_STEP                    | Nested blob scale step (scale=minscale + step) (default=1)"
		echo "nestedBlobKernFactor = $NESTED_BLOB_KERN_FACTOR                  | Nested blob curvature/LoG kernel size factor f (kern size=f x sigma) (default=1)"
		echo '###'
		echo '###'
		echo '//=================================='
		echo '//==  SOURCE FITTING OPTIONS      =='
		echo '//=================================='
		echo "fitSources = $FIT_SOURCES                             | Deblend and fit point-like sources with multi-component gaus fit (T/F)"
		echo "fitUseThreads = $FIT_USETHREADS                       | Use multithread in source fitting (default=false) (T/F)"
		echo "fitScaleDataToMax = $FIT_SCALE_DATA_TO_MAX            | Scale source flux data to max peak flux if true, otherwise scale to mJy units (default=false)"
		echo "fitMinimizer = $FIT_MINIMIZER                         | Minimizer {Minuit,Minuit2} (default=Minuit) (T/F)"
		echo "fitMinimizerAlgo = $FIT_MINIMIZER_ALGO         				| Minimizer algorithm: {migrad,simplex,scan,minimize,fumili} (default=minimize)"
		echo "fitPrintLevel = $FIT_PRINTLEVEL                				| Minimizer print level (default=1)"
		echo "fitStrategy = $FIT_STRATEGY                    				| Minimizer strategy (higher means more accurate but more fcn calls) (default=2)"
		echo "nBeamsMaxToFit = $FIT_MAX_NBEAMS							 				| Maximum number of beams in compact source for fitting (if above thr fitting not performed)"
		echo "fitUseNestedAsComponents = $FIT_USE_NESTED_AS_COMPONENTS  | If true use nested sources (if any) to estimate fitted components, otherwise estimate blended blobs (default=false)"
		echo "fitMaxNComponents = $FIT_MAX_COMPONENTS               | Maximum number of components fitted in a blob"
		echo "fitWithCentroidLimits = $FIT_WITH_POS_LIMITS          | Use limits when fitting gaussian centroid (T/F)"
		echo "fitCentroidLimit = $FIT_POS_LIMIT                     | Source centroid limits in pixels "
		echo "fixCentroidInPreFit = $PREFIT_FIX_POS                 | Fix centroid pars in prefit (T/F)"
		echo "fitWithBkgLimits = $FIT_WITH_BKG_LIMITS			          | Use limits when fitting bkg offset (T/F)"
		echo "fitWithFixedBkg = $FIT_WITH_FIXED_BKG                 | Fix bkg level parameter in fit (T/F)"
		echo "fitUseEstimatedBkgLevel = $FIT_WITH_ESTIMATED_BKG     | Use estimated (avg bkg) as bkg level par in fit (T/F)"
		echo "fitUseBkgBoxEstimate = $FIT_USE_BKGBOX_ESTIMATE       | Use bkg estimated in a box around source (if available) as bkg level par in fit (T/F)"
		echo "fitBkgLevel = $FIT_BKG                                | Fixed bkg level used when fitWithFixedBkg=true"
		echo "fitWithAmplLimits = $FIT_WITH_AMPL_LIMITS             | Use limits when fitting gaussian amplitude (T/F)"
		echo "fixAmplInPreFit = $PREFIT_FIX_AMPL                    | Fix amplitude par in pre-fit (T/F)"
		echo "fitAmplLimit = $FIT_AMPL_LIMIT                        | Flux amplitude limit around source peak (e.g. Speak*(1+-fitAmplLimit))"
		echo "fixSigmaInPreFit = $PREFIT_FIX_SIGMA                  | Fix sigma in prefit (T/F)"
		echo "fitWithFixedSigma = $FIT_WITH_SIGMA_FIXED             | Fix sigmas in fit (T/F)"
		echo "fitWithSigmaLimits = $FIT_WITH_SIGMA_LIMITS           | Use limits when fitting gaussian sigmas (T/F)"
		echo "fitSigmaLimit = $FIT_SIGMA_LIMIT                      | Gaussian sigma limit around psf or beam (e.g. Bmaj*(1+-fitSigmaLimit))"
		echo "fitWithFixedTheta = $FIT_WITH_THETA_FIXED             | Fix gaussian ellipse theta par in fit (T/F)"
		echo "fixThetaInPreFit = $PREFIT_FIX_THETA                  | Fix theta in prefit (T/F)"
		echo "fitWithThetaLimits = $FIT_WITH_THETA_LIMITS           | Use limits when fitting gaussian theta par (T/F)"
		echo "fitThetaLimit = $FIT_THETA_LIMIT                      | Gaussian theta limit around psf or beam in degrees (e.g. Bpa +- fitThetaLimit)"
		echo 'useFluxZCutInFit = false                              | Include in fit only source pixels above a given flux significance level (T/F)'
		echo 'fitZCutMin = 2.5                                      | Flux significance below which source pixels are not included in the fit'
		echo "peakMinKernelSize = $PEAK_MIN_KERN_SIZE               | Minimum dilation kernel size (in pixels) used to detect peaks (default=3)"
		echo "peakMaxKernelSize = $PEAK_MAX_KERN_SIZE               | Maximum dilation kernel size (in pixels) used to detect peaks (default=7)"
		echo "peakKernelMultiplicityThr = $PEAK_KERNEL_MULTIPLICITY_THR  | Requested peak multiplicity across different dilation kernels (-1=peak found in all given kernels,1=only in one kernel, etc)"
		echo "peakShiftTolerance = $PEAK_SHIFT_TOLERANCE            | Shift tolerance (in pixels) used to compare peaks in different dilation kernels (default=1 pixel)"
		echo "peakZThrMin = $PEAK_ZTHR_MIN                          | Minimum peak flux significance (in nsigmas above avg source bkg & noise) below which peak is skipped (default=1)"
		echo "fitFcnTolerance = $FIT_FCNTOL								          | Fit function minimization tolerance (default=1.e-5)"
		echo "fitMaxIters = $FIT_MAXNITERS                          | Fit max number of iterations (default=100000)"
    echo "fitImproveConvergence = $FIT_IMPROVE_CONVERGENCE      | Try to improve convergence by iterating fit if not converged or converged with pars at limits (default=true)"
		echo "fitNRetries = $FIT_NRETRIES                           | Number of times fit is repeated (with enlarged limits) if improve convergence flag is enabled (default=1000)"
		echo "fitDoFinalMinimizerStep = $FIT_IMPROVE_ERRORS         | Switch on/off running of final minimizer step after fit convergence with MIGRAD (default=true)"
		echo "###fitFinalMinimizer = 2                              | Final minimizer (1=MIGRAD,2=HESS,3=MINOS) (default=2)"
		echo "fitChi2RegPar = 1															 				| Chi2 regularization par chi2=chi2_signal + regpar*chi2_bkg (default=1)"	
		echo "fitParBoundIncreaseStepSize = $FIT_PARBOUNDINCREASE_STEPSIZE    | Par bound increase step size (e.g. parmax= parmax_old+(1+nretry)*fitParBoundIncreaseStepSize*0.5*|max-min| (default=0.1)"
		echo "fitRetryWithLessComponents = $FIT_RETRY_WITH_LESS_COMPONENTS    | If fit does not converge repeat it iteratively with one component less at each cycle (default=true)"
		echo "fitApplyRedChi2Cut = $FIT_USE_CHI2_CUT	              | Apply fit Chi2/NDF cut. Used to set fit quality flag. If Chi2/NDF>cut the good fit cut is not passed (default=true)"
		echo "fitRedChi2Cut = $FIT_CHI2_CUT									        | Chi2/NDF cut value (default=5)"
		echo "fitApplyFitEllipseCuts = $FIT_USE_ELLIPSE_CUTS        | Apply fit ellipse selection cuts. If not passed fit component is tagged as fake (default=false)"
		echo "fitEllipseEccentricityRatioMinCut = 0.5               | Ellipse eccentricity ratio (fit/beam) min cut value (default=0.5)"
		echo "fitEllipseEccentricityRatioMaxCut = 1.5               | Ellipse eccentricity ratio (fit/beam) max cut value (default=1.5)"
		echo "fitEllipseAreaRatioMinCut = 0.01 							        | Ellipse area ratio (fit/beam) min cut value (default=0.01)"
		echo "fitEllipseAreaRatioMaxCut = 10                        | Ellipse area ratio (fit/beam) min cut value (default=10)"
		echo "fitEllipseRotAngleCut = 45                            | Ellipse rot angle diff (|fit-beam|) cut value in degrees (default=45)"
		echo '###'
		echo '###'
		echo '//================================'
		echo '//==  SOURCE SELECTION OPTIONS  =='
		echo '//================================'
		echo "applySourceSelection = $SELECT_SOURCES  						    | Apply selection cuts to sources (T/F)"
		echo "useMinBoundingBoxCut = $USE_BOUNDING_BOX_CUT				    | Use bounding box cut (T/F)"
		echo "sourceMinBoundingBox = $MIN_BOUNDING_BOX  					    | Minimum bounding box cut (source tagged as bad if below this threshold)"
		echo "useCircRatioCut = $USE_CIRC_RATIO_CUT 							    | Use circularity ratio cut (T/F)"
		echo "psCircRatioThr = $CIRC_RATIO_THR										    | Circular ratio threshold (source passes point-like cut if above this threshold)"
		echo "useElongCut = $USE_ELONGATION_CUT										    | Use elongation cut (T/F)"
		echo "psElongThr = $ELONGATION_THR												    | Elongation threshold (source passes point-like cut if below this threshold"
		echo "useEllipseAreaRatioCut = $USE_ELLIPSE_AREA_RATIO_CUT	  | Use Ellipse area ratio cut (T/F)"
		echo "psEllipseAreaRatioMinThr = $ELLIPSE_AREA_RATIO_MIN_THR	| Ellipse area ratio min threshold"
		echo "psEllipseAreaRatioMaxThr = $ELLIPSE_AREA_RATIO_MAX_THR  | Ellipse area ratio max threshold"
		echo "useMaxNPixCut = $USE_NMAXPIX_CUT  										  | Use max npixels cut (T/F)"
		echo "psMaxNPix = $MAX_NPIX   															  | Max number of pixels for point-like sources (source passes point-like cut if below this threshold)"
		echo "useNBeamsCut = $USE_NBEAMS_CUT                          | Use nBeams cut (T/F)"
		echo "psNBeamsThr = $NBEAMS_THR                               | nBeams threshold (sources passes point-like cut if nBeams<thr)"
		echo '###'
		echo '###'
		echo '//================================'
		echo '//==  SOURCE RESIDUAL OPTIONS   =='
		echo '//================================'
		echo "removeNestedSources = $RESIDUAL_REMOVE_NESTED				| Dilate sources nested inside bright sources (T/F)"
		echo "residualZThr = $RESIDUAL_ZTHR                      | Significance threshold (in sigmas) above which sources of selected type are dilated"
		echo "residualZHighThr = $RESIDUAL_ZHIGHTHR               | Significance threshold (in sigmas) above which sources are always dilated (even if they have nested or different type)"
		echo "dilateKernelSize = $DILATE_KERNEL_SIZE							| Size of kernel (odd) to be used in dilation operation"
		echo "removedSourceType = $RESIDUAL_REMOVED_SOURCE_TYPE   | Type of bright sources to be dilated from the input image (-1=ALL,1=COMPACT,2=POINT-LIKE,3=EXTENDED)"
		echo 'residualModel = 1																    | Model used to replace residual pixel values (1=bkg,2=source median)'
		echo 'residualModelRandomize = false											| Randomize pixel values used to replace residual pixels (T/F)'
		echo "psSubtractionMethod = $PS_SUBTRACTION_METHOD        | Point-source subtraction method (1=dilation, 2=model subtraction (default=1)"
		echo '###'
		echo '###'
		echo '//============================================='
		echo '//==  EXTENDED SOURCE FINDING OPTIONS        =='
		echo '//============================================='
		echo "searchExtendedSources = $SEARCH_EXTENDED_SOURCES		| Search extended sources after bright source removal (T/F)"
		echo "extendedSearchMethod = $EXT_SFINDER_METHOD					| Extended source search method (1=WT-thresholding,2=SPSegmentation,3=ActiveContour,4=Saliency thresholding)"
		echo 'useResidualInExtendedSearch = true									| Use residual image (with selected sources dilated) as input for extended source search'
		echo '###'
		echo '###'
		echo '//================================================'
		echo '//==  WAVELET TRANSFORM ALGORITHM MAIN OPTIONS  =='
		echo '//================================================'
		echo "wtScaleSearchMin = $WTSCALE_MIN									    | Minimum Wavelet scale to be used for extended source search"
		echo "wtScaleSearchMax = $WTSCALE_MAX											| Maximum Wavelet scale to be used for extended source search"		
		echo '###'
		echo '###'
		echo '//=============================================='
		echo '//==  ACTIVE CONTOUR ALGORITHMS MAIN OPTIONS  =='
		echo '//=============================================='
		echo "acMethod = $AC_METHOD										            | Active contour method (1=Chanvese, 2=LRAC)"
		echo "acNIters = $AC_NITERS																| Maximum number of iterations"
		echo "acInitLevelSetMethod = $AC_LEVELSET_METHOD          | Level set initialization method (1=circle,2=checkerboard,3=saliency)"
		echo "acInitLevelSetSizePar = $AC_LEVELSET_SIZE           | Level set size fraction wrt to minimum image size (e.g. circle radius=fraction x image size)"
		echo "acTolerance = $AC_TOLERANCE                         | Tolerance parameter to stop main iteration loop"
		echo '###'
		echo '###'
		echo '//==================================='
		echo '//==  CHAN-VESE ALGORITHM OPTIONS  =='
		echo '//==================================='
		echo "cvNItersInner = $CV_NITERS_INNER                    | Chan-Vese maximum number of inner iterations"
		echo "cvNItersReInit = $CV_NITERS_REINIT                  | Chan-Vese maximum number of iterations performed in re-initialization step"
		echo "cvTimeStepPar = $CV_TIMESTEP												| Chan-Vese time step par"
		echo "cvWindowSizePar = $CV_WINDOWSIZE										| Chan-Vese window size par"
		echo "cvLambda1Par = $CV_LAMBDA1													| Chan-Vese lambda1 par"
		echo "cvLambda2Par = $CV_LAMBDA2													| Chan-Vese lambda2 par"
		echo "cvMuPar = $CV_MU																		| Chan-Vese mu par"
		echo "cvNuPar = $CV_NU																		|	Chan-Vese nu par"
		echo "cvPPar = $CV_P																			| Chan-Vese p par"
		##echo "cvInitContourToSaliencyMap = false                  | Init contour to binarized saliency map"
		echo '###'
		echo '###'
		echo '//==================================='
		echo '//==  LRAC ALGORITHM OPTIONS       =='
		echo '//==================================='
		##echo 'lracNIters = 1000 																	| Number of iterations'
		echo 'lracLambdaPar = 0.1																  | Regularization par'
		echo 'lracRadiusPar = 1																	  | Radius of locatization ball par'
		echo 'lracEpsPar = 0.1																	  | Convergence par'
		##echo 'lracInitContourToSaliencyMap = false                  | Init contour to binarized saliency map'
		echo '###'
		echo '###'
		echo '//==============================='
		echo '//==  SALIENCY FILTER OPTIONS  =='
		echo '//==============================='
		echo "saliencyUseOptimalThr = $USE_OPTIMAL_THR_IN_SALIENCY     | Use optimal threshold in multiscale saliency thresholding (T/F)"
		echo "saliencyThrFactor = $SALIENCY_THR										     | Saliency threshold factor for tagging signal regions (thr=<saliency>*factor)"
		echo 'saliencyBkgThrFactor = 1														     | Saliency threshold factor for tagging bkg regions (thr=<saliency>*factor)'
		echo 'saliencyImgThrFactor = 1														     | Threshold factor to consider a region as significant (thr=<img>*factor)'
		echo "saliencyResoMin = $SALIENCY_MIN_RESO								     | Saliency min reso par"
		echo "saliencyResoMax = $SALIENCY_MAX_RESO								     | Saliency max reso par"
		echo "saliencyResoStep = $SALIENCY_RESO_STEP							     | Saliency reso step par"
		echo 'saliencyUseCurvInDiss = false 											     | Use curvature parameter in dissimilarity estimation (T/F)'
		echo "saliencyUseRobustPars = $SALIENCY_USE_ROBUST_PARS        | Use robust pars in saliency map computation (T/F)"
		echo "saliencyUseBkgMap = $USE_BKG_MAP_IN_SALIENCY    		     | Use bkg map in saliency map computation (T/F)"
		echo "saliencyUseNoiseMap = $USE_RMS_MAP_IN_SALIENCY					 | Use noise map in saliency map computation (T/F)"
		echo "saliencyNNFactor = $SALIENCY_NN_PAR   							     | Fraction of most similar neighbors used in saliency map computation"
		echo '###saliencySpatialRegFactor = 6												     | Spatial regularization factor (ruling exp decay in saliency spatial weighting)'
		echo 'saliencyMultiResoCombThrFactor = 0.7								     | Fraction of combined salient multi-resolution maps to consider global saliency'
		echo 'saliencyDissExpFalloffPar = 100                          | Dissimilarity exponential cutoff parameter (value)'
		echo 'saliencySpatialDistRegPar = 1                            | Spatial-color distance regularization par (value, 1=equal weights)'
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
		echo '###spHierMergingPars = 1  0.3  0.5  3         | Hierarch algo pars: MIN_SEGMENTS, MERGING_RATIO, DIST_REGULARIZATION, DIST_THRESHOLD es 0.25 (value,value,value,value)'
		echo '###spHierMergingMaxDissRatio = 1000  1.05     | Maximum mutual dissimilarity among regions for merging 1st and 2nd neighbors es. 1.15 (value)'
		echo '###spMergingEdgeModel = 2                     | Edge model (1=Kirsch,2=Chan-Vese) (value)'
		echo '###use2ndNeighborsInSPMerging  T              | Use 2nd-order neighbors in superpixel merging (T/F)'
		echo '###useCurvatureInSPMerging  T                 | Use curvature params in superpixel merging (T/F)'
		echo '###useLogContrastInSPGeneration  F            | Use logarithmic contrast to generate initial partition (T/F)'
		echo '###usePixelRatioCut  T  0.4                   | Use pixel ratio to tag background regions (if npix/npixtot>cut)  (T/F,value)'
		echo '###tagSignificantSP  T  1  0.5                | Tag significant superpixels if fraction of significant subregions is above cut (1=saliency,2=Mahalanobis) (T/F,method,cut)'
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
	local logfile=$5
	
	echo "INFO: Creating sh file $shfile (jobindex=$jobindex, exe=$exe, exe_args=$exe_args)..."
	( 
			echo "#!/bin/bash"
			#echo "#PBS -N SFinderJob$jobindex"
			#echo "#PBS -j oe"
  		#echo "#PBS -o $BASEDIR"
			#echo "#PBS -l select=1:ncpus=1:mpiprocs=$NPROC:mem=$JOB_MEMORY"'GB'
    	#echo "#PBS -l walltime=$JOB_WALLTIME"
			#echo "#PBS -l place=scatter"
    	#echo '#PBS -r n'
      #echo '#PBS -S /bin/bash' 
      #echo '#PBS -p 1'
			#echo "$JOB_USER_GROUP_OPTION"

			#echo "$BATCH_JOB_NAME_DIRECTIVE SFinderJob$jobindex"
			#echo "$BATCH_JOB_OUTFILE_DIRECTIVE $BASEDIR"
			#echo "$BATCH_JOB_ERRFILE_DIRECTIVE $BASEDIR"
			#echo "$BATCH_JOB_JOINOUTERR_DIRECTIVE"
			#echo "$BATCH_JOB_WALLTIME_DIRECTIVE$JOB_WALLTIME"
			#echo "$BATCH_JOB_SHELL_DIRECTIVE"
			#echo "$BATCH_JOB_USERGRP_DIRECTIVE $JOB_USER_GROUP"
			#echo "$BATCH_JOB_PRIORITY"
			#echo "$BATCH_JOB_NOREQUEUE_DIRECTIVE"
			#echo "$BATCH_JOB_SCATTER_DIRECTIVE"
			#echo "$BATCH_JOB_NNODES_DIRECTIVE$JOB_NNODES"
			#echo "$BATCH_JOB_NPROC_DIRECTIVE$NPROC"
			#echo "$BATCH_JOB_MEM_DIRECTIVE$JOB_MEMORY"'gb'
			#echo "$BATCH_JOB_NCORE_DIRECTIVE$JOB_NCPUS"

			echo "$BATCH_JOB_NAME_DIRECTIVE SFinderJob$jobindex"
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

      echo " "

      echo 'echo ""'
      echo 'echo "INFO: Source the software environment vars in file '"$ENV_FILE"' ..."'
      echo "source $ENV_FILE"
      echo 'echo ""'
      
      echo "JOBDIR=$BASEDIR"
     
      echo " "
      echo " "

      echo " "
      echo 'echo "*************************************************"'
      echo 'echo "****         RUN SOURCE FINDER               ****"'
      echo 'echo "*************************************************"'
      echo 'echo ""'
      echo '  cd $JOBDIR'

      echo "  $exe $exe_args >& $logfile"
      
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
		outputfile="out-$filename_base_noext"'.root'
 		ds9region_file="ds9-$filename_base_noext"'.reg'
		ds9fitregion_file="ds9_fitcomp-$filename_base_noext"'.reg'
		catalog_file="catalog-$filename_base_noext"'.dat'
		catalog_fitcomp_file="catalog_fitcomp-$filename_base_noext"'.dat'

		## Define output log filename
		logfile="output_$filename_base_noext"'.log'

		## Define config & run script file names 
		if [ "$APPEND_RUN_INDEX" = true ]; then
			configfile="config_$filename_base_noext"'_'"$index.cfg"
			shfile="SFinderRun_$filename_base_noext"'_'"$index.sh"
		else
			configfile="config_$filename_base_noext"'.cfg'
			shfile="SFinderRun_$filename_base_noext"'.sh'
		fi

		
		## Generate config file
  	echo "INFO: Creating config file $configfile for input file: $inputfile ..."
		generate_config $configfile


		## Generate script
		echo "INFO: Creating script file $shfile for input file: $inputfile ..."

		CMD=""
		if [ "$MPI_ENABLED" = true ]; then	
			CMD="mpirun -np $NPROC_TOT $MPI_OPTIONS "
			if [ "$HOSTFILE_GIVEN" = true ] ; then
				CMD="$CMD -f $HOSTFILE "
			fi
		fi
		if [ "$RUN_IN_CONTAINER" = true ] ; then
			EXE="$CMD singularity run $CONTAINER_OPTIONS --app sfinder $CONTAINER_IMG"		
		else
			EXE="$CMD $CAESAR_DIR/bin/FindSourceMPI"
		fi
		EXE_ARGS="--config=$configfile $MPI_ENABLED_OPTION"

		generate_exec_script "$shfile" "$index" "$EXE" "$EXE_ARGS" "$logfile"


		# Submits the job to batch system
		if [ "$SUBMIT" = true ] ; then
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
	catalog_file="catalog-$filename_base_noext"'.dat'
	catalog_fitcomp_file="catalog_fitcomp-$filename_base_noext"'.dat'


	## Define output log filename
	logfile="output_$filename_base_noext"'.log'

	## Define and generate config file
	configfile="config_$filename_base_noext"'.cfg'
  echo "INFO: Creating config file $configfile for input file: $inputfile ..."
	generate_config $configfile


	## Define executable & args variables and generate script
	shfile="Run_$filename_base_noext"'.sh'
	
	##EXE="$CAESAR_DIR/scripts/RunSFinderMPI.sh"
	##EXE_ARGS="--nproc=$NPROC_TOT --config=$configfile $RUN_IN_CONTAINER_FLAG $CONTAINER_IMG_FLAG"
	##if [ "$HOSTFILE_GIVEN" = true ] ; then
	##	EXE_ARGS="$EXE_ARGS --hostfile=$HOSTFILE"
	##fi

	CMD=""
	if [ "$MPI_ENABLED" = true ]; then	
		CMD="mpirun -np $NPROC_TOT $MPI_OPTIONS "
		if [ "$HOSTFILE_GIVEN" = true ] ; then
			CMD="$CMD -f $HOSTFILE "
		fi
	fi

	if [ "$RUN_IN_CONTAINER" = true ] ; then
		EXE="$CMD singularity run $CONTAINER_OPTIONS --app sfinder $CONTAINER_IMG"		
	else
		EXE="$CMD $CAESAR_DIR/bin/FindSourceMPI"
	fi
	EXE_ARGS="--config=$configfile $MPI_ENABLED_OPTION"


	echo "INFO: Creating script file $shfile for input file: $inputfile ..."
	jobId=" "
	generate_exec_script "$shfile" "$jobId" "$EXE" "$EXE_ARGS" "$logfile"

	# Submits the job to batch system
	if [ "$SUBMIT" = true ] ; then
		echo "INFO: Submitting script $shfile to QUEUE $BATCH_QUEUE using $BATCH_SYSTEM batch system ..."
		$BATCH_SUB_CMD $BATCH_QUEUE_NAME_OPTION $BATCH_QUEUE $CURRENTJOBDIR/$shfile
	fi

fi

echo "*** END SUBMISSION ***"

