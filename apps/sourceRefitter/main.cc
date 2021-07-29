// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************

#include <Image.h>
#include <Source.h>
#include <SourceFitter.h>
#include <SourceExporter.h>

#include <ConfigParser.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>

//ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TKey.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>


using namespace std;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input=[INPUT_FILE] \t Input file in ROOT format with Caesar source collection to be refitted"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (ROOT format) where to store selected sources (default=sources.root)"<<endl;
	cout<<"-c, --config \t Config file containing option settings"<<endl;
	cout<<"-f, --fitUnfitted \t Perform fitting of sources found without fit information (default=false)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "config", required_argument, 0, 'c' },
	{ "input", required_argument, 0, 'i' },
	{ "output", required_argument, 0, 'o' },	
	{ "fitUnfitted", no_argument, 0, 'f' },
	{ "verbosity", required_argument, 0, 'v'},
  {(char*)0, (int)0, (int*)0, (int)0}
};

struct SourceIdData
{
	int sourceIndex;
	int nestedSourceIndex;
	SourceIdData(int sid,int nid=-1)
		: sourceIndex(sid), nestedSourceIndex(nid)
	{}
};


//Options
std::string configFileName= "";
std::string fileName= "";
std::string outputFileName= "sources.root";
std::string regionOutputFileName= "sources.reg";
std::string regionOutputFileName_refitted= "sources_refitted.reg";
std::string regionComponentsOutputFileName_refitted= "sources_fitcomp_refitted.reg";
std::string regionComponentsOutputFileName= "sources_fitcomp.reg";
std::string catalogOutputFileName= "catalog.dat";
std::string catalogComponentsOutputFileName= "catalog_fitcomp.dat";
int verbosity= 4;//INFO level
bool fitUnfitted= false;

//Variables
TFile* outputFile= 0;
TTree* outputTree= 0;
TFile* inputFile= 0;
Source* m_source= 0;
std::vector<Source*> m_sources;
std::vector<Source*> m_fittedSources;
std::vector<SourceIdData> m_fittedSourceIds;
int m_ds9WCSType= 0;//use original WCS type to save catalog
bool m_useSimpleWCSEllipseConversion;
double m_nestedBlobMinScale;
double m_nestedBlobMaxScale;
double m_nestedBlobScaleStep;
double m_nestedBlobKernFactor;
double m_NestedBlobThrFactor;
int m_NMinPix;
SourceFitOptions m_fitOptions;
double m_nBeamsMaxToFit;
double m_beamBmaj;
double m_beamBmin;
double m_beamBpa;
double m_pixSize;
double m_pixSizeX;
double m_pixSizeY;
int m_fitMaxNComponents;
bool m_fitWithCentroidLimits;
bool m_fixCentroidInPreFit;
double m_fitCentroidLimit;
bool m_fitWithFixedBkg;
bool m_fitWithBkgLimits;
double m_fitBkgLevel;	
bool m_fitUseEstimatedBkgLevel;
bool m_fitWithAmplLimits;
bool m_fixAmplInPreFit;
double m_fitAmplLimit;
bool m_fixSigmaInPreFit;
bool m_fitWithSigmaLimits;
double m_fitSigmaLimit;
bool m_fitWithFixedSigma;
bool m_fitWithFixedTheta;
bool m_fitWithThetaLimits;
bool m_fixThetaInPreFit;
double m_fitThetaLimit;
bool m_useFluxZCutInFit;
double m_fitZCutMin;
int m_peakMinKernelSize;
int m_peakMaxKernelSize;
int m_peakKernelMultiplicityThr;
int m_peakShiftTolerance;
double m_peakZThrMin;
double m_fitFcnTolerance;
long int m_fitMaxIters;
bool m_fitImproveConvergence;
long int m_fitNRetries;
bool m_fitDoFinalMinimizerStep;
int m_fitFinalMinimizer;
bool m_fitUseNestedAsComponents;
double m_fitChi2RegPar;
bool m_fitRetryWithLessComponents;
std::string m_fitMinimizer;		
std::string m_fitMinimizerAlgo;
int m_fitStrategy;
int m_fitPrintLevel;	
double m_fitParBoundIncreaseStepSize;
bool m_fitUseThreads;
bool fitInMultithread;
bool m_fitScaleDataToMax;
int m_sourceBkgBoxBorderSize;
bool m_fitUseBkgBoxEstimate;
bool m_fitApplyRedChi2Cut;
double m_fitRedChi2Cut;
bool m_fitApplyFitEllipseCuts;
double m_fitEllipseEccentricityRatioMinCut;
double m_fitEllipseEccentricityRatioMaxCut;
double m_fitEllipseAreaRatioMinCut;
double m_fitEllipseAreaRatioMaxCut;
double m_fitEllipseRotAngleCut;


//Functions
int ParseOptions(int argc, char *argv[]);
int Init();
std::string GetStringLogLevel(int verbosity);
int ReadData(std::string filename);
bool IsFittableSource(Source* source);
int FitSources(std::vector<Source*>& sources);
int FitSource(std::vector<SourceIdData>& fittedSourceIds,Source* source,SourceFitOptions& fitOptions,int sindex,long int nindex=-1);
int CloneObjectsInFile(std::vector<std::string> excludedObjNames);
void Save();
int SaveSources();
int SaveDS9Regions();
int SaveCatalog();

int main(int argc, char *argv[])
{
	//================================
	//== PRINT LOGO
	//================================
	SysUtils::PrintAsciiLogo();

	//================================
	//== PARSE OPTIONS
	//================================
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to parse command line options!");
		#endif
		return -1;
	}
	
	//=======================
	//== INIT
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Initializing data...");
	#endif
	if(Init()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to initialize data!");
		#endif
		return -1;
	}

	//=======================
	//==  READ DATA
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading source data in file "<<fileName<<" ...");
	#endif
	if(ReadData(fileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of source data failed!");
		#endif
		return -1;
	}


	//=======================
	//==  FIT SOURCES
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Fitting sources in file "<<fileName<<" ...");
	#endif
	if(FitSources(m_sources)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Source re-fitting failed!");
		#endif
		return -1;
	}

	//=======================
	//== SAVE
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving data to file ...");
	#endif
	Save();

	#ifdef LOGGING_ENABLED
		INFO_LOG("End source refitter");
	#endif

	return 0;

}//close main



int ParseOptions(int argc, char *argv[])
{
	//## Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments...see macro usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Parse options
	std::string configFileName= "";
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hc:i:o:v:f",options_tab, &option_index)) != -1) {
    
    switch (c) {
			case 0 : 
			{
				break;
			}
			case 'h':
			{
      	Usage(argv[0]);	
				exit(0);
			}
			case 'i':	
			{
				fileName= std::string(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
    	case 'c':	
			{
				configFileName= std::string(optarg);	
				break;	
			}
			case 'f':	
			{
				fitUnfitted= true;
				break;
		 	}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			default:
			{
      	Usage(argv[0]);	
				exit(0);
			}
    }//close switch
	}//close while
 
	//=======================
	//== CHECK ARGS 
	//=======================
	//Check input file name
	if(fileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty input file name given!");
		#endif
		return -1;
	}

	//Check output file name
	if(outputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output file name given!");
		#endif
		return -1;
	}

	//Check region output file name
	if(regionOutputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty region output file name given!");
		#endif
		return -1;
	}

	//Check catalog output file name
	if(catalogOutputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty catalog output file name given!");
		#endif
		return -1;
	}

	//Set DS9 region component file name
	regionComponentsOutputFileName= CodeUtils::ExtractFileNameFromPath(regionOutputFileName,true);
	regionComponentsOutputFileName+= "_fitcomp.reg";

	//Set catalog component file name
	catalogComponentsOutputFileName= CodeUtils::ExtractFileNameFromPath(catalogOutputFileName,true);
	catalogComponentsOutputFileName+= "_fitcomp.dat";


	//## Read config options 
	if(ConfigParser::Instance().Parse(configFileName)<0){
		cerr<<"ERROR: Failed to parse config options!"<<endl;
		return -1;
	}
	cout<<endl;
	cout<<"== RUN CONFIG OPTIONS =="<<endl;
	PRINT_OPTIONS();
	cout<<"========================"<<endl;
	cout<<endl;


	//=======================
	//== Init Logger 
	//=======================
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	//=======================
	//== Init thread numbers 
	//=======================
	int nThreads;
	if(GET_OPTION_VALUE(nThreads,nThreads)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get nThreads option!");
		#endif
		return -1;
	}
	if(nThreads>0) SysUtils::SetOMPThreads(nThreads);

	return 0;

}//close ParseOptions()

int Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");
	
	//Init output source TTree
	if(!outputTree) outputTree= new TTree("SourceInfo","SourceInfo");
	m_source= 0;
	m_sources.clear();
	outputTree->Branch("Source",&m_source);
	
	//## Set options
	//Get user-supplied map & beam options
	GET_OPTION_VALUE(pixSize,m_pixSize);
	GET_OPTION_VALUE(beamBmaj,m_beamBmaj);
	GET_OPTION_VALUE(beamBmin,m_beamBmin);
	GET_OPTION_VALUE(beamTheta,m_beamBpa);
	m_pixSizeX= m_pixSize;
	m_pixSizeY= m_pixSize;

	//Get output file
	GET_OPTION_VALUE(ds9WCSType,m_ds9WCSType);
	GET_OPTION_VALUE(useSimpleWCSEllipseConversion,m_useSimpleWCSEllipseConversion);
	
	//Get nested source options
	GET_OPTION_VALUE(nestedBlobMinScale,m_nestedBlobMinScale);
	GET_OPTION_VALUE(nestedBlobMaxScale,m_nestedBlobMaxScale);
	GET_OPTION_VALUE(nestedBlobScaleStep,m_nestedBlobScaleStep);
	GET_OPTION_VALUE(nestedBlobKernFactor,m_nestedBlobKernFactor);
	GET_OPTION_VALUE(nestedBlobThrFactor,m_NestedBlobThrFactor);

	if(m_nestedBlobMinScale>m_nestedBlobMaxScale){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid nested blob search scales given (hint: min scale cannot be larger than max scale)!");
		#endif
		return -1;
	}
	if(m_nestedBlobScaleStep<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid nested blob scale step given (hint: step must be >0)!");
		#endif
		return -1;
	}

	//Get source search options
	GET_OPTION_VALUE(minNPix,m_NMinPix);

	//Get source fitting options
	GET_OPTION_VALUE(nBeamsMaxToFit,m_nBeamsMaxToFit);
	GET_OPTION_VALUE(fitMaxNComponents,m_fitMaxNComponents);
	GET_OPTION_VALUE(fitWithCentroidLimits,m_fitWithCentroidLimits);
	GET_OPTION_VALUE(fixCentroidInPreFit,m_fixCentroidInPreFit);
	GET_OPTION_VALUE(fitCentroidLimit,m_fitCentroidLimit);
	GET_OPTION_VALUE(fitWithBkgLimits,m_fitWithBkgLimits);
	GET_OPTION_VALUE(fitWithFixedBkg,m_fitWithFixedBkg);
	GET_OPTION_VALUE(fitUseEstimatedBkgLevel,m_fitUseEstimatedBkgLevel);
	GET_OPTION_VALUE(fitBkgLevel,m_fitBkgLevel);
	GET_OPTION_VALUE(fitWithAmplLimits,m_fitWithAmplLimits);
	GET_OPTION_VALUE(fixAmplInPreFit,m_fixAmplInPreFit);
	GET_OPTION_VALUE(fitAmplLimit,m_fitAmplLimit);
	GET_OPTION_VALUE(fitWithSigmaLimits,m_fitWithSigmaLimits);
	GET_OPTION_VALUE(fixSigmaInPreFit,m_fixSigmaInPreFit);
	GET_OPTION_VALUE(fitSigmaLimit,m_fitSigmaLimit);
	GET_OPTION_VALUE(fitWithFixedSigma,m_fitWithFixedSigma);
	GET_OPTION_VALUE(fitWithThetaLimits,m_fitWithThetaLimits);
	GET_OPTION_VALUE(fixThetaInPreFit,m_fixThetaInPreFit);
	GET_OPTION_VALUE(fitWithFixedTheta,m_fitWithFixedTheta);
	GET_OPTION_VALUE(fitThetaLimit,m_fitThetaLimit);
	GET_OPTION_VALUE(useFluxZCutInFit,m_useFluxZCutInFit);
	GET_OPTION_VALUE(fitZCutMin,m_fitZCutMin);
	GET_OPTION_VALUE(peakMinKernelSize,m_peakMinKernelSize);
	GET_OPTION_VALUE(peakMaxKernelSize,m_peakMaxKernelSize);
	GET_OPTION_VALUE(peakKernelMultiplicityThr,m_peakKernelMultiplicityThr);
	GET_OPTION_VALUE(peakShiftTolerance,m_peakShiftTolerance);	
	GET_OPTION_VALUE(peakZThrMin,m_peakZThrMin);
	
	GET_OPTION_VALUE(fitFcnTolerance,m_fitFcnTolerance);
	GET_OPTION_VALUE(fitMaxIters,m_fitMaxIters);
	GET_OPTION_VALUE(fitImproveConvergence,m_fitImproveConvergence);
	GET_OPTION_VALUE(fitNRetries,m_fitNRetries);
	GET_OPTION_VALUE(fitDoFinalMinimizerStep,m_fitDoFinalMinimizerStep);
	GET_OPTION_VALUE(fitUseNestedAsComponents,m_fitUseNestedAsComponents);
	GET_OPTION_VALUE(fitChi2RegPar,m_fitChi2RegPar);
		
	GET_OPTION_VALUE(fitMinimizer,m_fitMinimizer);
	GET_OPTION_VALUE(fitMinimizerAlgo,m_fitMinimizerAlgo);
	GET_OPTION_VALUE(fitStrategy,m_fitStrategy);
	GET_OPTION_VALUE(fitPrintLevel,m_fitPrintLevel);
	GET_OPTION_VALUE(fitParBoundIncreaseStepSize,m_fitParBoundIncreaseStepSize);
	GET_OPTION_VALUE(fitUseThreads,m_fitUseThreads);
	GET_OPTION_VALUE(fitScaleDataToMax,m_fitScaleDataToMax);

	GET_OPTION_VALUE(sourceBkgBoxBorderSize,m_sourceBkgBoxBorderSize);
	GET_OPTION_VALUE(fitUseBkgBoxEstimate,m_fitUseBkgBoxEstimate);
	GET_OPTION_VALUE(fitRetryWithLessComponents,m_fitRetryWithLessComponents);

	GET_OPTION_VALUE(fitApplyRedChi2Cut,m_fitApplyRedChi2Cut);
	GET_OPTION_VALUE(fitRedChi2Cut,m_fitRedChi2Cut);
	GET_OPTION_VALUE(fitApplyFitEllipseCuts,m_fitApplyFitEllipseCuts);
	GET_OPTION_VALUE(fitEllipseEccentricityRatioMinCut,m_fitEllipseEccentricityRatioMinCut);
	GET_OPTION_VALUE(fitEllipseEccentricityRatioMaxCut,m_fitEllipseEccentricityRatioMaxCut);
	GET_OPTION_VALUE(fitEllipseAreaRatioMinCut,m_fitEllipseAreaRatioMinCut);
	GET_OPTION_VALUE(fitEllipseAreaRatioMaxCut,m_fitEllipseAreaRatioMaxCut);
	GET_OPTION_VALUE(fitEllipseRotAngleCut,m_fitEllipseRotAngleCut);
	if(m_fitApplyFitEllipseCuts && m_fitEllipseEccentricityRatioMinCut>=m_fitEllipseEccentricityRatioMaxCut){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid fit ellipse eccentricity ratio cut option given (hint: min cut value must be smaller than max cut value)!");
		#endif
		return -1;
	}
	if(m_fitApplyFitEllipseCuts && m_fitEllipseAreaRatioMinCut>=m_fitEllipseAreaRatioMaxCut){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid fit ellipse area ratio cut option given (hint: min cut value must be smaller than max cut value)!");
		#endif
		return -1;
	}	

	if(m_peakMinKernelSize>m_peakMaxKernelSize){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid peak kernel size option given (hint: min kernel must be larger or equal to max kernel size)!");
		#endif
		return -1;
	}
	if(m_peakMinKernelSize<=0 || m_peakMinKernelSize%2==0 || m_peakMaxKernelSize<=0 || m_peakMaxKernelSize%2==0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid peak kernel sizes given (hint: kernel size must be positive and odd)!");
		#endif
		return -1;
	}

	//Set fit options
	m_fitOptions.bmaj= fabs(m_beamBmaj/m_pixSizeX);//converted in pixels
	m_fitOptions.bmin= fabs(m_beamBmin/m_pixSizeY);//converted in pixels
	m_fitOptions.bpa= m_beamBpa + 90.;//beam pos angle is computed from North. We are assuming angles computed from x axis
	m_fitOptions.nMaxComponents= m_fitMaxNComponents;
	m_fitOptions.limitCentroidInFit= m_fitWithCentroidLimits;
	m_fitOptions.fixCentroidInPreFit= m_fixCentroidInPreFit;
	m_fitOptions.centroidLimit= m_fitCentroidLimit;
	m_fitOptions.limitBkgInFit= m_fitWithBkgLimits;
	m_fitOptions.fixBkg= m_fitWithFixedBkg;
	m_fitOptions.useEstimatedBkgLevel= m_fitUseEstimatedBkgLevel;
	m_fitOptions.useBkgBoxEstimate= m_fitUseBkgBoxEstimate;
	m_fitOptions.fixedBkgLevel= m_fitBkgLevel;
	m_fitOptions.limitAmplInFit= m_fitWithAmplLimits;
	m_fitOptions.fixAmplInPreFit= m_fixAmplInPreFit;
	m_fitOptions.amplLimit= m_fitAmplLimit;
	m_fitOptions.fixSigmaInPreFit= m_fixSigmaInPreFit;
	m_fitOptions.limitSigmaInFit= m_fitWithSigmaLimits;
	m_fitOptions.sigmaLimit= m_fitSigmaLimit;
	m_fitOptions.fixSigma= m_fitWithFixedSigma;
	m_fitOptions.limitThetaInFit= m_fitWithThetaLimits;
	m_fitOptions.fixThetaInPreFit= m_fixThetaInPreFit;
	m_fitOptions.fixTheta= m_fitWithFixedTheta;
	m_fitOptions.thetaLimit= m_fitThetaLimit;
	m_fitOptions.useFluxZCut= m_useFluxZCutInFit;
	m_fitOptions.fluxZThrMin= m_fitZCutMin;
	m_fitOptions.peakMinKernelSize= m_peakMinKernelSize;
	m_fitOptions.peakMaxKernelSize= m_peakMaxKernelSize;
	m_fitOptions.peakZThrMin= m_peakZThrMin;
	m_fitOptions.peakKernelMultiplicityThr= m_peakKernelMultiplicityThr;
	m_fitOptions.peakShiftTolerance= m_peakShiftTolerance;

	m_fitOptions.fitFcnTolerance= m_fitFcnTolerance;
	m_fitOptions.fitMaxIters= m_fitMaxIters;
	m_fitOptions.fitImproveConvergence= m_fitImproveConvergence;
	m_fitOptions.fitNRetries= m_fitNRetries;
	m_fitOptions.fitDoFinalMinimizerStep= m_fitDoFinalMinimizerStep;
	//m_fitOptions.fitFinalMinimizer= m_fitFinalMinimizer;
	m_fitOptions.useNestedAsComponents= m_fitUseNestedAsComponents;
	m_fitOptions.chi2RegPar= m_fitChi2RegPar;
	m_fitOptions.fitRetryWithLessComponents= m_fitRetryWithLessComponents;

	m_fitOptions.useRedChi2Cut= m_fitApplyRedChi2Cut;
	m_fitOptions.fitRedChi2Cut= m_fitRedChi2Cut;
	m_fitOptions.useFitEllipseCuts= m_fitApplyFitEllipseCuts;
	m_fitOptions.fitEllipseEccentricityRatioMinCut= m_fitEllipseEccentricityRatioMinCut;
	m_fitOptions.fitEllipseEccentricityRatioMaxCut= m_fitEllipseEccentricityRatioMaxCut;
	m_fitOptions.fitEllipseAreaRatioMinCut= m_fitEllipseAreaRatioMinCut;
	m_fitOptions.fitEllipseAreaRatioMaxCut= m_fitEllipseAreaRatioMaxCut;
	m_fitOptions.fitEllipseRotAngleCut= m_fitEllipseRotAngleCut;

	m_fitOptions.fitMinimizer= m_fitMinimizer;		
	m_fitOptions.fitMinimizerAlgo= m_fitMinimizerAlgo;
	m_fitOptions.fitStrategy= m_fitStrategy;
	m_fitOptions.fitPrintLevel= m_fitPrintLevel;
	m_fitOptions.fitParBoundIncreaseStepSize= m_fitParBoundIncreaseStepSize;
	m_fitOptions.fitScaleDataToMax= m_fitScaleDataToMax;
	m_fitOptions.wcsType= m_ds9WCSType;
	m_fitOptions.useSimpleWCSEllipseConversion= m_useSimpleWCSEllipseConversion;

	//## Check minimizer support
	if(m_fitOptions.fitMinimizer=="Minuit2" || m_fitOptions.fitMinimizer=="minuit2"){
		#ifndef MINUIT2_ENABLED
			#ifdef LOGGING_ENABLED
				WARN_LOG("Minuit2 minimizer was selected as option but not available/found in the system, switching to Minuit+Migrad as fallback.");
			#endif
			m_fitOptions.fitMinimizer= "Minuit";
			m_fitOptions.fitMinimizerAlgo= "Migrad";
		#endif
	}
	if(m_fitOptions.fitMinimizer=="R" || m_fitOptions.fitMinimizer=="r"){
		#ifndef ROOTR_ENABLED
			#ifdef LOGGING_ENABLED
				WARN_LOG("R minimizer was selected as option but not available/found in the system, switching to Minuit+Migrad as fallback.");
			#endif
			m_fitOptions.fitMinimizer= "Minuit";
			m_fitOptions.fitMinimizerAlgo= "Migrad";
		#endif
	}

	//## Check fit minimizer multithread support
	bool fitInMultithread= m_fitUseThreads;
	if(m_fitUseThreads && (m_fitOptions.fitMinimizer=="Minuit" || m_fitOptions.fitMinimizer=="minuit")){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Selected Minuit minimizer is not thread-safe, switching off source fit multithread.");
		#endif
		fitInMultithread= false;
	}
	if(m_fitUseThreads && (m_fitOptions.fitMinimizer=="R" || m_fitOptions.fitMinimizer=="r")){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Selected R minimizer is not thread-safe, switching off source fit multithread.");
		#endif
		fitInMultithread= false;
	}
	#ifdef LOGGING_ENABLED
		INFO_LOG("Fitting sources in multithread? "<<fitInMultithread);
	#endif

	//## NB: Convert scale pars in pixels assuming they represent multiple of beam width (Bmin)	
	//double pixSize= fabs(std::min(m_pixSizeX,m_pixSizeY));
	double pixSize= std::min(fabs(m_pixSizeX),fabs(m_pixSizeY));
	//double beamWidth= fabs(std::min(m_beamBmaj,m_beamBmin));	
	double beamWidth= std::min(fabs(m_beamBmaj),fabs(m_beamBmin));
	double beamPixSize= beamWidth/pixSize;
	double sigmaMin= m_nestedBlobMinScale*beamPixSize/GausSigma2FWHM;//convert from FWHM to sigma
	double sigmaMax= m_nestedBlobMaxScale*beamPixSize/GausSigma2FWHM;//convert from FWHM to sigma
	double sigmaStep= m_nestedBlobScaleStep;	
	m_fitOptions.scaleMin= sigmaMin;
	m_fitOptions.scaleMax= sigmaMax;
	m_fitOptions.scaleStep= sigmaStep;
	m_fitOptions.minBlobSize= m_NMinPix;
	m_fitOptions.blobMapThrFactor= m_NestedBlobThrFactor;
	m_fitOptions.blobMapKernelFactor= m_nestedBlobKernFactor; 
	
	return 0;

}//close Init()

bool IsFittableSource(Source* aSource)
{
	//Check if not point-like or compact
	int sourceType= aSource->Type;
	bool isCompact= (sourceType==ePointLike || sourceType==eCompact);
	if(!isCompact) return false;

	//If compact source check nbeams (if too large do not perform fit)
	if(sourceType==eCompact){
		double NPix= aSource->GetNPixels();
		double beamArea= aSource->GetBeamFluxIntegral();
		double nBeams= 0;
		if(beamArea>0) nBeams= NPix/beamArea;
		if(nBeams>m_nBeamsMaxToFit) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source "<<aSource->GetName()<<" not fittable as a whole (nBeams="<<nBeams<<">"<<m_nBeamsMaxToFit<<")");
			#endif
			return false;
		}
	}

	return true;

}//close IsFittableSource()

int FitSources(std::vector<Source*>& sources)
{
	//Check given source list
	if(sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty source list, nothing to be fitted!");
		#endif
		return 0;
	}

	//## Loop over image sources and perform fitting stage for non-extended sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Loop over image sources and perform fitting stage for non-extended sources...");
	#endif

	m_fittedSourceIds.clear();
	
	#ifdef OPENMP_ENABLED
		#pragma omp parallel for if(fitInMultithread)
	#endif
	for(size_t i=0;i<sources.size();i++)
	{
		int sindex= (int)(i);
		if(FitSource(m_fittedSourceIds,sources[i],m_fitOptions,sindex,-1)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to fit source "<<sources[i]->GetName()<<", skip to next ...");
			#endif
			continue;
		}
	
	}//end loop sources
	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Fitted #"<<m_fittedSourceIds.size()<<"/"<<sources.size()<<" sources at this stage...");
	#endif

	
	for(size_t i=0;i<m_fittedSourceIds.size();i++){
		int sid= m_fittedSourceIds[i].sourceIndex;	
		int nid= m_fittedSourceIds[i].nestedSourceIndex;	
		Source* aSource= sources[sid];
		if(nid!=-1) aSource= sources[sid]->GetNestedSource(nid);
		m_fittedSources.push_back(aSource);
	}

	//== DEBUG ==
	std::stringstream ss;
	ss<<"fitted sources {";
	for(size_t i=0;i<m_fittedSourceIds.size();i++){	
		int sid= m_fittedSourceIds[i].sourceIndex;	
		int nid= m_fittedSourceIds[i].nestedSourceIndex;	
		Source* aSource= sources[sid];
		if(nid!=-1) aSource= sources[sid]->GetNestedSource(nid);
		std::string sname= aSource->GetName();
		ss<<sname<<",";
	}
	ss<<"}";
	#ifdef LOGGING_ENABLED
		INFO_LOG(ss.str());
	#endif

	return 0;

}//close FitSources()


int FitSource(std::vector<SourceIdData>& fittedSourceIds, Source* source, SourceFitOptions& fitOptions,int sindex,long int nindex)
{
	//Retrieve existing fit information
	bool hasFitInfo= source->HasFitInfo();
	int fitStatus= source->GetFitStatus();

	//Check if source is fittable
	bool isFittable= IsFittableSource(source);

	//If no fit information is available, try to fit if fittable otherwise go to nested
	if(!hasFitInfo){
		if(isFittable && fitUnfitted) {
			//Fit mother source
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source "<<source->GetName()<<") is fittable but is found not fitted in input catalog, it will be fitted now ...");
			#endif

			fittedSourceIds.push_back(SourceIdData(sindex,nindex));

			if(source->Fit(fitOptions)<0) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to fit source "<<source->GetName()<<"), skip to next...");
				#endif
				return -1;
			}
		}//close if fittable
		else{//fit nested

			//Fit nested sources
			std::vector<Source*> nestedSources= source->GetNestedSources();	
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source "<<source->GetName()<<" not fittable as a whole (extended or large compact), fitting nested components individually (#"<<nestedSources.size()<<" components present) ...");
			#endif

			for(size_t j=0;j<nestedSources.size();j++){
				int nestedSourceIndex= (int)(j);
				if(!nestedSources[j]) continue;

				if(FitSource(fittedSourceIds,nestedSources[j],fitOptions,sindex,nestedSourceIndex)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to fit nested source no. "<<j<<" of source "<<source->GetName()<<"), skip to next nested...");
					#endif
					continue;
				}
			}//end loop nested sources
		}//close else fit nested

	}//close if !hasFitInfo
	else{

		//Fit information available, retrieve number of components
		int nComponents= source->GetNFitComponents();
		int nSelComponents= source->GetNSelFitComponents();

		//Refit source if there is a mismatch between nComponents and nSelComponents
		//due to some source components tagged as bad/spurious in the selection process
		if(nComponents!=nSelComponents){
			
			if(isFittable) {
				//Update fit options forcing max components to number of selected components
				SourceFitOptions fitOptions_new= fitOptions;
				fitOptions_new.nMaxComponents= nSelComponents;

				//Fit mother source
				#ifdef LOGGING_ENABLED
					INFO_LOG("Source "<<source->GetName()<<") will be refitted as nComponents="<<nComponents<<"!="<<"nSelComponents="<<nSelComponents<<" ...");
				#endif

				fittedSourceIds.push_back(SourceIdData(sindex,nindex));

				if(source->Fit(fitOptions_new)<0) {
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to fit source "<<source->GetName()<<"), skip to next...");
					#endif
					return -1;
				}

			}//close if fittable
			else{

				//Fit nested sources
				std::vector<Source*> nestedSources= source->GetNestedSources();	
				#ifdef LOGGING_ENABLED
					INFO_LOG("Source "<<source->GetName()<<" not fittable as a whole (extended or large compact), fitting nested components individually (#"<<nestedSources.size()<<" components present) ...");
				#endif

				for(size_t j=0;j<nestedSources.size();j++){
					int nestedSourceIndex= (int)(j);
					if(!nestedSources[j]) continue;

					if(FitSource(fittedSourceIds,nestedSources[j],fitOptions,sindex,nestedSourceIndex)<0){
						#ifdef LOGGING_ENABLED
							WARN_LOG("Failed to fit nested source no. "<<j<<" of source "<<source->GetName()<<"), skip to next nested...");
						#endif
						continue;
					}
				}//end loop nested sources

			}//close else !fittable
		}//close if
	}//close else Refit source
	
	
	return 0;

}//close FitSource()


int ReadData(std::string filename)
{
	//Open file
	inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filename<<"!");
		#endif
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)inputFile->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get access to source tree in file "<<filename<<"!");	
		#endif
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	
	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Found #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif

	for(int i=0;i<sourceTree->GetEntries();i++)
	{
		sourceTree->GetEntry(i);
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif

		//Copy source
		Source* source= new Source;
		*source= *aSource;

		//Add sources to list
		m_sources.push_back(source);

	}//end loop sources

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<m_sources.size()<<" sources read ...");
	#endif

	return 0;

}//close ReadData()


int SaveDS9Regions()
{
	//Save DS9 regions for islands
	bool convertDS9RegionsToWCS= false;
	int status= SourceExporter::WriteToDS9(regionOutputFileName,m_sources,convertDS9RegionsToWCS);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to write source island DS9 regions!");
		#endif
		return -1;
	}

	//Save DS9 regions for components
	status= SourceExporter::WriteComponentsToDS9(regionComponentsOutputFileName,m_sources,convertDS9RegionsToWCS);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to write source component DS9 regions!");
		#endif
		return -1;
	}

	//Save DS9 region for re-fitted sources only
	status= SourceExporter::WriteToDS9(regionOutputFileName_refitted,m_fittedSources,convertDS9RegionsToWCS);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to write refitted source island DS9 regions!");
		#endif
		return -1;
	}

	//Save DS9 region for re-fitted source components only
	status= SourceExporter::WriteComponentsToDS9(regionComponentsOutputFileName_refitted,m_fittedSources,convertDS9RegionsToWCS);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to write refitted source component DS9 regions!");
		#endif
		return -1;
	}

	return 0;

}//close SaveDS9Regions()

int SaveSources()
{
	//Check output file is open
	if(!outputFile || !outputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Output ROOT file not allocated or open!");
		#endif
		return -1;
	}
	if(!outputTree){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Output ROOT source TTree not allocated!");
		#endif
		return -1;
	}

	//Loop over selected sources and write TTree to output file
	for(size_t i=0;i<m_sources.size();i++){
		m_source= m_sources[i];
		outputTree->Fill();
	}

	outputTree->Write();

	return 0;

}//close SaveSources()


int SaveCatalog()
{
	//Return if no sources are found
	if(m_sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No sources selected, no catalog file will be written!");
		#endif
		return 0;
	}

	//Retrieve source WCS
	WCS* wcs= m_sources[0]->GetWCS(m_ds9WCSType);
	if(!wcs) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute WCS from sources!");
		#endif
	}	

	//Saving island/blob catalog to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogOutputFileName<<" ...");
	#endif
	bool dumpNestedSourceInfo= true;
	int status= SourceExporter::WriteToAscii(catalogOutputFileName,m_sources,dumpNestedSourceInfo,m_ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source catalog to file "<<catalogOutputFileName<<" failed!");
		#endif
	}
	
	//Saving source fitted components to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogComponentsOutputFileName<<" ...");
	#endif
	status= SourceExporter::WriteComponentsToAscii(catalogComponentsOutputFileName,m_sources,dumpNestedSourceInfo,m_ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source fitted component catalog to file "<<catalogComponentsOutputFileName<<" failed!");
		#endif
	}

	return 0;

}//close SaveCatalog()

void Save()
{
	//Save TTree to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		
		//Clone all objects present in input file but the Source TTree
		CloneObjectsInFile({"SourceInfo"});

		//Write selected source TTree
		SaveSources();

		//Close file
		outputFile->Close();

	}//close if

	//Save DS9 regions with selected sources	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving DS9 regions with selected sources...");
	#endif
	SaveDS9Regions();

	//Save ascii catalog with selected sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving ascii catalog with selected sources...");
	#endif
	SaveCatalog();

}//close Save()

int CloneObjectsInFile(std::vector<std::string> excludedObjNames)
{	
	//Check if input file is open
	if(!inputFile || !inputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Input file is not open!");
		#endif
		return -1;
	}

	//Check if output file is open
	if(!outputFile || !outputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Output file is not open!");
		#endif
		return -1;
	}

	//Loop on all entries present in input file
	TKey* key;
  TIter nextkey(inputFile->GetListOfKeys());
	
	while ((key = (TKey*)nextkey())) 
	{
		std::string className= key->GetClassName();
		std::string keyName= key->GetName();
    TClass* cl = gROOT->GetClass(className.c_str());
    if (!cl) continue;
	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Processing object key "<<keyName<<" (name="<<className<<")...");
		#endif

		bool excludeObj= false;
		for(size_t i=0;i<excludedObjNames.size();i++){
			if(keyName==excludedObjNames[i]){
				excludeObj= true;
				break;
			}
		}
		if(excludeObj) {
			INFO_LOG("Object "<<keyName<<" exluded from the list of objects that will be saved ...");
			continue;
		}

		if (cl->InheritsFrom(TTree::Class())) {
    	TTree* T= (TTree*)inputFile->Get(key->GetName());

			//Write to output file
			outputFile->cd();
      TTree* newT= T->CloneTree(-1,"fast");
      newT->Write();
			
		}//close if TTree object
		else{
			TObject* obj= key->ReadObj();
      outputFile->cd();
      obj->Write();
      delete obj;			
		}

	}//end loop objects

	return 0;

}//close CloneObjectsInFile()

std::string GetStringLogLevel(int verbosity){

	std::string slevel= "";
	if(verbosity<=0) slevel= "FATAL";
	else if(verbosity==1) slevel= "FATAL";
	else if(verbosity==2) slevel= "ERROR";
	else if(verbosity==3) slevel= "WARN";
	else if(verbosity==4) slevel= "INFO";
	else if(verbosity>=5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

