#include <Consts.h>
#include <Logger.h>
#include <Source.h>
#include <SourceFitter.h>

#include <TFile.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>

#include <iostream>

using namespace Caesar;
using namespace std;
	
double gScaleMin= 1;//x beam
double gScaleMax= 2;//x beam
double gScaleStep= 1;
double gBlobThrFactor= 0;
int gBlobKernelFactor= 6;
double gPeakZThr= 1;
int maxNPeaks= 10;
int gMinBlobSize= 5;
int gPeakShiftTolerance= 2;
double gParBoundIncreaseStepSize= 0.1;
int gMaxFitComponents= 10;
bool gUseNestedAsComponents= false;
bool gFitImproveConvergence= true;
int gFitNRetries= 100;
bool gDoFinalMinimizerStep= false;

bool gFixAmplInPreFit= true;
bool gFixSigmaInPreFit= false;
bool gFixThetaInPreFit= false;
bool gFixCentroidInPreFit= false;
double gAmplLimit= 0.5;
double gSigmaLimit= 0.5;
double gCentroidLimit= 5;
double gThetaLimit= 10;

bool gFitScaleDataToMax= false;
int gFitPrintLevel= 100;
double gFitFcnTolerance= 1.e-2;
long int gFitMaxIters= 1000000;
std::string gFitMinimizer= "Minuit2";//"Minuit", "GSL"
std::string gFitMinimizerAlgo= "minimizer";//Migrad, BFGS
int gImgPixMargin= 5;
bool gFitRetryWithLessComponents= false;

int TestFitSource(std::string filename,int sourceIndex,std::string minimizer="Minuit2",std::string minimizerAlgo="minimizer",bool improveFit=true,double amplLimit=0.5,double sigmaLimit=0.5,double centroidLimit=5,double thetaLimit=10,double fitTol=1.e-2,bool scaleDataToMax=false,bool retryWithLessComponents=false)
{
	//#####################################
	//##         SET ARGS
	//#####################################
	gAmplLimit= amplLimit;
	gSigmaLimit= sigmaLimit;
	gCentroidLimit= centroidLimit;
	gThetaLimit= thetaLimit;
	gFitMinimizer= minimizer;
	gFitMinimizerAlgo= minimizerAlgo;
	gFitImproveConvergence= improveFit;
	gFitFcnTolerance= fitTol;
	gFitScaleDataToMax= scaleDataToMax;
	gFitRetryWithLessComponents= retryWithLessComponents;

	//#####################################
	//##         READ FILE & IMAGE
	//#####################################
	//Read file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile){
		cerr<<"ERROR: Failed to open file "<<filename<<"!"<<endl;
		return -1;
	}

	//REad Tree with source list
	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");
	if(!SourceInfo){
		cerr<<"ERROR: Failed to read source data tree!"<<endl;
		return -1;
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);
	SourceInfo->GetEntry(sourceIndex);

	//Compute source image
	Image* sourceImg= aSource->GetImage(eFluxMap,gImgPixMargin);	
	double bkgAvg= aSource->GetBkgSum()/aSource->GetNPixels();
	double bkgRMSAvg= aSource->GetBkgRMSSum()/aSource->GetNPixels();


	//#####################################
	//##         SOURCE PARS
	//#####################################	
	//Get metadata pars
	ImgMetaData* metadata= aSource->GetImageMetaData();
	double beamBmaj= metadata->Bmaj;
 	double beamBmin= metadata->Bmin;
	double beamBpa= metadata->Bpa;
	double pixSizeX= metadata->dX; 
	double pixSizeY= metadata->dY;
	double pixSize= fabs(std::min(pixSizeX,pixSizeY));
	double beamWidth= fabs(std::min(beamBmaj,beamBmin));
	double beamPixSize= beamWidth/pixSize;
	double sigmaMin= gScaleMin*beamPixSize/GausSigma2FWHM;//convert from FWHM to sigma
	double sigmaMax= gScaleMax*beamPixSize/GausSigma2FWHM;//convert from FWHM to sigma
	double sigmaStep= gScaleStep;
	
	cout<<"INFO: Beam pars ("<<beamBmaj<<","<<beamBmin<<","<<beamBpa<<"), beamWidth="<<beamWidth<<", beamPixSize="<<beamPixSize<<", pixSizeXY("<<pixSizeX<<","<<pixSizeY<<"), pixSize="<<pixSize<<endl;
	cout<<"INFO: sigma min/max/step="<<sigmaMin<<"/"<<sigmaMax<<"/"<<sigmaStep<<endl;
	
	

	//#######################################
	//##     FIT SOURCE
	//#######################################
	//Set fit options
	SourceFitOptions fitOptions;	
	fitOptions.bmaj= fabs(beamBmaj/pixSizeX);//converted in pixels
	fitOptions.bmin= fabs(beamBmin/pixSizeY);//converted in pixels
	fitOptions.bpa= beamBpa + 90.;//beam pos angle is computed from North. We are assuming angles computed from x axis;
	fitOptions.nMaxComponents= gMaxFitComponents;
	fitOptions.limitCentroidInFit= true;
	fitOptions.fixCentroidInPreFit= gFixCentroidInPreFit;
	fitOptions.centroidLimit= gCentroidLimit;
	fitOptions.limitBkgInFit= true;
	fitOptions.fixBkg= true;
	fitOptions.useEstimatedBkgLevel= true;
	fitOptions.fixedBkgLevel= 0;
	fitOptions.limitAmplInFit= true;
	fitOptions.fixAmplInPreFit= gFixAmplInPreFit;
	fitOptions.amplLimit= gAmplLimit;
	fitOptions.fixSigmaInPreFit= gFixSigmaInPreFit;
	fitOptions.limitSigmaInFit= true;
	fitOptions.sigmaLimit= gSigmaLimit;
	fitOptions.fixSigma= false;
	fitOptions.limitThetaInFit= true;
	fitOptions.fixThetaInPreFit= gFixThetaInPreFit;
	fitOptions.fixTheta= false;
	fitOptions.thetaLimit= gThetaLimit;
	fitOptions.useFluxZCut= true;
	fitOptions.fluxZThrMin= 2.5;
	fitOptions.peakMinKernelSize= 3;
	fitOptions.peakMaxKernelSize= 7;
	fitOptions.peakZThrMin= gPeakZThr;
	fitOptions.peakKernelMultiplicityThr= 1;
	fitOptions.peakShiftTolerance= gPeakShiftTolerance;
	fitOptions.scaleMin= sigmaMin;
	fitOptions.scaleMax= sigmaMax;
	fitOptions.scaleStep= sigmaStep;
	fitOptions.minBlobSize= gMinBlobSize;
	fitOptions.blobMapThrFactor= gBlobThrFactor;
	fitOptions.blobMapKernelFactor= gBlobKernelFactor; 
	fitOptions.useNestedAsComponents= gUseNestedAsComponents;
	fitOptions.chi2RegPar= 1;

	//New test options
	fitOptions.fitFcnTolerance= gFitFcnTolerance;
	fitOptions.fitMaxIters= gFitMaxIters;
	fitOptions.fitPrintLevel= gFitPrintLevel;
	fitOptions.fitStrategy= 2;
	fitOptions.fitImproveConvergence= gFitImproveConvergence;
	fitOptions.fitNRetries= gFitNRetries;
	fitOptions.fitParBoundIncreaseStepSize= gParBoundIncreaseStepSize;
	fitOptions.fitMinimizer= gFitMinimizer;
	fitOptions.fitMinimizerAlgo= gFitMinimizerAlgo;
	//fitOptions.fitMinimizer= "GSL";
	//fitOptions.fitMinimizerAlgo= "BFGS";
	fitOptions.fitDoFinalMinimizerStep= gDoFinalMinimizerStep;
	fitOptions.fitScaleDataToMax= gFitScaleDataToMax;
	fitOptions.fitRetryWithLessComponents= gFitRetryWithLessComponents;

	if(aSource->Fit(fitOptions)<0) {
		WARN_LOG("Failed to fit source!");
	}

	
	//#######################################
	//##         DRAW PLOTS
	//#######################################
	cout<<"INFO: Drawing plots ..."<<endl;

	//Draw source map
	cout<<"INFO: Draw source map ..."<<endl;
	TCanvas* SourceMapPlot= new TCanvas("SourceMapPlot","SourceMapPlot");
	SourceMapPlot->cd();
	
	TH2D* sourceHisto= sourceImg->GetHisto2D("sourceHisto");
	sourceHisto->SetStats(0);
	sourceHisto->Draw("COLZ");

	
	//Draw source & fit	
	cout<<"INFO: Draw source & fit..."<<endl;
	int pixMargin= gImgPixMargin;
	bool drawImg= true;
	bool drawContours= true;
	bool drawNested= false;
	bool drawFitComponents= true;
	aSource->Draw(pixMargin,eFluxMap,drawImg,drawContours,drawNested,drawFitComponents,kRed,kSolid);
	
	cout<<"INFO: End macro"<<endl;

	return 0;

}//close macro



