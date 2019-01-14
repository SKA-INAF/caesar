#include <TROOT.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <Math/QuantFuncMathCore.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
using namespace TMVA;

#include <Image.h>
#include <BkgData.h>
#include <Logger.h>
#include <Source.h>
#include <MathUtils.h>
#include <AstroUtils.h>
#include <StatsUtils.h>


using namespace Caesar;

//## Define global vars
std::vector<std::string> FileNames;

//#############################################
struct MacroOptions {

	MacroOptions(){
		SetDefaults();
	};

	//## Set default options
	void SetDefaults()
	{	
		//- Output filename
		outputFileName= "DrawOutput.root";
		
		//- Classifier options
		applySourceClassifier= true;
		applyUserRecCuts= true;
		recCutsClassifierWeights= "dataset/weights/SourceCutsClassification_Cuts.weights.xml";
		signalCutEff= 0.9;
		nnClassifierWeights= "dataset/weights/SourceNNClassification_MLP.weights.xml";
		nnCut= 0.5;
		addSourceSNRToNNVars= false;

		//- Beam & map info
		//bkgLevel_true= 10.e-6;
		noiseLevel_true= 400.e-6;
		Bmaj= 13.347706794699599;//arcsec
		Bmin= 8.35635757446;//arcsec
		Bpa= 6.416473388670e-01;//deg
		pixSize= 1;//in arcsec

		//Source selection cuts
		selectResolvedTrueSources= true;
		mutualTrueSourceDistThr= 1;//arcsec
		selectFitComponents= true;
		excludeSourceAtEdge= true;
		sourceRegionXMin= 200;//in pixels
		sourceRegionXMax= 2360;//in pixels
		sourceRegionYMin= 200;//in pixels
		sourceRegionYMax= 2360;//in pixels

		applyFitEllipseCuts= false;
		deltaFitThetaMinThr= -45;//in deg
		deltaFitThetaMaxThr= 45;//in deg
		eccentricityRatioMinThr= 0.5;//0.7;
		eccentricityRatioMaxThr= 1.5;//1.3;
		areaToBeamRatioMinThr= 0.01;
		areaToBeamRatioMaxThr= 10;

		applyChi2Cut= true;
		chi2Cut= 10;
		//applyFitComponentFlagCut= true;
		applyFitQualityCut= true;
		fitQualityCut= 3;

		//- Drawing options
		drawErrorX= false;
	
	}//close SetDefaults()

	//- Output options 
	std::string outputFileName;

	//- Classifier options
	bool applySourceClassifier;
	std::string recCutsClassifierWeights;
	bool applyUserRecCuts;
	double signalCutEff;
	std::string nnClassifierWeights;
	double nnCut;
	bool addSourceSNRToNNVars;

	//- Beam & map info
	//double bkgLevel_true;//in Jy/beam
	double noiseLevel_true;//in Jy/beam
	double Bmaj;//in arcsec
	double Bmin;//in arcsec
	double Bpa;//in deg
	double pixSize;//in arcsec

	//- Source selection cuts
	bool selectResolvedTrueSources;
	double mutualTrueSourceDistThr;//in arcsec
	bool selectFitComponents;
	bool excludeSourceAtEdge;
	double sourceRegionXMin;//in pixels
	double sourceRegionXMax;//in pixels
	double sourceRegionYMin;//in pixels
	double sourceRegionYMax;//in pixels
	
	bool applyFitEllipseCuts;
	double deltaFitThetaMinThr;//in deg
	double deltaFitThetaMaxThr;//in deg
	double eccentricityRatioMinThr;
	double eccentricityRatioMaxThr;
	double areaToBeamRatioMinThr;
	double areaToBeamRatioMaxThr;

	bool applyChi2Cut;
	double chi2Cut;

	bool applyFitQualityCut;
	int fitQualityCut;

	//- Draw options
	bool drawErrorX;

};
MacroOptions opt;
//##############################################




//====================================
//        VARIABLES
//===================================
//## Beam & map info
double beamArea= 0;
double beamTheta= 0;
TEllipse* beamEllipse= 0;
double beamEllipseEccentricity= 0;
double beamEllipseArea= 0;


//Draw info
std::vector<double> Zbins= {0.1,0.5,1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000,2500,5000,10000};
std::vector<double> LgFluxBins= {
	-4.5,-4,-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.5,0,1,2
};
std::vector<double> LgFluxBins_ext= {
	-4.5,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,1,2
};

//- Output file
TFile* outputFile= 0;

//- Output data
TTree* FitEllipseParTree= 0;
TTree* PointSourcePerfInfo= 0;
TTree* RecPointSourcePerfInfo= 0;
std::string sourceFileName;
std::string sourceName;
float sourcePosX;
float sourcePosY;
float sourceTruePosX;
float sourceTruePosY;
float sourceTrueFlux;
float sourceFlux;
float ellMaj;
float ellMin;
float ellTheta;
float ellArea;
float ellEccentricity;
float dTheta;
float eccentricityRatio;
float areaToBeamRatio;
float sourceSNR;
float sourceTrueSNR;
int isTrueSource;
float fitChi2;
float fitNDF;
int fitQuality;
int fitComponentFlag;

//- Classifier data
TMVA::Reader* reader= 0;
TMVA::Reader* reader_NN= 0;

//- Compact source histos
TH1D* NTrueSourceHisto_compact= 0;
TH1D* NRecSourceHisto_compact_fit= 0;
TH1D* NRecSourceHisto_compact_fit_presel= 0;
TH1D* NRecSourceHisto_compact_fit_cutsel= 0;
TH1D* NRecSourceHisto_compact_fit_nnsel= 0;

TH1D* NTrueSourceHisto_reliability_compact_fit= 0;
TH1D* NRecSourceHisto_reliability_compact_fit= 0;
TH1D* NTrueSourceHisto_reliability_compact_fit_presel= 0;
TH1D* NRecSourceHisto_reliability_compact_fit_presel= 0;
TH1D* NTrueSourceHisto_reliability_compact_fit_cutsel= 0;
TH1D* NRecSourceHisto_reliability_compact_fit_cutsel= 0;
TH1D* NTrueSourceHisto_reliability_compact_fit_nnsel= 0;
TH1D* NRecSourceHisto_reliability_compact_fit_nnsel= 0;

//- Compact source completeness/reliability
TEfficiency::EStatOption gEfficiencyErrModel= TEfficiency::kFCP;
//TEfficiency::EStatOption gEfficiencyErrModel= TEfficiency::kFFC;

TEfficiency* Efficiency_compact_fit= 0;
TEfficiency* Efficiency_compact_fit_presel= 0;
TEfficiency* Efficiency_compact_fit_cutsel= 0;
TEfficiency* Efficiency_compact_fit_nnsel= 0;
TEfficiency* Reliability_compact_fit= 0;
TEfficiency* Reliability_compact_fit_presel= 0;
TEfficiency* Reliability_compact_fit_cutsel= 0;
TEfficiency* Reliability_compact_fit_nnsel= 0;

//- Compact source graphs
TGraphAsymmErrors* xPosAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* xPosResolutionGraph_compact_fit= 0;
TGraphAsymmErrors* yPosResolutionGraph_compact_fit= 0;
TGraphAsymmErrors* xPosAccuracyGraph_compact_fit_presel= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact_fit_presel= 0;
TGraphAsymmErrors* xPosResolutionGraph_compact_fit_presel= 0;
TGraphAsymmErrors* yPosResolutionGraph_compact_fit_presel= 0;
TGraphAsymmErrors* xPosAccuracyGraph_compact_fit_cutsel= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact_fit_cutsel= 0;
TGraphAsymmErrors* xPosResolutionGraph_compact_fit_cutsel= 0;
TGraphAsymmErrors* yPosResolutionGraph_compact_fit_cutsel= 0;
TGraphAsymmErrors* xPosAccuracyGraph_compact_fit_nnsel= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact_fit_nnsel= 0;
TGraphAsymmErrors* xPosResolutionGraph_compact_fit_nnsel= 0;
TGraphAsymmErrors* yPosResolutionGraph_compact_fit_nnsel= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* FluxDensityResolutionGraph_compact_fit= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact_fit_presel= 0;
TGraphAsymmErrors* FluxDensityResolutionGraph_compact_fit_presel= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact_fit_cutsel= 0;
TGraphAsymmErrors* FluxDensityResolutionGraph_compact_fit_cutsel= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact_fit_nnsel= 0;
TGraphAsymmErrors* FluxDensityResolutionGraph_compact_fit_nnsel= 0;
TGraph* FluxDensityAccuracyGraphPoints_compact= 0;


//- Compact source data vectors
std::vector< std::vector<double> > FluxList_compact_fit;
std::vector< std::vector<double> > FluxList_compact_fit_presel;
std::vector< std::vector<double> > FluxList_compact_fit_cutsel;
std::vector< std::vector<double> > FluxList_compact_fit_nnsel;
std::vector< std::vector<double> > xPosPullList_compact_fit;
std::vector< std::vector<double> > yPosPullList_compact_fit;
std::vector< std::vector<double> > xPosPullList_compact_fit_presel;
std::vector< std::vector<double> > yPosPullList_compact_fit_presel;
std::vector< std::vector<double> > xPosPullList_compact_fit_cutsel;
std::vector< std::vector<double> > yPosPullList_compact_fit_cutsel;
std::vector< std::vector<double> > xPosPullList_compact_fit_nnsel;
std::vector< std::vector<double> > yPosPullList_compact_fit_nnsel;

std::vector< std::vector<double> > FluxDensityPullList_compact_fit;
std::vector< std::vector<double> > FluxDensityPullList_compact_fit_presel;
std::vector< std::vector<double> > FluxDensityPullList_compact_fit_cutsel;
std::vector< std::vector<double> > FluxDensityPullList_compact_fit_nnsel;

std::vector<TH1D*> FluxDensityPullHistos_compact_fit_presel;

//- Extended source histos/graphs
TH1D* NTrueSourceHisto_ext= 0;
TH1D* NRecSourceHisto_ext= 0;
TH1D* NTrueSourceHisto_reliability_ext= 0;
TH1D* NRecSourceHisto_reliability_ext= 0;
TEfficiency* Efficiency_ext= 0;
TEfficiency* Reliability_ext= 0;

const int nSimTypes= 6;
std::vector<std::string> SimTypeLabels {"ring","bubble","ellipse","disk","blob","composite"};
std::vector<int> SimTypeColors {kRed,kGreen,kBlue,kMagenta,kOrange,kGray+1};
std::map<int,int> SimTypeToIndexMap {
	{1,0}, {15,0},
	{2,1}, {25,1},
	{3,2}, {35,2},
	{4,3}, {45,3},
	{5,4}, {55,4},
	{135,5}, {235,5}, {125,5}, {145,5}, {345,5}, {1235,5}
};
std::vector<TH1D*> NTrueSourceHisto_simtypes_ext;
std::vector<TH1D*> NRecSourceHisto_simtypes_ext;
std::vector<TEfficiency*> Efficiency_simtypes_ext;

TGraphAsymmErrors* xPosAccuracyGraph_ext= 0;
TGraphAsymmErrors* yPosAccuracyGraph_ext= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_ext= 0;

std::vector< std::vector<double> > FluxList_ext;
std::vector< std::vector<double> > xPosPullList_ext;
std::vector< std::vector<double> > yPosPullList_ext;
std::vector< std::vector<double> > FluxDensityPullList_ext;
std::vector< std::vector<double> > FluxSignificanceList_ext;

//====================================
//        FUNCTIONS
//===================================
//## Define functions
void Init();
void Draw();
void Save();
int AnalyzeData(std::string inputFileName);
int FillAnalysisHisto();
void ComputeAnalysisHistos();
double GausPosSigmaFcn(double* x,double* par);
double FitPosErrFcn(double* x,double* par);
double GetPosAngleInRange(double pa)
{
	while (pa <= -90) pa += 180;
  while (pa > 90) pa -= 180;
  return pa;
}

//====================================
//        MACRO
//===================================
int DrawSFinderPerformances(std::string _fileName,bool isFileList=false,MacroOptions macroOptions=MacroOptions()){

	//## Set logging level
	LoggerManager::Instance().CreateConsoleLogger("INFO","logger","System.out");

	//Set vars
	opt= macroOptions;
	
	//## Get list of files to be processed
	FileNames.clear();

	if(isFileList){//file list

		//Check filenames
		std::ifstream fileStream(_fileName);	
		std::string line;
		if (fileStream.fail() || !fileStream.is_open()){
			ERROR_LOG("Failed to open file "<<_fileName<<" for reading...");
			return -1;
		}
		
		//Store filenames present in lists
		std::string filename= "";
		while (std::getline(fileStream, line)) {
    	std::istringstream iss(line);
    	if (!(iss >> filename)) { 
				ERROR_LOG("Failed to read line from file "<<_fileName<<"!");
				return -1; 
			}
    	FileNames.push_back(filename);
		}//end file read
				
		
	}//close if
	else{//single file
		FileNames.push_back(_fileName);
	}//close if

	if(FileNames.empty()){
		ERROR_LOG("Empty filelist, nothing to be processed!");
		return -1;
	}


	//## Init data
	Init();

	//## Analyze data
	for(size_t i=0;i<FileNames.size();i++){
		if(AnalyzeData(FileNames[i])<0){
			ERROR_LOG("Failed to analyze data for file no. "<<i+1<<", skip to next...");
			continue;
		}
	}//end loop files

	//## Fill analysis histos
	FillAnalysisHisto();

	//## Compute analysis histo
	ComputeAnalysisHistos();

	//## Draw 
	Draw();

	//Save data
	Save();

	return 0;
	
}//close macro

void Save()
{
	//Save data to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();
	
		//Store tree
		if(FitEllipseParTree) FitEllipseParTree->Write();
		if(PointSourcePerfInfo) PointSourcePerfInfo->Write();
		if(RecPointSourcePerfInfo) RecPointSourcePerfInfo->Write();

		//Close file
		outputFile->Close();
	}

}//close Save()

void Init()
{
	//Open output ROOT file
	outputFile= new TFile(opt.outputFileName.c_str(),"RECREATE");
	outputFile->cd();

	//Create output trees
	FitEllipseParTree= new TTree("FitEllipseInfo","FitEllipseInfo");
	FitEllipseParTree->Branch("isTrueSource",&isTrueSource,"isTrueSource/I");
	FitEllipseParTree->Branch("beamBmaj",&opt.Bmaj,"beamBmaj/D");
	FitEllipseParTree->Branch("beamBmin",&opt.Bmin,"beamBmin/D");
	FitEllipseParTree->Branch("beamTheta",&beamTheta,"beamTheta/D");
	FitEllipseParTree->Branch("beamArea",&beamEllipseArea,"beamArea/D");
	FitEllipseParTree->Branch("beamEccentricity",&beamEllipseEccentricity,"beamEccentricity/D");
	FitEllipseParTree->Branch("fitEllipseBmaj",&ellMaj,"fitEllipseBmaj/F");
	FitEllipseParTree->Branch("fitEllipseBmin",&ellMin,"fitEllipseBmin/F");
	FitEllipseParTree->Branch("fitEllipseTheta",&ellTheta,"fitEllipseTheta/F");
	FitEllipseParTree->Branch("fitEllipseArea",&ellArea,"fitEllipseArea/F");
	FitEllipseParTree->Branch("fitEllipseEccentricity",&ellEccentricity,"fitEllipseEccentricity/F");	
	FitEllipseParTree->Branch("thetaDiff",&dTheta,"thetaDiff/F");
	FitEllipseParTree->Branch("eccentricityRatio",&eccentricityRatio,"eccentricityRatio/F");
	FitEllipseParTree->Branch("sourceToBeamRatio",&areaToBeamRatio,"sourceToBeamRatio/F");
	FitEllipseParTree->Branch("sourceSNR",&sourceSNR,"sourceSNR/F");

	PointSourcePerfInfo= new TTree("PointSourcePerfInfo","PointSourcePerfInfo");
	PointSourcePerfInfo->Branch("filename",&sourceFileName);	
	PointSourcePerfInfo->Branch("sourceName",&sourceName);	
	PointSourcePerfInfo->Branch("sourcePosX",&sourcePosX,"sourcePosX/F");
	PointSourcePerfInfo->Branch("sourcePosY",&sourcePosY,"sourcePosY/F");
	PointSourcePerfInfo->Branch("sourceTruePosX",&sourceTruePosX,"sourceTruePosX/F");
	PointSourcePerfInfo->Branch("sourceTruePosY",&sourceTruePosY,"sourceTruePosY/F");
	PointSourcePerfInfo->Branch("sourceFlux",&sourceFlux,"sourceFlux/F");
	PointSourcePerfInfo->Branch("sourceTrueFlux",&sourceTrueFlux,"sourceTrueFlux/F");
	PointSourcePerfInfo->Branch("thetaDiff",&dTheta,"thetaDiff/F");
	PointSourcePerfInfo->Branch("eccentricityRatio",&eccentricityRatio,"eccentricityRatio/F");
	PointSourcePerfInfo->Branch("sourceToBeamRatio",&areaToBeamRatio,"sourceToBeamRatio/F");
	PointSourcePerfInfo->Branch("sourceSNR",&sourceSNR,"sourceSNR/F");
	PointSourcePerfInfo->Branch("sourceTrueSNR",&sourceTrueSNR,"sourceTrueSNR/F");
	PointSourcePerfInfo->Branch("fitChi2",&fitChi2,"fitChi2/F");
	PointSourcePerfInfo->Branch("fitNDF",&fitNDF,"fitNDF/F");
	PointSourcePerfInfo->Branch("fitQuality",&fitQuality,"fitQuality/I");
	PointSourcePerfInfo->Branch("fitComponentFlag",&fitComponentFlag,"fitComponentFlag/I");

	RecPointSourcePerfInfo= new TTree("RecPointSourcePerfInfo","RecPointSourcePerfInfo");
	RecPointSourcePerfInfo->Branch("filename",&sourceFileName);
	RecPointSourcePerfInfo->Branch("sourceName",&sourceName);	
	RecPointSourcePerfInfo->Branch("sourcePosX",&sourcePosX,"sourcePosX/F");
	RecPointSourcePerfInfo->Branch("sourcePosY",&sourcePosY,"sourcePosY/F");
	RecPointSourcePerfInfo->Branch("sourceFlux",&sourceFlux,"sourceFlux/F");
	RecPointSourcePerfInfo->Branch("isTrueSource",&isTrueSource,"isTrueSource/I");
	RecPointSourcePerfInfo->Branch("thetaDiff",&dTheta,"thetaDiff/F");
	RecPointSourcePerfInfo->Branch("eccentricityRatio",&eccentricityRatio,"eccentricityRatio/F");
	RecPointSourcePerfInfo->Branch("sourceToBeamRatio",&areaToBeamRatio,"sourceToBeamRatio/F");
	RecPointSourcePerfInfo->Branch("sourceSNR",&sourceSNR,"sourceSNR/F");
	RecPointSourcePerfInfo->Branch("fitChi2",&fitChi2,"fitChi2/F");
	RecPointSourcePerfInfo->Branch("fitNDF",&fitNDF,"fitNDF/F");
	RecPointSourcePerfInfo->Branch("fitQuality",&fitQuality,"fitQuality/I");
	RecPointSourcePerfInfo->Branch("fitComponentFlag",&fitComponentFlag,"fitComponentFlag/I");

	// Create the MVA reader object	
	if(opt.applySourceClassifier){
		//Initialize the cut reader
		if(!opt.applyUserRecCuts){
			INFO_LOG("Initializing the MVA cut reader...");
  		reader = new TMVA::Reader( "!Color:!Silent" );
			reader->AddVariable("thetaDiff",&dTheta);
			reader->AddVariable("eccentricityRatio",&eccentricityRatio);
			reader->AddVariable("sourceToBeamRatio",&areaToBeamRatio);

			TString weightFile_cut= opt.recCutsClassifierWeights.c_str();
			INFO_LOG("Booking MVA method using weight file "<<weightFile_cut.Data()<<" ...");
			reader->BookMVA("Cuts",weightFile_cut);
		}


		//Initialize the NN reader
		INFO_LOG("Initializing the MVA NN reader...");
  	reader_NN = new TMVA::Reader( "!Color:!Silent" );
		reader_NN->AddVariable("thetaDiff",&dTheta);
		reader_NN->AddVariable("eccentricityRatio",&eccentricityRatio);
		reader_NN->AddVariable("sourceToBeamRatio",&areaToBeamRatio);
		if(opt.addSourceSNRToNNVars) reader_NN->AddVariable("sourceSNR",&sourceSNR);
	
		TString weightFile_NN= opt.nnClassifierWeights.c_str();
		INFO_LOG("Booking MVA NN method using weight file "<<weightFile_NN.Data()<<" ...");
		reader_NN->BookMVA("MLP",weightFile_NN);
	}

	//Compute map variables
	//beamTheta= MathUtils::Mod(opt.Bpa+90.,180.);
	beamTheta= opt.Bpa;//wrt to North
	beamArea= AstroUtils::GetBeamAreaInPixels(opt.Bmaj,opt.Bmin,opt.pixSize,opt.pixSize);
	beamEllipse= new TEllipse(0.,0.,opt.Bmaj/GausSigma2FWHM,opt.Bmin/GausSigma2FWHM,0,360,beamTheta);
	//beamEllipseEccentricity= MathUtils::ComputeEllipseEccentricity(beamEllipse);
	beamEllipseEccentricity= MathUtils::ComputeEllipseEccentricity(opt.Bmaj,opt.Bmin);
	//beamEllipseArea= MathUtils::ComputeEllipseArea(beamEllipse);
	beamEllipseArea= MathUtils::ComputeEllipseArea(opt.Bmaj,opt.Bmin);
	INFO_LOG("Beam ellipse: bmaj/bmin/bpa="<<opt.Bmaj<<"/"<<opt.Bmin<<"/"<<opt.Bpa<<", a="<<std::min(beamEllipse->GetR1(),beamEllipse->GetR2())<<", b="<<std::max(beamEllipse->GetR1(),beamEllipse->GetR2())<<", theta="<<beamEllipse->GetTheta()<<", A="<<beamEllipseArea<<", E="<<beamEllipseEccentricity);
	

	int nBins= (int)(LgFluxBins.size()-1);
	int nBins_ext= (int)(LgFluxBins_ext.size()-1);
	int nBins_Z= (int)(Zbins.size()-1);
	
	//## Init data vector
	for(int i=0;i<nBins;i++){
		FluxList_compact_fit.push_back( std::vector<double>() );
		FluxList_compact_fit_presel.push_back( std::vector<double>() );	
		FluxList_compact_fit_cutsel.push_back( std::vector<double>() );	
		FluxList_compact_fit_nnsel.push_back( std::vector<double>() );
		xPosPullList_compact_fit.push_back( std::vector<double>() );
		yPosPullList_compact_fit.push_back( std::vector<double>() );
		xPosPullList_compact_fit_presel.push_back( std::vector<double>() );
		yPosPullList_compact_fit_presel.push_back( std::vector<double>() );
		xPosPullList_compact_fit_cutsel.push_back( std::vector<double>() );
		yPosPullList_compact_fit_cutsel.push_back( std::vector<double>() );
		xPosPullList_compact_fit_nnsel.push_back( std::vector<double>() );
		yPosPullList_compact_fit_nnsel.push_back( std::vector<double>() );
		FluxDensityPullList_compact_fit.push_back( std::vector<double>() );
		FluxDensityPullList_compact_fit_presel.push_back( std::vector<double>() );
		FluxDensityPullList_compact_fit_cutsel.push_back( std::vector<double>() );
		FluxDensityPullList_compact_fit_nnsel.push_back( std::vector<double>() );
	}

	for(int i=0;i<nBins_ext;i++){
		FluxList_ext.push_back( std::vector<double>() );
		xPosPullList_ext.push_back( std::vector<double>() );
		yPosPullList_ext.push_back( std::vector<double>() );
		FluxDensityPullList_ext.push_back( std::vector<double>() );
		FluxSignificanceList_ext.push_back( std::vector<double>() );
	}

	//## Init histos	
	//- Completeness histos
	TString histoName= "NTrueSourceHisto_compact";
	NTrueSourceHisto_compact= new TH1D(histoName,"True compact source distribution",nBins,LgFluxBins.data());
	NTrueSourceHisto_compact->Sumw2();

	histoName= "NRecSourceHisto_compact_fit";
	NRecSourceHisto_compact_fit= new TH1D(histoName,"Fitted compact source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_compact_fit->Sumw2();
	
	histoName= "NRecSourceHisto_compact_fit_presel";
	NRecSourceHisto_compact_fit_presel= new TH1D(histoName,"Fitted compact source distribution (pre selection)",nBins,LgFluxBins.data());
	NRecSourceHisto_compact_fit_presel->Sumw2();

	histoName= "NRecSourceHisto_compact_fit_cutsel";
	NRecSourceHisto_compact_fit_cutsel= new TH1D(histoName,"Fitted compact source distribution (cut selection)",nBins,LgFluxBins.data());
	NRecSourceHisto_compact_fit_cutsel->Sumw2();

	histoName= "NRecSourceHisto_compact_fit_nnsel";
	NRecSourceHisto_compact_fit_nnsel= new TH1D(histoName,"Fitted compact source distribution (NN selection)",nBins,LgFluxBins.data());
	NRecSourceHisto_compact_fit_nnsel->Sumw2();

	histoName= "NTrueSourceHisto_ext";
	NTrueSourceHisto_ext= new TH1D(histoName,"True extended source distribution",nBins_ext,LgFluxBins_ext.data());
	NTrueSourceHisto_ext->Sumw2();

	histoName= "NRecSourceHisto_ext";
	NRecSourceHisto_ext= new TH1D(histoName,"Rec extended source distribution",nBins_ext,LgFluxBins_ext.data());
	NRecSourceHisto_ext->Sumw2();
		
	TH1D* histo= 0;
	for(int k=0;k<nSimTypes;k++){
		histoName= Form("NTrueSourceHisto_ext_simtype%d",k+1);
		histo= new TH1D(histoName,histoName,nBins_ext,LgFluxBins_ext.data());
		histo->Sumw2();
		NTrueSourceHisto_simtypes_ext.push_back(histo);
	
		histoName= Form("NRecSourceHisto_ext_simtype%d",k+1);
		histo= new TH1D(histoName,histoName,nBins_ext,LgFluxBins_ext.data());
		histo->Sumw2();
		NRecSourceHisto_simtypes_ext.push_back(histo);		
	}//end loop simtypes

	//- Reliability histos
	histoName= "NTrueSourceHisto_reliability_compact_fit";
	NTrueSourceHisto_reliability_compact_fit= new TH1D(histoName,"True compact fitted source distribution ",nBins,LgFluxBins.data());
	NTrueSourceHisto_reliability_compact_fit->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_compact_fit";
	NRecSourceHisto_reliability_compact_fit= new TH1D(histoName,"Rec compact fitted source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_reliability_compact_fit->Sumw2();

	histoName= "NTrueSourceHisto_reliability_compact_fit_presel";
	NTrueSourceHisto_reliability_compact_fit_presel= new TH1D(histoName,"True compact fitted source distribution (pre selection)",nBins,LgFluxBins.data());
	NTrueSourceHisto_reliability_compact_fit_presel->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_compact_fit_presel";
	NRecSourceHisto_reliability_compact_fit_presel= new TH1D(histoName,"Rec compact fitted source distribution (pre selection)",nBins,LgFluxBins.data());
	NRecSourceHisto_reliability_compact_fit_presel->Sumw2();

	histoName= "NTrueSourceHisto_reliability_compact_fit_cutsel";
	NTrueSourceHisto_reliability_compact_fit_cutsel= new TH1D(histoName,"True compact fitted source distribution (cut selection)",nBins,LgFluxBins.data());
	NTrueSourceHisto_reliability_compact_fit_cutsel->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_compact_fit_cutsel";
	NRecSourceHisto_reliability_compact_fit_cutsel= new TH1D(histoName,"Rec compact fitted source distribution (cut selection)",nBins,LgFluxBins.data());
	NRecSourceHisto_reliability_compact_fit_cutsel->Sumw2();

	histoName= "NTrueSourceHisto_reliability_compact_fit_nnsel";
	NTrueSourceHisto_reliability_compact_fit_nnsel= new TH1D(histoName,"True compact fitted source distribution (NN selection)",nBins,LgFluxBins.data());
	NTrueSourceHisto_reliability_compact_fit_nnsel->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_compact_fit_nnsel";
	NRecSourceHisto_reliability_compact_fit_nnsel= new TH1D(histoName,"Rec compact fitted source distribution (NN selection)",nBins,LgFluxBins.data());
	NRecSourceHisto_reliability_compact_fit_nnsel->Sumw2();


	histoName= "NTrueSourceHisto_reliability_ext";
	NTrueSourceHisto_reliability_ext= new TH1D(histoName,"True extended source distribution",nBins_ext,LgFluxBins_ext.data());
	NTrueSourceHisto_reliability_ext->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_ext";
	NRecSourceHisto_reliability_ext= new TH1D(histoName,"Rec extended source distribution",nBins_ext,LgFluxBins_ext.data());
	NRecSourceHisto_reliability_ext->Sumw2();

	//- Pos accuracy graphs
	TString graphName= "";

	graphName= "xPosAccuracyGraph_compact_fit";
	xPosAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact_fit->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact_fit";
	yPosAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact_fit->SetName(graphName);

	graphName= "xPosResolutionGraph_compact_fit";
	xPosResolutionGraph_compact_fit= new TGraphAsymmErrors;
	xPosResolutionGraph_compact_fit->SetName(graphName);
	
	graphName= "yPosResolutionGraph_compact_fit";
	yPosResolutionGraph_compact_fit= new TGraphAsymmErrors;
	yPosResolutionGraph_compact_fit->SetName(graphName);

	graphName= "xPosAccuracyGraph_compact_fit_presel";
	xPosAccuracyGraph_compact_fit_presel= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact_fit_presel->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact_fit_presel";
	yPosAccuracyGraph_compact_fit_presel= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact_fit_presel->SetName(graphName);

	graphName= "xPosResolutionGraph_compact_fit_presel";
	xPosResolutionGraph_compact_fit_presel= new TGraphAsymmErrors;
	xPosResolutionGraph_compact_fit_presel->SetName(graphName);
	
	graphName= "yPosResolutionGraph_compact_fit_presel";
	yPosResolutionGraph_compact_fit_presel= new TGraphAsymmErrors;
	yPosResolutionGraph_compact_fit_presel->SetName(graphName);


	graphName= "xPosAccuracyGraph_compact_fit_cutsel";
	xPosAccuracyGraph_compact_fit_cutsel= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact_fit_cutsel->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact_fit_cutsel";
	yPosAccuracyGraph_compact_fit_cutsel= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact_fit_cutsel->SetName(graphName);

	graphName= "xPosResolutionGraph_compact_fit_cutsel";
	xPosResolutionGraph_compact_fit_cutsel= new TGraphAsymmErrors;
	xPosResolutionGraph_compact_fit_cutsel->SetName(graphName);
	
	graphName= "yPosResolutionGraph_compact_fit_cutsel";
	yPosResolutionGraph_compact_fit_cutsel= new TGraphAsymmErrors;
	yPosResolutionGraph_compact_fit_cutsel->SetName(graphName);


	graphName= "xPosAccuracyGraph_compact_fit_nnsel";
	xPosAccuracyGraph_compact_fit_nnsel= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact_fit_nnsel->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact_fit_nnsel";
	yPosAccuracyGraph_compact_fit_nnsel= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact_fit_nnsel->SetName(graphName);

	graphName= "xPosResolutionGraph_compact_fit_nnsel";
	xPosResolutionGraph_compact_fit_nnsel= new TGraphAsymmErrors;
	xPosResolutionGraph_compact_fit_nnsel->SetName(graphName);
	
	graphName= "yPosResolutionGraph_compact_fit_nnsel";
	yPosResolutionGraph_compact_fit_nnsel= new TGraphAsymmErrors;
	yPosResolutionGraph_compact_fit_nnsel->SetName(graphName);

	graphName= "xPosAccuracyGraph_ext";
	xPosAccuracyGraph_ext= new TGraphAsymmErrors;
	xPosAccuracyGraph_ext->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_ext";
	yPosAccuracyGraph_ext= new TGraphAsymmErrors;
	yPosAccuracyGraph_ext->SetName(graphName);

	//- Flux accuracy graphs
	graphName= "FluxDensityAccuracyGraph_compact_fit";
	FluxDensityAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact_fit->SetName(graphName);

	graphName= "FluxDensityAccuracyGraph_compact_fit_presel";
	FluxDensityAccuracyGraph_compact_fit_presel= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact_fit_presel->SetName(graphName);

	graphName= "FluxDensityAccuracyGraph_compact_fit_cutsel";
	FluxDensityAccuracyGraph_compact_fit_cutsel= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact_fit_cutsel->SetName(graphName);

	graphName= "FluxDensityAccuracyGraph_compact_fit_nnsel";
	FluxDensityAccuracyGraph_compact_fit_nnsel= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact_fit_nnsel->SetName(graphName);


	graphName= "FluxDensityResolutionGraph_compact_fit";
	FluxDensityResolutionGraph_compact_fit= new TGraphAsymmErrors;
	FluxDensityResolutionGraph_compact_fit->SetName(graphName);

	graphName= "FluxDensityResolutionGraph_compact_fit_presel";
	FluxDensityResolutionGraph_compact_fit_presel= new TGraphAsymmErrors;
	FluxDensityResolutionGraph_compact_fit_presel->SetName(graphName);

	graphName= "FluxDensityResolutionGraph_compact_fit_cutsel";
	FluxDensityResolutionGraph_compact_fit_cutsel= new TGraphAsymmErrors;
	FluxDensityResolutionGraph_compact_fit_cutsel->SetName(graphName);

	graphName= "FluxDensityResolutionGraph_compact_fit_nnsel";
	FluxDensityResolutionGraph_compact_fit_nnsel= new TGraphAsymmErrors;
	FluxDensityResolutionGraph_compact_fit_nnsel->SetName(graphName);

	graphName= "FluxDensityAccuracyGraphPoints_compact";
	FluxDensityAccuracyGraphPoints_compact= new TGraph;
	FluxDensityAccuracyGraphPoints_compact->SetName(graphName);
	

	graphName= "FluxDensityAccuracyGraph_ext";
	FluxDensityAccuracyGraph_ext= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_ext->SetName(graphName);

	
	
}//close Init()



void ComputeAnalysisHistos()
{
	//##Compute detection efficiency
	//- Compact sources
	Efficiency_compact_fit = new TEfficiency(*NRecSourceHisto_compact_fit,*NTrueSourceHisto_compact); 
	Efficiency_compact_fit->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Efficiency_compact_fit->SetConfidenceLevel(0.68);

	Efficiency_compact_fit_presel = new TEfficiency(*NRecSourceHisto_compact_fit_presel,*NTrueSourceHisto_compact); 
	Efficiency_compact_fit_presel->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Efficiency_compact_fit_presel->SetConfidenceLevel(0.68);

	Efficiency_compact_fit_cutsel = new TEfficiency(*NRecSourceHisto_compact_fit_cutsel,*NTrueSourceHisto_compact); 
	Efficiency_compact_fit_cutsel->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Efficiency_compact_fit_cutsel->SetConfidenceLevel(0.68);

	Efficiency_compact_fit_nnsel = new TEfficiency(*NRecSourceHisto_compact_fit_nnsel,*NTrueSourceHisto_compact); 
	Efficiency_compact_fit_nnsel->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Efficiency_compact_fit_nnsel->SetConfidenceLevel(0.68);

	//- Extended sources
	Efficiency_ext = new TEfficiency(*NRecSourceHisto_ext,*NTrueSourceHisto_ext); 
	Efficiency_ext->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Efficiency_ext->SetConfidenceLevel(0.68);

	INFO_LOG("#"<<NTrueSourceHisto_ext->GetEntries()<<" extended sources analyzed");
	for(int i=0;i<NTrueSourceHisto_ext->GetNbinsX();i++){
		INFO_LOG("Bin "<<i+1<<": #"<<NTrueSourceHisto_ext->GetBinContent(i+1)<<" ext sources");
	}
	
	TEfficiency* effhisto= 0;
	for(int k=0;k<nSimTypes;k++){
		effhisto = new TEfficiency(*NRecSourceHisto_simtypes_ext[k],*NTrueSourceHisto_simtypes_ext[k]); 
		effhisto->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
		effhisto->SetConfidenceLevel(0.68);
		Efficiency_simtypes_ext.push_back(effhisto);
	}

	for(int k=0;k<nSimTypes;k++){
		INFO_LOG("Source type="<<SimTypeLabels[k]<<": found/tot="<<NRecSourceHisto_simtypes_ext[k]->GetEntries()<<"/"<<NTrueSourceHisto_simtypes_ext[k]->GetEntries()<<", eff="<<(double)NRecSourceHisto_simtypes_ext[k]->GetEntries()/(double)NTrueSourceHisto_simtypes_ext[k]->GetEntries());
		for(int i=0;i<NTrueSourceHisto_simtypes_ext[k]->GetNbinsX();i++){
			INFO_LOG("Bin "<<i+1<<": #ext sources="<<NTrueSourceHisto_simtypes_ext[k]->GetBinContent(i+1)<<", eff="<<(double)NRecSourceHisto_simtypes_ext[k]->GetBinContent(i+1)/(double)NTrueSourceHisto_simtypes_ext[k]->GetBinContent(i+1));
		}
	}
	

	//## Compute reliability 
	//- Compact sources
	Reliability_compact_fit= new TEfficiency(*NTrueSourceHisto_reliability_compact_fit,*NRecSourceHisto_reliability_compact_fit); 
	Reliability_compact_fit->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Reliability_compact_fit->SetConfidenceLevel(0.68);

	Reliability_compact_fit_presel= new TEfficiency(*NTrueSourceHisto_reliability_compact_fit_presel,*NRecSourceHisto_reliability_compact_fit_presel); 
	Reliability_compact_fit_presel->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Reliability_compact_fit_presel->SetConfidenceLevel(0.68);

	Reliability_compact_fit_cutsel= new TEfficiency(*NTrueSourceHisto_reliability_compact_fit_cutsel,*NRecSourceHisto_reliability_compact_fit_cutsel); 
	Reliability_compact_fit_cutsel->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Reliability_compact_fit_cutsel->SetConfidenceLevel(0.68);
	
	Reliability_compact_fit_nnsel= new TEfficiency(*NTrueSourceHisto_reliability_compact_fit_nnsel,*NRecSourceHisto_reliability_compact_fit_nnsel); 
	Reliability_compact_fit_nnsel->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Reliability_compact_fit_nnsel->SetConfidenceLevel(0.68);

	//- Extended sources
	Reliability_ext= new TEfficiency(*NTrueSourceHisto_reliability_ext,*NRecSourceHisto_reliability_ext); 
	Reliability_ext->SetStatisticOption(gEfficiencyErrModel);  // to set option for errors (see ref doc)
	Reliability_ext->SetConfidenceLevel(0.68);


	//Compute data list stats
	int nMinPointsToDraw= 2;
	int nBins= (int)(LgFluxBins.size()-1);
	int nBins_ext= (int)(LgFluxBins_ext.size()-1);
	int nPoints_xPull_fit= 0;
	int nPoints_yPull_fit= 0;	
	int nPoints_xPull_fit_presel= 0;
	int nPoints_yPull_fit_presel= 0;
	int nPoints_xPull_fit_cutsel= 0;
	int nPoints_yPull_fit_cutsel= 0;
	int nPoints_xPull_fit_nnsel= 0;
	int nPoints_yPull_fit_nnsel= 0;
	
	int nPoints_FluxPull_all= 0;
	int nPoints_FluxPull_fit= 0;	
	int nPoints_FluxPull_fit_presel= 0;
	int nPoints_FluxPull_fit_cutsel= 0;
	int nPoints_FluxPull_fit_nnsel= 0;
	int nPoints_OffsetPull= 0;
	int nPoints_xPull_ext= 0;
	int nPoints_yPull_ext= 0;	
	int nPoints_FluxPull_ext= 0;
	int nPoints_FluxPullVSSignificance_ext= 0;

	double pointOffsetX= 0;
	//double pointOffsetX= 0.04;

	for(int i=0;i<nBins;i++){
		double x_avg= LgFluxBins[i] + 0.5*(LgFluxBins[i+1]-LgFluxBins[i]);

		//Compute xpos accuracy	fit	
		if(xPosPullList_compact_fit[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(xPosPullList_compact_fit[i].size());
			Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(xPosPullList_compact_fit[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			x-= pointOffsetX;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
           
			double posx_mean= 0;
			double posx_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,xPosPullList_compact_fit[i]);
			double posx_median= stats_posx.median;
			double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
			double posx_iqr= stats_posx.Q3-stats_posx.Q1;
			double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
			double posx_iqr_half= posx_iqr/2.;
			double posx_iqr_half_err= 0.5*posx_iqr_err;

			//Bias graph
			double y= posx_median; 
			//double yerr_low= posx_median-stats_posx.Q1;
			//double yerr_up= stats_posx.Q3-posx_median;
			double yerr_low= posx_median_err;
			double yerr_up= posx_median_err;
			xPosAccuracyGraph_compact_fit->SetPoint(nPoints_xPull_fit,x,y);
			xPosAccuracyGraph_compact_fit->SetPointError(nPoints_xPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posx_iqr_half;
			yerr_low= posx_iqr_half_err;
			yerr_up= posx_iqr_half_err;
			xPosResolutionGraph_compact_fit->SetPoint(nPoints_xPull_fit,x,y);
			xPosResolutionGraph_compact_fit->SetPointError(nPoints_xPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);

			nPoints_xPull_fit++;
		}//close if

		//Compute ypos accuracy fit
		if(yPosPullList_compact_fit[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(yPosPullList_compact_fit[i].size());
			Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(yPosPullList_compact_fit[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			x+= pointOffsetX;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}

			double posy_mean= 0;
			double posy_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,yPosPullList_compact_fit[i]);
			double posy_median= stats_posy.median;
			double posy_err= posy_rms/sqrt(nEntries);
			double posy_median_err= 1.253*posy_rms/sqrt(nEntries);
			double posy_iqr= stats_posy.Q3-stats_posy.Q1;
			double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
			double posy_iqr_half= posy_iqr/2.;
			double posy_iqr_half_err= 0.5*posy_iqr_err;

			//Bias graph
			double y= posy_median;
			//double yerr_low= y-stats_posy.Q1;
			//double yerr_up= stats_posy.Q3-y;
			double yerr_low= posy_median_err;
			double yerr_up= posy_median_err;
			yPosAccuracyGraph_compact_fit->SetPoint(nPoints_yPull_fit,x,y);
			yPosAccuracyGraph_compact_fit->SetPointError(nPoints_yPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posy_iqr_half;
			yerr_low= posy_iqr_half_err;
			yerr_up= posy_iqr_half_err;
			yPosResolutionGraph_compact_fit->SetPoint(nPoints_yPull_fit,x,y);
			yPosResolutionGraph_compact_fit->SetPointError(nPoints_yPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);

			nPoints_yPull_fit++;
		}	


		//Compute xpos accuracy	fit	(PRE-SELECTION)
		if(xPosPullList_compact_fit_presel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(xPosPullList_compact_fit_presel[i].size());
			Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(xPosPullList_compact_fit_presel[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_presel[i]);

			std::map<std::string,double> statsBootstrapErr_posx;			
			int nBootstrapSamples= 1000;
			Caesar::StatsUtils::ComputeStatsBootstrapError(statsBootstrapErr_posx,xPosPullList_compact_fit_presel[i],nBootstrapSamples);

			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}

			double posx_mean= 0;
			double posx_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,xPosPullList_compact_fit_presel[i]);
			double posx_median= stats_posx.median;
			double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
			double posx_median_robust_err= statsBootstrapErr_posx["median"];
			double posx_iqr= stats_posx.Q3-stats_posx.Q1;
			double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
			double posx_iqr_robust_err= statsBootstrapErr_posx["iqr"];
			double posx_iqr_half= posx_iqr/2.;
			double posx_iqr_half_err= 0.5*posx_iqr_err;
			double posx_iqr_half_robust_err= 0.5*posx_iqr_robust_err;

			//Bias graph
			double y= posx_median;
			//double yerr_low= y-stats_posx.Q1;
			//double yerr_up= stats_posx.Q3-y;
			//double yerr_low= posx_median_err;
			//double yerr_up= posx_median_err;
			double yerr_low= posx_median_robust_err;
			double yerr_up= posx_median_robust_err;
			xPosAccuracyGraph_compact_fit_presel->SetPoint(nPoints_xPull_fit_presel,x,y);
			xPosAccuracyGraph_compact_fit_presel->SetPointError(nPoints_xPull_fit_presel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posx_iqr_half;
			//yerr_low= posx_iqr_half_err;
			//yerr_up= posx_iqr_half_err;
			yerr_low= posx_iqr_half_robust_err;
			yerr_up= posx_iqr_half_robust_err;
			xPosResolutionGraph_compact_fit_presel->SetPoint(nPoints_xPull_fit_presel,x,y);
			xPosResolutionGraph_compact_fit_presel->SetPointError(nPoints_xPull_fit_presel,xerr_low,xerr_up,yerr_low,yerr_up);

			nPoints_xPull_fit_presel++;
		}

		//Compute ypos accuracy fit (PRE-SELECTION)
		if(yPosPullList_compact_fit_presel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(yPosPullList_compact_fit_presel[i].size());
			Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(yPosPullList_compact_fit_presel[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_presel[i]);

			std::map<std::string,double> statsBootstrapErr_posy;			
			int nBootstrapSamples= 1000;
			Caesar::StatsUtils::ComputeStatsBootstrapError(statsBootstrapErr_posy,yPosPullList_compact_fit_presel[i],nBootstrapSamples);

			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
		
			double posy_mean= 0;
			double posy_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,yPosPullList_compact_fit_presel[i]);
			double posy_median= stats_posy.median;
			double posy_err= posy_rms/sqrt(nEntries);
			double posy_median_err= 1.253*posy_rms/sqrt(nEntries);	
			double posy_median_robust_err= statsBootstrapErr_posy["median"];
			double posy_iqr= stats_posy.Q3-stats_posy.Q1;
			double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
			double posy_iqr_robust_err= statsBootstrapErr_posy["iqr"];
			double posy_iqr_half= posy_iqr/2.;
			double posy_iqr_half_err= 0.5*posy_iqr_err;
			double posy_iqr_half_robust_err= 0.5*posy_iqr_robust_err;

			//Bias graph
			double y= posy_median;
			//double yerr_low= y-stats_posy.Q1;
			//double yerr_up= stats_posy.Q3-y;
			//double yerr_low= posy_median_err;
			//double yerr_up= posy_median_err;
			double yerr_low= posy_median_robust_err;
			double yerr_up= posy_median_robust_err;
			yPosAccuracyGraph_compact_fit_presel->SetPoint(nPoints_yPull_fit_presel,x,y);
			yPosAccuracyGraph_compact_fit_presel->SetPointError(nPoints_yPull_fit_presel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posy_iqr_half;
			//yerr_low= posy_iqr_half_err;
			//yerr_up= posy_iqr_half_err;
			yerr_low= posy_iqr_half_robust_err;
			yerr_up= posy_iqr_half_robust_err;
			yPosResolutionGraph_compact_fit_presel->SetPoint(nPoints_yPull_fit_presel,x,y);
			yPosResolutionGraph_compact_fit_presel->SetPointError(nPoints_yPull_fit_presel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_yPull_fit_presel++;
		}	

		//Compute xpos accuracy	fit	(CUT SELECTION)
		if(xPosPullList_compact_fit_cutsel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(xPosPullList_compact_fit_cutsel[i].size());
			Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(xPosPullList_compact_fit_cutsel[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_cutsel[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}

			double posx_mean= 0;
			double posx_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,xPosPullList_compact_fit_cutsel[i]);
			double posx_median= stats_posx.median;
			double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
			double posx_iqr= stats_posx.Q3-stats_posx.Q1;
			double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
			double posx_iqr_half= posx_iqr/2.;
			double posx_iqr_half_err= 0.5*posx_iqr_err;

			//Bias graph
			double y= posx_median;
			//double yerr_low= y-stats_posx.Q1;
			//double yerr_up= stats_posx.Q3-y;
			double yerr_low= posx_median_err;
			double yerr_up= posx_median_err;
			xPosAccuracyGraph_compact_fit_cutsel->SetPoint(nPoints_xPull_fit_cutsel,x,y);
			xPosAccuracyGraph_compact_fit_cutsel->SetPointError(nPoints_xPull_fit_cutsel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posx_iqr_half;
			yerr_low= posx_iqr_half_err;
			yerr_up= posx_iqr_half_err;
			yPosResolutionGraph_compact_fit_cutsel->SetPoint(nPoints_yPull_fit_cutsel,x,y);
			yPosResolutionGraph_compact_fit_cutsel->SetPointError(nPoints_yPull_fit_cutsel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_xPull_fit_cutsel++;
		}

		//Compute ypos accuracy fit (CUT SELECTION)
		if(yPosPullList_compact_fit_cutsel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(yPosPullList_compact_fit_cutsel[i].size());
			Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(yPosPullList_compact_fit_cutsel[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_cutsel[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
		
			double posy_mean= 0;
			double posy_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,yPosPullList_compact_fit_cutsel[i]);
			double posy_median= stats_posy.median;
			double posy_err= posy_rms/sqrt(nEntries);
			double posy_median_err= 1.253*posy_rms/sqrt(nEntries);
			double posy_iqr= stats_posy.Q3-stats_posy.Q1;
			double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
			double posy_iqr_half= posy_iqr/2.;
			double posy_iqr_half_err= 0.5*posy_iqr_err;

			//Bias graph
			double y= posy_median;
			//double yerr_low= y-stats_posy.Q1;
			//double yerr_up= stats_posy.Q3-y;
			double yerr_low= posy_median_err;
			double yerr_up= posy_median_err;
			yPosAccuracyGraph_compact_fit_cutsel->SetPoint(nPoints_yPull_fit_cutsel,x,y);
			yPosAccuracyGraph_compact_fit_cutsel->SetPointError(nPoints_yPull_fit_cutsel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posy_iqr_half;
			yerr_low= posy_iqr_half_err;
			yerr_up= posy_iqr_half_err;
			yPosResolutionGraph_compact_fit_cutsel->SetPoint(nPoints_yPull_fit_cutsel,x,y);
			yPosResolutionGraph_compact_fit_cutsel->SetPointError(nPoints_yPull_fit_cutsel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_yPull_fit_cutsel++;
		}	

		//Compute xpos accuracy	fit	(NN SELECTION)
		if(xPosPullList_compact_fit_nnsel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(xPosPullList_compact_fit_nnsel[i].size());
			Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(xPosPullList_compact_fit_nnsel[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_nnsel[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}

			double posx_mean= 0;
			double posx_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,xPosPullList_compact_fit_nnsel[i]);
			double posx_median= stats_posx.median;
			double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
			double posx_iqr= stats_posx.Q3-stats_posx.Q1;
			double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
			double posx_iqr_half= posx_iqr/2.;
			double posx_iqr_half_err= 0.5*posx_iqr_err;

			//Bias graph
			double y= posx_median;
			//double yerr_low= y-stats_posx.Q1;
			//double yerr_up= stats_posx.Q3-y;
			double yerr_low= posx_median_err;
			double yerr_up= posx_median_err;
			xPosAccuracyGraph_compact_fit_nnsel->SetPoint(nPoints_xPull_fit_nnsel,x,y);
			xPosAccuracyGraph_compact_fit_nnsel->SetPointError(nPoints_xPull_fit_nnsel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posx_iqr_half;
			yerr_low= posx_iqr_half_err;
			yerr_up= posx_iqr_half_err;
			yPosResolutionGraph_compact_fit_nnsel->SetPoint(nPoints_yPull_fit_nnsel,x,y);
			yPosResolutionGraph_compact_fit_nnsel->SetPointError(nPoints_yPull_fit_nnsel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			
			nPoints_xPull_fit_nnsel++;
		}

		//Compute ypos accuracy fit (NN SELECTION)
		if(yPosPullList_compact_fit_nnsel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(yPosPullList_compact_fit_nnsel[i].size());
			Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(yPosPullList_compact_fit_nnsel[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_nnsel[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
			
			double posy_mean= 0;
			double posy_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,yPosPullList_compact_fit_nnsel[i]);
			double posy_median= stats_posy.median;
			double posy_err= posy_rms/sqrt(nEntries);
			double posy_median_err= 1.253*posy_rms/sqrt(nEntries);
			double posy_iqr= stats_posy.Q3-stats_posy.Q1;
			double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
			double posy_iqr_half= posy_iqr/2.;
			double posy_iqr_half_err= 0.5*posy_iqr_err;

			//Bias graph
			double y= posy_median;
			//double yerr_low= y-stats_posy.Q1;
			//double yerr_up= stats_posy.Q3-y;
			double yerr_low= posy_median_err;
			double yerr_up= posy_median_err;
			yPosAccuracyGraph_compact_fit_nnsel->SetPoint(nPoints_yPull_fit_nnsel,x,y);
			yPosAccuracyGraph_compact_fit_nnsel->SetPointError(nPoints_yPull_fit_nnsel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= posy_iqr_half;
			yerr_low= posy_iqr_half_err;
			yerr_up= posy_iqr_half_err;
			yPosResolutionGraph_compact_fit_nnsel->SetPoint(nPoints_yPull_fit_nnsel,x,y);
			yPosResolutionGraph_compact_fit_nnsel->SetPointError(nPoints_yPull_fit_nnsel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_yPull_fit_nnsel++;
		}	


		//Compute Flux accuracy
		for(size_t j=0;j<FluxList_compact_fit[i].size();j++){
			FluxDensityAccuracyGraphPoints_compact->SetPoint(nPoints_FluxPull_all,FluxList_compact_fit[i][j],FluxDensityPullList_compact_fit[i][j]);
			nPoints_FluxPull_all++;
		}
		
		//Compute Flux accuracy for fitted sources
		if(FluxDensityPullList_compact_fit[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(FluxDensityPullList_compact_fit[i].size());
			Caesar::BoxStats<double> stats_fluxpull= StatsUtils::ComputeBoxStats(FluxDensityPullList_compact_fit[i]);	
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
		
			double fluxpull_mean= 0;
			double fluxpull_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(fluxpull_mean,fluxpull_rms,FluxDensityPullList_compact_fit[i]);
			double fluxpull_median= stats_fluxpull.median;
			double fluxpull_median_err= 1.253*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr= stats_fluxpull.Q3-stats_fluxpull.Q1;
			double fluxpull_iqr_err= 1.573*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr_half= fluxpull_iqr/2.;
			double fluxpull_iqr_half_err= 0.5*fluxpull_iqr_err;

			//Bias graph
			double y= fluxpull_median;
			//double yerr_low= y-stats_fluxpull.Q1;
			//double yerr_up= stats_fluxpull.Q3-y;
			double yerr_low= fluxpull_median_err;
			double yerr_up= fluxpull_median_err;
			FluxDensityAccuracyGraph_compact_fit->SetPoint(nPoints_FluxPull_fit,x,y);
			FluxDensityAccuracyGraph_compact_fit->SetPointError(nPoints_FluxPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= fluxpull_iqr_half;
			yerr_low= fluxpull_iqr_half_err;
			yerr_up= fluxpull_iqr_half_err;
			FluxDensityResolutionGraph_compact_fit->SetPoint(nPoints_FluxPull_fit,x,y);
			FluxDensityResolutionGraph_compact_fit->SetPointError(nPoints_FluxPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_FluxPull_fit++;
		}	


		//Compute Flux accuracy for fitted sources (PRE-SELECTION)
		if(FluxDensityPullList_compact_fit_presel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(FluxDensityPullList_compact_fit_presel[i].size());
			Caesar::BoxStats<double> stats_fluxpull= StatsUtils::ComputeBoxStats(FluxDensityPullList_compact_fit_presel[i]);	
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_presel[i]);
			
			std::map<std::string,double> statsBootstrapErr_fluxpull;			
			int nBootstrapSamples= 1000;
			Caesar::StatsUtils::ComputeStatsBootstrapError(statsBootstrapErr_fluxpull,FluxDensityPullList_compact_fit_presel[i],nBootstrapSamples);

			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;	
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
		
			double fluxpull_mean= 0;
			double fluxpull_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(fluxpull_mean,fluxpull_rms,FluxDensityPullList_compact_fit_presel[i]);
			double fluxpull_median= stats_fluxpull.median;
			double fluxpull_median_err= 1.253*fluxpull_rms/sqrt(nEntries);
			double fluxpull_median_robust_err= statsBootstrapErr_fluxpull["median"];
			double fluxpull_iqr= stats_fluxpull.Q3-stats_fluxpull.Q1;
			double fluxpull_iqr_err= 1.573*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr_robust_err= statsBootstrapErr_fluxpull["iqr"];
			double fluxpull_iqr_half= fluxpull_iqr/2.;
			double fluxpull_iqr_half_err= 0.5*fluxpull_iqr_err;
			double fluxpull_iqr_half_robust_err= 0.5*fluxpull_iqr_robust_err;
			double fluxpull_min= stats_fluxpull.minVal;
			double fluxpull_max= stats_fluxpull.maxVal;

			INFO_LOG("Bin "<<i+1<<": N="<<nEntries<<", <x>="<<x<<", fluxpull_mean="<<fluxpull_mean<<", fluxpull_rms="<<fluxpull_rms<<", fluxpull_median="<<fluxpull_median<<", fluxpull_median_err="<<fluxpull_median_err<<", fluxpull_median_robust_err="<<fluxpull_median_robust_err<<", iqr_half="<<fluxpull_iqr_half<<", fluxpull_iqr_half_err="<<fluxpull_iqr_half_err<<", fluxpull_iqr_half_robust_err="<<fluxpull_iqr_half_robust_err<<", fluxpull_min="<<fluxpull_min<<", fluxpull_max="<<fluxpull_max);

			//Bias graph
			double y= fluxpull_median;
			//double yerr_low= y-stats_fluxpull.Q1;
			//double yerr_up= stats_fluxpull.Q3-y;
			//double yerr_low= fluxpull_median_err;
			//double yerr_up= fluxpull_median_err;
			double yerr_low= fluxpull_median_robust_err;
			double yerr_up= fluxpull_median_robust_err;
			FluxDensityAccuracyGraph_compact_fit_presel->SetPoint(nPoints_FluxPull_fit_presel,x,y);
			FluxDensityAccuracyGraph_compact_fit_presel->SetPointError(nPoints_FluxPull_fit_presel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= fluxpull_iqr_half;
			//yerr_low= fluxpull_iqr_half_err;
			//yerr_up= fluxpull_iqr_half_err;
			yerr_low= fluxpull_iqr_half_robust_err;
			yerr_up= fluxpull_iqr_half_robust_err;
			FluxDensityResolutionGraph_compact_fit_presel->SetPoint(nPoints_FluxPull_fit_presel,x,y);
			FluxDensityResolutionGraph_compact_fit_presel->SetPointError(nPoints_FluxPull_fit_presel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			//Histo
			TString histoName= Form("FluxDensityPullHisto%d__compact_fit_presel",i+1);
			TH1D* h= new TH1D(histoName,histoName,500,fluxpull_min-10,fluxpull_max+10);
			h->Sumw2();
			for(size_t k=0;k<FluxDensityPullList_compact_fit_presel[i].size();k++){
				h->Fill(FluxDensityPullList_compact_fit_presel[i][k]);
			}
			FluxDensityPullHistos_compact_fit_presel.push_back(h);

			/*
			std::vector<double> fluxpull_data= FluxDensityPullList_compact_fit_presel[i];
			Caesar::StatsUtils::RemoveOutliers(fluxpull_data);
			nEntries= (int)(fluxpull_data.size());
			Caesar::BoxStats<double> stats_fluxpull_nooutliers= StatsUtils::ComputeBoxStats(fluxpull_data);
			fluxpull_mean= 0;
			fluxpull_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(fluxpull_mean,fluxpull_rms,fluxpull_data);
			fluxpull_median= stats_fluxpull_nooutliers.median;
			fluxpull_median_err= 1.253*fluxpull_rms/sqrt(nEntries);
			fluxpull_iqr= stats_fluxpull_nooutliers.Q3-stats_fluxpull_nooutliers.Q1;
			fluxpull_iqr_err= 1.573*fluxpull_rms/sqrt(nEntries);
			fluxpull_iqr_half= fluxpull_iqr/2.;
			fluxpull_iqr_half_err= 0.5*fluxpull_iqr_err;
			fluxpull_min= stats_fluxpull_nooutliers.minVal;
			fluxpull_max= stats_fluxpull_nooutliers.maxVal;

			INFO_LOG("Bin "<<i+1<<": N="<<nEntries<<", <x>="<<x<<", fluxpull_mean="<<fluxpull_mean<<", fluxpull_rms="<<fluxpull_rms<<", fluxpull_median="<<fluxpull_median<<", fluxpull_median_err="<<fluxpull_median_err<<", iqr_half="<<fluxpull_iqr_half<<", fluxpull_iqr_half_err="<<fluxpull_iqr_half_err<<", fluxpull_min="<<fluxpull_min<<", fluxpull_max="<<fluxpull_max);	
			*/
			

			nPoints_FluxPull_fit_presel++;
		}	
		
		//Compute Flux accuracy for fitted sources (CUT SELECTION)
		if(FluxDensityPullList_compact_fit_cutsel[i].size()>=nMinPointsToDraw){
			int nEntries= (int)(FluxDensityPullList_compact_fit_cutsel[i].size());
			Caesar::BoxStats<double> stats_fluxpull= StatsUtils::ComputeBoxStats(FluxDensityPullList_compact_fit_cutsel[i]);	
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_cutsel[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;	
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
		
			double fluxpull_mean= 0;
			double fluxpull_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(fluxpull_mean,fluxpull_rms,FluxDensityPullList_compact_fit_cutsel[i]);
			double fluxpull_median= stats_fluxpull.median;
			double fluxpull_median_err= 1.253*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr= stats_fluxpull.Q3-stats_fluxpull.Q1;
			double fluxpull_iqr_err= 1.573*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr_half= fluxpull_iqr/2.;
			double fluxpull_iqr_half_err= 0.5*fluxpull_iqr_err;

			//Bias graph
			double y= fluxpull_median;
			//double yerr_low= y-stats_fluxpull.Q1;
			//double yerr_up= stats_fluxpull.Q3-y;
			double yerr_low= fluxpull_median_err;
			double yerr_up= fluxpull_median_err;
			FluxDensityAccuracyGraph_compact_fit_cutsel->SetPoint(nPoints_FluxPull_fit_cutsel,x,y);
			FluxDensityAccuracyGraph_compact_fit_cutsel->SetPointError(nPoints_FluxPull_fit_cutsel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= fluxpull_iqr_half;
			yerr_low= fluxpull_iqr_half_err;
			yerr_up= fluxpull_iqr_half_err;
			FluxDensityResolutionGraph_compact_fit_cutsel->SetPoint(nPoints_FluxPull_fit_cutsel,x,y);
			FluxDensityResolutionGraph_compact_fit_cutsel->SetPointError(nPoints_FluxPull_fit_cutsel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_FluxPull_fit_cutsel++;
		}	

		//Compute Flux accuracy for fitted sources (NN SELECTION)
		if(FluxDensityPullList_compact_fit_nnsel[i].size()>=nMinPointsToDraw){	
			int nEntries= (int)(FluxDensityPullList_compact_fit_nnsel[i].size());
			Caesar::BoxStats<double> stats_fluxpull= StatsUtils::ComputeBoxStats(FluxDensityPullList_compact_fit_nnsel[i]);	
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit_nnsel[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}

			double fluxpull_mean= 0;
			double fluxpull_rms= 0;
			Caesar::StatsUtils::ComputeMeanAndRMS(fluxpull_mean,fluxpull_rms,FluxDensityPullList_compact_fit_nnsel[i]);
			double fluxpull_median= stats_fluxpull.median;
			double fluxpull_median_err= 1.253*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr= stats_fluxpull.Q3-stats_fluxpull.Q1;
			double fluxpull_iqr_err= 1.573*fluxpull_rms/sqrt(nEntries);
			double fluxpull_iqr_half= fluxpull_iqr/2.;
			double fluxpull_iqr_half_err= 0.5*fluxpull_iqr_err;

			//Bias graph
			double y= fluxpull_median;
			//double yerr_low= y-stats_fluxpull.Q1;
			//double yerr_up= stats_fluxpull.Q3-y;
			double yerr_low= fluxpull_median_err;
			double yerr_up= fluxpull_median_err;
			FluxDensityAccuracyGraph_compact_fit_nnsel->SetPoint(nPoints_FluxPull_fit_nnsel,x,y);
			FluxDensityAccuracyGraph_compact_fit_nnsel->SetPointError(nPoints_FluxPull_fit_nnsel,xerr_low,xerr_up,yerr_low,yerr_up);

			//Reso graph
			y= fluxpull_iqr_half;
			yerr_low= fluxpull_iqr_half_err;
			yerr_up= fluxpull_iqr_half_err;
			FluxDensityResolutionGraph_compact_fit_nnsel->SetPoint(nPoints_FluxPull_fit_nnsel,x,y);
			FluxDensityResolutionGraph_compact_fit_nnsel->SetPointError(nPoints_FluxPull_fit_nnsel,xerr_low,xerr_up,yerr_low,yerr_up);
			
			nPoints_FluxPull_fit_nnsel++;
		}	
		

		

		
	}//end loop bins

	//Fill graph for ext
	for(int i=0;i<nBins_ext;i++){
		double x_avg= LgFluxBins_ext[i] + 0.5*(LgFluxBins_ext[i+1]-LgFluxBins_ext[i]);



		//Compute Flux accuracy for extended sources
		if(FluxDensityPullList_ext[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(FluxDensityPullList_ext[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_ext[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			FluxDensityAccuracyGraph_ext->SetPoint(nPoints_FluxPull_ext,x,y);
			FluxDensityAccuracyGraph_ext->SetPointError(nPoints_FluxPull_ext,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_FluxPull_ext++;
		}	

		

		//Compute xpos accuracy	extended sources	
		if(xPosPullList_ext[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(xPosPullList_ext[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_ext[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			xPosAccuracyGraph_ext->SetPoint(nPoints_xPull_ext,x,y);
			xPosAccuracyGraph_ext->SetPointError(nPoints_xPull_ext,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_xPull_ext++;
		}

		//Compute ypos accuracy extended sources
		if(yPosPullList_ext[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(yPosPullList_ext[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_ext[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			if(!opt.drawErrorX){
				xerr_low= 0;
				xerr_up= 0;
			}
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			yPosAccuracyGraph_ext->SetPoint(nPoints_yPull_ext,x,y);
			yPosAccuracyGraph_ext->SetPointError(nPoints_yPull_ext,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_yPull_ext++;
		}	

	}//end loop bin ext 

}//close ComputeAnalysisHistos()





int FillAnalysisHisto()
{
	//Set branch addresses in trees	
	PointSourcePerfInfo->SetBranchAddress("sourcePosX",&sourcePosX);
	PointSourcePerfInfo->SetBranchAddress("sourcePosY",&sourcePosY);
	PointSourcePerfInfo->SetBranchAddress("sourceTruePosX",&sourceTruePosX);
	PointSourcePerfInfo->SetBranchAddress("sourceTruePosY",&sourceTruePosY);
	PointSourcePerfInfo->SetBranchAddress("sourceFlux",&sourceFlux);
	PointSourcePerfInfo->SetBranchAddress("sourceTrueFlux",&sourceTrueFlux);
	PointSourcePerfInfo->SetBranchAddress("thetaDiff",&dTheta);
	PointSourcePerfInfo->SetBranchAddress("eccentricityRatio",&eccentricityRatio);
	PointSourcePerfInfo->SetBranchAddress("sourceToBeamRatio",&areaToBeamRatio);
	PointSourcePerfInfo->SetBranchAddress("sourceSNR",&sourceSNR);
	PointSourcePerfInfo->SetBranchAddress("sourceTrueSNR",&sourceTrueSNR);
	PointSourcePerfInfo->SetBranchAddress("fitChi2",&fitChi2);
	PointSourcePerfInfo->SetBranchAddress("fitNDF",&fitNDF);
	PointSourcePerfInfo->SetBranchAddress("fitQuality",&fitQuality);
	PointSourcePerfInfo->SetBranchAddress("fitComponentFlag",&fitComponentFlag);

	RecPointSourcePerfInfo->SetBranchAddress("sourcePosX",&sourcePosX);
	RecPointSourcePerfInfo->SetBranchAddress("sourcePosY",&sourcePosY);
	RecPointSourcePerfInfo->SetBranchAddress("sourceFlux",&sourceFlux);
	RecPointSourcePerfInfo->SetBranchAddress("thetaDiff",&dTheta);
	RecPointSourcePerfInfo->SetBranchAddress("eccentricityRatio",&eccentricityRatio);
	RecPointSourcePerfInfo->SetBranchAddress("sourceToBeamRatio",&areaToBeamRatio);
	RecPointSourcePerfInfo->SetBranchAddress("sourceSNR",&sourceSNR);
	RecPointSourcePerfInfo->SetBranchAddress("isTrueSource",&isTrueSource);
	RecPointSourcePerfInfo->SetBranchAddress("fitChi2",&fitChi2);
	RecPointSourcePerfInfo->SetBranchAddress("fitNDF",&fitNDF);
	RecPointSourcePerfInfo->SetBranchAddress("fitQuality",&fitQuality);
	RecPointSourcePerfInfo->SetBranchAddress("fitComponentFlag",&fitComponentFlag);

	//Loop events and fill histos
	INFO_LOG("Filling histos ...");
	long int nSources= PointSourcePerfInfo->GetEntries();
	long int nSources_presel= 0;
	long int nSources_cutsel= 0;
	long int nSources_nnsel= 0;
	for(int i=0;i<PointSourcePerfInfo->GetEntries();i++){
		PointSourcePerfInfo->GetEntry(i);

		//Get data 
		double lgFlux_true= log10(sourceTrueFlux);
		double xOffset= sourcePosX-sourceTruePosX;
		double yOffset= sourcePosY-sourceTruePosY;
		double FluxPull= sourceFlux/sourceTrueFlux-1;
		double redChi2= fitChi2/fitNDF;
	
		//## Apply minimal cuts
		//- Select flux bin
		int gBin= NTrueSourceHisto_compact->FindBin(lgFlux_true);
		if(NTrueSourceHisto_compact->IsBinUnderflow(gBin) || NTrueSourceHisto_compact->IsBinOverflow(gBin)) continue;


		//Fill histos & vectors for minimal selection data
		NRecSourceHisto_compact_fit->Fill(lgFlux_true,1);
		//NRecSourceVSSignificanceHisto_compact_fit->Fill(sourceTrueSNR,1);
		FluxList_compact_fit[gBin-1].push_back(lgFlux_true);
		xPosPullList_compact_fit[gBin-1].push_back(xOffset);
		yPosPullList_compact_fit[gBin-1].push_back(yOffset);
		FluxDensityPullList_compact_fit[gBin-1].push_back(FluxPull);
		
	
		//## Apply cut & NN selection?
		if(!opt.selectFitComponents) continue;

		//## Apply pre-selection cuts
		//- Apply chi2 cut	
		if(opt.applyChi2Cut && redChi2>=opt.chi2Cut) {
			INFO_LOG("Skipping source with lgFlux_true="<<lgFlux_true<<" as not passing chi2 cut (redChi2="<<redChi2<<">="<<opt.chi2Cut<<")");
			continue;
		}

		//- Apply fit quality cut
		if(opt.applyFitQualityCut && fitQuality<opt.fitQualityCut){
			INFO_LOG("Skipping source with lgFlux_true="<<lgFlux_true<<" as not passing fitQuality cut (fitQuality="<<fitQuality<<"<"<<opt.fitQualityCut<<")");
			continue;
		}

		//- Apply fit ellipse cuts
		if(opt.applyFitEllipseCuts){
			bool isGoodTheta= (dTheta>opt.deltaFitThetaMinThr && dTheta<opt.deltaFitThetaMaxThr);
			bool isGoodEccentricity= (eccentricityRatio>opt.eccentricityRatioMinThr && eccentricityRatio<opt.eccentricityRatioMaxThr);
			bool isGoodAreaToBeam= (areaToBeamRatio>opt.areaToBeamRatioMinThr && areaToBeamRatio<opt.areaToBeamRatioMaxThr);
			bool passed= (isGoodTheta && isGoodEccentricity && isGoodAreaToBeam);
			if(!passed) {
				INFO_LOG("Skipping source with lgFlux_true="<<lgFlux_true<<" as not passing fit ellipse cut");
				continue;		
			}
		}	

		nSources_presel++;
		NRecSourceHisto_compact_fit_presel->Fill(lgFlux_true,1);
		FluxList_compact_fit_presel[gBin-1].push_back(lgFlux_true);
		xPosPullList_compact_fit_presel[gBin-1].push_back(xOffset);
		yPosPullList_compact_fit_presel[gBin-1].push_back(yOffset);
		FluxDensityPullList_compact_fit_presel[gBin-1].push_back(FluxPull);

		//## Apply advanced cut selection
		if(opt.applySourceClassifier){

			//- Apply rectangular cut selection
			bool passed= true;
			if(opt.applyUserRecCuts){
				bool isGoodTheta= (dTheta>opt.deltaFitThetaMinThr && dTheta<opt.deltaFitThetaMaxThr);
				bool isGoodEccentricity= (eccentricityRatio>opt.eccentricityRatioMinThr && eccentricityRatio<opt.eccentricityRatioMaxThr);
				bool isGoodAreaToBeam= (areaToBeamRatio>opt.areaToBeamRatioMinThr && areaToBeamRatio<opt.areaToBeamRatioMaxThr);
				passed= (isGoodTheta && isGoodEccentricity && isGoodAreaToBeam);
			}
			else{
				passed = reader->EvaluateMVA("Cuts",opt.signalCutEff);
			}

      if(passed) {
				nSources_cutsel++;
				NRecSourceHisto_compact_fit_cutsel->Fill(lgFlux_true,1);
				//NRecSourceVSSignificanceHisto_compact_fit_cutsel->Fill(sourceTrueSNR,1);
	
				FluxList_compact_fit_cutsel[gBin-1].push_back(lgFlux_true);
				xPosPullList_compact_fit_cutsel[gBin-1].push_back(xOffset);
				yPosPullList_compact_fit_cutsel[gBin-1].push_back(yOffset);
				FluxDensityPullList_compact_fit_cutsel[gBin-1].push_back(FluxPull);
			}

			//- Apply NN selection
			double NNOut= reader_NN->EvaluateMVA("MLP");
			passed= (NNOut>=opt.nnCut);
      if(passed) {
				nSources_nnsel++;
				NRecSourceHisto_compact_fit_nnsel->Fill(lgFlux_true,1);
				//NRecSourceVSSignificanceHisto_compact_fit_nnsel->Fill(sourceTrueSNR,1);
	
				FluxList_compact_fit_nnsel[gBin-1].push_back(lgFlux_true);
				xPosPullList_compact_fit_nnsel[gBin-1].push_back(xOffset);
				yPosPullList_compact_fit_nnsel[gBin-1].push_back(yOffset);
				FluxDensityPullList_compact_fit_nnsel[gBin-1].push_back(FluxPull);
			}
		}//close if 

		
	}//end loop events

	INFO_LOG("True detected sources after pre-selection cuts: sel/tot="<<nSources_presel<<"/"<<nSources);
	INFO_LOG("True detected sources after cut selection: sel/tot="<<nSources_cutsel<<"/"<<nSources);
	INFO_LOG("True detected sources after NN selection: sel/tot="<<nSources_nnsel<<"/"<<nSources);

	INFO_LOG("Filling rec histos ...");
	long int nRecSources= RecPointSourcePerfInfo->GetEntries();
	long int nRecSources_presel= 0;
	long int nRecSources_cutsel= 0;
	long int nRecSources_nnsel= 0;
	long int nFitSources= 0;
	long int nFitSources_fake= 0;
	
	for(int i=0;i<RecPointSourcePerfInfo->GetEntries();i++){
		RecPointSourcePerfInfo->GetEntry(i);

		//Get data 
		double lgFlux_rec= log10(sourceFlux);
		double Z_rec= sourceFlux/opt.noiseLevel_true;
		double redChi2= fitChi2/fitNDF;

		//## Apply minimal cuts

		//Fill histos
		NRecSourceHisto_reliability_compact_fit->Fill(lgFlux_rec,1);
		if(isTrueSource==1){
			NTrueSourceHisto_reliability_compact_fit->Fill(lgFlux_rec,1);	
			nFitSources++;
		}
		else{
			nFitSources_fake++;
		}
							
		//Apply selection?
		if(!opt.selectFitComponents) continue;

		//## Apply pre-selection cuts
		//- Apply chi2 cut	
		if(opt.applyChi2Cut && redChi2>=opt.chi2Cut) continue;

		//- Apply fit quality cut
		if(opt.applyFitQualityCut && fitQuality<opt.fitQualityCut){
			continue;
		}

		//- Apply fit ellipse cuts
		if(opt.applyFitEllipseCuts){
			bool isGoodTheta= (dTheta>opt.deltaFitThetaMinThr && dTheta<opt.deltaFitThetaMaxThr);
			bool isGoodEccentricity= (eccentricityRatio>opt.eccentricityRatioMinThr && eccentricityRatio<opt.eccentricityRatioMaxThr);
			bool isGoodAreaToBeam= (areaToBeamRatio>opt.areaToBeamRatioMinThr && areaToBeamRatio<opt.areaToBeamRatioMaxThr);
			bool passed= (isGoodTheta && isGoodEccentricity && isGoodAreaToBeam);
			if(!passed) continue;
		}	

		nRecSources_presel++;
		NRecSourceHisto_reliability_compact_fit_presel->Fill(lgFlux_rec,1);
		if(isTrueSource==1){
			NTrueSourceHisto_reliability_compact_fit_presel->Fill(lgFlux_rec,1);	
		}

		//## Apply advanced cut selection
		if(opt.applySourceClassifier){
			//Apply cut selection
			bool passed= true;
			if(opt.applyUserRecCuts){
				bool isGoodTheta= (dTheta>opt.deltaFitThetaMinThr && dTheta<opt.deltaFitThetaMaxThr);
				bool isGoodEccentricity= (eccentricityRatio>opt.eccentricityRatioMinThr && eccentricityRatio<opt.eccentricityRatioMaxThr);
				bool isGoodAreaToBeam= (areaToBeamRatio>opt.areaToBeamRatioMinThr && areaToBeamRatio<opt.areaToBeamRatioMaxThr);
				passed= (isGoodTheta && isGoodEccentricity && isGoodAreaToBeam);
			}
			else{
				passed = reader->EvaluateMVA("Cuts",opt.signalCutEff);
			}
			
      if(passed) {
				nRecSources_cutsel++;
				NRecSourceHisto_reliability_compact_fit_cutsel->Fill(lgFlux_rec,1);
				if(isTrueSource==1){
					NTrueSourceHisto_reliability_compact_fit_cutsel->Fill(lgFlux_rec,1);	
				}
			}


			//Apply NN selection
			double NNOut= reader_NN->EvaluateMVA("MLP");
			passed= (NNOut>=opt.nnCut);
      if(passed) {
				nRecSources_nnsel++;
				NRecSourceHisto_reliability_compact_fit_nnsel->Fill(lgFlux_rec,1);
				if(isTrueSource==1){
					NTrueSourceHisto_reliability_compact_fit_nnsel->Fill(lgFlux_rec,1);	
				}
			}
		}//close if 

		
	}//end loop events

	INFO_LOG("Detected sources after preselection cuts: sel/tot="<<nRecSources_presel<<"/"<<nRecSources);
	INFO_LOG("Detected sources cut selection: sel/tot="<<nRecSources_cutsel<<"/"<<nRecSources);
	INFO_LOG("Detected sources NN selection: sel/tot="<<nRecSources_nnsel<<"/"<<nRecSources);


	return 0;

}//close FillAnalysisHisto()


int AnalyzeData(std::string filename){

	INFO_LOG("Analyzing input data "<<filename<<" ...");

	sourceFileName= CodeUtils::ExtractFileNameFromPath(filename,false);

	//Define vars
	int found;
	std::string* Name_true= 0;
	int Type_true;
	int SimType_true;
	long long NPix;
	double X0_true;//true pos x
	double Y0_true;//true pos y
	double S_true;//true peak amplitude
	double X0;//pos x
	double Y0;//pos y
	double S;
	double beamArea_true;
	double X0_sweighted;
	double Y0_sweighted;

	std::string* Name_rec= 0;
	long long NPix_rec;
	int Type_rec;
	double S_rec;
	double Smax_rec;//rec peak amplitude
	double X0_rec;//rec pos x
	double Y0_rec;//rec pos y
	double fluxDensity_rec;
	double beamArea_rec;
	double X0_sweighted_rec;
	double Y0_sweighted_rec;
	double S_bkg;
	double AvgBkg;
	double AvgRMS;

	int HasFitInfo;
	double X0_fit;//fitted pos x
 	double Y0_fit;//fitted pos y
	double S_fit;//fitted peak amplitude	
	double fluxDensity_fit;//flux density reconstructed using fit pars
	double offset_fit;//fitted offset	
	double sigmaX_fit;//fitted sigmaX
	double sigmaY_fit;//fitted sigmaY
	double theta_fit;//fitted theta
	double chi2_fit;
	double ndf_fit;
	
	int nTrueMatchedSources;
	std::vector<std::string>* Names_true= 0;
	std::vector<double>* XPos_true= 0;
	std::vector<double>* YPos_true= 0;

	//Open file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile){
		ERROR_LOG("Failed to open file "<<filename<<"!");
		return -1;
	}


	//Get point/compact source match tree
	TTree* SourceMatchInfo= (TTree*)inputFile->Get("SourceMatchInfo");
	if(!SourceMatchInfo){
		ERROR_LOG("Failed to get SourceMatchInfo tree from file "<<filename<<"!");
		return -1;
	}

	SourceMatchInfo->SetBranchAddress("name",&Name_true);
	SourceMatchInfo->SetBranchAddress("type",&Type_true);
	SourceMatchInfo->SetBranchAddress("simtype",&SimType_true);
	SourceMatchInfo->SetBranchAddress("X0_true",&X0_true);
 	SourceMatchInfo->SetBranchAddress("Y0_true",&Y0_true);
	SourceMatchInfo->SetBranchAddress("X0_sweighted",&X0_sweighted);
 	SourceMatchInfo->SetBranchAddress("Y0_sweighted",&Y0_sweighted);
 	SourceMatchInfo->SetBranchAddress("S_true",&S_true);
	SourceMatchInfo->SetBranchAddress("S",&S);	
	SourceMatchInfo->SetBranchAddress("beamArea_true",&beamArea_true);
 	SourceMatchInfo->SetBranchAddress("found",&found);
	SourceMatchInfo->SetBranchAddress("name_rec",&Name_rec);	
	SourceMatchInfo->SetBranchAddress("NPix_rec",&NPix_rec);
	SourceMatchInfo->SetBranchAddress("S_rec",&S_rec);	
 	SourceMatchInfo->SetBranchAddress("Smax_rec",&Smax_rec);
 	SourceMatchInfo->SetBranchAddress("X0_rec",&X0_rec);
 	SourceMatchInfo->SetBranchAddress("Y0_rec",&Y0_rec);
	SourceMatchInfo->SetBranchAddress("X0_sweighted_rec",&X0_sweighted_rec);
 	SourceMatchInfo->SetBranchAddress("Y0_sweighted_rec",&Y0_sweighted_rec);
	SourceMatchInfo->SetBranchAddress("beamArea_rec",&beamArea_rec);
	SourceMatchInfo->SetBranchAddress("S_bkg",&S_bkg);
	SourceMatchInfo->SetBranchAddress("AvgBkg",&AvgBkg);
	SourceMatchInfo->SetBranchAddress("AvgRMS",&AvgRMS);
  SourceMatchInfo->SetBranchAddress("HasFitInfo",&HasFitInfo);
 	SourceMatchInfo->SetBranchAddress("X0_fit",&X0_fit);
 	SourceMatchInfo->SetBranchAddress("Y0_fit",&Y0_fit);
 	SourceMatchInfo->SetBranchAddress("S_fit",&S_fit);
	SourceMatchInfo->SetBranchAddress("fluxDensity_fit",&fluxDensity_fit);
	SourceMatchInfo->SetBranchAddress("offset_fit",&offset_fit);
	SourceMatchInfo->SetBranchAddress("sigmaX_fit",&sigmaX_fit);
	SourceMatchInfo->SetBranchAddress("sigmaY_fit",&sigmaY_fit);
	SourceMatchInfo->SetBranchAddress("theta_fit",&theta_fit);
	SourceMatchInfo->SetBranchAddress("chi2_fit",&chi2_fit);
	SourceMatchInfo->SetBranchAddress("ndf_fit",&ndf_fit);
	SourceMatchInfo->SetBranchAddress("fitQuality",&fitQuality);
	SourceMatchInfo->SetBranchAddress("fitComponentFlag",&fitComponentFlag);


	//Get point/compact source reliability tree
	TTree* RecSourceInfo= (TTree*)inputFile->Get("RecSourceInfo");
	if(!RecSourceInfo){
		ERROR_LOG("Failed to get RecSourceInfo tree from file "<<filename<<"!");
		return -1;
	}

	RecSourceInfo->SetBranchAddress("name_rec",&Name_rec);
	RecSourceInfo->SetBranchAddress("type_rec",&Type_rec);	
	RecSourceInfo->SetBranchAddress("HasFitInfo",&HasFitInfo);
	RecSourceInfo->SetBranchAddress("sigmaX_fit",&sigmaX_fit);
	RecSourceInfo->SetBranchAddress("sigmaY_fit",&sigmaY_fit);
	RecSourceInfo->SetBranchAddress("theta_fit",&theta_fit);	
	RecSourceInfo->SetBranchAddress("X0_rec",&X0_rec);
	RecSourceInfo->SetBranchAddress("Y0_rec",&Y0_rec);
	RecSourceInfo->SetBranchAddress("Smax_rec",&Smax_rec);
	RecSourceInfo->SetBranchAddress("fluxDensity_rec",&fluxDensity_rec);
	RecSourceInfo->SetBranchAddress("beamArea_rec",&beamArea_rec);
	RecSourceInfo->SetBranchAddress("NPix_rec",&NPix_rec);
	RecSourceInfo->SetBranchAddress("nTrueMatchedSources",&nTrueMatchedSources);	
	RecSourceInfo->SetBranchAddress("name_true",&Names_true);	
	RecSourceInfo->SetBranchAddress("X0_true",&XPos_true);	
	RecSourceInfo->SetBranchAddress("Y0_true",&YPos_true);	
	RecSourceInfo->SetBranchAddress("chi2_fit",&chi2_fit);
	RecSourceInfo->SetBranchAddress("ndf_fit",&ndf_fit);
	RecSourceInfo->SetBranchAddress("fitQuality",&fitQuality);
	RecSourceInfo->SetBranchAddress("fitComponentFlag",&fitComponentFlag);

	//Get extended source match tree
	TTree* ExtSourceMatchInfo= (TTree*)inputFile->Get("ExtSourceMatchInfo");
	if(!ExtSourceMatchInfo){
		ERROR_LOG("Failed to get ExtSourceMatchInfo tree from file "<<filename<<"!");
		return -1;
	}

	ExtSourceMatchInfo->SetBranchAddress("name",&Name_true);
	ExtSourceMatchInfo->SetBranchAddress("type",&Type_true);
	ExtSourceMatchInfo->SetBranchAddress("simtype",&SimType_true);
	ExtSourceMatchInfo->SetBranchAddress("NPix",&NPix);
	ExtSourceMatchInfo->SetBranchAddress("X0_true",&X0_true);
 	ExtSourceMatchInfo->SetBranchAddress("Y0_true",&Y0_true);
 	ExtSourceMatchInfo->SetBranchAddress("S_true",&S_true);	
	ExtSourceMatchInfo->SetBranchAddress("X0",&X0);
 	ExtSourceMatchInfo->SetBranchAddress("Y0",&Y0);
	ExtSourceMatchInfo->SetBranchAddress("X0_sweighted",&X0_sweighted);
 	ExtSourceMatchInfo->SetBranchAddress("Y0_sweighted",&Y0_sweighted);
	ExtSourceMatchInfo->SetBranchAddress("S",&S);
	ExtSourceMatchInfo->SetBranchAddress("beamArea_true",&beamArea_true);
 	ExtSourceMatchInfo->SetBranchAddress("found",&found);
	ExtSourceMatchInfo->SetBranchAddress("name_rec",&Name_rec);	
	ExtSourceMatchInfo->SetBranchAddress("NPix_rec",&NPix_rec);
	ExtSourceMatchInfo->SetBranchAddress("S_rec",&S_rec);	
 	ExtSourceMatchInfo->SetBranchAddress("Smax_rec",&Smax_rec);
 	ExtSourceMatchInfo->SetBranchAddress("X0_rec",&X0_rec);
 	ExtSourceMatchInfo->SetBranchAddress("Y0_rec",&Y0_rec);
	ExtSourceMatchInfo->SetBranchAddress("beamArea_rec",&beamArea_rec);
	ExtSourceMatchInfo->SetBranchAddress("X0_sweighted_rec",&X0_sweighted_rec);
 	ExtSourceMatchInfo->SetBranchAddress("Y0_sweighted_rec",&Y0_sweighted_rec);
	ExtSourceMatchInfo->SetBranchAddress("S_bkg",&S_bkg);
	ExtSourceMatchInfo->SetBranchAddress("AvgBkg",&AvgBkg);
	ExtSourceMatchInfo->SetBranchAddress("AvgRMS",&AvgRMS);

	/*
	//Get point/compact source reliability tree
	TTree* RecExtSourceInfo= (TTree*)inputFile->Get("RecExtSourceInfo");
	if(!RecSourceInfo){
		ERROR_LOG("Failed to get RecExtSourceInfo tree from file "<<filename<<"!");
		return -1;
	}

	RecExtSourceInfo->SetBranchAddress("name_rec",&Name_rec);
	RecExtSourceInfo->SetBranchAddress("type_rec",&Type_rec);	
	RecExtSourceInfo->SetBranchAddress("HasFitInfo",&HasFitInfo);	
	RecExtSourceInfo->SetBranchAddress("X0_rec",&X0_rec);
	RecExtSourceInfo->SetBranchAddress("Y0_rec",&Y0_rec);
	RecExtSourceInfo->SetBranchAddress("Smax_rec",&Smax_rec);
	RecExtSourceInfo->SetBranchAddress("fluxDensity_rec",&fluxDensity_rec);
	RecExtSourceInfo->SetBranchAddress("beamArea_rec",&beamArea_rec);
	RecExtSourceInfo->SetBranchAddress("NPix_rec",&NPix_rec);
	RecExtSourceInfo->SetBranchAddress("nTrueMatchedSources",&nTrueMatchedSources);	
	RecExtSourceInfo->SetBranchAddress("name_true",&Names_true);	
	RecExtSourceInfo->SetBranchAddress("X0_true",&XPos_true);	
	RecExtSourceInfo->SetBranchAddress("Y0_true",&YPos_true);	
	*/

	//============================================================
	//===       SELECT TRUE SOURCES
	//============================================================	
	//Select true point-like sources to be included in the analysis?
	std::vector<bool> trueSourceSelectionFlags(SourceMatchInfo->GetEntries(),true);
	std::map<std::string,bool> trueSourceSelectionFlagMap;
	
	//Init flags
	for(int i=0;i<SourceMatchInfo->GetEntries();i++){
		SourceMatchInfo->GetEntry(i);		
		std::string sourceName_i= *Name_true;
	
		trueSourceSelectionFlags[i]= true;
		trueSourceSelectionFlagMap[sourceName_i]= true;
	}//end loop sources

	//Set flags
	for(int i=0;i<SourceMatchInfo->GetEntries();i++){
		//Get true source data
		SourceMatchInfo->GetEntry(i);		
		std::string sourceName_i= *Name_true;
		int sourceType_i= Type_true;
		double X_i= X0_true;
		double Y_i= Y0_true;
		

		//Skip if source is at edge?	
		bool atEdgeX= (X_i<opt.sourceRegionXMin || X_i>opt.sourceRegionXMax);
		bool atEdgeY= (Y_i<opt.sourceRegionYMin || Y_i>opt.sourceRegionYMax);
		bool atEdge= (atEdgeX || atEdgeY);
		if(opt.excludeSourceAtEdge && atEdge) {
			INFO_LOG("Skipping true source "<<sourceName_i<<" at pos("<<X_i<<","<<Y_i<<") as found at border...");
			trueSourceSelectionFlags[i]= false;
			trueSourceSelectionFlagMap[sourceName_i]= false;
			continue;		
		}

		//Search for sources too close to be resolved
		if(opt.selectResolvedTrueSources){
			for(int j=i+1;j<SourceMatchInfo->GetEntries();j++){
				SourceMatchInfo->GetEntry(j);
				std::string sourceName_j= *Name_true;
				int sourceType_j= Type_true;
				double X_j= X0_true;
				double Y_j= Y0_true;
				double distX= fabs(X_i-X_j);
				double distY= fabs(Y_i-Y_j);
				double dist= sqrt( (X_i-X_j)*(X_i-X_j) + (Y_i-Y_j)*(Y_i-Y_j) );

				//Check if they are both point-like sources
				if(sourceType_i!=ePointLike || sourceType_j!=ePointLike) continue;

				//if(distX<=opt.mutualTrueSourceDistThr/opt.pixSize && distY<=opt.mutualTrueSourceDistThr/opt.pixSize){//mark both sources as not selected (not resolved)
				if(dist<opt.mutualTrueSourceDistThr/opt.pixSize){//mark both sources as not selected (not resolved)
					INFO_LOG("Skipping true sources "<<sourceName_i<<"-"<<sourceName_j<<" at positions pos("<<X_i<<","<<Y_i<<")-("<<X_j<<","<<Y_j<<") as too close given tolerance ("<<opt.mutualTrueSourceDistThr/opt.pixSize<<" pix) ...");
					trueSourceSelectionFlags[i]= false;
					trueSourceSelectionFlags[j]= false;
					trueSourceSelectionFlagMap[sourceName_i]= false;
					trueSourceSelectionFlagMap[sourceName_j]= false;
				}
			}//end loop next sources
		}//close if selectResolvedTrueSources
	}//end loop sources
	



	//============================================================
	//===    COMPUTE POINT/COMPACT SOURCE FINDING EFFICIENCY
	//============================================================	
	//Loop over source entries
	INFO_LOG("Loop over source entries");
	int nTrueSources= 0;
	int nFitSources= 0;
	TEllipse* fitEllipse= new TEllipse(0,0,1,1,0,360,0);

	for(int i=0;i<SourceMatchInfo->GetEntries();i++){
		SourceMatchInfo->GetEntry(i);		

		
		//Fill true info
		std::string sourceName_true= *Name_true;
		double fluxDensity_true= S_true;
		double lgFlux_true= log10(fluxDensity_true);
		double Z_true= fluxDensity_true/opt.noiseLevel_true;
		//double Z_true= fluxDensity_true/AvgRMS;
		
		//Apply selection cuts to true source
		//- Skip true source if overlapping with another true source
		if(!trueSourceSelectionFlags[i]) {
			INFO_LOG("Skipping true source "<<sourceName_true<<" as overlapping with another true source within chosen tolerance...");
			continue;
		}		

		//Fill true histo
		int gBin= NTrueSourceHisto_compact->FindBin(lgFlux_true);
		if(NTrueSourceHisto_compact->IsBinUnderflow(gBin) || NTrueSourceHisto_compact->IsBinOverflow(gBin)) continue;
		NTrueSourceHisto_compact->Fill(lgFlux_true,1);
		//NTrueSourceVSSignificanceHisto_compact->Fill(Z_true,1);
		nTrueSources++;

		//Skip if true source is not detected
		if(found==0) continue;

		//Skip if no fit source are present
		if(!HasFitInfo) continue;


		//## Set Tree pars
		double fluxDensity= fluxDensity_fit/beamArea_rec;	
		double Z_rec= fluxDensity/opt.noiseLevel_true;
		sourceName= *Name_rec;
		sourceFlux= fluxDensity;
		sourcePosX= X0_fit;
		sourcePosY= Y0_fit;
		sourceSNR= Z_rec;
		sourceTrueFlux= fluxDensity_true;
		sourceTrueSNR= Z_true;
		sourceTruePosX= X0_true;
		sourceTruePosY= Y0_true;
		fitChi2= chi2_fit;
		fitNDF= ndf_fit;	
			

		//Get fit ellipse pars
		//ellMaj= std::max(sigmaX_fit,sigmaY_fit);
		//ellMin= std::min(sigmaX_fit,sigmaY_fit);
		ellMaj= std::max(sigmaX_fit,sigmaY_fit)*GausSigma2FWHM*opt.pixSize;//in arcsec
		ellMin= std::min(sigmaX_fit,sigmaY_fit)*GausSigma2FWHM*opt.pixSize;//in arcsec
		//ellTheta= MathUtils::Mod(theta_fit,180.);
		ellTheta= theta_fit-90.;
		//if(sigmaY_fit>sigmaX_fit) ellTheta= MathUtils::Mod(theta_fit-90.,180.);
		if(sigmaY_fit>sigmaX_fit) ellTheta+= 90.;
		ellTheta= GetPosAngleInRange(ellTheta);
		fitEllipse->SetR1(ellMaj);
		fitEllipse->SetR2(ellMin);
		fitEllipse->SetTheta(ellTheta);
		//ellArea= MathUtils::ComputeEllipseArea(fitEllipse);
		ellArea= MathUtils::ComputeEllipseArea(ellMaj,ellMin);
		//ellEccentricity= MathUtils::ComputeEllipseEccentricity(fitEllipse);
		ellEccentricity= MathUtils::ComputeEllipseEccentricity(ellMaj,ellMin);
		//dTheta= fitEllipse->GetTheta()-beamEllipse->GetTheta();
		dTheta= ellTheta-beamTheta;	
		eccentricityRatio= ellEccentricity/beamEllipseEccentricity;
		areaToBeamRatio= ellArea/beamEllipseArea;
		isTrueSource= 1;	
		
		//## Fill tree
		PointSourcePerfInfo->Fill();
	
	}//end loop sources

	INFO_LOG("Selected "<<nTrueSources<<"/"<<SourceMatchInfo->GetEntries()<<" true compact sources for analysis...");

	//============================================================
	//===    COMPUTE COMPACT SOURCE FINDING RELIABILITY
	//============================================================	
	//## Loop over rec sources
	nFitSources= 0;
	int nFitSources_fake= 0;
	for(int i=0;i<RecSourceInfo->GetEntries();i++){
		RecSourceInfo->GetEntry(i);
		
		//## Apply selection cuts
		//- Skip if not compact source
		bool isPointSource= (Type_rec==ePointLike);
		bool isCompactSource= (Type_rec==eCompact);
		//if(!isPointSource) continue;	
		if(!isPointSource && !isCompactSource) continue;

		//- Skip if source is at edge?	
		bool atEdgeX= (X0_rec<opt.sourceRegionXMin || X0_rec>opt.sourceRegionXMax);
		bool atEdgeY= (Y0_rec<opt.sourceRegionYMin || Y0_rec>opt.sourceRegionYMax);
		bool atEdge= (atEdgeX || atEdgeY);
		if(opt.excludeSourceAtEdge && atEdge) continue;

		//- Skip if source has no fit info
		if(!HasFitInfo) continue;

		//- Skip if source has matches with excluded true sources
		if(nTrueMatchedSources>0){
			int nTrueSourceMatch_sel= 0;
			std::vector<std::string> sourceNames_true_removed;
			for(size_t k=0;k<Names_true->size();k++){
				std::string sourceName_true= Names_true->at(k);
				bool isTrueSourceSelected= trueSourceSelectionFlagMap[sourceName_true];
				if(isTrueSourceSelected) nTrueSourceMatch_sel++;
				else sourceNames_true_removed.push_back(sourceName_true);
			}//end loop true source match

			if(nTrueSourceMatch_sel<=0){
				INFO_LOG("Skipping rec source "<<*Name_rec<<" as it matches with true sources that have been excluded from the analysis...");
				continue;
			}
			if(nTrueSourceMatch_sel!=nTrueMatchedSources){
				DEBUG_LOG("Removed "<<nTrueMatchedSources-nTrueSourceMatch_sel<<"/"<<nTrueMatchedSources<<" true sources from rec source "<<*Name_rec<<" as it matches with true sources that have been excluded from the analysis...");
				for(size_t k=0;k<sourceNames_true_removed.size();k++){
					DEBUG_LOG("True source removed="<<sourceNames_true_removed[k]);
				}
			}

		}//close if


		//## Set tree data
		double lgFlux_rec= log10(fluxDensity_rec/beamArea_rec);
		double Z_rec= (fluxDensity_rec/beamArea_rec)/opt.noiseLevel_true;
		sourceName= *Name_rec;
		sourceFlux= fluxDensity_rec/beamArea_rec;
		sourceSNR= Z_rec;
		sourcePosX= X0_rec;
		sourcePosY= Y0_rec;
		isTrueSource= 0;
		if(nTrueMatchedSources>0) isTrueSource= 1;

		fitChi2= chi2_fit;
		fitNDF= ndf_fit;

		//Compute ellipse pars
		//ellMaj= std::max(sigmaX_fit,sigmaY_fit);
		//ellMin= std::min(sigmaX_fit,sigmaY_fit);
		ellMaj= std::max(sigmaX_fit,sigmaY_fit)*GausSigma2FWHM*opt.pixSize;//in arcsec
		ellMin= std::min(sigmaX_fit,sigmaY_fit)*GausSigma2FWHM*opt.pixSize;//in arcsec
		//ellTheta= MathUtils::Mod(theta_fit,180.);
		ellTheta= theta_fit-90.;
		//if(sigmaY_fit>sigmaX_fit) ellTheta= MathUtils::Mod(theta_fit-90.,180.);
		if(sigmaY_fit>sigmaX_fit) ellTheta+= 90.;
		ellTheta= GetPosAngleInRange(ellTheta);
		fitEllipse->SetR1(ellMaj);
		fitEllipse->SetR2(ellMin);
		fitEllipse->SetTheta(ellTheta);
		//ellArea= MathUtils::ComputeEllipseArea(fitEllipse);
		ellArea= MathUtils::ComputeEllipseArea(ellMaj,ellMin);
		//ellEccentricity= MathUtils::ComputeEllipseEccentricity(fitEllipse);
		ellEccentricity= MathUtils::ComputeEllipseEccentricity(ellMaj,ellMin);
		eccentricityRatio= ellEccentricity/beamEllipseEccentricity;
		//dTheta= fitEllipse->GetTheta()-beamEllipse->GetTheta();	
		dTheta= ellTheta-beamTheta;		
		areaToBeamRatio= ellArea/beamEllipseArea;

		//Fill tree
		FitEllipseParTree->Fill();
		RecPointSourcePerfInfo->Fill();

	}//end loop rec sources


	//============================================================
	//===    COMPUTE EXTENDED SOURCE FINDING EFFICIENCY
	//============================================================	
	//Loop over source entries
	int nTrueSources_ext= 0;
	for(int i=0;i<ExtSourceMatchInfo->GetEntries();i++){
		ExtSourceMatchInfo->GetEntry(i);		

		//Fill true info
		double fluxDensity_true= S/beamArea_true;
		double lgFlux_true= log10(fluxDensity_true);

		double nBeams_true= (double)(NPix)/beamArea_true;
		double Z_true= fluxDensity_true/(opt.noiseLevel_true*sqrt(nBeams_true));
		//double Z_true= fluxDensity_true/(AvgRMS*sqrt(nBeams_true));
		
		int gBin= NTrueSourceHisto_ext->FindBin(lgFlux_true);
		if(NTrueSourceHisto_ext->IsBinUnderflow(gBin) || NTrueSourceHisto_ext->IsBinOverflow(gBin)) continue;	
		NTrueSourceHisto_ext->Fill(lgFlux_true,1);
		//NTrueSourceVSSignificanceHisto_ext->Fill(Z_true,1);
		
		int simtype_index= SimTypeToIndexMap[SimType_true];
		NTrueSourceHisto_simtypes_ext[simtype_index]->Fill(lgFlux_true,1);
		
		nTrueSources_ext++;

		//Skip if true source is not detected
		if(found==0) continue;

		//Fill rec info	
		NRecSourceHisto_ext->Fill(lgFlux_true,1);
		//NRecSourceVSSignificanceHisto_ext->Fill(Z_true,1);
		NRecSourceHisto_simtypes_ext[simtype_index]->Fill(lgFlux_true,1);
		
		//double xOffset= fabs(X0_rec-X0);
		//double yOffset= fabs(Y0_rec-Y0);
		//double xOffset= fabs(X0_sweighted_rec-X0_sweighted);
		//double yOffset= fabs(Y0_sweighted_rec-Y0_sweighted);
		double xOffset= X0_sweighted_rec-X0_sweighted;
		double yOffset= Y0_sweighted_rec-Y0_sweighted;
		double fluxDensity= S_rec;
		double fluxPull= S_rec/S-1;

		FluxList_ext[gBin-1].push_back(lgFlux_true);
		xPosPullList_ext[gBin-1].push_back(xOffset);
		yPosPullList_ext[gBin-1].push_back(yOffset);
		FluxDensityPullList_ext[gBin-1].push_back(fluxPull);
		FluxSignificanceList_ext[gBin-1].push_back(Z_true);
						
	}//end loop sources

	INFO_LOG("Selected "<<nTrueSources_ext<<"/"<<ExtSourceMatchInfo->GetEntries()<<" true extended sources for analysis...");

	//============================================================
	//===    COMPUTE EXTENDED SOURCE FINDING RELIABILITY
	//============================================================	
	for(int i=0;i<RecSourceInfo->GetEntries();i++){
		RecSourceInfo->GetEntry(i);
		
		//Skip if not extended source
		bool isExtendedSource= (Type_rec==eExtended);
		bool isCompactPlusExtendedSource= (Type_rec==eCompactPlusExtended);
		if(!isExtendedSource && !isCompactPlusExtendedSource) continue;

		//Skip if source is at edge?	
		bool atEdgeX= (X0_rec<opt.sourceRegionXMin || X0_rec>opt.sourceRegionXMax);
		bool atEdgeY= (Y0_rec<opt.sourceRegionYMin || Y0_rec>opt.sourceRegionYMax);
		bool atEdge= (atEdgeX || atEdgeY);
		//if(opt.excludeSourceAtEdge && atEdge) continue;

		//Fill rec info
		double lgFlux_rec= log10(fluxDensity_rec/beamArea_rec);
		//double nBeams_rec= (double)(NPix_rec)/beamArea_rec;
		double nBeams_rec= 1;//ADD NPIXRec to ROOT TTRee
		double Z_rec= (fluxDensity_rec/beamArea_rec)/(opt.noiseLevel_true*sqrt(nBeams_rec));
		NRecSourceHisto_reliability_ext->Fill(lgFlux_rec,1);
		//NRecSourceVSSignificanceHisto_reliability_ext->Fill(Z_rec,1);

		//Fill true info
		if(nTrueMatchedSources>0){
			NTrueSourceHisto_reliability_ext->Fill(lgFlux_rec,1);
			//NTrueSourceVSSignificanceHisto_reliability_ext->Fill(Z_rec,1);
		}

	}//end loop rec sources

	return 0;

}//close AnalyzeData()


void Draw(){

	gROOT->SetStyle("myStyle2");
	//gStyle->SetOptLogx();
	//gStyle->SetErrorX(0.0001);

	//double xMin_draw= LgFluxBins[0];
	//double xMax_draw= LgFluxBins[LgFluxBins.size()-1];
	double xMin_draw= -4;
	double xMax_draw= 0.5;
	double ZMin_draw= Zbins[0];
	double ZMax_draw= Zbins[Zbins.size()-1];

	//double zThr= 5;//(S-bkg)/noise= S/noise-bkg/noise= 5
	//double SNThr= zThr + bkgLevel_true/noiseLevel_true;
	//TLine* detectionThrLine= new TLine(SNThr,0,SNThr,1);
	//detectionThrLine->SetLineColor(kBlack);
	//detectionThrLine->SetLineStyle(kDashed);

	std::vector<double> zRefThrList {5,10,50};
	std::vector<TLine*> detectionThrLines;
	TLine* detectionThrLine= 0;
	for(size_t i=0;i<zRefThrList.size();i++){
		double zRefThr= zRefThrList[i];
		double fluxThrRef= zRefThr*opt.noiseLevel_true;
		double lgFluxThrRef= log10(fluxThrRef);
		detectionThrLine= new TLine(lgFluxThrRef,0,lgFluxThrRef,1);
		detectionThrLine->SetLineColor(kBlack);
		detectionThrLine->SetLineStyle(kDashed);
		detectionThrLines.push_back(detectionThrLine);
	}

	//===============================================
	//==          DRAW COMPLETENESS
	//===============================================	
	//- Efficiency plot compact source
	TCanvas* EffPlot= new TCanvas("EffPlot","EffPlot");
	EffPlot->cd();

	TH2D* EffPlotBkg= new TH2D("EffPlotBkg","",100,xMin_draw,xMax_draw,100,0,1.1);
	EffPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
	EffPlotBkg->GetYaxis()->SetTitle("Completeness");
	EffPlotBkg->Draw();

	double Eff_detectionAreaX[]= {xMin_draw,log10(5*opt.noiseLevel_true)};
	double Eff_detectionAreaY[]= {0,0};
	TGraph* Eff_sigmaDetectionArea = new TGraph(2,Eff_detectionAreaX,Eff_detectionAreaY);
  Eff_sigmaDetectionArea->SetLineColor(kGray+2);
  Eff_sigmaDetectionArea->SetLineWidth(15001);
  Eff_sigmaDetectionArea->SetFillStyle(3005);	
	Eff_sigmaDetectionArea->SetFillColor(kGray);
	Eff_sigmaDetectionArea->Draw("C");

	TGraphAsymmErrors* EfficiencyGraph_compact_fit= Efficiency_compact_fit->CreateGraph(); 
	TGraphAsymmErrors* EfficiencyGraph_compact_fit_presel= Efficiency_compact_fit_presel->CreateGraph(); 
	TGraphAsymmErrors* EfficiencyGraph_compact_fit_cutsel= Efficiency_compact_fit_cutsel->CreateGraph(); 
	TGraphAsymmErrors* EfficiencyGraph_compact_fit_nnsel= Efficiency_compact_fit_nnsel->CreateGraph(); 
	if(!opt.drawErrorX){
		for(int i=0;i<EfficiencyGraph_compact_fit->GetN();i++){
			EfficiencyGraph_compact_fit->SetPointEXhigh(i,0);
			EfficiencyGraph_compact_fit->SetPointEXlow(i,0);		
		}//end loop graph points
		for(int i=0;i<EfficiencyGraph_compact_fit_presel->GetN();i++){
			EfficiencyGraph_compact_fit_presel->SetPointEXhigh(i,0);
			EfficiencyGraph_compact_fit_presel->SetPointEXlow(i,0);		
		}//end loop graph points
		for(int i=0;i<EfficiencyGraph_compact_fit_cutsel->GetN();i++){
			EfficiencyGraph_compact_fit_cutsel->SetPointEXhigh(i,0);
			EfficiencyGraph_compact_fit_cutsel->SetPointEXlow(i,0);		
		}//end loop graph points
		for(int i=0;i<EfficiencyGraph_compact_fit_nnsel->GetN();i++){
			EfficiencyGraph_compact_fit_nnsel->SetPointEXhigh(i,0);
			EfficiencyGraph_compact_fit_nnsel->SetPointEXlow(i,0);		
		}//end loop graph points
	}//close if


	EfficiencyGraph_compact_fit->SetMarkerStyle(8);
	EfficiencyGraph_compact_fit->SetMarkerSize(1.3);
	EfficiencyGraph_compact_fit->SetMarkerColor(kBlack);
	EfficiencyGraph_compact_fit->SetLineColor(kBlack);
	EfficiencyGraph_compact_fit->Draw("PLZ same");
	
	EfficiencyGraph_compact_fit_presel->SetMarkerStyle(21);
	EfficiencyGraph_compact_fit_presel->SetMarkerSize(1.3);
	EfficiencyGraph_compact_fit_presel->SetMarkerColor(kRed);
	EfficiencyGraph_compact_fit_presel->SetLineColor(kRed);
	EfficiencyGraph_compact_fit_presel->Draw("PLZ same");

	EfficiencyGraph_compact_fit_cutsel->SetMarkerStyle(33);
	EfficiencyGraph_compact_fit_cutsel->SetMarkerSize(1.7);
	EfficiencyGraph_compact_fit_cutsel->SetMarkerColor(kBlue);
	EfficiencyGraph_compact_fit_cutsel->SetLineColor(kBlue);
	EfficiencyGraph_compact_fit_cutsel->Draw("PLZ same");

	EfficiencyGraph_compact_fit_nnsel->SetMarkerStyle(22);
	EfficiencyGraph_compact_fit_nnsel->SetMarkerSize(1.7);
	EfficiencyGraph_compact_fit_nnsel->SetMarkerColor(kGreen+1);
	EfficiencyGraph_compact_fit_nnsel->SetLineColor(kGreen+1);
	EfficiencyGraph_compact_fit_nnsel->Draw("PLZ same");

	TLegend* EffPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	EffPlotLegend->SetFillColor(0);
	EffPlotLegend->SetTextSize(0.045);
	EffPlotLegend->SetTextFont(52);
	EffPlotLegend->AddEntry(EfficiencyGraph_compact_fit,"fit","PL");
	EffPlotLegend->AddEntry(EfficiencyGraph_compact_fit_presel,"presel","PL");	
	EffPlotLegend->AddEntry(EfficiencyGraph_compact_fit_cutsel,"cut sel","PL");
	EffPlotLegend->AddEntry(EfficiencyGraph_compact_fit_nnsel,"NN sel","PL");
	EffPlotLegend->AddEntry(Eff_sigmaDetectionArea,"5#sigma detection limit","F");
	EffPlotLegend->Draw("same");


	//- Efficiency plot extended source
	TCanvas* EffPlot_ext= new TCanvas("EffPlot_ext","EffPlot_ext");
	EffPlot_ext->cd();

	TH2D* EffPlotBkg_ext= new TH2D("EffPlotBkg_ext","",100,xMin_draw,xMax_draw,100,0,1.2);
	EffPlotBkg_ext->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
	EffPlotBkg_ext->GetYaxis()->SetTitle("Completeness");
	EffPlotBkg_ext->Draw();

	double Eff_detectionAreaX_ext[]= {xMin_draw,log10(5*opt.noiseLevel_true)};
	double Eff_detectionAreaY_ext[]= {0,0};

	TGraph* Eff_sigmaDetectionArea_ext = new TGraph(2,Eff_detectionAreaX,Eff_detectionAreaY);
  Eff_sigmaDetectionArea_ext->SetLineColor(kGray+2);
  Eff_sigmaDetectionArea_ext->SetLineWidth(15001);
  Eff_sigmaDetectionArea_ext->SetFillStyle(3005);	
	Eff_sigmaDetectionArea_ext->SetFillColor(kGray);
	Eff_sigmaDetectionArea_ext->Draw("C");

	TGraphAsymmErrors* EfficiencyGraph_ext= Efficiency_ext->CreateGraph(); 
	std::vector<TGraphAsymmErrors*> EfficiencyGraph_simtypes_ext;
	for(int k=0;k<nSimTypes;k++){
		TGraphAsymmErrors* EfficiencyGraph_simtype= Efficiency_simtypes_ext[k]->CreateGraph(); 
		EfficiencyGraph_simtypes_ext.push_back(EfficiencyGraph_simtype);
	}

	if(!opt.drawErrorX){
		for(int i=0;i<EfficiencyGraph_ext->GetN();i++){
			EfficiencyGraph_ext->SetPointEXhigh(i,0);
			EfficiencyGraph_ext->SetPointEXlow(i,0);		
		}//end loop graph points

		for(int k=0;k<nSimTypes;k++){	
			for(int i=0;i<EfficiencyGraph_simtypes_ext[k]->GetN();i++){
				EfficiencyGraph_simtypes_ext[k]->SetPointEXhigh(i,0);
				EfficiencyGraph_simtypes_ext[k]->SetPointEXlow(i,0);		
			}//end loop graph points
		}//end loop sim types
	}//close if

	EfficiencyGraph_ext->SetMarkerStyle(8);
	EfficiencyGraph_ext->SetMarkerSize(1.3);
	EfficiencyGraph_ext->SetMarkerColor(kBlack);
	EfficiencyGraph_ext->SetLineColor(kBlack);
	EfficiencyGraph_ext->Draw("PLZ same");

	for(int k=0;k<nSimTypes;k++){	
		EfficiencyGraph_simtypes_ext[k]->SetMarkerStyle(8);
		EfficiencyGraph_simtypes_ext[k]->SetMarkerSize(1.3);
		EfficiencyGraph_simtypes_ext[k]->SetMarkerColor(SimTypeColors[k]);
		EfficiencyGraph_simtypes_ext[k]->SetLineColor(SimTypeColors[k]);
		//EfficiencyGraph_simtypes_ext[k]->Draw("PLZ same");
	}
	
	TLegend* EffPlotLegend_ext= new TLegend(0.6,0.7,0.7,0.8);
	EffPlotLegend_ext->SetFillColor(0);
	EffPlotLegend_ext->SetTextSize(0.045);
	EffPlotLegend_ext->SetTextFont(52);
	EffPlotLegend_ext->AddEntry(EfficiencyGraph_ext,"all","PL");
	//for(int k=0;k<nSimTypes;k++) EffPlotLegend_ext->AddEntry(EfficiencyGraph_simtypes_ext[k],SimTypeLabels[k].c_str(),"PL");
	EffPlotLegend_ext->AddEntry(Eff_sigmaDetectionArea_ext,"5#sigma detection limit","F");
	EffPlotLegend_ext->Draw("same");

	
	//===============================================
	//==          DRAW RELIABILITY
	//===============================================
	double minReliabilityX_draw= -4;
	double maxReliabilityX_draw= 0.1;
	double maxReliabilityErr= 0.5;

	//Compact
	TCanvas* ReliabilityPlot= new TCanvas("ReliabilityPlot","ReliabilityPlot");
	ReliabilityPlot->cd();

	TH2D* ReliabilityPlotBkg= new TH2D("ReliabilityPlotBkg","",100,minReliabilityX_draw,maxReliabilityX_draw,100,0,1.1);
	ReliabilityPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{meas}/Jy)");
	ReliabilityPlotBkg->GetYaxis()->SetTitle("Reliability");
	ReliabilityPlotBkg->Draw();

	Eff_sigmaDetectionArea->Draw("C");

	TGraphAsymmErrors* ReliabilityGraph_compact_fit= Reliability_compact_fit->CreateGraph(); 
	TGraphAsymmErrors* ReliabilityGraph_compact_fit_presel= Reliability_compact_fit_presel->CreateGraph(); 
	TGraphAsymmErrors* ReliabilityGraph_compact_fit_cutsel= Reliability_compact_fit_cutsel->CreateGraph(); 
	TGraphAsymmErrors* ReliabilityGraph_compact_fit_nnsel= Reliability_compact_fit_nnsel->CreateGraph();
	if(!opt.drawErrorX){
		int pointCounter= 0;
		int N= ReliabilityGraph_compact_fit->GetN();
		for(int i=0;i<N;i++){
			
			double x, y;
			ReliabilityGraph_compact_fit->GetPoint(pointCounter,x,y);
			cout<<"Bin "<<i<<", pointCounter="<<pointCounter<<", x="<<x<<", Reliability(@fit level)="<<y<<endl;
			if(x<minReliabilityX_draw || x>maxReliabilityX_draw) {
				cout<<"Remove bin "<<i<<", x="<<x<<endl;
				ReliabilityGraph_compact_fit->RemovePoint(pointCounter);
				continue;
			}

			double yerr_low= ReliabilityGraph_compact_fit->GetErrorYlow(pointCounter);
			double yerr_up= ReliabilityGraph_compact_fit->GetErrorYhigh(pointCounter);
			if(yerr_low>maxReliabilityErr){
				ReliabilityGraph_compact_fit->SetPointEYlow(pointCounter,0);
			}			
			if(yerr_up>maxReliabilityErr){
				ReliabilityGraph_compact_fit->SetPointEYhigh(pointCounter,0);
			}
	
			ReliabilityGraph_compact_fit->SetPointEXhigh(pointCounter,0);
			ReliabilityGraph_compact_fit->SetPointEXlow(pointCounter,0);	
			pointCounter++;	

		}//end loop graph points

		pointCounter= 0;
		N= ReliabilityGraph_compact_fit_presel->GetN();
		cout<<"N="<<N<<", ReliabilityGraph_compact_fit_presel->GetN()="<<ReliabilityGraph_compact_fit_presel->GetN()<<endl;
		for(int i=0;i<N;i++){
			double x, y;
			ReliabilityGraph_compact_fit_presel->GetPoint(pointCounter,x,y);
			cout<<"Bin "<<i<<", pointCounter="<<pointCounter<<", x="<<x<<", Reliability(@presel level)="<<y<<endl;
			
			if(x<minReliabilityX_draw || x>maxReliabilityX_draw) {
				ReliabilityGraph_compact_fit_presel->RemovePoint(pointCounter);
				continue;
			}
			double yerr_low= ReliabilityGraph_compact_fit_presel->GetErrorYlow(pointCounter);
			double yerr_up= ReliabilityGraph_compact_fit_presel->GetErrorYhigh(pointCounter);
			if(yerr_low>maxReliabilityErr){
				ReliabilityGraph_compact_fit_presel->SetPointEYlow(pointCounter,0);
			}			
			if(yerr_up>maxReliabilityErr){
				ReliabilityGraph_compact_fit_presel->SetPointEYhigh(pointCounter,0);
			}
			ReliabilityGraph_compact_fit_presel->SetPointEXhigh(pointCounter,0);
			ReliabilityGraph_compact_fit_presel->SetPointEXlow(pointCounter,0);	
			pointCounter++;	
		}//end loop graph points

		pointCounter= 0;
		N= ReliabilityGraph_compact_fit_cutsel->GetN();
		for(int i=0;i<N;i++){
			double x, y;
			ReliabilityGraph_compact_fit_cutsel->GetPoint(pointCounter,x,y);
			if(x<minReliabilityX_draw || x>maxReliabilityX_draw) {
				ReliabilityGraph_compact_fit_cutsel->RemovePoint(pointCounter);
				continue;
			}

			double yerr_low= ReliabilityGraph_compact_fit_cutsel->GetErrorYlow(pointCounter);
			double yerr_up= ReliabilityGraph_compact_fit_cutsel->GetErrorYhigh(pointCounter);
			if(yerr_low>maxReliabilityErr){
				ReliabilityGraph_compact_fit_cutsel->SetPointEYlow(pointCounter,0);
			}			
			if(yerr_up>maxReliabilityErr){
				ReliabilityGraph_compact_fit_cutsel->SetPointEYhigh(pointCounter,0);
			}

			ReliabilityGraph_compact_fit_cutsel->SetPointEXhigh(pointCounter,0);
			ReliabilityGraph_compact_fit_cutsel->SetPointEXlow(pointCounter,0);	
			pointCounter++;	
		}//end loop graph points

		pointCounter= 0;
		N= ReliabilityGraph_compact_fit_nnsel->GetN();
		for(int i=0;i<N;i++){
			double x, y;
			ReliabilityGraph_compact_fit_nnsel->GetPoint(pointCounter,x,y);
			if(x<minReliabilityX_draw || x>maxReliabilityX_draw) {
				ReliabilityGraph_compact_fit_nnsel->RemovePoint(pointCounter);
				continue;
			}
			double yerr_low= ReliabilityGraph_compact_fit_nnsel->GetErrorYlow(pointCounter);
			double yerr_up= ReliabilityGraph_compact_fit_nnsel->GetErrorYhigh(pointCounter);
			if(yerr_low>maxReliabilityErr){
				ReliabilityGraph_compact_fit_nnsel->SetPointEYlow(pointCounter,0);
			}			
			if(yerr_up>maxReliabilityErr){
				ReliabilityGraph_compact_fit_nnsel->SetPointEYhigh(pointCounter,0);
			}
			ReliabilityGraph_compact_fit_nnsel->SetPointEXhigh(pointCounter,0);
			ReliabilityGraph_compact_fit_nnsel->SetPointEXlow(pointCounter,0);	
			pointCounter++;	
		}//end loop graph points
	}//close if


	ReliabilityGraph_compact_fit->SetMarkerStyle(8);
	ReliabilityGraph_compact_fit->SetMarkerSize(1.3);
	ReliabilityGraph_compact_fit->SetMarkerColor(kBlack);
	ReliabilityGraph_compact_fit->SetLineColor(kBlack);
	ReliabilityGraph_compact_fit->Draw("PZL same");

	ReliabilityGraph_compact_fit_presel->SetMarkerStyle(21);
	ReliabilityGraph_compact_fit_presel->SetMarkerSize(1.3);
	ReliabilityGraph_compact_fit_presel->SetMarkerColor(kRed);
	ReliabilityGraph_compact_fit_presel->SetLineColor(kRed);
	ReliabilityGraph_compact_fit_presel->Draw("PZL same");

	ReliabilityGraph_compact_fit_cutsel->SetMarkerStyle(33);
	ReliabilityGraph_compact_fit_cutsel->SetMarkerSize(1.7);
	ReliabilityGraph_compact_fit_cutsel->SetMarkerColor(kBlue);
	ReliabilityGraph_compact_fit_cutsel->SetLineColor(kBlue);
	ReliabilityGraph_compact_fit_cutsel->Draw("PZL same");

	ReliabilityGraph_compact_fit_nnsel->SetMarkerStyle(22);
	ReliabilityGraph_compact_fit_nnsel->SetMarkerSize(1.7);
	ReliabilityGraph_compact_fit_nnsel->SetMarkerColor(kGreen+1);
	ReliabilityGraph_compact_fit_nnsel->SetLineColor(kGreen+1);
	ReliabilityGraph_compact_fit_nnsel->Draw("PZL same");
	
	TLegend* ReliabilityPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	ReliabilityPlotLegend->SetFillColor(0);
	ReliabilityPlotLegend->SetTextSize(0.045);
	ReliabilityPlotLegend->SetTextFont(52);
	ReliabilityPlotLegend->AddEntry(ReliabilityGraph_compact_fit,"fit","PL");
	ReliabilityPlotLegend->AddEntry(ReliabilityGraph_compact_fit_presel,"presel","PL");
	ReliabilityPlotLegend->AddEntry(ReliabilityGraph_compact_fit_cutsel,"cut sel","PL");
	ReliabilityPlotLegend->AddEntry(ReliabilityGraph_compact_fit_nnsel,"NN sel","PL");
	ReliabilityPlotLegend->AddEntry(Eff_sigmaDetectionArea,"5#sigma detection limit","F");
	ReliabilityPlotLegend->Draw("same");


	//Extended
	TCanvas* ReliabilityPlot_ext= new TCanvas("ReliabilityPlot_ext","ReliabilityPlot_ext");
	ReliabilityPlot_ext->cd();

	TH2D* ReliabilityPlotBkg_ext= new TH2D("ReliabilityPlotBkg_ext","",100,xMin_draw,xMax_draw,100,0,1.2);
	ReliabilityPlotBkg_ext->GetXaxis()->SetTitle("log_{10}(S_{rec}/Jy)");
	ReliabilityPlotBkg_ext->GetYaxis()->SetTitle("Reliability");
	ReliabilityPlotBkg_ext->Draw();

	Reliability_ext->SetMarkerStyle(8);
	Reliability_ext->SetMarkerColor(kBlack);
	Reliability_ext->SetLineColor(kBlack);
	Reliability_ext->Draw("p same");
	
	detectionThrLine->Draw("l");

	

	//===============================================
	//==          DRAW POSITIONAL ACCURACY
	//===============================================
	double minPosAccuracy_draw= -4;
	double maxPosAccuracy_draw= 4;
	double minPosAccuracyX_draw= -4;
	double maxPosAccuracyX_draw= 0.4;
	double minPosResolution_draw= 0;
	double maxPosResolution_draw= 5.;
	double PosAccuracy_detectionAreaX[]= {minPosAccuracyX_draw,log10(5*opt.noiseLevel_true)};
	double PosAccuracy_detectionAreaY[]= {minPosAccuracy_draw,minPosAccuracy_draw};
	double PosResolution_detectionAreaX[]= {minPosAccuracyX_draw,log10(5*opt.noiseLevel_true)};
	double PosResolution_detectionAreaY[]= {minPosResolution_draw,minPosResolution_draw};

	
	TCanvas* PosAccuracyPlot= new TCanvas("PosAccuracyPlot","PosAccuracyPlot",550,950);
	PosAccuracyPlot->cd();

	//- Draw pos bias
	TPad* PosAccuracyPlotPad1= new TPad("PosAccuracyPlotPad1","PosAccuracyPlotPad1",0,0.5,1,1);
	PosAccuracyPlotPad1->SetBottomMargin(0);
	PosAccuracyPlotPad1->SetLeftMargin(0.15);
  PosAccuracyPlotPad1->SetRightMargin(0.05);
	PosAccuracyPlotPad1->cd();

	TH2D* PosAccuracyPlotBkg= new TH2D("PosAccuracyPlotBkg","",100,minPosAccuracyX_draw,maxPosAccuracyX_draw,100,minPosAccuracy_draw,maxPosAccuracy_draw);
	PosAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");	
	PosAccuracyPlotBkg->GetXaxis()->SetTitleSize(0.08);
	PosAccuracyPlotBkg->GetXaxis()->SetTitleOffset(0.85);
	PosAccuracyPlotBkg->GetXaxis()->SetLabelSize(0.06);
	PosAccuracyPlotBkg->GetYaxis()->SetTitle("<RA>, <Dec> ('')");
	PosAccuracyPlotBkg->GetYaxis()->SetTitleSize(0.08);
	PosAccuracyPlotBkg->GetYaxis()->SetTitleOffset(0.8);
	PosAccuracyPlotBkg->GetYaxis()->SetLabelSize(0.06);
	PosAccuracyPlotBkg->Draw();

	TGraph* PosAccuracy_sigmaDetectionArea = new TGraph(2,PosAccuracy_detectionAreaX,PosAccuracy_detectionAreaY);
  PosAccuracy_sigmaDetectionArea->SetLineColor(kGray+2);
  PosAccuracy_sigmaDetectionArea->SetLineWidth(15001);
  PosAccuracy_sigmaDetectionArea->SetFillStyle(3005);	
	PosAccuracy_sigmaDetectionArea->SetFillColor(kGray);
	PosAccuracy_sigmaDetectionArea->Draw("C");

	xPosAccuracyGraph_compact_fit->SetMarkerStyle(8);
	xPosAccuracyGraph_compact_fit->SetMarkerSize(1.3);
	xPosAccuracyGraph_compact_fit->SetMarkerColor(kBlack);
	xPosAccuracyGraph_compact_fit->SetLineColor(kBlack);
	//xPosAccuracyGraph_compact_fit->Draw("PZ same");

	xPosAccuracyGraph_compact_fit_presel->SetMarkerStyle(8);
	xPosAccuracyGraph_compact_fit_presel->SetMarkerSize(1.3);
	xPosAccuracyGraph_compact_fit_presel->SetMarkerColor(kBlack);
	xPosAccuracyGraph_compact_fit_presel->SetLineColor(kBlack);
	xPosAccuracyGraph_compact_fit_presel->Draw("PZ same");
	
	xPosAccuracyGraph_compact_fit_cutsel->SetMarkerStyle(8);
	xPosAccuracyGraph_compact_fit_cutsel->SetMarkerSize(1.3);
	xPosAccuracyGraph_compact_fit_cutsel->SetMarkerColor(kBlack);
	xPosAccuracyGraph_compact_fit_cutsel->SetLineColor(kBlack);
	//xPosAccuracyGraph_compact_fit_cutsel->Draw("PZ same");
	
	xPosAccuracyGraph_compact_fit_nnsel->SetMarkerStyle(8);
	xPosAccuracyGraph_compact_fit_nnsel->SetMarkerSize(1.3);
	xPosAccuracyGraph_compact_fit_nnsel->SetMarkerColor(kBlack);
	xPosAccuracyGraph_compact_fit_nnsel->SetLineColor(kBlack);
	//xPosAccuracyGraph_compact_fit_nnsel->Draw("PZ same");

	yPosAccuracyGraph_compact_fit->SetMarkerStyle(21);
	yPosAccuracyGraph_compact_fit->SetMarkerSize(1.3);
	yPosAccuracyGraph_compact_fit->SetMarkerColor(kRed);
	yPosAccuracyGraph_compact_fit->SetLineColor(kRed);
	//yPosAccuracyGraph_compact_fit->Draw("PZ same");

	yPosAccuracyGraph_compact_fit_presel->SetMarkerStyle(21);
	yPosAccuracyGraph_compact_fit_presel->SetMarkerSize(1.3);
	yPosAccuracyGraph_compact_fit_presel->SetMarkerColor(kRed);
	yPosAccuracyGraph_compact_fit_presel->SetLineColor(kRed);
	yPosAccuracyGraph_compact_fit_presel->Draw("PZ same");

	yPosAccuracyGraph_compact_fit_cutsel->SetMarkerStyle(21);
	yPosAccuracyGraph_compact_fit_cutsel->SetMarkerSize(1.3);
	yPosAccuracyGraph_compact_fit_cutsel->SetMarkerColor(kRed);
	yPosAccuracyGraph_compact_fit_cutsel->SetLineColor(kRed);
	//yPosAccuracyGraph_compact_fit_cutsel->Draw("ep same");

	yPosAccuracyGraph_compact_fit_nnsel->SetMarkerStyle(21);
	yPosAccuracyGraph_compact_fit_nnsel->SetMarkerSize(1.3);
	yPosAccuracyGraph_compact_fit_nnsel->SetMarkerColor(kRed);
	yPosAccuracyGraph_compact_fit_nnsel->SetLineColor(kRed);
	//yPosAccuracyGraph_compact_fit_nnsel->Draw("ep same");
	
	
	TLine* refLine_posAccuracy= new TLine(minPosAccuracyX_draw,0,maxPosAccuracyX_draw,0);
	refLine_posAccuracy->SetLineColor(kBlack);
	refLine_posAccuracy->SetLineStyle(kDashed);
	refLine_posAccuracy->Draw("same");


	TLegend* PosAccuracyPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	PosAccuracyPlotLegend->SetFillColor(0);
	PosAccuracyPlotLegend->SetTextSize(0.05);
	PosAccuracyPlotLegend->SetTextFont(52);
	PosAccuracyPlotLegend->AddEntry(xPosAccuracyGraph_compact_fit,"<RA>","PL");
	PosAccuracyPlotLegend->AddEntry(yPosAccuracyGraph_compact_fit,"<Dec>","PL");
	PosAccuracyPlotLegend->AddEntry(PosAccuracy_sigmaDetectionArea,"5#sigma detection limit","F");
	PosAccuracyPlotLegend->Draw("same");

	//- Draw pos resolution
	TPad* PosAccuracyPlotPad2= new TPad("PosAccuracyPlotPad2","PosAccuracyPlotPad2",0,0.01,1,0.5);
	PosAccuracyPlotPad2->SetTopMargin(0);
	PosAccuracyPlotPad2->SetLeftMargin(0.15);
	PosAccuracyPlotPad2->SetBottomMargin(0.12);
  PosAccuracyPlotPad2->SetRightMargin(0.05);
	PosAccuracyPlotPad2->cd();

	TH2D* PosAccuracyPlotBkg2= new TH2D("PosAccuracyPlotBkg2","",100,minPosAccuracyX_draw,maxPosAccuracyX_draw,100,minPosResolution_draw,maxPosResolution_draw);
	PosAccuracyPlotBkg2->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
	PosAccuracyPlotBkg2->GetXaxis()->SetTitleSize(0.08);
	PosAccuracyPlotBkg2->GetXaxis()->SetTitleOffset(0.85);
	PosAccuracyPlotBkg2->GetXaxis()->SetLabelSize(0.06);
	PosAccuracyPlotBkg2->GetYaxis()->SetTitle("#sigma(RA), #sigma(Dec)  ('')");
	PosAccuracyPlotBkg2->GetYaxis()->SetTitleSize(0.08);
	PosAccuracyPlotBkg2->GetYaxis()->SetTitleOffset(0.8);
	PosAccuracyPlotBkg2->GetYaxis()->SetLabelSize(0.06);
	PosAccuracyPlotBkg2->Draw();

	TGraph* PosResolution_sigmaDetectionArea = new TGraph(2,PosResolution_detectionAreaX,PosResolution_detectionAreaY);
  PosResolution_sigmaDetectionArea->SetLineColor(kGray+2);
  PosResolution_sigmaDetectionArea->SetLineWidth(15001);
  PosResolution_sigmaDetectionArea->SetFillStyle(3005);	
	PosResolution_sigmaDetectionArea->SetFillColor(kGray);
	PosResolution_sigmaDetectionArea->Draw("C");

	xPosResolutionGraph_compact_fit->SetMarkerStyle(8);
	xPosResolutionGraph_compact_fit->SetMarkerSize(1.3);
	xPosResolutionGraph_compact_fit->SetMarkerColor(kBlack);
	xPosResolutionGraph_compact_fit->SetLineColor(kBlack);
	xPosResolutionGraph_compact_fit->Draw("PZ same");

	yPosResolutionGraph_compact_fit->SetMarkerStyle(21);
	yPosResolutionGraph_compact_fit->SetMarkerSize(1.3);
	yPosResolutionGraph_compact_fit->SetMarkerColor(kRed);
	yPosResolutionGraph_compact_fit->SetLineColor(kRed);
	yPosResolutionGraph_compact_fit->Draw("PZ same");

	double sigmaToQ1= ROOT::Math::gaussian_quantile(0.75,1);
	double sigmaToIQR= 2*ROOT::Math::gaussian_quantile(0.75,1);
	double sigmaToHalfIQR= ROOT::Math::gaussian_quantile(0.75,1);
	
	TF1* expPosXSigmaFcn_minus= new TF1("expPosXSigmaFcn_minus",FitPosErrFcn,xMin_draw,xMax_draw,5);
	expPosXSigmaFcn_minus->SetNpx(50);
	expPosXSigmaFcn_minus->SetParameters(-1,opt.Bmaj,opt.Bmin,opt.pixSize,opt.noiseLevel_true);
	expPosXSigmaFcn_minus->SetLineColor(kBlack);
	expPosXSigmaFcn_minus->SetLineStyle(9);
	//expPosXSigmaFcn_minus->Draw("C same");

	TF1* expPosXSigmaFcn_plus= new TF1("expPosXSigmaFcn_plus",FitPosErrFcn,xMin_draw,xMax_draw,5);
	expPosXSigmaFcn_plus->SetNpx(50);
	expPosXSigmaFcn_plus->SetParameters(sigmaToHalfIQR,opt.Bmaj,opt.Bmin,opt.pixSize,opt.noiseLevel_true);
	expPosXSigmaFcn_plus->SetLineColor(kBlack);
	expPosXSigmaFcn_plus->SetLineStyle(9);
	expPosXSigmaFcn_plus->Draw("C same");

	TF1* expPosYSigmaFcn_minus= new TF1("expPosYSigmaFcn_minus",FitPosErrFcn,xMin_draw,xMax_draw,5);
	expPosYSigmaFcn_minus->SetNpx(50);
	expPosYSigmaFcn_minus->SetParameters(-1,opt.Bmin,opt.Bmaj,opt.pixSize,opt.noiseLevel_true);
	expPosYSigmaFcn_minus->SetLineColor(kRed);
	expPosYSigmaFcn_minus->SetLineStyle(kDotted);
	//expPosYSigmaFcn_minus->Draw("C same");

	TF1* expPosYSigmaFcn_plus= new TF1("expPosYSigmaFcn_plus",FitPosErrFcn,xMin_draw,xMax_draw,5);
	expPosYSigmaFcn_plus->SetNpx(50);
	expPosYSigmaFcn_plus->SetParameters(sigmaToHalfIQR,opt.Bmin,opt.Bmaj,opt.pixSize,opt.noiseLevel_true);
	expPosYSigmaFcn_plus->SetLineColor(kRed);
	expPosYSigmaFcn_plus->SetLineStyle(kDotted);
	expPosYSigmaFcn_plus->Draw("C same");


	TLegend* PosAccuracyPlotLegend2= new TLegend(0.5,0.5,0.7,0.7);
	PosAccuracyPlotLegend2->SetFillColor(0);
	PosAccuracyPlotLegend2->SetTextSize(0.05);
	PosAccuracyPlotLegend2->SetTextFont(52);
	PosAccuracyPlotLegend2->AddEntry(xPosResolutionGraph_compact_fit,"#sigma_{RA}","PL");
	PosAccuracyPlotLegend2->AddEntry(yPosResolutionGraph_compact_fit,"#sigma_{Dec}","PL");
	PosAccuracyPlotLegend2->AddEntry(PosResolution_sigmaDetectionArea,"5#sigma detection limit","F");
	PosAccuracyPlotLegend2->AddEntry(expPosXSigmaFcn_plus,"#sigma_{RA}^{ideal}","L");
	PosAccuracyPlotLegend2->AddEntry(expPosYSigmaFcn_plus,"#sigma_{Dec}^{ideal}","L");
	PosAccuracyPlotLegend2->Draw("same");

	TPad* clearLabelPad = new TPad("clearLabelPad", "clearLabelPad",0.05201342,0.479021,0.08557047,0.5314685);
  clearLabelPad->SetBorderMode(0);

	PosAccuracyPlot->cd();
	PosAccuracyPlotPad1->Draw();
	PosAccuracyPlotPad2->Draw();
	clearLabelPad->Draw();
	PosAccuracyPlot->Update();

	
	//- Pos accuracy for extended sources
	TCanvas* PosAccuracyPlot_ext= new TCanvas("PosAccuracyPlot_ext","PosAccuracyPlot_ext");
	PosAccuracyPlot_ext->cd();

	TH2D* PosAccuracyPlotBkg_ext= new TH2D("PosAccuracyPlotBkg_ext","",100,xMin_draw,xMax_draw,100,-20,20);
	PosAccuracyPlotBkg_ext->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	PosAccuracyPlotBkg_ext->GetYaxis()->SetTitle("<#deltaRA>, <#deltaDec> ('')");
	PosAccuracyPlotBkg_ext->Draw();

	xPosAccuracyGraph_ext->SetMarkerStyle(8);
	xPosAccuracyGraph_ext->SetMarkerColor(kBlack);
	xPosAccuracyGraph_ext->SetLineColor(kBlack);
	xPosAccuracyGraph_ext->Draw("ep same");
	
	yPosAccuracyGraph_ext->SetMarkerStyle(21);
	yPosAccuracyGraph_ext->SetMarkerColor(kRed);
	yPosAccuracyGraph_ext->SetLineColor(kRed);
	yPosAccuracyGraph_ext->Draw("ep same");
	
	TLine* refLine_posAccuracy_ext= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_posAccuracy_ext->SetLineColor(kBlack);
	refLine_posAccuracy_ext->SetLineStyle(kDashed);
	refLine_posAccuracy_ext->Draw("same");

	TLegend* PosAccuracyPlotLegend_ext= new TLegend(0.6,0.7,0.7,0.8);
	PosAccuracyPlotLegend_ext->SetFillColor(0);
	PosAccuracyPlotLegend_ext->SetTextSize(0.045);
	PosAccuracyPlotLegend_ext->SetTextFont(52);
	PosAccuracyPlotLegend_ext->AddEntry(xPosAccuracyGraph_ext,"<#deltaRA>","PL");
	PosAccuracyPlotLegend_ext->AddEntry(yPosAccuracyGraph_ext,"<#deltaDec>","PL");
	PosAccuracyPlotLegend_ext->Draw("same");

	
	//===============================================
	//==          DRAW FLUX ACCURACY
	//===============================================
	//- Flux accuracy
	double minFluxAccuracy_draw= -0.5;
	double maxFluxAccuracy_draw= 1.;
	double minFluxResolution_draw= -0.2;
	double maxFluxResolution_draw= 1.;
	double minFluxAccuracyX_draw= -4;
	double maxFluxAccuracyX_draw= 0.5;
	

	double detectionAreaX[]= {minFluxAccuracyX_draw,log10(5*opt.noiseLevel_true)};
	double detectionAreaY[]= {minFluxAccuracy_draw,minFluxAccuracy_draw};
	double detectionAreaX_reso[]= {minFluxAccuracyX_draw,log10(5*opt.noiseLevel_true)};
	double detectionAreaY_reso[]= {minFluxResolution_draw,minFluxResolution_draw};

	TCanvas* FluxAccuracyPlot= new TCanvas("FluxAccuracyPlot","FluxAccuracyPlot",550,950);
	FluxAccuracyPlot->cd();

	//Flux bias
	TPad* FluxAccuracyPlotPad1= new TPad("FluxAccuracyPlotPad1","FluxAccuracyPlotPad1",0,0.5,1,1);
	FluxAccuracyPlotPad1->SetBottomMargin(0);
	FluxAccuracyPlotPad1->SetLeftMargin(0.15);
  FluxAccuracyPlotPad1->SetRightMargin(0.05);
	FluxAccuracyPlotPad1->cd();

	TH2D* FluxAccuracyPlotBkg= new TH2D("FluxAccuracyPlotBkg","",100,minFluxAccuracyX_draw,maxFluxAccuracyX_draw,100,minFluxAccuracy_draw,maxFluxAccuracy_draw);
	FluxAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
	FluxAccuracyPlotBkg->GetYaxis()->SetTitle("<#frac{S_{meas}-S_{gen}}{S_{gen}}>");
	FluxAccuracyPlotBkg->GetXaxis()->SetTitleSize(0.08);
	FluxAccuracyPlotBkg->GetXaxis()->SetTitleOffset(0.85);
	FluxAccuracyPlotBkg->GetXaxis()->SetLabelSize(0.06);
	FluxAccuracyPlotBkg->GetYaxis()->SetTitleSize(0.08);
	FluxAccuracyPlotBkg->GetYaxis()->SetTitleOffset(0.8);
	FluxAccuracyPlotBkg->GetYaxis()->SetLabelSize(0.06);
	FluxAccuracyPlotBkg->Draw();
	
	
	TGraph* sigmaDetectionArea = new TGraph(2,detectionAreaX,detectionAreaY);
  sigmaDetectionArea->SetLineColor(kGray+2);
  sigmaDetectionArea->SetLineWidth(15001);
  sigmaDetectionArea->SetFillStyle(3005);	
	sigmaDetectionArea->SetFillColor(kGray);
	sigmaDetectionArea->Draw("C");

	FluxDensityAccuracyGraphPoints_compact->SetMarkerColor(kGray);
	FluxDensityAccuracyGraphPoints_compact->SetMarkerStyle(8);
	FluxDensityAccuracyGraphPoints_compact->SetMarkerSize(0.3);
	//FluxDensityAccuracyGraphPoints_compact->Draw("P");

	FluxDensityAccuracyGraph_compact_fit->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_compact_fit->SetMarkerSize(1.3);
	FluxDensityAccuracyGraph_compact_fit->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_compact_fit->SetLineColor(kBlack);
	//FluxDensityAccuracyGraph_compact_fit->Draw("EPZ same");

	FluxDensityAccuracyGraph_compact_fit_presel->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_compact_fit_presel->SetMarkerSize(1.3);
	FluxDensityAccuracyGraph_compact_fit_presel->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_compact_fit_presel->SetLineColor(kBlack);
	FluxDensityAccuracyGraph_compact_fit_presel->Draw("EPZ same");

	FluxDensityAccuracyGraph_compact_fit_cutsel->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_compact_fit_cutsel->SetMarkerSize(1.3);
	FluxDensityAccuracyGraph_compact_fit_cutsel->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_compact_fit_cutsel->SetLineColor(kBlack);
	//FluxDensityAccuracyGraph_compact_fit_cutsel->Draw("EPZ same");

	FluxDensityAccuracyGraph_compact_fit_nnsel->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_compact_fit_nnsel->SetMarkerSize(1.3);
	FluxDensityAccuracyGraph_compact_fit_nnsel->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_compact_fit_nnsel->SetLineColor(kBlack);
	//FluxDensityAccuracyGraph_compact_fit_nnsel->Draw("EPZ same");

	TLine* refLine_fluxAccuracy= new TLine(minFluxAccuracyX_draw,0,maxFluxAccuracyX_draw,0);
	refLine_fluxAccuracy->SetLineColor(kBlack);
	refLine_fluxAccuracy->SetLineStyle(kDashed);
	refLine_fluxAccuracy->Draw("same");

	TLegend* FluxAccuracyPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	FluxAccuracyPlotLegend->SetFillColor(0);
	FluxAccuracyPlotLegend->SetTextSize(0.045);
	FluxAccuracyPlotLegend->SetTextFont(52);
	FluxAccuracyPlotLegend->AddEntry(sigmaDetectionArea,"5#sigma detection limit","F");
	FluxAccuracyPlotLegend->Draw("same");

	//Flux resolution
	TPad* FluxAccuracyPlotPad2= new TPad("FluxAccuracyPlotPad2","FluxAccuracyPlotPad2",0,0.01,1,0.5);
	FluxAccuracyPlotPad2->SetTopMargin(0);
	FluxAccuracyPlotPad2->SetLeftMargin(0.15);
	FluxAccuracyPlotPad2->SetBottomMargin(0.12);
  FluxAccuracyPlotPad2->SetRightMargin(0.05);
	FluxAccuracyPlotPad2->cd();

	TH2D* FluxAccuracyPlotBkg2= new TH2D("FluxAccuracyPlotBkg2","",100,minFluxAccuracyX_draw,maxFluxAccuracyX_draw,100,minFluxResolution_draw,maxFluxResolution_draw);
	FluxAccuracyPlotBkg2->GetXaxis()->SetTitle("log_{10}(S_{gen}/Jy)");
	FluxAccuracyPlotBkg2->GetYaxis()->SetTitle("#sigma(#frac{S_{meas}-S_{gen}}{S_{gen}})");
	FluxAccuracyPlotBkg2->GetXaxis()->SetTitleSize(0.08);
	FluxAccuracyPlotBkg2->GetXaxis()->SetTitleOffset(0.85);
	FluxAccuracyPlotBkg2->GetXaxis()->SetLabelSize(0.06);
	FluxAccuracyPlotBkg2->GetYaxis()->SetTitleSize(0.08);
	FluxAccuracyPlotBkg2->GetYaxis()->SetTitleOffset(0.8);
	FluxAccuracyPlotBkg2->GetYaxis()->SetLabelSize(0.06);
	FluxAccuracyPlotBkg2->Draw();

	TGraph* sigmaDetectionArea2 = new TGraph(2,detectionAreaX_reso,detectionAreaY_reso);
  sigmaDetectionArea2->SetLineColor(kGray+2);
  sigmaDetectionArea2->SetLineWidth(15001);
  sigmaDetectionArea2->SetFillStyle(3005);	
	sigmaDetectionArea2->SetFillColor(kGray);
	sigmaDetectionArea2->Draw("C");

	FluxDensityResolutionGraph_compact_fit->SetMarkerStyle(8);
	FluxDensityResolutionGraph_compact_fit->SetMarkerSize(1.3);
	FluxDensityResolutionGraph_compact_fit->SetMarkerColor(kBlack);
	FluxDensityResolutionGraph_compact_fit->SetLineColor(kBlack);
	//FluxDensityResolutionGraph_compact_fit->Draw("EPZ same");

	FluxDensityResolutionGraph_compact_fit_presel->SetMarkerStyle(8);
	FluxDensityResolutionGraph_compact_fit_presel->SetMarkerSize(1.3);
	FluxDensityResolutionGraph_compact_fit_presel->SetMarkerColor(kBlack);
	FluxDensityResolutionGraph_compact_fit_presel->SetLineColor(kBlack);
	FluxDensityResolutionGraph_compact_fit_presel->Draw("EPZ same");

	FluxDensityResolutionGraph_compact_fit_cutsel->SetMarkerStyle(8);
	FluxDensityResolutionGraph_compact_fit_cutsel->SetMarkerSize(1.3);
	FluxDensityResolutionGraph_compact_fit_cutsel->SetMarkerColor(kBlack);
	FluxDensityResolutionGraph_compact_fit_cutsel->SetLineColor(kBlack);
	//FluxDensityResolutionGraph_compact_fit_cutsel->Draw("EPZ same");

	FluxDensityResolutionGraph_compact_fit_nnsel->SetMarkerStyle(8);
	FluxDensityResolutionGraph_compact_fit_nnsel->SetMarkerSize(1.3);
	FluxDensityResolutionGraph_compact_fit_nnsel->SetMarkerColor(kBlack);
	FluxDensityResolutionGraph_compact_fit_nnsel->SetLineColor(kBlack);
	//FluxDensityResolutionGraph_compact_fit_nnsel->Draw("EPZ same");

	TF1* fluxDensityErrLine_1sigma_plus= new TF1("fluxDensityErrLine_1sigma_plus","[0]*[1]/pow(10,x)",minFluxAccuracyX_draw,maxFluxAccuracyX_draw);
	fluxDensityErrLine_1sigma_plus->SetParameters(opt.noiseLevel_true,1);
	fluxDensityErrLine_1sigma_plus->SetLineColor(kBlack);
	fluxDensityErrLine_1sigma_plus->SetLineStyle(kDashed);
	fluxDensityErrLine_1sigma_plus->SetLineWidth(2);
	fluxDensityErrLine_1sigma_plus->Draw("l same");
	
	TF1* fluxDensityErrLine_1sigma_minus= new TF1("fluxDensityErrLine_1sigma_minus","[0]*[1]/pow(10,x)",minFluxAccuracyX_draw,maxFluxAccuracyX_draw);
	fluxDensityErrLine_1sigma_minus->SetParameters(opt.noiseLevel_true,-1);
	fluxDensityErrLine_1sigma_minus->SetLineColor(kBlack);
	fluxDensityErrLine_1sigma_minus->SetLineStyle(kDashed);
	fluxDensityErrLine_1sigma_minus->SetLineWidth(2);
	//fluxDensityErrLine_1sigma_minus->Draw("l same");

	TF1* fluxDensityErrLine_3sigma_plus= new TF1("fluxDensityErrLine_3sigma_plus","[0]*[1]/pow(10,x)",minFluxAccuracyX_draw,maxFluxAccuracyX_draw);
	fluxDensityErrLine_3sigma_plus->SetParameters(opt.noiseLevel_true,3);
	fluxDensityErrLine_3sigma_plus->SetLineColor(kBlack);
	fluxDensityErrLine_3sigma_plus->SetLineStyle(kDotted);
	fluxDensityErrLine_3sigma_plus->SetLineWidth(2);
	fluxDensityErrLine_3sigma_plus->Draw("l same");
	
	TF1* fluxDensityErrLine_3sigma_minus= new TF1("fluxDensityErrLine_3sigma_minus","[0]*[1]/pow(10,x)",minFluxAccuracyX_draw,maxFluxAccuracyX_draw);
	fluxDensityErrLine_3sigma_minus->SetParameters(opt.noiseLevel_true,-3);
	fluxDensityErrLine_3sigma_minus->SetLineColor(kBlack);
	fluxDensityErrLine_3sigma_minus->SetLineStyle(kDotted);
	fluxDensityErrLine_3sigma_minus->SetLineWidth(2);
	//fluxDensityErrLine_3sigma_minus->Draw("l same");

	TLegend* FluxAccuracyPlotLegend2= new TLegend(0.6,0.7,0.7,0.8);
	FluxAccuracyPlotLegend2->SetFillColor(0);
	FluxAccuracyPlotLegend2->SetTextSize(0.045);
	FluxAccuracyPlotLegend2->SetTextFont(52);
	FluxAccuracyPlotLegend2->AddEntry(fluxDensityErrLine_1sigma_plus,"1#sigma","L");
	FluxAccuracyPlotLegend2->AddEntry(fluxDensityErrLine_3sigma_plus,"3#sigma","L");
	FluxAccuracyPlotLegend2->AddEntry(sigmaDetectionArea,"5#sigma detection limit","F");
	FluxAccuracyPlotLegend2->Draw("same");

	TPad* clearLabelPad2 = new TPad("clearLabelPad2", "clearLabelPad",0.05201342,0.479021,0.08557047,0.5314685);
  clearLabelPad2->SetBorderMode(0);

	FluxAccuracyPlot->cd();
	FluxAccuracyPlotPad1->Draw();
	FluxAccuracyPlotPad2->Draw();
	clearLabelPad2->Draw();
	FluxAccuracyPlot->Update();


	//Flux histos
	TCanvas* FluxPullHistoPlot= 0;
	
	for(size_t i=0;i<FluxDensityPullHistos_compact_fit_presel.size();i++){
		TString canvasName= Form("FluxPullHistoPlot_bin%d",(int)(i+1));
		FluxPullHistoPlot= new TCanvas(canvasName,canvasName);
		FluxPullHistoPlot->cd();

		FluxDensityPullHistos_compact_fit_presel[i]->Draw("hist");
	}



	//- FLux accuracy for extended sources 
	TCanvas* FluxAccuracyPlot_ext= new TCanvas("FluxAccuracyPlot_ext","FluxAccuracyPlot_ext");
	FluxAccuracyPlot_ext->cd();

	TH2D* FluxAccuracyPlotBkg_ext= new TH2D("FluxAccuracyPlotBkg_ext","",100,xMin_draw,xMax_draw,100,-5,5);
	FluxAccuracyPlotBkg_ext->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	FluxAccuracyPlotBkg_ext->GetYaxis()->SetTitle("<S_{rec}/S_{true}-1>");
	FluxAccuracyPlotBkg_ext->Draw();

	FluxDensityAccuracyGraph_ext->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_ext->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_ext->SetLineColor(kBlack);
	FluxDensityAccuracyGraph_ext->Draw("ep same");

	/*
	TLine* refLine_fluxDensityAccuracy_ext= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_fluxDensityAccuracy_ext->SetLineColor(kBlack);
	refLine_fluxDensityAccuracy_ext->SetLineStyle(kDashed);
	refLine_fluxDensityAccuracy_ext->Draw("same");

	TLine* detectionThrLine_fluxDensityAccuracy_ext= new TLine(SNThr,-5,SNThr,5);
	detectionThrLine_fluxDensityAccuracy_ext->SetLineColor(kBlack);
	detectionThrLine_fluxDensityAccuracy_ext->SetLineStyle(kDashed);
	detectionThrLine_fluxDensityAccuracy_ext->Draw("same");
	*/

	
}//close Draw()

double FitPosErrFcn(double* x,double* par){

	//double fluxDensity= pow(10,x[0]);
	double A= pow(10,x[0]);
	
	double scaleFactor= par[0];
	double bmaj= par[1];
	double bmin= par[2];
	//double sigmaX= bmaj/GausSigma2FWHM;
	//double sigmaY= bmin/GausSigma2FWHM;
	//double A= fluxDensity/(2*TMath::Pi()*sigmaX*sigmaY);

	double pixelSize= par[3];
	double noiseRMS= par[4];
	double fcn= scaleFactor*sqrt(2./TMath::Pi())*sqrt(bmaj/bmin)*pixelSize*noiseRMS/A;
	return fcn;
}


double GausPosSigmaFcn(double* x,double* par){

	double X= pow(10,x[0]);//SN
	double scaleFactor= par[0];
	double pixelSize= par[1];
	double sigmaX= par[2];
	double sigmaY= par[3];

	double fcn= scaleFactor*sqrt(2*sigmaX/(TMath::Pi()*sigmaY))*pixelSize/X;
	//cout<<"X="<<X<<" fcn="<<fcn<<endl;	
	
	return fcn;

}//close GausPosVarianceFcn()

