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

//## Beam & map info
double bkgLevel_true= 10.e-6;
double noiseLevel_true= 400.e-6;
double Bmaj= 13.347706794699599;//arcsec
double Bmin= 8.35635757446;//arcsec
double Bpa= 6.416473388670e-01;//deg
double pixSize= 1;//in arcsec
double beamArea= 0;
TEllipse* beamEllipse= 0;
double beamEllipseEccentricity= 0;
double beamEllipseArea= 0;

//### Selection cuts
bool selectResolvedTrueSources= false;
double mutualTrueSourceDistThr= 20;//in arcsec 
bool selectFitComponents= true;
int maxFitComponentsThr= 1;
bool excludeSourceAtEdge= true;
double sourceRegionXMin= 200;//in pixels
double sourceRegionXMax= 2360;//in pixels
double sourceRegionYMin= 200;//in pixels
double sourceRegionYMax= 2360;//in pixels

//Draw info
//std::vector<double> Zbins= {0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};
std::vector<double> Zbins= {0.1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,10000};

std::vector<double> LgFluxBins= {
	-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.5,0,1,2
};

//- Compact source histos/graphs
TH1D* NTrueSourceHisto_compact= 0;
TH1D* NRecSourceHisto_compact= 0;
TH1D* NRecSourceHisto_compact_fit= 0;
TH1D* NTrueSourceVSSignificanceHisto_compact= 0;
TH1D* NRecSourceVSSignificanceHisto_compact= 0;
TH1D* NRecSourceVSSignificanceHisto_compact_fit= 0;
TH1D* NTrueSourceHisto_reliability_compact= 0;
TH1D* NRecSourceHisto_reliability_compact= 0;
TH1D* NTrueSourceVSSignificanceHisto_reliability_compact= 0;
TH1D* NRecSourceVSSignificanceHisto_reliability_compact= 0;
TEfficiency* Efficiency_compact= 0;
TEfficiency* Efficiency_compact_fit= 0;
TEfficiency* Reliability_compact= 0;
TEfficiency* EfficiencyVSSignificance_compact= 0;
TEfficiency* EfficiencyVSSignificance_compact_fit= 0;
TEfficiency* ReliabilityVSSignificance_compact= 0;

TGraphAsymmErrors* xPosAccuracyGraph_compact_centroid= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact_centroid= 0;
TGraphAsymmErrors* xPosAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* xPosAccuracyGraph_compact= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact= 0;
TGraphAsymmErrors* xPosTrueAccuracyGraph_compact= 0;
TGraphAsymmErrors* yPosTrueAccuracyGraph_compact= 0;
TGraphAsymmErrors* xPosTrueAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* yPosTrueAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact_fit= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact_peak= 0;
TGraph* FluxDensityAccuracyGraphPoints_compact= 0;
TGraph* FitEllipseVSBeam_compact_real= 0;
TGraph* FitEllipseVSBeam_compact_fake= 0;


std::vector< std::vector<double> > FluxList_compact;
std::vector< std::vector<double> > FluxList_compact_fit;
std::vector< std::vector<double> > xPosPullList_compact_centroid;
std::vector< std::vector<double> > yPosPullList_compact_centroid;
std::vector< std::vector<double> > xPosPullList_compact_fit;
std::vector< std::vector<double> > yPosPullList_compact_fit;
std::vector< std::vector<double> > xPosPullList_compact;
std::vector< std::vector<double> > yPosPullList_compact;
std::vector< std::vector<double> > xPosTruePullList_compact;
std::vector< std::vector<double> > yPosTruePullList_compact;
std::vector< std::vector<double> > xPosTruePullList_compact_fit;
std::vector< std::vector<double> > yPosTruePullList_compact_fit;
std::vector< std::vector<double> > FluxDensityPullList_compact;
std::vector< std::vector<double> > FluxDensityPullList_compact_peak;
std::vector< std::vector<double> > FluxDensityPullList_compact_fit;
std::vector< std::vector<double> > FluxDensityTruePullList_compact;

//- Extended source histos/graphs
TH1D* NTrueSourceHisto_ext= 0;
TH1D* NRecSourceHisto_ext= 0;
TH1D* NTrueSourceHisto_reliability_ext= 0;
TH1D* NRecSourceHisto_reliability_ext= 0;
TEfficiency* Efficiency_ext= 0;
TEfficiency* Reliability_ext= 0;

TH1D* NTrueSourceVSSignificanceHisto_ext= 0;
TH1D* NRecSourceVSSignificanceHisto_ext= 0;
TH1D* NTrueSourceVSSignificanceHisto_reliability_ext= 0;
TH1D* NRecSourceVSSignificanceHisto_reliability_ext= 0;
TEfficiency* EfficiencyVSSignificance_ext= 0;
TEfficiency* ReliabilityVSSignificance_ext= 0;

TGraphAsymmErrors* xPosAccuracyGraph_ext= 0;
TGraphAsymmErrors* yPosAccuracyGraph_ext= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_ext= 0;

std::vector< std::vector<double> > FluxList_ext;
std::vector< std::vector<double> > xPosPullList_ext;
std::vector< std::vector<double> > yPosPullList_ext;
std::vector< std::vector<double> > FluxDensityPullList_ext;


//## Define functions
void Init();
void Draw();
int AnalyzeData(std::string inputFileName);
void ComputeAnalysisHistos();
double GausPosSigmaFcn(double* x,double* par);
double FitPosErrFcn(double* x,double* par);

//====================================
//        MACRO
//===================================
int DrawSFinderPerformances(std::string _fileName,bool isFileList=false){

	//## Set logging level
	LoggerManager::Instance().CreateConsoleLogger("INFO","logger","System.out");

	
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
			ERROR_LOG("Failed to analyze data for file no. "<<i+1);
			return -1;
		}
	}//end loop files

	//## Compute analysis histo
	ComputeAnalysisHistos();

	//## Draw 
	Draw();

	//Save data
	//...

	return 0;
	
}//close macro


void Init(){

	//Compute map variables
	double theta_beam= MathUtils::Mod(Bpa+90.,180.);
	beamArea= AstroUtils::GetBeamAreaInPixels(Bmaj,Bmin,pixSize,pixSize);
	beamEllipse= new TEllipse(0.,0.,Bmaj/GausSigma2FWHM,Bmin/GausSigma2FWHM,0,360,theta_beam);
	beamEllipseEccentricity= MathUtils::ComputeEllipseEccentricity(beamEllipse);
	beamEllipseArea= MathUtils::ComputeEllipseArea(beamEllipse);
	INFO_LOG("Beam ellipse: bmaj/bmin/bpa="<<Bmaj<<"/"<<Bmin<<"/"<<Bpa<<", a="<<std::min(beamEllipse->GetR1(),beamEllipse->GetR2())<<", b="<<std::max(beamEllipse->GetR1(),beamEllipse->GetR2())<<", theta="<<beamEllipse->GetTheta()<<", A="<<beamEllipseArea<<", E="<<beamEllipseEccentricity);
	
	int nBins= (int)(LgFluxBins.size()-1);
	int nBins_Z= (int)(Zbins.size()-1);
	
	//## Init data vector
	for(int i=0;i<nBins;i++){
		FluxList_compact.push_back( std::vector<double>() );
		FluxList_compact_fit.push_back( std::vector<double>() );
		xPosPullList_compact.push_back( std::vector<double>() );
		yPosPullList_compact.push_back( std::vector<double>() );
		xPosPullList_compact_fit.push_back( std::vector<double>() );
		yPosPullList_compact_fit.push_back( std::vector<double>() );
		xPosPullList_compact_centroid.push_back( std::vector<double>() );
		yPosPullList_compact_centroid.push_back( std::vector<double>() );
		xPosTruePullList_compact.push_back( std::vector<double>() );
		yPosTruePullList_compact.push_back( std::vector<double>() );
		xPosTruePullList_compact_fit.push_back( std::vector<double>() );
		yPosTruePullList_compact_fit.push_back( std::vector<double>() );
		FluxDensityPullList_compact.push_back( std::vector<double>() );
		FluxDensityPullList_compact_fit.push_back( std::vector<double>() );
		FluxDensityPullList_compact_peak.push_back( std::vector<double>() );
		FluxDensityTruePullList_compact.push_back( std::vector<double>() );
	
		FluxList_ext.push_back( std::vector<double>() );
		xPosPullList_ext.push_back( std::vector<double>() );
		yPosPullList_ext.push_back( std::vector<double>() );
		FluxDensityPullList_ext.push_back( std::vector<double>() );
	}

	//## Init histos	
	//- Completeness histos
	TString histoName= "NTrueSourceHisto_compact";
	NTrueSourceHisto_compact= new TH1D(histoName,"True compact source distribution",nBins,LgFluxBins.data());
	NTrueSourceHisto_compact->Sumw2();
	
	histoName= "NRecSourceHisto_compact";
	NRecSourceHisto_compact= new TH1D(histoName,"Rec compact source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_compact->Sumw2();
	
	histoName= "NRecSourceHisto_compact_fit";
	NRecSourceHisto_compact_fit= new TH1D(histoName,"Fitted compact source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_compact_fit->Sumw2();

	histoName= "NTrueSourceVSSignificanceHisto_compact";
	NTrueSourceVSSignificanceHisto_compact= new TH1D(histoName,"True compact source distribution",nBins,Zbins.data());
	NTrueSourceVSSignificanceHisto_compact->Sumw2();
	
	histoName= "NRecSourceVSSignificanceHisto_compact";
	NRecSourceVSSignificanceHisto_compact= new TH1D(histoName,"Rec compact source distribution",nBins,Zbins.data());
	NRecSourceVSSignificanceHisto_compact->Sumw2();
	
	histoName= "NRecSourceVSSignificanceHisto_compact_fit";
	NRecSourceVSSignificanceHisto_compact_fit= new TH1D(histoName,"Fitted compact source distribution",nBins,Zbins.data());
	NRecSourceVSSignificanceHisto_compact_fit->Sumw2();

	histoName= "NTrueSourceHisto_ext";
	NTrueSourceHisto_ext= new TH1D(histoName,"True extended source distribution",nBins,LgFluxBins.data());
	NTrueSourceHisto_ext->Sumw2();
	
	histoName= "NRecSourceHisto_ext";
	NRecSourceHisto_ext= new TH1D(histoName,"Rec extended source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_ext->Sumw2();

	histoName= "NTrueSourceVSSignificanceHisto_ext";
	NTrueSourceVSSignificanceHisto_ext= new TH1D(histoName,"True extended source distribution",nBins_Z,Zbins.data());
	NTrueSourceVSSignificanceHisto_ext->Sumw2();
	
	histoName= "NRecSourceVSSignificanceHisto_ext";
	NRecSourceVSSignificanceHisto_ext= new TH1D(histoName,"Rec extended source distribution",nBins_Z,Zbins.data());
	NRecSourceVSSignificanceHisto_ext->Sumw2();
		
	//- Reliability histos
	histoName= "NTrueSourceHisto_reliability_compact";
	NTrueSourceHisto_reliability_compact= new TH1D(histoName,"True compact source distribution",nBins,LgFluxBins.data());
	NTrueSourceHisto_reliability_compact->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_compact";
	NRecSourceHisto_reliability_compact= new TH1D(histoName,"Rec compact source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_reliability_compact->Sumw2();

	histoName= "NTrueSourceVSSignificanceHisto_reliability_compact";
	NTrueSourceVSSignificanceHisto_reliability_compact= new TH1D(histoName,"True compact source distribution",nBins,Zbins.data());
	NTrueSourceVSSignificanceHisto_reliability_compact->Sumw2();
		
	histoName= "NRecSourceVSSignificanceHisto_reliability_compact";
	NRecSourceVSSignificanceHisto_reliability_compact= new TH1D(histoName,"Rec compact source distribution",nBins,Zbins.data());
	NRecSourceVSSignificanceHisto_reliability_compact->Sumw2();
	
	histoName= "NTrueSourceHisto_reliability_ext";
	NTrueSourceHisto_reliability_ext= new TH1D(histoName,"True extended source distribution",nBins,LgFluxBins.data());
	NTrueSourceHisto_reliability_ext->Sumw2();
		
	histoName= "NRecSourceHisto_reliability_ext";
	NRecSourceHisto_reliability_ext= new TH1D(histoName,"Rec extended source distribution",nBins,LgFluxBins.data());
	NRecSourceHisto_reliability_ext->Sumw2();

	histoName= "NTrueSourceVSSignificanceHisto_reliability_ext";
	NTrueSourceVSSignificanceHisto_reliability_ext= new TH1D(histoName,"True extended source distribution",nBins,Zbins.data());
	NTrueSourceVSSignificanceHisto_reliability_ext->Sumw2();
		
	histoName= "NRecSourceVSSignificanceHisto_reliability_ext";
	NRecSourceVSSignificanceHisto_reliability_ext= new TH1D(histoName,"Rec extended source distribution",nBins,Zbins.data());
	NRecSourceVSSignificanceHisto_reliability_ext->Sumw2();

	//- Pos accuracy graphs
	TString graphName= "xPosAccuracyGraph_compact";
	xPosAccuracyGraph_compact= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact";
	yPosAccuracyGraph_compact= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact->SetName(graphName);
	
	graphName= "xPosAccuracyGraph_compact_fit";
	xPosAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact_fit->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact_fit";
	yPosAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact_fit->SetName(graphName);

	graphName= "xPosAccuracyGraph_compact_centroid";
	xPosAccuracyGraph_compact_centroid= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact_centroid->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact_centroid";
	yPosAccuracyGraph_compact_centroid= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact_centroid->SetName(graphName);

	graphName= "xPosTrueAccuracyGraph_compact";
	xPosTrueAccuracyGraph_compact= new TGraphAsymmErrors;
	xPosTrueAccuracyGraph_compact->SetName(graphName);
	
	graphName= "yPosTrueAccuracyGraph_compact";
	yPosTrueAccuracyGraph_compact= new TGraphAsymmErrors;
	yPosTrueAccuracyGraph_compact->SetName(graphName);

	graphName= "xPosTrueAccuracyGraph_compact_fit";
	xPosTrueAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	xPosTrueAccuracyGraph_compact_fit->SetName(graphName);
	
	graphName= "yPosTrueAccuracyGraph_compact_fit";
	yPosTrueAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	yPosTrueAccuracyGraph_compact_fit->SetName(graphName);

	graphName= "xPosAccuracyGraph_ext";
	xPosAccuracyGraph_ext= new TGraphAsymmErrors;
	xPosAccuracyGraph_ext->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_ext";
	yPosAccuracyGraph_ext= new TGraphAsymmErrors;
	yPosAccuracyGraph_ext->SetName(graphName);

	//- Flux accuracy graphs
	graphName= "FluxDensityAccuracyGraph_compact";
	FluxDensityAccuracyGraph_compact= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact->SetName(graphName);
	
	graphName= "FluxDensityAccuracyGraph_compact_fit";
	FluxDensityAccuracyGraph_compact_fit= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact_fit->SetName(graphName);

	graphName= "FluxDensityAccuracyGraph_compact_peak";
	FluxDensityAccuracyGraph_compact_peak= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact_peak->SetName(graphName);

	graphName= "FluxDensityAccuracyGraphPoints_compact";
	FluxDensityAccuracyGraphPoints_compact= new TGraph;
	FluxDensityAccuracyGraphPoints_compact->SetName(graphName);
	

	graphName= "FluxDensityAccuracyGraph_ext";
	FluxDensityAccuracyGraph_ext= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_ext->SetName(graphName);

	//- Fit ellipse graph
	graphName= "FitEllipseVSBeam_compact_real";
	FitEllipseVSBeam_compact_real= new TGraph;
	FitEllipseVSBeam_compact_real->SetName(graphName);

	graphName= "FitEllipseVSBeam_compact_fake";
	FitEllipseVSBeam_compact_fake= new TGraph;
	FitEllipseVSBeam_compact_fake->SetName(graphName);

}//close Init()



void ComputeAnalysisHistos(){

	//Compute detection efficiency
	Efficiency_compact = new TEfficiency(*NRecSourceHisto_compact,*NTrueSourceHisto_compact); 
	Efficiency_compact->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Efficiency_compact->SetConfidenceLevel(0.68);
	
	Efficiency_compact_fit = new TEfficiency(*NRecSourceHisto_compact_fit,*NTrueSourceHisto_compact); 
	Efficiency_compact_fit->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Efficiency_compact_fit->SetConfidenceLevel(0.68);

	EfficiencyVSSignificance_compact = new TEfficiency(*NRecSourceVSSignificanceHisto_compact,*NTrueSourceVSSignificanceHisto_compact); 
	EfficiencyVSSignificance_compact->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	EfficiencyVSSignificance_compact->SetConfidenceLevel(0.68);
	
	EfficiencyVSSignificance_compact_fit = new TEfficiency(*NRecSourceVSSignificanceHisto_compact_fit,*NTrueSourceVSSignificanceHisto_compact); 
	EfficiencyVSSignificance_compact_fit->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	EfficiencyVSSignificance_compact_fit->SetConfidenceLevel(0.68);
	
	Efficiency_ext = new TEfficiency(*NRecSourceHisto_ext,*NTrueSourceHisto_ext); 
	Efficiency_ext->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Efficiency_ext->SetConfidenceLevel(0.68);

	EfficiencyVSSignificance_ext = new TEfficiency(*NRecSourceVSSignificanceHisto_ext,*NTrueSourceVSSignificanceHisto_ext); 
	EfficiencyVSSignificance_ext->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	EfficiencyVSSignificance_ext->SetConfidenceLevel(0.68);

	//Compute reliability 
	Reliability_compact= new TEfficiency(*NTrueSourceHisto_reliability_compact,*NRecSourceHisto_reliability_compact); 
	Reliability_compact->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Reliability_compact->SetConfidenceLevel(0.68);
	
	Reliability_ext= new TEfficiency(*NTrueSourceHisto_reliability_ext,*NRecSourceHisto_reliability_ext); 
	Reliability_ext->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Reliability_ext->SetConfidenceLevel(0.68);

	ReliabilityVSSignificance_compact= new TEfficiency(*NTrueSourceVSSignificanceHisto_reliability_compact,*NRecSourceVSSignificanceHisto_reliability_compact); 
	ReliabilityVSSignificance_compact->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	ReliabilityVSSignificance_compact->SetConfidenceLevel(0.68);
	
	ReliabilityVSSignificance_ext= new TEfficiency(*NTrueSourceVSSignificanceHisto_reliability_ext,*NRecSourceVSSignificanceHisto_reliability_ext); 
	ReliabilityVSSignificance_ext->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	ReliabilityVSSignificance_ext->SetConfidenceLevel(0.68);
	
	//Compute data list stats
	int nMinPointsToDraw= 2;
	int nBins= (int)(LgFluxBins.size()-1);
	int nPoints_xPull= 0;
	int nPoints_yPull= 0;	
	int nPoints_xPull_fit= 0;
	int nPoints_yPull_fit= 0;
	int nPoints_xPull_true= 0;
	int nPoints_yPull_true= 0;
	int nPoints_xPull_true_fit= 0;
	int nPoints_yPull_true_fit= 0;		
	int nPoints_SpeakPull= 0;
	int nPoints_FluxPull= 0;
	int nPoints_FluxPull_all= 0;
	int nPoints_OffsetPull= 0;
	int nPoints_xPull_ext= 0;
	int nPoints_yPull_ext= 0;	
	int nPoints_FluxPull_ext= 0;

	for(int i=0;i<nBins;i++){
		double x_avg= LgFluxBins[i] + 0.5*(LgFluxBins[i+1]-LgFluxBins[i]);

		//Compute xpos accuracy		
		if(xPosPullList_compact[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(xPosPullList_compact[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			xPosAccuracyGraph_compact->SetPoint(nPoints_xPull,x,y);
			xPosAccuracyGraph_compact->SetPointError(nPoints_xPull,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_xPull++;
		}

		//Compute ypos accuracy
		if(yPosPullList_compact[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(yPosPullList_compact[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			yPosAccuracyGraph_compact->SetPoint(nPoints_yPull,x,y);
			yPosAccuracyGraph_compact->SetPointError(nPoints_yPull,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_yPull++;
		}	

		//Compute xpos accuracy	fit	
		if(xPosPullList_compact_fit[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(xPosPullList_compact_fit[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			xPosAccuracyGraph_compact_fit->SetPoint(nPoints_xPull_fit,x,y);
			xPosAccuracyGraph_compact_fit->SetPointError(nPoints_xPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_xPull_fit++;
		}

		//Compute ypos accuracy fit
		if(yPosPullList_compact_fit[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(yPosPullList_compact_fit[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			yPosAccuracyGraph_compact_fit->SetPoint(nPoints_yPull_fit,x,y);
			yPosAccuracyGraph_compact_fit->SetPointError(nPoints_yPull_fit,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_yPull_fit++;
		}	

		//Compute xpos accuracy (wrt true pos)
		if(xPosPullList_compact[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(xPosTruePullList_compact[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			xPosTrueAccuracyGraph_compact->SetPoint(nPoints_xPull_true,x,y);
			xPosTrueAccuracyGraph_compact->SetPointError(nPoints_xPull_true,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_xPull_true++;
		}

		//Compute ypos accuracy (wrt true pos)
		if(yPosPullList_compact[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(yPosTruePullList_compact[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			yPosTrueAccuracyGraph_compact->SetPoint(nPoints_yPull_true,x,y);
			yPosTrueAccuracyGraph_compact->SetPointError(nPoints_yPull_true,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_yPull_true++;
		}	
			
		//Compute xpos accuracy fit (wrt true pos)
		if(xPosPullList_compact_fit[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(xPosTruePullList_compact_fit[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			xPosTrueAccuracyGraph_compact_fit->SetPoint(nPoints_xPull_true_fit,x,y);
			xPosTrueAccuracyGraph_compact_fit->SetPointError(nPoints_xPull_true_fit,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_xPull_true_fit++;
		}

		//Compute ypos accuracy fit (wrt true pos)
		if(yPosPullList_compact_fit[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(yPosTruePullList_compact_fit[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact_fit[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			yPosTrueAccuracyGraph_compact_fit->SetPoint(nPoints_yPull_true_fit,x,y);
			yPosTrueAccuracyGraph_compact_fit->SetPointError(nPoints_yPull_true_fit,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_yPull_true_fit++;
		}	
			

		//Compute Flux accuracy
		for(size_t j=0;j<FluxList_compact[i].size();j++){
			FluxDensityAccuracyGraphPoints_compact->SetPoint(nPoints_FluxPull_all,FluxList_compact[i][j],FluxDensityPullList_compact[i][j]);
			nPoints_FluxPull_all++;
		}
		if(FluxDensityPullList_compact[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(FluxDensityPullList_compact[i]);	
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			FluxDensityAccuracyGraph_compact->SetPoint(nPoints_FluxPull,x,y);
			FluxDensityAccuracyGraph_compact->SetPointError(nPoints_FluxPull,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_FluxPull++;
		}	
		

		//Compute Flux accuracy for extended sources
		if(FluxDensityPullList_ext[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(FluxDensityPullList_ext[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_ext[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
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
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			yPosAccuracyGraph_ext->SetPoint(nPoints_yPull_ext,x,y);
			yPosAccuracyGraph_ext->SetPointError(nPoints_yPull_ext,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_yPull_ext++;
		}	


	}//end loop bins

}//close ComputeAnalysisHistos()


void Draw(){

	gROOT->SetStyle("myStyle2");
	//gStyle->SetOptLogx();

	//double xMin_draw= LgFluxBins[0];
	//double xMax_draw= LgFluxBins[LgFluxBins.size()-1];
	double xMin_draw= -4;
	double xMax_draw= 1;
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
		double fluxThrRef= zRefThr*noiseLevel_true;
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
	EffPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	EffPlotBkg->GetYaxis()->SetTitle("Completeness");
	EffPlotBkg->Draw();

	double Eff_detectionAreaX[]= {xMin_draw,log10(5*noiseLevel_true)};
	double Eff_detectionAreaY[]= {0,0};
	TGraph* Eff_sigmaDetectionArea = new TGraph(2,Eff_detectionAreaX,Eff_detectionAreaY);
  Eff_sigmaDetectionArea->SetLineColor(kGray+2);
  Eff_sigmaDetectionArea->SetLineWidth(15001);
  Eff_sigmaDetectionArea->SetFillStyle(3005);	
	Eff_sigmaDetectionArea->SetFillColor(kGray);
	Eff_sigmaDetectionArea->Draw("C");

	Efficiency_compact->SetMarkerStyle(8);
	Efficiency_compact->SetMarkerColor(kBlack);
	Efficiency_compact->SetLineColor(kBlack);
	Efficiency_compact->Draw("p same");

	Efficiency_compact_fit->SetMarkerStyle(21);
	Efficiency_compact_fit->SetMarkerColor(kRed);
	Efficiency_compact_fit->SetLineColor(kRed);
	Efficiency_compact_fit->Draw("p same");

	//for(size_t i=0;i<detectionThrLines.size();i++) {
	//	detectionThrLines[i]->Draw("l same");
	//}

	TLegend* EffPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	EffPlotLegend->SetFillColor(0);
	EffPlotLegend->SetTextSize(0.045);
	EffPlotLegend->SetTextFont(52);
	EffPlotLegend->AddEntry(Efficiency_compact,"detected","PL");
	EffPlotLegend->AddEntry(Efficiency_compact_fit,"fitted","PL");
	EffPlotLegend->AddEntry(Eff_sigmaDetectionArea,"S_{true}<5#sigma","F");
	EffPlotLegend->Draw("same");


	TCanvas* EffVSSignificancePlot= new TCanvas("EffVSSignificancePlot","EffVSSignificancePlot");
	EffVSSignificancePlot->cd();

	TH2D* EffVSSignificancePlotBkg= new TH2D("EffVSSignificancePlotBkg","",100,ZMin_draw,ZMax_draw,100,0,1.1);
	EffVSSignificancePlotBkg->GetXaxis()->SetTitle("S/N");
	EffVSSignificancePlotBkg->GetYaxis()->SetTitle("Completeness");
	EffVSSignificancePlotBkg->Draw();

	EfficiencyVSSignificance_compact->SetMarkerStyle(8);
	EfficiencyVSSignificance_compact->SetMarkerColor(kBlack);
	EfficiencyVSSignificance_compact->SetLineColor(kBlack);
	EfficiencyVSSignificance_compact->Draw("p same");

	EfficiencyVSSignificance_compact_fit->SetMarkerStyle(21);
	EfficiencyVSSignificance_compact_fit->SetMarkerColor(kRed);
	EfficiencyVSSignificance_compact_fit->SetLineColor(kRed);
	EfficiencyVSSignificance_compact_fit->Draw("p same");

	TLegend* EffVSSignificancePlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	EffVSSignificancePlotLegend->SetFillColor(0);
	EffVSSignificancePlotLegend->SetTextSize(0.045);
	EffVSSignificancePlotLegend->SetTextFont(52);
	EffVSSignificancePlotLegend->AddEntry(EfficiencyVSSignificance_compact,"detected","PL");
	EffVSSignificancePlotLegend->AddEntry(EfficiencyVSSignificance_compact_fit,"fitted","PL");
	EffVSSignificancePlotLegend->Draw("same");



	//- Efficiency plot extended source
	TCanvas* EffPlot_ext= new TCanvas("EffPlot_ext","EffPlot_ext");
	EffPlot_ext->cd();

	TH2D* EffPlotBkg_ext= new TH2D("EffPlotBkg_ext","",100,xMin_draw,xMax_draw,100,0,1.2);
	EffPlotBkg_ext->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	EffPlotBkg_ext->GetYaxis()->SetTitle("Completeness");
	EffPlotBkg_ext->Draw();

	Efficiency_ext->SetMarkerStyle(8);
	Efficiency_ext->SetMarkerColor(kBlack);
	Efficiency_ext->SetLineColor(kBlack);
	Efficiency_ext->Draw("p same");
	
	

	TCanvas* EffVSSignificancePlot_ext= new TCanvas("EffVSSignificancePlot_ext","EffVSSignificancePlot_ext");
	EffVSSignificancePlot_ext->cd();

	TH2D* EffVSSignificancePlotBkg_ext= new TH2D("EffVSSignificancePlotBkg_ext","",100,ZMin_draw,ZMax_draw,100,0,1.2);
	EffVSSignificancePlotBkg_ext->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	EffVSSignificancePlotBkg_ext->GetYaxis()->SetTitle("Completeness");
	EffVSSignificancePlotBkg_ext->Draw();

	EfficiencyVSSignificance_ext->SetMarkerStyle(8);
	EfficiencyVSSignificance_ext->SetMarkerColor(kBlack);
	EfficiencyVSSignificance_ext->SetLineColor(kBlack);
	EfficiencyVSSignificance_ext->Draw("p same");
	
	//===============================================
	//==          DRAW RELIABILITY
	//===============================================
	
	//Compact
	TCanvas* ReliabilityPlot= new TCanvas("ReliabilityPlot","ReliabilityPlot");
	ReliabilityPlot->cd();

	TH2D* ReliabilityPlotBkg= new TH2D("ReliabilityPlotBkg","",100,xMin_draw,xMax_draw,100,0,1.2);
	ReliabilityPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{rec}/Jy)");
	ReliabilityPlotBkg->GetYaxis()->SetTitle("Reliability");
	ReliabilityPlotBkg->Draw();

	Eff_sigmaDetectionArea->Draw("C");

	Reliability_compact->SetMarkerStyle(8);
	Reliability_compact->SetMarkerColor(kBlack);
	Reliability_compact->SetLineColor(kBlack);
	Reliability_compact->Draw("p same");
	
	TLegend* ReliabilityPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	ReliabilityPlotLegend->SetFillColor(0);
	ReliabilityPlotLegend->SetTextSize(0.045);
	ReliabilityPlotLegend->SetTextFont(52);
	ReliabilityPlotLegend->AddEntry(Reliability_compact,"point-source reliability","PL");
	ReliabilityPlotLegend->AddEntry(Eff_sigmaDetectionArea,"S_{rec}<5#sigma","F");
	ReliabilityPlotLegend->Draw("same");

	TCanvas* ReliabilityVSSignificancePlot= new TCanvas("ReliabilityVSSignificancePlot","ReliabilityVSSignificancePlot");
	ReliabilityVSSignificancePlot->cd();

	TH2D* ReliabilityVSSignificancePlotBkg= new TH2D("ReliabilityVSSignificancePlotBkg","",100,ZMin_draw,ZMax_draw,100,0,1.2);
	ReliabilityVSSignificancePlotBkg->GetXaxis()->SetTitle("S/N");
	ReliabilityVSSignificancePlotBkg->GetYaxis()->SetTitle("Reliability");
	ReliabilityVSSignificancePlotBkg->Draw();

	ReliabilityVSSignificance_compact->SetMarkerStyle(8);
	ReliabilityVSSignificance_compact->SetMarkerColor(kBlack);
	ReliabilityVSSignificance_compact->SetLineColor(kBlack);
	ReliabilityVSSignificance_compact->Draw("p same");
	

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

	TCanvas* ReliabilityVSSignificancePlot_ext= new TCanvas("ReliabilityVSSignificancePlot_ext","ReliabilityVSSignificancePlot_ext");
	ReliabilityVSSignificancePlot_ext->cd();

	TH2D* ReliabilityVSSignificancePlotBkg_ext= new TH2D("ReliabilityVSSignificancePlotBkg_ext","",100,ZMin_draw,ZMax_draw,100,0,1.2);
	ReliabilityVSSignificancePlotBkg_ext->GetXaxis()->SetTitle("S/N");
	ReliabilityVSSignificancePlotBkg_ext->GetYaxis()->SetTitle("Reliability");
	ReliabilityVSSignificancePlotBkg_ext->Draw();

	ReliabilityVSSignificance_ext->SetMarkerStyle(8);
	ReliabilityVSSignificance_ext->SetMarkerColor(kBlack);
	ReliabilityVSSignificance_ext->SetLineColor(kBlack);
	ReliabilityVSSignificance_ext->Draw("p same");

	//===============================================
	//==          DRAW POSITIONAL ACCURACY
	//===============================================
	TCanvas* PosAccuracyPlot= new TCanvas("PosAccuracyPlot","PosAccuracyPlot");
	PosAccuracyPlot->cd();

	TH2D* PosAccuracyPlotBkg= new TH2D("PosAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-20,20);
	PosAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	PosAccuracyPlotBkg->GetYaxis()->SetTitle("<#deltaRA>, <#deltaDec> ('')");
	PosAccuracyPlotBkg->Draw();

	xPosAccuracyGraph_compact->SetMarkerStyle(8);
	xPosAccuracyGraph_compact->SetMarkerColor(kBlack);
	xPosAccuracyGraph_compact->SetLineColor(kBlack);
	xPosAccuracyGraph_compact->Draw("ep same");
	
	yPosAccuracyGraph_compact->SetMarkerStyle(24);
	yPosAccuracyGraph_compact->SetMarkerColor(kBlack);
	yPosAccuracyGraph_compact->SetLineColor(kBlack);
	yPosAccuracyGraph_compact->Draw("ep same");

	xPosAccuracyGraph_compact_fit->SetMarkerStyle(21);
	xPosAccuracyGraph_compact_fit->SetMarkerColor(kRed);
	xPosAccuracyGraph_compact_fit->SetLineColor(kRed);
	xPosAccuracyGraph_compact_fit->Draw("ep same");
	
	yPosAccuracyGraph_compact_fit->SetMarkerStyle(25);
	yPosAccuracyGraph_compact_fit->SetMarkerColor(kRed);
	yPosAccuracyGraph_compact_fit->SetLineColor(kRed);
	yPosAccuracyGraph_compact_fit->Draw("ep same");
	
	TLine* refLine_posAccuracy= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_posAccuracy->SetLineColor(kBlack);
	refLine_posAccuracy->SetLineStyle(kDashed);
	refLine_posAccuracy->Draw("same");

	double sigmaToFirstQuantile= ROOT::Math::gaussian_quantile(0.75,1);
	TF1* expPosSigmaFcn_plus= new TF1("expPosSigmaFcn_plus",GausPosSigmaFcn,xMin_draw,xMax_draw,4);
	expPosSigmaFcn_plus->SetParameters(sigmaToFirstQuantile,pixSize,Bmaj,Bmin);
	expPosSigmaFcn_plus->SetLineColor(kBlack);
	expPosSigmaFcn_plus->SetLineStyle(9);
	//expPosSigmaFcn_plus->Draw("l same");

	TF1* expPosSigmaFcn_minus= new TF1("expPosSigmaFcn_minus",GausPosSigmaFcn,xMin_draw,xMax_draw,4);
	expPosSigmaFcn_minus->SetParameters(-sigmaToFirstQuantile,pixSize,Bmaj,Bmin);
	expPosSigmaFcn_minus->SetLineColor(kBlack);
	expPosSigmaFcn_minus->SetLineStyle(9);
	//expPosSigmaFcn_minus->Draw("l same");


	
	TLegend* PosAccuracyPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	PosAccuracyPlotLegend->SetFillColor(0);
	PosAccuracyPlotLegend->SetTextSize(0.045);
	PosAccuracyPlotLegend->SetTextFont(52);
	PosAccuracyPlotLegend->AddEntry(xPosAccuracyGraph_compact,"<#deltaRA>","PL");
	PosAccuracyPlotLegend->AddEntry(yPosAccuracyGraph_compact,"<#deltaDec>","PL");
	PosAccuracyPlotLegend->Draw("same");

	TLegend* PosAccuracyPlotLegend2= new TLegend(0.2,0.7,0.3,0.8);
	PosAccuracyPlotLegend2->SetFillColor(0);
	PosAccuracyPlotLegend2->SetTextSize(0.045);
	PosAccuracyPlotLegend2->SetTextFont(52);
	PosAccuracyPlotLegend2->SetHeader("Fitted");
	PosAccuracyPlotLegend2->AddEntry(xPosAccuracyGraph_compact_fit,"<#deltaRA>","PL");
	PosAccuracyPlotLegend2->AddEntry(yPosAccuracyGraph_compact_fit,"<#deltaDec>","PL");
	//PosAccuracyPlotLegend2->Draw("same");

	//True pos accuracy
	TCanvas* TruePosAccuracyPlot= new TCanvas("TruePosAccuracyPlot","TruePosAccuracyPlot");
	TruePosAccuracyPlot->cd();

	TH2D* TruePosAccuracyPlotBkg= new TH2D("TruePosAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-6,6);
	TruePosAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	TruePosAccuracyPlotBkg->GetYaxis()->SetTitle("<#deltaRA>, <#deltaDec> ('')");
	TruePosAccuracyPlotBkg->Draw();

	double detectionAreaX_pos[]= {xMin_draw,log10(5*noiseLevel_true)};
	double detectionAreaY_pos[]= {-6,-6};
	TGraph* sigmaDetectionArea_pos = new TGraph(2,detectionAreaX_pos,detectionAreaY_pos);
  sigmaDetectionArea_pos->SetLineColor(kGray+2);
  sigmaDetectionArea_pos->SetLineWidth(15001);
  sigmaDetectionArea_pos->SetFillStyle(3005);	
	sigmaDetectionArea_pos->SetFillColor(kGray);
	sigmaDetectionArea_pos->Draw("C");
	
	//expPosSigmaFcn_plus->Draw("l same");
	//expPosSigmaFcn_minus->Draw("l same");

	//TF1* posErr_minus= new TF1("posErr_minus","[0]*sqrt(2/TMath::Pi())*sqrt([1]/[2])*[3]*[4]/pow(10,x)",xMin_draw,xMax_draw,5);
	TF1* posErr_minus= new TF1("posErr_minus",FitPosErrFcn,xMin_draw,xMax_draw,5);
	posErr_minus->SetNpx(50);
	posErr_minus->SetParameters(-1,Bmaj,Bmin,pixSize,noiseLevel_true);
	posErr_minus->SetLineColor(kBlack);
	posErr_minus->SetLineStyle(9);
	posErr_minus->Draw("C same");

	//TF1* posErr_plus= new TF1("posErr_plus","sqrt(2/TMath::Pi())*[0]*[1]/pow(10,x)",xMin_draw,xMax_draw,3);
	TF1* posErr_plus= new TF1("posErr_plus",FitPosErrFcn,xMin_draw,xMax_draw,5);
	posErr_plus->SetNpx(50);
	posErr_plus->SetParameters(1,Bmaj,Bmin,pixSize,noiseLevel_true);
	posErr_plus->SetLineColor(kBlack);
	posErr_plus->SetLineStyle(9);
	posErr_plus->Draw("C same");
	

	xPosTrueAccuracyGraph_compact->SetMarkerStyle(8);
	xPosTrueAccuracyGraph_compact->SetMarkerColor(kBlack);
	xPosTrueAccuracyGraph_compact->SetLineColor(kBlack);
	xPosTrueAccuracyGraph_compact->Draw("ep same");
	
	yPosTrueAccuracyGraph_compact->SetMarkerStyle(24);
	yPosTrueAccuracyGraph_compact->SetMarkerColor(kBlack);
	yPosTrueAccuracyGraph_compact->SetLineColor(kBlack);
	yPosTrueAccuracyGraph_compact->Draw("ep same");

	xPosTrueAccuracyGraph_compact_fit->SetMarkerStyle(21);
	xPosTrueAccuracyGraph_compact_fit->SetMarkerColor(kRed);
	xPosTrueAccuracyGraph_compact_fit->SetLineColor(kRed);
	//xPosTrueAccuracyGraph_compact_fit->Draw("ep same");
	
	yPosTrueAccuracyGraph_compact_fit->SetMarkerStyle(25);
	yPosTrueAccuracyGraph_compact_fit->SetMarkerColor(kRed);
	yPosTrueAccuracyGraph_compact_fit->SetLineColor(kRed);
	//yPosTrueAccuracyGraph_compact_fit->Draw("ep same");
	
	refLine_posAccuracy->Draw("same");
	//expPosSigmaFcn_plus->Draw("l same");
	//expPosSigmaFcn_minus->Draw("l same");

	TLegend* TruePosAccuracyPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	TruePosAccuracyPlotLegend->SetFillColor(0);
	TruePosAccuracyPlotLegend->SetTextSize(0.045);
	TruePosAccuracyPlotLegend->SetTextFont(52);	
	//TruePosAccuracyPlotLegend->SetHeader("Detected");
	TruePosAccuracyPlotLegend->AddEntry(xPosTrueAccuracyGraph_compact,"<#deltaRA>","PL");
	TruePosAccuracyPlotLegend->AddEntry(yPosTrueAccuracyGraph_compact,"<#deltaDec>","PL");
	TruePosAccuracyPlotLegend->AddEntry(posErr_plus,"#pm #sigma_{pos}^{exp}","L");
	TruePosAccuracyPlotLegend->AddEntry(sigmaDetectionArea_pos,"S_{true}<5#sigma","F");
	TruePosAccuracyPlotLegend->Draw("same");

	TLegend* TruePosAccuracyPlotLegend2= new TLegend(0.6,0.7,0.7,0.8);
	TruePosAccuracyPlotLegend2->SetFillColor(0);
	TruePosAccuracyPlotLegend2->SetTextSize(0.045);
	TruePosAccuracyPlotLegend2->SetTextFont(52);
	TruePosAccuracyPlotLegend2->SetHeader("Fitted");
	TruePosAccuracyPlotLegend2->AddEntry(xPosTrueAccuracyGraph_compact_fit,"<#deltaRA>","PL");
	TruePosAccuracyPlotLegend2->AddEntry(yPosTrueAccuracyGraph_compact_fit,"<#deltaDec>","PL");
	//TruePosAccuracyPlotLegend2->Draw("same");
	
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
	//- FLux accuracy
	TCanvas* FluxAccuracyPlot= new TCanvas("FluxAccuracyPlot","FluxAccuracyPlot");
	FluxAccuracyPlot->cd();

	TH2D* FluxAccuracyPlotBkg= new TH2D("FluxAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-0.5,0.5);
	FluxAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	FluxAccuracyPlotBkg->GetYaxis()->SetTitle("<S_{rec}/S_{true}-1>");
	FluxAccuracyPlotBkg->Draw();
	
	double detectionAreaX[]= {xMin_draw,log10(5*noiseLevel_true)};
	double detectionAreaY[]= {-0.5,-0.5};
	TGraph* sigmaDetectionArea = new TGraph(2,detectionAreaX,detectionAreaY);
  sigmaDetectionArea->SetLineColor(kGray+2);
  sigmaDetectionArea->SetLineWidth(15001);
  sigmaDetectionArea->SetFillStyle(3005);	
	sigmaDetectionArea->SetFillColor(kGray);
	sigmaDetectionArea->Draw("C");

	FluxDensityAccuracyGraphPoints_compact->SetMarkerStyle(1);
	FluxDensityAccuracyGraphPoints_compact->SetMarkerColor(kGray);
	FluxDensityAccuracyGraphPoints_compact->SetLineColor(kGray);
	FluxDensityAccuracyGraphPoints_compact->Draw("p same");

	FluxDensityAccuracyGraph_compact->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_compact->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_compact->SetLineColor(kBlack);
	FluxDensityAccuracyGraph_compact->Draw("ep same");

	TLine* refLine_fluxAccuracy= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_fluxAccuracy->SetLineColor(kBlack);
	refLine_fluxAccuracy->SetLineStyle(kDashed);
	refLine_fluxAccuracy->Draw("same");


	TF1* fluxDensityErrLine_1sigma_plus= new TF1("fluxDensityErrLine_1sigma_plus","[0]*[1]/pow(10,x)",xMin_draw,xMax_draw);
	fluxDensityErrLine_1sigma_plus->SetParameters(noiseLevel_true,1);
	fluxDensityErrLine_1sigma_plus->SetLineColor(kBlack);
	fluxDensityErrLine_1sigma_plus->SetLineStyle(kDashed);
	fluxDensityErrLine_1sigma_plus->SetLineWidth(2);
	fluxDensityErrLine_1sigma_plus->Draw("l same");
	
	TF1* fluxDensityErrLine_1sigma_minus= new TF1("fluxDensityErrLine_1sigma_minus","[0]*[1]/pow(10,x)",xMin_draw,xMax_draw);
	fluxDensityErrLine_1sigma_minus->SetParameters(noiseLevel_true,-1);
	fluxDensityErrLine_1sigma_minus->SetLineColor(kBlack);
	fluxDensityErrLine_1sigma_minus->SetLineStyle(kDashed);
	fluxDensityErrLine_1sigma_minus->SetLineWidth(2);
	fluxDensityErrLine_1sigma_minus->Draw("l same");

	TF1* fluxDensityErrLine_3sigma_plus= new TF1("fluxDensityErrLine_3sigma_plus","[0]*[1]/pow(10,x)",xMin_draw,xMax_draw);
	fluxDensityErrLine_3sigma_plus->SetParameters(noiseLevel_true,3);
	fluxDensityErrLine_3sigma_plus->SetLineColor(kBlack);
	fluxDensityErrLine_3sigma_plus->SetLineStyle(kDotted);
	fluxDensityErrLine_3sigma_plus->SetLineWidth(2);
	fluxDensityErrLine_3sigma_plus->Draw("l same");
	
	TF1* fluxDensityErrLine_3sigma_minus= new TF1("fluxDensityErrLine_3sigma_minus","[0]*[1]/pow(10,x)",xMin_draw,xMax_draw);
	fluxDensityErrLine_3sigma_minus->SetParameters(noiseLevel_true,-3);
	fluxDensityErrLine_3sigma_minus->SetLineColor(kBlack);
	fluxDensityErrLine_3sigma_minus->SetLineStyle(kDotted);
	fluxDensityErrLine_3sigma_minus->SetLineWidth(2);
	fluxDensityErrLine_3sigma_minus->Draw("l same");


	TLegend* FluxAccuracyPlotLegend= new TLegend(0.6,0.7,0.7,0.8);
	FluxAccuracyPlotLegend->SetFillColor(0);
	FluxAccuracyPlotLegend->SetTextSize(0.045);
	FluxAccuracyPlotLegend->SetTextFont(52);
	FluxAccuracyPlotLegend->AddEntry(FluxDensityAccuracyGraph_compact,"<S_{rec}/S_{true}-1>","PL");
	FluxAccuracyPlotLegend->AddEntry(fluxDensityErrLine_1sigma_plus,"#pm 1#sigma","L");
	FluxAccuracyPlotLegend->AddEntry(fluxDensityErrLine_3sigma_plus,"#pm 3#sigma","L");
	FluxAccuracyPlotLegend->AddEntry(sigmaDetectionArea,"S_{true}<5#sigma","F");
	FluxAccuracyPlotLegend->Draw("same");


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

	//===============================================
	//==          DRAW FIT ELLIPSE VS BEAM PARS
	//===============================================
	TCanvas* FitEllipseVSBeamPlot_compact= new TCanvas("FitEllipseVSBeamPlot_compact","FitEllipseVSBeamPlot_compact");
	FitEllipseVSBeamPlot_compact->cd();
	
	FitEllipseVSBeam_compact_real->GetXaxis()->SetTitle("d#theta (deg)");
	FitEllipseVSBeam_compact_real->GetYaxis()->SetTitle("E_{fit}/E_{beam}");
	FitEllipseVSBeam_compact_real->SetMarkerStyle(8);
	FitEllipseVSBeam_compact_real->SetMarkerSize(0.7);
	FitEllipseVSBeam_compact_real->Draw("AP");

	FitEllipseVSBeam_compact_fake->SetMarkerColor(kRed);
	FitEllipseVSBeam_compact_fake->SetMarkerStyle(8);
	FitEllipseVSBeam_compact_fake->SetMarkerSize(0.7);
	FitEllipseVSBeam_compact_fake->Draw("P");

}//close Draw()


int AnalyzeData(std::string filename){

	INFO_LOG("Analyzing input data "<<filename<<" ...");

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

	//============================================================
	//===       SELECT TRUE SOURCES
	//============================================================	
	//Select true point-like sources to be included in the analysis?
	std::vector<bool> trueSourceSelectionFlags(SourceMatchInfo->GetEntries(),true);
	if(selectResolvedTrueSources){
		for(int i=0;i<SourceMatchInfo->GetEntries();i++){
			SourceMatchInfo->GetEntry(i);		
			double X_i= X0_true;
			double Y_i= Y0_true;
			for(int j=i+1;j<SourceMatchInfo->GetEntries();j++){
				SourceMatchInfo->GetEntry(j);
				double X_j= X0_true;
				double Y_j= Y0_true;
				double dist= sqrt( (X_i-X_j)*(X_i-X_j) + (Y_i-Y_j)*(Y_i-Y_j) );
				if(dist<mutualTrueSourceDistThr/pixSize){//mark both sources as not selected (not resolved)
					trueSourceSelectionFlags[i]= false;
					trueSourceSelectionFlags[j]= false;
				}
			}//end loop next sources
		}//end loop sources
	}//close if selectResolvedTrueSources



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

		//Apply selection cuts to true source
		//if(!trueSourceSelectionFlags[i]) continue;//skip true source

		//Fill true info
		double fluxDensity_true= S_true;
		//double fluxDensity_true= S/beamArea_true;
		double lgFlux_true= log10(fluxDensity_true);
		double Z_true= fluxDensity_true/noiseLevel_true;
		//double Z_true= fluxDensity_true/AvgRMS;

		//INFO_LOG("S="<<S<<" fluxDensity_true="<<fluxDensity_true<<", lgFlux_true="<<lgFlux_true<<", AvgBkg="<<AvgBkg<<", AvgRMS="<<AvgRMS<<", Z_true="<<Z_true);
		
		int gBin= NTrueSourceHisto_compact->FindBin(lgFlux_true);
		if(NTrueSourceHisto_compact->IsBinUnderflow(gBin) || NTrueSourceHisto_compact->IsBinOverflow(gBin)) continue;
		NTrueSourceHisto_compact->Fill(lgFlux_true,1);
		NTrueSourceVSSignificanceHisto_compact->Fill(Z_true,1);
		nTrueSources++;

		//Skip if not found
		if(found==0) continue;			

		//Fill rec info
		NRecSourceHisto_compact->Fill(lgFlux_true,1);
		NRecSourceVSSignificanceHisto_compact->Fill(Z_true,1);

		double xOffset_true= X0_sweighted_rec-X0_true;
		double yOffset_true= Y0_sweighted_rec-Y0_true;

		double xOffset= X0_sweighted_rec-X0_sweighted;
		double yOffset= Y0_sweighted_rec-Y0_sweighted;
		double Speak= Smax_rec;
		double fluxDensity= S_rec/beamArea_rec;
		double SpeakPull= Speak/S_true-1;
		double FluxPull= fluxDensity/fluxDensity_true-1;
		double FluxPull_peak= Speak/fluxDensity_true-1;
			
		FluxList_compact[gBin-1].push_back(lgFlux_true);
		xPosPullList_compact_centroid[gBin-1].push_back(xOffset);
		yPosPullList_compact_centroid[gBin-1].push_back(yOffset);	
		FluxDensityPullList_compact_peak[gBin-1].push_back(FluxPull_peak);					

		if(HasFitInfo){
			
			//Get fit info
			xOffset_true= X0_fit-X0_true;
			yOffset_true= Y0_fit-Y0_true;
			xOffset= X0_fit-X0_sweighted;
			yOffset= Y0_fit-Y0_sweighted;

			Speak= S_fit;
			fluxDensity= fluxDensity_fit/beamArea_rec;	
			//fluxDensity= fluxDensity_fit;	
			SpeakPull= Speak/S_true-1;
			FluxPull= fluxDensity/fluxDensity_true-1;
				
			//Fill histos
			NRecSourceHisto_compact_fit->Fill(lgFlux_true,1);
			NRecSourceVSSignificanceHisto_compact_fit->Fill(Z_true,1);
		
			//Fill vector fit info
			FluxList_compact_fit[gBin-1].push_back(lgFlux_true);
			xPosPullList_compact_fit[gBin-1].push_back(xOffset);
			yPosPullList_compact_fit[gBin-1].push_back(yOffset);
			xPosTruePullList_compact_fit[gBin-1].push_back(xOffset_true);
			yPosTruePullList_compact_fit[gBin-1].push_back(yOffset_true);
			FluxDensityPullList_compact_fit[gBin-1].push_back(FluxPull);
							
			//Fill vectors cumulated info
			xPosPullList_compact[gBin-1].push_back(xOffset);
			yPosPullList_compact[gBin-1].push_back(yOffset);
			xPosTruePullList_compact[gBin-1].push_back(xOffset_true);
			yPosTruePullList_compact[gBin-1].push_back(yOffset_true);
			FluxDensityPullList_compact[gBin-1].push_back(FluxPull);

			//Fill fit ellipse pars
			double ellMaj= std::max(sigmaX_fit,sigmaY_fit);
			double ellMin= std::min(sigmaX_fit,sigmaY_fit);
			double ellAngle= MathUtils::Mod(theta_fit,180.);
			if(sigmaY_fit>sigmaX_fit) ellAngle= MathUtils::Mod(theta_fit-90.,180.);
			fitEllipse->SetR1(ellMaj);
			fitEllipse->SetR2(ellMin);
			fitEllipse->SetTheta(ellAngle);
			double ellipseArea= MathUtils::ComputeEllipseArea(fitEllipse);
			double ellipseEccentricity= MathUtils::ComputeEllipseEccentricity(fitEllipse);
			double dTheta= fitEllipse->GetTheta()-beamEllipse->GetTheta();	
			//double dTheta= theta_fit-(bpa);
			FitEllipseVSBeam_compact_real->SetPoint(nFitSources,dTheta,ellipseEccentricity/beamEllipseEccentricity);
			INFO_LOG("FitEllipse (sigmaX/sigmaY/theta="<<sigmaX_fit<<"/"<<sigmaY_fit<<"/"<<theta_fit<<", A="<<ellipseArea<<", E="<<ellipseEccentricity<<", R1="<<fitEllipse->GetR1()<<", R2="<<fitEllipse->GetR2()<<", theta="<<fitEllipse->GetTheta()<<"), A/A_beam="<<ellipseArea/beamEllipseArea<<", E/E_beam="<<ellipseEccentricity/beamEllipseEccentricity<<", dTheta(deg)="<<dTheta);
			nFitSources++;

			//INFO_LOG("S="<<S<<" fluxDensity_true="<<fluxDensity_true<<", lgFlux_true="<<lgFlux_true<<", AvgBkg="<<AvgBkg<<", AvgRMS="<<AvgRMS<<", Z_true="<<Z_true<<", fluxDensity_fit="<<fluxDensity_fit<<", FluxPull="<<FluxPull);
				
		}//close if has fit info
		else{
			//Fill vectors cumulated info
			xPosPullList_compact[gBin-1].push_back(xOffset);
			yPosPullList_compact[gBin-1].push_back(yOffset);
			xPosTruePullList_compact[gBin-1].push_back(xOffset_true);
			yPosTruePullList_compact[gBin-1].push_back(yOffset_true);
			FluxDensityPullList_compact[gBin-1].push_back(FluxPull_peak);

			//INFO_LOG("S="<<S<<" fluxDensity_true="<<fluxDensity_true<<", lgFlux_true="<<lgFlux_true<<", AvgBkg="<<AvgBkg<<", AvgRMS="<<AvgRMS<<", Z_true="<<Z_true<<", Speak="<<Speak<<",  FluxPull="<<FluxPull_peak);
		
		}
			
	}//end loop sources

	INFO_LOG("Selected "<<nTrueSources<<"/"<<SourceMatchInfo->GetEntries()<<" true compact sources for analysis...");

	//============================================================
	//===    COMPUTE COMPACT SOURCE FINDING RELIABILITY
	//============================================================	
	//## Loop over rec sources
	nFitSources= 0;
	for(int i=0;i<RecSourceInfo->GetEntries();i++){
		RecSourceInfo->GetEntry(i);
		
		//Skip if not compact source
		bool isPointSource= (Type_rec==Source::ePointLike);
		if(!isPointSource) continue;

		//Skip if source is at edge?	
		bool atEdgeX= (X0_rec<sourceRegionXMin || X0_rec>sourceRegionXMax);
		bool atEdgeY= (Y0_rec<sourceRegionYMin || Y0_rec>sourceRegionYMax);
		bool atEdge= (atEdgeX || atEdgeY);
		if(excludeSourceAtEdge && atEdge) continue;


		//Fill rec info
		double lgFlux_rec= log10(fluxDensity_rec/beamArea_rec);
		double Z_rec= (fluxDensity_rec/beamArea_rec)/noiseLevel_true;
		//INFO_LOG("fluxDensity_rec="<<fluxDensity_rec<<", lgFlux_rec="<<lgFlux_rec<<", Z_rec="<<Z_rec);
		NRecSourceHisto_reliability_compact->Fill(lgFlux_rec,1);
		NRecSourceVSSignificanceHisto_reliability_compact->Fill(Z_rec,1);

		//Fill true info
		if(nTrueMatchedSources>0){
			NTrueSourceHisto_reliability_compact->Fill(lgFlux_rec,1);	
			NTrueSourceVSSignificanceHisto_reliability_compact->Fill(Z_rec,1);
		}
		else{
			//Fill fit ellipse pars
			double ellMaj= std::max(sigmaX_fit,sigmaY_fit);
			double ellMin= std::min(sigmaX_fit,sigmaY_fit);
			double ellAngle= MathUtils::Mod(theta_fit,180.);
			if(sigmaY_fit>sigmaX_fit) ellAngle= MathUtils::Mod(theta_fit-90.,180.);
			fitEllipse->SetR1(ellMaj);
			fitEllipse->SetR2(ellMin);
			fitEllipse->SetTheta(ellAngle);
			double ellipseArea= MathUtils::ComputeEllipseArea(fitEllipse);
			double ellipseEccentricity= MathUtils::ComputeEllipseEccentricity(fitEllipse);
			double dTheta= fitEllipse->GetTheta()-beamEllipse->GetTheta();	
			//double dTheta= theta_fit-(bpa);
			FitEllipseVSBeam_compact_fake->SetPoint(nFitSources,dTheta,ellipseEccentricity/beamEllipseEccentricity);
			INFO_LOG("FALSE SOURCE: FitEllipse (sigmaX/sigmaY/theta="<<sigmaX_fit<<"/"<<sigmaY_fit<<"/"<<theta_fit<<", A="<<ellipseArea<<", E="<<ellipseEccentricity<<", R1="<<fitEllipse->GetR1()<<", R2="<<fitEllipse->GetR2()<<", theta="<<fitEllipse->GetTheta()<<"), A/A_beam="<<ellipseArea/beamEllipseArea<<", E/E_beam="<<ellipseEccentricity/beamEllipseEccentricity<<", dTheta(deg)="<<dTheta);
			nFitSources++;
		}

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
		double Z_true= fluxDensity_true/(noiseLevel_true*sqrt(nBeams_true));
		//double Z_true= fluxDensity_true/(AvgRMS*sqrt(nBeams_true));
		//INFO_LOG("nBeams_true="<<nBeams_true<<", NPix="<<NPix<<", AvgBkg="<<AvgBkg<<", AvgRMS="<<AvgRMS<<", Z_true="<<Z_true);

		int gBin= NTrueSourceHisto_ext->FindBin(lgFlux_true);
		if(NTrueSourceHisto_ext->IsBinUnderflow(gBin) || NTrueSourceHisto_ext->IsBinOverflow(gBin)) continue;	
		NTrueSourceHisto_ext->Fill(lgFlux_true,1);
		NTrueSourceVSSignificanceHisto_ext->Fill(Z_true,1);
		nTrueSources_ext++;

		//Fill rec info
		if(found==1) {	
			NRecSourceHisto_ext->Fill(lgFlux_true,1);
			NRecSourceVSSignificanceHisto_ext->Fill(Z_true,1);

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
	
		}//close if
						
	}//end loop sources

	INFO_LOG("Selected "<<nTrueSources_ext<<"/"<<ExtSourceMatchInfo->GetEntries()<<" true extended sources for analysis...");

	//============================================================
	//===    COMPUTE EXTENDED SOURCE FINDING RELIABILITY
	//============================================================	
	//## Loop over rec sources
	for(int i=0;i<RecExtSourceInfo->GetEntries();i++){
		RecExtSourceInfo->GetEntry(i);
		
		//Fill rec info
		double S_rec= fluxDensity_rec;
		fluxDensity_rec= S_rec/beamArea_rec;
		double lgFlux_rec= log10(fluxDensity_rec);
		double nBeams_rec= (double)(NPix_rec)/beamArea_rec;
		double Z_rec= fluxDensity_rec/(noiseLevel_true*sqrt(nBeams_rec));
		NRecSourceHisto_reliability_ext->Fill(lgFlux_rec,1);
		NRecSourceVSSignificanceHisto_reliability_ext->Fill(Z_rec,1);
		

		//Fill true info
		if(nTrueMatchedSources>0){
			NTrueSourceHisto_reliability_ext->Fill(lgFlux_rec,1);
			NTrueSourceVSSignificanceHisto_reliability_ext->Fill(Z_rec,1);
		}

	}//end loop rec sources

	return 0;

}//close AnalyzeData()

double FitPosErrFcn(double* x,double* par){

	double A= pow(10,x[0]);
	double scaleFactor= par[0];
	double bmaj= par[1];
	double bmin= par[2];
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

