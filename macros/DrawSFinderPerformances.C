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
double Bmaj= 9.8;//arcsec
double Bmin= 5.8;//arcsec
double pixSize= 1;//in arcsec
double beamArea= 0;

//### Selection cuts
bool selectResolvedTrueSources= false;
double mutualTrueSourceDistThr= 20;//in arcsec 
bool selectFitComponents= true;
int maxFitComponentsThr= 1;

//Draw info
//std::vector<double> Zbins= {0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};
std::vector<double> Zbins= {0.1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100};

std::vector<double> LgFluxBins= {
	-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.5,0,1,2
};

//- Compact source histos/graphs
TH1D* NTrueSourceHisto_compact= 0;
TH1D* NRecSourceHisto_compact= 0;
TH1D* NFitSourceHisto_compact= 0;
TH1D* NTrueSourceVSSignificanceHisto_compact= 0;
TH1D* NRecSourceVSSignificanceHisto_compact= 0;
TH1D* NFitSourceVSSignificanceHisto_compact= 0;
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

TGraphAsymmErrors* xPosAccuracyGraph_compact= 0;
TGraphAsymmErrors* yPosAccuracyGraph_compact= 0;
TGraphAsymmErrors* SpeakAccuracyGraph_compact= 0;
TGraphAsymmErrors* FluxDensityAccuracyGraph_compact= 0;

std::vector< std::vector<double> > FluxList_compact;
std::vector< std::vector<double> > xPosPullList_compact;
std::vector< std::vector<double> > yPosPullList_compact;
std::vector< std::vector<double> > SpeakPullList_compact;
std::vector< std::vector<double> > FluxDensityPullList_compact;

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
	beamArea= AstroUtils::GetBeamAreaInPixels(Bmaj,Bmin,pixSize,pixSize);
	
	int nBins= (int)(LgFluxBins.size()-1);
	int nBins_Z= (int)(Zbins.size()-1);
	
	//## Init data vector
	for(int i=0;i<nBins;i++){
		FluxList_compact.push_back( std::vector<double>() );
		xPosPullList_compact.push_back( std::vector<double>() );
		yPosPullList_compact.push_back( std::vector<double>() );
		SpeakPullList_compact.push_back( std::vector<double>() );
		FluxDensityPullList_compact.push_back( std::vector<double>() );
		
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
	
	histoName= "NFitSourceHisto_compact";
	NFitSourceHisto_compact= new TH1D(histoName,"Fitted compact source distribution",nBins,LgFluxBins.data());
	NFitSourceHisto_compact->Sumw2();

	histoName= "NTrueSourceVSSignificanceHisto_compact";
	NTrueSourceVSSignificanceHisto_compact= new TH1D(histoName,"True compact source distribution",nBins,Zbins.data());
	NTrueSourceVSSignificanceHisto_compact->Sumw2();
	
	histoName= "NRecSourceVSSignificanceHisto_compact";
	NRecSourceVSSignificanceHisto_compact= new TH1D(histoName,"Rec compact source distribution",nBins,Zbins.data());
	NRecSourceVSSignificanceHisto_compact->Sumw2();
	
	histoName= "NFitSourceVSSignificanceHisto_compact";
	NFitSourceVSSignificanceHisto_compact= new TH1D(histoName,"Fitted compact source distribution",nBins,Zbins.data());
	NFitSourceVSSignificanceHisto_compact->Sumw2();

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

	/*
	//- Positional accuracy histos
	histoName= "xPosAccuracyHisto_compact";
	xPosAccuracyHisto_compact= new TProfile(histoName,"x position accuracy",nBins,LgFluxBins.data(),"");
	xPosAccuracyHisto_compact->Sumw2();
	
	histoName= "yPosAccuracyHisto_compact";
	yPosAccuracyHisto_compact= new TProfile(histoName,"y position accuracy",nBins,LgFluxBins.data(),"");
	yPosAccuracyHisto_compact->Sumw2();
	
	//- Flux accuracy histos
	histoName= "SpeakAccuracyHisto_compact";
	SpeakAccuracyHisto_compact= new TProfile(histoName,"Speak accuracy",nBins,LgFluxBins.data(),"");
	SpeakAccuracyHisto_compact->Sumw2();
	
	histoName= "FluxDensityAccuracyHisto";
	FluxDensityAccuracyHisto= new TProfile(histoName,"Flux density accuracy",nBins,LgFluxBins.data(),"");
	FluxDensityAccuracyHisto->Sumw2();

	histoName= "FluxDensityAccuracyHisto_ext";
	FluxDensityAccuracyHisto_ext= new TProfile(histoName,"Extended source flux density accuracy",nBins,LgFluxBins.data(),"");
	FluxDensityAccuracyHisto_ext->Sumw2();
	
	//- Offset accuracy histos
	histoName= "BkgAccuracyHisto";
	BkgAccuracyHisto= new TProfile(histoName,"Offset accuracy",nBins,LgFluxBins.data(),"");
	BkgAccuracyHisto->Sumw2();
	*/

	//- Pos accuracy graphs
	TString graphName= "xPosAccuracyGraph_compact";
	xPosAccuracyGraph_compact= new TGraphAsymmErrors;
	xPosAccuracyGraph_compact->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_compact";
	yPosAccuracyGraph_compact= new TGraphAsymmErrors;
	yPosAccuracyGraph_compact->SetName(graphName);
		
	graphName= "xPosAccuracyGraph_ext";
	xPosAccuracyGraph_ext= new TGraphAsymmErrors;
	xPosAccuracyGraph_ext->SetName(graphName);
	
	graphName= "yPosAccuracyGraph_ext";
	yPosAccuracyGraph_ext= new TGraphAsymmErrors;
	yPosAccuracyGraph_ext->SetName(graphName);

	//- Flux accuracy graphs
	graphName= "SpeakAccuracyGraph_compact";
	SpeakAccuracyGraph_compact= new TGraphAsymmErrors;
	SpeakAccuracyGraph_compact->SetName(graphName);
	
	graphName= "FluxDensityAccuracyGraph_compact";
	FluxDensityAccuracyGraph_compact= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_compact->SetName(graphName);

	graphName= "FluxDensityAccuracyGraph_ext";
	FluxDensityAccuracyGraph_ext= new TGraphAsymmErrors;
	FluxDensityAccuracyGraph_ext->SetName(graphName);
	
	//graphName= "BkgAccuracyGraph";
	//BkgAccuracyGraph= new TGraphAsymmErrors;
	//BkgAccuracyGraph->SetName(graphName);
		
	
}//close Init()



void ComputeAnalysisHistos(){

	//Compute detection efficiency
	Efficiency_compact = new TEfficiency(*NRecSourceHisto_compact,*NTrueSourceHisto_compact); 
	Efficiency_compact->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Efficiency_compact->SetConfidenceLevel(0.68);
	
	Efficiency_compact_fit = new TEfficiency(*NFitSourceHisto_compact,*NTrueSourceHisto_compact); 
	Efficiency_compact_fit->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	Efficiency_compact_fit->SetConfidenceLevel(0.68);

	EfficiencyVSSignificance_compact = new TEfficiency(*NRecSourceVSSignificanceHisto_compact,*NTrueSourceVSSignificanceHisto_compact); 
	EfficiencyVSSignificance_compact->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	EfficiencyVSSignificance_compact->SetConfidenceLevel(0.68);
	
	EfficiencyVSSignificance_compact_fit = new TEfficiency(*NFitSourceVSSignificanceHisto_compact,*NTrueSourceVSSignificanceHisto_compact); 
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
	int nPoints_SpeakPull= 0;
	int nPoints_FluxPull= 0;
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

		//Compute Speak accuracy
		if(SpeakPullList_compact[i].size()>=nMinPointsToDraw){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(SpeakPullList_compact[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(FluxList_compact[i]);
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			SpeakAccuracyGraph_compact->SetPoint(nPoints_SpeakPull,x,y);
			SpeakAccuracyGraph_compact->SetPointError(nPoints_SpeakPull,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_SpeakPull++;
		}	

		//Compute Flux accuracy
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

	//===============================================
	//==          DRAW HISTO
	//===============================================
	TCanvas* TrueSourceDistributionPlot= new TCanvas("TrueSourceDistributionPlot","TrueSourceDistributionPlot");
	TrueSourceDistributionPlot->cd();
	
	NTrueSourceHisto_compact->Draw("ep");

	TCanvas* FittedSourceDistributionPlot= new TCanvas("FittedSourceDistributionPlot","FittedSourceDistributionPlot");
	FittedSourceDistributionPlot->cd();
	
	NFitSourceHisto_compact->Draw("ep");

	//===============================================
	//==          DRAW COMPLETENESS
	//===============================================
	double xMin_draw= LgFluxBins[0];
	double xMax_draw= LgFluxBins[LgFluxBins.size()-1];
	double ZMin_draw= Zbins[0];
	double ZMax_draw= Zbins[Zbins.size()-1];

	double zThr= 5;//(S-bkg)/noise= S/noise-bkg/noise= 5
	double SNThr= zThr + bkgLevel_true/noiseLevel_true;
	TLine* detectionThrLine= new TLine(SNThr,0,SNThr,1);
	detectionThrLine->SetLineColor(kBlack);
	detectionThrLine->SetLineStyle(kDashed);
	
	//- Efficiency plot compact source
	TCanvas* EffPlot= new TCanvas("EffPlot","EffPlot");
	EffPlot->cd();

	TH2D* EffPlotBkg= new TH2D("EffPlotBkg","",100,xMin_draw,xMax_draw,100,0,1.2);
	EffPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	EffPlotBkg->GetYaxis()->SetTitle("Completeness");
	EffPlotBkg->Draw();

	Efficiency_compact->SetMarkerStyle(8);
	Efficiency_compact->SetMarkerColor(kBlack);
	Efficiency_compact->SetLineColor(kBlack);
	Efficiency_compact->Draw("p same");
	
	detectionThrLine->Draw("l");

	TCanvas* EffVSSignificancePlot= new TCanvas("EffVSSignificancePlot","EffVSSignificancePlot");
	EffVSSignificancePlot->cd();

	TH2D* EffVSSignificancePlotBkg= new TH2D("EffVSSignificancePlotBkg","",100,ZMin_draw,ZMax_draw,100,0,1.2);
	EffVSSignificancePlotBkg->GetXaxis()->SetTitle("S/N");
	EffVSSignificancePlotBkg->GetYaxis()->SetTitle("Completeness");
	EffVSSignificancePlotBkg->Draw();

	EfficiencyVSSignificance_compact->SetMarkerStyle(8);
	EfficiencyVSSignificance_compact->SetMarkerColor(kBlack);
	EfficiencyVSSignificance_compact->SetLineColor(kBlack);
	EfficiencyVSSignificance_compact->Draw("p same");

	//- Efficiency plot compact source (fitted data)
	TCanvas* EffPlot_fit= new TCanvas("EffPlot_fit","EffPlot_fit");
	EffPlot_fit->cd();

	TH2D* EffPlotBkg_fit= new TH2D("EffPlotBkg_fit","",100,xMin_draw,xMax_draw,100,0,1.2);
	EffPlotBkg_fit->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	EffPlotBkg_fit->GetYaxis()->SetTitle("Completeness");
	EffPlotBkg_fit->Draw();

	Efficiency_compact_fit->SetMarkerStyle(8);
	Efficiency_compact_fit->SetMarkerColor(kBlack);
	Efficiency_compact_fit->SetLineColor(kBlack);
	Efficiency_compact_fit->Draw("p same");
	
	detectionThrLine->Draw("l");

	TCanvas* EffVSSignificancePlot_fit= new TCanvas("EffVSSignificancePlot_fit","EffVSSignificancePlot_fit");
	EffVSSignificancePlot_fit->cd();

	TH2D* EffVSSignificancePlotBkg_fit= new TH2D("EffVSSignificancePlotBkg_fit","",100,ZMin_draw,ZMax_draw,100,0,1.2);
	EffVSSignificancePlotBkg_fit->GetXaxis()->SetTitle("S/N");
	EffVSSignificancePlotBkg_fit->GetYaxis()->SetTitle("Completeness");
	EffVSSignificancePlotBkg_fit->Draw();

	EfficiencyVSSignificance_compact_fit->SetMarkerStyle(8);
	EfficiencyVSSignificance_compact_fit->SetMarkerColor(kBlack);
	EfficiencyVSSignificance_compact_fit->SetLineColor(kBlack);
	EfficiencyVSSignificance_compact_fit->Draw("p same");

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
	
	detectionThrLine->Draw("l");


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

	Reliability_compact->SetMarkerStyle(8);
	Reliability_compact->SetMarkerColor(kBlack);
	Reliability_compact->SetLineColor(kBlack);
	Reliability_compact->Draw("p same");
	
	detectionThrLine->Draw("l");

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
	//- X pos accuracy
	TCanvas* xPosAccuracyPlot= new TCanvas("xPosAccuracyPlot","xPosAccuracyPlot");
	xPosAccuracyPlot->cd();

	TH2D* xPosAccuracyPlotBkg= new TH2D("xPosAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-5,5);
	xPosAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	xPosAccuracyPlotBkg->GetYaxis()->SetTitle("<X0_{rec}-X0_{true}> ('')");
	xPosAccuracyPlotBkg->Draw();

	xPosAccuracyGraph_compact->SetMarkerStyle(8);
	xPosAccuracyGraph_compact->SetMarkerColor(kBlack);
	xPosAccuracyGraph_compact->SetLineColor(kBlack);
	xPosAccuracyGraph_compact->Draw("ep same");
	
	TLine* refLine_xposAccuracy= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_xposAccuracy->SetLineColor(kBlack);
	refLine_xposAccuracy->SetLineStyle(kDashed);
	refLine_xposAccuracy->Draw("same");

	TLine* detectionThrLine_xposAccuracy= new TLine(SNThr,-5,SNThr,5);
	detectionThrLine_xposAccuracy->SetLineColor(kBlack);
	detectionThrLine_xposAccuracy->SetLineStyle(kDashed);
	detectionThrLine_xposAccuracy->Draw("same");

	double sigmaToFirstQuantile= ROOT::Math::gaussian_quantile(0.75,1);
	TF1* expPosXSigmaFcn_plus= new TF1("expPosXSigmaFcn_plus",GausPosSigmaFcn,xMin_draw,xMax_draw,4);
	expPosXSigmaFcn_plus->SetParameters(sigmaToFirstQuantile,pixSize,Bmaj,Bmin);
	expPosXSigmaFcn_plus->SetLineColor(kBlack);
	expPosXSigmaFcn_plus->SetLineStyle(9);
	expPosXSigmaFcn_plus->Draw("l same");

	TF1* expPosXSigmaFcn_minus= new TF1("expPosXSigmaFcn_minus",GausPosSigmaFcn,xMin_draw,xMax_draw,4);
	expPosXSigmaFcn_minus->SetParameters(-sigmaToFirstQuantile,pixSize,Bmaj,Bmin);
	expPosXSigmaFcn_minus->SetLineColor(kBlack);
	expPosXSigmaFcn_minus->SetLineStyle(9);
	expPosXSigmaFcn_minus->Draw("l same");

	//- Y pos accuracy
	TCanvas* yPosAccuracyPlot= new TCanvas("yPosAccuracyPlot","yPosAccuracyPlot");
	yPosAccuracyPlot->cd();

	TH2D* yPosAccuracyPlotBkg= new TH2D("yPosAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-5,5);
	yPosAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	yPosAccuracyPlotBkg->GetYaxis()->SetTitle("<Y0_{rec}-Y0_{true}> ('')");
	yPosAccuracyPlotBkg->Draw();

	yPosAccuracyGraph_compact->SetMarkerStyle(8);
	yPosAccuracyGraph_compact->SetMarkerColor(kBlack);
	yPosAccuracyGraph_compact->SetLineColor(kBlack);
	yPosAccuracyGraph_compact->Draw("ep same");

	TLine* refLine_yposAccuracy= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_yposAccuracy->SetLineColor(kBlack);
	refLine_yposAccuracy->SetLineStyle(kDashed);
	refLine_yposAccuracy->Draw("same");

	TLine* detectionThrLine_yposAccuracy= new TLine(SNThr,-5,SNThr,5);
	detectionThrLine_yposAccuracy->SetLineColor(kBlack);
	detectionThrLine_yposAccuracy->SetLineStyle(kDashed);
	detectionThrLine_yposAccuracy->Draw("same");

	TF1* expPosYSigmaFcn_plus= new TF1("expPosYSigmaFcn_plus",GausPosSigmaFcn,xMin_draw,xMax_draw,4);
	expPosYSigmaFcn_plus->SetParameters(sigmaToFirstQuantile,pixSize,Bmaj,Bmin);
	expPosYSigmaFcn_plus->SetLineColor(kBlack);
	expPosYSigmaFcn_plus->SetLineStyle(9);
	expPosYSigmaFcn_plus->Draw("l same");

	TF1* expPosYSigmaFcn_minus= new TF1("expPosYSigmaFcn_minus",GausPosSigmaFcn,xMin_draw,xMax_draw,4);
	expPosYSigmaFcn_minus->SetParameters(-sigmaToFirstQuantile,pixSize,Bmaj,Bmin);
	expPosYSigmaFcn_minus->SetLineColor(kBlack);
	expPosYSigmaFcn_minus->SetLineStyle(9);
	expPosYSigmaFcn_minus->Draw("l same");


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
	//- Peak accuracy
	TCanvas* SPeakAccuracyPlot= new TCanvas("SPeakAccuracyPlot","SPeakAccuracyPlot");
	SPeakAccuracyPlot->cd();

	TH2D* SPeakAccuracyPlotBkg= new TH2D("SPeakAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-5,5);
	SPeakAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	SPeakAccuracyPlotBkg->GetYaxis()->SetTitle("<Speak_{rec}/Speak_{true}-1>");
	SPeakAccuracyPlotBkg->Draw();

	SpeakAccuracyGraph_compact->SetMarkerStyle(8);
	SpeakAccuracyGraph_compact->SetMarkerColor(kBlack);
	SpeakAccuracyGraph_compact->SetLineColor(kBlack);
	SpeakAccuracyGraph_compact->Draw("ep same");
	
	TLine* refLine_SpeakAccuracy= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_SpeakAccuracy->SetLineColor(kBlack);
	refLine_SpeakAccuracy->SetLineStyle(kDashed);
	refLine_SpeakAccuracy->Draw("same");

	TLine* detectionThrLine_SpeakAccuracy= new TLine(SNThr,-5,SNThr,5);
	detectionThrLine_SpeakAccuracy->SetLineColor(kBlack);
	detectionThrLine_SpeakAccuracy->SetLineStyle(kDashed);
	detectionThrLine_SpeakAccuracy->Draw("same");

	//- FLux accuracy
	TCanvas* FluxAccuracyPlot= new TCanvas("FluxAccuracyPlot","FluxAccuracyPlot");
	FluxAccuracyPlot->cd();

	TH2D* FluxAccuracyPlotBkg= new TH2D("FluxAccuracyPlotBkg","",100,xMin_draw,xMax_draw,100,-5,5);
	FluxAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S/Jy)");
	FluxAccuracyPlotBkg->GetYaxis()->SetTitle("<S_{rec}/S_{true}-1>");
	FluxAccuracyPlotBkg->Draw();

	FluxDensityAccuracyGraph_compact->SetMarkerStyle(8);
	FluxDensityAccuracyGraph_compact->SetMarkerColor(kBlack);
	FluxDensityAccuracyGraph_compact->SetLineColor(kBlack);
	FluxDensityAccuracyGraph_compact->Draw("ep same");

	TLine* refLine_fluxDensityAccuracy= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_fluxDensityAccuracy->SetLineColor(kBlack);
	refLine_fluxDensityAccuracy->SetLineStyle(kDashed);
	refLine_fluxDensityAccuracy->Draw("same");

	TLine* detectionThrLine_fluxDensityAccuracy= new TLine(SNThr,-5,SNThr,5);
	detectionThrLine_fluxDensityAccuracy->SetLineColor(kBlack);
	detectionThrLine_fluxDensityAccuracy->SetLineStyle(kDashed);
	detectionThrLine_fluxDensityAccuracy->Draw("same");

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

	TLine* refLine_fluxDensityAccuracy_ext= new TLine(xMin_draw,0,xMax_draw,0);
	refLine_fluxDensityAccuracy_ext->SetLineColor(kBlack);
	refLine_fluxDensityAccuracy_ext->SetLineStyle(kDashed);
	refLine_fluxDensityAccuracy_ext->Draw("same");

	TLine* detectionThrLine_fluxDensityAccuracy_ext= new TLine(SNThr,-5,SNThr,5);
	detectionThrLine_fluxDensityAccuracy_ext->SetLineColor(kBlack);
	detectionThrLine_fluxDensityAccuracy_ext->SetLineStyle(kDashed);
	detectionThrLine_fluxDensityAccuracy_ext->Draw("same");

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

	//Get point/compact source reliability tree
	TTree* RecSourceInfo= (TTree*)inputFile->Get("RecSourceInfo");
	if(!RecSourceInfo){
		ERROR_LOG("Failed to get RecSourceInfo tree from file "<<filename<<"!");
		return -1;
	}

	RecSourceInfo->SetBranchAddress("name_rec",&Name_rec);
	RecSourceInfo->SetBranchAddress("type_rec",&Type_rec);	
	RecSourceInfo->SetBranchAddress("HasFitInfo",&HasFitInfo);	
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

		INFO_LOG("S="<<S<<" fluxDensity_true="<<fluxDensity_true<<", lgFlux_true="<<lgFlux_true<<", AvgBkg="<<AvgBkg<<", AvgRMS="<<AvgRMS<<", Z_true="<<Z_true);
		
		int gBin= NTrueSourceHisto_compact->FindBin(lgFlux_true);
		if(NTrueSourceHisto_compact->IsBinUnderflow(gBin) || NTrueSourceHisto_compact->IsBinOverflow(gBin)) continue;
		NTrueSourceHisto_compact->Fill(lgFlux_true,1);
		NTrueSourceVSSignificanceHisto_compact->Fill(Z_true,1);
		nTrueSources++;

			
		//Fill rec info
		if(found==1) {	
			NRecSourceHisto_compact->Fill(lgFlux_true,1);
			NRecSourceVSSignificanceHisto_compact->Fill(Z_true,1);

			//double xOffset= fabs(X0_rec-X0_true);
			//double yOffset= fabs(Y0_rec-Y0_true);
			double xOffset= X0_sweighted_rec-X0_true;
			double yOffset= Y0_sweighted_rec-Y0_true;
			double Speak= Smax_rec;
			double fluxDensity= S_rec/beamArea_rec;
			double SpeakPull= Speak/S_true-1;
			double FluxPull= fluxDensity/fluxDensity_true-1;
			
			if(HasFitInfo){

				//Get fit info
				xOffset= X0_fit-X0_true;
				yOffset= Y0_fit-Y0_true;
				Speak= S_fit;
				//fluxDensity= fluxDensity_fit/beamArea_rec;	
				fluxDensity= fluxDensity_fit;	
				SpeakPull= Speak/S_true-1;
				FluxPull= fluxDensity/fluxDensity_true-1;
						
				NFitSourceHisto_compact->Fill(lgFlux_true,1);	
				
				FluxList_compact[gBin-1].push_back(lgFlux_true);
				xPosPullList_compact[gBin-1].push_back(xOffset);
				yPosPullList_compact[gBin-1].push_back(yOffset);
				SpeakPullList_compact[gBin-1].push_back(SpeakPull);
				FluxDensityPullList_compact[gBin-1].push_back(FluxPull);
							
			}//close if has fit info			
						
		}//close if found
		
	}//end loop sources

	INFO_LOG("Selected "<<nTrueSources<<"/"<<SourceMatchInfo->GetEntries()<<" true compact sources for analysis...");

	//============================================================
	//===    COMPUTE SOURCE FINDING RELIABILITY
	//============================================================	
	//## Loop over rec sources
	for(int i=0;i<RecSourceInfo->GetEntries();i++){
		RecSourceInfo->GetEntry(i);
		
		//Fill rec info
		double lgFlux_rec= log10(fluxDensity_rec/beamArea_rec);
		double Z_rec= fluxDensity_rec/noiseLevel_true;
		NRecSourceHisto_reliability_compact->Fill(lgFlux_rec,1);
		NRecSourceVSSignificanceHisto_reliability_compact->Fill(Z_rec,1);

		//Fill true info
		if(nTrueMatchedSources>0){
			NTrueSourceHisto_reliability_compact->Fill(lgFlux_rec);	
			NTrueSourceVSSignificanceHisto_reliability_compact->Fill(Z_rec,1);
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
		INFO_LOG("nBeams_true="<<nBeams_true<<", NPix="<<NPix<<", AvgBkg="<<AvgBkg<<", AvgRMS="<<AvgRMS<<", Z_true="<<Z_true);

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


double GausPosSigmaFcn(double* x,double* par){

	double X= x[0];//SN
	double scaleFactor= par[0];
	double pixelSize= par[1];
	double sigmaX= par[2];
	double sigmaY= par[3];

	double fcn= scaleFactor*sqrt(2*sigmaX/(TMath::Pi()*sigmaY))*pixelSize/X;
	//cout<<"X="<<X<<" fcn="<<fcn<<endl;	
	
	return fcn;

}//close GausPosVarianceFcn()

