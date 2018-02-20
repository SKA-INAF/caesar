#include <TROOT.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLine.h>

#include <Logger.h>
#include <Source.h>
#include <MathUtils.h>
#include <AstroUtils.h>
#include <StatsUtils.h>

using namespace Caesar;

//Vars
std::vector<double> LgFluxBins= {
	-6,-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.5,0,1,1.5
};

std::map<int,int> typeToHistoIdMap;


int DrawImagingPerformances(std::string fileName){

	//## Set logging level
	LoggerManager::Instance().CreateConsoleLogger("INFO","logger","System.out");
	
	//## Read file
	TFile* inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile){
		ERROR_LOG("Failed to read file "<<fileName<<"!");
		return -1;		
	}

	//## Get tree
	TTree* data= (TTree*)inputFile->Get("data");
	if(!data){
		ERROR_LOG("Failed to read tree from file "<<fileName<<"!");
		return -1;
	}
	
	//# Set tree branches
	std::string name;
	int type;
	int simtype;
	double simmaxscale;
	double S;
	double S_true;
	double S_rec;
	double BkgAvg;
	double S_bkg;
	double S_bkg_mean;
	double beamArea;
	double Npix_true;
	double Npix_rec;

	//data->SetBranchAddress("name",&name);
	data->SetBranchAddress("type",&type);
	data->SetBranchAddress("simtype",&simtype);
	data->SetBranchAddress("simmaxscale",&simmaxscale);
	data->SetBranchAddress("S",&S);
	data->SetBranchAddress("S_true",&S_true);
	data->SetBranchAddress("S_rec",&S_rec);
	data->SetBranchAddress("beamArea",&beamArea);
	data->SetBranchAddress("BkgAvg",&BkgAvg);	
	data->SetBranchAddress("S_bkg",&S_bkg);
	data->SetBranchAddress("S_bkg_mean",&S_bkg_mean);
	data->SetBranchAddress("Npix_true",&Npix_true);
	data->SetBranchAddress("Npix_rec",&Npix_rec);

	//## Define histos
	int nBins= (int)(LgFluxBins.size()-1);
	INFO_LOG("nBins="<<nBins);
	TH1D* dummyHisto= new TH1D("dummyHisto","dummyHisto",nBins,LgFluxBins.data());
	
	TGraph* fluxPoints_compact= new TGraph;
	TGraphAsymmErrors* fluxAccuracy_compact= new TGraphAsymmErrors;
	int nPoints_compact= 0;

	TGraph* fluxPoints_extended= new TGraph;
	TGraphAsymmErrors* fluxAccuracy_extended= new TGraphAsymmErrors;
	int nPoints_extended= 0;

	TGraph* fluxPoints_compactextended= new TGraph;
	TGraphAsymmErrors* fluxAccuracy_compactextended= new TGraphAsymmErrors;
	int nPoints_compactextended= 0;
		
	//## Init data vector
	std::vector<std::vector<double>> fluxList_compact;
	std::vector<std::vector<double>> fluxList_extended;
	std::vector<std::vector<double>> fluxList_compactextended;
	std::vector<std::vector<double>> fluxRatioList_compact;
	std::vector<std::vector<double>> fluxRatioList_extended;
	std::vector<std::vector<double>> fluxRatioList_compactextended;
	for(int i=0;i<nBins;i++){
		fluxList_compact.push_back( std::vector<double>() );
		fluxList_extended.push_back( std::vector<double>() );
		fluxList_compactextended.push_back( std::vector<double>() );
		fluxRatioList_compact.push_back( std::vector<double>() );
		fluxRatioList_extended.push_back( std::vector<double>() );
		fluxRatioList_compactextended.push_back( std::vector<double>() );
	}

	//## Loop over data and fill histos
	INFO_LOG("#"<<data->GetEntries()<<" sources to be read...");
	for(int i=0;i<data->GetEntries();i++){
		data->GetEntry(i);

		//Get info
		double fluxDensity_true= S/beamArea;
		double lgFlux_true= log10(fluxDensity_true);
		
		double bkg= 0;
		if(type==Source::eCompact || type==Source::ePointLike) {
			bkg= S_bkg;
		}
		else if(type==Source::eExtended){
			bkg= S_bkg;
		}
		else if(type==Source::eCompactPlusExtended){
			bkg= Npix_rec*BkgAvg;
		}	
		else{
			bkg= 0;
		}
			
		//double fluxRatio= (S_rec-bkg)/S;
		double fluxRatio= (S_rec-bkg-S)/S;				
		int gBin= dummyHisto->FindBin(lgFlux_true);
		if(dummyHisto->IsBinOverflow(gBin) || dummyHisto->IsBinUnderflow(gBin)){
			WARN_LOG("Overflow/underflow flux bin (S_true="<<S_true<<", S="<<S<<", S_rec="<<S_rec<<", fluxDensity_true="<<fluxDensity_true<<"), skip entry...");
			continue;
		}

		//Fill histo/graph
		if(type==Source::eCompact || type==Source::ePointLike){
			fluxRatioList_compact[gBin-1].push_back(fluxRatio);
			fluxList_compact[gBin-1].push_back(lgFlux_true);
			fluxPoints_compact->SetPoint(nPoints_compact,lgFlux_true,fluxRatio);	
			nPoints_compact++;
		}	
		else if(type==Source::eExtended){
			fluxRatioList_extended[gBin-1].push_back(fluxRatio);
			fluxList_extended[gBin-1].push_back(lgFlux_true);
			fluxPoints_extended->SetPoint(nPoints_extended,lgFlux_true,fluxRatio);	
			nPoints_extended++;
		}
		else if(type==Source::eCompactPlusExtended){
			fluxRatioList_compactextended[gBin-1].push_back(fluxRatio);
			fluxList_compactextended[gBin-1].push_back(lgFlux_true);
			fluxPoints_compactextended->SetPoint(nPoints_compactextended,lgFlux_true,fluxRatio);	
			nPoints_compactextended++;
		}
		
	}//end loop sources

	//===============================================
	//==          FILL STATS HISTO
	//===============================================	
	//Compute data list stats
	nPoints_compact= 0;
	nPoints_extended= 0;
	nPoints_compactextended= 0;
	int nMinEntriesToPlot= 10;

	for(int i=0;i<nBins;i++){
		//double x= LgFluxBins[i] + 0.5*(LgFluxBins[i+1]-LgFluxBins[i]);

		//Compute flux accuracy
		if(fluxList_compact[i].size()>=nMinEntriesToPlot){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(fluxRatioList_compact[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(fluxList_compact[i]);
			double xmin= *(std::min_element(fluxList_compact[i].begin(),fluxList_compact[i].end()));
			double xmax= *(std::max_element(fluxList_compact[i].begin(),fluxList_compact[i].end()));
			double ymin= *(std::min_element(fluxRatioList_compact[i].begin(),fluxRatioList_compact[i].end()));
			double ymax= *(std::max_element(fluxRatioList_compact[i].begin(),fluxRatioList_compact[i].end()));
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			fluxAccuracy_compact->SetPoint(nPoints_compact,x,y);
			fluxAccuracy_compact->SetPointError(nPoints_compact,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_compact++;
			INFO_LOG("Compact flux bin "<<i<<": x="<<x<<", y="<<y<<" (min/max="<<ymin<<"/"<<ymax<<", N="<<fluxList_compact[i].size()<<")");
		}
		if(fluxList_extended[i].size()>=nMinEntriesToPlot){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(fluxRatioList_extended[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(fluxList_extended[i]);
			double xmin= *(std::min_element(fluxList_extended[i].begin(),fluxList_extended[i].end()));
			double xmax= *(std::max_element(fluxList_extended[i].begin(),fluxList_extended[i].end()));
			double ymin= *(std::min_element(fluxRatioList_extended[i].begin(),fluxRatioList_extended[i].end()));
			double ymax= *(std::max_element(fluxRatioList_extended[i].begin(),fluxRatioList_extended[i].end()));
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			fluxAccuracy_extended->SetPoint(nPoints_extended,x,y);
			fluxAccuracy_extended->SetPointError(nPoints_extended,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_extended++;
			INFO_LOG("Extended flux bin "<<i<<": x="<<x<<", y="<<y<<" (min/max="<<ymin<<"/"<<ymax<<", N="<<fluxList_extended[i].size()<<")");
		}
		if(fluxList_compactextended[i].size()>=nMinEntriesToPlot){
			Caesar::BoxStats<double> stats= StatsUtils::ComputeBoxStats(fluxRatioList_compactextended[i]);
			Caesar::BoxStats<double> stats_x= StatsUtils::ComputeBoxStats(fluxList_compactextended[i]);
			double xmin= *(std::min_element(fluxList_compactextended[i].begin(),fluxList_compactextended[i].end()));
			double xmax= *(std::max_element(fluxList_compactextended[i].begin(),fluxList_compactextended[i].end()));
			double ymin= *(std::min_element(fluxRatioList_compactextended[i].begin(),fluxRatioList_compactextended[i].end()));
			double ymax= *(std::max_element(fluxRatioList_compactextended[i].begin(),fluxRatioList_compactextended[i].end()));
			double y= stats.median;
			double yerr_low= y-stats.Q1;
			double yerr_up= stats.Q3-y;
			double x= stats_x.median;
			double xerr_low= x-stats_x.Q1;
			double xerr_up= stats_x.Q3-x;
			fluxAccuracy_compactextended->SetPoint(nPoints_compactextended,x,y);
			fluxAccuracy_compactextended->SetPointError(nPoints_compactextended,xerr_low,xerr_up,yerr_low,yerr_up);
			nPoints_compactextended++;
			INFO_LOG("Compact+extended flux bin "<<i<<": x="<<x<<", y="<<y<<" (min/max="<<ymin<<"/"<<ymax<<", N="<<fluxList_compactextended[i].size()<<")");
		}

	}//end loop bins

	//===============================================
	//==          DRAW PLOTS
	//===============================================	
	gROOT->SetStyle("myStyle2");
	double fluxRatio_min= -10;
	double fluxRatio_max= 10;

	//Reference line
	TLine* refLine_fluxAccuracy= new TLine(LgFluxBins[0],0,LgFluxBins[LgFluxBins.size()-1],0);
	refLine_fluxAccuracy->SetLineColor(kBlack);
	refLine_fluxAccuracy->SetLineStyle(kDashed);

	TCanvas* FluxAccuracyPlot= new TCanvas("FluxAccuracyPlot","FluxAccuracyPlot");
	FluxAccuracyPlot->cd();

	TH2D* FluxAccuracyPlotBkg= new TH2D("FluxAccuracyPlotBkg","",100,LgFluxBins[0]-0.5,LgFluxBins[LgFluxBins.size()-1]+0.5,100,fluxRatio_min,fluxRatio_max);
	FluxAccuracyPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	FluxAccuracyPlotBkg->GetYaxis()->SetTitle("S_{rec}/S_{true}-1");
	FluxAccuracyPlotBkg->SetStats(0);
	FluxAccuracyPlotBkg->Draw();

	fluxPoints_compact->SetMarkerSize(1);
	fluxPoints_compact->SetMarkerStyle(1);
	fluxPoints_compact->SetMarkerColor(kRed-10);
	fluxPoints_compact->SetLineColor(kRed-10);
	fluxPoints_compact->Draw("P same");
	
	fluxAccuracy_compact->SetMarkerSize(1.3);
	fluxAccuracy_compact->SetMarkerStyle(8);
	fluxAccuracy_compact->SetMarkerColor(kRed);
	fluxAccuracy_compact->SetLineColor(kRed);
	fluxAccuracy_compact->Draw("ep same");

	refLine_fluxAccuracy->Draw("same");

	//###  EXTENDED ###
	TCanvas* FluxAccuracyPlot_extended= new TCanvas("FluxAccuracyPlot_extended","FluxAccuracyPlot_extended");
	FluxAccuracyPlot_extended->cd();

	TH2D* FluxAccuracyPlotBkg_extended= new TH2D("FluxAccuracyPlotBkg_extended","",100,LgFluxBins[0]-0.5,LgFluxBins[LgFluxBins.size()-1]+0.5,100,fluxRatio_min,fluxRatio_max);
	FluxAccuracyPlotBkg_extended->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	FluxAccuracyPlotBkg_extended->GetYaxis()->SetTitle("S_{rec}/S_{true}-1");
	FluxAccuracyPlotBkg_extended->SetStats(0);
	FluxAccuracyPlotBkg_extended->Draw();

	fluxPoints_extended->SetMarkerSize(1);
	fluxPoints_extended->SetMarkerStyle(1);
	fluxPoints_extended->SetMarkerColor(kGreen+1);
	fluxPoints_extended->SetLineColor(kGreen+1);
	fluxPoints_extended->Draw("P same");

	fluxAccuracy_extended->SetMarkerSize(1.3);
	fluxAccuracy_extended->SetMarkerStyle(8);
	fluxAccuracy_extended->SetMarkerColor(kGreen+1);
	fluxAccuracy_extended->SetLineColor(kGreen+1);
	fluxAccuracy_extended->Draw("ep same");

	refLine_fluxAccuracy->Draw("same");

	//###  EXTENDED + COMPACT ###
	TCanvas* FluxAccuracyPlot_compactextended= new TCanvas("FluxAccuracyPlot_compactextended","FluxAccuracyPlot_compactextended");
	FluxAccuracyPlot_compactextended->cd();

	TH2D* FluxAccuracyPlotBkg_compactextended= new TH2D("FluxAccuracyPlotBkg_compactextended","",100,LgFluxBins[0]-0.5,LgFluxBins[LgFluxBins.size()-1]+0.5,100,fluxRatio_min,fluxRatio_max);
	FluxAccuracyPlotBkg_compactextended->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	FluxAccuracyPlotBkg_compactextended->GetYaxis()->SetTitle("S_{rec}/S_{true}-1");
	FluxAccuracyPlotBkg_compactextended->SetStats(0);
	FluxAccuracyPlotBkg_compactextended->Draw();

	fluxPoints_compactextended->SetMarkerSize(1);
	fluxPoints_compactextended->SetMarkerStyle(1);
	fluxPoints_compactextended->SetMarkerColor(kBlue);
	fluxPoints_compactextended->SetLineColor(kBlue);
	fluxPoints_compactextended->Draw("P same");

	fluxAccuracy_compactextended->SetMarkerSize(1.3);
	fluxAccuracy_compactextended->SetMarkerStyle(8);
	fluxAccuracy_compactextended->SetMarkerColor(kBlue);
	fluxAccuracy_compactextended->SetLineColor(kBlue);
	fluxAccuracy_compactextended->Draw("ep same");

	refLine_fluxAccuracy->Draw("same");

	//### All sources
	TCanvas* FluxAccuracyPlot_all= new TCanvas("FluxAccuracyPlot_all","FluxAccuracyPlot_all");
	FluxAccuracyPlot_all->cd();

	TH2D* FluxAccuracyPlotBkg_all= new TH2D("FluxAccuracyPlotBkg_all","",100,LgFluxBins[0]-0.5,LgFluxBins[LgFluxBins.size()-1]+0.5,100,fluxRatio_min,fluxRatio_max);
	FluxAccuracyPlotBkg_all->GetXaxis()->SetTitle("log_{10}(S_{true}/Jy)");
	FluxAccuracyPlotBkg_all->GetYaxis()->SetTitle("S_{rec}/S_{true}");
	FluxAccuracyPlotBkg_all->SetStats(0);
	FluxAccuracyPlotBkg_all->Draw();

	
	fluxPoints_compact->SetMarkerSize(1);
	fluxPoints_compact->SetMarkerStyle(1);
	fluxPoints_compact->SetMarkerColor(kRed-10);
	fluxPoints_compact->SetLineColor(kRed-10);
	//fluxPoints_compact->Draw("P same");

	fluxPoints_extended->SetMarkerSize(1);
	fluxPoints_extended->SetMarkerStyle(1);
	fluxPoints_extended->SetMarkerColor(kGreen+1);
	fluxPoints_extended->SetLineColor(kGreen+1);
	//fluxPoints_extended->Draw("P same");

	fluxPoints_compactextended->SetMarkerSize(1);
	fluxPoints_compactextended->SetMarkerStyle(1);
	fluxPoints_compactextended->SetMarkerColor(kBlue);
	fluxPoints_compactextended->SetLineColor(kBlue);
	//fluxPoints_compactextended->Draw("P same");
	

	fluxAccuracy_compact->SetMarkerSize(1.3);
	fluxAccuracy_compact->SetMarkerStyle(8);
	fluxAccuracy_compact->SetMarkerColor(kRed);
	fluxAccuracy_compact->SetLineColor(kRed);
	fluxAccuracy_compact->Draw("ep same");


	fluxAccuracy_extended->SetMarkerSize(1.3);
	fluxAccuracy_extended->SetMarkerStyle(21);
	fluxAccuracy_extended->SetMarkerColor(kGreen+1);
	fluxAccuracy_extended->SetLineColor(kGreen+1);
	fluxAccuracy_extended->Draw("ep same");

	
	fluxAccuracy_compactextended->SetMarkerSize(1.3);
	fluxAccuracy_compactextended->SetMarkerStyle(22);
	fluxAccuracy_compactextended->SetMarkerColor(kBlue);
	fluxAccuracy_compactextended->SetLineColor(kBlue);
	fluxAccuracy_compactextended->Draw("ep same");

	
	refLine_fluxAccuracy->Draw("same");
	

	return 0;

}//close DrawImagingPerformances()
