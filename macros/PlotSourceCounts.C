
//=================================
//==   FUNCTIONS
//=================================
TH1D* ComputeSourceCounts(std::string filename,double A_sr,bool normalize=false,double beta=2.5,bool logyscale=false);


//=================================
//==   VARIABLES
//=================================
/*
std::vector<double> LgFluxBins= {
-4,-3.75,-3.5,-3.25,-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.75,-0.5,-0.25,0,0.5,1
};
*/
std::vector<double> LgFluxBins= {
-4,-3.8,-3.6,-3.4,-3.2,-3,-2.8,-2.6,-2.4,-2.2,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.5,1,1.25
};


//=================================
//==   MACRO
//=================================

void PlotSourceCounts(std::string filename,double A_sr=1,bool normalize=false,double beta=2.5,bool logyscale=false)
{
	
	//# - Compute source count histo
	TH1D* SourceCounts= ComputeSourceCounts(filename,A_sr,normalize,beta,logyscale);

	
	//# - Draw plots
	if(!logyscale){
		gStyle->SetOptLogx();
	}
	gStyle->SetOptLogy();

	TCanvas* SourceCountsPlot= new TCanvas("SourceCountsPlot","SourceCountsPlot");
	SourceCountsPlot->cd();

	TString yTitle= "dN/dS (Jy^{-1}sr^{-1})";
	if(normalize) yTitle= Form("S^{%1.1f}*dN/dS (Jy^{-%1.1f}sr^{-1})",beta,beta-1);
	TString xTitle= "S/Jy";
	if(logyscale) xTitle= "log_{10}(S/Jy)";

	int nBins= (int)(LgFluxBins.size()-1);
	double xmin= LgFluxBins[0]-0.5;
	double xmax= LgFluxBins[LgFluxBins.size()-1]+0.5;
	double ymin= 1.e-7;
	double ymax= SourceCounts->GetMaximum()*10;

	if(!logyscale){
		xmin= pow(10,LgFluxBins[0]-0.5);
		xmax= pow(10,LgFluxBins[LgFluxBins.size()-1]+0.5);
	}

	
	TH2D* SourceCountsPlotBkg= new TH2D("SourceCountsPlotBkg","",100,xmin,xmax,100,ymin,ymax);
	SourceCountsPlotBkg->SetStats(0);
	SourceCountsPlotBkg->GetXaxis()->SetTitle(xTitle);
	SourceCountsPlotBkg->GetYaxis()->SetTitle(yTitle);
	SourceCountsPlotBkg->Draw();

	SourceCounts->SetMarkerStyle(8);
	SourceCounts->SetMarkerSize(1.3);
	SourceCounts->SetMarkerColor(kBlack);
	SourceCounts->SetLineColor(kBlack);
	SourceCounts->Draw("EPZ same");
	
	
}//close macro



TH1D* ComputeSourceCounts(std::string filename,double A_sr,bool normalize=false,double beta=2.5,bool logyscale=false)
{
	// - Adjust spectral index. If lgS is used as x axis gamma=beta-1
	double gamma= beta;
	if(logyscale) gamma= beta-1;

	//- Init histo
	std::vector<double> FluxBins;
	for(size_t i=0;i<LgFluxBins.size();i++){
		double Flux= pow(10,LgFluxBins[i]);
		FluxBins.push_back(Flux);
	}
	int nBins= (int)(FluxBins.size()-1);

	TH1D* dNdSHisto= 0;
	if(logyscale) dNdSHisto= new TH1D("dNdSHisto","dNdSHisto",nBins,LgFluxBins.data());
	else dNdSHisto= new TH1D("dNdSHisto","dNdSHisto",nBins,FluxBins.data());
	dNdSHisto->Sumw2();

	//- Read input file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()) {
		cerr<<"ERROR: Cannot open file "<<filename<<"!"<<endl;
		return nullptr;
	}
	
	//- Get source tree
	TTree* data= (TTree*)inputFile->Get("data");	
	if(!data){
		cerr<<"ERROR: Cannot get source data from file "<<filename<<"!"<<endl;
		return nullptr;
	}

	double fittedFlux;
	double beamArea;

	data->SetBranchAddress("fittedFlux",&fittedFlux);
	data->SetBranchAddress("beamArea",&beamArea);

	cout<<"# "<<data->GetEntries()<<" sources in file "<<filename<<" ..."<<endl;
	
	//- Fill histo
	for(int i=0;i<data->GetEntries();i++){
		data->GetEntry(i);

		double S= fittedFlux/beamArea;
		double lgS= log10(S);

		if(logyscale) dNdSHisto->Fill(lgS);	
		else dNdSHisto->Fill(S);	

	}//end loop sources

	
	//- Compute source counts histo	
	TH1D* SourceCounts= (TH1D*)dNdSHisto->Clone("SourceCounts");
	/*
	TH1D* SourceCounts= 0;
	if(logyscale) SourceCounts= new TH1D("SourceCounts","SourceCounts",nBins,LgFluxBins.data());
	else SourceCounts= (TH1D*)dNdSHisto->Clone("SourceCounts");
	SourceCounts->Sumw2();
	*/

	if(logyscale){
		for(int i=0;i<dNdSHisto->GetNbinsX();i++)
		{
			double N= dNdSHisto->GetBinContent(i+1);
			double NErr= dNdSHisto->GetBinError(i+1); 
			if(N<=0) continue;
			double lgS= dNdSHisto->GetBinCenter(i+1);		
			double dlgS= dNdSHisto->GetBinWidth(i+1);
			double lgS_min= dNdSHisto->GetBinLowEdge(i+1);
			double lgS_max= lgS_min + dlgS;
			double dS= pow(10,lgS_max-lgS_min);

			double S= pow(10,lgS);		
			double NormFactor= pow(S,-gamma);
		
			double counts= N/(dS*A_sr);//in units: Jy^-1 sr^-1 
			double counts_norm= counts/NormFactor;//in units: Jy^-2.5 sr^-1
			double countsErr= NErr/(dS*A_sr);
			double countsErr_norm= countsErr/NormFactor;

			if(normalize){
				SourceCounts->SetBinContent(i+1,counts_norm);
				SourceCounts->SetBinError(i+1,countsErr_norm);
			}
			else{
				SourceCounts->SetBinContent(i+1,counts);
				SourceCounts->SetBinError(i+1,countsErr);
			}
		
			cout<<"INFO: bin "<<i+1<<": N="<<N<<", S(Jy)="<<S<<", dS(Jy)="<<dS<<", normFactor="<<NormFactor<<", counts="<<counts<<" +- "<<countsErr<<", counts(norm)="<<counts_norm<<" +- "<<countsErr_norm<<endl;

		}//end loop bins
	}
	else{

		for(int i=0;i<dNdSHisto->GetNbinsX();i++)
		{
			double N= dNdSHisto->GetBinContent(i+1);
			double NErr= dNdSHisto->GetBinError(i+1); 
			if(N<=0) continue;
			double S= dNdSHisto->GetBinCenter(i+1);		
			double dS= dNdSHisto->GetBinWidth(i+1);		
			double NormFactor= pow(S,-gamma);
		
			double counts= N/(dS*A_sr);//in units: Jy^-1 sr^-1 
			double counts_norm= counts/NormFactor;//in units: Jy^-2.5 sr^-1
			double countsErr= NErr/(dS*A_sr);
			double countsErr_norm= countsErr/NormFactor;

			if(normalize){
				SourceCounts->SetBinContent(i+1,counts_norm);
				SourceCounts->SetBinError(i+1,countsErr_norm);
			}
			else{
				SourceCounts->SetBinContent(i+1,counts);
				SourceCounts->SetBinError(i+1,countsErr);
			}
		
			cout<<"INFO: bin "<<i+1<<": N="<<N<<", S(Jy)="<<S<<", dS(Jy)="<<dS<<", normFactor="<<NormFactor<<", counts="<<counts<<" +- "<<countsErr<<", counts(norm)="<<counts_norm<<" +- "<<countsErr_norm<<endl;

		}//end loop bins
	}//close else

	// - Clear histo
	if(dNdSHisto){
		dNdSHisto->Delete();
	}

	return SourceCounts;

}//close ComputeSourceCounts()



