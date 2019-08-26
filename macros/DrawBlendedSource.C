#include <Logger.h>

int peakShiftTolerance= 2;
int peakKernelMultiplicityThr= 1;
int maxNComponents= -1;
//std::vector<int> kernels= {3,5,7};
std::vector<int> kernels= {3};
double thrFactor= 1;
bool skipBorders= true;

int DrawBlendedSource(std::string fileName,int sourceIndex,int nestedSourceIndex=-1,double peakZThr=1,double curvThr=0)
{
	//Open file with source collection
	TFile* inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()){
		ERROR_LOG("Failed to open file "<<fileName<<"!");
		return -1;
	}

	//Get source tree
	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");
	if(!SourceInfo){
		ERROR_LOG("Failed to get source info tree from file "<<fileName<<"!");
		return -1;
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	//Read source
	SourceInfo->GetEntry(sourceIndex);
	if(!aSource){
		ERROR_LOG("Failed to read source with index "<<sourceIndex<<" from file "<<fileName<<"!");
		return -1;
	}
	INFO_LOG("Reading source "<<aSource->GetName()<<" ...");

	//If nested source index is given access to nested source
	Source* source= aSource;
	if(nestedSourceIndex!=-1 && nestedSourceIndex>=0){
		source= aSource->GetNestedSource(nestedSourceIndex);
		if(!source){
			ERROR_LOG("Failed to get access to nested source with index "<<nestedSourceIndex<<"!");
			return -1;
		}
	}

	//## Get source map & histo
	Image* sourceImg= source->GetImage(eFluxMap);
	Image* sourceImg_norm= sourceImg->GetNormalizedImage();
	TH2D* sourceHisto= sourceImg_norm->GetHisto2D("sourceHisto");

	//std::string pixelArrayStr= sourceImg->GetPixelNumpyArrayStr();
	//cout<<"== PIXEL ARRAY =="<<endl;
	//cout<<pixelArrayStr<<endl;
	//cout<<"================="<<endl;

	Image* sourceImg_inv= source->GetImage(eFluxMap);
	sourceImg_inv->Scale(-1);
	Image* sourceImg_inv_norm= sourceImg_inv->GetNormalizedImage();
	TH2D* sourceHisto_inv= sourceImg_inv_norm->GetHisto2D("sourceHisto_inv");

	
	//## Get source curvature map
	bool useRange= true;
	double minRangeThr= 0;
	bool computeRobustStats= true;
	bool forceRecomputing= true;

	//Laplacian filter
	Image* curvImg_laplFilt= sourceImg->GetLaplacianImage(true);
	Image* curvImg_laplFilt_norm= curvImg_laplFilt->GetNormalizedImage();
	//curvImg_laplFilt->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
	//double median= (curvImg_laplFilt->GetPixelStats())->median;
	//double medianThr= thrFactor*median;
	//double thrLevel= medianThr;
	//curvImg_laplFilt->ApplyThreshold(thrLevel);
	//TH2D* curvHisto_laplFilt= curvImg_laplFilt->GetHisto2D("curvHisto_laplFilt");
	TH2D* curvHisto_laplFilt= curvImg_laplFilt_norm->GetHisto2D("curvHisto_laplFilt");

	Image* curvImg_laplFilt_inv= sourceImg->GetLaplacianImage(false);
	Image* curvImg_laplFilt_inv_norm= curvImg_laplFilt_inv->GetNormalizedImage();

	//LoG filter
	Image* curvImg_logFilt= sourceImg->GetLoGImage(true);
	//curvImg_logFilt->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
	//median= (curvImg_logFilt->GetPixelStats())->median;
	//medianThr= thrFactor*median;
	//thrLevel= medianThr;
	//curvImg_logFilt->ApplyThreshold(thrLevel);
	TH2D* curvHisto_logFilt= curvImg_logFilt->GetHisto2D("curvHisto_logFilt");

	//TopHat filter
	Image* filtMap_tophat= sourceImg->GetTopHatImage();
	TH2D* filtHisto_tophat= filtMap_tophat->GetHisto2D("filtHisto_tophat");

	//Reconstruction filter mask
	Image* filtMap_morphReco= sourceImg->GetTopHatImage();
	TH2D* filtHisto_tophat= filtMap_tophat->GetHisto2D("filtHisto_tophat");



	//Get source peaks
	INFO_LOG("Finding source peaks points...");
	/*
	std::vector<TVector2> peaks;
	source->FindComponentPeaks(peaks,peakZThr,maxNComponents);

	TGraph* peakGraph= new TGraph;
	for(int i=0;i<peaks.size();i++){
		peakGraph->SetPoint(i,peaks[i].X(),peaks[i].Y());
	}
	*/

	TGraph* peakGraph= curvImg_laplFilt_norm->ComputePeakGraph(kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr);
	

	//Get source valleys
	INFO_LOG("Finding source valley points...");

	/*
	std::vector<TVector2> valleys;
	bool invertSearch= true;
	source->FindComponentPeaks(valleys,peakZThr,maxNComponents,peakShiftTolerance,kernels,peakKernelMultiplicityThr,invertSearch);

	TGraph* valleyGraph= new TGraph;
	for(int i=0;i<valleys.size();i++){
		valleyGraph->SetPoint(i,valleys[i].X(),valleys[i].Y());
	}
	*/
	

	
	//TGraph* valleyGraph= sourceImg_inv_norm->ComputePeakGraph(kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr);
	TGraph* valleyGraph= curvImg_laplFilt_inv_norm->ComputePeakGraph(kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr);
	

	//Draw source & detected peaks
	TCanvas* Plot= new TCanvas("Plot","Plot");
	Plot->cd();

	sourceHisto->SetStats(0);
	sourceHisto->Draw("COLZ");

	peakGraph->SetMarkerSize(1.3);
	peakGraph->SetMarkerColor(kWhite);
	peakGraph->Draw("P");
	
	valleyGraph->SetMarkerStyle(24);
	valleyGraph->SetMarkerSize(1.3);
	valleyGraph->SetMarkerColor(kBlack);
	valleyGraph->Draw("P");

	TCanvas* Plot2= new TCanvas("Plot2","Plot2");
	Plot2->cd();
	sourceHisto_inv->SetStats(0);
	sourceHisto_inv->Draw("COLZ");

	//Draw curvature map
	TCanvas* CurvMapPlot= new TCanvas("CurvMapPlot","CurvMapPlot");
	CurvMapPlot->cd();

	curvHisto_laplFilt->SetStats(0);
	curvHisto_laplFilt->Draw("COLZ");

	TCanvas* CurvMapPlot2= new TCanvas("CurvMapPlot2","CurvMapPlot2");
	CurvMapPlot2->cd();

	curvHisto_logFilt->SetStats(0);
	curvHisto_logFilt->Draw("COLZ");

	//Draw curvature map
	TCanvas* FiltMapPlot= new TCanvas("FiltMapPlot","FiltMapPlot");
	FiltMapPlot->cd();

	filtHisto_tophat->SetStats(0);
	filtHisto_tophat->Draw("COLZ");

	return 0;

}//close macro

