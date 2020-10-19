#include <Image.h>
#include <MathUtils.h>
#include <StatsUtils.h>
#include <AstroUtils.h>
#include <DS9Region.h>
#include <DS9RegionParser.h>

#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace Caesar;

//==========================================
//==     CLASS DEFINITIONS
//==========================================
struct SourceInfo 
{
	std::string name;
	double x;
	double y;
	double ra;
	double dec;
	double l;
	double b;
	double flux;
	double rms;
	double snr;
	double chi2;
	double ndf;
	double A;
	double bmaj_pix;
	double bmin_pix;
	double pa_pix;
	double bmaj_wcs;
	double bmin_wcs;
	double pa_wcs;
	double eccentricityRatio;
	double areaRatio;
	double rotAngleVSBeam;
	int fitQuality;
	int flag;

	SourceInfo()
	{
		name= "";
		x= 0;
		y= 0;
		ra= 0;
		dec= 0;
		l= 0;
		b= 0;
		flux= 0;
		rms= 0;
		snr= 0;
		chi2= 0;
		ndf= 0;
		A= 0;
		eccentricityRatio= 0;
		areaRatio= 0;
		rotAngleVSBeam= 0;
		fitQuality= 0;
		flag= 0;
		bmaj_pix= 0;
		bmin_pix= 0;
		pa_pix= 0;
		bmaj_wcs= 0;
		bmin_wcs= 0;
		pa_wcs= 0;
	}

	SourceInfo(std::string _name,double _x,double _y,double _flux)
		: name(_name), x(_x), y(_y), flux(_flux)
	{}

};//close struct SourceInfo

enum CatalogType {eCAESAR_CATALOG=1,eSELAVY_CATALOG=2,eAEGEAN_CATALOG=3};

struct MacroOptions
{
	double redChi2Thr;
	double mapPixSize;//in arcsec
	double beamBmaj;//in arcsec
	double beamBmin;//in arcsec
	int catalogType; 
	int catalogVersion;
	double matchPosThr;//in arcsec
	bool saveToFile;
	bool applySelection;	
	bool useWCS;
	bool selectSourcesInRegion;
	std::string regionFileName;
	std::string mapFileName;

	double galMinLong;
	double galMaxLong;
	double galLongStep; 
	double galMinLat;
	double galMaxLat;
	double galLatStep;

	MacroOptions()
	{
		Init();
	}
	~MacroOptions(){}

	void Init()
	{
		redChi2Thr= 10.;
		mapPixSize= 4;//in arcsec
		beamBmaj= 24.007906;//in arcsec
		beamBmin= 20.992698;//in arcsec	
		catalogType= 1;
		catalogVersion= 3;
		matchPosThr= 8;//arcsec
		saveToFile= true;
		applySelection= true;
		selectSourcesInRegion= false;
		regionFileName= "";
		useWCS= true;
		mapFileName= "";

		galMinLong= 335; 
		galMaxLong= 355;
		galLongStep= 1; 
		galMinLat= -5;
		galMaxLat= 5;
		galLatStep= 0.25;
	}
	
};//close struct MacroOptions()

//==========================================
//==     METHODS
//==========================================
int ReadData(std::vector<SourceInfo>& sources,std::string filename,int (*fcn)(SourceInfo&, const std::vector<std::string>&), std::vector<std::string> excludedPatterns= {});
static int ParseData(SourceInfo& info,const std::vector<std::string>& fields);
static int ParseData_caesar(SourceInfo& info,const std::vector<std::string>& fields);
static int ParseData_selavy(SourceInfo& info,const std::vector<std::string>& fields);
static int ParseData_aegean(SourceInfo& info,const std::vector<std::string>& fields);
bool HasPattern(std::string,std::string);
int FindSourceMatch(const std::vector<SourceInfo>& sources_ref,const std::vector<SourceInfo>& sources_rec);
void Init();
void Draw();
void Save();
int AnalyzeData();

//==========================================
//==     VARS
//==========================================
MacroOptions gOpts;
double gMatchPosThr= 10;//in arcsec
int gCatalogType= eCAESAR_CATALOG;
int gCatalogVersion= 2;
bool gApplySelection= true;
double gRedChi2Thr= 10.;
double gMapPixSize= 4;//in arcsec
double gBeamBmaj= 24.007906;//in arcsec
double gBeamBmin= 20.992698;//in arcsec
double gBeamArea_wcs= MathUtils::ComputeEllipseArea(gBeamBmaj,gBeamBmin);
bool gSelectRegion= false;
bool gUseWCS= true;
std::string gRegionFileName= "";
std::vector<DS9Region*> gRegions;
std::string gMapFileName;
Image* gImage= 0;
WCS* gGalWCS= 0;
WCS* gWCS= 0;
bool gSaveToFile= true;

//- REFERENCE CATALOG
const int gNDataCols= 6;
std::vector<std::string> gFileNames {};
std::vector<std::string> gSkipLinePatterns {};


//- REC CATALOG
const int gNDataCols_caesar= 42;
const int gNDataCols_caesar_v2= 47;
const int gNDataCols_caesar_v3= 47;
const int gNDataCols_selavy= 37;
const int gNDataCols_aegean= 27;
std::vector<std::string> gSkipLinePatterns_caesar {};
std::vector<std::string> gSkipLinePatterns_selavy {};
std::vector<std::string> gSkipLinePatterns_aegean {"island"};
std::vector<std::string> gFileNames_rec {};

//- OUTPUT FILE
TFile* gOutputFile= 0;
TTree* gMatchDataTree= 0;
std::string gSourceName_true;
double gSourcePosX_true;
double gSourcePosY_true;
double gSourceGalPosX_true;
double gSourceGalPosY_true;
double gSourceJ2000PosX_true;
double gSourceJ2000PosY_true;
double gSourceFlux_true;
std::string gSourceName_rec;
double gSourcePosX_rec;
double gSourcePosY_rec;
double gSourceFlux_rec;
double gSourceGalPosX_rec;
double gSourceGalPosY_rec;
double gSourceJ2000PosX_rec;
double gSourceJ2000PosY_rec;
bool gIsFound;

//- HISTOS
double gLgFluxMin_draw= -5;
double gLgFluxMax_draw= 1;
//std::vector<double> gLgFluxBins {-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
//std::vector<double> gLgRecFluxBins {-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
std::vector<double> gLgFluxBins {-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
std::vector<double> gLgRecFluxBins {-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
//std::vector<double> gSNRBins {0,1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,1000,10000};
std::vector<double> gSNRBins {0,2,4,5,7,10,15,20,25,50,100,1000,10000};
std::vector<double> gGalLatBins {-5,-4,-3,-2,-1,0,1,2,3,4,5};
TProfile* gPosXDiffHisto= 0;
TProfile* gPosYDiffHisto= 0;
TProfile* gFluxDiffHisto= 0;
TProfile* gFluxDiffVSRecFluxHisto= 0;
TH1D* gFluxHisto= 0;
TH1D* gFluxHisto_reliability= 0;
TH1D* gFluxHisto_fdr= 0;
TH2D* gFluxGalLatHisto2D= 0;
TH2D* gFluxGalLatHisto2D_reliability= 0;
TH2D* gFluxHisto2D= 0;
TH2D* gCompleteness2D= 0;
TH2D* gCompletenessVSFluxVSGalLat= 0;
TH2D* gReliabilityVSFluxVSGalLat= 0;

TH2D* gSourceCountsHisto2D_galcoord= 0;
TH2D* gSourceCountsHisto2D_galcoord_rec= 0;
TH2D* gCompleteness2D_galcoord= 0;

std::vector<TH2D*> gSourceCountsHistos2D_galcoord;
std::vector<TH2D*> gSourceCountsHistos2D_galcoord_rec;
std::vector<TH2D*> gCompletenessList2D_galcoord;

TH1D* gEccentricityRatioHisto= 0;
TH1D* gAreaRatioHisto= 0;
TH1D* gRotAngleVSBeamHisto= 0;


//- Rec vectors & graphs
std::vector<std::vector<double>> gFluxBinData;
std::vector<std::vector<double>> gFluxDiffBinData;
std::vector<std::vector<double>> gLgFluxDiffBinData;
std::vector<std::vector<double>> gPosXDiffBinData;
std::vector<std::vector<double>> gPosYDiffBinData;

std::vector<std::vector<double>> gSNRBinData;
std::vector<std::vector<double>> gFluxDiffVSSNRBinData;
std::vector<std::vector<double>> gPosXDiffVSSNRBinData;
std::vector<std::vector<double>> gPosYDiffVSSNRBinData;

std::vector<std::vector<double>> gRecFluxBinData;
std::vector<std::vector<double>> gLgFluxDiffVSRecFluxBinData;

std::vector<TH1D*> gFluxDiffHistos;
std::vector<TF1*> gFluxDiffFitFcns;
TH1D* gFluxHisto_rec= 0;
TH1D* gFluxHisto_rec_reliability= 0;
TH1D* gFluxHisto_rec_fdr= 0;
TH2D* gFluxHisto2D_rec= 0;
TH2D* gFluxGalLatHisto2D_rec= 0;
TH2D* gFluxGalLatHisto2D_rec_reliability= 0;
std::vector<TH2D*> gSourceCountsHistos2D_galcoord_reliability;
std::vector<TH2D*> gSourceCountsHistos2D_galcoord_rec_reliability;
std::vector<TH2D*> gReliabilityList2D_galcoord;

TEfficiency* gCompleteness= 0;
TGraphAsymmErrors* gCompletenessGraph= 0;
TEfficiency* gReliability= 0;
TGraphAsymmErrors* gReliabilityGraph= 0;
TEfficiency* gFalseDetectionRate= 0;
TGraphAsymmErrors* gFalseDetectionRateGraph= 0;
TGraphErrors* gFluxBiasGraph= 0;
TGraphErrors* gFluxResoGraph= 0; 
TGraphErrors* gLgFluxBiasGraph= 0;
TGraphErrors* gLgFluxResoGraph= 0;
TGraphErrors* gLgFluxBiasVSRecFluxGraph= 0;
TGraphErrors* gLgFluxResoVSRecFluxGraph= 0; 
TGraphErrors* gPosXBiasGraph= 0; 
TGraphErrors* gPosYBiasGraph= 0; 
TGraphErrors* gPosXResoGraph= 0; 
TGraphErrors* gPosYResoGraph= 0;

TH1D* gSNRHisto_rec= 0;
TGraphErrors* gFluxBiasVSSNRGraph= 0;
TGraphErrors* gFluxResoVSSNRGraph= 0; 
TGraphErrors* gPosXBiasVSSNRGraph= 0; 
TGraphErrors* gPosYBiasVSSNRGraph= 0; 
TGraphErrors* gPosXResoVSSNRGraph= 0; 
TGraphErrors* gPosYResoVSSNRGraph= 0; 

//- Selection vectors & graphs
std::vector<std::vector<double>> gFluxBinData_sel;
std::vector<std::vector<double>> gLgFluxDiffBinData_sel;
std::vector<std::vector<double>> gFluxDiffBinData_sel;
std::vector<std::vector<double>> gPosXDiffBinData_sel;
std::vector<std::vector<double>> gPosYDiffBinData_sel;

std::vector<std::vector<double>> gRecFluxBinData_sel;
std::vector<std::vector<double>> gLgFluxDiffVSRecFluxBinData_sel;

TH1D* gFluxHisto_sel= 0;
TH1D* gFluxHisto_sel_reliability= 0;
TH1D* gFluxHisto_sel_fdr= 0;
TEfficiency* gCompleteness_sel= 0;
TGraphAsymmErrors* gCompletenessGraph_sel= 0;
TEfficiency* gReliability_sel= 0;
TGraphAsymmErrors* gReliabilityGraph_sel= 0;
TEfficiency* gFalseDetectionRate_sel= 0;
TGraphAsymmErrors* gFalseDetectionRateGraph_sel= 0;
TGraphErrors* gFluxBiasGraph_sel= 0;
TGraphErrors* gFluxResoGraph_sel= 0; 
TGraphErrors* gLgFluxBiasGraph_sel= 0;
TGraphErrors* gLgFluxResoGraph_sel= 0; 
TGraphErrors* gLgFluxBiasVSRecFluxGraph_sel= 0;
TGraphErrors* gLgFluxResoVSRecFluxGraph_sel= 0; 
TGraphErrors* gPosXBiasGraph_sel= 0; 
TGraphErrors* gPosYBiasGraph_sel= 0; 
TGraphErrors* gPosXResoGraph_sel= 0; 
TGraphErrors* gPosYResoGraph_sel= 0; 


int CrossMatchSources(std::string fileName_rec,std::string fileName,bool isFileList=false,MacroOptions options=MacroOptions()) 
{
	//================================
	//==      CHECK ARGS
	//================================
	cout<<"INFO: Check macro args ..."<<endl;
	gOpts= options;
	gSaveToFile= options.saveToFile;
	gMatchPosThr= options.matchPosThr;
	gCatalogType= options.catalogType;
	gApplySelection= options.applySelection;
	gCatalogVersion= options.catalogVersion;
	gSelectRegion= options.selectSourcesInRegion;
	gRegionFileName= options.regionFileName;
	gMapFileName= options.mapFileName;

	gRedChi2Thr= options.redChi2Thr;
	gMapPixSize= options.mapPixSize;
	gBeamBmaj= options.beamBmaj; 
	gBeamBmin= options.beamBmin;
	gBeamArea_wcs= MathUtils::ComputeEllipseArea(gBeamBmaj,gBeamBmin);
	gUseWCS= options.useWCS;
	
	if( gCatalogType!=eCAESAR_CATALOG && gCatalogType!=eSELAVY_CATALOG && gCatalogType!=eAEGEAN_CATALOG )
	{
		cerr<<"ERROR: Unknown catalog type given ("<<gCatalogType<<")!"<<endl;
		return -1;
	}

	//Fill filelist	
	if(fileName=="" || fileName_rec==""){
		cerr<<"ERROR: Empty filename/filelist given!"<<endl;
		return -1;
	}

	//================================
	//==      INIT
	//================================
	cout<<"INFO: Init macro data ..."<<endl;
	Init();


	//================================
	//==      READ FILE LISTS
	//================================
	gFileNames_rec.clear();
	gFileNames.clear();

	if(isFileList)
	{
		//Check filenames
		std::ifstream fileStream(fileName);	
		std::ifstream fileStream_rec(fileName_rec);
		std::string line;
		std::string filename= "";
		if (fileStream.fail() || !fileStream.is_open()){
			cerr<<"ERROR: Failed to open file "<<fileName<<" for reading..."<<endl;
			return -1;
		}
		if (fileStream_rec.fail() || !fileStream_rec.is_open()){
			cerr<<"ERROR: Failed to open file "<<fileName_rec<<" for reading..."<<endl;
			return -1;
		}

		//Store filenames present in lists
		std::vector<std::string> fileNames;
		while (std::getline(fileStream, line)) {
    	std::istringstream iss(line);
    	if (!(iss >> filename)) { 
				cerr<<"ERROR: Failed to read line from file "<<fileName<<"!"<<endl;
				return -1; 
			}
    	gFileNames.push_back(filename);
		}//end file read
				
		while (std::getline(fileStream_rec, line)) {
    	std::istringstream iss(line);
    	if (!(iss >> filename)) { 
				cerr<<"ERROR: Failed to read line from file "<<fileName<<"!"<<endl;
				return -1; 
			}
    	gFileNames_rec.push_back(filename);
		}//end file read

		//Check files have same number of entries
		if(gFileNames.size()!=gFileNames_rec.size()){
			cerr<<"ERROR: Input filelist must have the same number of entries!"<<endl;
			return -1;
		}
	}//close if
	else{
		gFileNames_rec.push_back(fileName_rec);
		gFileNames.push_back(fileName);
	}

	//================================
	//==      SELECT REGION
	//================================
	if(gSelectRegion){
		//Read and parse DS9 regions
		if (DS9RegionParser::Parse(gRegions,gRegionFileName)<0){
			cerr<<"ERROR: Failed to read and parse DS9 region file "<<gRegionFileName<<"!"<<endl;
			return -1;
		}

		//Check if empty
		if(gRegions.empty()){
			cerr<<"No regions read from file "<<gRegionFileName<<"!"<<endl;
			return -1;
		}
		else{
			cout<<"INFO: #"<<gRegions.size()<<" regions read from file "<<gRegionFileName<<"..."<<endl;
		}

		//Print regions
		for(size_t i=0;i<gRegions.size();i++){
			cout<<"--> DS9 Region no. "<<i+1<<endl;
			gRegions[i]->Print();
		}
	}//close if gSelectRegion()


	//================================
	//==      READ MAP
	//================================
	if(gMapFileName!="")
	{
		//Read image
		cout<<"INFO: Reading map from file "<<gMapFileName<<" ..."<<endl;
		gImage= new Image;
		if(gImage->ReadFITS(gMapFileName)<0){
			cerr<<"ERROR: Failed to read map from file "<<gMapFileName<<"!"<<endl;
			return -1;
		}

		//Get metadata
		ImgMetaData* metadata= gImage->GetMetaData();
		if(!metadata){
			cerr<<"ERROR: Failed to get image metadata!"<<endl;
			return -1;
		}
			
		//Get WCS
		gGalWCS= metadata->GetWCS(eGALACTIC);
		if(!gGalWCS){
			cerr<<"ERROR: Failed to get galactic WCS!"<<endl;
			return -1;
		}
	
		gWCS= metadata->GetWCS(eJ2000);
		if(!gWCS){
			cerr<<"ERROR: Failed to get galactic WCS!"<<endl;
			return -1;
		}
	}

	

	//================================
	//==      READ CATALOG DATA
	//================================
	//## Read catalog data
	int nFiles= (int)(gFileNames.size());
	for(int i=0;i<nFiles;i++){
		//Read reference catalog data
		std::vector<SourceInfo> sources;
		ReadData(sources,gFileNames[i],ParseData,gSkipLinePatterns);
		cout<<"INFO: #"<<sources.size()<<" sources found in reference catalog file "<<gFileNames[i]<<endl;
		
		//Read rec catalog data
		std::vector<SourceInfo> sources_rec;
		if(gCatalogType==eCAESAR_CATALOG){
			ReadData(sources_rec,gFileNames_rec[i],ParseData_caesar,gSkipLinePatterns_caesar);
		}
		else if(gCatalogType==eSELAVY_CATALOG){
			ReadData(sources_rec,gFileNames_rec[i],ParseData_selavy,gSkipLinePatterns_selavy);
		}
		else if(gCatalogType==eAEGEAN_CATALOG){
			ReadData(sources_rec,gFileNames_rec[i],ParseData_aegean,gSkipLinePatterns_aegean);
		}
		else{
			cerr<<"ERROR: Unknown/invalid catalog type ("<<gCatalogType<<")!"<<endl;
			return -1;
		}
		cout<<"INFO: #"<<sources_rec.size()<<" ref sources found in catalog file "<<gFileNames_rec[i]<<endl;

		//Find matches
		FindSourceMatch(sources,sources_rec);

	}//end loop files

	//================================
	//==      ANALYZE DATA
	//================================
	AnalyzeData();

	//================================
	//==      DRAW PLOTS
	//================================
	//Draw();

	//================================
	//==      SAVE
	//================================
	if(gSaveToFile) Save();

	return 0;

}//close macro

void Save()
{
	//Save file
	if(gOutputFile && gOutputFile->IsOpen())
	{
		if(gMatchDataTree) gMatchDataTree->Write();

		if(gPosXDiffHisto) gPosXDiffHisto->Write();
 		if(gPosYDiffHisto) gPosYDiffHisto->Write();
		if(gFluxDiffHisto) gFluxDiffHisto->Write();
		if(gFluxDiffVSRecFluxHisto) gFluxDiffVSRecFluxHisto->Write();
		if(gFluxHisto) gFluxHisto->Write();
		if(gFluxHisto_reliability) gFluxHisto_reliability->Write();			
		if(gFluxHisto_fdr) gFluxHisto_fdr->Write();	
		if(gFluxHisto2D) gFluxHisto2D->Write();
		if(gFluxGalLatHisto2D) gFluxGalLatHisto2D->Write();	
		if(gFluxGalLatHisto2D_reliability) gFluxGalLatHisto2D_reliability->Write();	
		if(gCompleteness2D) gCompleteness2D->Write();
		if(gCompletenessVSFluxVSGalLat) gCompletenessVSFluxVSGalLat->Write();
		if(gSourceCountsHisto2D_galcoord) gSourceCountsHisto2D_galcoord->Write();
		if(gSourceCountsHisto2D_galcoord_rec) gSourceCountsHisto2D_galcoord_rec->Write(); 
		if(gCompleteness2D_galcoord) gCompleteness2D_galcoord->Write();
		for(size_t i=0;i<gCompletenessList2D_galcoord.size();i++) gCompletenessList2D_galcoord[i]->Write();
		for(size_t i=0;i<gReliabilityList2D_galcoord.size();i++) gReliabilityList2D_galcoord[i]->Write();

		if(gEccentricityRatioHisto) gEccentricityRatioHisto->Write();
		if(gAreaRatioHisto) gAreaRatioHisto->Write();
		if(gRotAngleVSBeamHisto) gRotAngleVSBeamHisto->Write();

		//Write rec histos & graphs
		if(gFluxHisto_rec) gFluxHisto_rec->Write();
		if(gFluxHisto_rec_reliability) gFluxHisto_rec_reliability->Write();
		if(gFluxHisto2D_rec) gFluxHisto2D_rec->Write();		
		if(gFluxGalLatHisto2D_rec) gFluxGalLatHisto2D_rec->Write();		
		if(gFluxGalLatHisto2D_rec_reliability) gFluxGalLatHisto2D_rec_reliability->Write();
		if(gFluxHisto_rec_fdr) gFluxHisto_rec_fdr->Write();
		if(gCompleteness) gCompleteness->Write();
		if(gCompletenessGraph) gCompletenessGraph->Write();
		if(gReliability) gReliability->Write();
		if(gReliabilityGraph) gReliabilityGraph->Write();
		if(gReliabilityVSFluxVSGalLat) gReliabilityVSFluxVSGalLat->Write();
		if(gFalseDetectionRate) gFalseDetectionRate->Write();
		if(gFalseDetectionRateGraph) gFalseDetectionRateGraph->Write();
		if(gFluxBiasGraph) gFluxBiasGraph->Write();
		if(gFluxResoGraph) gFluxResoGraph->Write();
		if(gLgFluxBiasGraph) gLgFluxBiasGraph->Write();
		if(gLgFluxResoGraph) gLgFluxResoGraph->Write();
		if(gLgFluxBiasVSRecFluxGraph) gLgFluxBiasVSRecFluxGraph->Write();
		if(gLgFluxResoVSRecFluxGraph) gLgFluxResoVSRecFluxGraph->Write();
		if(gPosXBiasGraph) gPosXBiasGraph->Write();
		if(gPosYBiasGraph) gPosYBiasGraph->Write();
		if(gPosXResoGraph) gPosXResoGraph->Write();
		if(gPosYResoGraph) gPosYResoGraph->Write();

		if(gSNRHisto_rec) gSNRHisto_rec->Write();
		if(gFluxBiasVSSNRGraph) gFluxBiasVSSNRGraph->Write();
		if(gFluxResoVSSNRGraph) gFluxResoVSSNRGraph->Write();
		if(gPosXBiasVSSNRGraph) gPosXBiasVSSNRGraph->Write();
		if(gPosYBiasVSSNRGraph) gPosYBiasVSSNRGraph->Write();
		if(gPosXResoVSSNRGraph) gPosXResoVSSNRGraph->Write();
		if(gPosYResoVSSNRGraph) gPosYResoVSSNRGraph->Write();

		//Write sel histo & graphs
		if(gFluxHisto_sel) gFluxHisto_sel->Write();
		if(gCompleteness_sel) gCompleteness_sel->Write();
		if(gCompletenessGraph_sel) gCompletenessGraph_sel->Write();
		if(gFluxHisto_sel_reliability) gFluxHisto_sel_reliability->Write();
		if(gFluxHisto_sel_fdr) gFluxHisto_sel_fdr->Write();
		if(gReliability_sel) gReliability_sel->Write();
		if(gReliabilityGraph_sel) gReliabilityGraph_sel->Write();
		if(gFalseDetectionRate_sel) gFalseDetectionRate_sel->Write();
		if(gFalseDetectionRateGraph_sel) gFalseDetectionRateGraph_sel->Write();
		if(gFluxBiasGraph_sel) gFluxBiasGraph_sel->Write();
		if(gFluxResoGraph_sel) gFluxResoGraph_sel->Write();
		if(gLgFluxBiasGraph_sel) gLgFluxBiasGraph_sel->Write();
		if(gLgFluxResoGraph_sel) gLgFluxResoGraph_sel->Write();
		if(gLgFluxBiasVSRecFluxGraph_sel) gLgFluxBiasVSRecFluxGraph_sel->Write();
		if(gLgFluxResoVSRecFluxGraph_sel) gLgFluxResoVSRecFluxGraph_sel->Write();
		if(gPosXBiasGraph_sel) gPosXBiasGraph_sel->Write();
		if(gPosYBiasGraph_sel) gPosYBiasGraph_sel->Write();
		if(gPosXResoGraph_sel) gPosXResoGraph_sel->Write();
		if(gPosYResoGraph_sel) gPosYResoGraph_sel->Write();		

		gOutputFile->Close();
	}

}//close Save()

void Init()
{
	//Create output file
	gOutputFile= new TFile("output.root","RECREATE");
	gOutputFile->cd();

	gMatchDataTree= new TTree("MatchInfo","MatchInfo");
	//gMatchDataTree->Branch("name",&gSourceName_true);
	gMatchDataTree->Branch("x",&gSourcePosX_true);
	gMatchDataTree->Branch("y",&gSourcePosY_true);
	gMatchDataTree->Branch("l",&gSourceGalPosX_true);
	gMatchDataTree->Branch("b",&gSourceGalPosY_true);
	gMatchDataTree->Branch("ra",&gSourceJ2000PosX_true);
	gMatchDataTree->Branch("dec",&gSourceJ2000PosY_true);
	gMatchDataTree->Branch("flux",&gSourceFlux_true);
	gMatchDataTree->Branch("isFound",&gIsFound);
	//gMatchDataTree->Branch("name_rec",&gSourceName_rec);
	gMatchDataTree->Branch("x_rec",&gSourcePosX_rec);
	gMatchDataTree->Branch("y_rec",&gSourcePosY_rec);
	gMatchDataTree->Branch("l_rec",&gSourceGalPosX_rec);
	gMatchDataTree->Branch("b_rec",&gSourceGalPosY_rec);
	gMatchDataTree->Branch("ra_rec",&gSourceJ2000PosX_rec);
	gMatchDataTree->Branch("dec_rec",&gSourceJ2000PosY_rec);
	gMatchDataTree->Branch("flux_rec",&gSourceFlux_rec);
	


	//Init histos
	int nBins= (int)(gLgFluxBins.size()-1);
	int nBins_snr= (int)(gSNRBins.size()-1);
	int nBins_rec= (int)(gLgRecFluxBins.size()-1);
	int nBins_galLat= (int)(gGalLatBins.size()-1);
	
	gPosXDiffHisto= new TProfile("posXDiffHisto","",nBins,gLgFluxBins.data(),"S");
	gPosXDiffHisto->Sumw2();

	gPosYDiffHisto= new TProfile("posYDiffHisto","",nBins,gLgFluxBins.data(),"S");
	gPosYDiffHisto->Sumw2();

	gFluxDiffHisto= new TProfile("fluxDiffHisto","",nBins,gLgFluxBins.data(),"S");
	gFluxDiffHisto->Sumw2();

	gFluxDiffVSRecFluxHisto= new TProfile("fluxDiffVSRecFluxHisto","",nBins,gLgRecFluxBins.data(),"S");
	gFluxDiffVSRecFluxHisto->Sumw2();

	gFluxHisto= new TH1D("fluxHisto","",nBins,gLgFluxBins.data());
	gFluxHisto->Sumw2();


	gFluxHisto_reliability= new TH1D("fluxHisto_reliability","",nBins,gLgFluxBins.data());
	gFluxHisto_reliability->Sumw2();

	gFluxHisto_fdr= new TH1D("fluxHisto_fdr","",nBins,gLgFluxBins.data());
	gFluxHisto_fdr->Sumw2();

	gFluxGalLatHisto2D= new TH2D("fluxGalLatHisto2D","",nBins,gLgFluxBins.data(),nBins_galLat,gGalLatBins.data());
	gFluxGalLatHisto2D->Sumw2();
	
	gFluxGalLatHisto2D_reliability= new TH2D("fluxGalLatHisto2D_reliability","",nBins,gLgFluxBins.data(),nBins_galLat,gGalLatBins.data());
	gFluxGalLatHisto2D_reliability->Sumw2();

	gEccentricityRatioHisto= new TH1D("eccentricityRatio","",100,0,10);
	gEccentricityRatioHisto->Sumw2();

	gAreaRatioHisto= new TH1D("areaRatio","",100,0,20);
	gAreaRatioHisto->Sumw2();

	gRotAngleVSBeamHisto= new TH1D("rotAngleVSBeam","",100,-180,180);
	gRotAngleVSBeamHisto->Sumw2();

	//- Rec graphs
	gFluxHisto_rec= new TH1D("fluxHisto_rec","",nBins,gLgFluxBins.data());
	gFluxHisto_rec->Sumw2();

	gFluxHisto_rec_reliability= new TH1D("fluxHisto_rec_reliability","",nBins,gLgFluxBins.data());
	gFluxHisto_rec_reliability->Sumw2();

	gFluxHisto_rec_fdr= new TH1D("fluxHisto_rec_fdr","",nBins,gLgFluxBins.data());
	gFluxHisto_rec_fdr->Sumw2();

	gFluxGalLatHisto2D_rec= new TH2D("fluxGalLatHisto2D_rec","",nBins,gLgFluxBins.data(),nBins_galLat,gGalLatBins.data());
	gFluxGalLatHisto2D_rec->Sumw2();	

	gFluxGalLatHisto2D_rec_reliability= new TH2D("fluxGalLatHisto2D_rec_reliability","",nBins,gLgFluxBins.data(),nBins_galLat,gGalLatBins.data());
	gFluxGalLatHisto2D_rec_reliability->Sumw2();

	gSNRHisto_rec= new TH1D("snrHisto","",nBins_snr,gSNRBins.data());
	gSNRHisto_rec->Sumw2();

	gFluxBiasGraph= new TGraphErrors;
	gFluxBiasGraph->SetName("fluxBiasGraph");

	gFluxResoGraph= new TGraphErrors;
	gFluxResoGraph->SetName("fluxResoGraph");

	gLgFluxBiasGraph= new TGraphErrors;
	gLgFluxBiasGraph->SetName("lgfluxBiasGraph");

	gLgFluxResoGraph= new TGraphErrors;
	gLgFluxResoGraph->SetName("lgfluxResoGraph");

	gLgFluxBiasVSRecFluxGraph= new TGraphErrors;
	gLgFluxBiasVSRecFluxGraph->SetName("lgfluxBiasVSRecFluxGraph");

	gLgFluxResoVSRecFluxGraph= new TGraphErrors;
	gLgFluxResoVSRecFluxGraph->SetName("lgfluxResoVSRecFluxGraph");

	gPosXBiasGraph= new TGraphErrors;
	gPosXBiasGraph->SetName("posXBiasGraph");

	gPosYBiasGraph= new TGraphErrors;
	gPosYBiasGraph->SetName("posYBiasGraph");

	gPosXResoGraph= new TGraphErrors;
	gPosXResoGraph->SetName("posXResoGraph");

	gPosYResoGraph= new TGraphErrors;
	gPosYResoGraph->SetName("posYResoGraph");

	gFluxBiasVSSNRGraph= new TGraphErrors;
	gFluxBiasVSSNRGraph->SetName("fluxBiasVSSNRGraph");

	gFluxResoVSSNRGraph= new TGraphErrors;
	gFluxResoVSSNRGraph->SetName("fluxResoVSSNRGraph");

	gPosXBiasVSSNRGraph= new TGraphErrors;
	gPosXBiasVSSNRGraph->SetName("posXBiasVSSNRGraph");

	gPosYBiasVSSNRGraph= new TGraphErrors;
	gPosYBiasVSSNRGraph->SetName("posYBiasVSSNRGraph");

	gPosXResoVSSNRGraph= new TGraphErrors;
	gPosXResoVSSNRGraph->SetName("posXResoVSSNRGraph");

	gPosYResoVSSNRGraph= new TGraphErrors;
	gPosYResoVSSNRGraph->SetName("posYResoVSSNRGraph");

	for(int i=0;i<nBins;i++){
		gFluxBinData.push_back( std::vector<double>() );
		gFluxDiffBinData.push_back( std::vector<double>() );
		gLgFluxDiffBinData.push_back( std::vector<double>() );		
		gPosXDiffBinData.push_back( std::vector<double>() );
		gPosYDiffBinData.push_back( std::vector<double>() );
	}
	for(int i=0;i<nBins_snr;i++){
		gSNRBinData.push_back( std::vector<double>() );
		gFluxDiffVSSNRBinData.push_back( std::vector<double>() );
		gPosXDiffVSSNRBinData.push_back( std::vector<double>() );
		gPosYDiffVSSNRBinData.push_back( std::vector<double>() );
	}
	for(int i=0;i<nBins_rec;i++){
		gRecFluxBinData.push_back( std::vector<double>() );
		gLgFluxDiffVSRecFluxBinData.push_back( std::vector<double>() );
	}

	//- Selection graphs
	gFluxHisto_sel= new TH1D("fluxHisto_sel","",nBins,gLgFluxBins.data());
	gFluxHisto_sel->Sumw2();

	gFluxHisto_sel_reliability= new TH1D("fluxHisto_sel_reliability","",nBins,gLgFluxBins.data());
	gFluxHisto_sel_reliability->Sumw2();

	gFluxHisto_sel_fdr= new TH1D("fluxHisto_sel_fdr","",nBins,gLgFluxBins.data());
	gFluxHisto_sel_fdr->Sumw2();

	gFluxBiasGraph_sel= new TGraphErrors;
	gFluxBiasGraph_sel->SetName("fluxBiasGraph_sel");

	gFluxResoGraph_sel= new TGraphErrors;
	gFluxResoGraph_sel->SetName("fluxResoGraph_sel");

	gLgFluxBiasGraph_sel= new TGraphErrors;
	gLgFluxBiasGraph_sel->SetName("lgfluxBiasGraph_sel");

	gLgFluxResoGraph_sel= new TGraphErrors;
	gLgFluxResoGraph_sel->SetName("lgfluxResoGraph_sel");

	gLgFluxBiasVSRecFluxGraph_sel= new TGraphErrors;
	gLgFluxBiasVSRecFluxGraph_sel->SetName("lgfluxBiasVSRecFluxGraph_sel");

	gLgFluxResoVSRecFluxGraph_sel= new TGraphErrors;
	gLgFluxResoVSRecFluxGraph_sel->SetName("lgfluxResoVSRecFluxGraph_sel");

	gPosXBiasGraph_sel= new TGraphErrors;
	gPosXBiasGraph_sel->SetName("posXBiasGraph_sel");

	gPosYBiasGraph_sel= new TGraphErrors;
	gPosYBiasGraph_sel->SetName("posYBiasGraph_sel");

	gPosXResoGraph_sel= new TGraphErrors;
	gPosXResoGraph_sel->SetName("posXResoGraph_sel");

	gPosYResoGraph_sel= new TGraphErrors;
	gPosYResoGraph_sel->SetName("posYResoGraph_sel");

	for(int i=0;i<nBins;i++){
		gFluxBinData_sel.push_back( std::vector<double>() );
		gFluxDiffBinData_sel.push_back( std::vector<double>() );	
		gLgFluxDiffBinData_sel.push_back( std::vector<double>() );
		gPosXDiffBinData_sel.push_back( std::vector<double>() );
		gPosYDiffBinData_sel.push_back( std::vector<double>() );
	}
	
	for(int i=0;i<nBins_rec;i++){
		gRecFluxBinData_sel.push_back( std::vector<double>() );
		gLgFluxDiffVSRecFluxBinData_sel.push_back( std::vector<double>() );	
	}
	
}//close Init()

void Draw()
{
	/*
	//Compute completeness
	cout<<"INFO: Compute Completeness ..."<<endl;
	gCompleteness= new TEfficiency(*gFluxHisto_rec,*gFluxHisto); 
	gCompleteness->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gCompleteness->SetConfidenceLevel(0.68);

	gCompletenessGraph= gCompleteness->CreateGraph(); 
	gCompletenessGraph->SetNameTitle("CompletenessGraph","CompletenessGraph");

	gCompleteness_sel= new TEfficiency(*gFluxHisto_sel,*gFluxHisto); 
	gCompleteness_sel->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gCompleteness_sel->SetConfidenceLevel(0.68);

	gCompletenessGraph_sel= gCompleteness_sel->CreateGraph(); 
	gCompletenessGraph_sel->SetNameTitle("CompletenessGraph_sel","CompletenessGraph_sel");

	//Compute completeness
	cout<<"INFO: Compute Completeness 2D map ..."<<endl;
	gCompleteness2D= (TH2D*)gFluxHisto2D_rec->Clone("Completeness2D");
	gCompleteness2D->Divide(gFluxHisto2D);

	cout<<"INFO: Compute Completeness 2D map in gal coord ..."<<endl;
	gCompleteness2D_galcoord= (TH2D*)gSourceCountsHisto2D_galcoord_rec->Clone("Completeness2D_galcoord");
	gCompleteness2D_galcoord->Divide(gSourceCountsHisto2D_galcoord);

	//Compute completeness vs flux vs gal lat
	gCompletenessVSFluxVSGalLat= (TH2D*)gFluxGalLatHisto2D_rec->Clone("CompletenessVSFluxVSGalLat");
	gCompletenessVSFluxVSGalLat->Divide(gFluxGalLatHisto2D);

	TH2D* hh= 0;
	for(size_t k=0;k<gSourceCountsHistos2D_galcoord.size();k++)
	{
		TString histoName= Form("Completeness2D_galcoord_fluxBin%d",(int)(k+1));
		hh= (TH2D*)gSourceCountsHistos2D_galcoord_rec[k]->Clone(histoName);
		hh->Divide(gSourceCountsHistos2D_galcoord[k]);
		gCompletenessList2D_galcoord.push_back(hh);
	}

	//Compute reliability
	cout<<"INFO: Compute Reliability ..."<<endl;
	gReliability= new TEfficiency(*gFluxHisto_rec_reliability,*gFluxHisto_reliability); 
	gReliability->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gReliability->SetConfidenceLevel(0.68);

	gReliabilityGraph= gReliability->CreateGraph(); 
	gReliabilityGraph->SetNameTitle("ReliabilityGraph","ReliabilityGraph");

	gReliability_sel= new TEfficiency(*gFluxHisto_sel_reliability,*gFluxHisto_reliability); 
	gReliability_sel->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gReliability_sel->SetConfidenceLevel(0.68);

	gReliabilityGraph_sel= gReliability_sel->CreateGraph(); 
	gReliabilityGraph_sel->SetNameTitle("ReliabilityGraph_sel","ReliabilityGraph_sel");

	gReliabilityVSFluxVSGalLat= (TH2D*)gFluxGalLatHisto2D_rec_reliability->Clone("ReliabilityVSFluxVSGalLat");
	gReliabilityVSFluxVSGalLat->Divide(gFluxGalLatHisto2D_reliability);

	for(size_t k=0;k<gSourceCountsHistos2D_galcoord_reliability.size();k++)
	{
		TString histoName= Form("Reliability2D_galcoord_fluxBin%d",(int)(k+1));
		hh= (TH2D*)gSourceCountsHistos2D_galcoord_rec_reliability[k]->Clone(histoName);
		hh->Divide(gSourceCountsHistos2D_galcoord_reliability[k]);
		gReliabilityList2D_galcoord.push_back(hh);
	}

	//Compute false detection rate
	cout<<"INFO: Compute False Detection Rate ..."<<endl;
	gFalseDetectionRate= new TEfficiency(*gFluxHisto_rec_fdr,*gFluxHisto_fdr); 
	gFalseDetectionRate->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gFalseDetectionRate->SetConfidenceLevel(0.68);

	gFalseDetectionRateGraph= gFalseDetectionRate->CreateGraph(); 
	gFalseDetectionRateGraph->SetNameTitle("FalseDetectionRateGraph","FalseDetectionRateGraph");

	gFalseDetectionRate_sel= new TEfficiency(*gFluxHisto_sel_fdr,*gFluxHisto_fdr); 
	gFalseDetectionRate_sel->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gFalseDetectionRate_sel->SetConfidenceLevel(0.68);

	gFalseDetectionRateGraph_sel= gFalseDetectionRate_sel->CreateGraph(); 
	gFalseDetectionRateGraph_sel->SetNameTitle("FalseDetectionRateGraph_sel","FalseDetectionRateGraph_sel");
	*/

	//Draw source number plots
	TCanvas* SourceNumberPlot= new TCanvas("SourceNumberPlot","SourceNumberPlot");
	SourceNumberPlot->cd();

	gFluxHisto->GetXaxis()->SetTitle("log_{10}(flux/Jy)");
	gFluxHisto->GetYaxis()->SetTitle("#sources");
	gFluxHisto->SetMarkerStyle(8);
	gFluxHisto->SetMarkerSize(1.3);
	gFluxHisto->SetMarkerColor(kBlack);
	gFluxHisto->SetLineColor(kBlack);
	gFluxHisto->Draw("hist");

	gFluxHisto_rec->SetMarkerStyle(8);
	gFluxHisto_rec->SetMarkerSize(1.3);
	gFluxHisto_rec->SetMarkerColor(kRed);
	gFluxHisto_rec->SetLineColor(kRed);
	gFluxHisto_rec->Draw("hist same");

	gFluxHisto_sel->SetMarkerStyle(8);
	gFluxHisto_sel->SetMarkerSize(1.3);
	gFluxHisto_sel->SetMarkerColor(kGreen);
	gFluxHisto_sel->SetLineColor(kGreen);
	gFluxHisto_sel->Draw("hist same");

	TLegend* SourceNumberPlotLegend= new TLegend(0.7,0.7,0.8,0.8);
	SourceNumberPlotLegend->AddEntry(gFluxHisto,"gen sources","L");
	SourceNumberPlotLegend->AddEntry(gFluxHisto_rec,"det sources","L");
	SourceNumberPlotLegend->AddEntry(gFluxHisto_sel,"sel sources","L");
	SourceNumberPlotLegend->Draw("same");


	//Draw completeness
	TCanvas* CompletenessPlot= new TCanvas("CompletenessPlot","CompletenessPlot");
	CompletenessPlot->cd();

	TH2D* CompletenessPlotBkg= new TH2D("CompletenessPlotBkg","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,0,1.05);
	CompletenessPlotBkg->GetXaxis()->SetTitle("log_{10}(flux/Jy)");
	CompletenessPlotBkg->GetYaxis()->SetTitle("completeness");
	CompletenessPlotBkg->SetStats(0);
	CompletenessPlotBkg->Draw();
	
	gCompletenessGraph->SetMarkerStyle(8);
	gCompletenessGraph->SetMarkerSize(1.3);
	gCompletenessGraph->SetMarkerColor(kBlack);
	gCompletenessGraph->SetLineColor(kBlack);
	gCompletenessGraph->Draw("PLZ same");

	if(gApplySelection){
		gCompletenessGraph_sel->SetMarkerStyle(21);
		gCompletenessGraph_sel->SetMarkerSize(1.3);
		gCompletenessGraph_sel->SetMarkerColor(kRed);
		gCompletenessGraph_sel->SetLineColor(kRed);
		gCompletenessGraph_sel->Draw("PLZ same");
	}

	TLegend* CompletenessPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	CompletenessPlotLegend->SetTextSize(0.045);
	CompletenessPlotLegend->SetTextFont(52);
	CompletenessPlotLegend->AddEntry(gCompletenessGraph,"det sources","PL"); 
	if(gApplySelection) CompletenessPlotLegend->AddEntry(gCompletenessGraph_sel,"sel sources","PL");
	CompletenessPlotLegend->Draw("same");


	//Draw efficiency 2D
	TCanvas* Completeness2DPlot= new TCanvas("Completeness2DPlot","Completeness2DPlot");
	Completeness2DPlot->cd();

	gCompleteness2D->GetXaxis()->SetTitle("x");
	gCompleteness2D->GetYaxis()->SetTitle("y");
	gCompleteness2D->Draw("COLZ");

	//Draw reliability
	TCanvas* ReliabilityPlot= new TCanvas("ReliabilityPlot","ReliabilityPlot");
	ReliabilityPlot->cd();

	TH2D* ReliabilityPlotBkg= new TH2D("ReliabilityPlotBkg","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,0,1.05);
	ReliabilityPlotBkg->GetXaxis()->SetTitle("log_{10}(S_{meas}/Jy)");
	ReliabilityPlotBkg->GetYaxis()->SetTitle("reliability");
	ReliabilityPlotBkg->SetStats(0);
	ReliabilityPlotBkg->Draw();
	
	gReliabilityGraph->SetMarkerStyle(8);
	gReliabilityGraph->SetMarkerSize(1.3);
	gReliabilityGraph->SetMarkerColor(kBlack);
	gReliabilityGraph->SetLineColor(kBlack);
	gReliabilityGraph->Draw("PLZ same");

	if(gApplySelection){
		gReliabilityGraph_sel->SetMarkerStyle(21);
		gReliabilityGraph_sel->SetMarkerSize(1.3);
		gReliabilityGraph_sel->SetMarkerColor(kRed);
		gReliabilityGraph_sel->SetLineColor(kRed);
		gReliabilityGraph_sel->Draw("PLZ same");
	}

	TLegend* ReliabilityPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	ReliabilityPlotLegend->SetTextSize(0.045);
	ReliabilityPlotLegend->SetTextFont(52);
	ReliabilityPlotLegend->AddEntry(gReliabilityGraph,"det sources","PL"); 
	if(gApplySelection) ReliabilityPlotLegend->AddEntry(gReliabilityGraph_sel,"sel sources","PL");
	ReliabilityPlotLegend->Draw("same");


	//Draw false detection rate
	TCanvas* FalseDetectionRatePlot= new TCanvas("FalseDetectionRatePlot","FalseDetectionRatePlot");
	FalseDetectionRatePlot->cd();

	TH2D* FalseDetectionRatePlotBkg= new TH2D("FalseDetectionRatePlotBkg","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,0,1.05);
	FalseDetectionRatePlotBkg->GetXaxis()->SetTitle("log_{10}(S_{meas}/Jy)");
	FalseDetectionRatePlotBkg->GetYaxis()->SetTitle("reliability");
	FalseDetectionRatePlotBkg->SetStats(0);
	FalseDetectionRatePlotBkg->Draw();
	
	gFalseDetectionRateGraph->SetMarkerStyle(8);
	gFalseDetectionRateGraph->SetMarkerSize(1.3);
	gFalseDetectionRateGraph->SetMarkerColor(kBlack);
	gFalseDetectionRateGraph->SetLineColor(kBlack);
	gFalseDetectionRateGraph->Draw("PLZ same");

	if(gApplySelection){
		gFalseDetectionRateGraph_sel->SetMarkerStyle(21);
		gFalseDetectionRateGraph_sel->SetMarkerSize(1.3);
		gFalseDetectionRateGraph_sel->SetMarkerColor(kRed);
		gFalseDetectionRateGraph_sel->SetLineColor(kRed);
		gFalseDetectionRateGraph_sel->Draw("PLZ same");
	}

	TLegend* FalseDetectionRatePlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	FalseDetectionRatePlotLegend->SetTextSize(0.045);
	FalseDetectionRatePlotLegend->SetTextFont(52);
	FalseDetectionRatePlotLegend->AddEntry(gFalseDetectionRateGraph,"det sources","PL"); 
	if(gApplySelection) FalseDetectionRatePlotLegend->AddEntry(gFalseDetectionRateGraph_sel,"sel sources","PL");
	FalseDetectionRatePlotLegend->Draw("same");

	//Draw pos bias
	TCanvas* PosDiffPlot= new TCanvas("PosDiffPlot","PosDiffPlot");
	PosDiffPlot->cd();

	TH2D* PosDiffPlotBkg= new TH2D("","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,-10,10);
	PosDiffPlotBkg->GetXaxis()->SetTitle("log_{10}(F/Jy)");
	PosDiffPlotBkg->GetYaxis()->SetTitle("<#Delta x>, <#Delta y> ('')");
	PosDiffPlotBkg->SetStats(0);
	PosDiffPlotBkg->Draw();

	gPosXDiffHisto->SetMarkerStyle(8);
	gPosXDiffHisto->SetMarkerSize(1.3);
	gPosXDiffHisto->SetStats(0);
	//gPosXDiffHisto->Draw("ep same");

	gPosYDiffHisto->SetMarkerStyle(8);
	gPosYDiffHisto->SetMarkerSize(1.3);
	gPosYDiffHisto->SetStats(0);
	//gPosYDiffHisto->Draw("ep same");
	
	gPosXBiasGraph->SetMarkerStyle(8);
	gPosXBiasGraph->SetMarkerSize(1.3);
	gPosXBiasGraph->SetMarkerColor(kBlack);
	gPosXBiasGraph->SetLineColor(kBlack);
	gPosXBiasGraph->Draw("PZ");

	gPosYBiasGraph->SetMarkerStyle(21);
	gPosYBiasGraph->SetMarkerSize(1.3);
	gPosYBiasGraph->SetMarkerColor(kRed);
	gPosYBiasGraph->SetLineColor(kRed);
	gPosYBiasGraph->Draw("PZ");

	if(gApplySelection){
		gPosXBiasGraph_sel->SetMarkerStyle(24);
		gPosXBiasGraph_sel->SetMarkerSize(1.3);
		gPosXBiasGraph_sel->SetMarkerColor(kBlack);
		gPosXBiasGraph_sel->SetLineColor(kBlack);
		gPosXBiasGraph_sel->Draw("PZ");

		gPosYBiasGraph_sel->SetMarkerStyle(25);
		gPosYBiasGraph_sel->SetMarkerSize(1.3);
		gPosYBiasGraph_sel->SetMarkerColor(kRed);
		gPosYBiasGraph_sel->SetLineColor(kRed);
		gPosYBiasGraph_sel->Draw("PZ");
	}

	TLegend* PosDiffPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	PosDiffPlotLegend->SetTextSize(0.045);
	PosDiffPlotLegend->SetTextFont(52);
	PosDiffPlotLegend->AddEntry(gPosXBiasGraph,"x pos bias (det)","PL"); 
	PosDiffPlotLegend->AddEntry(gPosYBiasGraph,"y pos bias (det)","PL"); 
	if(gApplySelection) {
		PosDiffPlotLegend->AddEntry(gPosXBiasGraph_sel,"x pos bias (sel)","PL");
		PosDiffPlotLegend->AddEntry(gPosYBiasGraph_sel,"y pos bias (sel)","PL");
	}
	PosDiffPlotLegend->Draw("same");

	//Draw pos reso
	TCanvas* PosDiffResoPlot= new TCanvas("PosDiffResoPlot","PosDiffResoPlot");
	PosDiffResoPlot->cd();

	TH2D* PosDiffResoPlotBkg= new TH2D("","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,-1,10);
	PosDiffResoPlotBkg->GetXaxis()->SetTitle("log_{10}(F/Jy)");
	PosDiffResoPlotBkg->GetYaxis()->SetTitle("#sigma(#Delta x), #sigma(#Delta y) ('')");
	PosDiffResoPlotBkg->SetStats(0);
	PosDiffResoPlotBkg->Draw();

	gPosXResoGraph->SetMarkerStyle(8);
	gPosXResoGraph->SetMarkerSize(1.3);
	gPosXResoGraph->SetMarkerColor(kBlack);
	gPosXResoGraph->SetLineColor(kBlack);
	gPosXResoGraph->Draw("PZ");

	gPosYResoGraph->SetMarkerStyle(21);
	gPosYResoGraph->SetMarkerSize(1.3);
	gPosYResoGraph->SetMarkerColor(kRed);
	gPosYResoGraph->SetLineColor(kRed);
	gPosYResoGraph->Draw("PZ");

	if(gApplySelection){
		gPosXResoGraph_sel->SetMarkerStyle(24);
		gPosXResoGraph_sel->SetMarkerSize(1.3);
		gPosXResoGraph_sel->SetMarkerColor(kBlack);
		gPosXResoGraph_sel->SetLineColor(kBlack);
		gPosXResoGraph_sel->Draw("PZ");

		gPosYResoGraph_sel->SetMarkerStyle(25);
		gPosYResoGraph_sel->SetMarkerSize(1.3);
		gPosYResoGraph_sel->SetMarkerColor(kRed);
		gPosYResoGraph_sel->SetLineColor(kRed);
		gPosYResoGraph_sel->Draw("PZ");
	}

	TLegend* PosDiffResoPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	PosDiffResoPlotLegend->SetTextSize(0.045);
	PosDiffResoPlotLegend->SetTextFont(52);
	PosDiffResoPlotLegend->AddEntry(gPosXResoGraph,"x pos reso (det)","PL"); 
	PosDiffResoPlotLegend->AddEntry(gPosYResoGraph,"y pos reso (det)","PL"); 
	if(gApplySelection) {
		PosDiffResoPlotLegend->AddEntry(gPosXResoGraph_sel,"x pos reso (sel)","PL");
		PosDiffResoPlotLegend->AddEntry(gPosYResoGraph_sel,"y pos reso (sel)","PL");
	}
	PosDiffResoPlotLegend->Draw("same");

	//Draw pos bias VS SNR
	TCanvas* PosDiffVSSNRPlot= new TCanvas("PosDiffVSSNRPlot","PosDiffVSSNRPlot");
	PosDiffVSSNRPlot->cd();

	TH2D* PosDiffVSSNRPlotBkg= new TH2D("PosDiffVSSNRPlotBkg","",100,0,100,100,-10,10);
	PosDiffVSSNRPlotBkg->GetXaxis()->SetTitle("S/N");
	PosDiffVSSNRPlotBkg->GetYaxis()->SetTitle("<#Delta x>, <#Delta y> ('')");
	PosDiffVSSNRPlotBkg->SetStats(0);
	PosDiffVSSNRPlotBkg->Draw();

	gPosXBiasVSSNRGraph->SetMarkerStyle(8);
	gPosXBiasVSSNRGraph->SetMarkerSize(1.3);
	gPosXBiasVSSNRGraph->SetMarkerColor(kBlack);
	gPosXBiasVSSNRGraph->SetLineColor(kBlack);
	gPosXBiasVSSNRGraph->Draw("PZ");

	gPosYBiasVSSNRGraph->SetMarkerStyle(21);
	gPosYBiasVSSNRGraph->SetMarkerSize(1.3);
	gPosYBiasVSSNRGraph->SetMarkerColor(kRed);
	gPosYBiasVSSNRGraph->SetLineColor(kRed);
	gPosYBiasVSSNRGraph->Draw("PZ");

	/*
	if(gApplySelection){
		gPosXBiasVSSNRGraph_sel->SetMarkerStyle(24);
		gPosXBiasVSSNRGraph_sel->SetMarkerSize(1.3);
		gPosXBiasVSSNRGraph_sel->SetMarkerColor(kBlack);
		gPosXBiasVSSNRGraph_sel->SetLineColor(kBlack);
		gPosXBiasVSSNRGraph_sel->Draw("PZ");

		gPosYBiasVSSNRGraph_sel->SetMarkerStyle(25);
		gPosYBiasVSSNRGraph_sel->SetMarkerSize(1.3);
		gPosYBiasVSSNRGraph_sel->SetMarkerColor(kRed);
		gPosYBiasVSSNRGraph_sel->SetLineColor(kRed);
		gPosYBiasVSSNRGraph_sel->Draw("PZ");
	}
	*/

	TLegend* PosDiffVSSNRPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	PosDiffVSSNRPlotLegend->SetTextSize(0.045);
	PosDiffVSSNRPlotLegend->SetTextFont(52);
	PosDiffVSSNRPlotLegend->AddEntry(gPosXBiasVSSNRGraph,"x pos bias (det)","PL"); 
	PosDiffVSSNRPlotLegend->AddEntry(gPosYBiasVSSNRGraph,"y pos bias (det)","PL"); 
	/*
	if(gApplySelection) {
		PosDiffVSSNRPlotLegend->AddEntry(gPosXBiasVSSNRGraph_sel,"x pos bias (sel)","PL");
		PosDiffVSSNRPlotLegend->AddEntry(gPosYBiasVSSNRGraph_sel,"y pos bias (sel)","PL");
	}
	*/
	PosDiffVSSNRPlotLegend->Draw("same");

	//Draw pos reso VS SNR
	TCanvas* PosDiffResoVSSNRPlot= new TCanvas("PosDiffResoVSSNRPlot","PosDiffResoVSSNRPlot");
	PosDiffResoVSSNRPlot->cd();

	TH2D* PosDiffResoVSSNRPlotBkg= new TH2D("PosDiffResoVSSNRPlotBkg","",100,0,100,100,-1,10);
	PosDiffResoVSSNRPlotBkg->GetXaxis()->SetTitle("S/N");
	PosDiffResoVSSNRPlotBkg->GetYaxis()->SetTitle("#sigma(#Delta x), #sigma(#Delta y) ('')");
	PosDiffResoVSSNRPlotBkg->SetStats(0);
	PosDiffResoVSSNRPlotBkg->Draw();

	gPosXResoVSSNRGraph->SetMarkerStyle(8);
	gPosXResoVSSNRGraph->SetMarkerSize(1.3);
	gPosXResoVSSNRGraph->SetMarkerColor(kBlack);
	gPosXResoVSSNRGraph->SetLineColor(kBlack);
	gPosXResoVSSNRGraph->Draw("PZ");

	gPosYResoVSSNRGraph->SetMarkerStyle(21);
	gPosYResoVSSNRGraph->SetMarkerSize(1.3);
	gPosYResoVSSNRGraph->SetMarkerColor(kRed);
	gPosYResoVSSNRGraph->SetLineColor(kRed);
	gPosYResoVSSNRGraph->Draw("PZ");

	/*
	if(gApplySelection){
		gPosXResoVSSNRGraph_sel->SetMarkerStyle(24);
		gPosXResoVSSNRGraph_sel->SetMarkerSize(1.3);
		gPosXResoVSSNRGraph_sel->SetMarkerColor(kBlack);
		gPosXResoVSSNRGraph_sel->SetLineColor(kBlack);
		gPosXResoVSSNRGraph_sel->Draw("PZ");

		gPosYResoVSSNRGraph_sel->SetMarkerStyle(25);
		gPosYResoVSSNRGraph_sel->SetMarkerSize(1.3);
		gPosYResoVSSNRGraph_sel->SetMarkerColor(kRed);
		gPosYResoVSSNRGraph_sel->SetLineColor(kRed);
		gPosYResoVSSNRGraph_sel->Draw("PZ");
	}
	*/

	TLegend* PosDiffResoVSSNRPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	PosDiffResoVSSNRPlotLegend->SetTextSize(0.045);
	PosDiffResoVSSNRPlotLegend->SetTextFont(52);
	PosDiffResoVSSNRPlotLegend->AddEntry(gPosXResoVSSNRGraph,"x pos reso (det)","PL"); 
	PosDiffResoVSSNRPlotLegend->AddEntry(gPosYResoVSSNRGraph,"y pos reso (det)","PL"); 
	/*
	if(gApplySelection) {
		PosDiffResoVSSNRPlotLegend->AddEntry(gPosXResoVSSNRGraph_sel,"x pos reso (sel)","PL");
		PosDiffResoVSSNRPlotLegend->AddEntry(gPosYResoVSSNRGraph_sel,"y pos reso (sel)","PL");
	}
	*/
	PosDiffResoVSSNRPlotLegend->Draw("same");

	//Draw flux bias
	TCanvas* FluxDiffPlot= new TCanvas("FluxDiffPlot","FluxDiffPlot");
	FluxDiffPlot->cd();

	TH2D* FluxDiffPlotBkg= new TH2D("FluxDiffPlotBkg","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,-10,10);
	FluxDiffPlotBkg->GetXaxis()->SetTitle("log_{10}(F/Jy)");
	FluxDiffPlotBkg->GetYaxis()->SetTitle("<#Delta F/F>");
	FluxDiffPlotBkg->SetStats(0);
	FluxDiffPlotBkg->Draw();

	gFluxDiffHisto->SetMarkerStyle(8);
	gFluxDiffHisto->SetMarkerSize(1.3);
	//gFluxDiffHisto->Draw("ep same");
	
	gFluxBiasGraph->SetMarkerStyle(8);
	gFluxBiasGraph->SetMarkerSize(1.3);
	gFluxBiasGraph->SetLineColor(kBlack);
	gFluxBiasGraph->Draw("PZ");

	if(gApplySelection){
		gFluxBiasGraph_sel->SetMarkerStyle(21);
		gFluxBiasGraph_sel->SetMarkerSize(1.3);
		gFluxBiasGraph_sel->SetLineColor(kRed);
		gFluxBiasGraph_sel->Draw("PZ");
	}

	TLegend* FluxDiffPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	FluxDiffPlotLegend->SetTextSize(0.045);
	FluxDiffPlotLegend->SetTextFont(52);
	FluxDiffPlotLegend->AddEntry(gFluxBiasGraph,"flux bias (det)","PL"); 
	if(gApplySelection) {
		FluxDiffPlotLegend->AddEntry(gFluxBiasGraph_sel,"flux bias (sel)","PL");
	}
	FluxDiffPlotLegend->Draw("same");


	

	//Draw flux reso
	TCanvas* FluxDiffResoPlot= new TCanvas("FluxDiffResoPlot","FluxDiffResoPlot");
	FluxDiffResoPlot->cd();

	TH2D* FluxDiffResoPlotBkg= new TH2D("FluxDiffResoPlotBkg","",100,gLgFluxMin_draw,gLgFluxMax_draw,100,-1,15);
	FluxDiffResoPlotBkg->GetXaxis()->SetTitle("log_{10}(F/Jy)");
	FluxDiffResoPlotBkg->GetYaxis()->SetTitle("#sigma(#Delta F/F)");
	FluxDiffResoPlotBkg->SetStats(0);
	FluxDiffResoPlotBkg->Draw();

	gFluxResoGraph->SetMarkerStyle(8);
	gFluxResoGraph->SetMarkerSize(1.3);
	gFluxResoGraph->SetLineColor(kBlack);
	gFluxResoGraph->Draw("PZ");

	if(gApplySelection){
		gFluxResoGraph_sel->SetMarkerStyle(21);
		gFluxResoGraph_sel->SetMarkerSize(1.3);
		gFluxResoGraph_sel->SetLineColor(kRed);
		gFluxResoGraph_sel->Draw("PZ");
	}

	TLegend* FluxDiffResoPlotLegend= new TLegend(0.4,0.2,0.7,0.3);
	FluxDiffResoPlotLegend->SetTextSize(0.045);
	FluxDiffResoPlotLegend->SetTextFont(52);
	FluxDiffResoPlotLegend->AddEntry(gFluxResoGraph,"flux reso (det)","PL"); 
	if(gApplySelection) {
		FluxDiffResoPlotLegend->AddEntry(gFluxResoGraph_sel,"flux reso (sel)","PL");
	}
	FluxDiffResoPlotLegend->Draw("same");


	//Draw flux bias vs SNR
	TCanvas* FluxDiffVSSNRPlot= new TCanvas("FluxDiffVSSNRPlot","FluxDiffVSSNRPlot");
	FluxDiffVSSNRPlot->cd();

	TH2D* FluxDiffVSSNRPlotBkg= new TH2D("FluxDiffVSSNRPlotBkg","",100,0,100,100,-10,10);
	FluxDiffVSSNRPlotBkg->GetXaxis()->SetTitle("S/N");
	FluxDiffVSSNRPlotBkg->GetYaxis()->SetTitle("<#Delta F/F>");
	FluxDiffVSSNRPlotBkg->SetStats(0);
	FluxDiffVSSNRPlotBkg->Draw();

	gFluxBiasVSSNRGraph->SetMarkerStyle(8);
	gFluxBiasVSSNRGraph->SetMarkerSize(1.3);
	gFluxBiasVSSNRGraph->SetLineColor(kBlack);
	gFluxBiasVSSNRGraph->Draw("PZ");

	/*
	if(gApplySelection){
		gFluxBiasVSSNRGraph_sel->SetMarkerStyle(21);
		gFluxBiasVSSNRGraph_sel->SetMarkerSize(1.3);
		gFluxBiasVSSNRGraph_sel->SetLineColor(kRed);
		gFluxBiasVSSNRGraph_sel->Draw("PZ");
	}
	*/

	TLegend* FluxDiffVSSNRPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	FluxDiffVSSNRPlotLegend->SetTextSize(0.045);
	FluxDiffVSSNRPlotLegend->SetTextFont(52);
	FluxDiffVSSNRPlotLegend->AddEntry(gFluxBiasVSSNRGraph,"flux bias (det)","PL"); 
	if(gApplySelection) {
		//FluxDiffVSSNRPlotLegend->AddEntry(gFluxBiasVSSNRGraph_sel,"flux bias (sel)","PL");
	}
	FluxDiffVSSNRPlotLegend->Draw("same");


	//Draw flux reso VS SNR
	TCanvas* FluxDiffResoVSSNRPlot= new TCanvas("FluxDiffResoVSSNRPlot","FluxDiffResoVSSNRPlot");
	FluxDiffResoVSSNRPlot->cd();

	TH2D* FluxDiffResoVSSNRPlotBkg= new TH2D("FluxDiffResoVSSNRPlotBkg","",100,0,100,100,-1,15);
	FluxDiffResoVSSNRPlotBkg->GetXaxis()->SetTitle("S/N");
	FluxDiffResoVSSNRPlotBkg->GetYaxis()->SetTitle("#sigma(#Delta F/F)");
	FluxDiffResoVSSNRPlotBkg->SetStats(0);
	FluxDiffResoVSSNRPlotBkg->Draw();

	gFluxResoVSSNRGraph->SetMarkerStyle(8);
	gFluxResoVSSNRGraph->SetMarkerSize(1.3);
	gFluxResoVSSNRGraph->SetLineColor(kBlack);
	gFluxResoVSSNRGraph->Draw("PZ");

	/*
	if(gApplySelection){
		gFluxResoVSSNRGraph_sel->SetMarkerStyle(21);
		gFluxResoVSSNRGraph_sel->SetMarkerSize(1.3);
		gFluxResoVSSNRGraph_sel->SetLineColor(kRed);
		gFluxResoVSSNRGraph_sel->Draw("PZ");
	}
	*/

	TLegend* FluxDiffResoVSSNRPlotLegend= new TLegend(0.6,0.2,0.8,0.3);
	FluxDiffResoVSSNRPlotLegend->SetTextSize(0.045);
	FluxDiffResoVSSNRPlotLegend->SetTextFont(52);
	FluxDiffResoVSSNRPlotLegend->AddEntry(gFluxResoVSSNRGraph,"flux reso (det)","PL"); 
	if(gApplySelection) {
		//FluxDiffResoVSSNRPlotLegend->AddEntry(gFluxResoVSSNRGraph_sel,"flux reso (sel)","PL");
	}
	FluxDiffResoVSSNRPlotLegend->Draw("same");

	//Draw fit histos
	TCanvas* c= 0;
	TH2D* c_bkg= 0;
	for(size_t i=0;i<gFluxDiffHistos.size();i++){
		TString canvasName= Form("FluxDiffHistoPlot_bin%d",(int)(i+1));
		c= new TCanvas(canvasName,canvasName);
		c->cd();

		gFluxDiffHistos[i]->GetXaxis()->SetTitle("#delta F/F");
		gFluxDiffHistos[i]->SetMarkerStyle(8);
		gFluxDiffHistos[i]->SetMarkerSize(1.3);
		gFluxDiffHistos[i]->GetYaxis()->SetRangeUser(0.,1.5);
		gFluxDiffHistos[i]->Draw("epz");
	
		gFluxDiffFitFcns[i]->Draw("l same");
	}//end loop histos	

	//## Draw eccentricity ratio
	TCanvas* EccentricityRatioPlot= new TCanvas("EccentricityRatioPlot","EccentricityRatioPlot");
	EccentricityRatioPlot->cd();

	gEccentricityRatioHisto->GetXaxis()->SetTitle("eccentricity ratio");
	gEccentricityRatioHisto->Draw("hist");

	//## Draw area ratio ratio
	TCanvas* AreaRatioPlot= new TCanvas("AreaRatioPlot","AreaRatioPlot");
	AreaRatioPlot->cd();

	gAreaRatioHisto->GetXaxis()->SetTitle("area ratio");
	gAreaRatioHisto->Draw("hist");

	//## Draw rotangle vs beam
	TCanvas* RotAnglePlot= new TCanvas("RotAnglePlot","RotAnglePlot");
	RotAnglePlot->cd();

	gRotAngleVSBeamHisto->GetXaxis()->SetTitle("rot angle wrt beam");
	gRotAngleVSBeamHisto->Draw("hist");

}//close Draw()


int FindSourceMatch(const std::vector<SourceInfo>& sources,const std::vector<SourceInfo>& sources_rec)
{
	//## Find map min & max
	auto xrange_it = std::minmax_element( sources.begin(), sources.end(),
                             []( const SourceInfo& l, const SourceInfo& r)
                             {
                                 return l.x < r.x;
                             } );
	auto yrange_it = std::minmax_element( sources.begin(), sources.end(),
                             []( const SourceInfo& l, const SourceInfo& r)
                             {
                                 return l.y < r.y;
                             } );
	double x_min= (*(xrange_it.first)).x;
	double x_max= (*(xrange_it.second)).x;
	double y_min= (*(yrange_it.first)).y;
	double y_max= (*(yrange_it.second)).y;
	
	//cout<<"=== DEBUG ===="<<endl;
	//for(size_t i=0;i<sources.size();i++){
	//	cout<<"(l,b)=("<<sources[i].l<<","<<sources[i].b<<")"<<endl;
	//}	
	//cout<<"=============="<<endl;

	auto l_range_it = std::minmax_element( sources.begin(), sources.end(),
                             []( const SourceInfo& l, const SourceInfo& r)
                             {
                                 return l.l < r.l;
                             } );
	auto b_range_it = std::minmax_element( sources.begin(), sources.end(),
                             []( const SourceInfo& l, const SourceInfo& r)
                             {
                                 return l.b < r.b;
                             } );
	
	double l_min= (*(l_range_it.first)).l;
	double l_max= (*(l_range_it.second)).l;
	double b_min= (*(b_range_it.first)).b;
	double b_max= (*(b_range_it.second)).b;

	//cout<<"INFO: Source ref catalog spatial range: x("<<x_min<<","<<x_max<<"), y("<<y_min<<","<<y_max<<"), l("<<l_min<<","<<l_max<<"), b("<<b_min<<","<<b_max<<")"<<endl;

	
	//## Creating map histos
	double offset= 0.5;
	int nBins= 20;
	if(!gFluxHisto2D){
		gFluxHisto2D= new TH2D("fluxHisto2D","",nBins,x_min-offset,x_max+offset,nBins,y_min-offset,y_max+offset);
		gFluxHisto2D->Sumw2();
	}

	if(!gFluxHisto2D_rec){
		gFluxHisto2D_rec= new TH2D("fluxHisto2D_rec","",nBins,x_min-offset,x_max+offset,nBins,y_min-offset,y_max+offset);
		gFluxHisto2D_rec->Sumw2();
	}

	int nBinsX= (gOpts.galMaxLong-gOpts.galMinLong)/gOpts.galLongStep;
	int nBinsY= (gOpts.galMaxLat-gOpts.galMinLat)/gOpts.galLatStep;

	if(!gSourceCountsHisto2D_galcoord){
		gSourceCountsHisto2D_galcoord= new TH2D("sourceCounts2D_galcoord","",nBinsX,gOpts.galMinLong,gOpts.galMaxLong,nBinsY,gOpts.galMinLat,gOpts.galMaxLat);
		gSourceCountsHisto2D_galcoord->Sumw2();
	}

	if(!gSourceCountsHisto2D_galcoord_rec){
		gSourceCountsHisto2D_galcoord_rec= new TH2D("sourceCounts2D_galcoord_rec","",nBinsX,gOpts.galMinLong,gOpts.galMaxLong,nBinsY,gOpts.galMinLat,gOpts.galMaxLat);
		gSourceCountsHisto2D_galcoord_rec->Sumw2();
	}

	//std::vector<double> gLgFluxBins {-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
	//std::vector<double> gLgRecFluxBins {-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};

	
	TH2D* hh= 0;
	if(gSourceCountsHistos2D_galcoord.empty()){
		for(size_t i=0;i<gLgFluxBins.size();i++)
		{
			TString histoName= Form("sourceCounts2D_galcoord_fluxbin%d",(int)(i+1));
			hh= new TH2D(histoName,"",nBinsX,gOpts.galMinLong,gOpts.galMaxLong,nBinsY,gOpts.galMinLat,gOpts.galMaxLat);
			hh->Sumw2();
			gSourceCountsHistos2D_galcoord.push_back(hh);

			histoName= Form("sourceCounts2D_galcoord_rec_fluxbin%d",(int)(i+1));
			hh= new TH2D(histoName,"",nBinsX,gOpts.galMinLong,gOpts.galMaxLong,nBinsY,gOpts.galMinLat,gOpts.galMaxLat);
			hh->Sumw2();
			gSourceCountsHistos2D_galcoord_rec.push_back(hh);
		}
	}
	if(gSourceCountsHistos2D_galcoord_reliability.empty()){
		for(size_t i=0;i<gLgRecFluxBins.size();i++)
		{
			TString histoName= Form("sourceCounts2D_galcoord_reliability_fluxbin%d",(int)(i+1));
			hh= new TH2D(histoName,"",nBinsX,gOpts.galMinLong,gOpts.galMaxLong,nBinsY,gOpts.galMinLat,gOpts.galMaxLat);
			hh->Sumw2();
			gSourceCountsHistos2D_galcoord_reliability.push_back(hh);

			histoName= Form("sourceCounts2D_galcoord_rec_reliability_fluxbin%d",(int)(i+1));
			hh= new TH2D(histoName,"",nBinsX,gOpts.galMinLong,gOpts.galMaxLong,nBinsY,gOpts.galMinLat,gOpts.galMaxLat);
			hh->Sumw2();
			gSourceCountsHistos2D_galcoord_rec_reliability.push_back(hh);
		}
	}
	
	//## Find if ref source if found in rec catalog
	long int nMatch= 0;
	long int nMatch_sel= 0;
	long int nSources= 0;
	std::vector<bool> isRecSourceFound(sources_rec.size(),false);
	std::vector<bool> isRecSourceFound_sel(sources_rec.size(),false);

	for(size_t i=0;i<sources.size();i++)
	{
		if(i%1000==0) cout<<"--> "<<i+1<<"/"<<sources.size()<<" sources in reference catalog processed..."<<endl;

		std::string name= sources[i].name;
		double x= sources[i].x;
		double y= sources[i].y;
		double ra= sources[i].ra;
		double dec= sources[i].dec;
		double l= sources[i].l;
		double b= sources[i].b;
		double flux= sources[i].flux;
		double lgFlux= log10(flux);

		gIsFound= false;
		gSourceName_true= name;
		gSourcePosX_true= x;
		gSourcePosY_true= y;
		gSourceFlux_true= flux;
		gSourceJ2000PosX_true= ra;
		gSourceJ2000PosY_true= dec;
		gSourceGalPosX_true= l;
		gSourceGalPosY_true= b;

		//Skip source if outside mosaic region
		if(gSelectRegion){
			bool isInside= false;
			for(size_t k=0;k<gRegions.size();k++){
				//if(gRegions[k]->IsPointInsideRegion(x,y)){
				if(gRegions[k]->IsPointInsideRegion(ra,dec)){
					isInside= true;
					break;
				}
			}
			if(!isInside) continue;
		}//close if gSelectRegion

		//Find flux bin
		int fluxBinIndex= -1;
		for(size_t k=0;k<gLgFluxBins.size()-1;k++){
			if(lgFlux>=gLgFluxBins[k] && lgFlux<=gLgFluxBins[k+1]){
				fluxBinIndex= k;
				break;
			}
		}

	
		gFluxHisto->Fill(lgFlux,1);
		//gFluxHisto2D->Fill(x,y,1);
		gFluxHisto2D->Fill(ra,dec,1);
		gFluxGalLatHisto2D->Fill(lgFlux,b,1);
		gSourceCountsHisto2D_galcoord->Fill(l,b,1);
		if(fluxBinIndex!=-1) gSourceCountsHistos2D_galcoord[fluxBinIndex]->Fill(l,b,1);
		
		nSources++;
	
		bool found= false;
		bool selected= false;
		std::vector<std::string> name_list;
		std::vector<double> xrec_list;
		std::vector<double> yrec_list;
		std::vector<double> ra_rec_list;
		std::vector<double> dec_rec_list;
		std::vector<double> l_rec_list;
		std::vector<double> b_rec_list;
		std::vector<double> dx_list;
		std::vector<double> dy_list;
		std::vector<double> dra_list;
		std::vector<double> ddec_list;
		std::vector<double> reclgflux_list;
		std::vector<double> dflux_list;
		std::vector<double> dlgflux_list;
		std::vector<double> snr_list;
		std::vector<double> redchi2_list;
		std::vector<double> Arec_list;
		std::vector<double> eratio_list;
		std::vector<double> arearatio_list;
		std::vector<double> rotangle_list;
		std::vector<double> ellipseArea_list;
		std::vector<double> ellipseArea_wcs_list;
		std::vector<size_t> recindex_list;	


		for(size_t j=0;j<sources_rec.size();j++)
		{	
			std::string name_rec= sources_rec[j].name;
			double x_rec= sources_rec[j].x;
			double y_rec= sources_rec[j].y;
			double ra_rec= sources_rec[j].ra;
			double dec_rec= sources_rec[j].dec;
			double l_rec= sources_rec[j].l;
			double b_rec= sources_rec[j].b;
			double flux_rec= sources_rec[j].flux;
			double lgflux_rec= log10(flux_rec);
			double rms_rec= sources_rec[j].rms;
			double snr_rec= sources_rec[j].snr;
			double chi2_rec= sources_rec[j].chi2;
			double ndf_rec= sources_rec[j].ndf;
			double redchi2_rec= chi2_rec/ndf_rec;	
			double chi2prob= TMath::Prob(chi2_rec,ndf_rec);		
			
			double bmaj= sources_rec[j].bmaj_pix;//in pix
			double bmin= sources_rec[j].bmin_pix;	//in pix
			//double bmaj_wcs= sources_rec[j].bmaj_wcs;	
			//double bmin_wcs= sources_rec[j].bmin_wcs;	
			double bmaj_wcs= bmaj*gMapPixSize;//converted in arcsec
			double bmin_wcs= bmin*gMapPixSize;//converted in arcsec
	 		double ellipseArea= MathUtils::ComputeEllipseArea(bmaj,bmin);
			double ellipseArea_wcs= MathUtils::ComputeEllipseArea(bmaj_wcs,bmin_wcs);
			double ellipseAreaRatio_wcs= ellipseArea_wcs/gBeamArea_wcs;

			double A_rec= sources_rec[j].A;
			double eccentricityRatio= sources_rec[j].eccentricityRatio;
			//double areaRatio= sources_rec[j].areaRatio;
			double areaRatio= ellipseAreaRatio_wcs;
			double rotAngle= sources_rec[j].rotAngleVSBeam;

			//Check if sources match given position
			double dx= (x_rec-x);//in pixels
			double dy= (y_rec-y);//in pixels
			double dra= (ra_rec-ra)*3600.;//in arcsec
			double ddec= (dec_rec-dec)*3600.;//in arcsec
			/*
			double dx= 0;		
			double dy= 0;
			if(gUseWCS){
				//dx= (x_rec-x)*3600.;//in arcsec
				//dy= (y_rec-y)*3600.;//in arcsec
				dx= (ra_rec-ra)*3600.;//in arcsec
				dy= (dec_rec-dec)*3600.;//in arcsec
			}
			else{
				dx= (x_rec-x);//in pixels
				dy= (y_rec-y);//in pixels
			}
			*/

			double dpos= sqrt(dx*dx + dy*dy);
			double dpos_wcs= sqrt(dra*dra + ddec*ddec);
			double dflux= flux_rec-flux;
			double dflux_rel= flux_rec/flux-1;
			double dlgflux= log10(flux_rec) - log10(flux);

			//Skip source if outside mosaic region
			if(gSelectRegion){
				bool isInside= false;
				for(size_t k=0;k<gRegions.size();k++){
					//if(gRegions[k]->IsPointInsideRegion(x_rec,y_rec)){
					if(gRegions[k]->IsPointInsideRegion(ra_rec,dec_rec)){
						isInside= true;
						break;
					}
				}
				if(!isInside) continue;
			}//close if gSelectRegion

			//if(dpos>gMatchPosThr) continue;
			if(dpos_wcs>gMatchPosThr) continue;

			//cout<<"Match source (name="<<name<<", x="<<x<<", y="<<y<<") with source (name="<<name_rec<<", x="<<x_rec<<", y="<<y_rec<<"), dx="<<dx<<", dy="<<dy<<", snr="<<snr_rec<<endl;
			found= true;

			//Fill vectors for matched sources
			recindex_list.push_back(j);
			name_list.push_back(name_rec);
			xrec_list.push_back(x_rec);
			yrec_list.push_back(y_rec);
			ra_rec_list.push_back(ra_rec);
			dec_rec_list.push_back(dec_rec);
			l_rec_list.push_back(l_rec);
			b_rec_list.push_back(b_rec);
			dx_list.push_back(dx);
			dy_list.push_back(dy);
			dra_list.push_back(dra);
			ddec_list.push_back(ddec);	
			reclgflux_list.push_back(lgflux_rec);
			dflux_list.push_back(dflux_rel);
			dlgflux_list.push_back(dlgflux);
			snr_list.push_back(snr_rec);
			redchi2_list.push_back(redchi2_rec);
			Arec_list.push_back(A_rec);

			eratio_list.push_back(eccentricityRatio);
			arearatio_list.push_back(areaRatio);
			rotangle_list.push_back(rotAngle);

		}//end loop rec sources

		

		//Fill histos & vector
		gIsFound= false;
		gSourceName_rec= "";
		gSourcePosX_rec= -999;
		gSourcePosY_rec= -999;
		gSourceFlux_rec= -999;
		gSourceJ2000PosX_rec= -999;
		gSourceJ2000PosY_rec= -999;
		gSourceGalPosX_rec= -999;
		gSourceGalPosY_rec= -999;

		if(found){
			nMatch++;
			
			//Check for multiple detections and select the closest
			int nFound= (int)(dx_list.size());
			std::string name_rec= name_list[0];
			double x_rec= xrec_list[0];
			double y_rec= yrec_list[0];
			double ra_rec= ra_rec_list[0];
			double dec_rec= dec_rec_list[0];
			double l_rec= l_rec_list[0];
			double b_rec= b_rec_list[0];
			double dx= dx_list[0];
			double dy= dy_list[0];
			double dra= dra_list[0];
			double ddec= ddec_list[0];
			double dpos= sqrt(dx*dx + dy*dy); 
			double dpos_wcs= sqrt(dra*dra + ddec*ddec); 	
			double lgFlux_rec= reclgflux_list[0];
			double dflux_rel= dflux_list[0];
			double dlgflux= dlgflux_list[0];
			double snr_rec= snr_list[0];
			double redchi2_rec= redchi2_list[0];
			double A_rec= Arec_list[0];
			double eccentricityRatio= eratio_list[0];
			double areaRatio= arearatio_list[0];
			double rotAngle= rotangle_list[0];
			size_t recIndex= recindex_list[0];

			if(nFound>1){
				int best_index= 0;
				for(int k=1;k<nFound;k++){
					/*
					double dpos_k= sqrt(dx_list[k]*dx_list[k] + dy_list[k]*dy_list[k]);
					if(dpos_k<dpos){
						best_index= k;
						dpos= dpos_k;
					}
					*/
					double dpos_k= sqrt(dra_list[k]*dra_list[k] + ddec_list[k]*ddec_list[k]);
					if(dpos_k<dpos_wcs){
						best_index= k;
						dpos_wcs= dpos_k;
					}
				}//end loop match

				name_rec= name_list[best_index];
				x_rec= xrec_list[best_index];
				y_rec= yrec_list[best_index];
				ra_rec= ra_rec_list[best_index];
				dec_rec= dec_rec_list[best_index];
				l_rec= l_rec_list[best_index];
				b_rec= b_rec_list[best_index];
				dx= dx_list[best_index];
				dy= dy_list[best_index];
				dpos= sqrt(dx*dx + dy*dy); 
				dra= dra_list[best_index];
				ddec= ddec_list[best_index];
				dpos_wcs= sqrt(dra*dra + ddec*ddec); 
				lgFlux_rec= reclgflux_list[best_index];
				dflux_rel= dflux_list[best_index];
				dlgflux= dlgflux_list[best_index];
				snr_rec= snr_list[best_index];
				redchi2_rec= redchi2_list[best_index];
				A_rec= Arec_list[best_index];
				eccentricityRatio= eratio_list[best_index];
				areaRatio= arearatio_list[best_index];
				rotAngle= rotangle_list[best_index];
				recIndex= recindex_list[best_index];
	
				cout<<"INFO: #"<<nFound<<" match found, match no. "<<best_index<<" selected as the best..."<<endl;

			}//close nFound>1

			gIsFound= true;
			gSourceName_rec= name_rec;
			gSourcePosX_rec= x_rec;
			gSourcePosY_rec= y_rec;
			gSourceFlux_rec= pow(10,lgFlux_rec);
			gSourceJ2000PosX_rec= ra_rec;
			gSourceJ2000PosY_rec= dec_rec;
			gSourceGalPosX_rec= l_rec;
			gSourceGalPosY_rec= b_rec;

			isRecSourceFound[recIndex]= true;

			//Apply selection
			if(gApplySelection){
				selected= (
					redchi2_rec<gRedChi2Thr &&
					A_rec>0
				);
			}
			else{
				nMatch_sel++;
				selected= true;
			}

			//Fill histos
			//gPosXDiffHisto->Fill(lgFlux,dx);
			//gPosYDiffHisto->Fill(lgFlux,dy);
			gPosXDiffHisto->Fill(lgFlux,dra);
			gPosYDiffHisto->Fill(lgFlux,ddec);
			gFluxDiffHisto->Fill(lgFlux,dflux_rel);
			gFluxDiffVSRecFluxHisto->Fill(lgFlux_rec,dflux_rel);
			gFluxHisto_rec->Fill(lgFlux,1);
			//gFluxHisto2D_rec->Fill(x,y,1);
			gFluxHisto2D_rec->Fill(ra,dec,1);
			gFluxGalLatHisto2D_rec->Fill(lgFlux,b,1);
			gSourceCountsHisto2D_galcoord_rec->Fill(l,b,1);
			if(fluxBinIndex!=-1) gSourceCountsHistos2D_galcoord_rec[fluxBinIndex]->Fill(l,b,1);
			gSNRHisto_rec->Fill(snr_rec);
			gEccentricityRatioHisto->Fill(eccentricityRatio);
			gAreaRatioHisto->Fill(areaRatio);
			gRotAngleVSBeamHisto->Fill(rotAngle);

			if(selected){
				gFluxHisto_sel->Fill(lgFlux,1);
				nMatch_sel++;
				isRecSourceFound_sel[recIndex]= true;
			}

			int binId= gFluxDiffHisto->FindBin(lgFlux);
			if(binId>0 && binId<gFluxBinData.size()) {
				gFluxBinData[binId-1].push_back(lgFlux);
				gFluxDiffBinData[binId-1].push_back(dflux_rel);
				gLgFluxDiffBinData[binId-1].push_back(dlgflux);
				//gPosXDiffBinData[binId-1].push_back(dx);
				//gPosYDiffBinData[binId-1].push_back(dy);
				gPosXDiffBinData[binId-1].push_back(dra);
				gPosYDiffBinData[binId-1].push_back(ddec);

				if(selected){
					gFluxBinData_sel[binId-1].push_back(lgFlux);
					gFluxDiffBinData_sel[binId-1].push_back(dflux_rel);
					gLgFluxDiffBinData_sel[binId-1].push_back(dlgflux);
					//gPosXDiffBinData_sel[binId-1].push_back(dx);
					//gPosYDiffBinData_sel[binId-1].push_back(dy);	
					gPosXDiffBinData_sel[binId-1].push_back(dra);
					gPosYDiffBinData_sel[binId-1].push_back(ddec);
				}
			}//close if

			int binId_rec= gFluxDiffVSRecFluxHisto->FindBin(lgFlux_rec);
			if(binId_rec>0 && binId_rec<gRecFluxBinData.size()) {
				gRecFluxBinData[binId_rec-1].push_back(lgFlux_rec);
				gLgFluxDiffVSRecFluxBinData[binId_rec-1].push_back(dlgflux);
				if(selected){
					gRecFluxBinData_sel[binId_rec-1].push_back(lgFlux_rec);
					gLgFluxDiffVSRecFluxBinData_sel[binId_rec-1].push_back(dlgflux);
				}
			}	

			int binId_snr= gSNRHisto_rec->FindBin(snr_rec);
			if(binId_snr>0 && binId_snr<gSNRBinData.size()) {
				gSNRBinData[binId_snr-1].push_back(snr_rec);
				gFluxDiffVSSNRBinData[binId_snr-1].push_back(dflux_rel);
				//gPosXDiffVSSNRBinData[binId_snr-1].push_back(dx);
				//gPosYDiffVSSNRBinData[binId_snr-1].push_back(dy);
				gPosXDiffVSSNRBinData[binId_snr-1].push_back(dra);
				gPosYDiffVSSNRBinData[binId_snr-1].push_back(ddec);
			}//close if

		}//close if found
		else{
			if(flux>1.e-1){
				cout<<"Bright source not found (name="<<name<<", pos("<<x<<","<<y<<"), skypos("<<ra<<","<<dec<<"), flux="<<flux<<endl;
			}			
		}

		gMatchDataTree->Fill();

	}//end loop ref sources

	cout<<"INFO: #"<<nMatch<<"/"<<nSources<<" source match (#"<<nMatch_sel<<" selected)"<<endl;


	//======================================
	//==      RELIABILITY
	//======================================	
	int nRecSources_true= 0;
	int nSelSources_true= 0;

	for(size_t j=0;j<sources_rec.size();j++)
	{	
		if(j%100==0) cout<<"--> "<<j+1<<"/"<<sources_rec.size()<<" sources in rec catalog processed..."<<endl;
		double x_rec= sources_rec[j].x;
		double y_rec= sources_rec[j].y;
		double ra_rec= sources_rec[j].ra;
		double dec_rec= sources_rec[j].dec;
		double flux_rec= sources_rec[j].flux;
		double lgFlux_rec= log10(flux_rec);
		double l_rec= sources_rec[j].l;
		double b_rec= sources_rec[j].b;

		//Skip source if outside mosaic region
		if(gSelectRegion){
			bool isInside= false;
			for(size_t k=0;k<gRegions.size();k++){
				//if(gRegions[k]->IsPointInsideRegion(x_rec,y_rec)){
				if(gRegions[k]->IsPointInsideRegion(ra_rec,dec_rec)){
					isInside= true;
					break;
				}
			}
			if(!isInside) continue;
		}//close if gSelectRegion

		//Find flux bin
		int fluxBinIndex= -1;
		for(size_t k=0;k<gLgRecFluxBins.size()-1;k++){
			if(lgFlux_rec>=gLgRecFluxBins[k] && lgFlux_rec<=gLgRecFluxBins[k+1]){
				fluxBinIndex= k;
				break;
			}
		}

		gFluxHisto_reliability->Fill(lgFlux_rec,1);
		gFluxGalLatHisto2D_reliability->Fill(lgFlux_rec,b_rec,1);
		if(fluxBinIndex!=-1) gSourceCountsHistos2D_galcoord_reliability[fluxBinIndex]->Fill(l_rec,b_rec,1);
		gFluxHisto_fdr->Fill(lgFlux_rec,1);
		if(isRecSourceFound[j]) {
			gFluxHisto_rec_reliability->Fill(lgFlux_rec,1);
			gFluxGalLatHisto2D_rec_reliability->Fill(lgFlux_rec,b_rec,1);
			if(fluxBinIndex!=-1) gSourceCountsHistos2D_galcoord_rec_reliability[fluxBinIndex]->Fill(l_rec,b_rec,1);
			nRecSources_true++;
		}
		else {
			gFluxHisto_rec_fdr->Fill(lgFlux_rec,1);
		}

		if(isRecSourceFound_sel[j]) {
			gFluxHisto_sel_reliability->Fill(lgFlux_rec,1);
			nSelSources_true++;
		}
		else {
			gFluxHisto_sel_fdr->Fill(lgFlux_rec,1);
		}

	}//end loop rec sources

	cout<<"INFO: #"<<nRecSources_true<<"/"<<sources_rec.size()<<" good detection (#"<<nSelSources_true<<" selected)"<<endl;
	
	return 0;

}//close FindSourceMatch()


int AnalyzeData()
{
	//Init histo & fcn 
	TH1D* h= 0;
	TF1* fcn= 0;
	int nHistoBins= 40;
	int pntCounter= 0;

	//Loop over data and find stats
	for(size_t i=0;i<gFluxBinData.size();i++){
		int nEntries= (int)(gFluxBinData[i].size());
		cout<<"INFO: #"<<nEntries<<" data available in flux bin "<<i+1<<" ..."<<endl;

		if(nEntries<=0) continue;
	
		//Compute x pos stats
		Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(gPosXDiffBinData[i]);
		double posx_mean= 0;
		double posx_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,gPosXDiffBinData[i]);
		double posx_median= stats_posx.median;
		double posx_err= posx_rms/sqrt(nEntries);
		double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
		double posx_iqr= stats_posx.Q3-stats_posx.Q1;
		double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
		double posx_iqr_half= posx_iqr/2.;
		double posx_iqr_half_err= 0.5*posx_iqr_err;

		//Compute y pos stats
		Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(gPosYDiffBinData[i]);
		double posy_mean= 0;
		double posy_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,gPosYDiffBinData[i]);
		double posy_median= stats_posy.median;
		double posy_err= posy_rms/sqrt(nEntries);
		double posy_median_err= 1.253*posy_rms/sqrt(nEntries);
		double posy_iqr= stats_posy.Q3-stats_posy.Q1;
		double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
		double posy_iqr_half= posy_iqr/2.;
		double posy_iqr_half_err= 0.5*posy_iqr_err;

		//Compute fluxdiff stats	
		double fluxdiff_mean= 0;
		double fluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(fluxdiff_mean,fluxdiff_rms,gFluxDiffBinData[i]);
		Caesar::BoxStats<double> stats_fluxdiff= StatsUtils::ComputeBoxStats(gFluxDiffBinData[i]);
		double fluxdiff_median= stats_fluxdiff.median;
		double fluxdiff_err= fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_median_err= 1.253*fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_iqr= stats_fluxdiff.Q3-stats_fluxdiff.Q1;
		double fluxdiff_iqr_err= 1.573*fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_iqr_half= fluxdiff_iqr/2.;
		double fluxdiff_iqr_half_err= 0.5*fluxdiff_iqr_err;
		double fluxdiff_err_low= fluxdiff_median-stats_fluxdiff.Q1;
		double fluxdiff_err_up= stats_fluxdiff.Q3-fluxdiff_median;
		double fluxdiff_min= stats_fluxdiff.minVal;
		double fluxdiff_max= stats_fluxdiff.maxVal;
		double dfluxdiff= fabs(fluxdiff_max-fluxdiff_min);

		//Compute flux stats
		double flux_mean= 0;
		double flux_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(flux_mean,flux_rms,gFluxBinData[i]);

		cout<<"INFO: bin "<<i+1<<" (lgFlux="<<flux_mean<<") flux diff stats {mean="<<fluxdiff_mean<<", rms="<<fluxdiff_rms<<", min/max="<<fluxdiff_min<<"/"<<fluxdiff_max<<", median="<<fluxdiff_median<<", IQR="<<fluxdiff_iqr<<"}"<<endl;

		//Compute lgfluxdiff stats	
		double lgfluxdiff_mean= 0;
		double lgfluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(lgfluxdiff_mean,lgfluxdiff_rms,gLgFluxDiffBinData[i]);
		Caesar::BoxStats<double> stats_lgfluxdiff= StatsUtils::ComputeBoxStats(gLgFluxDiffBinData[i]);
		double lgfluxdiff_median= stats_lgfluxdiff.median;
		double lgfluxdiff_err= lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_median_err= 1.253*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr= stats_lgfluxdiff.Q3-stats_lgfluxdiff.Q1;
		double lgfluxdiff_iqr_err= 1.573*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr_half= lgfluxdiff_iqr/2.;
		double lgfluxdiff_iqr_half_err= 0.5*lgfluxdiff_iqr_err;
		double lgfluxdiff_err_low= lgfluxdiff_median-stats_lgfluxdiff.Q1;
		double lgfluxdiff_err_up= stats_lgfluxdiff.Q3-lgfluxdiff_median;
		double lgfluxdiff_min= stats_lgfluxdiff.minVal;
		double lgfluxdiff_max= stats_lgfluxdiff.maxVal;
		double dlgfluxdiff= fabs(lgfluxdiff_max-lgfluxdiff_min);
		
		//Fill graph
		//gFluxBiasGraph->SetPoint(pntCounter,flux_mean,fluxdiff_mean);
		//gFluxBiasGraph->SetPointError(pntCounter,0,fluxdiff_err);
		gFluxBiasGraph->SetPoint(pntCounter,flux_mean,fluxdiff_median);
		gFluxBiasGraph->SetPointError(pntCounter,0,fluxdiff_median_err);
		gFluxResoGraph->SetPoint(pntCounter,flux_mean,fluxdiff_iqr_half);
		gFluxResoGraph->SetPointError(pntCounter,0,fluxdiff_iqr_half_err);

		gLgFluxBiasGraph->SetPoint(pntCounter,flux_mean,lgfluxdiff_median);
		gLgFluxBiasGraph->SetPointError(pntCounter,0,lgfluxdiff_median_err);
		gLgFluxResoGraph->SetPoint(pntCounter,flux_mean,lgfluxdiff_iqr_half);
		gLgFluxResoGraph->SetPointError(pntCounter,0,lgfluxdiff_iqr_half_err);
	
		gPosXBiasGraph->SetPoint(pntCounter,flux_mean,posx_median);
		gPosXBiasGraph->SetPointError(pntCounter,0,posx_median_err);
		gPosYBiasGraph->SetPoint(pntCounter,flux_mean,posy_median);
		gPosYBiasGraph->SetPointError(pntCounter,0,posy_median_err);

		gPosXResoGraph->SetPoint(pntCounter,flux_mean,posx_iqr_half);
		gPosXResoGraph->SetPointError(pntCounter,0,posx_iqr_half_err);
	
		gPosYResoGraph->SetPoint(pntCounter,flux_mean,posy_iqr_half);
		gPosYResoGraph->SetPointError(pntCounter,0,posy_iqr_half_err);

		//Fill histos
		TString histoName= Form("fluxDiffHist_bin%d",(int)(i+1));
		h= new TH1D(histoName,histoName,nHistoBins,fluxdiff_min-0.1*dfluxdiff,fluxdiff_max+0.1*dfluxdiff);
		h->Sumw2();
		for(size_t j=0;j<gFluxDiffBinData[i].size();j++){
			double fluxDiff= gFluxDiffBinData[i][j];
			h->Fill(fluxDiff);	
		}		
		h->Scale(1./h->Integral());
		
		//Create fit fcn
		TString fcnName= Form("fluxDiffFitFcn_bin%d",(int)(i+1));
		fcn= new TF1(fcnName,"[0]*TMath::Gaus(x,[1],[2])",fluxdiff_min-0.1*dfluxdiff,fluxdiff_max+0.1*dfluxdiff);
		fcn->SetParNames("A","mean","sigma");
		fcn->SetParameters(1,fluxdiff_mean,fluxdiff_rms);
		fcn->SetNpx(500);
		
		//Fit histos
		h->Fit(fcn,"RS");

		//Add histos & fcn to list
		gFluxDiffHistos.push_back(h);
		gFluxDiffFitFcns.push_back(fcn);


		pntCounter++;

	}//end loop 


	
	//Loop over selected data and find stats
	pntCounter= 0;

	for(size_t i=0;i<gFluxBinData_sel.size();i++){
		int nEntries= (int)(gFluxBinData_sel[i].size());
		cout<<"INFO: #"<<nEntries<<" data available in flux bin "<<i+1<<" ..."<<endl;

		if(nEntries<=0) continue;
	
		//Compute x pos stats
		Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(gPosXDiffBinData_sel[i]);
		double posx_mean= 0;
		double posx_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,gPosXDiffBinData_sel[i]);
		double posx_median= stats_posx.median;
		double posx_err= posx_rms/sqrt(nEntries);
		double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
		double posx_iqr= stats_posx.Q3-stats_posx.Q1;
		double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
		double posx_iqr_half= posx_iqr/2.;
		double posx_iqr_half_err= 0.5*posx_iqr_err;

		//Compute y pos stats
		Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(gPosYDiffBinData_sel[i]);
		double posy_mean= 0;
		double posy_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,gPosYDiffBinData_sel[i]);
		double posy_median= stats_posy.median;
		double posy_err= posy_rms/sqrt(nEntries);
		double posy_median_err= 1.253*posy_rms/sqrt(nEntries);
		double posy_iqr= stats_posy.Q3-stats_posy.Q1;
		double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
		double posy_iqr_half= posy_iqr/2.;
		double posy_iqr_half_err= 0.5*posy_iqr_err;

		//Compute fluxdiff stats	
		double fluxdiff_mean= 0;
		double fluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(fluxdiff_mean,fluxdiff_rms,gFluxDiffBinData_sel[i]);
		Caesar::BoxStats<double> stats_fluxdiff= StatsUtils::ComputeBoxStats(gFluxDiffBinData_sel[i]);
		double fluxdiff_median= stats_fluxdiff.median;
		double fluxdiff_err= fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_median_err= 1.253*fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_iqr= stats_fluxdiff.Q3-stats_fluxdiff.Q1;
		double fluxdiff_iqr_err= 1.573*fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_iqr_half= fluxdiff_iqr/2.;
		double fluxdiff_iqr_half_err= 0.5*fluxdiff_iqr_err;
		double fluxdiff_err_low= fluxdiff_median-stats_fluxdiff.Q1;
		double fluxdiff_err_up= stats_fluxdiff.Q3-fluxdiff_median;
		double fluxdiff_min= stats_fluxdiff.minVal;
		double fluxdiff_max= stats_fluxdiff.maxVal;
		double dfluxdiff= fabs(fluxdiff_max-fluxdiff_min);

		//Compute flux stats
		double flux_mean= 0;
		double flux_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(flux_mean,flux_rms,gFluxBinData_sel[i]);

		cout<<"INFO: bin "<<i+1<<" (lgFlux="<<flux_mean<<") flux diff stats {mean="<<fluxdiff_mean<<", rms="<<fluxdiff_rms<<", min/max="<<fluxdiff_min<<"/"<<fluxdiff_max<<", median="<<fluxdiff_median<<", IQR="<<fluxdiff_iqr<<"}"<<endl;

		//Compute lgfluxdiff stats	
		double lgfluxdiff_mean= 0;
		double lgfluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(lgfluxdiff_mean,lgfluxdiff_rms,gLgFluxDiffBinData_sel[i]);
		Caesar::BoxStats<double> stats_lgfluxdiff= StatsUtils::ComputeBoxStats(gLgFluxDiffBinData_sel[i]);
		double lgfluxdiff_median= stats_lgfluxdiff.median;
		double lgfluxdiff_err= lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_median_err= 1.253*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr= stats_lgfluxdiff.Q3-stats_lgfluxdiff.Q1;
		double lgfluxdiff_iqr_err= 1.573*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr_half= lgfluxdiff_iqr/2.;
		double lgfluxdiff_iqr_half_err= 0.5*lgfluxdiff_iqr_err;
		double lgfluxdiff_err_low= lgfluxdiff_median-stats_lgfluxdiff.Q1;
		double lgfluxdiff_err_up= stats_lgfluxdiff.Q3-lgfluxdiff_median;
		double lgfluxdiff_min= stats_lgfluxdiff.minVal;
		double lgfluxdiff_max= stats_lgfluxdiff.maxVal;
		double dlgfluxdiff= fabs(lgfluxdiff_max-lgfluxdiff_min);

		//Fill graph
		//gFluxBiasGraph_sel->SetPoint(pntCounter,flux_mean,fluxdiff_mean);
		//gFluxBiasGraph_sel->SetPointError(pntCounter,0,fluxdiff_err);
		gFluxBiasGraph_sel->SetPoint(pntCounter,flux_mean,fluxdiff_median);
		gFluxBiasGraph_sel->SetPointError(pntCounter,0,fluxdiff_median_err);
		gFluxResoGraph_sel->SetPoint(pntCounter,flux_mean,fluxdiff_iqr_half);
		gFluxResoGraph_sel->SetPointError(pntCounter,0,fluxdiff_iqr_half_err);

		gLgFluxBiasGraph_sel->SetPoint(pntCounter,flux_mean,lgfluxdiff_median);
		gLgFluxBiasGraph_sel->SetPointError(pntCounter,0,lgfluxdiff_median_err);
		gLgFluxResoGraph_sel->SetPoint(pntCounter,flux_mean,lgfluxdiff_iqr_half);
		gLgFluxResoGraph_sel->SetPointError(pntCounter,0,lgfluxdiff_iqr_half_err);
	
		gPosXBiasGraph_sel->SetPoint(pntCounter,flux_mean,posx_median);
		gPosXBiasGraph_sel->SetPointError(pntCounter,0,posx_median_err);
		gPosYBiasGraph_sel->SetPoint(pntCounter,flux_mean,posy_median);
		gPosYBiasGraph_sel->SetPointError(pntCounter,0,posy_median_err);

		gPosXResoGraph_sel->SetPoint(pntCounter,flux_mean,posx_iqr_half);
		gPosXResoGraph_sel->SetPointError(pntCounter,0,posx_iqr_half_err);
	
		gPosYResoGraph_sel->SetPoint(pntCounter,flux_mean,posy_iqr_half);
		gPosYResoGraph_sel->SetPointError(pntCounter,0,posy_iqr_half_err);

		pntCounter++;

	}//end loop 



	//Loop over rec flux data	
	pntCounter= 0;

	for(size_t i=0;i<gRecFluxBinData.size();i++){
		int nEntries= (int)(gRecFluxBinData[i].size());
		cout<<"INFO: #"<<nEntries<<" data available in flux bin "<<i+1<<" ..."<<endl;

		if(nEntries<=0) continue;
	
		//Compute flux stats
		double flux_mean= 0;
		double flux_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(flux_mean,flux_rms,gRecFluxBinData[i]);

		//Compute lgfluxdiff stats	
		double lgfluxdiff_mean= 0;
		double lgfluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(lgfluxdiff_mean,lgfluxdiff_rms,gLgFluxDiffVSRecFluxBinData[i]);
		Caesar::BoxStats<double> stats_lgfluxdiff= StatsUtils::ComputeBoxStats(gLgFluxDiffVSRecFluxBinData[i]);
		double lgfluxdiff_median= stats_lgfluxdiff.median;
		double lgfluxdiff_err= lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_median_err= 1.253*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr= stats_lgfluxdiff.Q3-stats_lgfluxdiff.Q1;
		double lgfluxdiff_iqr_err= 1.573*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr_half= lgfluxdiff_iqr/2.;
		double lgfluxdiff_iqr_half_err= 0.5*lgfluxdiff_iqr_err;
		double lgfluxdiff_err_low= lgfluxdiff_median-stats_lgfluxdiff.Q1;
		double lgfluxdiff_err_up= stats_lgfluxdiff.Q3-lgfluxdiff_median;
		double lgfluxdiff_min= stats_lgfluxdiff.minVal;
		double lgfluxdiff_max= stats_lgfluxdiff.maxVal;
		double dlgfluxdiff= fabs(lgfluxdiff_max-lgfluxdiff_min);
		
		//Fill graph
		gLgFluxBiasVSRecFluxGraph->SetPoint(pntCounter,flux_mean,lgfluxdiff_median);
		gLgFluxBiasVSRecFluxGraph->SetPointError(pntCounter,0,lgfluxdiff_median_err);
		gLgFluxResoVSRecFluxGraph->SetPoint(pntCounter,flux_mean,lgfluxdiff_iqr_half);
		gLgFluxResoVSRecFluxGraph->SetPointError(pntCounter,0,lgfluxdiff_iqr_half_err);
		
		pntCounter++;

	}//end loop rec bin data

	//Loop over rec flux selected data	
	pntCounter= 0;

	for(size_t i=0;i<gRecFluxBinData_sel.size();i++){
		int nEntries= (int)(gRecFluxBinData_sel[i].size());
		cout<<"INFO: #"<<nEntries<<" data available in flux bin "<<i+1<<" ..."<<endl;

		if(nEntries<=0) continue;
	
		//Compute flux stats
		double flux_mean= 0;
		double flux_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(flux_mean,flux_rms,gRecFluxBinData_sel[i]);

		//Compute lgfluxdiff stats	
		double lgfluxdiff_mean= 0;
		double lgfluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(lgfluxdiff_mean,lgfluxdiff_rms,gLgFluxDiffVSRecFluxBinData_sel[i]);
		Caesar::BoxStats<double> stats_lgfluxdiff= StatsUtils::ComputeBoxStats(gLgFluxDiffVSRecFluxBinData_sel[i]);
		double lgfluxdiff_median= stats_lgfluxdiff.median;
		double lgfluxdiff_err= lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_median_err= 1.253*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr= stats_lgfluxdiff.Q3-stats_lgfluxdiff.Q1;
		double lgfluxdiff_iqr_err= 1.573*lgfluxdiff_rms/sqrt(nEntries);
		double lgfluxdiff_iqr_half= lgfluxdiff_iqr/2.;
		double lgfluxdiff_iqr_half_err= 0.5*lgfluxdiff_iqr_err;
		double lgfluxdiff_err_low= lgfluxdiff_median-stats_lgfluxdiff.Q1;
		double lgfluxdiff_err_up= stats_lgfluxdiff.Q3-lgfluxdiff_median;
		double lgfluxdiff_min= stats_lgfluxdiff.minVal;
		double lgfluxdiff_max= stats_lgfluxdiff.maxVal;
		double dlgfluxdiff= fabs(lgfluxdiff_max-lgfluxdiff_min);

		//Fill graph
		gLgFluxBiasVSRecFluxGraph_sel->SetPoint(pntCounter,flux_mean,lgfluxdiff_median);
		gLgFluxBiasVSRecFluxGraph_sel->SetPointError(pntCounter,0,lgfluxdiff_median_err);
		gLgFluxResoVSRecFluxGraph_sel->SetPoint(pntCounter,flux_mean,lgfluxdiff_iqr_half);
		gLgFluxResoVSRecFluxGraph_sel->SetPointError(pntCounter,0,lgfluxdiff_iqr_half_err);
		
	
		pntCounter++;

	}//end loop rec bin selected data




	//Loop over SNR data and find stats
	pntCounter= 0;

	for(size_t i=0;i<gSNRBinData.size();i++){
		int nEntries= (int)(gSNRBinData[i].size());
		cout<<"INFO: #"<<nEntries<<" data available in SNR bin "<<i+1<<" ..."<<endl;

		if(nEntries<=0) continue;
	
		//Compute x pos stats
		Caesar::BoxStats<double> stats_posx= StatsUtils::ComputeBoxStats(gPosXDiffVSSNRBinData[i]);
		double posx_mean= 0;
		double posx_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(posx_mean,posx_rms,gPosXDiffVSSNRBinData[i]);
		double posx_median= stats_posx.median;
		double posx_err= posx_rms/sqrt(nEntries);
		double posx_median_err= 1.253*posx_rms/sqrt(nEntries);
		double posx_iqr= stats_posx.Q3-stats_posx.Q1;
		double posx_iqr_err= 1.573*posx_rms/sqrt(nEntries);
		double posx_iqr_half= posx_iqr/2.;
		double posx_iqr_half_err= 0.5*posx_iqr_err;

		//Compute y pos stats
		Caesar::BoxStats<double> stats_posy= StatsUtils::ComputeBoxStats(gPosYDiffVSSNRBinData[i]);
		double posy_mean= 0;
		double posy_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(posy_mean,posy_rms,gPosYDiffVSSNRBinData[i]);
		double posy_median= stats_posy.median;
		double posy_err= posy_rms/sqrt(nEntries);
		double posy_median_err= 1.253*posy_rms/sqrt(nEntries);
		double posy_iqr= stats_posy.Q3-stats_posy.Q1;
		double posy_iqr_err= 1.573*posy_rms/sqrt(nEntries);
		double posy_iqr_half= posy_iqr/2.;
		double posy_iqr_half_err= 0.5*posy_iqr_err;

		//Compute fluxdiff stats	
		double fluxdiff_mean= 0;
		double fluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(fluxdiff_mean,fluxdiff_rms,gFluxDiffVSSNRBinData[i]);
		Caesar::BoxStats<double> stats_fluxdiff= StatsUtils::ComputeBoxStats(gFluxDiffVSSNRBinData[i]);
		double fluxdiff_median= stats_fluxdiff.median;
		double fluxdiff_err= fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_median_err= 1.253*fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_iqr= stats_fluxdiff.Q3-stats_fluxdiff.Q1;
		double fluxdiff_iqr_err= 1.573*fluxdiff_rms/sqrt(nEntries);
		double fluxdiff_iqr_half= fluxdiff_iqr/2.;
		double fluxdiff_iqr_half_err= 0.5*fluxdiff_iqr_err;
		double fluxdiff_err_low= fluxdiff_median-stats_fluxdiff.Q1;
		double fluxdiff_err_up= stats_fluxdiff.Q3-fluxdiff_median;
		double fluxdiff_min= stats_fluxdiff.minVal;
		double fluxdiff_max= stats_fluxdiff.maxVal;
		double dfluxdiff= fabs(fluxdiff_max-fluxdiff_min);

		//Compute SNR stats
		double snr_mean= 0;
		double snr_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(snr_mean,snr_rms,gSNRBinData[i]);

		cout<<"INFO: bin "<<i+1<<" (snr="<<snr_mean<<") flux diff stats {mean="<<fluxdiff_mean<<", rms="<<fluxdiff_rms<<", min/max="<<fluxdiff_min<<"/"<<fluxdiff_max<<", median="<<fluxdiff_median<<", IQR="<<fluxdiff_iqr<<"}"<<endl;

		//Fill graph
		//gFluxBiasVSSNRGraph->SetPoint(pntCounter,snr_mean,fluxdiff_mean);
		//gFluxBiasVSSNRGraph->SetPointError(pntCounter,0,fluxdiff_err);
		gFluxBiasVSSNRGraph->SetPoint(pntCounter,snr_mean,fluxdiff_median);
		gFluxBiasVSSNRGraph->SetPointError(pntCounter,0,fluxdiff_median_err);
		gFluxResoVSSNRGraph->SetPoint(pntCounter,snr_mean,fluxdiff_iqr_half);
		gFluxResoVSSNRGraph->SetPointError(pntCounter,0,fluxdiff_iqr_half_err);
	
		gPosXBiasVSSNRGraph->SetPoint(pntCounter,snr_mean,posx_median);
		gPosXBiasVSSNRGraph->SetPointError(pntCounter,0,posx_median_err);
		gPosYBiasVSSNRGraph->SetPoint(pntCounter,snr_mean,posy_median);
		gPosYBiasVSSNRGraph->SetPointError(pntCounter,0,posy_median_err);

		gPosXResoVSSNRGraph->SetPoint(pntCounter,snr_mean,posx_iqr_half);
		gPosXResoVSSNRGraph->SetPointError(pntCounter,0,posx_iqr_half_err);
	
		gPosYResoVSSNRGraph->SetPoint(pntCounter,snr_mean,posy_iqr_half);
		gPosYResoVSSNRGraph->SetPointError(pntCounter,0,posy_iqr_half_err);

		pntCounter++;

	}//end loop 


	//Compute completeness
	cout<<"INFO: Compute Completeness ..."<<endl;
	gCompleteness= new TEfficiency(*gFluxHisto_rec,*gFluxHisto); 
	gCompleteness->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gCompleteness->SetConfidenceLevel(0.68);

	gCompletenessGraph= gCompleteness->CreateGraph(); 
	gCompletenessGraph->SetNameTitle("CompletenessGraph","CompletenessGraph");

	gCompleteness_sel= new TEfficiency(*gFluxHisto_sel,*gFluxHisto); 
	gCompleteness_sel->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gCompleteness_sel->SetConfidenceLevel(0.68);

	gCompletenessGraph_sel= gCompleteness_sel->CreateGraph(); 
	gCompletenessGraph_sel->SetNameTitle("CompletenessGraph_sel","CompletenessGraph_sel");

	//Compute completeness
	cout<<"INFO: Compute Completeness 2D map ..."<<endl;
	gCompleteness2D= (TH2D*)gFluxHisto2D_rec->Clone("Completeness2D");
	gCompleteness2D->Divide(gFluxHisto2D);

	cout<<"INFO: Compute Completeness 2D map in gal coord ..."<<endl;
	gCompleteness2D_galcoord= (TH2D*)gSourceCountsHisto2D_galcoord_rec->Clone("Completeness2D_galcoord");
	gCompleteness2D_galcoord->Divide(gSourceCountsHisto2D_galcoord);

	//Compute completeness vs flux vs gal lat
	gCompletenessVSFluxVSGalLat= (TH2D*)gFluxGalLatHisto2D_rec->Clone("CompletenessVSFluxVSGalLat");
	gCompletenessVSFluxVSGalLat->Divide(gFluxGalLatHisto2D);

	TH2D* hh= 0;
	for(size_t k=0;k<gSourceCountsHistos2D_galcoord.size();k++)
	{
		TString histoName= Form("Completeness2D_galcoord_fluxBin%d",(int)(k+1));
		hh= (TH2D*)gSourceCountsHistos2D_galcoord_rec[k]->Clone(histoName);
		hh->Divide(gSourceCountsHistos2D_galcoord[k]);
		gCompletenessList2D_galcoord.push_back(hh);
	}

	//Compute reliability
	cout<<"INFO: Compute Reliability ..."<<endl;
	gReliability= new TEfficiency(*gFluxHisto_rec_reliability,*gFluxHisto_reliability); 
	gReliability->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gReliability->SetConfidenceLevel(0.68);

	gReliabilityGraph= gReliability->CreateGraph(); 
	gReliabilityGraph->SetNameTitle("ReliabilityGraph","ReliabilityGraph");

	gReliability_sel= new TEfficiency(*gFluxHisto_sel_reliability,*gFluxHisto_reliability); 
	gReliability_sel->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gReliability_sel->SetConfidenceLevel(0.68);

	gReliabilityGraph_sel= gReliability_sel->CreateGraph(); 
	gReliabilityGraph_sel->SetNameTitle("ReliabilityGraph_sel","ReliabilityGraph_sel");

	gReliabilityVSFluxVSGalLat= (TH2D*)gFluxGalLatHisto2D_rec_reliability->Clone("ReliabilityVSFluxVSGalLat");
	gReliabilityVSFluxVSGalLat->Divide(gFluxGalLatHisto2D_reliability);

	for(size_t k=0;k<gSourceCountsHistos2D_galcoord_reliability.size();k++)
	{
		TString histoName= Form("Reliability2D_galcoord_fluxBin%d",(int)(k+1));
		hh= (TH2D*)gSourceCountsHistos2D_galcoord_rec_reliability[k]->Clone(histoName);
		hh->Divide(gSourceCountsHistos2D_galcoord_reliability[k]);
		gReliabilityList2D_galcoord.push_back(hh);
	}

	//Compute false detection rate
	cout<<"INFO: Compute False Detection Rate ..."<<endl;
	gFalseDetectionRate= new TEfficiency(*gFluxHisto_rec_fdr,*gFluxHisto_fdr); 
	gFalseDetectionRate->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gFalseDetectionRate->SetConfidenceLevel(0.68);

	gFalseDetectionRateGraph= gFalseDetectionRate->CreateGraph(); 
	gFalseDetectionRateGraph->SetNameTitle("FalseDetectionRateGraph","FalseDetectionRateGraph");

	gFalseDetectionRate_sel= new TEfficiency(*gFluxHisto_sel_fdr,*gFluxHisto_fdr); 
	gFalseDetectionRate_sel->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gFalseDetectionRate_sel->SetConfidenceLevel(0.68);

	gFalseDetectionRateGraph_sel= gFalseDetectionRate_sel->CreateGraph(); 
	gFalseDetectionRateGraph_sel->SetNameTitle("FalseDetectionRateGraph_sel","FalseDetectionRateGraph_sel");

	return 0;

}//close AnalyzeData()

int ParseData_caesar(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check catalog version
	int NDataCols= 0;
	if(gCatalogVersion==1) NDataCols= gNDataCols_caesar;
	else if(gCatalogVersion==2) NDataCols= gNDataCols_caesar_v2; 

	//Check number of fields present
	if(fields.size()<NDataCols){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<NDataCols<<" expected!"<<endl;
		return -1;
	}
	
	//Set var values
	std::string blobName= std::string(fields[0]);
	long int nPix= atol(fields[1].c_str());
	std::string compId= std::string(fields[2]);
	std::string name= blobName + std::string("_fitcomp") + compId;
	double x= atof(fields[4].c_str());
	double y= atof(fields[5].c_str());
	double x_wcs= atof(fields[8].c_str());//wcs x
	double y_wcs= atof(fields[9].c_str());//wcs y
	double A= atof(fields[13].c_str());//in Jy/beam
	double beamArea= atof(fields[19].c_str());
	double fluxDensity= atof(fields[15].c_str());	
	double flux= fluxDensity/beamArea;
	if(gCatalogVersion==3){
		flux= atof(fields[15].c_str());
		fluxDensity= flux*beamArea;		
	}

	double bmaj_pix= atof(fields[20].c_str());
	double bmin_pix= atof(fields[21].c_str());
	double pa_pix= atof(fields[22].c_str());
	double bmaj_wcs= atof(fields[26].c_str());
	double bmin_wcs= atof(fields[27].c_str());
	double pa_wcs= atof(fields[28].c_str());

	double bkgSum= 0;
	double rmsSum= 0;
	double chi2= 0;
	double ndf= 0;
	double eccentricityRatio= 0;
	double areaRatio= 0;
	double rotAngleVSBeam= 0;
	int fitQuality= 0;
	int flag= 0;
	if(gCatalogVersion==1){
		bkgSum= atof(fields[38].c_str());
		rmsSum= atof(fields[39].c_str());
		chi2= atof(fields[NDataCols-2].c_str());
		ndf= atof(fields[NDataCols-1].c_str());
		eccentricityRatio= 0;//not defined in v1
		areaRatio= 0;//not defined in v1
		rotAngleVSBeam= 0;//not defined in v1
		fitQuality= 0;//not defined in v1
		flag= 0;//not defined in v1
	}
	else if(gCatalogVersion==2){
		bkgSum= atof(fields[41].c_str());
		rmsSum= atof(fields[42].c_str());
		chi2= atof(fields[NDataCols-4].c_str());
		ndf= atof(fields[NDataCols-3].c_str());
		eccentricityRatio= atof(fields[38].c_str());
		areaRatio= atof(fields[39].c_str());
		rotAngleVSBeam= atof(fields[40].c_str());
		fitQuality= atoi(fields[NDataCols-2].c_str());
		flag= atoi(fields[NDataCols-1].c_str());
	}

	double avgBkg= bkgSum/(double)(nPix);
	double avgRMS= rmsSum/(double)(nPix);
	double snr= A/avgRMS;
	
	double ra= 0;
	double dec= 0;
	AstroUtils::PixelToWCSCoords(ra,dec,gWCS,x,y); 

	double l= 0;
	double b= 0;
	AstroUtils::PixelToWCSCoords(l,b,gGalWCS,x,y); 
	//cout<<"(x,y)=("<<x<<","<<y<<"), (l,b)=("<<l<<","<<b<<")"<<endl;

	info.name= name;
	/*
	if(gUseWCS){
		info.x= x_wcs;
		info.y= y_wcs;
	}
	else{
		info.x= x;
		info.y= y;
	}
	*/
	info.x= x;
	info.y= y;

	info.ra= ra;
	info.dec= dec;
	info.l= l;
	info.b= b;
	info.flux= flux;
	info.rms= avgRMS;
	info.snr= snr;
	info.chi2= chi2;
	info.ndf= ndf;
	info.A= A;
	info.bmaj_pix= bmaj_pix;
	info.bmin_pix= bmin_pix;
	info.pa_pix= pa_pix;
	info.bmaj_wcs= bmaj_wcs;
	info.bmin_wcs= bmin_wcs;
	info.pa_wcs= pa_wcs;
	

	info.eccentricityRatio= eccentricityRatio;
	info.areaRatio= areaRatio;
	info.rotAngleVSBeam= rotAngleVSBeam;
	info.fitQuality= fitQuality;
	info.flag= flag;

	//cout<<"Source {name="<<info.name<<", x="<<info.x<<", y="<<info.y<<", flux="<<info.flux<<", beam="<<beamArea<<", fluxDensity="<<fluxDensity<<", rmsSum="<<rmsSum<<", avgRMS="<<avgRMS<<", snr="<<snr<<"}"<<endl;
	
	return 0;

}//close ParseData_caesar()


int ParseData_selavy(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols_selavy){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols_selavy<<" expected!"<<endl;
		return -1;
	}
	
	//Set var values
	std::string blobName= std::string(fields[0]);
	std::string compId= std::string(fields[1]);
	std::string name= blobName + std::string("_fitcomp") + compId;
	double x_wcs= atof(fields[5].c_str());
	double y_wcs= atof(fields[6].c_str());
	double A= atof(fields[10].c_str())/1000.;//in Jy/beam
	double flux= atof(fields[12].c_str())/1000;//in Jy
	double rms= atof(fields[32].c_str())/1000;//in Jy/beam
	double snr= A/rms;
	double chi2= atof(fields[gNDataCols_selavy-8].c_str());
	double ndf= atof(fields[gNDataCols_selavy-4].c_str());

	double bmaj_pix= 0;//not defined in catalog
	double bmin_pix= 0;//not defined in catalog
	double pa_pix= 0;//not defined in catalog
	double bmaj_wcs= atof(fields[10].c_str());
	double bmin_wcs= atof(fields[11].c_str());
	double pa_wcs= atof(fields[12].c_str());

	double x= 0;
	double y= 0;
	AstroUtils::WCSToPixelCoords(x,y,gWCS,x_wcs,y_wcs); 
	
	double ra= x_wcs;
	double dec= y_wcs;
	
	double l= 0;
	double b= 0;
	AstroUtils::PixelToWCSCoords(l,b,gGalWCS,x,y); 


	info.name= name;
	/*
	if(gUseWCS){
		info.x= x_wcs;
		info.y= y_wcs;
	}
	else{
		//cerr<<"ERROR: Selavy does not provide source centroid in pixel coordinates!"<<endl;
		//return -1;
		info.x= x;
		info.y= y;
	}
	*/
	info.x= x;
	info.y= y;

	info.ra= ra;
	info.dec= dec;
	info.l= l;
	info.b= b;

	info.flux= flux;
	info.rms= rms;
	info.snr= snr;
	info.chi2= chi2;
	info.ndf= ndf;
	info.A= A;

	info.bmaj_pix= bmaj_pix;
	info.bmin_pix= bmin_pix;
	info.pa_pix= pa_pix;
	info.bmaj_wcs= bmaj_wcs;
	info.bmin_wcs= bmin_wcs;
	info.pa_wcs= pa_wcs;

	//cout<<"Source {name="<<info.name<<", x="<<info.x<<", y="<<info.y<<", flux="<<info.flux<<", beam="<<beamArea<<", fluxDensity="<<fluxDensity<<", rmsSum="<<rmsSum<<", avgRMS="<<avgRMS<<", snr="<<snr<<"}"<<endl;
	
	return 0;

}//close ParseData_selavy()




int ParseData_aegean(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols_aegean){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols_aegean<<" expected!"<<endl;
		return -1;
	}
	
	//#   island	source	background	local_rms	ra_str	dec_str	ra	err_ra	dec	err_dec	peak_flux	err_peak_flux	int_flux	err_int_flux	a	err_a	b	err_b	pa	err_pa	flags	residual_mean	residual_std	uuid	psf_a	psf_b	psf_pa

	//Set var values
	std::string blobName= std::string(fields[0]);
	std::string compId= std::string(fields[1]);
	std::string name= std::string("S") + blobName + std::string("_fitcomp") + compId;
	double x_wcs= atof(fields[6].c_str());
	double y_wcs= atof(fields[8].c_str());
	double A= atof(fields[10].c_str());//in Jy/beam
	double flux= atof(fields[12].c_str());//in Jy
	double rms= atof(fields[3].c_str());//in Jy/beam
	double snr= A/rms;
	double chi2= 0;//not present in catalogue
	double ndf= 0;//not present in catalogue

	double bmaj_pix= 0;//not defined in catalog
	double bmin_pix= 0;//not defined in catalog
	double pa_pix= 0;//not defined in catalog
	double bmaj_wcs= atof(fields[14].c_str());
	double bmin_wcs= atof(fields[16].c_str());
	double pa_wcs= atof(fields[18].c_str());

	double x= 0;
	double y= 0;
	AstroUtils::WCSToPixelCoords(x,y,gWCS,x_wcs,y_wcs); 
	
	double ra= x_wcs;
	double dec= y_wcs;
	
	double l= 0;
	double b= 0;
	AstroUtils::PixelToWCSCoords(l,b,gGalWCS,x,y); 


	info.name= name;

	/*
	if(gUseWCS){
		info.x= x_wcs;
		info.y= y_wcs;
	}
	else{
		//cerr<<"ERROR: Aegean does not provide source centroid in pixel coordinates!"<<endl;
		//return -1;
		info.x= x;
		info.y= y;
	}
	*/
	info.x= x;
	info.y= y;

	info.ra= ra;
	info.dec= dec;
	info.l= l;
	info.b= b;

	info.flux= flux;
	info.rms= rms;
	info.snr= snr;
	info.chi2= chi2;
	info.ndf= ndf;
	info.A= A;

	info.bmaj_pix= bmaj_pix;
	info.bmin_pix= bmin_pix;
	info.pa_pix= pa_pix;
	info.bmaj_wcs= bmaj_wcs;
	info.bmin_wcs= bmin_wcs;
	info.pa_wcs= pa_wcs;

	//cout<<"Source {name="<<info.name<<", x="<<info.x<<", y="<<info.y<<", flux="<<info.flux<<", beam="<<beamArea<<", fluxDensity="<<fluxDensity<<", rmsSum="<<rmsSum<<", avgRMS="<<avgRMS<<", snr="<<snr<<"}"<<endl;
	
	return 0;

}//close ParseData_aegean()

int ParseData(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols){
		cerr<<"ERROR: Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols<<" expected!"<<endl;
		return -1;
	}
	
	//Parse data
	std::string name= std::string(fields[0]);
	double x= atof(fields[1].c_str());
	double y= atof(fields[2].c_str());
	double x_wcs= atof(fields[3].c_str());//wcs x
	double y_wcs= atof(fields[4].c_str());//wcs y
	double flux= atof(fields[5].c_str());
	
	double ra= 0;
	double dec= 0;
	AstroUtils::PixelToWCSCoords(ra,dec,gWCS,x,y); 

	double l= 0;
	double b= 0;
	AstroUtils::PixelToWCSCoords(l,b,gGalWCS,x,y); 

	//Set var values
	info.name= name;
	 
	/*
	if(gUseWCS){
		info.x= x_wcs;
		info.y= y_wcs;
	}
	else{
		info.x= x;
		info.y= y;
	}
	*/
	info.x= x;
	info.y= y;

	info.flux= flux;

	info.ra= ra;
	info.dec= dec;
	info.l= l;
	info.b= b;

	
	return 0;

}//close ParseData()


int ReadData(std::vector<SourceInfo>& sources,std::string filename,int (*ParseFcn)(SourceInfo&, const std::vector<std::string>&),std::vector<std::string> excludedPatterns)
{
	//Clear collection
	sources.clear();

	//Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		cerr<<"ERROR: Failed to open config file "<<filename<<" for reading!"<<endl;
		return -1;
  }

	//Parsing file
	cout<<"INFO: Reading and parsing file: "<<filename<<endl;

	std::string parsedline= "";	
	int line_counter= 0;

	while(std::getline(in,parsedline)) {
		line_counter++;

		//Check file
		if (!in.good()) break;

		//Parse line fields
		istringstream ss(parsedline);
		std::string field= "";
		bool skipLine= false;
		std::vector<std::string> fields;

		while(ss >> field)
    {
			//Skip first line
			char first_char= *(field.c_str());
			if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char=='\t' || first_char==' '){
				skipLine= true;
				break;
			}

			//Search for pattern and skip if present
			for(size_t i=0;i<excludedPatterns.size();i++){
				bool hasPattern= HasPattern(field,excludedPatterns[i]);
				if(hasPattern) {
					skipLine= true;
					break;
				}
			}
			if(skipLine) break;
	
			//cout<<field<<"  ";
			fields.push_back(field);

		}//end loop line
		
		if(skipLine) continue;
		//else cout<<endl;

		
		//Process line
		SourceInfo info;
 		int status= (*ParseFcn)(info,fields);
		if(status==0){
			sources.push_back(info);
		}
		else{
			cerr<<"ERROR: Failed to parse line "<<line_counter<<" of file "<<filename<<"!"<<endl;
		}	

		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	return 0;

}//close ReadData()


bool HasPattern(std::string str,std::string pattern)
{
	std::size_t found = str.find(pattern);
  if (found!=std::string::npos) return true;

	return false;

}//close HasPattern()
