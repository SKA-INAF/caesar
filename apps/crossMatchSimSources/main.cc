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
#include <AstroUtils.h>
#include <StatsUtils.h>
#include <DS9Region.h>
#include <DS9RegionParser.h>


#include <ConfigParser.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>
#include <MathUtils.h>
#include <EllipseUtils.h>
#include <Contour.h>
#include <WCSUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <numeric>

using namespace std;
using namespace Caesar;

void Usage(char* exeName)
{
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --inputfile=[INPUT_FILE] \t Input file or filelist containing sim files"<<endl;
	cout<<"-I, --inputfile_rec=[INPUT_FILE_REC] \t Input file or filelist containing files to be crosmatched"<<endl;
	cout<<"-m, --mapfile=[INPUT_FILE] \t Input map file in FITS format"<<endl;
	cout<<"-r, --regionfile=[INPUT_FILE] \t Input DS9 region file used to select sources (default=empty)"<<endl;
	cout<<"-c, --catalogType=[CATALOG_TYPE] \t Catalog type id {1=CAESAR,SELAVY,AEGEAN} (default=CAESAR)"<<endl;
	cout<<"-C, --catalogVersion=[CATALOG_VERSION] \t Catalog version (default=3)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (root format) in which to store source match info (default=output.root)"<<endl;
	cout<<"-T, --matchPosThr=[POS_THRESHOLD] \t Source island centroid distance in arcsec below which we have a match (default=8)"<<endl;
	
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "inputfile", required_argument, 0, 'i' },
	{ "inputfile_rec", required_argument, 0, 'I' },
	{ "mapfile", required_argument, 0, 'm' },
	{ "regionfile", required_argument, 0, 'r' },
	{ "catalogVersion", required_argument, 0, 'C' },
	{ "catalogType", required_argument, 0, 'c' },	
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", required_argument, 0, 'o' },	
	{ "matchPosThr", required_argument, 0, 't'},	
  {(char*)0, (int)0, (int*)0, (int)0}
};



//Options
std::string fileName= "";
std::string fileName_rec= "";
std::string gOutputFileName= "output.root";
bool gSelectRegion= false;
std::string gRegionFileName= "";
std::string gMapFileName= "";
float gMatchPosThr= 8;//dist in arcsec below which two sources can be matched
int verbosity= 4;//INFO level

enum CatalogType {
	eCAESAR_CATALOG=1,
	eSELAVY_CATALOG=2,
	eAEGEAN_CATALOG=3
};
int gCatalogType= 1;
int gCatalogVersion= 3;


//==========================================
//==     GLOBAL VARS
//==========================================
//- INPUT DATA
std::vector<std::string> gFileNames {};
std::vector<std::string> gFileNames_rec {};
std::vector<DS9Region*> gRegions;
WCS* gGalWCS= 0;
WCS* gWCS= 0;

// - CATALOG READ
const int gNDataCols= 6;
const int gNDataCols_caesar= 42;
const int gNDataCols_caesar_v2= 47;
const int gNDataCols_caesar_v3= 47;
const int gNDataCols_selavy= 37;
const int gNDataCols_aegean= 27;
std::vector<std::string> gSkipLinePatterns {};
std::vector<std::string> gSkipLinePatterns_caesar {};
std::vector<std::string> gSkipLinePatterns_selavy {};
std::vector<std::string> gSkipLinePatterns_aegean {"island"};

// - GRAPH/HISTO DATA
std::vector<double> gLgFluxBins {-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
std::vector<double> gLgFluxBins_min {-5.0,-4.5,-4.0,-3.5,-3.00,-2.75,-2.50,-2.25,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0};
std::vector<double> gLgFluxBins_max {-4.5,-4.0,-3.5,-3.0,-2.75,-2.50,-2.25,-2.00,-1.5,-1.0,-0.5, 0.0,0.5,1.0,5.0};

//std::vector<double> gLgRecFluxBins {-5,-4.5,-4,-3.5,-3,-2.75,-2.5,-2.25,-2,-1.5,-1,-0.5,0,0.5,1};
std::vector<double> gSNRBins {0,2,4,5,7,10,15,20,25,50,100,1000,10000};
std::vector<double> gSNRBins_min {0,2,4,5,7, 10,15,20,25,50, 100,  1000};
std::vector<double> gSNRBins_max {2,4,5,7,10,15,20,25,50,100,1000,10000};

std::vector<double> gGalLatBins {-5,-4,-3,-2,-1,0,1,2,3,4,5};
std::vector<double> gGalLatBins_min {-5,-4,-3,-2,-1,0,1,2,3,4};
std::vector<double> gGalLatBins_max {-4,-3,-2,-1, 0,1,2,3,4,5};

TH1D* gFluxHisto= 0;
TH1D* gFluxHisto_rec= 0;
TH1D* gFluxHisto_reliability= 0;
TH1D* gFluxHisto_fdr= 0;
TH1D* gFluxHisto_rec_reliability= 0;
TH1D* gFluxHisto_rec_fdr= 0;
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
TGraphErrors* gFluxBiasVSRecFluxGraph= 0;
TGraphErrors* gFluxResoVSRecFluxGraph= 0; 
TGraphErrors* gPosXBiasGraph= 0; 
TGraphErrors* gPosYBiasGraph= 0; 
TGraphErrors* gPosXResoGraph= 0; 
TGraphErrors* gPosYResoGraph= 0;

TGraphErrors* gFluxBiasVSSNRGraph= 0;
TGraphErrors* gFluxResoVSSNRGraph= 0; 
TGraphErrors* gPeakFluxBiasVSSNRGraph= 0;
TGraphErrors* gPeakFluxResoVSSNRGraph= 0; 
TGraphErrors* gPosXBiasVSSNRGraph= 0; 
TGraphErrors* gPosYBiasVSSNRGraph= 0; 
TGraphErrors* gPosXResoVSSNRGraph= 0; 
TGraphErrors* gPosYResoVSSNRGraph= 0; 

//- Rec vectors & graphs
std::vector<std::vector<double>> gFluxBinData;
std::vector<std::vector<double>> gFluxDiffBinData;
std::vector<std::vector<double>> gLgFluxDiffBinData;
std::vector<std::vector<double>> gPosXDiffBinData;
std::vector<std::vector<double>> gPosYDiffBinData;

std::vector<std::vector<double>> gSNRBinData;
std::vector<std::vector<double>> gFluxDiffVSSNRBinData;
std::vector<std::vector<double>> gPeakFluxDiffVSSNRBinData;
std::vector<std::vector<double>> gPosXDiffVSSNRBinData;
std::vector<std::vector<double>> gPosYDiffVSSNRBinData;

std::vector<std::vector<double>> gRecFluxBinData;
std::vector<std::vector<double>> gFluxDiffVSRecFluxBinData;
std::vector<std::vector<double>> gLgFluxDiffVSRecFluxBinData;

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
double gSourceSNR_true;
std::string gSourceName_rec;
double gSourcePosX_rec;
double gSourcePosY_rec;
double gSourceFlux_rec;
double gSourceGalPosX_rec;
double gSourceGalPosY_rec;
double gSourceJ2000PosX_rec;
double gSourceJ2000PosY_rec;
double gSourceSNR_rec;
bool gIsFound;

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




//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int Init();
int ReadFileLists();
int ReadMap();
int ReadRegions();
int ReadCatalogData();
int AnalyzeData();
int ReadData(std::vector<SourceInfo>& sources,std::string filename,int (*fcn)(SourceInfo&, const std::vector<std::string>&), std::vector<std::string> excludedPatterns= {});
static int ParseData(SourceInfo& info,const std::vector<std::string>& fields);
static int ParseData_caesar(SourceInfo& info,const std::vector<std::string>& fields);
static int ParseData_selavy(SourceInfo& info,const std::vector<std::string>& fields);
static int ParseData_aegean(SourceInfo& info,const std::vector<std::string>& fields);
bool HasPattern(std::string,std::string);
int FindSourceMatch(const std::vector<SourceInfo>& sources_ref,const std::vector<SourceInfo>& sources_rec);
void Save();


int main(int argc, char *argv[])
{
	//================================
	//== PARSE CMD LINE OPTIONS
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

	//================================
	//==      READ FILE LISTS
	//================================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading input file lists...");
	#endif
	if(ReadFileLists()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read input file lists!");
		#endif
		return -1;
	}

	
	//=======================
	//== READ REGIONS
	//=======================
	if(gSelectRegion){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading regions in file "<<gRegionFileName<<" ...");
		#endif
		if(ReadRegions()<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Region reading failed!");
			#endif
			return -1;
		}
	}//close if


	//=======================
	//== READ MAP
	//=======================
	if(ReadMap()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading input map failed!");
		#endif
		return -1;
	}
	
	//================================
	//==      READ CATALOG DATA
	//================================
	if(ReadCatalogData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading catalog data failed!");
		#endif
		return -1;
	}
	
	//================================
	//==      ANALYZE DATA
	//================================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Analyzing data...");
	#endif
	if(AnalyzeData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Analyze data failed!");
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
		INFO_LOG("End macro");
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
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hi:I:C:c:o:v:m:r:t:",options_tab, &option_index)) != -1) {
    
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
			case 'I':	
			{
				fileName_rec= std::string(optarg);	
				break;	
			}
			case 'C':	
			{
				gCatalogVersion= atoi(optarg);
				break;	
			}
			case 'c':	
			{
				gCatalogType= atoi(optarg);	
				break;	
			}
			case 'o':	
			{
				gOutputFileName= std::string(optarg);	
				break;	
			}
			case 'r':	
			{
				gRegionFileName= std::string(optarg);	
				if(gRegionFileName!="") gSelectRegion= true;
				break;	
			}
			case 'm':	
			{
				gMapFileName= std::string(optarg);	
				break;	
			}
			case 't':
			{
				gMatchPosThr= atof(optarg);
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
	//== INIT LOGGER
	//=======================
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	//=======================
	//== CHECK OPTIONS
	//=======================
	if( gCatalogType!=eCAESAR_CATALOG &&
			gCatalogType!=eSELAVY_CATALOG &&
			gCatalogType!=eAEGEAN_CATALOG
	){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unsupported catalog type given (type="<<gCatalogType<<")!");
		#endif
		return -1;
	}	

	if(gMapFileName=="")
	{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty map filename given!");
		#endif
		return -1;
	}

	if(gOutputFileName=="")
	{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output filename given!");
		#endif
		return -1;
	}

	//=======================
	//== PRINT OPTIONS
	//=======================
	#ifdef LOGGING_ENABLED	
		INFO_LOG("catalogType="<<gCatalogType<<", catalogVersion="<<gCatalogVersion);
		INFO_LOG("matchPosThr(arcsec)="<<gMatchPosThr);
	#endif

	return 0;

}//close ParseOptions()

int Init()
{
	//Create output file
	gOutputFile= new TFile(gOutputFileName.c_str(),"RECREATE");
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
	gMatchDataTree->Branch("snr_rec",&gSourceSNR_rec);

	//Init histos
	int nBins= (int)(gLgFluxBins.size()-1);
	int nBins_snr= (int)(gSNRBins.size()-1);
	//int nBins_rec= (int)(gLgRecFluxBins.size()-1);
	int nBins_rec= (int)(gLgFluxBins.size()-1);
	int nBins_galLat= (int)(gGalLatBins.size()-1);
	
	gFluxHisto= new TH1D("fluxHisto","",nBins,gLgFluxBins.data());
	gFluxHisto->Sumw2();

	gFluxHisto_reliability= new TH1D("fluxHisto_reliability","",nBins,gLgFluxBins.data());
	gFluxHisto_reliability->Sumw2();

	gFluxHisto_fdr= new TH1D("fluxHisto_fdr","",nBins,gLgFluxBins.data());
	gFluxHisto_fdr->Sumw2();

	gFluxHisto_rec= new TH1D("fluxHisto_rec","",nBins,gLgFluxBins.data());
	gFluxHisto_rec->Sumw2();

	gFluxHisto_rec_reliability= new TH1D("fluxHisto_rec_reliability","",nBins,gLgFluxBins.data());
	gFluxHisto_rec_reliability->Sumw2();

	gFluxHisto_rec_fdr= new TH1D("fluxHisto_rec_fdr","",nBins,gLgFluxBins.data());
	gFluxHisto_rec_fdr->Sumw2();

	//Init graphs
	gFluxBiasGraph= new TGraphErrors;
	gFluxBiasGraph->SetName("fluxBiasGraph");
	gFluxBiasGraph->SetMarkerStyle(8);
	gFluxBiasGraph->SetMarkerSize(1.3);
	gFluxBiasGraph->SetMarkerColor(kBlack);
	gFluxBiasGraph->SetLineColor(kBlack);

	gFluxResoGraph= new TGraphErrors;
	gFluxResoGraph->SetName("fluxResoGraph");
	gFluxResoGraph->SetMarkerStyle(8);
	gFluxResoGraph->SetMarkerSize(1.3);
	gFluxResoGraph->SetMarkerColor(kBlack);
	gFluxResoGraph->SetLineColor(kBlack);

	gLgFluxBiasGraph= new TGraphErrors;
	gLgFluxBiasGraph->SetName("lgfluxBiasGraph");
	gLgFluxBiasGraph->SetMarkerStyle(8);
	gLgFluxBiasGraph->SetMarkerSize(1.3);
	gLgFluxBiasGraph->SetMarkerColor(kBlack);
	gLgFluxBiasGraph->SetLineColor(kBlack);

	gLgFluxResoGraph= new TGraphErrors;
	gLgFluxResoGraph->SetName("lgfluxResoGraph");
	gLgFluxResoGraph->SetMarkerStyle(8);
	gLgFluxResoGraph->SetMarkerSize(1.3);
	gLgFluxResoGraph->SetMarkerColor(kBlack);
	gLgFluxResoGraph->SetLineColor(kBlack);

	gFluxBiasVSRecFluxGraph= new TGraphErrors;
	gFluxBiasVSRecFluxGraph->SetName("fluxBiasVSRecFluxGraph");
	gFluxBiasVSRecFluxGraph->SetMarkerStyle(8);
	gFluxBiasVSRecFluxGraph->SetMarkerSize(1.3);
	gFluxBiasVSRecFluxGraph->SetMarkerColor(kBlack);
	gFluxBiasVSRecFluxGraph->SetLineColor(kBlack);

	gFluxResoVSRecFluxGraph= new TGraphErrors;
	gFluxResoVSRecFluxGraph->SetName("fluxResoVSRecFluxGraph");
	gFluxResoVSRecFluxGraph->SetMarkerStyle(8);
	gFluxResoVSRecFluxGraph->SetMarkerSize(1.3);
	gFluxResoVSRecFluxGraph->SetMarkerColor(kBlack);
	gFluxResoVSRecFluxGraph->SetLineColor(kBlack);

	gLgFluxBiasVSRecFluxGraph= new TGraphErrors;
	gLgFluxBiasVSRecFluxGraph->SetName("lgfluxBiasVSRecFluxGraph");
	gLgFluxBiasVSRecFluxGraph->SetMarkerStyle(8);
	gLgFluxBiasVSRecFluxGraph->SetMarkerSize(1.3);
	gLgFluxBiasVSRecFluxGraph->SetMarkerColor(kBlack);
	gLgFluxBiasVSRecFluxGraph->SetLineColor(kBlack);

	gLgFluxResoVSRecFluxGraph= new TGraphErrors;
	gLgFluxResoVSRecFluxGraph->SetName("lgfluxResoVSRecFluxGraph");
	gLgFluxResoVSRecFluxGraph->SetMarkerStyle(8);
	gLgFluxResoVSRecFluxGraph->SetMarkerSize(1.3);
	gLgFluxResoVSRecFluxGraph->SetMarkerColor(kBlack);
	gLgFluxResoVSRecFluxGraph->SetLineColor(kBlack);

	gPosXBiasGraph= new TGraphErrors;
	gPosXBiasGraph->SetName("posXBiasGraph");
	gPosXBiasGraph->SetMarkerStyle(8);
	gPosXBiasGraph->SetMarkerSize(1.3);
	gPosXBiasGraph->SetMarkerColor(kBlack);
	gPosXBiasGraph->SetLineColor(kBlack);

	gPosYBiasGraph= new TGraphErrors;
	gPosYBiasGraph->SetName("posYBiasGraph");
	gPosYBiasGraph->SetMarkerStyle(8);
	gPosYBiasGraph->SetMarkerSize(1.3);
	gPosYBiasGraph->SetMarkerColor(kBlack);
	gPosYBiasGraph->SetLineColor(kBlack);

	gPosXResoGraph= new TGraphErrors;
	gPosXResoGraph->SetName("posXResoGraph");
	gPosXResoGraph->SetMarkerStyle(8);
	gPosXResoGraph->SetMarkerSize(1.3);
	gPosXResoGraph->SetMarkerColor(kBlack);
	gPosXResoGraph->SetLineColor(kBlack);

	gPosYResoGraph= new TGraphErrors;
	gPosYResoGraph->SetName("posYResoGraph");
	gPosYResoGraph->SetMarkerStyle(8);
	gPosYResoGraph->SetMarkerSize(1.3);
	gPosYResoGraph->SetMarkerColor(kBlack);
	gPosYResoGraph->SetLineColor(kBlack);

	gFluxBiasVSSNRGraph= new TGraphErrors;
	gFluxBiasVSSNRGraph->SetName("fluxBiasVSSNRGraph");
	gFluxBiasVSSNRGraph->SetMarkerStyle(8);
	gFluxBiasVSSNRGraph->SetMarkerSize(1.3);
	gFluxBiasVSSNRGraph->SetMarkerColor(kBlack);
	gFluxBiasVSSNRGraph->SetLineColor(kBlack);

	gFluxResoVSSNRGraph= new TGraphErrors;
	gFluxResoVSSNRGraph->SetName("fluxResoVSSNRGraph");
	gFluxResoVSSNRGraph->SetMarkerStyle(8);
	gFluxResoVSSNRGraph->SetMarkerSize(1.3);
	gFluxResoVSSNRGraph->SetMarkerColor(kBlack);
	gFluxResoVSSNRGraph->SetLineColor(kBlack);


	gPeakFluxBiasVSSNRGraph= new TGraphErrors;
	gPeakFluxBiasVSSNRGraph->SetName("peakFluxBiasVSSNRGraph");
	gPeakFluxBiasVSSNRGraph->SetMarkerStyle(8);
	gPeakFluxBiasVSSNRGraph->SetMarkerSize(1.3);
	gPeakFluxBiasVSSNRGraph->SetMarkerColor(kBlack);
	gPeakFluxBiasVSSNRGraph->SetLineColor(kBlack);

	gPeakFluxResoVSSNRGraph= new TGraphErrors;
	gPeakFluxResoVSSNRGraph->SetName("peakFluxResoVSSNRGraph");
	gPeakFluxResoVSSNRGraph->SetMarkerStyle(8);
	gPeakFluxResoVSSNRGraph->SetMarkerSize(1.3);
	gPeakFluxResoVSSNRGraph->SetMarkerColor(kBlack);
	gPeakFluxResoVSSNRGraph->SetLineColor(kBlack);

	gPosXBiasVSSNRGraph= new TGraphErrors;
	gPosXBiasVSSNRGraph->SetName("posXBiasVSSNRGraph");
	gPosXBiasVSSNRGraph->SetMarkerStyle(8);
	gPosXBiasVSSNRGraph->SetMarkerSize(1.3);
	gPosXBiasVSSNRGraph->SetMarkerColor(kBlack);
	gPosXBiasVSSNRGraph->SetLineColor(kBlack);

	gPosYBiasVSSNRGraph= new TGraphErrors;
	gPosYBiasVSSNRGraph->SetName("posYBiasVSSNRGraph");
	gPosYBiasVSSNRGraph->SetMarkerStyle(8);
	gPosYBiasVSSNRGraph->SetMarkerSize(1.3);
	gPosYBiasVSSNRGraph->SetMarkerColor(kBlack);
	gPosYBiasVSSNRGraph->SetLineColor(kBlack);

	gPosXResoVSSNRGraph= new TGraphErrors;
	gPosXResoVSSNRGraph->SetName("posXResoVSSNRGraph");
	gPosXResoVSSNRGraph->SetMarkerStyle(8);
	gPosXResoVSSNRGraph->SetMarkerSize(1.3);
	gPosXResoVSSNRGraph->SetMarkerColor(kBlack);
	gPosXResoVSSNRGraph->SetLineColor(kBlack);

	gPosYResoVSSNRGraph= new TGraphErrors;
	gPosYResoVSSNRGraph->SetName("posYResoVSSNRGraph");
	gPosYResoVSSNRGraph->SetMarkerStyle(8);
	gPosYResoVSSNRGraph->SetMarkerSize(1.3);
	gPosYResoVSSNRGraph->SetMarkerColor(kBlack);
	gPosYResoVSSNRGraph->SetLineColor(kBlack);

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
		gPeakFluxDiffVSSNRBinData.push_back( std::vector<double>() );
		gPosXDiffVSSNRBinData.push_back( std::vector<double>() );
		gPosYDiffVSSNRBinData.push_back( std::vector<double>() );
	}
	for(int i=0;i<nBins_rec;i++){
		gRecFluxBinData.push_back( std::vector<double>() );
		gLgFluxDiffVSRecFluxBinData.push_back( std::vector<double>() );
		gFluxDiffVSRecFluxBinData.push_back( std::vector<double>() );
	}


	return 0;

}//close Init()


int FindSourceMatch(const std::vector<SourceInfo>& sources,const std::vector<SourceInfo>& sources_rec)
{
	
	//- Find if ref source if found in rec catalog
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
		double peakFlux= sources[i].flux;
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
				if(gRegions[k]->IsPointInsideRegion(ra,dec)){
					isInside= true;
					break;
				}
			}
			if(!isInside) continue;
		}//close if gSelectRegion

		
		gFluxHisto->Fill(lgFlux,1);
		nSources++;
	
		bool found= false;
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
		std::vector<double> dpeakflux_list;
		std::vector<double> snr_list;
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
			double peakFlux_rec= sources_rec[j].A;
			double rms_rec= sources_rec[j].rms;
			double snr_rec= sources_rec[j].snr;
			
			//Check if sources match given position
			double dx= (x_rec-x);//in pixels
			double dy= (y_rec-y);//in pixels
			double dra= (ra_rec-ra)*3600.;//in arcsec
			double ddec= (dec_rec-dec)*3600.;//in arcsec
			double dpos= sqrt(dx*dx + dy*dy);
			double dpos_wcs= sqrt(dra*dra + ddec*ddec);
			double dflux= flux_rec-flux;
			double dflux_rel= flux_rec/flux-1;
			double dlgflux= log10(flux_rec) - log10(flux);
			double dpeakflux_rel= peakFlux_rec/peakFlux-1;

			//Skip source if outside mosaic region
			if(gSelectRegion){
				bool isInside= false;
				for(size_t k=0;k<gRegions.size();k++){
					if(gRegions[k]->IsPointInsideRegion(ra_rec,dec_rec)){
						isInside= true;
						break;
					}
				}
				if(!isInside) continue;
			}//close if gSelectRegion

			if(dpos_wcs>gMatchPosThr) continue;
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
			dpeakflux_list.push_back(dpeakflux_rel);	
			snr_list.push_back(snr_rec);

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
		gSourceSNR_rec= -999;

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
			double dpeakflux_rel= dpeakflux_list[0];
			double snr_rec= snr_list[0];
			size_t recIndex= recindex_list[0];

			if(nFound>1){
				int best_index= 0;
				for(int k=1;k<nFound;k++){
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
				dpeakflux_rel= dpeakflux_list[best_index];
				recIndex= recindex_list[best_index];
	
				//cout<<"INFO: #"<<nFound<<" match found, match no. "<<best_index<<" selected as the best..."<<endl;

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
			gSourceSNR_rec= snr_rec;

			isRecSourceFound[recIndex]= true;

			//Fill histos
			//gPosXDiffHisto->Fill(lgFlux,dra);
			//gPosYDiffHisto->Fill(lgFlux,ddec);
			//gFluxDiffHisto->Fill(lgFlux,dflux_rel);
			//gFluxDiffVSRecFluxHisto->Fill(lgFlux_rec,dflux_rel);
			gFluxHisto_rec->Fill(lgFlux,1);
			//gSNRHisto_rec->Fill(snr_rec);
			
			//Find true flux bin
			int binId= -1;
			for(size_t k=0;k<gLgFluxBins_min.size();k++){
				if(lgFlux>=gLgFluxBins_min[k] && lgFlux<=gLgFluxBins_max[k]){
					binId= k;
					break;
				}
			}
			//Find rec flux bin
			int binId_rec= -1;
			for(size_t k=0;k<gLgFluxBins_min.size();k++){
				if(lgFlux_rec>=gLgFluxBins_min[k] && lgFlux_rec<=gLgFluxBins_max[k]){
					binId_rec= k;
					break;
				}
			}

			//Find rec SNR bin
			int binId_snr= -1;
			for(size_t k=0;k<gSNRBins_min.size();k++){
				if(snr_rec>=gSNRBins_min[k] && snr_rec<=gSNRBins_max[k]){
					binId_snr= k;
					break;
				}
			}
			//cout<<"snr_rec="<<snr_rec<<", binId_snr="<<binId_snr<<endl;
			
			if(binId>=0 && binId<gFluxBinData.size()) {
				gFluxBinData[binId].push_back(lgFlux);
				gFluxDiffBinData[binId].push_back(dflux_rel);
				gLgFluxDiffBinData[binId].push_back(dlgflux);
				gPosXDiffBinData[binId].push_back(dra);
				gPosYDiffBinData[binId].push_back(ddec);
			}

			if(binId_rec>=0 && binId_rec<gRecFluxBinData.size()) {
				gRecFluxBinData[binId_rec].push_back(lgFlux_rec);
				gLgFluxDiffVSRecFluxBinData[binId_rec].push_back(dlgflux);
				gFluxDiffVSRecFluxBinData[binId_rec].push_back(dflux_rel);
			}	

			if(binId_snr>=0 && binId_snr<gSNRBinData.size()) {
				gSNRBinData[binId_snr].push_back(snr_rec);
				gFluxDiffVSSNRBinData[binId_snr].push_back(dflux_rel);
				gPeakFluxDiffVSSNRBinData[binId_snr].push_back(dpeakflux_rel);
				gPosXDiffVSSNRBinData[binId_snr].push_back(dra);
				gPosYDiffVSSNRBinData[binId_snr].push_back(ddec);
			}

		}//close if found
		

		gMatchDataTree->Fill();

	}//end loop ref sources

	#ifdef LOGGING_ENABLED	
		INFO_LOG("#"<<nMatch<<"/"<<nSources<<" source match ...");
	#endif

	//======================================
	//==      RELIABILITY
	//======================================	
	int nRecSources_true= 0;
	int nSelSources_true= 0;

	for(size_t j=0;j<sources_rec.size();j++)
	{	
		if(j%1000==0) cout<<"--> "<<j+1<<"/"<<sources_rec.size()<<" sources in rec catalog processed..."<<endl;
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
				if(gRegions[k]->IsPointInsideRegion(ra_rec,dec_rec)){
					isInside= true;
					break;
				}
			}
			if(!isInside) continue;
		}//close if gSelectRegion

		//Fill histos
		gFluxHisto_reliability->Fill(lgFlux_rec,1);
		gFluxHisto_fdr->Fill(lgFlux_rec,1);
		if(isRecSourceFound[j]) {
			gFluxHisto_rec_reliability->Fill(lgFlux_rec,1);
			nRecSources_true++;
		}
		else {
			gFluxHisto_rec_fdr->Fill(lgFlux_rec,1);
		}

	}//end loop rec sources

	#ifdef LOGGING_ENABLED	
		INFO_LOG("#"<<nRecSources_true<<"/"<<sources_rec.size()<<" good detection ...");
	#endif

	return 0;

}//close FindSourceMatch()



int AnalyzeData()
{
	//Init histo & fcn 
	int pntCounter= 0;

	//Loop over data and find stats
	for(size_t i=0;i<gFluxBinData.size();i++)
	{
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

		//Compute fluxdiff stats	
		double fluxdiff_mean= 0;
		double fluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(fluxdiff_mean,fluxdiff_rms,gFluxDiffVSRecFluxBinData[i]);
		Caesar::BoxStats<double> stats_fluxdiff= StatsUtils::ComputeBoxStats(gFluxDiffVSRecFluxBinData[i]);
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

		
		//Fill graph
		gLgFluxBiasVSRecFluxGraph->SetPoint(pntCounter,flux_mean,lgfluxdiff_median);
		gLgFluxBiasVSRecFluxGraph->SetPointError(pntCounter,0,lgfluxdiff_median_err);
		gLgFluxResoVSRecFluxGraph->SetPoint(pntCounter,flux_mean,lgfluxdiff_iqr_half);
		gLgFluxResoVSRecFluxGraph->SetPointError(pntCounter,0,lgfluxdiff_iqr_half_err);

		gFluxBiasVSRecFluxGraph->SetPoint(pntCounter,flux_mean,fluxdiff_median);
		gFluxBiasVSRecFluxGraph->SetPointError(pntCounter,0,fluxdiff_median_err);
		gFluxResoVSRecFluxGraph->SetPoint(pntCounter,flux_mean,fluxdiff_iqr_half);
		gFluxResoVSRecFluxGraph->SetPointError(pntCounter,0,fluxdiff_iqr_half_err);
		
		pntCounter++;

	}//end loop rec bin data

	
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

		//Compute peak fluxdiff stats	
		double peakfluxdiff_mean= 0;
		double peakfluxdiff_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(peakfluxdiff_mean,peakfluxdiff_rms,gPeakFluxDiffVSSNRBinData[i]);
		Caesar::BoxStats<double> stats_peakfluxdiff= StatsUtils::ComputeBoxStats(gPeakFluxDiffVSSNRBinData[i]);
		double peakfluxdiff_median= stats_peakfluxdiff.median;
		double peakfluxdiff_err= peakfluxdiff_rms/sqrt(nEntries);
		double peakfluxdiff_median_err= 1.253*peakfluxdiff_rms/sqrt(nEntries);
		double peakfluxdiff_iqr= stats_peakfluxdiff.Q3-stats_peakfluxdiff.Q1;
		double peakfluxdiff_iqr_err= 1.573*peakfluxdiff_rms/sqrt(nEntries);
		double peakfluxdiff_iqr_half= peakfluxdiff_iqr/2.;
		double peakfluxdiff_iqr_half_err= 0.5*peakfluxdiff_iqr_err;
		double peakfluxdiff_err_low= peakfluxdiff_median-stats_peakfluxdiff.Q1;
		double peakfluxdiff_err_up= stats_peakfluxdiff.Q3-peakfluxdiff_median;
		double peakfluxdiff_min= stats_peakfluxdiff.minVal;
		double peakfluxdiff_max= stats_peakfluxdiff.maxVal;
		double dpeakfluxdiff= fabs(peakfluxdiff_max-peakfluxdiff_min);

		//Compute SNR stats
		double snr_mean= 0;
		double snr_rms= 0;
		Caesar::StatsUtils::ComputeMeanAndRMS(snr_mean,snr_rms,gSNRBinData[i]);

		cout<<"INFO: bin "<<i+1<<" (snr="<<snr_mean<<") flux diff stats {mean="<<fluxdiff_mean<<", rms="<<fluxdiff_rms<<", min/max="<<fluxdiff_min<<"/"<<fluxdiff_max<<", median="<<fluxdiff_median<<", IQR="<<fluxdiff_iqr<<"}"<<endl;

		//Fill graph
		gFluxBiasVSSNRGraph->SetPoint(pntCounter,snr_mean,fluxdiff_median);
		gFluxBiasVSSNRGraph->SetPointError(pntCounter,0,fluxdiff_median_err);
		gFluxResoVSSNRGraph->SetPoint(pntCounter,snr_mean,fluxdiff_iqr_half);
		gFluxResoVSSNRGraph->SetPointError(pntCounter,0,fluxdiff_iqr_half_err);
	
		gPeakFluxBiasVSSNRGraph->SetPoint(pntCounter,snr_mean,peakfluxdiff_median);
		gPeakFluxBiasVSSNRGraph->SetPointError(pntCounter,0,peakfluxdiff_median_err);
		gPeakFluxResoVSSNRGraph->SetPoint(pntCounter,snr_mean,peakfluxdiff_iqr_half);
		gPeakFluxResoVSSNRGraph->SetPointError(pntCounter,0,peakfluxdiff_iqr_half_err);
	
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

	//Compute reliability
	cout<<"INFO: Compute Reliability ..."<<endl;
	gReliability= new TEfficiency(*gFluxHisto_rec_reliability,*gFluxHisto_reliability); 
	gReliability->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gReliability->SetConfidenceLevel(0.68);

	gReliabilityGraph= gReliability->CreateGraph(); 
	gReliabilityGraph->SetNameTitle("ReliabilityGraph","ReliabilityGraph");

	//Compute false detection rate
	cout<<"INFO: Compute False Detection Rate ..."<<endl;
	gFalseDetectionRate= new TEfficiency(*gFluxHisto_rec_fdr,*gFluxHisto_fdr); 
	gFalseDetectionRate->SetStatisticOption(TEfficiency::kFCP);  // to set option for errors (see ref doc)
	gFalseDetectionRate->SetConfidenceLevel(0.68);

	gFalseDetectionRateGraph= gFalseDetectionRate->CreateGraph(); 
	gFalseDetectionRateGraph->SetNameTitle("FalseDetectionRateGraph","FalseDetectionRateGraph");

	return 0;

}//close AnalyzeData()



int ReadCatalogData()
{
	int nFiles= (int)(gFileNames.size());

	for(int i=0;i<nFiles;i++)
	{
		//- Read reference catalog data
		std::vector<SourceInfo> sources;
		ReadData(sources,gFileNames[i],ParseData,gSkipLinePatterns);
		#ifdef LOGGING_ENABLED	
			INFO_LOG("#"<<sources.size()<<" sources found in reference catalog file "<<gFileNames[i]<<" ...");
		#endif
		
		//- Read rec catalog data
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

		//- Find matches
		FindSourceMatch(sources,sources_rec);

	}//end loop files

	return 0;

}//close ReadCatalogData()



int ReadFileLists()
{
	//Clear lists
	gFileNames_rec.clear();
	gFileNames.clear();

	//Check filenames
	std::ifstream fileStream(fileName);	
	std::ifstream fileStream_rec(fileName_rec);
	std::string line;
	std::string filename= "";
	if (fileStream.fail() || !fileStream.is_open()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<fileName<<" for reading!");
		#endif
		return -1;
	}
	if (fileStream_rec.fail() || !fileStream_rec.is_open()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<fileName_rec<<" for reading!");
		#endif
		return -1;
	}

	//Store filenames present in lists
	std::vector<std::string> fileNames;
	while (std::getline(fileStream, line)) {
   	std::istringstream iss(line);
   	if (!(iss >> filename)) { 
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read line from file "<<fileName<<"!");
			#endif
			return -1; 
		}
   	gFileNames.push_back(filename);
	}//end file read
				
	while (std::getline(fileStream_rec, line)) {
   	std::istringstream iss(line);
   	if (!(iss >> filename)) { 
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read line from file "<<fileName_rec<<"!");
			#endif
			return -1; 
		}
   	gFileNames_rec.push_back(filename);
	}//end file read

	//Check files have same number of entries
	if(gFileNames.size()!=gFileNames_rec.size()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Input filelist must have the same number of entries!");
		#endif
		return -1;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<gFileNames.size()<<" files to be processed...");
	#endif
	
	return 0;

}//close ReadFileLists()

int ReadRegions()
{
	//Read and parse DS9 regions
	if (DS9RegionParser::Parse(gRegions,gRegionFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse DS9 region file "<<gRegionFileName<<"!");
		#endif
		return -1;
	}

	//Check if empty
	if(gRegions.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No regions read from file "<<gRegionFileName<<"!");
		#endif
		return -1;
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<gRegions.size()<<" regions read from file "<<gRegionFileName<<"...");
		#endif
	}

	//Print regions
	for(size_t i=0;i<gRegions.size();i++){
		cout<<"--> DS9 Region no. "<<i+1<<endl;
		gRegions[i]->Print();
	}

	return 0;

}//close ReadRegions()


int ReadMap()
{
	//Read image
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading map from file "<<gMapFileName<<" ...");
	#endif
	Image* img= new Image;
	if(img->ReadFITS(gMapFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read map from file "<<gMapFileName<<"!");
		#endif
		delete img;
		img= 0;
		return -1;
	}

	//Get metadata
	ImgMetaData* metadata= img->GetMetaData();
	if(!metadata){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get image metadata!");
		#endif
		delete img;
		img= 0;
		return -1;
	}
			
	//Get WCS
	gGalWCS= metadata->GetWCS(eGALACTIC);
	if(!gGalWCS){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get galactic WCS!");
		#endif
		return -1;
	}
	
	gWCS= metadata->GetWCS(eJ2000);
	if(!gWCS){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get J2000 WCS!");
		#endif
		return -1;
	}

	return 0;

}//close ReadMap()


void Save()
{
	//Save file
	if(gOutputFile)
	{
		gOutputFile->cd();
		if(gMatchDataTree) gMatchDataTree->Write();

		if(gFluxHisto) gFluxHisto->Write();
		if(gFluxHisto_reliability) gFluxHisto_reliability->Write();			
		if(gFluxHisto_fdr) gFluxHisto_fdr->Write();	
		if(gFluxHisto_rec) gFluxHisto_rec->Write();
		if(gFluxHisto_rec_reliability) gFluxHisto_rec_reliability->Write();
		if(gFluxHisto_rec_fdr) gFluxHisto_rec_fdr->Write();
		//if(gCompleteness) gCompleteness->Write();
		if(gCompletenessGraph) gCompletenessGraph->Write();
		//if(gReliability) gReliability->Write();
		if(gReliabilityGraph) gReliabilityGraph->Write();
		if(gFalseDetectionRate) gFalseDetectionRate->Write();
		if(gFalseDetectionRateGraph) gFalseDetectionRateGraph->Write();
		
		if(gFluxBiasGraph) gFluxBiasGraph->Write();
		if(gFluxResoGraph) gFluxResoGraph->Write();
		if(gLgFluxBiasGraph) gLgFluxBiasGraph->Write();
		if(gLgFluxResoGraph) gLgFluxResoGraph->Write();
		if(gFluxBiasVSRecFluxGraph) gFluxBiasVSRecFluxGraph->Write();
		if(gFluxResoVSRecFluxGraph) gFluxResoVSRecFluxGraph->Write();
		if(gLgFluxBiasVSRecFluxGraph) gLgFluxBiasVSRecFluxGraph->Write();
		if(gLgFluxResoVSRecFluxGraph) gLgFluxResoVSRecFluxGraph->Write();
		if(gPosXBiasGraph) gPosXBiasGraph->Write();
		if(gPosYBiasGraph) gPosYBiasGraph->Write();
		if(gPosXResoGraph) gPosXResoGraph->Write();
		if(gPosYResoGraph) gPosYResoGraph->Write();

		if(gFluxBiasVSSNRGraph) gFluxBiasVSSNRGraph->Write();
		if(gFluxResoVSSNRGraph) gFluxResoVSSNRGraph->Write();
		if(gPeakFluxBiasVSSNRGraph) gPeakFluxBiasVSSNRGraph->Write();
		if(gPeakFluxResoVSSNRGraph) gPeakFluxResoVSSNRGraph->Write();
		if(gPosXBiasVSSNRGraph) gPosXBiasVSSNRGraph->Write();
		if(gPosYBiasVSSNRGraph) gPosYBiasVSSNRGraph->Write();
		if(gPosXResoVSSNRGraph) gPosXResoVSSNRGraph->Write();
		if(gPosYResoVSSNRGraph) gPosYResoVSSNRGraph->Write();

		gOutputFile->Close();
	}
	

}//close Save()



int ParseData_caesar(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check catalog version
	int NDataCols= 0;
	if(gCatalogVersion==1) NDataCols= gNDataCols_caesar;
	else if(gCatalogVersion==2) NDataCols= gNDataCols_caesar_v2; 
	else if(gCatalogVersion==3) NDataCols= gNDataCols_caesar_v3; 

	//Check number of fields present
	if(fields.size()<NDataCols){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of fields found in line ("<<fields.size()<<") when "<<NDataCols<<" expected!");
		#endif
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
	else if(gCatalogVersion==2 || gCatalogVersion==3){
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
	//cout<<"A="<<A<<", avgRMS="<<avgRMS<<", snr="<<snr<<endl;
	
	double ra= 0;
	double dec= 0;
	AstroUtils::PixelToWCSCoords(ra,dec,gWCS,x,y); 

	double l= 0;
	double b= 0;
	AstroUtils::PixelToWCSCoords(l,b,gGalWCS,x,y); 
	
	info.name= name;
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

	return 0;

}//close ParseData_caesar()


int ParseData_selavy(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols_selavy){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols_selavy<<" expected!");
		#endif
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

	return 0;

}//close ParseData_selavy()




int ParseData_aegean(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols_aegean){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols_aegean<<" expected!");
		#endif
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
	
	return 0;

}//close ParseData_aegean()

int ParseData(SourceInfo& info,const std::vector<std::string>& fields)
{
	//Check number of fields present
	if(fields.size()<gNDataCols){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid number of fields found in line ("<<fields.size()<<") when "<<gNDataCols<<" expected!");
		#endif
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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open config file "<<filename<<" for reading!");
		#endif
		return -1;
  }

	//Parsing file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading and parsing file: "<<filename);
	#endif

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
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to parse line "<<line_counter<<" of file "<<filename<<"!");
			#endif
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

