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
#include <SourceMatchData.h>

#include <ConfigParser.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>
#include <MathUtils.h>
#include <EllipseUtils.h>
#include <Contour.h>
#include <WCSUtils.h>
#include <DS9Region.h>
#include <DS9RegionParser.h>
#include <TCanvas.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TGraphErrors.h>


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
	cout<<"-i, --input=[INPUT_FILE] \t Input file in ROOT format to be cross-matched with catalogs"<<endl;
	cout<<"-C, --catalogs=[CATALOG_FILE] \t Input file list name (ascii format) containing all catalog files in ROOT format to be cross-matched each other"<<endl;
	cout<<"-B, --regions=[REGION_FILE] \t Input file list name (ascii format) containing all region files in DS9 format used as shuffle boundaries"<<endl;
	cout<<"-o, --output=[MATCH_OUTPUT_FILE] \t Output file name (root format) in which to store source match info"<<endl;
	cout<<"-J, --soutput=[OUTPUT_FILE] \t Output file name (root format) in which to store source catalog with spectral index info added"<<endl;
	cout<<"-j, --saveAllSources \t Save all sources (default: save only source matches)"<<endl;
	cout<<"-n, --nx=[NX] \t Number of divisions along X (default=4)"<<endl;
	cout<<"-N, --ny=[NY] \t Number of divisions along Y (default=4)"<<endl;
	cout<<"-m, --matchSourcesByFlux \t Match sources by flux (default=false)"<<endl;
	cout<<"-M, --matchFluxRelDiffThr \t Flux relative difference threshold (default=0.5)"<<endl;
	cout<<"-O, --matchByOverlap \t Match source islands by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-A, --matchOverlapThr \t Source island contour overlap fraction (wrt to total area) threshold (default=0.5)"<<endl;
	cout<<"-L, --applyAreaRatioThr \t If enabled and if source island is fully enclosed inside the region (or viceversa) a match requires that the source/region area ratio is higher than the overlap threshold. If disabled, assume a match whenever a source is fully enclosed in a region (or viceversa) (default=false)"<<endl;
	cout<<"-P, --matchByPos \t Match source islands by position centroid distance (default=false)"<<endl;
	cout<<"-T, --matchPosThr=[POS_THESHOLD] \t Source island centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	cout<<"-c, --matchSourceComponent \t Match source fit components by position centroid (if enabled) and ellipse overlap (if enabled) (default=false)"<<endl;
	cout<<"-p, --matchComponentByPos \t Match source fit components by position centroid distance (default=false)"<<endl;
	cout<<"-b, --matchComponentByOverlap \t Match source fit component by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-t, --compMatchPosThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	cout<<"-e, --compMatchPosHighThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=5)"<<endl;
	cout<<"-a, --compMatchOverlapThr \t Source fit component contour overlap fraction (wrt to total area) threshold (default=0.8)"<<endl;
	cout<<"-d, --compMatchOverlapLowThr \t Source fit component contour overlap fraction (wrt to total area) low threshold (default=0.2)"<<endl;
	cout<<"-l, --applyComponentAreaRatioThr \t If enabled and if source fit component is fully enclosed inside another (or viceversa) a match requires that the source1/source2 area ratio is higher than the component overlap threshold. If disabled, assume a match whenever a source component is fully enclosed in a source (or viceversa) (default=false)"<<endl;
	cout<<"-f, --filterByType \t Consider only true sources with given type when searching the match (default=no)"<<endl;
	cout<<"-s, --selectedType=[TYPE] \t True source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED) (default=-1)"<<endl;
	cout<<"-F, --filterBySimType \t Consider only true sources with given sim type when searching the match (default=no)"<<endl;
	cout<<"-S, --selectedSimType=[TYPE] \t True source sim types to be crossmatched (eRingLike=1,eBubbleLike=2,eEllipseLike=3,eDiskLike=4,eBlobLike=5) (default=-1)"<<endl;
	cout<<"-g, --shuffleSources \t Randomize sources and catalog sources in an annulus centred on their original position (default=no)"<<endl;
	cout<<"-r, --shuffleRmin=[Rmin] \t Minimum annulus radius for source shuffling (default=30)"<<endl;
	cout<<"-R, --shuffleRmax=[Rmax] \t Maximum annulus radius for source shuffling (default=50)"<<endl;
	cout<<"-u, --shiftFlux \t Shift flux density of main catalog (not that of the matching catalogs) (default=false)"<<endl;
	cout<<"-U, --fluxShift=[FLUX_SHIFT] \t Flux density shift (default=0)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "catalogs", required_argument, 0, 'C' },
	{ "regions", required_argument, 0, 'B' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", required_argument, 0, 'o' },
	{ "soutput", required_argument, 0, 'J' },	
	{ "saveAllSources", no_argument, 0, 'j' },	
	{ "nx", required_argument, 0, 'n' },
	{ "ny", required_argument, 0, 'N' },
	{ "matchSourcesByFlux", no_argument, 0, 'm'},
	{ "matchFluxRelDiffThr", required_argument, 0, 'M'},
	{ "matchByOverlap", no_argument, 0, 'O'},	
	{ "matchOverlapThr", required_argument, 0, 'A'},	
	{ "applyAreaRatioThr", no_argument, 0, 'L'},
	{ "matchByPos", no_argument, 0, 'P'},
	{ "matchPosThr", required_argument, 0, 'T'},
	{ "matchSourceComponent", no_argument, 0, 'c'},		
	{ "matchComponentByPos", no_argument, 0, 'p'},
	{ "matchComponentByOverlap", no_argument, 0, 'b'},	
	{ "compMatchPosThr", required_argument, 0, 't'},
	{ "compMatchPosHighThr", required_argument, 0, 'e'},
	{ "compMatchOverlapThr", required_argument, 0, 'a'},
	{ "compMatchOverlapLowThr", required_argument, 0, 'd'},
	{ "applyComponentAreaRatioThr", no_argument, 0, 'l'},	
	{ "minSourceMatchClusterSize",required_argument,0,'n'},
	{ "filterByType", no_argument, 0, 'f'},
	{ "filterBySimType", no_argument, 0, 'F'},
	{ "selectedType", required_argument, 0, 's'},	
	{ "selectedSimType", required_argument, 0, 'S'},
	{ "shuffleSources", no_argument, 0, 'g'},	
	{ "shuffleRmin", required_argument, 0, 'r'},
	{ "shuffleRmax", required_argument, 0, 'R'},
	{ "shiftFlux", no_argument, 0, 'u'},
	{ "fluxShift", required_argument, 0, 'U'},
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string catalogFileName= "";
std::string regionBoundaryFileName= "";
int verbosity= 4;//INFO level
bool matchSourcesByFlux= false;
float matchFluxRelDiffThr= 0.5;//50%
bool matchSourcesByOverlap= false;
float matchOverlapThr= 0.5;//fraction of overlap above which two sources are matched
bool matchSourcesByPos= false;
bool useWCSSimpleGeometry= true;
float matchPosThr= 2.5;//dist in arcsec below which two sources can be matched
bool matchSourceComponent= false;
bool matchSourceComponentsByPos= false;
float compMatchPosThr= 2.5;//dist in arcsec below which source fit component and region can be matched
float compMatchPosHighThr= 5;//dist in arcsec below which source fit component and region can be matched
bool matchSourceComponentsByOverlap= false;
float compMatchOverlapThr= 0.8;//fraction of overlap above which source fit component and region are matched
float compMatchOverlapLowThr= 0.2;//fraction of overlap above which source fit component and region are matched
bool applySourceComponentAreaRatioThr= false;
bool correctFlux= false;
bool shiftFlux= false;
double fluxShift= 0;
bool saveAllSources= false;
bool selectSourceByType= false;//default=all true sources searched 
std::vector<int> stypes;
std::vector<int> ssimtypes;
bool selectSourceBySimType= false;//default=all true sources searched 
bool applyFluxOverlapThreshold= false;
bool applySourceAreaRatioThr= false;
double fluxOverlapThr= 0.5;
float tileXmin= 0;
float tileXmax= 0;
float tileXstep= 0;
float tileYmin= 0;
float tileYmax= 0;
float tileYstep= 0;
long int nTilesX= 0;
long int nTilesY= 0;
bool shuffleSources= false;
double shuffleRMin= 30;
double shuffleRMax= 50;
bool readBoundaryRegions= false;


//Globar vars
TFile* inputFile= 0;
TFile* outputFile= 0;
std::string outputFileName= "sources.root";
TTree* outputTree= 0;
TFile* matchOutputFile= 0;
std::string matchOutputFileName= "MatchOutput.root";
TTree* matchedSourceInfo= 0;
SourceMatchData* sourceMatchData= 0;
Source* matchedSource= 0;
Source* inputSource= 0;
Source* m_source= 0;
std::vector<std::vector<DS9Region*>> m_regions;
std::vector<Contour*> m_regionContours;
std::vector<SourceMatchData*> sourceMatchDataCollection;
TTree* matchOptionTree= 0;
bool cloneSource= false;
SourceMatchPars* sourceMatchPars= 0;
SourceComponentMatchPars* sourceComponentMatchPars= 0;
int nTotSources= 0;
int nMatchedSources= 0;
int nTotSourceComponents= 0;
int nMatchedSourceComponents= 0;
int nCatalogs= 0;
std::vector<int> nMatchedSourcesPerCatalog;
std::vector<int> nMatchedSourceCatalogMultiplicity;
std::vector<int> nMatchedSourceComponentsPerCatalog;
std::vector<int> nMatchedSourceComponentCatalogMultiplicity;
std::vector<TCanvas*> m_SEDPlots;
std::vector<TGraphErrors*> m_SEDGraphs;
std::vector<TGraphErrors*> m_SEDComponentGraphs;


//====================================
//     ComponentPars struct
//====================================
struct ComponentPars {
	ComponentPars(){	
		fitComponentIndex= -1;
		sname= "";
		catalogIndex= -1;
		fitEllipse= 0;
		A= 0;
		fluxDensity= 0;
		isSelected= true;
	}
	~ComponentPars(){
		CodeUtils::DeletePtr<TEllipse>(fitEllipse);
	}
	int fitComponentIndex;
	int catalogIndex;
	std::string sname;
	double A;//fitted amplitude
	double fluxDensity;//fitted flux density
	TEllipse* fitEllipse;
	bool isSelected;//is selected component
};


//====================================
//     SourcePars struct
//====================================
struct SourcePars 
{
	SourcePars(int _catalogIndex,int _sourceIndex,int _nestedSourceIndex=-1)
		: catalogIndex(_catalogIndex), sourceIndex(_sourceIndex), nestedSourceIndex(_nestedSourceIndex)
	{
		S= 0;
		fluxDensity= 0;
		contour= 0;
		X0_wcs= 0;
		Y0_wcs= 0;
		componentPars.clear();
	}
	~SourcePars(){
		CodeUtils::DeletePtr<Contour>(contour);
		CodeUtils::DeletePtrCollection<ComponentPars>(componentPars);
	}
	void AddComponentPars(ComponentPars* pars){
		componentPars.push_back(pars);
	}

	int catalogIndex;
	int sourceIndex;
	int nestedSourceIndex;
	std::string sname;
	double X0_wcs;
	double Y0_wcs;
	Contour* contour;
	double S;//pixel-integrated flux
	double fluxDensity;//fitted flux density
	std::vector<ComponentPars*> componentPars;
};

//====================================
//     TileData struct
//====================================
struct TileData 
{
	//Constructor
	TileData(float _xmin,float _xmax,float _ymin,float _ymax)
		: xmin(_xmin),xmax(_xmax),ymin(_ymin),ymax(_ymax)
	{
		neighborTileIndexes.clear();
		sourcePars.clear();
		sourceCatalogPars.clear();
	}

	//Destructor 
	~TileData()
	{
		CodeUtils::DeletePtrCollection<SourcePars>(sourcePars);
		for(size_t i=0;i<sourceCatalogPars.size();i++){
			for(size_t j=0;j<sourceCatalogPars[i].size();j++){
				if(sourceCatalogPars[i][j]){
					delete sourceCatalogPars[i][j];
					sourceCatalogPars[i][j]= 0;
				}
			}
		}	
	}

	//Check if tile 
	bool IsNeighbor(TileData* aTile){
		if(!aTile) return false;
		float xmin_N= aTile->xmin;
		float xmax_N= aTile->xmax;
		float ymin_N= aTile->ymin;
		float ymax_N= aTile->ymax;

		//Consider 3x3 box
		bool isLeftRight= (xmin_N==xmax || xmax_N==xmin) && (ymin_N==ymin && ymax==ymax);
		bool isTopBottom= (ymin_N==ymax || ymax_N==ymin) && (xmin_N==xmin && xmax==xmax);
		bool isMajorDiag= (xmin_N==xmax && ymax_N==ymin) || (xmax_N==xmin && ymin_N==ymax);
		bool isMinorDiag= (xmin_N==xmax && ymin_N==ymax) || (xmax_N==xmax && ymax_N==ymin);
		bool isAdjacent= isLeftRight || isTopBottom || isMajorDiag || isMinorDiag;
		return isAdjacent;
	}

	//Add neighbor
	void AddNeighbor(size_t index){
		neighborTileIndexes.push_back(index);
	}

	//Add source pars
	int AddSourcePars(SourcePars* pars){
		if(!pars) return -1;
		sourcePars.push_back(pars);
		return 0;
	}

	//Get number of sources pars in this tile
	long int GetNSourcePars(){
		return static_cast<long int>(sourcePars.size());
	}

	//Init source catalog pars
	int InitCatalogSourcePars(int NCatalogs){
		if(NCatalogs<=0) return -1;
		sourceCatalogPars.clear();
		for(int i=0;i<NCatalogs;i++) sourceCatalogPars.push_back( std::vector<SourcePars*>() );
		return 0;
	}

	//Add source catalog pars
	int AddCatalogSourcePars(int catalogIndex,SourcePars* pars){
		if(!pars) return -1;
		if(catalogIndex<0 || catalogIndex>=(int)(sourceCatalogPars.size())) return -1;
		sourceCatalogPars[catalogIndex].push_back(pars);
		return 0;
	}

	//Get number of sources pars in this tile
	long int GetNSourceCatalogs(){
		return static_cast<long int>(sourceCatalogPars.size());
	}

	//Get number of sources catalog pars in this tile
	long int GetNSourceParsInCatalog(int catalogIndex){
		if(catalogIndex<0 || catalogIndex>=(int)(sourceCatalogPars.size())) return 0;
		return static_cast<long int>(sourceCatalogPars[catalogIndex].size());
	}

	//Get source catalog pars
	SourcePars* GetSourceCatalogPars(int catalogIndex,int cindex)
	{
		if(catalogIndex<0 || catalogIndex>=(int)(sourceCatalogPars.size())) return nullptr;
		if(cindex<0 || cindex>=(int)(sourceCatalogPars[catalogIndex].size())) return nullptr;	
		return sourceCatalogPars[catalogIndex][cindex];
	}

	//Tile coordinates
	float xmin;
	float xmax;
	float ymin;
	float ymax;

	//List of neighbors tiles
	std::vector<size_t> neighborTileIndexes;

	//Source pars for this tile
	std::vector<SourcePars*> sourcePars;

	//Source catalog pars for this tile
	std::vector<std::vector<SourcePars*>> sourceCatalogPars;

};//close TileData()

std::vector<TileData*> tileDataList;

std::vector<Source*> m_sources;
std::vector<SourcePars*> m_source_pars;

std::vector<std::vector<Source*>> m_catalog_sources;
std::vector<std::vector<SourcePars*>> m_catalog_source_pars;




//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int InitGrid(float xmin,float xmax,float xstep,float ymin,float ymax,float ystep);
int FindSourceMatchesInTiles();
int ComputeSpectralIndices();
int ReadData(std::string filename);
int ReadCatalogData(std::string filelist);
int ReadSourceData(std::string filename,int collectionIndex);
int FillSourcePars(std::vector<SourcePars*>& pars,Source* aSource,int collectionIndex,int sourceIndex,int nestedSourceIndex=-1,WCS* wcs=0,int coordSystem=0,bool shuffleData=false,double offsetX=0,double offsetY=0);
int ReadRegionData(std::string filelist);
int ReadCatalogRegionData(std::string filename,int catalogIndex);
bool HaveSourceComponentMatch(SourceComponentMatchPars& scompmatchpars,ComponentPars* pars1,ComponentPars* pars2);
bool HaveSourceIslandMatch(SourceMatchPars& smatchpars,SourcePars* pars1,SourcePars* pars2);
int GetRandPosInAnnulus(double& offsetX,double& offsetY);
int GetRandPosInRegion(double& offsetX,double& offsetY,Contour* contour,double xmin,double xmax,double ymin,double ymax,int nMaxGen=1000);
void Save();
int SaveSED();
int SaveSources();
int CloneObjectsInFile(std::vector<std::string> excludedObjNames);
int Init();

int main(int argc, char *argv[])
{
	//================================
	//== Parse command line options
	//================================
	if(ParseOptions(argc,argv)<0){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Failed to parse command line options!");
		#endif
		return -1;
	}
	
	//=======================
	//== Init
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
	//== Read data
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


	//==================================
	//== Read catalogs boundary regions
	//==================================
	if(readBoundaryRegions){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading catalog region boundary data from given files...");
		#endif
		if(ReadRegionData(regionBoundaryFileName)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Reading of catalog region data failed!");
			#endif
			return -1;
		}
	}

	//=======================
	//== Read catalogs data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading catalog source data from given files...");
	#endif
	if(ReadCatalogData(catalogFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of source catalog data failed!");
		#endif
		return -1;
	}


	//=======================
	//== Find source matches
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Correlating source catalogs to find matches...");
	#endif
	if(FindSourceMatchesInTiles()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Correlating source data to find matches failed!");
		#endif
		return -1;
	}

	//=============================
	//== Compute Spectral Indices
	//=============================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing spectral indices of matched sources ...");
	#endif
	if(ComputeSpectralIndices()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Computing spectral indices of matched sources failed!");
		#endif
		return -1;
	}
	

	//=======================
	//== Save
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving data to file ...");
	#endif
	Save();

	#ifdef LOGGING_ENABLED
		INFO_LOG("End source match finder");
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

	while((c = getopt_long(argc, argv, "hi:C:o:v:n:N:mM:POLlA:T:t:e:pba:d:cfFs:S:x:X:w:y:Y:k:gr:R:jJ:B:uU:",options_tab, &option_index)) != -1) {
    
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
			case 'C':	
			{
				catalogFileName= std::string(optarg);	
				break;	
			}
			case 'B':	
			{
				regionBoundaryFileName= std::string(optarg);	
				readBoundaryRegions= true;
				break;	
			}
			case 'o':	
			{
				matchOutputFileName= std::string(optarg);	
				break;	
			}
			case 'J':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'j':
			{
				saveAllSources= true;
				break;
			}
			case 'n':	
			{
				nTilesX= atol(optarg);
				break;	
			}
			case 'N':	
			{
				nTilesY= atol(optarg);
				break;	
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			case 'm':
			{
				matchSourcesByFlux= true;
				break;
			}
			case 'M':
			{
				matchFluxRelDiffThr= atof(optarg);
				break;
			}
			case 'O':
			{
				matchSourcesByOverlap= true;
				break;
			}
			case 'A':
			{
				matchOverlapThr= atof(optarg);
				break;
			}	
			case 'P':
			{
				matchSourcesByPos= true;
				break;
			}
			case 'T':
			{
				matchPosThr= atof(optarg);
				break;
			}
			case 'c':
			{
				matchSourceComponent= true;
				break;
			}
			case 'p':
			{
				matchSourceComponentsByPos= true;
				break;
			}	
			case 't':
			{
				compMatchPosThr= atof(optarg);
				break;
			}
			case 'e':
			{
				compMatchPosHighThr= atof(optarg);
				break;
			}	
			case 'b':
			{
				matchSourceComponentsByOverlap= true;
				break;
			}
			case 'a':
			{
				compMatchOverlapThr= atof(optarg);
				break;
			}	
			case 'd':
			{
				compMatchOverlapLowThr= atof(optarg);
				break;
			}	
			case 'f':
			{
				selectSourceByType= true;
				break;
			}
			case 's':	
			{
				int stype= atoi(optarg);
				stypes.push_back(stype);	
				break;	
			}	
			case 'F':
			{
				selectSourceBySimType= true;
				break;
			}
			case 'S':
			{
				int ssimtype= atoi(optarg);
				ssimtypes.push_back(ssimtype);	
				break;
			}	
			
			case 'L':
			{
				applySourceAreaRatioThr= true;
				break;
			}	
			case 'l':
			{
				applySourceComponentAreaRatioThr= true;
				break;
			}
			case 'g':
			{
				shuffleSources= true;
				break;
			}
			case 'r':
			{
				shuffleRMin= atof(optarg);
				break;
			}
			case 'R':
			{
				shuffleRMax= atof(optarg);
				break;
			}	
			case 'u':
			{
				shiftFlux= true;
				break;
			}
			case 'U':
			{
				fluxShift= atof(optarg);
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
	//== Init Logger 
	//=======================
	//## Set logging level
	std::string sloglevel= GetStringLogLevel(verbosity);
	#ifdef LOGGING_ENABLED
		LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	#endif

	
	//=======================
	//== Print options
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("matchSourcesByPos? "<<matchSourcesByPos<<", matchPosThr(arcsec)="<<matchPosThr);
		INFO_LOG("matchSourcesByOverlap? "<<matchSourcesByOverlap<<", matchOverlapThr="<<matchOverlapThr);
		INFO_LOG("applySourceAreaRatioThr? "<<applySourceAreaRatioThr);
		INFO_LOG("matchSourcesByFlux? "<<matchSourcesByFlux<<", matchFluxRelDiffThr="<<matchFluxRelDiffThr);
		INFO_LOG("matchSourceComponent? "<<matchSourceComponent);
		INFO_LOG("matchSourceComponentsByPos? "<<matchSourceComponentsByPos<<", compMatchPosThr(arcsec)="<<compMatchPosThr<<", compMatchPosHighThr(arcsec)="<<compMatchPosHighThr);
		INFO_LOG("matchSourceComponentsByOverlap? "<<matchSourceComponentsByOverlap<<", compMatchOverlapThr="<<compMatchOverlapThr<<", compMatchOverlapLowThr="<<compMatchOverlapLowThr);
		INFO_LOG("applySourceComponentAreaRatioThr? "<<applySourceComponentAreaRatioThr);
		INFO_LOG("shiftFlux? "<<shiftFlux<<", shift="<<fluxShift);
	#endif

	return 0;

}//close ParseOptions()

int Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");
	
	//Init output source TTree
	if(!outputTree) outputTree= new TTree("SourceInfo","SourceInfo");
	m_source= 0;
	outputTree->Branch("Source",&m_source);

	//Open match output file
	if(!matchOutputFile) matchOutputFile= new TFile(matchOutputFileName.c_str(),"RECREATE");
	
	//Init source cube TTree
	if(!matchedSourceInfo) matchedSourceInfo= new TTree("SourceMatchInfo","SourceMatchInfo");
	matchedSourceInfo->Branch("SourceMatchData",&sourceMatchData);
	sourceMatchDataCollection.clear();

	//Init match option TTree
	if(!matchOptionTree) matchOptionTree= new TTree("MatchSummary","MatchSummary");
	matchOptionTree->Branch("matchSourcesByFlux",&matchSourcesByFlux);
	matchOptionTree->Branch("matchFluxRelDiffThr",&matchFluxRelDiffThr);
	matchOptionTree->Branch("matchSourcesByOverlap",&matchSourcesByOverlap);
	matchOptionTree->Branch("matchOverlapThr",&matchOverlapThr);
	matchOptionTree->Branch("matchSourcesByPos",&matchSourcesByPos);
	matchOptionTree->Branch("matchPosThr",&matchPosThr);
	matchOptionTree->Branch("matchSourceComponent",&matchSourceComponent);
	matchOptionTree->Branch("matchSourceComponentsByPos",&matchSourceComponentsByPos);
	matchOptionTree->Branch("compMatchPosThr",&compMatchPosThr);
	matchOptionTree->Branch("compMatchPosHighThr",&compMatchPosHighThr);
	matchOptionTree->Branch("matchSourceComponentsByOverlap",&matchSourceComponentsByOverlap);
	matchOptionTree->Branch("compMatchOverlapThr",&compMatchOverlapThr);
	matchOptionTree->Branch("compMatchOverlapLowThr",&compMatchOverlapLowThr);
	matchOptionTree->Branch("selectSourceByType",&selectSourceByType);
	matchOptionTree->Branch("selectSourceBySimType",&selectSourceBySimType);
	matchOptionTree->Branch("applySourceAreaRatioThr",&applySourceAreaRatioThr);
	matchOptionTree->Branch("applySourceComponentAreaRatioThr",&applySourceComponentAreaRatioThr);

	matchOptionTree->Branch("nTotSources",&nTotSources);
	matchOptionTree->Branch("nMatchedSources",&nMatchedSources);
	matchOptionTree->Branch("nTotSourceComponents",&nTotSourceComponents);
	matchOptionTree->Branch("nMatchedSourceComponents",&nMatchedSourceComponents);
	matchOptionTree->Branch("nCatalogs",&nCatalogs);
	matchOptionTree->Branch("nMatchedSourcesPerCatalog",&nMatchedSourcesPerCatalog);
	matchOptionTree->Branch("nMatchedSourceCatalogMultiplicity",&nMatchedSourceCatalogMultiplicity);
	matchOptionTree->Branch("nMatchedSourceComponentsPerCatalog",&nMatchedSourceComponentsPerCatalog);
	matchOptionTree->Branch("nMatchedSourceComponentCatalogMultiplicity",&nMatchedSourceComponentCatalogMultiplicity);
	matchOptionTree->Branch("shiftFlux",&shiftFlux);
	matchOptionTree->Branch("fluxShift",&fluxShift);

	//Set random seed
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	
	return 0;

}//close Init()


int InitGrid(float xmin,float xmax,float xstep,float ymin,float ymax,float ystep)
{
	//Create 2d grid
	#ifdef LOGGING_ENABLED
		INFO_LOG("Creating 2d tile grid (x min/max/step="<<xmin<<"/"<<xmax<<"/"<<xstep<<", y min/max/step="<<ymin<<"/"<<ymax<<"/"<<ystep);
	#endif
	std::vector<float> ix_min;
	std::vector<float> ix_max;
	std::vector<float> iy_min;
	std::vector<float> iy_max;
	MathUtils::Compute2DFloatGrid(ix_min,ix_max,iy_min,iy_max,xmin,xmax,xstep,ymin,ymax,ystep);
	nTilesX= static_cast<long int>(ix_min.size());
	nTilesY= static_cast<long int>(iy_min.size());
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<nTilesX<<" x "<<nTilesY<<" tiles created...");
	#endif

	TileData* aTileData= 0;
	for(size_t j=0;j<iy_min.size();j++){
		float ymin= iy_min[j];
		float ymax= iy_max[j];

		for(size_t i=0;i<ix_min.size();i++){
			float xmin= ix_min[i];
			float xmax= ix_max[i];
	
			aTileData= new TileData(xmin,xmax,ymin,ymax);
			tileDataList.push_back(aTileData);
		}//end loop x
	}//end loop y

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<tileDataList.size()<<" tile initialized...");
	#endif

	//Fill neighbors id per tile
	for(size_t i=0;i<tileDataList.size()-1;i++){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Tile no. "<<i<<": xrange["<<tileDataList[i]->xmin<<","<<tileDataList[i]->xmax<<"], yrange=["<<tileDataList[i]->ymin<<","<<tileDataList[i]->ymax<<"]");
		#endif

		for(size_t j=i+1;j<tileDataList.size();j++){
			bool areNeighbors= tileDataList[i]->IsNeighbor(tileDataList[j]);
			if(areNeighbors){
				#ifdef LOGGING_ENABLED
					INFO_LOG("Tiles ("<<i<<","<<j<<") are neighbors");
				#endif
				tileDataList[i]->AddNeighbor(j);
				tileDataList[j]->AddNeighbor(i);		
			}
		}//end loop tiles next
	}//end loop tiles

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Tile no. "<<tileDataList.size()-1<<": xrange["<<tileDataList[tileDataList.size()-1]->xmin<<","<<tileDataList[tileDataList.size()-1]->xmax<<"], yrange="<<tileDataList[tileDataList.size()-1]->ymin<<","<<tileDataList[tileDataList.size()-1]->ymax<<"]");
	#endif

	return 0;

}//close InitGrid()


int ComputeSpectralIndices()
{
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing spectral indices...");
	#endif

	TCanvas* Plot= 0;

	for(size_t i=0;i<sourceMatchDataCollection.size();i++)
	{
		//Get source
		Source* source= sourceMatchDataCollection[i]->GetSource();

		//Skip if no match found
		if(!sourceMatchDataCollection[i]->HasSourceMatch()) continue;

		//Compute source SED and spectral index
		int status= sourceMatchDataCollection[i]->ComputeSourceSEDs();
		if(status==0){//Set spectral index in source
			SpectralIndexData spectralIndexData= sourceMatchDataCollection[i]->GetSpectralIndexData();
			source->SetSpectralIndexData(spectralIndexData);

			//Save sed graph
			TGraphErrors* sed= sourceMatchDataCollection[i]->GetSEDGraph();
			m_SEDGraphs.push_back(sed);
		}
		else {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute spectral indices for source match data no. "<<i+1<<"!");
			#endif
		}

		//Compute source component SEDs and relative spectral indices
		status= sourceMatchDataCollection[i]->ComputeSourceComponentSEDs();
		if(status==0){
			std::vector<SpectralIndexData> componentSpectralIndexData= sourceMatchDataCollection[i]->GetSourceComponentSpectralIndexData();
			source->SetComponentSpectralIndexData(componentSpectralIndexData);

			//Save sed component graph
			std::vector<TGraphErrors*> componentSEDs= sourceMatchDataCollection[i]->GetSEDComponentGraphs();
			for(size_t k=0;k<componentSEDs.size();k++) m_SEDComponentGraphs.push_back(componentSEDs[k]);

		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute source component spectral indices for source match data no. "<<i+1<<"!");
			#endif
		}

		//Get SED+fit plot
		//Plot= sourceMatchDataCollection[i]->DrawSED();
		//if(Plot) m_SEDPlots.push_back(Plot);

	}//end loop match data

	return 0;


}//close ComputeSpectralIndices()




int FindSourceMatchesInTiles()
{
	//####################################
	//##     METHOD
	//####################################
	//1) Compare source by contour (centroid, overlap area)
	//2) Compare source component ellipses (centroid, overlap)
	//####################################
	nTotSources= 0;
	nMatchedSources= 0;
	nTotSourceComponents= 0;
	nMatchedSourceComponents= 0;
	

	//Loop over tiles and find matches with regions in the same tile and in neighbor tiles
	for(size_t i=0;i<tileDataList.size();i++)
	{
		long int NSourcePars= tileDataList[i]->GetNSourcePars();
		long int NCatalogs= tileDataList[i]->GetNSourceCatalogs();
		std::vector<size_t> neighbors= tileDataList[i]->neighborTileIndexes;

		//Init counters
		if(nMatchedSourcesPerCatalog.empty()){
			nMatchedSourcesPerCatalog= std::vector<int>(NCatalogs,0);
		}
		if(nMatchedSourceCatalogMultiplicity.empty()){
			nMatchedSourceCatalogMultiplicity= std::vector<int>(NCatalogs+1,0);
		}
		if(nMatchedSourceComponentsPerCatalog.empty()){
			nMatchedSourceComponentsPerCatalog= std::vector<int>(NCatalogs,0);
		}
		if(nMatchedSourceComponentCatalogMultiplicity.empty()){
			nMatchedSourceComponentCatalogMultiplicity= std::vector<int>(NCatalogs+1,0);
		}

		for(long int j=0;j<NSourcePars;j++) 
		{		
			//Get source info
			SourcePars* sourcePars= (tileDataList[i]->sourcePars)[j];
			std::vector<ComponentPars*> componentPars= sourcePars->componentPars;
			long int sourceIndex= sourcePars->sourceIndex;
			long int nestedSourceIndex= sourcePars->nestedSourceIndex;
			Source* source= 0;
			if(nestedSourceIndex==-1) source= m_sources[sourceIndex];
			else source= m_sources[sourceIndex]->GetNestedSource(nestedSourceIndex);
			
			if(!source){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<") found, this should not occur, skip source!");
				#endif
				continue;
			}

			#ifdef LOGGING_ENABLED
			if(j%100==0) INFO_LOG("#"<<j+1<<"/"<<NSourcePars<<" sources to be matched in tile no. "<<i+1<<" against "<<NCatalogs<<" catalogs ...");
			#endif

			//Create source match data for this source
			inputSource= new Source;
			*inputSource= *source;
			sourceMatchData= new SourceMatchData(inputSource,NCatalogs,cloneSource);
			sourceMatchData->ShiftFlux(shiftFlux,fluxShift);
	
			//Get source fit info
			bool hasFitInfo= source->HasFitInfo();
			int nFitComponents= 0;
			int nSelFitComponents= 0;
			if(hasFitInfo) {
				SourceFitPars fitPars= source->GetFitPars();
				nFitComponents= fitPars.GetNComponents();
				nSelFitComponents= fitPars.GetNSelComponents();
			}

			//Count total number of sources & components					
			nTotSources++;
			nTotSourceComponents+= nSelFitComponents;
			std::vector<int> sourceMatchCounter(NCatalogs,0);
			std::vector<std::vector<int>> compMatchCounter;
			for(size_t c=0;c<componentPars.size();c++){
				compMatchCounter.push_back(std::vector<int>());
				for(int p=0;p<NCatalogs;p++){
					compMatchCounter[c].push_back(0);
				}
			}

			//==============================================
			//===        CATALOG LOOP
			//==============================================
			//Loop over source catalogs in the same tile
			for(long int s=0;s<NCatalogs;s++)
			{
				int nMatches= 0;	
									
				long int NSourceParsInCatalog= tileDataList[i]->GetNSourceParsInCatalog(s);
				if(NSourceParsInCatalog<=0) continue;

				#ifdef LOGGING_ENABLED
					DEBUG_LOG("#"<<NSourcePars<<" sources to be cross-matched against #"<<NSourceParsInCatalog<<" sources (catalog #"<<s+1<<") in tile no. "<<i+1<<" ...");
				#endif

				//==============================================
				//===      CATALOG SOURCE LOOP (SAME TILE)
				//==============================================	
				for(long int k=0;k<NSourceParsInCatalog;k++)
				{
					//Get catalog source
					SourcePars* sourceCatalogPars= tileDataList[i]->GetSourceCatalogPars(s,k);
					if(!sourceCatalogPars){
						#ifdef LOGGING_ENABLED
							ERROR_LOG("Null ptr to source catalog pars (catalog no. "<<s+1<<", source pars no. "<<k+1<<")!");
						#endif
						return -1;
					}
					long int sourceIndex_catalog= sourceCatalogPars->sourceIndex;
					long int nestedSourceIndex_catalog= sourceCatalogPars->nestedSourceIndex;
					Source* source_catalog= 0;
					if(nestedSourceIndex_catalog==-1) source_catalog= m_catalog_sources[s][sourceIndex_catalog];
					else source_catalog= m_catalog_sources[s][sourceIndex_catalog]->GetNestedSource(nestedSourceIndex_catalog);
			
					if(!source_catalog){
						#ifdef LOGGING_ENABLED
							WARN_LOG("No catalog source (catalog="<<s+1<<", index="<<sourceIndex_catalog<<", nestedIndex="<<nestedSourceIndex_catalog<<") found in tile no. "<<i+1<<", this should not occur, skip source!");
						#endif
						continue;
					}
		
					//Check source match
					SourceMatchPars smatchpars;
					bool hasMatch= HaveSourceIslandMatch(smatchpars,sourcePars,sourceCatalogPars);
					if(!hasMatch) {
						continue;
					}
						
					//Add matched source 
					matchedSource= new Source;
					*matchedSource= *source_catalog;
					int status= sourceMatchData->AddMatchedSourceToGroup(s,matchedSource,cloneSource);
					if(status<0){
						#ifdef LOGGING_ENABLED
							WARN_LOG("Failed to add matched source to group!");
						#endif		
						continue;
					}

					sourceMatchPars= new SourceMatchPars;
					*sourceMatchPars= smatchpars;
					if(sourceMatchData->AddMatchedSourcePars(s,sourceMatchPars)<0){
						#ifdef LOGGING_ENABLED
							WARN_LOG("Failed to add matched source pars!");
						#endif
						continue;
					}
						

					int matchedSourceIndex= nMatches; 
					nMatches++;
					sourceMatchCounter[s]++;
					
					
					//Search component matches	
					std::vector<ComponentPars*> sourceCatalogComponentPars= sourceCatalogPars->componentPars;

					for(size_t l=0;l<componentPars.size();l++)
					{
						ComponentPars* cpars= componentPars[l];
						int fitComponentIndex= cpars->fitComponentIndex;
				
						for(size_t t=0;t<sourceCatalogComponentPars.size();t++)
						{
							ComponentPars* cpars_catalog= sourceCatalogComponentPars[t];
							int fitComponentIndex_catalog= cpars_catalog->fitComponentIndex;
							
							SourceComponentMatchPars cmatchpars;
							bool hasComponentMatch= HaveSourceComponentMatch(cmatchpars,cpars,cpars_catalog);
							if(!hasComponentMatch) continue;
								
							#ifdef LOGGING_ENABLED
								INFO_LOG("SOURCE "<<inputSource->GetName()<<" (COMP"<<fitComponentIndex<<"): Matched source component (cat="<<s<<", sindex="<<matchedSourceIndex<<", cindex="<<fitComponentIndex_catalog<<")");
							#endif

							status= sourceMatchData->AddMatchedComponentIndex(fitComponentIndex,s,matchedSourceIndex,fitComponentIndex_catalog);
							if(status<0){
								#ifdef LOGGING_ENABLED
									WARN_LOG("Failed to add matched source component to group!");
								#endif	
								continue;
							}

							sourceComponentMatchPars= new SourceComponentMatchPars;
							*sourceComponentMatchPars= cmatchpars;
							status= sourceMatchData->AddMatchedSourceComponentPars(fitComponentIndex,s,sourceComponentMatchPars);
							if(status<0){
								#ifdef LOGGING_ENABLED
									WARN_LOG("Failed to add matched source component pars to group!");
								#endif	
								continue;
							}

							compMatchCounter[l][s]++;			
						}//end (t index) loop components in this catalog source 
					}//end (l index) loop source components
					

				}//end (k index) loop catalog sources in this tile


				//==============================================
				//===     CATALOG SOURCE LOOP (NEIGHBOUR TILES)
				//==============================================	
				for(size_t n=0;n<neighbors.size();n++)	
				{
					size_t neighborIndex= neighbors[n];
					if(i==neighborIndex) continue;//skip same tile

					long int NSourceParsInCatalog_neighbor= tileDataList[neighborIndex]->GetNSourceParsInCatalog(s);
					if(NSourceParsInCatalog_neighbor<=0) continue;

					#ifdef LOGGING_ENABLED
						DEBUG_LOG("#"<<NSourcePars<<" sources to be cross-matched against #"<<NSourceParsInCatalog_neighbor<<" sources (catalog #"<<s+1<<") in tile no. "<<neighborIndex<<" ...");
					#endif	
				
					for(long int k=0;k<NSourceParsInCatalog_neighbor;k++)
					{
						//Get catalog source
						SourcePars* sourceCatalogPars= tileDataList[neighborIndex]->GetSourceCatalogPars(s,k);

						if(!sourceCatalogPars){
							#ifdef LOGGING_ENABLED
								ERROR_LOG("Null ptr to source catalog pars (catalog no. "<<s+1<<", source pars no. "<<k+1<<"!");
							#endif
							return -1;
						}

						long int sourceIndex_catalog= sourceCatalogPars->sourceIndex;
						long int nestedSourceIndex_catalog= sourceCatalogPars->nestedSourceIndex;
						Source* source_catalog= 0;
						if(nestedSourceIndex_catalog==-1) source_catalog= m_catalog_sources[s][sourceIndex_catalog];
						else source_catalog= m_catalog_sources[s][sourceIndex_catalog]->GetNestedSource(nestedSourceIndex_catalog);
			
						if(!source_catalog){
							#ifdef LOGGING_ENABLED
								WARN_LOG("No catalog source (catalog="<<s+1<<", index="<<sourceIndex_catalog<<", nestedIndex="<<nestedSourceIndex_catalog<<") found in tile no. "<<i+1<<", this should not occur, skip source!");
							#endif
							continue;
						}

						//Check source match
						SourceMatchPars smatchpars;
						bool hasMatch= HaveSourceIslandMatch(smatchpars,sourcePars,sourceCatalogPars);
						if(!hasMatch) {
							continue;
						}
						
						//Add matched source & pars
						matchedSource= new Source;
						*matchedSource= *source_catalog;
						int status= sourceMatchData->AddMatchedSourceToGroup(s,matchedSource,cloneSource);
						if(status<0){
							#ifdef LOGGING_ENABLED
								WARN_LOG("Failed to add matched source to group!");
							#endif
							continue;
						}
		
						sourceMatchPars= new SourceMatchPars;
						*sourceMatchPars= smatchpars;
						if(sourceMatchData->AddMatchedSourcePars(s,sourceMatchPars)<0){
							#ifdef LOGGING_ENABLED
								WARN_LOG("Failed to add matched source pars!");
							#endif
							continue;
						}
						

						int matchedSourceIndex= nMatches; 
						nMatches++;	
						sourceMatchCounter[s]++;
		
						
						//Search component matches
						std::vector<ComponentPars*> sourceCatalogComponentPars= sourceCatalogPars->componentPars;

						for(size_t l=0;l<componentPars.size();l++){
							ComponentPars* cpars= componentPars[l];
							int fitComponentIndex= cpars->fitComponentIndex;
							
							for(size_t t=0;t<sourceCatalogComponentPars.size();t++){
								ComponentPars* cpars_catalog= sourceCatalogComponentPars[t];
								int fitComponentIndex_catalog= cpars_catalog->fitComponentIndex;
							
								SourceComponentMatchPars cmatchpars;
								bool hasComponentMatch= HaveSourceComponentMatch(cmatchpars,cpars,cpars_catalog);
								if(!hasComponentMatch) continue;
								
								#ifdef LOGGING_ENABLED
									INFO_LOG("SOURCE "<<inputSource->GetName()<<" (COMP"<<fitComponentIndex<<"): Matched source component (cat="<<s<<", sindex="<<matchedSourceIndex<<", cindex="<<fitComponentIndex_catalog<<")");
								#endif

								status= sourceMatchData->AddMatchedComponentIndex(fitComponentIndex,s,matchedSourceIndex,fitComponentIndex_catalog);
								if(status<0){
									#ifdef LOGGING_ENABLED
										WARN_LOG("Failed to add matched source component to group!");
									#endif	
									continue;
								}

								sourceComponentMatchPars= new SourceComponentMatchPars;
								*sourceComponentMatchPars= cmatchpars;
								status= sourceMatchData->AddMatchedSourceComponentPars(fitComponentIndex,s,sourceComponentMatchPars);
								if(status<0){
									#ifdef LOGGING_ENABLED
										WARN_LOG("Failed to add matched source component pars to group!");
									#endif	
									continue;
								}

								compMatchCounter[l][s]++;

							}//end loop components in this catalog source 
						}//end loop source components
						

					}//end (k index) loop sources in neighbor tiles (s catalog)
				}//end (n index) loop neighbor tiles

				#ifdef LOGGING_ENABLED
					INFO_LOG("Source "<<inputSource->GetName()<<": #"<<nMatches<<" matches found for catalog no. "<<s+1<<" ...");
				#endif

			}//end (s index) loop catalogs 


			//Add source match data to collection
			if(sourceMatchData->HasSourceMatch()){
				nMatchedSources++;
			}
			
			if(saveAllSources || sourceMatchData->HasSourceMatch()){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Add source match data for source "<<inputSource->GetName()<<" to collection...");
				#endif
				sourceMatchDataCollection.push_back(sourceMatchData);
			}

			//Update counters
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Update source match counters...");
			#endif
			std::vector<int> sourceMatchMultiplicityCounter(NCatalogs,0);
			for(size_t c=0;c<sourceMatchCounter.size();c++){
				if(sourceMatchCounter[c]>0) {
					nMatchedSourcesPerCatalog[c]++;
					sourceMatchMultiplicityCounter[c]++;
				}
			}
			for(size_t c=0;c<compMatchCounter.size();c++){
				std::vector<int> catalogMultiplicityCounter(compMatchCounter[c].size(),0);
				for(size_t p=0;p<compMatchCounter[c].size();p++){
					if(compMatchCounter[c][p]>0) {
						nMatchedSourceComponentsPerCatalog[p]++;
						catalogMultiplicityCounter[p]++;
					}
				}
				int compMultiplicity= std::accumulate(catalogMultiplicityCounter.begin(), catalogMultiplicityCounter.end(), 0);
				if(compMultiplicity>0) nMatchedSourceComponents++;
				nMatchedSourceComponentCatalogMultiplicity[compMultiplicity]++;
			}
			
			//Compute catalog match multiplicity
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Compute catalog match multiplicity...");
			#endif
			int multiplicity= std::accumulate(sourceMatchMultiplicityCounter.begin(), sourceMatchMultiplicityCounter.end(), 0);
			nMatchedSourceCatalogMultiplicity[multiplicity]++;
	
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Source match multiplicity="<<multiplicity);
			#endif

			
		}//end loop source pars in this tile
	}//end loop tiles

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sourceMatchDataCollection.size()<<" source match data added to collection...");
	#endif
	
	cout<<"== MATCH RESULTS =="<<endl;
	cout<<"#matched/tot sources="<<nMatchedSources<<"/"<<nTotSources<<endl;
	for(size_t i=0;i<nMatchedSourcesPerCatalog.size();i++){
		cout<<"--> CAT"<<i+1<<": #match/tot="<<nMatchedSourcesPerCatalog[i]<<"/"<<nTotSources<<endl;
	}
	cout<<"source match multiplicity {";
	for(size_t i=0;i<nMatchedSourceCatalogMultiplicity.size()-1;i++){
		cout<<nMatchedSourceCatalogMultiplicity[i]<<",";
	}
	cout<<nMatchedSourceCatalogMultiplicity[nMatchedSourceCatalogMultiplicity.size()-1]<<"}"<<endl;
	cout<<"#matched/tot source components="<<nMatchedSourceComponents<<"/"<<nTotSourceComponents<<endl;
	for(size_t i=0;i<nMatchedSourceComponentsPerCatalog.size();i++){
		cout<<"--> CAT"<<i+1<<": #match/tot="<<nMatchedSourceComponentsPerCatalog[i]<<"/"<<nTotSourceComponents<<endl;
	}
	cout<<"source component match multiplicity {";
	for(size_t i=0;i<nMatchedSourceComponentCatalogMultiplicity.size()-1;i++){
		cout<<nMatchedSourceComponentCatalogMultiplicity[i]<<",";
	}
	cout<<nMatchedSourceComponentCatalogMultiplicity[nMatchedSourceComponentCatalogMultiplicity.size()-1]<<"}"<<endl;
	cout<<"===================="<<endl;


	//Fill Tree
	matchOptionTree->Fill();


	return 0;

}//close FindSourceMatchesInTiles()




bool HaveSourceIslandMatch(SourceMatchPars& smatchpars,SourcePars* pars1,SourcePars* pars2)
{
	//## Check input data
	if(!pars1 || !pars2){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given for contour!");
		#endif
		return false;
	}

	//## If sources are from the same collection return NO MATCH
	if(pars1->catalogIndex==pars2->catalogIndex) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Skip match as sources are from same collection (catalogIndex_1="<<pars1->catalogIndex<<", catalogIndex_2="<<pars2->catalogIndex<<")...");
		#endif
		return false;
	}

	//## Check contour
	Contour* cont1= pars1->contour;
	Contour* cont2= pars2->contour;
	if(!cont1 || !cont2){
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/both contour are nullptr!");
		#endif
		return false;
	}
	
	//## Check if contour pars were computed
	if(!cont1->HasParameters || !cont2->HasParameters){
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/both contours have no computed parameters!");
		#endif
		return false;
	}

	//## Fill match parameters
	// - Centroid dstance
	TVector2 centroid_1= cont1->Centroid;
	TVector2 centroid_2= cont2->Centroid;
	double posX_1= centroid_1.X();
	double posX_2= centroid_2.X();
	double posY_1= centroid_1.Y();
	double posY_2= centroid_2.Y();
	double posThr= matchPosThr/3600.;//convert threshold in deg (contour are given in deg)
	double posDist_Euclidean= MathUtils::GetEuclideanDist(posX_1,posY_1,posX_2,posY_2);
	double posDist_Haversine= AstroUtils::GetWCSPointDist_Haversine(posX_1,posY_1,posX_2,posY_2);
	double posDist= 0;
	if(useWCSSimpleGeometry) posDist= posDist_Euclidean;
	else posDist= posDist_Haversine;

	//- Contour overlaps
	double contArea_1= 0;
	double contArea_2= 0;
	int status= MathUtils::ComputeContourArea(contArea_1,cont1);
	if(status<0 && matchSourcesByOverlap){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute contour 1 area!");
		#endif
		return false;
	}
	status= MathUtils::ComputeContourArea(contArea_2,cont2);
	if(status<0 && matchSourcesByOverlap){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute contour 2 area!");
		#endif
		return false;
	} 
		
	double overlapArea= 0;
	int overlapFlag;
	status= MathUtils::ComputeContourOverlapArea(overlapArea,overlapFlag,cont1,cont2);
	if(status<0 && matchSourcesByOverlap){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute contour overlap area!");
		#endif
		return false;
	}
	double overlapAreaFraction_1= overlapArea/contArea_1;
	double overlapAreaFraction_2= overlapArea/contArea_2;

	// - Flux diff
	double flux_1= pars1->fluxDensity;
	double flux_2= pars2->fluxDensity;
	double fluxRelDiff= (flux_1-flux_2)/flux_2;
	
	smatchpars.posDist_Euclidean= posDist_Euclidean;
	smatchpars.posDist_Haversine= posDist_Haversine;
	smatchpars.contourOverlapArea= overlapArea;
	smatchpars.contourArea1= contArea_1;
	smatchpars.contourArea2= contArea_2;
	smatchpars.contourOverlapFlag= overlapFlag;
	smatchpars.contourOverlapAreaRatio1= overlapAreaFraction_1;
	smatchpars.contourOverlapAreaRatio2= overlapAreaFraction_2;
	smatchpars.fluxRelDiff= fluxRelDiff;



	//## Compare contour centroids 
	if(matchSourcesByPos){	
		
		if(fabs(posDist)>posThr) {	
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("NO POS MATCH: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<"), dist(arcsec)="<<posDist*3600.);
			#endif
			return false;
		}
		#ifdef LOGGING_ENABLED
			INFO_LOG("POS MATCH: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<"), dist(arcsec)="<<posDist*3600.);
		#endif
	}

	//## Compute contour overlap
	if(matchSourcesByOverlap){
		/*
		double contArea_1= 0;
		double contArea_2= 0;
		if(MathUtils::ComputeContourArea(contArea_1,cont1)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute contour 1 area!");
			#endif
			return false;
		}
		if(MathUtils::ComputeContourArea(contArea_2,cont2)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute contour 2 area!");
			#endif
			return false;
		} 
		double overlapArea= 0;
		int overlapFlag;
		if(MathUtils::ComputeContourOverlapArea(overlapArea,overlapFlag,cont1,cont2)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute contour overlap area!");
			#endif
			return false;
		}
		//Compute relative overlap areas
		double overlapAreaFraction_1= overlapArea/contArea_1;
		double overlapAreaFraction_2= overlapArea/contArea_2;

		*/
		
		//Check cases
		//1) Disjoint contours
		//2) Contour 1 inside contour 2 (overlap=Area1)
		//3) Contour 2 inside contour 1 (overlap=Area2)
		
		if(overlapFlag==eCONT_NOT_OVERLAPPING){
			#ifdef LOGGING_ENABLED	
				DEBUG_LOG("NO CONTOUR OVERLAP: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
			#endif
			return false;
		}

		else if(overlapFlag==eCONT1_INSIDE_CONT2){
			if(applySourceAreaRatioThr && overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT1 INSIDE CONT2 (OVERLAP BELOW THR): Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		else if(overlapFlag==eCONT2_INSIDE_CONT1){
			if(applySourceAreaRatioThr && overlapAreaFraction_1<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT2 INSIDE CONT1 (OVERLAP BELOW THR): Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		else if(overlapFlag==eCONT_OVERLAPPING){
			if(overlapAreaFraction_1<matchOverlapThr || overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT OVERLAP BELOW THR: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("CONTOUR OVERLAP MATCH: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
		#endif

	}//close if

	//## Check flux rel difference (if requested)	
	if(matchSourcesByFlux){
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}
	

	return true;

}//close HaveSourceIslandMatch()


bool HaveSourceComponentMatch(SourceComponentMatchPars& scompmatchpars,ComponentPars* pars1,ComponentPars* pars2)
{
	//################################################	
	//##        METHOD 
	//################################################
	//1) Check overlap first.
	//      - No overlap --> no match	
	//      - Complete overlap (e.g. one inside the other) --> check pos match	
	//          - Small pos match (e.g. within high pos thr) --> match
	//      - Small overlap (e.g. above low overlap thr) --> check pos match
	//           - High pos match (e.g. within low pos thr) --> match 	
	//      - Large overlap (e.g. above high overlap thr) --> check pos match	
	//           - Small pos match (e.g. within high pos thr) --> match
	//################################################

	//## Check data
	if(!pars1 || !pars2){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given for component pars!");
		#endif
		return false;
	}

	//## If sources are from the same collection return NO MATCH
	if(pars1->catalogIndex==pars2->catalogIndex) return false;

	//## Check fit ellipses
	TEllipse* ellipse1= pars1->fitEllipse;
	TEllipse* ellipse2= pars2->fitEllipse;
	if(!ellipse1 || !ellipse2){		
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/both fit ellipses are nullptr!");
		#endif
		return false;
	}

	//## Check ellipse centroid sky distance (if requested)
	double posThr= compMatchPosThr/3600.;//convert in deg as ellipse coords are in deg
	double posHighThr= compMatchPosHighThr/3600.;//convert in deg as ellipse coords are in deg
	double Xc_1= ellipse1->GetX1();
	double Yc_1= ellipse1->GetY1();
	double R1_1= ellipse1->GetR1();	
	double R2_1= ellipse1->GetR2();
	double Theta_1= ellipse1->GetTheta();
	double Xc_2= ellipse2->GetX1();
	double Yc_2= ellipse2->GetY1();
	double R1_2= ellipse2->GetR1();	
	double R2_2= ellipse2->GetR2();
	double Theta_2= ellipse2->GetTheta();

	//## Compute ellipse area, distance and overlap
	double posDist_Euclidean= MathUtils::GetEuclideanDist(Xc_1,Yc_1,Xc_2,Yc_2);
	double posDist_Haversine= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);
	double posDist= 0.;
	if(useWCSSimpleGeometry) posDist= posDist_Euclidean;
	else posDist= posDist_Haversine;

	double ellipseArea_1= MathUtils::ComputeEllipseArea(ellipse1);
	double ellipseArea_2= MathUtils::ComputeEllipseArea(ellipse2);
	double ellipseOverlapArea= -1;
	double err= 0;
	int rtn= 0;
	if(MathUtils::ComputeEllipseOverlapArea(ellipseOverlapArea,err,rtn,ellipse1,ellipse2)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute ellipse overlap area (err status="<<rtn<<"), return no match!");
		#endif
		return false;
	}	
	double overlapAreaFraction_1= ellipseOverlapArea/ellipseArea_1;
	double overlapAreaFraction_2= ellipseOverlapArea/ellipseArea_2;

	//## Check overlap
	bool hasPosLargeMatch= (fabs(posDist)<=posThr);
	bool hasPosSmallMatch= (fabs(posDist)<=posHighThr);
	bool noOverlap= (rtn==EllipseUtils::DISJOINT_ELLIPSES);
	bool hasOverlapCompleteMatch= (
		rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 || rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1
	);
	bool hasOverlapLargeMatch= (
		(!noOverlap && !hasOverlapCompleteMatch && overlapAreaFraction_1>=compMatchOverlapThr && overlapAreaFraction_2>=compMatchOverlapThr)
	);
	bool hasOverlapSmallMatch= (
		(!noOverlap && !hasOverlapCompleteMatch && overlapAreaFraction_1>=compMatchOverlapLowThr && overlapAreaFraction_2>=compMatchOverlapLowThr)
	);
	
	//## Compute ellipse flux diff
	double flux_1= pars1->fluxDensity;
	double flux_2= pars2->fluxDensity;
	double fluxRelDiff= (flux_1-flux_2)/flux_2;

	//## Fill match pars
	scompmatchpars.posDist_Euclidean= posDist_Euclidean;
	scompmatchpars.posDist_Haversine= posDist_Haversine;
	scompmatchpars.ellipseOverlapArea= ellipseOverlapArea;
	scompmatchpars.ellipseArea1= ellipseArea_1;
	scompmatchpars.ellipseArea2= ellipseArea_2;
	scompmatchpars.ellipseOverlapFlag= rtn;
	scompmatchpars.ellipseOverlapAreaRatio1= overlapAreaFraction_1;
	scompmatchpars.ellipseOverlapAreaRatio2= overlapAreaFraction_2;
	scompmatchpars.fluxRelDiff= fluxRelDiff;


	//Check overlap
	if(matchSourceComponentsByOverlap)
	{
		//- No overlap: Disjoint ellipses
		if(rtn==EllipseUtils::DISJOINT_ELLIPSES){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("DISJOINT ELLIPSES: Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
			#endif
			return false;
		}

		//- Full overlap: Ellipse 1 inside ellipse 2 (overlap=Area1)
		else if(rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 ){
			//Check area overlap ratio
			if(applySourceComponentAreaRatioThr && overlapAreaFraction_2<compMatchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE1 INSIDE ELLIPSE2 (OVERLAP BELOW THR): overlapArea2="<<overlapAreaFraction_2<<"<"<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			//Check small pos match
			if(matchSourceComponentsByPos && !hasPosSmallMatch){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}
		}

		//- Full overlap: Ellipse 2 inside ellipse 1 (overlap=Area2)
		else if(rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1){
			//Check area overlap ratio
			if(applySourceComponentAreaRatioThr && overlapAreaFraction_1<compMatchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE2 INSIDE ELLIPSE1 (OVERLAP BELOW THR): overlapArea1="<<overlapAreaFraction_1<<"<"<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			//Check small pos match
			if(matchSourceComponentsByPos && !hasPosSmallMatch){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}
		}

		//- Partial overlap
		else {
			//- Large overlap
			if(hasOverlapLargeMatch){	
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE LARGE OVERLAP: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				
				//Check small pos match
				if(matchSourceComponentsByPos && !hasPosSmallMatch){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
					#endif
					return false;
				}
			}//close if

			//- Small overlap
			else if(hasOverlapSmallMatch){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE SMALL OVERLAP: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapLowThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapLowThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif

				//Check large pos match
				if(matchSourceComponentsByPos && !hasPosLargeMatch){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("NO ELLIPSE POS LARGE MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
					#endif
					return false;
				}
			}//close else if
			
			//- Partial overlap not sufficient
			else {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE OVERLAP BELOW SMALL & HIGH THR: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			#ifdef LOGGING_ENABLED
				DEBUG_LOG("ELLIPSE OVERLAP MATCH: overlapArea1="<<overlapAreaFraction_1<<", overlapArea2="<<overlapAreaFraction_2<<", compMatchOverlapThr="<<compMatchOverlapThr<<", compMatchOverlapLowThr="<<compMatchOverlapLowThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
			#endif
		}//close else

	}//close if check overlap
	else{
		//## Check ellipse centroid sky distance (if requested)
		if(matchSourceComponentsByPos && !hasPosLargeMatch){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("NO ELLIPSE POS MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}");
			#endif
			return false;
		}

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("ELLIPSE POS MATCH: posDist(arcsec)="<<posDist*3600<<"<"<<posThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}");
		#endif
	}


	//## Check flux rel difference (if requested)	
	if(matchSourcesByFlux){
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("SOURCE MATCH: Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, overlapArea1="<<overlapAreaFraction_1<<", overlapArea2="<<overlapAreaFraction_2<<", posDist(arcsec)="<<posDist*3600<<", ");
	#endif

	return true;

}//close HaveSourceComponentMatch()


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

	//Compute WCS?
	WCS* wcs= 0;
	//int wcsType= 0;
	int wcsType= eJ2000;
	double sourceXmin= 1.e+99;
	double sourceXmax= -1.e+99;
	double sourceYmin= 1.e+99;
	double sourceYmax= -1.e+99;
	std::vector<SourcePars*> pars;
	int catalogIndex= -1;
	int sourceIndex= 0;
	
	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Found #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif

	for(int i=0;i<sourceTree->GetEntries();i++)
	{
		sourceTree->GetEntry(i);
		int type= aSource->Type;
		int simType= aSource->SimType;
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif


		//Select source by type?
		if(selectSourceByType){
			bool skipSource= true;
			for(size_t j=0;j<stypes.size();j++){
				if( stypes[j]==-1 || type==stypes[j]) {
					skipSource= false;
					break;
				}
			}
			if(skipSource) continue;
		}

		//Select source by sim type?
		if(selectSourceBySimType){
			bool skipSource= true;
			for(size_t j=0;j<ssimtypes.size();j++){
				if( ssimtypes[j]==-1 || type==ssimtypes[j]) {
					skipSource= false;
					break;
				}
			}
			if(skipSource) continue;
		}

		//Copy source
		Source* source= new Source;
		*source= *aSource;

		//Compute wcs for this source collection if not done
		if(!wcs){
			wcs= source->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to compute WCS from source no. "<<i+1<<" (name="<<source->GetName()<<")!");
				#endif
				return -1;
			}
		}

		//Shuffle source?
		
		double offsetX= 0;
		double offsetY= 0;
		/*
		if(shuffleSources){
			double A = 2/(shuffleRMax*shuffleRMax - shuffleRMin*shuffleRMin);
			double r= sqrt(2*gRandom->Uniform()/A + shuffleRMin*shuffleRMin);
			double theta= gRandom->Uniform(0,2*TMath::Pi());
			offsetX= r*cos(theta);
			offsetY= r*sin(theta);
			offsetX/= 3600;//convert in deg
			offsetY/= 3600;//convert in deg
		}
		*/

		//Compute source position in WCS coords (needed to compute source min/max coords)
		double X0_wcs= 0;
		double Y0_wcs= 0;
		source->GetWCSPos(X0_wcs,Y0_wcs,wcs,wcsType);

		if(shuffleSources){
			X0_wcs+= offsetX;
			Y0_wcs+= offsetY;
		}

		//Find min & max coordinates
		if(X0_wcs>sourceXmax) sourceXmax= X0_wcs;
		if(X0_wcs<sourceXmin) sourceXmin= X0_wcs;
		if(Y0_wcs>sourceYmax) sourceYmax= Y0_wcs;
		if(Y0_wcs<sourceYmin) sourceYmin= Y0_wcs;	
		
		//Fill source pars
		std::vector<SourcePars*> thisPars;
		//if(FillSourcePars(thisPars,source,catalogIndex,sourceIndex,-1,wcs,wcsType,shuffleSources,offsetX,offsetY)<0){
		if(FillSourcePars(thisPars,source,catalogIndex,sourceIndex,-1,wcs,wcsType,false,offsetX,offsetY)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill pars for source no. "<<i+1<<" (name="<<source->GetName()<<", collectionIndex="<<catalogIndex<<")!");
			#endif
			return -1;
		}

		//Add pars to pars
		pars.insert(pars.end(),thisPars.begin(),thisPars.end());

		//Add sources to list
		m_sources.push_back(source);
		sourceIndex++;

	}//end loop sources

	//Append source pars to main collection
	m_source_pars.insert(m_source_pars.end(),pars.begin(),pars.end());

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source min/max WCS coordinate range: (Xmin,Xmax)=("<<sourceXmin<<","<<sourceXmax<<"), (Ymin,Ymax)=("<<sourceYmin<<","<<sourceYmax<<")");
		INFO_LOG("#"<<m_sources.size()<<" sources to be cross-matched (#"<<m_source_pars.size()<<" source pars added) ...");
	#endif

	//Delete WCS for this collection
	WCSUtils::DeleteWCS(&wcs);

	//Define 2D grid with available sources	
	double borderTol= 0.1;
	double dX= fabs(sourceXmax-sourceXmin);
	double dY= fabs(sourceYmax-sourceYmin);
	tileXmin= sourceXmin - borderTol*dX;
	tileXmax= sourceXmax + borderTol*dX;
	tileYmin= sourceYmin - borderTol*dY;
	tileYmax= sourceYmax + borderTol*dY;
	tileXstep= (tileXmax-tileXmin)/nTilesX;
	tileYstep= (tileYmax-tileYmin)/nTilesY;
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing 2D "<<nTilesX<<" x "<<nTilesY<<" grid with source coordinates: (Xmin,Xmax,Xstep)=("<<tileXmin<<","<<tileXmax<<","<<tileXstep<<"), (Ymin,Ymax,Ystep)=("<<tileYmin<<","<<tileYmax<<","<<tileYstep<<")");
	#endif
	
	if(InitGrid(tileXmin,tileXmax,tileXstep,tileYmin,tileYmax,tileYstep)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to initialize tile data grid!");
		#endif
		return -1;
	}


	//## Add source data to tile
	for(size_t i=0;i<m_source_pars.size();i++)
	{
		double X0_wcs= m_source_pars[i]->X0_wcs;
		double Y0_wcs= m_source_pars[i]->Y0_wcs;
		std::string sourceName= m_source_pars[i]->sname;

		long int tileIndex= MathUtils::FindGrid2DBin(
			X0_wcs,Y0_wcs,
			nTilesX,tileXmin,tileXmax,tileXstep,
			nTilesY,tileYmin,tileYmax,tileYstep
		);

		if(tileIndex>=0 && tileIndex<(long int)(tileDataList.size())){//Add to tile data
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding source (name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<") to list...");
			#endif
			tileDataList[tileIndex]->AddSourcePars(m_source_pars[i]);
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot find tile index or invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<") for source (name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<"), check tile index calculation!");
			#endif
			continue;
		}
	}//end loop pars

	//Printing num sources in each tile
	for(size_t i=0;i<tileDataList.size();i++){
		long int nSourcePars= tileDataList[i]->GetNSourcePars();
		#ifdef LOGGING_ENABLED
			INFO_LOG("Tile no. "<<i+1<<": #"<<nSourcePars<<" sources to be matched...");
		#endif
	}

	return 0;

}//close ReadData()


int ReadCatalogData(std::string filelist)
{
	//Check filenames
	std::ifstream fileStream(filelist);	
	std::string line;
	if (fileStream.fail() || !fileStream.is_open()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filelist<<" for reading...");
		#endif
		return -1;
	}

	//Store filenames present in lists
	std::vector<std::string> fileNames;
	std::string filename= "";
	while (std::getline(fileStream, line)) {
  	std::istringstream iss(line);
    if (!(iss >> filename)) { 
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read line from file "<<filelist<<"!");
			#endif
			return -1; 
		}
    fileNames.push_back(filename);
	}//end file read
				
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<fileNames.size()<<" files present...");
	#endif

	//Initialize tile source catalog pars
	nCatalogs= static_cast<int>(fileNames.size());
	for(size_t i=0;i<tileDataList.size();i++){
		tileDataList[i]->InitCatalogSourcePars(nCatalogs);
	}
	for(int i=0;i<nCatalogs;i++){
		m_catalog_sources.push_back(std::vector<Source*>());
		m_catalog_source_pars.push_back(std::vector<SourcePars*>());
	}

	//Check nCatalogs is equal to region boundary size (if given as program input)
	if(!m_regions.empty()){
		int nCatalogRegions= static_cast<int>(m_regions.size());
		if(nCatalogRegions!=nCatalogs){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Mismatch between number of boundary regions ("<<nCatalogRegions<<") and number of source catalogs ("<<nCatalogs<<")!");
			#endif
			return -1;
		}
	}

	//Finally reading source data
	for(size_t i=0;i<fileNames.size();i++){
		if(ReadSourceData(fileNames[i],i)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read source data for file no. "<<i+1<<"!");
			#endif
			return -1;
		}
	}//end loop files

	//Printing num sources in each tile
	for(size_t i=0;i<tileDataList.size();i++){
		long int nSourceCatalogs= tileDataList[i]->GetNSourceCatalogs();
		std::stringstream ss;
		ss<<"Tile no. "<<i+1<<": {";
		for(long int j=0;j<nSourceCatalogs-1;j++){
			ss<<"#"<<tileDataList[i]->GetNSourceParsInCatalog(j)<<",";
		}
		ss<<"#"<<tileDataList[i]->GetNSourceParsInCatalog(nSourceCatalogs-1)<<"} catalog sources ...";
		#ifdef LOGGING_ENABLED
			INFO_LOG(ss.str());
		#endif
	}

	return 0;

}//close ReadCatalogData()


int ReadSourceData(std::string filename,int catalogIndex)
{
	//Open files
	TFile* f= new TFile(filename.c_str(),"READ");
	if(!f){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filename<<"!");
		#endif
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)f->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get access to source tree in file "<<filename<<"!");	
		#endif
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	//Initialize WCS & source pars
	WCS* wcs= 0;
	//int wcsType= 0;
	int wcsType= eJ2000;
	int sourceIndex= 0;

	//Find region boundary vertex
	double regionXmin= 1.e+99;
	double regionXmax= -1.e+99;
	double regionYmin= 1.e+99;
	double regionYmax= -1.e+99;

	if(catalogIndex!=-1 && readBoundaryRegions){
		std::vector<TVector2> bb= m_regionContours[catalogIndex]->BoundingBoxVertex;
		for(size_t k=0;k<bb.size();k++){
			double x= bb[k].X();
			double y= bb[k].Y();
			if(x>=regionXmax) regionXmax= x;
			if(x<=regionXmin) regionXmin= x;
			if(y>=regionYmax) regionYmax= y;
			if(y<=regionYmin) regionYmin= y;
		}
		if(shuffleSources){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source shuffle in region with boundary x["<<regionXmin<<","<<regionXmax<<"], y["<<regionYmin<<","<<regionYmax<<"] ...");
			#endif
		}
	}

	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);
		int type= aSource->Type;
		int simType= aSource->SimType;

		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source "<<i+1<<"/"<<sourceTree->GetEntries()<<" (name="<<aSource->GetName()<<") ...");
		#endif


		//Select source by type?
		if(selectSourceByType){
			bool skipSource= true;
			for(size_t j=0;j<stypes.size();j++){
				if( stypes[j]==-1 || type==stypes[j]) {
					skipSource= false;
					break;
				}
			}
			if(skipSource) continue;
		}

		//Select source by sim type?
		if(selectSourceBySimType){
			bool skipSource= true;
			for(size_t j=0;j<ssimtypes.size();j++){
				if( ssimtypes[j]==-1 || type==ssimtypes[j]) {
					skipSource= false;
					break;
				}
			}
			if(skipSource) continue;
		}

		
		//Copy source
		Source* source= new Source;
		*source= *aSource;


		//Compute wcs for this source collection if not done
		if(!wcs){
			wcs= source->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to compute WCS from source no. "<<i+1<<" (name="<<source->GetName()<<")!");
				#endif
				return -1;
			}
		}

		//Compute wcs source pos 
		double x0= source->X0;
		double y0= source->Y0;
		double x0_wcs= 0;
		double y0_wcs= 0;
		AstroUtils::PixelToWCSCoords(x0_wcs,y0_wcs,wcs,x0,y0);

		//Shuffle source?
		//NB: Generate uniform in annulus (e.g. see https://stackoverflow.com/questions/9048095/create-random-number-within-an-annulus)
		double offsetX= 0;
		double offsetY= 0;
		if(shuffleSources){
			if(readBoundaryRegions){
				double x0_wcs_rand= 0;
				double y0_wcs_rand= 0;
				if(GetRandPosInRegion(x0_wcs_rand,y0_wcs_rand,m_regionContours[catalogIndex],regionXmin,regionXmax,regionYmin,regionYmax)<0){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Failed to compute rand pos in region for source no. "<<sourceIndex<<"!");
					#endif
					return -1;
				}
				offsetX= x0_wcs_rand - x0_wcs;
				offsetY= y0_wcs_rand - y0_wcs;
				#ifdef LOGGING_ENABLED
					INFO_LOG("Shuffle source no. "<<i+1<<" (name="<<source->GetName()<<"), pos("<<x0_wcs<<","<<y0_wcs<<"), shuffle pos("<<x0_wcs_rand<<","<<y0_wcs_rand<<"), offset("<<offsetX<<","<<offsetY<<"), offset/arcsec("<<offsetX*3600.<<","<<offsetY*3600<<")");
				#endif
			}
			else{
				GetRandPosInAnnulus(offsetX,offsetY);
				offsetX/= 3600;//convert to deg
				offsetY/= 3600;//convert to deg
			}	

		}//close if shuffleSources

		/*
		if(shuffleSources){
			double A = 2/(shuffleRMax*shuffleRMax - shuffleRMin*shuffleRMin);
			double r= sqrt(2*gRandom->Uniform()/A + shuffleRMin*shuffleRMin);
			double theta= gRandom->Uniform(0,2*TMath::Pi());
			offsetX= r*cos(theta);
			offsetY= r*sin(theta);
			offsetX/= 3600;//convert in deg
			offsetY/= 3600;//convert in deg
		}
		*/

		//Fill source pars
		if(FillSourcePars(m_catalog_source_pars[catalogIndex],source,catalogIndex,sourceIndex,-1,wcs,wcsType,shuffleSources,offsetX,offsetY)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill pars for source no. "<<i+1<<" (name="<<source->GetName()<<", collectioIndex="<<catalogIndex<<")!");
			#endif
			return -1;
		}


		//Add source in collection	
		m_catalog_sources[catalogIndex].push_back(source);

		sourceIndex++;
	}//end loop sources

	

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<m_catalog_sources[catalogIndex].size()<<" sources to be cross-matched in catalog no. "<<catalogIndex<<" (#"<<m_catalog_source_pars[catalogIndex].size()<<" source pars added) ...");
	#endif

	
	//## Add source catalog data to tile
	#ifdef LOGGING_ENABLED
		INFO_LOG("Adding #"<<m_catalog_source_pars[catalogIndex].size()<<" source pars from catalog "<<catalogIndex<<" to list ...");
	#endif
	for(size_t i=0;i<m_catalog_source_pars[catalogIndex].size();i++)
	{
		if(!m_catalog_source_pars[catalogIndex][i]) cerr<<"null ptr!"<<endl;
		double X0_wcs= m_catalog_source_pars[catalogIndex][i]->X0_wcs;
		double Y0_wcs= m_catalog_source_pars[catalogIndex][i]->Y0_wcs;
		std::string sourceName= m_catalog_source_pars[catalogIndex][i]->sname;
		
		long int tileIndex= MathUtils::FindGrid2DBin(
			X0_wcs,Y0_wcs,
			nTilesX,tileXmin,tileXmax,tileXstep,
			nTilesY,tileYmin,tileYmax,tileYstep
		);

		if(tileIndex>=0 && tileIndex<(long int)(tileDataList.size())){//Add to tile data
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding source (catalogId="<<catalogIndex<<", name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<") to list...");
			#endif
			int status= tileDataList[tileIndex]->AddCatalogSourcePars(catalogIndex,m_catalog_source_pars[catalogIndex][i]);
			if(status<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to add source par no. "<<i+1<<" from catalog "<<catalogIndex<<" to tile!");
				#endif
				continue;
			}
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot find tile index or invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<") for source (catalogId="<<catalogIndex<<", name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<"), check tile index calculation!");
			#endif
			continue;
		}

	}//end loop pars

	//Delete WCS for this collection
	WCSUtils::DeleteWCS(&wcs);

	return 0;

}//close ReadSourceData()

int FillSourcePars(std::vector<SourcePars*>& pars,Source* aSource,int catalogIndex,int sourceIndex,int nestedSourceIndex,WCS* wcs,int coordSystem,bool shuffleData,double offsetX,double offsetY)
{
	//####  METHOD ##############################
	//Distinguish different source cases
	//1) Single-source (point/compact) with fit info (1 component)
	//     1.1) No nested sources
	//     1.2) 1 nested source: ignore in collection
	//2) Multi-component source (point/compact) with fit info (N component)
	//     2.1) No nested source
	//     2.2) N nested sources: ignore in collection
	//3) Multi-component source (extended) WITHOUT fit info
	//     3.1) No nested source --> skip it
	//     3.2) Nested source without fit info --> skip it
	//     3.3) Nested source with fit info --> add to collection with its components
	//###########################################

	//Check input source
	if(!aSource) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Nullptr to input source given!");
		#endif
		return -1;	
	}

	//Get source data
	std::string sourceName= std::string(aSource->GetName());
	bool hasFitInfo= aSource->HasFitInfo();
	bool hasNestedSources= aSource->HasNestedSources();
	double beamArea= aSource->GetBeamFluxIntegral();
	double S= aSource->GetS();	
	if(correctFlux) S/= beamArea;

	//Get 0th contour converted to WCS & relative centroid
	bool useFWHM= true;
	bool convertToWCS= true;
	int pixOffset= 0;
	bool computeContourPars= true;
	Contour* contour= aSource->GetWCSContour(0,wcs,coordSystem,pixOffset,computeContourPars);
	if(!contour){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute WCS contour for source "<<sourceName<<"!");
		#endif
		return -1;
	}

	//Shuffle contour using offset
	if(shuffleData){
		contour->ApplyOffset(offsetX,offsetY);
	}
	
	TVector2 contourCentroid= contour->Centroid;

	
	//## Fill source pars for island first
	SourcePars* sourcePars= new SourcePars(catalogIndex,sourceIndex,nestedSourceIndex);
	sourcePars->sname= sourceName;
	sourcePars->contour= contour;
	sourcePars->S= S;
	sourcePars->X0_wcs= contourCentroid.X();
	sourcePars->Y0_wcs= contourCentroid.Y();
	

	//## Fill source fit components
	if(hasFitInfo)
	{
		//Get fitted pars & ellipse converted to WCS
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= fitPars.GetNComponents();
		std::vector<TEllipse*> fitEllipses;
		if(aSource->GetFitEllipses(fitEllipses,useFWHM,convertToWCS,wcs,coordSystem,pixOffset,useWCSSimpleGeometry)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute WCS ellipse for fitted components of source "<<sourceName<<"!");
			#endif
			CodeUtils::DeletePtr<Contour>(contour);
			CodeUtils::DeletePtr<SourcePars>(sourcePars);
			return -1;
		}

		//Check ellipses and pars size
		if(nComponents!=(int)(fitEllipses.size())){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Number of fit components shall be equal to fitted ellipses (this should not occur)!");
			#endif
			CodeUtils::DeletePtr<Contour>(contour);
			CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
			CodeUtils::DeletePtr<SourcePars>(sourcePars);
			return -1;
		}


		//Set fitted flux density
		double fluxDensity= 0;
		aSource->GetFluxDensity(fluxDensity);
		if(correctFlux) fluxDensity/= beamArea;
		sourcePars->fluxDensity= fluxDensity;


		//Fill component pars & ellipse
		ComponentPars* thisComponentPars= 0;
		for(int k=0;k<nComponents;k++){
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));
			double A= fitPars.GetParValue(k,"A");
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			if(correctFlux) fluxDensity/= beamArea;

			//Apply offset to ellipse (if shuffle source is enabled)
			if(shuffleData){
				double Cx= fitEllipses[k]->GetX1();
				double Cy= fitEllipses[k]->GetY1();
				fitEllipses[k]->SetX1(Cx + offsetX);
				fitEllipses[k]->SetY1(Cy + offsetY);
			}

			//Fill component par
			thisComponentPars= new ComponentPars;
			thisComponentPars->sname= sname;
			thisComponentPars->catalogIndex= catalogIndex;	
			thisComponentPars->fitComponentIndex= k;
			thisComponentPars->fitEllipse= fitEllipses[k];
			thisComponentPars->A= A;
			thisComponentPars->fluxDensity= fluxDensity;
			thisComponentPars->isSelected= isSelected;

			//Add component to source pars
			sourcePars->AddComponentPars(thisComponentPars);	
		}//end loop components

	}//close if has fit info

	
	//Append this source pars to collection
	pars.push_back(sourcePars);

	//## Fill pars for nested sources (if any)
	if(hasNestedSources){
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			if(FillSourcePars(pars,nestedSources[j],catalogIndex,sourceIndex,j,wcs,coordSystem,shuffleData,offsetX,offsetY)<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
				#endif
				return -1;
			}
		}//end loop sources
	}//close if nested sources
	
	return 0;

}//close FillSourcePars()


int ReadRegionData(std::string filelist)
{
	//Check filenames
	std::ifstream fileStream(filelist);	
	std::string line;
	if (fileStream.fail() || !fileStream.is_open()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filelist<<" for reading...");
		#endif
		return -1;
	}

	//Store filenames present in lists
	std::vector<std::string> fileNames;
	std::string filename= "";
	while (std::getline(fileStream, line)) {
  	std::istringstream iss(line);
    if (!(iss >> filename)) { 
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read line from file "<<filelist<<"!");
			#endif
			return -1; 
		}
    fileNames.push_back(filename);
	}//end file read
				
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<fileNames.size()<<" files present...");
	#endif


	//Initialize region boundary 
	int nFiles= static_cast<int>(fileNames.size());
	for(int i=0;i<nFiles;i++){
		m_regions.push_back(std::vector<DS9Region*>());
	}

	//Finally reading region data
	for(size_t i=0;i<fileNames.size();i++){
		if(ReadCatalogRegionData(fileNames[i],i)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read region data for file no. "<<i+1<<"!");
			#endif
			return -1;
		}
	}//end loop files

	//Compute region contours and check that only one region is given. Multiple region not supported.
	for(size_t i=0;i<m_regions.size();i++)
	{
		if(m_regions[i].size()!=1){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("More than one region given for catalog file no. "<<i+1<<" (hint: multi-region not supported)!");
			#endif
			return -1;
		}
		
		//Compute contour
		Contour* contour= m_regions[i][0]->GetContour(true);
		if(!contour){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to get contour for catalog region "<<i+1<<"!");
			#endif
			return -1;
		}
		m_regionContours.push_back(contour);

	}//end loop catalog regions

	return 0;

}//close ReadRegionData()

int ReadCatalogRegionData(std::string filename,int catalogIndex)
{
	//Read and parse DS9 regions
	if (DS9RegionParser::Parse(m_regions[catalogIndex],filename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse DS9 region file "<<filename<<"!");
		#endif
		return -1;
	}

	//Check if empty
	if(m_regions[catalogIndex].empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No regions read from file "<<filename<<"!");
		#endif
		return -1;
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<m_regions[catalogIndex].size()<<" regions read from file "<<filename<<"...");
		#endif
	}

	return 0;

}//close ReadCatalogRegionData()

int GetRandPosInAnnulus(double& offsetX,double& offsetY)
{
	//NB: Generate uniform in annulus (e.g. see https://stackoverflow.com/questions/9048095/create-random-number-within-an-annulus)
	//    Offset returned in arcsec
	double A = 2/(shuffleRMax*shuffleRMax - shuffleRMin*shuffleRMin);
	double r= sqrt(2*gRandom->Uniform()/A + shuffleRMin*shuffleRMin);
	double theta= gRandom->Uniform(0,2*TMath::Pi());
	offsetX= r*cos(theta);
	offsetY= r*sin(theta);
	
	return 0;

}//close GetRandPosInAnnulus()


int GetRandPosInRegion(double& x0,double& y0,Contour* contour,double xmin,double xmax,double ymin,double ymax,int nMaxGen)
{
	//NB: Position returned in deg
	int nGen= 0;
	x0= 0;
	y0= 0;
	bool isGoodGen= false;

	while(nGen<nMaxGen) 
	{
		double xrand= gRandom->Uniform(xmin,xmax);				
		double yrand= gRandom->Uniform(ymin,ymax);				
		bool isInsideContour= contour->IsPointInsideContour(xrand,yrand);
		if(isInsideContour){
			x0= xrand;	
			y0= yrand;		
			isGoodGen= true;
			break;
		}
		nGen++;

	}//end while

	if(!isGoodGen){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to generate pos in region (nGen="<<nGen<<">="<<nMaxGen<<")!");
		#endif
		return -1;
	}

	return 0;

}//close GetRandPosInRegion()

void Save()
{
	//Save match TTree to file
	if(matchOutputFile && matchOutputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving #"<<sourceMatchDataCollection.size()<<" source match data to file...");
		#endif
		matchOutputFile->cd();
		if(matchedSourceInfo){		
			for(size_t i=0;i<sourceMatchDataCollection.size();i++){
				sourceMatchData= sourceMatchDataCollection[i];
				matchedSourceInfo->Fill();
			}
			matchedSourceInfo->Write();
		}
		if(matchOptionTree){
			matchOptionTree->Write();
		}

		//Write SED plots to file
		//SaveSED();

		matchOutputFile->Close();
	}

	//Save source TTree to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		
		//Clone all objects present in input file but the Source TTree
		CloneObjectsInFile({"SourceInfo"});

		//Write selected source TTree
		SaveSources();

		//Close file
		outputFile->Close();

	}//close if

}//close Save()

int SaveSED()
{
	//Save SED graphs
	for(size_t i=0;i<m_SEDGraphs.size();i++){
		if(m_SEDGraphs[i]) m_SEDGraphs[i]->Write();
	}

	//Save SED component graphs
	for(size_t i=0;i<m_SEDComponentGraphs.size();i++){
		if(m_SEDComponentGraphs[i]) m_SEDComponentGraphs[i]->Write();
	}

	//Save sed plots to file
	for(size_t i=0;i<m_SEDPlots.size();i++){
		if(m_SEDPlots[i]) m_SEDPlots[i]->Write();
	}

	return 0;

}//close SavePlots()

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
	for(size_t i=0;i<sourceMatchDataCollection.size();i++){
		m_source= sourceMatchDataCollection[i]->GetSource();
		outputTree->Fill();
	}

	outputTree->Write();

	return 0;

}//close SaveSources()

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

