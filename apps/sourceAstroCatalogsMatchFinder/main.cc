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
#include <AstroObject.h>
#include <AstroObjectParser.h>
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

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

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
	cout<<"-c, --catalogType=[CATALOG_TYPE] \t Catalog type id {1=SIMBAD, 2=NED, 3=MGPS, 4=HASH, 5=WISE-HII, 6=ATNF-PSR, 7=WOLF-RAYET, 8=SELAVY} (default=SIMBAD)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (root format) in which to store source match info"<<endl;
	cout<<"-n, --nx=[NX] \t Number of divisions along X (default=4)"<<endl;
	cout<<"-N, --ny=[NY] \t Number of divisions along Y (default=4)"<<endl;
	cout<<"-q, --selectObjectsInRegions \t Select catalog object inside given DS9 regions (default=no)"<<endl;
	cout<<"-Q, --regionFile \t DS9 region filename used to select catalog objects"<<endl;
	cout<<"-m, --matchSourcesByFlux \t Match sources by flux (default=false)"<<endl;
	cout<<"-M, --matchFluxRelDiffThr \t Flux relative difference threshold (default=0.5)"<<endl;
	cout<<"-p, --matchComponentByPos \t Match source fit components by position centroid distance (default=false)"<<endl;
	cout<<"-b, --matchComponentByOverlap \t Match source fit component by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-t, --compMatchPosThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	cout<<"-e, --compMatchPosHighThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=5)"<<endl;
	cout<<"-a, --compMatchOverlapThr \t Source fit component contour overlap fraction (wrt to total area) threshold (default=0.8)"<<endl;
	cout<<"-d, --compMatchOverlapLowThr \t Source fit component contour overlap fraction (wrt to total area) low threshold (default=0.2)"<<endl;
	cout<<"-l, --applyComponentAreaRatioThr \t If enabled and if source fit component is fully enclosed inside another (or viceversa) a match requires that the source1/source2 area ratio is higher than the component overlap threshold. If disabled, assume a match whenever a source component is fully enclosed in a source (or viceversa) (default=false)"<<endl;
	cout<<"-f, --filterByType \t Consider only sources with given type when searching the match (default=no)"<<endl;
	cout<<"-s, --selectedType=[TYPE] \t Source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED) (default=-1)"<<endl;
	cout<<"-F, --filterObjectByType \t Consider only object with given type when searching the match (default=no)"<<endl;
	cout<<"-S, --selectedObjectType=[TYPE] \t Object types to be crossmatched (default=-1)"<<endl;
	cout<<"-g, --shuffleSources \t Randomize sources and catalog sources in an annulus centred on their original position (default=no)"<<endl;
	cout<<"-r, --shuffleRmin=[Rmin] \t Minimum annulus radius in arcsec for source shuffling (default=30 arcsec)"<<endl;
	cout<<"-R, --shuffleRmax=[Rmax] \t Maximum annulus radius in arcsec for source shuffling (default=50 arcsec)"<<endl;
	
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "catalogs", required_argument, 0, 'C' },
	{ "catalogType", required_argument, 0, 'c' },	
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", required_argument, 0, 'o' },	
	{ "selectObjectsInRegions", no_argument, 0, 'q' },	
	{ "regionFile", required_argument, 0, 'Q' },	
	{ "nx", required_argument, 0, 'n' },
	{ "ny", required_argument, 0, 'N' },
	{ "matchSourcesByFlux", no_argument, 0, 'm'},
	{ "matchFluxRelDiffThr", required_argument, 0, 'M'},
	{ "matchComponentByPos", no_argument, 0, 'p'},
	{ "matchComponentByOverlap", no_argument, 0, 'b'},	
	{ "compMatchPosThr", required_argument, 0, 't'},
	{ "compMatchPosHighThr", required_argument, 0, 'e'},
	{ "compMatchOverlapThr", required_argument, 0, 'a'},
	{ "compMatchOverlapLowThr", required_argument, 0, 'd'},
	{ "applyComponentAreaRatioThr", no_argument, 0, 'l'},	
	{ "filterByType", no_argument, 0, 'f'},
	{ "selectedType", required_argument, 0, 's'},
	{ "filterObjectByType", no_argument, 0, 'F'},
	{ "selectedObjectType", required_argument, 0, 'S'},
	{ "shuffleSources", no_argument, 0, 'g'},	
	{ "shuffleRmin", required_argument, 0, 'r'},
	{ "shuffleRmax", required_argument, 0, 'R'},
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string catalogFileName= "";
bool selectObjectsInRegions= false;
std::string regionFileName= "";
int verbosity= 4;//INFO level
bool useWCSSimpleGeometry= true;
bool matchSourcesByFlux= false;
float matchFluxRelDiffThr= 0.5;//50%
bool matchSourceComponentsByPos= false;
float compMatchPosThr= 2.5;//dist in arcsec below which source fit component and region can be matched
float compMatchPosHighThr= 5;//dist in arcsec below which source fit component and region can be matched
bool matchSourceComponentsByOverlap= false;
float compMatchOverlapThr= 0.8;//fraction of overlap above which source fit component and region are matched
float compMatchOverlapLowThr= 0.2;//fraction of overlap above which source fit component and region are matched
bool applySourceComponentAreaRatioThr= false;
bool correctFlux= false;
bool selectSourceByType= false;//default=all true sources searched 
std::vector<int> stypes;
bool selectObjectByType= false;//default=all object searched 
std::vector<int> otypes;
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

enum CatalogType 
{
	eSIMBAD_CATALOG=1,
	eNED_CATALOG=2,
	eMGPS_CATALOG=3,
	eHASH_CATALOG=4,
	eWISE_HII_CATALOG=5,
	eATNF_PSR_CATALOG=6,
	eWOLF_RAYET_CATALOG=7,
	eSELAVY_CATALOG=8,
};
int catalogType= 1;


//Globar vars
TFile* outputFile= 0;
std::string outputFileName= "MatchOutput.root";
TTree* matchedSourceInfo= 0;
TTree* matchOptionTree= 0;
std::string sourceName= "";
double sourcePosX= 0;
double sourcePosY= 0;
double sourcePosX_wcs= 0;
double sourcePosY_wcs= 0;
double sourceFluxDensity= 0;
double sourceFlux= 0;
double sourceFluxDensityErr= 0;
double sourceFluxErr= 0;
int nMatchedObjects= 0;
std::vector<double> fluxDiff;
std::vector<int> matchedObjectMultiplicity;
std::vector<AstroObject*> matchedObjects;
int nTotSourceComponents= 0;
int nMatchedSourceComponents= 0;
int nCatalogs= 0;
std::vector<int> nMatchedSourceComponentsPerCatalog;
std::vector<int> nMatchedSourceComponentCatalogMultiplicity;

//====================================
//     ComponentPars struct
//====================================
struct ComponentPars 
{
	ComponentPars(int _componentIndex,int _sourceIndex,int _nestedSourceIndex=-1,int _catalogIndex=-1)
		: catalogIndex(_catalogIndex), sourceIndex(_sourceIndex), nestedSourceIndex(_nestedSourceIndex), fitComponentIndex(_componentIndex)
	{
		sname= "";
		fitEllipse= nullptr;
		fluxDensity= 0;
		flux= 0;
		fluxDensityErr= 0;
		fluxErr= 0;
		X0= 0;
		Y0= 0;
		X0_wcs= 0;
		Y0_wcs= 0;
		isSelected= true;
	}
	~ComponentPars(){
		CodeUtils::DeletePtr<TEllipse>(fitEllipse);
	}
	
	int catalogIndex;
	int sourceIndex;
	int nestedSourceIndex;
	int fitComponentIndex;
	std::string sname;
	double X0_wcs;
	double Y0_wcs;
	double X0;
	double Y0;
	double fluxDensity;//fitted flux density
	double flux;//fitted flux
	double fluxDensityErr;//fitted flux density err
	double fluxErr;//fitted flux
	TEllipse* fitEllipse;
	bool isSelected;//is selected component
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
		sourceComponentPars.clear();
		catalogObjects.clear();
	}

	//Destructor 
	~TileData()
	{
		CodeUtils::DeletePtrCollection<ComponentPars>(sourceComponentPars);

		for(size_t i=0;i<catalogObjects.size();i++){
			for(size_t j=0;j<catalogObjects[i].size();j++){
				if(catalogObjects[i][j]){
					delete catalogObjects[i][j];
					catalogObjects[i][j]= 0;
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
	int AddSourceComponentPars(ComponentPars* pars){
		if(!pars) return -1;
		sourceComponentPars.push_back(pars);
		return 0;
	}

	//Get number of sources pars in this tile
	long int GetNSourceComponentPars(){
		return static_cast<long int>(sourceComponentPars.size());
	}

	//Init source catalog pars
	int InitCatalogObjects(int NCatalogs){
		if(NCatalogs<=0) return -1;
		catalogObjects.clear();
		for(int i=0;i<NCatalogs;i++) catalogObjects.push_back( std::vector<AstroObject*>() );
		return 0;
	}

	//Add source catalog pars
	int AddCatalogObject(int catalogIndex,AstroObject* obj){
		if(!obj) return -1;
		if(catalogIndex<0 || catalogIndex>=(int)(catalogObjects.size())) return -1;
		catalogObjects[catalogIndex].push_back(obj);
		return 0;
	}

	//Get number of sources pars in this tile
	long int GetNCatalogObjects(){
		return static_cast<long int>(catalogObjects.size());
	}

	//Get number of sources catalog pars in this tile
	long int GetNObjectsInCatalog(int catalogIndex){
		if(catalogIndex<0 || catalogIndex>=(int)(catalogObjects.size())) return 0;
		return static_cast<long int>(catalogObjects[catalogIndex].size());
	}

	//Get source catalog pars
	AstroObject* GetCatalogObject(int catalogIndex,int cindex)
	{
		if(catalogIndex<0 || catalogIndex>=(int)(catalogObjects.size())) return nullptr;
		if(cindex<0 || cindex>=(int)(catalogObjects[catalogIndex].size())) return nullptr;	
		return catalogObjects[catalogIndex][cindex];
	}

	//Tile coordinates
	float xmin;
	float xmax;
	float ymin;
	float ymax;

	//List of neighbors tiles
	std::vector<size_t> neighborTileIndexes;

	//Source component pars for this tile
	std::vector<ComponentPars*> sourceComponentPars;

	//Catalog objects for this tile
	std::vector<std::vector<AstroObject*>> catalogObjects;

};//close TileData()

std::vector<TileData*> tileDataList;

std::vector<Source*> m_sources;
std::vector<ComponentPars*> m_sourceComponentPars;
std::vector<std::vector<AstroObject*>> m_catalogObjects;
std::vector<DS9Region*> m_regions;



//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int InitGrid(float xmin,float xmax,float xstep,float ymin,float ymax,float ystep);
int FindSourceMatchesInTiles();
int ComputeSpectralIndices();
int ReadRegions();
int ReadData(std::string filename);
int ReadCatalogData(std::string filelist);
int ReadAstroObjectData(std::string filename,int collectionIndex);
int FillSourceComponentPars(std::vector<ComponentPars*>& pars,Source* aSource,int collectionIndex,int sourceIndex,int nestedSourceIndex=-1,WCS* wcs=0,int coordSystem=0,double offsetX=0,double offsetY=0);
bool HaveSourceComponentMatch(SourceComponentMatchPars& scompmatchpars,ComponentPars* cpars,AstroObject* obj);
void Save();
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
	//== Read regions
	//=======================
	if(selectObjectsInRegions){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading regions in file "<<regionFileName<<" ...");
		#endif
		if(ReadRegions()<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Region reading failed!");
			#endif
			return -1;
		}
	}//close if

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

	while((c = getopt_long(argc, argv, "hi:C:c:o:v:qQ:n:N:mM:Llt:e:pba:d:fs:FS:gr:R:",options_tab, &option_index)) != -1) {
    
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
			case 'c':	
			{
				catalogType= atoi(optarg);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'q':	
			{
				selectObjectsInRegions= true;
				break;	
			}
			case 'Q':	
			{
				regionFileName= std::string(optarg);	
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
				selectObjectByType= true;
				break;
			}
			case 'S':	
			{
				int type= atoi(optarg);
				otypes.push_back(type);	
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
	//== Check options
	//=======================
	if( catalogType!=eSIMBAD_CATALOG &&
			catalogType!=eNED_CATALOG &&
			catalogType!=eMGPS_CATALOG &&
			catalogType!=eHASH_CATALOG &&
			catalogType!=eWISE_HII_CATALOG &&
			catalogType!=eATNF_PSR_CATALOG &&
			catalogType!=eWOLF_RAYET_CATALOG &&
			catalogType!=eSELAVY_CATALOG
	){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unsupported catalog type given (type="<<catalogType<<")!");
		#endif
		return -1;
	}	

	//=======================
	//== Print options
	//=======================
	#ifdef LOGGING_ENABLED	
		INFO_LOG("catalogType="<<catalogType);
		INFO_LOG("matchSourcesByFlux? "<<matchSourcesByFlux<<", matchFluxRelDiffThr="<<matchFluxRelDiffThr);
		INFO_LOG("matchSourceComponentsByPos? "<<matchSourceComponentsByPos<<", compMatchPosThr(arcsec)="<<compMatchPosThr<<", compMatchPosHighThr(arcsec)="<<compMatchPosHighThr);
		INFO_LOG("matchSourceComponentsByOverlap? "<<matchSourceComponentsByOverlap<<", compMatchOverlapThr="<<compMatchOverlapThr<<", compMatchOverlapLowThr="<<compMatchOverlapLowThr);
		INFO_LOG("applySourceComponentAreaRatioThr? "<<applySourceComponentAreaRatioThr);
	#endif

	return 0;

}//close ParseOptions()

int Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//Init source cube TTree
	if(!matchedSourceInfo) matchedSourceInfo= new TTree("SourceMatchInfo","SourceMatchInfo");
	matchedSourceInfo->Branch("name",&sourceName);
	matchedSourceInfo->Branch("x",&sourcePosX);
	matchedSourceInfo->Branch("y",&sourcePosY);
	matchedSourceInfo->Branch("x_wcs",&sourcePosX_wcs);
	matchedSourceInfo->Branch("y_wcs",&sourcePosY_wcs);
	matchedSourceInfo->Branch("fluxDensity",&sourceFluxDensity);
	matchedSourceInfo->Branch("fluxDensityErr",&sourceFluxDensityErr);	
	matchedSourceInfo->Branch("flux",&sourceFlux);
	matchedSourceInfo->Branch("fluxErr",&sourceFluxErr);
	matchedSourceInfo->Branch("nMatchedObjects",&nMatchedObjects);
	matchedSourceInfo->Branch("matchedObjects",&matchedObjects);
	//matchedSourceInfo->Branch("fluxDiff",&fluxDiff);
	matchedSourceInfo->Branch("matchedObjectMultiplicity",&matchedObjectMultiplicity);
	

	//Init match option TTree
	if(!matchOptionTree) matchOptionTree= new TTree("MatchSummary","MatchSummary");
	matchOptionTree->Branch("matchSourceComponentsByPos",&matchSourceComponentsByPos);
	matchOptionTree->Branch("compMatchPosThr",&compMatchPosThr);
	matchOptionTree->Branch("compMatchPosHighThr",&compMatchPosHighThr);
	matchOptionTree->Branch("matchSourceComponentsByOverlap",&matchSourceComponentsByOverlap);
	matchOptionTree->Branch("compMatchOverlapThr",&compMatchOverlapThr);
	matchOptionTree->Branch("compMatchOverlapLowThr",&compMatchOverlapLowThr);
	matchOptionTree->Branch("selectSourceByType",&selectSourceByType);
	matchOptionTree->Branch("applySourceAreaRatioThr",&applySourceAreaRatioThr);
	matchOptionTree->Branch("applySourceComponentAreaRatioThr",&applySourceComponentAreaRatioThr);

	matchOptionTree->Branch("nTotSourceComponents",&nTotSourceComponents);
	matchOptionTree->Branch("nMatchedSourceComponents",&nMatchedSourceComponents);
	matchOptionTree->Branch("nCatalogs",&nCatalogs);
	matchOptionTree->Branch("nMatchedSourceComponentsPerCatalog",&nMatchedSourceComponentsPerCatalog);
	matchOptionTree->Branch("nMatchedSourceComponentCatalogMultiplicity",&nMatchedSourceComponentCatalogMultiplicity);

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





int FindSourceMatchesInTiles()
{
	//Init counters	
	nTotSourceComponents= 0;
	nMatchedSourceComponents= 0;
	std::map<std::string,int> matchedSourceNMatchMap;

	struct MatchInfo 
	{
		MatchInfo(){
			cpars= 0;
			objs.clear();
		}
		ComponentPars* cpars;
		std::vector<AstroObject*> objs;
	};
	std::vector<MatchInfo*> matchInfoCollection;
	MatchInfo* matchInfo= 0;

	//Loop over tiles and find matches with regions in the same tile and in neighbor tiles
	for(size_t i=0;i<tileDataList.size();i++)
	{
		long int NSourceComponentPars= tileDataList[i]->GetNSourceComponentPars();
		long int NCatalogs= tileDataList[i]->GetNCatalogObjects();
		std::vector<size_t> neighbors= tileDataList[i]->neighborTileIndexes;
		

		//Init counters
		if(nMatchedSourceComponentsPerCatalog.empty()){
			nMatchedSourceComponentsPerCatalog= std::vector<int>(NCatalogs,0);
		}
		if(nMatchedSourceComponentCatalogMultiplicity.empty()){
			nMatchedSourceComponentCatalogMultiplicity= std::vector<int>(NCatalogs+1,0);
		}

		for(long int j=0;j<NSourceComponentPars;j++) 
		{		
			//Get source component 
			ComponentPars* componentPars= (tileDataList[i]->sourceComponentPars)[j];
			
			#ifdef LOGGING_ENABLED
			if(j%100==0) INFO_LOG("#"<<j+1<<"/"<<NSourceComponentPars<<" source components to be matched in tile no. "<<i+1<<" against "<<NCatalogs<<" catalogs ...");
			#endif

			//Count total number of sources & components					
			nTotSourceComponents++;
			std::vector<int> matchCounter(NCatalogs,0);
			bool matchFound= false;
			nMatchedObjects= 0;

			//Create MatchInfo 
			matchInfo= new MatchInfo;
			matchInfo->cpars= componentPars;

			//==============================================
			//===        CATALOG LOOP
			//==============================================
			//Loop over source catalogs in the same tile
			for(long int s=0;s<NCatalogs;s++)
			{
				int nMatches= 0;	
									
				long int NObjectsInCatalog= tileDataList[i]->GetNObjectsInCatalog(s);
				if(NObjectsInCatalog<=0) continue;

				#ifdef LOGGING_ENABLED
					DEBUG_LOG("#"<<NSourceComponentPars<<" sources to be cross-matched against #"<<NObjectsInCatalog<<" objects (catalog #"<<s+1<<") in tile no. "<<i+1<<" ...");
				#endif

				//==============================================
				//===      CATALOG SOURCE LOOP (SAME TILE)
				//==============================================	
				for(long int k=0;k<NObjectsInCatalog;k++)
				{
					//Get catalog object
					AstroObject* catalogObject= tileDataList[i]->GetCatalogObject(s,k);
					if(!catalogObject){
						#ifdef LOGGING_ENABLED
							ERROR_LOG("Null ptr to catalog object (catalog no. "<<s+1<<", object no. "<<k+1<<")!");
						#endif
						return -1;
					}
					
					//Check match
					SourceComponentMatchPars cmatchpars;
					bool hasMatch= HaveSourceComponentMatch(cmatchpars,componentPars,catalogObject);
					if(!hasMatch) {
						continue;
					}
					matchFound= true;
					matchCounter[s]++;	
					matchedSourceNMatchMap[catalogObject->name]++;
					//matchedObjects.push_back(catalogObject);
					nMatchedObjects++;

					//double dflux= sourceFlux- catalogObject->flux;
					//fluxDiff.push_back(dflux);

					(matchInfo->objs).push_back(catalogObject);

				}//end (k index) loop objects in catalog


				//==============================================
				//===     CATALOG SOURCE LOOP (NEIGHBOUR TILES)
				//==============================================	
				for(size_t n=0;n<neighbors.size();n++)	
				{
					size_t neighborIndex= neighbors[n];
					if(i==neighborIndex) continue;//skip same tile

					long int NObjectsInCatalog_neighbor= tileDataList[neighborIndex]->GetNObjectsInCatalog(s);
					if(NObjectsInCatalog_neighbor<=0) continue;

					#ifdef LOGGING_ENABLED
						DEBUG_LOG("#"<<NSourceComponentPars<<" sources to be cross-matched against #"<<NObjectsInCatalog_neighbor<<" catalog objects (catalog #"<<s+1<<") in tile no. "<<neighborIndex<<" ...");
					#endif	
				
					for(long int k=0;k<NObjectsInCatalog_neighbor;k++)
					{
						//Get catalog object
						AstroObject* catalogObject= tileDataList[neighborIndex]->GetCatalogObject(s,k);

						if(!catalogObject){
							#ifdef LOGGING_ENABLED
								ERROR_LOG("Null ptr to catalog object (catalog no. "<<s+1<<", object no. "<<k+1<<"!");
							#endif
							return -1;
						}
						
						//Check source match
						SourceComponentMatchPars cmatchpars;
						bool hasMatch= HaveSourceComponentMatch(cmatchpars,componentPars,catalogObject);
						if(!hasMatch) {
							continue;
						}
						matchFound= true;
						matchCounter[s]++;
						matchedSourceNMatchMap[catalogObject->name]++;
						//matchedObjects.push_back(catalogObject);
						nMatchedObjects++;

						//double dflux= sourceFlux- catalogObject->flux;
						//fluxDiff.push_back(dflux);

						(matchInfo->objs).push_back(catalogObject);

					}//end (k index) loop objects in catalog
				}//end (n index) loop neighbor tiles
	
			}//end (s index) loop catalogs


			//Update counters
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Update source match counters...");
			#endif
			if(matchFound){
				nMatchedSourceComponents++;
			}
			std::vector<int> sourceComponentMatchMultiplicityCounter(NCatalogs,0);
			for(size_t c=0;c<matchCounter.size();c++){
				if(matchCounter[c]>0) {
					nMatchedSourceComponentsPerCatalog[c]++;
					sourceComponentMatchMultiplicityCounter[c]++;
				}
			}

			//Compute catalog match multiplicity
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Compute catalog match multiplicity...");
			#endif
			int multiplicity= std::accumulate(sourceComponentMatchMultiplicityCounter.begin(), sourceComponentMatchMultiplicityCounter.end(), 0);
			nMatchedSourceComponentCatalogMultiplicity[multiplicity]++;
	
			//Fill Tree
			if(nMatchedObjects>0){
				//matchedSourceInfo->Fill();
				matchInfoCollection.push_back(matchInfo);
			}

		}//end (j index) loop component pars in this tile

	}//end (i index) loop tiles


	cout<<"== MATCH RESULTS =="<<endl;
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

	//Fill TTree
	for(size_t i=0;i<matchInfoCollection.size();i++)
	{
		ComponentPars* componentPars= matchInfoCollection[i]->cpars;
		
		//Set TTree data
		matchedObjects.clear();
		fluxDiff.clear();
		matchedObjectMultiplicity.clear();
		nMatchedObjects= 0;
		sourceName= componentPars->sname;
		sourcePosX= componentPars->X0;
		sourcePosY= componentPars->Y0;
		sourcePosX_wcs= componentPars->X0_wcs;
		sourcePosY_wcs= componentPars->Y0_wcs;
		sourceFluxDensity= componentPars->fluxDensity;
		sourceFlux= componentPars->flux;
		sourceFluxDensityErr= componentPars->fluxDensityErr;
		sourceFluxErr= componentPars->fluxErr;

		for(size_t j=0;j<(matchInfoCollection[i]->objs).size();j++)
		{
			AstroObject* obj= (matchInfoCollection[i]->objs)[j];
			std::string name= obj->name;
			int mult= matchedSourceNMatchMap[name];

			matchedObjects.push_back(obj);
			nMatchedObjects++;
			double dflux= sourceFlux- obj->flux;
			fluxDiff.push_back(dflux);
			matchedObjectMultiplicity.push_back(mult);

		}//end loop matched objects

		//Fill TTree 
		matchedSourceInfo->Fill();

	}//end loop matched components

	
	//Fill Tree
	matchOptionTree->Fill();

	return 0;

}//close FindSourceMatchesInTiles()







bool HaveSourceComponentMatch(SourceComponentMatchPars& scompmatchpars,ComponentPars* cpars,AstroObject* obj)
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
	if(!cpars || !obj){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given for component pars and/or catalog object!");
		#endif
		return false;
	}

	//## Check fit ellipses
	TEllipse* ellipse1= cpars->fitEllipse;
	TEllipse* ellipse2= obj->GetFitEllipse();
	if(matchSourceComponentsByOverlap && (!ellipse1 || !ellipse2)){		
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/both fit ellipses are nullptr!");
		#endif
		return false;
	}

	//## Check if flux information is available on Astro object
	if(matchSourcesByFlux && !obj->hasFluxInfo){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Flux info is not available for catalog object "<<obj->name<<", returning no match!");
		#endif
		return false;
	}

	//## Check ellipse centroid sky distance (if requested)
	double posThr= compMatchPosThr/3600.;//convert in deg as ellipse coords are in deg
	double posHighThr= compMatchPosHighThr/3600.;//convert in deg as ellipse coords are in deg
	double Xc_1= cpars->X0_wcs;//ellipse1->GetX1();
	double Yc_1= cpars->Y0_wcs;//ellipse1->GetY1();
	double Xc_2= obj->x;//ellipse2->GetX1();
	double Yc_2= obj->y;//ellipse2->GetY1();
	double posDist_Euclidean= MathUtils::GetEuclideanDist(Xc_1,Yc_1,Xc_2,Yc_2);
	double posDist_Haversine= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);
	double posDist= 0.;
	if(useWCSSimpleGeometry) posDist= posDist_Euclidean;
	else posDist= posDist_Haversine;


	//## Compute ellipse area, distance and overlap
	double R1_1= 0;
	double R2_1= 0;
	double Theta_1= 0;
	double ellipseArea_1= 0;
	if(ellipse1){
		R1_1= ellipse1->GetR1();	
		R2_1= ellipse1->GetR2();
		Theta_1= ellipse1->GetTheta();
		ellipseArea_1= MathUtils::ComputeEllipseArea(ellipse1);
	}
	double R1_2= 0;
	double R2_2= 0;
	double Theta_2= 0;
	double ellipseArea_2= 0;
	if(ellipse2){
		R1_2= ellipse2->GetR1();	
		R2_2= ellipse2->GetR2();
		Theta_2= ellipse2->GetTheta();
		ellipseArea_2= MathUtils::ComputeEllipseArea(ellipse2);
	}
	
	double ellipseOverlapArea= -1;
	double err= 0;
	int rtn= 0;
	double overlapAreaFraction_1= 0;
	double overlapAreaFraction_2= 0;
	
	if(ellipse1 && ellipse2){
		if(MathUtils::ComputeEllipseOverlapArea(ellipseOverlapArea,err,rtn,ellipse1,ellipse2)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute ellipse overlap area (err status="<<rtn<<"), return no match!");
			#endif
			return false;
		}
		overlapAreaFraction_1= ellipseOverlapArea/ellipseArea_1;
		overlapAreaFraction_2= ellipseOverlapArea/ellipseArea_2;
	}	
	
	//## Check overlap
	bool hasPosLargeMatch= (fabs(posDist)<=posThr);
	bool hasPosSmallMatch= (fabs(posDist)<=posHighThr);
	bool noPosMatch= (hasPosSmallMatch || hasPosLargeMatch);
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
	double flux_1= cpars->fluxDensity;
	double flux_2= obj->fluxDensity;
	double fluxRelDiff= 0;
	if(flux_2>0){
		fluxRelDiff= (flux_1-flux_2)/flux_2;
	}

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
				DEBUG_LOG("DISJOINT ELLIPSES: Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
			#endif
			return false;
		}

		//- Full overlap: Ellipse 1 inside ellipse 2 (overlap=Area1)
		else if(rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 ){
			//Check area overlap ratio
			if(applySourceComponentAreaRatioThr && overlapAreaFraction_2<compMatchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE1 INSIDE ELLIPSE2 (OVERLAP BELOW THR): overlapArea2="<<overlapAreaFraction_2<<"<"<<compMatchOverlapThr<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			//Check small pos match
			if(matchSourceComponentsByPos && !hasPosSmallMatch){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}
		}

		//- Full overlap: Ellipse 2 inside ellipse 1 (overlap=Area2)
		else if(rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1){
			//Check area overlap ratio
			if(applySourceComponentAreaRatioThr && overlapAreaFraction_1<compMatchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE2 INSIDE ELLIPSE1 (OVERLAP BELOW THR): overlapArea1="<<overlapAreaFraction_1<<"<"<<compMatchOverlapThr<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			//Check small pos match
			if(matchSourceComponentsByPos && !hasPosSmallMatch){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}
		}

		//- Partial overlap
		else {
			//- Large overlap
			if(hasOverlapLargeMatch){	
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE LARGE OVERLAP: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapThr<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				
				//Check small pos match
				if(matchSourceComponentsByPos && !hasPosSmallMatch){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
					#endif
					return false;
				}
			}//close if

			//- Small overlap
			else if(hasOverlapSmallMatch){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE SMALL OVERLAP: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapLowThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapLowThr<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif

				//Check large pos match
				if(matchSourceComponentsByPos && !hasPosLargeMatch){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("NO ELLIPSE POS LARGE MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posThr*3600<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
					#endif
					return false;
				}
			}//close else if
			
			//- Partial overlap not sufficient
			else {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("ELLIPSE OVERLAP BELOW SMALL & HIGH THR: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapThr<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			#ifdef LOGGING_ENABLED
				DEBUG_LOG("ELLIPSE OVERLAP MATCH: overlapArea1="<<overlapAreaFraction_1<<", overlapArea2="<<overlapAreaFraction_2<<", compMatchOverlapThr="<<compMatchOverlapThr<<", compMatchOverlapLowThr="<<compMatchOverlapLowThr<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
			#endif
		}//close else

	}//close if check overlap
	else{
		//## Check ellipse centroid sky distance (if requested)
		if(matchSourceComponentsByPos && !hasPosLargeMatch){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("NO ELLIPSE POS MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posThr*3600<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}");
			#endif
			return false;
		}

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("ELLIPSE POS MATCH: posDist(arcsec)="<<posDist*3600<<"<"<<posThr*3600<<", Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}");
		#endif
	}


	//## Check flux rel difference (if requested)	
	if(matchSourcesByFlux){
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("SOURCE MATCH: Source {name="<<cpars->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, CatalogObject {name="<<obj->name<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, overlapArea1="<<overlapAreaFraction_1<<", overlapArea2="<<overlapAreaFraction_2<<", posDist(arcsec)="<<posDist*3600<<", ");
	#endif

	return true;

}//close HaveSourceComponentMatch()


int ReadRegions()
{
	//Read and parse DS9 regions
	if (DS9RegionParser::Parse(m_regions,regionFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse DS9 region file "<<regionFileName<<"!");
		#endif
		return -1;
	}

	//Check if empty
	if(m_regions.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No regions read from file "<<regionFileName<<"!");
		#endif
		return -1;
	}
	
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<m_regions.size()<<" regions read from file "<<regionFileName<<"...");
	#endif	

	return 0;

}//close ReadRegions()


int ReadData(std::string filename)
{
	//Open files
	TFile* inputFile= new TFile(filename.c_str(),"READ");
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
	std::vector<ComponentPars*> pars;
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
		//NB: Generate uniform in annulus (e.g. see https://stackoverflow.com/questions/9048095/create-random-number-within-an-annulus)
		double offsetX= 0;
		double offsetY= 0;
		if(shuffleSources){
			double A = 2/(shuffleRMax*shuffleRMax - shuffleRMin*shuffleRMin);
			double r= sqrt(2*gRandom->Uniform()/A + shuffleRMin*shuffleRMin);
			double theta= gRandom->Uniform(0,2*TMath::Pi());
			offsetX= r*cos(theta);
			offsetY= r*sin(theta);
			offsetX/= 3600;//convert in deg
			offsetY/= 3600;//convert in deg
		}

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
		std::vector<ComponentPars*> thisPars;
		if(FillSourceComponentPars(thisPars,source,catalogIndex,sourceIndex,-1,wcs,wcsType,offsetX,offsetY)<0){
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
	m_sourceComponentPars.insert(m_sourceComponentPars.end(),pars.begin(),pars.end());

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source min/max WCS coordinate range: (Xmin,Xmax)=("<<sourceXmin<<","<<sourceXmax<<"), (Ymin,Ymax)=("<<sourceYmin<<","<<sourceYmax<<")");
		INFO_LOG("#"<<m_sources.size()<<" sources to be cross-matched (#"<<m_sourceComponentPars.size()<<" source pars added) ...");
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
	for(size_t i=0;i<m_sourceComponentPars.size();i++)
	{
		double X0_wcs= m_sourceComponentPars[i]->X0_wcs;
		double Y0_wcs= m_sourceComponentPars[i]->Y0_wcs;
		std::string sourceName= m_sourceComponentPars[i]->sname;

		long int tileIndex= MathUtils::FindGrid2DBin(
			X0_wcs,Y0_wcs,
			nTilesX,tileXmin,tileXmax,tileXstep,
			nTilesY,tileYmin,tileYmax,tileYstep
		);

		if(tileIndex>=0 && tileIndex<(long int)(tileDataList.size())){//Add to tile data
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding source (name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<") to list...");
			#endif
			tileDataList[tileIndex]->AddSourceComponentPars(m_sourceComponentPars[i]);
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
		long int nSourceComponentPars= tileDataList[i]->GetNSourceComponentPars();
		#ifdef LOGGING_ENABLED
			INFO_LOG("Tile no. "<<i+1<<": #"<<nSourceComponentPars<<" source components to be matched...");
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
		tileDataList[i]->InitCatalogObjects(nCatalogs);
	}
	for(int i=0;i<nCatalogs;i++){
		m_catalogObjects.push_back(std::vector<AstroObject*>());
	}

	//Finally reading source data
	for(size_t i=0;i<fileNames.size();i++){
		if(ReadAstroObjectData(fileNames[i],i)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read source data for file no. "<<i+1<<"!");
			#endif
			return -1;
		}
	}//end loop files

	//Printing num sources in each tile
	for(size_t i=0;i<tileDataList.size();i++){
		long int nCatalogObjects= tileDataList[i]->GetNCatalogObjects();
		std::stringstream ss;
		ss<<"Tile no. "<<i+1<<": {";
		for(long int j=0;j<nCatalogObjects-1;j++){
			ss<<"#"<<tileDataList[i]->GetNObjectsInCatalog(j)<<",";
		}
		ss<<"#"<<tileDataList[i]->GetNObjectsInCatalog(nCatalogObjects-1)<<"} catalog objects ...";
		#ifdef LOGGING_ENABLED
			INFO_LOG(ss.str());
		#endif
	}

	return 0;

}//close ReadCatalogData()


int ReadAstroObjectData(std::string filename,int catalogIndex)
{
	//## Read catalog data
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading objects from catalog file "<<filename<<" ...");
	#endif
	int status= 0;
	std::vector<AstroObject*> astroObjects;
	if(catalogType==eMGPS_CATALOG){
		AstroObjectParser::ParseMGPSData(astroObjects,filename,'|');
	}
	else if(catalogType==eSELAVY_CATALOG){
		AstroObjectParser::ParseSelavyData(astroObjects,filename,' ');
	}
	else if(catalogType==eSIMBAD_CATALOG){
		AstroObjectParser::ParseSimbadData(astroObjects,filename,'|');
	}
	else if(catalogType==eNED_CATALOG){
		AstroObjectParser::ParseNedData(astroObjects,filename,'|');
	}
	else if(catalogType==eHASH_CATALOG){
		AstroObjectParser::ParseHASHData(astroObjects,filename,'|');
	}
	else if(catalogType==eWISE_HII_CATALOG){
		AstroObjectParser::ParseWiseHIIData(astroObjects,filename,'|');
	}
	else if(catalogType==eATNF_PSR_CATALOG){
		AstroObjectParser::ParseATNFPsrData(astroObjects,filename,'|');
	}
	else if(catalogType==eWOLF_RAYET_CATALOG){
		AstroObjectParser::ParseWRCatData(astroObjects,filename,'|');
	}
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/unsupported catalog type (type="<<catalogType<<")!");
		#endif
		return -1;
	}

	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read objects from catalog file "<<filename<<"!");
		#endif
		return -1;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<astroObjects.size()<<" objects read from catalog file "<<filename<<" ...");
	#endif

	//## Select catalog data
	m_catalogObjects[catalogIndex].clear();

	//Loop over objects
	for(size_t i=0;i<astroObjects.size();i++)
	{
		int id= astroObjects[i]->id;
		double x= astroObjects[i]->x;
		double y= astroObjects[i]->y;

		//Select by object identifier
		if(selectObjectByType){
			bool skipObject= true;
			for(size_t j=0;j<otypes.size();j++){
				if( otypes[j]==-1 || id==otypes[j]) {
					skipObject= false;
					break;
				}
			}
			if(skipObject) continue;
		}

		//Select objects in regions
		if(selectObjectsInRegions)
		{
			//Loop over regions
			bool isInsideRegion= false;
			for(size_t j=0;j<m_regions.size();j++)
			{
				if(m_regions[j]->IsPointInsideRegion(x,y)){
					isInsideRegion= true;
					break;
				}
			}//end loop regions

			if(!isInsideRegion) continue;

		}//close if

		//Add to collection 
		m_catalogObjects[catalogIndex].push_back(astroObjects[i]);
		
	}//end loop objects
	

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<m_catalogObjects[catalogIndex].size()<<" objects selected to be cross-matched in catalog no. "<<catalogIndex<<" ...");
	#endif


	//## Add catalog object data to tile
	for(size_t i=0;i<m_catalogObjects[catalogIndex].size();i++)
	{
		double X0_wcs= m_catalogObjects[catalogIndex][i]->x;
		double Y0_wcs= m_catalogObjects[catalogIndex][i]->y;
		std::string name= m_catalogObjects[catalogIndex][i]->name;
		
		long int tileIndex= MathUtils::FindGrid2DBin(
			X0_wcs,Y0_wcs,
			nTilesX,tileXmin,tileXmax,tileXstep,
			nTilesY,tileYmin,tileYmax,tileYstep
		);

		if(tileIndex>=0 && tileIndex<(long int)(tileDataList.size())){//Add to tile data
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding catalog object (catalogId="<<catalogIndex<<", name="<<name<<", pos("<<X0_wcs<<","<<Y0_wcs<<") to tile...");
			#endif
			int status= tileDataList[tileIndex]->AddCatalogObject(catalogIndex,m_catalogObjects[catalogIndex][i]);
			if(status<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to add catalog object no. "<<i+1<<" from catalog "<<catalogIndex<<" to tile!");
				#endif
				continue;
			}
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot find tile index or invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<") for catalog object (catalogId="<<catalogIndex<<", name="<<name<<", pos("<<X0_wcs<<","<<Y0_wcs<<"), check tile index calculation!");
			#endif
			continue;
		}

	}//end loop pars

	
	return 0;

}//close ReadAstroObjectData()

int FillSourceComponentPars(std::vector<ComponentPars*>& pars,Source* aSource,int catalogIndex,int sourceIndex,int nestedSourceIndex,WCS* wcs,int coordSystem,double offsetX,double offsetY)
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
	bool useFWHM= true;
	bool convertToWCS= true;
	int pixOffset= 0;

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
			CodeUtils::DeletePtrCollection<ComponentPars>(pars);
			return -1;
		}

		//Check ellipses and pars size
		if(nComponents!=(int)(fitEllipses.size())){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Number of fit components shall be equal to fitted ellipses (this should not occur)!");
			#endif
			CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
			CodeUtils::DeletePtrCollection<ComponentPars>(pars);
			return -1;
		}


		//Fill component pars & ellipse
		ComponentPars* thisComponentPars= 0;

		for(int k=0;k<nComponents;k++){
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));
			double X0= 0;
			double Y0= 0;
			fitPars.GetComponentPosition(X0,Y0,k);
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			double fluxDensityErr= fitPars.GetComponentFluxDensityErr(k);
			double flux= fluxDensity/beamArea;
			double fluxErr= fluxDensityErr/beamArea;

			//Apply offset to ellipse (if shuffle source is enabled)
			double Cx= fitEllipses[k]->GetX1();
			double Cy= fitEllipses[k]->GetY1();
			if(shuffleSources){
				Cx+= offsetX;
				Cy+= offsetY;
				fitEllipses[k]->SetX1(Cx);
				fitEllipses[k]->SetY1(Cy);
			}

			//Fill component par
			thisComponentPars= new ComponentPars(k,sourceIndex,nestedSourceIndex,catalogIndex);
			thisComponentPars->sname= sname;
			thisComponentPars->X0= X0;
			thisComponentPars->Y0= Y0;
			thisComponentPars->X0_wcs= Cx;
			thisComponentPars->Y0_wcs= Cy;
			thisComponentPars->fitEllipse= fitEllipses[k];
			thisComponentPars->fluxDensity= fluxDensity;	
			thisComponentPars->flux= flux;
			thisComponentPars->fluxDensityErr= fluxDensityErr;	
			thisComponentPars->fluxErr= fluxErr;
			thisComponentPars->isSelected= isSelected;

			//Append this source pars to collection
			pars.push_back(thisComponentPars);

		}//end loop components
	}//close if has fit info

	
	
	//## Fill pars for nested sources (if any)
	if(hasNestedSources){
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			if(FillSourceComponentPars(pars,nestedSources[j],catalogIndex,sourceIndex,j,wcs,coordSystem,offsetX,offsetY)<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
				#endif
				return -1;
			}
		}//end loop sources
	}//close if nested sources
	
	return 0;

}//close FillSourceComponentPars()

void Save()
{
	//Save TTree to file
	if(outputFile && outputFile->IsOpen())
	{		
		outputFile->cd();
		if(matchedSourceInfo){
			matchedSourceInfo->Write();
		}
		if(matchOptionTree){
			matchOptionTree->Write();
		}
		outputFile->Close();
	}

}//close Save()

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

