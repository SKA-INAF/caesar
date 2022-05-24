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
#include <TKey.h>
#include <TROOT.h>

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
	cout<<"-c, --catalogType=[CATALOG_TYPE] \t Catalog type id {1=SIMBAD, 2=NED, 3=MGPS, 4=HASH, 5=WISE-HII, 6=ATNF-PSR, 7=WOLF-RAYET, 8=SELAVY, 9=MASH, 10=GAIA} (default=SIMBAD)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (root format) in which to store source match info"<<endl;
	cout<<"-J, --soutput=[OUTPUT_FILE] \t Output file name (root format) in which to store source catalog with spectral index info added"<<endl;
	cout<<"-n, --nx=[NX] \t Number of divisions along X (default=4)"<<endl;
	cout<<"-N, --ny=[NY] \t Number of divisions along Y (default=4)"<<endl;
	cout<<"-q, --selectObjectsInRegions \t Select catalog object inside given DS9 regions (default=no)"<<endl;
	cout<<"-Q, --regionFile \t DS9 region filename used to select catalog objects"<<endl;
	cout<<"-m, --matchSourcesByFlux \t Match sources by flux (default=false)"<<endl;
	cout<<"-M, --matchFluxRelDiffThr \t Flux relative difference threshold (default=0.5)"<<endl;

	cout<<"-O, --matchByOverlap \t Match source islands by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-A, --matchOverlapThr \t Source island contour overlap fraction (wrt to total area) threshold (default=0.5)"<<endl;
	cout<<"-P, --matchByPos \t Match source islands by position centroid distance (default=false)"<<endl;
	cout<<"-T, --matchPosThr=[POS_THESHOLD] \t Source island centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	
	cout<<"-p, --matchComponentByPos \t Match source fit components by position centroid distance (default=false)"<<endl;
	cout<<"-b, --matchComponentByOverlap \t Match source fit component by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-t, --compMatchPosThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	cout<<"-e, --compMatchPosHighThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=5)"<<endl;
	cout<<"-a, --compMatchOverlapThr \t Source fit component contour overlap fraction (wrt to total area) threshold (default=0.8)"<<endl;
	cout<<"-d, --compMatchOverlapLowThr \t Source fit component contour overlap fraction (wrt to total area) low threshold (default=0.2)"<<endl;
	cout<<"-l, --applyComponentAreaRatioThr \t If enabled and if source fit component is fully enclosed inside another (or viceversa) a match requires that the source1/source2 area ratio is higher than the component overlap threshold. If disabled, assume a match whenever a source component is fully enclosed in a source (or viceversa) (default=false)"<<endl;
	cout<<"-f, --filterByMorphId \t Consider only sources with given morph id when searching the match (default=no)"<<endl;
	cout<<"-s, --selectedMorphId=[MORPH_ID] \t Source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED, 5=DIFFUSE) (default=-1)"<<endl;
	cout<<"-F, --filterObjectByType \t Consider only object with given type when searching the match (default=no)"<<endl;
	cout<<"-S, --selectedObjectType=[TYPE] \t Object types to be crossmatched (default=-1)"<<endl;
	cout<<"-g, --shuffleSources \t Randomize sources and catalog sources in an annulus centred on their original position (default=no)"<<endl;
	cout<<"-r, --shuffleRmin=[Rmin] \t Minimum annulus radius in arcsec for source shuffling (default=30 arcsec)"<<endl;
	cout<<"-R, --shuffleRmax=[Rmax] \t Maximum annulus radius in arcsec for source shuffling (default=50 arcsec)"<<endl;
	cout<<"-B, --shuffleInRegion \t Randomize sources in boundary region (rather than in annulus) (default=no)"<<endl;
	
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
	{ "soutput", required_argument, 0, 'J' },	
	{ "selectObjectsInRegions", no_argument, 0, 'q' },	
	{ "regionFile", required_argument, 0, 'Q' },	
	{ "nx", required_argument, 0, 'n' },
	{ "ny", required_argument, 0, 'N' },
	{ "matchSourcesByFlux", no_argument, 0, 'm'},
	{ "matchFluxRelDiffThr", required_argument, 0, 'M'},
	{ "matchByOverlap", no_argument, 0, 'O'},	
	{ "matchOverlapThr", required_argument, 0, 'A'},
	{ "matchByPos", no_argument, 0, 'P'},
	{ "matchPosThr", required_argument, 0, 'T'},
	{ "matchComponentByPos", no_argument, 0, 'p'},
	{ "matchComponentByOverlap", no_argument, 0, 'b'},	
	{ "compMatchPosThr", required_argument, 0, 't'},
	{ "compMatchPosHighThr", required_argument, 0, 'e'},
	{ "compMatchOverlapThr", required_argument, 0, 'a'},
	{ "compMatchOverlapLowThr", required_argument, 0, 'd'},
	{ "applyComponentAreaRatioThr", no_argument, 0, 'l'},	
	{ "filterByMorphId", no_argument, 0, 'f'},
	{ "selectedMorphId", required_argument, 0, 's'},
	{ "filterObjectByType", no_argument, 0, 'F'},
	{ "selectedObjectType", required_argument, 0, 'S'},
	{ "shuffleSources", no_argument, 0, 'g'},	
	{ "shuffleRmin", required_argument, 0, 'r'},
	{ "shuffleRmax", required_argument, 0, 'R'},
	{ "shuffleInRegion", no_argument, 0, 'B'},		
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string catalogFileName= "";
bool selectObjectsInRegions= false;
std::string regionFileName= "";
int verbosity= 4;//INFO level
bool useWCSSimpleGeometry= true;

bool matchSourcesByOverlap= false;
float matchOverlapThr= 0.5;//fraction of overlap above which two sources are matched
bool matchSourcesByPos= false;
float matchPosThr= 2.5;//dist in arcsec below which two sources can be matched

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
bool selectSourceByMorphId= false;//default=all true sources searched 
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
bool shuffleInRegion= false;

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
	eMASH_CATALOG=9,
	eGAIA_CATALOG=10,
};
int catalogType= 1;


//Globar vars
TFile* inputFile= 0;
TFile* outputFile= 0;
std::string outputFileName= "sources.root";
TTree* outputTree= 0;
TFile* matchOutputFile= 0;
std::string matchOutputFileName= "MatchOutput.root";
TTree* matchedSourceInfo= 0;
TTree* matchedSourceComponentInfo= 0;
TTree* matchSummaryTree= 0;
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
int nTotSources= 0;
int nMatchedSources= 0;
int nTotSourceComponents= 0;
int nMatchedSourceComponents= 0;
int nCatalogs= 0;
std::vector<int> nMatchedSourcesPerCatalog;
std::vector<int> nMatchedSourceCatalogMultiplicity;
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
//     SourcePars struct
//====================================
struct SourcePars 
{
	SourcePars(int _catalogIndex,int _sourceIndex,int _nestedSourceIndex=-1)
		: catalogIndex(_catalogIndex), sourceIndex(_sourceIndex), nestedSourceIndex(_nestedSourceIndex)
	{
		sname= "";	
		contour= nullptr;
		S= 0;
		fluxDensity= 0;
		flux= 0;
		fluxDensityErr= 0;
		fluxErr= 0;
		X0= 0;
		Y0= 0;
		X0_wcs= 0;
		Y0_wcs= 0;
		//componentPars.clear();
	}
	~SourcePars(){
		CodeUtils::DeletePtr<Contour>(contour);
		//CodeUtils::DeletePtrCollection<ComponentPars>(componentPars);
	}
	/*
	void AddComponentPars(ComponentPars* pars){
		componentPars.push_back(pars);
	}
	*/
	int catalogIndex;
	int sourceIndex;
	int nestedSourceIndex;
	std::string sname;
	double X0;
	double Y0;
	double X0_wcs;
	double Y0_wcs;
	Contour* contour;
	double S;//pixel-integrated flux
	double fluxDensity;//fitted flux density
	double flux;//fitted flux
	double fluxDensityErr;//fitted flux density err
	double fluxErr;//fitted flux
	//std::vector<ComponentPars*> componentPars;
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
		sourcePars.clear();
	}

	//Destructor 
	~TileData()
	{
		CodeUtils::DeletePtrCollection<ComponentPars>(sourceComponentPars);
		CodeUtils::DeletePtrCollection<SourcePars>(sourcePars);

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
	int AddSourcePars(SourcePars* pars){
		if(!pars) return -1;
		sourcePars.push_back(pars);
		return 0;
	}

	//Get number of source pars in this tile
	long int GetNSourcePars(){
		return static_cast<long int>(sourcePars.size());
	}

	//Add source component pars
	int AddSourceComponentPars(ComponentPars* pars){
		if(!pars) return -1;
		sourceComponentPars.push_back(pars);
		return 0;
	}

	//Get number of sources component pars in this tile
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

	//Source island pars for this tile
	std::vector<SourcePars*> sourcePars;

	//Source component pars for this tile
	std::vector<ComponentPars*> sourceComponentPars;

	//Catalog objects for this tile
	std::vector<std::vector<AstroObject*>> catalogObjects;

};//close TileData()

std::vector<TileData*> tileDataList;

std::vector<Source*> m_sources;
Source* m_source= 0;
std::vector<SourcePars*> m_sourcePars;
std::vector<ComponentPars*> m_sourceComponentPars;
std::vector<std::vector<AstroObject*>> m_catalogObjects;
std::vector<DS9Region*> m_regions;
Contour* m_regionContour= 0;


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int InitGrid(float xmin,float xmax,float xstep,float ymin,float ymax,float ystep);
int FindSourceMatchesInTiles();
int FindSourceComponentMatchesInTiles();
int ComputeSpectralIndices();
int ReadRegions();
int ReadData(std::string filename);
int ReadCatalogData(std::string filelist);
int ReadAstroObjectData(std::string filename,int collectionIndex);
int FillSourcePars(std::vector<SourcePars*>& spars,std::vector<ComponentPars*>& cpars,Source* aSource,int collectionIndex,int sourceIndex,int nestedSourceIndex=-1,WCS* wcs=0,int coordSystem=0,double offsetX=0,double offsetY=0);
bool HaveSourceComponentMatch(SourceComponentMatchPars& scompmatchpars,ComponentPars* cpars,AstroObject* obj);
bool HaveSourceIslandMatch(SourceMatchPars& smatchpars,SourcePars* spars,AstroObject* obj);
int GetRandPosInAnnulus(double& offsetX,double& offsetY);
int GetRandPosInRegion(double& offsetX,double& offsetY,Contour* contour,double xmin,double xmax,double ymin,double ymax,int nMaxGen=1000);
void Save();
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
		INFO_LOG("Finding source island matches with astro catalog...");
	#endif
	if(FindSourceMatchesInTiles()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Finding source island matches with astro catalog failed!");
		#endif
		return -1;
	}

	//==================================
	//== Find source component matches
	//==================================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Finding source component matches with astro catalog");
	#endif
	if(FindSourceComponentMatchesInTiles()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Finding source component matches with astro catalog failed!");
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

	while((c = getopt_long(argc, argv, "hi:C:c:o:v:qQ:n:N:mM:OA:PT:Llt:e:pba:d:fs:FS:gr:R:J:B",options_tab, &option_index)) != -1) {
    
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
				matchOutputFileName= std::string(optarg);	
				break;	
			}
			case 'J':	
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
				selectSourceByMorphId= true;
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
			case 'B':
			{
				shuffleInRegion= true;
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
			catalogType!=eSELAVY_CATALOG &&
			catalogType!=eMASH_CATALOG && 
			catalogType!=eGAIA_CATALOG
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
		INFO_LOG("matchSourcesByPos? "<<matchSourcesByPos<<", matchPosThr(arcsec)="<<matchPosThr);
		INFO_LOG("matchSourcesByOverlap? "<<matchSourcesByOverlap<<", matchOverlapThr="<<matchOverlapThr);
		INFO_LOG("matchSourceComponentsByPos? "<<matchSourceComponentsByPos<<", compMatchPosThr(arcsec)="<<compMatchPosThr<<", compMatchPosHighThr(arcsec)="<<compMatchPosHighThr);
		INFO_LOG("matchSourceComponentsByOverlap? "<<matchSourceComponentsByOverlap<<", compMatchOverlapThr="<<compMatchOverlapThr<<", compMatchOverlapLowThr="<<compMatchOverlapLowThr);
		INFO_LOG("applySourceComponentAreaRatioThr? "<<applySourceComponentAreaRatioThr);
	#endif

	return 0;

}//close ParseOptions()

int Init()
{
	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");
	
	//Init output source TTree
	if(!outputTree) outputTree= new TTree("SourceInfo","SourceInfo");
	m_source= 0;
	outputTree->Branch("Source",&m_source);

	//Open match output file
	if(!matchOutputFile) matchOutputFile= new TFile(matchOutputFileName.c_str(),"RECREATE");
	
	//Init source match TTree
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

	//Init source component match TTree
	if(!matchedSourceComponentInfo) matchedSourceComponentInfo= new TTree("SourceComponentMatchInfo","SourceComponentMatchInfo");
	matchedSourceComponentInfo->Branch("name",&sourceName);
	matchedSourceComponentInfo->Branch("x",&sourcePosX);
	matchedSourceComponentInfo->Branch("y",&sourcePosY);
	matchedSourceComponentInfo->Branch("x_wcs",&sourcePosX_wcs);
	matchedSourceComponentInfo->Branch("y_wcs",&sourcePosY_wcs);
	matchedSourceComponentInfo->Branch("fluxDensity",&sourceFluxDensity);
	matchedSourceComponentInfo->Branch("fluxDensityErr",&sourceFluxDensityErr);	
	matchedSourceComponentInfo->Branch("flux",&sourceFlux);
	matchedSourceComponentInfo->Branch("fluxErr",&sourceFluxErr);
	matchedSourceComponentInfo->Branch("nMatchedObjects",&nMatchedObjects);
	matchedSourceComponentInfo->Branch("matchedObjects",&matchedObjects);
	//matchedSourceComponentInfo->Branch("fluxDiff",&fluxDiff);
	matchedSourceComponentInfo->Branch("matchedObjectMultiplicity",&matchedObjectMultiplicity);
	

	//Init match option TTree
	if(!matchSummaryTree) matchSummaryTree= new TTree("MatchSummary","MatchSummary");
	matchSummaryTree->Branch("matchSourcesByPos",&matchSourcesByPos);
	matchSummaryTree->Branch("matchPosThr",&matchPosThr);
	matchSummaryTree->Branch("matchSourcesByOverlap",&matchSourcesByOverlap);
	matchSummaryTree->Branch("matchOverlapThr",&matchOverlapThr);
	matchSummaryTree->Branch("matchSourceComponentsByPos",&matchSourceComponentsByPos);
	matchSummaryTree->Branch("compMatchPosThr",&compMatchPosThr);
	matchSummaryTree->Branch("compMatchPosHighThr",&compMatchPosHighThr);
	matchSummaryTree->Branch("matchSourceComponentsByOverlap",&matchSourceComponentsByOverlap);
	matchSummaryTree->Branch("compMatchOverlapThr",&compMatchOverlapThr);
	matchSummaryTree->Branch("compMatchOverlapLowThr",&compMatchOverlapLowThr);
	matchSummaryTree->Branch("selectSourceByMorphId",&selectSourceByMorphId);
	matchSummaryTree->Branch("applySourceAreaRatioThr",&applySourceAreaRatioThr);
	matchSummaryTree->Branch("applySourceComponentAreaRatioThr",&applySourceComponentAreaRatioThr);
	matchSummaryTree->Branch("nCatalogs",&nCatalogs);
	matchSummaryTree->Branch("nTotSources",&nTotSources);
	matchSummaryTree->Branch("nMatchedSources",&nMatchedSources);
	matchSummaryTree->Branch("nMatchedSourcesPerCatalog",&nMatchedSourcesPerCatalog);
	matchSummaryTree->Branch("nMatchedSourceCatalogMultiplicity",&nMatchedSourceCatalogMultiplicity);
	matchSummaryTree->Branch("nTotSourceComponents",&nTotSourceComponents);
	matchSummaryTree->Branch("nMatchedSourceComponents",&nMatchedSourceComponents);
	matchSummaryTree->Branch("nMatchedSourceComponentsPerCatalog",&nMatchedSourceComponentsPerCatalog);
	matchSummaryTree->Branch("nMatchedSourceComponentCatalogMultiplicity",&nMatchedSourceComponentCatalogMultiplicity);

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
	nTotSources= 0;
	nMatchedSources= 0;
	std::map<std::string,int> matchedSourceNMatchMap;

	struct MatchInfo 
	{
		MatchInfo(){
			spars= 0;
			objs.clear();
		}
		SourcePars* spars;
		std::vector<AstroObject*> objs;
	};
	std::vector<MatchInfo*> matchInfoCollection;
	MatchInfo* matchInfo= 0;

	//Loop over tiles and find matches with regions in the same tile and in neighbor tiles
	for(size_t i=0;i<tileDataList.size();i++)
	{
		long int NSourcePars= tileDataList[i]->GetNSourcePars();
		long int NCatalogs= tileDataList[i]->GetNCatalogObjects();
		std::vector<size_t> neighbors= tileDataList[i]->neighborTileIndexes;
		
		//Init counters
		if(nMatchedSourcesPerCatalog.empty()){
			nMatchedSourcesPerCatalog= std::vector<int>(NCatalogs,0);
		}
		if(nMatchedSourceCatalogMultiplicity.empty()){
			nMatchedSourceCatalogMultiplicity= std::vector<int>(NCatalogs+1,0);
		}

		for(long int j=0;j<NSourcePars;j++) 
		{		
			//Get source component 
			SourcePars* sourcePars= (tileDataList[i]->sourcePars)[j];
			
			#ifdef LOGGING_ENABLED
			if(j%100==0) INFO_LOG("#"<<j+1<<"/"<<NSourcePars<<" source islands to be matched in tile no. "<<i+1<<" against "<<NCatalogs<<" catalogs ...");
			#endif

			//Count total number of sources					
			nTotSources++;
			std::vector<int> matchCounter(NCatalogs,0);
			bool matchFound= false;
			nMatchedObjects= 0;

			//Create MatchInfo 
			matchInfo= new MatchInfo;
			matchInfo->spars= sourcePars;

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
					DEBUG_LOG("#"<<NSourcePars<<" sources to be cross-matched against #"<<NObjectsInCatalog<<" objects (catalog #"<<s+1<<") in tile no. "<<i+1<<" ...");
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
					SourceMatchPars smatchpars;
					bool hasMatch= HaveSourceIslandMatch(smatchpars,sourcePars,catalogObject);
					if(!hasMatch) {
						continue;
					}
					matchFound= true;
					matchCounter[s]++;	
					matchedSourceNMatchMap[catalogObject->name]++;
					nMatchedObjects++;

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
						DEBUG_LOG("#"<<NSourcePars<<" sources to be cross-matched against #"<<NObjectsInCatalog_neighbor<<" catalog objects (catalog #"<<s+1<<") in tile no. "<<neighborIndex<<" ...");
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
						SourceMatchPars smatchpars;
						bool hasMatch= HaveSourceIslandMatch(smatchpars,sourcePars,catalogObject);
						if(!hasMatch) {
							continue;
						}
						matchFound= true;
						matchCounter[s]++;
						matchedSourceNMatchMap[catalogObject->name]++;
						nMatchedObjects++;

						(matchInfo->objs).push_back(catalogObject);

					}//end (k index) loop objects in catalog
				}//end (n index) loop neighbor tiles
	
			}//end (s index) loop catalogs

			//Update counters
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Update source match counters...");
			#endif
			if(matchFound){
				nMatchedSources++;
			}
			std::vector<int> sourceMatchMultiplicityCounter(NCatalogs,0);
			for(size_t c=0;c<matchCounter.size();c++){
				if(matchCounter[c]>0) {
					nMatchedSourcesPerCatalog[c]++;
					sourceMatchMultiplicityCounter[c]++;
				}
			}

			//Compute catalog match multiplicity
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Compute catalog match multiplicity...");
			#endif
			int multiplicity= std::accumulate(sourceMatchMultiplicityCounter.begin(), sourceMatchMultiplicityCounter.end(), 0);
			nMatchedSourceCatalogMultiplicity[multiplicity]++;
	
			//Fill Tree
			if(nMatchedObjects>0){
				matchInfoCollection.push_back(matchInfo);
			}

		}//end (j index) loop source in tile

	}//end (i index) loop tiles

	cout<<"== ISLAND MATCH RESULTS =="<<endl;
	cout<<"#matched/tot sources="<<nMatchedSources<<"/"<<nTotSources<<endl;
	for(size_t i=0;i<nMatchedSourcesPerCatalog.size();i++){
		cout<<"--> CAT"<<i+1<<": #match/tot="<<nMatchedSourcesPerCatalog[i]<<"/"<<nTotSources<<endl;
	}
	cout<<"source match multiplicity {";
	for(size_t i=0;i<nMatchedSourceCatalogMultiplicity.size()-1;i++){
		cout<<nMatchedSourceCatalogMultiplicity[i]<<",";
	}
	cout<<nMatchedSourceCatalogMultiplicity[nMatchedSourceCatalogMultiplicity.size()-1]<<"}"<<endl;
	cout<<"===================="<<endl;

	
	//Fill TTree
	for(size_t i=0;i<matchInfoCollection.size();i++)
	{
		SourcePars* sourcePars= matchInfoCollection[i]->spars;
		int sourceIndex= sourcePars->sourceIndex;
		int nestedSourceIndex= sourcePars->nestedSourceIndex;
		std::vector<AstroObject> objs;
		
		//Set TTree data
		matchedObjects.clear();
		fluxDiff.clear();
		matchedObjectMultiplicity.clear();
		nMatchedObjects= 0;
		sourceName= sourcePars->sname;
		sourcePosX= sourcePars->X0;
		sourcePosY= sourcePars->Y0;
		sourcePosX_wcs= sourcePars->X0_wcs;
		sourcePosY_wcs= sourcePars->Y0_wcs;
		sourceFluxDensity= sourcePars->fluxDensity;
		sourceFlux= sourcePars->flux;
		sourceFluxDensityErr= sourcePars->fluxDensityErr;
		sourceFluxErr= sourcePars->fluxErr;

		for(size_t j=0;j<(matchInfoCollection[i]->objs).size();j++)
		{
			AstroObject* obj= (matchInfoCollection[i]->objs)[j];
			std::string name= obj->name;
			int mult= matchedSourceNMatchMap[name];

			matchedObjects.push_back(obj);
			objs.push_back(*obj);

			nMatchedObjects++;
			double dflux= sourceFlux- obj->flux;
			fluxDiff.push_back(dflux);
			matchedObjectMultiplicity.push_back(mult);

		}//end loop matched objects

		//Fill TTree 
		matchedSourceInfo->Fill();

		//Set match object in source
		Source* source= 0;
		if(nestedSourceIndex==-1) {
			source= m_sources[sourceIndex];
		}
		else {
			source= m_sources[sourceIndex]->GetNestedSource(nestedSourceIndex);
		}
		if(source){
			for(size_t j=0;j<objs.size();j++){
				source->AddAstroObject(objs[j]);
			}
		}

	}//end loop matched sources

	return 0;

}//close FindSourceMatchesInTiles()


int FindSourceComponentMatchesInTiles()
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
					DEBUG_LOG("#"<<NSourceComponentPars<<" source components to be cross-matched against #"<<NObjectsInCatalog<<" objects (catalog #"<<s+1<<") in tile no. "<<i+1<<" ...");
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
				DEBUG_LOG("Update source component match counters...");
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
				DEBUG_LOG("Compute catalog component match multiplicity...");
			#endif
			int multiplicity= std::accumulate(sourceComponentMatchMultiplicityCounter.begin(), sourceComponentMatchMultiplicityCounter.end(), 0);
			nMatchedSourceComponentCatalogMultiplicity[multiplicity]++;
	
			//Fill Tree
			if(nMatchedObjects>0){
				matchInfoCollection.push_back(matchInfo);
			}

		}//end (j index) loop component pars in this tile

	}//end (i index) loop tiles


	cout<<"== COMPONENT MATCH RESULTS =="<<endl;
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
		int sourceIndex= componentPars->sourceIndex;
		int nestedSourceIndex= componentPars->nestedSourceIndex;
		int componentIndex= componentPars->fitComponentIndex;
		std::vector<AstroObject> objs;

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
			objs.push_back(*obj);

		}//end loop matched objects

		//Fill TTree 
		matchedSourceComponentInfo->Fill();

		
		//Set match object in source
		Source* source= 0;
		if(nestedSourceIndex==-1) {
			source= m_sources[sourceIndex];
		}
		else {
			source= m_sources[sourceIndex]->GetNestedSource(nestedSourceIndex);
		}
		if(source){
			for(size_t j=0;j<objs.size();j++){
				source->AddComponentAstroObject(componentIndex,objs[j]);
			}
		}
		

	}//end loop matched components

	return 0;

}//close FindSourceComponentMatchesInTiles()


bool HaveSourceIslandMatch(SourceMatchPars& smatchpars,SourcePars* spars,AstroObject* obj)
{
	//## Check data
	if(!spars || !obj){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given for source pars and/or catalog object!");
		#endif
		return false;
	}

	//## Check contours
	Contour* cont1= spars->contour;
	Contour* cont2= obj->GetContour(true);
	if(matchSourcesByOverlap)
	{
		if(!cont1 || !cont2){		
			#ifdef LOGGING_ENABLED
				WARN_LOG("One/both contour are nullptr!");
			#endif
			return false;
		}
		// Check if contour pars were computed
		if(!cont1->HasParameters || !cont2->HasParameters){
			#ifdef LOGGING_ENABLED
				WARN_LOG("One/both contours have no computed parameters!");
			#endif
			return false;
		}
	}//close if check contour

	//## Fill match parameters
	// - Centroid dstance
	double Xc_1= spars->X0_wcs;
	double Yc_1= spars->Y0_wcs;
	double Xc_2= obj->x;
	double Yc_2= obj->y;
	double posThr= matchPosThr/3600.;//convert threshold in deg (contour are given in deg)
	double posDist_Euclidean= MathUtils::GetEuclideanDist(Xc_1,Yc_1,Xc_2,Yc_2);
	double posDist_Haversine= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);
	double posDist= 0.;
	if(useWCSSimpleGeometry) posDist= posDist_Euclidean;
	else posDist= posDist_Haversine;
	
	//- Contour overlaps
	double contArea_1= 0;
	double contArea_2= 0;
	double overlapArea= 0;
	int overlapFlag;
	double overlapAreaFraction_1= 0;
	double overlapAreaFraction_2= 0;

	if(cont1 && cont2){
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
		
		status= MathUtils::ComputeContourOverlapArea(overlapArea,overlapFlag,cont1,cont2);
		if(status<0 && matchSourcesByOverlap){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute contour overlap area!");
			#endif
			return false;
		}
		overlapAreaFraction_1= overlapArea/contArea_1;
		overlapAreaFraction_2= overlapArea/contArea_2;

	}//close if contours

	// - Flux diff
	double flux_1= spars->flux;
	double flux_2= obj->flux;
	double fluxRelDiff= 0;
	if(flux_2>0){
		fluxRelDiff= (flux_1-flux_2)/flux_2;
	}
	
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
				DEBUG_LOG("NO POS MATCH: Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<")}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<")}, dist(arcsec)="<<posDist*3600.);
			#endif
			return false;
		}
		#ifdef LOGGING_ENABLED
			INFO_LOG("POS MATCH: Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<")}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<")}, dist(arcsec)="<<posDist*3600.);
		#endif
	}

	//## Compute contour overlap
	if(matchSourcesByOverlap){
		
		//Check cases
		//1) Disjoint contours
		//2) Contour 1 inside contour 2 (overlap=Area1)
		//3) Contour 2 inside contour 1 (overlap=Area2)
		
		if(overlapFlag==eCONT_NOT_OVERLAPPING){
			#ifdef LOGGING_ENABLED	
				DEBUG_LOG("NO CONTOUR OVERLAP: Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<"), area="<<contArea_1<<"}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<"), area="<<contArea_2<<"}, overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
			#endif
			return false;
		}

		
		else if(overlapFlag==eCONT1_INSIDE_CONT2){
			if(applySourceAreaRatioThr && overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT1 INSIDE CONT2 (OVERLAP BELOW THR): Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<"), area="<<contArea_1<<"}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<"), area="<<contArea_2<<"}, overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		
		else if(overlapFlag==eCONT2_INSIDE_CONT1){
			if(applySourceAreaRatioThr && overlapAreaFraction_1<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT2 INSIDE CONT1 (OVERLAP BELOW THR): Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<"), area="<<contArea_1<<"}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<"), area="<<contArea_2<<"}, overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}
		
		else if(overlapFlag==eCONT_OVERLAPPING){
			if(overlapAreaFraction_1<matchOverlapThr || overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT OVERLAP BELOW THR: Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<"), area="<<contArea_1<<"}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<"), area="<<contArea_2<<"}, overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("CONTOUR OVERLAP MATCH: Source {name="<<spars->sname<<", pos("<<Xc_1<<","<<Yc_1<<"), area="<<contArea_1<<"}, CatalogObject {name="<<obj->name<<", pos("<<Xc_2<<","<<Yc_2<<"), area="<<contArea_2<<"}, overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
		#endif
		
	}//close if

	//## Check flux rel difference (if requested)	
	if(matchSourcesByFlux){
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}
	

	return true;

}//close HaveSourceIslandMatch()




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
	double flux_1= cpars->flux;
	double flux_2= obj->flux;
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

	//Compute region contours and check that only one region is given. Multiple region not supported
	if(m_regions.size()==1){
		//Compute contour
		m_regionContour= m_regions[0]->GetContour(true);
		if(!m_regionContour){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to get region contour, will not be able to shuffle sources in region!");
			#endif
		}
	}//close if
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Shuffling not supported with more than 1 region");
		#endif
	}

	return 0;

}//close ReadRegions()


int ReadData(std::string filename)
{
	//Open files
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
	std::vector<SourcePars*> spars;
	std::vector<ComponentPars*> cpars;
	int catalogIndex= -1;
	int sourceIndex= 0;


	//Find region boundary vertex
	double regionXmin= 1.e+99;
	double regionXmax= -1.e+99;
	double regionYmin= 1.e+99;
	double regionYmax= -1.e+99;

	if(m_regionContour){
		std::vector<TVector2> bb= m_regionContour->BoundingBoxVertex;
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
		INFO_LOG("Found #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif

	for(int i=0;i<sourceTree->GetEntries();i++)
	{
		sourceTree->GetEntry(i);
		int morphId= aSource->MorphId;
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif


		//Select source by morph id?
		if(selectSourceByMorphId){
			bool skipSource= true;
			for(size_t j=0;j<stypes.size();j++){
				if( stypes[j]==-1 || morphId==stypes[j]) {
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

		//Compute source position in WCS coords (needed to compute source min/max coords)
		double x0_wcs= 0;
		double y0_wcs= 0;
		source->GetWCSPos(x0_wcs,y0_wcs,wcs,wcsType);

		//Shuffle source?
		//NB: Generate uniform in annulus (e.g. see https://stackoverflow.com/questions/9048095/create-random-number-within-an-annulus)
		double offsetX= 0;
		double offsetY= 0;
		if(shuffleSources){
			if(shuffleInRegion){
				double x0_wcs_rand= 0;
				double y0_wcs_rand= 0;
				if(GetRandPosInRegion(x0_wcs_rand,y0_wcs_rand,m_regionContour,regionXmin,regionXmax,regionYmin,regionYmax)<0){
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

			x0_wcs+= offsetX;
			y0_wcs+= offsetY;

		}//close if shuffleSources

		/*
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
		
		if(shuffleSources){
			X0_wcs+= offsetX;
			Y0_wcs+= offsetY;
		}
		*/


		//Find min & max coordinates
		if(x0_wcs>sourceXmax) sourceXmax= x0_wcs;
		if(x0_wcs<sourceXmin) sourceXmin= x0_wcs;
		if(y0_wcs>sourceYmax) sourceYmax= y0_wcs;
		if(y0_wcs<sourceYmin) sourceYmin= y0_wcs;	
		
		//Fill source pars
		std::vector<ComponentPars*> thisComponentPars;
		std::vector<SourcePars*> thisSourcePars;
		if(FillSourcePars(thisSourcePars,thisComponentPars,source,catalogIndex,sourceIndex,-1,wcs,wcsType,offsetX,offsetY)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill pars for source no. "<<i+1<<" (name="<<source->GetName()<<", collectionIndex="<<catalogIndex<<")!");
			#endif
			return -1;
		}

		//Add pars to pars
		cpars.insert(cpars.end(),thisComponentPars.begin(),thisComponentPars.end());
		spars.insert(spars.end(),thisSourcePars.begin(),thisSourcePars.end());

		//Add sources to list
		m_sources.push_back(source);
		sourceIndex++;

	}//end loop sources

	//Append source pars to main collection
	m_sourceComponentPars.insert(m_sourceComponentPars.end(),cpars.begin(),cpars.end());
	m_sourcePars.insert(m_sourcePars.end(),spars.begin(),spars.end());

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source min/max WCS coordinate range: (Xmin,Xmax)=("<<sourceXmin<<","<<sourceXmax<<"), (Ymin,Ymax)=("<<sourceYmin<<","<<sourceYmax<<")");
		INFO_LOG("#"<<m_sources.size()<<" sources to be cross-matched (#"<<m_sourceComponentPars.size()<<" source components) ...");
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


	//## Add source component data to tile
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
				DEBUG_LOG("Adding source component (name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<") to list...");
			#endif
			tileDataList[tileIndex]->AddSourceComponentPars(m_sourceComponentPars[i]);
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot find tile index or invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<") for source component (name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<"), check tile index calculation!");
			#endif
			continue;
		}
	}//end loop pars

	//Printing num source components in each tile
	for(size_t i=0;i<tileDataList.size();i++){
		long int nSourceComponentPars= tileDataList[i]->GetNSourceComponentPars();
		#ifdef LOGGING_ENABLED
			INFO_LOG("Tile no. "<<i+1<<": #"<<nSourceComponentPars<<" source components to be matched...");
		#endif
	}


	//## Add source island data to tile
	for(size_t i=0;i<m_sourcePars.size();i++)
	{
		double X0_wcs= m_sourcePars[i]->X0_wcs;
		double Y0_wcs= m_sourcePars[i]->Y0_wcs;
		std::string sourceName= m_sourcePars[i]->sname;

		long int tileIndex= MathUtils::FindGrid2DBin(
			X0_wcs,Y0_wcs,
			nTilesX,tileXmin,tileXmax,tileXstep,
			nTilesY,tileYmin,tileYmax,tileYstep
		);

		if(tileIndex>=0 && tileIndex<(long int)(tileDataList.size())){//Add to tile data
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding source (name="<<sourceName<<", pos("<<X0_wcs<<","<<Y0_wcs<<") to list...");
			#endif
			tileDataList[tileIndex]->AddSourcePars(m_sourcePars[i]);
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
			INFO_LOG("Tile no. "<<i+1<<": #"<<nSourcePars<<" source islands to be matched...");
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
	else if(catalogType==eGAIA_CATALOG){
		AstroObjectParser::ParseGaiaData(astroObjects,filename,'|');
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
	else if(catalogType==eMASH_CATALOG){
		AstroObjectParser::ParseMASHData(astroObjects,filename,'|');
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

int FillSourcePars(std::vector<SourcePars*>& spars,std::vector<ComponentPars*>& cpars,Source* aSource,int catalogIndex,int sourceIndex,int nestedSourceIndex,WCS* wcs,int coordSystem,double offsetX,double offsetY)
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
	bool useFWHM= true;
	bool convertToWCS= true;
	int pixOffset= 0;

	//Get 0th contour converted to WCS & relative centroid
	bool computeContourPars= true;
	Contour* contour= aSource->GetWCSContour(0,wcs,coordSystem,pixOffset,computeContourPars);
	if(!contour){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute WCS contour for source "<<sourceName<<"!");
		#endif
		return -1;
	}

	//Shuffle contour using offset
	if(shuffleSources){
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
		//Set fitted flux density
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= fitPars.GetNComponents();
		double fluxDensity= fitPars.GetFluxDensity();
		double fluxDensityErr= fitPars.GetFluxDensityErr();
		double flux= fluxDensity/beamArea;
		double fluxErr= fluxDensityErr/beamArea;
		sourcePars->fluxDensity= fluxDensity;
		sourcePars->fluxDensityErr= fluxDensityErr;
		sourcePars->flux= flux;
		sourcePars->fluxErr= fluxErr;

		//Get fitted pars & ellipse converted to WCS
		std::vector<TEllipse*> fitEllipses;
		if(aSource->GetFitEllipses(fitEllipses,useFWHM,convertToWCS,wcs,coordSystem,pixOffset,useWCSSimpleGeometry)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute WCS ellipse for fitted components of source "<<sourceName<<"!");
			#endif
			CodeUtils::DeletePtr<Contour>(contour);
			CodeUtils::DeletePtrCollection<SourcePars>(spars);
			CodeUtils::DeletePtrCollection<ComponentPars>(cpars);
			return -1;
		}

		//Check ellipses and pars size
		if(nComponents!=(int)(fitEllipses.size())){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Number of fit components shall be equal to fitted ellipses (this should not occur)!");
			#endif
			CodeUtils::DeletePtr<Contour>(contour);
			CodeUtils::DeletePtrCollection<SourcePars>(spars);
			CodeUtils::DeletePtrCollection<ComponentPars>(cpars);
			CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
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
			fluxDensity= fitPars.GetComponentFluxDensity(k);
			fluxDensityErr= fitPars.GetComponentFluxDensityErr(k);
			flux= fluxDensity/beamArea;
			fluxErr= fluxDensityErr/beamArea;

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
			cpars.push_back(thisComponentPars);

		}//end loop components
	}//close if has fit info

	//Append source pars to collection
	spars.push_back(sourcePars);
	
	//## Fill pars for nested sources (if any)
	if(hasNestedSources){
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			if(FillSourcePars(spars,cpars,nestedSources[j],catalogIndex,sourceIndex,j,wcs,coordSystem,offsetX,offsetY)<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
				#endif
				return -1;
			}
		}//end loop sources
	}//close if nested sources
	
	return 0;

}//close FillSourceComponentPars()


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
	//Fill Summary Tree
	if(matchSummaryTree) matchSummaryTree->Fill();

	//Save match TTree to file
	if(matchOutputFile && matchOutputFile->IsOpen())
	{		
		matchOutputFile->cd();
		if(matchedSourceInfo){
			matchedSourceInfo->Write();
		}
		if(matchedSourceComponentInfo){
			matchedSourceComponentInfo->Write();
		}
		if(matchSummaryTree){
			matchSummaryTree->Write();
		}
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
	for(size_t i=0;i<m_sources.size();i++){
		m_source= m_sources[i];
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

