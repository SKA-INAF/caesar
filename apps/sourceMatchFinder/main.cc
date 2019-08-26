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
#include <SourceCube.h>

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

//WCS headers (TO BE DEPRECATED)
//#include <wcs.h>

//BOOST headers
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "boost/graph/graph_traits.hpp"

typedef boost::adjacency_list<
	boost::vecS, 
	boost::vecS, 
	boost::undirectedS
> UndirectedGraph;//undirected graph


#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace boost;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input=[INPUT_FILE] \t Input file list name (ascii format) containing all catalog files in ROOT format to be cross-matched each other"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (root format) in which to store source match info"<<endl;
	
	cout<<"-x, --xmin=[XMIN] \t Grid RA min to be searched (in deg)"<<endl; 
	cout<<"-X, --xmax=[XMAX] \t Grid RA max to be searched (in deg)"<<endl;
 	cout<<"-w, --xstep=[XSTEP] \t Grid RA step (in arcmin)"<<endl;
 	cout<<"-y, --ymin=[XMIN] \t Grid DEC min to be searched (in deg)"<<endl; 
	cout<<"-Y, --ymax=[XMAX] \t Grid DEC max to be searched (in deg)"<<endl;
 	cout<<"-k, --ystep=[XSTEP] \t Grid DEC step (in arcmin)"<<endl;

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
	cout<<"-n, --minSourceMatchClusterSize=[MIN_SOURCE_MATCH_CLUSTER_SIZE] \t Minimum number of matched sources to be retained (e.g. skip source match cluster below this threshold) (default=2)"<<endl;
	
	cout<<"-f, --filterByType \t Consider only true sources with given type when searching the match (default=no)"<<endl;
	cout<<"-s, --selectedType=[TYPE] \t True source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED) (default=-1)"<<endl;
	cout<<"-F, --filterBySimType \t Consider only true sources with given sim type when searching the match (default=no)"<<endl;
	cout<<"-S, --selectedSimType=[TYPE] \t True source sim types to be crossmatched (eRingLike=1,eBubbleLike=2,eEllipseLike=3,eDiskLike=4,eBlobLike=5) (default=-1)"<<endl;
	cout<<"-j, --no-compactSourceCorrelation \t Disable correlation search for compact sources (default=enabled)"<<endl;
	cout<<"-J, --no-extendedSourceCorrelation \t Disable correlation search for extended sources (default=enabled)"<<endl;
	//cout<<"-c, --correctFlux \t Correct rec integrated flux by beam area (default=no correction)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", optional_argument, 0, 'o' },
	{ "xmin", required_argument, 0, 'x' },
	{ "xmax", required_argument, 0, 'X' },
	{ "xstep", required_argument, 0, 'w' },
	{ "ymin", required_argument, 0, 'y' },
	{ "ymax", required_argument, 0, 'Y' },
	{ "ystep", required_argument, 0, 'k' },

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

	//{ "correctFlux", no_argument, 0, 'c'},
	{ "filterByType", no_argument, 0, 'f'},
	{ "filterBySimType", no_argument, 0, 'F'},
	{ "selectedType", required_argument, 0, 's'},	
	{ "selectedSimType", required_argument, 0, 'S'},
	{ "no-compactSourceCorrelation", no_argument, 0, 'j'},
	{ "no-extendedSourceCorrelation", no_argument, 0, 'J'},	
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
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
int minSourceMatchClusterSize= 2;
bool correctFlux= false;
bool selectTrueSourceByType= false;//default=all true sources searched 
//int selectedTrueSourceType= -1;
std::vector<int> stypes;
std::vector<int> ssimtypes;
bool selectTrueSourceBySimType= false;//default=all true sources searched 
int selectedTrueSourceSimType= -1;
bool enableCompactSourceCorrelation= true;
bool enableExtendedSourceCorrelation= true;
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

//Globar vars
TFile* outputFile= 0;
std::string outputFileName= "MatchOutput.root";
TTree* matchedSourceInfo= 0;
SourceCube* sourceCube= 0;
std::vector<SourceCube*> sourceCubes;

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
};


//====================================
//     SourcePars struct
//====================================
struct SourcePars {
	//SourcePars(){}
	SourcePars(int _catalogIndex,int _sourceIndex,int _nestedSourceIndex=-1)
		: catalogIndex(_catalogIndex), sourceIndex(_sourceIndex), nestedSourceIndex(_nestedSourceIndex)
	{
		S= 0;
		fluxDensity= 0;
		contour= 0;
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
	Contour* contour;
	double S;//pixel-integrated flux
	double fluxDensity;//fitted flux density
	std::vector<ComponentPars*> componentPars;
};

//====================================
//     TileData struct
//====================================
struct TileData {

	//Constructor
	TileData(float _xmin,float _xmax,float _ymin,float _ymax)
		: xmin(_xmin),xmax(_xmax),ymin(_ymin),ymax(_ymax)
	{
		neighborTileIndexes.clear();
		sourcePars.clear();
	}

	//Destructor 
	~TileData(){
		CodeUtils::DeletePtrCollection<SourcePars>(sourcePars);
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

	//Tile coordinates
	float xmin;
	float xmax;
	float ymin;
	float ymax;

	//List of neighbors tiles
	std::vector<size_t> neighborTileIndexes;

	//Source pars for this tile
	std::vector<SourcePars*> sourcePars;

};//close TileData()

std::vector<TileData*> tileDataList;
std::vector<std::vector<Source*>> sources;
std::vector<SourcePars*> source_pars;

struct clique_store {
	//Constructor
	clique_store(std::vector<std::vector<size_t>>& _clusterIds)
		: clusterIds(_clusterIds)
	{}

  template <typename Clique, typename GraphType>
  void clique(const Clique& c, const GraphType& g)
  {
		//Create a new cluster
		std::vector<size_t> newClusterIds;

  	// Iterate over the clique and print each vertex within it
    typename Clique::const_iterator it, end = c.end();
    for(it = c.begin(); it != end; ++it) {
			newClusterIds.push_back(*it);
			//cout << *it << " ";
    }
    //cout << endl;
	
		clusterIds.push_back(newClusterIds);
  }//close clique()
  
	std::vector<std::vector<size_t>>& clusterIds;
};


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int InitGrid(float xmin,float xmax,float xstep,float ymin,float ymax,float ystep);
int FindSourceMatchesInTiles();
int ReadData();
int ReadSourceData(std::string filename,int collectionIndex);
int FillSourcePars(std::vector<SourcePars*>& pars,Source* aSource,int collectionIndex,int sourceIndex,int nestedSourceIndex=-1,WCS* wcs=0,int coordSystem=0);
bool HaveSourceComponentMatch(ComponentPars* pars1,ComponentPars* pars2);
bool HaveSourceIslandMatch(SourcePars* pars1,SourcePars* pars2);
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
	//== Read data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading source data from given files...");
	#endif
	if(ReadData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of source data failed!");
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

	while((c = getopt_long(argc, argv, "hi:o:v:mM:POLlA:T:t:e:pba:d:n:cfFs:S:jJx:X:w:y:Y:k:",options_tab, &option_index)) != -1) {
    
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
			case 'o':	
			{
				outputFileName= std::string(optarg);	
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
			case 'n':
			{
				minSourceMatchClusterSize= atoi(optarg);
				break;
			}
			/*
			case 'c':
			{
				correctFlux= true;
				break;
			}
			*/
	
			case 'f':
			{
				selectTrueSourceByType= true;
				break;
			}
			/*
			case 's':
			{
				selectedTrueSourceType= atoi(optarg);
				break;
			}
			*/
			case 's':	
			{
				int stype= atoi(optarg);
				stypes.push_back(stype);	
				break;	
			}	
			case 'F':
			{
				selectTrueSourceBySimType= true;
				break;
			}
			case 'S':
			{
				//selectedTrueSourceSimType= atoi(optarg);
				int ssimtype= atoi(optarg);
				ssimtypes.push_back(ssimtype);	
				break;
			}	
			case 'j':
			{
				enableCompactSourceCorrelation= false;
				break;
			}	
			case 'J':
			{
				enableExtendedSourceCorrelation= false;
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
			case 'x':
			{
				tileXmin= atof(optarg);
				break;
			}
			case 'X':
			{
				tileXmax= atof(optarg);
				break;
			}
			case 'w':
			{
				tileXstep= atof(optarg)/60.;//convert in deg
				break;
			}
			case 'y':
			{
				tileYmin= atof(optarg);
				break;
			}
			case 'Y':
			{
				tileYmax= atof(optarg);
				break;
			}
			case 'k':
			{
				tileYstep= atof(optarg)/60.;//convert in deg
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
	if( tileXmin==tileXmax || (tileXmax==0 && tileXmin==0) || tileXstep==0){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Invalid or missing tile x range options!");
		#endif
		return -1;
	}
	if( tileYmin==tileYmax || (tileYmax==0 && tileYmin==0) || tileYstep==0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid or missing tile y range options!");
		#endif
		return -1;
	}	

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
	#endif

	return 0;

}//close ParseOptions()

int Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//Init source cube TTree
	if(!matchedSourceInfo) matchedSourceInfo= new TTree("SourceMatchInfo","SourceMatchInfo");
	matchedSourceInfo->Branch("SourceCube",&sourceCube);
	sourceCubes.clear();

	//Initialize tile data grid
	if(InitGrid(tileXmin,tileXmax,tileXstep,tileYmin,tileYmax,tileYstep)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to initialize tile data grid!");
		#endif
		return -1;
	}

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
	//Init looked tile links
	std::map<std::pair<size_t,size_t>,bool> tileLinkSearched;
	for(size_t i=0;i<tileDataList.size()-1;i++){
		tileLinkSearched.insert( std::make_pair(std::make_pair(i,i),false) );
		for(size_t j=i+1;j<tileDataList.size();j++){
			tileLinkSearched.insert( std::make_pair(std::make_pair(i,j),false) );
			tileLinkSearched.insert( std::make_pair(std::make_pair(j,i),false) );
		}//end loop tile next
	}//end loop 

	//Init source cubes
	SourceCube* aSourceCube= 0;
	long int scubeCounter= 0;

	//Define match component pars
	struct MatchComponentPars {
		MatchComponentPars(){
			componentPars= 0;
			sourceClusterIndex= -1;
		}
		~MatchComponentPars(){

		}
		long int sourceParIndex;
		size_t sourceClusterIndex;
		size_t componentIndex;
		ComponentPars* componentPars;
	};

	//Loop over tiles and find matches over neighbors	
	for(size_t i=0;i<tileDataList.size();i++){
		size_t index= i;
		long int NSourcePars= tileDataList[i]->GetNSourcePars();
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<NSourcePars<<" sources to be cross-matched...");
		#endif

		//Fill list of source pars to be cross-matched 
		//NB: Include only this tile and neighbor tiles not searched before
		std::vector<SourcePars*> sourceParsToBeMatched;
		sourceParsToBeMatched.insert(sourceParsToBeMatched.end(),(tileDataList[i]->sourcePars).begin(),(tileDataList[i]->sourcePars).end());

		std::vector<size_t> neighbors= tileDataList[i]->neighborTileIndexes;

		//=== DEBUG
		std::stringstream sstream;
		sstream<<"Tile no. "<<i+1<<": neighbors{";
		for(size_t j=0;j<neighbors.size();j++){
			sstream<<neighbors[j]<<",";
		}
		sstream<<"}";
		#ifdef LOGGING_ENABLED	
			DEBUG_LOG(sstream.str());
		#endif
		//====DEBUG
		

		for(size_t j=0;j<neighbors.size();j++){
			size_t neighborIndex= neighbors[j];
			std::pair<size_t,size_t> indexLink(index,neighborIndex);
			bool linkSearched= tileLinkSearched[indexLink];
			if(linkSearched){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Skip tile link ("<<index<<"-"<<neighborIndex<<") as already searched before...");
				#endif
				continue;
			}
			sourceParsToBeMatched.insert(sourceParsToBeMatched.end(),(tileDataList[neighborIndex]->sourcePars).begin(),(tileDataList[neighborIndex]->sourcePars).end());

			//Set tile link as visited
			tileLinkSearched[indexLink]= true;

			//Set inverse link as visited
			std::pair<size_t,size_t> indexLink_inv(neighborIndex,index);	
			tileLinkSearched[indexLink_inv]= true;
		}//end loop neighbors


		//Check if there are sources to be crossmatched
		if(sourceParsToBeMatched.empty()){
			#ifdef LOGGING_ENABLED	
				DEBUG_LOG("No sources to be crossmatched for tile no. "<<i+1<<" and its neighbors, skip to next tile...");
			#endif
			continue;
		}

		//Init graph
		UndirectedGraph undirGraph(sourceParsToBeMatched.size());
		std::vector<std::vector<size_t>> clusterIds;
  	clique_store vis(clusterIds);

		//Loop over source parameters and build adjacency list
		for(size_t k=0;k<sourceParsToBeMatched.size()-1;k++){
			#ifdef LOGGING_ENABLED
				if(k%100==0) INFO_LOG("#"<<k+1<<"/"<<sourceParsToBeMatched.size()<<" source scanned for adjacency for tile no. "<<i+1<<" ...");
			#endif

			for(size_t l=k+1;l<sourceParsToBeMatched.size();l++){
				bool hasMatch= HaveSourceIslandMatch(sourceParsToBeMatched[k],sourceParsToBeMatched[l]);
				if(hasMatch){//Set link in graph between two sources
					boost::add_edge(k,l,undirGraph);
				}
			}//end loop next source pars
		}//end loop source pars

		// Use the Bron-Kerbosch algorithm to find all cliques, printing them as they are found
		#ifdef LOGGING_ENABLED
			INFO_LOG("Using the Bron-Kerbosch algorithm to find all source matches (min clust size="<<minSourceMatchClusterSize<<") ...");
		#endif
		std::size_t minCliqueSize= minSourceMatchClusterSize;
		bron_kerbosch_all_cliques(undirGraph, vis, minCliqueSize);

		//Get cluster ids from clique and print them
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<clusterIds.size()<<" source match clusters found with multiplicity >"<<minCliqueSize<<" in tile no. "<<i+1<<" and its neighbors ...");
		#endif
		std::stringstream ss;
	

		for(size_t k=0;k<clusterIds.size();k++){
			//Print source cluster founds
			ss.str("");
			ss<<"SOURCE CLUST MATCH no. "<<k+1<<": {";
			for(size_t l=0;l<clusterIds[k].size();l++){
				size_t index= clusterIds[k][l];
				std::string sname= sourceParsToBeMatched[index]->sname;
				Contour* contour= sourceParsToBeMatched[index]->contour;
				int catalogIndex= sourceParsToBeMatched[index]->catalogIndex;
				double posX= 0;
				double posY= 0;
				if(contour){
					posX= (contour->Centroid).X();
					posY= (contour->Centroid).Y();
				}
				if(l==clusterIds[k].size()-1) ss<<"(CAT"<<catalogIndex<<", index="<<index<<", pos("<<posX<<","<<posY<<"), sname="<<sname<<")}";
				else ss<<"(CAT"<<catalogIndex<<", index="<<index<<", pos("<<posX<<","<<posY<<"), sname="<<sname<<"), ";
			}//end loop cluster members in this cluster

			#ifdef LOGGING_ENABLED
				INFO_LOG(ss.str());		
			#endif

			//Fill list of components for graph search
			std::vector<MatchComponentPars*> matchComponentParList;
			MatchComponentPars* matchPars= 0;
			for(size_t l=0;l<clusterIds[k].size();l++){
				size_t sourceParIndex= clusterIds[k][l];
				for(size_t t=0;t<(sourceParsToBeMatched[sourceParIndex]->componentPars).size();t++){
					matchPars= new MatchComponentPars;	
					matchPars->sourceParIndex= sourceParIndex;
					matchPars->sourceClusterIndex= l;
					matchPars->componentIndex= t;
					matchPars->componentPars= (sourceParsToBeMatched[sourceParIndex]->componentPars)[t];
					matchComponentParList.push_back(matchPars);	
				}//end loop components for this source
			}//end loop source contour cluster

			//Search for match in components for this contour cluster
			#ifdef LOGGING_ENABLED
				INFO_LOG("Searching for cluster components in source contour cluster no. "<<k+1<<" (#"<<clusterIds[k].size()<<" sources present) for tile no. "<<i+1<<" ...");
			#endif
			UndirectedGraph undirGraph_component(matchComponentParList.size());
			std::vector<std::vector<size_t>> clusterComponentIds;
 	 		clique_store vis_component(clusterComponentIds);

			for(size_t l=0;l<matchComponentParList.size()-1;l++){
				for(size_t t=l+1;t<matchComponentParList.size();t++){
					bool hasComponentMatch= HaveSourceComponentMatch(matchComponentParList[l]->componentPars,matchComponentParList[t]->componentPars);
					if(hasComponentMatch){//Set link in graph between two sources
						boost::add_edge(l,t,undirGraph_component);
					}
				}//end loop next component pars
			}//end loop component pars
			
			#ifdef LOGGING_ENABLED
				INFO_LOG("Using the Bron-Kerbosch algorithm to find all source component matches (min clust size="<<minSourceMatchClusterSize<<") ...");
			#endif
			std::size_t minCliqueSize_component= minSourceMatchClusterSize;
			bron_kerbosch_all_cliques(undirGraph_component, vis_component, minCliqueSize_component);

			//Print component cluster founds
			ss.str("");
			for(size_t l=0;l<clusterComponentIds.size();l++){
				ss<<"COMPONENT CLUST MATCH no. "<<l+1<<": {";
				for(size_t t=0;t<clusterComponentIds[l].size();t++){	
					size_t index= clusterComponentIds[l][t];
					int sourceParIndex= matchComponentParList[index]->sourceParIndex;
					ComponentPars* cpars= matchComponentParList[index]->componentPars;
					int catalogIndex= cpars->catalogIndex;
					std::string sname= cpars->sname;
					TEllipse* ellipse= cpars->fitEllipse;
					double posX= 0;
					double posY= 0;
					if(ellipse){
						posX= ellipse->GetX1();
						posY= ellipse->GetY1();
					}
				
					if(t==clusterComponentIds[l].size()-1) ss<<"(index="<<index<<", pos("<<posX<<","<<posY<<", sname="<<sname<<")}";
					else ss<<"(index="<<index<<", pos("<<posX<<","<<posY<<", sname="<<sname<<"), "; 
				}
				#ifdef LOGGING_ENABLED
					INFO_LOG(ss.str());		
				#endif
			}//end loop cluster components

			//Check if any cluster component is found
			if(clusterComponentIds.empty()){
				#ifdef LOGGING_ENABLED
					INFO_LOG("No cluster found with components, skip to next cluster...");
				#endif
				continue;
			}

			//Check if a number of matching components less than available sources is found
			//if(){
			//
			//}

			//Store source cube found
			#ifdef LOGGING_ENABLED
				INFO_LOG("Filling source cube info...");
			#endif
			TString scubeName= Form("SCube%d",scubeCounter);
			aSourceCube= new SourceCube(std::string(scubeName));
			scubeCounter++;

			for(size_t l=0;l<clusterIds[k].size();l++){
				//Retrieve source to be stored from catalog collection
				size_t index= clusterIds[k][l];
				SourcePars* spars= sourceParsToBeMatched[index];
				int catalogIndex= spars->catalogIndex;
				int sourceIndex= spars->sourceIndex;
				int nestedSourceIndex= spars->nestedSourceIndex;
	
				Source* sourceToBeStored= 0;
				if(nestedSourceIndex==-1){
					sourceToBeStored= sources[catalogIndex][sourceIndex];
				}
				else{
					sourceToBeStored= sources[catalogIndex][sourceIndex]->GetNestedSource(nestedSourceIndex);
				}

				if(!sourceToBeStored){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Cannot retrieve source to be stored (catalogIndex="<<catalogIndex<<", sindex="<<sourceIndex<<", nestedSourceIndex="<<nestedSourceIndex<<")");	
					#endif
					continue;
				}

				//Add source to cube
				aSourceCube->AddSource(sourceToBeStored);
				
			}//end loop cluster members in this cluster
			
			//Add components to cube
			for(size_t l=0;l<clusterComponentIds.size();l++){
				if(clusterComponentIds[l].size()>0) aSourceCube->AddComponentMatch();
				for(size_t t=0;t<clusterComponentIds[l].size();t++){	
					size_t index= clusterComponentIds[l][t];
					MatchComponentPars* matchCPars= matchComponentParList[index];
					size_t sindex= matchCPars->sourceClusterIndex;
					size_t cindex= matchCPars->componentIndex;
					if(aSourceCube->AddIndexToComponent(l,sindex,cindex)<0){	
						#ifdef LOGGING_ENABLED
							ERROR_LOG("Failed to add component index to source cube!");
						#endif
						continue;
					}
				}//end loop components in cluster
			}//end loop component clusters
			
			//Add cube to global collection
			sourceCubes.push_back(aSourceCube);
	
		}//end loop clusters


	}//end loop tiles

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sourceCubes.size()<<" source cubes found and added to collection...");
	#endif

	return 0;

}//close FindSourceMatchesInTiles()

/*
int FindSourceMatches()
{
	INFO_LOG("#"<<source_pars.size()<<" sources to be cross-matched...");

	//Build undirected graph and initialize cluster store
	UndirectedGraph undirGraph(source_pars.size());
	std::vector<std::vector<size_t>> clusterIds;
  clique_store vis(clusterIds);

	//Loop over source parameters and build adjacency list
	for(size_t i=0;i<source_pars.size()-1;i++){
		if(i%100==0) INFO_LOG("#"<<i+1<<"/"<<source_pars.size()<<" source scanned for adjacency...");

		for(size_t j=i+1;j<source_pars.size();j++){
			bool hasMatch= HaveSourceIslandMatch(source_pars[i],source_pars[j]);
			if(hasMatch){//Set link in graph between two sources
				boost::add_edge(i,j,undirGraph);
			}
		}//end loop next source pars
	}//end loop source pars

	
	// Use the Bron-Kerbosch algorithm to find all cliques, printing them as they are found
	INFO_LOG("Using the Bron-Kerbosch algorithm to find all source matches (min clust size="<<minSourceMatchClusterSize<<") ...");
	std::size_t minCliqueSize= minSourceMatchClusterSize;
	bron_kerbosch_all_cliques(undirGraph, vis, minCliqueSize);

	//Get cluster ids from clique and print them
	INFO_LOG("#"<<clusterIds.size()<<" source match clusters found with multiplicity >"<<minCliqueSize<<" ...");
	std::stringstream ss;

	struct MatchComponentPars {
		long int sourceParIndex;
		ComponentPars* componentPars;
	};

	for(size_t i=0;i<clusterIds.size();i++){
		//Print source cluster founds
		ss.str("");
		ss<<"SOURCE CLUST MATCH no. "<<i+1<<": {";
		for(size_t j=0;j<clusterIds[i].size();j++){
			int index= clusterIds[i][j];
			std::string sname= source_pars[index]->sname;
			Contour* contour= source_pars[index]->contour;
			double posX= 0;
			double posY= 0;
			if(contour){
				posX= (contour->Centroid).X();
				posY= (contour->Centroid).Y();
			}
			if(j==clusterIds[i].size()-1) ss<<"(index="<<index<<", pos("<<posX<<","<<posY<<"), sname="<<sname<<")}";
			else ss<<"(index="<<index<<", pos("<<posX<<","<<posY<<"), sname="<<sname<<"), ";
		}
		INFO_LOG(ss.str());		

		//Fill list of components for graph search
		std::vector<MatchComponentPars> matchComponentParList;
		for(size_t j=0;j<clusterIds[i].size();j++){
			int sourceParIndex= clusterIds[i][j];
			for(size_t k=0;k<(source_pars[sourceParIndex]->componentPars).size();k++){
				MatchComponentPars matchPars;	
				matchPars.sourceParIndex= sourceParIndex;
				matchPars.componentPars= (source_pars[sourceParIndex]->componentPars)[k];
				matchComponentParList.push_back(matchPars);	
			}//end loop components for this source
		}//end loop source contour cluster

		//Search for match in components for this contour cluster
		INFO_LOG("Searching for cluster components in source contour cluster no. "<<i+1<<" (#"<<clusterIds[i].size()<<" sources present)...");
		UndirectedGraph undirGraph_component(matchComponentParList.size());
		std::vector<std::vector<size_t>> clusterComponentIds;
 	 	clique_store vis_component(clusterComponentIds);

		for(size_t j=0;j<matchComponentParList.size()-1;j++){
			for(size_t k=j+1;k<matchComponentParList.size();k++){
				bool hasComponentMatch= HaveSourceComponentMatch(matchComponentParList[j].componentPars,matchComponentParList[k].componentPars);
				if(hasComponentMatch){//Set link in graph between two sources
					boost::add_edge(j,k,undirGraph_component);
				}
			}//end loop next component pars
		}//end loop component pars
		
		INFO_LOG("Using the Bron-Kerbosch algorithm to find all source component matches (min clust size="<<minSourceMatchClusterSize<<") ...");
		std::size_t minCliqueSize_component= minSourceMatchClusterSize;
		bron_kerbosch_all_cliques(undirGraph_component, vis_component, minCliqueSize_component);

		//Print component cluster founds
		ss.str("");
		for(size_t j=0;j<clusterComponentIds.size();j++){
			ss<<"COMPONENT CLUST MATCH no. "<<j+1<<": {";
			for(size_t k=0;k<clusterComponentIds[j].size()-1;k++){	
				int index= clusterComponentIds[j][k];
				int sourceParIndex= matchComponentParList[index].sourceParIndex;
				ComponentPars* cpars= matchComponentParList[index].componentPars;
				int catalogIndex= cpars->catalogIndex;
				std::string sname= cpars->sname;
				TEllipse* ellipse= cpars->fitEllipse;
				double posX= 0;
				double posY= 0;
				if(ellipse){
					posX= ellipse->GetX1();
					posY= ellipse->GetY1();
				}
				
				if(k==clusterComponentIds[j].size()-1) ss<<"(index="<<index<<", pos("<<posX<<","<<posY<<", sname="<<sname<<")}";
				else ss<<"(index="<<index<<", pos("<<posX<<","<<posY<<", sname="<<sname<<"), "; 
			}
			INFO_LOG(ss.str());		

		}//end loop cluster components

		//Store source cube found
		//...
		//...

	}//end loop clusters


	return 0;

}//close FindSourceMatches()
*/

bool HaveSourceIslandMatch(SourcePars* pars1,SourcePars* pars2)
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

	//## Compare contour centroids 
	TVector2 centroid_1= cont1->Centroid;
	TVector2 centroid_2= cont2->Centroid;
	double posX_1= centroid_1.X();
	double posX_2= centroid_2.X();
	double posY_1= centroid_1.Y();
	double posY_2= centroid_2.Y();
	double posThr= matchPosThr/3600.;//convert threshold in deg (contour are given in deg)

	if(matchSourcesByPos){	
		double posDist= 0;
		if(useWCSSimpleGeometry) posDist= MathUtils::GetEuclideanDist(posX_1,posY_1,posX_2,posY_2);
		else posDist= AstroUtils::GetWCSPointDist_Haversine(posX_1,posY_1,posX_2,posY_2);

		if(fabs(posDist)>posThr) {	
			#ifdef LOGGING_ENABLED
				INFO_LOG("NO POS MATCH: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<"), dist(arcsec)="<<posDist*3600.);
			#endif
			return false;
		}
		#ifdef LOGGING_ENABLED
			INFO_LOG("POS MATCH: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<"), dist(arcsec)="<<posDist*3600.);
		#endif
	}

	//## Compute contour overlap
	if(matchSourcesByOverlap){
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

		//Check cases
		//1) Disjoint contours
		//2) Contour 1 inside contour 2 (overlap=Area1)
		//3) Contour 2 inside contour 1 (overlap=Area2)
		//Compute relative overlap areas
		double overlapAreaFraction_1= overlapArea/contArea_1;
		double overlapAreaFraction_2= overlapArea/contArea_2;
		if(overlapFlag==eCONT_NOT_OVERLAPPING){
			#ifdef LOGGING_ENABLED	
				INFO_LOG("NO CONTOUR OVERLAP: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
			#endif
			return false;
		}

		else if(overlapFlag==eCONT1_INSIDE_CONT2){
			if(applySourceAreaRatioThr && overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					INFO_LOG("CONT1 INSIDE CONT2 (OVERLAP BELOW THR): Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		else if(overlapFlag==eCONT2_INSIDE_CONT1){
			if(applySourceAreaRatioThr && overlapAreaFraction_1<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					INFO_LOG("CONT2 INSIDE CONT1 (OVERLAP BELOW THR): Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}

		else if(overlapFlag==eCONT_OVERLAPPING){
			if(overlapAreaFraction_1<matchOverlapThr || overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					INFO_LOG("CONT OVERLAP BELOW THR: Source (CAT"<<pars1->catalogIndex<<", name="<<pars1->sname<<", pos="<<posX_1<<","<<posY_1<<", area="<<contArea_1<<"), Source (CAT"<<pars2->catalogIndex<<", name="<<pars2->sname<<", pos("<<posX_2<<","<<posY_2<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
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
		double flux_1= pars1->fluxDensity;
		double flux_2= pars2->fluxDensity;
		double fluxRelDiff= (flux_1-flux_2)/flux_2;
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}
	

	return true;

}//close HaveSourceIslandMatch()


bool HaveSourceComponentMatch(ComponentPars* pars1,ComponentPars* pars2)
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
	double posDist= 0.;
	if(useWCSSimpleGeometry) posDist= MathUtils::GetEuclideanDist(Xc_1,Yc_1,Xc_2,Yc_2);
	else posDist= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);

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
	
	

	if(matchSourceComponentsByOverlap)
	{
		//- No overlap: Disjoint ellipses
		if(rtn==EllipseUtils::DISJOINT_ELLIPSES){
			#ifdef LOGGING_ENABLED
				INFO_LOG("DISJOINT ELLIPSES: Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
			#endif
			return false;
		}

		//- Full overlap: Ellipse 1 inside ellipse 2 (overlap=Area1)
		else if(rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 ){
			//Check area overlap ratio
			if(applySourceComponentAreaRatioThr && overlapAreaFraction_2<compMatchOverlapThr){
				#ifdef LOGGING_ENABLED
					INFO_LOG("ELLIPSE1 INSIDE ELLIPSE2 (OVERLAP BELOW THR): overlapArea2="<<overlapAreaFraction_2<<"<"<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			//Check small pos match
			if(matchSourceComponentsByPos && !hasPosSmallMatch){
				#ifdef LOGGING_ENABLED
					INFO_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}
		}

		//- Full overlap: Ellipse 2 inside ellipse 1 (overlap=Area2)
		else if(rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1){
			//Check area overlap ratio
			if(applySourceComponentAreaRatioThr && overlapAreaFraction_1<compMatchOverlapThr){
				#ifdef LOGGING_ENABLED
					INFO_LOG("ELLIPSE2 INSIDE ELLIPSE1 (OVERLAP BELOW THR): overlapArea1="<<overlapAreaFraction_1<<"<"<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			//Check small pos match
			if(matchSourceComponentsByPos && !hasPosSmallMatch){
				#ifdef LOGGING_ENABLED
					INFO_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}
		}

		//- Partial overlap
		else {
			//- Large overlap
			if(hasOverlapLargeMatch){	
				#ifdef LOGGING_ENABLED
					INFO_LOG("ELLIPSE LARGE OVERLAP: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				
				//Check small pos match
				if(matchSourceComponentsByPos && !hasPosSmallMatch){
					#ifdef LOGGING_ENABLED
						INFO_LOG("NO ELLIPSE POS SMALL MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posHighThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
					#endif
					return false;
				}
			}//close if

			//- Small overlap
			else if(hasOverlapSmallMatch){
				#ifdef LOGGING_ENABLED
					INFO_LOG("ELLIPSE SMALL OVERLAP: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapLowThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapLowThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif

				//Check large pos match
				if(matchSourceComponentsByPos && !hasPosLargeMatch){
					#ifdef LOGGING_ENABLED
						INFO_LOG("NO ELLIPSE POS LARGE MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
					#endif
					return false;
				}
			}//close else if
			
			//- Partial overlap not sufficient
			else {
				#ifdef LOGGING_ENABLED
					INFO_LOG("ELLIPSE OVERLAP BELOW SMALL & HIGH THR: overlapArea1="<<overlapAreaFraction_1<<">="<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">="<<compMatchOverlapThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
				#endif
				return false;
			}

			#ifdef LOGGING_ENABLED
				INFO_LOG("ELLIPSE OVERLAP MATCH: overlapArea1="<<overlapAreaFraction_1<<", overlapArea2="<<overlapAreaFraction_2<<", compMatchOverlapThr="<<compMatchOverlapThr<<", compMatchOverlapLowThr="<<compMatchOverlapLowThr<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}, rtn="<<rtn);
			#endif
		}//close else

	}//close if check overlap
	else{
		//## Check ellipse centroid sky distance (if requested)
		if(matchSourceComponentsByPos && !hasPosLargeMatch){
			#ifdef LOGGING_ENABLED
				INFO_LOG("NO ELLIPSE POS MATCH: posDist(arcsec)="<<posDist*3600<<">"<<posThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}");
			#endif
			return false;
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("ELLIPSE POS MATCH: posDist(arcsec)="<<posDist*3600<<"<"<<posThr*3600<<", Source1 {name="<<pars1->sname<<", ellipse("<<Xc_1<<","<<Yc_1<<","<<R1_1<<","<<R2_1<<","<<Theta_1<<")}, Source2 {name="<<pars2->sname<<", ellipse("<<Xc_2<<","<<Yc_2<<","<<R1_2<<","<<R2_2<<","<<Theta_2<<")}");
		#endif
	}


	/*
	if(matchSourcesByPos){
		double Xc_1= ellipse1->GetX1();
		double Yc_1= ellipse1->GetY1();
		double Xc_2= ellipse2->GetX1();
		double Yc_2= ellipse2->GetY1();
		double posDist= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);
		if(fabs(posDist)>posThr) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Ellipse pos diff above thr ("<<posDist<<">"<<posThr<<", pos1("<<Xc_1<<","<<Yc_1<<"), pos2("<<Xc_2<<","<<Yc_2<<"))");
			#endif
			return false;
		}
	}

	//## Check ellipse overlap (if requested)
	if(matchSourcesByOverlap){
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

		//Check cases
		//1) Disjoint ellipses
		//2) Ellipse 1 inside ellipse 2 (overlap=Area1)
		//3) Ellipse 2 inside ellipse 1 (overlap=Area2)
		//Compute relative overlap areas
		double overlapAreaFraction_1= ellipseOverlapArea/ellipseArea_1;
		double overlapAreaFraction_2= ellipseOverlapArea/ellipseArea_2;
		if(rtn==EllipseUtils::DISJOINT_ELLIPSES){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Disjoint ellipses");
			#endif
			return false;
		}
		else if(rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 && overlapAreaFraction_2<matchOverlapThr){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Ellipse 1 inside 2 below overlap thr ("<<overlapAreaFraction_2<<"<"<<matchOverlapThr<<")");
			#endif
			return false;
		}
		else if(rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1 && overlapAreaFraction_1<matchOverlapThr){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Ellipse 2 inside 1 below overlap thr ("<<overlapAreaFraction_1<<"<"<<matchOverlapThr<<")");
			#endif
			return false;
		}
		else {
			if(overlapAreaFraction_1<matchOverlapThr || overlapAreaFraction_2<matchOverlapThr){	
				#ifdef LOGGING_ENABLED
					INFO_LOG("Ellipse overlap below thr ("<<overlapAreaFraction_1<<"<"<<matchOverlapThr<<", "<<overlapAreaFraction_2<<"<"<<matchOverlapThr<<")");
				#endif
				return false;
			}
		}
	}//close if matchSourcesByOverlap
	*/

	//## Check flux rel difference (if requested)	
	if(matchSourcesByFlux){
		double flux_1= pars1->fluxDensity;
		double flux_2= pars2->fluxDensity;
		double fluxRelDiff= (flux_1-flux_2)/flux_2;
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}

	return true;

}//close HaveSourceComponentMatch()



int ReadData()
{
	//Check filenames
	std::ifstream fileStream(fileName);	
	std::string line;
	if (fileStream.fail() || !fileStream.is_open()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<fileName<<" for reading...");
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
				ERROR_LOG("Failed to read line from file "<<fileName<<"!");
			#endif
			return -1; 
		}
    fileNames.push_back(filename);
	}//end file read
				
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<fileNames.size()<<" files present...");
	#endif

	//Finally reading source data
	source_pars.clear();
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
		long int nSourcePars= tileDataList[i]->GetNSourcePars();
		#ifdef LOGGING_ENABLED
			INFO_LOG("Tile no. "<<i+1<<": #"<<nSourcePars<<" sources to be matched...");
		#endif
	}

	return 0;

}//close ReadData()


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
	int wcsType= 0;
	WCS* wcs= 0;
	int coordSys= 0;//eJ2000
	std::vector<SourcePars*> pars;
	int sourceIndex= 0;
	sources.push_back(std::vector<Source*>());

	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);
		int type= aSource->Type;
		int simType= aSource->SimType;

		/*
		//Select source by type?
		if( selectTrueSourceByType && type!=selectedTrueSourceType && type!=-1) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Skip true source type "<<type<<" (selected type="<<selectedTrueSourceType<<")...");
			#endif
			continue;
		}
		*/

		//Select source by type?
		if(selectTrueSourceByType){
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
		/*
		if( selectTrueSourceBySimType && simType!=selectedTrueSourceSimType && simType!=-1) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Skip true source type "<<simType<<" (selected type="<<selectedTrueSourceSimType<<")...");
			#endif
			continue;
		}
		*/
		if(selectTrueSourceBySimType){
			bool skipSource= true;
			for(size_t j=0;j<ssimtypes.size();j++){
				if( ssimtypes[j]==-1 || type==ssimtypes[j]) {
					skipSource= false;
					break;
				}
			}
			if(skipSource) continue;
		}

		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading source (CAT"<<catalogIndex<<", name="<<aSource->GetName()<<")...");
		#endif

		//Copy source
		Source* source= new Source;
		*source= *aSource;
		
		//Compute wcs for this source collection if not done
		if(!wcs){
			wcs= source->GetWCS(coordSys);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to compute WCS from source no. "<<i+1<<" (name="<<source->GetName()<<")!");
				#endif
				return -1;
			}
		}

		//Fill source pars
		std::vector<SourcePars*> thisPars;
		if(FillSourcePars(thisPars,source,catalogIndex,sourceIndex,-1,wcs,wcsType)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill pars for source no. "<<i+1<<" (name="<<source->GetName()<<", collectioIndex="<<catalogIndex<<")!");
			#endif
			return -1;
		}

		//Add pars to pars
		pars.insert(pars.end(),thisPars.begin(),thisPars.end());

		//Add source in collection	
		sources[catalogIndex].push_back(source);

		sourceIndex++;
	}//end loop sources

	//Append source pars to main collection
	source_pars.insert(source_pars.end(),pars.begin(),pars.end());

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sources[catalogIndex].size()<<" sources to be cross-matched in this collection (#"<<source_pars.size()<<" source pars added from catalog no. "<<catalogIndex<<")...");
	#endif

	//Delete WCS for this collection
	WCSUtils::DeleteWCS(&wcs);

	return 0;

}//close ReadSourceData()

int FillSourcePars(std::vector<SourcePars*>& pars,Source* aSource,int catalogIndex,int sourceIndex,int nestedSourceIndex,WCS* wcs,int coordSystem)
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
	TVector2 contourCentroid= contour->Centroid;

	//## Fill source pars for island first
	SourcePars* sourcePars= new SourcePars(catalogIndex,sourceIndex,nestedSourceIndex);
	sourcePars->sname= sourceName;
	sourcePars->contour= contour;
	sourcePars->S= S;

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

			//Fill component par
			thisComponentPars= new ComponentPars;
			thisComponentPars->sname= sname;
			thisComponentPars->catalogIndex= catalogIndex;	
			thisComponentPars->fitComponentIndex= k;
			thisComponentPars->fitEllipse= fitEllipses[k];
			thisComponentPars->A= A;
			thisComponentPars->fluxDensity= fluxDensity;

			//Add component to source pars
			sourcePars->AddComponentPars(thisComponentPars);	
		}//end loop components

	}//close if has fit info


	//## Add source data to tile
	//Compute tile data index on the basis of contour centroid
	//If index is found, add to corresponding tile data
	long int tileIndex= MathUtils::FindGrid2DBin(
		contourCentroid.X(),contourCentroid.Y(),
		nTilesX,tileXmin,tileXmax,tileXstep,
		nTilesY,tileYmin,tileYmax,tileYstep
	);

	if(tileIndex>=0 && tileIndex<(long int)(tileDataList.size())){//Add to tile data
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Adding source (CAT"<<catalogIndex<<", name="<<sourceName<<", pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<") to list...");
		#endif
		tileDataList[tileIndex]->AddSourcePars(sourcePars);
	}
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Cannot find tile index or invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<") for source (CAT"<<catalogIndex<<", name="<<sourceName<<", pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<"), check tile index calculation!");
		#endif
		CodeUtils::DeletePtr<SourcePars>(sourcePars);
	}

	//Append this source pars to collection
	pars.push_back(sourcePars);



	/*
	//Fill source pars for fit componentd first
	bool useFWHM= true;
	bool convertToWCS= true;

	if(hasFitInfo){	
		//Get 0th contour converted to WCS
		int pixOffset= 0;
		bool computeContourPars= true;
		Contour* contour= aSource->GetWCSContour(0,wcs,coordSystem,pixOffset,computeContourPars);
		if(!contour){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute WCS contour for source "<<sourceName<<"!");
			#endif
			return -1;
		}

		//Get contour centroid
		TVector2 contourCentroid= contour->Centroid;

		//Get fitted pars & ellipse converted to WCS
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= fitPars.GetNComponents();
		std::vector<TEllipse*> fitEllipses;
		if(aSource->GetFitEllipses(fitEllipses,useFWHM,convertToWCS,wcs,coordSystem)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute WCS ellipse for fitted components of source "<<sourceName<<"!");
			#endif
			CodeUtils::DeletePtr<Contour>(contour);
			return -1;
		}

		//Check ellipses and pars size
		if(nComponents!=(int)(fitEllipses.size())){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Number of fit components shall be equal to fitted ellipses (this should not occur)!");
			#endif
			CodeUtils::DeletePtr<Contour>(contour);
			CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
			return -1;
		}

		//Get integrated & fitted flux density
		double S= aSource->GetS();
		double fluxDensity= 0;
		aSource->GetFluxDensity(fluxDensity);
		if(correctFlux){
			S/= beamArea;
			fluxDensity/= beamArea;
		}
		
		//Loop over fitted pars and ellipses and fill pars
		SourcePars* sourcePars= new SourcePars(catalogIndex,sourceIndex,nestedSourceIndex);
		sourcePars->sname= sourceName;
		sourcePars->contour= contour;
		sourcePars->fluxDensity= fluxDensity;
		sourcePars->S= S;
	
		//Fill component pars & ellipse
		ComponentPars* thisComponentPars= 0;
		for(int k=0;k<nComponents;k++){
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));
			double A= fitPars.GetParValue(k,"A");
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			if(correctFlux) fluxDensity/= beamArea;

			//Fill component par
			thisComponentPars= new ComponentPars;
			thisComponentPars->sname= sname;
			thisComponentPars->catalogIndex= catalogIndex;	
			thisComponentPars->fitComponentIndex= k;
			thisComponentPars->fitEllipse= fitEllipses[k];
			thisComponentPars->A= A;
			thisComponentPars->fluxDensity= fluxDensity;

			//Add component to source pars
			sourcePars->AddComponentPars(thisComponentPars);
			
		}//close if

		//Compute tile data index on the basis of contour centroid
		//If index is found, add to corresponding tile data
		long int tileIndex= MathUtils::FindGrid2DBin(
			contourCentroid.X(),contourCentroid.Y(),
			nTilesX,tileXmin,tileXmax,tileXstep,
			nTilesY,tileYmin,tileYmax,tileYstep
		);
		if(tileIndex<0){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Cannot find tile index of source (CAT"<<catalogIndex<<", name="<<sourceName<<", pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<")!");
			#endif
		}
		else if(tileIndex>=(long int)(tileDataList.size())){		
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<"), check tile index calculation!");
			#endif
		}
		else{//Add to tile data
			#ifdef LOGGING_ENABLED
				INFO_LOG("Adding source (CAT"<<catalogIndex<<", name="<<sourceName<<", pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<") to list...");
			#endif
			tileDataList[tileIndex]->AddSourcePars(sourcePars);
		}

		//Append this source pars to collection
		pars.push_back(sourcePars);

	}//close has fit info
	else{
	
		//## Fill pars for nested sources (if any)
		if(hasNestedSources){
			std::vector<Source*> nestedSources= aSource->GetNestedSources();
			for(size_t j=0;j<nestedSources.size();j++){
				if(FillSourcePars(pars,nestedSources[j],catalogIndex,sourceIndex,j,wcs,coordSystem)<0){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
					#endif
					return -1;
				}
			}//end loop sources
		}//close if nested sources

	}//close else no fit info
	*/


	//## Fill pars for nested sources (if any)
	if(hasNestedSources){
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			if(FillSourcePars(pars,nestedSources[j],catalogIndex,sourceIndex,j,wcs,coordSystem)<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
				#endif
				return -1;
			}
		}//end loop sources
	}//close if nested sources
	
	return 0;

}//close FillSourcePars()

void Save()
{
	//Save TTree to file
	if(outputFile && outputFile->IsOpen()){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving #"<<sourceCubes.size()<<" source cubes to file...");
		#endif
		outputFile->cd();		
		for(size_t i=0;i<sourceCubes.size();i++){
			sourceCube= sourceCubes[i];
			matchedSourceInfo->Fill();
		}
		matchedSourceInfo->Write();
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

