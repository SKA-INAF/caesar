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

#include <Source.h>
#include <DS9Region.h>
#include <DS9RegionParser.h>
#include <SourceExporter.h>

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
	cout<<"-i, --input=[INPUT_FILE] \t Input ROOT file produced by CAESAR containing the source collection to be matched with regions"<<endl;
	cout<<"-r, --region=[REGION_FILE] \t Input DS9 region file containing the region(s) to be used to find matches"<<endl;
	cout<<"-s, --stype=[SOURCE_TYPE] \t Source type to be cross-match (-1=all,1=compact,2=point-like,3=extended,4=ext+point) (default=-1)"<<endl;
	cout<<"-n, --nx=[NX] \t Number of divisions along X (default=4)"<<endl;
	cout<<"-m, --ny=[NY] \t Number of divisions along Y (default=4)"<<endl;
	cout<<"-P, --matchByPos \t Match source island and region by position centroid distance (default=false)"<<endl;
	cout<<"-O, --matchByOverlap \t Match source island and region by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-p, --matchComponentByPos \t Match source fit component and region by position centroid distance (default=false)"<<endl;
	cout<<"-b, --matchComponentByOverlap \t Match source fit component and region by contour overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-T, --matchPosThr=[POS_THESHOLD] \t Source island-region centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	cout<<"-A, --matchOverlapThr \t Source island-region contour overlap fraction (wrt to total area) threshold (default=0.5)"<<endl;
	cout<<"-c, --matchSourceComponent \t Match source fit component and region by position centroid (if enabled) and ellipse overlap (if enabled) (default=false)"<<endl;
	cout<<"-t, --compMatchPosThr=[POS_THESHOLD] \t Source fit component-region centroid distance in arcsec below which we have a match (default=2.5)"<<endl;
	cout<<"-a, --compMatchOverlapThr \t Source fit component-region contour overlap fraction (wrt to total area) threshold (default=0.5)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name (ROOT format) where to store matched sources (default=sources.root)"<<endl;
	cout<<"-R, --region-output=[REGION_OUTPUT_FILE] \t Output DS9 region file name where to store selected sources (default=sources.reg)"<<endl;
	cout<<"-C, --catalog-output=[CATALOG_OUTPUT_FILE] \t Output catalog file name where to store selected sources (default=catalog.dat)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "region", required_argument, 0, 'r' },
	{ "stype", required_argument, 0, 's' },
	{ "nx", required_argument, 0, 'n' },
	{ "ny", required_argument, 0, 'm' },
	{ "matchByPos", no_argument, 0, 'P'},
	{ "matchByOverlap", no_argument, 0, 'O'},	
	{ "matchPosThr", required_argument, 0, 'T'},
	{ "matchOverlapThr", required_argument, 0, 'A'},	
	{ "matchSourceComponent", no_argument, 0, 'c'},	
	{ "matchComponentByPos", no_argument, 0, 'p'},
	{ "matchComponentByOverlap", no_argument, 0, 'b'},	
	{ "compMatchPosThr", required_argument, 0, 't'},
	{ "compMatchOverlapThr", required_argument, 0, 'a'},	
	{ "output", required_argument, 0, 'o' },
	{ "region-output", required_argument, 0, 'R' },
	{ "catalog-output", required_argument, 0, 'C' },
	{ "verbosity", required_argument, 0, 'v'},
  {(char*)0, (int)0, (int*)0, (int)0}
};

//Options
std::string fileName= "";
std::string regionFileName= "";
std::string outputFileName= "sources.root";
std::string regionOutputFileName= "sources.reg";
std::string regionComponentsOutputFileName= "sources_fitcomp.reg";
std::string catalogOutputFileName= "catalog.dat";
std::string catalogComponentsOutputFileName= "catalog_fitcomp.dat";
std::vector<int> stypes;
bool matchSourcesByPos= false;
bool matchSourcesByOverlap= false;
bool matchSourceComponent= false;
bool matchSourceComponentsByPos= false;
bool matchSourceComponentsByOverlap= false;
float matchPosThr= 2.5;//dist in arcsec below which source island and region can be matched
float matchOverlapThr= 0.5;//fraction of overlap above which source island-region are matched
float compMatchPosThr= 2.5;//dist in arcsec below which source fit component and region can be matched
float compMatchOverlapThr= 0.5;//fraction of overlap above which source fit component and region are matched
int ds9WCSType= 0;//use original WCS type to save catalog
int verbosity= 4;//INFO level
TFile* inputFile= 0;
TTree* sourceTree= 0;
TTree* perfTree= 0;
TTree* configTree= 0;

//====================================
//     ComponentPars struct
//====================================
struct ComponentPars 
{
	//Constructor
	ComponentPars(){	
		fitComponentIndex= -1;
		sname= "";
		fitEllipse= 0;
		A= 0;
		fluxDensity= 0;
	}

	//Destructor
	~ComponentPars(){
		CodeUtils::DeletePtr<TEllipse>(fitEllipse);
	}

	//Vars
	int fitComponentIndex;
	std::string sname;
	double A;//fitted amplitude
	double fluxDensity;//fitted flux density
	TEllipse* fitEllipse;
};

//====================================
//     SourcePars struct
//====================================
struct SourcePars 
{
	//Constructor
	SourcePars(long int _sourceIndex,long int _nestedSourceIndex=-1)
		: sourceIndex(_sourceIndex), nestedSourceIndex(_nestedSourceIndex)
	{
		S= 0;
		fluxDensity= 0;
		contour= 0;
		componentPars.clear();
	}

	//Destructor
	~SourcePars(){
		CodeUtils::DeletePtr<Contour>(contour);
		CodeUtils::DeletePtrCollection<ComponentPars>(componentPars);
	}

	//Methods
	void AddComponentPars(ComponentPars* pars){
		componentPars.push_back(pars);
	}

	//Vars
	long int sourceIndex;
	long int nestedSourceIndex;
	std::string sname;
	Contour* contour;
	double S;//pixel-integrated flux
	double fluxDensity;//fitted flux density
	std::vector<ComponentPars*> componentPars;
};


//====================================
//     RegionPars struct
//====================================
struct RegionPars 
{
	//Constructor
	RegionPars(long int index)
		: regionIndex(index)
	{	
		contour= 0;
		ellipse= 0;
	}

	//Destructor
	~RegionPars(){
		CodeUtils::DeletePtr<TEllipse>(ellipse);
		CodeUtils::DeletePtr<Contour>(contour);
	}

	//Vars
	long int regionIndex;
	Contour* contour;
	TEllipse* ellipse;

};//close RegionPars()


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
		regionPars.clear();
	}

	//Destructor 
	~TileData(){
		CodeUtils::DeletePtrCollection<SourcePars>(sourcePars);
		CodeUtils::DeletePtrCollection<RegionPars>(regionPars);
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

	//Add region pars
	int AddRegionPars(RegionPars* pars){
		if(!pars) return -1;
		regionPars.push_back(pars);
		return 0;
	}

	//Get number of sources pars in this tile
	long int GetNSourcePars(){
		return static_cast<long int>(sourcePars.size());
	}

	//Get number of region pars in this tile
	long int GetNRegionPars(){
		return static_cast<long int>(regionPars.size());
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

	//Region pars for this tile
	std::vector<RegionPars*> regionPars;

};//close TileData()

//Globar vars
TFile* outputFile= 0;
TTree* outputTree= 0;

Source* m_source= 0;
std::vector<Source*> m_sources;
std::vector<DS9Region*> m_regions;
std::vector<TileData*> tileDataList;
std::vector<SourcePars*> source_pars;
std::vector<RegionPars*> region_pars;
WCS* wcs= 0;
int wcsType= -1;
bool useWCS= false;
float tileXmin= 0;
float tileXmax= 0;
float tileXstep= 0;
float tileYmin= 0;
float tileYmax= 0;
float tileYstep= 0;
long int nTilesX= 0;
long int nTilesY= 0;
bool correctFlux= false;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadRegionData(std::string filename);
int ReadSourceData(std::string filename);
int FillTileData();
int FillSourcePars(Source* aSource,int sourceIndex,int nestedSourceIndex=-1);
int FillRegionPars(DS9Region* aRegion,int regionIndex);
int FindSourceMatchesInTiles();
bool HaveSourceIslandMatch(SourcePars* sourcePars,RegionPars* regionPars);
bool HaveSourceComponentMatch(ComponentPars* componentPars,RegionPars* regionPars);
int CloneObjectsInFile(std::vector<std::string> excludedObjNames);
void Save();
int SaveSources();
int SaveDS9Regions();
int SaveCatalog();
int Init();
int InitGrid(float xmin,float xmax,float xstep,float ymin,float ymax,float ystep);
void ClearData();

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
		ClearData();
		return -1;
	}

	//=======================
	//== Read DS9 region data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading DS9 region file "<<regionFileName<<" ...");
	#endif
	if(ReadRegionData(regionFileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of DS9 region failed!");
		#endif
		ClearData();
		return -1;
	}

	//=======================
	//== Read source data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading source data from file "<<fileName<<" ...");
	#endif
	if(ReadSourceData(fileName)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of source data failed!");
		#endif
		ClearData();
		return -1;
	}

	//=================================
	//== Fill region/source tile data
	//=================================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Filling source/region tile data ...");
	#endif
	if(FillTileData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to fill source/region tile data!");
		#endif
		ClearData();
		return -1;
	}

	//=======================
	//== Find source matches
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Finding source/region matches ...");
	#endif
	if(FindSourceMatchesInTiles()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to find source/region matches!");
		#endif
		return -1;
	}

	//...
	//...

	//=======================
	//== Save
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving data to file ...");
	#endif
	Save();

	#ifdef LOGGING_ENABLED
		INFO_LOG("End source region selector");
	#endif


}//close main()


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

	while((c = getopt_long(argc, argv, "hi:r:s:n:m:cPOpbT:A:t:a:o:R:C:v:",options_tab, &option_index)) != -1) {
    
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
			case 'r':	
			{
				regionFileName= std::string(optarg);	
				break;	
			}
			case 's':	
			{
				int stype= atoi(optarg);
				stypes.push_back(stype);	
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				break;	
			}
			case 'n':	
			{
				nTilesX= atol(optarg);
				break;	
			}
			case 'm':	
			{
				nTilesY= atol(optarg);
				break;	
			}
			case 'c':
			{
				matchSourceComponent= true;
				break;
			}
			case 'P':
			{
				matchSourcesByPos= true;
				break;
			}
			
			case 'O':
			{
				matchSourcesByOverlap= true;
				break;
			}	
			case 'p':
			{
				matchSourceComponentsByPos= true;
				break;
			}
			case 'b':
			{
				matchSourceComponentsByOverlap= true;
				break;
			}	
			case 'T':
			{
				matchPosThr= atof(optarg);
				break;
			}	
			case 'A':
			{
				matchOverlapThr= atof(optarg);
				break;
			}	
			case 't':
			{
				compMatchPosThr= atof(optarg);
				break;
			}	
			case 'a':
			{
				compMatchOverlapThr= atof(optarg);
				break;
			}	
			case 'R':	
			{
				regionOutputFileName= std::string(optarg);	
				break;	
			}
			case 'C':	
			{
				catalogOutputFileName= std::string(optarg);	
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
 
	//for(size_t i=0;i<stypes.size();i++) cout<<"stype["<<i<<"]="<<stypes[i]<<endl;	

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
		INFO_LOG("matchSourceComponent? "<<matchSourceComponent);
		INFO_LOG("matchSourceComponentsByPos? "<<matchSourceComponentsByPos<<", compMatchPosThr(arcsec)="<<compMatchPosThr);
		INFO_LOG("matchSourceComponentsByOverlap? "<<matchSourceComponentsByOverlap<<", compMatchOverlapThr="<<compMatchOverlapThr);
	#endif

	//=======================
	//== CHECK ARGS 
	//=======================
	//Check input file name
	if(fileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty input file name given!");
		#endif
		return -1;
	}

	//Check input region file name
	if(regionFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty input region file name given!");
		#endif
		return -1;
	}

	//Check output file name
	if(outputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output file name given!");
		#endif
		return -1;
	}

	//Check region output file name
	if(regionOutputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty region output file name given!");
		#endif
		return -1;
	}

	//Check catalog output file name
	if(catalogOutputFileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty catalog output file name given!");
		#endif
		return -1;
	}

	//Set DS9 region component file name
	regionComponentsOutputFileName= regionOutputFileName;
	CodeUtils::ExtractFileNameFromPath(regionComponentsOutputFileName,true);
	regionComponentsOutputFileName+= "_fitcomp.reg";

	//Set catalog component file name
	catalogComponentsOutputFileName= catalogOutputFileName;
	CodeUtils::ExtractFileNameFromPath(catalogComponentsOutputFileName,true);
	catalogComponentsOutputFileName+= "_fitcomp.dat";

	return 0;

}//close ParseOptions()

int Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//Init output source TTree
	if(!outputTree) outputTree= new TTree("SourceInfo","SourceInfo");
	m_source= 0;
	m_sources.clear();
	outputTree->Branch("Source",&m_source);

	
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
					DEBUG_LOG("Tiles ("<<i<<","<<j<<") are neighbors");
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



void ClearData()
{
	//Delete WCS
	WCSUtils::DeleteWCS(&wcs);

	//Delete TTree
	if(outputTree){
		delete outputTree;
		outputTree= 0;
	}

	//Close file
	if(outputFile && outputFile->IsOpen()){
		outputFile->Close();
	}

}//close ClearData()


int ReadRegionData(std::string filename)
{
	//Read and parse DS9 regions
	if (DS9RegionParser::Parse(m_regions,filename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse DS9 region file "<<filename<<"!");
		#endif
		return -1;
	}

	//Check if empty
	if(m_regions.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No regions read from file "<<filename<<"!");
		#endif
		return -1;
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<m_regions.size()<<" regions read from file "<<filename<<"...");
		#endif
	}

	return 0;

}//close ReadRegionData()


int ReadSourceData(std::string filename)
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
	wcsType= m_regions[0]->csType;
	useWCS= (wcsType!=eIMG_CS && wcsType!=eUNKNOWN_CS);
	double sourceXmin= 1.e+99;
	double sourceXmax= -1.e+99;
	double sourceYmin= 1.e+99;
	double sourceYmax= -1.e+99;
	
	
	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Found #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	#endif

	
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);
		
		#ifdef LOGGING_ENABLED
			if(i%1000==0) INFO_LOG("Reading source no. "<<i+1<<"/"<<sourceTree->GetEntries()<<"...");
		#endif

		//Select source by type?
		int type= aSource->Type;
		bool skipSource= true;
		for(size_t j=0;j<stypes.size();j++){
			if( stypes[j]==-1 || type==stypes[j]) {
				skipSource= false;
				break;
			}
		}
		if(skipSource) continue;

		//Copy source
		Source* source= new Source;
		*source= *aSource;

		double X0= source->X0;
		double Y0= source->Y0;

		//Convert coordinates to WCS?
		if(useWCS){
			//Compute wcs for this source collection if not done
			if(!wcs){
				wcs= source->GetWCS(wcsType);
			}
			
			//Check WCS
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to compute WCS from source no. "<<i+1<<" (name="<<source->GetName()<<")!");
				#endif
				return -1;
			}

			//If using WCS for the matching, use WCS centroid to compute source min/max coords
			source->GetWCSPos(X0,Y0,wcs,wcsType);

		}//close if


		//Find min & max coordinates
		if(X0>sourceXmax) sourceXmax= X0;
		if(X0<sourceXmin) sourceXmin= X0;
		if(Y0>sourceYmax) sourceYmax= Y0;
		if(Y0<sourceYmin) sourceYmin= Y0;	

		
		//Add sources to list
		m_sources.push_back(source);
		
	}//end loop sources

	#ifdef LOGGING_ENABLED
		INFO_LOG("Source min/max coordinate range: (Xmin,Xmax)=("<<sourceXmin<<","<<sourceXmax<<"), (Ymin,Ymax)=("<<sourceYmin<<","<<sourceYmax<<")");
	#endif

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


	return 0;

}//close ReadSourceData()


int FillTileData()
{
	//Compute source pars and add to tile data
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing source pars and adding to tile data...");
	#endif

	for(size_t i=0;i<m_sources.size();i++)
	{
		if(FillSourcePars(m_sources[i],i,-1)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill pars for source no. "<<i+1<<" (name="<<m_sources[i]->GetName()<<")!");
			#endif
			return -1;
		}
	}//end loop sources

	//Compute region pars and add to tile data
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing region pars and adding to tile data...");
	#endif

	for(size_t i=0;i<m_regions.size();i++)
	{
		if(FillRegionPars(m_regions[i],i)<0){
			#ifdef LOGGING_ENABLED	
				ERROR_LOG("Failed to fill pars for region no. "<<i+1<<")!");
			#endif
			return -1;
		}
	}//end loop regions

	return 0;

}//close FillTileData()

int FillRegionPars(DS9Region* aRegion,int regionIndex)
{
	//####  METHOD ##############################
	//- Compute region ellipse
	//- Compute region contour
	//###########################################

	//Check input region
	if(!aRegion) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Nullptr to input region given!");
		#endif
		return -1;	
	}

	//Compute region ellipse
	TEllipse* ellipse= aRegion->GetEllipse();
	if(!ellipse){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute ellipse from region no. "<<regionIndex<<"!");
		#endif
		return -1;
	}

	//Compute region contour
	bool computePars= true;
	Contour* contour= aRegion->GetContour(computePars);
	if(!contour){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute contour from region no. "<<regionIndex<<"!");
		#endif
		delete ellipse;
		ellipse= 0;
		return -1;
	}

	//Create and fill region pars
	RegionPars* regionPars= new RegionPars(regionIndex);
	regionPars->contour= contour;
	regionPars->ellipse= ellipse;	

	//Compute tile data index on the basis of contour centroid
	//If index is found, add to corresponding tile data
	TVector2 contourCentroid= contour->Centroid;

	long int tileIndex= MathUtils::FindGrid2DBin(
		contourCentroid.X(),contourCentroid.Y(),
		nTilesX,tileXmin,tileXmax,tileXstep,
		nTilesY,tileYmin,tileYmax,tileYstep
	);
	if(tileIndex<0){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Cannot find tile index of region no. "<<regionIndex<<" (pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<")!");
		#endif
	}
	else if(tileIndex>=(long int)(tileDataList.size())){		
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<"), check tile index calculation!");
		#endif
	}
	else{//Add to tile data
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Adding region no. "<<regionIndex<<" (pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<") to list...");
		#endif
		tileDataList[tileIndex]->AddRegionPars(regionPars);
	}

	return 0;

}//close FillRegionPars()


int FillSourcePars(Source* aSource,int sourceIndex,int nestedSourceIndex)
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
		

	//Fill source pars for fit componentd first
	bool useFWHM= true;
	bool convertToWCS= useWCS;

	if(hasFitInfo){	
		//Get 0th contour converted to WCS
		int pixOffset= 0;
		bool computeContourPars= true;
		Contour* contour= aSource->GetWCSContour(0,wcs,wcsType,pixOffset,computeContourPars);
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
		if(aSource->GetFitEllipses(fitEllipses,useFWHM,convertToWCS,wcs,wcsType)<0){
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
		SourcePars* sourcePars= new SourcePars(sourceIndex,nestedSourceIndex);
		sourcePars->sname= sourceName;
		sourcePars->contour= contour;
		sourcePars->fluxDensity= fluxDensity;
		sourcePars->S= S;
	
		//Fill component pars & ellipse
		ComponentPars* thisComponentPars= 0;
		for(int k=0;k<nComponents;k++){
			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));
			double A= fitPars.GetParValue(k,"A");
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			if(correctFlux) fluxDensity/= beamArea;

			//Fill component par
			thisComponentPars= new ComponentPars;
			thisComponentPars->sname= sname;
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
				WARN_LOG("Cannot find tile index of source (name="<<sourceName<<", pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<")!");
			#endif
		}
		else if(tileIndex>=(long int)(tileDataList.size())){		
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid tile index found (index="<<tileIndex<<", tilesize="<<tileDataList.size()<<"), check tile index calculation!");
			#endif
		}
		else{//Add to tile data
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding source (name="<<sourceName<<", pos("<<contourCentroid.X()<<","<<contourCentroid.Y()<<") to list...");
			#endif
			tileDataList[tileIndex]->AddSourcePars(sourcePars);
		}

	}//close has fit info
	else{
	
		//## Fill pars for nested sources (if any)
		if(hasNestedSources){
			std::vector<Source*> nestedSources= aSource->GetNestedSources();
			for(size_t j=0;j<nestedSources.size();j++){
				if(FillSourcePars(nestedSources[j],sourceIndex,j)<0){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
					#endif
					return -1;
				}
			}//end loop sources
		}//close if nested sources

	}//close else no fit info

	return 0;

}//close FillSourcePars()


bool HaveSourceIslandMatch(SourcePars* sourcePars,RegionPars* regionPars)
{
	//## Check input data
	if(!sourcePars || !regionPars){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given for source/region pars!");
		#endif
		return false;
	}

	//## Check contour
	Contour* cont1= sourcePars->contour;
	Contour* cont2= regionPars->contour;
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
	double posThr= matchPosThr/3600.;//convert threshold in deg (contour are given in deg)
	if(matchSourcesByPos){	
		double posDist= AstroUtils::GetWCSPointDist_Haversine(centroid_1.X(),centroid_1.Y(),centroid_2.X(),centroid_2.Y());
		
		if(fabs(posDist)>posThr) {	
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("NO POS MATCH: Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<"), Region (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<"), dist(arcsec)="<<posDist*3600.);
			#endif
			return false;
		}
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("POS MATCH: Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<"), Source (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<"), dist(arcsec)="<<posDist*3600.);
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
				DEBUG_LOG("NO CONTOUR OVERLAP: Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<", area="<<contArea_1<<"), Region (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
			#endif
			return false;
		}
		else if(overlapFlag==eCONT1_INSIDE_CONT2 && overlapAreaFraction_2<matchOverlapThr){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("CONT1 INSIDE CONT2 (OVERLAP BELOW THR): Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<", area="<<contArea_1<<"), Region (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
			#endif
			return false;
		}
		else if(overlapFlag==eCONT2_INSIDE_CONT1 && overlapAreaFraction_1<matchOverlapThr){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("CONT2 INSIDE CONT1 (OVERLAP BELOW THR): Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<", area="<<contArea_1<<"), Region (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
			#endif
			return false;
		}
		else if(overlapFlag==eCONT_OVERLAPPING){
			if(overlapAreaFraction_1<matchOverlapThr || overlapAreaFraction_2<matchOverlapThr){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("CONT OVERLAP BELOW THR: Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<", area="<<contArea_1<<"), Region (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
				#endif
				return false;
			}
		}
	
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("CONTOUR OVERLAP MATCH: Source (index="<<sourcePars->sourceIndex<<", nestedIndex="<<sourcePars->nestedSourceIndex<<", name="<<sourcePars->sname<<", pos="<<centroid_1.X()<<","<<centroid_1.Y()<<", area="<<contArea_1<<"), Region (index="<<regionPars->regionIndex<<", pos("<<centroid_2.X()<<","<<centroid_2.Y()<<", area="<<contArea_2<<"), overlapArea="<<overlapArea<<", overlapAreaFraction_1="<<overlapAreaFraction_1<<", overlapAreaFraction_2="<<overlapAreaFraction_2);
		#endif

	}//close if

	return true;

}//close HaveSourceIslandMatch()


bool HaveSourceComponentMatch(ComponentPars* componentPars,RegionPars* regionPars)
{
	//## Check data
	if(!componentPars || !regionPars){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given for component/region pars!");
		#endif
		return false;
	}

	//## Check fit ellipses
	TEllipse* ellipse1= componentPars->fitEllipse;
	TEllipse* ellipse2= regionPars->ellipse;
	if(!ellipse1 || !ellipse2){		
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/both fit ellipses are nullptr!");
		#endif
		return false;
	}

	//## Check ellipse centroid sky distance (if requested)
	double posThr= compMatchPosThr/3600.;//convert in deg as ellipse coords are in deg
	if(matchSourceComponentsByPos){
		double Xc_1= ellipse1->GetX1();
		double Yc_1= ellipse1->GetY1();
		double Xc_2= ellipse2->GetX1();
		double Yc_2= ellipse2->GetY1();
		double posDist= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);
		if(fabs(posDist)>posThr) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("NO ELLIPSE POS MATCH: Ellipse pos diff above thr ("<<posDist<<">"<<posThr<<", pos1("<<Xc_1<<","<<Yc_1<<"), pos2("<<Xc_2<<","<<Yc_2<<"))");
			#endif
			return false;
		}
		#ifdef LOGGING_ENABLED
			INFO_LOG("ELLIPSE POS MATCH: Ellipse pos diff below thr ("<<posDist<<"<"<<posThr<<", pos1("<<Xc_1<<","<<Yc_1<<"), pos2("<<Xc_2<<","<<Yc_2<<"))");
		#endif
	}

	//## Check ellipse overlap (if requested)
	if(matchSourceComponentsByOverlap){
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
				DEBUG_LOG("Disjoint ellipses");
			#endif
			return false;
		}
		else if(rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 && overlapAreaFraction_2<compMatchOverlapThr){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Ellipse 1 inside 2 below overlap thr ("<<overlapAreaFraction_2<<"<"<<compMatchOverlapThr<<")");
			#endif
			return false;
		}
		else if(rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1 && overlapAreaFraction_1<compMatchOverlapThr){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Ellipse 2 inside 1 below overlap thr ("<<overlapAreaFraction_1<<"<"<<compMatchOverlapThr<<")");
			#endif
			return false;
		}
		else {
			if(overlapAreaFraction_1<compMatchOverlapThr || overlapAreaFraction_2<compMatchOverlapThr){	
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Ellipse overlap below thr ("<<overlapAreaFraction_1<<"<"<<compMatchOverlapThr<<", "<<overlapAreaFraction_2<<"<"<<compMatchOverlapThr<<")");
				#endif
				return false;
			}
		}
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("ELLIPSE OVERLAP MATCH: overlapArea1="<<overlapAreaFraction_1<<">"<<compMatchOverlapThr<<", overlapArea2="<<overlapAreaFraction_2<<">"<<compMatchOverlapThr<<")");
		#endif

	}//close if matchSourcesByOverlap

	return true;

}//close HaveSourceComponentMatch()


int FindSourceMatchesInTiles()
{
	//####################################
	//##     METHOD
	//####################################
	//1) Compare source by contour (centroid, overlap area)
	//2) Compare source component ellipses (centroid, overlap)
	//####################################

	//Loop over tiles and find matches with regions in the same tile and in neighbor tiles
	long int nTotSources= 0;
	long int nMatchedSources= 0;
	long int nTotSourceComponents= 0;
	long int nMatchedSourceComponents= 0;

	for(size_t i=0;i<tileDataList.size();i++)
	{
		long int NSourcePars= tileDataList[i]->GetNSourcePars();
		long int NRegionPars= tileDataList[i]->GetNRegionPars();
		std::vector<size_t> neighbors= tileDataList[i]->neighborTileIndexes;

		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<NSourcePars<<" sources to be cross-matched against #"<<NRegionPars<<" regions in tile no. "<<i+1<<" ...");
		#endif

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
			

			//Get source fit info
			bool hasFitInfo= source->HasFitInfo();
			int nFitComponents= 0;
			if(hasFitInfo) {
				SourceFitPars fitPars= source->GetFitPars();
				nFitComponents= fitPars.GetNComponents();
			}
			
			//Count total number of sources
			nTotSources++;
			nTotSourceComponents+= nFitComponents;
	
			int nMatches= 0;

			//Loop over regions in the same tile
			std::vector<std::pair<long int,long int>> matchedRegionIndexes;

			for(long int k=0;k<NRegionPars;k++){
				RegionPars* regionPars= (tileDataList[i]->regionPars)[k];
				bool hasMatch= HaveSourceIslandMatch(sourcePars,regionPars);
				if(!hasMatch) continue;

				//Record match
				nMatches++;
				matchedRegionIndexes.push_back(std::make_pair(i,k));
			}//end loop regions in this tile

			

			//Loop over regions in neighbor tiles
			for(size_t l=0;l<neighbors.size();l++){
				size_t neighborIndex= neighbors[l];
				if(i==neighborIndex) continue;//skip same tile

				long int NRegionPars_neighbor= tileDataList[neighborIndex]->GetNRegionPars();

				for(long int k=0;k<NRegionPars_neighbor;k++){
					RegionPars* regionPars= (tileDataList[neighborIndex]->regionPars)[k];
					bool hasMatch= HaveSourceIslandMatch(sourcePars,regionPars);
					if(!hasMatch) continue;

					//Record match
					nMatches++;
					matchedRegionIndexes.push_back(std::make_pair(neighborIndex,k));
				}//end loop regions in this neighbor tile

			}//end loop neighbor tiles

			//## Tag source according to match results
			//## Set to unknown flag if not assigned to any region
			if(nMatches<=0){
				source->SetFlag(eFake);
				if(hasFitInfo){
					SourceFitPars fitPars= source->GetFitPars();
					for(int l=0;l<nFitComponents;l++) fitPars.SetComponentFlag(l,eFake);
					source->SetFitPars(fitPars);
				}

				#ifdef LOGGING_ENABLED
					INFO_LOG("Tagged source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<", name="<<source->GetName()<<") and its component as fake (no match found)");
				#endif

				//Skip further processing
				continue;
			}
			else if(nMatches==1){
				nMatchedSources++;
				source->SetFlag(eReal);
				source->SetType(ePointLike);
				#ifdef LOGGING_ENABLED
					INFO_LOG("Tagged source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<", name="<<source->GetName()<<") as real and point-like (1 match found)");
				#endif
			}
			else if(nMatches>1){
				nMatchedSources++;
				source->SetFlag(eReal);
				source->SetType(eCompact);
				#ifdef LOGGING_ENABLED
					INFO_LOG("Tagged source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<", name="<<source->GetName()<<") as real and compact ("<<nMatches<<" match found)");
				#endif
			}

			//Match source components
			if(matchSourceComponent){
				

				for(size_t l=0;l<componentPars.size();l++){
					ComponentPars* cpars= componentPars[l];
					int nComponentMatches= 0;

					for(size_t k=0;k<matchedRegionIndexes.size();k++){
						long int tileIndex= matchedRegionIndexes[k].first;
						long int regionIndex= matchedRegionIndexes[k].second;
						RegionPars* regionPars= (tileDataList[tileIndex]->regionPars)[regionIndex];
				
						bool hasMatch= HaveSourceComponentMatch(cpars,regionPars);
						if(!hasMatch) continue;
								
						nComponentMatches++;
					}//end loop regions

					//Tag source component
					SourceFitPars fitPars= source->GetFitPars();
					if(nComponentMatches<=0){
						fitPars.SetComponentFlag(l,eFake);
						fitPars.SetComponentType(l,eUnknownType);

						#ifdef LOGGING_ENABLED
							INFO_LOG("Tagged component no. "<<l+1<<" of source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<", name="<<source->GetName()<<") as fake (no match found)");
						#endif
					}
					else if(nComponentMatches==1){
						nMatchedSourceComponents++;
						fitPars.SetComponentFlag(l,eReal);
						fitPars.SetComponentType(l,ePointLike);
						#ifdef LOGGING_ENABLED
							INFO_LOG("Tagged component no. "<<l+1<<" of source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<", name="<<source->GetName()<<") as real and point-like (1 match found)");
						#endif
					}
					else if(nComponentMatches>1){
						nMatchedSourceComponents++;
						fitPars.SetComponentFlag(l,eReal);
						fitPars.SetComponentType(l,eCompact);
						#ifdef LOGGING_ENABLED
							INFO_LOG("Tagged component no. "<<l+1<<" of source (index="<<sourceIndex<<", nestedIndex="<<nestedSourceIndex<<", name="<<source->GetName()<<") as real and compact ("<<nComponentMatches<<" match found)");
						#endif
					}
					source->SetFitPars(fitPars);

				}//end loop component pars
			}//close if matchSourceComponent

		}//end loop sources in this tile
	}//end loop tiles

	#ifdef LOGGING_ENABLED
		INFO_LOG("#matched/tot sources: "<<nMatchedSources<<"/"<<nTotSources);
		INFO_LOG("#matched/tot source fit components: "<<nMatchedSourceComponents<<"/"<<nTotSourceComponents);
	#endif

	return 0;

}//close FindSourceMatchesInTiles()


int SaveDS9Regions()
{
	//Save DS9 regions for islands
	bool convertDS9RegionsToWCS= false;
	int status= SourceExporter::WriteToDS9(regionOutputFileName,m_sources,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	//Save DS9 regions for components
	status= SourceExporter::WriteComponentsToDS9(regionComponentsOutputFileName,m_sources,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	return 0;

}//close SaveDS9Regions()


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


int SaveCatalog()
{
	//Return if no sources are found
	if(m_sources.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No sources selected, no catalog file will be written!");
		#endif
		return 0;
	}

	//Retrieve source WCS
	WCS* wcs= m_sources[0]->GetWCS(ds9WCSType);
	if(!wcs) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute WCS from sources!");
		#endif
	}	

	//Saving island/blob catalog to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogOutputFileName<<" ...");
	#endif
	bool dumpNestedSourceInfo= true;
	int status= SourceExporter::WriteToAscii(catalogOutputFileName,m_sources,dumpNestedSourceInfo,ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source catalog to file "<<catalogOutputFileName<<" failed!");
		#endif
	}
	
	//Saving source fitted components to ascii file
	#ifdef LOGGING_ENABLED
		INFO_LOG("Writing source catalog to file "<<catalogComponentsOutputFileName<<" ...");
	#endif
	status= SourceExporter::WriteComponentsToAscii(catalogComponentsOutputFileName,m_sources,dumpNestedSourceInfo,ds9WCSType,wcs);
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Writing source fitted component catalog to file "<<catalogComponentsOutputFileName<<" failed!");
		#endif
	}

	return 0;

}//close SaveCatalog()


void Save()
{
	//Save TTree to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		
		//Clone all objects present in input file but the Source TTree
		CloneObjectsInFile({"SourceInfo"});

		//Write selected source TTree
		SaveSources();

		//Close file
		outputFile->Close();

	}//close if

	//Save DS9 regions with selected sources	
	#ifdef LOGGING_ENABLED
		ERROR_LOG("Saving DS9 regions with selected sources...");
	#endif
	SaveDS9Regions();

	//Save ascii catalog with selected sources
	#ifdef LOGGING_ENABLED
		ERROR_LOG("Saving ascii catalog with selected sources...");
	#endif
	SaveCatalog();

}//close Save()

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

std::string GetStringLogLevel(int verbosity)
{
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
