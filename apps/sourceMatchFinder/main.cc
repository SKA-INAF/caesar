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

#include <ConfigParser.h>
#include <Logger.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <EllipseUtils.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>

//WCS headers
#include <wcs.h>

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
	
	cout<<"-m, --matchSourcesByFlux \t Match sources by flux (default=false)"<<endl;
	cout<<"-M, --matchFluxRelDiffThr \t Flux relative difference threshold (default=0.5)"<<endl;
	cout<<"-a, --matchSourcesByOverlap \t Match sources by contour/ellipse overlap fraction (NB: skip if below thr) (default=false)"<<endl;
	cout<<"-A, --matchOverlapThr \t Overlap fraction (wrt to total area) threshold (default=0.5)"<<endl;
	cout<<"-t, --matchSourcesByPos \t Match sources by position (default=false)"<<endl;
	cout<<"-T, --matchPosThr=[POS_THESHOLD] \t Source centroid distance in arcsec below which two sources are matched (default=2.5)"<<endl;
	cout<<"-n, --minSourceMatchClusterSize=[MIN_SOURCE_MATCH_CLUSTER_SIZE] \t Minimum number of matched sources to be retained (e.g. skip source match cluster below this threshold) (default=2)"<<endl;

	cout<<"-f, --filterByType \t Consider only true sources with given type when searching the match (default=no)"<<endl;
	cout<<"-s, --selectedType=[TYPE] \t True source types to be crossmatched (1=COMPACT, 2=POINT-LIKE, 3=EXTENDED, 4=COMPACT_WITH_EXTENDED) (default=-1)"<<endl;
	cout<<"-F, --filterBySimType \t Consider only true sources with given sim type when searching the match (default=no)"<<endl;
	cout<<"-S, --selectedSimType=[TYPE] \t True source sim types to be crossmatched (eRingLike=1,eBubbleLike=2,eEllipseLike=3,eDiskLike=4,eBlobLike=5) (default=-1)"<<endl;
	cout<<"-e, --no-compactSourceCorrelation \t Disable correlation search for compact sources (default=enabled)"<<endl;
	cout<<"-E, --no-extendedSourceCorrelation \t Disable correlation search for extended sources (default=enabled)"<<endl;
	cout<<"-c, --correctFlux \t Correct rec integrated flux by beam area (default=no correction)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "input", required_argument, 0, 'i' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", optional_argument, 0, 'o' },
	{ "matchSourcesByFlux", no_argument, 0, 'm'},
	{ "matchFluxRelDiffThr", required_argument, 0, 'M'},
	{ "matchSourcesByOverlap", no_argument, 0, 'a'},	
	{ "matchOverlapThr", required_argument, 0, 'A'},	
	{ "matchSourcesByPos", no_argument, 0, 't'},
	{ "matchPosThr", required_argument, 0, 'T'},
	{ "minSourceMatchClusterSize",required_argument,0,'n'},

	{ "correctFlux", no_argument, 0, 'c'},
	{ "filterByType", no_argument, 0, 'f'},
	{ "filterBySimType", no_argument, 0, 'F'},
	{ "selectedType", required_argument, 0, 's'},	
	{ "selectedSimType", required_argument, 0, 'S'},
	{ "no-compactSourceCorrelation", no_argument, 0, 'e'},
	{ "no-extendedSourceCorrelation", no_argument, 0, 'E'},	
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
float matchPosThr= 2.5;//dist in arcsec below which two sources can be matched
int minSourceMatchClusterSize= 2;
bool correctFlux= false;
bool selectTrueSourceByType= false;//default=all true sources searched 
int selectedTrueSourceType= -1;
bool selectTrueSourceBySimType= false;//default=all true sources searched 
int selectedTrueSourceSimType= -1;
bool enableCompactSourceCorrelation= true;
bool enableExtendedSourceCorrelation= true;
bool applyFluxOverlapThreshold= false;
double fluxOverlapThr= 0.5;

//Globar vars
TFile* outputFile= 0;
std::string outputFileName= "MatchOutput.root";
TTree* matchedSourceInfo= 0;


struct SourceContourPars {
	int catalogIndex;	
	int sourceIndex;
	int nestedSourceIndex;
	std::string sname;
	double fluxDensity;
	Contour* contour;
};

struct SourcePars {
	int catalogIndex;
	int sourceIndex;
	int nestedSourceIndex;
	int fitComponentIndex;
	std::string sname;
	double Smax;
	double fluxDensity;
	TEllipse* fitEllipse;
};
std::vector<std::vector<Source*>> sources;
std::vector<SourcePars> source_pars;
std::vector<SourceContourPars> sourceContourPars;

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
			cout << *it << " ";
    }
    cout << endl;
	
		clusterIds.push_back(newClusterIds);
  }//close clique()
  
	std::vector<std::vector<size_t>>& clusterIds;
};


//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int FindSourceMatches();
int ReadData();
int ReadSourceData(std::string filename,int collectionIndex);
int FillSourcePars(std::vector<SourcePars>& pars,Source* aSource,int collectionIndex,int sourceIndex,int nestedSourceIndex=-1,WorldCoor* wcs=0,int coordSystem=0);
bool HaveSourceMatch(SourcePars& pars1,SourcePars& pars2);
void Save();
void Init();

int main(int argc, char *argv[])
{
	//================================
	//== Parse command line options
	//================================
	if(ParseOptions(argc,argv)<0){
		ERROR_LOG("Failed to parse command line options!");
		return -1;
	}
	
	//=======================
	//== Init
	//=======================
	INFO_LOG("Initializing data...");
	Init();

	//=======================
	//== Read data
	//=======================
	INFO_LOG("Reading source data from given files...");
	if(ReadData()<0){
		ERROR_LOG("Reading of source data failed!");
		return -1;
	}

	//=======================
	//== Find source matches
	//=======================
	INFO_LOG("Correlating source catalogs to find matches...");
	if(FindSourceMatches()<0){
		ERROR_LOG("Correlating source data to find matches failed!");
		return -1;
	}

	//=======================
	//== Save
	//=======================
	INFO_LOG("Saving data to file ...");
	Save();

	INFO_LOG("End source match finder");

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

	while((c = getopt_long(argc, argv, "hi:o:v:mM:aA:tT:n:cfFs:S:eE",options_tab, &option_index)) != -1) {
    
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
			case 'a':
			{
				matchSourcesByOverlap= true;
				break;
			}	
			case 'A':
			{
				matchOverlapThr= atof(optarg);
				break;
			}	
			case 't':
			{
				matchSourcesByPos= true;
				break;
			}
			case 'T':
			{
				matchPosThr= atof(optarg);
				break;
			}	
			case 'n':
			{
				minSourceMatchClusterSize= atoi(optarg);
				break;
			}
			case 'c':
			{
				correctFlux= true;
				break;
			}
			case 'f':
			{
				selectTrueSourceByType= true;
				break;
			}
			case 's':
			{
				selectedTrueSourceType= atoi(optarg);
				break;
			}
			case 'F':
			{
				selectTrueSourceBySimType= true;
				break;
			}
			case 'S':
			{
				selectedTrueSourceSimType= atoi(optarg);
				break;
			}	
			case 'e':
			{
				enableCompactSourceCorrelation= false;
				break;
			}	
			case 'E':
			{
				enableExtendedSourceCorrelation= false;
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
	LoggerManager::Instance().CreateConsoleLogger(sloglevel,"logger","System.out");
	
	return 0;

}//close ParseOptions()

void Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	//...
	//...

}//close Init()


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
			bool hasMatch= HaveSourceMatch(source_pars[i],source_pars[j]);
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
	for(size_t i=0;i<clusterIds.size();i++){
		ss<<"CLUST no. "<<i+1<<" {";
		for(size_t j=0;j<clusterIds[i].size()-1;j++){
			ss<<clusterIds[i][j]<<",";
		}
		ss<<clusterIds[i][clusterIds[i].size()-1]<<"}";
		INFO_LOG(ss.str());		
	}//end loop clusters


	return 0;

}//close FindSourceMatches()


bool HaveContourMatch(SourceContourPars& pars1,SourceContourPars& pars2)
{
	//## If sources are from the same collection return NO MATCH
	if(pars1.catalogIndex==pars2.catalogIndex) return false;

	//## Check fit ellipses
	Contour* cont1= pars1.contour;
	Contour* cont2= pars2.contour;
	if(!cont1 || !cont2){
		WARN_LOG("One/both fit ellipses are nullptr!");
		return false;
	}

}//close HaveContourMatch()


bool HaveSourceMatch(SourcePars& pars1,SourcePars& pars2)
{
	//## If sources are from the same collection return NO MATCH
	if(pars1.catalogIndex==pars2.catalogIndex) return false;

	//## Check fit ellipses
	TEllipse* ellipse1= pars1.fitEllipse;
	TEllipse* ellipse2= pars2.fitEllipse;
	if(!ellipse1 || !ellipse2){
		WARN_LOG("One/both fit ellipses are nullptr!");
		return false;
	}

	//## Check ellipse centroid sky distance (if requested)
	if(matchSourcesByPos){
		double Xc_1= ellipse1->GetX1();
		double Yc_1= ellipse1->GetY1();
		double Xc_2= ellipse2->GetX1();
		double Yc_2= ellipse2->GetY1();
		double posDist= AstroUtils::GetWCSPointDist_Haversine(Xc_1,Yc_1,Xc_2,Yc_2);
		if(fabs(posDist)>matchPosThr) return false;
	}

	//## Check ellipse overlap (if requested)
	if(matchSourcesByOverlap){
		double ellipseArea_1= MathUtils::ComputeEllipseArea(ellipse1);
		double ellipseArea_2= MathUtils::ComputeEllipseArea(ellipse2);
		double ellipseOverlapArea= -1;
		double err= 0;
		int rtn= 0;
		if(MathUtils::ComputeEllipseOverlapArea(ellipseOverlapArea,err,rtn,ellipse1,ellipse2)<0){
			WARN_LOG("Failed to compute ellipse overlap area (err status="<<rtn<<"), return no match!");
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
			return false;
		}
		else if(rtn==EllipseUtils::ELLIPSE1_INSIDE_ELLIPSE2 && overlapAreaFraction_2<matchOverlapThr){
			return false;
		}
		else if(rtn==EllipseUtils::ELLIPSE2_INSIDE_ELLIPSE1 && overlapAreaFraction_1<matchOverlapThr){
			return false;
		}
		else {
			if(overlapAreaFraction_1<matchOverlapThr || overlapAreaFraction_2<matchOverlapThr){
				return false;
			}
		}
	}//close if matchSourcesByOverlap

	//## Check flux rel difference (if requested)	
	if(matchSourcesByFlux){
		double flux_1= pars1.fluxDensity;
		double flux_2= pars2.fluxDensity;
		double fluxRelDiff= (flux_1-flux_2)/flux_2;
		if( fabs(fluxRelDiff)>matchFluxRelDiffThr ) return false;
	}

	return true;

}//close HaveSourceMatch()



int ReadData()
{

	//Check filenames
	std::ifstream fileStream(fileName);	
	std::string line;
	if (fileStream.fail() || !fileStream.is_open()){
		ERROR_LOG("Failed to open file "<<fileName<<" for reading...");
		return -1;
	}

	//Store filenames present in lists
	std::vector<std::string> fileNames;
	std::string filename= "";
	while (std::getline(fileStream, line)) {
  	std::istringstream iss(line);
    if (!(iss >> filename)) { 
			ERROR_LOG("Failed to read line from file "<<fileName<<"!");
			return -1; 
		}
    fileNames.push_back(filename);
	}//end file read
				
	INFO_LOG("#"<<fileNames.size()<<" files present...");

	//Finally reading source data
	source_pars.clear();
	for(size_t i=0;i<fileNames.size();i++){
		if(ReadSourceData(fileNames[i],i)<0){
			ERROR_LOG("Failed to read source data for file no. "<<i+1<<"!");
			return -1;
		}
	}//end loop files

	return 0;

}//close ReadData()


int ReadSourceData(std::string filename,int catalogIndex)
{
	//Open files
	TFile* f= new TFile(filename.c_str(),"READ");
	if(!f){
		ERROR_LOG("Failed to open file "<<filename<<"!");
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)f->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<filename<<"!");	
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	//Initialize WCS & source pars
	int wcsType= 0;
	WorldCoor* wcs= 0;
	int coordSys= 0;//eJ2000
	std::vector<SourcePars> pars;
	int sourceIndex= 0;

	//Read sources
	INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<filename<<"...");
	for(int i=0;i<sourceTree->GetEntries();i++){
		sourceTree->GetEntry(i);
		int type= aSource->Type;
		int simType= aSource->SimType;

		//Select source by type?
		if( selectTrueSourceByType && type!=selectedTrueSourceType && type!=-1) {
			DEBUG_LOG("Skip true source type "<<type<<" (selected type="<<selectedTrueSourceType<<")...");
			continue;
		}
		//Select source by sim type?
		if( selectTrueSourceBySimType && simType!=selectedTrueSourceSimType && simType!=-1) {
			DEBUG_LOG("Skip true source type "<<simType<<" (selected type="<<selectedTrueSourceSimType<<")...");
			continue;
		}

		//Copy source
		Source* source= new Source;
		*source= *aSource;
		
		//Compute wcs for this source collection if not done
		if(!wcs){
			wcs= source->GetWCS(coordSys);
			if(!wcs){
				ERROR_LOG("Failed to compute WCS from source no. "<<i+1<<" (name="<<source->GetName()<<")!");
				return -1;
			}
		}

		//Fill source pars
		std::vector<SourcePars> thisPars;
		if(FillSourcePars(thisPars,source,catalogIndex,sourceIndex,-1,wcs,wcsType)<0){
			ERROR_LOG("Failed to fill pars for source no. "<<i+1<<" (name="<<source->GetName()<<", collectioIndex="<<catalogIndex<<")!");
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

	INFO_LOG("#"<<sources[catalogIndex].size()<<" sources to be cross-matched in this collection (#"<<source_pars.size()<<" source pars added)...");

	return 0;

}//close ReadSourceData()

int FillSourcePars(std::vector<SourcePars>& pars,Source* aSource,int catalogIndex,int sourceIndex,int nestedSourceIndex,WorldCoor* wcs,int coordSystem)
{
	//Check input source
	if(!aSource) return -1;

	//Get source data
	std::string sourceName= std::string(aSource->GetName());
	bool hasFitInfo= aSource->HasFitInfo();
	bool hasNestedSources= aSource->HasNestedSources();
	double beamArea= aSource->GetBeamFluxIntegral();
	
	//Fill source pars for fit componentd first
	bool useFWHM= true;
	bool convertToWCS= true;

	if(hasFitInfo){	
		//Get fitted pars & ellipse converted to WCS
		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= fitPars.GetNComponents();
		std::vector<TEllipse*> fitEllipses;
		if(aSource->GetFitEllipses(fitEllipses,useFWHM,convertToWCS,wcs,coordSystem)<0){
			ERROR_LOG("Failed to compute WCS ellipse for fitted components of source "<<sourceName<<"!");
			return -1;
		}

		//Check ellipses and pars size
		if(nComponents!=(int)(fitEllipses.size())){
			ERROR_LOG("Number of fit components shall be equal to fitted ellipses (this should not occur)!");
			return -1;
		}

		//Loop over fitted pars and ellipses and fill pars
		for(int k=0;k<nComponents;k++){
			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));
			double Smax= fitPars.GetParValue(k,"A");
			double fluxDensity= fitPars.GetComponentFluxDensity(k);
			if(correctFlux) fluxDensity/= beamArea;

			//Fill pars
			SourcePars thisPars;
			thisPars.catalogIndex= catalogIndex;
			thisPars.sname= sname;
			thisPars.sourceIndex= sourceIndex;
			thisPars.nestedSourceIndex= nestedSourceIndex;
			thisPars.fitComponentIndex= k;
			thisPars.fitEllipse= fitEllipses[k];
			thisPars.Smax= Smax;
			thisPars.fluxDensity= fluxDensity;
			pars.push_back(thisPars);
		}//close if


	}//close has fit info
	
	//## Fill pars for nested sources (if any)
	if(hasNestedSources){
		std::vector<Source*> nestedSources= aSource->GetNestedSources();
		for(size_t j=0;j<nestedSources.size();j++){
			if(FillSourcePars(pars,nestedSources[j],catalogIndex,sourceIndex,j,wcs,coordSystem)<0){
				ERROR_LOG("Failed to get nested source pars for source "<<sourceName<<"!");
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
		outputFile->cd();		
		//...
		//...
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
	else if(verbosity>5) slevel= "DEBUG";
	else slevel= "OFF";

	return slevel;

}//close GetStringLogLevel()

