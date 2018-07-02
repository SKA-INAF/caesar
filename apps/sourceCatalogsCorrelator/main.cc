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

//ROOT headers
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;
using namespace Caesar;

void Usage(char* exeName){
	cout<<"=========== USAGE ==========="<<endl;
	cout<<"Usage: "<<exeName<<" [options]"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
  cout<<"-h, --help \t Show help message and exit"<<endl;
	cout<<"-i, --input=[INPUT_FILE] \t Input file name containing sources to be read in ROOT TTree (.root)"<<endl;
	cout<<"-I, --input-rec=[INPUT_REC_FILE] \t Input file name containing reconstructed/detected sources to be read in ROOT TTree (.root)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name "<<endl;
	cout<<"-t, --overlapThr=[THESHOLD] \t Fraction of matching pixels to consider sources equal (default=0.9)"<<endl;
	cout<<"-T, --posThr=[POS_THESHOLD] \t Pixel dustance below which two sources are matched (default=2.5)"<<endl;
	cout<<"-a, --applyFluxOverlapThr \t Apply flux overlap threshold in extended source matching (default=not applied)"<<endl;
	cout<<"-A, --fluxOverlapThr=[FLUX_THRESHOLD] \t Flux overlap threshold in extended source matching (default=0.5)"<<endl;
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
	{ "input-rec", required_argument, 0, 'I' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "output", optional_argument, 0, 'o' },
	{ "overlapThr", required_argument, 0, 't'},
	{ "posThr", required_argument, 0, 'T'},
	{ "correctFlux", no_argument, 0, 'c'},
	{ "filterByType", no_argument, 0, 'f'},
	{ "filterBySimType", no_argument, 0, 'F'},
	{ "selectedType", required_argument, 0, 's'},	
	{ "selectedSimType", required_argument, 0, 'S'},
	{ "applyFluxOverlapThr", no_argument, 0, 'a'},	
	{ "fluxOverlapThr", required_argument, 0, 'A'},
	{ "no-compactSourceCorrelation", no_argument, 0, 'e'},
	{ "no-extendedSourceCorrelation", no_argument, 0, 'E'},	
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
std::string fileName= "";
std::string fileName_rec= "";
int verbosity= 4;//INFO level
float matchOverlapThr= 0.9;//fraction of overlap above which two sources are matched
bool correctFlux= false;
float matchPosThr= 2.5;//#pixels below which two sources are matched
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
TTree* matchedExtSourceInfo= 0;
TTree* recSourceInfo= 0;
TTree* recExtSourceInfo= 0;
std::vector<Source*> sources;
std::vector<Source*> sources_rec;
int SourceFoundFlag;
std::string SourceName;
int SourceType;
int SourceSimType;
long int NPix;
double X0;
double Y0;
double X0_sweighted;
double Y0_sweighted;
double S;
double Smax;
double S_true;
double X0_true;
double Y0_true;
double fluxDensity_true;
double beamArea_true;
std::string SourceName_rec;
int SourceType_rec;
int IsNestedSource_rec;
long int NPix_rec;
double S_rec;
double X0_rec;
double Y0_rec;
double X0_sweighted_rec;
double Y0_sweighted_rec;
double Smax_rec;
double fluxDensity_rec;
double beamArea_rec;
double MatchFraction;
double MatchFraction_rec;
double Sratio;
double Sratio_rec;
double dX;
double dY;
double S_bkg;
double AvgBkg;
double AvgRMS;
int HasFitInfo;
int FitStatus;
double S_fit;
double X0_fit;
double Y0_fit;
double sigmaX_fit;
double sigmaY_fit;
double theta_fit;
double offset_fit;
double fluxDensity_fit;
double residualMin_fit;
double residualMax_fit;
double residualMean_fit;
double residualRMS_fit;
double residualMedian_fit;
double residualMAD_fit;
double chi2_fit;
double ndf_fit;
double ncomponents_fit;
const int MAX_NTRUE_MATCH_SOURCES= 1000;
int nTrueMatchedSources;
std::vector<std::string> SourceNameList_true;
std::vector<double> PosXList_true;
std::vector<double> PosYList_true;


struct MatchingSourceInfo {
	//Define struct fields		
	size_t sourceIndex;
	int fitComponentIndex;
	int nestedSourceIndex;

	MatchingSourceInfo(){
		sourceIndex= -1;
		fitComponentIndex= -1;
		nestedSourceIndex= -1;
	}
	MatchingSourceInfo(size_t index,int comp_index=-1,int nested_index=-1)
		: sourceIndex(index), fitComponentIndex(comp_index), nestedSourceIndex(nested_index)
	{}	

	bool operator < (const MatchingSourceInfo& obj) const {
    return std::tie(sourceIndex,fitComponentIndex,nestedSourceIndex) < std::tie(obj.sourceIndex,obj.fitComponentIndex,obj.nestedSourceIndex);
	}
};//close struct MatchingSourceInfo

std::map<MatchingSourceInfo,std::vector<int>> RecSourceAssociationMap;
std::map<MatchingSourceInfo,std::vector<int>> RecSourceAssociationMap_ext;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int CorrelateSourceCatalogs();
bool FindPointSourceMatch(int source_true_index);
bool FindExtendedSourceMatch(int source_true_index);
int FindSourceMatches();
int FindRealAndFakeSources();
int FindRecSourceMatch(Source* source_rec,int sourceIndex,int nestedIndex);
int FindPointSourceRealAndFakeSources(int source_rec_index);
int FindExtendedRealAndFakeSources(int source_rec_index);
int ReadSourceData();
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
	//== Read source data
	//=======================
	INFO_LOG("Reading source data from given files...");
	if(ReadSourceData()<0){
		ERROR_LOG("Reading of source data failed!");
		return -1;
	}

	//=======================
	//== Correlate catalogs
	//=======================
	INFO_LOG("Correlating source catalogs...");
	if(CorrelateSourceCatalogs()<0){
		ERROR_LOG("Correlating source data failed!");
		return -1;
	}

	//=======================
	//== Save
	//=======================
	INFO_LOG("Saving data to file ...");
	Save();

	INFO_LOG("End source catalog correlator");

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

	while((c = getopt_long(argc, argv, "hi:I:o:v:t:T:cfFs:S:eEaA:",options_tab, &option_index)) != -1) {
    
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
			case 't':
			{
				matchOverlapThr= atof(optarg);
				break;
			}
			case 'T':
			{
				matchPosThr= atof(optarg);
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
			case 'a':
			{
				applyFluxOverlapThreshold= true;
				break;
			}	
			case 'A':
			{
				fluxOverlapThr= atof(optarg);
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




int CorrelateSourceCatalogs()
{
	//Check source catalogs (true & rec)
	if(sources.empty() || sources_rec.empty()){
		WARN_LOG("One or both source collections are empty, nothing to correlate...");
		return 0;
	}

	//For each true source find best match in rec source catalogue
	if(FindSourceMatches()<0){
		ERROR_LOG("Failed to find source matches!");
		return -1;
	}

	//For each rec source store associated true sources
	if(FindRealAndFakeSources()<0){
		ERROR_LOG("Failed to store rec-true matches!");
		return -1;
	}

	return 0;

}//close CorrelateSourceCatalogs()


int FindSourceMatches()
{
	//Check source sizes
	if(sources.empty() || sources_rec.empty()){
		WARN_LOG("One or both source collections are empty, nothing to correlate...");
		return -1;
	}

	INFO_LOG("Correlating ("<<sources.size()<<","<<sources_rec.size()<<") sources...");
	long int nFoundSources= 0;
	long int nCompactSources= 0;
	long int nCompactSources_found= 0;	
	long int nExtendedSources= 0;
	long int nExtendedSources_found= 0;	

	//Loop over first catalogue (e.g. simulated/true sources)
	for(size_t i=0;i<sources.size();i++){
		
		if(i%100==0) INFO_LOG("Finding cross-match for source no. "<<i<<"/"<<sources.size()<<"...");

		//Find source match according to source type (e.g. compact/point-like vs extended)
		int sourceType= sources[i]->Type;
		if(enableCompactSourceCorrelation && (sourceType==Source::eCompact || sourceType==Source::ePointLike) ){
			nCompactSources++;
			bool isCompactSourceFound= FindPointSourceMatch(i);
			if(isCompactSourceFound) nCompactSources_found++;
		}
		if(enableExtendedSourceCorrelation && (sourceType==Source::eExtended || sourceType==Source::eCompactPlusExtended) ){
			nExtendedSources++;
			bool isExtendedSourceFound= FindExtendedSourceMatch(i);
			if(isExtendedSourceFound) nExtendedSources_found++;
		}

	}//end loop collection

	INFO_LOG("#"<<nCompactSources_found<<"/"<<nCompactSources<<" compact sources found (search enabled? "<<enableCompactSourceCorrelation<<") ...");
	INFO_LOG("#"<<nExtendedSources_found<<"/"<<nExtendedSources<<" extended sources found (search enabled? "<<enableExtendedSourceCorrelation<<") ...");

  return 0;

}//close FindSourceMatches()



bool FindPointSourceMatch(int source_true_index)
{
	//## Set true source info
	Source* source_true= sources[source_true_index]; 
	NPix= source_true->GetNPixels();
	SourceFoundFlag= 0;
	SourceName= std::string(source_true->GetName());
	SourceType= source_true->Type;
	SourceSimType= source_true->SimType;
	X0= source_true->X0;
	Y0= source_true->Y0;
	S= source_true->GetS();
	Smax= source_true->GetSmax();
	S_true= S;
	X0_true= X0;
	Y0_true= Y0;
	X0_sweighted= source_true->GetSx();
	Y0_sweighted= source_true->GetSy();
	beamArea_true= source_true->GetBeamFluxIntegral();
	bool hasTrueInfo= source_true->HasTrueInfo();
	if(hasTrueInfo){
		S_true= source_true->GetTrueFlux();
		source_true->GetTruePos(X0_true,Y0_true);
		fluxDensity_true= -1;//to be fixed
	}		

	//## Init stored rec source match info
	SourceFoundFlag= 0;
	SourceName_rec= "";
	NPix_rec= 0;
	S_rec= -1;
	Smax_rec= -1;
	X0_rec= -1;
	Y0_rec= -1;
	X0_sweighted_rec= -1;
	Y0_sweighted_rec= -1;
	SourceType_rec= -1;
	MatchFraction= 0.;
	MatchFraction_rec= 0.;
	S_bkg= 0;
	AvgBkg= 0;
	AvgRMS= 0;
	HasFitInfo= 0;
	FitStatus= SourceFitter::eFitUnknownStatus;
	S_fit= -1;
	X0_fit= -1;
	Y0_fit= -1;
	sigmaX_fit= -1;
	sigmaY_fit= -1;
	offset_fit= -1;
	theta_fit= -1;
	fluxDensity_fit= -1;
	residualMin_fit= -1;
	residualMax_fit= -1;
	residualMean_fit= -1;
	residualRMS_fit= -1;
	residualMedian_fit= -1;
	residualMAD_fit= -1;
	chi2_fit= -1;
	ndf_fit= -1;
	ncomponents_fit= -1;


	//## Search for extended source associations by matching pixels
	//SourcePosMatchPars match_info;
	std::vector<SourcePosMatchPars> match_info_list;
	//bool foundSource= source_true->FindSourceMatchByPos(match_info,sources_rec,matchPosThr);
	bool foundSource= source_true->FindSourceMatchByPos(match_info_list,sources_rec,matchPosThr);

	//## Store tree info
	if(foundSource){
		SourceFoundFlag= 1;	
		SourcePosMatchPars match_info= match_info_list[0];//in case of multiple match pick the first (they are sorted by posDiff)

		long int match_source_index= match_info.index;
		int componentIndex= match_info.fitComponentIndex;
		int nestedIndex= match_info.nestedIndex;
		float posDiff= match_info.posDiff;

		Source* match_source= sources_rec[match_source_index];
		if(nestedIndex!=-1) match_source= sources_rec[match_source_index]->GetNestedSource(nestedIndex);

		NPix_rec= match_source->GetNPixels();
		SourceName_rec= std::string(match_source->GetName());
		SourceType_rec= match_source->Type;
		X0_rec= match_source->X0;
		Y0_rec= match_source->Y0;
		X0_sweighted_rec= match_source->GetSx();
		Y0_sweighted_rec= match_source->GetSy();
		beamArea_rec= match_source->GetBeamFluxIntegral();
		S_rec= match_source->GetS();
		Smax_rec= match_source->GetSmax();
		long int NMatchingPixels= source_true->GetNMatchingPixels(match_source);
		MatchFraction= (double)(NMatchingPixels)/(double)(NPix);
		MatchFraction_rec= (double)(NMatchingPixels)/(double)(NPix_rec);
			
		//Correct flux from Jy/beam to Jy
		if(correctFlux){
			S_rec/= beamArea_rec;
			//Smax_rec/= beamArea_rec;
		}

		//Compute average bkg/rms
		std::vector<Pixel*> pixels= match_source->GetPixels();
		std::vector<double> bkgValues;
		std::vector<double> rmsValues;
		S_bkg= 0;
		for(size_t k=0;k<pixels.size();k++){	
			double bkgLevel= pixels[k]->GetBkg().first;
			double bkgRMS= pixels[k]->GetBkg().second;
			S_bkg+= bkgLevel;
			bkgValues.push_back(bkgLevel);
			rmsValues.push_back(bkgRMS);
		}
		AvgBkg= StatsUtils::GetMedianFast(bkgValues);
		AvgRMS= StatsUtils::GetMedianFast(rmsValues);

		//Store fit info
		if(componentIndex!=-1){
			HasFitInfo= 1;
			SourceFitPars fitPars= match_source->GetFitPars();
			FitStatus= fitPars.GetStatus();
			X0_fit= fitPars.GetParValue(componentIndex,"x0");	
			Y0_fit= fitPars.GetParValue(componentIndex,"y0");	
			S_fit= fitPars.GetParValue(componentIndex,"A");
			sigmaX_fit= fitPars.GetParValue(componentIndex,"sigmaX");
			sigmaY_fit= fitPars.GetParValue(componentIndex,"sigmaY");
			theta_fit= fitPars.GetParValue(componentIndex,"theta");
			offset_fit= fitPars.GetOffsetPar();
			fluxDensity_fit= fitPars.GetComponentFluxDensity(componentIndex);
			residualMin_fit= fitPars.GetResidualMin();
			residualMax_fit= fitPars.GetResidualMax();
			residualMean_fit= fitPars.GetResidualMean();
			residualRMS_fit= fitPars.GetResidualRMS();
			residualMedian_fit= fitPars.GetResidualMedian();
			residualMAD_fit= fitPars.GetResidualMAD();
			chi2_fit= fitPars.GetChi2();
			ndf_fit= fitPars.GetNDF();
			ncomponents_fit= fitPars.GetNComponents();

			//Correct flux from Jy/beam to Jy
			if(correctFlux){
				//S_fit/= beamArea_rec;
				fluxDensity_fit/= beamArea_rec;
			}
		}//close if has component index match


		INFO_LOG("True source "<<SourceName<<" (index="<<source_true_index<<", X0="<<X0<<", Y0="<<Y0<<", N="<<NPix<<") reconstructed by source "<<SourceName_rec<<" (index="<<match_source_index<<", X0="<<X0_rec<<", Y0="<<Y0_rec<<", N="<<NPix_rec<<"), NMatchingPixels="<<MatchFraction*NPix<<" f="<<MatchFraction<<" f_rec="<<MatchFraction_rec<<" (t="<<matchOverlapThr<<"), posDiff="<<posDiff<<" (posThr="<<matchPosThr<<")");


		//## Store rec-true association map
		//## NB: Find if this rec source was already associated to other true sources	
		//## Use all matched found previously	
		for(size_t k=0;k<match_info_list.size();k++){		
			
			long int match_source_index= match_info_list[k].index;
			int componentIndex= match_info_list[k].fitComponentIndex;
			int nestedIndex= match_info_list[k].nestedIndex;
			MatchingSourceInfo info(match_source_index,componentIndex,nestedIndex);

			std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap.find(info);
			if(RecSourceAssociationMap.empty() || it==RecSourceAssociationMap.end()){//item not found
				INFO_LOG("Match rec source (name="<<SourceName_rec<<", index="<<match_source_index<<", nestedIndex="<<nestedIndex<<", componentIndex="<<componentIndex<<") not found in map, adding it...");
				RecSourceAssociationMap[info].push_back(source_true_index);
			}
			else{
				//Find if true source was already associated to this source
				std::vector<int>::iterator vIt= std::find(RecSourceAssociationMap[info].begin(),RecSourceAssociationMap[info].end(),source_true_index);
				bool itemAlreadyPresent= (
					!RecSourceAssociationMap[info].empty() && 
					vIt!=RecSourceAssociationMap[info].end()
				);
				INFO_LOG("Match rec source (name="<<SourceName_rec<<", index="<<match_source_index<<", nestedIndex="<<nestedIndex<<", componentIndex="<<componentIndex<<") found in map, appending to it...");
				if(!itemAlreadyPresent){
					RecSourceAssociationMap[info].push_back(source_true_index);
				}
			}//close else
		}//end loop all matches

	
	}//close if match found
	

	//Fill ROOT tree with match info
	matchedSourceInfo->Fill();

	return foundSource;

}//close FindPointSourceMatch()

bool FindExtendedSourceMatch(int source_true_index)
{
	//## Set true source info
	Source* source_true= sources[source_true_index]; 
	NPix= source_true->GetNPixels();
	SourceName= std::string(source_true->GetName());
	SourceType= source_true->Type;
	SourceSimType= source_true->SimType;
	X0= source_true->X0;
	Y0= source_true->Y0;
	X0_sweighted= source_true->GetSx();
	Y0_sweighted= source_true->GetSy();
	S= source_true->GetS();
	Smax= source_true->GetSmax();
	S_true= S;
	X0_true= X0;
	Y0_true= Y0;
	beamArea_true= source_true->GetBeamFluxIntegral();
	bool hasTrueInfo= source_true->HasTrueInfo();
	if(hasTrueInfo){
		S_true= source_true->GetTrueFlux();
		source_true->GetTruePos(X0_true,Y0_true);
		fluxDensity_true= -1;//to be fixed
	}		

	//## Init stored rec source match info
	SourceFoundFlag= 0;
	SourceName_rec= "";
	NPix_rec= 0;
	S_rec= -1;
	Smax_rec= -1;
	X0_rec= -1;
	Y0_rec= -1;
	X0_sweighted_rec= -1;
	Y0_sweighted_rec= -1;
	SourceType_rec= -1;
	MatchFraction= 0.;
	MatchFraction_rec= 0.;
	Sratio= 0.;
	Sratio_rec= 0.;
	dX= 0;
	dY= 0;
	S_bkg= 0;
	AvgBkg= 0;
	AvgRMS= 0;

	//## Search for extended source associations by matching pixels
	SourceOverlapMatchPars match_info;
	bool foundSource= source_true->FindSourceMatchByOverlapArea(match_info,sources_rec,matchOverlapThr);

	//## Check flux overlap threshold
	if(applyFluxOverlapThreshold){
		double fluxOverlapThr_max= 1+fluxOverlapThr;
		if(fluxOverlapThr>0.5) fluxOverlapThr_max= 1 + (1-fluxOverlapThr);
		foundSource= (
			match_info.Sratio>fluxOverlapThr && match_info.Sratio<fluxOverlapThr_max &&
			match_info.Sratio_rec>fluxOverlapThr && match_info.Sratio_rec<fluxOverlapThr_max
		); 
	}

	//## Store tree info
	if(foundSource){
		SourceFoundFlag= 1;	

		long int match_source_index= match_info.index;
		MatchFraction= static_cast<double>(match_info.overlapFraction);
		MatchFraction_rec= static_cast<double>(match_info.overlapFraction_rec);
		Sratio= static_cast<double>(match_info.Sratio);
		Sratio_rec= static_cast<double>(match_info.Sratio_rec);
		dX= static_cast<double>(match_info.dX);
		dY= static_cast<double>(match_info.dY);
		NPix_rec= sources_rec[match_source_index]->GetNPixels();
		SourceName_rec= std::string(sources_rec[match_source_index]->GetName());
		SourceType_rec= sources_rec[match_source_index]->Type;
		X0_rec= sources_rec[match_source_index]->X0;
		Y0_rec= sources_rec[match_source_index]->Y0;
		X0_sweighted_rec= sources_rec[match_source_index]->GetSx();
		Y0_sweighted_rec= sources_rec[match_source_index]->GetSy();
		beamArea_rec= sources_rec[match_source_index]->GetBeamFluxIntegral();
		S_rec= sources_rec[match_source_index]->GetS();
		Smax_rec= sources_rec[match_source_index]->GetSmax();
		
		//Remove compact sources from total flux S_rec?
		//...

		//Correct flux from Jy/beam to Jy
		if(correctFlux){
			S_rec/= beamArea_rec;
			//Smax_rec/= beamArea_rec;
		}

		//Compute average bkg/rms
		std::vector<Pixel*> pixels= sources_rec[match_source_index]->GetPixels();
		std::vector<double> bkgValues;
		std::vector<double> rmsValues;
		S_bkg= 0;
		for(size_t k=0;k<pixels.size();k++){	
			double bkgLevel= pixels[k]->GetBkg().first;
			double bkgRMS= pixels[k]->GetBkg().second;
			S_bkg+= bkgLevel;
			bkgValues.push_back(bkgLevel);
			rmsValues.push_back(bkgRMS);
		}
		AvgBkg= StatsUtils::GetMedianFast(bkgValues);
		AvgRMS= StatsUtils::GetMedianFast(rmsValues);

		INFO_LOG("True source "<<SourceName<<" (index="<<source_true_index<<", X0="<<X0<<", Y0="<<Y0<<", N="<<NPix<<") reconstructed by source "<<SourceName_rec<<" (index="<<match_source_index<<", X0="<<X0_rec<<", Y0="<<Y0_rec<<", N="<<NPix_rec<<"), NMatchingPixels="<<MatchFraction*NPix<<" f="<<MatchFraction<<" f_rec="<<MatchFraction_rec<<" (t="<<matchOverlapThr<<")");

		
		//## Store rec-true association map
		//## NB: Find if this rec source was already associated to other true sources
		MatchingSourceInfo info(match_source_index,-1);
		std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap.find(info);
		if(RecSourceAssociationMap.empty() || it==RecSourceAssociationMap.end()){//item not found
			RecSourceAssociationMap[info].push_back(source_true_index);
		}
		else{
			//Find if true source was already associated to this source
			std::vector<int>::iterator vIt= std::find(RecSourceAssociationMap[info].begin(),RecSourceAssociationMap[info].end(),source_true_index);
			bool itemAlreadyPresent= (
				!RecSourceAssociationMap[info].empty() && 
				vIt!=RecSourceAssociationMap[info].end()
			);
			if(!itemAlreadyPresent){
				RecSourceAssociationMap[info].push_back(source_true_index);
			}
		}
		/*
		std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap_ext.find(info);
		if(RecSourceAssociationMap_ext.empty() || it==RecSourceAssociationMap_ext.end()){//item not found
			RecSourceAssociationMap_ext[info].push_back(source_true_index);
		}
		else{
			//Find if true source was already associated to this source
			std::vector<int>::iterator vIt= std::find(RecSourceAssociationMap_ext[info].begin(),RecSourceAssociationMap_ext[info].end(),source_true_index);
			bool itemAlreadyPresent= (
				!RecSourceAssociationMap_ext[info].empty() && 
				vIt!=RecSourceAssociationMap_ext[info].end()
			);
			if(!itemAlreadyPresent){
				RecSourceAssociationMap_ext[info].push_back(source_true_index);
			}
		}
		*/
	}//close if match found
	
	//## Fill ROOT tree with match info
	matchedExtSourceInfo->Fill();

	return foundSource;

}//close FindExtendedSourceMatch()

int FindRealAndFakeSources()
{
	//Loop over rec sources and find associations to true
	for(size_t i=0;i<sources_rec.size();i++){

		if(i%100==0) INFO_LOG("Finding associations for rec source no. "<<i<<"/"<<sources_rec.size()<<"...");

		int sourceType= sources_rec[i]->Type;

		//Store mother rec source associations to true sources
		INFO_LOG("Store mother rec source no. "<<i+1<<" (name="<<sources_rec[i]->GetName()<<") associations to true sources");
		if(FindRecSourceMatch(sources_rec[i],i,-1)<0){
			WARN_LOG("Failed to find true associations to rec source "<<i<<" (type="<<sourceType<<")!");
		}

		//Loop over nested sources
		std::vector<Source*> nestedSources= sources_rec[i]->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++){
			INFO_LOG("Store nested source no. "<<k+1<<" (name="<<nestedSources[k]->GetName()<<") of mother rec source no. "<<i+1<<" (name="<<sources_rec[i]->GetName()<<") associations to true sources");
		
			if(FindRecSourceMatch(nestedSources[k],i,k)<0){
				WARN_LOG("Failed to find true associations to nested source no. "<<k+1<<" of rec source "<<i<<" (type="<<sourceType<<")!");
			}
		}//end loop nested

		/*
		if(enableCompactSourceCorrelation && (sourceType==Source::eCompact || sourceType==Source::ePointLike) ){
			if(FindPointSourceRealAndFakeSources(i)<0){
				WARN_LOG("Failed to find true associations to rec source "<<i<<"!");
			}
		}
		if(enableExtendedSourceCorrelation && (sourceType==Source::eExtended || sourceType==Source::eCompactPlusExtended) ){
			if(FindExtendedRealAndFakeSources(i)<0){
				WARN_LOG("Failed to find true associations to rec source "<<i<<"!");
			}
		}
		*/

	}//end loop rec sources


	return 0;

}//close FindRealAndFakeSources()


int FindRecSourceMatch(Source* source_rec,int sourceIndex,int nestedIndex)
{
	//Get rec source info
	if(!source_rec){
		WARN_LOG("Null ptr to source given!");
		return -1;
	}
	std::string sname= std::string(source_rec->GetName());
	SourceName_rec= std::string(source_rec->GetName());
	IsNestedSource_rec= 0;
	if(nestedIndex!=-1) IsNestedSource_rec= 1;
	SourceType_rec= source_rec->Type;
	NPix_rec= source_rec->GetNPixels();			
	X0_rec= source_rec->X0;
	Y0_rec= source_rec->Y0;
	Smax_rec= source_rec->GetSmax();
	fluxDensity_rec= source_rec->GetS();
	beamArea_rec= source_rec->GetBeamFluxIntegral();
	if(correctFlux){
		//Smax_rec/= beamArea_rec;
		fluxDensity_rec/= beamArea_rec;
	}
	HasFitInfo= source_rec->HasFitInfo();
	sigmaX_fit= -1;
	sigmaY_fit= -1;
	theta_fit= -1;

	//Init association info
	nTrueMatchedSources= 0;
	SourceNameList_true.clear();
	PosXList_true.clear();
	PosYList_true.clear();

	//Find association to true sources
	if(HasFitInfo){
		SourceFitPars fitPars= source_rec->GetFitPars();
		for(int k=0;k<fitPars.GetNComponents();k++){
			SourceName_rec= sname + std::string(Form("_fitcomp%d",k+1));
			SourceType_rec= Source::ePointLike;
			X0_rec= fitPars.GetParValue(k,"x0");	
			Y0_rec= fitPars.GetParValue(k,"y0");
			Smax_rec= fitPars.GetParValue(k,"A");	
			sigmaX_fit= fitPars.GetParValue(k,"sigmaX");
			sigmaY_fit= fitPars.GetParValue(k,"sigmaY");
			theta_fit= fitPars.GetParValue(k,"theta");
			fluxDensity_rec= fitPars.GetComponentFluxDensity(k);
			if(correctFlux){
				//Smax_rec/= beamArea_rec;
				fluxDensity_rec/= beamArea_rec;
			}

			nTrueMatchedSources= 0;
			SourceNameList_true.clear();
			PosXList_true.clear();
			PosYList_true.clear();

			//Find true sources associated to rec source
			std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap.find(MatchingSourceInfo(sourceIndex,k,nestedIndex));
			if(it!=RecSourceAssociationMap.end()){
				std::vector<int> true_source_indexes= it->second;
				nTrueMatchedSources= static_cast<int>(true_source_indexes.size());
				for(size_t j=0;j<true_source_indexes.size();j++){
					int index= true_source_indexes[j];
					SourceName= std::string(sources[index]->GetName());
					X0_true= sources[index]->X0;
					Y0_true= sources[index]->Y0;
					bool hasTrueInfo= sources[index]->HasTrueInfo();
					if(hasTrueInfo) sources[index]->GetTruePos(X0_true,Y0_true);
					SourceNameList_true.push_back(SourceName);
					PosXList_true.push_back(X0_true);
					PosYList_true.push_back(Y0_true);
				}//end loop associated true sources
			}//close if found true source association

			//Fill tree
			recSourceInfo->Fill();
				
		}//end loop fitted components	
	}//close has fit info
	else{
			
		//Find true sources associated to rec source
		std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap.find(MatchingSourceInfo(sourceIndex,-1,nestedIndex));
		if(it!=RecSourceAssociationMap.end()){
			std::vector<int> true_source_indexes= it->second;
			nTrueMatchedSources= static_cast<int>(true_source_indexes.size());
			for(size_t j=0;j<true_source_indexes.size();j++){
				int index= true_source_indexes[j];
				SourceName= std::string(sources[index]->GetName());
				X0_true= sources[index]->X0;
				Y0_true= sources[index]->Y0;
				bool hasTrueInfo= sources[index]->HasTrueInfo();
				if(hasTrueInfo) sources[index]->GetTruePos(X0_true,Y0_true);
				SourceNameList_true.push_back(SourceName);
				PosXList_true.push_back(X0_true);
				PosYList_true.push_back(Y0_true);
			}//end loop associated true sources
		}//close if found true source association

		//Fill tree
		recSourceInfo->Fill();	

	}//close !hasFitInfo


	return 0;

}//close FindRealAndFakeSources()

int FindPointSourceRealAndFakeSources(int source_rec_index)
{
	//Get rec source info
	int sourceIndex= source_rec_index;
	Source* source_rec= sources_rec[source_rec_index]; 
	std::string sname= std::string(source_rec->GetName());
	SourceName_rec= std::string(source_rec->GetName());
	SourceType_rec= source_rec->Type;			
	X0_rec= source_rec->X0;
	Y0_rec= source_rec->Y0;
	Smax_rec= source_rec->GetSmax();
	fluxDensity_rec= source_rec->GetS();
	beamArea_rec= source_rec->GetBeamFluxIntegral();
	if(correctFlux){
		//Smax_rec/= beamArea_rec;
		fluxDensity_rec/= beamArea_rec;
	}
	HasFitInfo= source_rec->HasFitInfo();
	sigmaX_fit= -1;
	sigmaY_fit= -1;
	theta_fit= -1;

	//Init association info
	nTrueMatchedSources= 0;
	SourceNameList_true.clear();
	PosXList_true.clear();
	PosYList_true.clear();

	//Find association to true sources
	if(HasFitInfo){
		SourceFitPars fitPars= source_rec->GetFitPars();
		for(int k=0;k<fitPars.GetNComponents();k++){
			SourceName_rec= sname + std::string(Form("_%d",k+1));
			X0_rec= fitPars.GetParValue(k,"x0");	
			Y0_rec= fitPars.GetParValue(k,"y0");
			Smax_rec= fitPars.GetParValue(k,"A");
			fluxDensity_rec= fitPars.GetComponentFluxDensity(k);
			sigmaX_fit= fitPars.GetParValue(k,"sigmaX");
			sigmaY_fit= fitPars.GetParValue(k,"sigmaY");
			theta_fit= fitPars.GetParValue(k,"theta");
			if(correctFlux){
				//Smax_rec/= beamArea_rec;
				fluxDensity_rec/= beamArea_rec;
			}

			nTrueMatchedSources= 0;
			SourceNameList_true.clear();
			PosXList_true.clear();
			PosYList_true.clear();

			//Find true sources associated to rec source
			std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap.find(MatchingSourceInfo(sourceIndex,k));
			if(it!=RecSourceAssociationMap.end()){
				std::vector<int> true_source_indexes= it->second;
				nTrueMatchedSources= static_cast<int>(true_source_indexes.size());
				for(size_t j=0;j<true_source_indexes.size();j++){
					int index= true_source_indexes[j];
					SourceName= std::string(sources[index]->GetName());
					X0_true= sources[index]->X0;
					Y0_true= sources[index]->Y0;
					bool hasTrueInfo= sources[index]->HasTrueInfo();
					if(hasTrueInfo) sources[index]->GetTruePos(X0_true,Y0_true);
					SourceNameList_true.push_back(SourceName);
					PosXList_true.push_back(X0_true);
					PosYList_true.push_back(Y0_true);
				}//end loop associated true sources
			}//close if found true source association

			//Fill tree
			recSourceInfo->Fill();
				
		}//end loop fitted components	
	}//close has fit info
	else{
			
		//Find true sources associated to rec source
		std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap.find(MatchingSourceInfo(sourceIndex,-1));
		if(it!=RecSourceAssociationMap.end()){
			std::vector<int> true_source_indexes= it->second;
			nTrueMatchedSources= static_cast<int>(true_source_indexes.size());
			for(size_t j=0;j<true_source_indexes.size();j++){
				int index= true_source_indexes[j];
				SourceName= std::string(sources[index]->GetName());
				X0_true= sources[index]->X0;
				Y0_true= sources[index]->Y0;
				bool hasTrueInfo= sources[index]->HasTrueInfo();
				if(hasTrueInfo) sources[index]->GetTruePos(X0_true,Y0_true);
				SourceNameList_true.push_back(SourceName);
				PosXList_true.push_back(X0_true);
				PosYList_true.push_back(Y0_true);
			}//end loop associated true sources
		}//close if found true source association

		//Fill tree
		recSourceInfo->Fill();	

	}//close !hasFitInfo

	return 0;

}//close FindPointSourceRealAndFakeSources()

int FindExtendedRealAndFakeSources(int source_rec_index)
{
	//Get rec source info
	int sourceIndex= source_rec_index;
	Source* source_rec= sources_rec[source_rec_index]; 
	SourceName_rec= std::string(source_rec->GetName());
	SourceType_rec= source_rec->Type;			
	X0_rec= source_rec->X0;
	Y0_rec= source_rec->Y0;
	Smax_rec= source_rec->GetSmax();
	fluxDensity_rec= source_rec->GetS();
	beamArea_rec= source_rec->GetBeamFluxIntegral();
	if(correctFlux){
		//Smax_rec/= beamArea_rec;
		fluxDensity_rec/= beamArea_rec;
	}
			
	//Find true sources associated to rec source
	//NB: This should loop over all sources (independently of their type), as an extended source can reconstruct a compact source (for example when residual is not working properly)
	nTrueMatchedSources= 0;
	SourceNameList_true.clear();
	PosXList_true.clear();
	PosYList_true.clear();
		
	std::map<MatchingSourceInfo,std::vector<int>>::iterator it= RecSourceAssociationMap_ext.find(MatchingSourceInfo(sourceIndex,-1));
	if(it!=RecSourceAssociationMap_ext.end()){
		std::vector<int> true_source_indexes= it->second;
		nTrueMatchedSources= static_cast<int>(true_source_indexes.size());
		for(size_t j=0;j<true_source_indexes.size();j++){
			int index= true_source_indexes[j];
			SourceName= std::string(sources[index]->GetName());
			X0_true= sources[index]->X0;
			Y0_true= sources[index]->Y0;
			bool hasTrueInfo= sources[index]->HasTrueInfo();
			if(hasTrueInfo) sources[index]->GetTruePos(X0_true,Y0_true);
			SourceNameList_true.push_back(SourceName);
			PosXList_true.push_back(X0_true);
			PosYList_true.push_back(Y0_true);
		}//end loop associated true sources
	}//close if found true source association

	//Fill tree
	recExtSourceInfo->Fill();	

	return 0;

}//close FindExtendedRealAndFakeSources()


void Init(){

	//Open output file
	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	if(!matchedSourceInfo) matchedSourceInfo= new TTree("SourceMatchInfo","SourceMatchInfo");
	matchedSourceInfo->Branch("found",&SourceFoundFlag,"found/I");
	matchedSourceInfo->Branch("name",&SourceName);
	matchedSourceInfo->Branch("type",&SourceType,"type/I");
	matchedSourceInfo->Branch("simtype",&SourceSimType,"simtype/I");
	matchedSourceInfo->Branch("NPix",&NPix,"NPix/L");	
	matchedSourceInfo->Branch("S",&S,"S/D");
	matchedSourceInfo->Branch("Smax",&Smax,"Smax/D");	
	matchedSourceInfo->Branch("X0",&X0,"X0/D");
	matchedSourceInfo->Branch("Y0",&Y0,"Y0/D");
	matchedSourceInfo->Branch("S_true",&S_true,"S_true/D");
	matchedSourceInfo->Branch("X0_true",&X0_true,"X0_true/D");
	matchedSourceInfo->Branch("Y0_true",&Y0_true,"Y0_true/D");
	matchedSourceInfo->Branch("X0_sweighted",&X0_sweighted,"X0_sweighted/D");
	matchedSourceInfo->Branch("Y0_sweighted",&Y0_sweighted,"Y0_sweighted/D");
	matchedSourceInfo->Branch("fluxDensity_true",&fluxDensity_true,"fluxDensity_true/D");
	matchedSourceInfo->Branch("beamArea_true",&beamArea_true,"beamArea_true/D");
	matchedSourceInfo->Branch("name_rec",&SourceName_rec);	
	matchedSourceInfo->Branch("type_rec",&SourceType_rec,"type_rec/I");
	matchedSourceInfo->Branch("NPix_rec",&NPix_rec,"NPix_rec/L");	
	matchedSourceInfo->Branch("S_rec",&S_rec,"S_rec/D");
	matchedSourceInfo->Branch("Smax_rec",&Smax_rec,"Smax_rec/D");
	matchedSourceInfo->Branch("X0_rec",&X0_rec,"X0_rec/D");
	matchedSourceInfo->Branch("Y0_rec",&Y0_rec,"Y0_rec/D");
	matchedSourceInfo->Branch("X0_sweighted_rec",&X0_sweighted_rec,"X0_sweighted_rec/D");
	matchedSourceInfo->Branch("Y0_sweighted_rec",&Y0_sweighted_rec,"Y0_sweighted_rec/D");
	matchedSourceInfo->Branch("beamArea_rec",&beamArea_rec,"beamArea_rec/D");
	matchedSourceInfo->Branch("MatchFraction",&MatchFraction,"MatchFraction/D");
	matchedSourceInfo->Branch("MatchFraction_rec",&MatchFraction_rec,"MatchFraction_rec/D");
	matchedSourceInfo->Branch("S_bkg",&S_bkg,"S_bkg/D");	
	matchedSourceInfo->Branch("AvgBkg",&AvgBkg,"AvgBkg/D");	
	matchedSourceInfo->Branch("AvgRMS",&AvgRMS,"AvgRMS/D");	
	matchedSourceInfo->Branch("HasFitInfo",&HasFitInfo,"HasFitInfo/I");
	matchedSourceInfo->Branch("FitStatus",&FitStatus,"FitStatus/I");
	matchedSourceInfo->Branch("S_fit",&S_fit,"S_fit/D");
	matchedSourceInfo->Branch("X0_fit",&X0_fit,"X0_fit/D");
	matchedSourceInfo->Branch("Y0_fit",&Y0_fit,"Y0_fit/D");
	matchedSourceInfo->Branch("sigmaX_fit",&sigmaX_fit,"sigmaX_fit/D");
	matchedSourceInfo->Branch("sigmaY_fit",&sigmaY_fit,"sigmaY_fit/D");
	matchedSourceInfo->Branch("theta_fit",&theta_fit,"theta_fit/D");
	matchedSourceInfo->Branch("offset_fit",&offset_fit,"offset_fit/D");
	matchedSourceInfo->Branch("fluxDensity_fit",&fluxDensity_fit,"fluxDensity_fit/D");
	matchedSourceInfo->Branch("residualMin_fit",&residualMin_fit,"residualMin_fit/D");
	matchedSourceInfo->Branch("residualMax_fit",&residualMax_fit,"residualMax_fit/D");
	matchedSourceInfo->Branch("residualMean_fit",&residualMean_fit,"residualMean_fit/D");
	matchedSourceInfo->Branch("residualRMS_fit",&residualRMS_fit,"residualRMS_fit/D");
	matchedSourceInfo->Branch("residualMedian_fit",&residualMedian_fit,"residualMedian_fit/D");
	matchedSourceInfo->Branch("residualMAD_fit",&residualMAD_fit,"residualMAD_fit/D");
	matchedSourceInfo->Branch("chi2_fit",&chi2_fit,"chi2_fit/D");
	matchedSourceInfo->Branch("ndf_fit",&ndf_fit,"ndf_fit/D");
	matchedSourceInfo->Branch("ncomponents_fit",&ncomponents_fit,"ncomponents_fit/D");
	
	if(!matchedExtSourceInfo) matchedExtSourceInfo= new TTree("ExtSourceMatchInfo","ExtSourceMatchInfo");
	matchedExtSourceInfo->Branch("found",&SourceFoundFlag,"found/I");
	matchedExtSourceInfo->Branch("name",&SourceName);
	matchedExtSourceInfo->Branch("type",&SourceType,"type/I");
	matchedExtSourceInfo->Branch("simtype",&SourceSimType,"simtype/I");
	matchedExtSourceInfo->Branch("NPix",&NPix,"NPix/L");	
	matchedExtSourceInfo->Branch("S",&S,"S/D");
	matchedExtSourceInfo->Branch("Smax",&Smax,"Smax/D");	
	matchedExtSourceInfo->Branch("X0",&X0,"X0/D");
	matchedExtSourceInfo->Branch("Y0",&Y0,"Y0/D");
	matchedExtSourceInfo->Branch("X0_sweighted",&X0_sweighted,"X0_sweighted/D");
	matchedExtSourceInfo->Branch("Y0_sweighted",&Y0_sweighted,"Y0_sweighted/D");
	matchedExtSourceInfo->Branch("S_true",&S_true,"S_true/D");
	matchedExtSourceInfo->Branch("X0_true",&X0_true,"X0_true/D");
	matchedExtSourceInfo->Branch("Y0_true",&Y0_true,"Y0_true/D");
	matchedExtSourceInfo->Branch("fluxDensity_true",&fluxDensity_true,"fluxDensity_true/D");
	matchedExtSourceInfo->Branch("beamArea_true",&beamArea_true,"beamArea_true/D");
	matchedExtSourceInfo->Branch("name_rec",&SourceName_rec);	
	matchedExtSourceInfo->Branch("type_rec",&SourceType_rec,"type_rec/I");
	matchedExtSourceInfo->Branch("NPix_rec",&NPix_rec,"NPix_rec/L");	
	matchedExtSourceInfo->Branch("S_rec",&S_rec,"S_rec/D");
	matchedExtSourceInfo->Branch("Smax_rec",&Smax_rec,"Smax_rec/D");
	matchedExtSourceInfo->Branch("X0_rec",&X0_rec,"X0_rec/D");
	matchedExtSourceInfo->Branch("Y0_rec",&Y0_rec,"Y0_rec/D");
	matchedExtSourceInfo->Branch("X0_sweighted_rec",&X0_sweighted_rec,"X0_sweighted_rec/D");
	matchedExtSourceInfo->Branch("Y0_sweighted_rec",&Y0_sweighted_rec,"Y0_sweighted_rec/D");
	matchedExtSourceInfo->Branch("beamArea_rec",&beamArea_rec,"beamArea_rec/D");
	matchedExtSourceInfo->Branch("MatchFraction",&MatchFraction,"MatchFraction/D");
	matchedExtSourceInfo->Branch("MatchFraction_rec",&MatchFraction_rec,"MatchFraction_rec/D");
	matchedExtSourceInfo->Branch("Sratio",&Sratio,"Sratio/D");
	matchedExtSourceInfo->Branch("Sratio_rec",&Sratio_rec,"Sratio_rec/D");
	matchedExtSourceInfo->Branch("dX",&dX,"dX/D");
	matchedExtSourceInfo->Branch("dY",&dY,"dY/D");
	matchedExtSourceInfo->Branch("S_bkg",&S_bkg,"S_bkg/D");	
	matchedExtSourceInfo->Branch("AvgBkg",&AvgBkg,"AvgBkg/D");	
	matchedExtSourceInfo->Branch("AvgRMS",&AvgRMS,"AvgRMS/D");
	
	//Create rec source info (useful for reliability estimation)
	if(!recSourceInfo) recSourceInfo= new TTree("RecSourceInfo","RecSourceInfo");
	recSourceInfo->Branch("name_rec",&SourceName_rec);
	recSourceInfo->Branch("type_rec",&SourceType_rec,"type_rec/I");
	recSourceInfo->Branch("HasFitInfo",&HasFitInfo,"HasFitInfo/I");	
	recSourceInfo->Branch("IsNestedSource",&IsNestedSource_rec,"IsNestedSource/I");
	recSourceInfo->Branch("NPix_rec",&NPix_rec,"NPix_rec/L");	
	recSourceInfo->Branch("X0_rec",&X0_rec,"X0_rec/D");
	recSourceInfo->Branch("Y0_rec",&Y0_rec,"Y0_rec/D");
	recSourceInfo->Branch("Smax_rec",&Smax_rec,"Smax_rec/D");	
	recSourceInfo->Branch("sigmaX_fit",&sigmaX_fit,"sigmaX_fit/D");
	recSourceInfo->Branch("sigmaY_fit",&sigmaY_fit,"sigmaY_fit/D");
	recSourceInfo->Branch("theta_fit",&theta_fit,"theta_fit/D");	
	recSourceInfo->Branch("fluxDensity_rec",&fluxDensity_rec,"fluxDensity_rec/D");	
	recSourceInfo->Branch("beamArea_rec",&beamArea_rec,"beamArea_rec/D");	
	recSourceInfo->Branch("nTrueMatchedSources",&nTrueMatchedSources,"nTrueMatchedSources/I");
	recSourceInfo->Branch("name_true",&SourceNameList_true);
	recSourceInfo->Branch("X0_true",&PosXList_true);
	recSourceInfo->Branch("Y0_true",&PosYList_true);

	/*
	if(!recExtSourceInfo) recExtSourceInfo= new TTree("RecExtSourceInfo","RecExtSourceInfo");
	recExtSourceInfo->Branch("name_rec",&SourceName_rec);
	recExtSourceInfo->Branch("type_rec",&SourceType_rec,"type_rec/I");
	recExtSourceInfo->Branch("X0_rec",&X0_rec,"X0_rec/D");
	recExtSourceInfo->Branch("Y0_rec",&Y0_rec,"Y0_rec/D");
	recExtSourceInfo->Branch("Smax_rec",&Smax_rec,"Smax_rec/D");	
	recExtSourceInfo->Branch("fluxDensity_rec",&fluxDensity_rec,"fluxDensity_rec/D");	
	recExtSourceInfo->Branch("beamArea_rec",&beamArea_rec,"beamArea_rec/D");		
	recExtSourceInfo->Branch("nTrueMatchedSources",&nTrueMatchedSources,"nTrueMatchedSources/I");
	recExtSourceInfo->Branch("name_true",&SourceNameList_true);
	recExtSourceInfo->Branch("X0_true",&PosXList_true);
	recExtSourceInfo->Branch("Y0_true",&PosYList_true);
	*/

}//close Init()


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

int ReadSourceData()
{
	//Open files
	TFile* f= new TFile(fileName.c_str(),"READ");
	if(!f){
		ERROR_LOG("Failed to open file "<<fileName<<"!");
		return -1;
	}
	
	TFile* f_rec= new TFile(fileName_rec.c_str(),"READ");
	if(!f_rec){
		ERROR_LOG("Failed to open file "<<fileName_rec<<"!");
		return -1;
	}

	//Get access to source trees
	Source* aSource= 0;

	TTree* sourceTree= (TTree*)f->Get("SourceInfo");
	if(!sourceTree || sourceTree->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName<<"!");	
		return -1;
	}
	sourceTree->SetBranchAddress("Source",&aSource);

	TTree* sourceTree_rec= (TTree*)f_rec->Get("SourceInfo");
	if(!sourceTree_rec || sourceTree_rec->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName_rec<<"!");	
		return -1;
	}
	sourceTree_rec->SetBranchAddress("Source",&aSource);

	//Read sources
	INFO_LOG("Reading #"<<sourceTree->GetEntries()<<" sources in file "<<fileName<<"...");
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

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	INFO_LOG("#"<<sources.size()<<" true sources to be cross-matched...");

	INFO_LOG("Reading #"<<sourceTree_rec->GetEntries()<<" sources in file "<<fileName_rec<<"...");
	for(int i=0;i<sourceTree_rec->GetEntries();i++){
		sourceTree_rec->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources_rec.push_back(source);
	}//end loop sources

	INFO_LOG("#"<<sources_rec.size()<<" rec sources selected in the match catalogue...");

	return 0;

}//close ReadSourceData()

void Save()
{
	//Save TTRees to file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		
		if(matchedSourceInfo) matchedSourceInfo->Write();
		if(matchedExtSourceInfo) matchedExtSourceInfo->Write();
		if(recSourceInfo) recSourceInfo->Write();
		//if(recExtSourceInfo) recExtSourceInfo->Write();
		outputFile->Close();
	}
}//close Save()



