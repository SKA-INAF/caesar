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
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>
#include <Graph.h>
#include <Consts.h>

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
	cout<<"-I, --input-rec=[INPUT_FILE_REC] \t Input FITS file name with rec map (.fits)"<<endl;
	cout<<"-o, --output=[OUTPUT_FILE] \t Output file name containing convolved sources (.root)"<<endl;
	cout<<"-O, --output-map=[OUTPUT_FILE_MAP] \t Output FITS file name containing convolved map (.fits)"<<endl;
	cout<<"-r, --output-ds9=[OUTPUT_FILE] \t Output DS9 region file name "<<endl;
	cout<<"-B, --bmaj=[BMAJ] \t Bmaj in arcsec (NB: Overridden if --input-rec is given)"<<endl;
	cout<<"-b, --bmin=[BMAJ] \t Bmin in arcsec (NB: Overridden if --input-rec is given)"<<endl;
	cout<<"-a, --bpa=[BMAJ] \t Bpa in degrees (NB: Overridden if --input-rec is given)"<<endl;
	cout<<"-f, --fluxtruncthr=[FLUX_TRUNC_THRESHOLD] \t Flux loss threshold for source truncation (default=0.001)"<<endl;
	cout<<"-u, --userthreshold \t Use fixed threshold (provided in option -t/-T) (default=no)"<<endl;
	cout<<"-t, --threshold=[THRESHOLD] \t Flux threshold below which pixels are removed from sources (default=0)"<<endl;
	cout<<"-T, --threshold-ext=[THRESHOLD_EXT] \t Flux threshold below which pixels are removed from extended sources (default=0)"<<endl;
	cout<<"-m, --mergesources \t Merge overlapping thresholds after convolution (default=no)"<<endl;
	cout<<"-e, --mergecompactsources \t Enable merging of compact sources (default=no)"<<endl;
	cout<<"-s, --nsigmas=[NSIGMAS] \t Number of sigmas used in convolution gaussian kernel (default=10)"<<endl;
	cout<<"-v, --verbosity=[LEVEL] \t Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)"<<endl;
	cout<<"=============================="<<endl;
}//close Usage()

static const struct option options_tab[] = {
  /* name, has_arg, &flag, val */
  { "help", no_argument, 0, 'h' },
	{ "verbosity", required_argument, 0, 'v'},
	{ "input", required_argument, 0, 'i' },
	{ "input-rec", required_argument, 0, 'I' },	
	{ "output", required_argument, 0, 'o' },
	{ "output-map", required_argument, 0, 'O' },
	{ "output-ds9", required_argument, 0, 'r' },
	{ "bmaj", required_argument, 0, 'B'}, //Bmaj
	{ "bmin", required_argument, 0, 'b'}, //Bmin
	{ "bpa", required_argument, 0, 'a'}, //Pos angle
	{ "fluxtruncthr", required_argument, 0, 'f'},
	{ "threshold", required_argument, 0, 't'}, //Flux threshold below which pixels are removed from convolved sources
	{ "threshold-ext", required_argument, 0, 'T'}, //Flux threshold below which pixels are removed from convolved extended sources
	{ "mergesources", no_argument, 0, 'm'},
	{ "userthreshold", no_argument, 0, 'u'},
	{ "mergecompactsources", no_argument, 0, 'e'},
	{ "mergeextsources", no_argument, 0, 'E'},
	{ "nsigmas", required_argument, 0, 's'}, //Gaus conv kernel nsigmas (default=10)
  {(char*)0, (int)0, (int*)0, (int)0}
};


//Options
int verbosity= 4;//INFO level
std::string fileName= "";
std::string fileName_recmap= "";
bool recMapGiven= false;
std::string outputFileName= "skymodel_conv.root";
bool outputFileNameGiven= false;
std::string outputFileName_fits= "skymodel_conv.fits";
bool outputFileNameGiven_fits= false;
std::string outputFileName_ds9regions= "ds9regions.reg";
bool outputFileNameGiven_ds9regions= false;

bool bmajGiven= false;
bool bminGiven= false;
bool bpaGiven= false;
double Bmaj= -1;
double Bmin= -1;
double Bpa= -1;
double fluxTruncThr= 0.001;//max 0.1% flux loss when truncating source
bool useUserThreshold= false;
double fluxThr= -1;
double fluxThr_ext= -1;
bool mergeOverlappingSources= false;
bool enableExtendedSourceMerging= false;
bool enableCompactSourceMerging= false;
int nSigmas= 10;
int minPixels= 5;

//Vars
Image* img= 0;
std::vector<Source*> sources;
Image* img_conv= 0;
Image* img_rec= 0;
std::vector<Source*> sources_conv;//true source convolved with beam
std::vector<Source*> sources_conv_merged;//true source convolved with beam and merged if overlapping
std::vector<Source*> sources_rec;//rec source formed from true convolved source mask
TFile* outputFile= 0;
TTree* sourceTree= 0;
Source* aSource= 0;	
TFile* outputFile_rec= 0;
TTree* sourceTree_rec= 0;

//Functions
int ParseOptions(int argc, char *argv[]);
std::string GetStringLogLevel(int verbosity);
int ReadData();
int RunConvolver();
int MergeSources();
void Save();
void SaveDS9RegionFile();
void Init();

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
	Init();

	//=======================
	//== Read data
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading source data ...");
	#endif
	if(ReadData()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Reading of source data failed!");
		#endif
		return -1;
	}

	//=======================
	//== Convolver
	//=======================
	#ifdef LOGGING_ENABLED
		INFO_LOG("Convolve skymodel sources ...");
	#endif
	if(RunConvolver()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Skymodel source convolver failed!");
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
		INFO_LOG("End skymodel convolver");
	#endif

	return 0;

}//close main



int ParseOptions(int argc, char *argv[])
{
	//## Check args
	if(argc<2){
		cerr<<"ERROR: Invalid number of arguments (n="<<argc<<") ...see program usage!"<<endl;
		Usage(argv[0]);
		exit(1);
	}

	//## Parse options
	int c = 0;
  int option_index = 0;

	while((c = getopt_long(argc, argv, "hi:I:o:O:r:v:f:t:T:B:b:a:s:meEu",options_tab, &option_index)) != -1) {
    
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
				fileName_recmap= std::string(optarg);	
				recMapGiven= true;
				break;	
			}
			case 'o':	
			{
				outputFileName= std::string(optarg);	
				outputFileNameGiven= true;
				break;	
			}
			case 'O':	
			{
				outputFileName_fits= std::string(optarg);
				outputFileNameGiven_fits= true;
				break;
			}
			case 'r':	
			{
				outputFileName_ds9regions= std::string(optarg);
				outputFileNameGiven_ds9regions= true;
				break;
			}
			case 'v':	
			{
				verbosity= atoi(optarg);	
				break;	
			}
			case 'u':	
			{
				useUserThreshold= true;
				break;	
			}
			case 'f':
			{
				fluxTruncThr= atof(optarg);
				break;
			}
			case 't':
			{
				fluxThr= atof(optarg);
				break;
			}
			case 'T':
			{
				fluxThr_ext= atof(optarg);
				break;
			}
			case 'B':
			{
				Bmaj= atof(optarg);
				bmajGiven= true;
				break;
			}
			case 'b':
			{
				Bmin= atof(optarg);
				bminGiven= true;
				break;
			}
			case 'a':
			{
				Bpa= atof(optarg);
				bpaGiven= true;
				break;
			}
			case 'm':
			{
				mergeOverlappingSources= true;
				break;
			}
			case 'e':
			{
				enableCompactSourceMerging= true;
				break;
			}
			case 'E':
			{
				enableExtendedSourceMerging= true;
				break;
			}
			case 's':
			{
				nSigmas= atoi(optarg);
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
	//== Check args 
	//=======================
	//Check input file
	if(fileName==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid/empty input file name given!");
		#endif
		return -1;
	}

	//Check rec map (if given)
	if(recMapGiven){
		if(fileName_recmap==""){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Empty rec map filename given!");	
			#endif
			return -1;
		}
	}
	else{
		bool userBeamGiven= (bmajGiven && bminGiven && bpaGiven);
		if(!userBeamGiven){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("No or incomplete beam parameters given (hint: when recmap is not given as argument all beam pars must be specified!)");
			#endif
			return -1;
		}
		if(Bmin<0 || Bmaj<0 || Bpa<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Invalid beam parameters given (hint: beam pars are mandatory and shall be >0)");
			#endif
			return -1;
		}
	}//close else

	//Check nsigmas
	if(nSigmas<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid nsigma arg given (hint: shall be >0)!");
		#endif
		return -1;
	}

	return 0;

}//close ParseOptions()


int RunConvolver()
{
	//================================================
	// ==> Convolver steps
	//1) Read skymodel image (to get image size)	
	//2) Loop over true sources
	//    - Read source i-th
	//    - Create mask image with source i-th
	//    - Convolve mask image with given beam
	//    - Apply threshold to convolved image (e.g. set to zero all pixels below a given flux threshold)
	//    - Find source in thresholded convolved image and add to convolved source list
	//3) Make a mask image with all convolved sources --> write to fits & root
	//4) Merge sources sharing pixels (if enabled) --> write to root
	//=================================================


	//Copy input map
	img_conv= img->GetCloned("",true,true);
	img_conv->Reset();

	//Compute flux correction factor
	ImgMetaData* metadata= img->GetMetaData();	
	double dX= 1;
	double dY= 1;
	if(metadata){
		dX= fabs(metadata->dX*3600);//convert to arcsec
		dY= fabs(metadata->dY*3600);//convert to arcsec
	}
	else{
		#ifdef LOGGING_ENABLED
			WARN_LOG("Input map has no metadata, assuming pixel sizes=1...");
		#endif
	}
	double beamArea= AstroUtils::GetBeamAreaInPixels(Bmaj,Bmin,dX,dY);

	int source_counter= 0;

	for(size_t i=0;i<sources.size();i++){
		#ifdef LOGGING_ENABLED
			if(i%100==0) INFO_LOG("#"<<i+1<<"/"<<sources.size()<<" sources convolved...");
		#endif

		//Get true source info
		int type= sources[i]->Type;
		int simtype= sources[i]->SimType;
		double simmaxscale= sources[i]->SimMaxScale;
		int flag= sources[i]->Flag;
		bool hasTrueInfo= sources[i]->HasTrueInfo();
		if(!hasTrueInfo){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Source no. "<<i<<" has no true info stored, skip it!");
			#endif
			continue;
		}
		double S_true= sources[i]->GetTrueFlux(); 
		double X0_true, Y0_true;
		sources[i]->GetTruePos(X0_true,Y0_true);

		//Fill image source mask
		//NB: For compact source use true info and fill pixel, for extended source use pixel list
		Image* sourceImg= img->GetCloned("",true,true);
		if(!sourceImg){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to create image mask for source "<<i+1<<"!");
			#endif
			return -1;
		}
		sourceImg->Reset();

		if(type==eCompact || type==ePointLike){
			sourceImg->Fill(X0_true, Y0_true, S_true);
		}
		else if(type==eExtended || type==eCompactPlusExtended){
			std::vector<Pixel*> pixels= sources[i]->GetPixels();
			for(size_t k=0;k<pixels.size();k++){	
				long int gBin= pixels[k]->id;
				double S= pixels[k]->S;
				sourceImg->SetPixelValue(gBin,S);
			}//end loop pixels
		}//close else
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Unknown type for source no. "<<i+1<<", skip it...");
			#endif
			continue;
		}

		/*
		//This is wrong!!!!
		bool isBinary= false;	
		bool invert= false;
		Image* sourceImg= img->GetSourceMask({sources[i]},isBinary,invert);
		if(!sourceImg){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute image mask for source "<<i+1<<"!");
			#endif
			return -1;
		}
		*/

		//Convolve source image with beam
		Image* sourceImg_conv= sourceImg->GetBeamConvolvedImage(Bmaj,Bmin,Bpa,nSigmas);
		if(!sourceImg_conv){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to convolve source image "<<i+1<<"!");
			#endif
			delete sourceImg;	
			sourceImg= 0;
			return -1;
		}

		//Delete source image
		delete sourceImg;
		sourceImg= 0;

		//Find truncation threshold (user-supplied or found from image)
		double thr= 0;//no threshold
		if(useUserThreshold){
			if(type==eCompact || type==ePointLike){
				thr= fluxThr;
			}
			else if(type==eExtended || type==eCompactPlusExtended){
				thr= fluxThr_ext;
			}
			else{
				thr= fluxThr;
			}
		}//close if
		else{
			//Find threshold
			bool skipNegativePixels= true;
			thr= sourceImg_conv->FindCumulativeSumThr(fluxTruncThr,skipNegativePixels);
		}

		//Find convolved source
		std::vector<Source*> csources;
		if(sourceImg_conv->FindCompactSource(csources,thr,minPixels)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to find convolved source "<<i+1<<"!");
			#endif
			delete sourceImg_conv;	
			sourceImg_conv= 0;
			return -1;
		}
	
		if(csources.empty()){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Source "<<i+1<<" not found after convolution (below npix/flux threshold?) will be removed...");
			#endif
			delete sourceImg_conv;
			sourceImg_conv= 0;
			continue;
		}

		int csourceIndex= 0;
		if(csources.size()>1){
			#ifdef LOGGING_ENABLED
				WARN_LOG("More than 1 source found in convolved image for source no. "<<i+1<<" (name="<<sources[i]->GetName()<<", id="<<sources[i]->Id<<", type="<<sources[i]->Type<<"), this should not occur normally (could be one extended source broke up at image edge), will take the larger one...");
			#endif
			long int nPix_max= -999;
			for(size_t k=0;k<csources.size();k++){
				long int nPix= csources[k]->NPix;
				if(nPix>=nPix_max){
					nPix_max= nPix;
					csourceIndex= k;
				}
			}
			//delete sourceImg_conv;
			//sourceImg_conv= 0;
			//return -1;
		}

		//Add convolved source to list
		source_counter++;
		TString sourceName= Form("S%d",source_counter);
		csources[csourceIndex]->SetName(std::string(sourceName));	
		csources[csourceIndex]->SetId(source_counter);
		csources[csourceIndex]->Type= type;
		csources[csourceIndex]->SimType= simtype;
		csources[csourceIndex]->SimMaxScale= simmaxscale;
		csources[csourceIndex]->Flag= flag;
		csources[csourceIndex]->SetTrueInfo(S_true,X0_true,Y0_true);
		csources[csourceIndex]->SetBeamFluxIntegral(beamArea);
		sources_conv.push_back(csources[csourceIndex]);		

		//Add convolved image to skymodel
		bool computeStats= false;
		img_conv->Add(sourceImg_conv,1,computeStats);

		//Delete convolved source image
		delete sourceImg_conv;
		sourceImg_conv= 0;

	}//end loop sources

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sources_conv.size()<<" sources present after convolution...");
	#endif

	//Set metadata in skymodel convolved image
	ImgMetaData* metadata_conv= img_conv->GetMetaData();
	if(metadata_conv){
		metadata_conv->BUnit= "Jy/beam"; 
		metadata_conv->Bmaj= Bmaj/3600.;
		metadata_conv->Bmin= Bmin/3600.;
		metadata_conv->Bpa= Bpa;
	}

	//Compute stats of skymodel convolved image
	#ifdef LOGGING_ENABLED
		INFO_LOG("Computing stats of skymodel convolved image...");
	#endif
	img_conv->ComputeStats(true);

	//Merge convolved sources
	if(mergeOverlappingSources && MergeSources()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to merge convolved sources!");
		#endif
		return -1;
	}

	//If rec map is given create rec sources from conv source mask
	if(recMapGiven && img_rec){

		//Compute flux correction factor for rec sources
		ImgMetaData* metadata_rec= img_rec->GetMetaData();	
		double dX_rec= 1;
		double dY_rec= 1;
		if(metadata_rec){
			dX_rec= fabs(metadata_rec->dX*3600);//convert to arcsec
			dY_rec= fabs(metadata_rec->dY*3600);//convert to arcsec
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("Rec map has no metadata, assuming pixel sizes=1...");
			#endif
		}
		double beamArea_rec= AstroUtils::GetBeamAreaInPixels(Bmaj,Bmin,dX_rec,dY_rec);

		//Get true convolved source collection to be iterated
		std::vector<Source*>::iterator it_start= sources_conv.begin();
		std::vector<Source*>::iterator it_end= sources_conv.end();
		if(mergeOverlappingSources){
			it_start= sources_conv_merged.begin();	
			it_end= sources_conv_merged.end();	
		}
		
		//Iterate over true convolved sources
		for(std::vector<Source*>::iterator it= it_start; it != it_end; ++it) {
			Source* aSource= *it;
    	
			//Get source data
			std::string sourceName= aSource->GetName();
			long int id= aSource->Id;
			int type= aSource->Type;
			int simtype= aSource->SimType;
			int flag= aSource->Flag;
			double simmaxscale= aSource->SimMaxScale;
			double X0= aSource->X0;
			double Y0= aSource->Y0;
			double S= aSource->GetS();//in Jy/beam
			double beamIntegral= aSource->GetBeamFluxIntegral(); 
			double S_true_estimated= S/beamIntegral;//convert to Jy/pixel
			double S_true= S_true_estimated;
			if(aSource->HasTrueInfo()){
				aSource->GetTruePos(X0,Y0);
				S_true= aSource->GetTrueFlux();
			}
			std::vector<Pixel*> pixels= aSource->GetPixels();
		
			//Make rec source using pixels present in true convolved source
			Source* rec_source= new Source;
			rec_source->SetName(sourceName);
			rec_source->Id= id;
			rec_source->Type= type;	
			rec_source->Flag= flag;	
			rec_source->SimType= simtype;	
			rec_source->SimMaxScale= simmaxscale;
			//rec_source->SetTrueInfo(S_true_estimated,X0,Y0);
			rec_source->SetTrueInfo(S_true,X0,Y0);
			rec_source->SetBeamFluxIntegral(beamArea_rec);

			for(size_t k=0;k<pixels.size();k++){	
				long int gBin= pixels[k]->id;
				long int ix= pixels[k]->ix;
				long int iy= pixels[k]->iy;
				double x= pixels[k]->x;
				double y= pixels[k]->y;
				double flux= img_rec->GetPixelValue(gBin);
				rec_source->AddPixel( new Pixel(gBin,ix,iy,x,y,flux) );
			}//end loop pixels
		
			// Compute stats
			rec_source->ComputeStats();
		
			//Compute morphology parameters
			rec_source->ComputeMorphologyParams();
	
			//Add rec source to collection
			sources_rec.push_back(rec_source);

		}//end loop sources

	}//close if

	return 0;

}//close RunConvolver()


int MergeSources()
{
	//## Return if there are no sources to be merged
	if(sources_conv.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No sources to be merged, nothing to be done...");
		#endif
		return 0;
	}

	//## Fill source graph
	#ifdef LOGGING_ENABLED
		INFO_LOG("Fill list of sources to be merged and fill corresponding graph data struct...");
	#endif
	Graph mergedSourceGraph;
	std::vector<bool> isMergeableSource;
	for(size_t i=0;i<sources_conv.size();i++){
		mergedSourceGraph.AddVertex();
		isMergeableSource.push_back(false);
	}

	//## Find adjacent sources	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Finding adjacent/overlapping sources (#"<<mergedSourceGraph.GetNVertexes()<<") ...");
	#endif

	for(size_t i=0;i<sources_conv.size()-1;i++){
		Source* source= sources_conv[i];
		int type= source->Type;
		bool isCompactSource= (type==eCompact || type==ePointLike);
		bool isExtendedSource= (type==eExtended || type==eCompactPlusExtended);

		//Loop neighbors
		for(size_t j=i+1;j<sources_conv.size();j++){	
			Source* source_neighbor= sources_conv[j];
			int type_neighbor= source_neighbor->Type;
			bool isCompactSource_neighbor= (type_neighbor==eCompact || type_neighbor==ePointLike);
			bool isExtendedSource_neighbor= (type_neighbor==eExtended || type_neighbor==eCompactPlusExtended);

			//Check if both sources are compact and if they are allowed to be merged
			if(isCompactSource && isCompactSource_neighbor && !enableCompactSourceMerging){	
				#ifdef LOGGING_ENABLED		
					DEBUG_LOG("Skip merging as both sources (i,j)=("<<i<<","<<j<<") are compact and merging among compact sources is disabled...");
				#endif
				continue;
			}
	
			//Check if both sources are extended or if one is extended and the other compact and if they are allowed to be merged
			if(isExtendedSource && isExtendedSource_neighbor && !enableExtendedSourceMerging){
				#ifdef LOGGING_ENABLED			
					DEBUG_LOG("Skip merging as both sources (i,j)=("<<i<<","<<j<<") are extended and merging among extended sources is disabled...");
				#endif
				continue;
			}
			if( ((isCompactSource && isExtendedSource_neighbor) || (isExtendedSource && isCompactSource_neighbor)) && !enableExtendedSourceMerging){		
				#ifdef LOGGING_ENABLED	
					DEBUG_LOG("Skip merging between sources (i,j)=("<<i<<","<<j<<") as merging among extended and compact sources is disabled...");
				#endif
				continue;
			}
		
			//Check is sources are adjacent
			//NB: This is time-consuming (N1xN2 more or less)!!!
			bool areAdjacentSources= source->IsAdjacentSource(source_neighbor);
			if(!areAdjacentSources) continue;

			//If they are adjacent add linking in graph
			#ifdef LOGGING_ENABLED
				INFO_LOG("Sources (i,j)=("<<i<<","<<j<<") are adjacent and selected for merging...");
			#endif
			mergedSourceGraph.AddEdge(i,j);
			isMergeableSource[i]= true;
			isMergeableSource[j]= true;
		}//end loop sources
	}//end loop sources


	//## Add to merged collection all sources not mergeable
	//## NB: If sources are not to be merged each other, add to merged collection 
	sources_conv_merged.clear();

	for(size_t i=0;i<sources_conv.size();i++){
		int type= sources_conv[i]->Type;
		bool isMergeable= isMergeableSource[i];
		bool isCompactSource= (type==eCompact || type==ePointLike);
		bool isExtendedSource= (type==eExtended || type==eCompactPlusExtended);		
		if(isMergeable) {	
			if(enableCompactSourceMerging && isCompactSource) continue;
			if(enableExtendedSourceMerging && isExtendedSource) continue;
		}
		Source* merged_source= new Source;
		*merged_source= *(sources_conv[i]);
		sources_conv_merged.push_back(merged_source);
	}


	//## Find all connected components in graph corresponding to sources to be merged
	#ifdef LOGGING_ENABLED
		INFO_LOG("Find all connected components in graph corresponding to sources to be merged...");
	#endif
	std::vector<std::vector<int>> connected_source_indexes;
	mergedSourceGraph.GetConnectedComponents(connected_source_indexes);
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<connected_source_indexes.size()<<"/"<<sources_conv.size()<<" selected for merging...");
	#endif

	//## Now merge the sources
	bool copyPixels= true;//create memory for new pixels
	bool checkIfAdjacent= false;//already done before
	bool sumMatchingPixels= true;
	bool computeStatPars= false;//do not compute stats& pars at each merging
	bool computeMorphPars= false;
	bool computeRobustStats= true;
	bool forceRecomputing= true;//need to re-compute moments because pixel flux of merged sources are summed up

	#ifdef LOGGING_ENABLED	
		INFO_LOG("Merging sources and adding them to collection...");
	#endif
	for(size_t i=0;i<connected_source_indexes.size();i++){
		//Skip empty or single sources
		if(connected_source_indexes[i].size()<=1) continue;

		//Get source id=0 of this component
		int index= connected_source_indexes[i][0];
		Source* source= sources_conv[index];
		
		//Create a new source which merges the two
		Source* merged_source= new Source;
		*merged_source= *source;

		//Merge other sources in the group if any 
		int nMerged= 0;
		
		for(size_t j=1;j<connected_source_indexes[i].size();j++){
			int index_adj= connected_source_indexes[i][j];
			Source* source_adj= sources_conv[index_adj];
				
			int status= merged_source->MergeSource(source_adj,copyPixels,checkIfAdjacent,computeStatPars,computeMorphPars,sumMatchingPixels);
			if(status<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to merge sources (i,j)=("<<index<<","<<index_adj<<"), skip to next...");
				#endif
				continue;
			}
			nMerged++;

		}//end loop of sources to be merged in this component

		//If at least one was merged recompute stats & pars of merged source
		if(nMerged>0) {
			//Set name
			TString sname= Form("Smerg%d",i+1);
			merged_source->SetId(i+1);
			merged_source->SetName(std::string(sname));

			//Compute stats
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Recomputing stats & moments of merged source in merge group "<<i<<" after #"<<nMerged<<" merged source...");
			#endif
			if(merged_source->ComputeStats(computeRobustStats,forceRecomputing)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to compute stats for merged source in merge group "<<i<<"...");
				#endif
				continue;
			}
	
			//Compute morph params
			if(merged_source->ComputeMorphologyParams()<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to compute morph pars for merged source in merge group "<<i<<"...");
				#endif
				continue;
			}
		}//close if

		//Add merged source to collection
		sources_conv_merged.push_back(merged_source);

	}//end loop number of components

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sources_conv_merged.size()<<" sources present in merged collection...");
	#endif

	return 0;

}//close MergeSources()



void Init()
{
	//Open output file
	if(!outputFile) {	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Opening ROOT file "<<outputFileName<<" ...");
		#endif
		outputFile= new TFile(outputFileName.c_str(),"RECREATE");
	}

	//Create source Tree
	if(!sourceTree) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("Creating ROOT Tree SourceInfo ...");
		#endif
		sourceTree= new TTree("SourceInfo","SourceInfo");	
	}
	aSource= 0;
	sourceTree->Branch("Source",&aSource);

	//Define output FITS file name
	std::string outputFileName_base= CodeUtils::ExtractSubString(outputFileName,".");
	if(!outputFileNameGiven_fits){		
		outputFileName_fits= outputFileName_base + std::string(".fits");
		#ifdef LOGGING_ENABLED
			INFO_LOG("Set skymodel convolved FITS output to "<<outputFileName_fits<<" ...");
		#endif
	}

	//Define output FITS file name
	if(!outputFileNameGiven_ds9regions){
		outputFileName_ds9regions= outputFileName_base + std::string(".reg");	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Set DS9 output regions file to "<<outputFileName_ds9regions<<" ...");
		#endif
	}

	//Create file & source Tree with rec sources (if rec map given)
	if(recMapGiven){
		std::string outputFileName_rec= outputFileName_base + std::string("_rec.root");
		if(!outputFile_rec) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Opening ROOT file "<<outputFileName_rec<<" ...");
			#endif
			outputFile_rec= new TFile(outputFileName_rec.c_str(),"RECREATE");
		}

		if(!sourceTree_rec) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Creating ROOT Tree SourceInfo for rec sources...");
			#endif
			sourceTree_rec= new TTree("SourceInfo","SourceInfo");	
		}
		aSource= 0;
		sourceTree_rec->Branch("Source",&aSource);
	}


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

int ReadData()
{
	//Read recmap and get beam info (if option is given)
	if(recMapGiven){
		img_rec= new Image();
		if(img_rec->ReadFITS(fileName_recmap)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to read rec map FITS file "<<fileName_recmap<<"!");
			#endif
			return -1;
		}
		
		//Get beam info
		ImgMetaData* metadata= img_rec->GetMetaData();	
		if(!metadata){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Rec map has no metadata!");
			#endif
			delete img_rec;
			img_rec= 0;
			return -1;
		}
		Bmaj= metadata->Bmaj*3600;
		Bmin= metadata->Bmin*3600;
		Bpa= metadata->Bpa;
		#ifdef LOGGING_ENABLED
			INFO_LOG("Read beam info from recmap: {Bmaj,Bmin,Bpa}={"<<Bmaj<<","<<Bmin<<","<<Bpa<<"}");
		#endif
	}//close if

	
	//Open file with source collection
	TFile* inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<fileName<<"!");
		#endif
		return -1;
	}

	//Read skymodel image
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading skymodel image from file "<<fileName<<"...");
	#endif
	img= (Image*)inputFile->Get("img");
	if(!img){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read skymodel image from input file "<<fileName<<"!");
		#endif
		return -1;
	}
	
	//Get access to source trees
	#ifdef LOGGING_ENABLED
		INFO_LOG("Get access to source tree from file "<<fileName<<"...");
	#endif

	TTree* sourceDataTree= (TTree*)inputFile->Get("SourceInfo");
	if(!sourceDataTree || sourceDataTree->IsZombie()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get access to source tree in file "<<fileName<<"!");	
		#endif
		return -1;
	}
	sourceDataTree->SetBranchAddress("Source",&aSource);

	//Read sources
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading #"<<sourceDataTree->GetEntries()<<" sources in file "<<fileName<<"...");
	#endif
	for(int i=0;i<sourceDataTree->GetEntries();i++){
		sourceDataTree->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<sources.size()<<" sources read...");
	#endif

	
	return 0;

}//close ReadData()

void Save()
{
	//Write DS9 regions
	SaveDS9RegionFile();
	
	//Write fits image
	if(img_conv){	
		#ifdef LOGGING_ENABLED
			INFO_LOG("Saving convolved skymodel to FITS file...");
		#endif
		if(img_conv->WriteFITS(outputFileName_fits)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to write skymodel convolved to FITS file "<<outputFileName_fits<<"!");
			#endif
		}
	}

	//Save data to ROOT file
	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();		

		//Save convolved skymodel image
		if(img_conv) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Saving convolved skymodel to ROOT file...");
			#endif
			img_conv->SetName("img");
			img_conv->Write(); 
		}

		//Save source tree?
		if(sourceTree){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Filling source ROOT TTree...");
			#endif
			if(mergeOverlappingSources){
				for(size_t k=0;k<sources_conv_merged.size();k++){
					aSource= sources_conv_merged[k];
					sourceTree->Fill();
				}
			}//close if
			else{
				for(size_t k=0;k<sources_conv.size();k++){
					aSource= sources_conv[k];
					sourceTree->Fill();
				}
			}
			#ifdef LOGGING_ENABLED
				INFO_LOG("Writing tree to file...");
			#endif
			sourceTree->Write();
		}//close if save source tree

		outputFile->Close();
	}//close if file open

	//Save rec source data to ROOT file
	if(recMapGiven && outputFile_rec && outputFile_rec->IsOpen()){
		outputFile_rec->cd();		

		//Save source tree?
		if(sourceTree_rec){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Filling rec source ROOT TTree...");
			#endif
			for(size_t k=0;k<sources_rec.size();k++){
				aSource= sources_rec[k];
				sourceTree_rec->Fill();
			}
			#ifdef LOGGING_ENABLED
				INFO_LOG("Writing source rec tree to file...");
			#endif
			sourceTree_rec->Write();
		}//close if save source tree

		outputFile_rec->Close();
	}//close if file open

}//close Save()


void SaveDS9RegionFile()
{
	#ifdef LOGGING_ENABLED
		INFO_LOG("Saving "<<sources_conv_merged.size()<<" sources to DS9 region file "<<outputFileName_ds9regions<<" ...");
	#endif

	//Open file
	FILE* fout= fopen(outputFileName_ds9regions.c_str(),"w");

	//Writing header
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Saving DS9 region header...");
	#endif
	fprintf(fout,"global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"image\n");
	
	//Writing source regions	
	std::string colorStr_last= "white";
	int nSources= static_cast<int>(sources_conv.size());
	if(mergeOverlappingSources) nSources= static_cast<int>(sources_conv_merged.size());

	for(int k=0;k<nSources;k++){
		Source* aSource= 0;	
		if(mergeOverlappingSources) aSource= sources_conv_merged[k];
		else aSource= sources_conv[k];

		int source_type= aSource->Type;
		bool isAtEdge= aSource->IsAtEdge();

		//Set source color
		std::string colorStr= "white";
		if(source_type==eExtended) colorStr= "green";
		else if(source_type==eCompactPlusExtended) colorStr= "orange";
		else if(source_type==ePointLike) colorStr= "red";
		else if(source_type==eCompact) colorStr= "blue";
		else colorStr= "magenta";
		if(colorStr!=colorStr_last){
			colorStr_last= colorStr;
			fprintf(fout,"global color=%s font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n",colorStr.c_str());
		}

		//Get region and write to file
		std::string regionInfo= aSource->GetDS9Region(true);
		fprintf(fout,"%s\n",regionInfo.c_str());
	  	
	}//end loop sources
		
	//Close file
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Closing DS9 file region...");
	#endif
	fclose(fout);

}//close SaveDS9RegionFile()


