#include <TROOT.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLine.h>

#include <Image.h>
#include <BkgData.h>
#include <Logger.h>
#include <Source.h>
#include <MathUtils.h>
#include <AstroUtils.h>
#include <StatsUtils.h>

using namespace Caesar;

//Bkg options
//int bkgEstimator= eMedianClippedBkg;
int bkgEstimator= eMedianBkg;
int bkgBoxBeamFactor= 30;
double bkgGridStepSize= 0.2;
bool use2ndPass= true;
bool skipOutliers= true;
double seedThr= 5;
double mergeThr= 2.5;
int minPixels= 5;

//Vars
std::string outputFileName= "Output.root";
TFile* outputFile= 0;
TTree* outputTree= 0;
std::string filename= "";
std::string filename_rec= "";
std::string SourceName;
double Npix;
int type;
int simtype;
double simmaxscale;
double X0;
double Y0;
double S;
double Smax;
double S_true;
double S_true_estimated;
double Npix_rec;
double X0_rec;
double Y0_rec;
double S_rec;
double Smax_rec;
bool correctRecFlux=false;
double beamArea;
double BkgAvg;
double S_bkg;
double S_bkg_mean;

//Methods
void Init();
int AnalyzeData(std::string,std::string);
void CompareData();
void Save();

int ComputeImagingPerformances(std::string _fileName,std::string _fileName_rec,bool isFileList=false,bool correctRecFluxByBeamArea=false){

	//## Set args
	correctRecFlux= correctRecFluxByBeamArea;

	//## Set logging level
	LoggerManager::Instance().CreateConsoleLogger("INFO","logger","System.out");
	
	//## Init data
	Init();


	//## Analyze data
	if(isFileList){//file list

		//Check filenames
		std::ifstream fileStream(_fileName);	
		std::ifstream fileStream_rec(_fileName_rec);
		std::string line;
		if (fileStream.fail() || !fileStream.is_open()){
			ERROR_LOG("Failed to open file "<<_fileName<<" for reading...");
			return -1;
		}
		if (fileStream_rec.fail() || !fileStream_rec.is_open()){
			ERROR_LOG("Failed to open file "<<_fileName_rec<<" for reading...");
			return -1;
		}

		//Store filenames present in lists
		std::vector<std::string> fileNames;
		while (std::getline(fileStream, line)) {
    	std::istringstream iss(line);
    	if (!(iss >> filename)) { 
				ERROR_LOG("Failed to read line from file "<<_fileName<<"!");
				return -1; 
			}
    	fileNames.push_back(filename);
		}//end file read
				
		std::vector<std::string> fileNames_rec;
		while (std::getline(fileStream_rec, line)) {
    	std::istringstream iss(line);
    	if (!(iss >> filename)) { 
				ERROR_LOG("Failed to read line from file "<<_fileName<<"!");
				return -1; 
			}
    	fileNames_rec.push_back(filename);
		}//end file read

		//Check files have same number of entries
		if(fileNames.size()!=fileNames_rec.size()){
			ERROR_LOG("Input filelist must have the same number of entries!");
			return -1;
		}

		//Finally analyze data
		for(size_t i=0;i<fileNames.size();i++){
			if(AnalyzeData(fileNames[i],fileNames_rec[i])<0){
				ERROR_LOG("Failed to analyze data for file no. "<<i+1);
				return -1;
			}
		}//end loop files

	}//close if
	else{//single file
		//Analyze data
		if(AnalyzeData(_fileName,_fileName_rec)<0){
			ERROR_LOG("Failed to analyze data!");
			return -1;
		}
	}//close if

	
	//## Save data
	Save();

	return 0;

}//close macro


int AnalyzeData(std::string fileName,std::string fileName_rec){

	filename= fileName;
	filename_rec= fileName_rec;

	//## Read rec image
	INFO_LOG("Reading rec image from file "<<fileName_rec<<"...");
	Image* img_rec= new Image();
	if(img_rec->ReadFITS(fileName_rec)<0){
		ERROR_LOG("Failed to read rec image from input file "<<fileName_rec<<"!");
		return -1;
	}

	//Get beam area
	ImgMetaData* metadata= img_rec->GetMetaData();
	int pixelWidthInBeam= 1;
	if(metadata){
		beamArea= metadata->GetBeamFluxIntegral();
		pixelWidthInBeam= metadata->GetBeamWidthInPixel();
	}
	else{
		ERROR_LOG("Rec map has no metadata stored, cannot get beam information needed for correction!");
		beamArea= 1;
		delete img_rec;
		img_rec= 0;
		return -1;
	}
	INFO_LOG("Beam area="<<beamArea<<", pixelWidthInBeam="<<pixelWidthInBeam);


	//## Read true source data
	//Open file with source collection
	TFile* inputFile= new TFile(fileName.c_str(),"READ");
	if(!inputFile){
		ERROR_LOG("Failed to open file "<<fileName<<"!");
		return -1;
	}

	//Get access to source trees
	INFO_LOG("Get access to source tree from file "<<fileName<<"...");
	
	TTree* sourceDataTree= (TTree*)inputFile->Get("SourceInfo");
	if(!sourceDataTree || sourceDataTree->IsZombie()){
		ERROR_LOG("Failed to get access to source tree in file "<<fileName<<"!");	
		if(img_rec){
			delete img_rec;
			img_rec= 0;
		}
		return -1;
	}

	Source* aSource= 0;
	sourceDataTree->SetBranchAddress("Source",&aSource);

	//Read sources
	std::vector<Source*> sources;
	INFO_LOG("Reading #"<<sourceDataTree->GetEntries()<<" sources in file "<<fileName<<"...");
	for(int i=0;i<sourceDataTree->GetEntries();i++){
		sourceDataTree->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		sources.push_back(source);
	}//end loop sources

	INFO_LOG("#"<<sources.size()<<" true sources read...");

	//## Compute bkg map in reconstructed data
	INFO_LOG("Computing image stats...");
	img_rec->ComputeStats(true);

	INFO_LOG("Computing image bkg");
	double bkgBoxSize= pixelWidthInBeam*bkgBoxBeamFactor;
	double bkgGridSize= bkgGridStepSize*bkgBoxSize;
	ImgBkgData* bkgData= img_rec->ComputeBkg(
		bkgEstimator,true,bkgBoxSize,bkgBoxSize,bkgGridSize,bkgGridSize,
		use2ndPass,skipOutliers,seedThr,mergeThr,minPixels
	);
	if(!bkgData){
		ERROR_LOG("Failed to compute rec image bkg!");
		if(img_rec){
			delete img_rec;
			img_rec= 0;
		}
		for(size_t i=0;i<sources.size();i++){
			if(sources[i]){
				delete sources[i];		
				sources[i]= 0;
			}
		}
		sources.clear();
		return -1;
	}

	//## Compute bkg stats
	(bkgData->BkgMap)->ComputeStats(true);
	BkgAvg= (bkgData->BkgMap)->GetPixelStats()->median;

	//## Compare sources
	for(size_t i=0;i<sources.size();i++){
		if(i%100==0) INFO_LOG("#"<<i+1<<"/"<<sources.size()<<" true sources read...");

		//Get source data
		Npix= sources[i]->GetNPixels();
		SourceName= sources[i]->GetName();	
		type= sources[i]->Type;
		simtype= sources[i]->SimType;
		int flag= sources[i]->Flag;
		simmaxscale= sources[i]->SimMaxScale;
		X0= sources[i]->X0;
		Y0= sources[i]->Y0;
		S= sources[i]->GetS();//in Jy/beam
		Smax= sources[i]->GetSmax();
		S_true= sources[i]->GetTrueFlux();
		double beamIntegral= sources[i]->GetBeamFluxIntegral(); 
		S_true_estimated= S/beamIntegral;//convert to Jy/pixel
		std::vector<Pixel*> pixels= sources[i]->GetPixels();
	

		//Get rec source 
		Source* rec_source= new Source;
		Pixel* aPixel= 0;
		std::vector<double> bkgValues;
		S_bkg= 0;
		for(size_t k=0;k<pixels.size();k++){	
			long int gBin= pixels[k]->id;
			long int ix= pixels[k]->ix;
			long int iy= pixels[k]->iy;
			double x= pixels[k]->x;
			double y= pixels[k]->y;
			double flux= img_rec->GetPixelValue(gBin);
			double bkgLevel= (bkgData->BkgMap)->GetPixelValue(gBin);
			double bkgRMS= (bkgData->NoiseMap)->GetPixelValue(gBin);
			S_bkg+= bkgLevel;

			bkgValues.push_back(bkgLevel);

			aPixel= new Pixel(gBin,ix,iy,x,y,flux);
			aPixel->SetBkg(bkgLevel,bkgRMS);
			rec_source->AddPixel(aPixel);
		}

		//## Compute stats
		DEBUG_LOG("Computing source stats...");
		rec_source->ComputeStats();
		
		//## Compute morphology parameters
		DEBUG_LOG("Computing blob morphology params...");
		rec_source->ComputeMorphologyParams();

		//Get rec source pars
		Npix_rec= rec_source->GetNPixels();
		X0_rec= rec_source->X0;
		Y0_rec= rec_source->Y0;
		S_rec= rec_source->GetS();
		Smax_rec= rec_source->GetSmax();
		S_bkg_mean= StatsUtils::GetMedianFast(bkgValues);		

		//Correct flux from Jy/beam to Jy
		if(correctRecFlux){
			S_rec/= beamArea;
			//Smax_rec/= beamArea;
		}
		
		//Store/compare data
		INFO_LOG("True source no. "<<i+1<<" (type="<<type<<", pos("<<X0<<","<<Y0<<"), S="<<S<<"), Rec source: pos("<<X0_rec<<","<<Y0_rec<<"), S_rec="<<S_rec<<", BkgAvg="<<BkgAvg<<", S_bkg="<<S_bkg<<")");
		outputTree->Fill();
	
	}//end loop sources

	//## Clear stuff
	if(img_rec){
		delete img_rec;
		img_rec= 0;
	}
	for(size_t i=0;i<sources.size();i++){
		if(sources[i]){
			delete sources[i];		
			sources[i]= 0;
		}
	}
	sources.clear();
	
	if(bkgData){
		delete bkgData;
		bkgData= 0;
	}

	return 0;

}//close AnalyzeData()


void Init(){

	if(!outputFile) outputFile= new TFile(outputFileName.c_str(),"RECREATE");

	if(!outputTree) outputTree= new TTree("data","data");
	outputTree->Branch("filename",&filename);
	outputTree->Branch("Npix",&Npix,"Npix/D");
	outputTree->Branch("name",&SourceName);
	outputTree->Branch("type",&type,"type/I");
	outputTree->Branch("simtype",&simtype,"simtype/I");
	outputTree->Branch("simmaxscale",&simmaxscale,"simmaxscale/D");
	outputTree->Branch("S_true",&S_true,"S_true/D");
	outputTree->Branch("S",&S,"S/D");
	outputTree->Branch("Smax",&Smax,"Smax/D");	
	outputTree->Branch("X0",&X0,"X0/D");
	outputTree->Branch("Y0",&Y0,"Y0/D");

	outputTree->Branch("filename_rec",&filename_rec);
	outputTree->Branch("Npix_rec",&Npix_rec,"Npix_rec/D");
	outputTree->Branch("S_rec",&S_rec,"S_rec/D");
	outputTree->Branch("Smax_rec",&Smax_rec,"Smax_rec/D");
	outputTree->Branch("X0_rec",&X0_rec,"X0_rec/D");
	outputTree->Branch("Y0_rec",&Y0_rec,"Y0_rec/D");
	outputTree->Branch("beamArea",&beamArea,"beamArea/D");
	outputTree->Branch("BkgAvg",&BkgAvg,"BkgAvg/D");
	outputTree->Branch("S_bkg",&S_bkg,"S_bkg/D");
	outputTree->Branch("S_bkg_mean",&S_bkg_mean,"S_bkg_mean/D");

}//close Init()

void Save(){

	if(outputFile && outputFile->IsOpen()){
		outputFile->cd();	
		if(outputTree) outputTree->Write();	
		outputFile->Close();
	}

}//close Save()


