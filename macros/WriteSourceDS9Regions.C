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

//Vars
bool useDashLines= false;
bool useColors= true;

//Functions
int WriteRegions(std::string);

int WriteSourceDS9Regions(std::string _fileName,bool isFileList=false,bool _useDashLines=false,bool _useColors=true){

	//## Set args
	useDashLines= _useDashLines;
	useColors= _useColors;

	//## Set logging level
	LoggerManager::Instance().CreateConsoleLogger("INFO","logger","System.out");
	
	
	//## Analyze data
	if(isFileList){//file list

		//Check filenames
		std::ifstream fileStream(_fileName);	
		std::string line;
		if (fileStream.fail() || !fileStream.is_open()){
			ERROR_LOG("Failed to open file "<<_fileName<<" for reading...");
			return -1;
		}

		//Store filenames present in lists
		std::string filename= "";
		std::vector<std::string> fileNames;
		while (std::getline(fileStream, line)) {
    	std::istringstream iss(line);
    	if (!(iss >> filename)) { 
				ERROR_LOG("Failed to read line from file "<<_fileName<<"!");
				return -1; 
			}
    	fileNames.push_back(filename);
		}//end file read
				
		//Write regions
		for(size_t i=0;i<fileNames.size();i++){
			if(WriteRegions(fileNames[i])<0){
				ERROR_LOG("Failed to write region for file no. "<<i+1);
				return -1;
			}
		}//end loop files

	}//close if
	else{//single file
		//WriteRegions
		if(WriteRegions(_fileName)<0){
			ERROR_LOG("Failed to write region file!");
			return -1;
		}
	}//close if

	
	return 0;

}//close macro


int WriteRegions(std::string fileName) 
{
	//Check file
	FileInfo info;
	bool match_extension= true;
	std::string extension= ".root";
	if(!SysUtils::CheckFile(fileName,info,match_extension,extension)){
		ERROR_LOG("Failed to check file "<<fileName<<" on filesystem!");
		return -1;
	}

	std::string fileName_noext= CodeUtils::ExtractFileNameFromPath(fileName,true);
	std::string DS9CatalogFileName= std::string("ds9regions_") + fileName_noext + std::string(".reg");
	INFO_LOG("Writing DS9 region file "<<DS9CatalogFileName<<"...");
	

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
		return -1;
	}

	Source* aSource= 0;
	sourceDataTree->SetBranchAddress("Source",&aSource);

	//Open DS9 file region and write header
	

	FILE* fout= fopen(DS9CatalogFileName.c_str(),"w");
	fprintf(fout,"global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"image\n");

	INFO_LOG("Reading #"<<sourceDataTree->GetEntries()<<" sources in file "<<fileName<<"...");
	std::string colorStr_last= "white";
	
	for(int i=0;i<sourceDataTree->GetEntries();i++){
		sourceDataTree->GetEntry(i);

		//Get source info
		int source_type= aSource->Type;
		bool isAtEdge= aSource->IsAtEdge();

		//Set source color/tag/...
		std::string colorStr= "white";
		std::string tagStr= "unknown";
		if(source_type==Source::eExtended) {
			colorStr= "green";
			tagStr= "extended";
		}
		else if(source_type==Source::eCompactPlusExtended) {
			colorStr= "magenta";
			tagStr= "extended-compact";
		}
		else if(source_type==Source::ePointLike) {
			colorStr= "red";
			tagStr= "point-like";
		}
		else if(source_type==Source::eCompact) {
			colorStr= "blue";
			tagStr= "compact";
		}
		else {
			colorStr= "white";
			tagStr= "unknown";
		}
		if(!useColors) colorStr= "white";

		//Get source region
		std::string regionInfo= aSource->GetDS9Region(true);
		regionInfo+= std::string(" color=") + colorStr;

		//Set source color/tag and other properties
		regionInfo+= std::string(" color=") + colorStr;
		regionInfo+= std::string(" tag={") + tagStr + std::string("}");
		if(useDashLines) regionInfo+= std::string(" dash=1");

		//Write source region to file
		fprintf(fout,"%s\n",regionInfo.c_str());
	  	
	}//end loop sources

	//Close file		
	fclose(fout);

	return 0;

}//close WriteRegions()

