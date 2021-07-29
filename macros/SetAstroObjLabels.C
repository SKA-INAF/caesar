#include <AstroObject.h>
#include <AstroObjectParser.h>

using namespace std;
using namespace Caesar;

//Classes
struct LabelData 
{
	std::string sname;
	int componentId;
	std::string objName;
	int objType;
	std::string objTypeStr;
	bool isConfirmed;
	double ra;
	double dec;
};

//Variables
std::string gInputFileName= "";
std::string gLabelFileName= "";
Source* gSource= 0;
std::vector<Source*> gSources;
std::map<std::string,std::vector<LabelData>> gLabelDataCollection;
TFile* gInputFile= 0;
TTree* gOutputTree= 0;
TFile* gOutputFile= 0;
std::string gOutputFileName= "output.root";
std::string gRegionOutFileName= "ds9.reg";
std::string gCatalogOutFileName= "catalog.dat";
std::string gRegionComponentsOutFileName= "";
std::string gCatalogComponentsOutFileName= "";
char gDelimiter= '\t';

//Functions
int Init();
int ReadSources();
int ReadLabelData();
int LabelSources();
int LabelSource(Source* source,std::vector<LabelData>& labelDataList);
int SaveDS9Regions();
int SaveSources();
int SaveCatalog();
int CloneObjectsInFile(std::vector<std::string> excludedObjNames);
void Save();

int SetAstroObjLabels(std::string inputFileName,std::string labelFileName="",std::string outputFileName="output.root",std::string regionOutFileName="ds9.reg",std::string catalogOutFileName="catalog.dat",char delimiter='|')
{
	//===================================
	//==       GET ARGS
	//===================================
	gInputFileName= inputFileName;
	gLabelFileName= labelFileName;
	gOutputFileName= outputFileName;		
	gRegionOutFileName= regionOutFileName;
	gCatalogOutFileName= catalogOutFileName;
	gDelimiter= delimiter;
	
	//===================================
	//==       INIT
	//===================================
	cout<<"INFO: Initializing macro data..."<<endl;
	Init();

	//===================================
	//==       READ LABEL FILE
	//===================================
	if(ReadLabelData()<0){
		cerr<<"ERROR: Failed to read labels!"<<endl;
		return -1;
	}

	//===================================
	//==       READ CATALOG SOURCES
	//===================================
	if(ReadSources()<0){
		cerr<<"ERROR: Failed to read sources from catalog!"<<endl;
		return -1;
	}

	//===================================
	//==       LABEL CATALOG SOURCES
	//===================================
	if(LabelSources()<0){
		cerr<<"ERROR: Failed to label sources!"<<endl;
		return -1;
	}

	//===================================
	//==      SAVE
	//===================================
	Save();

	return 0;

}//close macro


int Init()
{
	//Check output file name
	if(gOutputFileName==""){
		cerr<<"ERROR: Empty output file name given!"<<endl;
		return -1;
	}

	//Check region output file name
	if(gRegionOutFileName==""){
		cerr<<"ERROR: Empty region output file name given!"<<endl;
		return -1;
	}

	//Check catalog output file name
	if(gCatalogOutFileName==""){
		cerr<<"ERROR: Empty catalog output file name given!"<<endl;
		return -1;
	}

	//Set DS9 region component file name
	gRegionComponentsOutFileName= CodeUtils::ExtractFileNameFromPath(gRegionOutFileName,true);
	gRegionComponentsOutFileName+= "_fitcomp.reg";

	//Set catalog component file name
	gCatalogComponentsOutFileName= CodeUtils::ExtractFileNameFromPath(gCatalogOutFileName,true);
	gCatalogComponentsOutFileName+= "_fitcomp.dat";

	
	//Initialize file
	gOutputFile= new TFile(gOutputFileName.c_str(),"RECREATE");
	gOutputFile->cd();
	
	//Init output source TTree
	gOutputTree= new TTree("SourceInfo","SourceInfo");
	gSource= 0;
	gSources.clear();
	gOutputTree->Branch("Source",&gSource);

	return 0;

}//close Init()

int ReadLabelData()
{
	//Read label data
	cout<<"INFO: Reading label data "<<gLabelFileName<<" ..."<<endl;
	char delimiter= '|';
	TTree* labelData= new TTree;
	labelData->ReadFile(gLabelFileName.c_str(),"sname/C:componentId/I:objType/C:objName/C:isConfirmed/I:objRA/D:objDec/D",delimiter);

	char sname[100];
	int componentId;
	char objType[100];
	char objName[100];
	int isConfirmed;
	double ra;
	double dec;

	labelData->SetBranchAddress("sname",sname);
	labelData->SetBranchAddress("objType",objType);
	labelData->SetBranchAddress("objName",objName);
	labelData->SetBranchAddress("componentId",&componentId);
	labelData->SetBranchAddress("isConfirmed",&isConfirmed);
	labelData->SetBranchAddress("objRA",&ra);
	labelData->SetBranchAddress("objDec",&dec);

	cout<<"# "<<labelData->GetEntries()<<" astro object labels read from file "<<gLabelFileName<<" ..."<<endl;

	//Create astro objects
	LabelData* ld= 0;

	for(int i=0;i<labelData->GetEntries();i++)
	{
		labelData->GetEntry(i);

		std::string objTypeStr= std::string(objType);
		std::string sourceName= std::string(sname);
		std::string objNameStr= std::string(objName);

		ld= new LabelData;
		ld->sname= sourceName;
		ld->componentId= componentId;
		ld->objName= objNameStr;
		ld->objType= GetAstroObjectType(objTypeStr);
		ld->objTypeStr= objTypeStr;
		ld->isConfirmed= isConfirmed;
		ld->ra= ra;
		ld->dec= dec;

		gLabelDataCollection[sourceName].push_back(*ld);

	}//end loop

	return 0;

}//close ReadLabelData()

int ReadSources()
{
	//- Read input file
	gInputFile= new TFile(gInputFileName.c_str(),"READ");
	if(!gInputFile || gInputFile->IsZombie()) {
		cerr<<"ERROR: Cannot open file "<<gInputFileName<<"!"<<endl;
		return -1;
	}
	
	//- Get source tree
	TTree* SourceInfo= (TTree*)gInputFile->Get("SourceInfo");	
	if(!SourceInfo){
		cerr<<"ERROR: Cannot get sources from file "<<gInputFileName<<"!"<<endl;
		return -1;
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	cout<<"# "<<SourceInfo->GetEntries()<<" sources to be read..."<<endl;
	
	//- Store source list 
	for(int i=0;i<SourceInfo->GetEntries();i++)
	{
		SourceInfo->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;
		gSources.push_back(source);

	}//end loop sources

	return 0;

}//close ReadSources()


int LabelSources()
{
	//Loop over sources
	for(size_t i=0;i<gSources.size();i++)
	{
		Source* source= gSources[i];
		std::string sname= source->GetName();
		std::vector<Source*> nestedSources= source->GetNestedSources();

		//Label other source
		auto it = gLabelDataCollection.find(sname);
		bool hasLabelData= (it!=gLabelDataCollection.end());

		if(hasLabelData){
			std::vector<LabelData> labelDataList= it->second;
			if(LabelSource(source,labelDataList)<0){
				cerr<<"WARN: Failed to label source "<<sname<<"!"<<endl;
			}	
		}

		//Label nested sources
		for(size_t j=0;j<nestedSources.size();j++)
		{
			Source* nestedSource= nestedSources[j];
			std::string sname_nested= nestedSource->GetName();
			auto it_nested = gLabelDataCollection.find(sname_nested);
			bool hasLabelData_nested= (it_nested!=gLabelDataCollection.end());

			if(!hasLabelData_nested) continue;

			std::vector<LabelData> labelDataList= it_nested->second;
			if(LabelSource(nestedSource,labelDataList)<0){
				cerr<<"WARN: Failed to label source "<<sname_nested<<"!"<<endl;
			}

		}//end loop

	}//end loop sources

	return 0;

}//close LabelSources()

int LabelSource(Source* source,std::vector<LabelData>& labelDataList)
{
	//Check ptr
	if(!source) return -1;

	cout<<"INFO: Labeling source "<<source->GetName()<<" ..."<<endl;

	AstroObject* astroObj= 0;
	std::vector<AstroObject*> astroObjs;

	for(size_t i=0;i<labelDataList.size();i++)
	{	
		bool hasFitInfo= source->HasFitInfo();
		int nFitComponents= source->GetNFitComponents();

		std::string objName= labelDataList[i].objName;
		//CodeUtils::StripBlankSpaces(objName);

		//Create astro object
		astroObj= new AstroObject;
		astroObj->index= i;
		astroObj->name= objName;
		astroObj->id_str= labelDataList[i].objTypeStr;
		astroObj->id= labelDataList[i].objType;
		astroObj->subid= labelDataList[i].objType;
		astroObj->x= labelDataList[i].ra;
		astroObj->y= labelDataList[i].dec;
		astroObj->confirmed= labelDataList[i].isConfirmed;

		//Label source
		source->AddAstroObject(*astroObj);

		int compId= labelDataList[i].componentId-1;
		if(labelDataList[i].componentId>=0){
			if(source->AddComponentAstroObject(compId,*astroObj)<0){
				cerr<<"WARN: Failed to add astro obj for component "<<compId<<" of source "<<source->GetName()<<" (nFitComponents="<<nFitComponents<<")"<<endl;
			}
		}
		else{
			for(int k=0;k<nFitComponents;k++){
				if(source->AddComponentAstroObject(k,*astroObj)<0){
					cerr<<"WARN: Failed to add astro obj for component "<<k<<" of source "<<source->GetName()<<" (nFitComponents="<<nFitComponents<<")"<<endl;
				}
			}
		}

	}//end loop labels

	return 0;

}//close LabelSource()



int SaveDS9Regions()
{
	//Save DS9 regions for islands
	bool convertDS9RegionsToWCS= false;
	int status= SourceExporter::WriteToDS9(gRegionOutFileName,gSources,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	//Save DS9 regions for components
	status= SourceExporter::WriteComponentsToDS9(gRegionComponentsOutFileName,gSources,convertDS9RegionsToWCS);
	if(status<0){
		return -1;
	}

	return 0;

}//close SaveDS9Regions()


int SaveSources()
{
	//Check output file is open
	if(!gOutputFile || !gOutputFile->IsOpen()){
		cerr<<"ERROR: Output ROOT file not allocated or open!"<<endl;
		return -1;
	}
	if(!gOutputTree){
		cerr<<"ERROR: Output ROOT source TTree not allocated!"<<endl;
		return -1;
	}

	//Loop over selected sources and write TTree to output file
	for(size_t i=0;i<gSources.size();i++){
		gSource= gSources[i];
		gOutputTree->Fill();
	}

	gOutputTree->Write();

	return 0;

}//close SaveSources()



int SaveCatalog()
{
	//Return if no sources are found
	if(gSources.empty()){
		cerr<<"WARN: No sources selected, no catalog file will be written!"<<endl;
		return 0;
	}	

	//Saving island/blob catalog to ascii file
	cout<<"INFO: Writing source catalog to file "<<gCatalogOutFileName<<" ..."<<endl;
	bool dumpNestedSourceInfo= true;
	int wcsType= eJ2000;
	WCS* wcs= 0;
	bool writeAdditionalSourceInfo= true;
	bool convertBrightenessToFlux= true;

	int status= SourceExporter::WriteToAscii(gCatalogOutFileName,gSources,dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightenessToFlux,gDelimiter);
	if(status<0){
		cerr<<"WARN: Writing source catalog to file "<<gCatalogOutFileName<<" failed!"<<endl;
	}
	
	//Saving source fitted components to ascii file
	cout<<"INFO: Writing source catalog to file "<<gCatalogComponentsOutFileName<<" ..."<<endl;
	
	status= SourceExporter::WriteComponentsToAscii(gCatalogComponentsOutFileName,gSources,dumpNestedSourceInfo,wcsType,wcs,writeAdditionalSourceInfo,convertBrightenessToFlux,gDelimiter);
	if(status<0){
		cerr<<"WARN: Writing source fitted component catalog to file "<<gCatalogComponentsOutFileName<<" failed!"<<endl;
	}

	return 0;

}//close SaveCatalog()


void Save()
{
	//Save TTree to file
	if(gOutputFile && gOutputFile->IsOpen()){
		gOutputFile->cd();		
		
		//Clone all objects present in input file but the Source TTree
		CloneObjectsInFile({"SourceInfo"});

		//Write selected source TTree
		SaveSources();

		//Close file
		gOutputFile->Close();

	}//close if

	//Save DS9 regions with selected sources	
	cout<<"INFO: Saving DS9 regions with selected sources..."<<endl;
	SaveDS9Regions();

	//Save ascii catalog with selected sources
	cout<<"INFO: Saving ascii catalog with selected sources..."<<endl;
	SaveCatalog();

}//close Save()

int CloneObjectsInFile(std::vector<std::string> excludedObjNames)
{	
	//Check if input file is open
	if(!gInputFile || !gInputFile->IsOpen()){
		cerr<<"ERROR: Input file is not open!"<<endl;
		return -1;
	}

	//Check if output file is open
	if(!gOutputFile || !gOutputFile->IsOpen()){
		cerr<<"ERROR: Output file is not open!"<<endl;
		return -1;
	}

	//Loop on all entries present in input file
	TKey* key;
  TIter nextkey(gInputFile->GetListOfKeys());
	
	while ((key = (TKey*)nextkey())) 
	{
		std::string className= key->GetClassName();
		std::string keyName= key->GetName();
    TClass* cl = gROOT->GetClass(className.c_str());
    if (!cl) continue;
	
		cout<<"INFO: Processing object key "<<keyName<<" (name="<<className<<") ..."<<endl;
	

		bool excludeObj= false;
		for(size_t i=0;i<excludedObjNames.size();i++){
			if(keyName==excludedObjNames[i]){
				excludeObj= true;
				break;
			}
		}
		if(excludeObj) {
			cout<<"INFO: Object "<<keyName<<" exluded from the list of objects that will be saved ..."<<endl;
			continue;
		}

		if (cl->InheritsFrom(TTree::Class())) {
    	TTree* T= (TTree*)gInputFile->Get(key->GetName());

			//Write to output file
			gOutputFile->cd();
      TTree* newT= T->CloneTree(-1,"fast");
      newT->Write();
			
		}//close if TTree object
		else{
			TObject* obj= key->ReadObj();
      gOutputFile->cd();
      obj->Write();
      delete obj;			
		}

	}//end loop objects

	return 0;

}//close CloneObjectsInFile()



