#include <AstroObject.h>
#include <AstroObjectParser.h>

using namespace std;
using namespace Caesar;

//Variables
bool gReadRegions= false;
std::string gRegionFileName= "";
std::vector<DS9Region*> gRegions;
std::string gCatalogFileName= "";
std::vector<AstroObject*> gAstroObjects;
std::vector<AstroObject*> gAstroObjects_sel;
AstroObject* gAstroObject= 0;
TTree* gOutputTree= 0;
TFile* gOutputFile= 0;
std::string gOutputFileName= "output.root";
std::string gRegionOutFileName= "ds9.reg";

//Functions
int Init();
int ReadRegions();
int ReadObjects();
int SelectObjects();
int WriteDS9Regions();
int Save();

int ReadSimbadObjects(std::string simbadFileName,std::string regionFileName="",std::string outputFileName="output.root",std::string regionOutFileName="ds9.reg")
{
	//===================================
	//==       GET ARGS
	//===================================
	gCatalogFileName= simbadFileName;
	gOutputFileName= outputFileName;		
	gRegionFileName= regionFileName;
	gRegionOutFileName= regionOutFileName;
	gReadRegions= (gRegionFileName!="");

	//===================================
	//==       INIT
	//===================================
	cout<<"INFO: Initializing macro data..."<<endl;
	Init();

	//===================================
	//==       READ REGIONS
	//===================================
	if(gReadRegions)
	{
		cout<<"INFO: Reading region file "<<gRegionFileName<<" ..."<<endl;
		if(ReadRegions()<0){
			cerr<<"ERROR: Failed to read region file "<<gRegionFileName<<"!"<<endl;	
			return -1;
		}
	}

	
	//===================================
	//==       READ CATALOG OBJECTS
	//===================================
	if(ReadObjects()<0){
		cerr<<"ERROR: Failed to read objects from catalog!"<<endl;
		return -1;
	}

	//===================================
	//==       SELECT CATALOG OBJECTS
	//===================================
	if(SelectObjects()<0){
		cerr<<"ERROR: Failed to select objects in regions!"<<endl;
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
	//Initialize file
	gOutputFile= new TFile(gOutputFileName.c_str(),"RECREATE");
	gOutputFile->cd();

	//Initialize TTree
	gOutputTree= new TTree("data","data");
	gOutputTree->Branch("AstroObject",&gAstroObject);

	return 0;

}//close Init()


int ReadRegions()
{
	//Read and parse DS9 regions
	if (DS9RegionParser::Parse(gRegions,gRegionFileName)<0){
		cerr<<"ERROR: Failed to read and parse DS9 region file "<<gRegionFileName<<"!"<<endl;
		return -1;
	}

	//Check if empty
	if(gRegions.empty()){
		cerr<<"No regions read from file "<<gRegionFileName<<"!"<<endl;
		return -1;
	}
	else{
		cout<<"INFO: #"<<gRegions.size()<<" regions read from file "<<gRegionFileName<<"..."<<endl;
	}

	return 0;

}//close ReadRegions()

int ReadObjects()
{
	cout<<"INFO: Reading objects from catalog file "<<gCatalogFileName<<" ..."<<endl;
	int status= AstroObjectParser::ParseSimbadData(gAstroObjects,gCatalogFileName,'|');
	if(status<0){
		cerr<<"ERROR: Failed to read objects from catalog file "<<gCatalogFileName<<"!"<<endl;
		return -1;
	}
	cout<<"#"<<gAstroObjects.size()<<" objects read from catalog file..."<<endl;

	return 0;

}//close ReadObjects()

int SelectObjects()
{
	//If no regions available set selected to all collection
	if(gRegions.empty()){
		gAstroObjects_sel= gAstroObjects;
		return 0;		
	}

	//Loop over sources
	for(size_t i=0;i<gAstroObjects.size();i++)
	{
		double x= gAstroObjects[i]->x;
		double y= gAstroObjects[i]->y;

		//Loop over regions
		bool isInsideRegion= false;
		for(size_t j=0;j<gRegions.size();j++)
		{
			if(gRegions[j]->IsPointInsideRegion(x,y)){
				isInsideRegion= true;
				break;
			}
		}//end loop regions

		//Add to collection if inside region
		if(isInsideRegion){
			gAstroObjects_sel.push_back(gAstroObjects[i]);
		}

	}//end loop objects

	cout<<"INFO: #"<<gAstroObjects_sel.size()<<" objects selected in regions..."<<endl;

	return 0;

}//close SelectObjects()

int WriteDS9Regions()
{
	//Create region file
	FILE* fout= fopen(gRegionOutFileName.c_str(),"w");
	
	//Write header
	fprintf(fout,"global color=green font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"fk5\n");

	//Loop over sources
	for(size_t i=0;i<gAstroObjects_sel.size();i++)
	{
		std::string regionText= gAstroObjects_sel[i]->GetDS9Region();
		fprintf(fout,"%s\n",regionText.c_str());
	}

	//Close file
	fclose(fout);
	
	return 0;

}//close WriteDS9Regions()

int Save()
{
	if(gOutputFile)
	{
		cout<<"INFO: Saving data to file..."<<endl;
		gOutputFile->cd();

		for(size_t i=0;i<gAstroObjects_sel.size();i++)
		{
			gAstroObject= gAstroObjects_sel[i];
			if(gAstroObject && gOutputTree){
				gOutputTree->Fill();
			}
		}//end loop objects

		if(gOutputTree) gOutputTree->Write();
		gOutputFile->Close();

	}//close if

	//Save DS9 regions
	cout<<"INFO: Saving objects to DS9 region..."<<endl;
	WriteDS9Regions();

	return 0;

}//close Save()



