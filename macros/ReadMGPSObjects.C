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
std::string gRegionOutFileName= "sources.reg";

//Functions
int Init();
int ReadRegions();
int ReadObjects();
int SelectObjects();
int SaveDS9Regions(std::string filename);
int Save();

int ReadMGPSObjects(std::string simbadFileName,std::string regionFileName="",std::string outputFileName="output.root")
{
	//===================================
	//==       GET ARGS
	//===================================
	gCatalogFileName= simbadFileName;
	gOutputFileName= outputFileName;		
	gRegionFileName= regionFileName;
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
	SaveDS9Regions(gRegionOutFileName);

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
	int status= AstroObjectParser::ParseMGPSData(gAstroObjects,gCatalogFileName,'|');
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

	//Loop over objects
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


int SaveDS9Regions(std::string filename)
{
	// Set DS9 header
	std::string header= AstroUtils::GetDS9WCSTypeHeader(eJ2000);
	
	// Open file
	FILE* fout= fopen(filename.c_str(),"w");
	fprintf(fout,"global color=blue font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n");
	fprintf(fout,"%s\n",header.c_str());

	//Loop over objects
	bool useImageCoords= false;
	
	for(size_t i=0;i<gAstroObjects_sel.size();i++)
	{
		double x0= gAstroObjects_sel[i]->x;
		double y0= gAstroObjects_sel[i]->y;
		std::string sname= gAstroObjects_sel[i]->name;

		//Get ellipse
		TEllipse* ellipse= gAstroObjects_sel[i]->GetFitEllipse();
		if(!ellipse) continue;

		double R1= ellipse->GetR1();
		double R2= ellipse->GetR2();
		double theta= ellipse->GetTheta();
		double pa= theta;

		//Get encoded string region
		std::string regionText(Form("%s",sname.c_str()));
		
		std::stringstream sstream;
		sstream<<"ellipse "<<x0<<" "<<y0<<" "<<R1<<"\" "<<R2<<"\" "<<pa<<" ";
		sstream<<"# ";
		sstream<<"text={"<<regionText<<"} ";
		std::string region= sstream.str();

		//Write region to file
		fprintf(fout,"%s\n",region.c_str());

	}//end loop objects

	fclose(fout);

	return 0;

}//close SaveDS9Regions()

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

	return 0;

}//close Save()



