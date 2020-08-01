
#include <Source.h>
#include <SourceExporter.h>
#include <iostream>

using namespace std;
using namespace Caesar;

//Variables
std::string gRegionFileName= "";
std::vector<DS9Region*> gRegions;
Contour* gRegionContour= 0;

//Functions
int ReadRegions(std::string);
int GetRandPosInRegion(double& x0,double& y0,Contour* contour,double xmin,double xmax,double ymin,double ymax,int nMaxGen=1000);


void ShuffleSourceCatalog(std::string filename="",std::string regionFileName="",std::string outfilename="shuffledSourceComp.dat")
{
	//===================================
	//==       INIT
	//===================================
	// - Init random engine
	delete gRandom;
	gRandom= new TRandom3(0);

	// - Init output file with shuffled sources
	FILE* fout= fopen(outfilename.c_str(),"w");

	//===================================
	//==       READ REGION DATA
	//===================================
	// - Read input file
	if(ReadRegions(regionFileName)<0)
	{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read region file!");
		#endif
		return -1;
	}

	//===================================
	//==       READ DATA
	//===================================
	// - Read input file
	TFile* inputFile= new TFile(filename.c_str(),"READ");
	if(!inputFile || inputFile->IsZombie()) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cannot open file "<<filename<<"!");
		#endif
		return -1;
	}
	
	//## Get source tree
	TTree* SourceInfo= (TTree*)inputFile->Get("SourceInfo");	
	if(!SourceInfo){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Cannot get sources from file "<<filename<<"!");
		#endif
		return -1;
	}

	Source* aSource= 0;
	SourceInfo->SetBranchAddress("Source",&aSource);

	cout<<"# "<<SourceInfo->GetEntries()<<" sources to be drawn..."<<endl;
	

	// - Find region boundary vertex
	double regionXmin= 1.e+99;
	double regionXmax= -1.e+99;
	double regionYmin= 1.e+99;
	double regionYmax= -1.e+99;

	std::vector<TVector2> bb= gRegionContour->BoundingBoxVertex;
	for(size_t k=0;k<bb.size();k++){
		double x= bb[k].X();
		double y= bb[k].Y();
		if(x>=regionXmax) regionXmax= x;
		if(x<=regionXmin) regionXmin= x;
		if(y>=regionYmax) regionYmax= y;
		if(y<=regionYmin) regionYmin= y;
	}
	#ifdef LOGGING_ENABLED
		INFO_LOG("Source shuffle in region with boundary x["<<regionXmin<<","<<regionXmax<<"], y["<<regionYmin<<","<<regionYmax<<"] ...");
	#endif
	
	//===================================
	//==       RANDOMIZE SOURCES
	//===================================
	//## Randomize sources
	WCS* wcs= 0;
	int wcsType= eJ2000;
	bool useFWHM= true;
	bool convertToWCS= true;
	int pixOffset= 0;
	bool useWCSSimpleGeometry= true;
	
	for(int i=0;i<SourceInfo->GetEntries();i++){
		SourceInfo->GetEntry(i);

		Source* source= new Source;
		*source= *aSource;

		//Compute wcs for this source collection if not done
		if(!wcs){
			wcs= source->GetWCS(wcsType);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to compute WCS from source no. "<<i+1<<" (name="<<source->GetName()<<")!");
				#endif
				return -1;
			}
		}

		//Compute wcs source pos 
		std::string sourceName= source->GetName();
		double x0= source->X0;
		double y0= source->Y0;
		double x0_wcs= 0;
		double y0_wcs= 0;
		AstroUtils::PixelToWCSCoords(x0_wcs,y0_wcs,wcs,x0,y0);

		//Shuffle source
		double offsetX= 0;
		double offsetY= 0;
		double x0_wcs_rand= 0;
		double y0_wcs_rand= 0;
		if(GetRandPosInRegion(x0_wcs_rand,y0_wcs_rand,gRegionContour,regionXmin,regionXmax,regionYmin,regionYmax)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute rand pos in region for source no. "<<sourceIndex<<"!");
			#endif
			return -1;
		}
		offsetX= x0_wcs_rand - x0_wcs;
		offsetY= y0_wcs_rand - y0_wcs;
		#ifdef LOGGING_ENABLED
			INFO_LOG("Shuffle source no. "<<i+1<<" (name="<<source->GetName()<<"), pos("<<x0_wcs<<","<<y0_wcs<<"), shuffle pos("<<x0_wcs_rand<<","<<y0_wcs_rand<<"), offset("<<offsetX<<","<<offsetY<<"), offset/arcsec("<<offsetX*3600.<<","<<offsetY*3600<<")");
		#endif

		//Get fitted pars
		bool hasFitInfo= source->HasFitInfo();
		if(!hasFitInfo) continue;

		SourceFitPars fitPars= aSource->GetFitPars();
		int nComponents= fitPars.GetNComponents();	
		std::vector<TEllipse*> fitEllipses;
		if(source->GetFitEllipses(fitEllipses,useFWHM,convertToWCS,wcs,wcsType,pixOffset,useWCSSimpleGeometry)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute WCS ellipse for fitted components of source "<<sourceName<<"!");
			#endif
			return -1;
		}

		for(int k=0;k<nComponents;k++){
			//Skip source component if not selected
			bool isSelected= fitPars.IsSelectedComponent(k);
			if(!isSelected) continue;

			std::string sname= sourceName + std::string(Form("_fitcomp%d",k+1));

			//Apply offset to ellipse 
			double Cx= fitEllipses[k]->GetX1() + offsetX;
			double Cy= fitEllipses[k]->GetY1() + offsetY;
			
			//Write to file
			fprintf(fout,"%s %f %f\n",sname.c_str(),Cx,Cy);
			
		}//end loop components
	}//end loop sources

	//===================================
	//==       SAVE
	//===================================
	// - Close file
	fclose(fout);

	return 0;
	
}//close macro




int ReadRegions(std::string filename)
{
	//Read and parse DS9 regions
	std::vector<DS9Region*> regions;
	if (DS9RegionParser::Parse(regions,filename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse DS9 region file "<<filename<<"!");
		#endif
		return -1;
	}

	//Check if empty
	if(regions.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No regions read from file "<<filename<<"!");
		#endif
		return -1;
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<regions.size()<<"  regions read from file "<<filename<<" ...");
		#endif
	}

	//Compute region contours and check that only one region is given. Multiple region not supported.
	if(regions.size()!=1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("More than one region given (hint: multi-region not supported)!");
		#endif
		return -1;
	}

	//Compute contour
	gRegionContour= regions[0]->GetContour(true);
	if(!gRegionContour){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get contour for catalog region "<<i+1<<"!");
		#endif
		return -1;
	}
	

	return 0;

}//close ReadRegions()


int GetRandPosInRegion(double& x0,double& y0,Contour* contour,double xmin,double xmax,double ymin,double ymax,int nMaxGen)
{
	//NB: Position returned in deg
	int nGen= 0;
	x0= 0;
	y0= 0;
	bool isGoodGen= false;

	while(nGen<nMaxGen) 
	{
		double xrand= gRandom->Uniform(xmin,xmax);				
		double yrand= gRandom->Uniform(ymin,ymax);				
		bool isInsideContour= contour->IsPointInsideContour(xrand,yrand);
		if(isInsideContour){
			x0= xrand;	
			y0= yrand;		
			isGoodGen= true;
			break;
		}
		nGen++;

	}//end while

	if(!isGoodGen){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to generate pos in region (nGen="<<nGen<<">="<<nMaxGen<<")!");
		#endif
		return -1;
	}

	return 0;

}//close GetRandPosInRegion()

