#include <SFinderThread.h>
#include <SFinder.h>

//Caesar headers
#include <Img.h>
#include <BkgData.h>
#include <Contour.h>
#include <Logger.h>
#include <WorkerData.h>
#include <Serializer.h>

#include <tango.h>

//ROOT headers
#include <TFile.h>

//## Standard headers
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <ctime>
#include <stdexcept>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>

#include <map>
#include <vector>
#include <thread>
#include <memory>
#include <functional>
#include <chrono>
#include <regex>
#include <exception>


using namespace std;
using namespace Caesar;

namespace SFinder_ns {


SFinder* SFinderThread::m_device;
log4tango::Logger* SFinderThread::m_logger;

SFinderThread::SFinderThread (SFinder* dev) 
	: Tango::LogAdapter(dev) 
{
  m_device = dev;
  m_logger= Tango::LogAdapter(dev).get_logger();
 	m_stopThread = false;

}//close constructor
  
SFinderThread::~SFinderThread(){
	
	DEBUG_LOG("Shutting down existing thread...");
	Stop();

}//close destructor



void SFinderThread::Run(const std::vector<Caesar::WorkerTask*> tasks) {

	//## Set RUNNING state
  INFO_LOG("Starting sfinder thread...");	
	m_device->UpdateState(Tango::RUNNING,"Source finding started");
	

	//## Source finding tasks
	int nTasks= (int)tasks.size();
	INFO_LOG("Starting source finding job (#"<<nTasks<<" tasks) ...");
	
	bool hasFailed= false;

	for(long int i=0;i<nTasks;i++){
		
		INFO_LOG("Starting source finding task no. "<<i+1<<" ...");
		
		if(RunTask(*(tasks[i]))<0){
			ERROR_LOG("Source finding task failed for reco task no. "<<i+1<<"!");
			hasFailed= true;
			break;
		}
	}//end loop reco

	if(hasFailed){
		ERROR_LOG("Source finding job failed!");
	}
	else {
		INFO_LOG("Source finding job terminated with success...");
	}
	
	//Set error reply and free resource			
	m_device->UpdateState(Tango::ON,"Source finding task completed, freeing resource for usage");

}//close Run()


int SFinderThread::RunTask(WorkerTask& task){

	/*
	//## Check task
	if(!task){
		ERROR_LOG("Null ptr to given task data!");
		return -1;
	}
	*/
	std::string filename= task.filename;
	long int tileMinX= task.ix_min;
	long int tileMaxX= task.ix_max;
	long int tileMinY= task.iy_min;
	long int tileMaxY= task.iy_max;
	

	//## Read input image
	INFO_LOG("Reading input image "<<filename<<" X["<<tileMinX<<","<<tileMaxX<<"] Y["<<tileMinY<<","<<tileMaxY<<"]...");
	Img* inputImg= ReadImage(filename,tileMinX,tileMaxX,tileMinY,tileMaxY);
	if(!inputImg){
		ERROR_LOG("Reading of input image failed!");
		return -1;
	}	

	//Check stop thread
	if(m_stopThread){
		INFO_LOG("Stopped computation");
		inputImg->Delete();
		return -1;
	}

	//## Compute input image stats & bkg
	INFO_LOG("Computing input image stats & bkg...");	
	BkgData* bkgData= ComputeStatsAndBkg(inputImg);	
	if(!bkgData){
		ERROR_LOG("ERROR: Failed to compute stats/bkg info!");
		inputImg->Delete();
		return -1;
	}

	//Check stop thread
	if(m_stopThread){
		INFO_LOG("Stopped computation");
		inputImg->Delete();
		delete bkgData;
		bkgData= 0;	
		return -1;
	}

	
	//## Find compact sources
	INFO_LOG("Searching compact sources...");

	WorkerData* compactSourceData= new WorkerData;
	compactSourceData->info= task;
	compactSourceData->data_type= WorkerData::eCOMPACT_SOURCE_DATA;
	if(m_device->attr_searchCompactSources_write && FindCompactSources(*compactSourceData,inputImg,false,bkgData)<0){
		ERROR_LOG("Compact source search failed!");
		return -1;
	}
	INFO_LOG("#"<<(compactSourceData->sources).size()<<" sources found (#"<<(compactSourceData->edge_sources).size()<<"@ edge) ...");

	
	//--> Push compact source data
	INFO_LOG("Pushing compact source data event...");
	if(PushWorkerDataEvent(compactSourceData)<0){
		ERROR_LOG("Failed to push compact source data event!");
		return -1;
	}

	//Check stop thread
	if(m_stopThread){
		INFO_LOG("Stopped computation");
		inputImg->Delete();
		delete bkgData;
		bkgData= 0;	
		return -1;
	}

	
	

	/*
	//## Push compact source pipe
	INFO_LOG("Encoding compact source collection to pipe blob...");
	if(Serializer::SourceCollectionToTangoPipe(m_device->m_compactSourcesPipeBlob,compact_sources)<0){
		ERROR_LOG("Failed to encode compact sources to pipe!");
		inputImg->Delete();
		delete bkgData;
		bkgData= 0;
		return -1;
	}

	INFO_LOG("Pushing compact source pipe event...");
	(m_device->m_compactSourcesPipeBlob).set_name("compactSourcesPipeBlob");
	m_device->push_pipe_event("compactSourcesPipe",&(m_device->m_compactSourcesPipeBlob),true);
	*/


	//...
	//...

	//## Clear data
	inputImg->Delete();
	delete bkgData;
	bkgData= 0;

	return 0;

}//close RunTask()


int SFinderThread::FindCompactSources(WorkerData& compactSourceData,Img* inputImg,bool computeStatsAndBkg,BkgData* inputBkgData){

	//## Check input image
	if(!inputImg) {
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//## Find sources
	INFO_LOG("Finding compact sources ...");	
	std::vector<Source*> sources;	
	int status= FindSources(sources,inputImg,computeStatsAndBkg,inputBkgData);
	if(status<0) {
		ERROR_LOG("Compact source finding failed!");
		return -1;
	}

	//## Retrieve found sources 
	int nSources= (int)sources.size();
	INFO_LOG("#"<<nSources<<" compact sources detected in input image...");
	if(nSources<=0) return 0;

	//## Apply source selection?
	int nSelSources= nSources;
	if(m_device->attr_selectCompactSources_write){
		if(SelectSources(sources)<0){
			ERROR_LOG("Failed to select sources!");
			return -1;
		}
		nSelSources= sources.size();
	}//close if source selection

	//## Add detected sources to the list	
	INFO_LOG("#"<<nSelSources<<" compact sources selected after cuts...");
	for(unsigned int i=0;i<sources.size();i++){
		if(sources[i]->IsAtEdge()){
			(compactSourceData.edge_sources).push_back(sources[i]);
		}
		else{
			(compactSourceData.sources).push_back(sources[i]);
		}
	}//end loop 

	return 0;

}//close FindCompactSources()


int SFinderThread::FindSources(std::vector<Source*>& sources,Img* inputImg,bool computeStatsAndBkg,BkgData* inputBkgData){

	//## Check input image
	if(!inputImg) {
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//## Compute stats and bkg?
	BkgData* bkgData= inputBkgData;
	if(computeStatsAndBkg || !bkgData){
		INFO_LOG("Computing image stats/bkg...");
		bkgData= ComputeStatsAndBkg(inputImg);	
		if(!bkgData){
			ERROR_LOG("Failed to compute stats/bkg info!");
			return -1;
		}
	}

	//Check stop thread
	if(m_stopThread){
		INFO_LOG("Stopped computation");
		if(!inputBkgData){
			delete bkgData;
			bkgData= 0;
		}	
		return -1;
	}

	//## Compute significance map
	Img* significanceMap= inputImg->GetSignificanceMap(bkgData,m_device->attr_useLocalBkg_write);
	if(!significanceMap){
		ERROR_LOG("Failed to compute significance map!");
		if(!inputBkgData){//delete only if inputBkgData is null (e.g. bkg is computed here)
			delete bkgData;
			bkgData= 0;
		}
		return -1;
	}

	//Check stop cnditions
	if(m_stopThread){
		INFO_LOG("Stopped computation");
		if(!inputBkgData){
			delete bkgData;
			bkgData= 0;
		}	
		significanceMap->Delete();
		return -1;
	}

	//## Find sources
	INFO_LOG("Finding sources...");	
	int status= inputImg->FindCompactSource(
		sources,significanceMap,bkgData,
		m_device->attr_seedThr_write,m_device->attr_mergeThr_write,m_device->attr_minNPix_write,m_device->attr_searchNegativeExcess_write,m_device->attr_mergeBelowSeed_write,
		m_device->attr_searchNestedSources_write,m_device->attr_nestedBlobThrFactor_write
	);

	if(status<0) {
		ERROR_LOG("Source finding failed!");
		if(!inputBkgData){
			delete bkgData;
			bkgData= 0;
		}
		significanceMap->Delete();
		return -1;
	}
	int nSources= (int)sources.size();
	INFO_LOG(nSources<<" sources detected in input image...");	
	
	//## Clear allocated data
	if(!inputBkgData){
		delete bkgData;
		bkgData= 0;
	}
	significanceMap->Delete();
		
	return 0;

}//close FindSources()


int SFinderThread::SelectSources(std::vector<Source*>& sources){

	//## Apply source selection?
	int nSources= (int)sources.size();
	if(nSources<=0) return 0;
	
	int nSelSources= 0;
	std::vector<Source*> sources_sel;

	for(int i=0;i<nSources;i++){	
		std::string sourceName= sources[i]->Name;
		int sourceId= sources[i]->Id;
		long int NPix= sources[i]->NPix;
		double X0= sources[i]->X0;
		double Y0= sources[i]->Y0;

		//Is bad source (i.e. line-like blob, etc...)?
		if(!IsGoodSource(sources[i])) {
			DEBUG_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as bad source, skipped!");
			sources[i]->SetGoodSourceFlag(false);
			continue;
		}
			
		//Is point-like source?
		if( IsPointLikeSource(sources[i]) ){
			DEBUG_LOG("Source no. "<<i<<" (name="<<sourceName<<",id="<<sourceId<<", n="<<NPix<<"("<<X0<<","<<Y0<<")) tagged as a point-like source ...");
			sources[i]->SetType(Source::ePointLike);
		}

		//Tag nested sources
		std::vector<Source*> nestedSources= sources[i]->GetNestedSources();
		for(unsigned int j=0;j<nestedSources.size();j++){
			std::string nestedSourceName= nestedSources[j]->Name;
			int nestedSourceId= nestedSources[j]->Id;
			long int nestedNPix= nestedSources[j]->NPix;
			double nestedX0= nestedSources[j]->X0;
			double nestedY0= nestedSources[j]->Y0;

			if(!IsGoodSource(nestedSources[j])) {
				DEBUG_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as bad source, skipped!");
				nestedSources[j]->SetGoodSourceFlag(false);
			}
			if( IsPointLikeSource(nestedSources[j]) ){
				DEBUG_LOG("Source no. "<<i<<": nested source no. "<<j<<" (name="<<nestedSourceName<<",id="<<nestedSourceId<<", n="<<nestedNPix<<"("<<nestedX0<<","<<nestedY0<<")) tagged as a point-like source ...");
				nestedSources[j]->SetType(Source::ePointLike);
			}
		}//end loop nested sources
			
		//Add source to the list	
		sources_sel.push_back(sources[i]);
		nSelSources++;
	}//end loop sources

	INFO_LOG("Added "<<nSelSources<<" compact sources to the selected list...");

	//Clear initial vector (DO NOT CLEAR MEMORY!) and fill with selection (then reset selection)
	sources.clear();
	sources.insert(sources.end(),sources_sel.begin(),sources_sel.end());
	sources_sel.clear();

	return 0;

}//close SelectSources()

bool SFinderThread::IsGoodSource(Source* aSource){
	
	if(!aSource) return false;

	//## Check for pixels 	
	if(aSource->NPix<=0 || (aSource->GetPixels()).size()<=0) return false;

	//## Check for line-like source
	if( (aSource->GetContours()).size()<=0) {
		WARN_LOG("No contour stored for this source, cannot perform check!");
		return true;
	}

	double BoundingBoxMin= ((aSource->GetContours())[0])->BoundingBoxMin;
	if(m_device->attr_useBoundingBoxCut_write && BoundingBoxMin<m_device->attr_minBoundingBoxThr_write) {
		DEBUG_LOG("BoundingBox cut not passed (BoundingBoxMin="<<BoundingBoxMin<<"<"<<m_device->attr_minBoundingBoxThr_write<<")");
		return false;
	}

	//## Add other check here ...
	//...

	return true;

}//close IsGoodSource()

bool SFinderThread::IsPointLikeSource(Source* aSource){

	if(!aSource) return false;
	if(!aSource->HasParameters()) {
		WARN_LOG("No parameters are available for this source (did you compute them?)...point-like check cannot be performed!");
		return true;
	}

	std::string sourceName= aSource->Name;
	int sourceId= aSource->Id;

	//Loop over contours and check if all of them have circular features
	bool isPointLike= true;
	std::vector<Contour*> contours= aSource->GetContours();

	for(unsigned int i=0;i<contours.size();i++){
		Contour* thisContour= contours[i];

		//Test circularity ratio: 1= circle
		if(m_device->attr_useCircRatioCut_write && thisContour->CircularityRatio<m_device->attr_psCircRatioThr_write) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass CircularityRatio cut (CR="<<thisContour->CircularityRatio<<"<"<<m_device->attr_psCircRatioThr_write<<")");
			isPointLike= false;
			break;
		}

		//Test elongation (how symmetrical is the shape): 0=circle,square
		if(m_device->attr_useElongCut_write && thisContour->Elongation>m_device->attr_psElongThr_write) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass Elongation cut (ELONG="<<thisContour->CircularityRatio<<">"<<m_device->attr_psElongThr_write<<")");
			isPointLike= false;
			break;	
		}

		//Test ellipse fit
		if(m_device->attr_useEllipseAreaRatioCut_write && (thisContour->EllipseAreaRatio<m_device->attr_psEllipseAreaRatioMinThr_write || thisContour->EllipseAreaRatio>m_device->attr_psEllipseAreaRatioMaxThr_write)) {
			DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass EllipseAreaRatio cut (EAR="<<thisContour->EllipseAreaRatio<<" outside range ["<<m_device->attr_psEllipseAreaRatioMinThr_write<<","<<m_device->attr_psEllipseAreaRatioMaxThr_write<<"])");
			isPointLike= false;
			break;	
		}

	}//end contour loop
	
	//Check number of pixels
	if(m_device->attr_useMaxNPixCut_write && aSource->NPix>m_device->attr_psMaxNPix_write){
		DEBUG_LOG("Source (name="<<sourceName<<","<<"id="<<sourceId<<") does not pass nMaxPix cut (NPix="<<aSource->NPix<<">"<<m_device->attr_psMaxNPix_write<<")");
		isPointLike= false;
	}

	if(!isPointLike) return false;

	return true;

}//close IsPointLikeSource()

Img* SFinderThread::ReadImage(const std::string& filename,long int tileMinX,long int tileMaxX,long int tileMinY,long int tileMaxY){

	//## Check file
	FileInfo info;
	bool match_extension= false;
	if(!SysUtils::CheckFile(filename,info,match_extension,"")){
		ERROR_LOG("Invalid input file name specified (invalid file path?)!");
		return 0;
	}
	
	//=== ROOT reading ===
	Img* inputImg= 0;
	if(info.extension==".root"){// Read image from ROOT file
		TFile* inputFile = new TFile(filename.c_str(),"READ");
		if(!inputFile || inputFile->IsZombie()){
			ERROR_LOG("Cannot open input file "<<filename<<"!");
			return 0;
		}
		inputImg=  (Img*)inputFile->Get(filename.c_str());
		if(!inputImg){
			ERROR_LOG("Cannot get image from input file "<<filename<<"!");
			return 0;
		}
	}//close if

	//=== FITS reading ===
	else if(info.extension==".fits"){// Read image from FITS file
		inputImg= new Img;
		int status= 0;
		if(tileMinX!=-1 && tileMaxX!=-1 && tileMinY!=-1 && tileMaxY!=-1){
		 	status= inputImg->ReadFITS(filename,tileMinX+1,tileMaxX+1,tileMinY+1,tileMaxY+1);//CFITSIO filter assumes image from 1 to N (not from 0 to N-1)
		}
		else {
			status= inputImg->ReadFITS(filename,tileMinX,tileMaxX,tileMinY,tileMaxY);
		}

		if(status<0){
			ERROR_LOG("Failed to read image from input file "<<filename<<"!");
			if(inputImg) inputImg->Delete();
			return 0;
		}
	}//close else if

	//== Invalid extension ==
	else{
		ERROR_LOG("Invalid file extension detected (ext="<<info.extension<<")!");
		return 0;
	}
	inputImg->SetNameTitle("img","img");
	
	return inputImg;

}//close ReadImage()


BkgData* SFinderThread::ComputeStatsAndBkg(Img* img){

	//## Check input img
	if(!img){
		ERROR_LOG("Null ptr to input image given!");
		return 0;
	}

	//## Compute stats
	INFO_LOG("Computing image stats...");
	bool computeRobustStats= true;
	bool skipNegativePix= m_device->attr_skipNegativePixels_write;//false
	bool forceRecomputing= false;
	if(!img->ComputeStats(computeRobustStats,skipNegativePix,forceRecomputing)<0){
		ERROR_LOG("Stats computing failed!");
		return 0;
	}
	img->LogStats("INFO");

	//## Set local bkg grid/box
	//## If MetaData & beam info are available, interpret grid&box options as multiple of beam
	//## If no info is available (or use of beam info is off) interpret grid&box options as fractions wrt image size
	double boxSizeX= m_device->attr_localBkgBoxSizeX_write;
	double boxSizeY= m_device->attr_localBkgBoxSizeY_write;
	bool useBeamInfoInBkg= m_device->attr_useBeamInfoInBkg_write;
	long int nPixelsInBeam= 0;
	if(useBeamInfoInBkg && img->HasMetaData()){
		nPixelsInBeam= img->GetMetaData()->GetBeamSizeInPixel();	
	}
	
	if(useBeamInfoInBkg && nPixelsInBeam>0){
		INFO_LOG("Setting bkg boxes as ("<<m_device->attr_localBkgBoxSizeX_write<<","<<m_device->attr_localBkgBoxSizeY_write<<") x beam (beam="<<nPixelsInBeam<<" pixels) ...");
		boxSizeX= nPixelsInBeam*m_device->attr_localBkgBoxSizeX_write;
		boxSizeY= nPixelsInBeam*m_device->attr_localBkgBoxSizeY_write;
	}
	else{
		WARN_LOG("Beam information is not available or its usage has been turned off, using image fractions...");
		double Nx= img->GetNbinsX();
		double Ny= img->GetNbinsY();
		boxSizeX= m_device->attr_localBkgBoxSizeX_write*Nx;
		boxSizeY= m_device->attr_localBkgBoxSizeY_write*Ny;
	}

	double gridSizeX= m_device->attr_localBkgGridStepSizeX_write*boxSizeX;
	double gridSizeY= m_device->attr_localBkgGridStepSizeY_write*boxSizeY;
	INFO_LOG("Bkg box ("<<boxSizeX<<","<<boxSizeY<<") pixels, Grid step ("<<gridSizeX<<","<<gridSizeY<<") pixels ...");
		
	//## Compute Bkg
	BkgData* bkgData= img->ComputeBkg(
		m_device->attr_bkgEstimator_write,m_device->attr_useLocalBkg_write,
		boxSizeX,boxSizeY,gridSizeX,gridSizeY,
		m_device->attr_use2ndPassInLocalBkg_write,m_device->attr_skipOutliersInLocalBkg_write,
		m_device->attr_seedThr_write,m_device->attr_mergeThr_write,m_device->attr_minNPix_write
	);

	if(!bkgData) {
		ERROR_LOG("Bkg computing failed!");
		return 0;
	}
		
	return bkgData;

}//close ComputeStatsAndBkg()


int SFinderThread::PushWorkerProgressEvent(){

	//...
	//...

	return 0;

}//close PushWorkerProgressEvent()


int SFinderThread::PushWorkerDataEvent(WorkerData* workerData){

	//Check workerData
	if(!workerData){
		ERROR_LOG("Null ptr to given worker data!");
		return -1;
	}

	/*
	//Encode workerData to buffer
	SBuffer buffer;	
	INFO_LOG("Serializing workerData to buffer ...");
	if(Serializer::WorkerDataToBuffer(buffer,workerData)<0){
		ERROR_LOG("Failed to encode worker data to buffer!");
		return -1;
	}
	INFO_LOG("Buffer size: "<<buffer.size<<", string size="<<(buffer.data).size());
	*/

		

	//Set attribute and push event
	bool release= false;//false by default
	
	(m_device->m_mutex)->lock();

	
	INFO_LOG("Set sourceData attribute ...");
	long int buffer_size= 0;
	try {
		//if(*m_device->attr_encodedSourceData_read) CORBA::string_free(*m_device->attr_encodedSourceData_read);
		//INFO_LOG("Allocating size for buffer: n="<<buffer.size);
		
		//Encode workerData to char array
			
		INFO_LOG("Serializing workerData to buffer ...");
		if(Serializer::WorkerDataToCharArray(m_device->attr_encodedSourceData_read,buffer_size,workerData)<0){
			ERROR_LOG("Failed to encode worker data to char array!");
			throw std::runtime_error("Failed to encode worker data to char array");
		}
		INFO_LOG("Buffer size: "<<buffer_size);
	

		//std::vector<unsigned char> workerDataCharArray(buffer.data.begin(), buffer.data.end());
		//INFO_LOG("workerDataCharArray size="<<workerDataCharArray.size());
			
		//m_device->attr_encodedSourceData_read= m_device->create_DevVarCharArray( (unsigned char*)(buffer.data.c_str()),buffer.size);
		//m_device->attr_encodedSourceData_read= m_device->create_DevVarCharArray(buffer.data.c_str(),buffer.size);

		/*
		if(*m_device->attr_sourceData_read) CORBA::string_free(*m_device->attr_sourceData_read);
		INFO_LOG("Allocating size for buffer: n="<<buffer.size);
		*(m_device->attr_sourceData_read)= CORBA::string_alloc(buffer.size);
		//*(m_device->attr_sourceData_read)= CORBA::string_dup((buffer.data).c_str());
		INFO_LOG("Copying data: n="<<buffer.size);
		memcpy ( *(m_device->attr_sourceData_read), (buffer.data).c_str(), buffer.size+1);	
		//INFO_LOG("Set null termination...");
		//(m_device->attr_sourceData_read)[buffer.size]= 0;
		*/
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Failed to set sourceData attribute!");
		(m_device->m_mutex)->unlock();
		return -1;
	}
	catch(const std::exception& e){
		ERROR_LOG("C++ exception occurred when setting sourceData attribute (err="<<e.what()<<")!");
		(m_device->m_mutex)->unlock();
		return -1;
	}
	catch(...){
		ERROR_LOG("Unknown exception occurred when setting sourceData attribute!");
		(m_device->m_mutex)->unlock();
		return -1;
	}

	INFO_LOG("Push sourceData attribute event...");
	try {
		//m_device->push_change_event("sourceData",m_device->attr_sourceData_read, 1, 0,release);
		m_device->push_change_event("encodedSourceData",m_device->attr_encodedSourceData_read, buffer_size, 0,release);
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Failed to push sourceData event!");
		(m_device->m_mutex)->unlock();
		return -1;
	}	
	(m_device->m_mutex)->unlock();

	return 0;

}//close PushWorkerDataEvent()

}//close namespace
