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
/**
* @file BlobFinder.cc
* @class BlobFinder
* @brief Blob finder class
*
* Class to perform blob finding 
* @author S. Riggi
* @date 20/01/2015
*/


#include <BlobFinder.h>
#include <Image.h>
#include <BkgData.h>
#include <Blob.h>
#include <Source.h>
#include <Logger.h>
#include <Graph.h>
#include <GausFilter.h>

#include <TObject.h>
#include <TMatrixD.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;

ClassImp(Caesar::BlobFinder)

namespace Caesar {

BlobFinder::BlobFinder() {

	
}//close costructor

BlobFinder::~BlobFinder(){
	

}//close destructor


template <class T>
int BlobFinder::FindBlobs(Image* inputImg,std::vector<T*>& blobs,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image* curvMap){

	//## Check input img
	if(!inputImg){
		ERROR_LOG("Null ptr to given input image!");
		return -1;
	}

	//Get metadata if available
	ImgMetaData* metadata= inputImg->GetMetaData();

	//Set img range data
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	long int Ntot= Nx*Ny;
	float Xmin= inputImg->GetXmin();
	float Ymin= inputImg->GetYmin();
	float Xmax= inputImg->GetXmax();
	float Ymax= inputImg->GetYmax();
	DEBUG_LOG("Image size ("<<Nx<<","<<Ny<<"), Image range(x["<<Xmin<<","<<Xmax<<") y["<<Ymin<<","<<Ymax<<"])");

	//## Check if the flood map is provided otherwise set to the input map
	//## NB: In source search it should be the significance map
	//## NB2: It could be used to fill blobs in the input map, conditional to another map (i.e. a binary mask)
	if(!floodImg) floodImg= inputImg;

	//## Check if bkg data are provided and if local bkg is available
	//## If local bkg is available use it, otherwise use global bkg
	bool hasBkgData= false;
	bool hasLocalBkg= false;
	if(bkgData){
		hasBkgData= true;
		hasLocalBkg= bkgData->HasLocalBkg();
	}

	//## Check curvature data 
	if(curvMap && !curvMap->HasSameBinning(inputImg)){
		ERROR_LOG("Given curvature map has different bnnning wrt input map!");
		return -1;
	}

	//## Set flood threshold
	//Example: merge=4, seed=5  [4,5] o [4,+inf]
	// merge=-4 seed=-5         [-5,-4] o [-inf,-4]            [-inf,4]
	double floodMinThr= mergeThr;
	double floodMaxThr= std::numeric_limits<double>::infinity();
	double floodMinThr_inv= -std::numeric_limits<double>::infinity();
	double floodMaxThr_inv= -mergeThr;
	if(mergeBelowSeed) {
		floodMaxThr= seedThr;
		floodMinThr_inv= seedThr;
	}

	DEBUG_LOG("Flood thr("<<floodMinThr<<","<<floodMaxThr<<") Flood inv thr("<<floodMinThr_inv<<","<<floodMaxThr_inv<<")");
	
	//## Find seed pixels (above seed threshold)	
	std::vector<long int> pixelSeeds;	
	std::vector<bool> isNegativeExcessSeed;
	
	
	for(long int i=0;i<Nx;i++){
		for(long int j=0;j<Ny;j++){	
			long int gBin= inputImg->GetBin(i,j);	
			double Z= floodImg->GetPixelValue(i,j);
			bool isNegative= false;
			if(fabs(Z)>=seedThr) {
				if(Z<0) isNegative= true;
				pixelSeeds.push_back(gBin);	
				isNegativeExcessSeed.push_back(isNegative);
			}
		}//end loop y
	}//end loop x
	
	DEBUG_LOG("#"<<pixelSeeds.size()<<" seeds found ...");

	//## Perform cluster finding starting from detected seeds
	int nBlobs= 0;
	T* aBlob= 0;
	Pixel* aPixel= 0;
	std::vector<bool> isAddedInCluster(Ntot,false);
	
	for(size_t k=0;k<pixelSeeds.size();k++){
		long int seedPixelId= pixelSeeds[k];	
		long int binX= inputImg->GetBinX(seedPixelId);
		long int binY= inputImg->GetBinY(seedPixelId);
		
		//Check if this seed bin has been already assigned to a cluster
		if(isAddedInCluster[seedPixelId]) {
			DEBUG_LOG("Skip pixel seed "<<seedPixelId<<" as was already assigned to a previous blob...");		
			continue;
		}
		
		//Skip negative excess seed if not requested
		if(!findNegativeExcess && isNegativeExcessSeed[k]) {
			DEBUG_LOG("Skip negative excess pixel seed "<<seedPixelId<<"...");
			continue;
		}
		
		//Compute flooded pixels
		DEBUG_LOG("Computing flood-fill around seed pixel "<<seedPixelId<<"...");
		std::vector<long int> clusterPixelIds;
		int status= 0;
		if(isNegativeExcessSeed[k]){
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr_inv,floodMaxThr_inv);
		}
		else {
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr);
		}
		if(status<0) {
			WARN_LOG("Flood fill failed, skip seed!");
			continue;
		}

		//Append cluster pixels to a blob object
		size_t nClusterPixels= clusterPixelIds.size();
		if(nClusterPixels==0 || (int)nClusterPixels<minPixels) {//skip small blobs
			DEBUG_LOG("Blob pixels found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<") below npix threshold (thr="<<minPixels<<"), skip blob!");
			continue;
		}
		DEBUG_LOG("Blob found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<")");
		
		nBlobs++;	
		
		DEBUG_LOG("Adding new blob (# "<<nBlobs<<") to list...");
		//TString blobName= Form("%s_blobId%d",inputImg->GetName().c_str(),nBlobs);
		TString blobName= Form("S%d",nBlobs);
		aBlob= new T();
		aBlob->SetId(nBlobs);	
		aBlob->SetName(std::string(blobName));
		
		for(size_t l=0;l<nClusterPixels;l++){
			long int clusterPixelId= clusterPixelIds[l];	
			if(isAddedInCluster[clusterPixelId]) continue;
			isAddedInCluster[clusterPixelId]= true;//do not forget to add to list of taken pixels!

			long int clusterPixelIdX= inputImg->GetBinX(clusterPixelId);
			long int clusterPixelIdY= inputImg->GetBinY(clusterPixelId);
			double S= inputImg->GetPixelValue(clusterPixelId);			
			double Z= floodImg->GetPixelValue(clusterPixelId);

			double x= inputImg->GetX(clusterPixelIdX);
			double y= inputImg->GetY(clusterPixelIdY);
			long int ix= clusterPixelIdX;
			long int iy= clusterPixelIdY;

			DEBUG_LOG("Adding pixel id="<<clusterPixelId<<", (x,y)=("<<x<<","<<y<<"), (ix,iy)=("<<ix<<","<<iy<<")");
			
			aPixel= new Pixel;
			aPixel->S= S;
			if(fabs(Z)>=seedThr) aPixel->type= Pixel::eSeed;
			else aPixel->type= Pixel::eNormal;
			aPixel->id= clusterPixelId;
			aPixel->SetPhysCoords(x,y);
			aPixel->SetCoords(ix,iy);
				
			//Check if this pixel is at edge, if so mark blob as at edge
			if( inputImg->IsEdgeBin(ix,iy) ) {
				aBlob->SetEdgeFlag(true);
			}
			
			//Set bkg data if available
			if(hasBkgData){
				double bkgLevel= bkgData->gBkg;
				double noiseLevel= bkgData->gNoise;
				if(hasLocalBkg){
					bkgLevel= (bkgData->BkgMap)->GetPixelValue(clusterPixelId);
					noiseLevel= (bkgData->NoiseMap)->GetPixelValue(clusterPixelId);
				}
				aPixel->SetBkg(bkgLevel,noiseLevel);
			}//close if

			//Set curvature data if available
			if(curvMap) {
				double Scurv= curvMap->GetPixelValue(clusterPixelId);
				aPixel->SetCurv(Scurv);
			}
	
			aBlob->AddPixel(aPixel);
		}//end loop cluster pixels

		//## Check if blobs has pixels
		if(!aBlob->HasPixels()){//no pixel...delete blob!
			delete aBlob;
			aBlob= 0;
			continue;
		}

		//## Compute stats
		DEBUG_LOG("Computing blob stats...");
		aBlob->ComputeStats();
		
		//## Compute morphology parameters
		DEBUG_LOG("Computing blob morphology params...");
		aBlob->ComputeMorphologyParams();

		//## Adding image metadata to image (needed for WCS)
		aBlob->SetImageMetaData(metadata,Xmin,Ymin);

		//## Add blob to list
		blobs.push_back(aBlob);
		
	}//end loop seeds

	DEBUG_LOG("#"<<blobs.size()<<" blobs found!");

	return 0;

}//close BlobFinder::FindBlobs()
template int BlobFinder::FindBlobs<Blob>(Image* img,std::vector<Blob*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);
template int BlobFinder::FindBlobs<Source>(Image* img,std::vector<Source*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);


int BlobFinder::FloodFill(Image* img,std::vector<long int>& clusterPixelIds,long int seedPixelId,double floodMinThr,double floodMaxThr){
	
	//Init
	clusterPixelIds.clear();

	//Check image and given seed id
	if(!img){
		ERROR_LOG("Null ptr to image given!");
		return -1;
	}
	if(!img->HasBin(seedPixelId)){//check if given seed actually exists
		ERROR_LOG("Given seed id is outside image range!");
		return -1;
	}

	//Check given flood range
	double seedSignal= img->GetPixelValue(seedPixelId);
	if(seedSignal<floodMinThr || seedSignal>floodMaxThr){
		WARN_LOG("Given flood threshold range does not contain seed, no blobs detected!");
		return -1;
	}
	
	//Add seed to queue and loop over queue
	DEBUG_LOG("Starting flood-fill from seed pixel "<<seedPixelId<<" (floodMinThr="<<floodMinThr<<", floodMaxThr="<<floodMaxThr<<")");
	std::queue<long int> pixelQueue;
	pixelQueue.push(seedPixelId);
	
	int Ntot= img->GetNPixels();
	std::vector<bool> isAddedInQueue(Ntot,false);	
	std::vector<bool> isAddedInCluster(Ntot,false);

	while(!pixelQueue.empty()){

		//Take first pixel in queue, process it and then remove from the queue
		long int gBinId= pixelQueue.front();
		if(!img->HasBin(gBinId)) {
			WARN_LOG("Invalid bin ("<<gBinId<<") put to queue (this should not occur, check!)");
			pixelQueue.pop();
			continue;
		}
		long int binIdX= img->GetBinX(gBinId);
		long int binIdY= img->GetBinY(gBinId);
		DEBUG_LOG("Processing top item in queue (id="<<gBinId<<", (ix,iy)=("<<binIdX<<","<<binIdY<<")");
		pixelQueue.pop();
		

		//Loop on row pixels above threshold
		while (img->IsBinContentInRange(binIdX-1,binIdY,floodMinThr,floodMaxThr)){
    	binIdX--;
    }//close while loop
    
		bool spanUp = false;
    bool spanDown = false;

		DEBUG_LOG("Start flood-fill spanning from (ix,iy)=("<<binIdX<<","<<binIdY<<")");
		 
		while (img->IsBinContentInRange(binIdX,binIdY,floodMinThr,floodMaxThr)) {
   		long int gBinId_cluster= img->GetBin(binIdX,binIdY);
			if(img->HasBin(binIdX,binIdY) && !isAddedInCluster[gBinId_cluster]) {
				DEBUG_LOG("Adding pixel to blob (id="<<gBinId_cluster<<", (ix,iy)=("<<binIdX<<","<<binIdY<<")");
				clusterPixelIds.push_back(gBinId_cluster);
				isAddedInCluster[gBinId_cluster]= true;
			}
			
			//Search up pixels
			long int gBinId_up= img->GetBin(binIdX,binIdY+1);

			if (!spanUp && img->IsBinContentInRange(binIdX,binIdY+1,floodMinThr,floodMaxThr)) {
      	if(img->HasBin(binIdX,binIdY+1) && !isAddedInQueue[gBinId_up]  ) {
					pixelQueue.push(gBinId_up);
					isAddedInQueue[gBinId_up]= true;
					spanUp = true;
				} 
			}//close if
			else if (spanUp && !img->IsBinContentInRange(binIdX,binIdY+1,floodMinThr,floodMaxThr)){
				spanUp = false;
      }

			//Search down pixel
			long int gBinId_down= img->GetBin(binIdX,binIdY-1);

			if (!spanDown && img->IsBinContentInRange(binIdX,binIdY-1,floodMinThr,floodMaxThr)) {
     		if(img->HasBin(binIdX,binIdY-1) && !isAddedInQueue[gBinId_down] ) {
					pixelQueue.push(gBinId_down);
					isAddedInQueue[gBinId_down]= true;
					spanDown = true;
				} 
      }//close if 
			else if (spanDown && !img->IsBinContentInRange(binIdX,binIdY-1,floodMinThr,floodMaxThr)) {
				spanDown = false;
      }
      binIdX++;
		}//end while loop
	}//end queue loop
	
	//Append cluster pixels to a source object
	DEBUG_LOG("#"<<clusterPixelIds.size()<<" cluster pixels found around given seed "<<seedPixelId);
	
	return 0;

}//close BlobFinder::FloodFill()


/*
Image* BlobFinder::ComputeMultiScaleBlobMap(Image* img,double sigmaMin,double sigmaMax,double sigmaStep,double thrFactor,int kernelFactor,double multiplicityThrFactor)
{
	//## Check image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return 0;
	}

	//## Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;	
	int multiplicityThr= static_cast<int>(std::round(multiplicityThrFactor*nScales));
	std::vector<Image*> filterMaps;
	std::vector<double> thresholdLevels;
	
	for(int i=0;i<nScales;i++){
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) kernelSize++;
		INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");

		//Compute LoG filter
		bool invert= true;
		Image* filterMap= img->GetNormLoGImage(kernelSize,sigma,invert);

		//Compute stats
		//NB: Skip negative pixels
		bool useRange= true;
		double minRangeThr= 0;
		//bool skipNegativePixels= true;
		bool computeRobustStats= true;
		bool forceRecomputing= true;
		//filterMap->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
		filterMap->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
		filterMaps.push_back(filterMap);

		//Compute threshold levels
		ImgStats* imgStats= filterMap->GetPixelStats();	
		double median= imgStats->median;
		double medianThr= thrFactor*median;
		double thrLevel= medianThr;
		thresholdLevels.push_back(thrLevel);	
	}//end loop reso

	//Find blobs across scales
	Image* blobMask= img->GetCloned("",true,true);
	blobMask->Reset();

	for(long int i=0;i<blobMask->GetNx();i++){
		for(long int j=0;j<blobMask->GetNy();j++){
			double imgBinContent= img->GetBinContent(i,j);	
			if(imgBinContent==0) continue;
			double wmin= TMath::Infinity();
			double wmax= -TMath::Infinity();
			int multiplicity= 0;
			
			for(size_t k=0;k<filterMaps.size();k++){
				double w= filterMaps[k]->GetPixelValue(i,j);
				if(w>thresholdLevels[k]) multiplicity++;
				if(w<wmin) wmin= w;
				if(w>=wmax) wmax= w;
			}//end loop scales

			if(multiplicity>=multiplicityThr){
				blobMask->FillPixel(i,j,wmax);
			}
			else {
				blobMask->FillPixel(i,j,wmin);
			}

		}//end loop y
	}//end loop x

	//Clear 
	for(unsigned int k=0;k<filterMaps.size();k++){
		if(filterMaps[k]) filterMaps[k]->Delete();		
	}
	filterMaps.clear();

	return blobMask;

}//close ComputeMultiScaleBlobMap()
*/


Image* BlobFinder::ComputeBlobMask(Image* img,double Bmaj,double Bmin,double Bpa,double kernNSigmaSize,double peakZThr,double peakZMergeThr,int minBlobSize,double thrFactor,int bkgEstimator,int bkgBox,double bkgGridStepSize)
{
	//## Check image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return nullptr;
	}

	//## Smooth image with elliptical kernel
	INFO_LOG("Computing smoothed map with elliptical gaussian kernel (bmaj/bmin/bpa="<<Bmaj<<","<<Bmin<<","<<Bpa<<" pixels) ...");
	double kernelScaleFactor= 1;
	Image* filtMap= img->GetBeamConvolvedImage(Bmaj,Bmin,Bpa,kernNSigmaSize,kernelScaleFactor);
	if(!filtMap){
		ERROR_LOG("Failed to compute smoothed map!");
		return nullptr;
	}

	//## Compute curvature map to be thresholded
	INFO_LOG("Computing curvature map ...");
	bool invert= true;
	Image* curvMap= filtMap->GetLaplacianImage(invert);
	if(!curvMap){
		ERROR_LOG("Failed to compute curvature filtered map...");
		CodeUtils::DeletePtr<Image>(filtMap);
		return nullptr;
	}

	//Clear filt map
	CodeUtils::DeletePtr<Image>(filtMap);

	//Compute stats
	//NB: Skip negative and zero pixels
	bool useRange= true;
	double minRangeThr= 0;
	bool computeRobustStats= true;
	bool forceRecomputing= true;
	curvMap->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
		
	//Compute threshold levels
	ImgStats* imgStats= curvMap->GetPixelStats();	
	double median= imgStats->median;
	double medianThr= thrFactor*median;
	double thrLevel= medianThr;
	
	//Zero-threshold filtered map
	INFO_LOG("Thresholding curvature map ...");
	curvMap->ApplyThreshold(thrLevel);

	//Compute bkg map
	INFO_LOG("Computing bkg/rms of curvature map...");
	bool useLocalBkg= true;
	bool use2ndPass= true;
	bool skipOutliers= false;
	bool useRangeInBkg= true;
	ImgBkgData* bkgData= curvMap->ComputeBkg(
		bkgEstimator,
		useLocalBkg,bkgBox,bkgBox,bkgGridStepSize,bkgGridStepSize,	
		use2ndPass,
		skipOutliers,peakZThr,peakZMergeThr,minBlobSize,
		useRangeInBkg,thrLevel
	);
	if(!bkgData){
		ERROR_LOG("Failed to compute bkg map of curvature image!");
		CodeUtils::DeletePtr<Image>(curvMap);	
		return nullptr;
	}

	//Compute significance map
	INFO_LOG("Computing curvature significance map ...");
	Image* significanceMap= curvMap->GetSignificanceMap(bkgData,useLocalBkg);
	if(!significanceMap){
		ERROR_LOG("Failed to compute curvature significance map!");
		CodeUtils::DeletePtr<Image>(curvMap);
		return nullptr;
	}

	//Clear bkg & curv data
	CodeUtils::DeletePtr<Image>(curvMap);
	CodeUtils::DeletePtr<ImgBkgData>(bkgData);	

	//Find peaks in curvature map
	INFO_LOG("Finding peaks in curvature significance map ...");
	//std::vector<TVector2> peakPoints;
	std::vector<ImgPeak> peakPoints;
	bool skipBorders= true;
	double peakKernelMultiplicityThr= 1;	
	int peakShiftTolerance= 2;
	std::vector<int> kernels {3};
	if(significanceMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
		ERROR_LOG("Failed to find peaks in curvature map!");
		CodeUtils::DeletePtr<Image>(significanceMap);
		return nullptr;		
	}
	INFO_LOG("#"<<peakPoints.size()<<" peaks found in curvature map ...");
		
	//## Select peaks (skip peaks at boundary or faint peaks)
	Image* blobMask= img->GetCloned("",true,true);
	blobMask->Reset();
	
	int nBlobs= 0;
	double floodMinThr= std::max(0.,peakZMergeThr);
	double floodMaxThr= std::numeric_limits<double>::infinity();

	for(size_t k=0;k<peakPoints.size();k++){
		//double x= peakPoints[k].X();
		//double y= peakPoints[k].Y();
		double x= peakPoints[k].x;
		double y= peakPoints[k].y;
		long int seedPixelId= significanceMap->FindBin(x,y);
		if(seedPixelId<0){
			ERROR_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
			CodeUtils::DeletePtr<Image>(blobMask);	
			CodeUtils::DeletePtr<Image>(significanceMap);
			return nullptr;			
		}
		long int ix= significanceMap->GetBinX(seedPixelId);
		long int iy= significanceMap->GetBinY(seedPixelId);
		
		//Skip peaks below threshold
		double Zpeak= significanceMap->GetBinContent(seedPixelId);
		if(Zpeak<peakZThr) {
			DEBUG_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak="<<Zpeak<<"<"<<peakZThr<<")");
			continue;
		}
			
		//Find blobs given current seed peak
		std::vector<long int> clusterPixelIds;
		if(FloodFill(significanceMap,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr)<0){
			WARN_LOG("Failed to find blobs in curvature map (seed pix="<<seedPixelId<<"), skip to next...");
			continue;
		}
			
		//Check blob size agaist min size required
		long int nPixInBlob= (long int)(clusterPixelIds.size());
		if(nPixInBlob<(long int)(minBlobSize)){
			INFO_LOG("Skip blob (id="<<seedPixelId<<") as below min size threshold (npix="<<nPixInBlob<<"<"<<minBlobSize<<")");
			continue;
		}
		nBlobs++;
		DEBUG_LOG("Blob found @ (x,y)=("<<ix<<","<<iy<<") (N="<<nPixInBlob<<")");
			
		//Mask blob pixels
		for(size_t k=0;k<clusterPixelIds.size();k++){
			long int gbin= clusterPixelIds[k];
			blobMask->SetPixelValue(gbin,1);
		}
	
	}//end loop peaks
		
	INFO_LOG("#"<<nBlobs<<" blobs found  ...");
		
	//Clear data
	CodeUtils::DeletePtr<Image>(significanceMap);

	return blobMask;

}//close ComputeBlobMask()



Image* BlobFinder::ComputeMultiScaleBlobMask(Image* img,double sigmaMin,double sigmaMax,double sigmaStep,double peakZThr,double peakZMergeThr,int minBlobSize,double thrFactor,int kernelFactor,bool useLocalBkg,int bkgEstimator,int bkgBox,double bkgGridStepSize)
{
	//## Check image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return nullptr;
	}

	//## Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;	
	std::vector<Image*> filterMaps;
	std::vector<Image*> filterSignificanceMaps;
	std::vector<double> thresholdLevels;
	int peakShiftTolerance= 2;
	struct PeakInfo {
		long int gbin;
		long int ix;
		long int iy;
		double S;
		int scale;
		double Z;
		PeakInfo(long int _gbin,long int _ix,long int _iy,double _S, int _scale)
			: gbin(_gbin), ix(_ix), iy(_iy), S(_S), scale(_scale)
		{}
	};
	std::vector<PeakInfo> peaks;
	
	for(int i=0;i<nScales;i++){
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) kernelSize++;
		
		//Compute LoG filter
		INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");
		bool invert= true;
		Image* filterMap= img->GetNormLoGImage(kernelSize,sigma,invert);
		filterMaps.push_back(filterMap);
	
		//Compute stats
		//NB: Skip negative pixels
		//bool skipNegativePixels= true;
		bool useRange= true;
		double minRangeThr= 0;
		bool computeRobustStats= true;
		bool forceRecomputing= true;
		//filterMap->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
		filterMap->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
		
		//Compute threshold levels
		ImgStats* imgStats= filterMap->GetPixelStats();	
		double median= imgStats->median;
		double medianThr= thrFactor*median;
		double thrLevel= medianThr;
		thresholdLevels.push_back(thrLevel);

		//Zero-threshold filtered map
		INFO_LOG("Zero-thresholding LoG map @ scale "<<sigma<<" ...");
		filterMap->ApplyThreshold(thrLevel);

		//Compute bkg map
		INFO_LOG("Computing bkg map @ scale "<<sigma<<"...");
		//bool useLocalBkg= true;
		bool use2ndPass= true;
		bool skipOutliers= false;
		bool useRangeInBkg= true;
		ImgBkgData* bkgData= filterMap->ComputeBkg(
			bkgEstimator,
			useLocalBkg,bkgBox,bkgBox,bkgGridStepSize,bkgGridStepSize,	
			use2ndPass,
			skipOutliers,peakZThr,peakZMergeThr,minBlobSize,
			useRangeInBkg,thrLevel
		);
		if(!bkgData){
			ERROR_LOG("Failed to compute bkg map @ scale "<<sigma<<"!");
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
			return nullptr;
		}

		//Compute significance map
		INFO_LOG("Computing significance map @ scale "<<sigma<<"...");
		Image* filterSignificanceMap= filterMap->GetSignificanceMap(bkgData,useLocalBkg);
		if(!filterSignificanceMap){
			ERROR_LOG("Failed to compute significance map @ scale "<<sigma<<"!");
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
			return nullptr;
		}
		filterSignificanceMaps.push_back(filterSignificanceMap);
		delete bkgData;
		bkgData= 0;

		//Find peaks in filter map
		INFO_LOG("Finding peaks in significance map @ scale "<<sigma<<" ...");
		//std::vector<TVector2> peakPoints;
		std::vector<ImgPeak> peakPoints;
		bool skipBorders= true;
		double peakKernelMultiplicityThr= 1;
		std::vector<int> kernels {3};
		if(filterSignificanceMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
			ERROR_LOG("Failed to find peaks @ scale "<<sigma<<"!");
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
			return nullptr;		
		}
		INFO_LOG("#"<<peakPoints.size()<<" peaks found @ scale "<<sigma<<" ...");
		
		//## Select peaks (skip peaks at boundary or faint peaks)
		std::vector<PeakInfo> peaks_scale;
		for(size_t k=0;k<peakPoints.size();k++){
			//double x= peakPoints[k].X();
			//double y= peakPoints[k].Y();
			double x= peakPoints[k].x;
			double y= peakPoints[k].y;
			long int gbin= filterSignificanceMap->FindBin(x,y);
			if(gbin<0){
				ERROR_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
				CodeUtils::DeletePtrCollection<Image>(filterMaps);	
				CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
				return nullptr;			
			}
			long int ix= filterSignificanceMap->GetBinX(gbin);
			long int iy= filterSignificanceMap->GetBinY(gbin);
		
			//Skip peaks below threshold
			double Zpeak= filterSignificanceMap->GetBinContent(gbin);
			if(Zpeak<peakZThr) {
				DEBUG_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak="<<Zpeak<<"<"<<peakZThr<<")");
				continue;
			}
			double Speak= filterMap->GetBinContent(gbin);
			
			//Add peak to selected peak
			PeakInfo peakInfo(gbin,ix,iy,Speak,(int)(i));
			peakInfo.Z= Zpeak;
			peaks_scale.push_back(peakInfo);
		}//end loop peaks
		peaks.insert(peaks.end(),peaks_scale.begin(),peaks_scale.end());		
		INFO_LOG("#"<<peaks_scale.size()<<" peaks selected @ scale "<<sigma<<" ...");
		
	}//end loop scales

	//## Clear significance maps
	//CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);

	//## Select peak scale according to max peak across scales
	std::vector<PeakInfo> peaks_best;
	std::vector<std::vector<long int>> peakIds;
	for(int i=0;i<nScales;i++){
		peakIds.push_back( std::vector<long int>() );
	}

	if(nScales>1){

		// Create graph with "matching/adjacent" peaks
		Graph linkGraph(peaks.size());
		for(size_t i=0;i<peaks.size()-1;i++) {
			for(size_t j=i+1;j<peaks.size();j++) {	
				long int distX= peaks[i].ix - peaks[j].ix;
				long int distY= peaks[i].iy - peaks[j].iy;
				bool areAdjacent= (fabs(distX)<=peakShiftTolerance && fabs(distY)<=peakShiftTolerance);
				if(!areAdjacent) continue;
				linkGraph.AddEdge(i,j);
			}//end loop peaks
		}//end loop peaks
	
		//Find connected peaks
		std::vector<std::vector<int>> connected_indexes;
		linkGraph.GetConnectedComponents(connected_indexes);

		//Find best scale according to max peak across scales
		for(size_t i=0;i<connected_indexes.size();i++){
			
			double Speak_max= -1.e+99;
			//int bestScaleIndex= 0;
			int index_best= connected_indexes[i][0]; 
			
			for(size_t j=1;j<connected_indexes[i].size();j++){
				int index= connected_indexes[i][j];
				double Speak= peaks[index].S;
				if(Speak>Speak_max){
					Speak_max= Speak;
					//bestScaleIndex= index;
					index_best= index;
				}
			}//end loop items in cluster
		
			int scale_best= peaks[index_best].scale;
			long int gbin_best= peaks[index_best].gbin;
			peaks_best.push_back(peaks[index_best]);
			peakIds[scale_best].push_back(gbin_best);
		}//end loop clusters	

	}//close if
	else{
		for(size_t i=0;i<peaks.size();i++) {
			peaks_best.push_back(peaks[i]);
			peakIds[0].push_back(peaks[i].gbin);
		}
	}//close else

	INFO_LOG("#"<<peaks_best.size()<<" best peaks selected across scales ...");

	
	//Find blobs across scales
	Image* blobMask= img->GetCloned("",true,true);
	blobMask->Reset();

	int nBlobs= 0;

	for(size_t i=0;i<filterMaps.size();i++){
		INFO_LOG("Finding blobs across scale no. "<<i+1<<" (#"<<peakIds[i].size()<<" peaks present) ...");		
		//double floodMinThr= thresholdLevels[i];
		double floodMinThr= std::max(0.,peakZMergeThr);
		double floodMaxThr= std::numeric_limits<double>::infinity();
		
		for(size_t j=0;j<peakIds[i].size();j++){

			//Find blobs given current seed peak
			long int seedPixelId= peakIds[i][j];	
			std::vector<long int> clusterPixelIds;
			//if(FloodFill(filterMaps[i],clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr)<0){
			if(FloodFill(filterSignificanceMaps[i],clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr)<0){
				WARN_LOG("Failed to find blobs by flood-fill @ scale "<<i+1<<" (seed pix="<<seedPixelId<<"), skip to next...");
				continue;
			}
			
			//Check blob size against min size required
			int nPixInBlob= (int)(clusterPixelIds.size());
			if(nPixInBlob<minBlobSize){
				INFO_LOG("Skip blob @ scale "<<i+1<<" (id="<<seedPixelId<<") as below min size threshold (npix="<<nPixInBlob<<"<"<<minBlobSize<<")");
				continue;
			}
			long int ix= filterSignificanceMaps[i]->GetBinX(seedPixelId);
			long int iy= filterSignificanceMaps[i]->GetBinY(seedPixelId);	
			nBlobs++;
			DEBUG_LOG("Blob found @ (x,y)=("<<ix<<","<<iy<<") (N="<<nPixInBlob<<")");
			
			//Mask blob pixels
			for(size_t k=0;k<clusterPixelIds.size();k++){
				long int gbin= clusterPixelIds[k];
				blobMask->SetPixelValue(gbin,1);
			}

		}//end loop blobs per scale
	}//end loop scales

	INFO_LOG("#"<<nBlobs<<" blobs present in mask");

	//## Clear memory 
	CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
	CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			
	return blobMask;

}//close ComputeMultiScaleBlobMask()



int BlobFinder::FindBlendedBlobs(std::vector<Source*>& blendedBlobs,std::vector<ImgPeak>& deblendedPeaks,Image* img,double sigmaMin,double sigmaMax,double sigmaStep,int minBlobSize,double thrFactor,int kernelFactor)
{
	//## Init
	blendedBlobs.clear();
	deblendedPeaks.clear();

	//## Check image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}

	//## Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;	
	std::vector<Image*> filterMaps;
	std::vector<double> thresholdLevels;
	int peakShiftTolerance= 2;
	struct PeakInfo {
		long int gbin;
		long int ix;
		long int iy;
		double x;
		double y;
		double Speak_img;
		double S;
		int scale;
		
		PeakInfo(long int _gbin,long int _ix,long int _iy,double _S, int _scale)
			: gbin(_gbin), ix(_ix), iy(_iy), S(_S), scale(_scale)
		{}
	};
	std::vector<PeakInfo> peaks;
	
	for(int i=0;i<nScales;i++){
		//Set kernel size for this scale
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) {	
			kernelSize++;
		}
		
		//Compute LoG filter
		INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");
		bool invert= true;
		Image* filterMap= img->GetNormLoGImage(kernelSize,sigma,invert);
		filterMaps.push_back(filterMap);
	
		//Compute stats
		//NB: Skip negative pixels
		bool useRange= true;
		double minRangeThr= 0;
		bool computeRobustStats= true;
		bool forceRecomputing= true;
		filterMap->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
		
		//Compute threshold levels
		ImgStats* imgStats= filterMap->GetPixelStats();	
		double median= imgStats->median;
		double medianThr= thrFactor*median;
		double thrLevel= medianThr;
		thresholdLevels.push_back(thrLevel);

		//Zero-threshold filtered map
		INFO_LOG("Zero-thresholding LoG map @ scale "<<sigma<<" (thr="<<thrLevel<<") ...");
		filterMap->ApplyThreshold(thrLevel);

		//Find peaks in filter map
		INFO_LOG("Finding peaks in significance map @ scale "<<sigma<<" ...");
		std::vector<ImgPeak> peakPoints;
		bool skipBorders= true;
		double peakKernelMultiplicityThr= 1;
		std::vector<int> kernels {3};
		if(filterMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
			ERROR_LOG("Failed to find peaks @ scale "<<sigma<<"!");
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			return -1;		
		}
		INFO_LOG("#"<<peakPoints.size()<<" peaks found @ scale "<<sigma<<" ...");
	
		//Skip to next scale if no peaks found at this scale
		if(peakPoints.empty()){
			continue;
		}
		
		//## Select peaks (skip peaks at boundary or faint peaks)
		std::vector<PeakInfo> peaks_scale;
		for(size_t k=0;k<peakPoints.size();k++){
			double x= peakPoints[k].x;
			double y= peakPoints[k].y;
			long int gbin= filterMap->FindBin(x,y);
			
			if(gbin<0){
				ERROR_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
				CodeUtils::DeletePtrCollection<Image>(filterMaps);	
				return -1;			
			}
			long int ix= filterMap->GetBinX(gbin);
			long int iy= filterMap->GetBinY(gbin);
			double Speak= filterMap->GetBinContent(gbin);
			double Speak_img= img->GetBinContent(gbin);
			INFO_LOG("Scale no. "<<i+1<<" (scale="<<sigma<<", peak no. "<<k+1<<", S="<<Speak<<", pos("<<x<<","<<y<<"), pixel pos("<<ix<<","<<iy<<")");


			//Add peak to selected peak
			PeakInfo peakInfo(gbin,ix,iy,Speak,(int)(i));
			peakInfo.x= x;
			peakInfo.y= y;
			peakInfo.Speak_img= Speak_img;
			peaks_scale.push_back(peakInfo);
		}//end loop peaks
		peaks.insert(peaks.end(),peaks_scale.begin(),peaks_scale.end());		
		INFO_LOG("#"<<peaks_scale.size()<<" peaks selected @ scale "<<sigma<<" ...");
		
	}//end loop scales

	//## Return if no peaks found
	if(peaks.empty()){	
		WARN_LOG("No peaks found at all searched scales (NB: this is strange, better check)!");
		CodeUtils::DeletePtrCollection<Image>(filterMaps);	
		return 0;
	}

	//## Select peak scale according to max peak across scales
	std::vector<PeakInfo> peaks_best;
	
	if(nScales>1){

		// Create graph with "matching/adjacent" peaks
		Graph linkGraph(peaks.size());
		for(size_t i=0;i<peaks.size()-1;i++) {
			for(size_t j=i+1;j<peaks.size();j++) {	
				long int distX= peaks[i].ix - peaks[j].ix;
				long int distY= peaks[i].iy - peaks[j].iy;
				bool areAdjacent= (fabs(distX)<=peakShiftTolerance && fabs(distY)<=peakShiftTolerance);
				if(!areAdjacent) continue;
				linkGraph.AddEdge(i,j);
			}//end loop peaks
		}//end loop peaks
	
		//Find connected peaks
		std::vector<std::vector<int>> connected_indexes;
		linkGraph.GetConnectedComponents(connected_indexes);

		//Find best scale according to max peak across scales
		for(size_t i=0;i<connected_indexes.size();i++){
			INFO_LOG("Peak no. "<<i+1<<" detected in "<<connected_indexes[i].size()<<" scales...");

			double Speak_max= -1.e+99;
			//int bestScaleIndex= 0;
			int index_best= connected_indexes[i][0];
			
			for(size_t j=1;j<connected_indexes[i].size();j++){
				int index= connected_indexes[i][j];
				double Speak= peaks[index].S;
				long int ix= peaks[index].ix;
				long int iy= peaks[index].iy;
				
				if(Speak>Speak_max){
					Speak_max= Speak;
					//bestScaleIndex= index;
					index_best= index;
				}
				INFO_LOG("Peak no. "<<i+1<<", scale="<<j<<": pos("<<ix<<","<<iy<<")");
			}//end loop items in cluster
		
			int scale_best= peaks[index_best].scale;
			long int ix_best= peaks[index_best].ix;
			long int iy_best= peaks[index_best].iy;
			INFO_LOG("Peak no. "<<i+1<<", best scale="<<scale_best<<": pos("<<ix_best<<","<<iy_best<<")");
			peaks_best.push_back(peaks[index_best]);
		}//end loop clusters	

	}//close if
	else{
		for(size_t i=0;i<peaks.size();i++) {
			peaks_best.push_back(peaks[i]);
		}
	}//close else

	INFO_LOG("#"<<peaks_best.size()<<" best peaks selected across scales ...");

	
	//Find blended blobs corresponding to best peaks
	std::vector<PeakInfo> peaks_final;
	for(size_t k=0;k<peaks_best.size();k++){
		int peakScale= peaks_best[k].scale;
		long int peakIx= peaks_best[k].ix;
		long int peakIy= peaks_best[k].iy;

		//Create a mask with these peak
		Image* markerImg= (Image*)filterMaps[peakScale]->GetCloned("",true,true);
		markerImg->Reset();
		for(int i=0;i<markerImg->GetNx();i++){
			for(int j=0;j<markerImg->GetNy();j++){
				double binContent= img->GetPixelValue(i,j);
				double w= filterMaps[peakScale]->GetPixelValue(i,j);
				if(binContent==0 || w==0) markerImg->SetPixelValue(i,j,1);//sure bkg
				else markerImg->SetPixelValue(i,j,0);//unknown
			}//end loop bins	
		}//end loop bins
			
		//markerImg->SetPixelValue(peakIx,peakIy,2);//set peak pixel to sure signal	
		int peakStepSize= 2;
		for(long int i=peakIx-peakStepSize;i<=peakIx+peakStepSize;i++){
			for(long int j=peakIy-peakStepSize;j<=peakIy+peakStepSize;j++){
				bool hasBin= markerImg->HasBin(i,j);
				if(!hasBin) continue;
				double binContent= img->GetPixelValue(i,j);
				if(binContent!=0) markerImg->SetPixelValue(i,j,2);
			}
		}
	
		//Extract blob mask by watershed transform
		Image* blobMask= MorphFilter::ComputeWatershedFilter(filterMaps[peakScale],markerImg);
		if(!blobMask){
			ERROR_LOG("Failed to extract blob mask for peak component no. "<<k+1<<"!");
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtr<Image>(markerImg);
			CodeUtils::DeletePtrCollection<Source>(blendedBlobs);
			return -1;
		}
		
		//Clear marker image
		CodeUtils::DeletePtr<Image>(markerImg);

		//Create flood image with mask + peak
		double peakMaskBinContent= blobMask->GetPixelValue(peakIx,peakIy);
		if(peakMaskBinContent<=0){
			WARN_LOG("No blobs extracted around peak "<<k+1<<", skip peak ...");
			CodeUtils::DeletePtr<Image>(blobMask);
			continue;
		}
		blobMask->SetPixelValue(peakIx,peakIy,peakMaskBinContent+1);

		//Extract blended blob using mask as flood image
		std::vector<Source*> blobs;
		double seedThr= 2;
		double mergeThr= 1;
		if(FindBlobs(img,blobs,blobMask,nullptr,seedThr,mergeThr,minBlobSize)<0){
			ERROR_LOG("Failed to find blended blobs from mask!");
			CodeUtils::DeletePtr<Image>(blobMask);
			CodeUtils::DeletePtrCollection<Source>(blendedBlobs);
			return -1;
		}

		//Clear blob mask
		CodeUtils::DeletePtr<Image>(blobMask);

		//Check if more than one blob is found
		//NB: Ideally only 1 blob around desired peak should be found
		if(blobs.empty()){
			WARN_LOG("No blended blob found for peak no. "<<k+1<<" (hint: current method was not able to extract blended blob or blob was below npix thr="<<minBlobSize<<"), go to next peak...");
			continue;
		}	
		else if(blobs.size()>1){
			WARN_LOG("More than one blended blob found for peak no. "<<k+1<<", this should not occur, so skip the peak...");
			continue;
		}
		else{//Add blob to blended blob collection
			blendedBlobs.push_back(blobs[0]);
			peaks_final.push_back(peaks_best[k]);
		}
		
	}//end loop best peaks

	INFO_LOG("#"<<blendedBlobs.size()<<" blended blobs found from #"<<peaks_best.size()<<" peaks...");

	for(size_t i=0;i<peaks_final.size();i++){	
		double x= peaks_final[i].x;
		double y= peaks_final[i].y;
		double S= peaks_final[i].Speak_img;
		long int ix= peaks_final[i].ix;
		long int iy= peaks_final[i].iy;
		deblendedPeaks.push_back( ImgPeak(x,y,S,ix,iy) );
	}

	//## Clear memory 
	CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			
	return 0;

}//close FindBlendedBlobs()

/*
Image* BlobFinder::GetMultiScaleBlobMask(Image* img,int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel,double thrFactor){

	//## Check image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return 0;
	}

	//## Init scales
	int nScales= (sigmaMax-sigmaMin)/sigmaStep + 1;
	
	TString imgName= Form("%s_blobMask",img->GetName().c_str());	
	Image* blobMask= img->GetCloned(std::string(imgName),true,true);
	blobMask->Reset();

	int nbins= 100;
	std::vector<Image*> filterMaps;
	std::vector<double> thresholdLevels;
	
	for(int i=0;i<nScales;i++){
		double sigma= sigmaMin + i*sigmaStep;
		int kernelSize= kernelFactor*sigma;	
		if(kernelSize%2==0) kernelSize++;
		INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");

		//Compute LoG filter
		bool invert= true;
		Image* filterMap= img->GetNormLoGImage(kernelSize,sigma,invert);

		//Compute stats
		//NB: Skip negative pixels
		//bool skipNegativePixels= true;
		bool useRange= true;	
		double minRangeThr= 0;
		bool computeRobustStats= true;
		bool forceRecomputing= true;
		//filterMap->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
		filterMap->ComputeStats(computeRobustStats,forceRecomputing,useRange,minRangeThr);
		filterMaps.push_back(filterMap);

		//Compute threshold levels
		ImgStats* imgStats= filterMap->GetPixelStats();	
		double median= imgStats->median;
		double medianRMS= imgStats->medianRMS;
		double medianThr= thrFactor*median;
		double medianRMSThr= thrFactor*medianRMS;
		double otsuThr= filterMap->FindOtsuThreshold(nbins);
		double valleyThr= filterMap->FindValleyThreshold(nbins,true);
		double optimalThr= std::max(std::min(otsuThr,valleyThr),medianThr);
		double thrLevel= medianRMSThr;
		if(thrModel==1) thrLevel= optimalThr;
		//else if(thrModel==2) thrLevel= medianRMSThr;
		//else thrLevel= medianRMSThr;
		else if(thrModel==2) thrLevel= medianThr;
		else thrLevel= medianThr;
		thresholdLevels.push_back(thrLevel);	
	}//end loop reso
	
	//Find blobs across scales
	for(long int i=0;i<blobMask->GetNx();i++){
		
		for(long int j=0;j<blobMask->GetNy();j++){
			double binContent= img->GetPixelValue(i,j);
			if(binContent==0) continue;

			double wsum= 0;
			int counter= 0;
			for(size_t k=0;k<filterMaps.size();k++){
				double w= filterMaps[k]->GetPixelValue(i,j);
				if(w<thresholdLevels[k]) continue;
				wsum+= w;
				counter++;
			}//end loop scales

			blobMask->FillPixel(i,j,counter);			
		}//end loop y
	}//end loop x

	//Clear 
	for(unsigned int k=0;k<filterMaps.size();k++){
		if(filterMaps[k]) filterMaps[k]->Delete();		
	}
	filterMaps.clear();
	
	return blobMask;

}//close GetMultiScaleBlobMask()
*/

}//close namespace
