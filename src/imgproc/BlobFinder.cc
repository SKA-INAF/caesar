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



Image* BlobFinder::GetMultiScaleBlobMask(Image* img,int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel,double thrFactor){

	//## Check imge
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
		Image* filterMap= img->GetNormLoGImage(kernelSize,sigma,true);
		filterMap->ComputeStats(true,false,true);
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
		else if(thrModel==2) thrLevel= medianRMSThr;
		else thrLevel= medianRMSThr;
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

}//close namespace
