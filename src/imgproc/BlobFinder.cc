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
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <Graph.h>
#include <GausFilter.h>
#include <Contour.h>

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
#include <chrono>

using namespace std;

ClassImp(Caesar::BlobFinder)

namespace Caesar {

BlobFinder::BlobFinder() {

	
}//close costructor

BlobFinder::~BlobFinder(){
	

}//close destructor

#ifdef OPENMP_ENABLED
template <class T>
int BlobFinder::FindBlobsMT(Image* inputImg,std::vector<T*>& blobs,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image* curvMap)
{
	//Start timer		
	auto start = chrono::steady_clock::now();

	//## Check input img
	if(!inputImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input image!");
		#endif
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

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image size ("<<Nx<<","<<Ny<<"), Image range(x["<<Xmin<<","<<Xmax<<") y["<<Ymin<<","<<Ymax<<"])");
	#endif

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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given curvature map has different bnnning wrt input map!");
		#endif
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

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Flood thr("<<floodMinThr<<","<<floodMaxThr<<") Flood inv thr("<<floodMinThr_inv<<","<<floodMaxThr_inv<<")");
	#endif

	//## Start blob search (a tile per thread)
	std::vector<T*> blobs_t;	
	std::vector< std::vector<T*> > blobs_per_tile;
	std::vector< std::vector<long int> > pixelSeeds_edge;
	std::vector<long int> tileIy_min;	
	std::vector<long int> tileIy_max;	
	

	#pragma omp parallel private(blobs_t) shared(blobs_per_tile) //reduction(merge: blobs_per_tile)
	{
		//Initialize data structures		
		int thread_id= omp_get_thread_num();
		int nthreads= SysUtils::GetOMPThreads();
		long int tileSize= std::round(Ny/nthreads);
	
		#pragma omp single
   	{
			//Compute image tile partition along y direction
			for(int i=0;i<nthreads;i++) {
				blobs_per_tile.push_back( std::vector<T*>() );
				pixelSeeds_edge.push_back( std::vector<long int>() );

				long int iy_min= i*tileSize;
				long int iy_max= iy_min + tileSize;
				if(iy_max>Ny) iy_max= Ny; 
				tileIy_min.push_back(iy_min);
				tileIy_max.push_back(iy_max);
			}//end loop threads

			#ifdef LOGGING_ENABLED
				DEBUG_LOG("thread_id="<<thread_id<<", nthreads="<<nthreads<<", tileSize="<<tileSize);
				for(size_t i=0;i<tileIy_min.size();i++){
					DEBUG_LOG("tile no. "<<i+1<<", iy=["<<tileIy_min[i]<<","<<tileIy_max[i]<<"]");
				}
			#endif
   	}//close single section

	
		//Search blobs in each tile
		std::vector<bool> isAddedInCluster_t(Ntot,false);
		std::vector<long int> pixelSeeds;
		T* aBlob_t= 0;
		Pixel* aPixel_t= 0;
		long int nBlobs_t= 0;

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Searching blobs in thread id "<<thread_id<<" iy=["<<tileIy_min[thread_id]<<","<<tileIy_max[thread_id]<<"]");
		#endif

		//#pragma omp task	
		//{
		for(long int j=tileIy_min[thread_id];j<tileIy_max[thread_id];j++){	
			
			for(long int i=0;i<Nx;i++){
				//Skip if pixel is below seed thr
				long int seedPixelId= inputImg->GetBin(i,j);	
				double Z= floodImg->GetPixelValue(i,j);
				bool isNegativeExcessSeed= false;
				if(fabs(Z)<seedThr) continue;//no seed pixel
				if(Z<0) isNegativeExcessSeed= true;
				
				//Check if this seed bin has been already assigned to a cluster	
				if(isAddedInCluster_t[seedPixelId]) {
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Skip pixel seed "<<seedPixelId<<" as was already assigned to a previous blob...");		
					#endif
					continue;
				}
		
				//Skip negative excess seed if not requested
				if(!findNegativeExcess && isNegativeExcessSeed) {	
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Skip negative excess pixel seed "<<seedPixelId<<"...");
					#endif
					continue;
				}

				//Add seed to collection
				pixelSeeds.push_back(seedPixelId);
	
				//Compute flooded pixels
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Computing flood-fill around seed pixel "<<seedPixelId<<"...");
				#endif
				std::vector<long int> clusterPixelIds;
				int status= 0;
				if(isNegativeExcessSeed){
					status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr_inv,floodMaxThr_inv);
				}
				else {
					status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr);
				}
				if(status<0) {
					#ifdef LOGGING_ENABLED
						WARN_LOG("Flood fill failed, skip seed!");
					#endif
					continue;
				}

				//CReate a blob from flooded pixels
				nBlobs_t++;
				TString blobName= Form("S%ld_thread%d",nBlobs_t,thread_id);
				aBlob_t= new T();
				aBlob_t->SetId(nBlobs_t);	
				aBlob_t->SetName(std::string(blobName));
		
				bool isBlobAtTileEdge= false;
				size_t nClusterPixels= clusterPixelIds.size();

				for(size_t l=0;l<nClusterPixels;l++){
					long int clusterPixelId= clusterPixelIds[l];	
					if(isAddedInCluster_t[clusterPixelId]) continue;
					isAddedInCluster_t[clusterPixelId]= true;//do not forget to add to list of taken pixels!

					long int ix= inputImg->GetBinX(clusterPixelId);
					long int iy= inputImg->GetBinY(clusterPixelId);
					double S= inputImg->GetPixelValue(clusterPixelId);			
					double Z= floodImg->GetPixelValue(clusterPixelId);
					double x= inputImg->GetX(ix);
					double y= inputImg->GetY(iy);
					if(iy==tileIy_min[thread_id] || iy==tileIy_max[thread_id]-1){
						isBlobAtTileEdge= true;
					}
					
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Adding pixel id="<<clusterPixelId<<", (x,y)=("<<x<<","<<y<<"), (ix,iy)=("<<ix<<","<<iy<<")");
					#endif
			
					aPixel_t= new Pixel;
					aPixel_t->S= S;
					if(fabs(Z)>=seedThr) aPixel_t->type= Pixel::eSeed;
					else aPixel_t->type= Pixel::eNormal;
					aPixel_t->id= clusterPixelId;
					aPixel_t->SetPhysCoords(x,y);
					aPixel_t->SetCoords(ix,iy);
				
					//Check if this pixel is at edge, if so mark blob as at edge
					if( inputImg->IsEdgeBin(ix,iy) ) {
						aBlob_t->SetEdgeFlag(true);
					}
			
					//Set bkg data if available
					if(hasBkgData){
						double bkgLevel= bkgData->gBkg;
						double noiseLevel= bkgData->gNoise;
						if(hasLocalBkg){
							bkgLevel= (bkgData->BkgMap)->GetPixelValue(clusterPixelId);
							noiseLevel= (bkgData->NoiseMap)->GetPixelValue(clusterPixelId);
						}
						aPixel_t->SetBkg(bkgLevel,noiseLevel);
					}//close if

					//Set curvature data if available
					if(curvMap) {
						double Scurv= curvMap->GetPixelValue(clusterPixelId);
						aPixel_t->SetCurv(Scurv);
					}
	
					aBlob_t->AddPixel(aPixel_t);
				}//end loop cluster pixels

				//## Check if blobs has pixels or has too few pixels
				long int nBlobPixels= aBlob_t->GetNPixels();
				if(nBlobPixels==0 || nBlobPixels<minPixels){//no or too few pixels...delete blob!
					delete aBlob_t;
					aBlob_t= 0;
					continue;
				}
	
				//## If blob is at edge delete and save seed for re-processing
				if(isBlobAtTileEdge){
					delete aBlob_t;
					aBlob_t= 0;
					pixelSeeds_edge[thread_id].push_back(seedPixelId);
					continue;			
				}

				//## Compute stats
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Computing blob stats...");
				#endif
				aBlob_t->ComputeStats();
		
				//## Compute morphology parameters
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Computing blob morphology params...");
				#endif
				aBlob_t->ComputeMorphologyParams();

				//## Adding image metadata to image (needed for WCS)
				aBlob_t->SetImageMetaData(metadata,Xmin,Ymin);

				//## Add blob to list
				blobs_t.push_back(aBlob_t);
				
			}//end loop x
		}//end loop y

		//blobs_per_tile[thread_id]= blobs_t; 
		#pragma omp critical
		{
    	blobs_per_tile[thread_id].insert(blobs_per_tile[thread_id].end(), blobs_t.begin(), blobs_t.end());
		}

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("thread_id="<<thread_id<<": #"<<pixelSeeds.size()<<" seeds found (#"<<pixelSeeds_edge[thread_id].size()<<" at edge), #"<<blobs_t.size()<<" blobs found ...");
		#endif

		//}//close task
	
	}//close parallel section

	//## Process "edge" blobs in a non-parallel way
	//First reset list
	std::vector<bool> isAddedInCluster(Ntot,false);
	long int nBlobs= 0;
	T* aBlob= 0;
	Pixel* aPixel= 0;

	//Flood fill around edge seed pixels
	for(size_t k=0;k<pixelSeeds_edge.size();k++){
		for(size_t t=0;t<pixelSeeds_edge[k].size();t++){
			long int seedPixelId= pixelSeeds_edge[k][t];
			double Z= floodImg->GetPixelValue(seedPixelId);
			bool isNegativeExcessSeed= false;
			if(Z<0) isNegativeExcessSeed= true;
				
			//Check if this seed bin has been already assigned to a cluster	
			if(isAddedInCluster[seedPixelId]) {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Skip pixel seed "<<seedPixelId<<" as was already assigned to a previous blob...");		
				#endif
				continue;
			}
		
			//Skip negative excess seed if not requested
			if(!findNegativeExcess && isNegativeExcessSeed) {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Skip negative excess pixel seed "<<seedPixelId<<"...");
				#endif
				continue;
			}

			//Compute flooded pixels
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Computing flood-fill around seed pixel "<<seedPixelId<<"...");
			#endif
			std::vector<long int> clusterPixelIds;
			int status= 0;
			if(isNegativeExcessSeed){
				status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr_inv,floodMaxThr_inv);
			}
			else {
				status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr);
			}
			if(status<0) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Flood fill failed, skip seed!");
				#endif
				continue;
			}

			//Create a blob from flooded pixels
			nBlobs++;
			TString blobName= Form("S%ld",nBlobs);
			aBlob= new T();
			aBlob->SetId(nBlobs);	
			aBlob->SetName(std::string(blobName));
		
			size_t nClusterPixels= clusterPixelIds.size();

			for(size_t l=0;l<nClusterPixels;l++){
				long int clusterPixelId= clusterPixelIds[l];	
				if(isAddedInCluster[clusterPixelId]) continue;
				isAddedInCluster[clusterPixelId]= true;//do not forget to add to list of taken pixels!

				long int ix= inputImg->GetBinX(clusterPixelId);
				long int iy= inputImg->GetBinY(clusterPixelId);
				double S= inputImg->GetPixelValue(clusterPixelId);			
				double Z= floodImg->GetPixelValue(clusterPixelId);
				double x= inputImg->GetX(ix);
				double y= inputImg->GetY(iy);	
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Adding pixel id="<<clusterPixelId<<", (x,y)=("<<x<<","<<y<<"), (ix,iy)=("<<ix<<","<<iy<<")");
				#endif

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

			//## Check if blobs has pixels or has too few pixels
			long int nBlobPixels= aBlob->GetNPixels();
			if(nBlobPixels==0 || nBlobPixels<minPixels){//no or too few pixels...delete blob!
				delete aBlob;
				aBlob= 0;
				continue;
			}

			//## Compute stats
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Computing blob stats...");
			#endif
			aBlob->ComputeStats();
		
			//## Compute morphology parameters
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Computing blob morphology params...");
			#endif
			aBlob->ComputeMorphologyParams();

			//## Adding image metadata to image (needed for WCS)
			aBlob->SetImageMetaData(metadata,Xmin,Ymin);
			
			//## Add blob to list
			blobs.push_back(aBlob);

		}//end loop edge seed pixels in tile
	}//end loop tiles

	//Append non-edge blobs (rename them)
	for(size_t k=0;k<blobs_per_tile.size();k++){
		for(size_t t=0;t<blobs_per_tile[k].size();t++){
			nBlobs++;
			TString blobName= Form("S%ld",nBlobs);
			blobs_per_tile[k][t]->SetId(nBlobs);	
			blobs_per_tile[k][t]->SetName(std::string(blobName));
			
			blobs.push_back(blobs_per_tile[k][t]);
		}//end loop blobs per tile
	}//end loop tiles

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<blobs.size()<<" blobs found!");
	#endif

	//Stop timer and print
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("FindBlobsMT completed in "<<dt<<" ms");	
	#endif

	return 0;

}//close FindBlobsMT()

#else
template <class T>
int BlobFinder::FindBlobsMT(Image* inputImg,std::vector<T*>& blobs,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image* curvMap)
{
	return FindBlobs(inputImg,blobs,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,curvMap);
}
#endif

template int BlobFinder::FindBlobsMT<Blob>(Image* img,std::vector<Blob*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);
template int BlobFinder::FindBlobsMT<Source>(Image* img,std::vector<Source*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);



template <class T>
int BlobFinder::FindBlobsST(Image* inputImg,std::vector<T*>& blobs,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image* curvMap)
{
	//Start timer		
	auto start = chrono::steady_clock::now();

	//## Check input img
	if(!inputImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input image!");
		#endif
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
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image size ("<<Nx<<","<<Ny<<"), Image range(x["<<Xmin<<","<<Xmax<<") y["<<Ymin<<","<<Ymax<<"])");
	#endif

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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given curvature map has different bnnning wrt input map!");
		#endif
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

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Flood thr("<<floodMinThr<<","<<floodMaxThr<<") Flood inv thr("<<floodMinThr_inv<<","<<floodMaxThr_inv<<")");
	#endif

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
	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<pixelSeeds.size()<<" seeds found ...");
	#endif

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
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Skip pixel seed "<<seedPixelId<<" as was already assigned to a previous blob...");		
			#endif
			continue;
		}
		
		//Skip negative excess seed if not requested
		if(!findNegativeExcess && isNegativeExcessSeed[k]) {
			#ifdef LOGGING_ENABLED	
				DEBUG_LOG("Skip negative excess pixel seed "<<seedPixelId<<"...");
			#endif
			continue;
		}
		
		//Compute flooded pixels
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing flood-fill around seed pixel "<<seedPixelId<<"...");
		#endif
		std::vector<long int> clusterPixelIds;
		int status= 0;
		if(isNegativeExcessSeed[k]){
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr_inv,floodMaxThr_inv);
		}
		else {
			status= FloodFill(floodImg,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr);
		}
		if(status<0) {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Flood fill failed, skip seed!");
			#endif
			continue;
		}

		//Append cluster pixels to a blob object
		size_t nClusterPixels= clusterPixelIds.size();
		if(nClusterPixels==0 || (int)nClusterPixels<minPixels) {//skip small blobs
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Blob pixels found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<") below npix threshold (thr="<<minPixels<<"), skip blob!");
			#endif
			continue;
		}
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Blob found @ (x,y)=("<<binX<<","<<binY<<") (N="<<nClusterPixels<<")");
		#endif
		
		nBlobs++;	
		
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Adding new blob (# "<<nBlobs<<") to list...");
		#endif
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

			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Adding pixel id="<<clusterPixelId<<", (x,y)=("<<x<<","<<y<<"), (ix,iy)=("<<ix<<","<<iy<<")");
			#endif

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
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing blob stats...");
		#endif
		aBlob->ComputeStats();
		
		//## Compute morphology parameters
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing blob morphology params...");
		#endif
		aBlob->ComputeMorphologyParams();

		//## Adding image metadata to image (needed for WCS)
		aBlob->SetImageMetaData(metadata,Xmin,Ymin);

		//## Add blob to list
		blobs.push_back(aBlob);
		
	}//end loop seeds

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<blobs.size()<<" blobs found!");
	#endif

	//Stop timer and print
	auto end = chrono::steady_clock::now();
	double dt= chrono::duration <double, milli> (end-start).count();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("FindBlobs completed in "<<dt<<" ms");	
	#endif

	return 0;

}//close BlobFinder::FindBlobsST()
template int BlobFinder::FindBlobsST<Blob>(Image* img,std::vector<Blob*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);
template int BlobFinder::FindBlobsST<Source>(Image* img,std::vector<Source*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);


template <class T>
int BlobFinder::FindBlobs(Image* inputImg,std::vector<T*>& blobs,Image* floodImg,ImgBkgData* bkgData,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image* curvMap)
{
	#ifdef OPENMP_ENABLED
		int nthreads_max= SysUtils::GetOMPMaxThreads();
		if(nthreads_max<=1){
			return FindBlobsST(inputImg,blobs,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,curvMap);
		}		
		else{
			return FindBlobsMT(inputImg,blobs,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,curvMap);
		}
	#else
		return FindBlobsST(inputImg,blobs,floodImg,bkgData,seedThr,mergeThr,minPixels,findNegativeExcess,mergeBelowSeed,curvMap);
	#endif

}//close FindBlobs()
template int BlobFinder::FindBlobs<Blob>(Image* img,std::vector<Blob*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);
template int BlobFinder::FindBlobs<Source>(Image* img,std::vector<Source*>& blobs,Image*,ImgBkgData*,double seedThr,double mergeThr,int minPixels,bool findNegativeExcess,bool mergeBelowSeed,Image*);


int BlobFinder::FloodFill(Image* img,std::vector<long int>& clusterPixelIds,long int seedPixelId,double floodMinThr,double floodMaxThr)
{	
	//Init
	clusterPixelIds.clear();

	//Check image and given seed id
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to image given!");
		#endif
		return -1;
	}
	if(!img->HasBin(seedPixelId)){//check if given seed actually exists
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given seed id is outside image range!");
		#endif
		return -1;
	}

	//Check given flood range
	double seedSignal= img->GetPixelValue(seedPixelId);
	if(seedSignal<floodMinThr || seedSignal>floodMaxThr){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Given flood threshold range does not contain seed, no blobs detected!");
		#endif
		return -1;
	}
	
	//Add seed to queue and loop over queue
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Starting flood-fill from seed pixel "<<seedPixelId<<" (floodMinThr="<<floodMinThr<<", floodMaxThr="<<floodMaxThr<<")");	
	#endif
	std::queue<long int> pixelQueue;
	pixelQueue.push(seedPixelId);
	
	int Ntot= img->GetNPixels();
	std::vector<bool> isAddedInQueue(Ntot,false);	
	std::vector<bool> isAddedInCluster(Ntot,false);

	while(!pixelQueue.empty()){

		//Take first pixel in queue, process it and then remove from the queue
		long int gBinId= pixelQueue.front();
		if(!img->HasBin(gBinId)) {
			#ifdef LOGGING_ENABLED
				WARN_LOG("Invalid bin ("<<gBinId<<") put to queue (this should not occur, check!)");
			#endif
			pixelQueue.pop();
			continue;
		}
		long int binIdX= img->GetBinX(gBinId);
		long int binIdY= img->GetBinY(gBinId);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Processing top item in queue (id="<<gBinId<<", (ix,iy)=("<<binIdX<<","<<binIdY<<")");
		#endif
		pixelQueue.pop();
		

		//Loop on row pixels above threshold
		while (img->IsBinContentInRange(binIdX-1,binIdY,floodMinThr,floodMaxThr)){
    	binIdX--;
    }//close while loop
    
		bool spanUp = false;
    bool spanDown = false;

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Start flood-fill spanning from (ix,iy)=("<<binIdX<<","<<binIdY<<")");
		#endif
		 
		while (img->IsBinContentInRange(binIdX,binIdY,floodMinThr,floodMaxThr)) {
   		long int gBinId_cluster= img->GetBin(binIdX,binIdY);
			if(img->HasBin(binIdX,binIdY) && !isAddedInCluster[gBinId_cluster]) {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Adding pixel to blob (id="<<gBinId_cluster<<", (ix,iy)=("<<binIdX<<","<<binIdY<<")");
				#endif
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
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<clusterPixelIds.size()<<" cluster pixels found around given seed "<<seedPixelId);
	#endif

	return 0;

}//close BlobFinder::FloodFill()


/*
Image* BlobFinder::ComputeMultiScaleBlobMap(Image* img,double sigmaMin,double sigmaMax,double sigmaStep,double thrFactor,int kernelFactor,double multiplicityThrFactor)
{
	//## Check image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
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
		#ifdef LOGGING_ENABLED
			INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");
		#endif

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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return nullptr;
	}

	//## Smooth image with elliptical kernel
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing smoothed map with elliptical gaussian kernel (bmaj/bmin/bpa="<<Bmaj<<","<<Bmin<<","<<Bpa<<" pixels) ...");
	#endif
	double kernelScaleFactor= 1;
	Image* filtMap= img->GetBeamConvolvedImage(Bmaj,Bmin,Bpa,kernNSigmaSize,kernelScaleFactor);
	if(!filtMap){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute smoothed map!");
		#endif
		return nullptr;
	}

	//## Compute curvature map to be thresholded
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing curvature map ...");
	#endif
	bool invert= true;
	Image* curvMap= filtMap->GetLaplacianImage(invert);
	if(!curvMap){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute curvature filtered map...");
		#endif
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
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Thresholding curvature map ...");
	#endif
	curvMap->ApplyThreshold(thrLevel);

	//Compute bkg map
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing bkg/rms of curvature map...");
	#endif
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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute bkg map of curvature image!");
		#endif
		CodeUtils::DeletePtr<Image>(curvMap);	
		return nullptr;
	}

	//Compute significance map
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing curvature significance map ...");
	#endif
	Image* significanceMap= curvMap->GetSignificanceMap(bkgData,useLocalBkg);
	if(!significanceMap){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute curvature significance map!");
		#endif
		CodeUtils::DeletePtr<Image>(curvMap);
		return nullptr;
	}

	//Clear bkg & curv data
	CodeUtils::DeletePtr<Image>(curvMap);
	CodeUtils::DeletePtr<ImgBkgData>(bkgData);	

	//Find peaks in curvature map
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Finding peaks in curvature significance map ...");
	#endif
	//std::vector<TVector2> peakPoints;
	std::vector<ImgPeak> peakPoints;
	bool skipBorders= true;
	double peakKernelMultiplicityThr= 1;	
	int peakShiftTolerance= 2;
	std::vector<int> kernels {3};
	if(significanceMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to find peaks in curvature map!");
		#endif
		CodeUtils::DeletePtr<Image>(significanceMap);
		return nullptr;		
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<peakPoints.size()<<" peaks found in curvature map ...");
	#endif

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
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
			#endif
			CodeUtils::DeletePtr<Image>(blobMask);	
			CodeUtils::DeletePtr<Image>(significanceMap);
			return nullptr;			
		}
		long int ix= significanceMap->GetBinX(seedPixelId);
		long int iy= significanceMap->GetBinY(seedPixelId);
		
		//Skip peaks below threshold
		double Zpeak= significanceMap->GetBinContent(seedPixelId);
		if(Zpeak<peakZThr) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak="<<Zpeak<<"<"<<peakZThr<<")");
			#endif
			continue;
		}
			
		//Find blobs given current seed peak
		std::vector<long int> clusterPixelIds;
		if(FloodFill(significanceMap,clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr)<0){
			#ifdef LOGGING_ENABLED	
				WARN_LOG("Failed to find blobs in curvature map (seed pix="<<seedPixelId<<"), skip to next...");
			#endif
			continue;
		}
			
		//Check blob size agaist min size required
		long int nPixInBlob= (long int)(clusterPixelIds.size());
		if(nPixInBlob<(long int)(minBlobSize)){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Skip blob (id="<<seedPixelId<<") as below min size threshold (npix="<<nPixInBlob<<"<"<<minBlobSize<<")");
			#endif
			continue;
		}
		nBlobs++;

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Blob found @ (x,y)=("<<ix<<","<<iy<<") (N="<<nPixInBlob<<")");
		#endif
			
		//Mask blob pixels
		for(size_t k=0;k<clusterPixelIds.size();k++){
			long int gbin= clusterPixelIds[k];
			blobMask->SetPixelValue(gbin,1);
		}
	
	}//end loop peaks
		
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<nBlobs<<" blobs found  ...");
	#endif

	//Clear data
	CodeUtils::DeletePtr<Image>(significanceMap);

	return blobMask;

}//close ComputeBlobMask()



Image* BlobFinder::ComputeMultiScaleBlobMask(Image* img,double sigmaMin,double sigmaMax,double sigmaStep,double peakZThr,double peakZMergeThr,int minBlobSize,double thrFactor,int kernelFactor,bool useLocalBkg,int bkgEstimator,int bkgBox,double bkgGridStepSize)
{
	//## Check image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
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
		#ifdef LOGGING_ENABLED
			INFO_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<") ...");
		#endif
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
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Zero-thresholding LoG map @ scale "<<sigma<<" ...");
		#endif
		filterMap->ApplyThreshold(thrLevel);

		//Compute bkg map
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing bkg map of LoG map @ scale "<<sigma<<"...");
		#endif
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
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute bkg map @ scale "<<sigma<<"!");
			#endif
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
			return nullptr;
		}

		//Compute significance map
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing significance map of LoG map @ scale "<<sigma<<"...");
		#endif
		Image* filterSignificanceMap= filterMap->GetSignificanceMap(bkgData,useLocalBkg);
		if(!filterSignificanceMap){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute significance map @ scale "<<sigma<<"!");
			#endif
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
			return nullptr;
		}
		filterSignificanceMaps.push_back(filterSignificanceMap);
		delete bkgData;
		bkgData= 0;

		//Find peaks in filter map
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Finding peaks in significance map of LoG map @ scale "<<sigma<<" ...");
		#endif
		std::vector<ImgPeak> peakPoints;
		bool skipBorders= true;
		double peakKernelMultiplicityThr= 1;
		std::vector<int> kernels {3};
		if(filterSignificanceMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to find peaks @ scale "<<sigma<<"!");
			#endif
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
			return nullptr;		
		}

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<peakPoints.size()<<" peaks found @ scale "<<sigma<<" ...");
		#endif
		
		//## Select peaks (skip peaks at boundary or faint peaks)
		std::vector<PeakInfo> peaks_scale;
		for(size_t k=0;k<peakPoints.size();k++){
			double x= peakPoints[k].x;
			double y= peakPoints[k].y;
			long int gbin= filterSignificanceMap->FindBin(x,y);
			if(gbin<0){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
				#endif
				CodeUtils::DeletePtrCollection<Image>(filterMaps);	
				CodeUtils::DeletePtrCollection<Image>(filterSignificanceMaps);
				return nullptr;			
			}
			long int ix= filterSignificanceMap->GetBinX(gbin);
			long int iy= filterSignificanceMap->GetBinY(gbin);
		
			//Skip peaks below threshold
			double Zpeak= filterSignificanceMap->GetBinContent(gbin);
			if(Zpeak<peakZThr) {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak="<<Zpeak<<"<"<<peakZThr<<")");
				#endif
				continue;
			}
			double Speak= filterMap->GetBinContent(gbin);
			
			//Add peak to selected peak
			PeakInfo peakInfo(gbin,ix,iy,Speak,(int)(i));
			peakInfo.Z= Zpeak;
			peaks_scale.push_back(peakInfo);
		}//end loop peaks
		peaks.insert(peaks.end(),peaks_scale.begin(),peaks_scale.end());		

		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<peaks_scale.size()<<" peaks selected @ scale "<<sigma<<" ...");
		#endif
		
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

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<peaks_best.size()<<" best peaks selected across scales ...");
	#endif
	
	//Find blobs across scales
	Image* blobMask= img->GetCloned("",true,true);
	blobMask->Reset();

	int nBlobs= 0;

	for(size_t i=0;i<filterMaps.size();i++){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Finding blobs across scale no. "<<i+1<<" (#"<<peakIds[i].size()<<" peaks present) ...");		
		#endif
		double floodMinThr= std::max(0.,peakZMergeThr);
		double floodMaxThr= std::numeric_limits<double>::infinity();
		
		for(size_t j=0;j<peakIds[i].size();j++){

			//Find blobs given current seed peak
			long int seedPixelId= peakIds[i][j];	
			std::vector<long int> clusterPixelIds;
			if(FloodFill(filterSignificanceMaps[i],clusterPixelIds,seedPixelId,floodMinThr,floodMaxThr)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to find blobs by flood-fill @ scale "<<i+1<<" (seed pix="<<seedPixelId<<"), skip to next...");
				#endif
				continue;
			}
			
			//Check blob size against min size required
			int nPixInBlob= (int)(clusterPixelIds.size());
			if(nPixInBlob<minBlobSize){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Skip blob @ scale "<<i+1<<" (id="<<seedPixelId<<") as below min size threshold (npix="<<nPixInBlob<<"<"<<minBlobSize<<")");
				#endif
				continue;
			}
			long int ix= filterSignificanceMaps[i]->GetBinX(seedPixelId);
			long int iy= filterSignificanceMaps[i]->GetBinY(seedPixelId);	
			nBlobs++;		
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Blob found @ (x,y)=("<<ix<<","<<iy<<") (N="<<nPixInBlob<<")");
			#endif
			
			//Mask blob pixels
			for(size_t k=0;k<clusterPixelIds.size();k++){
				long int gbin= clusterPixelIds[k];
				blobMask->SetPixelValue(gbin,1);
			}

		}//end loop blobs per scale
	}//end loop scales

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<nBlobs<<" blobs found in mask");
	#endif

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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
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
		#ifdef LOGGING_ENABLED	
			DEBUG_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");
		#endif
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
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Zero-thresholding LoG map @ scale "<<sigma<<" (thr="<<thrLevel<<") ...");
		#endif
		filterMap->ApplyThreshold(thrLevel);

		//Find peaks in filter map
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Finding peaks in LoG map @ scale "<<sigma<<" ...");
		#endif
		std::vector<ImgPeak> peakPoints;
		bool skipBorders= true;
		double peakKernelMultiplicityThr= 1;
		std::vector<int> kernels {3};
		if(filterMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to find peaks @ scale "<<sigma<<"!");
			#endif
			CodeUtils::DeletePtrCollection<Image>(filterMaps);	
			return -1;		
		}
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<peakPoints.size()<<" peaks found @ scale "<<sigma<<" ...");
		#endif
	
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
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
				#endif
				CodeUtils::DeletePtrCollection<Image>(filterMaps);	
				return -1;			
			}
			long int ix= filterMap->GetBinX(gbin);
			long int iy= filterMap->GetBinY(gbin);
			double Speak= filterMap->GetBinContent(gbin);
			double Speak_img= img->GetBinContent(gbin);
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Scale no. "<<i+1<<" (scale="<<sigma<<", peak no. "<<k+1<<", S="<<Speak<<", pos("<<x<<","<<y<<"), pixel pos("<<ix<<","<<iy<<")");
			#endif

			//Add peak to selected peak
			PeakInfo peakInfo(gbin,ix,iy,Speak,(int)(i));
			peakInfo.x= x;
			peakInfo.y= y;
			peakInfo.Speak_img= Speak_img;
			peaks_scale.push_back(peakInfo);
		}//end loop peaks
		peaks.insert(peaks.end(),peaks_scale.begin(),peaks_scale.end());		
	
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<peaks_scale.size()<<" peaks selected @ scale "<<sigma<<" ...");
		#endif
		
	}//end loop scales

	//## Return if no peaks found
	if(peaks.empty()){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("No peaks found at all searched scales (NB: this is strange, better check)!");
		#endif
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
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Peak no. "<<i+1<<" detected in "<<connected_indexes[i].size()<<" scales...");
			#endif

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
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Peak no. "<<i+1<<", scale="<<j<<": pos("<<ix<<","<<iy<<")");
				#endif
			}//end loop items in cluster
		
			int scale_best= peaks[index_best].scale;
			long int ix_best= peaks[index_best].ix;
			long int iy_best= peaks[index_best].iy;
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Peak no. "<<i+1<<", best scale="<<scale_best<<": pos("<<ix_best<<","<<iy_best<<")");
			#endif
			peaks_best.push_back(peaks[index_best]);
		}//end loop clusters	

	}//close if
	else{
		for(size_t i=0;i<peaks.size();i++) {
			peaks_best.push_back(peaks[i]);
		}
	}//close else

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<peaks_best.size()<<" best peaks selected across scales ...");
	#endif

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
				if(binContent!=0) markerImg->SetPixelValue(i,j,2);//signal
			}
		}
	
		//Extract blob mask by watershed transform
		Image* blobMask= MorphFilter::ComputeWatershedFilter(filterMaps[peakScale],markerImg);
		if(!blobMask){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to extract blob mask for peak component no. "<<k+1<<"!");
			#endif
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
			#ifdef LOGGING_ENABLED
				WARN_LOG("No blobs extracted around peak "<<k+1<<", skip peak ...");
			#endif
			CodeUtils::DeletePtr<Image>(blobMask);
			continue;
		}
		blobMask->SetPixelValue(peakIx,peakIy,peakMaskBinContent+1);

		//Extract blended blob using mask as flood image
		std::vector<Source*> blobs;
		double seedThr= 2;
		double mergeThr= 1;
		if(FindBlobs(img,blobs,blobMask,nullptr,seedThr,mergeThr,minBlobSize)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to find blended blobs from mask!");
			#endif
			CodeUtils::DeletePtr<Image>(blobMask);
			CodeUtils::DeletePtrCollection<Source>(blendedBlobs);
			return -1;
		}

		//Clear blob mask
		CodeUtils::DeletePtr<Image>(blobMask);

		//Check if more than one blob is found
		//NB: Ideally only 1 blob around desired peak should be found
		if(blobs.empty()){
			#ifdef LOGGING_ENABLED
				WARN_LOG("No blended blob found for peak no. "<<k+1<<" (hint: current method was not able to extract blended blob or blob was below npix thr="<<minBlobSize<<"), go to next peak...");
			#endif
			continue;
		}	
		else if(blobs.size()>1){
			#ifdef LOGGING_ENABLED
				WARN_LOG("More than one blended blob found for peak no. "<<k+1<<", this should not occur, so skip the peak...");
			#endif
			continue;
		}
		else{//Add blob to blended blob collection
			blendedBlobs.push_back(blobs[0]);
			peaks_final.push_back(peaks_best[k]);
		}
		
	}//end loop best peaks

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<blendedBlobs.size()<<" blended blobs found from #"<<peaks_best.size()<<" peaks...");
	#endif

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
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
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
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing LoG map @ scale "<<sigma<<" (step="<<sigmaStep<<", kernsize="<<kernelSize<<")");
		#endif

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
