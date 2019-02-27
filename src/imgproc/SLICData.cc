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
* @file SLICData.cc
* @class SLICData
* @brief SLIC data class
*
* Superpixel data
* @author S. Riggi
* @date 20/01/2015
*/

#include <SLICData.h>
#include <Region.h>
#include <SLIC.h>
#include <Image.h>
#include <CodeUtils.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

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

ClassImp(Caesar::SLICData)
ClassImp(Caesar::SLICContourData)
ClassImp(Caesar::SLICNeighborData)
ClassImp(Caesar::SLICSimilarityData)

namespace Caesar {

SLICData::SLICData() 
{
	//Init data
	inputImg= 0;
	edgeImg= 0;
	laplImg= 0;
	regions.clear();
	pixel_labels.clear();

}//close costructor

SLICData::~SLICData()
{
	//Clear data
	Clear();

}//close destructor


SLICData::SLICData(const SLICData& data)
{
	// Copy constructor
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy constuctor called...");
	#endif
  ((SLICData&)data).Copy(*this);

}//close copy constructor



void SLICData::Copy(TObject& obj) const
{
	//## Copying images
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying images...");
	#endif

	//Delete existing images
	if( ((SLICData&)obj).inputImg ){
		delete ((SLICData&)obj).inputImg;
		((SLICData&)obj).inputImg= 0;
	}
	if( ((SLICData&)obj).edgeImg ){
		delete ((SLICData&)obj).edgeImg;
		((SLICData&)obj).edgeImg= 0;
	}
	if( ((SLICData&)obj).laplImg ){
		delete ((SLICData&)obj).laplImg;
		((SLICData&)obj).laplImg= 0;
	}
	
	//Creating and copying images	
	if(inputImg){
		((SLICData&)obj).inputImg= new Image;
		*((SLICData&)obj).inputImg = *inputImg;
	}
	if(edgeImg){
		((SLICData&)obj).edgeImg= new Image;
		*((SLICData&)obj).edgeImg = *edgeImg;
	}
	if(laplImg){
		((SLICData&)obj).laplImg= new Image;
		*((SLICData&)obj).laplImg = *laplImg;
	}

	//## Copy region collection
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying region collection...");
	#endif

	//Delete first any existing collection
	for(unsigned int i=0;i<(((SLICData&)obj).regions).size();i++){
		if( (((SLICData&)obj).regions)[i] ){
			delete (((SLICData&)obj).regions)[i];
			(((SLICData&)obj).regions)[i]= 0;
		}
	}
	(((SLICData&)obj).regions).clear();
		
	//Allocate and copy
	Region* region= 0;
	for(unsigned int i=0;i<regions.size();i++){
		region= new Region;
		*region= *(regions[i]);
		(((SLICData&)obj).regions).push_back(region);
	}

	//## Copy labels
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copying pixel labels...");
	#endif

	((SLICData&)obj).pixel_labels= pixel_labels;
	
}//close Copy()


SLICData& SLICData::operator=(const SLICData& data) { 
	// Operator =
  if (this != &data) ((SLICData&)data).Copy(*this);
  return *this;
}

void SLICData::Clear()
{
	//Clear images
	ClearImages();
	
	//Clear regions
	ClearRegions();

}//close ClearAll() 


void SLICData::ClearImages(){

	//Delete images
	if(inputImg) inputImg->Delete();
	if(edgeImg) edgeImg->Delete();
	if(laplImg) laplImg->Delete();	

}//close ClearImages()

void SLICData::ClearRegions(){
	
	//Delete regions
	for(unsigned int i=0;i<regions.size();i++){
		if(regions[i]){
			delete regions[i];
			regions[i]= 0;
		}
	}
	regions.clear();
	pixel_labels.clear();

}//close ClearRegions()


int SLICData::Init(Image* img,bool useLogScaleMapping,Image* edgeImage)
{
	//## Clear existing data
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Clearing SLIC data...");
	#endif
	Clear();
	
	//## Compute image stats (needed to normalize it)
	if(!img->HasStats()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Input image has no stats computed, computing them...");
		#endif
		bool computeRobustStats= true;
		bool useRange= false;
		bool forceRecomputing= true;
		//if(img->ComputeStats(true,false,true)<0) {
		if(img->ComputeStats(computeRobustStats,forceRecomputing,useRange)<0) {
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute input image stats!");
			#endif
			return -1;
		}
	}
	
	//## Compute and set input image for SLIC generation (normalized to range [1-256])
	double normmin= 1;
	double normmax= 256;
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Normalize input image in range ["<<normmin<<","<<normmax<<"] ...");
	#endif
	if(useLogScaleMapping) {
		inputImg= img->GetNormalizedImage("LOG",normmin,normmax,true);
	}
	else {
		inputImg= img->GetNormalizedImage("LINEAR",normmin,normmax,true);
	}
	if(!inputImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute input norm image!");
		#endif
		return -1;
	}

	if( !inputImg->HasStats() ){	
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Computing norm image stats...");	
		#endif
		bool computeRobustStats= true;
		bool useRange= false;
		bool forceRecomputing= false;
		//if(inputImg->ComputeStats(true,false,false)<0) {
		if(inputImg->ComputeStats(computeRobustStats,forceRecomputing,useRange)<0) {
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute norm image stats!");
			#endif
			return -1;
		}
	}//close if hasStats
	

	//## Compute and set laplacian image
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing laplacian image...");
	#endif
	laplImg= img->GetLaplacianImage(true);
	
	if(!laplImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute lapl image!");
		#endif
		return -1;
	}
	if(!laplImg->HasStats()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Laplacian image has no stats computed, computing them...");
		#endif
		bool computeRobustStats= true;
		bool useRange= false;
		bool forceRecomputing= false;
		//if(laplImg->ComputeStats(true,false,false)<0) {
		if(laplImg->ComputeStats(computeRobustStats,forceRecomputing,useRange)<0) {
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute lapl image stats!");
			#endif
			return -1;
		}
	}
	
	//## Set edge image (if given in input)
	if(edgeImage){
		if(!edgeImage->HasStats()){
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Edge image has no stats computed, computing them...");
			#endif
			bool computeRobustStats= true;
			bool useRange= false;
			bool forceRecomputing= false;
			//if(edgeImage->ComputeStats(true,false,false)<0) {
			if(edgeImage->ComputeStats(computeRobustStats,forceRecomputing,useRange)<0) {
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to compute edge image stats!");
				#endif
				return -1;
			}
		}
		edgeImg= edgeImage;
	}//close if is edgeImg

	//## Init pixel labels
	long int npixels= img->GetNPixels();
	pixel_labels.assign(npixels,0);

	return 0;

}//close Init()

int SLICData::SetData(Image* img,Image* lapl_img,Image* edge_img)
{
	//Check given pointers			
	if(!img || !lapl_img) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given images, will not initialize data!");
		#endif
		return -1;
	}
	
	//Check image sizes
	if(!img->HasPixelData()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given input image has no data stored, will not initialize data!");
		#endif
		return -1;
	}

	if(!lapl_img->HasSameBinning(img,true)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given laplacian image has a different binning wrt input image, will not initialize data!");
		#endif
		return -1;
	}
				
	if(edge_img && !edge_img->HasSameBinning(img,true)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given edge image has a different binning wrt input image, will not initialize data!");
		#endif
		return -1;
	}

	//Set data
	inputImg= img;	
	laplImg= lapl_img;
	if(edge_img) edgeImg= edge_img;
		
	//Set pixel labels to 0	
	long int npixels= img->GetNPixels();
	pixel_labels.assign(npixels,0);

	return 0;

}//close SetData()


double SLICData::GetSedge(long int ix,long int iy)
{
	if(!edgeImg) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Trying to access to edge map, which is not allocated, returning zero!");
		#endif
		return 0;
	}
			
	return edgeImg->GetPixelValue(ix,iy);
		
}//close GetSedge()

int SLICData::SPGenerator(Image* img,int regionSize,double regParam, int minRegionSize, bool useLogScaleMapping, Image* edgeImage)
{
	//## Check inputs
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return -1;
	}
	long long Nx= img->GetNx();
	long long Ny= img->GetNy();
	if(Nx<=0 || Ny<=0 || regionSize<=0 || regParam<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid image size ("<<Nx<<"x"<<Ny<<") or input options (regionSize="<<regionSize<<", regParam="<<regParam<<")");
		#endif
		return -1;
	}

	//## Initialize SLIC storage data
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Initialize SLIC storage data...");
	#endif
	if(Init(img,useLogScaleMapping,edgeImage)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Initialization failed!");
		#endif
		return -1;
	}

	//## Init SLIC algo data
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Allocating SLIC algo data...");
	#endif
	Image* normImg= inputImg;
	const unsigned long long maxNumIterations = 100;
	const long long numRegionsX = (unsigned long long) ceil( (double)Nx/regionSize) ;
  const long long numRegionsY = (unsigned long long) ceil( (double)Ny/regionSize) ;
  const long long numRegions = numRegionsX * numRegionsY ;
  const long long numPixels = Nx * Ny;
	
	float* edgeMap = 0;
	edgeMap= (float*)calloc(numPixels, sizeof(float));
	float* centers = 0;
	centers= (float*)malloc(sizeof(float)*3*numRegions);
	unsigned int* masses = 0;
	masses= (unsigned int*)malloc(sizeof(unsigned int) * numPixels);
	if(!edgeMap || !centers || !masses){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to allocate memory for the algorithm!");
		#endif
		if(masses) free(masses);
  	if(centers) free(centers);
  	if(edgeMap) free(edgeMap);
		return -1;
	}

	//## Compute edge map (gradient strength)
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Compute edge map (gradient strength)...");	
	#endif

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for (long long j=1;j<Ny-1;++j) {
  	for (long long i=1;i<Nx-1;++i) {
			long int iy= j;
			long int ix= i;
			long int index= i+j*Nx;
			double w_left= normImg->GetPixelValue(ix-1,iy);
			double w_right= normImg->GetPixelValue(ix+1,iy);
			double w_up= normImg->GetPixelValue(ix,iy+1);
			double w_down= normImg->GetPixelValue(ix,iy-1);
			double grad= (w_left-w_right)*(w_left-w_right) + (w_up-w_down)*(w_up-w_down);
			edgeMap[index]= grad;
    }//end loop x
  }//end loop y
  
	//## Initialize K-means centers	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Initializing K-means centers...");
	#endif

 	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(long long index=0; index<(signed)numRegions; index++){
		long long u= index % numRegionsX;
		long long v= ((index-u)/numRegionsX) % numRegionsY;
		
		long long centerx = 0;
    long long centery = 0;
    double minEdgeValue = std::numeric_limits<double>::infinity();

		long long x = (long long)round(regionSize * (u + 0.5));
    long long y = (long long)round(regionSize * (v + 0.5));
    x = max(min(x, Nx-1),0LL);
    y = max(min(y, Ny-1),0LL);

		// Search in a 3x3 neighbourhood the smallest edge response
    for (long long yp = max(0LL, y-1); yp <= min(Ny-1, y+1) ; ++yp) {
      for (long long xp = max(0LL, x-1); xp <= min(Nx-1, x+1) ; ++xp) {
        double thisEdgeValue = edgeMap[xp+yp*Nx];
        if (thisEdgeValue < minEdgeValue) {
        	minEdgeValue = thisEdgeValue ;
          centerx = xp ;
          centery = yp ;
        }
      }//end loop xp
    }//end loop yp
		
		// Initialize the new center at this location
		long long array_index= 3*index;
		centers[array_index] = (float) centerx ;
    centers[array_index+1] = (float) centery ;
    centers[array_index+2] = normImg->GetPixelValue(centerx,centery);
	}
	
	

	//## Run k-means iterations
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Running "<<maxNumIterations<<" iterations...");
	#endif

 	double previousEnergy = std::numeric_limits<double>::infinity();
  double startingEnergy= 0;

  for (unsigned long long iter=0;iter < maxNumIterations ; ++iter) {
    double factor = regParam/(regionSize*regionSize);
    double energy = 0 ;

    // assign pixels to centers
		#ifdef OPENMP_ENABLED
		#pragma omp parallel for collapse(2) reduction(+: energy) 
		#endif
    for (long long y = 0 ; y < (signed)Ny ; ++y) {
      for (long long x = 0 ; x < (signed)Nx ; ++x) {
        long long u = floor((double)x/regionSize - 0.5);
        long long v = floor((double)y/regionSize - 0.5);
				float z = normImg->GetPixelValue(x,y);

        double minDistance = std::numeric_limits<double>::infinity();
        for (long long vp = max(0LL, v); vp <= min(numRegionsY-1, v+1) ; ++vp) {
          for (long long up = max(0LL, u); up <= min(numRegionsX-1, u+1) ; ++up) {
            long long region = up  + vp * numRegionsX;
            double centerx = centers[3*region + 0];
            double centery = centers[3*region + 1];
						double centerz = centers[3*region + 2];
						//float z = normImg->GetPixelValue(x,y);//moved out from nested loop

            double spatial = (x - centerx) * (x - centerx) + (y - centery) * (y - centery);
            double appearance = (z - centerz) * (z - centerz);
            double distance= appearance + factor * spatial;
            if (distance<minDistance) {
              minDistance = distance;
              //segmentation[x + y * Nx] = (unsigned int)region;
							pixel_labels[x + y * Nx] = static_cast<long int>(region);
            }
          }//end loop up
        }//end loop vp
        energy += minDistance;
      }//end loop X
    }//end loop Y

    // check energy termination conditions
    if (iter == 0) {
      startingEnergy = energy ;
    } 
		else {
      if ((previousEnergy - energy) < 1.e-5 * (startingEnergy - energy)) {
        break ;
      }
    }
    previousEnergy = energy;

    // Recompute centers
    memset(masses, 0, sizeof(unsigned int) * Nx * Ny);
    memset(centers, 0, sizeof(float) * 3*numRegions);

    for (long long y = 0 ; y<Ny ; ++y) {
      for (long long x = 0 ; x<Nx ; ++x) {
        long long pixel = x + y * Nx;
        long long region = pixel_labels[pixel] ;
        masses[region]++;
        centers[region * 3 + 0] += x;
        centers[region * 3 + 1] += y;
				centers[region * 3 + 2] += normImg->GetPixelValue(x,y);
      }//end loop x
    }//end loop y

    for (long long region = 0 ; region < (signed)numRegions ; ++region) {
      double mass = max((double)masses[region], 1e-8) ;
			for (long long i = 3*region; i<3*(region + 1); ++i) {
        centers[i]/= mass;
			}
    }//end loop region

  }//end loop iterations


	//Free stuff
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Freeing allocated memory for cluster centers & edge map...");
	#endif

  if(masses) free(masses);
  if(centers) free(centers);
  if(edgeMap) free(edgeMap);

	//## Eliminate small regions
	#ifdef LOGGING_ENABLED
 		DEBUG_LOG("Eliminating small regions...");
	#endif

  unsigned int* cleaned = (unsigned int*)calloc(numPixels, sizeof(unsigned int));
  unsigned long long* segment = (unsigned long long*)malloc(sizeof(unsigned long long) * numPixels);
	if(!cleaned || !segment){	
		if(cleaned) free(cleaned);
		if(segment) free(segment);
		return -1;	
	}

 	const long long dx [] = {+1, -1,  0,  0} ;
  const long long dy [] = { 0,  0, +1, -1} ;
 
  for (long long pixel=0; pixel < (signed)numPixels ; ++pixel) {
		if (cleaned[pixel]) continue;
		
    unsigned int label = pixel_labels[pixel];
    unsigned long long numExpanded = 0 ;
    unsigned long long segmentSize = 0 ;
    segment[segmentSize++] = pixel;

		// Find an adjacent label, for possible use later
    // Find cleanedLabel as the label of an already cleaned region neighbour of this pixel
		unsigned int cleanedLabel = label + 1;
    cleaned[pixel] = label + 1;
    long long x = pixel % Nx;
    long long y = pixel / Nx;
    for (long long direction = 0 ; direction < 4 ; ++direction) {
    	long long xp = x + dx[direction] ;
      long long yp = y + dy[direction] ;
      long long neighbor = xp + yp * Nx ;
      if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]) {
      	cleanedLabel = cleaned[neighbor] ;
      }
    }//end loop direction

    // Expand the segment	
		while (numExpanded < segmentSize) {
			long long open = segment[numExpanded++] ;
      x = open % Nx ;
      y = open / Nx ;
      for (long long direction = 0 ; direction < 4 ; ++direction) {
      	long long xp = x + dx[direction] ;
        long long yp = y + dy[direction] ;
        long long neighbor = xp + yp * Nx ;
        if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]==0 && pixel_labels[neighbor]==label) {
        	cleaned[neighbor] = label + 1 ;
          segment[segmentSize++] = neighbor ;
        }//close if
      }//end loop for
    }//end while loop

    //Change label to cleanedLabel if the segment is too small
		if ((signed)segmentSize < minRegionSize) {
			while (segmentSize > 0) {
      	cleaned[segment[--segmentSize]] = cleanedLabel ;
      }//end while loop
    }//close if 

  }//end pixel loop

	//## Restore base 0 indexing of the regions
  //for (long long pixel=0; pixel<(signed)numPixels; ++pixel) cleaned[pixel]--;
  //memcpy(segmentation, cleaned, numPixels * sizeof(unsigned int)) ;
		
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for (long long pixel=0; pixel<(signed)numPixels; ++pixel) {
		cleaned[pixel]--;
		pixel_labels[pixel]= cleaned[pixel];
	}

	//## Delete data
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Freeing allocated memory for cleaned and segment...");
	#endif
	if(cleaned) free(cleaned);
  if(segment) free(segment);

	
	//## Allocate regions
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Allocating "<<numRegions<<" regions in list...");
	#endif

	Region* aRegion= 0;
	for(long long k=0;k<numRegions;k++) {		
		aRegion= new Region();
		regions.push_back(aRegion);
	}

	//## Fill regions
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Filling regions...");
	#endif

	Pixel* aPixel= 0;

	for (long int ix=0; ix<Nx; ix++) {
  	for (long int iy=0;iy<Ny; iy++) {

			long int id= ix + iy * Nx;
			//long int id= img->GetBin(ix,iy);
			double x= img->GetX(ix);
			double y= img->GetY(iy);
			double S= img->GetPixelValue(ix,iy);
			double S_curv= laplImg->GetPixelValue(ix,iy);//(slicData->laplImg)->GetPixelValue(ix,iy);
			double S_edge= 0;
			if(edgeImg) edgeImg->GetPixelValue(ix,iy);

			long int label= pixel_labels[id];
			//labels[ix][iy]= label;
			
			if(label<0 || label>=(signed)numRegions){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Skip label "<<label<<" for pixel ("<<ix<<","<<iy<<")...");
				#endif
				continue;
			}
				
			//Create and fill a pixel
			aPixel= new Pixel;
			aPixel->id= id;
			aPixel->SetPhysCoords(x,y);
			aPixel->SetCoords(ix,iy);
			aPixel->S= S;
			aPixel->SetCurv(S_curv);
			aPixel->SetEdge(S_edge);
		
			//Set region id and add pixel to region
			regions[label]->Id= label;
			if(S!=0) regions[label]->AddPixel(aPixel);		
			
    }//end loop Ny
  }//end loop Nx

	
	//## Remove regions without pixels (important in case of empty image zones)
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Removing regions without pixels...");
	#endif
	RemoveEmptyRegions();

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<regions.size()<<" regions present after cleanup...");
	#endif

	//## Compute region parameters
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Computing region parameters...");
	#endif

	if(ComputeRegionParameters()<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/more errors occurred while computing region parameters!");
		#endif
		return -1;
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("End superpixel generation...");
	#endif

	return 0;

}//close SPGenerator()


Image* SLICData::GetSegmentedImage(Image* image,int selectedTag,bool normalize,bool binarize)
{
	//## Check if there are regions
	if(regions.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No regions available (did you generate the superpixel partition?), returning nullptr!");
		#endif
		return nullptr;
	}

	//## Compute segmented image
	return SLIC::GetSegmentedImage(image,regions,selectedTag,normalize,binarize);

}//close GetSegmentedImage()


int SLICData::ComputeRegionParameters(){

	//## Compute region parameters
	int status= 0;
	bool computeRobustStats= true;
	bool forceRecomputing= false;

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for(size_t k=0;k<regions.size();k++) {	
		long int nPix= regions[k]->NPix;
		if(nPix>0 && regions[k]->ComputeStats(computeRobustStats,forceRecomputing)<0){	
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute pars for region no. "<<k);
			#endif
			status= -1;
		}
	}//end loop regions

	return status;

}//close ComputeRegionParameters()


void SLICData::RemoveEmptyRegions(){

	//## Remove regions without pixels (important in case of empty image zones)
	std::vector<size_t> regionsToBeDeleted;
	std::map<long int,long int> goodRegionIdMap;

	//Find regions with empty pixels
	for(size_t k=0;k<regions.size();k++) {	
		long int regionId= regions[k]->Id;
		long int nPix= regions[k]->NPix;
		goodRegionIdMap.insert( std::pair<long int,long int>(regionId,1) );
		
		if(nPix>0) continue;

		regionsToBeDeleted.push_back(k);
		goodRegionIdMap[regionId]= -1;
		
	}//end loop regions

	//Check if there are regions to be deleted
	if(regionsToBeDeleted.empty()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("No empty region (without pixels) to be deleted.");
		#endif
		return;
	}

	//Set negative id for regions to be deleted
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for
	#endif
	for (size_t i=0; i<pixel_labels.size(); i++) {
		long int labelId= pixel_labels[i];		
		long int corrFactor= goodRegionIdMap[labelId];		
		pixel_labels[i]*= corrFactor;
		//labels[ix][iy]*= corrFactor; 
	}

	//Delete regions
	CodeUtils::DeletePtrItems(regions, regionsToBeDeleted);

}//close RemoveEmptyRegions()


void SLICData::GetRegionIdMap(std::map<long int,long int>& regionIdMap) const 
{			
	regionIdMap.clear();
	for(size_t k=0;k<regions.size();k++){
		long int regionId= regions[k]->Id;
		regionIdMap.insert( std::pair<long int,long int>(regionId,static_cast<long int>(k)) );
	}
		
}//close GetRegionIdMap()

void SLICData::GetRegionIndexMap(std::map<long int,long int>& regionIndexMap) const
{
	regionIndexMap.clear();
	for(size_t k=0;k<regions.size();k++){
		long int regionId= regions[k]->Id;
		regionIndexMap.insert( std::pair<long int,long int>(static_cast<long int>(k),regionId) );
	}
		
}//close GetRegionMapping()


int SLICData::SetPixelLabel(long int gBin,long int label,bool check)
{
	/*
	//Check inputs
	if( check && (gBin<0 || gBin>=pixel_labels.size()) ){
		#ifdef LOGGING_ENABLED		
			ERROR_LOG("Invalid pixel ("<<gBin<<") required to be set!");
		#endif
		return -1;
	}
	*/
	
	//Set pixel label
	pixel_labels[gBin]= label;

	return 0;
	
}//close SetPixelLabel()

int SLICData::SetPixelLabel(long int ix,long int iy,long int label,bool check)
{
	/*
	//Check inputs
	if( check && 
			(
				!inputImg || pixel_labels.empty() || (signed)(pixel_labels.size())!=inputImg->GetNPixels() ||
				ix<0 || ix>=inputImg->GetNx() || iy<0 || iy>=inputImg->GetNy()
			) 
	){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Image/pixel labels mismatch or invalid pixel ("<<ix<<","<<iy<<") required to be set!");
		#endif
		return -1;
	}
	*/

	//NB: No check is performed on argument and data
	long int Nx= inputImg->GetNx();
	long int gBin= ix + iy*Nx;
	pixel_labels[gBin]= label;
	//labels[ix][iy]= label; 
	
	return 0;

}//close SetPixelLabel()


//====================================
//==  SLICNeighborCollection class
//====================================
int SLICNeighborCollection::SetDtot(int index,double Dtot,double Dtot_n)
{			
	//Check index
	if(index<0 || index>=GetN()) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid neighbor index ("<<index<<") requested, valid values are [0,"<<m_neighbors.size()-1<<"]!");
		#endif
		return -1;
	}

	//Set Dtot
	m_neighbors[index].Dtot= Dtot;
	m_neighbors[index].Dtot_n= Dtot_n;

	return 0;
		
}//close SetDtot()


}//close namespace


