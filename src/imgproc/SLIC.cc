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
* @file SLIC.cc
* @class SLIC
* @brief SLIC generator class
*
* Superpixel generator
* @author S. Riggi
* @date 20/01/2015
*/

#include <SLIC.h>
#include <SLICData.h>
#include <Region.h>
#include <Image.h>
#include <CodeUtils.h>
#include <Logger.h>

#include <TObject.h>

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

ClassImp(Caesar::SLIC)


namespace Caesar {


SLIC::SLIC() 
{
	
}

SLIC::~SLIC()
{
	
}


SLICData* SLIC::SPGenerator(Image* img, int regionSize,double regParam, int minRegionSize, bool normalizeImage, bool useLogScaleMapping, Image* laplImg, Image* edgeImg)
{
	//==================================
	//==     Check inputs
	//==================================
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return nullptr;
	}

	long long Nx= img->GetNx();
	long long Ny= img->GetNy();
	const long long numPixels = Nx * Ny;
	if(Nx<=0 || Ny<=0 || regionSize<=0 || regParam<0){
		ERROR_LOG("Invalid image size ("<<Nx<<"x"<<Ny<<") or input options (regionSize="<<regionSize<<", regParam="<<regParam<<")");
		return nullptr;
	}

	//Check image binnings
	if(laplImg && !laplImg->HasSameBinning(img,true)){
		ERROR_LOG("Given laplacian image has a different binning wrt input image!");
		return nullptr;
	}	
	if(edgeImg && !edgeImg->HasSameBinning(img,true)){
		ERROR_LOG("Given edge image has a different binning wrt input image!");
		return nullptr;
	}

	//==================================
	//==     Initialize SLIC data
	//==================================
	INFO_LOG("Create and initialize SLIC data...");
	SLICData* slicData= new SLICData;
	(slicData->pixel_labels).assign(numPixels,0);//initialize with zeros to image size
	
	// Normalize input image for SLIC generation (normalized to range [1-256])?
	Image* normImg= 0;
	if(normalizeImage){	
		// Compute image stats (needed to normalize it)
		if(!img->HasStats()){
			INFO_LOG("Input image has no stats computed, computing them...");
			if(img->ComputeStats(true,false,true)<0) {
				ERROR_LOG("Failed to compute input image stats!");
				return nullptr;
			}
		}
		
		double normmin= 1;
		double normmax= 256;
		INFO_LOG("Normalize input image in range ["<<normmin<<","<<normmax<<"]...");

		if(useLogScaleMapping) {
			normImg= img->GetNormalizedImage("LOG",normmin,normmax,true);
		}
		else {
			normImg= img->GetNormalizedImage("LINEAR",normmin,normmax,true);
		}
		if(!normImg){
			ERROR_LOG("Failed to compute input norm image!");
			return nullptr;
		}	
	}//close if normalize image
	else{
		normImg= img->GetCloned("",true,false);
	}

	//Set input image in SLIC data (copy image)
	slicData->inputImg= normImg;

	//Set edge & lapl image in SLIC data (copy images)
	if(edgeImg) slicData->edgeImg= edgeImg->GetCloned("",true,false);
	if(laplImg) slicData->laplImg= laplImg->GetCloned("",true,false);
	

	//=====================================
	//==     Initialize SLIC algorithm
	//=====================================
	//## Init SLIC algo data
	INFO_LOG("Allocating SLIC algo data...");
	const unsigned long long maxNumIterations = 100;
	const long long numRegionsX = (unsigned long long) ceil( (double)Nx/regionSize) ;
  const long long numRegionsY = (unsigned long long) ceil( (double)Ny/regionSize) ;
  const long long numRegions = numRegionsX * numRegionsY ;
  
	float* edgeMap = 0;
	edgeMap= (float*)calloc(numPixels, sizeof(float));
	float* centers = 0;
	centers= (float*)malloc(sizeof(float)*3*numRegions);
	unsigned int* masses = 0;
	masses= (unsigned int*)malloc(sizeof(unsigned int) * numPixels);
	if(!edgeMap || !centers || !masses){
		ERROR_LOG("Failed to allocate memory for the algorithm!");
		if(masses) free(masses);
  	if(centers) free(centers);
  	if(edgeMap) free(edgeMap);
		delete slicData;
		slicData= 0;
		return nullptr;
	}

	//## Compute edge map (gradient strength)
	INFO_LOG("Compute edge map (gradient strength)...");	
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
	INFO_LOG("Initializing K-means centers...");

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
		//DEBUG_LOG("(u,v)("<<u<<","<<v<<") (x,y)=("<<x<<","<<y<<"), array_index="<<array_index);
	}
	
	
	//=====================================
	//==     Run SLIC algorithm
	//=====================================
	//## Run k-means iterations
	INFO_LOG("Running "<<maxNumIterations<<" iterations...");
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
							(slicData->pixel_labels)[x + y * Nx] = static_cast<long int>(region);
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
        long long region = (slicData->pixel_labels)[pixel] ;
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
	DEBUG_LOG("Freeing allocated memory for cluster centers & edge map...");
  if(masses) free(masses);
  if(centers) free(centers);
  if(edgeMap) free(edgeMap);
	
	//=====================================
	//==     Run SLIC post-processing 
	//=====================================
	//## Eliminate small regions
  INFO_LOG("Eliminating small regions...");
  unsigned int* cleaned = (unsigned int*)calloc(numPixels, sizeof(unsigned int));
  unsigned long long* segment = (unsigned long long*)malloc(sizeof(unsigned long long) * numPixels);
	if(!cleaned || !segment){	
		if(cleaned) free(cleaned);
		if(segment) free(segment);
		delete slicData;
		slicData= 0;
		return nullptr;
	}

 	const long long dx [] = {+1, -1,  0,  0} ;
  const long long dy [] = { 0,  0, +1, -1} ;
 
  for (long long pixel=0; pixel < (signed)numPixels ; ++pixel) {
		if (cleaned[pixel]) continue;
		
    unsigned int label = (slicData->pixel_labels)[pixel];
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
        if (0<=xp && xp<(signed)Nx && 0<=yp && yp<(signed)Ny && cleaned[neighbor]==0 && (slicData->pixel_labels)[neighbor]==label) {
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
		(slicData->pixel_labels)[pixel]= cleaned[pixel];
	}

	//## Delete data
	DEBUG_LOG("Freeing allocated memory for cleaned and segment...");
	if(cleaned) free(cleaned);
  if(segment) free(segment);

	//=====================================
	//==     Set SLIC region data
	//=====================================
	//## Allocate regions
	INFO_LOG("Allocating "<<numRegions<<" regions in list...");
	Region* aRegion= 0;
	for(long long k=0;k<numRegions;k++) {		
		aRegion= new Region();
		(slicData->regions).push_back(aRegion);
	}

	//## Fill regions
	INFO_LOG("Filling regions...");
	Pixel* aPixel= 0;

	for (long int ix=0; ix<Nx; ix++) {
  	for (long int iy=0;iy<Ny; iy++) {

			long int id= ix + iy * Nx;
			double x= img->GetX(ix);
			double y= img->GetY(iy);
			double S= img->GetPixelValue(ix,iy);
			double S_curv= 0;
			if(laplImg) laplImg->GetPixelValue(ix,iy);
			double S_edge= 0;
			if(edgeImg) edgeImg->GetPixelValue(ix,iy);

			//CHeck pixel label
			long int label= (slicData->pixel_labels)[id];		
			if(label<0 || label>=(signed)numRegions){
				WARN_LOG("Skip label "<<label<<" for pixel ("<<ix<<","<<iy<<")...");
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
			(slicData->regions)[label]->Id= label;
			if(S!=0) (slicData->regions)[label]->AddPixel(aPixel);		
			
    }//end loop Ny
  }//end loop Nx

	
	//## Remove regions without pixels (important in case of empty image zones)
	INFO_LOG("Removing regions without pixels...");
	slicData->RemoveEmptyRegions();
	INFO_LOG("#"<<slicData->GetNRegions()<<" regions present after cleanup...");

	//## Compute region parameters
	INFO_LOG("Computing region parameters...");
	if(slicData->ComputeRegionParameters()<0){
		WARN_LOG("One/more errors occurred while computing region parameters!");
	}

	INFO_LOG("End superpixel generation...");

	return slicData;

}//close SPGenerator()


Image* SLIC::GetSegmentedImage(Image* image,std::vector<Region*>const& regions,int selectedTag,bool normalize,bool binarize){

	//## Check input
	if(!image) {
		ERROR_LOG("Null ptr to input image given!");
		return 0;
	}

	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		ERROR_LOG("No regions available, nothing to be done!");
		return 0;
	}

	//## Use image normalization?
	double A= 0;
	double B= 1;
	if(image->HasStats() && normalize){
		ImgStats* stats= image->GetPixelStats();
		double NormMin= image->GetMinimum();
		double NormMax= image->GetMaximum();
		A= NormMin - (NormMax-NormMin)*stats->min/(stats->max-stats->min);
		B= (NormMax-NormMin)/(stats->max-stats->min);
	}	
	
	//## Fill image with region means
	TString imgName= Form("%s_segmented",image->GetName().c_str());
	Image* segmentedImg= (Image*)image->GetCloned(std::string(imgName),true,true);
	segmentedImg->Reset();
	
  for(int i=0;i<nRegions;i++){//loop on regions
		int regionTag= regions[i]->Tag;
		long int nPixelsInRegion= regions[i]->NPix;
		DEBUG_LOG("Region no. "<<i<<" N="<<nPixelsInRegion<<" regionTag="<<regionTag);
		
		if(selectedTag!=-1 && regionTag!=selectedTag) continue;
		
		double Mean= regions[i]->Mean;
		double colorValue= A+B*Mean;
		if(binarize) colorValue= 1;

		for(int j=0;j<nPixelsInRegion;j++){//loop on pixels inside region
			long int ix= (regions[i]->GetPixel(j))->ix;
			long int iy= (regions[i]->GetPixel(j))->iy;
			segmentedImg->FillPixel(ix,iy,colorValue);
			//double x= (regions[i]->GetPixel(j))->x;
			//double y= (regions[i]->GetPixel(j))->y;
			//segmentedImg->Fill(x,y,colorValue);
		}//end loop pixels in region	
	}//end loop regions
	
	return segmentedImg;

}//close GetSegmentedImage()


SLICContourData* SLIC::ComputeBoundaryContours(SLICData* slicData) {
  	
	//## Check input
	if(!slicData || !slicData->inputImg) {
		ERROR_LOG("Null ptr to given input image and/or slic data!");
		return nullptr;
	}
	long int nRegions= slicData->GetNRegions();
	if(nRegions<=0) return 0;

	//## Init data
	INFO_LOG("Computing contours from NR="<<nRegions<<" regions...");
	long int Nx= (slicData->inputImg)->GetNx();
	long int Ny= (slicData->inputImg)->GetNy();
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
		
	std::map<long int,long int> mapping;
	std::vector< std::vector<bool> > isNeighborLinkTaken;

	SLICContourData* contourData= new SLICContourData;
	(contourData->contour)= new Contour;
	
	for(long int k=0;k<nRegions;k++){
	  long int regionId= (slicData->regions)[k]->Id;
		mapping.insert( std::pair<long int,long int>(regionId,k) );
		contourData->connectedRegionIds.push_back( std::vector<int>() );
		contourData->boundaryData.push_back( SLICBoundaryPixMap() );
		isNeighborLinkTaken.push_back( std::vector<bool>() );
		for(int s=0;s<nRegions;s++) isNeighborLinkTaken[k].push_back(false);
	}//end loop regions 

  //## Loop over all the pixels
  for (long int i=0; i<Nx; i++) {
  	for (long int j=0; j<Ny; j++) {
			long int index= i+j*Nx;

			long int regionId= (slicData->pixel_labels)[index];
			if(regionId<0) continue;
			long int regionIndex= mapping[regionId];			
   		int nr_p = 0;
			std::map<int,int> connectedRegionCounterMap;
			connectedRegionCounterMap.clear();
			std::map<int,std::vector<int>> connectedRegionBoundaryPixelListMap;
			connectedRegionBoundaryPixelListMap.clear();
            
      // Compare the pixel to its 8 neighbours
      for (int k=0; k<8; k++) {
				int x = i + dx8[k];
				int y = j + dy8[k];
				         
        if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
					long int index_neighbor= x + y*Nx;
					long int regionId_neighbor= (slicData->pixel_labels)[index_neighbor];
       	
        	if (regionId!=regionId_neighbor && regionId_neighbor>=0) {	
          	nr_p++;
						connectedRegionCounterMap[regionId_neighbor]++;
						//long int gBin_neighbor= inputImg->GetBin(x,y);
						connectedRegionBoundaryPixelListMap[regionId_neighbor];
						//connectedRegionBoundaryPixelListMap[regionId_neighbor].push_back(gBin_neighbor);
						connectedRegionBoundaryPixelListMap[regionId_neighbor].push_back(index_neighbor);
          }
        }//close if
      }//end loop neighbours
            
      // Add the pixel to the contour list of corresponding region
     	if (nr_p>= 2) {
				//long int gBin= inputImg->GetBin(i,j);
				//double binX= (slicData->inputImg)->GetX(i);
				//double binY= (slicData->inputImg)->GetY(j);	
				//(contourData->contour)->AddPoint(TVector2(binX,binY));

				long int gBin= index;
				(contourData->contour)->AddPoint(TVector2(i,j));
      	
				//Add neighbor ids and boundary pixel ids
				std::map<int,std::vector<int>>::iterator counterListIt = connectedRegionBoundaryPixelListMap.begin();
				for (counterListIt=connectedRegionBoundaryPixelListMap.begin(); counterListIt!=connectedRegionBoundaryPixelListMap.end(); ++counterListIt){
					int neighborId= counterListIt->first;
					std::vector<int> neighborPixIds= counterListIt->second;
					int counts= (int)neighborPixIds.size();	
					int neighborIndex= mapping[neighborId];

					if(counts>=2) {		
						std::vector<int>::iterator it = std::find (contourData->connectedRegionIds[regionIndex].begin(), contourData->connectedRegionIds[regionIndex].end(), neighborIndex);
						
						if( contourData->connectedRegionIds[regionIndex].empty() || it==contourData->connectedRegionIds[regionIndex].end() ) {
							contourData->connectedRegionIds[regionIndex].push_back(neighborIndex);
						}//close if	

						(contourData->boundaryData[regionIndex])[neighborId];//add connection with neighbor if not existing in the map
						((contourData->boundaryData[regionIndex])[neighborId]).push_back(gBin);//add this contour pixel
						for(int t=0;t<counts;t++){
							int gBin_neighbor= neighborPixIds[t];
							it = std::find ( ((contourData->boundaryData[regionIndex])[neighborId]).begin(), ((contourData->boundaryData[regionIndex])[neighborId]).end(), gBin_neighbor);
							if( ((contourData->boundaryData[regionIndex])[neighborId]).empty() || it==((contourData->boundaryData[regionIndex])[neighborId]).end() ) {
								((contourData->boundaryData[regionIndex])[neighborId]).push_back(gBin_neighbor);
							}//close if
						}//end loop counts
						
					}//close if counts>2
				}//end loop map counter iterator
      }//close if nr_p>2
    }//end loop image Ny
  }//end loop image Nx
	
	return contourData;

}//close ComputeBoundaryContours()




int SLIC::FindNeighbors(std::vector<SLICNeighborCollection>& neighbors,SLICData* slicData,SLICContourData* contourData,bool get2ndNeighbors,int selectedTag,bool includeSpatialDist,bool normalizeParams,bool useRobustParams,bool addCurvDist){
	
	//## Check input
	if(!slicData || !contourData){
		ERROR_LOG("Null ptr to SLIC data given!");
		return -1;
	}
	
	//## Check slic data image 
	if(!slicData->edgeImg) {
		ERROR_LOG("Missing (nullptr) edge image in slic data!");
		return -1;
	}
	if(!slicData->inputImg || !slicData->laplImg) {
		ERROR_LOG("Missing (nullptr) input and/or laplacian image in slic data!");
		return -1;
	}

	
	//## Get image & curvature ranges
	double Xmin= (slicData->inputImg)->GetXmin();
	double Xmax= (slicData->inputImg)->GetXmax();
	double Ymin= (slicData->inputImg)->GetYmin();
	double Ymax= (slicData->inputImg)->GetYmax();
	double width= fabs(Xmax-Xmin);
	double height= fabs(Ymax-Ymin);
	double diagonal2= width*width + height*height;

	RegionDistNormData normPars;
	normPars.ImgDiagonal= sqrt(diagonal2);	
	normPars.Smin= (slicData->inputImg)->GetMinimum();
	normPars.Smax= (slicData->inputImg)->GetMaximum();
	normPars.Smin_curv= (slicData->laplImg)->GetMinimum();
	normPars.Smax_curv= (slicData->laplImg)->GetMaximum();
	normPars.NormMin= 0;
	normPars.NormMax= 1;

	//## Check number of regions
	int nRegions= slicData->GetNRegions();
	if(nRegions<=0){
		WARN_LOG("No regions available, nothing to be done!");
		return 0;
	}
	
	//## Fill list of 1st and 2-nd neighbors for each region
	//## Include only selected tags 	
	//SLICConnectedRegions connectedRegionIds= contourData->connectedRegionIds;
	std::vector<SLICBoundaryPixMap> boundaryData= contourData->boundaryData;
	std::vector< std::vector<long int> > neighborIndexList;
	double Emax= 1;

	for(int i=0; i<nRegions; i++) {
		long int regionId= (slicData->regions)[i]->Id;			
		neighborIndexList.push_back( std::vector<long int>() );
		SLICNeighborCollection aNeighborCollection;
		neighbors.push_back(aNeighborCollection);

		//Fill 1st-order neighbors
		for(unsigned int j=0;j<(contourData->connectedRegionIds)[i].size();j++){//loop over 1st neighbors
			int neighborIndex= (contourData->connectedRegionIds)[i][j];
			int neighborTag= (slicData->regions)[neighborIndex]->Tag;
			long int neighborId= (slicData->regions)[neighborIndex]->Id;

			if(selectedTag==-1 || neighborTag==selectedTag){

				//Find symmetric dissimilarity
				double Dsym= 0;
				double Dsym_space= 0;
				if(ComputeRegionDistance(Dsym,Dsym_space,(slicData->regions)[i],(slicData->regions)[neighborIndex],normPars,normalizeParams,useRobustParams,addCurvDist)<0){
					ERROR_LOG("Symm distance calculation failed!");
					return -1;
				}		

				//Find asymm dissimilarities
				double Diss= 0;
				double DissNeighbor= 0;
				if( ComputeRegionAsymmDistance(Diss, DissNeighbor, (slicData->regions)[i],(slicData->regions)[neighborIndex], normPars,normalizeParams,useRobustParams,addCurvDist,includeSpatialDist)<0)
				{
					ERROR_LOG("Asymmetric distance calculation failed!");
					return -1;
				}

				
				//Find edgeness
				SLICBoundaryPixMapIterator it= (contourData->boundaryData)[i].find(neighborId);
				double E= 0;
				if( it != (contourData->boundaryData)[i].end() ){//region is among neighbors, compute edgeness terms
					std::vector<int> sharedPixelIds= ((contourData->boundaryData)[i])[neighborId];
					int nBoundaryPixels= (int)sharedPixelIds.size();			
					for(int t=0;t<nBoundaryPixels;t++) {
						int gBin= sharedPixelIds[t];
						double S_edge= (slicData->edgeImg)->GetBinContent(gBin);
						E+= S_edge;
					}	
					if(nBoundaryPixels>0) E/= (double)nBoundaryPixels;  
				}//close if
	
				SLICBoundaryPixMapIterator it_neighbor= (contourData->boundaryData)[neighborIndex].find(regionId);
				double E_neighbor= 0;

				if( it_neighbor != (contourData->boundaryData)[neighborIndex].end()){//region is among neighbors, compute edgeness terms
					std::vector<int> sharedPixelIds= ((contourData->boundaryData)[neighborIndex])[regionId];
					int nBoundaryPixels= (int)sharedPixelIds.size();			
					for(int t=0;t<nBoundaryPixels;t++) {
						int gBin= sharedPixelIds[t];
						double S_edge= (slicData->edgeImg)->GetBinContent(gBin);
						E_neighbor+= S_edge;
					}	
					if(nBoundaryPixels>0) E_neighbor/= (double)nBoundaryPixels;  
				}//close if

				SLICNeighborData neighborData;
				neighborData.Order= 1;
				neighborData.Id= neighborId;
				neighborData.Index= neighborIndex;
				neighborData.Tag= neighborTag;
				neighborData.D= Diss;
				neighborData.D_n= DissNeighbor;
				neighborData.E= E;	
				neighborData.E_n= E_neighbor;
				neighborData.Dsym= Dsym;

				neighborIndexList[i].push_back(neighborIndex);
				neighbors[i].Add(neighborData);
			}//close if selected region tag

			//Fill 2nd-order neighbors?
			if(!get2ndNeighbors) continue;

			for(unsigned int t=0;t<(contourData->connectedRegionIds)[neighborIndex].size();t++){//loop over 2nd neighbors
				int neighborIndex_2nd= (contourData->connectedRegionIds)[neighborIndex][t];
				int neighborTag_2nd= (slicData->regions)[neighborIndex_2nd]->Tag;
				int neighborId_2nd= (slicData->regions)[neighborIndex_2nd]->Id;
			
				if( (selectedTag!=-1 || neighborTag_2nd==selectedTag) && neighborId_2nd!=regionId ){
					std::vector<long int>::iterator it= std::find(neighborIndexList[i].begin(),neighborIndexList[i].end(),neighborIndex_2nd);
						
					if( neighborIndexList[i].size()==0 || it==neighborIndexList[i].end() ){

						//Find symm diss for 2nd neighbors
						double Dsym_2nd= 0;
						double Dsym_space_2nd= 0;
						if(ComputeRegionDistance(Dsym_2nd,Dsym_space_2nd,(slicData->regions)[i],(slicData->regions)[neighborIndex_2nd],normPars,normalizeParams,useRobustParams,addCurvDist)<0){
							ERROR_LOG("Symm distance calculation failed!");
							return -1;
						}	
						
						//Find asymm diss for 2nd neighbors
						double Diss_2nd= 0; 
						double DissNeighbor_2nd= 0;
						if( ComputeRegionAsymmDistance(Diss_2nd, DissNeighbor_2nd, (slicData->regions)[i],(slicData->regions)[neighborIndex_2nd], normPars,normalizeParams,useRobustParams,addCurvDist,includeSpatialDist)<0)
						{
							ERROR_LOG("Asymmetric distance calculation failed!");
							return -1;
						}

						SLICNeighborData neighborData_2nd;
						neighborData_2nd.Order= 2;
						neighborData_2nd.Id= neighborId_2nd;
						neighborData_2nd.Index= neighborIndex_2nd;
						neighborData_2nd.Tag= neighborTag_2nd;
						neighborData_2nd.D= Diss_2nd;
						neighborData_2nd.D_n= DissNeighbor_2nd;
						neighborData_2nd.E= Emax;	
						neighborData_2nd.E_n= Emax;
						neighborData_2nd.Dsym= Dsym_2nd;
						neighbors[i].Add(neighborData_2nd);

						neighborIndexList[i].push_back(neighborIndex_2nd);
					}//close if item not found in collection	
				}//close if
			}//end loop 2nd-neighbors
		}//end loop 1st-neighbors
	}//end loop regions

	neighborIndexList.clear();

	return 0;

}//close FindNeighbors()

int SLIC::ComputeRegionDistance(double& Dsym,double& Dsym_space,Region* region_i,Region* region_j, RegionDistNormData normPars, bool normalizeParams,bool useRobustParams,bool addCurvDist)
{
	//Init
	Dsym= 0;
	Dsym_space= 0;	
	
	//Check regions
	if(!region_i || !region_j){
		ERROR_LOG("Null ptr to given regions!");
		return -1;
	}

	//Compute un-normalized distance squared
	DistPars distPars;
	if(region_i->GetDistance(distPars,region_j,useRobustParams)<0){
		ERROR_LOG("Symm distance calculation failed!");
		return -1;	
	}

	double dist2_color= distPars.dist2;//color distance squared
	double dist2_curv= distPars.dist2_curv;//curvature distance squared
	double dist2_spatial= distPars.dist2_spatial;//spatial distance squared
	
	//Normalize distances?
	if(normalizeParams){
		double dist2_color_norm= Caesar::StatsUtils::ComputeNormDiffSqr(dist2_color,normPars.Smin,normPars.Smax,normPars.NormMin,normPars.NormMax);
		double dist2_curv_norm= Caesar::StatsUtils::ComputeNormDiffSqr(dist2_curv,normPars.Smin_curv,normPars.Smax_curv,normPars.NormMin,normPars.NormMax);
		double dist2_spatial_norm= dist2_spatial/pow(normPars.ImgDiagonal,2);

		double dist2= dist2_color_norm;
		if(addCurvDist) dist2+= dist2_curv_norm;				
		
		Dsym= sqrt(dist2);
		Dsym_space= sqrt(dist2_spatial_norm);
	}
	else{
		double dist2= dist2_color;
		if(addCurvDist) dist2+= dist2_curv;
		
		Dsym= sqrt(dist2);
		Dsym_space= sqrt(dist2_spatial);
	}

	return 0;

}//close ComputeRegionDistance()



int SLIC::ComputeRegionAsymmDistance(double& Diss,double& DissNeighbor,Region* region_i,Region* region_j, RegionDistNormData normPars, bool normalizeParams,bool useRobustParams,bool addCurvDist,bool addSpatialDist)
{
	//Init
	Diss= 0;
	DissNeighbor= 0;	
	
	//Check regions
	if(!region_i || !region_j){
		ERROR_LOG("Null ptr to given regions!");
		return -1;
	}

	//Compute un-normalized distance squared
	DistPars asymmDistPars;
	DistPars asymmDistPars_neighbor;
	if(region_i->GetAsymmDistance(asymmDistPars,asymmDistPars_neighbor,region_j,useRobustParams)<0){
		ERROR_LOG("Asymm distance calculation failed!");
		return -1;						
	}

	double asymm_dist2_color= asymmDistPars.dist2;//color asymm distance squared
	double asymm_dist2_curv= asymmDistPars.dist2_curv;//curvature asymm distance squared
	double asymm_dist2_spatial= asymmDistPars.dist2_spatial;//spatial asymm distance squared
	
	double asymm_dist2_neighbor_color= asymmDistPars_neighbor.dist2;//color asymm distance squared
	double asymm_dist2_neighbor_curv= asymmDistPars_neighbor.dist2_curv;//curvature asymm distance squared
	double asymm_dist2_neighbor_spatial= asymmDistPars_neighbor.dist2_spatial;//spatial asymm distance squared
	
	//Normalize distances?
	if(normalizeParams){
		double asymm_dist2_color_norm= Caesar::StatsUtils::ComputeNormDiffSqr(asymm_dist2_color,normPars.Smin,normPars.Smax,normPars.NormMin,normPars.NormMax);
		double asymm_dist2_curv_norm= Caesar::StatsUtils::ComputeNormDiffSqr(asymm_dist2_curv,normPars.Smin,normPars.Smax,normPars.NormMin,normPars.NormMax);
		double asymm_dist2_spatial_norm= asymm_dist2_spatial/pow(normPars.ImgDiagonal,2);
		
		double asymm_dist2_neighbor_color_norm= Caesar::StatsUtils::ComputeNormDiffSqr(asymm_dist2_neighbor_color,normPars.Smin_curv,normPars.Smax_curv,normPars.NormMin,normPars.NormMax);
		double asymm_dist2_neighbor_curv_norm= Caesar::StatsUtils::ComputeNormDiffSqr(asymm_dist2_neighbor_curv,normPars.Smin_curv,normPars.Smax_curv,normPars.NormMin,normPars.NormMax);
		double asymm_dist2_neighbor_spatial_norm= asymm_dist2_neighbor_spatial/pow(normPars.ImgDiagonal,2);
		
		double asymm_dist2= asymm_dist2_color_norm;
		double asymm_dist2_neighbor= asymm_dist2_neighbor_color_norm;
		if(addCurvDist) {
			asymm_dist2+= asymm_dist2_curv_norm;
			asymm_dist2_neighbor+= asymm_dist2_neighbor_curv_norm;
		}
		if(addSpatialDist){
			asymm_dist2+= asymm_dist2_spatial_norm;
			asymm_dist2_neighbor+= asymm_dist2_neighbor_spatial_norm;
		}
		
		Diss= sqrt(asymm_dist2);
		DissNeighbor= sqrt(asymm_dist2_neighbor);

	}//close normalizePars
	else{
		double asymm_dist2= asymm_dist2_color;
		double asymm_dist2_neighbor= asymm_dist2_neighbor_color;

		if(addCurvDist) {
			asymm_dist2+= asymm_dist2_curv;
			asymm_dist2_neighbor+= asymm_dist2_neighbor_curv;	
		}
		if(addSpatialDist){
			asymm_dist2+= asymm_dist2_spatial;
			asymm_dist2_neighbor+= asymm_dist2_neighbor_spatial;
		}
		Diss= sqrt(asymm_dist2);
		DissNeighbor= sqrt(asymm_dist2_neighbor);

	}//close else

	return 0;

}//close ComputeRegionAsymmDistance()



SLICSimilarityData* SLIC::ComputeRegionSimilarity(SLICData* slicData,std::vector<SLICNeighborCollection>& neighbors,double beta){

	//## Check input data
	if(!slicData) {
		ERROR_LOG("Null ptr to given slic data!");
		return nullptr;
	}

	int nRegions= slicData->GetNRegions();
	if(nRegions<=0) {
		WARN_LOG("No regions present in slic data, returning nullptr!");
		return nullptr;
	}
	INFO_LOG("Compute region similarities (nRegions="<<nRegions<<")");

	if(!slicData->edgeImg) {
		ERROR_LOG("No edge image stored in slic data!");
		return nullptr;
	}

	double Emin_norm= 0;
	double Emax_norm= 1;
	double Dmin_norm= 0;
	double Dmax_norm= 1;
	
	//## Init matrix
	TMatrixD* AdjacencyMatrix= new TMatrixD(nRegions,nRegions);
	AdjacencyMatrix->Zero();

	//## Fill matrix
	double Dmin= 1.e+99;
	double Dmax= -1.e+99;
	double Emin= 1.e+99;
	double Emax= -1.e+99;
	std::vector<double> Dlist;

	for(size_t i=0;i<neighbors.size();i++){
		std::vector<SLICNeighborData> thisNeighbors= neighbors[i].GetNeighbors();
		for(size_t j=0;j<thisNeighbors.size();j++){//loop over neighbor list for this region
			
			int neighborIndex= thisNeighbors[j].Index;
			double D= thisNeighbors[j].D;
			double E= thisNeighbors[j].E;

			if(D>Dmax) Dmax= D;
			if(D<Dmin) Dmin= D;
			if(E>Emax) Emax= E;
			if(E<Emin) Emin= E;
			Dlist.push_back(D);
		}//end loop neighbor regions
	}//end loop regions

	//Check min & max
	if(Dmax<=Dmin || Emax<=Emin){
		WARN_LOG("Invalid normalization values for dissimilarity (min/max="<<Dmin<<","<<Dmax<<") and/or edgeness (min/max="<<Emin<<"/"<<Emax<<")!");
		delete AdjacencyMatrix;
		AdjacencyMatrix= 0;
		return nullptr;
	}
	
	//Compute diss median & mad
	double Dmedian= StatsUtils::GetMedianFast<double>(Dlist);
	double Dmad= StatsUtils::GetMADFast(Dlist,Dmedian);
	double Dmedianrms= Dmad*1.4826;
	

	//## Normalize matrix
	const double SMALL_NUMBER= 0.00000000001;
	
	for(unsigned int i=0;i<neighbors.size();i++){
		std::vector<SLICNeighborData> thisNeighbors= neighbors[i].GetNeighbors();
		for(unsigned int j=0;j<thisNeighbors.size();j++){//loop over neighbor list for this region
			int neighborIndex= thisNeighbors[j].Index;
			int neighborOrder= thisNeighbors[j].Order;
			double D= thisNeighbors[j].D;
			double D_n= thisNeighbors[j].D_n;
			double E= thisNeighbors[j].E;
			double E_n= thisNeighbors[j].E_n;

			double D_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D-Dmin)/(Dmax-Dmin);
			double D_n_norm= Dmin_norm + (Dmax_norm-Dmin_norm)*(D_n-Dmin)/(Dmax-Dmin);
			
			double E_norm= Emin_norm + (Emax_norm-Emin_norm)*(E-Emin)/(Emax-Emin);
			double E_n_norm= Emin_norm + (Emax_norm-Emin_norm)*(E_n-Emin)/(Emax-Emin);
			double Dtot= (1-beta)*D_norm + beta*E_norm;
			double Dtot_n= (1-beta)*D_n_norm + beta*E_n_norm;
			//double Dtot= D_norm + beta*E_norm;
			Dtot+= SMALL_NUMBER;//to avoid dividing by zero!
			Dtot_n+= SMALL_NUMBER;//to avoid dividing by zero!

			//(neighbors[i].GetNeighbor(j))->Dtot= Dtot;//Set Dtot in current neighbor
			//(neighbors[i].GetNeighbor(j))->Dtot_n= Dtot_n;//Set Dtot in current neighbor
			neighbors[i].SetDtot(j,Dtot,Dtot_n);

			//(*DissimilarityMatrix)(i,neighborIndex)= Dtot;
			if(neighborOrder>0) (*AdjacencyMatrix)(i,neighborIndex)= 1./Dtot;	

		}//end loop neighbor regions
	}//end loop regions
		
	//## Normalize adjacency matrix by rows
	for(int i=0;i<AdjacencyMatrix->GetNrows();i++){
		double sum= 0;	
		for(int j=0;j<AdjacencyMatrix->GetNcols();j++) {
			sum+= (*AdjacencyMatrix)(i,j);
		}
		if(sum!=0) {
			for(int j=0;j<AdjacencyMatrix->GetNcols();j++) (*AdjacencyMatrix)(i,j)/= sum;
		}
	}//end loop rows		
			
				
	//## Return data
	SLICSimilarityData* SimilarityData= new SLICSimilarityData;
	//SimilarityData->DissimilarityMatrix= DissimilarityMatrix;
	SimilarityData->AdjacencyMatrix= AdjacencyMatrix;  
	SimilarityData->Dmin= Dmin;
	SimilarityData->Dmax= Dmax;
	SimilarityData->Dmedian= Dmedian;
	SimilarityData->Dmedianrms= Dmedianrms;
	SimilarityData->Emin= Emin;
	SimilarityData->Emax= Emax;
	
	return SimilarityData;

}//close ComputeRegionSimilarity()




int SLIC::TagRegions(std::vector<Region*>& regions,Image* binaryMap_bkg,Image* binaryMap_signal){

	if(!binaryMap_bkg || !binaryMap_signal ) {
		ERROR_LOG("No binary maps provided!");
		return -1;
	}
	
	int nRegions= (int)regions.size();
	if(nRegions<=0) {
		ERROR_LOG("No regions available, nothing to be tagged!");
		return -1;
	}
	
	DEBUG_LOG("Tag regions...");
	int nBkgReg= 0;
	int nSignalReg= 0;
	int nUntaggedReg= 0;

	for(int i=0;i<nRegions;i++){	
		int nPix= regions[i]->NPix;
		int nBkg= 0;
		int nUntagged= 0;
		int nSignal= 0;
		if(nPix<=0) {
			regions[i]->Tag= Region::eUntagged;	
			//regions[i]->Saliency= 0;
			continue;
		}

		//double regionMean= regions[i]->Mean;
		//double averageSaliency= 0;

		for(int j=0;j<nPix;j++){//loop on pixels inside region
			long int pixId= (regions[i]->GetPixel(j))->id;
			//double saliency= saliencyMap->GetBinContent(pixId);	
			//averageSaliency+= saliency;

			double S_signal= binaryMap_signal->GetPixelValue(pixId);
			double S_bkg= binaryMap_bkg->GetPixelValue(pixId);
			bool isSignal= (S_signal>0);
			bool isBkg= (S_bkg>0);

			if(isBkg && !isSignal) nBkg++;
			else if(isSignal && !isBkg) nSignal++;
			else nUntagged++;
		}//end loop pixels in region
		
		//averageSaliency/= (double)nPix;
		//regions[i]->Saliency= averageSaliency;
	
		//Tag using a majority criterion
		regions[i]->Tag= Region::eUntagged;
		if(nSignal>nBkg && nSignal>nUntagged) {
			regions[i]->Tag= Region::eSignalTag;
			nSignalReg++;
		}
		else if(nBkg>nSignal && nBkg>nUntagged) {
			regions[i]->Tag= Region::eBkgTag;
			nBkgReg++;
		}
		else if(nUntagged>nSignal && nUntagged>nBkg) {
			regions[i]->Tag= Region::eUntagged;
			nUntaggedReg++;
		}

	}//end loop regions

	INFO_LOG("(nS,nB,nU)=("<<nSignalReg<<","<<nBkgReg<<","<<nUntaggedReg<<")");

	return 0;

}//close TagRegions()

int SLIC::CountTaggedRegions(std::vector<Region*>const& regions,int& nSig,int& nBkg,int& nUntagged){

	//Check regions
	int nRegions= (int)regions.size();
	if(nRegions<=0){
		nSig= 0;
		nBkg= 0;
		nUntagged= 0;
		ERROR_LOG("No regions given!");
		return -1;
	}

	nBkg= 0;
	nSig= 0;
	nUntagged= 0;
	for(unsigned int i=0;i<regions.size();i++) {
		int tag= regions[i]->Tag;
		if(tag==Region::eBkgTag) nBkg++;
		else if(tag==Region::eSignalTag) nSig++;
		else if(tag==Region::eUntagged) nUntagged++;
	}//end loop regions
			
	return 0;

}//close CountTaggedRegions()




}//close namespace

