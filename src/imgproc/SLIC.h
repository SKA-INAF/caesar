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
* @file SLIC.h
* @class SLIC
* @brief SLIC generator class
*
* Superpixel generator
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _SLIC_h
#define _SLIC_h 1

#include <SLICData.h>
#include <TObject.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>


namespace Caesar {

class Image;


struct RegionDistNormData {

	RegionDistNormData(){
		ImgDiagonal= 0;
		Smin= 0;
		Smax= 0;
		Smin_curv= 0;
		Smax_curv= 0;
		NormMin= 0;
		NormMax= 1;
	}

	double ImgDiagonal;
	double Smin;
	double Smax;
	double Smin_curv;
	double Smax_curv;
	double NormMin;
	double NormMax;

};//close RegionDistNormData()

class SLIC : public TObject {

  public:
		
    /** 
		\brief Class constructor: initialize structures.
 		*/
		SLIC();
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~SLIC();


	public:


		/** 
		\brief Generate superpixel partition
 		*/
		static SLICData* SPGenerator(Image* img, int regionSize=20,double regParam=1, int minRegionSize=10, bool normalizeImage=true, bool useLogScaleMapping=false, Image* laplImg=0, Image* edgeImg=0);

		/** 
		\brief Compute superpixel boundary contours
 		*/
		static SLICContourData* ComputeBoundaryContours(SLICData* slicData);
		
		/** 
		\brief Compute region similarities
 		*/
		static SLICSimilarityData* ComputeRegionSimilarity(SLICData* slicData,std::vector<SLICNeighborCollection>& neighbors,double beta=0.5);

		/** 
		\brief Find superpixel neighbors
 		*/
		static int FindNeighbors(std::vector<SLICNeighborCollection>& neighbors,SLICData* slicData,SLICContourData* contourData,bool get2ndNeighbors=true,int selectedTag=-1,bool includeSpatialDist=false,bool normalizeParams=true,bool useRobustParams=false,bool addCurvDist=true);

		/** 
		\brief Compute segmented image given a list of tagged regions
 		*/
		static Image* GetSegmentedImage(Image* img,std::vector<Region*>const& regions,int selectedTag=-1,bool normalize=false,bool binarize=false);

		/** 
		\brief Count number of regions per tag
 		*/
		static int CountTaggedRegions(std::vector<Region*>const& regions,int& NSig,int& NBkg,int& NUntagged);		

		/** 
		\brief Tag regions into signal/bkg according to signal & bkg marker images
 		*/
		static int TagRegions(std::vector<Region*>& regions,Image* binaryMap_bkg,Image* binaryMap_signal);

		/** 
		\brief Compute distance between regions
 		*/
		static int ComputeRegionDistance(double& dist,double& dist_spatial,Region* region_i,Region* region_j,RegionDistNormData normPars, bool normalizeParams=true,bool useRobustParams=false,bool addCurvDist=false);

		/** 
		\brief Compute asymmetric distance between regions
 		*/
		static int ComputeRegionAsymmDistance(double& dist,double& dist_neighbor,Region* region_i,Region* region_j,RegionDistNormData normPars, bool normalizeParams=true,bool useRobustParams=false,bool addCurvDist=false,bool addSpatialDist=false);

		

	ClassDef(SLIC,1)

};

#ifdef __MAKECINT__
#pragma link C++ class SLIC+;
#endif

}//close namespace

#endif

