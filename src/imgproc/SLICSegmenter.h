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
* @file SLICSegmenter.h
* @class SLICSegmenter
* @brief SLICSegmenter
*
* In this class, an over-segmentation is created of an image, provided by the
* step-size (distance between initial cluster locations) and the colour
* distance parameter.
* @author S. Riggi
* @date 15/06/2015
*/


#ifndef _SLIC_SEGMENTER_h
#define _SLIC_SEGMENTER_h 1


#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>

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

using namespace std;


namespace Caesar {

class SLICData;

class SLICSegmenter : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SLICSegmenter();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SLICSegmenter();


	public:
		/**
		\brief Merge the superpixels (main methods)
 		*/
		static int FindSegmentation(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,int minMergedSP=1,double SPMergingRatio=0.3,double SPMergingMaxDissRatio=1000,double SPMergingMaxDissRatio_2nd=1.05,double SPMergingDissThreshold=3,bool SPMergingIncludeSpatialPars=true,bool SPMergingUseRobustPars=false,bool SPMergingUseCurvDist=true);

		
	private:		
		/**
		* \brief Initialize and check data
		*/
		static int CheckData(SLICData const& slicData);
		/**
		\brief Adaptively merge the superpixels using max similarity
 		*/
		static int MultiStepSPMerger(SLICData const& slicData,SLICData& segmSlicData,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,bool SPMergingIncludeSpatialPars,bool SPMergingUseRobustPars,bool SPMergingUseCurvDist);
		/**
		\brief Hierarchical merge the superpixels 
 		*/
		static int SPHierarchicalMerger(SLICData& slicData,int mergerTag,int mergedTag,int minMergedSP,double SPMergingRatio,double SPMergingRegularization,bool use2ndNeighborsInSPMerging,double SPMergingMaxDissRatio,double SPMergingMaxDissRatio_2nd,double SPMergingDissThreshold,bool SPMergingIncludeSpatialPars,bool SPMergingUseRobustPars,bool SPMergingUseCurvDist);

		/**
		\brief Merge the superpixels according to max similarity
 		*/
		static int SPMaxSimilarityMerger(SLICData& segmSlicData,int mergerTag=1,int mergedTag=1,double SPMergingRegularization=0.5,bool use2ndNeighborsInSPMerging=false,bool SPMergingIncludeSpatialPars=false,bool SPMergingUseRobustPars=false,bool SPMergingUseCurvDist=true);		

		
	private:

	ClassDef(SLICSegmenter,1)

};//close SLICSegmenter()



#ifdef __MAKECINT__
#pragma link C++ class SLICSegmenter+;
#endif

}//close namespace

#endif

