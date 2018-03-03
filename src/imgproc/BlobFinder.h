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
* @file BlobFinder.h
* @class BlobFinder
* @brief Blob finder class
*
* Class to perform blob finding 
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _BLOB_FINDER_h
#define _BLOB_FINDER_h 1

#include <Consts.h>

#include <TObject.h>
#include <TMatrixD.h>

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
class ImgBkgData;
class Blob;
class Source;

class BlobFinder : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		BlobFinder();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~BlobFinder();

	public:

		/**
		* \brief Find blobs
		*/
		template <class T>
		static int FindBlobs(Image* inputImg,std::vector<T*>& blobs,Image* floodImg=0,ImgBkgData* bkgData=0,double seedThr=5,double mergeThr=2.6,int minPixels=10,bool findNegativeExcess=false,bool mergeBelowSeed=false,Image* curvMap=0);
		/**
		* \brief Flood fill algorithm
		*/
		static int FloodFill(Image* img,std::vector<long int>& clusterPixelIds,long int seedPixelId,double floodMinThr,double floodMaxThr);
		
		/**
		* \brief Get multiscale blob mask
		*/
		//static Image* GetMultiScaleBlobMask(Image* img,int kernelFactor,double sigmaMin,double sigmaMax,double sigmaStep,int thrModel=2,double thrFactor=1);

		/**
		* \brief Compute multiscale blob filtered map
		*/
		//static Image* ComputeMultiScaleBlobMap(Image* img,double sigmaMin,double sigmaMax,double sigmaStep,double thrFactor,int kernelFactor=1,double multiplicityThrFactor=0.7);

		/**
		* \brief Compute single-scale blob mask using elliptical gaus smoothing + laplacian filter map
		*/
		static Image* ComputeBlobMask(Image* img,double Bmaj,double Bmin,double Bpa=0,double kernNSigmaSize=2.5,double peakZThr=5,double peakZMergeThr=2.5,int minBlobSize=5,double thrFactor=0,int bkgEstimator=eMedianBkg,int bkgBox=100,double bkgGridStepSize=20);

		/**
		* \brief Compute multiscale blob mask using LoG filter maps
		*/
		static Image* ComputeMultiScaleBlobMask(Image* img,double sigmaMin,double sigmaMax,double sigmaStep,double peakZThr=5,double peakZMergeThr=2.5,int minBlobSize=5,double thrFactor=0,int kernelFactor=1,int bkgEstimator=eMedianBkg,int bkgBox=100,double bkgGridStepSize=20);

	public:
		
		ClassDef(BlobFinder,1)

};//close BlobFinder()



#ifdef __MAKECINT__
#pragma link C++ class BlobFinder+;
#pragma link C++ function FindBlob<Blob>(Caesar::Image*,std::vector<Blob*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool);
#pragma link C++ function FindBlob<Source>(Caesar::Image*,std::vector<Source*>&,Caesar::Image*,Caesar::ImgBkgData*,double,double,int,bool,bool);
#endif

}//close namespace

#endif

