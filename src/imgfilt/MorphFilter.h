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
* @file MorphFilter.h
* @class MorphFilter
* @brief Class implementing morphological filtering
*
* Morphological Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _MORPH_FILTER_h
#define _MORPH_FILTER_h 1

//ROOT
#include <TObject.h>
#include <TVector2.h>

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


namespace Caesar{

class Source;
class ImgBkgData;
class Img;
class Image;

class MorphFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MorphFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~MorphFilter();

		enum DilationModel {eDilateWithBkg=1,eDilateWithSourceMedian=2};

	public:
			
		/**
		* \brief Find peaks in image by combining local peaks found with 3 different dilation kernel sizes
		*/
		static int FindPeaks(std::vector<TVector2>& peakPoints,Image* img,int peakShiftTolerance=1,bool skipBorders=true);

		/**
		* \brief Dilate image with specified kernel
		*/
		static Image* Dilate(std::vector<long int>& peakPixelIds,Image* img,int KernSize,bool skipBorders=true);

		/**
		* \brief Dilate image around specified sources position
		*/
		static int DilateAroundSources(Image* img,std::vector<Source*>const& sources,int KernSize=5,int dilateModel=eDilateWithBkg,int dilateSourceType=-1,bool skipToNested=false,ImgBkgData* bkgData=0,bool useLocalBkg=false,bool randomize=false,double zThr=1.e+6);
	
		
	private:

		/**
		* \brief Dilate image around a specified source position
		*/
		static int DilateAroundSource(Image* img,Source* source,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr);
		/**
		* \brief Find pixels to be dilated
		*/
		static int FindDilatedSourcePixels(Image* img,Source* source,int KernSize,std::vector<long int>& pixelsToBeDilated);
		

	private:

	ClassDef(MorphFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class MorphFilter+;
#endif

}//close namespace

#endif
