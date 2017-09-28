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
* @file WTFilter.h
* @class WTFilter
* @brief Class implementing 2D Wavelet transform
*
* Wavelet Transform Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _WT_FILTER_h
#define _WT_FILTER_h 1


//ROOT
#include <TObject.h>

#include <vector>

namespace cv {
class Mat;
}

namespace Caesar{

class Image;

class WTFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    WTFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~WTFilter();

	public:
	
		//================================================
		//==         NEW IMAGE METHODS 
		//================================================
		/**
		* \brief Get Wavelet Transform decomposition of input image
		*/
		static std::vector<Image*> GetDecomposition(Image* image,int nScales);
		

	private:
		/**
		* \brief Build B3 spline kernel
		*/
		static cv::Mat BuildB3SlineKernel();
		
	private:

	ClassDef(WTFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class WTFilter+;
#endif

}//close namespace

#endif
