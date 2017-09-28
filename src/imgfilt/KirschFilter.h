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
* @file KirschFilter.h
* @class KirschFilter
* @brief Class implementing Kirsch filtering
*
* Kirsch Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _KIRSCH_FILTER_h
#define _KIRSCH_FILTER_h 1

#include <TObject.h>

#include <vector>

namespace cv {
	class Mat;
}

namespace Caesar{

class Image;

class KirschFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    KirschFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~KirschFilter();

	public:
	
		//===================================
		//==        NEW IMAGE METHODS
		//===================================
		/**
		* \brief Apply Kirsh filter to input image and return filtered image
		*/
		static Image* GetKirschFilter(Image* image);
		

	private:

		/**
		* \brief Build kirsch kernels
		*/
		static std::vector<cv::Mat> BuildKernels();

	private:

	ClassDef(KirschFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class KirschFilter+;
#endif

}//close namespace

#endif
