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
* @file GradientFilter.h
* @class GradientFilter
* @brief Class implementing gradient filtering
*
* Gradient Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _GRADIENT_FILTER_h
#define _GRADIENT_FILTER_h 1

//ROOT
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

namespace cv {
	class Mat;
}

namespace Caesar{

class Image;

class GradientFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    GradientFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~GradientFilter();

	public:

		//===================================
		//==        NEW IMAGE METHODS
		//===================================
		/**
		* \brief Apply a gradient filter to image and return filtered image
		*/	
		static Image* GetGradientFilter(Image* image);
	
		/**
		* \brief Apply a laplacian filter to image and return filtered image
		*/
		static Image* GetLaplaceFilter(Image* image);


	private:

		/**
		* \brief Build the gradient filter kernel
		*/	
		static std::vector<cv::Mat> BuildGradientKernels();
		/**
		* \brief Build the Laplacian filter kernel
		*/
		static cv::Mat BuildLaplaceKernel();
		
	private:

	ClassDef(GradientFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class GradientFilter+;
#endif

}//close namespace

#endif
