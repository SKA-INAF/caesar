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
* @file LoGFilter.h
* @class LoGFilter
* @brief Class implementing Logarithm of Gaussian (LoG) filtering
*
* Logarithm of Gaussian (LoG) Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _LOG_FILTER_h
#define _LOG_FILTER_h 1

//ROOT
#include <TObject.h>

#include <vector>


namespace cv {
	class Mat;
}

namespace Caesar{

class Image;

class LoGFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    LoGFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~LoGFilter();

	public:
	

		//===================================
		//==        NEW IMAGE METHODS
		//===================================
		static Image* GetLoGFilter(Image* image);
		static Image* GetNormLoGFilter(Image* image,int size,double scale);


	private:
		/**
		* \brief Build standard filter kernel
		*/
		static cv::Mat BuildStandardKernel();
		/**
		* \brief Build filter kernel
		*/
		static cv::Mat BuildKernel(int kernSize,double scale);
		/**
		* \brief Normalized LoG kernel definition
		*/
		static double NormLoGKernel(double x,double y,double sigma);

	private:

	ClassDef(LoGFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class LoGFilter+;
#endif

}//close namespace

#endif
