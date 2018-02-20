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
* @file GausFilter.h
* @class GausFilter
* @brief Class implementing elliptical gaussian filtering
*
* Elliptical gaussian Filter
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _GAUS_FILTER_h
#define _GAUS_FILTER_h 1

//ROOT
#include <TObject.h>

#include <vector>


namespace cv {
	class Mat;
}

namespace Caesar{

class Image;

class GausFilter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    GausFilter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~GausFilter();

	public:
	

		/**
		* \brief Return gaussian-filtered image
		*/
		static Image* GetGausFilter(Image* image,double bmaj,double bmin,double bpa,int nSigmas=5,double scale=1);
		

	private:
		
		/**
		* \brief Build filter kernel
		*/
		static cv::Mat BuildKernel(int kernSize,double sigmaX,double sigmaY,double theta,double scale);
		/**
		* \brief Normalized LoG kernel definition
		*/
		static double Gaus2DFcn(double A,double x,double y,double sigmaX,double sigmaY,double theta);

	private:

	ClassDef(GausFilter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class GausFilter+;
#endif

}//close namespace

#endif
