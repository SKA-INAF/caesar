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
* @file ZernikeMoments.h
* @class ZernikeMoments
* @brief ZernikeMoments
*
* @author S. Riggi
* @date 15/06/2015
*/


#ifndef _ZERNIKE_MOMENTS_H
#define _ZERNIKE_MOMENTS_H

#include <Image.h>


#include <TVector2.h>
#include <TGraph.h>
#include <TText.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>

namespace Caesar {

class Image;

class ZernikeMoments {

	public:
  	/**
		\brief Class constructor
		*/
    ZernikeMoments();
		/**
		\brief Class destructor
		*/
    ~ZernikeMoments();
		
		
	public:
    
		static std::vector<double> GetZernike2D_Direct(Image* img, double order, double radius);
		static std::vector<double> GetZernike2D (Image* img, double order, double rad);
		static std::vector<double> GetZernike2DOld(Image* img, double D, double R);
		static std::vector<double> mb_Znl(double *X, double *Y, double *P, int size, double D, double m10_m00, double m01_m00, double R, double psum);

	private:

};//close class ZernikeMoments

}//close namespace 

#endif

