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
* @file AstroUtils.h
* @class AstroUtils
* @brief Utility functions for astronomical tasks
*
* Utility functions for astronomical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef _ASTRO_UTILS_h
#define _ASTRO_UTILS_h 1

#include <wcs.h>

#include <TObject.h>
#include <TMath.h>

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
#include <time.h>
#include <ctime>

using namespace std;


namespace Caesar {

class Image;
class ImgMetaData;
class ImgStatsData;

class AstroUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    AstroUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~AstroUtils();

		
	public:

		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(double& xpos, double& ypos,WorldCoor* wcs,double ix,double iy);

		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(Caesar::Image* image,WorldCoor* wcs,double ix,double iy,double& xpos, double& ypos);
		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(Caesar::Image* image,double ix,double iy,double& xpos, double& ypos,int coordSystem=-1);
		
		/**
		* \brief Get beam area from BMAJ, BMIN 
		*/
		static double GetBeamArea(double Bmaj,double Bmin){
			double A= TMath::Pi()*Bmaj*Bmin/(4*log(2));//2d gaussian area with FWHM=fx,fy
			return A;
		}

		/**
		* \brief Get beam area in pixels given Bmaj, Bmin and pixel sizes (dx, dy) in deg
		*/
		static double GetBeamAreaInPixels(double Bmaj,double Bmin,double dX,double dY){	
			double beamArea= GetBeamArea(Bmaj,Bmin);
			double pixelArea= fabs(dX*dY);
			double A= beamArea/pixelArea;
			return A;
		}

		/**
		* \brief Get beam width in pixels given Bmaj, Bmin and pixel sizes (dx, dy) in deg
		*/
		static int GetBeamWidthInPixels(double Bmaj,double Bmin,double dX,double dY){	
			double beamWidth= sqrt(fabs(Bmaj*Bmin));
			double pixWidth= sqrt(fabs(dX*dY));
			int beamWidthInPixel= static_cast<int>( ceil(beamWidth/pixWidth) );
			return beamWidthInPixel;
		}
	
	private:
	
		ClassDef(AstroUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class AstroUtils+;
#endif	

}//close namespace


#endif 
