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

//WCSTOOLS (TO BE DEPRECATED)
//#include <wcs.h>

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

class TEllipse;

namespace Caesar {

class Image;
class ImgMetaData;
class ImgStatsData;
class Contour;
class WCS;

//========================================
//==   ASTRO UTILS
//========================================
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
		* \brief Get IAU naming from WCS string coordinates
		*/
		static int GetIAUCoords(std::string& iau,const std::string& s);
		

		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(double& xpos, double& ypos,WCS* wcs,double ix,double iy);
		/**
		* \brief Get WCS coordinates in string format corresponding to image coordinates 
		*/
		static int PixelToWCSStrCoords(std::string& wcs_str,WCS* wcs,double ix,double iy,int max_str_length=4096);


		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(Caesar::Image* image,WCS* wcs,double ix,double iy,double& xpos, double& ypos,bool useImageCoords=true);
		/**
		* \brief Get WCS coordinates in string format corresponding to image coordinates
		*/
		static int PixelToWCSStrCoords(Caesar::Image* image,WCS* wcs,double ix,double iy,std::string& wcs_str,bool useImageCoords=true,int max_str_length=4096); 

		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(Caesar::Image* image,double ix,double iy,double& xpos, double& ypos,int coordSystem=-1,bool useImageCoords=true);
		
		/**
		* \brief Get WCS coordinates in string format corresponding to image coordinates
		*/
		static int PixelToWCSStrCoords(Caesar::Image* image,double ix,double iy,std::string& wcs_pos,int coordSystem=-1,bool useImageCoords=true,int max_str_length=4096);
		

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

		
		/**
		* \brief Convert contour from pixel coordinates to sky coordinates
		*/
		static Contour* PixelToWCSContour(Contour* contour,WCS* wcs,int pixOffset=0);

		/**
		* \brief Convert contour list from pixel coordinates to sky coordinates
		*/
		static int PixelToWCSContours(std::vector<Contour*>& contours_wcs,std::vector<Contour*>const& contours,WCS* wcs,int pixOffset=0);
	
		/**
		* \brief Convert ellipse from pixel coordinates to sky coordinates
		*/
		static TEllipse* PixelToWCSEllipse(TEllipse* ellipse,WCS* wcs,int pixOffset=0);
	
		/**
		* \brief Convert ellipse from pixel coordinates to sky coordinates neglecting sky projection and assuming Euclidean distances
		*/
		static TEllipse* PixelToWCSEllipseSimple(TEllipse* ellipse,WCS* wcs,int pixOffset=0);
	
		/**
		* \brief Compute ellipse deconvolved from ellipse beam
		*/
		static TEllipse* GetBeamDeconvolvedEllipse(TEllipse* ellipse,TEllipse* beam);

		/**
		* \brief Compute ellipse pars deconvolved from ellipse beam
		*/
		static int GetBeamDeconvolvedEllipsePars(double& bmaj_deconv,double& bmin_deconv,double& bpa_deconv,double bmaj,double bmin,double bpa,double bmaj_beam,double bmin_beam,double bpa_beam);

		/**
		* \brief Compute WCS beam ellipse pars from pixel ellipse pars
		*/
		//static int GetBeamEllipsePars(double& x0_wcs,double& y0_wcs,double& bmaj_wcs,double& bmin_wcs,double& bpa_wcs,double x0,double y0,double bmaj,double bmin,double bpa);

		

		/**
		* \brief Get distance (in degrees) between two points on the sky using Haversine formula
		*/
		static double GetWCSPointDist_Haversine(double ra1,double dec1,double ra2,double dec2);

		/**
		* \brief Get distance (in degrees) between two points on the sky using Vincenty formula
		*/
		static double GetWCSPointDist_Vincenty(double ra1,double dec1,double ra2,double dec2);
		/**
		* \brief Get point bearing (in degrees) between two points on the sky 
		*/
		static double GetWCSPointBearing(double ra1,double dec1,double ra2,double dec2);

		/**
		* \brief Convert ellipse to DS9 format
		*/
		static std::string EllipseToDS9Region(TEllipse* ellipse,std::string text="",std::string color="white",std::vector<std::string> tags={},bool useImageCoords=true);
		
		/**
		* \brief Convert contour to DS9 format
		*/
		static std::string ContourToDS9Region(Contour* contour,std::string text="",std::string color="white",std::vector<std::string> tags={},bool useImageCoords=true);

		/**
		* \brief Returns DS9 WCS type header from flag
		*/
		static std::string GetDS9WCSTypeHeader(int coordSys);

	private:
	
		ClassDef(AstroUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class AstroUtils+;
#endif	

}//close namespace


#endif 
