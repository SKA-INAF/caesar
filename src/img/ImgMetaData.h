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
* @file ImgMetaData.h
* @class ImgMetaData
* @brief ImgMetaData
*
* Image metadata class
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _IMG_METADATA_h
#define _IMG_METADATA_h 1

#include <FITSReader.h>

//WCSTOOLS
#include <wcs.h>

#include <TObject.h>
#include <string>

namespace Caesar {

class ImgMetaData : public TObject {

	public:
		/**
		* \brief Constructor
		*/
		ImgMetaData();
		/**
		* \brief Destructor
		*/
		virtual ~ImgMetaData();

	public: 
		/**
		* \brief Set cards from FITS file
		*/
		void SetFITSCards(Caesar::FITSFileInfo& fits_info);

		/**
		* \brief Get current WCS type
		*/
		std::string GetWCSType(){return m_wcsType;}
		
		/**
		* \brief Get world coordinate system 
		*/
		WorldCoor* GetWorldCoord(int coordSystem=-1);

		/**
		* \brief Get pixel area in deg^2
		*/
		double GetPixelArea(){
			double pixelArea= fabs(dX*dY);
			return pixelArea;
		}

		/**
		* \brief Get width of synthetic beam
		*/
		double GetBeamWidth(){
			return sqrt(fabs(Bmaj*Bmin));
		}

		/**
		* \brief Get area of synthetic beam
		*/
		double GetBeamArea(){
			double fx= Bmaj;
			double fy= Bmin;
			double A= TMath::Pi()*fx*fy/(4*log(2));//2d gaussian area with FWHM=fx,fy
			return A;
		}
	
		/**
		* \brief Get pixel scale
		*/
		int GetBeamWidthInPixel(){
			double beamWidth= GetBeamWidth();
			double pixScale= sqrt(fabs(dX*dY));
			int beamWidthInPixel= int(ceil(beamWidth/pixScale));
			return beamWidthInPixel;
		}
		/**
		* \brief Get flux correction from beam
		*/
		double GetBeamFluxIntegral(){	
			double beamArea= GetBeamArea();
			double pixelArea= GetPixelArea();
			double f= beamArea/pixelArea;
			return f;
		}
	
	protected:
	
		/**
		* \brief Initialize fields
		*/
		void Init();

	public:
		//Image size
		int Nx;
		int Ny;

		//Reference pixel id
		int Cx;
		int Cy;

		//Reference pixel coords
		double Xc;
		double Yc;

		//Pixel size
		double dX;
		double dY;

		//System rotation info
		double RotX;
		double RotY;

		//Type of astro coords
		std::string CoordTypeX;
		std::string CoordTypeY;

		//Units
		std::string BUnit;

		//Beam info
		double Bmaj;
		double Bmin;
		double Bpa;

		//Obs Epoch
		double Epoch;

	private:
		//Type of WCS
		std::string m_wcsType;
		
	ClassDef(ImgMetaData,3)

};

#ifdef __MAKECINT__
#pragma link C++ class ImgMetaData+;
#endif


}//close namespace

#endif
