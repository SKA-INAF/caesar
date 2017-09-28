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
* @file FITSReader.h
* @class FITSReader
* @brief FITSReader
*
* Image Reader class for FITS files
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _FITS_READER_h
#define _FITS_READER_h 1

#include <SysUtils.h>
#include <StatsUtils.h>

#include <TObject.h>
#include <TFITS.h>
#include <TMath.h>


//CFITSIO headers
#include <fitsio.h>

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

using namespace std;


namespace Caesar {


class ImgMetaData;
class Image;

class FITSHeader : public TObject {

	public:
		FITSHeader(){};
		virtual ~FITSHeader(){};
		
	public:
		void Print(){
			cout<<"*** HEADER INFO ***"<<endl;
			cout<<"Image Size: "<<Nx<<"x"<<Ny<<" pixels, nRec="<<nRec<<endl;
			cout<<"Obs Coords: ("<<ObsRA<<","<<ObsDEC<<")"<<endl;
			cout<<"BUnit: "<<BUnit<<endl;
			cout<<"Coords Type: ("<<CoordTypeX<<","<<CoordTypeY<<")"<<endl;
			cout<<"PixelCoordCenter: ("<<Cx<<","<<Cy<<")"<<endl;
			cout<<"CoordCenter: ("<<Xc<<","<<Yc<<")"<<endl;
			cout<<"PixelStep: ("<<dX<<","<<dY<<")"<<endl;
			cout<<"BeamSize: ("<<Bmaj<<","<<Bmin<<","<<Bpa<<")"<<endl;
			cout<<"Rot: ("<<RotX<<","<<RotY<<")"<<endl;
			cout<<"Epoch: "<<Epoch<<endl;
			cout<<"***********************"<<endl;	
		}

	public:
		int Nx;
		int Ny;
		int nRec;
		double ObsRA;
		double ObsDEC;
		std::string BUnit;
		std::string CoordTypeX;
		std::string CoordTypeY;
		int Cx;
		int Cy;
		double Xc;
		double Yc;
		double dX;
		double dY;
		double Bmaj;
		double Bmin;
		double Bpa;
		double Epoch;
		double RotX;
		double RotY;

	ClassDef(FITSHeader,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class FITSHeader+;
#endif


class FITSFileInfo : public TObject {

	public:
		FITSFileInfo(){};
		virtual ~FITSFileInfo(){};
		
	public:
		void Print(bool printHeader=true){
			info.Print();	
			if(printHeader)	header.Print();
		}	
		void PrintHeader(){header.Print();}

	public:	
		FITSHeader header;
		Caesar::FileInfo info;

	ClassDef(FITSFileInfo,1)

};

#ifdef __MAKECINT__
#pragma link C++ class FITSFileInfo+;
#endif

class FITSReader : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    FITSReader();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~FITSReader();


	public:
	
		
		/**
		* \brief Read a FITS image & header and store it in Caesar img format (based on CFITSIO)
		*/
		static int Read(Caesar::Image& img,Caesar::FITSFileInfo& info,std::string filename,int hdu_id=1,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1,bool checkFile=true);
	
	private:

		
		/**
		* \brief Read header of currently open HDU (based on CFITSIO)
		*/		
		static int ReadHeader(Caesar::FITSFileInfo& fits_info,fitsfile* fp);
		/**
		* \brief Read image (based on CFITSIO)
		*/		
		static int ReadImage(Image& img,Caesar::FITSFileInfo& fits_info,fitsfile* fp,std::string filename,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1);
		
		#ifdef OPENMP_ENABLED
		/**
		* \brief Read image (multithread version) (based on CFITSIO)
		*/
		static int ReadImageMT(Image& img,Caesar::FITSFileInfo& fits_info,fitsfile* fp,std::string filename,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1);
		#endif

		/**
		* \brief Read image data
		*/	
		static int ReadAndFillImageData(Image& img,long int Nx,long int Ny,fitsfile* fp,int& read_data_status);
		static int ReadAndFillImageDataFast(Image& img,long int Nx,long int Ny,fitsfile* fp,int& read_data_status);
		
		#ifdef OPENMP_ENABLED
		/**
		* \brief Read image data (multithreaded version)
		*/	
		static int ReadAndFillImageDataMT(Image& img,long int Nx,long int Ny,fitsfile* fp,int& read_data_status,std::vector<Caesar::StatMoments<double>>& moments);
		static int ReadAndFillImageDataFastMT(Image& img,long int Nx,long int Ny,fitsfile* fp,int& read_data_status);
		#endif

		/**
		* \brief Internal method to handle errors occurred while processing cfitsio data structures
		*/
		static void HandleError(int& status,fitsfile* fp);

	private:

		ClassDef(FITSReader,1)

};

#ifdef __MAKECINT__
#pragma link C++ class FITSReader+; 
#endif

}//close namespace 


#endif


