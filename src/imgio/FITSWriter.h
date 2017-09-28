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
* @file FITSWriter.h
* @class FITSWriter
* @brief FITSWriter
*
* Image writer class for FITS files
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _FITS_WRITER_h
#define _FITS_WRITER_h 1


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

using namespace std;

namespace Caesar {


class Image;
class ImgMetaData;

class FITSWriter : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    FITSWriter();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~FITSWriter();


	public:
		/**
		* \brief Write image to FITS file
		*/
		static int WriteFITS(Image* img,std::string outfilename);
		
	private:

		/**
		* \brief Initialize python interface
		*/
		static int Init();
		
	private:

		ClassDef(FITSWriter,1)

};

#ifdef __MAKECINT__
#pragma link C++ class FITSWriter+; 
#endif


}//close namespace


#endif
