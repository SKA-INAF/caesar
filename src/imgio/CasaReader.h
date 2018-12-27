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
* @file CasaReader.h
* @class CasaReader
* @brief CasaReader
*
* Image Reader class for CASA files
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _CASA_READER_h
#define _CASA_READER_h 1

#include <SysUtils.h>
#include <StatsUtils.h>

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
#include <string>


namespace Caesar {


class ImgMetaData;
class Image;

class CasaReader : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    CasaReader();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~CasaReader();


	public:

		/**
		* \brief Read a CASA image & header and store it in Caesar img format
		*/
		static int Read(Caesar::Image& img,std::string filename,int ix_min=-1,int ix_max=-1,int iy_min=-1,int iy_max=-1,bool checkFile=true);
	

		ClassDef(CasaReader,1)

};

#ifdef __MAKECINT__
#pragma link C++ class CasaReader+; 
#endif

}//close namespace 


#endif
