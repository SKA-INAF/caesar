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
* @file CasaReader.cc
* @class CasaReader
* @brief CasaReader
*
* CASA Image Reader class
* @author S. Riggi
* @date 20/01/2015
*/

#include <CasaReader.h>
#include <Image.h>
#include <SysUtils.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <CodeUtils.h>

#include <TObject.h>

//- CASACORE headers
#include <casacore/images/Images/PagedImage.h>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <chrono>
using namespace std;

ClassImp(Caesar::CasaReader)

namespace Caesar {


CasaReader::CasaReader() {

}//close costructor


CasaReader::~CasaReader() {

}//close destructor


int CasaReader::Read(Caesar::Image& img,std::string filename,int ix_min,int ix_max,int iy_min,int iy_max,bool checkFile)
{
	//## Check dir
	if(checkFile && !SysUtils::CheckDir(filename)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("File path "<<filename<<" not found on filesystem!");
		#endif
		return -1;
	}

	//## Read CASA image and import it as PagedImage
	casa::PagedImage<casa::Float>* casa_img= 0;
	try { 
		casa_img= new casa::PagedImage<casa::Float>(filename.c_str());
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to import file "<<filename<<" as CASA PagedImage!");
		#endif
		return -1;
	}

	//...
	//...

	
	//## Delete paged image
	if(casa_img){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Deleting imported casa image...");
		#endif
		delete casa_img;
		casa_img= 0;
	}

	return 0;

}//close Read()

}//close namespace
