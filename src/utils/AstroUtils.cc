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
* @file AstroUtils.cc
* @class AstroUtils
* @brief Utility functions for astronomical tasks
*
* Utility functions for astronomical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <AstroUtils.h>
#include <Image.h>

#include <TObject.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::AstroUtils)

namespace Caesar {

AstroUtils::AstroUtils()
{

}

AstroUtils::~AstroUtils()
{

}


int AstroUtils::PixelToWCSCoords(Caesar::Image* image,WorldCoor* wcs,double ix,double iy,double& xpos, double& ypos) {

	//Check pixel values in input
	if(!image){
		cerr<<"AstroUtils::PixelToWCSCoords(): ERROR: Null image ptr given!"<<endl;
		return -1;	
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();

	if(ix<xmin || iy<ymin || ix>xmax || iy>ymax ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);

	return 0;

}//close PixelToWCSCoords()


int AstroUtils::PixelToWCSCoords(Caesar::Image* image,double ix,double iy,double& xpos, double& ypos,int coordSystem) {

	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();

	if(ix<xmin || iy<ymin || ix>xmax || iy>ymax ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check image meta-data
	if(!image->HasMetaData() ){
    ERROR_LOG("No metadata available in image!");
		return -1;
	}
	Caesar::ImgMetaData* metadata= image->GetMetaData();	
	
	//Get the coord system
	WorldCoor* wcs= metadata->GetWorldCoord(coordSystem);
	if(!wcs){
		ERROR_LOG("Failed to get WorldCoord system from metadata!");
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);

	//Clear up
	delete wcs;
	wcs= 0;

	return 0;
		
}//close PixelToWCSCoords()



}//close namespace



