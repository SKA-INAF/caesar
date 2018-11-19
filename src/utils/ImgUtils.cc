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
* @file ImgUtils.cc
* @class ImgUtils
* @brief ImgUtils
*
* Image utility class
* @author S. Riggi
* @date 12/07/2018
*/

#include <ImgUtils.h>
#include <Image.h>
#include <Logger.h>

//C++ headers
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>
#include <chrono>

using namespace std;


ClassImp(Caesar::ImgUtils)

namespace Caesar {

Image* ImgUtils::GetCircleLevelSetImage(long int nX,long int nY,double f)
{
	//Check inputs
	if(nX<=0 || nY<=0){
		ERROR_LOG("Invalid image dimensions given (must be >0)");
		return nullptr;
	}

	//Allocate image
	Image* img= 0;
	try{
		img= new Image(nX,nY);
	}
	catch(...){
		ERROR_LOG("Failed to allocate image of size ("<<nX<<" x "<<nY<<")!");
		return nullptr;
	}	
	
	//Fill image
	double R= std::min(nX,nY)*f;//circle radius
	long int Xc= nX/2;
	long int Yc= nY/2;
	
	#ifdef OPENMP_ENABLED
	#pragma omp for collapse(2)
	#endif
	for (long int i=0; i<nX; i++) {
    for (long int j=0; j<nY; j++) {
			double x= i - Xc;
			double y= j - Yc;
			double w = sqrt(x*x + y*y)-R;
			img->SetPixelValue(i,j,w);
    }//end loop cols
  }//end loop rows

	return img;

}//close GetCircleLevelSetImage()


Image* ImgUtils::GetCheckerBoardLevelSetImage(long int nX,long int nY,double f)
{
	//Check inputs
	if(nX<=0 || nY<=0){
		ERROR_LOG("Invalid image dimensions given (must be >0)");
		return nullptr;
	}

	//Allocate image
	Image* img= 0;
	try{
		img= new Image(nX,nY);
	}
	catch(...){
		ERROR_LOG("Failed to allocate image of size ("<<nX<<" x "<<nY<<")!");
		return nullptr;
	}	
	
	//Fill image
	double square_size= std::min(nX,nY)*f;//square size
	
	#ifdef OPENMP_ENABLED
	#pragma omp for collapse(2)
	#endif
	for (long int i=0; i<nX; i++) {
    for (long int j=0; j<nY; j++) {
			double xx= sin(TMath::Pi()/square_size*i);
			double yy= sin(TMath::Pi()/square_size*j);
			double w= xx*yy;
			img->SetPixelValue(i,j,w);
    }//end loop cols
  }//end loop rows

	return img;

}//close GetCheckerBoardLevelSetImage()




}//close namespace
