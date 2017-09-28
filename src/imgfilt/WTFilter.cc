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
* @file WTFilter.cc
* @class WTFilter
* @brief WTFilter
*
* Stationary Wavelet Transform filter class
* @author S. Riggi
* @date 20/01/2015
*/

#include <WTFilter.h>
#include <Image.h>
#include <MathUtils.h>
#include <WTFilter.h>

#include <TObject.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(Caesar::WTFilter)

namespace Caesar {

WTFilter::WTFilter() {

}//close costructor


WTFilter::~WTFilter(){

}//close destructor



//================================================
//==         NEW IMAGE METHODS 
//================================================
std::vector<Image*> WTFilter::GetDecomposition(Image* image,int nScales){
	
	//## Check image
	std::vector<Image*> imgCollection;
	imgCollection.clear();
	imgCollection.resize(0);
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return imgCollection;
	}

	//## Init kernels
	cv::Mat H= BuildB3SlineKernel();
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");
	int nrows= I.rows;
	int ncols= I.cols;
	
	//## Init filt mat
	std::vector<cv::Mat> F;
	std::vector<cv::Mat> W;//W(0) is the smoothed array
	for(int k=0;k<nScales+1;k++){
		F.push_back( cv::Mat::zeros(nrows,ncols,CV_64FC1) );
		W.push_back( cv::Mat::zeros(nrows,ncols,CV_64FC1) );
	}

	//## Compute convolution at all WT scales
	//Init F[0] with full reso input image
	F[0]= I;

	for(int n=1;n<nScales+1;n++){
		DEBUG_LOG("Compute convolution at scale "<<n<<" ...");
		F[n]= Caesar::MathUtils::GetATrousConvolution(F[n-1],H,n);
		W[n-1]= F[n-1]-F[n];
	}//end loop scales

	W[nScales]= F[nScales];

	//## Convert OpenCV mat list to Img
	Image* image_helper= 0;

	for(int n=0;n<nScales+1;n++){
		TString imgName= Form("%s_WT%d",image->GetName().c_str(),n);
		image_helper= image->GetCloned(std::string(imgName),true,true);
		image_helper->Reset();
		image_helper->FillFromMat(W[n]);
		imgCollection.push_back(image_helper);
	}//end loop scales
	
	return imgCollection;

}//close GetDecomposition()


cv::Mat WTFilter::BuildB3SlineKernel(){

	const int kernSize= 5;
	double kernValues[kernSize][kernSize]= { {1,4,6,4,1}, {4,16,24,16,4}, {6,24,36,24,6}, {4,16,24,16,4}, {1,4,6,4,1} };
	double kernScale= 256.;
	
	cv::Mat kernel(kernSize,kernSize,CV_64F);
	int rows = kernel.rows;
	int cols = kernel.cols;

	for (int i=0; i<rows;i++){
  	for (int j=0; j<cols;j++) { 
			kernel.at<double>(i,j)= kernValues[i][j]/kernScale;
		}
	}	
	return kernel;

}//close BuildB3SlineKernel()

}//close namespace 

