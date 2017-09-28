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
* @file KirschFilter.cc
* @class KirschFilter
* @brief Class implementing Kirsch filtering
*
* Kirsch Filter
* @author S. Riggi
* @date 20/01/2015
*/

#include <KirschFilter.h>
#include <Image.h>
#include <MathUtils.h>

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

ClassImp(Caesar::KirschFilter)

namespace Caesar {

KirschFilter::KirschFilter() {

}//close costructor


KirschFilter::~KirschFilter(){

}//close destructor



//===================================
//==        NEW IMAGE METHODS
//===================================
Image* KirschFilter::GetKirschFilter(Image* image)
{
	//## Check image
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return 0;
	}

	//## Init kernels
	std::vector<cv::Mat> kernelList= BuildKernels();
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");
	int nrows= I.rows;
	int ncols= I.cols;
	
	//## Compute convolutions for all kernel directions
	std::vector<cv::Mat> filteredEnsemble;
	cv::Mat filteredMat= cv::Mat::zeros(nrows,ncols,CV_64FC1);
	for(unsigned int k=0;k<kernelList.size();k++){
		filteredEnsemble.push_back( cv::Mat::zeros(nrows,ncols,CV_64FC1) );
		filteredEnsemble[k]= Caesar::MathUtils::GetConvolution(I,kernelList[k]);	
	}//end loop kernel directions

	//## Set final filtered image to maximum pixel value across all direction convolutions
	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for(int i=0;i<nrows;i++){
		for(int j=0;j<ncols;j++){
			double maxVal= -1.e+99;
			for(unsigned int k=0;k<kernelList.size();k++){
				double val= filteredEnsemble[k].at<double>(i,j);
				if(val>maxVal) maxVal= val;
			}	
			filteredMat.at<double>(i,j)= maxVal;
		}//end loop cols
	}//end loop rows


	//## Convert OpenCV mat list to Img
	TString imgName= Form("%s_Kirsch",image->GetName().c_str());
	Image* filteredImg= image->GetCloned(std::string(imgName),true,true);
	filteredImg->Reset();
	filteredImg->FillFromMat(filteredMat);

	return filteredImg;

}//close GetKirschFilter()



std::vector<cv::Mat> KirschFilter::BuildKernels(){

	//## Init Kirsch kernel values (8 kernels in different directions)	
	const int nkernels= 8;
	const int kernSize= 3;
	double kernValues[nkernels][kernSize][kernSize]= {
		{ {-3,-3,-3}, {-3,0,-3}, {5,5,5} },	//N
		{ {5,5,5}, {-3,0,-3}, {-3,-3,-3} }, //S
		{ {-3,-3,5}, {-3,0,5}, {-3,-3,5} }, //W
		{ {5,-3,-3}, {5,0,-3}, {5,-3,-3} }, //E
		{ {-3,-3,-3}, {-3,0,5}, {-3,5,5} }, //NW
		{ {-3,-3,-3}, {5,0,-3}, {5,5,-3} }, //NE
		{ {-3,5,5}, {-3,0,5}, {-3,-3,-3} }, //SW
		{ {5,5,-3}, {5,0,-3}, {-3,-3,-3} } //SE
	};

	std::vector<cv::Mat> kernelList;
	for(int k=0;k<nkernels;k++){
		kernelList.push_back( cv::Mat::zeros(kernSize,kernSize,CV_64FC1) );
		for (int i=0; i<kernelList[k].rows;i++){
  		for (int j=0; j<kernelList[k].cols;j++) {
				kernelList[k].at<double>(i,j)= kernValues[k][i][j];
			}
		}	
	}//end loop kernel list

	return kernelList;

}//close BuildKernels()


}//close namespace 

