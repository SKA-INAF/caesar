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
* @file GradientFilter.cc
* @class GradientFilter
* @brief Class implementing gradient filtering
*
* Gradient Filter
* @author S. Riggi
* @date 20/01/2015
*/


#include <GradientFilter.h>
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

ClassImp(Caesar::GradientFilter)

namespace Caesar {

GradientFilter::GradientFilter() {

}//close costructor


GradientFilter::~GradientFilter(){

}//close destructor

std::vector<cv::Mat> GradientFilter::BuildGradientKernels(){

	//## Init gradient kernel values (2 kernels in X & Y directions)	
	const int nkernels= 2;
	const int kernSize= 3;
	double kernValues[nkernels][kernSize][kernSize]= {
		{ {-3,0,+3}, {-10,0,+10}, {-3,0,+3} }, //X
		{ {-3,-10,-3}, {0,0,0}, {+3,+10,+3} } //Y
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

}//close BuildGradientKernels()


cv::Mat GradientFilter::BuildLaplaceKernel(){

	const int kernSize= 3;
	double kernValues[kernSize][kernSize]= { {1, 1, 1}, {1,-8, 1}, {1, 1, 1} };

	cv::Mat kernel(kernSize,kernSize,CV_64F);
	int nrows = kernel.rows;
	int ncols = kernel.cols;
	
	for (int i=0; i<nrows;i++){
  	for (int j=0; j<ncols;j++) { 
			kernel.at<double>(i,j)= kernValues[i][j];
		}
	}	
	return kernel;
	
}//close BuildLaplaceKernel()




//===================================
//==        NEW IMAGE METHODS
//===================================
Image* GradientFilter::GetGradientFilter(Image* image)
{
	//## Check image
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return 0;
	}

	//## Init kernels
	std::vector<cv::Mat> kernelList= BuildGradientKernels();
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");
	int nrows= I.rows;
	int ncols= I.cols;

	//## Compute convolutions for all kernel directions
	std::vector<cv::Mat> filteredEnsemble;
	for(size_t k=0;k<kernelList.size();k++){
		filteredEnsemble.push_back( cv::Mat::zeros(nrows,ncols,CV_64FC1) );
		filteredEnsemble[k]= Caesar::MathUtils::GetConvolution(I,kernelList[k]);	
	}//end loop kernel directions

	//## Compute gradient as sqrt(gradX^2+gradY^2)
	cv::Mat filteredMat= cv::Mat::zeros(nrows,ncols,CV_64FC1);

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for collapse(2)
	#endif
	for (int i=0; i<nrows;i++){
  	for (int j=0; j<ncols;j++) { 
			double sum2= 0;
			for(unsigned int k=0;k<filteredEnsemble.size();k++){
				double w= filteredEnsemble[k].at<double>(i,j);
				double w2= w*w;
				sum2+= w2;
			}//end loop kernels
			double val= sqrt(sum2);
			filteredMat.at<double>(i,j)= val;
		}//end loop cols	
	}//end loop rows

	//## Convert OpenCV mat list to Img
	TString imgName= Form("%s_Grad",image->GetName().c_str());
	Image* filteredImg= image->GetCloned(std::string(imgName),true,true);
	filteredImg->Reset();
	filteredImg->FillFromMat(filteredMat);

	return filteredImg;

}//close GetGradientFilter()

Image* GradientFilter::GetLaplaceFilter(Image* image)
{
	
	//## Check image
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return nullptr;
	}

	//## Init kernels
	cv::Mat kernel= BuildLaplaceKernel();
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");
	
	//## Compute convolution
	cv::Mat filteredMat= Caesar::MathUtils::GetConvolution(I,kernel);	

	//## Convert OpenCV mat list to Img
	TString imgName= Form("%s_Lapl",image->GetName().c_str());
	Image* filteredImg= image->GetCloned(std::string(imgName),true,true);
	if(!filteredImg){
		ERROR_LOG("Failed to clone image!");
		return nullptr;
	}

	//## Fill filtered image with laplacian Mat
	filteredImg->Reset();
	filteredImg->FillFromMat(filteredMat);

	return filteredImg;

}//close GetLaplaceFilter()

}//close namespace 

