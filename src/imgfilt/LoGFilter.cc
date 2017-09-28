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
* @file LoGFilter.cc
* @class LoGFilter
* @brief Class implementing Logarithm of Gaussian (LoG) filtering
*
* Logarithm of Gaussian (LoG) Filter
* @author S. Riggi
* @date 20/01/2015
*/

#include <LoGFilter.h>
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

ClassImp(Caesar::LoGFilter)

namespace Caesar {

LoGFilter::LoGFilter() {

}//close costructor


LoGFilter::~LoGFilter(){

}//close destructor


//===================================
//==        NEW IMAGE METHODS
//===================================
Image* LoGFilter::GetLoGFilter(Image* image){

	//## Check image
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return 0;
	}

	//## Init kernels
	cv::Mat kernel= BuildStandardKernel();
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");

	//## Compute convolution
	cv::Mat filteredMat= Caesar::MathUtils::GetConvolution(I,kernel);	

	//## Convert OpenCV mat list to Img
	TString imgName= Form("%s_NormLoG",image->GetName().c_str());
	Image* filteredImg= image->GetCloned(std::string(imgName),true,true);
	filteredImg->Reset();
	filteredImg->FillFromMat(filteredMat);

	return filteredImg;

}//close GetLoGFilter()


Image* LoGFilter::GetNormLoGFilter(Image* image,int size,double scale){
	
	//## Check image
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return 0;
	}

	//## Init kernels
	cv::Mat kernel= BuildKernel(size,scale);
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");
	
	//## Compute convolution
	cv::Mat filteredMat= Caesar::MathUtils::GetConvolution(I,kernel);	

	//## Convert OpenCV mat list to Img
	TString imgName= Form("%s_NormLoG",image->GetName().c_str());
	Image* filteredImg= image->GetCloned(std::string(imgName),true,true);
	filteredImg->Reset();
	filteredImg->FillFromMat(filteredMat);

	return filteredImg;

}//close GetNormLoGFilter()



cv::Mat LoGFilter::BuildStandardKernel(){

	const int kernSize= 7;
	double kernValues[kernSize][kernSize]= { 
		{0,0,1,1,1,0,0}, 
		{0,1,3,3,3,1,0}, 
		{1,3,0,-7,0,3,1}, 
		{1,3,-7,-24,-7,3,1}, 
		{0,0,1,1,1,0,0}, 
		{0,1,3,3,3,1,0}, 
		{1,3,0,-7,0,3,1} 
	};
	
	cv::Mat kernel(kernSize,kernSize,CV_64F);
	int rows = kernel.rows;
	int cols = kernel.cols;

	for (int i=0; i<rows;i++){
  	for (int j=0; j<cols;j++) { 
			kernel.at<double>(i,j)= kernValues[i][j];
		}
	}	
	return kernel;

}//close BuildStandardKernel()


cv::Mat LoGFilter::BuildKernel(int kernSize,double scale){

	cv::Mat kernel(kernSize,kernSize,CV_64F);
	int nrows = kernel.rows;
	int ncols = kernel.cols;
	int halfSize= (kernSize-1)/2;
		
	for (int i=0;i<nrows;i++){
  	for (int j=0;j<ncols;j++){
  		double x = j - halfSize;
  	  double y = i - halfSize;
			double kernValue= NormLoGKernel(x,y,scale);
			kernel.at<double>(i,j)= kernValue;
		}//end loop j
	}//end loop i
	
	cv::Scalar kernelMean= cv::mean(kernel);
	cv::subtract(kernel,kernelMean,kernel);

	return kernel;

}//close BuildKernel()


double LoGFilter::NormLoGKernel(double x,double y,double sigma){

	double arg= (x*x+y*y)/(2.*sigma*sigma);
	double norm= 1./(TMath::Pi()*pow(sigma,4));
	double fcn= -pow(sigma,2)*norm*(1-arg)*exp(-arg);

	return fcn;

}//close LoGFilter::NormLoGKernel()


}//close namespace 








