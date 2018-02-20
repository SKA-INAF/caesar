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
* @file GausFilter.cc
* @class GausFilter
* @brief Class implementing elliptical gaussian filtering
*
* Elliptical gaussian Filter
* @author S. Riggi
* @date 20/01/2015
*/

#include <GausFilter.h>
#include <Image.h>
#include <MathUtils.h>
#include <Consts.h>

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

ClassImp(Caesar::GausFilter)

namespace Caesar {

GausFilter::GausFilter() {

}//close costructor


GausFilter::~GausFilter(){

}//close destructor



Image* GausFilter::GetGausFilter(Image* image,double bmaj,double bmin,double bpa,int nSigmas,double scale){
	
	//## Check image
	if(!image){
		ERROR_LOG("Null prt to given image!");
		return nullptr;
	}

	//## Check if image has metadata and get 
	ImgMetaData* metadata= image->GetMetaData();
	if(!metadata){
		ERROR_LOG("Input image has no metadata (needed to convert given bmaj/bmin/bpa into pixels)!");
		return nullptr;
	}
	double dX= fabs(metadata->dX*3600);//convert to arcsec
	double dY= fabs(metadata->dY*3600);//convert to arcsec

	//## Check bmaj/bmin
	if(bmin>bmaj){
		ERROR_LOG("Invalid bmaj/bmin given (hint: bmaj shall be larger or equal to bmin)!");
		return nullptr;
	}
	
	//## Convert bmaj/bin to pixels 
	//## NB: The following is wrong, assume rectangular pixels (need sky conversions)!!!
	double pixSize= max(dX,dY);
	double bmaj_pix= bmaj/pixSize;
	double bmin_pix= bmin/pixSize;

	//## Compute sigmaX/sigmaY from bmaj/bmin
	double sigmaX= bmaj_pix/GausSigma2FWHM;
	double sigmaY= bmin_pix/GausSigma2FWHM;
	double theta= bpa + 90.;//# NB: BPA is the positional angle of the major axis measuring from North (up) counter clockwise, while theta is measured wrt to x axis
	double theta_rad= theta*TMath::DegToRad();//convert to radians

	//## Determine kernel size optimal. Assume nsigma volume (default 5 sigma volume)
	int kernSizeX= (static_cast<int>(nSigmas*sigmaX + 0.5) + 1) * 2;
	int kernSizeY= (static_cast<int>(nSigmas*sigmaY + 0.5) + 1) * 2;
	int kernSize= max(kernSizeX,kernSizeY);
	if (kernSize%2==0) kernSize++;
	DEBUG_LOG("kernSize="<<kernSize<<" pixSize="<<pixSize<<", sigmaX="<<sigmaX<<", sigmaY="<<sigmaY);

	//## Init kernels
	cv::Mat kernel= BuildKernel(kernSize,sigmaX,sigmaY,theta_rad,scale);
	
	//## Convert input image to OpenCV mat
	cv::Mat I= image->GetOpenCVMat("64");
	
	//## Compute convolution
	cv::Mat filteredMat= Caesar::MathUtils::GetConvolution(I,kernel);	

	//## Convert OpenCV mat list to Img
	Image* filteredImg= image->GetCloned("",true,true);
	filteredImg->Reset();
	filteredImg->FillFromMat(filteredMat);

	return filteredImg;

}//close GetGausFilter()


cv::Mat GausFilter::BuildKernel(int kernSize,double sigmaX,double sigmaY,double theta,double scale){

	cv::Mat kernel(kernSize,kernSize,CV_64F);
	int nrows = kernel.rows;
	int ncols = kernel.cols;
	int halfSize= (kernSize-1)/2;
	double ampl= scale;
	double maxValue= -1.e+99;

	//Find kernel maximum
	for (int i=0;i<nrows;i++){
  	for (int j=0;j<ncols;j++){
  		double x = j - halfSize;
  	  double y = i - halfSize;
			double kernValue= Gaus2DFcn(ampl,x,y,sigmaX,sigmaY,theta);
			if(kernValue>maxValue) maxValue= kernValue;
			kernel.at<double>(i,j)= kernValue;
		}//end loop j
	}//end loop i

	//Normalize kernel by maximum value
	for (int i=0;i<nrows;i++){
  	for (int j=0;j<ncols;j++){
  		double x = j - halfSize;
  	  double y = i - halfSize;
			double kernValue= Gaus2DFcn(ampl,x,y,sigmaX,sigmaY,theta);
			kernel.at<double>(i,j)= kernValue/maxValue;
		}//end loop j
	}//end loop i
	
	return kernel;

}//close BuildKernel()


double GausFilter::Gaus2DFcn(double ampl,double x,double y,double sigmaX,double sigmaY,double theta){

	double x_mean= 0;
	double y_mean= 0;
	double cost2 = cos(theta)*cos(theta);
	double sint2 = sin(theta)*sin(theta);
	double sin2t = sin(2. * theta);
	double xstd2 = sigmaX*sigmaX;
	double ystd2 = sigmaY*sigmaY;
	double xdiff = x - x_mean;
  double ydiff = y - y_mean;
	double a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2));
  double b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2));
 	double c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2));

	double fcn= ampl * exp(-((a * xdiff * xdiff) + (b * xdiff * ydiff) + (c * ydiff * ydiff)));
	
	return fcn;

}//close Gaus2DFcn()


}//close namespace 


