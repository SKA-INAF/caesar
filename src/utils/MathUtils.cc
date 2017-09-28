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
* @file MathUtils.cc
* @class MathUtils
* @brief Utility functions for math tasks
*
* Utility functions for math tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <MathUtils.h>
#include <CodeUtils.h>
#include <Logger.h>

#include <TMath.h>

#include <linterp.h>

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

ClassImp(Caesar::MathUtils)

namespace Caesar {

MathUtils::MathUtils(){

}

MathUtils::~MathUtils(){

}

int MathUtils::Compute2DGrid(std::vector<long int>& ix_min,std::vector<long int>& ix_max,std::vector<long int>& iy_min,std::vector<long int>& iy_max,long int Nx,long int Ny,long int boxSizeX,long int boxSizeY,float gridStepSizeX,float gridStepSizeY){

	//## Check given arguments
	if(Nx<=0 || Ny<=0){
		ERROR_LOG("Invalid Nx/Ny given (negative or zero)!");
		return -1;
	}
	if(boxSizeX<=0 || boxSizeY<=0) {
		ERROR_LOG("Invalid box size given!");
		return -1;
	}
	if(gridStepSizeX<=0 || gridStepSizeY<=0 || gridStepSizeX>1 || gridStepSizeY>1){
		ERROR_LOG("Invalid grid step size given (null or negative)!");
		return -1;
	}

	//## Check if image size is smaller than required box
	if(boxSizeX>Nx || boxSizeY>Ny) {
		WARN_LOG("Invalid box size given (too small or larger than image size)!");
		return -1;
	}

	long int stepSizeX= std::round(gridStepSizeX*boxSizeX);
	long int stepSizeY= std::round(gridStepSizeY*boxSizeY);
	long int indexX= 0;
	long int indexY= 0;
	ix_min.clear();
	ix_max.clear();
	iy_min.clear();
	iy_max.clear();
	

	while(indexY<=Ny){
		long int offsetY= min(boxSizeY-1,Ny-1-indexY);
		long int ymin= indexY;
		long int ymax= indexY+offsetY;
		if(ymin>=Ny || offsetY==0) break;	
		iy_min.push_back(ymin);
		iy_max.push_back(ymax);
		indexY+= stepSizeY;
	}//end while loop Y
		
	while(indexX<=Nx){
		long int offsetX= min(boxSizeX-1,Nx-1-indexX);
		long int xmin= indexX;
		long int xmax= indexX+offsetX;
		if(xmin>=Nx || offsetX==0) break;	
		ix_min.push_back(xmin);
		ix_max.push_back(xmax);
		indexX+= stepSizeX;
	}//end while loop Y

	
	return 0;

}//close Compute2DGrid()


int MathUtils::BiLinearInterpolation(std::vector<double>const& sampled_gridX, std::vector<double>const& sampled_gridY,std::vector<double>const& sampledZ,std::vector<double>const& interp_gridX,std::vector<double>const& interp_gridY,std::vector<double>& interpZ){

	//## Check args
	long int nSamplesX= (long int)sampled_gridX.size();
	long int nSamplesY= (long int)sampled_gridY.size();
	long int nSamplesZ= (long int)sampledZ.size();
	if(nSamplesX<=0 || nSamplesY<=0){
		ERROR_LOG("Invalid sample grid size given!");
		return -1;
	}
	long int num_elements = nSamplesX*nSamplesY;
	if(nSamplesZ!=num_elements){
		ERROR_LOG("Invalid sample Z given (it must be equal to Nx x Ny and indexed as ix*Ny+iy)!");
		return -1;
	}

	long int nInterpX= (long int)interp_gridX.size();
	long int nInterpY= (long int)interp_gridY.size();
  long int num_interp_elements = nInterpX*nInterpY;
	if(nInterpX<=0 || nInterpY<=0){
		ERROR_LOG("Invalid interpolation grid given (size must be >0 in both directions)!");
		return -1;
	}
	interpZ.clear();

	//## Perform the 2D interpolation
	try {		
	
		// Construct the grid in each dimension (note that we will pass in a sequence of iterators pointing to the beginning of each grid)
		DEBUG_LOG("Build 2D grid for interpolation (nSamplesX="<<nSamplesX<<", nSamplesY="<<nSamplesY<<")...");

  	std::vector< std::vector<double>::const_iterator > grid_iter_list;
  	grid_iter_list.push_back(sampled_gridX.begin());
  	grid_iter_list.push_back(sampled_gridY.begin());
  
  	// the size of the grid in each dimension
  	array<int,2> grid_sizes;
  	grid_sizes[0] = nSamplesX;
  	grid_sizes[1] = nSamplesY;
  
  	// construct the interpolator. the last two arguments are pointers to the underlying data
		DEBUG_LOG("Build the bkg interpolator...");
  	InterpMultilinear<2, double> interpolator_ML(grid_iter_list.begin(), grid_sizes.begin(), sampledZ.data(), sampledZ.data() + num_elements);
		
		
		//Construct interpolated grid
		// interpolate multiple values: create sequences for each coordinate
		DEBUG_LOG("Build the interpolated grid ("<<nInterpX<<","<<nInterpY<<") x("<<interp_gridX[0]<<","<<interp_gridX[nInterpX-1]<<") y("<<interp_gridY[0]<<","<<interp_gridY[nInterpY-1]<<")...");
  	
  	//std::vector<double> interp_gridX = CodeUtils::linspace(xlim[0],xlim[1], Nx);
		//std::vector<double> interp_gridY = CodeUtils::linspace(ylim[0],ylim[1], Ny);
		
  	std::vector<double> interp_x(num_interp_elements);
  	std::vector<double> interp_y(num_interp_elements);
  	for (unsigned int i=0; i<interp_gridX.size(); i++) {
    	for (unsigned int j=0; j<interp_gridY.size(); j++) {
				long int gBin= i*interp_gridY.size() + j;
	  		interp_x[gBin] = interp_gridX[i];
	  		interp_y[gBin] = interp_gridY[j];
			}
  	}
 	 	interpZ.assign(num_interp_elements,0);

		// pass in a sequence of iterators, one for each coordinate
  	std::vector< std::vector<double>::iterator > interp_list;
  	interp_list.push_back(interp_x.begin());
  	interp_list.push_back(interp_y.begin());
  
		//Interpolate sequence
		DEBUG_LOG("Run the interpolation on grid...");
		interpolator_ML.interp_vec(num_interp_elements, interp_list.begin(), interp_list.end(), interpZ.begin());
  	
	}//close try block
	catch( std::exception &ex ) {
		ERROR_LOG("Exception detected in interpolation (err=" << ex.what()<<")");
		return -1;
  } 
	catch(...) { 
		ERROR_LOG("Unknown exception caught in interpolation!");
		return -1;
  }		
	
	return 0;

}//close BiLinearInterpolation()


cv::Mat MathUtils::GetConvolution(cv::Mat I, cv::Mat kernel){

	// anchor Point
 	cv::Point anchor(kernel.cols - kernel.cols/2 - 1, kernel.rows - kernel.rows/2 - 1);

 	// flip the kernel
 	cv::Mat kernel_flipped;
 	int flip_code = -1; // -1-both axis, 0-x-axis, 1-y-axis
 	cv::flip(kernel, kernel_flipped, flip_code);

	//Apply convolution
	cv::Mat dst;
	filter2D(I, dst, I.depth(), kernel_flipped, anchor, cv::BORDER_CONSTANT);

	return dst;

}//close GetConvolution()

cv::Mat MathUtils::GetConvolution2(cv::Mat I, cv::Mat kernel){

	int nrows= I.rows;
	int ncols= I.cols;
	int nrows_kernel= kernel.rows;
	int ncols_kernel= kernel.cols;
	cv::Mat dst= cv::Mat::zeros(nrows,ncols,CV_64FC1);

	for(int i=0;i<nrows;i++){//loop rows
		for(int j=0;j<ncols;j++){//loop cols

			for(int k=0;k<nrows_kernel;k++){//start loops on filter box rows
				int index_x= i-k+nrows_kernel/2;
				int mirror_index_x= GetMirrorIndex(index_x,nrows);

				for(int l=0;l<ncols_kernel;l++){//start loops on filter box cols
					int index_y= j-l+ncols_kernel/2;
					int mirror_index_y= GetMirrorIndex(index_y,ncols);	

					double K= kernel.at<double>(k,l);
					double f= I.at<double>(mirror_index_x,mirror_index_y);	
					dst.at<double>(i,j)= K*f;
				}//end loop kernel cols
			}//end loop kernel rows
		}//end loop y
	}//end loop x

	return dst;
		
}//close GetConvolution2()

cv::Mat MathUtils::GetATrousConvolution(cv::Mat I, cv::Mat kernel,int scale){

	int nrows= I.rows;
	int ncols= I.cols;
	int nrows_kernel= kernel.rows;
	int ncols_kernel= kernel.cols;
	int rowIndex= (nrows_kernel-1)/2;
	int colIndex= (ncols_kernel-1)/2;
	int filterGap= pow(2,scale-1); 
	cv::Mat dst= cv::Mat::zeros(nrows,ncols,CV_64FC1);
			
	//Compute convolution
	for(int i=0;i<nrows;i++){
		for(int j=0;j<ncols;j++){

			for(int k=-rowIndex;k<rowIndex;k++){//start loops on filter box rows
				int index_x= i + filterGap*k;
				int mirror_index_x= GetMirrorIndex(index_x,nrows);
	
				for(int l=-colIndex;l<colIndex;l++){//start loops on filter box cols
					int index_y= j + filterGap*l;
					int mirror_index_y= GetMirrorIndex(index_y,ncols);
		
					double K= kernel.at<double>(k+rowIndex,l+colIndex);
					double f= I.at<double>(mirror_index_x,mirror_index_y); 
					dst.at<double>(i,j)+= K*f;
				}//end loop filter cols
			}//end loop filter rows
		}//end loop cols
	}//end loop rows

	return dst;

}//close GetATrousConvolution()


int MathUtils::GetMirrorIndex(int index,int N){
			
	int mirror_index= 0;
	if(index>=0 && index<=(N-1)) {
		mirror_index= index;
	}
	else if(index<0){
		mirror_index= -index;
	}
	else if(index>(N-1)){
		mirror_index= 2*(N-1)-index;
	}
	else{
		WARN_LOG("Invalid index of matrix size passed...exit!");
		return -1;
	}	
	return mirror_index;

}//close GetMirrorIndex()


std::vector< std::complex<double> > MathUtils::DFTShifted(std::vector< std::complex<double> > data, int n){

	int N= (int)data.size();
	if(n>N || n<0) n= N;//check size

	std::vector< std::complex<double> > Fn(n,0);
	for(int i=0; i<n; i++) {//loop over n
		Fn[i] = std::complex<double>(0.,0.);
		int s= -floor(N/2.) + i;
		//int s= i;
	
		for(int j=0;j<N;j++) {//loop over data size
			int k= j;
			//int k= -N/2 + j;
			double arg= 2.*TMath::Pi()*s*k/N;
			std::complex<double> prod= std::polar(1.,-arg);
			Fn[i]+= data[j]*prod;
		}//end loop data size
		Fn[i]/= (double)N;
	}//end loop n

	return Fn;

}//close DFTShifted()

std::vector< std::complex<double> > MathUtils::DFT(std::vector< std::complex<double> > data, int n){
			
	int N= (int)data.size();
	if(n>N || n<0) n= N;//check size

	std::vector< std::complex<double> > Fn(n,0);
	for(int i=0; i<n; i++) {//loop over n
		Fn[i] = std::complex<double>(0.,0.);
		//int s= -floor(N/2.) + i;
		int s= i;

		for(int j=0;j<N;j++) {//loop over data size
			int k= j;
			//int k= -floor(N/2.) + j;
			//double arg= 2.*TMath::Pi()*i*k/N;
			double arg= 2.*TMath::Pi()*s*k/N;
			std::complex<double> prod= std::polar(1.,-arg);
			Fn[i]+= data[j]*prod;
		}//end loop data size
	}//end loop n

	return Fn;

}//close DFT()


std::vector< std::complex<double> > MathUtils::IDFT(std::vector< std::complex<double> > data, int n){
			
	int N= (int)data.size();
	if(n>N || n<0) n= N;//check size

	std::vector< std::complex<double> > fn(n,0);
	for(int i=0; i<n; i++) {//loop over n
		fn[i] = std::complex<double>(0.,0.);
		int s= i;

		for(int j=0;j<N;j++) {//loop over data size
			int k= j;
			//int k= -floor(N/2.) + j;
				
			double arg= 2.*TMath::Pi()*s*k/N;
			std::complex<double> prod= std::polar(1.,arg);
			fn[i]+= data[j]*prod;
		}//end loop data size
		fn[i]/= (double)N;	
	}//end loop n

	return fn;

}//close IDFT()


std::vector<double> MathUtils::GetContourCurvature(std::vector< std::complex<double> > data,double sigma){
			
	int N= (int)data.size();
	std::vector< std::complex<double> > Us= DFT(data,N);

	//Compute Us smoothed with a gaussian
	std::vector< std::complex<double> > Us_firstDeriv; 
	std::vector< std::complex<double> > Us_secondDeriv; 
	std::vector< std::complex<double> > Us_firstDeriv_smoothed;
	std::vector< std::complex<double> > Us_secondDeriv_smoothed;
	for(int i=0;i<N;i++){
		//int s= -floor(N/2.) + i;
		//double Gs= sigma/sqrt(2*TMath::Pi())*exp(-pow(0.5*sigma*s,2));
		//double arg= 2*TMath::Pi()*s;
		int eta= EtaAuxiliaryFcn(i,N);
		double Gs= exp(-pow(sigma*eta,2));
		double arg= 2*TMath::Pi()*eta;
		std::complex<double> z(0,arg); 
		std::complex<double> thisUs_firstDeriv= z*Us[i];	
		std::complex<double> thisUs_secondDeriv= -pow(arg,2)*Us[i];
		std::complex<double> thisUs_firstDeriv_smoothed= thisUs_firstDeriv*Gs;
		std::complex<double> thisUs_secondDeriv_smoothed= thisUs_secondDeriv*Gs;
		Us_firstDeriv.push_back(thisUs_firstDeriv);
		Us_secondDeriv.push_back(thisUs_secondDeriv);
		Us_firstDeriv_smoothed.push_back(thisUs_firstDeriv_smoothed);
		Us_secondDeriv_smoothed.push_back(thisUs_secondDeriv_smoothed);
				
	}//end loop points

	//Compute ut'= IDST(Us')=IDST(i x 2pi x s x Us)	
	std::vector< std::complex<double> > ut_firstDeriv= IDFT(Us_firstDeriv,N);
	double L= 0;
	for(unsigned int i=0;i<ut_firstDeriv.size();i++) L+= std::abs(ut_firstDeriv[i]);	
	L*= 2.*TMath::Pi()/N;
		
	//Compute inverse transform of smoothed
	std::vector< std::complex<double> > ut_firstDeriv_smoothed= IDFT(Us_firstDeriv_smoothed, N);
	std::vector< std::complex<double> > ut_secondDeriv_smoothed= IDFT(Us_secondDeriv_smoothed, N);

	//Compute smoothed perymeter
	double L_smoothed= 0;
	for(int i=0;i<N;i++){
		std::complex<double> u1= ut_firstDeriv_smoothed[i];
		double u1_mod= std::abs(u1);  	
		L_smoothed+= u1_mod;
	}//end loop points
	L_smoothed*= 2.*TMath::Pi()/N;
		
	//Compute curvature
	double normFactor= L/L_smoothed;
	std::vector<double> curvature;
	for(int i=0;i<N;i++){
		//Normalize ut
		ut_firstDeriv_smoothed[i]*= normFactor;
		ut_secondDeriv_smoothed[i]*= normFactor;

		std::complex<double> u2_conj= std::conj(ut_secondDeriv_smoothed[i]); 
		std::complex<double> u1= ut_firstDeriv_smoothed[i];
		double u1_mod= std::abs(u1);  	
		double curv= -std::imag(u1*u2_conj)/pow(u1_mod,3);
		//curv*= L;//multiply by perymeter length
		//cout<<"Pnt no. "<<i<<" usmooth'="<<ut_firstDeriv_smoothed[i]<<" usmooth''="<<ut_secondDeriv_smoothed[i]<<" u1_mod="<<u1_mod<<" curv="<<curv<<endl;
		curvature.push_back(curv);
	}//end loop points

	return curvature;

}//close GetContourCurvature()


}//close namespace



