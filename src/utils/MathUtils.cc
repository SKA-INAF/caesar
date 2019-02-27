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
#include <EllipseUtils.h>
#include <Contour.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//ROOT headers
#include <TMath.h>
#include <TEllipse.h>
#include <TGraph.h>

//BOOST headers
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/scoped_array.hpp>

//linterp header
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
#include <deque>

using namespace std;

ClassImp(Caesar::MathUtils)

namespace Caesar {

MathUtils::MathUtils(){

}

MathUtils::~MathUtils(){

}

int MathUtils::Compute2DGrid(std::vector<long int>& ix_min,std::vector<long int>& ix_max,std::vector<long int>& iy_min,std::vector<long int>& iy_max,long int Nx,long int Ny,long int boxSizeX,long int boxSizeY,float gridStepSizeX,float gridStepSizeY)
{
	//## Check given arguments
	if(Nx<=0 || Ny<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid Nx/Ny given (negative or zero)!");
		#endif
		return -1;
	}
	if(boxSizeX<=0 || boxSizeY<=0) {
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Invalid box size given!");
		#endif
		return -1;
	}
	if(gridStepSizeX<=0 || gridStepSizeY<=0 || gridStepSizeX>1 || gridStepSizeY>1){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid grid step size given (null or negative)!");
		#endif
		return -1;
	}

	//## Check if image size is smaller than required box
	if(boxSizeX>Nx || boxSizeY>Ny) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid box size given (too small or larger than image size)!");		
		#endif
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


int MathUtils::FindGrid2DAxisBin(float x,long int nx,float xmin,float xmax,float xstep)
{
	if(x<xmin || x>xmax) return -1;
	long int index= static_cast<long int>(nx*(x-xmin)/(nx*xstep) );

	return index;

}//close FindGrid2DAxisBin()


long int MathUtils::FindGrid2DBin(float x,float y,long int nx,float xmin,float xmax,float xstep,long int ny,float ymin,float ymax,float ystep)
{
	//Find x index
	long int binX= FindGrid2DAxisBin(x,nx,xmin,xmax,xstep);
	if(binX<0) return -1;

	//Find y index
	long int binY= FindGrid2DAxisBin(y,ny,ymin,ymax,ystep);
	if(binY<0) return -1;

	//Find global bin
	long int binId= binX + binY*nx;

	return binId;

}//close FindTileIndex()

int MathUtils::Compute2DFloatGrid(std::vector<float>& ix_min,std::vector<float>& ix_max,std::vector<float>& iy_min,std::vector<float>& iy_max,float xmin,float xmax,float xstep,float ymin,float ymax,float ystep)
{
	long int nTilesX= static_cast<long int>(std::ceil(fabs(xmax-xmin)/xstep));
	long int nTilesY= static_cast<long int>(std::ceil(fabs(ymax-ymin)/ystep));
	float TileSizeX= xstep;
	float TileSizeY= ystep;	

	ix_min.clear();
	ix_max.clear();
	iy_min.clear();
	iy_max.clear();

	float y= ymin;
	for(int j=0;j<nTilesY;j++){
		float offset_y= fabs(std::min(ystep,ymax-y));
		float tile_ymin= y;
		float tile_ymax= y + offset_y;
		iy_min.push_back(tile_ymin);
		iy_max.push_back(tile_ymax);
		y+= offset_y;
	}//end loop tile y

	float x= xmin;
	for(int i=0;i<nTilesX;i++){
		float offset_x= fabs(std::min(xstep,xmax-x));
		float tile_xmin= x;
		float tile_xmax= x + offset_x;
		ix_min.push_back(tile_xmin);
		ix_max.push_back(tile_xmax);
		x+= offset_x;
	}//end loop tile y

	return 0;

}//close Compute2DGrid()

int MathUtils::BiLinearInterpolation(std::vector<double>const& sampled_gridX, std::vector<double>const& sampled_gridY,std::vector<double>const& sampledZ,std::vector<double>const& interp_gridX,std::vector<double>const& interp_gridY,std::vector<double>& interpZ){

	//## Check args
	long int nSamplesX= (long int)sampled_gridX.size();
	long int nSamplesY= (long int)sampled_gridY.size();
	long int nSamplesZ= (long int)sampledZ.size();
	if(nSamplesX<=0 || nSamplesY<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid sample grid size given!");
		#endif
		return -1;
	}
	long int num_elements = nSamplesX*nSamplesY;
	if(nSamplesZ!=num_elements){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid sample Z given (it must be equal to Nx x Ny and indexed as ix*Ny+iy)!");
		#endif
		return -1;
	}

	long int nInterpX= (long int)interp_gridX.size();
	long int nInterpY= (long int)interp_gridY.size();
  long int num_interp_elements = nInterpX*nInterpY;
	if(nInterpX<=0 || nInterpY<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid interpolation grid given (size must be >0 in both directions)!");
		#endif
		return -1;
	}
	interpZ.clear();

	//## Perform the 2D interpolation
	try {		
	
		// Construct the grid in each dimension (note that we will pass in a sequence of iterators pointing to the beginning of each grid)
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Build 2D grid for interpolation (nSamplesX="<<nSamplesX<<", nSamplesY="<<nSamplesY<<")...");
		#endif

  	std::vector< std::vector<double>::const_iterator > grid_iter_list;
  	grid_iter_list.push_back(sampled_gridX.begin());
  	grid_iter_list.push_back(sampled_gridY.begin());
  
  	// the size of the grid in each dimension
  	array<int,2> grid_sizes;
  	grid_sizes[0] = nSamplesX;
  	grid_sizes[1] = nSamplesY;
  
  	// construct the interpolator. the last two arguments are pointers to the underlying data
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Build the bkg interpolator...");
		#endif
  	InterpMultilinear<2, double> interpolator_ML(grid_iter_list.begin(), grid_sizes.begin(), sampledZ.data(), sampledZ.data() + num_elements);
		
		
		//Construct interpolated grid
		// interpolate multiple values: create sequences for each coordinate
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Build the interpolated grid ("<<nInterpX<<","<<nInterpY<<") x("<<interp_gridX[0]<<","<<interp_gridX[nInterpX-1]<<") y("<<interp_gridY[0]<<","<<interp_gridY[nInterpY-1]<<")...");
  	#endif

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
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Run the interpolation on grid...");
		#endif
		interpolator_ML.interp_vec(num_interp_elements, interp_list.begin(), interp_list.end(), interpZ.begin());
  	
	}//close try block
	catch( std::exception &ex ) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Exception detected in interpolation (err=" << ex.what()<<")");
		#endif
		return -1;
  } 
	catch(...) { 
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Unknown exception caught in interpolation!");
		#endif
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
	//filter2D(I, dst, I.depth(), kernel_flipped, anchor, cv::BORDER_REPLICATE);
	filter2D(I, dst, I.depth(), kernel_flipped, anchor, cv::BORDER_CONSTANT);
	//filter2D(I, dst, I.depth(), kernel_flipped, anchor, cv::BORDER_DEFAULT);

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

			double w= 0;
			for(int k=0;k<nrows_kernel;k++){//start loops on filter box rows
				int index_x= i-k+nrows_kernel/2;
				int mirror_index_x= GetMirrorIndex(index_x,nrows);

				for(int l=0;l<ncols_kernel;l++){//start loops on filter box cols
					int index_y= j-l+ncols_kernel/2;
					int mirror_index_y= GetMirrorIndex(index_y,ncols);	

					double K= kernel.at<double>(k,l);
					double f= I.at<double>(mirror_index_x,mirror_index_y);	
					//dst.at<double>(i,j)= K*f;
					if(std::isnormal(f)) w+= K*f;
				}//end loop kernel cols
			}//end loop kernel rows

			dst.at<double>(i,j)= w;

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
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid index of matrix size passed...exit!");
		#endif
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


int MathUtils::ComputeEllipseOverlapArea(double& overlapArea,double& err,int& rtn,TEllipse* ellipse1, TEllipse* ellipse2,int method,Contour* overlapContour)
{
	//Init area
	overlapArea= 0;

	//Check inputs
	if(!ellipse1 || !ellipse2) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input ellipses given!");
		#endif
		return -1;
	}
	
	//Get ellipse params
	double Theta_1= ellipse1->GetTheta();
	double R1_1= ellipse1->GetR1();
	double R2_1= ellipse1->GetR2(); 
	double Cx_1= ellipse1->GetX1();
	double Cy_1= ellipse1->GetY1();
	double Area_1= TMath::Pi()*R1_1*R2_1;

	double Theta_2= ellipse2->GetTheta();
	double R1_2= ellipse2->GetR1();
	double R2_2= ellipse2->GetR2(); 
	double Cx_2= ellipse2->GetX1();
	double Cy_2= ellipse2->GetY1(); 
	double Area_2= TMath::Pi()*R1_2*R2_2;

	double x[4], y[4];
  int nroots= 0;
	
	//Compute precise ellipse overlap area
	overlapArea= EllipseUtils::ellipse_ellipse_overlap (
			Theta_1, R1_1, R2_1, Cx_1, Cy_1, 
			Theta_2, R1_2, R2_2, Cx_2, Cy_2,
			x,y,&nroots, &rtn, method
	);
	if(overlapArea<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Ellipse overlap computation failed with status "<<rtn<<"!");	
		#endif
		return -1;
	}

	//Convert ellipse in polygons and find polygon overlapping areas
	polygon_2d poly, poly2;
	if(Ellipse2Polygon(poly,Cx_1, Cy_1,R1_1,R2_1,Theta_1)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to convert first ellipse to polygon!");
		#endif
		return -1;
	}
	if(Ellipse2Polygon(poly2,Cx_2,Cy_2,R1_2,R2_2,Theta_2)<0){	
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to convert second ellipse to polygon!");
		#endif
		return -1;
	}

	//Compute overlapping area between the two polygons
	double overlapArea_poly = 0;
	polygon_2d overlap_poly;
	if(ComputePolygonOverlapArea(overlapArea_poly,overlap_poly,poly,poly2)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute overlapping area between the two polygonized ellipses!");
		#endif
		return -1;
	}

	//Estimate errors between precise and approximate calculations
	double eps = 0.0001;
	err = (fabs(overlapArea_poly)<eps)? 0.0 : fabs(overlapArea_poly-overlapArea)/overlapArea_poly;

	//If graph given fill it with overlap area points
	if(overlapContour){
		//overlapAreaGraph->Set(0);
		for(int i=0; i<nroots; i++) {
			double theta= -Theta_1*TMath::DegToRad();	
			double Xroot= x[i]*cos(theta) + y[i]*sin(theta) + Cx_1;
			double Yroot= -x[i]*sin(theta) + y[i]*cos(theta) + Cy_1;
			//overlapAreaGraph->SetPoint(i,Xroot,Yroot);
			overlapContour->AddPoint(TVector2(Xroot,Yroot));
		}
		/*
  	overlapAreaGraph->SetMarkerColor(kBlue-6);
		overlapAreaGraph->SetMarkerStyle(24);
		overlapAreaGraph->SetMarkerSize(1.3);
		overlapAreaGraph->SetFillColor(kBlue-6);
		overlapAreaGraph->SetLineColor(kBlue-6);
		*/
		//Sort contour points counter clockwise
		overlapContour->SortPointsCounterClockWise();
	
	}//close if

	return 0;

}//close ComputeEllipseOverlapArea()

double MathUtils::ComputePolygonArea(polygon_2d& poly)
{
	double area = boost::geometry::area(poly);
	return area;

}//close ComputePolygonArea()

int MathUtils::ComputeContourArea(double& area,Contour* contour)
{
	//Init area
	area= -1;

	//Check contour
	if(!contour){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to contour given!");
		#endif
		return -1;
	}
	
	//Transform contour in polygon
	polygon_2d poly2d;
	if(Contour2Polygon(poly2d,contour)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to convert contour to polygon!");
		#endif
		return -1;
	}

	//Compute polygon area
	area= ComputePolygonArea(poly2d);

	return 0;

}//close ComputePolygonArea()


int MathUtils::ComputeContourOverlapArea(double& overlapArea,int& overlapFlag,Contour* contour,Contour* contour2,Contour* overlapContour)
{
	//Init area 
	overlapArea= -1;
	overlapFlag= eCONT_NOT_OVERLAPPING;
	
	//Convert contours to polygons
	polygon_2d poly, poly2;
	if(Contour2Polygon(poly,contour)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to convert first contour to polygon!");
		#endif
		return -1;
	}
	if(Contour2Polygon(poly2,contour2)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to convert second contour to polygon!");
		#endif
		return -1;
	}

	//Compute overlap area between polygons
	polygon_2d overlap_poly;
	if(ComputePolygonOverlapArea(overlapArea,overlap_poly,poly,poly2)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute polygon overlap area!");
		#endif
		overlapFlag= eCONT_OVERLAP_FAILED;
		return -1;
	}

	//Compute overlap flags
	if(overlapArea>0){
		//Set overlapping flag
		overlapFlag= eCONT_OVERLAPPING;

		//Check if contour are one inside the other
		bool isFirstInsideSecond= boost::geometry::within(poly, poly2);
		bool isSecondInsideFirst= false;
		if(isFirstInsideSecond) {
			overlapFlag= eCONT1_INSIDE_CONT2;
		}
		else{
			isSecondInsideFirst= boost::geometry::within(poly2, poly);
			if(isSecondInsideFirst) overlapFlag= eCONT2_INSIDE_CONT1;
		}
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Contour 1 inside 2? "<<isFirstInsideSecond<<", Contour 2 inside 1? "<<isSecondInsideFirst);
		#endif
	}
	else{
		//Set non-overlapping contours
		overlapFlag= eCONT_NOT_OVERLAPPING;
	}

	
	//If overlap contour given fill it
	if(overlapContour && overlapFlag!=eCONT_NOT_OVERLAPPING){
		std::vector<point_xy> const& points = overlap_poly.outer();
		std::vector<TVector2> points_noduplicates;
		int pos= -1;
		for (auto point : points) {
    	double x= point.x();
			double y= point.y();
			TVector2 pnt(x,y);
			std::vector<TVector2>::iterator it= std::find_if(points_noduplicates.begin(),points_noduplicates.end(),ContourPointComparator(pnt));
			if(it==points_noduplicates.end() || points_noduplicates.empty()){
				overlapContour->AddPoint(TVector2(x,y));
				points_noduplicates.push_back(pnt);
			}
    }//end loop polygon points
		
		//Compute contour pars (need centroid)
		int status= overlapContour->ComputeParameters();
		if(status<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to compute overlap contour pars");
			#endif
		}	
		else{
			//Sort contour points counter-clockwise
			overlapContour->SortPointsCounterClockWise();
		}
	}//close if

	
	return 0;

}//close ComputeContourOverlapArea()

int MathUtils::ComputePolygonOverlapArea(double& overlapArea,polygon_2d& overlap_poly,polygon_2d& poly, polygon_2d& poly2)
{
	//Init overlap area
	overlapArea = 0.;  

	//Find polygon intersection  
  std::deque<polygon_2d> output;
  bool ret = boost::geometry::intersection(poly, poly2, output);
  if(!ret) {
		#ifdef LOGGING_ENABLED
  		WARN_LOG("Polygons not intersecting or could not calculate the overlap, returning -1");
		#endif
		overlapArea= -1;
    return -1;
  }

	//Compute overlap area
	BOOST_FOREACH(overlap_poly, output) {
 		overlapArea = boost::geometry::area(overlap_poly);
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Polygon intersection: "<<boost::geometry::wkt(overlap_poly));
		#endif
  }

	//Check number of points
	size_t npoints= boost::geometry::num_points(overlap_poly);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<npoints<<" in polygon intersection");
  #endif

	//Correct overlap polygon
  boost::geometry::correct(overlap_poly);
	#ifdef LOGGING_ENABLED
 		DEBUG_LOG("Polygon intersection: "<<boost::geometry::wkt(overlap_poly));
	#endif

	return 0;

}//close ComputePolygonOverlapArea()

int MathUtils::Contour2Polygon(polygon_2d& poly,Contour* contour)
{
	//Check given contour
	if(!contour){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Null ptr to contour given!");
		#endif
		return -1;	
	}
	
	//Get contour points
	std::vector<TVector2> contourPoints= contour->GetPoints();

	//Loop over contour pointd and convert to polygon
	for(size_t i=0;i<contourPoints.size();i++){
		double x= contourPoints[i].X();
		double y= contourPoints[i].Y();
		boost::geometry::append(poly, boost::geometry::make<point_2d>(x,y));
	}//end loop points

	//Correct polygon
  boost::geometry::correct(poly);
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Polygon: "<<boost::geometry::wkt(poly));
	#endif

	return 0;

}//close Contour2Polygon()

int MathUtils::Ellipse2Polygon(polygon_2d& poly,double xc, double yc, double a, double b, double theta, int n)
{
  //Check n
	if(!n) {
		#ifdef LOGGING_ENABLED
 			ERROR_LOG("n should be >0");
		#endif
    return -1;
  }

	double w = theta*TMath::DegToRad();//convert theta to radians

	//Initialize polygon points
  double coor[n*2+1][2];		
	const std::size_t points_size = n*2+1;
	boost::scoped_array<point_xy> points(new point_xy[points_size]); 

	//Loop over n points and fill polygon
  double x= 0;
	double y= 0;
	double t = 0;
  double step = TMath::Pi()/n;
  double sinphi = sin(w);
  double cosphi = cos(w);
  for(int i=0; i<2*n+1; i++) {   
  	x = xc + a*cos(t)*cosphi - b*sin(t)*sinphi;
    y = yc + a*cos(t)*sinphi + b*sin(t)*cosphi;
    if(fabs(x) < 1e-4) x = 0;
    if(fabs(y) < 1e-4) y = 0;

    coor[i][0] = x;
    coor[i][1] = y;
		points[i]= point_xy(x,y);
    t += step;
  }//end loop points

	//Fill polygon
  boost::geometry::assign_points(poly,std::make_pair(&points[0], &points[0] + points_size));

	//Correct polygon
  boost::geometry::correct(poly);
    
	return 0;

}//close Ellipse2Polygon()

int MathUtils::Ellipse2Polygon(polygon_2d& poly,TEllipse* ellipse, int n)
{
	//Check ellipse
	if(!ellipse) {
		#ifdef LOGGING_ENABLED
 			ERROR_LOG("Null ptr to ellipse given!");
		#endif
    return -1;
  }

	//Get ellipse pars
	double Cx= ellipse->GetX1();
	double Cy= ellipse->GetY1();
	double theta= ellipse->GetTheta();
	double R1= ellipse->GetR1();
	double R2= ellipse->GetR2(); 
	
	return Ellipse2Polygon(poly,Cx,Cy,R1,R2,theta,n);

}//close Ellipse2Polygon()

double MathUtils::ComputeEllipseEccentricity(TEllipse* ellipse)
{
	//Check ellipse
	if(!ellipse) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given ellipse!");
		#endif
		return -1;
	}

	//Compute eccentricity
	double r1= ellipse->GetR1();
	double r2= ellipse->GetR2();
	double b= std::max(r1,r2);//semi-major axis
	double a= std::min(r1,r2);//semi-minor axis
	double a2= a*a;
	double b2= b*b;
	double e= 1-a2/b2; 

	return e;

}//close GetEllipseEccentricity()

double MathUtils::ComputeEllipseEccentricity(double bmaj,double bmin)
{	
	//NB: bmaj, bmin are major-axis
	//Check if they need to be switched
	double bmaj_half= bmaj/2.;
	double bmin_half= bmin/2.;
	double b= std::max(bmaj_half,bmin_half);//semi-major axis
	double a= std::min(bmaj_half,bmin_half);
	double a2= a*a;
	double b2= b*b;
	double e= 1-a2/b2;

	return e;

}//close ComputeEllipseEccentricity()


}//close namespace



