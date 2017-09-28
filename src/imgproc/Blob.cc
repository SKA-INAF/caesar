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
* @file Blob.cc
* @class Blob
* @brief Blob class
*
* Class representing an image blob
* @author S. Riggi
* @date 20/01/2015
*/


#include <Blob.h>
#include <Image.h>
#include <Pixel.h>
#include <Contour.h>
#include <StatsUtils.h>
#include <Logger.h>
#include <ZernikeMoments.h>

#include <TObject.h>
#include <TMatrixD.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

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
#include <climits>

using namespace std;

ClassImp(Caesar::Blob)

namespace Caesar {

Blob::Blob() : TNamed() 
{
	//Initialize
	Init();
	
}//close costructor

Blob::Blob(std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Init pars
	Init();

}//close constructor


Blob::Blob(std::vector<Pixel*>const& pixels,std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Init pars
  Init();

	//Fill pixels	
	bool makeCopy= false;
	for(size_t i=0;i<pixels.size();i++){
		if(pixels[i] && AddPixel(pixels[i],makeCopy)<0){
			ERROR_LOG("Failed to add pixel no. "<<i<<" to blob, skip to next!");
			continue;
		}
	}

}//close constructor

Blob::~Blob()
{
	//Clear pixels
	ClearPixels();
	ClearContours();	

}//close destructor


Blob::Blob(const Blob& blob) {
  // Contour copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((Blob&)blob).Copy(*this);
}


void Blob::Copy(TObject &obj) const {

	//Copy TNamed object
	DEBUG_LOG("Copying parent TNamed...");
	TNamed::Copy((Blob&)obj);

	// Copy this blob to blob
  ((Blob&)obj).HasPixelsAtEdge = HasPixelsAtEdge;
	((Blob&)obj).Id = Id;
	((Blob&)obj).NPix = NPix;
	((Blob&)obj).Mean = Mean;
	((Blob&)obj).RMS = RMS;
	((Blob&)obj).Skewness = Skewness;
	((Blob&)obj).Median = Median;
	((Blob&)obj).MedianRMS = MedianRMS;
	((Blob&)obj).X0 = X0;
	((Blob&)obj).Y0 = Y0;
	((Blob&)obj).Mean_curv = Mean_curv;
	((Blob&)obj).RMS_curv = RMS_curv;
	((Blob&)obj).Median_curv = Median_curv;
	((Blob&)obj).MedianRMS_curv = MedianRMS_curv;
	((Blob&)obj).Moments = Moments;
	((Blob&)obj).HuMoments = HuMoments;
	((Blob&)obj).ZMMoments = ZMMoments;

	((Blob&)obj).m_HasStats = m_HasStats;
	((Blob&)obj).m_HasParameters = m_HasParameters;
	((Blob&)obj).m_M1 = m_M1;
	((Blob&)obj).m_M2 = m_M2;
	((Blob&)obj).m_M3 = m_M3;
	((Blob&)obj).m_M4 = m_M4;
	((Blob&)obj).m_M1_curv = m_M1_curv;
	((Blob&)obj).m_M2_curv = m_M2_curv;
	((Blob&)obj).m_S = m_S;
	((Blob&)obj).m_Smax = m_Smax;
	((Blob&)obj).m_Smin = m_Smin;
	((Blob&)obj).m_Sxx = m_Sxx;
	((Blob&)obj).m_Syy = m_Syy;
	((Blob&)obj).m_Sxy = m_Sxy;
	((Blob&)obj).m_Sx = m_Sx;
	((Blob&)obj).m_Sy = m_Sy;
	((Blob&)obj).m_PixIdmax = m_PixIdmax;
	((Blob&)obj).m_PixIdmin = m_PixIdmin;
	((Blob&)obj).m_S_curv = m_S_curv;
	((Blob&)obj).m_S_edge = m_S_edge;

	//((Blob&)obj).m_ImgRange = m_ImgRange;
	//((Blob&)obj).m_ImageSizeX = m_ImageSizeX;
	//((Blob&)obj).m_ImageSizeY = m_ImageSizeY;
	//((Blob&)obj).m_ImageMinX = m_ImageMinX;
	//((Blob&)obj).m_ImageMaxX = m_ImageMaxX;
	//((Blob&)obj).m_ImageMinY = m_ImageMinY;
	//((Blob&)obj).m_ImageMaxY = m_ImageMaxY;
	//((Blob&)obj).m_ImageMinS = m_ImageMinS;
	//((Blob&)obj).m_ImageMaxS = m_ImageMaxS;
	//((Blob&)obj).m_ImageMinScurv = m_ImageMinScurv;
	//((Blob&)obj).m_ImageMaxScurv = m_ImageMaxScurv;
	//((Blob&)obj).m_ImageMinSedge = m_ImageMinSedge;
	//((Blob&)obj).m_ImageMaxSedge = m_ImageMaxSedge;
	//((Blob&)obj).m_ImageRMS = m_ImageRMS;
	
	((Blob&)obj).m_Xmin = m_Xmin;
	((Blob&)obj).m_Xmax = m_Xmax;
	((Blob&)obj).m_Ymin = m_Ymin;
	((Blob&)obj).m_Ymax = m_Ymax;
	((Blob&)obj).m_Ix_min = m_Ix_min;
	((Blob&)obj).m_Ix_max = m_Ix_max;
	((Blob&)obj).m_Iy_min = m_Iy_min;
	((Blob&)obj).m_Iy_max = m_Iy_max;
	((Blob&)obj).HasPixelsAtEdge = HasPixelsAtEdge;
	((Blob&)obj).HasPixelsAtEdge = HasPixelsAtEdge;
				
	//Copy pixel collection
	//Delete first any existing collection
	for(unsigned int i=0;i<(((Blob&)obj).m_Pixels).size();i++){
		if( (((Blob&)obj).m_Pixels)[i] ){
			delete (((Blob&)obj).m_Pixels)[i];
			(((Blob&)obj).m_Pixels)[i]= 0;
		}
	}
	(((Blob&)obj).m_Pixels).clear();

	Pixel* aPixel= 0;
	for(unsigned int i=0;i<m_Pixels.size();i++){
		aPixel= new Pixel;
		*aPixel= *(m_Pixels[i]);
		(((Blob&)obj).m_Pixels).push_back(aPixel);
	}

		
	//Copy contour collection
	//Delete first any existing collection
	for(unsigned int i=0;i<(((Blob&)obj).m_Contours).size();i++){
		if( (((Blob&)obj).m_Contours)[i] ){
			delete (((Blob&)obj).m_Contours)[i];
			(((Blob&)obj).m_Contours)[i]= 0;
		}
	}
	(((Blob&)obj).m_Contours).clear();

	Contour* aContour= 0;
	for(unsigned int i=0;i<m_Contours.size();i++){
		aContour= new Contour;
		*aContour= *(m_Contours[i]);
		(((Blob&)obj).m_Contours).push_back(aContour);
	}

}//close Copy()

Blob& Blob::operator=(const Blob& blob) { 
	// Operator =
  if (this != &blob)  ((Blob&)blob).Copy(*this);
  return *this;
}


void Blob::Init(){

	//Initialize parameters & data
	Id= -999;
	//this->SetNameTitle("","");
	HasPixelsAtEdge= false;
	ResetStats();
	ResetMoments();
	ClearPixels();
	ClearContours();	

}//close Init()

void Blob::ResetStats(){

	//Reset stats
	NPix= 0;
	Mean= 0;
	RMS= 0;
	Skewness= 0;
	Median= 0;
	MedianRMS= 0;

	X0= 0;
	Y0= 0;

	Mean_curv= 0;
	RMS_curv= 0;
	Median_curv= 0;
	MedianRMS_curv= 0;
	
	m_HasParameters= false;
	m_HasStats= false;
	
}//close ResetStats()

void Blob::ClearPixels(){

	for(size_t i=0;i<m_Pixels.size();i++){
		if(m_Pixels[i]){
			delete m_Pixels[i];
			m_Pixels[i]= 0;
		}
	}
	m_Pixels.clear();
	
}//close ClearPixels()


void Blob::ClearContours(){

	for(size_t i=0;i<m_Contours.size();i++){
		if(m_Contours[i]){
			delete m_Contours[i];
			m_Contours[i]= 0;
		}
	}
	m_Contours.clear();
	
}//close ClearContours()

void Blob::ResetMoments(){
	
	//Reset moments
	m_M1= 0;
  m_M2= 0;
	m_M3= 0;
	m_M4= 0;
	m_M1_curv= 0;
  m_M2_curv= 0;

	//Reset stat accumulator pars
	NPix= 0;	  
	X0= 0; 
	Y0= 0;
	m_S= 0; m_S_curv= 0.; m_S_edge= 0.;
	m_Smax= -1.e+99;
	m_Smin= 1.e+99;
	m_Sxx= 0; m_Syy= 0; m_Sxy= 0;
	m_Sx= 0; m_Sy= 0;
	m_PixIdmax= -1; m_PixIdmin= -1;
	
	//Reset image ranges
	m_Xmin= std::numeric_limits<float>::infinity();
	m_Xmax= -std::numeric_limits<float>::infinity();
	m_Ymin= std::numeric_limits<float>::infinity();
	m_Ymax= -std::numeric_limits<float>::infinity();

	m_Ix_min= LONG_MAX;
	m_Ix_max= LONG_MIN;
	m_Iy_min= LONG_MAX;
	m_Iy_max= LONG_MIN;
	
	m_HasStats= false;

}//close ResetMoments()


void Blob::UpdateMoments(Pixel* pixel){

	if(!pixel) return;

	//Get pixel data
	double w= pixel->S;
	double x= pixel->x;	
	double y= pixel->y;
	int ix= pixel->ix;
	int iy= pixel->iy;
	int id= pixel->id;
	double w_edge= pixel->S_edge;
	double w_curv= pixel->S_curv;
	
	//Update accumulator
	if(w<m_Smin) {
		m_Smin= w;
		m_PixIdmin= id;
	}
	if(w>m_Smax) {
		m_Smax= w;
		m_PixIdmax= id;
	}
	m_S+= w;
	m_Sx+= w*x;
	m_Sy+= w*y;
	m_Sxx+= w*x*x;
	m_Syy+= w*y*y;
	m_Sxy+= w*x*y;
	m_S_curv+= w_curv;
	m_S_edge+= w_edge;	
	X0+= x;
	Y0+= y;

	if(x<m_Xmin) m_Xmin= x;
	if(x>m_Xmax) m_Xmax= x;
	if(y<m_Ymin) m_Ymin= y;
	if(y>m_Ymax) m_Ymax= y;

	if(ix<m_Ix_min) m_Ix_min= ix;
	if(ix>m_Ix_max) m_Ix_max= ix;
	if(iy<m_Iy_min) m_Iy_min= iy;
	if(iy>m_Iy_max) m_Iy_max= iy;

	
	//Update moments
	NPix++;
  double delta = w - m_M1;
  double delta_n = delta/NPix;
  double delta_n2 = delta_n * delta_n;
  double f = delta * delta_n * (NPix-1);
  m_M1+= delta_n;
  m_M4+= f * delta_n2 * (NPix*NPix - 3*NPix + 3) + 6 * delta_n2 * m_M2 - 4 * delta_n * m_M3;
  m_M3+= f * delta_n * (NPix - 2) - 3 * delta_n * m_M2;
  m_M2+= f;	

	double delta_curv = w_curv - m_M1_curv;
  double delta_curv_n = delta_curv/NPix;
  double f_curv = delta_curv * delta_curv_n * (NPix-1);
  m_M1_curv+= delta_curv_n;
  m_M2_curv+= f_curv;	

}//close UpdateMoments()

int Blob::AddPixel(Pixel* aPixel,bool makeCopy){

	//Check pixel
	if(!aPixel){
		ERROR_LOG("Null ptr to given pixel, nothing will be added!");
		return -1;
	}

	//Append pixel to list
	Pixel* pix= aPixel;
	if(makeCopy){
		Pixel* clonedPixel= new Pixel;
		*clonedPixel= *aPixel;
		pix= clonedPixel;
	}
	m_Pixels.push_back(pix);

	//Update moment counts
	UpdateMoments(pix);

	return 0;

}//close AddPixel()




int Blob::ComputeStats(bool computeRobustStats,bool forceRecomputing){

	//## Compute region stats & shape parameters
	if(NPix<=0) return -1;

	//## If stats was already computed before return (unless forced)
	//## NB: Without this is stats are computed again the X0, Y0 etc are scaled again by Npix ==> WRONG!!!!
	//##     This was a bug in old version
	//if(m_HasStats && !forceRecomputing) {
	//	ERROR_LOG("Blob has already stats computed (you need to force recomputation)!");
	//	return -1;
	//}

	//## Recomputing moments?
	if(forceRecomputing){
		ResetMoments();//reset moments
		for(unsigned int k=0;k<m_Pixels.size();k++) UpdateMoments(m_Pixels[k]);//update moments
	}
		
	//## NB: Compute stats parameters only if not already done or if forced
	//##     This was a bug in older version affecting saliency calculation (the spatial component in region distance) and eventually the extended source extraction!!!
	//##     Need to introduce a rescaling of spatial vs color distance (a factor 100 more or less?) 
	if((signed)(m_Pixels.size())!=NPix){
		ERROR_LOG("Mismatch between number of pixels present in collection ("<<m_Pixels.size()<<") and NPix="<<NPix<<" (fix bug!!!)");
		return -1;
	}

	if(!m_HasStats || forceRecomputing){
		X0/= (double)(NPix);
		Y0/= (double)(NPix);
		m_Sx/= m_S;
		m_Sy/= m_S;
		m_Sxx/= m_S;
		m_Syy/= m_S;
		m_Sxy/= m_S;
	}

	Mean= m_M1;
	Mean_curv= m_M1_curv;
	RMS= 0;
	RMS_curv= 0;	
	if(NPix>1) {
		RMS= sqrt(m_M2/(NPix-1));
		RMS_curv= sqrt(m_M2_curv/(NPix-1));
	}

	Skewness= 0;
	if(m_M2!=0) {
  	Skewness= sqrt(NPix)*m_M3/pow(m_M2,1.5);//need to adjust for finite population?
	}
	
	//## Compute robust stats (median, MAD, Entropy, ...)
	if(computeRobustStats){	
		std::vector<double> SList;
		std::vector<double> SCurvList;
		
		for(size_t i=0;i<m_Pixels.size();i++){
			double S= m_Pixels[i]->S;		
			double S_curv= m_Pixels[i]->S_curv;
			SList.push_back(S);
			SCurvList.push_back(S_curv);
		}
		Median= StatsUtils::GetMedianFast<double>(SList);
		MedianRMS= 1.4826*StatsUtils::GetMADFast(SList,Median);
		Median_curv= StatsUtils::GetMedianFast<double>(SCurvList);
		MedianRMS_curv= 1.4826*StatsUtils::GetMADFast(SCurvList,Median);
	}//close if computeRobustStats

	m_HasStats= true;

	return 0;

}//close ComputeStats()

int Blob::ComputeMorphologyParams(){

	DEBUG_LOG("Computing blob morphology parameters...");
	if(NPix<=0 || m_Pixels.size()<=0) return -1;

	//## Reset existing pars
	m_HasParameters= false;
	DEBUG_LOG("Clear & reset existing contours (if any)...");
	for(size_t i=0;i<m_Contours.size();i++){
		if(m_Contours[i]){
			delete m_Contours[i];
			m_Contours[i]= 0;
		}
	}
	m_Contours.clear();
	Moments.clear();
	HuMoments.clear();
	ZMMoments.clear();
		
	//######################################
	//## Find the source bounding box
	//######################################
	float xRange[2]= {m_Xmin,m_Xmax};
	float yRange[2]= {m_Ymin,m_Ymax};	
	long int ixRange[2]= {m_Ix_min,m_Ix_max};
	long int iyRange[2]= {m_Iy_min,m_Iy_max};
	
	//Bounding box in (x,y) coordinates
	long int boundingBoxX[2];
	long int boundingBoxY[2];
	int deltaPix= 50;
	boundingBoxX[0]= xRange[0]-deltaPix;
	boundingBoxX[1]= xRange[1]+deltaPix;
	boundingBoxY[0]= yRange[0]-deltaPix;
	boundingBoxY[1]= yRange[1]+deltaPix;
	long int nBoxX= boundingBoxX[1]-boundingBoxX[0]+1;
	long int nBoxY= boundingBoxY[1]-boundingBoxY[0]+1;
	
	//Bounding box in (ix,iy) coordinates
	long int boundingBoxIX[2];
	long int boundingBoxIY[2];	
	boundingBoxIX[0]= ixRange[0]-deltaPix;
	boundingBoxIX[1]= ixRange[1]+deltaPix;
	boundingBoxIY[0]= iyRange[0]-deltaPix;
	boundingBoxIY[1]= iyRange[1]+deltaPix;
	long int nBoxIX= boundingBoxIX[1]-boundingBoxIX[0]+1;
	long int nBoxIY= boundingBoxIY[1]-boundingBoxIY[0]+1;

	DEBUG_LOG("xRange("<<xRange[0]<<","<<xRange[1]<<"), yRange("<<yRange[0]<<","<<yRange[1]<<")");
	DEBUG_LOG("ixRange("<<ixRange[0]<<","<<ixRange[1]<<"), iyRange("<<yRange[0]<<","<<iyRange[1]<<")");
	DEBUG_LOG("boundingBoxX("<<boundingBoxX[0]<<","<<boundingBoxX[1]<<"), boundingBoxY("<<boundingBoxY[0]<<","<<boundingBoxY[1]<<")"<<"  nBoxX="<<nBoxX<<", nBoxY="<<nBoxY);
	DEBUG_LOG("boundingBoxIX("<<boundingBoxIX[0]<<","<<boundingBoxIX[1]<<"), boundingBoxIY("<<boundingBoxIY[0]<<","<<boundingBoxIY[1]<<")"<<"  nBoxIX="<<nBoxIX<<", nBoxIY="<<nBoxIY);

	//## Fill image and binarized image
	//cv::Mat binarizedImg = cv::Mat::zeros(nBoxIY, nBoxIX, CV_8UC1);
	//cv::Mat rasterImg = cv::Mat::zeros(nBoxIY, nBoxIX, CV_64FC1);
	cv::Mat binarizedImg = cv::Mat::zeros(nBoxY, nBoxX, CV_8UC1);//last
	cv::Mat rasterImg = cv::Mat::zeros(nBoxY, nBoxX, CV_64FC1);//last

	for(unsigned int k=0;k<m_Pixels.size();k++){
		Pixel* thisPixel= m_Pixels[k];
		double thisS= thisPixel->S;	
		double x= thisPixel->x;
		double y= thisPixel->y;
		x-= boundingBoxX[0];
		y-= boundingBoxY[0];	
		
		long int ix= thisPixel->ix;
		long int iy= thisPixel->iy;	
		ix-= boundingBoxIX[0];
		iy-= boundingBoxIY[0];

		//long int rowId= nBoxIY-1-iy;
		//long int colId= nBoxIX-1-ix;
		long int rowId= nBoxY-1-y;//TEST (LAST)
		long int colId= nBoxX-1-x;//TEST (LAST)
		
		binarizedImg.at<uchar>(rowId, colId, 0) = 1;
		rasterImg.at<double>(rowId, colId, 0) = thisS;		
	}//end loop pixels


	//## Compute contour
	cv::Mat binarizedImg_clone= binarizedImg.clone();//image will be modified by findContours! 
	std::vector<std::vector<cv::Point>> contours; // Vector for storing contour
  std::vector<cv::Vec4i> hierarchy;
	cv::findContours( binarizedImg_clone, contours, hierarchy,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE, cv::Point(0,0) ); // Find only the external contours in the image
  //cv::findContours( binarizedImg_clone, contours, hierarchy,CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE ); // Find the contours in the image organizing in a 2-level hierarchy
	//cv::findContours( binarizedImg_clone, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE, Point(0,0) );
	
	Contour* aContour= 0;

	for(size_t i=0; i<contours.size(); i++){ // iterate through each contour
		int nContourPts= (int)contours[i].size();
		if(nContourPts<=0) continue;

		//Create and fill contour
		aContour= new Contour;

		std::stringstream sstream;
		sstream<<"Contour no. "<<i+1<<": (";
		for(int j=0;j<nContourPts;j++){
			int contx= contours[i][j].x;
			int conty= contours[i][j].y;

			/*		
			int rowId= nBoxIY-1-conty;
			int colId= nBoxIX-1-contx;
			int contx_transf= colId + boundingBoxX[0];
			int conty_transf= rowId + boundingBoxY[0];
			*/
			
			//LAST
			long int rowId= nBoxY-1-conty;
			long int colId= nBoxX-1-contx;
			long int contx_transf= colId + boundingBoxX[0];
			long int conty_transf= rowId + boundingBoxY[0];
			
			aContour->AddPoint(TVector2(contx_transf,conty_transf));
			sstream<<"("<<contx_transf<<","<<conty_transf<<"), ";
		}//end loop points in contour
		sstream<<")";
		DEBUG_LOG(sstream.str());
		

		//Compute contour parameters
		if(aContour->ComputeParameters()<0){
			WARN_LOG("One/more failures occurred while computing contour no. "<<i<<" parameters for blob id "<<Id<<"!");
		}
		
		//Add contour to the list
		m_Contours.push_back(aContour);	

	}//end loop contours


	//## Compute HuMoments
	cv::Moments moments= cv::moments(rasterImg, false);
	Moments.push_back(moments.m00);
	Moments.push_back(moments.m10);
	Moments.push_back(moments.m01);
	Moments.push_back(moments.m20);
	Moments.push_back(moments.m11);
	Moments.push_back(moments.m02);
	Moments.push_back(moments.m30);
	Moments.push_back(moments.m21);
	Moments.push_back(moments.m12);
	Moments.push_back(moments.m03);
	
	Moments.push_back(moments.mu20);
	Moments.push_back(moments.mu11);
	Moments.push_back(moments.mu02);
	Moments.push_back(moments.mu30);
	Moments.push_back(moments.mu21);
	Moments.push_back(moments.mu12);
	Moments.push_back(moments.mu03);
	
	Moments.push_back(moments.nu20);
	Moments.push_back(moments.nu11);
	Moments.push_back(moments.nu02);
	Moments.push_back(moments.nu30);
	Moments.push_back(moments.nu21);
	Moments.push_back(moments.nu12);
	Moments.push_back(moments.nu03);
		

	double hu_moments[7];
	cv::HuMoments(moments, hu_moments);
	for(int i=0;i<7;i++) HuMoments.push_back(hu_moments[i]);

	//## Compute zernike moments
	//int order= 6;
	//ComputeZernikeMoments(order);
		
	m_HasParameters= true;

	return 0;

}//close ComputeMorphologyParams()

Image* Blob::GetImage(ImgType mode,int pixMargin){

	//Bounding box in (x,y) coordinates
	double xRange[2]= {m_Xmin,m_Xmax};
	double yRange[2]= {m_Ymin,m_Ymax};	
	int boundingBoxX[2];
	int boundingBoxY[2];
	int deltaPix= pixMargin;//50;
	boundingBoxX[0]= xRange[0]-deltaPix;
	boundingBoxX[1]= xRange[1]+deltaPix;
	boundingBoxY[0]= yRange[0]-deltaPix;
	boundingBoxY[1]= yRange[1]+deltaPix;
	long int nBoxX= boundingBoxX[1]-boundingBoxX[0]+1;
	long int nBoxY= boundingBoxY[1]-boundingBoxY[0]+1;

	//## Fill image and binarized image
	TString imgName= Form("SourceImg_%s_mode%d",this->GetName(),mode);
	//Image* blobImg= new Image(nBoxX,boundingBoxX[0]-0.5,boundingBoxX[1]+0.5,nBoxY,boundingBoxY[0]-0.5,boundingBoxY[1]+0.5,imgName);
	Image* blobImg= new Image(nBoxX,nBoxY,boundingBoxX[0]-0.5,boundingBoxY[0]-0.5,std::string(imgName));
	
	if(mode==eBinaryMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			blobImg->Fill(thisX,thisY,1);
		}//end loop pixels
	}//close if
	
	else if(mode==eFluxMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			double thisS= thisPixel->S;	
			blobImg->Fill(thisX,thisY,thisS);
		}//end loop pixels
	}//close else
	
	else if(mode==eSignificanceMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			double thisS= thisPixel->S;		
			double thisBkg= thisPixel->bkgLevel;
			double thisNoise= thisPixel->noiseLevel;
			double thisZ= 0;
			if(thisNoise!=0) thisZ= (thisS-thisBkg)/thisNoise;	
			blobImg->Fill(thisX,thisY,thisZ);
		}//end loop pixels
	}//close else
	
	else if(mode==ePullMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			double thisS= thisPixel->S;	
			double thisPull= (thisS-Median)/MedianRMS;
			blobImg->Fill(thisX,thisY,thisPull);
		}//end loop pixels
	}//close else if
	
	else if(mode==eCurvatureMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			double thisCurv= thisPixel->S_curv;	
			blobImg->Fill(thisX,thisY,thisCurv);
		}//end loop pixels
	}//close else if
	
	
	else if(mode==eMeanFluxMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			blobImg->Fill(thisX,thisY,Mean);
		}//end loop pixels
	}

	else if(mode==eBkgMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			double thisBkg= thisPixel->bkgLevel;
			blobImg->Fill(thisX,thisY,thisBkg);
		}//end loop pixels
	}

	else if(mode==eNoiseMap){
		for(size_t k=0;k<m_Pixels.size();k++){
			Pixel* thisPixel= m_Pixels[k];
			double thisX= thisPixel->x;
			double thisY= thisPixel->y;
			double thisNoise= thisPixel->noiseLevel;
			blobImg->Fill(thisX,thisY,thisNoise);
		}//end loop pixels
	}

	return blobImg;

}//close GetImage()



int Blob::ComputeZernikeMoments(int order){

	//## Get source image
	Image* fluxMap= GetImage(eFluxMap);
	if(!fluxMap) return -1;
	
	//## Compute Zernike moments
	double radius= -1;
	ZMMoments.clear();
	ZMMoments= ZernikeMoments::GetZernike2D_Direct(fluxMap, order,radius);
	//ZMMoments= ZernikeMoments::GetZernike2D(fluxMap, order,radius);

	delete fluxMap;
	fluxMap= 0;

	return 0;

}//close ComputeZernikeMoments()


}//close namespace


