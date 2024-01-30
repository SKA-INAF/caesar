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
* @file MorphFilter.cc
* @class MorphFilter
* @brief Class implementing morphological filtering
*
* Morphological Filter
* @author S. Riggi
* @date 20/01/2015
*/


#include <MorphFilter.h>
#include <Image.h>
#include <MathUtils.h>
#include <Source.h>
#include <Pixel.h>
#include <BkgData.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <Graph.h>
#include <Contour.h>

#include <rtnorm.h>

//ROOT headers
#include <TObject.h>
#include <TVector2.h>
//#include <TRInterface.h>

//OpenCV headers
#include <opencv2/imgproc/types_c.h>

//C++ headers
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <cfloat>

using namespace std;

ClassImp(Caesar::MorphFilter)

namespace Caesar {

MorphFilter::MorphFilter() {

}//close costructor


MorphFilter::~MorphFilter(){

}//close destructor

Image* MorphFilter::ComputeWatershedFilter(std::vector<Contour*>& contours,Image* img,Image* markerImg)
{
	//## Init
	contours.clear();

	//## Compute watershed filter map
	Image* morphImg= ComputeWatershedFilter(img,markerImg);
	if(!morphImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute Watershed filter map!");
		#endif
		return nullptr;
	}

	//## Get image data
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	double Xmin= img->GetXmin();
	double Ymin= img->GetYmin();
	
	//## Fill binary map to compute contours
	//cv::Mat mark = cv::Mat::zeros(markers.size(), CV_8UC1);
 	cv::Mat mark = cv::Mat::zeros(Ny,Nx,CV_8UC1);

	for(long int j=0;j<Ny;j++){
		long int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			long int colId= i;
			double binContent= img->GetPixelValue(i,j);
			if(binContent==0) continue;
			double w= morphImg->GetPixelValue(i,j);
			if(w>1) mark.at<uchar>(rowId, colId, 0) = 1;
		}//end loop x
	}//end loop y

	//## Compute contours
	std::vector<std::vector<cv::Point>> markerContours; // Vector for storing contour
  std::vector<cv::Vec4i> hierarchy;
	cv::findContours(mark, markerContours, hierarchy,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE, cv::Point(0,0) );
	
	Contour* aContour= 0;

	for(size_t i=0; i<markerContours.size(); i++){ // iterate through each contour
		int nContourPts= (int)markerContours[i].size();
		if(nContourPts<=0) continue;

		//Create and fill contour
		aContour= new Contour;

		for(int j=0;j<nContourPts;j++){
			int contx= markerContours[i][j].x;
			int conty= markerContours[i][j].y;
			long int rowId= Ny-1-conty;
			//long int colId= Nx-1-contx;
			long int colId= contx;
			long int x= colId + Xmin;
			long int y= rowId + Ymin;

			//double x= img->GetX(contx);
			//double y= img->GetY(conty);
			aContour->AddPoint(TVector2(x,y));
		}//end loop points in contour
		
		//Compute contour parameters
		if(aContour->ComputeParameters()<0){
			#ifdef LOGGING_ENABLED	
				WARN_LOG("One/more failures occurred while computing contour no. "<<i<<"!");
			#endif
		}

		//Compute fitted ellipse
		if(aContour->ComputeFittedEllipse()<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Ellipse fitting to contour no. "<<i<<" failed!");
			#endif
		}
		
		//Add contour to the list
		contours.push_back(aContour);	

	}//end loop contours

	return morphImg;

}//close ComputeWatershedFilter()


Image* MorphFilter::ComputeWatershedFilter(Image* img,Image* markerImg)
{
	//Check input images
	if(!img || !markerImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr given to image and/or marker image!");
		#endif
		return nullptr;
	}

	//Check input image have same dimensions
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	long int Nx_marker= markerImg->GetNx();
	long int Ny_marker= markerImg->GetNy();
	if(Nx!=Nx_marker || Ny!=Ny_marker){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Input image and given marker image have different dimensions!");
		#endif
		return nullptr;
	}
	double Xmin= img->GetXmin();
	double Ymin= img->GetYmin();
	
	// Convert images to OpenCV format
	cv::Mat mat= img->GetOpenCVMat("32");
	cv::Mat mat_col;	
	cv::Mat src;
	//cv::cvtColor(mat,mat_col,CV_GRAY2BGR);//deprecated C API
	cv::cvtColor(mat,mat_col,cv::COLOR_GRAY2BGR);
	mat_col.convertTo(src,CV_8UC3); 
	
	cv::Mat mat_markers = markerImg->GetOpenCVMat("32I");
	cv::Mat markers;
	mat_markers.convertTo(markers,CV_32SC1); 
	long int nRows = mat.rows;
  long int nCols = mat.cols;

	// Perform the watershed algorithm
	cv::watershed(src, markers);

	// Create the result image
	Image* morphImg= (Image*)img->GetCloned("",true,true);
	morphImg->Reset();

	for(long int j=0;j<Ny;j++){
		long int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			long int colId= i;
			double binContent= img->GetPixelValue(i,j);
			if(binContent==0) continue;
			int clusterId = markers.at<int>(rowId,colId);
			if (clusterId>1) {
				//morphImg->FillPixel(i,j,clusterId);
				morphImg->FillPixel(i,j,1);
			}
			else{	
				morphImg->FillPixel(i,j,0);
			}
		}//end loop x
	}//end loop y

	return morphImg;

}//close ComputeWatershedFilter()


Image* MorphFilter::ComputeHDomeFilter(Image* img,double baseline,int kernSize)
{
	//Compute image morph reconstruction
	Image* morphRecoImg= ComputeMorphRecoFilter(img,baseline,kernSize);
	if(!morphRecoImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute morph reco filter map!");
		#endif
		return nullptr;
	}
	
	//Subtract morph reco from input map to get blob mask
	Image* blobMask= (Image*)img->GetCloned("",true,true);
	blobMask->Add(morphRecoImg,-1);

	//Clear morphRecoImg
	CodeUtils::DeletePtr<Image>(morphRecoImg);

	return blobMask;

}//close ComputeHDomeFilter()

/*
Image* MorphFilter::ComputeHDomeFilter(Image* img,int kernSize)
{
	//Check image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image given!");
		#endif
		return nullptr;
	}

	//Compute pixel histo image valley threshold
	if(!img->HasStats()){
		img->ComputeStats(true);
	}
	double Smin= img->GetMinimum();
	double Smax= img->GetMaximum();
	double valleyThr1D= img->FindValleyThreshold();

	//Copy input image and invert it to find valleys
	Image* tmpMap= (Image*)img->GetCloned("",true,true);
	tmpMap->Scale(-1);

	//Normalize image
	int normmin= 1;
	int normmax= 256;
	bool skipEmptyBins= true;
	Image* searchMap= tmpMap->GetNormalizedImage("LINEAR",normmin,normmax,skipEmptyBins);
	
	//Clear tmp map
	CodeUtils::DeletePtr<Image>(tmpMap);

	//Find valleys
	std::vector<TVector2> valleyPoints;
	int peakShiftTolerance= 2;
	bool skipBorders= true;
	int peakKernelMultiplicityThr= 1;
	if(FindPeaks(valleyPoints,searchMap,{kernSize},peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to find valleys in input map!");	
		#endif
		CodeUtils::DeletePtr<Image>(searchMap);
		return nullptr;
	}

	//Clear search map
	CodeUtils::DeletePtr<Image>(searchMap);

	//Initialize baseline thr to 1D valley thr. If 2D valleys are found get the smaller one
	double baseline= valleyThr1D;
	if(valleyPoints.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No valley points found, setting baseline to 1D valley thr ("<<valleyThr1D<<") ...");
		#endif
	}
	else{

		//Get flux corresponding to each 2D valley point
		std::vector<double> valleyThrList;
		for(size_t i=0;i<valleyPoints.size();i++){
			double x= valleyPoints[i].X();
			double y= valleyPoints[i].Y();
			long int gbin= img->FindBin(x,y);
			if(gbin<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to find gbin of valley point ("<<x<<","<<y<<"), this should not occur!");
				#endif
				return nullptr;
			}
			double S= img->GetBinContent(gbin);
			if(std::isnormal(S) && S!=Smin){
				valleyThrList.push_back(S);
			}
		}//end loop

		if(!valleyThrList.empty()){
			//Sort valley thr descending and set baseline thr to smaller value
			std::sort(valleyThrList.begin(),valleyThrList.end());
			baseline= valleyThrList[0];
			#ifdef LOGGING_ENABLED
				INFO_LOG("Setting baseline to smaller valley point ("<<baseline<<") ...");
			#endif
		}
		else{
			#ifdef LOGGING_ENABLED
				WARN_LOG("No good valley points found, setting baseline to 1D valley thr ("<<valleyThr1D<<") ...");
			#endif
		}

	}//close else

	//Compute image morph reconstruction
	Image* morphRecoImg= ComputeMorphRecoFilter(img,baseline,kernSize);
	if(!morphRecoImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute morph reco filter map!");
		#endif
		return nullptr;
	}
	
	//Subtract morph reco from input map to get blob mask
	Image* blobMask= (Image*)img->GetCloned("",true,true);
	blobMask->Add(morphRecoImg,-1);

	//Clear morphRecoImg
	CodeUtils::DeletePtr<Image>(morphRecoImg);

	return blobMask;

}//close ComputeHDomeFilter()
*/

Image* MorphFilter::ComputeMorphRecoFilter(Image* img,Image* markerImg,int kernSize,double tol)
{
	//Check input image
	if(!img || !markerImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image given!");
		#endif
		return nullptr;
	}

	//Check input image have same dimensions
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	long int Nx_marker= markerImg->GetNx();
	long int Ny_marker= markerImg->GetNy();
	if(Nx!=Nx_marker || Ny!=Ny_marker){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Input image and given marker image have different dimensions!");
		#endif
		return nullptr;
	}
	
	//Check kern size
	if(kernSize%2==0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Given kern size is even, adding +1 to make it odd...");
		#endif
		kernSize++;
	}

	// Convert images to OpenCV format
	cv::Mat mask= img->GetOpenCVMat("64");	
	cv::Mat marker = markerImg->GetOpenCVMat("64");
	long int nRows = mask.rows;
  long int nCols = mask.cols;

	//## Init dilation kernel
	cv::Size kernel_size(kernSize,kernSize);
	cv::Mat element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(-1,-1));
	
	//## Perform morphological reconstruction (i.e., geodesic dilation until stability)
	cv::Mat rec = marker;
	cv::Mat rec_old= cv::Mat::zeros(nRows,nCols,CV_64FC1);
	double eps= cv::sum(rec - rec_old)[0];
	int iter_counter= 0;
	int max_iters= 100;

	while(eps>tol){
		//Retain output of previous iteration
   	rec_old = rec;
	
		//Perform dilation
		int iterations= 1;
		cv::Mat rec_dilated;
		cv::dilate(rec, rec_dilated, element, cv::Point(-1,-1),iterations,cv::BORDER_CONSTANT);
   	rec = rec_dilated;
		
		//Restrict the dilated values using the mask
		for(long int i=0;i<nRows;++i) {
  		double* pixels_rec = rec.ptr<double>(i);
			double* pixels_mask = mask.ptr<double>(i);
    	for (long int j=0;j<nCols;++j){
    		if(pixels_rec[j]>pixels_mask[j]){
					pixels_rec[j]= pixels_mask[j];
				}
    	}//end loop cols
  	}//end loop rows

		//Update eps
		eps= cv::sum(rec - rec_old)[0];
		iter_counter++;
		if(iter_counter>=max_iters) break;
	}//end while loop

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Image reconstruction completed in #"<<iter_counter<<" iterations...");
	#endif

	//## Convert back to image 
	Image* morphImg= (Image*)img->GetCloned("",true,true);
	morphImg->Reset();

	for(long int j=0;j<Ny;j++){
		long int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			long int colId= i;
			double binContent= img->GetPixelValue(i,j);
			if(binContent==0) continue;
			double w= rec.at<double>(rowId,colId);			
			morphImg->FillPixel(i,j,w);
		}//end loop x
	}//end loop y

	return morphImg;
	
}//close ComputeMorphRecoFilter()

Image* MorphFilter::ComputeMorphRecoFilter(Image* img,double baseline,int kernSize,double tol)
{
	//Check input image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image given!");
		#endif
		return nullptr;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	//Check baseline is positive
	if(baseline<=0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Given baseline value is <=0 (hint: it must be >0)!");
		#endif
		return nullptr;
	}
	
	//Create marker image
	Image* markerImg= (Image*)img->GetCloned("",true,true);
	markerImg->Reset();
	for(int i=0;i<Nx;i++){
		for(long int j=0;j<Ny;j++){
			double w= img->GetPixelValue(i,j);
			double w_marker= w-baseline;
			markerImg->FillPixel(i,j,w_marker);
		}//end loop bins y
	}//end loop bins x

	//Compute morph filter
	return ComputeMorphRecoFilter(img,markerImg,kernSize,tol);
	
}//close ComputeMorphRecoFilter()


Image* MorphFilter::ComputeMorphFilter(Image* img,int morphOp,int KernSize,int structElementType,int niters,bool skipZeroPixels) 
{
	//Check input image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image given!");
		#endif
		return nullptr;
	}
	
	//Check kern size
	if(KernSize%2==0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Given kern size is even, adding +1 to make it odd...");
		#endif
		KernSize++;
	}

	//## Convert image to OpenCV format
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	cv::Mat mat= img->GetOpenCVMat("64");

	//## Init struct element
	int ptSize= -1;//(KernSize-1)/2;
	cv::Size kernel_size(KernSize,KernSize);
	cv::Mat element;
	if(structElementType==eMORPH_RECT) element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(ptSize,ptSize));
	else if(structElementType==eMORPH_ELLIPSE) element= cv::getStructuringElement(cv::MORPH_ELLIPSE, kernel_size, cv::Point(ptSize,ptSize));
	else if(structElementType==eMORPH_CROSS) element= cv::getStructuringElement(cv::MORPH_CROSS, kernel_size, cv::Point(ptSize,ptSize));	
	else element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(ptSize,ptSize));
	
	//## Morphology operation
	//MORPH_OPEN - an opening operation
	//MORPH_CLOSE - a closing operation
	//MORPH_GRADIENT - a morphological gradient
	//MORPH_TOPHAT - “top hat”
	//MORPH_BLACKHAT - “black hat”
	cv::Mat mat_morph;
	if(morphOp==eMORPH_OPENING) cv::morphologyEx(mat, mat_morph, cv::MORPH_OPEN, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_CLOSING) cv::morphologyEx(mat, mat_morph, cv::MORPH_CLOSE, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_GRADIENT) cv::morphologyEx(mat, mat_morph, cv::MORPH_GRADIENT, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_TOPHAT) cv::morphologyEx(mat, mat_morph, cv::MORPH_TOPHAT, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_BLACKHAT) cv::morphologyEx(mat, mat_morph, cv::MORPH_BLACKHAT, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_EROSION) cv::erode(mat,mat_morph,element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);
	else if(morphOp==eMORPH_DILATION) cv::dilate(mat,mat_morph,element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid morph operation given ("<<morphOp<<")!");
		#endif
		return nullptr;
	}


	//## Convert back to image 
	Image* morphImg= (Image*)img->GetCloned("",true,true);
	morphImg->Reset();

	for(long int j=0;j<Ny;j++){
		long int rowId= Ny-1-j;
		for(int i=0;i<Nx;i++){
			long int colId= i;
			double binContent= img->GetPixelValue(i,j);
			if(binContent==0 && skipZeroPixels) continue;
			double w= mat_morph.at<double>(rowId,colId);			
			morphImg->FillPixel(i,j,w);
		}//end loop x
	}//end loop y

	return morphImg;

}//close ComputeMorphFilter()


int MorphFilter::FindPeaks(std::vector<ImgPeak>& peakPoints,Image* img,std::vector<int> kernelSizes,int peakShiftTolerance,bool skipBorders,int peakKernelMultiplicityThr)
{
	//Check input image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image given!");
		#endif
		return -1;
	}

	//Check multiplicity thr
	int nKernels= static_cast<int>(kernelSizes.size());
	if(peakKernelMultiplicityThr==-1){
		peakKernelMultiplicityThr= nKernels;
	}
	if(peakKernelMultiplicityThr>nKernels){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Multiplicy thr ("<<peakKernelMultiplicityThr<<") larger than nkernels ("<<nKernels<<"), set it to "<<nKernels<<"...");
		#endif
		peakKernelMultiplicityThr= nKernels;
	}
	if(peakKernelMultiplicityThr==0 || (peakKernelMultiplicityThr<0 && peakKernelMultiplicityThr!=-1)){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid multiplicity thr ("<<peakKernelMultiplicityThr<<") given!");
		#endif
		return -1;
	}

	//Init data
	peakPoints.clear();
	struct PeakData{
		float x;
		float y;
		int id;
		PeakData(float _x,float _y,float _id): x(_x), y(_y), id(_id) {}
	};
	std::vector<PeakData> peaks;

	//## Find peaks with different kernel sizes
	bool hasPeaks= true;
	for(int k=0;k<nKernels;k++){

		//Find peaks for this kernel
		std::vector<long int> peakPixelIds;
		Image* dilatedImg= MorphFilter::Dilate(peakPixelIds,img,kernelSizes[k],skipBorders);
		if(!dilatedImg){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to compute dilated image for kernel no. "<<k<<" (size="<<kernelSizes[k]<<")!");
			#endif
			return -1;
		}
		delete dilatedImg;
		dilatedImg= 0;		
		
		//Stop if no peaks are detected
		if(peakPixelIds.empty()){
			hasPeaks= false;
			break;
		}
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("#"<<peakPixelIds.size()<<" peak pixels found with kernel no. "<<k<<" (size="<<kernelSizes[k]<<")!");
		#endif

		//Store peak points for current kernel
		for(size_t j=0;j<peakPixelIds.size();j++){
			long int gBin= peakPixelIds[j];
			long int binX= img->GetBinX(gBin); 
			long int binY= img->GetBinY(gBin); 
			double x= img->GetX(binX);
			double y= img->GetY(binY);
			peaks.push_back(PeakData(x,y,k));
		}//end loop peak points

	}//end loop kernels

	if(!hasPeaks) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("No peaks detected in one or all dilated kernel runs.");
		#endif
		return 0;
	}

	//## Create graph with "matching/adjacent" peaks
	Graph linkGraph(peaks.size());
	for(size_t i=0;i<peaks.size()-1;i++) {
		for(size_t j=i+1;j<peaks.size();j++) {	
			float distX= peaks[i].x - peaks[j].x;
			float distY= peaks[i].y - peaks[j].y;
			bool areAdjacent= (fabs(distX)<=peakShiftTolerance && fabs(distY)<=peakShiftTolerance);
			if(!areAdjacent) continue;
			linkGraph.AddEdge(i,j);
		}
	}
	
	//Find connected peaks
	std::vector<std::vector<int>> connected_indexes;
	linkGraph.GetConnectedComponents(connected_indexes);

	
	//Match peaks found with different kernels (given a tolerance)
	int npeaks= 0;
	
	for(size_t i=0;i<connected_indexes.size();i++){
		std::set<int> id_list;

		for(size_t j=0;j<connected_indexes[i].size();j++){
			int index= connected_indexes[i][j];
			id_list.insert(peaks[index].id);
		}//end loop items in cluster
		
		//Check multiplicity
		int multiplicity= static_cast<int>(id_list.size());
		if( multiplicity<peakKernelMultiplicityThr ) {
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Skip peak group no. "<<i+1<<" as below multiplicity thr ("<<multiplicity<<"<"<<peakKernelMultiplicityThr<<")");
			#endif
			continue;
		}

		//Compute peak mean and fill final peaks list
		TVector2 peakMean(0,0);
		for(size_t j=0;j<connected_indexes[i].size();j++){
			int index= connected_indexes[i][j];
			peakMean+= TVector2(peaks[index].x,peaks[index].y);
		}
		peakMean*= 1./(float)(connected_indexes[i].size());

		//Find image bin x y
		double x= peakMean.X();
		double y= peakMean.Y();
		long int gBin= img->FindBin(x,y);	
		if(gBin<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Cannot find image bin corresponding to peak("<<x<<","<<y<<")");
			#endif
			continue;
		}
		long int ix= img->GetBinX(gBin);
		long int iy= img->GetBinY(gBin);
		double S= img->GetPixelValue(ix,iy);
		
		
		peakPoints.push_back(ImgPeak(x,y,S,ix,iy));
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Peak no. "<<npeaks+1<<" C("<<peakMean.X()<<","<<peakMean.Y()<<")");
		#endif
		npeaks++;	

	}//end loop clusters

	if(npeaks<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No matching peaks across the three dilate kernels detected!");
		#endif
		return 0;
	}
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<npeaks<<" peaks detected!");
	#endif

	return 0;

}//close FindPeaks()


Image* MorphFilter::Dilate(std::vector<long int>& peakPixelIds,Image* img,int KernSize,bool skipBorders)
{
	//## Check input image		
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return nullptr;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	peakPixelIds.clear();

	//## Convert image to OpenCV format
	cv::Mat mat= img->GetOpenCVMat("64");

	//## Init dilation options
	cv::Size kernel_size(KernSize,KernSize);
	cv::Mat element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(-1,-1));
	
	//## Dilate image
	cv::Mat mat_dilated;
	int iterations= 1;
	cv::dilate(mat, mat_dilated, element, cv::Point(-1,-1),iterations,cv::BORDER_CONSTANT);

	//## Compare original and dilated image
	cv::Mat mat_cmp = cv::Mat::zeros(Ny,Nx,CV_8UC1);
	cv::compare(mat, mat_dilated, mat_cmp, cv::CMP_EQ);
	
	//## Convert back dilated image 
	Image* DilatedImg= img->GetCloned("",true,true);
	DilatedImg->Reset();
	
	int npeaks= 0;
	
	for(long int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		
		for(long int i=0;i<Nx;i++){
			int colId= i;
			double binContent= img->GetPixelValue(i,j);
			double matrixElement= mat_dilated.at<double>(rowId,colId);
				
			DilatedImg->FillPixel(i,j,matrixElement);

			float mat_comparison= (float)mat_cmp.at<unsigned char>(rowId,colId);
			bool borderCheck= (!skipBorders || (skipBorders && i>0 && j>0 && i<Nx-1 && j<Ny-1) );
			bool peakFound= (mat_comparison==255 && borderCheck && std::isnormal(binContent));

			if(peakFound){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Found local maximum pixel ("<<i<<","<<j<<"), check surrounding pixel to confirm...");
				#endif

				//## Check surrounding pixels (do not tag as peak if the surrounding 3x3 is flat)
				bool isFlatArea= true;
				for(int ix=i-1;ix<i+1;ix++){
					for(int iy=j-1;iy<j+1;iy++){
						if(ix!=i && iy!=j && ix>=0 && ix<Nx && iy>=0 && iy<Ny){
							double w= img->GetPixelValue(ix,iy);
							if(w!=binContent) {
								isFlatArea= false;
								break;
							}
						}
					}//end loop kernel y
				}//end loop kernel x

				if(!isFlatArea){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Peaks #"<<npeaks<<" detected @ ("<<i<<","<<j<<")");
					#endif
					long int gBin= img->GetBin(i,j);
					peakPixelIds.push_back(gBin);
					npeaks++;			
				}
				
			}//close if check local maximum
		}//end loop x
	}//end loop y
	
	return DilatedImg;

}//close Dilate()



Image* MorphFilter::GetFiltered(std::vector<long int>& peakPixelIds,Image* img,int KernSize,int morphOp,int structElementType,int niters,bool skipBorders) 
{
	//## Check input image		
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return nullptr;
	}
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	//## Convert image to OpenCV format
	cv::Mat mat= img->GetOpenCVMat("64");

	//## Init struct element
	int ptSize= -1;//(KernSize-1)/2;
	cv::Size kernel_size(KernSize,KernSize);
	cv::Mat element;
	if(structElementType==eMORPH_RECT) element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(ptSize,ptSize));
	else if(structElementType==eMORPH_ELLIPSE) element= cv::getStructuringElement(cv::MORPH_ELLIPSE, kernel_size, cv::Point(ptSize,ptSize));
	else if(structElementType==eMORPH_CROSS) element= cv::getStructuringElement(cv::MORPH_CROSS, kernel_size, cv::Point(ptSize,ptSize));	
	else element= cv::getStructuringElement(cv::MORPH_RECT, kernel_size, cv::Point(ptSize,ptSize));

	//## Morphology operation
	//MORPH_OPEN - an opening operation
	//MORPH_CLOSE - a closing operation
	//MORPH_GRADIENT - a morphological gradient
	//MORPH_TOPHAT - “top hat”
	//MORPH_BLACKHAT - “black hat”
	cv::Mat mat_morph;
	if(morphOp==eMORPH_OPENING) cv::morphologyEx(mat, mat_morph, cv::MORPH_OPEN, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_CLOSING) cv::morphologyEx(mat, mat_morph, cv::MORPH_CLOSE, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_GRADIENT) cv::morphologyEx(mat, mat_morph, cv::MORPH_GRADIENT, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_TOPHAT) cv::morphologyEx(mat, mat_morph, cv::MORPH_TOPHAT, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_BLACKHAT) cv::morphologyEx(mat, mat_morph, cv::MORPH_BLACKHAT, element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);	
	else if(morphOp==eMORPH_EROSION) cv::erode(mat,mat_morph,element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);
	else if(morphOp==eMORPH_DILATION) cv::dilate(mat,mat_morph,element,cv::Point(-1,-1),niters,cv::BORDER_CONSTANT);
	else{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Invalid morph operation given ("<<morphOp<<")!");
		#endif
		return nullptr;
	}

	//## Compare original and filtered image
	cv::Mat mat_cmp = cv::Mat::zeros(Ny,Nx,CV_8UC1);
	cv::compare(mat, mat_morph, mat_cmp, cv::CMP_EQ);

	//## Fill morhological filtered image
	Image* MorphImg= img->GetCloned("",true,true);
	MorphImg->Reset();
	
	int npeaks= 0;
	
	for(long int j=0;j<Ny;j++){
		int rowId= Ny-1-j;
		
		for(long int i=0;i<Nx;i++){
			int colId= i;
			double binContent= img->GetPixelValue(i,j);
			double matrixElement= mat_morph.at<double>(rowId,colId);
			//MorphImg->FillPixel(i,j,matrixElement);
				
			float mat_comparison= (float)mat_cmp.at<unsigned char>(rowId,colId);
			MorphImg->FillPixel(i,j,mat_comparison);//DEBUG

			bool borderCheck= (!skipBorders || (skipBorders && i>0 && j>0 && i<Nx-1 && j<Ny-1) );

			if( mat_comparison==255 && borderCheck && std::isnormal(binContent) ){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Found local maximum pixel ("<<i<<","<<j<<"), check surrounding pixel to confirm...");
				#endif

				//## Check surrounding pixels (do not tag as peak if the surrounding 3x3 is flat)
				bool isFlatArea= true;
				for(int ix=i-1;ix<i+1;ix++){
					for(int iy=j-1;iy<j+1;iy++){
						if(ix!=i && iy!=j && ix>=0 && ix<Nx && iy>=0 && iy<Ny){
							double w= img->GetPixelValue(ix,iy);
							if(w!=binContent) {
								isFlatArea= false;
								break;
							}
						}
					}//end loop kernel y
				}//end loop kernel x

				if(!isFlatArea){
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Peaks #"<<npeaks<<" detected @ ("<<i<<","<<j<<")");
					#endif
					long int gBin= img->GetBin(i,j);
					peakPixelIds.push_back(gBin);
					npeaks++;			
				}
				else{
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Local maximum pixel found ("<<i<<","<<j<<") not confirmed (flat surrounding)...");
					#endif
				}
			}//close if check local maximum
		}//end loop x
	}//end loop y
	
	return MorphImg;

}//close GetFiltered()


int MorphFilter::FindDilatedSourcePixels(Image* img,Source* source,int KernSize,std::vector<long int>& pixelsToBeDilated){

	//## Check input source
	if(!source) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given!");
		#endif
		return -1;
	}
	
	//## Init dilation kernel
	if(KernSize%2==0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("KernSize argument should be an odd number!");
		#endif
		return -1;
	}
	int dilateSize= KernSize/2;
	cv::Mat element= cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(KernSize,KernSize));
	
	//PixelCollection sourcePixels= source->GetPixels();
	int nDilatedPixels= 0;
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	for(int l=0;l<source->GetNPixels();l++){
		Pixel* thisPixel= source->GetPixel(l);
		long int id= thisPixel->id;
		long int binx= img->GetBinX(id);
		long int biny= img->GetBinY(id);
		long int ix= binx;
		long int iy= biny;			

		//Find if pixel was already added to the list of dilated pixels
		std::vector<long int>::iterator it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),id);
		if( (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) pixelsToBeDilated.push_back(id);
			
		for(int tx=-dilateSize;tx<=dilateSize;tx++){
			long int binx_next= binx + tx;
			long int colId= tx + dilateSize;

			for(int ty=-dilateSize;ty<=dilateSize;ty++){	
				long int biny_next= biny+ty;
				long int rowId= tx + dilateSize;

				if(ix+tx<Nx && ix+tx>=0 && iy+ty<Ny && iy+ty>=0){
					double kernValue= (double)element.at<char>(rowId,colId);
					long int gBinId= img->GetBin(binx_next,biny_next);
					it= std::find(pixelsToBeDilated.begin(),pixelsToBeDilated.end(),gBinId);
					if( kernValue>0 && (it==pixelsToBeDilated.end() || pixelsToBeDilated.empty()) ) {
						pixelsToBeDilated.push_back(gBinId);
						nDilatedPixels++;
					}
				}//close if
			}//end loop kernel
		}//end loop kernel
	}//end loop pixels		
	
	//## Replace selected pixels
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<nDilatedPixels<<" pixels to be dilated...");
	#endif

	return 0;

}//close FindDilatedSourcePixels()


int MorphFilter::FindDilatedSourcePixels(std::vector<long int>& pixelsToBeDilated,Image* img,Source* source,int kernSize)
{
	//## Check input source
	if(!source) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given!");
		#endif
		return -1;
	}
	
	//## Init dilation kernel
	if(kernSize%2==0){
		kernSize++;
		#ifdef LOGGING_ENABLED
			WARN_LOG("kern size argument given was converted to an odd number ("<<kernSize<<") ...");
		#endif
	}
	int dilateSize= kernSize/2;
	cv::Mat element= cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(kernSize,kernSize));
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	//## Add source pixels to dilated pixel list
	std::set<long int> dilatedPixelIds;
	for(long int l=0;l<source->GetNPixels();l++){
		Pixel* thisPixel= source->GetPixel(l);
		long int id= thisPixel->id;
		dilatedPixelIds.insert(id);
	}

	//## Get source contours
	std::vector<Contour*> contours= source->GetContours();
	if(contours.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No contours available for source "<<source->GetName()<<" (hint: check if contours were computed)!");
		#endif
	}

	for(size_t k=0;k<contours.size();k++){
		Contour* contour= contours[k];
		if(!contour) continue;
		double contx= 0;
		double conty= 0;
		
		for(int i=0;i<contour->GetN();i++){
			contour->GetPointXY(contx,conty,i);
			long int id= img->FindBin(contx,conty);
			if(id<0) continue;
			long int binx= img->GetBinX(id);
			long int biny= img->GetBinY(id);

			for(int tx=-dilateSize;tx<=dilateSize;tx++){
				long int binx_next= binx + tx;
				long int colId= tx + dilateSize;

				for(int ty=-dilateSize;ty<=dilateSize;ty++){	
					long int biny_next= biny+ty;
					long int rowId= tx + dilateSize;
					long int gBinId= img->GetBin(binx_next,biny_next);
					if(gBinId<0) continue;
					double kernValue= (double)element.at<char>(rowId,colId);
					if(kernValue>0) dilatedPixelIds.insert(gBinId);	
				}//end loop kernel
			}//end loop kernel
		}//end loop contour points	
	}//end loop contours

	//Add additional pixel to list
	std::copy(dilatedPixelIds.begin(), dilatedPixelIds.end(), std::back_inserter(pixelsToBeDilated));
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<pixelsToBeDilated.size()<<" pixels to be dilated...");
	#endif

	return 0;

}//close FindDilatedSourcePixels()


int MorphFilter::DilateAroundSource(Image* img,Source* source,int KernSize,int dilateModel,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,Image* mask,int bkgBoxThickness)
{
	//## Check input source
	if(!source) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given!");
		#endif
		return -1;
	}

	//## Check if source has stats computed, if not compute
	bool hasStats= source->HasStats();
	if(!hasStats){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No stats computed for input source...computing!");
		#endif
		source->ComputeStats(true,true);
	}
	double sourceMedian= source->Median;
	double sourceMedianRMS= source->MedianRMS;

	//## Compute background around source (used if no bkgdata are given)
	BkgSampleData bkgBoxData;
	int bkgEstimator= eMedianBkg;
	bool useParallelVersion= false;
	std::vector<float> maskedValues= {};
	int status= img->GetBkgInfoAroundSource(
		bkgBoxData,
		source,
		bkgBoxThickness,
		bkgEstimator,
		mask,
		useParallelVersion,
		maskedValues
	);
	if(status<0 && !bkgData){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to compute bkg around source (needed since no bkgdata were given in input)!");
		#endif
		return -1;
	}
	double boxBkg= bkgBoxData.bkgLevel;
	double boxBkgRMS= bkgBoxData.bkgRMS;
	
	//## Initialize GSL random init
	#ifdef LOGGING_ENABLED	
		DEBUG_LOG("Initialize GSL random engine...");
	#endif

  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng* rand_generator = gsl_rng_alloc (type); 

	//## Find pixels to be dilated
	std::vector<long int> pixelsToBeDilated;
	FindDilatedSourcePixels(pixelsToBeDilated,img,source,KernSize);

	//## Replace dilated pixels with model
	double sigmaTrunc= 1;//trunc random gaussian to +-sigmaTrunc	
	
	if(dilateModel==eDilateWithSourceMedian){
		double BkgRealization= sourceMedian;
		double BkgRMS= sourceMedianRMS;
		if(randomize){
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t l=0;l<pixelsToBeDilated.size();l++){
				long int id= pixelsToBeDilated[l];	
				double r= 0;
				RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
				double bkg= BkgRealization + r*BkgRMS;
				img->SetPixelValue(id,bkg);
			}//end loop pixels 	
		}
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t l=0;l<pixelsToBeDilated.size();l++){
				long int id= pixelsToBeDilated[l];	
				img->SetPixelValue(id,BkgRealization);
			}//end loop pixels 
		}
	}//close if
  else if(dilateModel==eDilateWithBkg)
	{
		if(bkgData)
		{

			if(useLocalBkg){
				if(randomize){
					#ifdef OPENMP_ENABLED
					#pragma omp parallel for
					#endif
					for(size_t l=0;l<pixelsToBeDilated.size();l++){
						long int id= pixelsToBeDilated[l];	
						double BkgRealization= (bkgData->BkgMap)->GetPixelValue(id);
						double BkgRMS= (bkgData->NoiseMap)->GetPixelValue(id);
						double r= 0;
						RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
						double bkg= BkgRealization + r*BkgRMS;
						img->SetPixelValue(id,bkg);
					}//end loop pixels
				}
				else{
					#ifdef OPENMP_ENABLED
					#pragma omp parallel for
					#endif
					for(size_t l=0;l<pixelsToBeDilated.size();l++){
						long int id= pixelsToBeDilated[l];	
						double BkgRealization= (bkgData->BkgMap)->GetPixelValue(id);
						img->SetPixelValue(id,BkgRealization);
					}//end loop pixels
				}
			}//close if
			else{
				double BkgRealization= bkgData->gBkg;
				double BkgRMS= bkgData->gNoise;	
				if(randomize){
					#ifdef OPENMP_ENABLED
					#pragma omp parallel for
					#endif
					for(size_t l=0;l<pixelsToBeDilated.size();l++){
						long int id= pixelsToBeDilated[l];	
						double r= 0;
						RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
						double bkg= BkgRealization + r*BkgRMS;
						img->SetPixelValue(id,bkg);
					}//end loop pixels
				}
				else{
					#ifdef OPENMP_ENABLED
					#pragma omp parallel for
					#endif
					for(size_t l=0;l<pixelsToBeDilated.size();l++){
						long int id= pixelsToBeDilated[l];	
						img->SetPixelValue(id,BkgRealization);
					}//end loop pixels
				}
			}//close else
		}//close else if has bkgData
		else{//use box bkg
			double BkgRealization= boxBkg;
			double BkgRMS= boxBkgRMS;	
				
			if(randomize){				
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];	
					double r= 0;
					RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
					double bkg= BkgRealization + r*BkgRMS;
					img->SetPixelValue(id,bkg);
				}//end loop pixels
			}//close if randomize
			else{
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];	
					img->SetPixelValue(id,BkgRealization);
				}//end loop pixels
			}//close else !randomize
		}//close else

	}//close else if

	// GSL rand generator deallocation
	gsl_rng_free(rand_generator); 

	return 0;

}//close DilateAroundSource()




int MorphFilter::DilateAroundSource(Image* img,Source* source,int KernSize,int dilateModel,int dilateSourceMorphId,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr,double zBrightThr)
{
	//## Check input source
	if(!source) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given!");
		#endif
		return -1;
	}
	bool hasNestedSources= source->HasNestedSources();
	bool hasStats= source->HasStats();
	if(!hasStats){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No stats computed for input source...computing!");
		#endif
		source->ComputeStats(true,true);
	}
	int sourceMorphId= source->MorphId;
	double sourceMedian= source->Median;
	double sourceMedianRMS= source->MedianRMS;
	
	//## Skip faint sources
	//Get pixel seeds
	bool isBrightSource= false;
	std::vector<int> seedPixelIndexes= source->GetSeedPixelIndexes();
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<seedPixelIndexes.size()<<" seed pixels...");
	#endif

	double Zmax= -1.e+99;

	#ifdef OPENMP_ENABLED
	#pragma omp parallel for reduction(max: Zmax)
	#endif
	for(size_t i=0;i<seedPixelIndexes.size();i++){
		int index= seedPixelIndexes[i];
		Pixel* aPixel= source->GetPixel(index);
		if(!aPixel) continue;
		long int id= aPixel->id;
		std::pair<double,double> bkgInfo= aPixel->GetBkg();
		double bkgLevel= bkgInfo.first;
		double noiseLevel= bkgInfo.second;
		double w= img->GetBinContent(id);
		double Z= (w-bkgLevel)/noiseLevel;
		if(fabs(Z)>Zmax) Zmax= Z;
	}

	if(Zmax<zThr){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Source is below significance threshold for dilation (Z="<<Zmax<<"<"<<zThr<<"), skip it!");
		#endif
		return 0;	
	}
	isBrightSource= (Zmax>=zBrightThr);
	
	//## Check R interface
	double sigmaTrunc= 1;//trunc random gaussian to +-sigmaTrunc	

	/*
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Retrieve RInterface instance...");
	#endif

	ROOT::R::TRInterface& fR= ROOT::R::TRInterface::Instance();
	std::string randomGenCmd= std::string("rtruncnorm(1, a=-sigmaTrunc, b=sigmaTrunc, mean = 0, sd = 1)");
	try{
		fR.Execute("library(\"truncnorm\");");
		fR["sigmaTrunc"]= sigmaTrunc;
	}
	catch( std::exception &ex ) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("C++ exception catched while loading R library truncnorm (err=" << ex.what() <<")");
		#endif
		return -1;
  } 
	catch(...) { 
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Unknown exception catched while loading R library truncnorm!");
		#endif
		return -1;
  }	
	*/

	//## Initialize GSL random init
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Initialize GSL random engine...");
	#endif
  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng* rand_generator = gsl_rng_alloc (type); 


	//## Find pixels to be dilated
	//## NB: If Z>Zthr_bright ==> dilate mother source always no matter what source morph, ignore all nested
	//##     If Zth<Z<Zthr_bright ==> dilate mother source is there are no nested and morph match, otherwise dilate nested

	std::vector<long int> pixelsToBeDilated;
	if(isBrightSource){//MOTHER SOURCE
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Selecting entire mother source for dilation (bright source) ...");
		#endif
		FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);
	}
	else{
		if(hasNestedSources && skipToNested){//NESTED SOURCES
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Dilating nested sources...");
			#endif
			std::vector<Source*> nestedSources= source->GetNestedSources();
			for(size_t k=0;k<nestedSources.size();k++){
				int nestedSourceMorphId= nestedSources[k]->MorphId;
				//if(dilateSourceMorphId==-1 || nestedSourceMorphId==sourceMorphId){	//possible bug!
				if(dilateSourceMorphId==-1 || nestedSourceMorphId==dilateSourceMorphId){
					FindDilatedSourcePixels(img,nestedSources[k],KernSize,pixelsToBeDilated);
				}
			}//end loop nested sources
		}//close if
		else{//MOTHER SOURCE
			if(dilateSourceMorphId==-1 || sourceMorphId==dilateSourceMorphId){
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Dilating entire mother source (no nested sources present + match source morph id) ...");
				#endif
				FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);
			}
		}
	}//close else


	//## Replace dilated pixels with model		
	if(dilateModel==eDilateWithSourceMedian){
		double BkgRealization= sourceMedian;
		double BkgRMS= sourceMedianRMS;
		if(randomize){
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t l=0;l<pixelsToBeDilated.size();l++){
				long int id= pixelsToBeDilated[l];			
				//double r= fR.Eval(randomGenCmd.c_str());
				double r= 0;
				RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
				double bkg= BkgRealization + r*BkgRMS;
				img->SetPixelValue(id,bkg);
			}//end loop pixels 	
		}
		else{
			#ifdef OPENMP_ENABLED
			#pragma omp parallel for
			#endif
			for(size_t l=0;l<pixelsToBeDilated.size();l++){
				long int id= pixelsToBeDilated[l];			
				img->SetPixelValue(id,BkgRealization);
			}//end loop pixels 
		}
	}//close if
  else if(dilateModel==eDilateWithBkg){
		if(useLocalBkg){
			if(randomize){
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					double BkgRealization= (bkgData->BkgMap)->GetPixelValue(id);
					double BkgRMS= (bkgData->NoiseMap)->GetPixelValue(id);
					double r= 0;
					RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
					//double r= fR.Eval(randomGenCmd.c_str());
					double bkg= BkgRealization + r*BkgRMS;
					img->SetPixelValue(id,bkg);
				}//end loop pixels
			}
			else{
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					double BkgRealization= (bkgData->BkgMap)->GetPixelValue(id);
					img->SetPixelValue(id,BkgRealization);
				}//end loop pixels
			}
		}//close if
		else{
			double BkgRealization= bkgData->gBkg;
			double BkgRMS= bkgData->gNoise;	
			if(randomize){
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					//double r= fR.Eval(randomGenCmd.c_str());
					double r= 0;
					RtNorm_ns::RtNorm::get_random(r,rand_generator,-sigmaTrunc,sigmaTrunc,0.,1.);
					double bkg= BkgRealization + r*BkgRMS;
					img->SetPixelValue(id,bkg);
				}//end loop pixels
			}
			else{
				#ifdef OPENMP_ENABLED
				#pragma omp parallel for
				#endif
				for(size_t l=0;l<pixelsToBeDilated.size();l++){
					long int id= pixelsToBeDilated[l];			
					img->SetPixelValue(id,BkgRealization);
				}//end loop pixels
			}
		}//close else
	}//close else if

	//## GSL rand generator deallocation
	gsl_rng_free(rand_generator); 

	//Force recomputation of stats if present, otherwise recompute only moments
	//bool skipNegativePixels= false;
	bool useRange= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	//if(img->HasStats()) status= img->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	//else status= img->ComputeMoments(skipNegativePixels);
	if(img->HasStats()) status= img->ComputeStats(computeRobustStats,forceRecomputing,useRange);
	else status= img->ComputeMoments(useRange);
		
	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to recompute moments/stats after source dilation!");
		#endif
		return -1;
	}

	return 0;

}//close DilateAroundSource()


int MorphFilter::DilateAroundSources(Image* img,std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceMorphId,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr,double zBrightThr)
{	
	//## Check input image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return -1;
	}
	
	//## Check bkg data
	if(dilateModel==eDilateWithBkg){
	 	if(!bkgData){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Selected to use bkg dilation but null ptr to bkg data!");
			#endif
			return -1;
		}
		if(useLocalBkg && !bkgData->HasLocalBkg()){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Selected to use local bkg but no local bkg data are available!");
			#endif
			return -1;
		}
	}//close if

	//## Check source list
	if(sources.size()<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source list empty, nothing to be dilated!");
		#endif
		return 0;
	}

	//## Start dilating sources
	for(size_t k=0;k<sources.size();k++){	
		int status= DilateAroundSource(img,sources[k],KernSize,dilateModel,dilateSourceMorphId,skipToNested,bkgData,useLocalBkg,randomize,zThr,zBrightThr);
		if(status<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Source dilation failed for source no. "<<k<<" ...");
			#endif
		}
	}//end loop sources

	return 0;

}//close DilateAroundSources()


}//close namespace
