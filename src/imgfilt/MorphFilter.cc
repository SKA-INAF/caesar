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
#include <Logger.h>
#include <Graph.h>

#include <rtnorm.h>

#include <TObject.h>
#include <TVector2.h>
//#include <TRInterface.h>

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


int MorphFilter::FindPeaks(std::vector<TVector2>& peakPoints,Image* img,std::vector<int> kernelSizes,int peakShiftTolerance,bool skipBorders,int peakKernelMultiplicityThr)
{
	//Check input image
	if(!img){
		ERROR_LOG("Null ptr to input image given!");
		return -1;
	}

	//Check multiplicity thr
	int nKernels= static_cast<int>(kernelSizes.size());
	if(peakKernelMultiplicityThr==-1){
		peakKernelMultiplicityThr= nKernels;
	}
	if(peakKernelMultiplicityThr>nKernels){
		WARN_LOG("Multiplicy thr ("<<peakKernelMultiplicityThr<<") larger than nkernels ("<<nKernels<<"), set it to "<<nKernels<<"...");
		peakKernelMultiplicityThr= nKernels;
	}
	if(peakKernelMultiplicityThr==0 || (peakKernelMultiplicityThr<0 && peakKernelMultiplicityThr!=-1)){
		ERROR_LOG("Invalid multiplicity thr ("<<peakKernelMultiplicityThr<<") given!");
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
			ERROR_LOG("Failed to compute dilated image for kernel no. "<<k<<" (size="<<kernelSizes[k]<<")!");
			return -1;
		}
		delete dilatedImg;
		dilatedImg= 0;		
		
		//Stop if no peaks are detected
		if(peakPixelIds.empty()){
			hasPeaks= false;
			break;
		}
		INFO_LOG("#"<<peakPixelIds.size()<<" peak pixels found with kernel no. "<<k<<" (size="<<kernelSizes[k]<<")!");

		//Store peak points for current kernel
		for(size_t j=0;j<peakPixelIds.size();j++){
			long int gBin= peakPixelIds[j];
			long int binX= img->GetBinX(gBin); 
			long int binY= img->GetBinY(gBin); 
			double x= img->GetX(binX);
			double y= img->GetY(binY);
			//points[k].push_back(TVector2(x,y));
			peaks.push_back(PeakData(x,y,k));
		}//end loop peak points

	}//end loop kernels

	if(!hasPeaks) {
		DEBUG_LOG("No peaks detected in one or all dilated kernel runs.");
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
			DEBUG_LOG("Skip peak group no. "<<i+1<<" as below multiplicity thr ("<<multiplicity<<"<"<<peakKernelMultiplicityThr<<")");
			continue;
		}

		//Compute peak mean and fill final peaks list
		TVector2 peakMean(0,0);
		for(size_t j=0;j<connected_indexes[i].size();j++){
			int index= connected_indexes[i][j];
			peakMean+= TVector2(peaks[index].x,peaks[index].y);
		}
		peakMean*= 1./(float)(connected_indexes[i].size());
		peakPoints.push_back(peakMean);
		INFO_LOG("Peak no. "<<npeaks+1<<" C("<<peakMean.X()<<","<<peakMean.Y()<<")");
		npeaks++;	

	}//end loop clusters

	if(npeaks<=0){
		WARN_LOG("No matching peaks across the three dilate kernels detected!");
		return 0;
	}
	INFO_LOG("#"<<npeaks<<" peaks detected!");

	return 0;

}//close FindPeaks()




Image* MorphFilter::Dilate(std::vector<long int>& peakPixelIds,Image* img,int KernSize,bool skipBorders)
{
	//## Check input image		
	if(!img){
		ERROR_LOG("Null ptr to given image!");
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

			if( mat_comparison==255 && borderCheck && std::isnormal(binContent) ){
				//INFO_LOG("Found local maximum pixel ("<<i<<","<<j<<"), check surrounding pixel to confirm...");

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
					INFO_LOG("Peaks #"<<npeaks<<" detected @ ("<<i<<","<<j<<")");
					long int gBin= img->GetBin(i,j);
					peakPixelIds.push_back(gBin);
					npeaks++;			
				}
				//else{
				//	INFO_LOG("Local maximum pixel found ("<<i<<","<<j<<") not confirmed (flat surrounding)...");
				//}
			}//close if check local maximum
		}//end loop x
	}//end loop y
	
	return DilatedImg;

}//close Dilate()



Image* MorphFilter::GetFiltered(std::vector<long int>& peakPixelIds,Image* img,int KernSize,int morphOp,int structElementType,int niters,bool skipBorders) 
{
	//## Check input image		
	if(!img){
		ERROR_LOG("Null ptr to given image!");
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
		ERROR_LOG("Invalid morph operation given ("<<morphOp<<")!");
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
				INFO_LOG("Found local maximum pixel ("<<i<<","<<j<<"), check surrounding pixel to confirm...");

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
					INFO_LOG("Peaks #"<<npeaks<<" detected @ ("<<i<<","<<j<<")");
					long int gBin= img->GetBin(i,j);
					peakPixelIds.push_back(gBin);
					npeaks++;			
				}
				else{
					INFO_LOG("Local maximum pixel found ("<<i<<","<<j<<") not confirmed (flat surrounding)...");
				}
			}//close if check local maximum
		}//end loop x
	}//end loop y
	
	return MorphImg;

}//close GetFiltered()


int MorphFilter::FindDilatedSourcePixels(Image* img,Source* source,int KernSize,std::vector<long int>& pixelsToBeDilated){

	//## Check input source
	if(!source) {
		ERROR_LOG("Null ptr to input source given!");
		return -1;
	}
	
	//## Init dilation kernel
	if(KernSize%2==0){
		ERROR_LOG("KernSize argument should be an odd number!");
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
	DEBUG_LOG("#"<<nDilatedPixels<<" pixels to be dilated...");
		
	return 0;

}//close FindDilatedSourcePixels()


int MorphFilter::DilateAroundSource(Image* img,Source* source,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr,double zBrightThr){

	//## Check input source
	if(!source) {
		ERROR_LOG("Null ptr to input source given!");
		return -1;
	}
	bool hasNestedSources= source->HasNestedSources();
	bool hasStats= source->HasStats();
	if(!hasStats){
		WARN_LOG("No stats computed for input source...computing!");
		source->ComputeStats(true,true);
	}
	int sourceType= source->Type;
	double sourceMedian= source->Median;
	double sourceMedianRMS= source->MedianRMS;
	
	//## Skip faint sources
	//Get pixel seeds
	bool isBrightSource= false;
	std::vector<int> seedPixelIndexes= source->GetSeedPixelIndexes();
	DEBUG_LOG("#"<<seedPixelIndexes.size()<<" seed pixels...");

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
		DEBUG_LOG("Source is below significance threshold for dilation (Z="<<Zmax<<"<"<<zThr<<"), skip it!");
		return 0;	
	}
	isBrightSource= (Zmax>=zBrightThr);
	
	//## Check R interface
	double sigmaTrunc= 1;//trunc random gaussian to +-sigmaTrunc	

	/*
	DEBUG_LOG("Retrieve RInterface instance...");
	ROOT::R::TRInterface& fR= ROOT::R::TRInterface::Instance();
	std::string randomGenCmd= std::string("rtruncnorm(1, a=-sigmaTrunc, b=sigmaTrunc, mean = 0, sd = 1)");
	try{
		fR.Execute("library(\"truncnorm\");");
		fR["sigmaTrunc"]= sigmaTrunc;
	}
	catch( std::exception &ex ) {
		ERROR_LOG("C++ exception catched while loading R library truncnorm (err=" << ex.what() <<")");
		return -1;
  } 
	catch(...) { 
		ERROR_LOG("Unknown exception catched while loading R library truncnorm!");
		return -1;
  }	
	*/

	//## Initialize GSL random init
	DEBUG_LOG("Initialize GSL random engine...");
  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng* rand_generator = gsl_rng_alloc (type); 


	//## Find pixels to be dilated
	//## NB: If Z>Zthr_bright ==> dilate mother source always no matter what source type, ignore all nested
	//##     If Zth<Z<Zthr_bright ==> dilate mother source is there are no nested and type match, otherwise dilate nested

	std::vector<long int> pixelsToBeDilated;
	if(isBrightSource){//MOTHER SOURCE
		DEBUG_LOG("Selecting entire mother source for dilation (bright source) ...");
		FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);
	}
	else{
		if(hasNestedSources && skipToNested){//NESTED SOURCES
			DEBUG_LOG("Dilating nested sources...");
			std::vector<Source*> nestedSources= source->GetNestedSources();
			for(size_t k=0;k<nestedSources.size();k++){
				int nestedSourceType= nestedSources[k]->Type;
				if(dilateSourceType==-1 || nestedSourceType==sourceType){
					FindDilatedSourcePixels(img,nestedSources[k],KernSize,pixelsToBeDilated);
				}
			}//end loop nested sources
		}//close if
		else{//MOTHER SOURCE
			if(dilateSourceType==-1 || sourceType==dilateSourceType){
				DEBUG_LOG("Dilating entire mother source (no nested sources present + match source type) ...");
				FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);
			}
		}
	}//close else



	/*
	//== MOTHER SOURCE
	std::vector<long int> pixelsToBeDilated;
	if(!hasNestedSources && (dilateSourceType==-1 || sourceType==dilateSourceType) ){
		DEBUG_LOG("Dilating mother sources...");
		FindDilatedSourcePixels(img,source,KernSize,pixelsToBeDilated);	
	}

	//== NESTED SOURCES
	if(skipToNested && hasNestedSources){
		DEBUG_LOG("Dilating nested sources...");
		std::vector<Source*> nestedSources= source->GetNestedSources();
		for(unsigned int k=0;k<nestedSources.size();k++){
			int nestedSourceType= nestedSources[k]->Type;
			if(dilateSourceType==-1 || nestedSourceType==sourceType){
				FindDilatedSourcePixels(img,nestedSources[k],KernSize,pixelsToBeDilated);
			}
		}//end loop nested sources
	}//close if dilate nested
	*/



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
	bool skipNegativePixels= false;
	bool computeRobustStats= true;	
	bool forceRecomputing= true;
	int status= 0;
	if(img->HasStats()) status= img->ComputeStats(computeRobustStats,skipNegativePixels,forceRecomputing);
	else status= img->ComputeMoments(skipNegativePixels);
		
	if(status<0){
		WARN_LOG("Failed to recompute moments/stats after source dilation!");
		return -1;
	}

	return 0;

}//close DilateAroundSource()


int MorphFilter::DilateAroundSources(Image* img,std::vector<Source*>const& sources,int KernSize,int dilateModel,int dilateSourceType,bool skipToNested,ImgBkgData* bkgData,bool useLocalBkg,bool randomize,double zThr,double zBrightThr){
	
	//## Check input image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}
	
	//## Check bkg data
	if(dilateModel==eDilateWithBkg){
	 	if(!bkgData){
			ERROR_LOG("Selected to use bkg dilation but null ptr to bkg data!");
			return -1;
		}
		if(useLocalBkg && !bkgData->HasLocalBkg()){
			ERROR_LOG("Selected to use local bkg but no local bkg data are available!");
			return -1;
		}
	}//close if

	//## Check source list
	if(sources.size()<=0){
		WARN_LOG("Source list empty, nothing to be dilated!");
		return 0;
	}

	//## Start dilating sources
	for(size_t k=0;k<sources.size();k++){	
		int status= DilateAroundSource(img,sources[k],KernSize,dilateModel,dilateSourceType,skipToNested,bkgData,useLocalBkg,randomize,zThr,zBrightThr);
		if(status<0){
			WARN_LOG("Source dilation failed for source no. "<<k<<" ...");
		}
	}//end loop sources

	return 0;

}//close DilateAroundSources()


}//close namespace
