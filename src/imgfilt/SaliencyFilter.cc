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
* @file SaliencyFilter.cc
* @class SaliencyFilter
* @brief Class implementing saliency filtering
*
* Saliency Filter
* @author S. Riggi
* @date 20/01/2015
*/

#include <SaliencyFilter.h>
#include <Image.h>
#include <Region.h>
#include <SLIC.h>
#include <SLICData.h>
#include <BkgData.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TVectorD.h>
#include <TGraph2D.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>

using namespace std;

ClassImp(Caesar::SaliencyFilter)

namespace Caesar {

SaliencyFilter::SaliencyFilter() {

}//close costructor


SaliencyFilter::~SaliencyFilter(){

}//close destructor

//===================================================
//==        NEW IMAGE METHODS
//===================================================
Image* SaliencyFilter::ComputeSaliencyMap(Image* img,int reso,double regFactor,int minRegionSize,double knnFactor,bool useRobust,double expFalloffPar,double distanceRegPar)
{
	//## Normalize image
	double NormMin= 1;
	double NormMax= 256;
	INFO_LOG("Normalize image to range ["<<NormMin<<","<<NormMax<<"]");
	Image* img_norm= img->GetNormalizedImage("LINEAR",NormMin,NormMax);
	if(!img_norm){
		ERROR_LOG("Failed to normalize input image!");
		return 0;
	}

	//## Compute segmentation in superpixels
	INFO_LOG("Generate superpixel partition @ reso="<<reso<<" (beta="<<regFactor<<", minRegionSize="<<minRegionSize<<")...");
	bool normalizeImg= true;
	bool useLogScaleMapping= false;
	SLICData* slicData= SLIC::SPGenerator(img_norm,reso,regFactor,minRegionSize,normalizeImg,useLogScaleMapping,0,0);	//pass null laplImg & edgeImg (not needed here)
	if(!slicData){
		ERROR_LOG("Superpixel segmentation failed!");
		return 0;
	}

	//Get results
	INFO_LOG("Getting access to generated superpixel list...");
	std::vector<Region*> regions= slicData->GetRegions();
	if(regions.empty()) {
		ERROR_LOG("Superpixel segmentation returned no regions!");
		delete slicData;
		slicData= 0;
		return 0;
	}

	//## Compute saliency map		
	INFO_LOG("Compute saliency map @ reso="<<reso<<" (knnFactor="<<knnFactor<<", useRobust="<<useRobust<<", expFalloffPar="<<expFalloffPar<<", distanceRegPar="<<distanceRegPar<<")...");
	Image* saliencyMap= ComputeSaliencyMap(img_norm,regions,knnFactor,useRobust,expFalloffPar,distanceRegPar);
	if(!saliencyMap){
		ERROR_LOG("Failed to compute saliency map!");
		delete slicData;
		slicData= 0;
		return 0;
	}

	//## Clear memory
	INFO_LOG("Clear allocated memory for saliency computation @ reso "<<reso);
	img_norm->Delete();
	delete slicData;
	slicData= 0;
	
	return saliencyMap;

}//close ComputeSaliencyMap()

Image* SaliencyFilter::ComputeSaliencyMap(Image* img,std::vector<Region*>const& regions,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar){

	//## Check regions
	size_t nRegions= regions.size();
	if(nRegions<=0) {
		ERROR_LOG("No regions given!");
		return 0;
	}

	//## Set number of nearest neighbor regions to be used in saliency computation
	long long int knn_min= 10;//in number of regions
	long long int knn_chosen= static_cast<long long int>(std::round(knnFactor*nRegions));//percentage of available regions
	long long int knn= std::max(knn_chosen,knn_min);
	size_t KNN= knn;
	if(knn>(signed)nRegions || knn<0) KNN= nRegions;
	
	//## Get image info
	double Xmin= img->GetXmin();
	double Xmax= img->GetXmax();
	double Ymin= img->GetYmin();
	double Ymax= img->GetYmax();
	double width= fabs(Xmax-Xmin);
	double height= fabs(Ymax-Ymin);
	double diagonal= sqrt(width*width + height*height);

	//## Compute region stats pars
	INFO_LOG("Compute region pars (#"<<nRegions<<" present) ...");
	for(size_t i=0;i<nRegions;i++) {
		if(!regions[i]->HasStats()) regions[i]->ComputeStats(true,false);
	}

	//## Compute dissimilarity matrix
	TMatrixD* ColorDistMatrix= new TMatrixD(nRegions,nRegions);
	ColorDistMatrix->Zero();

	TMatrixD* SpatialDistMatrix= new TMatrixD(nRegions,nRegions);
	SpatialDistMatrix->Zero();
	double dist_c_min= 1.e+99;
	double dist_c_max= -1.e+99;
	double dist_s_min= 1.e+99;
	double dist_s_max= -1.e+99;
	
	for(size_t i=0;i<nRegions-1;i++){

		//Region pars i-th
		double mu_i= regions[i]->Mean;
		double median_i= regions[i]->Median;
		double Xc_i= regions[i]->X0;
		double Yc_i= regions[i]->Y0;
		
		for(int j=i+1;j<nRegions;j++){

			//Region pars j-th
			double mu_j= regions[j]->Mean;
			double median_j= regions[j]->Median;
			double Xc_j= regions[j]->X0;
			double Yc_j= regions[j]->Y0;
					
			//Compute color & spatial distances
			double dist_c= fabs(mu_i-mu_j);
			if(useRobustPars) dist_c= fabs(median_i-median_j);
			double dist_s= sqrt( (Xc_i-Xc_j)*(Xc_i-Xc_j) + (Yc_i-Yc_j)*(Yc_i-Yc_j) ); 
			
			(*ColorDistMatrix)(i,j)= dist_c;
			(*ColorDistMatrix)(j,i)= dist_c;

			(*SpatialDistMatrix)(i,j)= dist_s;
			(*SpatialDistMatrix)(j,i)= dist_s;

			//Find min & max
			if(dist_c<dist_c_min) dist_c_min= dist_c;
			if(dist_c>dist_c_max) dist_c_max= dist_c;
			if(dist_s<dist_s_min) dist_s_min= dist_s;
			if(dist_s>dist_s_max) dist_s_max= dist_s;
	
		}//end loop regions
	}//end loop regions

	//## Normalize distances to [0,1]
	//## Color distances normalized to min & max
	//## Spatial distances normalized to image diagonal
	INFO_LOG("Color dist min/max: "<<dist_c_min<<"/"<<dist_c_max<<", Spatial dist min/max: "<<dist_s_min<<"/"<<dist_s_max<<" img size("<<width<<" x "<<height<<" (diagonal="<<diagonal<<")");
	
	double NormMin= 0;
	double NormMax= 1;

	for(size_t i=0;i<nRegions;i++){
		for(size_t j=i+1;j<nRegions;j++){
			
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_c_norm= NormMin + (NormMax-NormMin)*(dist_c-dist_c_min)/(dist_c_max-dist_c_min);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist_s_norm= NormMin + (NormMax-NormMin)*(dist_s-dist_s_min)/(dist_s_max-dist_s_min);	
			//double dist_s_norm= dist_s/diagonal;

			(*ColorDistMatrix)(i,j)= dist_c_norm;
			(*ColorDistMatrix)(j,i)= dist_c_norm;

			(*SpatialDistMatrix)(i,j)= dist_s_norm;
			(*SpatialDistMatrix)(j,i)= dist_s_norm;
		}//end loop regions
	}//end loop regions
	
	INFO_LOG("Color dist min/max: "<<ColorDistMatrix->Min()<<"/"<<ColorDistMatrix->Max()<<", Spatial dist min/max: "<<SpatialDistMatrix->Min()<<"/"<<SpatialDistMatrix->Max());

	//## Create saliency image
	TString imgName= Form("%s_saliency",img->GetName().c_str());
	Image* saliencyImg= (Image*)img->GetCloned(std::string(imgName),true,true);
	saliencyImg->SetName(std::string(imgName));
	saliencyImg->Reset();

	//## Compute saliency 
	INFO_LOG("Computing saliency map ...");
	double Smin= 1.e+99;
	double Smax= -1.e+99;
	std::vector<double> SList;

	for(int i=0;i<nRegions;i++){

		std::vector<double> dissList;

		for(int j=0;j<nRegions;j++){
			//if(i==j) continue;
			double dist_c= (*ColorDistMatrix)(i,j);
			double dist_s= (*SpatialDistMatrix)(i,j);
			double dist= dist_c/(1 + distanceRegPar*dist_s);
			double dissimilarity= exp(-expFalloffPar*dist);

			dissList.push_back(dissimilarity);
		}//end loop regions

		//Sort color dissimilarities for region i-th to use only K-th neighbors in color
		std::vector<double> sorted;
		std::vector<size_t> sort_index;//sorting index
		CodeUtils::sort(dissList,sorted,sort_index);	

		//Compute saliency over k-th neighbors
		double S= 0;
		for(size_t k=0;k<KNN;k++){
			size_t index= sort_index[k];
			double D= dissList[index];	
			S+= D;
		}
		S/= static_cast<double>(KNN);
		SList.push_back(S);

		
		if(S<Smin) Smin= S;
		if(S>Smax) Smax= S;

	}//end loop regions

	DEBUG_LOG("Saliency min/max="<<Smin<<"/"<<Smax);
	
	//## Delete matrix		
	if(ColorDistMatrix) ColorDistMatrix->Delete();
	if(SpatialDistMatrix) SpatialDistMatrix->Delete();

	
	//## Normalize saliency and fill maps
	for(size_t i=0;i<nRegions;i++){
			
		//Normalize Saliency color
		double S= SList[i];
		double Snorm= NormMin + (NormMax-NormMin)*(S-Smin)/(Smax-Smin);
		double Saliency= 1.-Snorm;
	
		//Fill image
		for(size_t j=0;j<regions[i]->GetNPixels();j++){//loop on pixels inside region
			Pixel* thisPixel= regions[i]->GetPixel(j);
			//double x= thisPixel->x;
			//double y= thisPixel->y;
			//saliencyImg->FillPixel(x,y,Saliency);
			long int ix= thisPixel->ix;
			long int iy= thisPixel->iy;
			saliencyImg->FillPixel(ix,iy,Saliency);
		}//end loop pixels in region

	}//end loop regions
	
	return saliencyImg;

}//close ComputeSaliencyMap()


Image* SaliencyFilter::ComputeMultiResoSaliencyMap(Image* img,int resoMin,int resoMax,int resoStep,double beta,int minRegionSize,double knnFactor,bool useRobustPars,double expFalloffPar,double distanceRegPar,double salientMultiplicityThrFactor,bool addBkgMap,bool addNoiseMap,ImgBkgData* bkgData,double saliencyThrFactor,double imgThrFactor)
{	
	//## Check input img
	if(!img){
		ERROR_LOG("Null ptr to given input image!");
		return 0;
	}
	
	//## Check bkg data	
	if(!bkgData && (addBkgMap || addNoiseMap)){
		ERROR_LOG("Selected to use bkgdata in saliency computation but no bkg data given!");
		return 0;
	}

	//## Compute img stats
	if(!img->HasStats()){
		img->ComputeStats(true,false,false);
	}
	ImgStats* imgStats= img->GetPixelStats();
	if(!imgStats){
		ERROR_LOG("No stats available for this image (hint: compute stats first)!");
		return 0;
	}
	double imgMedian= imgStats->median;
	double imgMedianThr= imgThrFactor*imgMedian;
	int nReso= (resoMax-resoMin)/resoStep + 1;
	
	TString imgName= Form("%s_saliencyMean",img->GetName().c_str());
	Image* saliencyImg_mean= (Image*)img->GetCloned(std::string(imgName),true,true);
	saliencyImg_mean->SetName(std::string(imgName));
	saliencyImg_mean->Reset();

	double NormMin= 1;
	double NormMax= 256;
	std::vector<Image*> salMaps;
	std::vector<double> salMapsThresholds;
	int nbins= 100;
	
	for(int i=0;i<nReso;i++){
		
		//Compute saliency map @ current reso
		int reso= resoMin + i*resoStep;
		INFO_LOG("Computing saliency map @ reso "<<reso<<" (step="<<resoStep<<"/"<<nReso<<")");
		Image* salMap= ComputeSaliencyMap(img,reso,beta,minRegionSize,knnFactor,useRobustPars,expFalloffPar,distanceRegPar);

		//Normalize saliency map
		INFO_LOG("Normalizing saliency map @ reso "<<reso<<" (step="<<resoStep<<"/"<<nReso<<")");
		Image* salMap_norm= salMap->GetNormalizedImage("LINEAR",NormMin,NormMax);
		if(salMap) {
			delete salMap;	
			salMap= 0;
		}

		INFO_LOG("Add saliency map @ reso "<<reso<<" to cumulative map...");
		saliencyImg_mean->Add(salMap_norm);
		salMaps.push_back(salMap_norm);

		//Compute stats		
		INFO_LOG("Computing stats for normalized saliency map @ reso "<<reso<<" (step="<<resoStep<<")");
		salMap_norm->ComputeStats(true,false,false);
		double salMedian= (salMap_norm->GetPixelStats())->median;
		double medianThr= saliencyThrFactor*salMedian;
		double salThr= medianThr;
		salMapsThresholds.push_back(salThr);

		INFO_LOG("Saliency map norm stats @ reso "<<reso<<" (median="<<salMedian<<", salThr="<<salThr<<")");
		//salMap_norm->PrintStats();

	}//end loop reso
	
	//Normalize final saliency
	if(nReso>1){
		saliencyImg_mean->Scale(1./(double)nReso);
	}
	double minSaliency= saliencyImg_mean->GetMinimum();
	
	imgName= Form("%s_saliencyCombined",img->GetName().c_str());
	Image* saliencyImg= (Image*)saliencyImg_mean->GetCloned(std::string(imgName),true,true);
	saliencyImg->SetName(std::string(imgName));
	saliencyImg->Reset();
	saliencyImg_mean->Delete();
	
	//## Normalize saliency (using adaptive threshold)
	int salientMultiplicityThr= static_cast<int>(std::round(salientMultiplicityThrFactor*nReso));

	INFO_LOG("Normalize saliency sum over reso (minSaliency="<<minSaliency<<", salientMultiplicityThr="<<salientMultiplicityThr<<")...");
	
	for(long int i=0;i<saliencyImg->GetNx();i++){
		for(long int j=0;j<saliencyImg->GetNy();j++){
			double imgBinContent= img->GetBinContent(i,j);	
			if(imgBinContent==0) continue;
			bool isPositiveExcess= (imgBinContent>imgMedianThr);
			double wmin= TMath::Infinity();//1.e+99;
			double wmax= -TMath::Infinity();//-1.e+99;
			int saliencyMultiplicity= 0;
			for(size_t k=0;k<salMaps.size();k++){
				double thisw= salMaps[k]->GetBinContent(i,j);
				double thisThreshold= salMapsThresholds[k];
				if(thisw>thisThreshold) saliencyMultiplicity++;
				if(thisw<wmin) wmin= thisw;
				if(thisw>=wmax) wmax= thisw;
			}//end loop multi reso
			
			if(fabs(wmin)==TMath::Infinity() || fabs(wmax)==TMath::Infinity()){
				WARN_LOG("Saliency min/max over scales for pixel ("<<i<<","<<j<<") is inf!");
			}

			if(saliencyMultiplicity>=salientMultiplicityThr){
				if(isPositiveExcess) saliencyImg->SetBinContent(i,j,wmax);
				else saliencyImg->SetBinContent(i,j,minSaliency);
			}
			else {
				saliencyImg->SetBinContent(i,j,wmin);
			}
		}//end loop bins Y
	}//end loop bins X
	

	//Clear map list
	for(unsigned int k=0;k<salMaps.size();k++){
		if(salMaps[k]) salMaps[k]->Delete();		
	}
	salMaps.clear();
	
	//Normalize bkg and noise maps
	if(addBkgMap && bkgData && bkgData->BkgMap){
		INFO_LOG("Normalize and add bkg map to saliency estimate...");
		Image* bkgImg= (bkgData->BkgMap)->GetNormalizedImage("LINEAR",NormMin,NormMax);
		if(bkgImg){
			INFO_LOG("Adding bkg map to saliency estimate...");
			saliencyImg->Add(bkgImg);	
			bkgImg->Delete();
		}
		else{
			WARN_LOG("Failed to normalize bkg map, cannot add it to saliency computation!");
		}
	}//close if

	if(addNoiseMap && bkgData && bkgData->NoiseMap){
		INFO_LOG("Normalize and add noise map to saliency estimate...");
		Image* noiseImg= (bkgData->NoiseMap)->GetNormalizedImage("LINEAR",NormMin,NormMax);
		if(noiseImg){
			saliencyImg->Add(noiseImg);
			noiseImg->Delete();	
		}
		else{
			WARN_LOG("Failed to normalize noise map, cannot add it to saliency computation!");
		}
	}//close if

	//Compute saliency stats
	INFO_LOG("Multi-reso combined saliency stats (before normalization and after bkg/noise sum)");
	saliencyImg->ComputeStats(true,false,true);
	//saliencyImg->PrintStats();

	//Normalize final map
	DEBUG_LOG("Normalize final maps...");
	imgName= Form("%s_saliencyMultiReso",img->GetName().c_str());
	Image* saliencyMap= saliencyImg->GetNormalizedImage("LINEAR",NormMin,NormMax);
	saliencyMap->SetName(std::string(imgName));
	if(saliencyImg) saliencyImg->Delete();

	for(long int i=0;i<saliencyMap->GetNx();i++){
		for(long int j=0;j<saliencyMap->GetNy();j++){
			double imgBinContent= img->GetBinContent(i,j);
			if(imgBinContent==0){
				saliencyMap->SetBinContent(i,j,0.);
			}
		}
	}		
	saliencyMap->ComputeStats(true,false,true);

	return saliencyMap;

}//close GetMultiResoSaliencyMap()
	


}//close namespace
