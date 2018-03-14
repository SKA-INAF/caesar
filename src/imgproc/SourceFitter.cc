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
* @file SourceFitter.cc
* @class SourceFitter
* @brief SourceFitter
*
* Class to fit a source image with a mixture of gaussian/skew normal/skew-t bivariate functions
* @author S. Riggi
* @date 01/09/2017
*/

#include <SourceFitter.h>

#include <Image.h>
#include <Source.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <SysUtils.h>
#include <Logger.h>
#include <Consts.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TFitResult.h>
#include <TBackCompFitter.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>

//C++ headers
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
#include <chrono>

using namespace std;

ClassImp(Caesar::SourceComponentPars)
ClassImp(Caesar::SourceFitPars)
ClassImp(Caesar::SourceFitter)

namespace Caesar {

//Static variables
int SourceFitter::m_NFitComponents= 0;
SourceFitPars SourceFitter::m_sourceFitPars;
int SourceFitter::m_fitStatus= eFitUnknownStatus;
TH2D* SourceFitter::m_fluxMapHisto= 0;

//Constructor
SourceFitter::SourceFitter()
	: TObject()
{
	m_fluxMapHisto= 0;
	m_bkgMean= 0.;
	m_rmsMean= 0.;

}//close costructor

//Destructor
SourceFitter::~SourceFitter()
{
	CodeUtils::DeletePtr<TH2D>(m_fluxMapHisto);
	
}//close destructor


int SourceFitter::InitData(Source* aSource,SourceFitOptions& fitOptions)
{
	//## Check if stats has been computed, otherwise compute them
	if(!aSource->HasStats()){
		INFO_LOG("Input source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") has no stats computed, computing them now...");
		aSource->ComputeStats();
	}

	//## Get source pars
	float Xmin, Xmax, Ymin, Ymax;
	aSource->GetSourceRange(Xmin,Xmax,Ymin,Ymax);
	long int nBoxX= Xmax - Xmin + 1;
	long int nBoxY= Ymax - Ymin + 1;
	std::vector<Pixel*> pixels= aSource->GetPixels();
		
	//## Initialize flux histo/map
	if(m_fluxMapHisto) CodeUtils::DeletePtr<TH2D>(m_fluxMapHisto);
	m_fluxMapHisto= new TH2D("","",nBoxX,Xmin-0.5,Xmax+0.5,nBoxY,Ymin-0.5,Ymax+0.5);
	m_fluxMapHisto->Sumw2();

	//## Fill source flux histo
	INFO_LOG("Filling flux histo for source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") ...");
	m_bkgMean= 0.;
	m_rmsMean= 0.;
	long int ndata= 0;

	if(fitOptions.useFluxZCut){
		for(size_t k=0;k<pixels.size();k++){
			double x= pixels[k]->x;
			double y= pixels[k]->y;
			double S= pixels[k]->S;
			//double Scurv= pixels[k]->GetCurv();	
			std::pair<double,double> bkgData= pixels[k]->GetBkg();	
			double bkg= bkgData.first;
			double rms= bkgData.second;
			if(rms>0) {
				double Z= (S-bkg)/rms;		
				if(Z>fitOptions.fluxZThrMin) {
					m_fluxMapHisto->Fill(x,y,S);
					m_bkgMean+= bkg;
					m_rmsMean+= rms;
					ndata++;
				}
			}
			else{		
				m_fluxMapHisto->Fill(x,y,S);
				m_bkgMean+= bkg;
				m_rmsMean+= rms;
				ndata++;
			}
		}//end loop pixels

		m_bkgMean/= (double)(ndata);
		m_rmsMean/= (double)(ndata);
	}//close if
	else{
		for(size_t k=0;k<pixels.size();k++){
			double x= pixels[k]->x;
			double y= pixels[k]->y;
			double S= pixels[k]->S;
			//double Scurv= pixels[k]->GetCurv();	
			std::pair<double,double> bkgData= pixels[k]->GetBkg();	
			double bkg= bkgData.first;
			double rms= bkgData.second;
			 
			m_fluxMapHisto->Fill(x,y,S);
			m_bkgMean+= bkg;
			m_rmsMean+= rms;
		}//end loop pixels
	
		m_bkgMean/= (double)(pixels.size());
		m_rmsMean/= (double)(pixels.size());
	}//close else
	
	INFO_LOG("Source (name="<<aSource->GetName()<<", N="<<pixels.size()<<", pos("<<aSource->X0<<","<<aSource->Y0<<")) bkg info: <bkg>="<<m_bkgMean<<", <rms>="<<m_rmsMean);

	//Check if histo has entries
	if(m_fluxMapHisto->GetSumOfWeights()<=0){
		ERROR_LOG("No entries present in flux map histogram (hint: check if properly filled, e.g. bin range problems)!");
		return -1;
	}

	return 0;

}//close InitData()

int SourceFitter::CheckFitOptions(SourceFitOptions& fitOptions)
{
	//## Check fit options
	if(fitOptions.peakMinKernelSize>fitOptions.peakMaxKernelSize){
		ERROR_LOG("Invalid peak kernel sizes given (hint: min kernel must be larger or equal to max kernel size)!");
		return -1;
	}
	if(fitOptions.peakMinKernelSize<=0 || fitOptions.peakMinKernelSize%2==0 || fitOptions.peakMaxKernelSize<=0 || fitOptions.peakMaxKernelSize%2==0){
		ERROR_LOG("Invalid peak kernel sizes given (hint: kernel size must be positive and odd)!");
		return -1;
	}

	return 0;

}//close CheckFitOptions()


int SourceFitter::EstimateFitComponents(std::vector<std::vector<double>>& fitPars_start,Source* aSource,SourceFitOptions& fitOptions)
{
	//## Init pars
	std::vector<int> kernels;
	for(int k=fitOptions.peakMinKernelSize;k<=fitOptions.peakMaxKernelSize;k+=2){
		kernels.push_back(k);
	}
	
	//## Estimate the number of components present in source from nested sources
	std::vector<Source*> nestedSources= aSource->GetNestedSources();
	int nComponents= static_cast<int>(nestedSources.size());

	if(nComponents==0){
		//Get centroid & peak
		double meanX= aSource->GetSx();
		double meanY= aSource->GetSy();
		double Smax= aSource->GetSmax();
		INFO_LOG("No nested components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<", pos("<<meanX<<","<<meanY<<")), will search for peak components ...");	
		
		//Compute standard deviations along axis
		double sigmaX_sample= 0;
		double sigmaY_sample= 0;
		double covXY_sample= 0;
		if(aSource->GetSampleStdDev(sigmaX_sample,sigmaY_sample,covXY_sample)<0){
			WARN_LOG("Failed to compute source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") standard deviations!");
			return -1;
		}
		double sigmaX= sigmaX_sample;
		double sigmaY= sigmaY_sample;
		double theta= StatsUtils::GetGaus2DThetaPar(sigmaX_sample,sigmaY_sample,covXY_sample);		
		StatsUtils::GetEllipseParsFromCovMatrix(sigmaX,sigmaY,theta,sigmaX_sample,sigmaY_sample,covXY_sample);

		//Estimate number of components from detected peaks
		std::vector<TVector2> peaks;
		int status= aSource->FindComponentPeaks(
			peaks,
			fitOptions.peakZThrMin, fitOptions.nMaxComponents,
			fitOptions.peakShiftTolerance,
			kernels,fitOptions.peakKernelMultiplicityThr
		);
		if(status<0){
			WARN_LOG("Failed to find component peaks in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")!");
			return -1;
		}

		//NB: If only one peak found initialize fit component to the entire source
		//otherwise use peaks found as centroid start values and beam as gaussian sigma pars
		if(peaks.empty()){
			INFO_LOG("No peaks found in this source (hint: could be a diffuse source and peaks are below desired threshold), nothing to be fit");
			return -1;
		}
		if(peaks.size()==1){
			fitPars_start.push_back( std::vector<double>() );
			fitPars_start[0].push_back(Smax);
			fitPars_start[0].push_back(meanX);
			fitPars_start[0].push_back(meanY);
			fitPars_start[0].push_back(sigmaX);
			fitPars_start[0].push_back(sigmaY);
			fitPars_start[0].push_back(theta);
		}
		else{
			for(size_t i=0;i<peaks.size();i++){
				double x= peaks[i].X();
				double y= peaks[i].Y();
				long int gbin= m_fluxMapHisto->FindBin(x,y);
				if(m_fluxMapHisto->IsBinOverflow(gbin) || m_fluxMapHisto->IsBinUnderflow(gbin)){
					WARN_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
					return -1;
				}
				double Speak= m_fluxMapHisto->GetBinContent(gbin);
				double sigmaX= fitOptions.bmaj/GausSigma2FWHM;
				double sigmaY= fitOptions.bmin/GausSigma2FWHM;
				double theta= fitOptions.bpa;

				fitPars_start.push_back( std::vector<double>() );
				fitPars_start[i].push_back(Speak);
				fitPars_start[i].push_back(x);
				fitPars_start[i].push_back(y);
				fitPars_start[i].push_back(sigmaX);
				fitPars_start[i].push_back(sigmaY);
				fitPars_start[i].push_back(theta);
			}//end loop peaks
		}//close else
		
		
	}//close if
	else{
		INFO_LOG("#"<<nComponents<<" nested components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<"), sorting and selecting best ones...");	
	
		//Sort nested sources by peak flux and select brighter
		std::sort(nestedSources.begin(),nestedSources.end(),SourceCompareByLargerPeakFlux());
		std::vector<size_t> selected_index;
		int nSelComponents= 0;
		int componentCounter= 0;

		for(size_t i=0;i<nestedSources.size();i++){
			//Get nested source pars
			double Smedian= nestedSources[i]->Median;
			double Smad= nestedSources[i]->MedianRMS;
			double x= nestedSources[i]->GetSx();
			double y= nestedSources[i]->GetSy();
			double Speak= nestedSources[i]->GetSmax();
			double Zpeak_imgbkg= 0;	
			double Zpeak_sourcebkg= 0;
			if(m_rmsMean!=0) Zpeak_imgbkg= (Speak-m_bkgMean)/m_rmsMean;
			if(Smad!=0) Zpeak_sourcebkg= (Speak-Smedian)/Smad;

			//Remove faint peaks
			//if(Zpeak_imgbkg<fitOptions.peakZThrMin) {
			if(Zpeak_sourcebkg<fitOptions.peakZThrMin) {
				INFO_LOG("Removing nested source no. "<<i<<" (name="<<nestedSources[i]->GetName()<<", pos("<<x<<","<<y<<")) from the component list as below peak significance thr (Zpeak(imgbkg)="<<Zpeak_imgbkg<<", Zpeak(sourcebkg)="<<Zpeak_sourcebkg<<"<"<<fitOptions.peakZThrMin<<")");
				continue;
			}

			//Add nested source index to list
			nSelComponents++;
			selected_index.push_back(i);

			//Check if number of maximum components has been reached
			if(nSelComponents>=fitOptions.nMaxComponents){
				INFO_LOG("Maximum number of fit components (N="<<fitOptions.nMaxComponents<<") reached, stop adding fit components...");
				break;
			}	

		}//end loop nested sources

		//Check if there are components left after selection
		if(selected_index.empty()){
			WARN_LOG("No components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") after selection (hint: check your fit component selection pars)!");	
			return -1;
		}

		//Loop over selected components and fill parameters
		for(size_t i=0;i<selected_index.size();i++){
			size_t index= selected_index[i];
			Source* nestedSource= nestedSources[index];

			//Get centroid & peak
			double meanX= nestedSource->GetSx();
			double meanY= nestedSource->GetSy();
			double Smax= nestedSource->GetSmax();

			//Compute sigma
			double sigmaX_sample= 0;
			double sigmaY_sample= 0;
			double covXY_sample= 0;
			if(nestedSource->GetSampleStdDev(sigmaX_sample,sigmaY_sample,covXY_sample)<0){
				WARN_LOG("Failed to compute standard deviations for nested component no. "<<i<<" of source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")!");
				return -1;
			}
			double sigmaX= sigmaX_sample;
			double sigmaY= sigmaY_sample;
			double theta= StatsUtils::GetGaus2DThetaPar(sigmaX_sample,sigmaY_sample,covXY_sample);		
			StatsUtils::GetEllipseParsFromCovMatrix(sigmaX,sigmaY,theta,sigmaX_sample,sigmaY_sample,covXY_sample);

			//Estimate number of components from detected peaks
			std::vector<TVector2> peaks;
			int status= nestedSource->FindComponentPeaks(
				peaks,
				fitOptions.peakZThrMin, fitOptions.nMaxComponents,
				fitOptions.peakShiftTolerance,
				kernels,fitOptions.peakKernelMultiplicityThr
			);
			if(status<0){
				WARN_LOG("Failed to find component peaks in nested component no. "<<i<<" of source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")!");
				return -1;
			}

			//NB: If only one peak found initialize fit component to the entire source
			//otherwise use peaks found as centroid start values and beam as gaussian sigma pars
			if(peaks.empty()){
				INFO_LOG("No peaks found in nested component no. "<<i<<" of source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") (hint: could be a diffuse source and peaks are below desired threshold), skip to next nested component...");
				continue;
			}			
			if(peaks.size()<=1){
				fitPars_start.push_back( std::vector<double>() );	
				fitPars_start[componentCounter].push_back(Smax);
				fitPars_start[componentCounter].push_back(meanX);
				fitPars_start[componentCounter].push_back(meanY);
				fitPars_start[componentCounter].push_back(sigmaX);
				fitPars_start[componentCounter].push_back(sigmaY);
				fitPars_start[componentCounter].push_back(theta);
				componentCounter++;
			}//close if

			else{
				//Loop over peaks and add component
				for(size_t j=0;j<peaks.size();j++){
					double x= peaks[j].X();
					double y= peaks[j].Y();
					long int gbin= m_fluxMapHisto->FindBin(x,y);
					if(m_fluxMapHisto->IsBinOverflow(gbin) || m_fluxMapHisto->IsBinUnderflow(gbin)){
						WARN_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
						return -1;
					}
					double Speak= m_fluxMapHisto->GetBinContent(gbin);
					double sigmaX= fitOptions.bmaj/GausSigma2FWHM;
					double sigmaY= fitOptions.bmin/GausSigma2FWHM;
					double theta= fitOptions.bpa;
	
					fitPars_start.push_back( std::vector<double>() );	
					fitPars_start[componentCounter].push_back(Speak);
					fitPars_start[componentCounter].push_back(x);
					fitPars_start[componentCounter].push_back(y);
					fitPars_start[componentCounter].push_back(sigmaX);
					fitPars_start[componentCounter].push_back(sigmaY);
					fitPars_start[componentCounter].push_back(theta);
					componentCounter++;
				}//end loop peaks
			}//close else
	
		}//end loop selected nested sources
		
	}//close else

	INFO_LOG("#"<<fitPars_start.size()<<" fit components will be fitted...");

	
	return 0;

}//close EstimateFitComponents()

int SourceFitter::DoFit(Source* aSource,SourceFitOptions& fitOptions,std::vector< std::vector<double> >& fitPars_start)
{
	//## Get source pars
	float Xmin, Xmax, Ymin, Ymax;
	aSource->GetSourceRange(Xmin,Xmax,Ymin,Ymax);
	double Smax= aSource->GetSmax();
	double Smin= aSource->GetSmin();

	//## Initialize fit function with start parameters
	m_NFitComponents= static_cast<int>(fitPars_start.size());
	int nComponentPars= 6;
	int nFitPars= nComponentPars*m_NFitComponents + 1;//fit components + constant offset
	double fitRangeXmin= m_fluxMapHisto->GetXaxis()->GetXmin();
	double fitRangeXmax= m_fluxMapHisto->GetXaxis()->GetXmax();
	double fitRangeYmin= m_fluxMapHisto->GetYaxis()->GetXmin();
	double fitRangeYmax= m_fluxMapHisto->GetYaxis()->GetXmax();
	
	std::string parNamePrefix[nComponentPars]= {"A","x0","y0","sigmaX","sigmaY","theta"};

	TF2* sourceFitFcn= new TF2("sourceFitFcn",
		SourceFitter::Gaus2DMixtureFcn,
		fitRangeXmin,fitRangeXmax,
		fitRangeYmin,fitRangeYmax,
		nFitPars
	);
	sourceFitFcn->SetNpx(1000);
	sourceFitFcn->SetNpy(1000);	
	INFO_LOG("Created source fit function (nPars="<<nFitPars<<") with range: X["<<fitRangeXmin<<","<<fitRangeXmax<<"], Y["<<fitRangeYmin<<","<<fitRangeYmax<<"]");

	//## Set start fit pars
	int par_counter= 0;
	std::vector<std::vector<double>> parLimits_min;
	std::vector<std::vector<double>> parLimits_max;

	//- Offset
	double offset= 0.;
	double offset_min= 0;
	double offset_max= 0;
	if(fitOptions.fixBkg){//fix bkg level fit par
		if(fitOptions.useEstimatedBkgLevel) offset= m_bkgMean;//use estimated avg bkg
		else offset= fitOptions.fixedBkgLevel;//use user-supplied bkg level
		sourceFitFcn->SetParameter(nFitPars-1,offset);
		sourceFitFcn->FixParameter(nFitPars-1,offset);	
	}
	else{
		if(fitOptions.useEstimatedBkgLevel) offset= m_bkgMean;//use estimated avg bkg
		else offset= fitOptions.fixedBkgLevel;//use user-supplied bkg level
		if(m_bkgMean<Smin || m_bkgMean>Smax) offset= Smin;
		offset_min= std::min(Smin,offset-fabs(m_rmsMean));
		offset_max= std::min(Smax,offset+fabs(m_rmsMean));
		sourceFitFcn->SetParameter(nFitPars-1,offset);
		if(fitOptions.limitBkgInFit) {
			sourceFitFcn->SetParLimits(nFitPars-1,offset_min,offset_max);
			INFO_LOG("offset="<<offset<<" ["<<offset_min<<","<<offset_max<<"]");	
		}
	}
	sourceFitFcn->SetParName(nFitPars-1,"offset");
	
	
	for(int i=0;i<m_NFitComponents;i++){
		parLimits_min.push_back(std::vector<double>(nComponentPars,0));
		parLimits_max.push_back(std::vector<double>(nComponentPars,0));
		double Speak= fitPars_start[i][0];
		double x= fitPars_start[i][1];
		double y= fitPars_start[i][2];
		double sigmaX= fitPars_start[i][3];
		double sigmaY= fitPars_start[i][4];
		double theta= fitPars_start[i][5];

		//Subtract offset from peak
		Speak-= offset;
		
		//## Set i-th component parameters
		//- Amplitude
		double Speak_min= std::min(Smin, Speak*(1 - fitOptions.amplLimit) );
		double Speak_max= std::max(Smax, Speak*(1 + fitOptions.amplLimit) );
		sourceFitFcn->SetParName(par_counter,Form("%s_%d",parNamePrefix[0].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter,Speak);
		if(fitOptions.limitAmplInFit){
			sourceFitFcn->SetParLimits(par_counter,Speak_min,Speak_max);
			parLimits_min[i][0]= Speak_min;
			parLimits_max[i][0]= Speak_max;
			INFO_LOG("Speak="<<Speak<<" ["<<Speak_min<<","<<Speak_max<<"]");
		}

		//- Centroids
		double centroidLimit= 0.5 * sqrt(pow(fitOptions.bmaj,2) + pow(fitOptions.bmin,2));
		double x0_min= x - centroidLimit;
		double x0_max= x + centroidLimit;
		double y0_min= y - centroidLimit;
		double y0_max= y + centroidLimit;
		sourceFitFcn->SetParName(par_counter+1,Form("%s_%d",parNamePrefix[1].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+1,x);
		sourceFitFcn->SetParName(par_counter+2,Form("%s_%d",parNamePrefix[2].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+2,y);
		if(fitOptions.limitCentroidInFit){
			sourceFitFcn->SetParLimits(par_counter+1,x0_min,x0_max);	
			sourceFitFcn->SetParLimits(par_counter+2,y0_min,y0_max);
			parLimits_min[i][1]= x0_min;
			parLimits_max[i][1]= x0_max;
			parLimits_min[i][2]= y0_min;
			parLimits_max[i][2]= y0_max;
			INFO_LOG("(x,y)=("<<x<<","<<y<<")"<<" bounds x("<<x0_min<<","<<x0_max<<") y("<<y0_min<<","<<y0_max<<")");
		}

		
		//- Sigmas
		double sourceSigmaMax_x= fabs(Xmax-Xmin)*sqrt(2.)/GausSigma2FWHM;
		double sourceSigmaMax_y= fabs(Ymax-Ymin)*sqrt(2.)/GausSigma2FWHM;
		double sigmaX_min= sigmaX*(1-fitOptions.sigmaLimit);
		double sigmaY_min= sigmaY*(1-fitOptions.sigmaLimit);
		double sigmaX_max= std::max(sourceSigmaMax_x,sigmaX*(1+fitOptions.sigmaLimit));
		double sigmaY_max= std::max(sourceSigmaMax_y,sigmaY*(1+fitOptions.sigmaLimit));
		sourceFitFcn->SetParName(par_counter+3,Form("%s_%d",parNamePrefix[3].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+3,sigmaX);		
		sourceFitFcn->SetParName(par_counter+4,Form("%s_%d",parNamePrefix[4].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+4,sigmaY);

		if(fitOptions.limitSigmaInFit){
			sourceFitFcn->SetParLimits(par_counter+3,sigmaX_min,sigmaX_max);
			sourceFitFcn->SetParLimits(par_counter+4,sigmaY_min,sigmaY_max);
			parLimits_min[i][3]= sigmaX_min;
			parLimits_max[i][3]= sigmaX_max;
			parLimits_min[i][4]= sigmaY_min;
			parLimits_max[i][4]= sigmaY_max;
			INFO_LOG("(sigmaX,sigmaY)=("<<sigmaX<<","<<sigmaY<<")"<<" bounds sigmaX("<<sigmaX_min<<","<<sigmaX_max<<") y("<<sigmaY_min<<","<<sigmaY_max<<"), bmaj="<<fitOptions.bmaj<<", bmin="<<fitOptions.bmin);
		}

		//- Theta		
		double theta_min= theta - fitOptions.thetaLimit;
		double theta_max= theta + fitOptions.thetaLimit;
		sourceFitFcn->SetParName(par_counter+5,Form("%s_%d",parNamePrefix[5].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+5,theta);
		if(fitOptions.limitThetaInFit){
			sourceFitFcn->SetParLimits(par_counter+5,theta_min,theta_max);
			parLimits_min[i][5]= theta_min;
			parLimits_max[i][5]= theta_max;
			INFO_LOG("theta="<<theta<<" bounds ("<<theta_min<<","<<theta_max<<")");
		}

		//Update par counter
		par_counter+= nComponentPars;
	
	}//end loop components

	
	//==============================================
	//==             PRE-FIT
	//==============================================
	//## Fix sigmas to blob size?
	if(fitOptions.fixSigmaInPreFit){
		for(int i=0;i<m_NFitComponents;i++){
			//Sigma X
			TString parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			double parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);

			//Sigma Y
			parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
			parNumber= sourceFitFcn->GetParNumber(parName);
			parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);
		}
	}//close if fix sigma

	//## Fix theta to beam bpa or to estimated blob theta 
	for(int i=0;i<m_NFitComponents;i++){
		TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
		int parNumber= sourceFitFcn->GetParNumber(parName);
		double parValue= sourceFitFcn->GetParameter(parName);
		//double parValue= fitOptions.bpa;
		sourceFitFcn->FixParameter(parNumber,parValue);
	}
	
	//## Fix offset 
	sourceFitFcn->FixParameter(nFitPars-1,offset);

	//## Perform pre-fit
	int prefitStatus= m_fluxMapHisto->Fit(sourceFitFcn,"RWN");
	if(prefitStatus!=0){
		WARN_LOG("Source pre-fit failed or did not converge.");
		m_fitStatus= eFitNotConverged;
	}

	//==============================================
	//==             FIT
	//==============================================
	//## Fix sigmas?
	if(fitOptions.fixSigma){
		for(int i=0;i<m_NFitComponents;i++){
			//Sigma X
			TString parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			double parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);

			//Sigma Y
			parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
			parNumber= sourceFitFcn->GetParNumber(parName);
			parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);
		}
	}//close if fix sigma
	else{
		for(int i=0;i<m_NFitComponents;i++){
			//Sigma X
			TString parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			sourceFitFcn->ReleaseParameter(parNumber);
			if(fitOptions.limitSigmaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sourceFitFcn->GetParameter(parNumber)<<") to ["<<parLimits_min[i][3]<<","<<parLimits_max[i][3]<<"]");
				sourceFitFcn->SetParLimits(parNumber,parLimits_min[i][3],parLimits_max[i][3]);
			}

			//Sigma Y
			parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
			parNumber= sourceFitFcn->GetParNumber(parName);
			sourceFitFcn->ReleaseParameter(parNumber);
			if(fitOptions.limitSigmaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sourceFitFcn->GetParameter(parNumber)<<") to ["<<parLimits_min[i][4]<<","<<parLimits_max[i][4]<<"]");
				sourceFitFcn->SetParLimits(parNumber,parLimits_min[i][4],parLimits_max[i][4]);	
			}
		}
	}//close else

	//## Fix theta?
	if(fitOptions.fixTheta){
		for(int i=0;i<m_NFitComponents;i++){
			TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			double parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);
		}
	}//close if fix theta
	else { 
		for(int i=0;i<m_NFitComponents;i++){
			TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			sourceFitFcn->ReleaseParameter(parNumber);
			if(fitOptions.limitThetaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sourceFitFcn->GetParameter(parNumber)<<") to ["<<parLimits_min[i][5]<<","<<parLimits_max[i][5]<<"]");
				sourceFitFcn->SetParLimits(parNumber,parLimits_min[i][5],parLimits_max[i][5]);
			}
		}
	}//close else
	
	//## Release offset?
	if(!fitOptions.fixBkg){
		sourceFitFcn->ReleaseParameter(nFitPars-1);
		if(fitOptions.limitBkgInFit) {
			INFO_LOG("Limiting offset par (parNumber="<<nFitPars-1<<", value="<<sourceFitFcn->GetParameter(nFitPars-1)<<") to ["<<offset_min<<","<<offset_max<<"]");
			sourceFitFcn->SetParLimits(nFitPars-1,offset_min,offset_max);
		}
	}

	//## Fit again
	m_fluxMapHisto->Fit(sourceFitFcn,"SRWN");
	TBackCompFitter* fitter = (TBackCompFitter *) TVirtualFitter::GetFitter();
	const ROOT::Fit::FitResult& fitRes = fitter->GetFitResult();
	
	//==============================================
	//==             FIT RESULTS
	//==============================================
	//## Retrieve fit results
	m_fitStatus= eFitNotConverged;
	int fitMinimizerStatus= -1;
	int fitStatus= -1;
	bool hasParsAtLimits= false;
	
	if(fitRes.IsValid()) {
		fitStatus= 0;
		m_fitStatus= eFitConverged;
	}
	
	//NB: Minuit fit status = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult (4: MIGRAD error, 40: MINOS error, 400: HESSE error, 4000: IMPROVE error)
	fitMinimizerStatus= fitRes.Status();

	//Check if has parameters at bounds
	hasParsAtLimits= HasFitParsAtLimit(fitRes);
	if(hasParsAtLimits) m_fitStatus= eFitConvergedWithWarns;
	

	double Chi2= sourceFitFcn->GetChisquare();//fitRes->Chi2();
	double NDF= sourceFitFcn->GetNDF();//fitRes->Ndf();
	int NFreePars= sourceFitFcn->GetNumberFreeParameters();//fitRes->NFreeParameters();
	int NFittedBins= sourceFitFcn->GetNumberFitPoints();//fitRes->NFreeParameters();
	INFO_LOG("Source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") fit info: status="<<fitStatus<<" fitterStatus="<<m_fitStatus<<" (hasParsAtLimits? "<<hasParsAtLimits<<"), NDF="<<NDF<<", Chi2="<<Chi2<<", NFreePars="<<NFreePars<<" NFittedBins="<<NFittedBins);

	m_sourceFitPars.SetNComponents(m_NFitComponents);
	m_sourceFitPars.SetChi2(Chi2);
	m_sourceFitPars.SetNDF(NDF);
	m_sourceFitPars.SetNFreePars(NFreePars);
	m_sourceFitPars.SetNFitPoints(NFittedBins);
	m_sourceFitPars.SetStatus(m_fitStatus);
	m_sourceFitPars.SetMinimizerStatus(fitMinimizerStatus);
	
	par_counter= 0;

	for(int i=0;i<m_NFitComponents;i++){
		for(int j=0;j<nComponentPars;j++){
			std::string parName_global= std::string(sourceFitFcn->GetParName(par_counter));
			std::string parName= CodeUtils::ExtractSubString(parName_global,"_");
			double parVal= sourceFitFcn->GetParameter(par_counter);
			double parErr= sourceFitFcn->GetParError(par_counter);
			if(m_sourceFitPars.SetParValueAndError(i,parName,parVal,parErr)<0){
				WARN_LOG("Failed to set par "<<parName<<" (parName_global="<<parName_global<<") value and error (check par name)!");
			}
			par_counter++;
		}//end loop pars per component
	}//end loop components
	
	double fittedOffset= sourceFitFcn->GetParameter(par_counter);
	double fittedOffsetErr= sourceFitFcn->GetParError(par_counter);
	m_sourceFitPars.SetOffsetPar(fittedOffset);
	m_sourceFitPars.SetOffsetParErr(fittedOffsetErr);

	//==============================================
	//==             FIT RESIDUALS
	//==============================================
	//Compute fit residuals	
	std::vector<double> residuals;
	double residualMax= -1.e+99;
	double residualMin= 1.e+99;
	for (int i=0;i<m_fluxMapHisto->GetNbinsX();i++) {
		double x= m_fluxMapHisto->GetXaxis()->GetBinCenter(i+1);
		for (int j=0;j<m_fluxMapHisto->GetNbinsY();j++) {
			double y= m_fluxMapHisto->GetYaxis()->GetBinCenter(j+1);
			double data= m_fluxMapHisto->GetBinContent(i+1,j+1);
			if(!std::isnormal(data) || data==0) continue;
			double model= sourceFitFcn->Eval(x,y);
  		double res = data - model;
			if(res<residualMin) residualMin= res;
			if(res>residualMax) residualMax= res;
			residuals.push_back(res);
		}
  }
	
	//Compute residual stats
	double residualMean= 0;
	double residualRMS= 0;
	StatsUtils::ComputeMeanAndRMS(residualMean,residualRMS,residuals);
	double residualMedian= StatsUtils::GetMedianFast<double>(residuals);
	double residualMAD= StatsUtils::GetMADFast(residuals,residualMedian);	
	
	m_sourceFitPars.SetResidualMin(residualMin);
	m_sourceFitPars.SetResidualMax(residualMax);
	m_sourceFitPars.SetResidualMean(residualMean);
	m_sourceFitPars.SetResidualRMS(residualRMS);
	m_sourceFitPars.SetResidualMedian(residualMedian);
	m_sourceFitPars.SetResidualMAD(residualMAD);
	
	INFO_LOG("Fit residual stats: min/max="<<residualMin<<"/"<<residualMax<<", mean="<<residualMean<<", rms="<<residualRMS<<", median="<<residualMedian<<", mad="<<residualMAD);
	
	//Clear up data
	CodeUtils::DeletePtr<TF2>(sourceFitFcn);
	
	return 0;

}//close DoFit()


int SourceFitter::FitSource(Source* aSource,SourceFitOptions& fitOptions)
{
	//## Check input source
	if(!aSource){
		ERROR_LOG("Null ptr to source given!");
		m_fitStatus= eFitAborted;
		return -1;
	}

	//## Check fit options
	if(CheckFitOptions(fitOptions)<0){
		ERROR_LOG("Invalid fit options given!");
		m_fitStatus= eFitAborted;
		return -1;
	}

	INFO_LOG("Fitting source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") assuming these blob pars: {Bmaj(pix)="<<fitOptions.bmaj<<", Bmin(pix)="<<fitOptions.bmin<<", Bpa(deg)="<<fitOptions.bpa<<"}");

	//## Initialize data histos
	if(InitData(aSource,fitOptions)<0){
		ERROR_LOG("Failed to initialize fit data histo!");
		m_fitStatus= eFitAborted;
		return -1;
	}

	//## Estimate number of components
	std::vector< std::vector<double> > fitPars_start;
	if(EstimateFitComponents(fitPars_start,aSource,fitOptions)<0){
		ERROR_LOG("Failed to estimate fit components!");
		m_fitStatus= eFitAborted;
		return -1;	
	}
	if(fitPars_start.empty()){
		ERROR_LOG("No fitted components estimated (this should not occur!)");
		m_fitStatus= eFitAborted;
		return -1;
	}

	//## Perform fit
	if(DoFit(aSource,fitOptions,fitPars_start)<0){
		ERROR_LOG("Failed to perform source fit!");
		m_fitStatus= eFitAborted;
		return -1;	
	}
	
	return 0;

}//close FitSource()



/*
int SourceFitter::FitSource(Source* aSource,SourceFitOptions& fitOptions)
{
	//## Check input source
	if(!aSource){
		ERROR_LOG("Null ptr to source given!");
		m_fitStatus= eFitAborted;
		return -1;
	}

	//## Check fit options
	INFO_LOG("Fitting source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") assuming these blob pars: {Bmaj(pix)="<<fitOptions.bmaj<<", Bmin(pix)="<<fitOptions.bmin<<", Bpa(deg)="<<fitOptions.bpa<<"}");

	if(fitOptions.peakMinKernelSize>fitOptions.peakMaxKernelSize){
		ERROR_LOG("Invalid peak kernel sizes given (hint: min kernel must be larger or equal to max kernel size)!");
		m_fitStatus= eFitAborted;
		return -1;
	}
	if(fitOptions.peakMinKernelSize<=0 || fitOptions.peakMinKernelSize%2==0 || fitOptions.peakMaxKernelSize<=0 || fitOptions.peakMaxKernelSize%2==0){
		ERROR_LOG("Invalid peak kernel sizes given (hint: kernel size must be positive and odd)!");
		m_fitStatus= eFitAborted;
		return -1;
	}

	//## Check if stats has been computed, otherwise compute them
	if(!aSource->HasStats()){
		INFO_LOG("Input source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") has no stats computed, computing them now...");
		aSource->ComputeStats();
	}

	//## Get source pars
	float Xmin, Xmax, Ymin, Ymax;
	aSource->GetSourceRange(Xmin,Xmax,Ymin,Ymax);
	double Smax= aSource->GetSmax();
	double Smin= aSource->GetSmin();
	double Srms=  aSource->RMS; 
	double Smedian= aSource->Median;
	double Smad= aSource->MedianRMS;
	INFO_LOG("Source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") range: X["<<Xmin<<","<<Xmax<<"], Y["<<Ymin<<","<<Ymax<<"], Smedian="<<Smedian<<", Smad="<<Smax<<", Srms="<<Srms);

	//## Get source flux image
	INFO_LOG("Get flux histo & curv map for source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")");	
	long int nBoxX= Xmax - Xmin + 1;
	long int nBoxY= Ymax - Ymin + 1;
	TH2D* fluxMapHisto= new TH2D("","",nBoxX,Xmin-0.5,Xmax+0.5,nBoxY,Ymin-0.5,Ymax+0.5);
	fluxMapHisto->Sumw2();
	Image* peakSearchMap= new Image(nBoxX,nBoxY,Xmin,Ymin,"");
	
	std::vector<Pixel*> pixels= aSource->GetPixels();
	double bkgMean= 0.;
	double rmsMean= 0.;
	long int ndata= 0;
	if(fitOptions.useFluxZCut){
		for(size_t k=0;k<pixels.size();k++){
			double x= pixels[k]->x;
			double y= pixels[k]->y;
			double S= pixels[k]->S;
			double Scurv= pixels[k]->GetCurv();	
			std::pair<double,double> bkgData= pixels[k]->GetBkg();	
			double bkg= bkgData.first;
			double rms= bkgData.second;
			if(rms>0) {
				double Z= (S-bkg)/rms;		
				if(Z>fitOptions.fluxZThrMin) {
					fluxMapHisto->Fill(x,y,S);
					bkgMean+= bkg;
					rmsMean+= rms;
					ndata++;
				}
			}
			else{		
				fluxMapHisto->Fill(x,y,S);
				bkgMean+= bkg;
				rmsMean+= rms;
				ndata++;
			}
			//peakSearchMap->Fill(x,y,Scurv);	
			peakSearchMap->Fill(x,y,S);	
		}//end loop pixels

		bkgMean/= (double)(ndata);
		rmsMean/= (double)(ndata);
	}//close if
	else{
		for(size_t k=0;k<pixels.size();k++){
			double x= pixels[k]->x;
			double y= pixels[k]->y;
			double S= pixels[k]->S;
			double Scurv= pixels[k]->GetCurv();	
			std::pair<double,double> bkgData= pixels[k]->GetBkg();	
			double bkg= bkgData.first;
			double rms= bkgData.second;
			 
			fluxMapHisto->Fill(x,y,S);
			//peakSearchMap->Fill(x,y,Scurv);	
			peakSearchMap->Fill(x,y,S);	
			bkgMean+= bkg;
			rmsMean+= rms;
		}//end loop pixels
	
		bkgMean/= (double)(pixels.size());
		rmsMean/= (double)(pixels.size());
	}//close else
	
	INFO_LOG("source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") bkg info: <bkg>="<<bkgMean<<", <rms>="<<rmsMean);

	//Check if histo has entries
	if(fluxMapHisto->GetSumOfWeights()<=0){
		ERROR_LOG("No entries present in flux map histogram (hint: check if properly filled, e.g. bin range problems)!");
		m_fitStatus= eFitAborted;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		delete peakSearchMap;
		peakSearchMap= 0;
		return -1;
	}

	//## Estimate the number of components present in source
	//## Method: Find peaks in curvature map
	INFO_LOG("Finding peaks in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")");	
	std::vector<TVector2> peakPoints;
	bool skipBorders= true;
	std::vector<int> kernels;
	for(int k=fitOptions.peakMinKernelSize;k<=fitOptions.peakMaxKernelSize;k+=2){
		kernels.push_back(k);
	}

	if(peakSearchMap->FindPeaks(peakPoints,kernels,fitOptions.peakShiftTolerance,skipBorders,fitOptions.peakKernelMultiplicityThr)<0){
		WARN_LOG("Failed to find peaks in source curvature map!");
		m_fitStatus= eFitAborted;
		delete peakSearchMap;
		peakSearchMap= 0;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		return -1;
	}
	
	//Clear peak search map
	if(peakSearchMap){
		delete peakSearchMap;	
		peakSearchMap= 0;
	}

	//## Select peaks (skip peaks at boundary or faint peaks)
	std::vector<TVector2> peaks_selected;
	for(size_t i=0;i<peakPoints.size();i++){
		double x= peakPoints[i].X();
		double y= peakPoints[i].Y();
		long int gbin= fluxMapHisto->FindBin(x,y);
		if(fluxMapHisto->IsBinOverflow(gbin) || fluxMapHisto->IsBinUnderflow(gbin)){
			WARN_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
			m_fitStatus= eFitAborted;
			delete fluxMapHisto;
			fluxMapHisto= 0;	
			return -1;
		}

		//Remove faint peaks in case more than one peak is found
		if(peakPoints.size()>1){
			double Speak= fluxMapHisto->GetBinContent(gbin);
			double Zpeak_imgbkg= 0;	
			double Zpeak_sourcebkg= 0;
			if(rmsMean!=0) Zpeak_imgbkg= (Speak-bkgMean)/rmsMean;
			if(Smad!=0) Zpeak_sourcebkg= (Speak-Smedian)/Smad;
			//if(Zpeak_imgbkg<fitOptions.peakZThrMin) {
			if(Zpeak_sourcebkg<fitOptions.peakZThrMin) {
				INFO_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak_imgbkg="<<Zpeak_imgbkg<<", Zpeak_sourcebkg="<<Zpeak_sourcebkg<<"<"<<fitOptions.peakZThrMin<<")");
				continue;
			}
		}//close if

		//Remove peaks lying on the source contour
		if(aSource->HasContours() && aSource->IsPointOnContour(x,y,0.5)) {
			INFO_LOG("Removing peak ("<<x<<","<<y<<") from the list as lying on source contour");
			continue;
		}

		//Add peak to selected peak
		peaks_selected.push_back(peakPoints[i]);
	}//end loop peaks
	peakPoints.clear();
	peakPoints= peaks_selected;


	int nComponents= static_cast<int>(peakPoints.size());
	if(nComponents<=0){
		WARN_LOG("No components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<"), this should not occur (at least one should be found)!");	
		m_fitStatus= eFitAborted;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		return -1;
	}

	INFO_LOG("#"<<nComponents<<" components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<"), #"<<fitOptions.nMaxComponents<<" components will be fitted at maximum...");
	
	//## Sort peak fluxes and get sort indexes 
	//## NB: If maximum number of components is exceeded, limit fit to brightest components
	std::vector<double> componentPeakFluxes;
	for(int i=0;i<nComponents;i++){
		double x= peakPoints[i].X();
		double y= peakPoints[i].Y();
		long int gbin= fluxMapHisto->FindBin(x,y);
		if(fluxMapHisto->IsBinOverflow(gbin) || fluxMapHisto->IsBinUnderflow(gbin)){
			WARN_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
			m_fitStatus= eFitAborted;
			delete fluxMapHisto;
			fluxMapHisto= 0;	
			return -1;
		}
		double Speak= fluxMapHisto->GetBinContent(gbin);
		componentPeakFluxes.push_back(Speak);
	}//end loop components

	INFO_LOG("Sorting component peak fluxes from largest to smallest...");
	std::vector<size_t> sort_index;//sorting index
	std::vector<double> componentPeakFluxes_sorted;
	CodeUtils::sort_descending(componentPeakFluxes,componentPeakFluxes_sorted,sort_index);

	//## Initialize fit function with start parameters
	m_NFitComponents= std::min(nComponents,fitOptions.nMaxComponents);
	int nComponentPars= 6;
	int nFitPars= nComponentPars*m_NFitComponents + 1;//fit components + constant offset
	double fitRangeXmin= fluxMapHisto->GetXaxis()->GetXmin();
	double fitRangeXmax= fluxMapHisto->GetXaxis()->GetXmax();
	double fitRangeYmin= fluxMapHisto->GetYaxis()->GetXmin();
	double fitRangeYmax= fluxMapHisto->GetYaxis()->GetXmax();
	
	std::string parNamePrefix[6]= {"A","x0","y0","sigmaX","sigmaY","theta"};

	TF2* sourceFitFcn= new TF2("sourceFitFcn",
		SourceFitter::Gaus2DMixtureFcn,
		fitRangeXmin,fitRangeXmax,
		fitRangeYmin,fitRangeYmax,
		nFitPars
	);	
	sourceFitFcn->SetNpx(1000);
	sourceFitFcn->SetNpy(1000);	
	INFO_LOG("Created source fit function (nPars="<<nFitPars<<") with range: X["<<fitRangeXmin<<","<<fitRangeXmax<<"], Y["<<fitRangeYmin<<","<<fitRangeYmax<<"]");

	//## Set start fit pars
	int par_counter= 0;

	//- Offset
	double offset= 0.;
	if(fitOptions.fixBkg){//fix bkg level fit par
		if(fitOptions.useEstimatedBkgLevel) offset= bkgMean;//use estimated avg bkg
		else offset= fitOptions.fixedBkgLevel;//use user-supplied bkg level
		sourceFitFcn->SetParameter(nFitPars-1,offset);
		sourceFitFcn->FixParameter(nFitPars-1,offset);	
	}
	else{
		if(fitOptions.useEstimatedBkgLevel) offset= bkgMean;//use estimated avg bkg
		else offset= fitOptions.fixedBkgLevel;//use user-supplied bkg level
		if(bkgMean<Smin || bkgMean>Smax) offset= Smin;
		double offset_min= std::min(Smin,offset-fabs(rmsMean));
		double offset_max= std::min(Smax,offset+fabs(rmsMean));
		sourceFitFcn->SetParameter(nFitPars-1,offset);
		if(fitOptions.limitBkgInFit) {
			sourceFitFcn->SetParLimits(nFitPars-1,offset_min,offset_max);
			INFO_LOG("offset="<<offset<<" ["<<offset_min<<","<<offset_max<<"]");	
		}
	}
	sourceFitFcn->SetParName(nFitPars-1,"offset");
	

	for(int i=0;i<m_NFitComponents;i++){
		size_t index= sort_index[i];

		double x= peakPoints[index].X();
		double y= peakPoints[index].Y();
		long int gbin= fluxMapHisto->FindBin(x,y);
		double Speak= fluxMapHisto->GetBinContent(gbin);

		//Subtract offset from peak
		Speak-= offset;
		
		//## Set i-th component parameters
		//- Amplitude
		double Speak_min= std::min(Smin, Speak*(1 - fitOptions.amplLimit) );
		double Speak_max= std::max(Smax, Speak*(1 + fitOptions.amplLimit) );
		sourceFitFcn->SetParName(par_counter,Form("%s_%d",parNamePrefix[0].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter,Speak);
		if(fitOptions.limitAmplInFit){
			sourceFitFcn->SetParLimits(par_counter,Speak_min,Speak_max);
			INFO_LOG("Speak="<<Speak<<" ["<<Speak_min<<","<<Speak_max<<"]");
		}

		//- Centroids
		double centroidLimit= 0.5 * sqrt(pow(fitOptions.bmaj,2) + pow(fitOptions.bmin,2));
		double x0_min= x - centroidLimit;
		double x0_max= x + centroidLimit;
		double y0_min= y - centroidLimit;
		double y0_max= y + centroidLimit;
		sourceFitFcn->SetParName(par_counter+1,Form("%s_%d",parNamePrefix[1].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+1,x);
		sourceFitFcn->SetParName(par_counter+2,Form("%s_%d",parNamePrefix[2].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+2,y);
		if(fitOptions.limitCentroidInFit){
			sourceFitFcn->SetParLimits(par_counter+1,x0_min,x0_max);	
			sourceFitFcn->SetParLimits(par_counter+2,y0_min,y0_max);
			INFO_LOG("(x,y)=("<<x<<","<<y<<")"<<" bounds x("<<x0_min<<","<<x0_max<<") y("<<y0_min<<","<<y0_max<<")");
		}

		
		//- Sigmas
		double sigmaX= fitOptions.bmaj/GausSigma2FWHM;
		double sigmaY= fitOptions.bmin/GausSigma2FWHM;
		double sourceSigmaMax_x= fabs(Xmax-Xmin)*sqrt(2.)/GausSigma2FWHM;
		double sourceSigmaMax_y= fabs(Ymax-Ymin)*sqrt(2.)/GausSigma2FWHM;
		double sigmaX_min= sigmaX*(1-fitOptions.sigmaLimit);
		double sigmaY_min= sigmaY*(1-fitOptions.sigmaLimit);
		double sigmaX_max= std::max(sourceSigmaMax_x,sigmaX*(1+fitOptions.sigmaLimit));
		double sigmaY_max= std::max(sourceSigmaMax_y,sigmaY*(1+fitOptions.sigmaLimit));
		sourceFitFcn->SetParName(par_counter+3,Form("%s_%d",parNamePrefix[3].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+3,sigmaX);		
		sourceFitFcn->SetParName(par_counter+4,Form("%s_%d",parNamePrefix[4].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+4,sigmaY);

		if(fitOptions.limitSigmaInFit){
			sourceFitFcn->SetParLimits(par_counter+3,sigmaX_min,sigmaX_max);
			sourceFitFcn->SetParLimits(par_counter+4,sigmaY_min,sigmaY_max);
			INFO_LOG("(sigmaX,sigmaY)=("<<sigmaX<<","<<sigmaY<<")"<<" bounds sigmaX("<<sigmaX_min<<","<<sigmaX_max<<") y("<<sigmaY_min<<","<<sigmaY_max<<"), bmaj="<<fitOptions.bmaj<<", bmin="<<fitOptions.bmin);
		}

		//- Theta		
		double theta= fitOptions.bpa;
		double theta_min= theta - fitOptions.thetaLimit;
		double theta_max= theta + fitOptions.thetaLimit;
		sourceFitFcn->SetParName(par_counter+5,Form("%s_%d",parNamePrefix[5].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+5,theta);
		if(fitOptions.limitThetaInFit){
			sourceFitFcn->SetParLimits(par_counter+5,theta_min,theta_max);
			INFO_LOG("theta="<<theta<<" bounds ("<<theta_min<<","<<theta_max<<")");
		}

		//Update par counter
		par_counter+= nComponentPars;
	
	}//end loop components

	
	//==============================================
	//==             PRE-FIT
	//==============================================
	//## Fix sigmas to blob size?
	if(fitOptions.fixSigmaInPreFit){
		for(int i=0;i<m_NFitComponents;i++){
			//Sigma X
			TString parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			double parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);

			//Sigma Y
			parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
			parNumber= sourceFitFcn->GetParNumber(parName);
			parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);
		}
	}//close if fix sigma

	//## Fix theta to beam bpa 
	for(int i=0;i<m_NFitComponents;i++){
		TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
		int parNumber= sourceFitFcn->GetParNumber(parName);
		double theta= fitOptions.bpa;
		sourceFitFcn->FixParameter(parNumber,theta);
	}
	
	//## Fix offset 
	sourceFitFcn->FixParameter(nFitPars-1,offset);

	//## Perform pre-fit
	int prefitStatus= fluxMapHisto->Fit(sourceFitFcn,"RWN");
	if(prefitStatus!=0){
		WARN_LOG("Source pre-fit failed or did not converge.");
		m_fitStatus= eFitNotConverged;
	}

	//==============================================
	//==             FIT
	//==============================================
	//## Fix sigmas?
	if(fitOptions.fixSigma){
		for(int i=0;i<m_NFitComponents;i++){
			//Sigma X
			TString parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			double parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);

			//Sigma Y
			parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
			parNumber= sourceFitFcn->GetParNumber(parName);
			parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);
		}
	}//close if fix sigma
	else{
		for(int i=0;i<m_NFitComponents;i++){
			//Sigma X
			TString parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			sourceFitFcn->ReleaseParameter(parNumber);

			//Sigma Y
			parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
			parNumber= sourceFitFcn->GetParNumber(parName);
			sourceFitFcn->ReleaseParameter(parNumber);
		}
	}//close else

	//## Fix theta?
	if(fitOptions.fixTheta){
		for(int i=0;i<m_NFitComponents;i++){
			TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			double parValue= sourceFitFcn->GetParameter(parName);
			sourceFitFcn->FixParameter(parNumber,parValue);
		}
	}//close if fix theta
	else { 
		for(int i=0;i<m_NFitComponents;i++){
			TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
			int parNumber= sourceFitFcn->GetParNumber(parName);
			sourceFitFcn->ReleaseParameter(parNumber);
		}
	}//close else
	
	//## Release offset?
	if(!fitOptions.fixBkg){
		sourceFitFcn->ReleaseParameter(nFitPars-1);
	}

	//## Fit again
	//TFitResultPtr fitRes= fluxMapHisto->Fit(sourceFitFcn,"S");
	fluxMapHisto->Fit(sourceFitFcn,"SRWN");
	TBackCompFitter* fitter = (TBackCompFitter *) TVirtualFitter::GetFitter();
	const ROOT::Fit::FitResult& fitRes = fitter->GetFitResult();
	
	//==============================================
	//==             FIT RESULTS
	//==============================================
	//## Retrieve fit results
	m_fitStatus= eFitNotConverged;
	int fitMinimizerStatus= -1;
	int fitStatus= -1;
	bool hasParsAtLimits= false;
	
	if(fitRes.IsValid()) {
		fitStatus= 0;
		m_fitStatus= eFitConverged;
	}
	
	//NB: Minuit fit status = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult (4: MIGRAD error, 40: MINOS error, 400: HESSE error, 4000: IMPROVE error)
	fitMinimizerStatus= fitRes.Status();

	//Check if has parameters at bounds
	hasParsAtLimits= HasFitParsAtLimit(fitRes);
	if(hasParsAtLimits) m_fitStatus= eFitConvergedWithWarns;
	

	double Chi2= sourceFitFcn->GetChisquare();//fitRes->Chi2();
	double NDF= sourceFitFcn->GetNDF();//fitRes->Ndf();
	int NFreePars= sourceFitFcn->GetNumberFreeParameters();//fitRes->NFreeParameters();
	int NFittedBins= sourceFitFcn->GetNumberFitPoints();//fitRes->NFreeParameters();
	INFO_LOG("Source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") fit info: status="<<fitStatus<<" fitterStatus="<<m_fitStatus<<" (hasParsAtLimits? "<<hasParsAtLimits<<"), NDF="<<NDF<<", Chi2="<<Chi2<<", NFreePars="<<NFreePars<<" NFittedBins="<<NFittedBins);

	m_sourceFitPars.SetNComponents(m_NFitComponents);
	m_sourceFitPars.SetChi2(Chi2);
	m_sourceFitPars.SetNDF(NDF);
	m_sourceFitPars.SetNFreePars(NFreePars);
	m_sourceFitPars.SetNFitPoints(NFittedBins);
	m_sourceFitPars.SetStatus(m_fitStatus);
	m_sourceFitPars.SetMinimizerStatus(fitMinimizerStatus);
	
	par_counter= 0;

	for(int i=0;i<m_NFitComponents;i++){
		for(int j=0;j<nComponentPars;j++){
			std::string parName_global= std::string(sourceFitFcn->GetParName(par_counter));
			std::string parName= CodeUtils::ExtractSubString(parName_global,"_");
			double parVal= sourceFitFcn->GetParameter(par_counter);
			double parErr= sourceFitFcn->GetParError(par_counter);
			if(m_sourceFitPars.SetParValueAndError(i,parName,parVal,parErr)<0){
				WARN_LOG("Failed to set par "<<parName<<" (parName_global="<<parName_global<<") value and error (check par name)!");
			}
			par_counter++;
		}//end loop pars per component
	}//end loop components
	
	double fittedOffset= sourceFitFcn->GetParameter(par_counter);
	double fittedOffsetErr= sourceFitFcn->GetParError(par_counter);
	m_sourceFitPars.SetOffsetPar(fittedOffset);
	m_sourceFitPars.SetOffsetParErr(fittedOffsetErr);

	//==============================================
	//==             FIT RESIDUALS
	//==============================================
	//Compute fit residuals	
	std::vector<double> residuals;
	double residualMax= -1.e+99;
	double residualMin= 1.e+99;
	for (int i=0;i<fluxMapHisto->GetNbinsX();i++) {
		double x= fluxMapHisto->GetXaxis()->GetBinCenter(i+1);
		for (int j=0;j<fluxMapHisto->GetNbinsY();j++) {
			double y= fluxMapHisto->GetYaxis()->GetBinCenter(j+1);
			double data= fluxMapHisto->GetBinContent(i+1,j+1);
			if(!std::isnormal(data) || data==0) continue;
			double model= sourceFitFcn->Eval(x,y);
  		double res = data - model;
			if(res<residualMin) residualMin= res;
			if(res>residualMax) residualMax= res;
			residuals.push_back(res);
		}
  }
	
	//Compute residual stats
	double residualMean= 0;
	double residualRMS= 0;
	StatsUtils::ComputeMeanAndRMS(residualMean,residualRMS,residuals);
	double residualMedian= StatsUtils::GetMedianFast<double>(residuals);
	double residualMAD= StatsUtils::GetMADFast(residuals,residualMedian);	
	
	m_sourceFitPars.SetResidualMin(residualMin);
	m_sourceFitPars.SetResidualMax(residualMax);
	m_sourceFitPars.SetResidualMean(residualMean);
	m_sourceFitPars.SetResidualRMS(residualRMS);
	m_sourceFitPars.SetResidualMedian(residualMedian);
	m_sourceFitPars.SetResidualMAD(residualMAD);
	
	INFO_LOG("Fit residual stats: min/max="<<residualMin<<"/"<<residualMax<<", mean="<<residualMean<<", rms="<<residualRMS<<", median="<<residualMedian<<", mad="<<residualMAD);
	
	//Clear up data
	if(fluxMapHisto){
		delete fluxMapHisto;
		fluxMapHisto= 0;
	}
	if(sourceFitFcn){
		delete sourceFitFcn;
		sourceFitFcn= 0;
	}
	return 0;

}//close FitSource()
*/

bool SourceFitter::HasFitParsAtLimit(const ROOT::Fit::FitResult& fitRes)
{
	//Loop over parameters and check if they converged at bounds
	bool hasParAtLimits= false;
	double parAtLimitThr= 1.e-6;//in percentage
	
	for(unsigned int i=0;i<fitRes.NPar();i++){
		if(!fitRes.IsParameterBound(i) || fitRes.IsParameterFixed(i) ) continue;
		
		std::string parName= fitRes.GetParameterName(i);

		//Check par limits
		double par_min= 0;
		double par_max= 0;
		fitRes.ParameterBounds(i,par_min,par_max);
		double par= fitRes.Parameter(i);		
		double parRelDiffMin= fabs(par/par_min-1);
		double parRelDiffMax= fabs(par/par_max-1);

		//if(par==par_min || par==par_max){
		if(parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr){	
			hasParAtLimits= true;
			INFO_LOG("Fit parameter "<<parName<<" is at limit (value="<<par<<", min/max="<<par_min<<"/"<<par_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
			break;
		}
	}//end loop parameters
	
	return hasParAtLimits;

}//close HasFitParsAtLimit()




double SourceFitter::Gaus2DMixtureFcn(double* x, double* p)
{
	//Define model pars
	//[1-7]: Gaussian component i-th pars
	//[last]: offset
	
	int nPars= 6;
	double pdfPars[nPars];
	double fcnValue= 0.;
	int par_counter= 0;
	
	for(int k=0;k<m_NFitComponents;k++){
		for(int j=0;j<nPars;j++){
			pdfPars[j]= p[par_counter];
			par_counter++;
		}		
		fcnValue+= Gaus2DFcn(x,pdfPars);
	}//end loop mixtures
	
	double offset= p[par_counter];
	fcnValue+= offset;

	return fcnValue;

}//close Gaus2DMixtureFcn()


double SourceFitter::Gaus2DFcn(double* x, double* par){

	double X= x[0];
	double Y= x[1];
	
	double A= par[0];
	double X0= par[1];
	double Y0= par[2];
	double sigmaX= par[3];
	double sigmaY= par[4];
	double theta= par[5];
	theta*= TMath::DegToRad();

	/*
	double sigmaX2= sigmaX*sigmaX;
	double sigmaY2= sigmaY*sigmaY;

	double cosTheta= cos(theta);
	double cos2Theta= cos(2*theta);
	double cosTheta2= cosTheta*cosTheta;
	double sinTheta= sin(theta);
	double sin2Theta= sin(2*theta);
	double sinTheta2= sinTheta*sinTheta;

	double a= cosTheta2/(2*sigmaX2) + sinTheta2/(2*sigmaY2);
	//double b= -sin2Theta/(4*sigmaX2) + sin2Theta/(4*sigmaY2);
	double b= sin2Theta/(4*sigmaX2) - sin2Theta/(4*sigmaY2);
	double c= sinTheta2/(2*sigmaX2) + cosTheta2/(2*sigmaY2);
	double argX= a*(X-X0)*(X-X0);
	double argY= c*(Y-Y0)*(Y-Y0);
	double argXY= 2*b*(X-X0)*(Y-Y0);
	double z= argX+argXY+argY;

	double fcn= A*exp(-z);
	*/

	
	double cost2= cos(theta)*cos(theta);
	double sint2= sin(theta)*sin(theta);
	double sin2t = sin(2. * theta);
	double xstd2= sigmaX*sigmaX;
	double ystd2= sigmaY*sigmaY;
	double xdiff = X - X0;
  double ydiff = Y - Y0;
	double a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2));
  double b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2));
  double c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2));
	double fcn= A * exp(-((a * xdiff*xdiff) + (b * xdiff * ydiff) + (c * ydiff*ydiff)));

	return fcn;

}//close Gaus2DFcn()


}//close namespace

