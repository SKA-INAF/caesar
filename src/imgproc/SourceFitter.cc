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

#include <TMinuitMinimizer.h>
#ifdef MINUIT2_ENABLED
	#include <Minuit2/Minuit2Minimizer.h>
#endif
#ifdef ROOTR_ENABLED
	#include <Math/RMinimizer.h>
#endif
#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/GSLNLSMinimizer.h>
#include <Math/Functor.h>
#include <Math/Factory.h>
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
using namespace std::placeholders;

ClassImp(Caesar::SourceComponentPars)
ClassImp(Caesar::SourceFitPars)
ClassImp(Caesar::SourceFitter)
ClassImp(Caesar::SourceFitter::SourceFitData)

namespace Caesar {

//Static variables
/*
int SourceFitter::m_NFitComponents= 0;
SourceFitPars SourceFitter::m_sourceFitPars;
int SourceFitter::m_fitStatus= eFitUnknownStatus;
//TH2D* SourceFitter::m_fluxMapHisto= 0;
std::vector<SourceFitter::SourceFitData> SourceFitter::m_fitData;
std::vector<SourceFitter::SourceFitData> SourceFitter::m_fitHaloData;
double SourceFitter::m_bkgMean= 0;
double SourceFitter::m_rmsMean= 0;
double SourceFitter::m_sourceX0= 0;
double SourceFitter::m_sourceY0= 0;
double SourceFitter::m_chi2RegPar= 0.1;
*/

//Constructor
SourceFitter::SourceFitter()
	: TObject()
{
	//m_fluxMapHisto= 0;
	m_bkgMean= 0.;
	m_rmsMean= 0.;
	m_fitData.clear();
	m_fitHaloData.clear();

}//close costructor

//Destructor
SourceFitter::~SourceFitter()
{
	//Delete flux map histo
	/*
	if(m_fluxMapHisto){
		delete m_fluxMapHisto;
		m_fluxMapHisto= 0;
	}
	*/

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
	m_sourceX0= aSource->GetSx();
	m_sourceY0= aSource->GetSy();

	
	//## Initialize flux histo/map
	/*
	if(m_fluxMapHisto) {
		delete m_fluxMapHisto;
		m_fluxMapHisto= 0;
	}
	m_fluxMapHisto= new TH2D("","",nBoxX,Xmin-0.5,Xmax+0.5,nBoxY,Ymin-0.5,Ymax+0.5);
	m_fluxMapHisto->SetDirectory(0);//disconnect from any existing file open
	m_fluxMapHisto->Sumw2();
	*/

	//## Fill source flux histo
	INFO_LOG("Filling flux histo for source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") ...");
	m_bkgMean= 0.;
	m_rmsMean= 0.;
	long int ndata= 0;
	long int ndata_rms= 0;
	//std::vector<int> filledBinIds;
	m_fitData.clear();
	m_fitHaloData.clear();

	for(size_t k=0;k<pixels.size();k++){
		double x= pixels[k]->x;
		double y= pixels[k]->y;
		double S= pixels[k]->S;
		if(!std::isnormal(S)) continue;
		std::pair<double,double> bkgData= pixels[k]->GetBkg();	
		double bkg= bkgData.first;
		double rms= bkgData.second;
	
		//Skip bin if below significance level?
		/*
		if(rms>0 && fitOptions.useFluxZCut){
			double Z= (S-bkg)/rms;		
			if(Z<fitOptions.fluxZThrMin) {
				DEBUG_LOG("Skipping pixel no. "<<k<<" (S="<<S<<", Z="<<Z<<"<"<<fitOptions.fluxZThrMin<<") as below desired significance level...");
				continue;
			}
		}
		*/

		//Find histo bin being filled	
		/*
		int binId= m_fluxMapHisto->FindBin(x,y);
		if(m_fluxMapHisto->IsBinUnderflow(binId) || m_fluxMapHisto->IsBinOverflow(binId)){
			WARN_LOG("Cannot find histo bin corresponding to pixel flux (S="<<S<<"), skip pixel...");
			continue;
		}
		*/

		//Compute pixel flux error
		//Set to 1 by default (e.g. ignored in fitting)
		double S_err= 1;
		if(fitOptions.setBinErrorsToMapRMS && rms>0){
			S_err= rms;
		}

		
		//Fill histogram bin and associated error
		/*
		m_fluxMapHisto->Fill(x,y,S);//increment sum of weights
		m_fluxMapHisto->SetBinError(binId,S_err);
		filledBinIds.push_back(binId);
		*/

		//Scale and fill fit data
		S*= 1.e+3;//convert to mJy/beam
		S_err*= 1.e+3;//convert to mJy/beam
		x-= m_sourceX0;
		y-= m_sourceY0;
		if(rms>0 && fitOptions.useFluxZCut){
			double Z= (S-bkg)/rms;		
			if(Z<fitOptions.fluxZThrMin) {
				m_fitHaloData.push_back(SourceFitData(x,y,S,S_err));
			}
			else{
				m_fitData.push_back(SourceFitData(x,y,S,S_err));
			}
		}//close if
		else{
			m_fitData.push_back(SourceFitData(x,y,S,S_err));
		}

		//Update number of fit entries & bkg/rms accumulators
		if(rms>0){
			m_bkgMean+= bkg;
			m_rmsMean+= rms;
			ndata_rms++;
		}
		ndata++;

	}//end loop pixels

	if(ndata_rms>0) {
		m_bkgMean/= (double)(ndata_rms);
		m_rmsMean/= (double)(ndata_rms);
	}

	//Set all bins error to average noise
	/*
	if(fitOptions.setBinErrorsToMapRMS && m_rmsMean>0){
		INFO_LOG("Setting fitted histo bin errors to average rms (<rms(Jy)>="<<m_rmsMean<<")");
		for(size_t i=0;i<filledBinIds.size();i++){
			int binId= filledBinIds[i];
			m_fluxMapHisto->SetBinError(binId,m_rmsMean);
		}
	}//close if
	*/

	//Convert bkg & rms to mJy
	m_bkgMean*= 1.e+3;//convert to mJy
	m_rmsMean*= 1.e+3;//convert to mJy

	INFO_LOG("Source (name="<<aSource->GetName()<<", N="<<pixels.size()<<", pos("<<m_sourceX0<<","<<m_sourceY0<<")) bkg info: <bkg(mJy)>="<<m_bkgMean<<", <rms(mJy)>="<<m_rmsMean);

	//Check if histo has entries
	/*
	if(m_fluxMapHisto->GetSumOfWeights()<=0){
		ERROR_LOG("No entries present in flux map histogram (hint: check if properly filled, e.g. bin range problems)!");
		return -1;
	}
	*/

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

	if(nComponents==0 || !fitOptions.useNestedAsComponents){
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


		//NB: Try to extract blended blobs to estimate component pars. 
		//    If this fails estimate number of components from detected peaks and assuming beam pars
		std::vector<Source*> deblendedBlobs;
		std::vector<ImgPeak> peaks;
		bool useDeblendedBlobPars= false;

		int status= aSource->FindBlendedComponents(	
			deblendedBlobs,peaks,
			fitOptions.peakZThrMin, fitOptions.nMaxComponents,
			fitOptions.scaleMin,fitOptions.scaleMax,fitOptions.scaleStep,
			fitOptions.minBlobSize,fitOptions.blobMapThrFactor,fitOptions.blobMapKernelFactor
		);
		if(status==0){
			INFO_LOG("#"<<deblendedBlobs.size()<<" blended component found in source "<<aSource->GetName()<<" ...");
			useDeblendedBlobPars= true;
		}//close if
		else{
			WARN_LOG("Failed to find blended blobs in source "<<aSource->GetName()<<", will estimate component peaks only and assume beam pars ...");
			status= aSource->FindComponentPeaks(
				peaks,
				fitOptions.peakZThrMin, fitOptions.nMaxComponents,
				fitOptions.peakShiftTolerance,
				kernels,fitOptions.peakKernelMultiplicityThr
			);
			if(status<0){
				ERROR_LOG("Failed to find component peaks in source (name="<<aSource->GetName()<<")!");
				return -1;
			}
		}//close else

	
		//NB: If only one peak found initialize fit component to the entire source
		//otherwise use peaks found as centroid start values and beam as gaussian sigma pars
		if(peaks.empty()){
			INFO_LOG("No peaks found in this source (hint: could be a diffuse source and peaks are below desired threshold), nothing to be fit");
			return -1;
		}
		if(peaks.size()==1){
			fitPars_start.push_back( std::vector<double>() );
			
			fitPars_start[0].push_back(Smax*1.e+3);//converted in mJy
			fitPars_start[0].push_back(meanX-m_sourceX0);//normalized to centroid
			fitPars_start[0].push_back(meanY-m_sourceY0);//normalized to centroid
			fitPars_start[0].push_back(sigmaX);
			fitPars_start[0].push_back(sigmaY);
			fitPars_start[0].push_back(theta*TMath::DegToRad());//converted to rad
			
		}
		else{
			for(size_t i=0;i<peaks.size();i++){
				double x= peaks[i].x;
				double y= peaks[i].y;
				double Speak= peaks[i].S;
				double sigmaX= fitOptions.bmaj/GausSigma2FWHM;
				double sigmaY= fitOptions.bmin/GausSigma2FWHM;
				double theta= fitOptions.bpa;
					
				if(useDeblendedBlobPars && deblendedBlobs[i]){
					INFO_LOG("Computing pars of blended component no. "<<i+1<<" (nPix="<<deblendedBlobs[i]->GetNPixels()<<") for source "<<aSource->GetName()<<" ...");
					double sigmaX_sample_blended= 0;
					double sigmaY_sample_blended= 0;
					double covXY_sample_blended= 0;
					if(deblendedBlobs[i]->GetSampleStdDev(sigmaX_sample_blended,sigmaY_sample_blended,covXY_sample_blended)==0){
						sigmaX= sigmaX_sample_blended;
						sigmaY= sigmaY_sample_blended;
						theta= StatsUtils::GetGaus2DThetaPar(sigmaX_sample_blended,sigmaY_sample_blended,covXY_sample_blended);		
						StatsUtils::GetEllipseParsFromCovMatrix(sigmaX,sigmaY,theta,sigmaX_sample_blended,sigmaY_sample_blended,covXY_sample_blended);
					}
					else{
						WARN_LOG("Failed to stddev pars of blended component no. "<<i+1<<" for source "<<aSource->GetName()<<", will use beam info...");
					}
				}//close if use deblended blob pars

				fitPars_start.push_back( std::vector<double>() );
				fitPars_start[i].push_back(Speak*1.e+3);//converted in mJy
				fitPars_start[i].push_back(x-m_sourceX0);//normalized to centroid
				fitPars_start[i].push_back(y-m_sourceY0);//normalized to centroid
				fitPars_start[i].push_back(sigmaX);
				fitPars_start[i].push_back(sigmaY);
				fitPars_start[i].push_back(theta*TMath::DegToRad());//converted to rad
				
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
			if(Zpeak_imgbkg<fitOptions.peakZThrMin) {
			//if(Zpeak_sourcebkg<fitOptions.peakZThrMin) {
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
			std::vector<ImgPeak> peaks;
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
				fitPars_start[componentCounter].push_back(Smax*1.e+3);//converted to mJy
				fitPars_start[componentCounter].push_back(meanX-m_sourceX0);//normalized to centroid
				fitPars_start[componentCounter].push_back(meanY-m_sourceY0);//normalized to centroid
				fitPars_start[componentCounter].push_back(sigmaX);
				fitPars_start[componentCounter].push_back(sigmaY);
				fitPars_start[componentCounter].push_back(theta*TMath::DegToRad());//converted to rad
				componentCounter++;
			}//close if

			else{
				//Loop over peaks and add component
				for(size_t j=0;j<peaks.size();j++){
					double x= peaks[j].x;
					double y= peaks[j].y;
					double Speak= peaks[j].S;
					/*
					long int gbin= m_fluxMapHisto->FindBin(x,y);
					if(m_fluxMapHisto->IsBinOverflow(gbin) || m_fluxMapHisto->IsBinUnderflow(gbin)){
						WARN_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
						return -1;
					}
					double Speak= m_fluxMapHisto->GetBinContent(gbin);
					*/
					double sigmaX= fitOptions.bmaj/GausSigma2FWHM;
					double sigmaY= fitOptions.bmin/GausSigma2FWHM;
					double theta= fitOptions.bpa;
	
					fitPars_start.push_back( std::vector<double>() );	
					fitPars_start[componentCounter].push_back(Speak*1.e+3);//converted to mJy
					fitPars_start[componentCounter].push_back(x-m_sourceX0);//normalized to centroid
					fitPars_start[componentCounter].push_back(y-m_sourceY0);//normalized to centroid
					fitPars_start[componentCounter].push_back(sigmaX);
					fitPars_start[componentCounter].push_back(sigmaY);
					fitPars_start[componentCounter].push_back(theta*TMath::DegToRad());//converted to rad	
					componentCounter++;
				}//end loop peaks
			}//close else
	
		}//end loop selected nested sources
		
	}//close else

	INFO_LOG("#"<<fitPars_start.size()<<" fit components will be fitted...");

	
	return 0;

}//close EstimateFitComponents()


/*
int SourceFitter::DoChi2Fit(Source* aSource,SourceFitOptions& fitOptions,std::vector< std::vector<double> >& fitPars_start)
{
	//## Get source pars
	float Xmin, Xmax, Ymin, Ymax;
	aSource->GetSourceRange(Xmin,Xmax,Ymin,Ymax);
	double Smax= aSource->GetSmax();
	double Smin= aSource->GetSmin();

	double Xmin_norm= Xmin-m_sourceX0;
	double Xmax_norm= Xmax-m_sourceX0;
	double Ymin_norm= Ymin-m_sourceY0;
	double Ymax_norm= Ymax-m_sourceY0;

	//Convert to mJy
	Smax*= 1.e+3;
	Smin*= 1.e+3;


	//==============================================
	//==             INITIALIZE FITTER
	//==============================================
	//## Initialize fitter
	m_NFitComponents= static_cast<int>(fitPars_start.size());
	int nComponentPars= 6;
	int nFitPars= nComponentPars*m_NFitComponents + 1;//fit components + constant offset

	//TVirtualFitter::SetDefaultFitter("Minuit");
	//TVirtualFitter* minuit= TVirtualFitter::Fitter(0,nFitPars);
	TFitter* minuit= new TFitter(nFitPars);
	minuit->SetDefaultFitter("Minuit");
	minuit->SetErrorDef(1);//set chi2 error definition (0.5 for likelihood fits)
	minuit->SetMaxIterations(5000);//default niters
	minuit->SetPrecision(1.e-6);//default tolerance

	
	//## Set start fit pars
	int par_counter= 0;
	std::string parNamePrefix[]= {"A","x0","y0","sigmaX","sigmaY","theta"};
	std::vector<std::vector<double>> parLimits_min;
	std::vector<std::vector<double>> parLimits_max;

	//- Offset
	double offset= 0.;
	double offset_min= 0;
	double offset_max= 0;
	double offset_step= 0.1;
	double offset_err= 0.;
	if(fitOptions.useEstimatedBkgLevel) offset= m_bkgMean;//use estimated avg bkg
	//else offset= fitOptions.fixedBkgLevel;//use user-supplied bkg level
	else offset= fitOptions.fixedBkgLevel*1.e+3;//use user-supplied bkg level converted in mJy
	if(m_bkgMean>Smax) {
		WARN_LOG("Offset start par given/estimated is below min or above max source flux, setting it to Smin...");
		offset= Smin;
	}
	offset_err= fabs(offset_step*offset);
	if(fitOptions.limitBkgInFit){
		offset_min= std::min(Smin,offset-fabs(m_rmsMean));
		offset_max= std::min(Smax,offset+fabs(m_rmsMean));
		INFO_LOG("Limiting offset par "<<offset<<" in range ["<<offset_min<<","<<offset_max<<"]");
	}
	
	minuit->SetParameter(nFitPars-1,"offset",offset,offset_err,offset_min,offset_max);
	minuit->FixParameter(nFitPars-1);//Fix offset par for pre-fit
	
	//- Component pars
	double Speak_step= 0.01;
	double centroid_step= 0.01;
	double sigma_step= 0.01;
	double theta_step= 0.01;
	std::map<std::string,int> parNameMap;

	for(int i=0;i<m_NFitComponents;i++){
		//Init min/max limits
		parLimits_min.push_back(std::vector<double>(nComponentPars,0));
		parLimits_max.push_back(std::vector<double>(nComponentPars,0));
		
		//## Set i-th component parameters
		//- Amplitude
		double Speak= fitPars_start[i][0] - offset;//subtract bkg offset from peak
		double Speak_min= 0;
		double Speak_max= 0;
		double Speak_err= fabs(Speak*Speak_step);
		TString parName= Form("%s_%d",parNamePrefix[0].c_str(),i+1);
		if(fitOptions.limitAmplInFit){
			Speak_min= std::min(Smin, Speak*(1 - fitOptions.amplLimit) );
			Speak_max= std::max(Smax, Speak*(1 + fitOptions.amplLimit) );
			parLimits_min[i][0]= Speak_min;
			parLimits_max[i][0]= Speak_max;
			INFO_LOG("Limiting amplitude par "<<Speak<<" in range ["<<Speak_min<<","<<Speak_max<<"]");
		}
		parNameMap.insert( std::make_pair(std::string(parName),par_counter) );
		minuit->SetParameter(par_counter,parName,Speak,Speak_err,Speak_min,Speak_max);
		if(fitOptions.fixAmplInPreFit){
			minuit->FixParameter(par_counter);
		}

		//- Centroids
		double x0= fitPars_start[i][1];
		double y0= fitPars_start[i][2];
		//double x0_err= fabs(x0*centroid_step);
		//double y0_err= fabs(y0*centroid_step);
		double x0_err= 1;//1 pixel around centroid
		double y0_err= 1;
		//double centroidLimit= 0.5 * sqrt(pow(fitOptions.bmaj,2) + pow(fitOptions.bmin,2));//in pixels
		//double centroidLimit= 0.5*fitOptions.bmaj;
		double centroidLimit= fitOptions.centroidLimit;
		double x0_min= 0;
		double x0_max= 0;
		double y0_min= 0;
		double y0_max= 0;
		if(fitOptions.limitCentroidInFit){
			//x0_min= x0 - centroidLimit;
			//x0_max= x0 + centroidLimit;
			//y0_min= y0 - centroidLimit;
			//y0_max= y0 + centroidLimit;
			x0_min= std::max(x0 - centroidLimit,Xmin_norm);
			x0_max= std::min(x0 + centroidLimit,Xmax_norm);
			y0_min= std::max(y0 - centroidLimit,Ymin_norm);
			y0_max= std::min(y0 + centroidLimit,Ymax_norm);
			parLimits_min[i][1]= x0_min;
			parLimits_max[i][1]= x0_max;
			parLimits_min[i][2]= y0_min;
			parLimits_max[i][2]= y0_max;
			INFO_LOG("Limiting centroid pars (x,y)=("<<x0<<","<<y0<<")"<<" to bounds x0("<<x0_min<<","<<x0_max<<"), y0("<<y0_min<<","<<y0_max<<")");
		}
		parName= Form("%s_%d",parNamePrefix[1].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+1) );
		minuit->SetParameter(par_counter+1,parName,x0,x0_err,x0_min,x0_max);
		parName= Form("%s_%d",parNamePrefix[2].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+2) );
		minuit->SetParameter(par_counter+2,parName,y0,y0_err,y0_min,y0_max);
		if(fitOptions.fixCentroidInPreFit){
			minuit->FixParameter(par_counter+1);		
			minuit->FixParameter(par_counter+2);
		}			
	

		//- Sigmas
		double sigmaX= fitPars_start[i][3];
		double sigmaY= fitPars_start[i][4];
		double sourceSigmaMax_x= fabs(Xmax-Xmin)*sqrt(2.)/GausSigma2FWHM;
		double sourceSigmaMax_y= fabs(Ymax-Ymin)*sqrt(2.)/GausSigma2FWHM;
		double sigmaX_err= fabs(sigmaX*sigma_step);
		double sigmaY_err= fabs(sigmaY*sigma_step);
		double sigmaX_min= 0;
		double sigmaY_min= 0;
		double sigmaX_max= 0;
		double sigmaY_max= 0;
		if(fitOptions.limitSigmaInFit){
			sigmaX_min= std::max(0.,sigmaX*(1-fitOptions.sigmaLimit));
			sigmaY_min= std::max(0.,sigmaY*(1-fitOptions.sigmaLimit));
			sigmaX_max= std::max(sourceSigmaMax_x,sigmaX*(1+fitOptions.sigmaLimit));
			sigmaY_max= std::max(sourceSigmaMax_y,sigmaY*(1+fitOptions.sigmaLimit));
			parLimits_min[i][3]= sigmaX_min;
			parLimits_max[i][3]= sigmaX_max;
			parLimits_min[i][4]= sigmaY_min;
			parLimits_max[i][4]= sigmaY_max;
			INFO_LOG("Limiting sigma pars (sigmaX,sigmaY)=("<<sigmaX<<","<<sigmaY<<")"<<" to bounds sigmaX("<<sigmaX_min<<","<<sigmaX_max<<"), sigmaY("<<sigmaY_min<<","<<sigmaY_max<<"), bmaj="<<fitOptions.bmaj<<", bmin="<<fitOptions.bmin);
		}
		parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+3) );	
		minuit->SetParameter(par_counter+3,parName,sigmaX,sigmaX_err,sigmaX_min,sigmaX_max);
		parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+4) );
		minuit->SetParameter(par_counter+4,parName,sigmaY,sigmaY_err,sigmaY_min,sigmaY_max);
		if(fitOptions.fixSigmaInPreFit){//fix sigma in pre-fit?
			minuit->FixParameter(par_counter+3);		
			minuit->FixParameter(par_counter+4);		
		}

		//- Theta		
		double theta= fitPars_start[i][5];
		double theta_min= 0;
		double theta_max= 0;
		double theta_err= fabs(theta*theta_step);
		double theta_limits= fitOptions.thetaLimit*TMath::DegToRad();
		//double theta_limits= fitOptions.thetaLimit;	
		if(fitOptions.limitThetaInFit){
			theta_min= theta - theta_limits;
			theta_max= theta + theta_limits;
			parLimits_min[i][5]= theta_min;
			parLimits_max[i][5]= theta_max;
			INFO_LOG("Limiting theta par "<<theta<<" to bounds ("<<theta_min<<","<<theta_max<<")");
		}
		parName= Form("%s_%d",parNamePrefix[5].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+5) );
		minuit->SetParameter(par_counter+5,parName,theta,theta_err,theta_min,theta_max);
		if(fitOptions.fixThetaInPreFit){//Fix theta par in pre-fit?
			minuit->FixParameter(par_counter+5);	
		}
		
		//Update par counter
		par_counter+= nComponentPars;
	
	}//end loop components

	//## Set minimization function
	minuit->SetFCN(SourceFitter::Chi2Fcn);
	//minuit->SetFCN( &std::bind(&SourceFitter::Chi2Fcn,this,_1,_2,_3,_4,_5) );

	ROOT::Math::Minimizer* minuitMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
	ROOT::Math::Functor fitFcn(std::bind(&SourceFitter::Chi2Fcn_test, this, _1), nFitPars); 
	minuitMinimizer->SetFunction(fitFcn);

	//==============================================
	//==             PRE-FIT
	//==============================================
	//## Perform pre-fit
	double arglist[100];
  arglist[0] = 0;
      
	// Set print level
  minuit->ExecuteCommand("SET PRINT",arglist,2);
	
	// Minimize
  arglist[0] = fitOptions.fitMaxIters; // number of function calls
  arglist[1] = fitOptions.fitFcnTolerance; // tolerance
  int prefitStatus= minuit->ExecuteCommand("MIGRAD",arglist,2);
	if(prefitStatus!=0){
		WARN_LOG("Source pre-fit failed or did not converge, trying to perform the full fit...");
		m_fitStatus= eFitNotConverged;
	}

	//==============================================
	//==             RELEASE PARS
	//==============================================
	
	//## Release pars and adjust par limits around pre-fitted estimates
	for(int i=0;i<m_NFitComponents;i++){
		//Amplitudes
		TString parName= Form("%s_%d",parNamePrefix[0].c_str(),i+1);
		int parNumber= parNameMap[std::string(parName)];
		double Ampl_fitted= minuit->GetParameter(parNumber);
		double Ampl_fitted_min= Ampl_fitted - fabs(fitOptions.amplLimit*Ampl_fitted);
		double Ampl_fitted_max= Ampl_fitted + fabs(fitOptions.amplLimit*Ampl_fitted);
		double Ampl_fitted_err= fabs(Ampl_fitted*Speak_step);
		if(minuit->IsFixed(parNumber)) minuit->ReleaseParameter(parNumber);	
		if(fitOptions.limitAmplInFit){
			INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<Ampl_fitted<<") to ["<<Ampl_fitted_min<<","<<Ampl_fitted_max<<"]");
			minuit->SetParameter(parNumber,parName,Ampl_fitted,Ampl_fitted_err,Ampl_fitted_min,Ampl_fitted_max);
		}

		//Centroid X
		parName= Form("%s_%d",parNamePrefix[1].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double X0_fitted= minuit->GetParameter(parNumber);
		double X0_fitted_min_default= X0_fitted - fabs(fitOptions.centroidLimit);
		double X0_fitted_max_default= X0_fitted + fabs(fitOptions.centroidLimit);
		double X0_fitted_min= std::max(Xmin_norm,X0_fitted_min_default);
		double X0_fitted_max= std::min(Xmax_norm,X0_fitted_max_default);
		double X0_fitted_err= 1;
		if(minuit->IsFixed(parNumber)) minuit->ReleaseParameter(parNumber);	
		if(fitOptions.limitCentroidInFit){
			INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<X0_fitted<<") to ["<<X0_fitted_min<<","<<X0_fitted_max<<"]");
			minuit->SetParameter(parNumber,parName,X0_fitted,X0_fitted_err,X0_fitted_min,X0_fitted_max);
		}

		//Centroid Y
		parName= Form("%s_%d",parNamePrefix[2].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double Y0_fitted= minuit->GetParameter(parNumber);
		double Y0_fitted_min_default= Y0_fitted - fabs(fitOptions.centroidLimit);
		double Y0_fitted_max_default= Y0_fitted + fabs(fitOptions.centroidLimit);
		double Y0_fitted_min= std::max(Ymin_norm,Y0_fitted_min_default);
		double Y0_fitted_max= std::min(Ymax_norm,Y0_fitted_max_default);
		double Y0_fitted_err= 1;
		if(minuit->IsFixed(parNumber)) minuit->ReleaseParameter(parNumber);	
		if(fitOptions.limitCentroidInFit){
			INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<Y0_fitted<<") to ["<<Y0_fitted_min<<","<<Y0_fitted_max<<"]");
			minuit->SetParameter(parNumber,parName,Y0_fitted,Y0_fitted_err,Y0_fitted_min,Y0_fitted_max);
		}

		//Sigma X
		parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double sigmaX_fitted= minuit->GetParameter(parNumber);
		double sigmaX_fitted_err= fabs(sigmaX_fitted*sigma_step);
		double sigmaX_fitted_min_default= sigmaX_fitted - fabs(fitOptions.sigmaLimit*sigmaX_fitted);
		double sigmaX_fitted_max_default= sigmaX_fitted + fabs(fitOptions.sigmaLimit*sigmaX_fitted);
		double sigmaX_fitted_min= std::max(0.,sigmaX_fitted_min_default);
		double sigmaX_fitted_max= sigmaX_fitted_max_default;

		if(fitOptions.fixSigma){
			if(!minuit->IsFixed(parNumber)) minuit->FixParameter(parNumber);
		}
		else{
			if(minuit->IsFixed(parNumber)) minuit->ReleaseParameter(parNumber);
			if(fitOptions.limitSigmaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sigmaX_fitted<<") to ["<<sigmaX_fitted_min<<","<<sigmaX_fitted_max<<"]");		
				minuit->SetParameter(parNumber,parName,sigmaX_fitted,sigmaX_fitted_err,sigmaX_fitted_min,sigmaX_fitted_max);
			}
		}

		//SigmaY
		parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double sigmaY_fitted= minuit->GetParameter(parNumber);
		double sigmaY_fitted_err= fabs(sigmaY_fitted*sigma_step);
		double sigmaY_fitted_min_default= sigmaY_fitted - fabs(fitOptions.sigmaLimit*sigmaY_fitted);
		double sigmaY_fitted_max_default= sigmaY_fitted + fabs(fitOptions.sigmaLimit*sigmaY_fitted);
		double sigmaY_fitted_min= std::max(0.,sigmaY_fitted_min_default);
		double sigmaY_fitted_max= sigmaY_fitted_max_default;
		if(fitOptions.fixSigma){
			if(!minuit->IsFixed(parNumber)) minuit->FixParameter(parNumber);
		}
		else{
			if(minuit->IsFixed(parNumber)) minuit->ReleaseParameter(parNumber);
			if(fitOptions.limitSigmaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sigmaY_fitted<<") to ["<<sigmaY_fitted_min<<","<<sigmaY_fitted_max<<"]");		
				minuit->SetParameter(parNumber,parName,sigmaY_fitted,sigmaY_fitted_err,sigmaY_fitted_min,sigmaY_fitted_max);
			}
		}

		//Theta
		parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double theta_fitted= minuit->GetParameter(parNumber);
		double theta_fitted_err= fabs(theta_fitted*theta_step);
		double theta_limits= fitOptions.thetaLimit*TMath::DegToRad();
		//double theta_limits= fitOptions.thetaLimit;
		double theta_fitted_min= theta_fitted - theta_limits;
		double theta_fitted_max= theta_fitted + theta_limits;
		if(fitOptions.fixTheta){
			if(!minuit->IsFixed(parNumber)) minuit->FixParameter(parNumber);
		} 
		else{
			if(minuit->IsFixed(parNumber)) minuit->ReleaseParameter(parNumber);
			if(fitOptions.limitThetaInFit){
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<theta_fitted<<") to ["<<theta_fitted_min<<","<<theta_fitted_max<<"]");
				minuit->SetParameter(parNumber,parName,theta_fitted,theta_fitted_err,theta_fitted_min,theta_fitted_max);
			}
		}
	}//end loop components

	
	//Bkg offset
	TString offset_name= minuit->GetParName(nFitPars-1);
	double offset_fitted= minuit->GetParameter(nFitPars-1);
	double offset_fitted_err= fabs(offset_fitted*offset_step);
	double offset_fitted_min= offset_fitted-fabs(m_rmsMean);
	double offset_fitted_max= offset_fitted+fabs(m_rmsMean);
	if(fitOptions.fixBkg){
		if(!minuit->IsFixed(nFitPars-1)) minuit->FixParameter(nFitPars-1);
	}
	else{
		if(minuit->IsFixed(nFitPars-1)) minuit->ReleaseParameter(nFitPars-1);
		if(fitOptions.limitBkgInFit){
			INFO_LOG("Limiting offset par (parNumber="<<nFitPars-1<<", value="<<offset_fitted<<") to ["<<offset_fitted_min<<","<<offset_fitted_max<<"]");
			minuit->SetParameter(nFitPars-1,offset_name,offset_fitted,offset_fitted_err,offset_fitted_min,offset_fitted_max);
		}
	}
	

	//==============================================
	//==             FULL FIT 
	//==============================================	
	// Set print level
	arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT",arglist,2);
	
	// Minimize
  arglist[0] = fitOptions.fitMaxIters; // number of function calls
  arglist[1] = fitOptions.fitFcnTolerance; // tolerance

	//NB: Minuit fit status = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult (4: MIGRAD error, 40: MINOS error, 400: HESSE error, 4000: IMPROVE error)
	//fitMinimizerStatus= fitRes.Status();

	//If any pars is at limit try to release it and fit again
	bool stopFit= false;
	int fitStatus= 0;	
	long int niters_fit= 0;
	long int niters_fit_max= fitOptions.fitNRetries;
	std::string defaultMinimizer= "MINIMIZE";//"MIGRAD"
	std::string minimizer= defaultMinimizer;
	std::string fallbackMinimizer= "SIMPLEX";
	std::string finalMinimizer= "HESS";
	if(fitOptions.fitFinalMinimizer==eHESS) finalMinimizer= "HESS";
	else if(fitOptions.fitFinalMinimizer==eMINOS) finalMinimizer= "MINOS";
	else {
		WARN_LOG("Invalid or redundant final minimizer step selected ("<<fitOptions.fitFinalMinimizer<<"), selected "<<finalMinimizer<<" by default (hint: you should choose HESS or MINOS)");
	}
	

	INFO_LOG("Start fit procedure with minimizer "<<minimizer<<" (final minimizer="<<finalMinimizer<<", nretries="<<niters_fit_max<<")");
	while(!stopFit){
		//Perform minimization
		INFO_LOG("Fitting source "<<aSource->GetName()<<" (#"<<niters_fit<<" iter cycle, minimizer="<<minimizer<<")...");
		fitStatus= minuit->ExecuteCommand(minimizer.c_str(),arglist,2);
		bool fitConverged= (fitStatus==0);
			
		//Check if any par is at limits. If so, enlarge limits and re-fit	
		std::vector<int> parsAtLimits;
		GetParsAtLimits(parsAtLimits,minuit);
		bool hasParsAtLimits= !parsAtLimits.empty();
		bool modifyParBounds= false;
		bool removeBounds= false;
		if(hasParsAtLimits){
			modifyParBounds= true;
			if(fitConverged){
				INFO_LOG("Fit converged with parameters at limits, will enlarge par bounds and retry...");
				m_fitStatus= eFitConvergedWithWarns;
				minimizer= defaultMinimizer;
			}
			else{
				INFO_LOG("Fit failed/not-converged and parameters at limits, will enlarge par bounds, switch to a simpler minimizer and retry...");
				m_fitStatus= eFitNotConverged;	
				minimizer= fallbackMinimizer;			
			}
		}//close if
		else{
			if(fitConverged){
				m_fitStatus= eFitConverged;
				if(minimizer==fallbackMinimizer){//If converged with SIMPLEX switch back to MIGRAD and fit again
					INFO_LOG("Fit converged without parameters at limits but with a simpler minimizer, switch back to default minimizer and fit again...");
					minimizer= defaultMinimizer;
				}
				else if(minimizer==defaultMinimizer && fitOptions.fitDoFinalMinimizerStep){//Run another fit iteration with HESS/MINOS for better error estimate
					INFO_LOG("Fit converged without parameters at limits, switching to final minimizer for better error estimation and fit again...");
					minimizer= finalMinimizer;
				}
				else{
					INFO_LOG("Fit converged without parameters at limits, will stop fit iterations and get results...");
					stopFit= true;
				}
			}
			else{
				INFO_LOG("Fit failed with no parameters at limits, will enlarge par bounds, switch to a simpler minimizer and retry...");
				modifyParBounds= true;	
				minimizer= fallbackMinimizer;
			}
		}//close else

		//Modify par bounds if requested
		if(modifyParBounds){
			double limitStep= (niters_fit+1)*0.1;
			for(size_t i=0;i<parsAtLimits.size();i++){
				int parIndex= parsAtLimits[i];
				TString parName;
				double parVal, parErr;
				double parVal_min, parVal_max;
				minuit->GetParameter(parIndex,(char*)(parName.Data()),parVal,parErr,parVal_min,parVal_max);
				
				
				double parBoundInterval= parVal_max-parVal_min;
				double parVal_min_new= parVal_min - limitStep*fabs(parBoundInterval/2.);
				double parVal_max_new= parVal_max + limitStep*fabs(parBoundInterval/2.);
				parVal_min= parVal_min_new;	
				parVal_max= parVal_max_new;
				std::string parBaseName= CodeUtils::ExtractSubString(std::string(parName),"_");
				if(parBaseName=="x0"){
					parVal_min= std::max(parVal_min_new,Xmin_norm);
					parVal_max= std::min(parVal_max_new,Xmax_norm);
				} 
				if(parBaseName=="y0"){
					parVal_min= std::max(parVal_min_new,Ymin_norm);
					parVal_max= std::min(parVal_max_new,Ymax_norm);
				} 
				if(parBaseName=="sigmaX" || parBaseName=="sigmaY"){
					parVal_min= std::max(parVal_min_new,0.);
				}

				if(niters_fit<niters_fit_max) minuit->SetParameter(parIndex,(char*)(parName.Data()),parVal,parErr,parVal_min,parVal_max);
				else minuit->SetParameter(parIndex,(char*)(parName.Data()),parVal,parErr,0,0);//if progressive bound enlarge does not work try remove the bound as the last attempt
				
				//minuit->ReleaseParameter(parIndex);
			}//end loop pars
		}//close if

		if(niters_fit>niters_fit_max || !fitOptions.fitImproveConvergence) {	
			WARN_LOG("Maximum number of fitting cycles reached ("<<niters_fit_max<<"), will stop fitting!");
			break;
		}
		niters_fit++;
			
	}//end loop fit

	INFO_LOG("Source fit ended with status code "<<m_fitStatus<<" ...");

	//==============================================
	//==             FIT RESULTS
	//==============================================
	
	//Retrieve covariance matrix and print its eigenvalues
	double* errMatrixValues= minuit->GetCovarianceMatrix();
	//minuit->ExecuteCommand("SHOw COVar",arglist,2);
	//minuit->ExecuteCommand("SHOw EIGenvalues",arglist,2);

	//Retrieve fitted component parameters
	m_sourceFitPars.SetNComponents(m_NFitComponents);
	if(fitOptions.fixBkg) m_sourceFitPars.SetOffsetFixed(true);
	if(fitOptions.fixSigma) m_sourceFitPars.SetSigmaFixed(true);
	if(fitOptions.fixTheta) m_sourceFitPars.SetThetaFixed(true);
	
	par_counter= 0;
	double parAtLimitThr= 1.e-6;//in percentage
	bool hasParsAtLimits= false;

	for(int i=0;i<m_NFitComponents;i++){
		for(int j=0;j<nComponentPars;j++){
			//Get par info
			TString parName_global;
			double parVal, parErr;
			double parVal_min, parVal_max;
			minuit->GetParameter(par_counter,(char*)(parName_global.Data()),parVal,parErr,parVal_min,parVal_max);
			
			//Check if at limits
			bool parFixed= minuit->IsFixed(par_counter);
			bool parBounded= (parVal_min!=parVal_max && parVal_min!=0 && parVal_max!=0);
			if(parBounded && !parFixed){
				double parRelDiffMin= fabs(parVal/parVal_min-1);
				double parRelDiffMax= fabs(parVal/parVal_max-1);
				bool parAtLimits= (parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr);
				if(parAtLimits) {
					hasParsAtLimits= true;
					INFO_LOG("Fit parameter "<<std::string(parName_global)<<" is at limit (value="<<parVal<<", min/max="<<parVal_min<<"/"<<parVal_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
				}
			}//close if

			//Set par values in source par struct
			std::string parName= CodeUtils::ExtractSubString(std::string(parName_global),"_");

			//Check if par needs unit conversion
			if(parName=="x0" ){
				parVal+= m_sourceX0;
			}
			if(parName=="y0"){
				parVal+= m_sourceY0;
			}
			if(parName=="A"){//convert to Jy
				parVal/= 1.e+3;
				parErr/= 1.e+3;
			}
			if(parName=="theta"){//convert to deg
				parVal*= TMath::RadToDeg();
				parErr*= TMath::RadToDeg();
			}

			//Set fitted parameter			
			if(m_sourceFitPars.SetParValueAndError(i,parName,parVal,parErr)<0){
				WARN_LOG("Failed to set par "<<parName<<" (parName_global="<<std::string(parName_global)<<") value and error (check par name)!");
			}
			par_counter++;
		}//end loop pars per component
	}//end loop components
	
	//Retrieve fitted offset
	TString fittedOffsetName;
	double fittedOffset, fittedOffsetErr;
	double fittedOffset_min, fittedOffset_max;
	minuit->GetParameter(par_counter,(char*)fittedOffsetName.Data(),fittedOffset,fittedOffsetErr,fittedOffset_min, fittedOffset_max);
	
	bool offsetFixed= minuit->IsFixed(par_counter);
	bool offsetBounded= (fittedOffset_min!=fittedOffset_max && fittedOffset_min!=0 && fittedOffset_max!=0);
	if(offsetBounded && !offsetFixed){
		double parRelDiffMin= fabs(fittedOffset/fittedOffset_min-1);
		double parRelDiffMax= fabs(fittedOffset/fittedOffset_max-1);
		bool parAtLimits= (parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr);
		if(parAtLimits) {
			hasParsAtLimits= true;
			INFO_LOG("Offset parameter is at limit (value="<<fittedOffset<<", min/max="<<fittedOffset_min<<"/"<<fittedOffset_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
		}
	}
	//m_sourceFitPars.SetOffsetPar(fittedOffset);
	//m_sourceFitPars.SetOffsetParErr(fittedOffsetErr);
	m_sourceFitPars.SetOffsetPar(fittedOffset/1.e+3);//convert back to Jy
	m_sourceFitPars.SetOffsetParErr(fittedOffsetErr/1.e+3);//convert back to Jy

	//Set fit status if any pars at limits
	if(m_fitStatus==eFitConverged && hasParsAtLimits) m_fitStatus= eFitConvergedWithWarns;	


	//Retrieve chi2 & NDF
	double Chi2, edm, errdef;
  int nvpar, nparx;
  minuit->GetStats(Chi2,edm,errdef,nvpar,nparx);
	int NFreePars= nvpar;
	int NPars= nparx;
	int NFittedBins= (int)(m_fitData.size());
	double NDF= NFittedBins-NFreePars;

	INFO_LOG("Source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") fit info: status="<<fitStatus<<" fitterStatus="<<m_fitStatus<<" (hasParsAtLimits? "<<hasParsAtLimits<<"), NDF="<<NDF<<", Chi2="<<Chi2<<", NPars="<<NPars<<", NFreePars="<<NFreePars<<" NFittedBins="<<NFittedBins);

	m_sourceFitPars.SetChi2(Chi2);
	m_sourceFitPars.SetNDF(NDF);
	m_sourceFitPars.SetNPars(NPars);
	m_sourceFitPars.SetNFreePars(NFreePars);
	m_sourceFitPars.SetNFitPoints(NFittedBins);
	m_sourceFitPars.SetStatus(m_fitStatus);
	
	//==============================================
	//==       FIT COVARIANCE & DERIVATIVE MATRIX
	//==============================================
	//NB: Check for NAN values in covariance matrix. Replace with 0 (TO BE INVESTIGATED)
	//for(int i=0;i<NFreePars*NFreePars;i++){
	for(int i=0;i<NPars*NPars;i++){
		if(std::isnan(errMatrixValues[i])){
			errMatrixValues[i]= 0.;
		}
	}
	for(int i=0;i<NFreePars;i++){
		INFO_LOG("Par no. "<<i<<" variance="<<minuit->GetCovarianceMatrixElement(i,i));
	}
	

	//Set covariance matrix
	//if(m_sourceFitPars.SetCovarianceMatrix(errMatrixValues,NFreePars)<0){
	if(m_sourceFitPars.SetCovarianceMatrix(errMatrixValues,NPars)<0){
		WARN_LOG("Failed to set covariance matrix (cannot retrieve it from fitter or array size issue), not able to compute flux density errors later!");
	}
	//else{
	//	cout<<"*** COVARIANCE MATRIX ***"<<endl;
	//	m_sourceFitPars.PrintCovarianceMatrix();
	//}

	//Compute flux density derivative matrix
	if(m_sourceFitPars.ComputeFluxDensityDerivMatrix()<0){
		WARN_LOG("Failed to compute flux density derivative matrix!");
	}
	//else{
		//cout<<"*** FLUX DENSITY DERIV MATRIX ***"<<endl;
		//m_sourceFitPars.PrintFluxDensityDerivMatrix();
	//}

	//==============================================
	//==       TOTAL FLUX DENSITY & ERROR
	//==============================================
	//Compute total flux density and error
	m_sourceFitPars.ComputeFluxDensity();
	if(m_sourceFitPars.ComputeFluxDensityError()<0){
		WARN_LOG("Failed to compute total flux density error!");
	}
	
	//==============================================
	//==             FIT RESIDUALS
	//==============================================
	//Get fitted parameters in array
	int nParsTot= minuit->GetNumberTotalParameters();
	double fittedPars[nParsTot];
	for(int i=0;i<nParsTot;i++){
		double parVal= minuit->GetParameter(i);
		fittedPars[i]= parVal;
	}

	//Compute fit residuals	
	std::vector<double> residuals;
	double residualMax= -1.e+99;
	double residualMin= 1.e+99;
	double X[2];
	for(size_t i=0;i<m_fitData.size();i++){
		double x= m_fitData[i].x;
		double y= m_fitData[i].y;
		double data= m_fitData[i].S;
		if(!std::isnormal(data) || data==0) continue;
		X[0]= x;
		X[1]= y;
		double model= Gaus2DMixtureFcn(X,fittedPars);
  	double res = data - model;
		if(res<residualMin) residualMin= res;
		if(res>residualMax) residualMax= res;
		residuals.push_back(res);

	}//end loop fitted data

	
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
	
	//==============================================
	//==             PRINT FIT RESULTS
	//==============================================
	m_sourceFitPars.Print();


	//Delete fitter
	//NB: Suggested not to do
	if(minuit){
		delete minuit;
		minuit= 0;
	}

	return 0;

}//close DoChi2Fit()
*/


ROOT::Math::Minimizer* SourceFitter::InitMinimizer(int nFitPars,SourceFitOptions& fitOptions)
{
	//Initialize fitter according to options
	ROOT::Math::Minimizer* fitter = 0;
	//fitter = ROOT::Math::Factory::CreateMinimizer(fitOptions.fitMinimizer,fitOptions.fitMinimizerAlgo);

	if(fitOptions.fitMinimizer=="Minuit2" || fitOptions.fitMinimizer=="minuit2"){
		#ifdef MINUIT2_ENABLED
			fitter= new ROOT::Minuit2::Minuit2Minimizer((fitOptions.fitMinimizerAlgo).c_str());
		#else	
			WARN_LOG("Minuit2 was selected but not available in the system, switching to MINUIT minimizer.");
			fitter= new TMinuitMinimizer((fitOptions.fitMinimizerAlgo).c_str(),nFitPars);
		#endif
	}
	else if(fitOptions.fitMinimizer=="Minuit" || fitOptions.fitMinimizer=="minuit"){
		fitter= new TMinuitMinimizer((fitOptions.fitMinimizerAlgo).c_str(),nFitPars);
	}
	else if(fitOptions.fitMinimizer=="R" || fitOptions.fitMinimizer=="r"){
		#ifdef ROOTR_ENABLED
			fitter= new ROOT::Math::RMinimizer((fitOptions.fitMinimizerAlgo).c_str());
		#else 
			WARN_LOG("RMinimizer was selected but not available in the system, switching to MINUIT minimizer.");
			fitter= new TMinuitMinimizer((fitOptions.fitMinimizerAlgo).c_str(),nFitPars);
		#endif
	}
	else if(fitOptions.fitMinimizer=="GSL" || fitOptions.fitMinimizer=="gsl"){
		fitter= new ROOT::Math::GSLNLSMinimizer();
	}
	else{
		ERROR_LOG("Invalid or unsupported minimizer ("<<fitOptions.fitMinimizer<<") given!");
		return nullptr;
	}

	if(!fitter){
		ERROR_LOG("Fitter creation failed (hint: chosen minimizer "<<fitOptions.fitMinimizer<<" not supported in the system?)");
		return nullptr;
	}

	//Set minimizer options
	//fitter->UseStaticMinuit(false);	
	fitter->SetMaxFunctionCalls(fitOptions.fitMaxIters);// for Minuit/Minuit2 
  fitter->SetMaxIterations(fitOptions.fitMaxIters);// for GSL 
  fitter->SetTolerance(fitOptions.fitFcnTolerance);//default tolerance
	fitter->SetPrecision(-1);//let minimizer choose the default
	fitter->SetErrorDef(1);//set chi2 error definition (0.5 for likelihood fits)
  fitter->SetPrintLevel(fitOptions.fitPrintLevel);//1=low
	fitter->SetStrategy(fitOptions.fitStrategy);//2=better error calculation
	fitter->SetValidError(fitOptions.fitDoFinalMinimizerStep);//run HESS for MINUIT if true

	return fitter;

}//close InitMinimizer()


int SourceFitter::DoChi2Fit(Source* aSource,SourceFitOptions& fitOptions,std::vector< std::vector<double> >& fitPars_start)
{
	//## Get source pars
	float Xmin, Xmax, Ymin, Ymax;
	aSource->GetSourceRange(Xmin,Xmax,Ymin,Ymax);
	double Smax= aSource->GetSmax();
	double Smin= aSource->GetSmin();

	double Xmin_norm= Xmin-m_sourceX0;
	double Xmax_norm= Xmax-m_sourceX0;
	double Ymin_norm= Ymin-m_sourceY0;
	double Ymax_norm= Ymax-m_sourceY0;

	//Convert to mJy
	Smax*= 1.e+3;
	Smin*= 1.e+3;

	//==============================================
	//==             INITIALIZE FITTER
	//==============================================
	//## Initialize fit parameters
	m_NFitComponents= static_cast<int>(fitPars_start.size());
	int nComponentPars= 6;
	int nFitPars= nComponentPars*m_NFitComponents + 1;//fit components + constant offset

	//## Initialize fitter
	ROOT::Math::Minimizer* fitter= InitMinimizer(nFitPars,fitOptions);
	if(!fitter){
		ERROR_LOG("Failed to initialize minimizer ("<<fitOptions.fitMinimizer<<")!");
		return -1;
	}

	/*
	//ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer(fitOptions.fitMinimizer,fitOptions.fitMinimizerAlgo);
	ROOT::Math::Minimizer* fitter = 0;
	if(fitOptions.fitMinimizer=="Minuit2" || fitOptions.fitMinimizer=="minuit2"){
		#ifdef MINUIT2_ENABLED
			fitter= new ROOT::Minuit2::Minuit2Minimizer((fitOptions.fitMinimizerAlgo).c_str());
		#else	
			WARN_LOG("Minuit2 was selected but not available in the system, switching to MINUIT minimizer.");
			fitter= new TMinuitMinimizer((fitOptions.fitMinimizerAlgo).c_str(),nFitPars);
		#endif
	}
	else if(fitOptions.fitMinimizer=="Minuit" || fitOptions.fitMinimizer=="minuit"){
		fitter= new TMinuitMinimizer((fitOptions.fitMinimizerAlgo).c_str(),nFitPars);
	}
	else if(fitOptions.fitMinimizer=="R" || fitOptions.fitMinimizer=="r"){
		#ifdef ROOTR_ENABLED
			fitter= new ROOT::Math::RMinimizer((fitOptions.fitMinimizerAlgo).c_str());
		#else 
			WARN_LOG("RMinimizer was selected but not available in the system, switching to MINUIT minimizer.");
			fitter= new TMinuitMinimizer((fitOptions.fitMinimizerAlgo).c_str(),nFitPars);
		#endif
	}
	else if(fitOptions.fitMinimizer=="GSL" || fitOptions.fitMinimizer=="gsl"){
		fitter= new ROOT::Math::GSLNLSMinimizer();
	}
	else{
		ERROR_LOG("Invalid or unsupported minimizer ("<<fitOptions.fitMinimizer<<") given!");
		return -1;
	}

	if(!fitter){
		ERROR_LOG("Fitter creation failed (hint: chosen minimizer "<<fitOptions.fitMinimizer<<" not supported in the system?)");
		return -1;
	}
	//fitter->UseStaticMinuit(false);	
	fitter->SetMaxFunctionCalls(fitOptions.fitMaxIters);// for Minuit/Minuit2 
  fitter->SetMaxIterations(fitOptions.fitMaxIters);// for GSL 
  fitter->SetTolerance(fitOptions.fitFcnTolerance);//default tolerance
	fitter->SetPrecision(-1);//let minimizer choose the default
	fitter->SetErrorDef(1);//set chi2 error definition (0.5 for likelihood fits)
  fitter->SetPrintLevel(fitOptions.fitPrintLevel);//1=low
	fitter->SetStrategy(fitOptions.fitStrategy);//2=better error calculation
	fitter->SetValidError(fitOptions.fitDoFinalMinimizerStep);//run HESS for MINUIT if true
	*/
	
	//## Set minimization function
	ROOT::Math::Functor fitFcn(std::bind(&SourceFitter::Chi2Fcn, this, _1), nFitPars);	
	fitter->SetFunction(fitFcn);
	

	//## Set start fit pars
	int par_counter= 0;
	std::string parNamePrefix[]= {"A","x0","y0","sigmaX","sigmaY","theta"};
	std::vector<std::vector<double>> parLimits_min;
	std::vector<std::vector<double>> parLimits_max;

	//- Offset
	double offset= 0.;
	double offset_min= 0;
	double offset_max= 0;
	double offset_step= 0.1;
	double offset_err= 0.;
	std::string offsetParName= "offset";
	if(fitOptions.useEstimatedBkgLevel) offset= m_bkgMean;//use estimated avg bkg
	else offset= fitOptions.fixedBkgLevel*1.e+3;//use user-supplied bkg level converted in mJy
	if(m_bkgMean>Smax) {
		WARN_LOG("Offset start par given/estimated is below min or above max source flux, setting it to Smin...");
		offset= Smin;
	}
	offset_err= fabs(offset_step*offset);
	if(fitOptions.limitBkgInFit){
		offset_min= std::min(Smin,offset-fabs(m_rmsMean));
		offset_max= std::min(Smax,offset+fabs(m_rmsMean));
		INFO_LOG("Limiting offset par "<<offset<<" in range ["<<offset_min<<","<<offset_max<<"]");
	}
	
	//- Component pars
	double Speak_step= 0.01;
	double centroid_step= 0.01;
	double sigma_step= 0.01;
	double theta_step= 0.01;
	std::map<std::string,int> parNameMap;

	for(int i=0;i<m_NFitComponents;i++){
		//Init min/max limits
		parLimits_min.push_back(std::vector<double>(nComponentPars,0));
		parLimits_max.push_back(std::vector<double>(nComponentPars,0));
		
		//## Set i-th component parameters
		//- Amplitude
		double Speak= fitPars_start[i][0] - offset;//subtract bkg offset from peak
		double Speak_min= 0;
		double Speak_max= 0;
		double Speak_err= fabs(Speak*Speak_step);
		std::string parName= Form("%s_%d",parNamePrefix[0].c_str(),i+1);
		if(fitOptions.limitAmplInFit){
			Speak_min= std::min(Smin, Speak*(1 - fitOptions.amplLimit) );
			Speak_max= std::max(Smax, Speak*(1 + fitOptions.amplLimit) );
			parLimits_min[i][0]= Speak_min;
			parLimits_max[i][0]= Speak_max;
			INFO_LOG("Limiting amplitude par "<<Speak<<" in range ["<<Speak_min<<","<<Speak_max<<"]");
		}
		parNameMap.insert( std::make_pair(std::string(parName),par_counter) );
		fitter->SetLimitedVariable(par_counter,parName,Speak,Speak_err,Speak_min,Speak_max);
		if(fitOptions.fixAmplInPreFit){
			fitter->FixVariable(par_counter);
		}

		//- Centroids
		double x0= fitPars_start[i][1];
		double y0= fitPars_start[i][2];
		double x0_err= 1;//1 pixel around centroid
		double y0_err= 1;
		double centroidLimit= fitOptions.centroidLimit;
		double x0_min= 0;
		double x0_max= 0;
		double y0_min= 0;
		double y0_max= 0;
		if(fitOptions.limitCentroidInFit){
			x0_min= std::max(x0 - centroidLimit,Xmin_norm);
			x0_max= std::min(x0 + centroidLimit,Xmax_norm);
			y0_min= std::max(y0 - centroidLimit,Ymin_norm);
			y0_max= std::min(y0 + centroidLimit,Ymax_norm);
			parLimits_min[i][1]= x0_min;
			parLimits_max[i][1]= x0_max;
			parLimits_min[i][2]= y0_min;
			parLimits_max[i][2]= y0_max;
			INFO_LOG("Limiting centroid pars (x,y)=("<<x0<<","<<y0<<")"<<" to bounds x0("<<x0_min<<","<<x0_max<<"), y0("<<y0_min<<","<<y0_max<<")");
		}
		parName= Form("%s_%d",parNamePrefix[1].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+1) );
		fitter->SetLimitedVariable(par_counter+1,parName,x0,x0_err,x0_min,x0_max);
		parName= Form("%s_%d",parNamePrefix[2].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+2) );
		fitter->SetLimitedVariable(par_counter+2,parName,y0,y0_err,y0_min,y0_max);
		if(fitOptions.fixCentroidInPreFit){
			fitter->FixVariable(par_counter+1);		
			fitter->FixVariable(par_counter+2);
		}			
	

		//- Sigmas
		double sigmaX= fitPars_start[i][3];
		double sigmaY= fitPars_start[i][4];
		double sourceSigmaMax_x= fabs(Xmax-Xmin)*sqrt(2.)/GausSigma2FWHM;
		double sourceSigmaMax_y= fabs(Ymax-Ymin)*sqrt(2.)/GausSigma2FWHM;
		double sigmaX_err= fabs(sigmaX*sigma_step);
		double sigmaY_err= fabs(sigmaY*sigma_step);
		double sigmaX_min= 0;
		double sigmaY_min= 0;
		double sigmaX_max= 0;
		double sigmaY_max= 0;
		if(fitOptions.limitSigmaInFit){
			sigmaX_min= std::max(0.,sigmaX*(1-fitOptions.sigmaLimit));
			sigmaY_min= std::max(0.,sigmaY*(1-fitOptions.sigmaLimit));
			sigmaX_max= std::max(sourceSigmaMax_x,sigmaX*(1+fitOptions.sigmaLimit));
			sigmaY_max= std::max(sourceSigmaMax_y,sigmaY*(1+fitOptions.sigmaLimit));
			parLimits_min[i][3]= sigmaX_min;
			parLimits_max[i][3]= sigmaX_max;
			parLimits_min[i][4]= sigmaY_min;
			parLimits_max[i][4]= sigmaY_max;
			INFO_LOG("Limiting sigma pars (sigmaX,sigmaY)=("<<sigmaX<<","<<sigmaY<<")"<<" to bounds sigmaX("<<sigmaX_min<<","<<sigmaX_max<<"), sigmaY("<<sigmaY_min<<","<<sigmaY_max<<"), bmaj="<<fitOptions.bmaj<<", bmin="<<fitOptions.bmin);
		}
		parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+3) );	
		fitter->SetLimitedVariable(par_counter+3,parName,sigmaX,sigmaX_err,sigmaX_min,sigmaX_max);
		parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+4) );
		fitter->SetLimitedVariable(par_counter+4,parName,sigmaY,sigmaY_err,sigmaY_min,sigmaY_max);
		if(fitOptions.fixSigmaInPreFit){//fix sigma in pre-fit?
			fitter->FixVariable(par_counter+3);		
			fitter->FixVariable(par_counter+4);		
		}

		//- Theta		
		double theta= fitPars_start[i][5];
		double theta_min= 0;
		double theta_max= 0;
		double theta_err= fabs(theta*theta_step);
		double theta_limits= fitOptions.thetaLimit*TMath::DegToRad();
		if(fitOptions.limitThetaInFit){
			theta_min= theta - theta_limits;
			theta_max= theta + theta_limits;
			parLimits_min[i][5]= theta_min;
			parLimits_max[i][5]= theta_max;
			INFO_LOG("Limiting theta par "<<theta<<" to bounds ("<<theta_min<<","<<theta_max<<")");
		}
		parName= Form("%s_%d",parNamePrefix[5].c_str(),i+1);
		parNameMap.insert( std::make_pair(std::string(parName),par_counter+5) );
		fitter->SetLimitedVariable(par_counter+5,parName,theta,theta_err,theta_min,theta_max);
		if(fitOptions.fixThetaInPreFit){//Fix theta par in pre-fit?
			fitter->FixVariable(par_counter+5);	
		}
		
		//Update par counter
		par_counter+= nComponentPars;
	
	}//end loop components

	//- Offset
	fitter->SetLimitedVariable(nFitPars-1,offsetParName,offset,offset_err,offset_min,offset_max);
	fitter->FixVariable(nFitPars-1);//Fix offset par for pre-fit
	

	//==============================================
	//==             PRE-FIT
	//==============================================
	INFO_LOG("Performing pre-fit (ndim="<<fitter->NDim()<<")");

	//## Perform pre-fit
	bool prefitStatus= fitter->Minimize();
	if(!prefitStatus){
		WARN_LOG("Source pre-fit failed or did not converge, trying to perform the full fit...");
		m_fitStatus= eFitNotConverged;
	}	

	//==============================================
	//==             RELEASE PARS
	//==============================================
	
	//## Release pars and adjust par limits around pre-fitted estimates
	const double* parValues= fitter->X();

	for(int i=0;i<m_NFitComponents;i++){
		//Amplitudes
		std::string parName= Form("%s_%d",parNamePrefix[0].c_str(),i+1);
		int parNumber= parNameMap[std::string(parName)];
		double Ampl_fitted= parValues[parNumber];
		double Ampl_fitted_min= Ampl_fitted - fabs(fitOptions.amplLimit*Ampl_fitted);
		double Ampl_fitted_max= Ampl_fitted + fabs(fitOptions.amplLimit*Ampl_fitted);
		double Ampl_fitted_err= fabs(Ampl_fitted*Speak_step);
		if(fitter->IsFixedVariable(parNumber)) fitter->ReleaseVariable(parNumber);	
		if(fitOptions.limitAmplInFit){
			INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<Ampl_fitted<<") to ["<<Ampl_fitted_min<<","<<Ampl_fitted_max<<"]");
			fitter->SetVariableLimits(parNumber,Ampl_fitted_min,Ampl_fitted_max);
			fitter->SetVariableStepSize(parNumber,Ampl_fitted_err);
		}

		//Centroid X
		parName= Form("%s_%d",parNamePrefix[1].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double X0_fitted= parValues[parNumber];
		double X0_fitted_min_default= X0_fitted - fabs(fitOptions.centroidLimit);
		double X0_fitted_max_default= X0_fitted + fabs(fitOptions.centroidLimit);
		double X0_fitted_min= std::max(Xmin_norm,X0_fitted_min_default);
		double X0_fitted_max= std::min(Xmax_norm,X0_fitted_max_default);
		double X0_fitted_err= 1;
		if(fitter->IsFixedVariable(parNumber)) fitter->ReleaseVariable(parNumber);	
		if(fitOptions.limitCentroidInFit){
			INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<X0_fitted<<") to ["<<X0_fitted_min<<","<<X0_fitted_max<<"]");
			fitter->SetVariableLimits(parNumber,X0_fitted_min,X0_fitted_max);
			fitter->SetVariableStepSize(parNumber,X0_fitted_err);
		}

		//Centroid Y
		parName= Form("%s_%d",parNamePrefix[2].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double Y0_fitted= parValues[parNumber];
		double Y0_fitted_min_default= Y0_fitted - fabs(fitOptions.centroidLimit);
		double Y0_fitted_max_default= Y0_fitted + fabs(fitOptions.centroidLimit);
		double Y0_fitted_min= std::max(Ymin_norm,Y0_fitted_min_default);
		double Y0_fitted_max= std::min(Ymax_norm,Y0_fitted_max_default);
		double Y0_fitted_err= 1;
		if(fitter->IsFixedVariable(parNumber)) fitter->ReleaseVariable(parNumber);
		if(fitOptions.limitCentroidInFit){
			INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<Y0_fitted<<") to ["<<Y0_fitted_min<<","<<Y0_fitted_max<<"]");
			fitter->SetVariableLimits(parNumber,Y0_fitted_min,Y0_fitted_max);
			fitter->SetVariableStepSize(parNumber,Y0_fitted_err);
		}

		//Sigma X
		parName= Form("%s_%d",parNamePrefix[3].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double sigmaX_fitted= parValues[parNumber];
		double sigmaX_fitted_err= fabs(sigmaX_fitted*sigma_step);
		double sigmaX_fitted_min_default= sigmaX_fitted - fabs(fitOptions.sigmaLimit*sigmaX_fitted);
		double sigmaX_fitted_max_default= sigmaX_fitted + fabs(fitOptions.sigmaLimit*sigmaX_fitted);
		double sigmaX_fitted_min= std::max(0.,sigmaX_fitted_min_default);
		double sigmaX_fitted_max= sigmaX_fitted_max_default;

		if(fitOptions.fixSigma){
			if(!fitter->IsFixedVariable(parNumber)) fitter->FixVariable(parNumber);	
		}
		else{
			if(fitter->IsFixedVariable(parNumber)) fitter->ReleaseVariable(parNumber);
			if(fitOptions.limitSigmaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sigmaX_fitted<<") to ["<<sigmaX_fitted_min<<","<<sigmaX_fitted_max<<"]");		
				fitter->SetVariableLimits(parNumber,sigmaX_fitted_min,sigmaX_fitted_max);
				fitter->SetVariableStepSize(parNumber,sigmaX_fitted_err);
			}
		}

		//SigmaY
		parName= Form("%s_%d",parNamePrefix[4].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double sigmaY_fitted= parValues[parNumber];
		double sigmaY_fitted_err= fabs(sigmaY_fitted*sigma_step);
		double sigmaY_fitted_min_default= sigmaY_fitted - fabs(fitOptions.sigmaLimit*sigmaY_fitted);
		double sigmaY_fitted_max_default= sigmaY_fitted + fabs(fitOptions.sigmaLimit*sigmaY_fitted);
		double sigmaY_fitted_min= std::max(0.,sigmaY_fitted_min_default);
		double sigmaY_fitted_max= sigmaY_fitted_max_default;
		if(fitOptions.fixSigma){
			if(!fitter->IsFixedVariable(parNumber)) fitter->FixVariable(parNumber);
		}
		else{
			if(fitter->IsFixedVariable(parNumber)) fitter->ReleaseVariable(parNumber);
			if(fitOptions.limitSigmaInFit) {
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<sigmaY_fitted<<") to ["<<sigmaY_fitted_min<<","<<sigmaY_fitted_max<<"]");		
				fitter->SetVariableLimits(parNumber,sigmaY_fitted_min,sigmaY_fitted_max);
				fitter->SetVariableStepSize(parNumber,sigmaY_fitted_err);
			}
		}

		//Theta
		parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
		parNumber= parNameMap[std::string(parName)];
		double theta_fitted= parValues[parNumber];
		double theta_fitted_err= fabs(theta_fitted*theta_step);
		double theta_limits= fitOptions.thetaLimit*TMath::DegToRad();
		double theta_fitted_min= theta_fitted - theta_limits;
		double theta_fitted_max= theta_fitted + theta_limits;
		if(fitOptions.fixTheta){
			if(!fitter->IsFixedVariable(parNumber)) fitter->FixVariable(parNumber);	
		} 
		else{
			if(fitter->IsFixedVariable(parNumber)) fitter->ReleaseVariable(parNumber);
			if(fitOptions.limitThetaInFit){
				INFO_LOG("Limiting par "<<parName<<" (parNumber="<<parNumber<<", value="<<theta_fitted<<") to ["<<theta_fitted_min<<","<<theta_fitted_max<<"]");
				fitter->SetVariableLimits(parNumber,theta_fitted_min,theta_fitted_max);
				fitter->SetVariableStepSize(parNumber,theta_fitted_err);
			}
		}
	}//end loop components

	
	//Bkg offset
	std::string offset_name= fitter->VariableName(nFitPars-1);
	double offset_fitted= parValues[nFitPars-1];
	double offset_fitted_err= fabs(offset_fitted*offset_step);
	double offset_fitted_min= offset_fitted-fabs(m_rmsMean);
	double offset_fitted_max= offset_fitted+fabs(m_rmsMean);
	if(fitOptions.fixBkg){
		if(!fitter->IsFixedVariable(nFitPars-1)) fitter->FixVariable(nFitPars-1);
	}
	else{
		if(fitter->IsFixedVariable(nFitPars-1)) fitter->ReleaseVariable(nFitPars-1);
		if(fitOptions.limitBkgInFit){
			INFO_LOG("Limiting offset par (parNumber="<<nFitPars-1<<", value="<<offset_fitted<<") to ["<<offset_fitted_min<<","<<offset_fitted_max<<"]");
			fitter->SetVariableLimits(nFitPars-1,offset_fitted_min,offset_fitted_max);
			fitter->SetVariableStepSize(nFitPars-1,offset_fitted_err);
		}
	}
	

	//==============================================
	//==             FULL FIT 
	//==============================================	
	//If any pars is at limit try to release it and fit again
	bool stopFit= false;
	int fitStatus= 0;	
	long int niters_fit= 0;
	long int niters_fit_max= fitOptions.fitNRetries;
	

	INFO_LOG("Start fit procedure (nretries="<<niters_fit_max<<")");
	while(!stopFit){
		//Perform minimization
		INFO_LOG("Fitting source "<<aSource->GetName()<<" (#"<<niters_fit<<" iter cycles performed) ...");
		bool fitConverged= fitter->Minimize();			

		//Check if any par is at limits. If so, enlarge limits and re-fit	
		std::vector<int> parsAtLimits;
		GetParsAtLimits(parsAtLimits,fitter);
		bool hasParsAtLimits= !parsAtLimits.empty();
		bool modifyParBounds= false;
		bool removeBounds= false;
		if(hasParsAtLimits){
			modifyParBounds= true;
			if(fitConverged){
				INFO_LOG("Fit converged with parameters at limits, will enlarge par bounds and retry...");
				m_fitStatus= eFitConvergedWithWarns;
			}
			else{
				INFO_LOG("Fit failed/not-converged and parameters at limits, will enlarge par bounds, switch to a simpler minimizer and retry...");
				m_fitStatus= eFitNotConverged;	
			}
		}//close if
		else{
			if(fitConverged){
				INFO_LOG("Fit converged without parameters at limits, will stop fit iterations and get results...");
				stopFit= true;
				m_fitStatus= eFitConverged;
			}
			else{
				INFO_LOG("Fit failed without parameters at limits, will set status to not converged and stop iterations ...");
				stopFit= true;
				m_fitStatus= eFitNotConverged;
			}
		}//close else

		//Modify par bounds if requested
		if(modifyParBounds){
			double limitStep= (niters_fit+1)*fitOptions.fitParBoundIncreaseStepSize;
			for(size_t i=0;i<parsAtLimits.size();i++){
				int parIndex= parsAtLimits[i];
					
				ROOT::Fit::ParameterSettings parInfo;
				fitter->GetVariableSettings(parIndex,parInfo);
				std::string parName= parInfo.Name(); 
				double parVal= parInfo.Value();
				double parErr= parInfo.StepSize();
				double parVal_min= parInfo.LowerLimit();
				double parVal_max= parInfo.UpperLimit();
			
				double parBoundInterval= parVal_max-parVal_min;
				double parVal_min_new= parVal_min - limitStep*fabs(parBoundInterval/2.);
				double parVal_max_new= parVal_max + limitStep*fabs(parBoundInterval/2.);
				parVal_min= parVal_min_new;	
				parVal_max= parVal_max_new;
				std::string parBaseName= CodeUtils::ExtractSubString(std::string(parName),"_");
				if(parBaseName=="x0"){
					parVal_min= std::max(parVal_min_new,Xmin_norm);
					parVal_max= std::min(parVal_max_new,Xmax_norm);
				} 
				if(parBaseName=="y0"){
					parVal_min= std::max(parVal_min_new,Ymin_norm);
					parVal_max= std::min(parVal_max_new,Ymax_norm);
				} 
				if(parBaseName=="sigmaX" || parBaseName=="sigmaY"){
					parVal_min= std::max(parVal_min_new,0.);
				}

				fitter->SetVariableStepSize(parIndex,parErr);
				if(niters_fit<niters_fit_max) {
					fitter->SetVariableLimits(parIndex,parVal_min,parVal_max);
				}
				else {
					//fitter->SetVariableLimits(parIndex,0,0);//valid only in MINUIT minimizer?
					fitter->SetVariableLimits(parIndex,2,1);//dummy value low>up to remove limits
				}
				
			}//end loop pars
		}//close if

		if(niters_fit>niters_fit_max || !fitOptions.fitImproveConvergence) {	
			WARN_LOG("Maximum number of fitting cycles reached ("<<niters_fit_max<<"), will stop fitting!");
			break;
		}
		niters_fit++;
			
	}//end loop fit

	INFO_LOG("Source fit ended with status code "<<m_fitStatus<<" ...");

	//==============================================
	//==             FIT RESULTS
	//==============================================
	//Retrieve chi2 & NDF
	double Chi2= fitter->MinValue();
	int NFreePars= fitter->NFree();
	int NPars= fitter->NDim();
	int NFittedBins= (int)(m_fitData.size());
	double NDF= NFittedBins-NFreePars;

	//Retrieve fitted pars & errors
	const double* fittedPars= fitter->X();
	const double* fittedParErrors= fitter->Errors();
	if(!fittedParErrors){
		WARN_LOG("Null ptr to fitted par errors (hint: minimizer provides errors? "<<fitter->ProvidesError()<<")");
	}

	//Retrieve covariance matrix and print its eigenvalues
	double errMatrixValues[NPars*NPars];	
	fitter->GetCovMatrix(errMatrixValues);
	
	//Retrieve fitted component parameters
	m_sourceFitPars.SetNComponents(m_NFitComponents);
	if(fitOptions.fixBkg) m_sourceFitPars.SetOffsetFixed(true);
	if(fitOptions.fixSigma) m_sourceFitPars.SetSigmaFixed(true);
	if(fitOptions.fixTheta) m_sourceFitPars.SetThetaFixed(true);
	
	par_counter= 0;
	double parAtLimitThr= 1.e-6;//in percentage
	bool hasParsAtLimits= false;

	for(int i=0;i<m_NFitComponents;i++){
		for(int j=0;j<nComponentPars;j++){
			//Get par info
			ROOT::Fit::ParameterSettings parInfo;
			fitter->GetVariableSettings(par_counter,parInfo);
			std::string parName_global= parInfo.Name(); 
			double parVal_min= parInfo.LowerLimit();
			double parVal_max= parInfo.UpperLimit();
			bool parFixed= parInfo.IsFixed();
			bool parBounded= parInfo.IsDoubleBound();

			double parVal= fittedPars[par_counter];
			double parErr= 0;
			if(fittedParErrors) parErr= fittedParErrors[par_counter];

			//Check if at limits
			if(parBounded && !parFixed){
				double parRelDiffMin= fabs(parVal/parVal_min-1);
				double parRelDiffMax= fabs(parVal/parVal_max-1);
				bool parAtLimits= (parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr);
				if(parAtLimits) {
					hasParsAtLimits= true;
					INFO_LOG("Fit parameter "<<parName_global<<" is at limit (value="<<parVal<<", min/max="<<parVal_min<<"/"<<parVal_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
				}
			}//close if

			//Set par values in source par struct
			std::string parName= CodeUtils::ExtractSubString(std::string(parName_global),"_");

			//Check if par needs unit conversion
			if(parName=="x0" ){
				parVal+= m_sourceX0;
			}
			if(parName=="y0"){
				parVal+= m_sourceY0;
			}
			if(parName=="A"){//convert to Jy
				parVal/= 1.e+3;
				parErr/= 1.e+3;
			}
			if(parName=="theta"){//convert to deg
				parVal*= TMath::RadToDeg();
				parErr*= TMath::RadToDeg();
			}

			//Set fitted parameter			
			if(m_sourceFitPars.SetParValueAndError(i,parName,parVal,parErr)<0){
				WARN_LOG("Failed to set par "<<parName<<" (parName_global="<<std::string(parName_global)<<") value and error (check par name)!");
			}
			par_counter++;
		}//end loop pars per component
	}//end loop components
	
	
	//Retrieve fitted offset
	ROOT::Fit::ParameterSettings fittedOffsetInfo;
	fitter->GetVariableSettings(par_counter,fittedOffsetInfo);
	std::string fittedOffsetName= fittedOffsetInfo.Name(); 
	double fittedOffset_min= fittedOffsetInfo.LowerLimit();
	double fittedOffset_max= fittedOffsetInfo.UpperLimit();
	bool offsetFixed= fittedOffsetInfo.IsFixed();
	bool offsetBounded= fittedOffsetInfo.IsDoubleBound();
	double fittedOffset= fittedPars[par_counter];
	double fittedOffsetErr= 0;
	if(fittedParErrors) fittedOffsetErr= fittedParErrors[par_counter];

	if(offsetBounded && !offsetFixed){
		double parRelDiffMin= fabs(fittedOffset/fittedOffset_min-1);
		double parRelDiffMax= fabs(fittedOffset/fittedOffset_max-1);
		bool parAtLimits= (parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr);
		if(parAtLimits) {
			hasParsAtLimits= true;
			INFO_LOG("Offset parameter is at limit (value="<<fittedOffset<<", min/max="<<fittedOffset_min<<"/"<<fittedOffset_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
		}
	}
	
	m_sourceFitPars.SetOffsetPar(fittedOffset/1.e+3);//convert back to Jy
	m_sourceFitPars.SetOffsetParErr(fittedOffsetErr/1.e+3);//convert back to Jy

	//Set fit status if any pars at limits
	if(m_fitStatus==eFitConverged && hasParsAtLimits) m_fitStatus= eFitConvergedWithWarns;	
	fitStatus= fitter->Status();

	
	INFO_LOG("Source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") fit info: status="<<fitStatus<<", fitterStatus="<<m_fitStatus<<", covMatrixStatus="<<fitter->CovMatrixStatus()<<", IsValidError? "<<fitter->IsValidError()<<", ProvidesError? "<<fitter->ProvidesError()<<", fit strategy="<<fitter->Strategy()<<" (hasParsAtLimits? "<<hasParsAtLimits<<"), NDF="<<NDF<<", Chi2="<<Chi2<<", NPars="<<NPars<<", NFreePars="<<NFreePars<<" NFittedBins="<<NFittedBins);

	m_sourceFitPars.SetChi2(Chi2);
	m_sourceFitPars.SetNDF(NDF);
	m_sourceFitPars.SetNPars(NPars);
	m_sourceFitPars.SetNFreePars(NFreePars);
	m_sourceFitPars.SetNFitPoints(NFittedBins);
	m_sourceFitPars.SetStatus(m_fitStatus);
	
	
	//==============================================
	//==       FIT COVARIANCE & DERIVATIVE MATRIX
	//==============================================
	//NB: Check for NAN values in covariance matrix. Replace with 0 (TO BE INVESTIGATED)
	//for(int i=0;i<NFreePars*NFreePars;i++){
	for(int i=0;i<NPars*NPars;i++){
		if(std::isnan(errMatrixValues[i])){
			errMatrixValues[i]= 0.;
		}
	}
	
	//Set covariance matrix
	if(m_sourceFitPars.SetCovarianceMatrix(errMatrixValues,NPars)<0){
		WARN_LOG("Failed to set covariance matrix (cannot retrieve it from fitter or array size issue), not able to compute flux density errors later!");
	}
	
	//Compute flux density derivative matrix
	if(m_sourceFitPars.ComputeFluxDensityDerivMatrix()<0){
		WARN_LOG("Failed to compute flux density derivative matrix!");
	}

	//==============================================
	//==       TOTAL FLUX DENSITY & ERROR
	//==============================================
	//Compute total flux density and error
	m_sourceFitPars.ComputeFluxDensity();
	if(m_sourceFitPars.ComputeFluxDensityError()<0){
		WARN_LOG("Failed to compute total flux density error!");
	}
	
	//==============================================
	//==             FIT RESIDUALS
	//==============================================
	//Compute fit residuals	
	std::vector<double> residuals;
	double residualMax= -1.e+99;
	double residualMin= 1.e+99;
	double X[2];
	for(size_t i=0;i<m_fitData.size();i++){
		double x= m_fitData[i].x;
		double y= m_fitData[i].y;
		double data= m_fitData[i].S;
		if(!std::isnormal(data) || data==0) continue;
		X[0]= x;
		X[1]= y;
		double model= Gaus2DMixtureFcn(X,fittedPars);
  	double res = data - model;
		if(res<residualMin) residualMin= res;
		if(res>residualMax) residualMax= res;
		residuals.push_back(res);

	}//end loop fitted data

	
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
	
	//==============================================
	//==             PRINT FIT RESULTS
	//==============================================
	m_sourceFitPars.Print();
	
	//==============================================
	//==             CLEAR DATA
	//==============================================
	//Delete minimizer
	if(fitter){
		delete fitter;
		fitter= 0;
	}

	return 0;

}//close DoChi2Fit()

int SourceFitter::GetParsAtLimits(std::vector<int>& parsAtLimits,ROOT::Math::Minimizer* fitter)
{
	//Check minimizer ptr
	if(!fitter){
		ERROR_LOG("Null ptr to minimizer given!");
		return -1;
	}

	//Loop over pars and found which are at limits
	double parAtLimitThr= 1.e-6;//in percentage
	
	//int nParsTot= fitter->GetNumberTotalParameters();
	int nParsTot= fitter->NDim();

	for(int i=0;i<nParsTot;i++){
		//Get par info
		ROOT::Fit::ParameterSettings parInfo;
		fitter->GetVariableSettings(i,parInfo);

		std::string parName= parInfo.Name();
		double parVal= parInfo.Value();
		double parVal_min= parInfo.LowerLimit();
		double parVal_max= parInfo.UpperLimit();
		bool parFixed= parInfo.IsFixed();
		bool parBounded= parInfo.IsDoubleBound();
		/*
		TString parName;
		double parVal, parErr;
		double parVal_min, parVal_max;
		fitter->GetParameter(i,(char*)(parName.Data()),parVal,parErr,parVal_min,parVal_max);
			
		//Check if at limits
		bool parFixed= fitter->IsFixed(i);
		bool parBounded= (parVal_min!=parVal_max && parVal_min!=0 && parVal_max!=0);
		*/

		if(parBounded && !parFixed){
			double parRelDiffMin= fabs(parVal/parVal_min-1);
			double parRelDiffMax= fabs(parVal/parVal_max-1);
			bool isParAtLimits= (parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr);
			if(isParAtLimits) {
				INFO_LOG("Fit parameter "<<parName<<" is at limit (value="<<parVal<<", min/max="<<parVal_min<<"/"<<parVal_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
				parsAtLimits.push_back(i);
			}
		}//close if
	}//end loop pars

	return 0;

}//close GetParsAtLimits()

/*
int SourceFitter::GetParsAtLimits(std::vector<int>& parsAtLimits,TFitter* fitter)
{
	//Check minimizer ptr
	if(!fitter){
		ERROR_LOG("Null ptr to MINUIT minimizer given!");
		return -1;
	}

	//Loop over pars and found which are at limits
	double parAtLimitThr= 1.e-6;//in percentage
	
	int nParsTot= fitter->GetNumberTotalParameters();
	for(int i=0;i<nParsTot;i++){
		//Get par info
		TString parName;
		double parVal, parErr;
		double parVal_min, parVal_max;
		fitter->GetParameter(i,(char*)(parName.Data()),parVal,parErr,parVal_min,parVal_max);
			
		//Check if at limits
		bool parFixed= fitter->IsFixed(i);
		bool parBounded= (parVal_min!=parVal_max && parVal_min!=0 && parVal_max!=0);
		if(parBounded && !parFixed){
			double parRelDiffMin= fabs(parVal/parVal_min-1);
			double parRelDiffMax= fabs(parVal/parVal_max-1);
			bool isParAtLimits= (parRelDiffMin<parAtLimitThr || parRelDiffMax<parAtLimitThr);
			if(isParAtLimits) {
				INFO_LOG("Fit parameter "<<std::string(parName)<<" is at limit (value="<<parVal<<", min/max="<<parVal_min<<"/"<<parVal_max<<" rel diff min/max="<<parRelDiffMin<<"/"<<parRelDiffMax<<", thr="<<parAtLimitThr<<")");
				parsAtLimits.push_back(i);
			}
		}//close if
	}//end loop pars

	return 0;

}//close GetParsAtLimits()
*/

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
	m_chi2RegPar= fitOptions.chi2RegPar;

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
	//if(DoFit(aSource,fitOptions,fitPars_start)<0){
	if(DoChi2Fit(aSource,fitOptions,fitPars_start)<0){		
		ERROR_LOG("Failed to perform source fit!");
		m_fitStatus= eFitAborted;
		return -1;	
	}
	
	return 0;

}//close FitSource()



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

/*
void SourceFitter::Chi2Fcn(int& nPar,double* grad,double& fval,double* p,int iflag)
{
	//Loop over fit data points
	double X[2]= {0,0};
	double Serr= 1;
	if(m_rmsMean>0) Serr= m_rmsMean;
	double res= 0;
	double chi2= 0;
	double chi2_signal= 0;
	double chi2_halo= 0;

	//Compute chi2 (signal pixels)
	for(size_t i=0;i<m_fitData.size();i++){
		double x= m_fitData[i].x;
		double y= m_fitData[i].y;
		double S= m_fitData[i].S;
		//double Serr= m_fitData[i].Serr;
		X[0]= x;
		X[1]= y;
		res= (S-Gaus2DMixtureFcn(X,p))/Serr;
		chi2_signal+= res*res;

	}//end loop data

	//Compute chi2 (halo pixels)
	for(size_t i=0;i<m_fitHaloData.size();i++){
		double x= m_fitHaloData[i].x;
		double y= m_fitHaloData[i].y;
		double S= m_fitHaloData[i].S;
		//double Serr= m_fitHaloData[i].Serr;
		X[0]= x;
		X[1]= y;
		res= (S-Gaus2DMixtureFcn(X,p))/Serr;
		chi2_halo+= res*res;

	}//end loop data

	chi2= chi2_signal + m_chi2RegPar*chi2_halo;
  
	fval = chi2;
	
}//close Chi2Fcn()
*/

double SourceFitter::Chi2Fcn(const double* p)
{
	//Loop over fit data points
	double X[2]= {0,0};
	double Serr= 1;
	if(m_rmsMean>0) Serr= m_rmsMean;
	double res= 0;
	double chi2= 0;
	double chi2_signal= 0;
	double chi2_halo= 0;

	//Compute chi2 (signal pixels)
	for(size_t i=0;i<m_fitData.size();i++){
		double x= m_fitData[i].x;
		double y= m_fitData[i].y;
		double S= m_fitData[i].S;
		//double Serr= m_fitData[i].Serr;
		X[0]= x;
		X[1]= y;
		res= (S-Gaus2DMixtureFcn(X,p))/Serr;
		chi2_signal+= res*res;

	}//end loop data

	//Compute chi2 (halo pixels)
	for(size_t i=0;i<m_fitHaloData.size();i++){
		double x= m_fitHaloData[i].x;
		double y= m_fitHaloData[i].y;
		double S= m_fitHaloData[i].S;
		//double Serr= m_fitHaloData[i].Serr;
		X[0]= x;
		X[1]= y;
		res= (S-Gaus2DMixtureFcn(X,p))/Serr;
		chi2_halo+= res*res;

	}//end loop data

	chi2= chi2_signal + m_chi2RegPar*chi2_halo;
  
	double fval = chi2;
	
	return fval;

}//close Chi2Fcn()


double SourceFitter::Gaus2DMixtureFcn(double* x, const double* p)
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


double SourceFitter::Gaus2DFcn(double* x, const double* par){

	double X= x[0];
	double Y= x[1];
	
	double A= par[0];
	double X0= par[1];
	double Y0= par[2];
	//double sigmaX= par[3];
	//double sigmaY= par[4];
	double sigmaX= fabs(par[3]);
	double sigmaY= fabs(par[4]);
	double theta= par[5];
	//theta*= TMath::DegToRad();

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

