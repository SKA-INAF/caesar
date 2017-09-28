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
//std::vector<TEllipse> SourceFitter::m_fitEllipses;
int SourceFitter::m_fitStatus= eFitUnknownStatus;

//Constructor
SourceFitter::SourceFitter()
	: TObject()
{
	
}//close costructor

//Destructor
SourceFitter::~SourceFitter()
{

}//close destructor


int SourceFitter::FitSource(Source* aSource,BlobPars blobPars,int nMaxComponents)
{
	//## Check input source
	if(!aSource){
		ERROR_LOG("Null ptr to source given!");
		m_fitStatus= eFitAborted;
		return -1;
	}
	INFO_LOG("Fitting source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") assuming these blob pars: {Bmaj(pix)="<<blobPars.bmaj<<", Bmin(pix)="<<blobPars.bmin<<", Bpa(deg)="<<blobPars.bpa<<"}");

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
	INFO_LOG("Source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") range: X["<<Xmin<<","<<Xmax<<"], Y["<<Ymin<<","<<Ymax<<"]");

	//## Get source flux image
	INFO_LOG("Get flux histo & curv map for source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")");	
	long int nBoxX= Xmax - Xmin + 1;
	long int nBoxY= Ymax - Ymin + 1;
	TH2D* fluxMapHisto= new TH2D("","",nBoxX,Xmin-0.5,Xmax+0.5,nBoxY,Ymin-0.5,Ymax+0.5);
	fluxMapHisto->Sumw2();
	Image* curvMap= new Image(nBoxX,nBoxY,Xmin,Ymin,"");
	
	std::vector<Pixel*> pixels= aSource->GetPixels();
	double bkgMean= 0.;
	double rmsMean= 0.;
	for(size_t k=0;k<pixels.size();k++){
		double x= pixels[k]->x;
		double y= pixels[k]->y;
		double S= pixels[k]->S;
		double Scurv= pixels[k]->GetCurv();	
		std::pair<double,double> bkgData= pixels[k]->GetBkg();	
		double bkg= bkgData.first;
		double rms= bkgData.second;
		fluxMapHisto->Fill(x,y,S);
		curvMap->Fill(x,y,Scurv);	
		bkgMean+= bkg;
		rmsMean+= rms;
	}//end loop pixels
	
	bkgMean/= (double)(pixels.size());
	rmsMean/= (double)(pixels.size());
	INFO_LOG("source (id="<<aSource->Id<<", name="<<aSource->GetName()<<") bkg info: <bkg>="<<bkgMean<<", <rms>="<<rmsMean);

	//Check if histo has entries
	if(fluxMapHisto->GetSumOfWeights()<=0){
		ERROR_LOG("No entries present in flux map histogram (hint: check if properly filled, e.g. bin range problems)!");
		m_fitStatus= eFitAborted;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		delete curvMap;
		curvMap= 0;
		return -1;
	}

	//## Estimate the number of components present in source
	//## Method: Find peaks in curvature map
	INFO_LOG("Finding peaks in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<")");	
	std::vector<TVector2> peakPoints;
	int peakShiftTolerance= 2;
	bool skipBorders= true;
	if(curvMap->FindPeaks(peakPoints,peakShiftTolerance,skipBorders)<0){
		WARN_LOG("Failed to find peaks in source curvature map!");
		m_fitStatus= eFitAborted;
		delete curvMap;
		curvMap= 0;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		return -1;
	}

	int nComponents= static_cast<int>(peakPoints.size());
	if(nComponents<=0){
		WARN_LOG("No components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<"), this should not occur (at least one should be found)!");	
		m_fitStatus= eFitAborted;
		delete curvMap;
		curvMap= 0;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		return -1;
	}

	INFO_LOG("#"<<nComponents<<" components found in source (id="<<aSource->Id<<", name="<<aSource->GetName()<<"), #"<<nMaxComponents<<" components will be fitted at maximum...");
	
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
			delete curvMap;
			curvMap= 0;
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
	m_NFitComponents= std::min(nComponents,nMaxComponents);
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

	int par_counter= 0;
	for(int i=0;i<m_NFitComponents;i++){
		size_t index= sort_index[i];

		double x= peakPoints[index].X();
		double y= peakPoints[index].Y();
		long int gbin= fluxMapHisto->FindBin(x,y);
		double Speak= fluxMapHisto->GetBinContent(gbin);
		
		//## Set i-th component parameters
		//- Amplitude
		double Speak_min= std::max(Smin, Speak*0.90);
		double Speak_max= std::min(Smax, Speak*1.10);
		sourceFitFcn->SetParName(par_counter,Form("%s_%d",parNamePrefix[0].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter,Speak);
		sourceFitFcn->SetParLimits(par_counter,Speak_min,Speak_max);
		INFO_LOG("Speak="<<Speak<<" ["<<Speak_min<<","<<Speak_max<<"]");

		//- Centroids
		double centroidLimit= 0.5 * sqrt(pow(blobPars.bmaj,2) + pow(blobPars.bmin,2));
		double x0_min= x - centroidLimit;
		double x0_max= x + centroidLimit;
		double y0_min= y - centroidLimit;
		double y0_max= y + centroidLimit;
		sourceFitFcn->SetParName(par_counter+1,Form("%s_%d",parNamePrefix[1].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+1,x);
		sourceFitFcn->SetParLimits(par_counter+1,x0_min,x0_max);

		sourceFitFcn->SetParName(par_counter+2,Form("%s_%d",parNamePrefix[2].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+2,y);
		sourceFitFcn->SetParLimits(par_counter+2,y0_min,y0_max);
		INFO_LOG("(x,y)=("<<x<<","<<y<<")"<<" bounds x("<<x0_min<<","<<x0_max<<") y("<<y0_min<<","<<y0_max<<")");

		//- Sigmas
		double sigmaX= blobPars.bmaj/GausSigma2FWHM;
		double sigmaY= blobPars.bmin/GausSigma2FWHM;
		double sourceSigmaMax_x= fabs(Xmax-Xmin)*sqrt(2.)/GausSigma2FWHM;
		double sourceSigmaMax_y= fabs(Ymax-Ymin)*sqrt(2.)/GausSigma2FWHM;
		double sigmaX_min= 0.8*sigmaX;
		double sigmaY_min= 0.8*sigmaY;
		double sigmaX_max= std::max(sourceSigmaMax_x,sigmaX*1.1);
		double sigmaY_max= std::max(sourceSigmaMax_y,sigmaY*1.1);
		sourceFitFcn->SetParName(par_counter+3,Form("%s_%d",parNamePrefix[3].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+3,sigmaX);
		sourceFitFcn->SetParLimits(par_counter+3,sigmaX_min,sigmaX_max);
		sourceFitFcn->SetParName(par_counter+4,Form("%s_%d",parNamePrefix[4].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+4,sigmaY);
		sourceFitFcn->SetParLimits(par_counter+4,sigmaY_min,sigmaY_max);
		INFO_LOG("(sigmaX,sigmaY)=("<<sigmaX<<","<<sigmaY<<")"<<" bounds sigmaX("<<sigmaX_min<<","<<sigmaX_max<<") y("<<sigmaY_min<<","<<sigmaY_max<<")");

		//- Theta		
		double theta= blobPars.bpa;
		sourceFitFcn->SetParName(par_counter+5,Form("%s_%d",parNamePrefix[5].c_str(),i+1));
		sourceFitFcn->SetParameter(par_counter+5,theta);
		INFO_LOG("theta="<<theta);

		//Update par counter
		par_counter+= nComponentPars;
	
	}//end loop components

	//- Offset
	double offset= bkgMean;
	if(bkgMean<Smin || bkgMean>Smax) offset= Smin;
	double offset_min= std::min(Smin,offset-fabs(rmsMean));
	double offset_max= std::min(Smax,offset+fabs(rmsMean));
	sourceFitFcn->SetParName(nFitPars-1,"offset");
	sourceFitFcn->SetParameter(nFitPars-1,offset);
	sourceFitFcn->SetParLimits(nFitPars-1,offset_min,offset_max);

	//## Pre-fit with theta fixed to beam bpa
	for(int i=0;i<m_NFitComponents;i++){
		TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
		int parNumber= sourceFitFcn->GetParNumber(parName);
		double theta= blobPars.bpa;
		sourceFitFcn->FixParameter(parNumber,theta);
	}
	int prefitStatus= fluxMapHisto->Fit(sourceFitFcn,"R");
	if(prefitStatus!=0){
		WARN_LOG("Source pre-fit failed or did not converge.");
		m_fitStatus= eFitNotConverged;
		delete curvMap;
		curvMap= 0;
		delete fluxMapHisto;
		fluxMapHisto= 0;
		delete sourceFitFcn;
		sourceFitFcn= 0;
		return 0;
	}

	//## Release theta and fit again
	for(int i=0;i<m_NFitComponents;i++){
		TString parName= Form("%s_%d",parNamePrefix[nComponentPars-1].c_str(),i+1);
		int parNumber= sourceFitFcn->GetParNumber(parName);
		sourceFitFcn->ReleaseParameter(parNumber);
	}
	//TFitResultPtr fitRes= fluxMapHisto->Fit(sourceFitFcn,"S");
	fluxMapHisto->Fit(sourceFitFcn,"S");
	TBackCompFitter* fitter = (TBackCompFitter *) TVirtualFitter::GetFitter();
	const ROOT::Fit::FitResult& fitRes = fitter->GetFitResult();
	
	//## Retrieve fit results
	m_fitStatus= eFitNotConverged;
	int fitMinimizerStatus= -1;
	int fitStatus= -1;
	bool hasParsAtLimits= false;
	//if(fitRes) {
		if(fitRes.IsValid()) {
			fitStatus= 0;
			m_fitStatus= eFitConverged;
		}
	
		//NB: Minuit fit status = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult (4: MIGRAD error, 40: MINOS error, 400: HESSE error, 4000: IMPROVE error)
		fitMinimizerStatus= fitRes.Status();

		//Check if has parameters at bounds
		hasParsAtLimits= HasFitParsAtLimit(fitRes);
		if(hasParsAtLimits) m_fitStatus= eFitConvergedWithWarns;
	//}
	//else{
	//	WARN_LOG("Null ptr to fit result pointer!");
	//}

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
	m_sourceFitPars.SetStatus(fitRes.IsValid());
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
	
	
	//Clear up data
	delete curvMap;
	curvMap= 0;
	delete fluxMapHisto;
	fluxMapHisto= 0;
	delete sourceFitFcn;
	sourceFitFcn= 0;
	
	return 0;

}//close FitSource()


//bool SourceFitter::HasFitParsAtLimit(TFitResultPtr fitRes){
bool SourceFitter::HasFitParsAtLimit(const ROOT::Fit::FitResult& fitRes){

	//if(!fitRes){
	//	WARN_LOG("Null ptr to fit results given, returning false!");
	//	return false;
	//}

	bool hasParAtLimits= false;
	
	for(int i=0;i<fitRes.NPar();i++){
		if(!fitRes.IsParameterBound(i) || fitRes.IsParameterFixed(i) ) continue;
	
		//Check par limits
		double par_min= 0;
		double par_max= 0;
		fitRes.ParameterBounds(i,par_min,par_max);
		double par= fitRes.Parameter(i);		
		
		if(par==par_min || par==par_max){
			hasParAtLimits= true;
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

	return fcn;

}//close Gaus2DFcn()


}//close namespace

