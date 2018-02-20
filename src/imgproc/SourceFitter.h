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
* @file SourceFitter.h
* @class SourceFitter
* @brief SourceFitter
*
* Class to fit a source image with a mixture of gaussian/skew normal/skew-t bivariate functions
* @author S. Riggi
* @date 01/09/2017
*/

#ifndef _SOURCE_FITTER_h
#define _SOURCE_FITTER_h 1

#include <SysUtils.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TFitResultPtr.h>
#include <Fit/FitResult.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>

namespace Caesar {

class Image;
class Source;


struct SourceFitOptions {
	//- Blob start fit pars
	double bmaj;//in pixels
	double bmin;//in pixels
	double bpa;//in degrees

	//- Max number of components to be fitted
	int nMaxComponents;

	//- Number of matching peaks across kernels (-1=peak detected in all kernels, ...) and distance tolerance
	int peakKernelMultiplicityThr;
	int peakShiftTolerance;

	//- Dilation kernels to be used when finding peaks
	int peakMinKernelSize;
	int peakMaxKernelSize;

	//- Peak flux significance min threshold (in nsigmas wrt to avg bkg & rms)
	double peakZThrMin;	

	//- Centroid options
	bool limitCentroidInFit;

	//- Bkg options
	bool fixBkg;
	bool limitBkgInFit;
	bool useEstimatedBkgLevel;
	double fixedBkgLevel;
		
	//- Amplitude fit par range (example +-20% around source peak)
	bool limitAmplInFit;
	double amplLimit;

	//- Sigma fit par range (example +-20% around source peak)
	bool fixSigmaInPreFit;
	bool fixSigma;	
	bool limitSigmaInFit;
	double sigmaLimit;

	//- Theta 
	bool fixTheta;
	bool limitThetaInFit;
	double thetaLimit;//in deg

	//- Flux significance min threshold (in nsigmas above the bkg)
	bool useFluxZCut;
	double fluxZThrMin;

	

	//Default constructor
	SourceFitOptions() 
	{
    bmaj= 5;//pix
    bmin= 5;//pix
		bpa= 0;
		nMaxComponents= 3;
		limitCentroidInFit= true;
		fixBkg= true;
		limitBkgInFit= true;
		useEstimatedBkgLevel= true;
		fixedBkgLevel= 0;
		limitAmplInFit= true;
		amplLimit= 0.2;
		limitSigmaInFit= true;
		sigmaLimit= 0.2;
		fixSigmaInPreFit= false;
		fixSigma= false;
		fixTheta= false;
		limitThetaInFit= true;
		thetaLimit= 5;//deg
		useFluxZCut= false;
		fluxZThrMin= 2.5;//in nsigmas
		peakZThrMin= 0;//in nsigmas
		peakMinKernelSize= 3;
		peakMaxKernelSize= 7;
		peakKernelMultiplicityThr= 1;
		peakShiftTolerance= 2;
	}//close constructor
};//close SourceFitOptions


//========================================
//==         SOURCE COMPONENT PARS
//========================================
class SourceComponentPars : public TObject {

	public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceComponentPars(){
			
			std::vector<std::string> parNames {"A","x0","y0","sigmaX","sigmaY","theta"};
			for(size_t i=0;i<parNames.size();i++){
				FitPars.insert( std::pair<std::string,double>(parNames[i],0.) );
				FitParsErr.insert( std::pair<std::string,double>(parNames[i],0.) );
			}		
		}//close contructor

		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceComponentPars(){
			FitPars.clear();
			FitParsErr.clear();
		}

	public:
		/** 
		\brief Set par value & error
 		*/
		int SetParValueAndError(std::string parName,double parVal,double parErr){
			if(!CodeUtils::HasMapKey(FitPars,parName)) {
				WARN_LOG("Invalid par name ("<<parName<<" given, cannot find par to be set!");
				return -1;
			}
			FitPars[parName]= parVal;
			FitParsErr[parName]= parErr;
			return 0;
		}

		/** 
		\brief Get par value
 		*/
		double GetParValue(std::string parName){
			if(!CodeUtils::HasMapKey(FitPars,parName)) return -999;
			return FitPars[parName];
		}

		/** 
		\brief Get par error
 		*/
		double GetParError(std::string parName){
			if(!CodeUtils::HasMapKey(FitParsErr,parName)) return -999;
			return FitParsErr[parName];
		}

		/** 
		\brief Get fit pars
 		*/
		std::map<std::string,double>const& GetFitPars() const {return FitPars;}		
		/** 
		\brief Get fit par errors
 		*/
		std::map<std::string,double>const& GetFitParErrors() const {return FitParsErr;}		

		/** 
		\brief Get fit ellipse
 		*/
		TEllipse* GetFitEllipse(bool useFWHM=true)
		{	
			//Check if has fit pars
			if(FitPars.empty()) {
				WARN_LOG("No fitted pars stored, returning nullptr ellipse!");
				return nullptr;
			}
		
			//Compute fit ellipse from pars
			double x0= FitPars["x0"];
			double y0= FitPars["y0"];
			double sigmaX= FitPars["sigmaX"];
			double sigmaY= FitPars["sigmaY"];
			double theta= FitPars["theta"];			
			//return StatsUtils::GetFitEllipse(x0,y0,sigmaX,sigmaY,theta,useFWHM);

			TEllipse* ellipse= new TEllipse(x0,y0,sigmaX,sigmaY,0.,360.,theta);
			ellipse->SetLineWidth(2);
			ellipse->SetFillColor(0);
			ellipse->SetFillStyle(0);

			return ellipse;
		}//close GetFitEllipse()

		/**
		* \brief Get flux density
		*/
		double GetFluxDensity(){
			//Gaussian Area= pi*bmaj*bmin/(4*log(2)) or Area=2*pi*sigmaX*sigmaY
			double ampl= FitPars["A"];
			double sigmaX= FitPars["sigmaX"];
			double sigmaY= FitPars["sigmaY"];
			double fluxDensity= 2*TMath::Pi()*ampl*sigmaX*sigmaY;
			return fluxDensity;
		}

	private:
		
	private:
		std::map<std::string,double> FitPars;
		std::map<std::string,double> FitParsErr;
	
	ClassDef(SourceComponentPars,1)

};//close SourceComponentPars()


//========================================
//==         SOURCE FIT PARS
//========================================
class SourceFitPars : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFitPars(){
			Init();
		}
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFitPars(int N)
		{
			Init();
			SetNComponents(N);
		}	

		/**
		* \brief Class destructor: free allocated memory
		*/		
		virtual ~SourceFitPars() {}
		
	public:
		/**
		* \brief Set offset par value 
		*/
		void SetOffsetPar(double value){offset=value;}
		/**
		* \brief Get offset par value 
		*/
		double GetOffsetPar(){return offset;}

		/**
		* \brief Set offset par error
		*/
		void SetOffsetParErr(double value){offset_err=value;}
		/**
		* \brief Get offset par error value 
		*/
		double GetOffsetParErr(){return offset_err;}

		/**
		* \brief Get fitted ellipses
		*/
		std::vector<TEllipse*> GetFittedEllipses(bool useFWHM=true){
			std::vector<TEllipse*> ellipses;
			for(size_t i=0;i<pars.size();i++){
				TEllipse* ellipse= pars[i].GetFitEllipse(useFWHM);
				ellipses.push_back(ellipse);
			}
			return ellipses;
		}

		/**
		* \brief Get pars
		*/
		std::vector<SourceComponentPars>const& GetPars() const {return pars;}

		/**
		* \brief Set pars for component i-th
		*/
		int SetComponentPars(int componentId,SourceComponentPars& componentPars){
			if(componentId<0 || componentId>=(signed)pars.size()) return -1;
			pars[componentId]= componentPars;
			return 0;
		}

		/**
		* \brief Set par value & error
		*/
		int SetParValueAndError(int componentId,std::string parName,double parValue,double parError){
			if(componentId<0 || componentId>=nComponents) {
				WARN_LOG("Invalid component id ("<<componentId<<" given, cannot find par to be set!");
				return -1;
			}
			return pars[componentId].SetParValueAndError(parName,parValue,parError);
		}

		/** 
		\brief Get par value
 		*/
		double GetParValue(int componentId,std::string parName){
			if(componentId<0 || componentId>=nComponents) return -999;
			return pars[componentId].GetParValue(parName);
		}

		/** 
		\brief Get par error
 		*/
		double GetParError(int componentId,std::string parName){
			if(componentId<0 || componentId>=nComponents) return -999;
			return pars[componentId].GetParError(parName);
		}

		/**
		* \brief Initialize component pars
		*/
		void SetNComponents(int N){
			pars.clear();
			for(int i=0;i<N;i++) pars.push_back(SourceComponentPars());
			nComponents= N;
		}
	
		/**
		* \brief Get number of components
		*/
		int GetNComponents(){return nComponents;}

		/**
		* \brief Get component flux density
		*/
		double GetComponentFluxDensity(int componentId){
			if(componentId<0 || componentId>=nComponents) return -999;
			return pars[componentId].GetFluxDensity();
		}

		/**
		* \brief Set chi2 
		*/
		void SetChi2(double value){chi2=value;}

		/**
		* \brief Get chi2 
		*/
		double GetChi2(){return chi2;}
		/**
		* \brief Set NDF 
		*/
		void SetNDF(double value){ndof=value;}
		/**
		* \brief Get NDF 
		*/
		double GetNDF(){return ndof;}

		/**
		* \brief Set Status 
		*/
		void SetStatus(int value){status=value;}

		/**
		* \brief Get Status 
		*/
		int GetStatus(){return status;}

		/**
		* \brief Set Minimizer Status 
		*/
		void SetMinimizerStatus(int value){minimizer_status=value;}

		/**
		* \brief Get Minimizer Status 
		*/
		int GetMinimizerStatus(){return minimizer_status;}

		

		/**
		* \brief Set number of free pars 
		*/
		void SetNFreePars(int value){npars_free=value;}
		/**
		* \brief Get number of free parameters 
		*/
		int GetNFreePars(){return npars_free;}

		/**
		* \brief Set number of free pars 
		*/
		void SetNFitPoints(int value){nfit_points=value;}

		/**
		* \brief Get number of fit points
		*/
		int GetNFitPoints(){return nfit_points;}

		/**
		* \brief Get residual mean
		*/
		double GetResidualMean(){return residualMean;}
		/**
		* \brief Set residual mean
		*/
		void SetResidualMean(double value){residualMean=value;}
		/**
		* \brief Get residual rms
		*/
		double GetResidualRMS(){return residualRMS;}
		/**
		* \brief Set residual rms
		*/
		void SetResidualRMS(double value){residualRMS=value;}
		/**
		* \brief Get residual median
		*/
		double GetResidualMedian(){return residualMedian;}
		/**
		* \brief Set residual median
		*/
		void SetResidualMedian(double value){residualMedian=value;}
		/**
		* \brief Get residual mad
		*/
		double GetResidualMAD(){return residualMAD;}
		/**
		* \brief Set residual mad
		*/
		void SetResidualMAD(double value){residualMAD=value;}
		/**
		* \brief Get residual min
		*/
		double GetResidualMin(){return residualMin;}
		/**
		* \brief Set residual min
		*/
		void SetResidualMin(double value){residualMin=value;}		
		/**
		* \brief Get residual max
		*/
		double GetResidualMax(){return residualMax;}
		/**
		* \brief Set residual max
		*/
		void SetResidualMax(double value){residualMax=value;}

		/**
		* \brief Print
		*/
		void Print(){
			cout<<"*** FIT RESULTS ***"<<endl;
			cout<<"status="<<status<<" (minimizer status="<<minimizer_status<<")"<<endl;
			cout<<"chi2="<<chi2<<", ndf="<<ndof<<", redchi2="<<chi2/ndof<<endl;
			for(int i=0;i<nComponents;i++){
				cout<<"--> Component "<<i+1<<endl;
				cout<<"A="<<pars[i].GetParValue("A")<<" +- "<<pars[i].GetParError("A")<<endl;
				cout<<"(x0,y0)=("<<	pars[i].GetParValue("x0")<<","<<pars[i].GetParValue("y0")<<") err("<<pars[i].GetParError("x0")<<","<<pars[i].GetParError("y0")<<")"<<endl;
				cout<<"(sigmaX,sigmaY)=("<<	pars[i].GetParValue("sigmaX")<<","<<pars[i].GetParValue("sigmaY")<<") err("<<pars[i].GetParError("sigmaX")<<","<<pars[i].GetParError("sigmaY")<<")"<<endl;
				cout<<"theta="<<pars[i].GetParValue("theta")<<" +- "<<pars[i].GetParError("theta")<<endl;
			}
			cout<<"****************"<<endl;
		}

	private:

		/**
		* \brief Initialize component pars
		*/
		void Init(){
			nComponents= 0;
			chi2= 0;
			ndof= 0;
			npars_free= 0;
			nfit_points= 0;
			status= -1;
			minimizer_status= -1;
			offset= 0;
			offset_err= 0;
			residualMean= 0;
			residualRMS= 0;
			residualMedian= 0;
			residualMAD= 0;
			residualMin= 0;
			residualMax= 0;
			pars.clear();
		}

	private:
		int nComponents;
		double chi2;
		double ndof;
		int npars_free;
		int nfit_points;
		int status;
		int minimizer_status;
		double offset;
		double offset_err;
		double residualMean;
		double residualRMS;	
		double residualMedian;
		double residualMAD;
		double residualMin;	
		double residualMax;
		std::vector<SourceComponentPars> pars;

	ClassDef(SourceFitPars,1)

};


//========================================
//==         SOURCE FITTER
//========================================

class SourceFitter : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceFitter();
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceFitter();

		/**
		* \brief Fit status enum flag
		*/
		enum FitStatusFlag {
			eFitUnknownStatus= 0,
			eFitAborted= 1,
			eFitNotConverged= 2,
			eFitConverged= 3,
			eFitConvergedWithWarns= 4
		};

	public:

		/**
		* \brief Fit source
		*/
		int FitSource(Source* source,SourceFitOptions& fitOptions);
		
		/**
		* \brief Get fit pars
		*/
		SourceFitPars GetFitPars(){return m_sourceFitPars;}

		/**
		* \brief Get fit status
		*/
		int GetFitStatus(){return m_fitStatus;}

		/**
		* \brief 2D Gaussian mixture model used for the fit
		*/
		static double Gaus2DMixtureFcn(double* x, double* p);

		/**
		* \brief 2D Gaussian model used for the fit
		*/
		static double Gaus2DFcn(double* x, double* p);

	protected:
		/**
		* \brief Check if fit has parameters converged at bounds
		*/
		bool HasFitParsAtLimit(const ROOT::Fit::FitResult& fitRes);

	private:
	
		static int m_NFitComponents;
		static int m_fitStatus;
		static SourceFitPars m_sourceFitPars;
		
	ClassDef(SourceFitter,1)

};//close SourceFitter class

#ifdef __MAKECINT__
#pragma link C++ class SourceComponentPars+;
#pragma link C++ class SourceFitPars+;
#pragma link C++ class SourceFitter+;
#pragma link C++ enum FitStatusFlag;
#pragma link C++ struct SourceFitOptions+;
#endif

}//close namespace

#endif
