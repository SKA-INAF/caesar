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

struct BlobPars {
	double bmaj;//in pixels
	double bmin;//in pixels
	double bpa;//in deg
};

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
			pars.clear();
			//for(int i=0;i<nComponents;i++) pars.push_back(SourceComponentPars());
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
		int FitSource(Source* source,BlobPars blobPars,int nMaxComponents=3);
		
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
		//bool HasFitParsAtLimit(TFitResultPtr fitRes);
		bool HasFitParsAtLimit(const ROOT::Fit::FitResult& fitRes);

	private:
	
		static int m_NFitComponents;
		static int m_fitStatus;
		static SourceFitPars m_sourceFitPars;
		//static std::vector<TEllipse> m_fitEllipses;

	ClassDef(SourceFitter,1)

};//close SourceFitter class

#ifdef __MAKECINT__
#pragma link C++ class SourceComponentPars+;
#pragma link C++ class SourceFitPars+;
#pragma link C++ class SourceFitter+;
#pragma link C++ enum FitStatusFlag;
#endif

}//close namespace

#endif
