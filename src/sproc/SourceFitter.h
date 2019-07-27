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
* Class to fit a source image with a mixture of gaussian functions
* @author S. Riggi
* @date 01/09/2017
*/

#ifndef _SOURCE_FITTER_h
#define _SOURCE_FITTER_h 1

#include <SysUtils.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <AstroUtils.h>
#include <WCSUtils.h>
#include <Consts.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TFitResultPtr.h>
#include <Fit/FitResult.h>
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TMinuitMinimizer.h>

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

class TH2D;

namespace Caesar {

class Image;
class Source;

struct SourceFitOptions {

	public:

	
	//Default constructor
	SourceFitOptions() 
	{
    bmaj= 5;//pix
    bmin= 5;//pix
		bpa= 0;
		nMaxComponents= 3;	
		fixCentroidInPreFit= false;
		limitCentroidInFit= true;
		centroidLimit= 5;//pixels
		fixBkg= true;
		limitBkgInFit= true;
		useEstimatedBkgLevel= true;
		useBkgBoxEstimate= false;
		fixedBkgLevel= 0;
		fixAmplInPreFit= false;
		limitAmplInFit= true;
		amplLimit= 0.2;
		limitSigmaInFit= true;
		sigmaLimit= 0.2;
		fixSigmaInPreFit= true;
		fixSigma= false;
		fixThetaInPreFit= true;
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
		setBinErrorsToMapRMS= true;
		
		scaleMin= 3;
		scaleMax= 3;
		scaleStep= 1;
		minBlobSize= 5;
		blobMapThrFactor= 0;
		blobMapKernelFactor= 6; 
		useNestedAsComponents= false;
	
		fitImproveConvergence= true;
		fitNRetries= 1000;
		fitDoFinalMinimizerStep= true;
		fitFinalMinimizer= eHESS;
		chi2RegPar= 0.1;
		fitRetryWithLessComponents= true;

		//Supported minimizers & algo
		//  - Minuit: {Migrad,Simplex,Scan,Minimize}
		//  - Minuit2: {Migrad,Simplex,Combined,Scan,Fumili}
		//  - GSLMultiMin: {conjugatefr,conjugatepr,bfgs,bfgs2,steepestdescent} NB: No error calculation supported so do not use it
		//  - GSLMultiFit: {} NB: Requiring a different chi2 function definition, so do not use it
		fitMinimizer= "Minuit2";
		fitMinimizerAlgo= "Migrad";
		fitPrintLevel= 1;
		fitStrategy= 2;
		fitFcnTolerance= 1.e-2;//default tolerance used in ROOT
		fitMaxIters= 100000;
		fitParBoundIncreaseStepSize= 0.1;

		wcsType= eJ2000;
		fitScaleDataToMax= false;

		//Selection cuts
		useRedChi2Cut= true;
		fitRedChi2Cut= 5.;
		useFitEllipseCuts= false;
		fitEllipseEccentricityRatioMinCut= 0.5;
		fitEllipseEccentricityRatioMaxCut= 1.5;
		fitEllipseAreaRatioMinCut= 0.01;
		fitEllipseAreaRatioMaxCut= 10;
		fitEllipseRotAngleCut= 45;
		
		useSimpleWCSEllipseConversion= true;

	}//close constructor

	public:
		//- Blob start fit pars
		double bmaj;//in pixels
		double bmin;//in pixels
		double bpa;//in degrees

		//- Max number of components to be fitted & options
		int nMaxComponents;
		bool useNestedAsComponents;
	
		//- Number of matching peaks across kernels (-1=peak detected in all kernels, ...) and distance tolerance
		int peakKernelMultiplicityThr;
		int peakShiftTolerance;

		//- Dilation kernels to be used when finding peaks
		int peakMinKernelSize;
		int peakMaxKernelSize;

		//- Multiscale blob finder options
		double scaleMin;
		double scaleMax;
		double scaleStep;
		int minBlobSize;
		double blobMapThrFactor;
		int blobMapKernelFactor;

		//- Peak flux significance min threshold (in nsigmas wrt to avg bkg & rms)
		double peakZThrMin;	

		//- Centroid options
		bool limitCentroidInFit;
		double centroidLimit;//in pixels
		bool fixCentroidInPreFit;

		//- Bkg options
		bool fixBkg;
		bool limitBkgInFit;
		bool useEstimatedBkgLevel;
		bool useBkgBoxEstimate;
		double fixedBkgLevel;
		
		//- Amplitude fit par range (example +-20% around source peak)
		bool limitAmplInFit;
		double amplLimit;
		bool fixAmplInPreFit;

		//- Sigma fit par range
		bool fixSigmaInPreFit;
		bool fixSigma;	
		bool limitSigmaInFit;
		double sigmaLimit;

		//- Theta 
		bool fixThetaInPreFit;
		bool fixTheta;
		bool limitThetaInFit;
		double thetaLimit;//in deg

		//- Flux significance min threshold (in nsigmas above the bkg)
		bool useFluxZCut;
		double fluxZThrMin;

		//- Fit data bin error (if true set bin errors to pixel noise rms, otherwise to 1)
		bool setBinErrorsToMapRMS;

		//- Fit minimization options
		double fitFcnTolerance;
		long int fitMaxIters;
		bool fitImproveConvergence;
		long int fitNRetries;
		bool fitDoFinalMinimizerStep;
		int fitFinalMinimizer;
		double chi2RegPar;
		bool fitRetryWithLessComponents;

		//- Fit minimizer
		int fitPrintLevel;
		int fitStrategy;
		std::string fitMinimizer;
		std::string fitMinimizerAlgo;
		double fitParBoundIncreaseStepSize;		

		//- Fit ellipse pars
		int wcsType;
	
		//- Scale data to max
		bool fitScaleDataToMax;

		//- Selection cuts
		bool useRedChi2Cut;
		double fitRedChi2Cut;
		bool useFitEllipseCuts;
		double fitEllipseEccentricityRatioMinCut;
		double fitEllipseEccentricityRatioMaxCut;
		double fitEllipseAreaRatioMinCut;
		double fitEllipseAreaRatioMaxCut;
		double fitEllipseRotAngleCut;

		//- Use simple WCS ellipse conversion
		bool useSimpleWCSEllipseConversion;
		
};//close SourceFitOptions




//===========================================
//==         SOURCE COMPONENT PARS
//===========================================
class SourceComponentPars : public TObject {

	public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceComponentPars()
		{
			//Init pars with 0
			std::vector<std::string> parNames {"A","x0","y0","sigmaX","sigmaY","theta"};
			for(size_t i=0;i<parNames.size();i++){
				FitPars.insert( std::pair<std::string,double>(parNames[i],0.) );
				FitParsErr.insert( std::pair<std::string,double>(parNames[i],0.) );
			}

			//Initialize ellipse pars
			InitEllipsePars();

			//Init other vars
			m_flag= eCandidate;
			m_type= ePointLike;
			m_selected= true;

		}//close contructor

		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceComponentPars(){
			FitPars.clear();
			FitParsErr.clear();
		}

		/**
		* \brief Copy constructor
		*/
		SourceComponentPars(const SourceComponentPars& pars)
		{
			// Copy constructor
  		Init();
  		((SourceComponentPars&)pars).Copy(*this);
		}
		
		/**
		* \brief Assignment Operator
		*/
		SourceComponentPars& operator=(const SourceComponentPars& pars)
		{
			// Operator =
  		if (this != &pars) ((SourceComponentPars&)pars).Copy(*this);
  		return *this;
		}

		/**
		* \brief Copy method
		*/
		void Copy(TObject& obj) const
		{
			//Copy vars
			((SourceComponentPars&)obj).m_pixSize = m_pixSize;

			((SourceComponentPars&)obj).m_hasBeamPars = m_hasBeamPars;
			((SourceComponentPars&)obj).m_beam_bmaj = m_beam_bmaj;
			((SourceComponentPars&)obj).m_beam_bmin = m_beam_bmin;
			((SourceComponentPars&)obj).m_beam_pa = m_beam_pa;
			((SourceComponentPars&)obj).m_beam_area = m_beam_area;
			((SourceComponentPars&)obj).m_beam_eccentricity = m_beam_eccentricity;

			((SourceComponentPars&)obj).m_hasEllipsePars = m_hasEllipsePars;
			((SourceComponentPars&)obj).m_x0 = m_x0;
			((SourceComponentPars&)obj).m_y0 = m_y0;
			((SourceComponentPars&)obj).m_bmaj = m_bmaj;
			((SourceComponentPars&)obj).m_bmin = m_bmin;
			((SourceComponentPars&)obj).m_pa = m_pa;
			((SourceComponentPars&)obj).m_x0_err = m_x0_err;
			((SourceComponentPars&)obj).m_y0_err = m_y0_err;
			((SourceComponentPars&)obj).m_bmaj_err = m_bmaj_err;		
			((SourceComponentPars&)obj).m_bmin_err = m_bmin_err;
			((SourceComponentPars&)obj).m_pa_err = m_pa_err;
			((SourceComponentPars&)obj).m_eccentricity = m_eccentricity;
			((SourceComponentPars&)obj).m_area = m_area;
			((SourceComponentPars&)obj).m_rotangle_vs_beam = m_rotangle_vs_beam;

			((SourceComponentPars&)obj).m_hasWCSEllipsePars = m_hasWCSEllipsePars;
			((SourceComponentPars&)obj).m_x0_wcs = m_x0_wcs;
			((SourceComponentPars&)obj).m_y0_wcs = m_y0_wcs;
			((SourceComponentPars&)obj).m_bmaj_wcs = m_bmaj_wcs;
			((SourceComponentPars&)obj).m_bmin_wcs = m_bmin_wcs;
			((SourceComponentPars&)obj).m_pa_wcs = m_pa_wcs;
			((SourceComponentPars&)obj).m_x0_err_wcs = m_x0_err_wcs;
			((SourceComponentPars&)obj).m_y0_err_wcs = m_y0_err_wcs;
			((SourceComponentPars&)obj).m_bmaj_err_wcs = m_bmaj_err_wcs;
			((SourceComponentPars&)obj).m_bmin_err_wcs = m_bmin_err_wcs;
			((SourceComponentPars&)obj).m_pa_err_wcs = m_pa_err_wcs;
	
			((SourceComponentPars&)obj).m_hasWCSDeconvolvedEllipsePars = m_hasWCSDeconvolvedEllipsePars;
			((SourceComponentPars&)obj).m_bmaj_deconv_wcs = m_bmaj_deconv_wcs;
			((SourceComponentPars&)obj).m_bmin_deconv_wcs = m_bmin_deconv_wcs;
			((SourceComponentPars&)obj).m_pa_deconv_wcs = m_pa_deconv_wcs;
				
			((SourceComponentPars&)obj).m_flag = m_flag;
			((SourceComponentPars&)obj).m_type = m_type;
			((SourceComponentPars&)obj).m_selected = m_selected;
						
			//Copy maps
			((SourceComponentPars&)obj).FitPars = FitPars;
			((SourceComponentPars&)obj).FitParsErr = FitParsErr;
			
		}//close Copy()


	public:
		/** 
		\brief Set par value & error
 		*/
		int SetParValueAndError(std::string parName,double parVal,double parErr){
			if(!CodeUtils::HasMapKey(FitPars,parName)) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Invalid par name ("<<parName<<" given, cannot find par to be set!");
				#endif
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
		\brief Get source component flag
 		*/	
		int GetFlag(){return m_flag;}

		/** 
		\brief Set source component flag
 		*/
		void SetFlag(int flag){m_flag=flag;}

		/** 
		\brief Get source component type
 		*/	
		int GetType(){return m_type;}

		/** 
		\brief Set source component type
 		*/
		void SetType(int type){m_type=type;}

		/** 
		\brief Get source component selected flag
 		*/	
		bool IsSelected(){return m_selected;}
		/** 
		\brief Set source component selected flag
 		*/
		void SetSelected(bool val){m_selected=val;}


		/** 
		\brief Get fit ellipse
 		*/
		TEllipse* GetFitEllipse(bool useFWHM=true)
		{	
			//Check if has fit pars
			if(FitPars.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No fitted pars stored, returning nullptr ellipse!");
				#endif
				return nullptr;
			}
		
			//Compute fit ellipse from pars
			double x0= FitPars["x0"];
			double y0= FitPars["y0"];
			double sigmaX= FitPars["sigmaX"];
			double sigmaY= FitPars["sigmaY"];
			double theta= FitPars["theta"];	
			if(useFWHM){
				sigmaX*= GausSigma2FWHM/2.;
				sigmaY*= GausSigma2FWHM/2.;
			}
		
			TEllipse* ellipse= new TEllipse(x0,y0,sigmaX,sigmaY,0.,360.,theta);
			ellipse->SetLineWidth(2);
			ellipse->SetFillColor(0);
			ellipse->SetFillStyle(0);

			return ellipse;

		}//close GetFitEllipse()

		/** 
		\brief Has ellipse pars?
 		*/
		bool HasEllipsePars(){return m_hasEllipsePars;}

		/** 
		\brief Compute ellipse pars in pixel coordinates
 		*/
		int ComputeEllipsePars()
		{
			//Check if has fit pars
			if(FitPars.empty() || FitParsErr.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No fitted pars and/or errors stored!");
				#endif
				return -1;
			}

			//Compute ellipse pars
			m_x0= FitPars["x0"];
			m_y0= FitPars["y0"];
			m_x0_err= FitParsErr["x0"];
			m_y0_err= FitParsErr["y0"];
			double sigmaX= FitPars["sigmaX"];
			double sigmaY= FitPars["sigmaY"];
			double theta= FitPars["theta"];	
			double sigmaX_err= FitParsErr["sigmaX"];
			double sigmaY_err= FitParsErr["sigmaY"];
			double theta_err= FitParsErr["theta"];	

			double bmaj= sigmaX*GausSigma2FWHM;
			double bmin= sigmaY*GausSigma2FWHM;
			double pa= theta-90.;//convention is to measure bpa from North CCW (theta is measured from x axis)
			//double pa= MathUtils::Mod(theta-90.,180.);//convention is to measure bpa from North CCW (theta is measured from x axis)
			double bmaj_err= sigmaX_err*GausSigma2FWHM;
			double bmin_err= sigmaY_err*GausSigma2FWHM;
			double pa_err= theta_err;
			
			m_bmaj= bmaj;
			m_bmin= bmin;
			m_pa= pa;
			m_bmaj_err= bmaj_err;
			m_bmin_err= bmin_err;
			m_pa_err= pa_err;

			//Force bmaj>bmin
			if(bmaj<bmin){
				m_bmin= bmaj;//swap bmaj/bmin
				m_bmaj= bmin;		
				m_pa= pa + 90.;//rotate pa
				//m_pa= MathUtils::Mod(pa+90.,180.); 

				m_bmin_err= bmaj_err;
				m_bmaj_err= bmin_err;
				m_pa_err= pa_err; 
			}

			//Limit pa in range [-90,90]
			m_pa= GetPosAngleInRange(m_pa);

			
			//Compute ellipse eccentricity & area
			//NB: Using simple conversion to arcsec rather than bmaj_wcs, bmin_wcs
			double bmaj_arcsec= m_bmaj*m_pixSize;
			double bmin_arcsec= m_bmin*m_pixSize;
			m_eccentricity= MathUtils::ComputeEllipseEccentricity(bmaj_arcsec,bmin_arcsec);
			m_area= MathUtils::ComputeEllipseArea(bmaj_arcsec,bmin_arcsec);			

			//Compute rotation angle vs beam (if beam info is available)
			if(m_hasBeamPars){
				double dtheta= m_pa-m_beam_pa;	
				//m_rotangle_vs_beam= MathUtils::Mod(dtheta,180.);
				m_rotangle_vs_beam= dtheta;
			}
			else{
				#ifdef LOGGING_ENABLED
					WARN_LOG("No beam information has been set, do not compute fit ellipse rot angle vs beam (set to 0 by default)!");
				#endif
			}
			

			//Set has ellipse par flag
			m_hasEllipsePars= true;

			return 0;

		}//close ComputeEllipsePars()


		/** 
		\brief Get fit ellipse pars.
 		*/
		int GetEllipsePars(double& x0,double& y0,double& bmaj,double& bmin,double& pa)
		{
			//Init pars
			x0= 0;
			y0= 0;
			bmaj= 0;
			bmin= 0;
			pa= 0;

			//Check if has pars computed
			if(!m_hasEllipsePars) {	
				#ifdef LOGGING_ENABLED
					WARN_LOG("Fitted ellipse pars not computed, return dummy values!");
				#endif
				return -1;
			}

			x0= m_x0;
			y0= m_y0;
			bmaj= m_bmaj;
			bmin= m_bmin;
			pa= m_pa;

			return 0;

		}//close GetEllipsePars()

		/** 
		\brief Set fit ellipse pars. Used for serialization.
 		*/
		void SetEllipsePars(double x0,double y0,double bmaj,double bmin,double pa)	
		{
			m_x0= x0;
			m_y0= y0;
			m_bmaj= bmaj;
			m_bmin= bmin;
			m_pa= pa;
			m_hasEllipsePars= true;
		}
		

		/** 
		\brief Get fit ellipse par errors
 		*/
		int GetEllipseParErrors(double& x0_err,double& y0_err,double& bmaj_err,double& bmin_err,double& pa_err)
		{
			//Init pars
			x0_err= 0;
			y0_err= 0;
			bmaj_err= 0;
			bmin_err= 0;
			pa_err= 0;

			//Check if has pars computed
			if(!m_hasEllipsePars) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Fitted ellipse pars not computed, return dummy values!");
				#endif
				return -1;
			}

			x0_err= m_x0_err;
			y0_err= m_y0_err;
			bmaj_err= m_bmaj_err;
			bmin_err= m_bmin_err;
			pa_err= m_pa_err;

			return 0;

		}//close GetEllipseParErrors()
		
		/** 
		\brief Set fit ellipse par errors. Used for serialization.
 		*/
		void SetEllipseParErrors(double x0_err,double y0_err,double bmaj_err,double bmin_err,double pa_err)	
		{
			m_x0_err= x0_err;
			m_y0_err= y0_err;
			m_bmaj_err= bmaj_err;
			m_bmin_err= bmin_err;
			m_pa_err= pa_err;
		}

		/** 
		\brief Get fitted ellipse eccentricity
 		*/
		double GetEllipseEccentricity(){return m_eccentricity;}
		/** 
		\brief Set fitted ellipse eccentricity
 		*/	
		void SetEllipseEccentricity(double val){m_eccentricity=val;}

		/** 
		\brief Get fitted ellipse area
 		*/
		double GetEllipseArea(){return m_area;}
		/** 
		\brief Set fitted ellipse area
 		*/	
		void SetEllipseArea(double val){m_area=val;}

		/** 
		\brief Get fitted ellipse rot angle wrt beam
 		*/
		double GetEllipseRotAngleVSBeam(){return m_rotangle_vs_beam;}
		/** 
		\brief Set fitted ellipse rot angle wrt beam
 		*/	
		void SetEllipseRotAngleVSBeam(double val){m_rotangle_vs_beam=val;}

		/**
		* \brief Compute ellipse pars in WCS coordinates
		*/
		int ComputeWCSEllipsePars(WCS* wcs)
		{
			//Check given WCS	
			if(!wcs){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Null ptr to WCS given!");	
				#endif
				return -1;
			}

			//Check if has fit pars
			if(FitPars.empty() || FitParsErr.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No fitted pars stored!");
				#endif
				return -1;
			}

			//Get fitted pars in pixel coordinates
			double x= FitPars["x0"];
			double y= FitPars["y0"];
			double sx= FitPars["sigmaX"];
			double sy= FitPars["sigmaY"];
			double theta= FitPars["theta"];	
			double theta_rad= theta*TMath::DegToRad();

			//Get fitted pars errors in pixel coordinates
			double x_err= FitParsErr["x0"];
			double y_err= FitParsErr["y0"]; 
			double sx_err= FitParsErr["sigmaX"];
			double sy_err= FitParsErr["sigmaY"];
			double theta_err= FitParsErr["theta"];
			double theta_err_rad= theta_err*TMath::DegToRad();

			//Compute ellipse centroid coords in WCS coords
			double x_wcs= 0;
			double y_wcs= 0;
			WCSUtils::pix2wcs (wcs,x,y,&x_wcs, &y_wcs);
			m_x0_wcs= x_wcs;
			m_y0_wcs= y_wcs;

			//Compute ellipse centroid errors in WCS coords
			double dx_wcs= 0;
			double dy_wcs= 0;
			WCSUtils::pix2wcs (wcs,x + x_err,y + y_err,&dx_wcs,&dy_wcs); 
			m_x0_err_wcs= AstroUtils::GetWCSPointDist_Haversine(x_wcs,y_wcs,dx_wcs,y_wcs);
			m_y0_err_wcs= AstroUtils::GetWCSPointDist_Haversine(x_wcs,y_wcs,x_wcs,dy_wcs);


			//Compute ellipse axis corner pixel coords
			double x1= x + sx * cos(theta_rad);
			double y1= y + sx * sin(theta_rad);
			double x2= x + sy * cos(theta_rad - TMath::Pi()/2.);
			double y2= y + sy * sin(theta_rad - TMath::Pi()/2.);


			//Convert ellipse axis corner WCS coords
			double x1_wcs= 0;
			double y1_wcs= 0;
			double x2_wcs= 0;
			double y2_wcs= 0;
			WCSUtils::pix2wcs (wcs,x1,y1,&x1_wcs, &y1_wcs);
			WCSUtils::pix2wcs (wcs,x2,y2,&x2_wcs, &y2_wcs);

			//Compute semi-axis in WCS coords
			double a_wcs= AstroUtils::GetWCSPointDist_Haversine(x_wcs,y_wcs,x1_wcs,y1_wcs);//semi-major axis
			double b_wcs= AstroUtils::GetWCSPointDist_Haversine(x_wcs,y_wcs,x2_wcs,y2_wcs);//semi-major axis		

			//Compute ellipse rot angle
			//NB: Bearing is returned from North to East (0 deg is north)
			//# From AEGEAN source finder: The a/b vectors are perpendicular in sky space, but not always in pixel space
      //# so we have to account for this by calculating the angle between the two vectors
      //# and modifying the minor axis length
			double pa_wcs = AstroUtils::GetWCSPointBearing(x_wcs,y_wcs,x1_wcs,y1_wcs);
			double pa2_wcs = AstroUtils::GetWCSPointBearing(x_wcs,y_wcs,x2_wcs,y2_wcs) - 90.;
			double defect= pa_wcs-pa2_wcs;
			
			//Correct ellipse minor axis
			b_wcs*= fabs(cos(defect*TMath::DegToRad()));
			
			//Compute axis errors
			//- Major axis
			double x1_ref= x + sx * cos(theta_rad); 
			double y1_ref= y + sy * sin(theta_rad); 
			double x1_err= x + (sx + sx_err) * cos(theta_rad);
			double y1_err= y + sy * sin(theta_rad);	
			double x1_ref_wcs= 0;
			double y1_ref_wcs= 0;
			double x1_err_wcs= 0;
			double y1_err_wcs= 0;
			WCSUtils::pix2wcs (wcs,x1_ref,y1_ref,&x1_ref_wcs, &y1_ref_wcs);
			WCSUtils::pix2wcs (wcs,x1_err,y1_err,&x1_err_wcs, &y1_err_wcs);
			double a_err_wcs= AstroUtils::GetWCSPointDist_Haversine(x1_ref_wcs,y1_ref_wcs,x1_err_wcs,y1_err_wcs);

			//- Minor axis
			double x2_ref= x + sx * cos(theta_rad + TMath::Pi()/2.);
			double y2_ref= y + sy * sin(theta_rad + TMath::Pi()/2.);
			double x2_err= x + sx * cos(theta_rad + TMath::Pi()/2.);
			double y2_err= y + (sy + sy_err) * sin(theta_rad + TMath::Pi()/2.);
			double x2_ref_wcs= 0;
			double y2_ref_wcs= 0;			
			double x2_err_wcs= 0;
			double y2_err_wcs= 0;
			WCSUtils::pix2wcs (wcs,x2_ref,y2_ref,&x2_ref_wcs, &y2_ref_wcs);
			WCSUtils::pix2wcs (wcs,x2_err,y2_err,&x2_err_wcs, &y2_err_wcs);
			double b_err_wcs= AstroUtils::GetWCSPointDist_Haversine(x2_ref_wcs,y2_ref_wcs,x2_err_wcs,y2_err_wcs);

			//Compute pos angle error
			double x1_thetaerr= x + sx * cos(theta_rad + theta_err_rad);
			double y1_thetaerr= y + sy * sin(theta_rad + theta_err_rad);
			double x1_thetaerr_wcs= 0;
			double y1_thetaerr_wcs= 0;
			WCSUtils::pix2wcs (wcs,x1_thetaerr,y1_thetaerr,&x1_thetaerr_wcs, &y1_thetaerr_wcs);
			double pa_err_wcs= fabs(AstroUtils::GetWCSPointBearing(x_wcs,y_wcs,x1_wcs,y1_wcs) - AstroUtils::GetWCSPointBearing(x_wcs,y_wcs,x1_thetaerr_wcs,y1_thetaerr_wcs));

			//Set bmaj/bmin (check at the end for shape)
			double bmaj= a_wcs*2;//x 2 to get axis 
			double bmin= b_wcs*2; 
			double pa= pa_wcs;
			double bmaj_err= a_err_wcs*2;
			double bmin_err= b_err_wcs*2;
			double pa_err= pa_err_wcs;
			m_bmaj_wcs= bmaj;
			m_bmin_wcs= bmin;	
			m_pa_wcs= pa;
			m_bmaj_err_wcs= bmaj_err;
			m_bmin_err_wcs= bmin_err;
			m_pa_err_wcs= pa_err;

			//Force bmaj>bmin (if bmaj<bmin, swap bmaj/bmin and increment pa by 90 deg)
			if(bmaj<bmin){
				m_bmaj_wcs= bmin;
				m_bmin_wcs= bmaj;
				m_pa_wcs+= 90.;
				m_bmaj_err_wcs= bmin_err;
				m_bmin_err_wcs= bmaj_err;
			}			

			//Convert in arcsec
			m_bmaj_wcs*= 3600;
			m_bmin_wcs*= 3600;
			m_bmaj_err_wcs*= 3600;
			m_bmin_err_wcs*= 3600;

			//Limit pa in [-90,90]
			m_pa_wcs= GetPosAngleInRange(m_pa_wcs);

			/*
			//Compute ellipse eccentricity & area
			m_eccentricity= MathUtils::ComputeEllipseEccentricity(m_bmaj_wcs,m_bmin_wcs);
			m_area= MathUtils::ComputeEllipseArea(m_bmaj_wcs,m_bmin_wcs);			

			//Compute rotation angle vs beam (if beam info is available)
			if(m_hasBeamPars){
				double dtheta= m_pa_wcs-m_beam_pa;	
				//m_rotangle_vs_beam= MathUtils::Mod(dtheta,180.);
				m_rotangle_vs_beam= dtheta;
			}
			else{
				#ifdef LOGGING_ENABLED
					WARN_LOG("No beam information has been set, do not compute fit ellipse rot angle vs beam (set to 0 by default)!");
				#endif
			}
			*/

			//Set has WCS ellipse par flag
			m_hasWCSEllipsePars= true;

			return 0;

		}//close ComputeWCSEllipsePars()

		/**
		* \brief Compute ellipse pars in WCS coordinates, neglecting sky projection and assuming Euclidean distances
		*/
		int ComputeWCSEllipseParsSimple(WCS* wcs)
		{
			//Check given WCS	
			if(!wcs){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Null ptr to WCS given!");	
				#endif
				return -1;
			}

			//Check if has fit pars
			if(FitPars.empty() || FitParsErr.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No fitted pars stored!");
				#endif
				return -1;
			}

			//Get fitted pars in pixel coordinates
			double x= FitPars["x0"];
			double y= FitPars["y0"];
			double sx= FitPars["sigmaX"];
			double sy= FitPars["sigmaY"];
			double theta= FitPars["theta"];	
			double theta_rad= theta*TMath::DegToRad();

			//Get fitted pars errors in pixel coordinates
			double x_err= FitParsErr["x0"];
			double y_err= FitParsErr["y0"]; 
			double sx_err= FitParsErr["sigmaX"];
			double sy_err= FitParsErr["sigmaY"];
			double theta_err= FitParsErr["theta"];
			double theta_err_rad= theta_err*TMath::DegToRad();

			//Compute ellipse centroid coords in WCS coords
			double x_wcs= 0;
			double y_wcs= 0;
			WCSUtils::pix2wcs (wcs,x,y,&x_wcs, &y_wcs);
			m_x0_wcs= x_wcs;
			m_y0_wcs= y_wcs;

			//Compute ellipse centroid errors in WCS coords
			double dx_wcs= 0;
			double dy_wcs= 0;
			WCSUtils::pix2wcs (wcs,x + x_err,y + y_err,&dx_wcs,&dy_wcs); 
			m_x0_err_wcs= fabs(x_wcs-dx_wcs);
			m_y0_err_wcs= fabs(y_wcs-dy_wcs);


			//Compute ellipse axis corner pixel coords
			double x1= x + sx * cos(theta_rad);
			double y1= y + sx * sin(theta_rad);
			double x2= x + sy * cos(theta_rad - TMath::Pi()/2.);
			double y2= y + sy * sin(theta_rad - TMath::Pi()/2.);


			//Convert ellipse axis corner WCS coords
			double x1_wcs= 0;
			double y1_wcs= 0;
			double x2_wcs= 0;
			double y2_wcs= 0;
			WCSUtils::pix2wcs (wcs,x1,y1,&x1_wcs, &y1_wcs);
			WCSUtils::pix2wcs (wcs,x2,y2,&x2_wcs, &y2_wcs);

			//Compute semi-axis in WCS coords
			double a_wcs= sqrt( (x_wcs-x1_wcs)*(x_wcs-x1_wcs) + (y_wcs-y1_wcs)*(y_wcs-y1_wcs) );//semi-major axis
			double b_wcs= sqrt( (x_wcs-x2_wcs)*(x_wcs-x2_wcs) + (y_wcs-y2_wcs)*(y_wcs-y2_wcs) );//semi-major axis		

			//Compute ellipse rot angle
			double pa_wcs= theta-90.;//convention is to measure bpa from North CCW (theta is measured from x axis)

			
			//Compute axis errors
			//- Major axis
			double x1_ref= x + sx * cos(theta_rad); 
			double y1_ref= y + sy * sin(theta_rad); 
			double x1_err= x + (sx + sx_err) * cos(theta_rad);
			double y1_err= y + sy * sin(theta_rad);	
			double x1_ref_wcs= 0;
			double y1_ref_wcs= 0;
			double x1_err_wcs= 0;
			double y1_err_wcs= 0;
			WCSUtils::pix2wcs (wcs,x1_ref,y1_ref,&x1_ref_wcs, &y1_ref_wcs);
			WCSUtils::pix2wcs (wcs,x1_err,y1_err,&x1_err_wcs, &y1_err_wcs);
			double a_err_wcs= sqrt( (x1_ref_wcs-x1_err_wcs)*(x1_ref_wcs-x1_err_wcs) + (y1_ref_wcs-y1_err_wcs)*(y1_ref_wcs-y1_err_wcs) );

			//- Minor axis
			double x2_ref= x + sx * cos(theta_rad + TMath::Pi()/2.);
			double y2_ref= y + sy * sin(theta_rad + TMath::Pi()/2.);
			double x2_err= x + sx * cos(theta_rad + TMath::Pi()/2.);
			double y2_err= y + (sy + sy_err) * sin(theta_rad + TMath::Pi()/2.);
			double x2_ref_wcs= 0;
			double y2_ref_wcs= 0;			
			double x2_err_wcs= 0;
			double y2_err_wcs= 0;
			WCSUtils::pix2wcs (wcs,x2_ref,y2_ref,&x2_ref_wcs, &y2_ref_wcs);
			WCSUtils::pix2wcs (wcs,x2_err,y2_err,&x2_err_wcs, &y2_err_wcs);
			//double b_err_wcs= AstroUtils::GetWCSPointDist_Haversine(x2_ref_wcs,y2_ref_wcs,x2_err_wcs,y2_err_wcs);
			double b_err_wcs= sqrt( (x2_ref_wcs-x2_err_wcs)*(x2_ref_wcs-x2_err_wcs) + (y2_ref_wcs-y2_err_wcs)*(y2_ref_wcs-y2_err_wcs) );

			//Compute pos angle error
			double pa_err_wcs= theta_err;

			//Set bmaj/bmin (check at the end for shape)
			double bmaj= a_wcs*2;//x 2 to get axis 
			double bmin= b_wcs*2; 
			double pa= pa_wcs;
			double bmaj_err= a_err_wcs*2;
			double bmin_err= b_err_wcs*2;
			double pa_err= pa_err_wcs;
			m_bmaj_wcs= bmaj;
			m_bmin_wcs= bmin;	
			m_pa_wcs= pa;
			m_bmaj_err_wcs= bmaj_err;
			m_bmin_err_wcs= bmin_err;
			m_pa_err_wcs= pa_err;

			//Force bmaj>bmin (if bmaj<bmin, swap bmaj/bmin and increment pa by 90 deg)
			if(bmaj<bmin){
				m_bmaj_wcs= bmin;
				m_bmin_wcs= bmaj;
				m_pa_wcs+= 90.;
				m_bmaj_err_wcs= bmin_err;
				m_bmin_err_wcs= bmaj_err;
			}			

			//Convert in arcsec
			m_bmaj_wcs*= 3600;
			m_bmin_wcs*= 3600;
			m_bmaj_err_wcs*= 3600;
			m_bmin_err_wcs*= 3600;

			//Limit pa in [-90,90]
			m_pa_wcs= GetPosAngleInRange(m_pa_wcs);

			//Set has WCS ellipse par flag
			m_hasWCSEllipsePars= true;

			return 0;

		}//close ComputeWCSEllipseParsSimple()


		/** 
		\brief Has WCS ellipse pars?
 		*/
		bool HasWCSEllipsePars(){return m_hasWCSEllipsePars;}

		/**
		* \brief Get ellipse pars in WCS coordinates
		*/
		int GetWCSEllipsePars(double& x0_wcs,double& y0_wcs,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
		{
			//Init 
			x0_wcs= 0;
			y0_wcs= 0;
			bmaj_wcs= 0;
			bmin_wcs= 0;
			pa_wcs= 0;
			
			//Check if has fit pars
			if(!m_hasWCSEllipsePars) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No WCS ellipse pars was computed!");
				#endif
				return -1;
			}

			x0_wcs= m_x0_wcs;
			y0_wcs= m_y0_wcs;
			bmaj_wcs= m_bmaj_wcs;
			bmin_wcs= m_bmin_wcs;
			pa_wcs= m_pa_wcs;

			return 0;

		}//close GetWCSEllipsePars()


		/**
		* \brief Set ellipse pars in WCS coordinates
		*/
		void SetWCSEllipsePars(double x0_wcs,double y0_wcs,double bmaj_wcs,double bmin_wcs,double pa_wcs)
		{
			m_x0_wcs= x0_wcs;
			m_y0_wcs= y0_wcs;
			m_bmaj_wcs= bmaj_wcs;
			m_bmin_wcs= bmin_wcs;
			m_pa_wcs= pa_wcs;
			m_hasWCSEllipsePars= true;
		}

		/** 
		\brief Get fit ellipse par errors
 		*/
		int GetWCSEllipseParErrors(double& x0_err_wcs,double& y0_err_wcs,double& bmaj_err_wcs,double& bmin_err_wcs,double& pa_err_wcs)
		{
			//Init pars
			x0_err_wcs= 0;
			y0_err_wcs= 0;
			bmaj_err_wcs= 0;
			bmin_err_wcs= 0;
			pa_err_wcs= 0;

			//Check if has pars computed
			if(!m_hasWCSEllipsePars) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("WCS ellipse pars not computed, return dummy values!");
				#endif
				return -1;
			}

			x0_err_wcs= m_x0_err_wcs;
			y0_err_wcs= m_y0_err_wcs;
			bmaj_err_wcs= m_bmaj_err_wcs;
			bmin_err_wcs= m_bmin_err_wcs;
			pa_err_wcs= m_pa_err_wcs;

			return 0;

		}//close GetWCSEllipseParErrors()
		
		/**
		* \brief Set ellipse par errors in WCS coordinates
		*/
		void SetWCSEllipseParErrors(double x0_err_wcs,double y0_err_wcs,double bmaj_err_wcs,double bmin_err_wcs,double pa_err_wcs)
		{
			m_x0_err_wcs= x0_err_wcs;
			m_y0_err_wcs= y0_err_wcs;
			m_bmaj_err_wcs= bmaj_err_wcs;
			m_bmin_err_wcs= bmin_err_wcs;
			m_pa_err_wcs= pa_err_wcs;
		}

		/**
		* \brief Set image pix size (it is assumed units are arcsec)
		*/
		void SetImagePixSize(double val){m_pixSize=val;}
		/**
		* \brief Get image pix size
		*/
		double GetImagePixSize(){return m_pixSize;}

		/**
		* \brief Set beam ellipse parameters (it is assumed units are arcsec)
		*/
		void SetBeamEllipsePars(double bmaj,double bmin,double pa)
		{
			m_hasBeamPars= true;
			m_beam_bmaj= bmaj;
			m_beam_bmin= bmin;
			m_beam_pa= pa;

			//Force bmaj>bmin
			if(bmaj<bmin){
				m_beam_bmin= bmaj;//swap bmaj/bmin
				m_beam_bmaj= bmin;		
				m_beam_pa= pa + 90.;//rotate pa
			}

			//Limit pa in range [-90,90]
			double pa_limited= GetPosAngleInRange(m_beam_pa);
			m_beam_pa= pa_limited;

			//Compute ellipse eccentricity & area
			m_beam_eccentricity= MathUtils::ComputeEllipseEccentricity(m_beam_bmaj,m_beam_bmin);
			m_beam_area= MathUtils::ComputeEllipseArea(m_beam_bmaj,m_beam_bmin);			
	
		}
	
		/**
		* \brief Has beam pars?
		*/
		bool HasBeamEllipsePars(){return m_hasBeamPars;}

		/**
		* \brief Get beam ellipse parameters
		*/
		int GetBeamEllipsePars(double& bmaj,double& bmin,double& pa)
		{
			//Check if has beam pars
			if(!m_hasBeamPars){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No beam pars stored!");
				#endif
				return -1;
			}
			bmaj= m_beam_bmaj;
			bmin= m_beam_bmin;
			pa= m_beam_pa;
			return 0;
		};

		/**
		* \brief Get beam ellipse eccentricity
		*/
		double GetBeamEllipseEccentricity(){return m_beam_eccentricity;}
		/**
		* \brief Set beam ellipse eccentricity
		*/
		void SetBeamEllipseEccentricity(double val){m_beam_eccentricity=val;}

		/**
		* \brief Get beam ellipse area
		*/
		double GetBeamEllipseArea(){return m_beam_area;}
		/**
		* \brief Set beam ellipse area
		*/
		void SetBeamEllipseArea(double val){m_beam_area=val;}
	
		/**
		* \brief Compute beam-deconvolved ellipse in WCS coords
		*/
		int ComputeWCSDeconvolvedEllipsePars()
		{
			int status= 0;
	
			//Check if has fit pars
			if(FitPars.empty() || FitParsErr.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No fitted pars stored!");
				#endif
				return -1;
			}
			//Check if has WCS ellipse pars
			if(!m_hasWCSEllipsePars) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No WCS ellipse pars was computed!");
				#endif
				return -1;
			}
			//Check if has ellipse beam pars
			if(!m_hasBeamPars) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("No beam pars are available!");
				#endif
				return -1;
			}
	
			//Compute deconvolved ellipse pars using formula of Wild 1970 (see also Mon. Not. R. Astron. Soc. 342, 1117–1130 (2003))			
			double sum2= m_bmaj_wcs*m_bmaj_wcs + m_bmin_wcs*m_bmin_wcs;
			double diff2= m_bmaj_wcs*m_bmaj_wcs - m_bmin_wcs*m_bmin_wcs;
			double sum2_beam= m_beam_bmaj*m_beam_bmaj + m_beam_bmin*m_beam_bmin;
			double diff2_beam= m_beam_bmaj*m_beam_bmaj - m_beam_bmin*m_beam_bmin;
			double pa_rad= m_pa_wcs*TMath::DegToRad();//in rad
			double pa_beam_rad= m_beam_pa*TMath::DegToRad();//in rad
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("ellipse pars("<<m_bmaj<<","<<m_bmin<<","<<m_pa<<"), wcs("<<m_bmaj_wcs<<","<<m_bmin_wcs<<","<<m_pa_wcs<<"), beam("<<m_beam_bmaj<<","<<m_beam_bmin<<","<<m_beam_pa<<")");
			#endif

			double beta2= pow(diff2,2) + pow(diff2_beam,2) - 2*diff2*diff2_beam*cos(2*(pa_rad-pa_beam_rad));
			if(beta2<0) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Numerical error (beta^2 is <0 in formula)!");
				#endif
				return -1;
			}
			double beta= sqrt(beta2);

			//Check if fit ellipse is smaller than beam.
			//If so, do not attempt too deconvolve. Set deconv ellipse to fitted ellipse.
			if(sum2<=sum2_beam){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Fitted ellipse beam is smaller than beam, do not deconvolve (set deconvolved ellipse to fitted ellipse)!");
				#endif
				m_bmaj_deconv_wcs= m_bmaj_wcs;
				m_bmin_deconv_wcs= m_bmin_wcs;
				m_pa_deconv_wcs= m_pa_wcs;
				m_hasWCSDeconvolvedEllipsePars= true;
				return 0;
			}

	
			double bmaj2_deconv= 0.5*(sum2 - sum2_beam + beta);
			double bmin2_deconv= 0.5*(sum2 - sum2_beam - beta);
			if(bmaj2_deconv<0 || bmin2_deconv) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Numerical error (bmaj^2/bmin^2 deconvolved are <0 in formula)!");	
				#endif
				return -1;
			}
			m_bmaj_deconv_wcs= sqrt(bmaj2_deconv);
			m_bmin_deconv_wcs= sqrt(bmin2_deconv);

			double arg= (diff2*sin(2*m_pa_wcs))/(diff2*cos(2*m_pa_wcs-diff2_beam));
			m_pa_deconv_wcs= 0.5*atan(arg)*TMath::RadToDeg();

			//Compute deconvolved ellipse par errors by error propagation
			//... WRITE ME ...
			//...
			//...
	

			m_hasWCSDeconvolvedEllipsePars= true;

			return 0;

		}//close ComputeWCSDeconvolvedEllipsePars()


		/**
		* \brief Has beam-deconvolved WCS ellipse pars?
		*/
		bool HasWCSDeconvolvedEllipsePars(){return m_hasWCSDeconvolvedEllipsePars;}


		/**
		* \brief Get ellipse pars in WCS coordinates
		*/
		int GetWCSDeconvolvedEllipsePars(double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
		{
			//Init 
			bmaj_wcs= 0;
			bmin_wcs= 0;
			pa_wcs= 0;
			
			//Check if has fit pars
			if(!m_hasWCSDeconvolvedEllipsePars) {
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("No WCS beam-deconvolved ellipse pars was computed!");
				#endif
				return -1;
			}

			bmaj_wcs= m_bmaj_deconv_wcs;
			bmin_wcs= m_bmin_deconv_wcs;
			pa_wcs= m_pa_deconv_wcs;

			return 0;

		}//close GetWCSDeconvolvedEllipsePars()
 
		void SetWCSDeconvolvedEllipsePars(double bmaj_wcs,double bmin_wcs,double pa_wcs)
		{
			m_bmaj_deconv_wcs= bmaj_wcs;
			m_bmin_deconv_wcs= bmin_wcs;
			m_pa_deconv_wcs= pa_wcs;
		}

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

		/**
		* \brief Limit position angle in [-90,90]. It is periodic with period 180 deg.
		*/
		double GetPosAngleInRange(double pa)
		{
    	while (pa <= -90) pa += 180;
    	while (pa > 90) pa -= 180;
    	return pa;
		}

		/**
		* \brief Init ellipse pars
		*/
		void InitEllipsePars()
		{
			m_hasEllipsePars= false;
			m_x0= -999;
			m_y0= -999;
			m_bmaj= -999;
			m_bmin= -999;
			m_pa= -999;
			m_x0_err= -999;
			m_y0_err= -999;
			m_bmaj_err= -999;
			m_bmin_err= -999;
			m_pa_err= -999;

			m_eccentricity= -999;
			m_area= -999;
			m_rotangle_vs_beam= 0;//NB: Set to 0 if beam info not available

			m_hasWCSEllipsePars= false;
			m_x0_wcs= -999;
			m_y0_wcs= -999;	
			m_bmaj_wcs= -999;
			m_bmin_wcs= -999;
			m_pa_wcs= -999;
			m_x0_err_wcs= -999;
			m_y0_err_wcs= -999;
			m_bmaj_err_wcs= -999;
			m_bmin_err_wcs= -999;
			m_pa_err_wcs= -999;

			//m_eccentricity_wcs= -999;
			//m_area_wcs= -999;
			//m_rotangle_vs_beam_wcs= 0;//NB: Set to 0 if beam info not available

			m_hasBeamPars= false;
			m_beam_bmaj= -999;
			m_beam_bmin= -999;
			m_beam_pa= -999;
			m_beam_area= -999;
			m_beam_eccentricity= -999;
			
			m_hasWCSDeconvolvedEllipsePars= false;
		 	m_bmaj_deconv_wcs= -999;
			m_bmin_deconv_wcs= -999;
			m_pa_deconv_wcs= -999;

			m_pixSize= 0;

		}//close InitEllipsePars()

		/**
		* \brief Init class data
		*/
		void Init()
		{
			//Init pars with 0
			std::vector<std::string> parNames {"A","x0","y0","sigmaX","sigmaY","theta"};
			for(size_t i=0;i<parNames.size();i++){
				FitPars.insert( std::pair<std::string,double>(parNames[i],0.) );
				FitParsErr.insert( std::pair<std::string,double>(parNames[i],0.) );
			}

			//Initialize ellipse pars
			InitEllipsePars();

			//Init other vars
			m_flag= eCandidate;
			m_type= ePointLike;
			m_selected= true;
		}
		
	private:

		std::map<std::string,double> FitPars;
		std::map<std::string,double> FitParsErr;
	
		//- Image pars
		double m_pixSize;//in arcsec

		//- Beam pars
		bool m_hasBeamPars;
		double m_beam_bmaj;
		double m_beam_bmin;
		double m_beam_pa;
		double m_beam_area;
		double m_beam_eccentricity;
				
		//- Ellipse pars
		bool m_hasEllipsePars;
		double m_x0;//ellipse x centroid in pixel coordinates
		double m_y0;//ellipse y centroid in pixel coordinates
		double m_bmaj;//ellipse major axis in pixel coordinates
		double m_bmin;//ellipse minor axis in pixel coordinates
		double m_pa;//ellipse position angle (CCW from North) in deg
		double m_x0_err;
		double m_y0_err;
		double m_bmaj_err;
		double m_bmin_err;
		double m_pa_err;

		double m_eccentricity;//ellipse eccentricity
		double m_area;//ellipse area
		double m_rotangle_vs_beam;//rotation angle vs beam

		//- WCS ellipse pars
		bool m_hasWCSEllipsePars;
		double m_x0_wcs;//ellipse x centroid in world coordinates (in deg)
		double m_y0_wcs;//ellipse y centroid in world coordinates (in deg)
		double m_bmaj_wcs;//ellipse major axis in world coordinates (in arcsec)
		double m_bmin_wcs;//ellipse minor axis in world coordinates (in arcsec)
		double m_pa_wcs;//ellipse position angle (CCW from North) in deg
		double m_x0_err_wcs;
		double m_y0_err_wcs;
		double m_bmaj_err_wcs;
		double m_bmin_err_wcs;
		double m_pa_err_wcs;

		//double m_eccentricity_wcs;//ellipse eccentricity
		//double m_area_wcs;//ellipse area
		//double m_rotangle_vs_beam_wcs;//rotation angle vs beam

		//- WCS beam deconvolved ellipse pars
		bool m_hasWCSDeconvolvedEllipsePars;
		double m_bmaj_deconv_wcs;
		double m_bmin_deconv_wcs;
		double m_pa_deconv_wcs;

		//- Source component flag
		int m_flag;

		//- Source component type
		int m_type;

		//- Source component selection flag
		bool m_selected;

	ClassDef(SourceComponentPars,7)

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
		* \brief Copy constructor
		*/
		SourceFitPars(const SourceFitPars& sourceFitPars){
			((SourceFitPars&)sourceFitPars).Copy(*this);
		}

		/**
		* \brief Class destructor: free allocated memory
		*/		
		virtual ~SourceFitPars() {}
		
		
		/**
		* \brief Assignment Operator
		*/
		SourceFitPars& operator=(const SourceFitPars& sourceFitPars){
			// Operator =
 		 	if (this != &sourceFitPars) ((SourceFitPars&)sourceFitPars).Copy(*this);
  		return *this;
		}

		/**
		* \brief Copy method
		*/
		void Copy(TObject& obj) const {
			// Copy this source to source obj	
  		((SourceFitPars&)obj).nComponents = nComponents;
			((SourceFitPars&)obj).chi2 = chi2;	
			((SourceFitPars&)obj).ndof = ndof;
			((SourceFitPars&)obj).npars = npars;	
			((SourceFitPars&)obj).npars_free = npars_free;	
			((SourceFitPars&)obj).npars_component= npars_component;
			((SourceFitPars&)obj).nfit_points = nfit_points;	
			((SourceFitPars&)obj).status = status;	
			((SourceFitPars&)obj).minimizer_status = minimizer_status;	
			((SourceFitPars&)obj).offset = offset;	
			((SourceFitPars&)obj).offset_err = offset_err;	
			((SourceFitPars&)obj).residualMean = residualMean;	
			((SourceFitPars&)obj).residualRMS = residualRMS;	
			((SourceFitPars&)obj).residualMedian = residualMedian;	
			((SourceFitPars&)obj).residualMAD = residualMAD;	
			((SourceFitPars&)obj).residualMin = residualMin;	
			((SourceFitPars&)obj).residualMax = residualMax;	

			((SourceFitPars&)obj).pars = pars;	
			((SourceFitPars&)obj).thetaFixed = thetaFixed;	
			((SourceFitPars&)obj).offsetFixed = offsetFixed;	
			((SourceFitPars&)obj).sigmaFixed = sigmaFixed;	
			((SourceFitPars&)obj).fluxDensity = fluxDensity;	
			((SourceFitPars&)obj).fluxDensityErr = fluxDensityErr;	

			((SourceFitPars&)obj).fitQuality = fitQuality;
		
			//Copy matrix
			int nRows= fitCovarianceMatrix.GetNrows();
			int nCols= fitCovarianceMatrix.GetNcols();
			((SourceFitPars&)obj).fitCovarianceMatrix.ResizeTo(nRows,nCols);
			for(int i=0;i<nRows;i++){
				for(int j=0;j<nCols;j++){
					double w= fitCovarianceMatrix(i,j);
					(((SourceFitPars&)obj).fitCovarianceMatrix)(i,j)= w;
				}
			}
		
			nRows= fluxDensityDerivMatrix.GetNrows();
			nCols= fluxDensityDerivMatrix.GetNcols();
			((SourceFitPars&)obj).fluxDensityDerivMatrix.ResizeTo(nRows,nCols);
			for(int i=0;i<nRows;i++){
				for(int j=0;j<nCols;j++){
					double w= fluxDensityDerivMatrix(i,j);
					(((SourceFitPars&)obj).fluxDensityDerivMatrix)(i,j)= w;
				}
			}
		
		}//close Copy()


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
				#ifdef LOGGING_ENABLED
					WARN_LOG("Invalid component id ("<<componentId<<" given, cannot find par to be set!");
				#endif
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
		* \brief Get number of selected components
		*/
		int GetNSelComponents()
		{
			int nSelComponents= 0;
			for(size_t i=0;i<pars.size();i++){
				bool isSelected= pars[i].IsSelected();
				if(isSelected) nSelComponents++;
			}
			return nSelComponents;
		}
	
		/**
		* \brief Check if component is selected
		*/
		bool IsSelectedComponent(int componentId){
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning false!");
				#endif
				return false;
			}
			return pars[componentId].IsSelected();
		}

		/**
		* \brief Set if component is selected
		*/
		int SetSelectedComponent(int componentId,bool selected)
		{
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning false!");
				#endif
				return -1;
			}
			pars[componentId].SetSelected(selected);
			return 0;
		}

		/**
		* \brief Remove fit components
		*/
		int RemoveComponents(std::vector<int> componentIds)
		{
			//Check components given
			int nComponentsToRemove= static_cast<int>(componentIds.size());
			if(nComponentsToRemove<=0 || nComponentsToRemove>nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("No components given to remove or too many (exceeding number of current components)!");
				#endif
				return -1;
			}
			for(size_t i=0;i<componentIds.size();i++){
				int componentId= componentIds[i];
				if(componentId<0 || componentId>=nComponents) {
					#ifdef LOGGING_ENABLED
						WARN_LOG("Component "<<componentId<<" does not exist!");
					#endif
					return -1;
				}
			}//end loop

			//Delete fit components from vector
			CodeUtils::DeleteItems(pars,componentIds);

			//Recompute nComponents
			nComponents= static_cast<int>(pars.size());

			return 0;
		}

		/**
		* \brief Get component flux density
		*/
		double GetComponentFluxDensity(int componentId){
			if(componentId<0 || componentId>=nComponents) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning zero flux!");
				#endif
				return 0;
			}
			return pars[componentId].GetFluxDensity();
		}

		/**
		* \brief Get component flux density error
		*/
		double GetComponentFluxDensityErr(int componentId){
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning zero error flux!");
				#endif
				return 0;
			}
			//Get component flux density deriv matrix
			TMatrixD D;
			if(GetComponentFluxDerivMatrix(D,componentId)<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to compute component flux derivative matrix, returning zero error flux!");
				#endif
				return 0;
			}

			TMatrixD D_t= TMatrixD(TMatrixD::kTransposed,D);
			TMatrixD VarMatrix= D*fitCovarianceMatrix*D_t;
			double Var= VarMatrix(0,0);
			if(Var<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Flux density variance for component "<<componentId<<" is negative (this should not occur, check for bugs or numerical roundoff errors!)");
				#endif
				return 0;
			}
			double Err= sqrt(Var);

			//Convert to Jy
			Err/= 1.e+3;

			return Err;			

		}//close GetComponentFluxDensityErr()

		/**
		* \brief Get component fit beam ellipse pars
		*/
		int GetComponentBeamEllipsePars(int componentId,double& bmaj,double& bmin,double& pa)
		{
			//Init values
			bmaj= 0;
			bmin= 0;
			pa= 0;
			if(componentId<0 || componentId>=nComponents){	
				#ifdef LOGGING_ENABLED	
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetBeamEllipsePars(bmaj,bmin,pa);
		}

		/**
		* \brief Get component beam ellipse eccentricity
		*/
		double GetComponentBeamEllipseEccentricity(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning E=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetBeamEllipseEccentricity();
		}

		/**
		* \brief Get component beam ellipse area
		*/
		double GetComponentBeamEllipseArea(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning A=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetBeamEllipseArea();
		}

		/**
		* \brief Has component beam pars stored?
		*/
		bool HasComponentBeamEllipsePars(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){	
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning false!");
				#endif
				return 0;
			}
			return pars[componentId].HasBeamEllipsePars();		
		}

		/**
		* \brief Get component fit ellipse pars
		*/
		int GetComponentFitEllipsePars(int componentId,double& x0,double& y0,double& bmaj,double& bmin,double& pa)
		{
			//Init values
			x0= 0;
			y0= 0;
			bmaj= 0;
			bmin= 0;
			pa= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetEllipsePars(x0,y0,bmaj,bmin,pa);
		}

		/**
		* \brief Get component fit ellipse eccentricity
		*/
		double GetComponentFitEllipseEccentricity(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning E=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetEllipseEccentricity();
		}

		/**
		* \brief Get component fit ellipse area
		*/
		double GetComponentFitEllipseArea(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist, returning A=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetEllipseArea();
		}

		/**
		* \brief Get component fit ellipse rot angle vs beam
		*/
		double GetComponentFitEllipseRotAngleVSBeam(int componentId)
		{
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED				
					WARN_LOG("Component "<<componentId<<" does not exist, returning A=0!");
				#endif
				return 0;
			}
			return pars[componentId].GetEllipseRotAngleVSBeam();
		}

		/**
		* \brief Get component fit ellipse pars
		*/
		int GetComponentFitEllipseParErrors(int componentId,double& x0_err,double& y0_err,double& bmaj_err,double& bmin_err,double& pa_err)
		{
			//Init values
			x0_err= 0;
			y0_err= 0;
			bmaj_err= 0;
			bmin_err= 0;
			pa_err= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetEllipseParErrors(x0_err,y0_err,bmaj_err,bmin_err,pa_err);
		}

		/**
		* \brief Get component WCS fit ellipse pars
		*/
		int GetComponentFitWCSEllipsePars(int componentId,double& x0_wcs,double& y0_wcs,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
		{
			//Init values
			x0_wcs= 0;
			y0_wcs= 0;
			bmaj_wcs= 0;
			bmin_wcs= 0;
			pa_wcs= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetWCSEllipsePars(x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs);
		}

		/**
		* \brief Get component WCS fit ellipse par errors
		*/
		int GetComponentFitWCSEllipseParErrors(int componentId,double& x0_wcs_err,double& y0_wcs_err,double& bmaj_wcs_err,double& bmin_wcs_err,double& pa_wcs_err)
		{
			//Init values
			x0_wcs_err= 0;
			y0_wcs_err= 0;
			bmaj_wcs_err= 0;
			bmin_wcs_err= 0;
			pa_wcs_err= 0;
			if(componentId<0 || componentId>=nComponents){	
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetWCSEllipseParErrors(x0_wcs_err,y0_wcs_err,bmaj_wcs_err,bmin_wcs_err,pa_wcs_err);
		}

		/**
		* \brief Get component WCS fit ellipse pars
		*/
		int GetComponentFitWCSDeconvolvedEllipsePars(int componentId,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
		{
			//Init values
			bmaj_wcs= 0;
			bmin_wcs= 0;
			pa_wcs= 0;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}

			//Get component fit ellipse pars
			return pars[componentId].GetWCSDeconvolvedEllipsePars(bmaj_wcs,bmin_wcs,pa_wcs);
		}

		/**
		* \brief Get flux density
		*/
		double GetFluxDensity(){return fluxDensity;}

		/**
		* \brief Get flux density error
		*/
		double GetFluxDensityErr(){return fluxDensityErr;}

		/**
		* \brief Set flux density
		*/
		void SetFluxDensity(double value){fluxDensity=value;}
		/**
		* \brief Set flux density error
		*/
		void SetFluxDensityErr(double value){fluxDensityErr=value;}

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
		* \brief Set number of component fit parameters 
		*/
		void SetNComponentPars(int nc){npars_component=nc;}
		/**
		* \brief Get number of component fit parameters 
		*/
		int GetNComponentPars(){return npars_component;}


		
		/**
		* \brief Set total number of pars 
		*/
		void SetNPars(int value){npars=value;}
		/**
		* \brief Get number of free parameters 
		*/
		int GetNPars(){return npars;}

		/**
		* \brief Set number of fitted data 
		*/
		void SetNFitPoints(int value){nfit_points=value;}

		/**
		* \brief Get number of fitted data 
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
		* \brief Set fit covariance matrix
		*/
		int SetCovarianceMatrix(double* errMatrixValues,int ndim)
		{
			if(!errMatrixValues || ndim<=0) return -1;
			fitCovarianceMatrix.ResizeTo(ndim,ndim);
			fitCovarianceMatrix.Zero();
			try{
				for(int i=0;i<ndim;i++){
					for(int j=0;j<ndim;j++){
						int index= i*ndim+j;
						fitCovarianceMatrix(i,j)= errMatrixValues[index];
					}//end loop dim
				}//end loop dim
			}//close try block
			catch(...){	
				#ifdef LOGGING_ENABLED
					ERROR_LOG("C++ exception occurred while filling fit covariance matrix (hint: array size and given dim are different?");
				#endif
				fitCovarianceMatrix.ResizeTo(0,0);
				return 0;
			}
			return 0;
		}//close SetCovarianceMatrix()

		/**
		* \brief Remove component(s) from covariance matrix
		*/
		/*
		int RemoveComponentsFromCovarianceMatrix(std::vector<int> componentIds)
		{	
			//Check input
			if(componentIds.empty()) return 0;
			int nComponents= (fitCovarianceMatrix.GetNrows()-1)/npars_component;
			
			for(size_t i=0;i<componentIds.size();i++){
				if(componentIds[i]<0 || componentIds[i]>=nComponents){
					#ifdef LOGGING_ENABLED
						ERROR_LOG("Invalid component id ("<<componentIds[i]<<") to be removed (#"<<nComponents<<" present in cov matrix)!");
					#endif
					return -1;
				}
			}
		
			//Set components to be kept
			std::vector<int> keptComponentIds;
			for(int i=0;i<nComponents;i++){
				bool keepComponent= true;
				for(size_t k=0;k<componentIds.size();k++){
					if(i==componentIds[k]){
						keepComponent= false;
						break;
					}
				}
				if(keepComponent) keptComponentIds.push_back(i);
			}

			//Sort component ids to be kept
			std::sort(keptComponentIds.begin(),keptComponentIds.end());

			//Check if there are components left
			int nComponents_left= static_cast<int>(keptComponentIds.size());	
			int nDim_new= nComponents_left*npars_component + 1;
			if(nComponents_left<=0){
				fitCovarianceMatrix.ResizeTo(0,0);
				return 0;
			}

			//Copy cov matrix
			TMatrixD C= fitCovarianceMatrix;
			fitCovarianceMatrix.ResizeTo(nDim_new,nDim_new);
			for(int i=0;i<fitCovarianceMatrix.GetNrows();i++){
				
				for(int j=0;j<fitCovarianceMatrix.GetNcols();j++){
					int componentIndex= j/npars_component;
					int parIndex= j%npars_component;
					int componentId= keptComponentIds[componentIndex];
					int colId= componentId*npars_component + parIndex;
					int rowId= componentId*npars_component + j;
					//fitCovarianceMatrix(i,j)= C(colId,colId);
				}
			}

			for(int k=0;k<nComponents_left;k++){
				int componentId= keptComponentIds[k];
				int row_start= componentId*npars_component;
				for(int i=0;i<npars_component;i++){
					for(int j=0;j<npars_component;j++){
						int index= i*npars_component + j;
						
						//fitCovarianceMatrix(i,j)= errMatrixValues[index];
					}
				}

				for(int i=0;i<fitCovarianceMatrix.GetNrows();i++){
					for(int j=0;j<fitCovarianceMatrix.GetNcols();j++){
						int componentId= 
						int rowId= i*npars_component + j;
						
						//fitCovarianceMatrix(i,j)= errMatrixValues[index];
					}
				}
			}

			return 0;

		}//close RemoveComponentsFromCovarianceMatrix()
		*/
		

		/**
		* \brief Get fit covariance matrix
		*/	
		TMatrixD& GetCovarianceMatrix(){return fitCovarianceMatrix;}

		/**
		* \brief Print fit covariance matrix
		*/
		void PrintCovarianceMatrix(){
			fitCovarianceMatrix.Print();
		}

		/**
		* \brief Get flux density derivative matrix
		*/	
		TMatrixD& GetFluxDensityDerivMatrix(){return fluxDensityDerivMatrix;}

		/**
		* \brief Print flux density derivative matrix
		*/
		void PrintFluxDensityDerivMatrix(){
			fluxDensityDerivMatrix.Print();
		}

		/**
		* \brief Compute flux density derivative matrix
		*/
		int ComputeFluxDensityDerivMatrix(){
			//Check size
			//if(npars_free<=0 || pars.empty()) {
			if(npars<=0 || pars.empty()) {
				#ifdef LOGGING_ENABLED
					WARN_LOG("Cannot compute derivative matrix as no fit pars are stored and/or number of free pars is not initialized!");
				#endif
				return -1;
			}

			//Init matrix to 1 x Nfree_pars
			//fluxDensityDerivMatrix.ResizeTo(1,npars_free);
			fluxDensityDerivMatrix.ResizeTo(1,npars);
			fluxDensityDerivMatrix.Zero();

			//Fill matrix 
			//NB: consider the cases in which sigma/offset/theta are fixed
			int parCounter= 0;
			for(int i=0;i<nComponents;i++){
				//Retrieve fitted pars for this component
				double A= pars[i].GetParValue("A");
				A*= 1.e+3;//convert to mJy
				double sigmaX= pars[i].GetParValue("sigmaX");
				double sigmaY= pars[i].GetParValue("sigmaY");

				//Fill derivative wrt to amplitude
				//deriv= 2 pi sigmaX sigmaY
				double deriv_ampl= 2.*TMath::Pi()*sigmaX*sigmaY;
				fluxDensityDerivMatrix(0,parCounter)= deriv_ampl;
				parCounter++;

				//Fill derivative wrt to centroids
				//deriv= 0 (for both x & y)
				fluxDensityDerivMatrix(0,parCounter)= 0;
				fluxDensityDerivMatrix(0,parCounter+1)= 0;
				parCounter+= 2;

				//Fill derivative wrt to sigmas (if not fixed)
				//deriv_x= 2 pi A sigma_y
				//deriv_y= 2 pi A sigma_x
				if(!sigmaFixed){
					double deriv_sigmaX= 2.*TMath::Pi()*A*sigmaY;
					double deriv_sigmaY= 2.*TMath::Pi()*A*sigmaX;
					fluxDensityDerivMatrix(0,parCounter)= deriv_sigmaX;
					fluxDensityDerivMatrix(0,parCounter+1)= deriv_sigmaY;
					parCounter+= 2;
				}

				//Fill derivative wrt to theta (if not fixed)
				//deriv_theta= 0;
				if(!thetaFixed){
					fluxDensityDerivMatrix(0,parCounter)= 0;
					parCounter++;
				}

			}//end loop components	

			//Fill deriv wrt to offset (if not fixed)
			//deriv_offset= 0;
			if(!offsetFixed){
				fluxDensityDerivMatrix(0,parCounter)= 0;
			}
	
			return 0;
		}//close ComputeFluxDensityDerivMatrix()
		
		/**
		* \brief Compute flux density
		*/
		void ComputeFluxDensity(){
			//Summing up flux for all components
			fluxDensity= 0;
			for(int i=0;i<nComponents;i++){
				double componentFluxDensity= pars[i].GetFluxDensity();
				fluxDensity+= componentFluxDensity;
			}//end loop components
		}//close ComputeFluxDensity()

		/**
		* \brief Compute flux density error
		*/
		int ComputeFluxDensityError(){
			//Do some checks on covariance & deriv matrix
			fluxDensityErr= 0;
			int nRows= fitCovarianceMatrix.GetNrows();
			int nCols= fluxDensityDerivMatrix.GetNcols();
			if(nCols<=0 || nRows<=0 || nRows!=nCols){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Fit covariance and/or deriv matrix were not computed or have invalid dimensions!");
				#endif
				return -1;
			}
	
			//Compute fluxDensityVariance= D x CovMatrix x D^t  (D=deriv matrix)
			TMatrixD fluxDensityDerivMatrix_t= TMatrixD(TMatrixD::kTransposed,fluxDensityDerivMatrix);
			TMatrixD fluxDensityVarianceMatrix= fluxDensityDerivMatrix*fitCovarianceMatrix*fluxDensityDerivMatrix_t;
			double fluxDensityVariance= fluxDensityVarianceMatrix(0,0);
			if(fluxDensityVariance<0){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Flux density variance is negative (this should not occur, check for bugs or numerical roundoff errors!)");
				#endif
				return -1;
			}
			fluxDensityErr= sqrt(fluxDensityVariance);

			//Convert back to Jy
			fluxDensityErr/= 1.e+3;
			
			return 0;
		}//close ComputeFluxDensityError()


		/**
		* \brief Compute component ellipse pars
		*/
		int ComputeComponentEllipsePars()
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeEllipsePars()<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentEllipsePars()


		/**
		* \brief Compute component WCS ellipse pars
		*/
		int ComputeComponentWCSEllipsePars(WCS* wcs)
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeWCSEllipsePars(wcs)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute WCS ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentWCSEllipsePars()

		/**
		* \brief Compute component WCS ellipse pars
		*/
		int ComputeComponentWCSEllipseParsSimple(WCS* wcs)
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeWCSEllipseParsSimple(wcs)<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute WCS ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentWCSEllipseParsSimple()

		
		/**
		* \brief Set component beam ellipse pars
		*/
		void SetComponentBeamEllipsePars(double bmaj,double bmin,double pa)
		{
			for(size_t i=0;i<pars.size();i++){
				pars[i].SetBeamEllipsePars(bmaj,bmin,pa);
			}
		}


		/**
		* \brief Compute component WCS ellipse pars
		*/
		int ComputeComponentWCSDeconvolvedEllipsePars()
		{
			int status= 0;
			for(size_t i=0;i<pars.size();i++){
				if(pars[i].ComputeWCSDeconvolvedEllipsePars()<0){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to compute WCS ellipse pars for fit component no. "<<i+1);
					#endif
					status= -1;
				}
			}
			return status;

		}//close ComputeComponentWCSDeconvolvedEllipsePars()

		/**
		* \brief Print
		*/
		void Print(){
			cout<<"*** FIT RESULTS ***"<<endl;
			cout<<"nPars="<<npars<<", nParsFree="<<npars_free<<", fitStatus="<<status<<" (minimizer status="<<minimizer_status<<"), fitQuality="<<fitQuality<<endl;
			cout<<"Chi2="<<chi2<<", ndf="<<ndof<<", Chi2/NDF="<<chi2/ndof<<endl;
			cout<<"fluxDensity="<<fluxDensity<<" +- "<<fluxDensityErr<<endl;
			for(int i=0;i<nComponents;i++){
				cout<<"--> Component "<<i+1<<endl;
				cout<<"fluxDensity="<<GetComponentFluxDensity(i)<<" +- "<<GetComponentFluxDensityErr(i)<<endl;
				cout<<"A="<<pars[i].GetParValue("A")<<" +- "<<pars[i].GetParError("A")<<endl;
				cout<<"(x0,y0)=("<<	pars[i].GetParValue("x0")<<","<<pars[i].GetParValue("y0")<<") err("<<pars[i].GetParError("x0")<<","<<pars[i].GetParError("y0")<<")"<<endl;
				cout<<"(sigmaX,sigmaY)=("<<	pars[i].GetParValue("sigmaX")<<","<<pars[i].GetParValue("sigmaY")<<") err("<<pars[i].GetParError("sigmaX")<<","<<pars[i].GetParError("sigmaY")<<")"<<endl;
				cout<<"theta="<<pars[i].GetParValue("theta")<<" +- "<<pars[i].GetParError("theta")<<endl;
			}
			cout<<"****************"<<endl;
		}

		/**
		* \brief Set theta fixed
		*/
		void SetThetaFixed(bool choice){thetaFixed=choice;}
		/**
		* \brief Is theta fixed?
		*/
		bool IsThetaFixed(){return thetaFixed;}
		/**
		* \brief Set theta fixed
		*/
		void SetSigmaFixed(bool choice){sigmaFixed=choice;}
		/**
		* \brief Is sigma fixed?
		*/
		bool IsSigmaFixed(){return sigmaFixed;}
		/**
		* \brief Set offset fixed
		*/
		void SetOffsetFixed(bool choice){offsetFixed=choice;}
		/**
		* \brief Is offset fixed?
		*/
		bool IsOffsetFixed(){return offsetFixed;}

		

		/**
		* \brief Get number of free parameters per component
		*/
		int GetFreeParsPerComponent(){
			int nFreeParsPerComponent= 3;//Ampl + centroids
			if(!sigmaFixed) nFreeParsPerComponent+= 2;//sigmas
			if(!thetaFixed) nFreeParsPerComponent++;
			return nFreeParsPerComponent;
		}

		/**
		* \brief Get component flux derivative matrix
		*/
		int GetComponentFluxDerivMatrix(TMatrixD& D,int componentId){
			//Keep only selected component
			//int nFreeParsPerComponent= GetFreeParsPerComponent();
			//int start_index= componentId*nFreeParsPerComponent;
			//int last_index= start_index + nFreeParsPerComponent-1;
			int start_index= componentId*npars_component;
			int last_index= start_index + npars_component-1;
			int nCols= fluxDensityDerivMatrix.GetNcols(); 
			if(nCols<=last_index){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Trying to access to an not-existing element (index="<<last_index<<") of derivative matrix (dim="<<nCols<<") (hint: derivative matrix not initialized)!");	
				#endif
				return -1;
			}
			
			//Init to flux derivative matrix
			D.ResizeTo(1,nCols);
			D.Zero();
			for(int i=start_index;i<=last_index;i++){
				D(0,i)= fluxDensityDerivMatrix(0,i);
			}

			//cout<<"*** DERIV MATRIX COMPONENT "<<componentId<<" ***"<<endl;
			//D.Print();

			return 0;

		}//close GetComponentFluxDerivMatrix()

		/**
		* \brief Set component flag
		*/
		int SetComponentFlag(int componentId,int flag)
		{		
			//Check component id	
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return 0;
			}

			//Get component fit ellipse pars
			pars[componentId].SetFlag(flag);

			return 0;
		}
	
		/**
		* \brief Get component flag
		*/
		int GetComponentFlag(int& flag,int componentId)
		{
			//Check component id
			flag= -1;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}
	
			//Retrieve flag
			flag= pars[componentId].GetFlag();

			return 0;
		}

		/**
		* \brief Set component flag
		*/
		int SetComponentType(int componentId,int type)
		{		
			//Check component id	
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return 0;
			}

			//Get component fit ellipse pars
			pars[componentId].SetType(type);

			return 0;
		}

		/**
		* \brief Get component type
		*/
		int GetComponentType(int& type,int componentId)
		{
			//Check component id
			type= -1;
			if(componentId<0 || componentId>=nComponents){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Component "<<componentId<<" does not exist!");
				#endif
				return -1;
			}
	
			//Retrieve flag
			type= pars[componentId].GetType();

			return 0;
		}

		/**
		* \brief Set component image pix size
		*/
		void SetComponentImagePixSize(double value)
		{		
			for(size_t i=0;i<pars.size();i++){
				pars[i].SetImagePixSize(value);
			}
		}

		/**
		* \brief Get is good fit flag
		*/
		int GetFitQuality(){return fitQuality;}

		/**
		* \brief Set is good fit flag
		*/
		void SetFitQuality(int flag){fitQuality=flag;}

		

	private:

		/**
		* \brief Initialize component pars
		*/
		void Init(){
			nComponents= 0;
			chi2= 0;
			ndof= 0;
			npars= 0;
			npars_free= 0;	
			npars_component= 6;
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
			thetaFixed= false;
			offsetFixed= false;
			sigmaFixed= false;	
			fluxDensity= 0;
			fluxDensityErr= 0;
			fitQuality= eBadFit;
		}

		
	private:

		int nComponents;
		double chi2;
		double ndof;
		int npars;
		int npars_free;
		int npars_component;
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
		TMatrixD fitCovarianceMatrix;
		TMatrixD fluxDensityDerivMatrix;
		bool thetaFixed;
		bool offsetFixed;
		bool sigmaFixed;
		double fluxDensity;
		double fluxDensityErr;

		int fitQuality;

	ClassDef(SourceFitPars,5)

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
		/*
		enum FitStatusFlag {
			eFitUnknownStatus= 0,
			eFitAborted= 1,
			eFitNotConverged= 2,
			eFitConverged= 3,
			eFitConvergedWithWarns= 4
		};
		*/

		/**
		* \brief Source fit data
		*/
		struct SourceFitData {
			SourceFitData(){
				x= 0;
				y= 0;
				S= 0;
				Serr= 0;	
			}

			SourceFitData(double _x,double _y,double _S,double _Serr)
				: x(_x),y(_y),S(_S),Serr(_Serr)
			{}
			double x;
			double y;
			double S;
			double Serr;
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
		//static double Gaus2DMixtureFcn(double* x, double* p);
		double Gaus2DMixtureFcn(double* x, const double* p);

		/**
		* \brief 2D Gaussian model used for the fit
		*/
		//static double Gaus2DFcn(double* x, double* p);
		double Gaus2DFcn(double* x, const double* p);

	protected:
		/**
		* \brief Check if fit has parameters converged at bounds
		*/
		bool HasFitParsAtLimit(const ROOT::Fit::FitResult& fitRes);
		/**
		* \brief Returns list of parameters at limits (version with Minimizer)
		*/
		int GetParsAtLimits(std::vector<int>& parsAtLimits,ROOT::Math::Minimizer* fitter);
		
		/**
		* \brief Returns list of parameters at limits (version with TFitter)
		*/
		//int GetParsAtLimits(std::vector<int>& parsAtLimits,TFitter* fitter);
		
		/**
		* \brief Estimate fit components
		*/
		int EstimateFitComponents(std::vector<std::vector<double>>& fitPars_start,Source* aSource,SourceFitOptions& fitOptions);
		
		/**
		* \brief Perform fit using given initial fit components pars
		*/
		int DoChi2Fit(Source* aSource,SourceFitOptions& fitOptions,std::vector< std::vector<double> >& fitPars_start);

		/**
		* \brief Init data
		*/
		int InitData(Source* aSource,SourceFitOptions& fitOptions);
		/**
		* \brief Check fit options
		*/
		int CheckFitOptions(SourceFitOptions& fitOptions);

		/**
		* \brief Chi2 fit function 
		*/
		double Chi2Fcn(const double* par);
			
		/**
		* \brief Chi2 fit function (old MINUIT version)
		*/
		//static void Chi2Fcn(int& nPar,double* grad,double& fval,double* p,int iflag);
		

		/**
		* \brief Initialize fitter
		*/
		ROOT::Math::Minimizer* InitMinimizer(int nFitPars,SourceFitOptions& fitOptions);
	
	private:
	
		//static int m_NFitComponents;
		//static int m_fitStatus;
		//static SourceFitPars m_sourceFitPars;
		//static TH2D* m_fluxMapHisto;
		//static std::vector<SourceFitData> m_fitData;
		//static std::vector<SourceFitData> m_fitHaloData;
		//static double m_chi2RegPar;
		//static double m_bkgMean;
		//static double m_rmsMean;
		//static double m_sourceX0;
		//static double m_sourceY0;

		int m_NFitComponents;
		int m_fitStatus;
		SourceFitPars m_sourceFitPars;
		
		//TH2D* m_fluxMapHisto;
		std::vector<SourceFitData> m_fitData;
		std::vector<SourceFitData> m_fitHaloData;
		double m_chi2RegPar;
		double m_bkgMean;
		double m_rmsMean;
		double m_sourceX0;
		double m_sourceY0;
		
	ClassDef(SourceFitter,1)

};//close SourceFitter class

#ifdef __MAKECINT__
#pragma link C++ class SourceComponentPars+;
#pragma link C++ class SourceFitPars+;
#pragma link C++ class SourceFitter+;
#pragma link C++ struct SourceFitter::SourceFitData+;
#pragma link C++ enum FitStatusFlag;
#pragma link C++ struct SourceFitOptions+;
#endif

}//close namespace

#endif
