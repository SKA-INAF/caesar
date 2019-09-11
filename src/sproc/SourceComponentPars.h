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
* @file SourceComponentPars.h
* @class SourceComponentPars
* @brief SourceComponentPars
*
* Source component parameters class
* @author S. Riggi
* @date 01/09/2017
*/

#ifndef _SOURCE_COMPONENT_PARS_h
#define _SOURCE_COMPONENT_PARS_h 1

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


namespace Caesar {


//===========================================
//==         SOURCE COMPONENT PARS
//===========================================
class SourceComponentPars : public TObject {

	public: 
		/** 
		\brief Class constructor: initialize structures.
 		*/
		SourceComponentPars();

		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~SourceComponentPars();

		/**
		* \brief Copy constructor
		*/
		SourceComponentPars(const SourceComponentPars& pars);
		
		/**
		* \brief Assignment Operator
		*/
		SourceComponentPars& operator=(const SourceComponentPars& pars);
		
		/**
		* \brief Copy method
		*/
		void Copy(TObject& obj) const;


	public:
		/** 
		\brief Set par value & error
 		*/
		int SetParValueAndError(std::string parName,double parVal,double parErr);
			

		/** 
		\brief Get par value
 		*/
		double GetParValue(std::string parName);

		/** 
		\brief Get par error
 		*/
		double GetParError(std::string parName);

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
		TEllipse* GetFitEllipse(bool useFWHM=true);

		/** 
		\brief Has ellipse pars?
 		*/
		bool HasEllipsePars(){return m_hasEllipsePars;}

		/** 
		\brief Compute ellipse pars in pixel coordinates
 		*/
		int ComputeEllipsePars();


		/** 
		\brief Get fit ellipse pars.
 		*/
		int GetEllipsePars(double& x0,double& y0,double& bmaj,double& bmin,double& pa);

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
		int GetEllipseParErrors(double& x0_err,double& y0_err,double& bmaj_err,double& bmin_err,double& pa_err);

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
		int ComputeWCSEllipsePars(WCS* wcs);

		/**
		* \brief Compute ellipse pars in WCS coordinates, neglecting sky projection and assuming Euclidean distances
		*/
		int ComputeWCSEllipseParsSimple(WCS* wcs);


		/** 
		\brief Has WCS ellipse pars?
 		*/
		bool HasWCSEllipsePars(){return m_hasWCSEllipsePars;}

		/**
		* \brief Get ellipse pars in WCS coordinates
		*/
		int GetWCSEllipsePars(double& x0_wcs,double& y0_wcs,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs);

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
		int GetWCSEllipseParErrors(double& x0_err_wcs,double& y0_err_wcs,double& bmaj_err_wcs,double& bmin_err_wcs,double& pa_err_wcs);
		
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
		void SetBeamEllipsePars(double bmaj,double bmin,double pa);
	
		/**
		* \brief Has beam pars?
		*/
		bool HasBeamEllipsePars(){return m_hasBeamPars;}

		/**
		* \brief Get beam ellipse parameters
		*/
		int GetBeamEllipsePars(double& bmaj,double& bmin,double& pa);

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
		int ComputeWCSDeconvolvedEllipsePars();

		/**
		* \brief Has beam-deconvolved WCS ellipse pars?
		*/
		bool HasWCSDeconvolvedEllipsePars(){return m_hasWCSDeconvolvedEllipsePars;}

		/**
		* \brief Get ellipse pars in WCS coordinates
		*/
		int GetWCSDeconvolvedEllipsePars(double& bmaj_wcs,double& bmin_wcs,double& pa_wcs);
 		/**
		* \brief Set WCS deconvolved ellipse pars
		*/
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

		/**
		* \brief Get position centroid
		*/
		void GetPosition(double& xpos,double& ypos)
		{
			xpos= FitPars["x0"];
			ypos= FitPars["y0"];
		}

		/**
		* \brief Get peak flux
		*/
		double GetPeakFlux(){return FitPars["A"];}

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
		void InitEllipsePars();

		/**
		* \brief Init class data
		*/
		void Init();
		
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


#ifdef __MAKECINT__
#pragma link C++ class SourceComponentPars+;
#endif

}//close namespace

#endif
