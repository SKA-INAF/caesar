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
* @file SourceComponentPars.cc
* @class SourceComponentPars
* @brief SourceComponentPars
*
* Source component parameters class
* @author S. Riggi
* @date 01/09/2017
*/

#include <SourceComponentPars.h>

#include <Image.h>
#include <Source.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <SysUtils.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif
#include <Consts.h>
#include <WCSUtils.h>

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

namespace Caesar {

//Constructor
SourceComponentPars::SourceComponentPars()
	: TObject()
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

}//close costructor

//Destructor
SourceComponentPars::~SourceComponentPars()
{
	FitPars.clear();
	FitParsErr.clear();

}//close destructor

// Copy constructor
SourceComponentPars::SourceComponentPars(const SourceComponentPars& pars) 
{
  Init();
  ((SourceComponentPars&)pars).Copy(*this);
}

//Source component pars
SourceComponentPars& SourceComponentPars::operator=(const SourceComponentPars& pars)
{
	// Operator =
  if (this != &pars) ((SourceComponentPars&)pars).Copy(*this);
  return *this;
}

//Copy method
void SourceComponentPars::Copy(TObject &obj) const 
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


void SourceComponentPars::Init()
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

}//close Init()
		

void SourceComponentPars::InitEllipsePars()
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


int SourceComponentPars::SetParValueAndError(std::string parName,double parVal,double parErr)
{			
	if(!CodeUtils::HasMapKey(FitPars,parName)) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid par name ("<<parName<<" given, cannot find par to be set!");
		#endif
		return -1;
	}
	FitPars[parName]= parVal;
	FitParsErr[parName]= parErr;
			
	return 0;
		
}//close SetParValueAndError()


double SourceComponentPars::GetParValue(std::string parName)
{
	if(!CodeUtils::HasMapKey(FitPars,parName)) return -999;
	return FitPars[parName];
}

double SourceComponentPars::GetParError(std::string parName)
{
	if(!CodeUtils::HasMapKey(FitParsErr,parName)) return -999;
	return FitParsErr[parName];
}


TEllipse* SourceComponentPars::GetFitEllipse(bool useFWHM)
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



int SourceComponentPars::ComputeEllipsePars()
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


int SourceComponentPars::GetEllipsePars(double& x0,double& y0,double& bmaj,double& bmin,double& pa)
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


int SourceComponentPars::GetEllipseParErrors(double& x0_err,double& y0_err,double& bmaj_err,double& bmin_err,double& pa_err)
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


int SourceComponentPars::ComputeWCSEllipsePars(WCS* wcs)
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


int SourceComponentPars::ComputeWCSEllipseParsSimple(WCS* wcs)
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


int SourceComponentPars::GetWCSEllipsePars(double& x0_wcs,double& y0_wcs,double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
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

int SourceComponentPars::GetWCSEllipseParErrors(double& x0_err_wcs,double& y0_err_wcs,double& bmaj_err_wcs,double& bmin_err_wcs,double& pa_err_wcs)
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


void SourceComponentPars::SetBeamEllipsePars(double bmaj,double bmin,double pa)
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
	
}//close SetBeamEllipsePars()


int SourceComponentPars::GetBeamEllipsePars(double& bmaj,double& bmin,double& pa)
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
		
}//close GetBeamEllipsePars()


int SourceComponentPars::ComputeWCSDeconvolvedEllipsePars()
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
			DEBUG_LOG("Numerical error (beta^2 is <0 in formula)!");
		#endif
		return -1;
	}
	double beta= sqrt(beta2);

	//Check if fit ellipse is smaller than beam.
	//If so, do not attempt too deconvolve. Set deconv ellipse to fitted ellipse.
	if(sum2<=sum2_beam){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Fitted ellipse beam is smaller than beam, do not deconvolve (set deconvolved ellipse to fitted ellipse)!");
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
			DEBUG_LOG("Numerical error (bmaj^2/bmin^2 deconvolved are <0 in formula)!");	
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


int SourceComponentPars::GetWCSDeconvolvedEllipsePars(double& bmaj_wcs,double& bmin_wcs,double& pa_wcs)
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

}//close namespace
