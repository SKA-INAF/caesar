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
* @file AstroObject.cc
* @class AstroObject
* @brief AstroObject class
*
* AstroObject class
* @author S. Riggi
* @date 06/08/2019
*/

#include <AstroObject.h>
#include <MathUtils.h>
#include <CodeUtils.h>
#include <AstroUtils.h>
#include <Contour.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//ROOT headers
#include <TMath.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>
#include <deque>

using namespace std;


ClassImp(Caesar::AstroObject)

namespace Caesar {

//===========================
//==    AstroObject CLASS
//===========================
AstroObject::AstroObject()
{
	//Initialize fields
	Init();

}//close constructor


AstroObject::~AstroObject()
{

}//close destructor


void AstroObject::Init()
{
	index= -1;
	name= "";
	id_str= "";
	id= eUNKNOWN_OBJECT;
	subid= eUNKNOWN_OBJECT;
	x= 0;
	y= 0;
	xerr= 0;
	yerr= 0;
	refs= "";
	confirmed= false;

	hasFluxInfo= false;
	peakFlux= 0;
	peakFluxErr= 0;
	fluxDensity= 0;
	fluxDensityErr= 0;
	flux= 0;
	fluxErr= 0;

	hasSizeInfo= false;
	radius= 0;

	hasEllipseInfo= false;
	bmaj= 0;
	bmin= 0;
	pa= 0;
		
	hasDeconvEllipseInfo= false;
	bmaj_deconv= 0;
	bmin_deconv= 0;
	pa_deconv= 0;

	hasFrequencyInfo= false;
	nu= 0;
	dnu= 0;

	hasSpectralIndexInfo= false;
	isMultiSourceMatchIndex= false;
	isSpectralIndexFit= false;
	spectralIndex= -999;
	spectralIndexErr= -999;

}//close Init()


TEllipse* AstroObject::GetFitEllipse()
{
	//Check if ellipse info are available
	if(!hasEllipseInfo){
		return nullptr;
	}

	//Compute ellipse pars
	double X0= x;
	double Y0= y;
	double R1= bmaj/2;
	double R2= bmin/2;
	double theta= pa+90;

	//Create ellipse
	TEllipse* ellipse= new TEllipse(X0,Y0,R1,R2,0.,360.,theta);
	ellipse->SetLineWidth(2);
	ellipse->SetFillColor(0);
	ellipse->SetFillStyle(0);

	return ellipse;

}//close GetFitEllipse()


Contour* AstroObject::GetContour(bool computePars)
{
	//Check if ellipse info are available
	if(!hasEllipseInfo){
		return nullptr;
	}
	
	//Compute ellipse pars
	double X0= x;
	double Y0= y;
	double a= bmaj/2;
	double b= bmin/2;
	double theta= pa+90;

	//Create contour
	Contour* contour= new Contour;

	//Fill contour
	double t= 0;
	double t_step= 2;
	double t_max= 360;
	double theta_rad= theta*TMath::DegToRad();

	while(true){
		if(t>=t_max) break;	
		double t_rad= t*TMath::DegToRad();
		double x= X0 + a*cos(t_rad)*cos(theta_rad) - b*sin(t_rad)*sin(theta_rad);
		double y= Y0 + a*cos(t_rad)*sin(theta_rad) + b*sin(t_rad)*cos(theta_rad);
		contour->AddPoint(TVector2(x,y));
		t+= t_step;		
	}//end loop 
	
	//Compute contour parameters
	if(computePars && contour->ComputeParameters()<0){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/more failures occurred while computing contour parameters!");
		#endif
	}

	return contour;

}//close GetContour()

std::string AstroObject::GetDS9Region(std::string text,std::string color,std::vector<std::string> tags)
{
	std::string regionText= "";
	std::stringstream ss;
	std::string regionName= text;

	//If no region text is given use region name
	if(regionName=="") regionName= name;

	//If no external tags are given use id
	std::vector<std::string> regionTags= tags;
	if(regionTags.empty()){
		std::stringstream sstream;
		sstream<<"id"<<id; 
		regionTags.push_back(sstream.str());
	}

	//Set region shape
	if(hasEllipseInfo && bmaj>0 && bmin>0) {
		//regionText= AstroUtils::EllipseToDS9Region(ellipse,regionName,color,regionTags,useImageCoords);

		double R1= bmaj/2;
		double R2= bmin/2;
		double theta= pa+90;
		ss<<"ellipse("<<x<<","<<y<<","<<R1<<"\","<<R2<<"\","<<theta<<") ";

	}
	else{
		//No ellipse pars available, check if radius is available	
		double r= 1;//dummy radius	
		if(hasSizeInfo && radius>0){
			ss<<"circle("<<x<<","<<y<<","<<radius<<"\""<<") ";
		}
		else{
			//No radius available use dummy radius
			ss<<"circle("<<x<<","<<y<<","<<r<<"\""<<") ";
		}
		
	}//close else

	ss<<"# ";
	ss<<"text={"<<regionName<<"} ";
	ss<<"color="<<color<<" ";
	for(size_t k=0;k<regionTags.size();k++){
		ss<<"tag={"<<regionTags[k]<<"} ";
	}


	regionText= ss.str();

	return regionText;

}//close GetDS9Region()


}//close namespace
