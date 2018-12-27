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
* @file AstroUtils.cc
* @class AstroUtils
* @brief Utility functions for astronomical tasks
*
* Utility functions for astronomical tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <AstroUtils.h>
#include <Image.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <MathUtils.h>
#include <WCSUtils.h>

#include <TObject.h>
#include <TEllipse.h>

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

using namespace std;

ClassImp(Caesar::AstroUtils)

namespace Caesar {

AstroUtils::AstroUtils()
{

}

AstroUtils::~AstroUtils()
{

}

int AstroUtils::GetIAUCoords(std::string& iau,const std::string& s)
{
	//Init string 
	iau= "";

	//Split string on whitespaces
	std::vector<std::string> fields= CodeUtils::SplitStringOnWhitespaces(s);
	if(fields.size()!=3){
		WARN_LOG("Invalid number of parsed fields (n="<<fields.size()<<") when 3 expected!");
		return -1;
	}

	std::string ra_field= fields[0];
	std::string dec_field= fields[1];
	std::string cs_field= fields[2];

	//Last field gives the coord system flag
	std::string cs= "";
	if(cs_field=="FK5" || cs_field=="fk5") cs= "J";
	else if(cs_field=="FK4" || cs_field=="fk4") cs= "B";
	else if(cs_field=="galactic") cs= "G";
	else {
		WARN_LOG("");
		return -1;
	}

	//Split ra & dec fields on :
	std::vector<std::string> ra_parsed_fields= CodeUtils::SplitStringOnPattern(ra_field,':');
	if(ra_parsed_fields.size()!=3){
		WARN_LOG("Invalid number of RA parsed fields (n="<<ra_parsed_fields.size()<<") when 3 expected!");
		return -1;
	}

	std::vector<std::string> dec_parsed_fields= CodeUtils::SplitStringOnPattern(dec_field,':');
	if(dec_parsed_fields.size()!=3){
		WARN_LOG("Invalid number of DEC parsed fields (n="<<dec_parsed_fields.size()<<") when 3 expected!");
		return -1;
	}

	//Create IAU coordinates
	std::stringstream ss;
	ss<<cs;
	//if(ra_parsed_fields[0]>0) ss<<"+";
	for(size_t i=0;i<ra_parsed_fields.size();i++) {
		int l= ra_parsed_fields[i].size();
		if(l==1) ss<<0<<ra_parsed_fields[i]; 
		else ss<<ra_parsed_fields[i];
	}
 	if(atol(dec_parsed_fields[0].c_str())>0) ss<<"+";
	for(size_t i=0;i<dec_parsed_fields.size();i++) {
		int l= dec_parsed_fields[i].size();
		if(l==1) ss<<0<<dec_parsed_fields[i]; 
		else ss<<dec_parsed_fields[i];
	}
 	
	iau= ss.str();	

	return 0;

}//close GetIAUCoords()


std::string AstroUtils::GetDS9WCSTypeHeader(int coordSys)
{
	std::string ds9CoordSys= "image";
	if(coordSys==eJ2000) ds9CoordSys= "fk5";
	else if(coordSys==eB1950) ds9CoordSys= "fk4";
	else if(coordSys==eGALACTIC) ds9CoordSys= "galactic";
	else{
		WARN_LOG("Unknown coord sys given ("<<coordSys<<"), setting to 'image'");
		ds9CoordSys= "image";
	}

	return ds9CoordSys;

}//close GetDS9WCSTypeHeader()

/*
int AstroUtils::PixelToWCSCoords(double& xpos, double& ypos,WorldCoor* wcs,double ix,double iy) 
{
	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);
	
	return 0;

}//close PixelToWCSCoords()
*/

int AstroUtils::PixelToWCSCoords(double& xpos, double& ypos,WCS* wcs,double ix,double iy) 
{
	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	WCSUtils::pix2wcs (wcs,ix,iy,&xpos, &ypos);
	
	return 0;

}//close PixelToWCSCoords()

/*
int AstroUtils::PixelToWCSStrCoords(std::string& wcs_str,WorldCoor* wcs,double ix,double iy,int max_str_length) 
{
	//Init str
	wcs_str= "";

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}
	
	if(max_str_length<=0){
		ERROR_LOG("Invalid max wcs string size given (must be >0)!");
		return -1;
	}

	//Convert coords
	char data[max_str_length];
	int status= pix2wcst (wcs,ix,iy,data,max_str_length);
	if(status==0){
		WARN_LOG("Failed to convert pixel coords ("<<ix<<","<<iy<<") to WCS string!");
		return -1;
	}	
	wcs_str= std::string(data);

	return 0;

}//close PixelToWCSStrCoords()
*/

int AstroUtils::PixelToWCSStrCoords(std::string& wcs_str,WCS* wcs,double ix,double iy,int max_str_length) 
{
	//Init str
	wcs_str= "";

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}
	
	if(max_str_length<=0){
		ERROR_LOG("Invalid max wcs string size given (must be >0)!");
		return -1;
	}

	//Convert coords
	char data[max_str_length];
	int status= WCSUtils::pix2wcst (wcs,ix,iy,data,max_str_length);
	if(status==0){
		WARN_LOG("Failed to convert pixel coords ("<<ix<<","<<iy<<") to WCS string!");
		return -1;
	}	
	wcs_str= std::string(data);

	return 0;

}//close PixelToWCSStrCoords()

/*
int AstroUtils::PixelToWCSCoords(Caesar::Image* image,WorldCoor* wcs,double ix,double iy,double& xpos, double& ypos,bool useImageCoords) {

	//## NB: ix & iy shall be the pixel coordinates (image matrix coordinates).
	//## When a tile is read in FITS the header keywords (crpix) are modified accordingly so you need to pass the new matrix coordinates (not the (x,y) corresponding to full image) 
	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();
	
	if(useImageCoords && (ix<xmin || iy<ymin || ix>xmax || iy>ymax) ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);
	
	return 0;

}//close PixelToWCSCoords()
*/



int AstroUtils::PixelToWCSCoords(Caesar::Image* image,WCS* wcs,double ix,double iy,double& xpos, double& ypos,bool useImageCoords) {

	//## NB: ix & iy shall be the pixel coordinates (image matrix coordinates).
	//## When a tile is read in FITS the header keywords (crpix) are modified accordingly so you need to pass the new matrix coordinates (not the (x,y) corresponding to full image) 
	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();
	
	if(useImageCoords && (ix<xmin || iy<ymin || ix>xmax || iy>ymax) ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	WCSUtils::pix2wcs (wcs,ix,iy,&xpos, &ypos);
	
	return 0;

}//close PixelToWCSCoords()

/*
int AstroUtils::PixelToWCSStrCoords(Caesar::Image* image,WorldCoor* wcs,double ix,double iy,std::string& wcs_str,bool useImageCoords,int max_str_length) 
{
	//Init str
	wcs_str= "";

	//## NB: ix & iy shall be the pixel coordinates (image matrix coordinates).
	//## When a tile is read in FITS the header keywords (crpix) are modified accordingly so you need to pass the new matrix coordinates (not the (x,y) corresponding to full image) 
	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}

	if(max_str_length<=0){
		ERROR_LOG("Invalid max wcs string size given (must be >0)!");
		return -1;
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();
	
	if(useImageCoords && (ix<xmin || iy<ymin || ix>xmax || iy>ymax) ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	char data[max_str_length];
	int status= pix2wcst(wcs,ix,iy,data,max_str_length);
	if(status==0){
		WARN_LOG("Failed to convert pixel coords ("<<ix<<","<<iy<<") to WCS string!");
		return -1;
	}	
	wcs_str= std::string(data);

	return 0;

}//close PixelToWCSStrCoords()
*/


int AstroUtils::PixelToWCSStrCoords(Caesar::Image* image,WCS* wcs,double ix,double iy,std::string& wcs_str,bool useImageCoords,int max_str_length) 
{
	//Init str
	wcs_str= "";

	//## NB: ix & iy shall be the pixel coordinates (image matrix coordinates).
	//## When a tile is read in FITS the header keywords (crpix) are modified accordingly so you need to pass the new matrix coordinates (not the (x,y) corresponding to full image) 
	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}

	if(max_str_length<=0){
		ERROR_LOG("Invalid max wcs string size given (must be >0)!");
		return -1;
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();
	
	if(useImageCoords && (ix<xmin || iy<ymin || ix>xmax || iy>ymax) ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	char data[max_str_length];
	int status= WCSUtils::pix2wcst(wcs,ix,iy,data,max_str_length);
	if(status==0){
		WARN_LOG("Failed to convert pixel coords ("<<ix<<","<<iy<<") to WCS string!");
		return -1;
	}	
	wcs_str= std::string(data);

	return 0;

}//close PixelToWCSStrCoords()


int AstroUtils::PixelToWCSCoords(Caesar::Image* image,double ix,double iy,double& xpos, double& ypos,int coordSystem,bool useImageCoords) 
{
	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();

	if(useImageCoords && (ix<xmin || iy<ymin || ix>xmax || iy>ymax) ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check image meta-data
	if(!image->HasMetaData() ){
    ERROR_LOG("No metadata available in image!");
		return -1;
	}
	Caesar::ImgMetaData* metadata= image->GetMetaData();	
	
	//Get the coord system
	//WorldCoor* wcs= metadata->GetWorldCoord(coordSystem);
	WCS* wcs= metadata->GetWCS(coordSystem);
	if(!wcs){
		ERROR_LOG("Failed to get WorldCoord system from metadata!");
		return -1;
	}

	//Convert coords
	//pix2wcs (wcs,ix,iy,&xpos, &ypos);
	WCSUtils::pix2wcs (wcs,ix,iy,&xpos, &ypos);
	
	//Clear up
	delete wcs;
	wcs= 0;

	return 0;
		
}//close PixelToWCSCoords()


int AstroUtils::PixelToWCSStrCoords(Caesar::Image* image,double ix,double iy,std::string& wcs_str,int coordSystem,bool useImageCoords,int max_str_length) 
{
	//Init str
	wcs_str= "";

	//Check pixel values in input
	if(!image){
		ERROR_LOG("Null image ptr given!");
		return -1;	
	}
	if(max_str_length<=0){
		ERROR_LOG("Invalid max wcs string size given (must be >0)!");
		return -1;
	}

	//Get image range
	double xmin= image->GetXmin();
	double ymin= image->GetYmin();
	double xmax= image->GetXmax();
	double ymax= image->GetYmax();

	if(useImageCoords && (ix<xmin || iy<ymin || ix>xmax || iy>ymax) ){
		ERROR_LOG("Invalid pix range selected (ix="<<ix<<", iy="<<iy<<")");
		return -1;	
	}

	//Check image meta-data
	if(!image->HasMetaData() ){
    ERROR_LOG("No metadata available in image!");
		return -1;
	}
	Caesar::ImgMetaData* metadata= image->GetMetaData();	
	
	//Get the coord system
	//WorldCoor* wcs= metadata->GetWorldCoord(coordSystem);
	WCS* wcs= metadata->GetWCS(coordSystem);
	if(!wcs){
		ERROR_LOG("Failed to get WorldCoord system from metadata!");
		return -1;
	}

	//Convert coords
	char data[max_str_length];
	//int status= pix2wcst(wcs,ix,iy,data,max_str_length);	
	int status= WCSUtils::pix2wcst(wcs,ix,iy,data,max_str_length);	
	if(status==0){
		WARN_LOG("Failed to convert pixel coords ("<<ix<<","<<iy<<") to WCS string!");
		return -1;
	}	
	wcs_str= std::string(data);

	//Clear up
	delete wcs;
	wcs= 0;

	return 0;
		
}//close PixelToWCSStrCoords()

/*
Contour* AstroUtils::PixelToWCSContour(Contour* contour,WorldCoor* wcs,int pixOffset)
{
	//Check input data
	if(!contour){
		WARN_LOG("Null ptr to input contour given, returning nullptr!");
		return nullptr;
	}
	if(!wcs){
		WARN_LOG("Null ptr to input coord system given, nothing can be done!");
		return nullptr;
	}

	//Loop over contour points and convert in sky coords
	double x_wcs= 0;
	double y_wcs= 0;
	Contour* contour_wcs= new Contour;
	
	for(int i=0;i<contour->GetN();i++){
		//Get pix contour point
		TVector2* contPnt= contour->GetPoint(i);
		if(!contPnt) {
			WARN_LOG("Null ptr point retrieved from contour, skip it!");
			continue;
		}
		double x= contPnt->X() + pixOffset; 
		double y= contPnt->Y() + pixOffset;
		pix2wcs (wcs,x,y,&x_wcs, &y_wcs);

		//Fill new contour with sky coords
		contour_wcs->AddPoint(TVector2(x_wcs,y_wcs));
	}//end loop points

	return contour_wcs;

}//close PixelToWCSContour()
*/

Contour* AstroUtils::PixelToWCSContour(Contour* contour,WCS* wcs,int pixOffset)
{
	//Check input data
	if(!contour){
		WARN_LOG("Null ptr to input contour given, returning nullptr!");
		return nullptr;
	}
	if(!wcs){
		WARN_LOG("Null ptr to input coord system given, nothing can be done!");
		return nullptr;
	}

	//Loop over contour points and convert in sky coords
	double x_wcs= 0;
	double y_wcs= 0;
	Contour* contour_wcs= new Contour;
	
	for(int i=0;i<contour->GetN();i++){
		//Get pix contour point
		TVector2* contPnt= contour->GetPoint(i);
		if(!contPnt) {
			WARN_LOG("Null ptr point retrieved from contour, skip it!");
			continue;
		}
		double x= contPnt->X() + pixOffset; 
		double y= contPnt->Y() + pixOffset;
		WCSUtils::pix2wcs (wcs,x,y,&x_wcs, &y_wcs);

		//Fill new contour with sky coords
		contour_wcs->AddPoint(TVector2(x_wcs,y_wcs));
	}//end loop points

	return contour_wcs;

}//close PixelToWCSContour()

/*
int AstroUtils::PixelToWCSContours(std::vector<Contour*>& contours_wcs,std::vector<Contour*>const& contours,WorldCoor* wcs,int pixOffset)
{
	//Check input data
	if(!wcs){
		WARN_LOG("Null ptr to input coord system given, nothing can be done!");
		return -1;
	}
	if(contours.empty()){
		WARN_LOG("Empty contour list, nothing to be done!");
		return -1;
	}

	//Clear existing collection
	contours_wcs.clear();

	//Loop over contours and convert them
	for(size_t i=0;i<contours.size();i++){
		Contour* contour_wcs= PixelToWCSContour(contours[i],wcs,pixOffset);
		if(!contour_wcs){
			ERROR_LOG("Failed to convert contour no. "<<i+1<<"!");
			CodeUtils::DeletePtrCollection<Contour>(contours_wcs);
			return -1;
		}
		contours_wcs.push_back(contour_wcs);
	}//end loop contours

	return 0;

}//close PixelToWCSContours()
*/

int AstroUtils::PixelToWCSContours(std::vector<Contour*>& contours_wcs,std::vector<Contour*>const& contours,WCS* wcs,int pixOffset)
{
	//Check input data
	if(!wcs){
		WARN_LOG("Null ptr to input coord system given, nothing can be done!");
		return -1;
	}
	if(contours.empty()){
		WARN_LOG("Empty contour list, nothing to be done!");
		return -1;
	}

	//Clear existing collection
	contours_wcs.clear();

	//Loop over contours and convert them
	for(size_t i=0;i<contours.size();i++){
		Contour* contour_wcs= PixelToWCSContour(contours[i],wcs,pixOffset);
		if(!contour_wcs){
			ERROR_LOG("Failed to convert contour no. "<<i+1<<"!");
			CodeUtils::DeletePtrCollection<Contour>(contours_wcs);
			return -1;
		}
		contours_wcs.push_back(contour_wcs);
	}//end loop contours

	return 0;

}//close PixelToWCSContours()

/*
TEllipse* AstroUtils::PixelToWCSEllipse(TEllipse* ellipse,WorldCoor* wcs,int pixOffset)
{
	//NB: See Aegean source finder wcs_helpers.py method
	//Check input data
	if(!ellipse){
		WARN_LOG("Null ptr to input contour given, returning nullptr!");
		return nullptr;
	}
	if(!wcs){
		WARN_LOG("Null ptr to input coord system given, nothing can be done!");
		return nullptr;
	}

	//Get ellipse pars
	double x= ellipse->GetX1();//ellipse centroid x
	double y= ellipse->GetY1();//ellipse centroid y
	x+= pixOffset;
	y+= pixOffset;
	double sx= ellipse->GetR1();//ellipse semi-major axis
	double sy= ellipse->GetR2();//ellipse semi-minor axis
	double theta= ellipse->GetTheta();//rotation angle (wrt x axis)
	double theta_rad= theta*TMath::DegToRad();

	//Compute ellipse axis coordinates
	double x1= x + sx * cos(theta_rad);
	double y1= y + sx * sin(theta_rad);
	double x2= x + sy * cos(theta_rad - TMath::Pi()/2.);
	double y2= y + sy * sin(theta_rad - TMath::Pi()/2.);

	//Convert ellipse centroid in sky coords
	double x_wcs= 0;
	double y_wcs= 0;
	pix2wcs (wcs,x,y,&x_wcs, &y_wcs);

	//Convert ellipse axis coord in sky coords
	double x1_wcs= 0;
	double y1_wcs= 0;
	double x2_wcs= 0;
	double y2_wcs= 0;
	pix2wcs (wcs,x1,y1,&x1_wcs, &y1_wcs);
	pix2wcs (wcs,x2,y2,&x2_wcs, &y2_wcs);

	//Compute ellipse axis points
	//double semimajor_wcs= GetWCSPointDist_Vincenty(x_wcs,y_wcs,x1_wcs,y1_wcs);
	//double semiminor_wcs= GetWCSPointDist_Vincenty(x_wcs,y_wcs,x2_wcs,y2_wcs);
	double semimajor_wcs= GetWCSPointDist_Haversine(x_wcs,y_wcs,x1_wcs,y1_wcs);
	double semiminor_wcs= GetWCSPointDist_Haversine(x_wcs,y_wcs,x2_wcs,y2_wcs);

	//Compute ellipse rot angle
	//NB: Bearing is returned from North to East (0 deg is north) but we want theta measured from x-axis so sum 90 deg
	double theta_wcs = GetWCSPointBearing(x_wcs,y_wcs,x1_wcs,y1_wcs);
	double theta2_wcs = GetWCSPointBearing(x_wcs,y_wcs,x2_wcs,y2_wcs) - 90;
	theta_wcs+= 90.;
	theta2_wcs+= 90.;
	double dtheta_wcs= theta_wcs-theta2_wcs; 

	//Correct ellipse minor axis
	semiminor_wcs*= fabs(cos(dtheta_wcs*TMath::DegToRad()));

	DEBUG_LOG("(x,y,sx,sy,theta)=("<<std::setprecision(8)<<x<<","<<y<<","<<sx<<","<<sy<<","<<theta<<"), (x1,y1)=("<<x1<<","<<y1<<"), (x2,y2)=("<<x2<<","<<y2<<"), (x1_wcs,y1_wcs)=("<<x1_wcs<<","<<y1_wcs<<"), (x2_wcs,y2_wcs)=("<<x2_wcs<<","<<y2_wcs<<"), (x_wcs,y_wcs,a,b,theta)=("<<x_wcs<<","<<y_wcs<<","<<semimajor_wcs<<","<<semiminor_wcs<<","<<theta_wcs<<"), theta2_wcs="<<theta2_wcs<<", dtheta_wcs="<<dtheta_wcs);

	TEllipse* ellipse_wcs= new TEllipse(x_wcs,y_wcs,semimajor_wcs,semiminor_wcs,0.,360.,theta_wcs);
	ellipse_wcs->SetLineWidth(2);
	ellipse_wcs->SetFillColor(0);
	ellipse_wcs->SetFillStyle(0);

	return ellipse_wcs;

}//close PixelToWCSEllipse()
*/


TEllipse* AstroUtils::PixelToWCSEllipse(TEllipse* ellipse,WCS* wcs,int pixOffset)
{
	//NB: See Aegean source finder wcs_helpers.py method
	//Check input data
	if(!ellipse){
		WARN_LOG("Null ptr to input contour given, returning nullptr!");
		return nullptr;
	}
	if(!wcs){
		WARN_LOG("Null ptr to input coord system given, nothing can be done!");
		return nullptr;
	}

	//Get ellipse pars
	double x= ellipse->GetX1();//ellipse centroid x
	double y= ellipse->GetY1();//ellipse centroid y
	x+= pixOffset;
	y+= pixOffset;
	double sx= ellipse->GetR1();//ellipse semi-major axis
	double sy= ellipse->GetR2();//ellipse semi-minor axis
	double theta= ellipse->GetTheta();//rotation angle (wrt x axis)
	double theta_rad= theta*TMath::DegToRad();

	//Compute ellipse axis coordinates
	double x1= x + sx * cos(theta_rad);
	double y1= y + sx * sin(theta_rad);
	double x2= x + sy * cos(theta_rad - TMath::Pi()/2.);
	double y2= y + sy * sin(theta_rad - TMath::Pi()/2.);

	//Convert ellipse centroid in sky coords
	double x_wcs= 0;
	double y_wcs= 0;
	WCSUtils::pix2wcs (wcs,x,y,&x_wcs, &y_wcs);

	//Convert ellipse axis coord in sky coords
	double x1_wcs= 0;
	double y1_wcs= 0;
	double x2_wcs= 0;
	double y2_wcs= 0;
	WCSUtils::pix2wcs (wcs,x1,y1,&x1_wcs, &y1_wcs);
	WCSUtils::pix2wcs (wcs,x2,y2,&x2_wcs, &y2_wcs);

	//Compute ellipse axis points
	//double semimajor_wcs= GetWCSPointDist_Vincenty(x_wcs,y_wcs,x1_wcs,y1_wcs);
	//double semiminor_wcs= GetWCSPointDist_Vincenty(x_wcs,y_wcs,x2_wcs,y2_wcs);
	double semimajor_wcs= GetWCSPointDist_Haversine(x_wcs,y_wcs,x1_wcs,y1_wcs);
	double semiminor_wcs= GetWCSPointDist_Haversine(x_wcs,y_wcs,x2_wcs,y2_wcs);

	//Compute ellipse rot angle
	//NB: Bearing is returned from North to East (0 deg is north) but we want theta measured from x-axis so sum 90 deg
	double theta_wcs = GetWCSPointBearing(x_wcs,y_wcs,x1_wcs,y1_wcs);
	double theta2_wcs = GetWCSPointBearing(x_wcs,y_wcs,x2_wcs,y2_wcs) - 90;
	theta_wcs+= 90.;
	theta2_wcs+= 90.;
	double dtheta_wcs= theta_wcs-theta2_wcs; 

	//Correct ellipse minor axis
	semiminor_wcs*= fabs(cos(dtheta_wcs*TMath::DegToRad()));

	DEBUG_LOG("(x,y,sx,sy,theta)=("<<std::setprecision(8)<<x<<","<<y<<","<<sx<<","<<sy<<","<<theta<<"), (x1,y1)=("<<x1<<","<<y1<<"), (x2,y2)=("<<x2<<","<<y2<<"), (x1_wcs,y1_wcs)=("<<x1_wcs<<","<<y1_wcs<<"), (x2_wcs,y2_wcs)=("<<x2_wcs<<","<<y2_wcs<<"), (x_wcs,y_wcs,a,b,theta)=("<<x_wcs<<","<<y_wcs<<","<<semimajor_wcs<<","<<semiminor_wcs<<","<<theta_wcs<<"), theta2_wcs="<<theta2_wcs<<", dtheta_wcs="<<dtheta_wcs);

	TEllipse* ellipse_wcs= new TEllipse(x_wcs,y_wcs,semimajor_wcs,semiminor_wcs,0.,360.,theta_wcs);
	ellipse_wcs->SetLineWidth(2);
	ellipse_wcs->SetFillColor(0);
	ellipse_wcs->SetFillStyle(0);

	return ellipse_wcs;

}//close PixelToWCSEllipse()


int AstroUtils::GetBeamDeconvolvedEllipsePars (
	double& bmaj_deconv,double& bmin_deconv,double& bpa_deconv,
	double bmaj,double bmin,double bpa,
	double bmaj_beam,double bmin_beam,double bpa_beam
)
{
	int status= 0;
	bmaj_deconv= 0;
	bmin_deconv= 0;
	bpa_deconv= 0;

	double sum2= bmaj*bmaj + bmin*bmin;
	double diff2= bmaj*bmaj - bmin*bmin;
	double sum2_beam= bmaj_beam*bmaj_beam + bmin_beam*bmin_beam;
	double diff2_beam= bmaj_beam*bmaj_beam - bmin_beam*bmin_beam;
	double bpa_rad= bpa*TMath::DegToRad();//in rad
	double bpa_beam_rad= bpa_beam*TMath::DegToRad();//in rad

	double beta2= pow(diff2,2) + pow(diff2_beam,2) - 2*diff2*diff2_beam*cos(2*(bpa_rad-bpa_beam_rad));
	if(beta2<0) {
		WARN_LOG("Numerical error (beta^2 is <0 in formula)!");
		status= -1;
	}
	double beta= sqrt(beta2);
	
	double bmaj2_deconv= 0.5*(sum2 - sum2_beam + beta);
	double bmin2_deconv= 0.5*(sum2 - sum2_beam - beta);
	if(bmaj2_deconv<0 || bmin2_deconv) {
		WARN_LOG("Numerical error (bmaj^2/bmin^2 deconvolved are <0 in formula)!");	
		status= -1;
	}
	bmaj_deconv= sqrt(bmaj2_deconv);
	bmin_deconv= sqrt(bmin2_deconv);

	double arg= (diff2*sin(2*bpa))/(diff2*cos(2*bpa-diff2_beam));
	bpa_deconv= 0.5*atan(arg)*TMath::RadToDeg();

	return status;

}//close GetBeamDeconvolvedEllipsePars()


TEllipse* AstroUtils::GetBeamDeconvolvedEllipse(TEllipse* ellipse,TEllipse* beam)
{
	//Check input data
	if(!ellipse || !beam){
		WARN_LOG("Null ptr to input ellipses given, returning nullptr!");
		return nullptr;
	}

	//Get ellipse axis pars
	double X0= ellipse->GetX1();
	double Y0= ellipse->GetY1();
	double D1= 2*ellipse->GetR1();//ellipse axis along x
	double D2= 2*ellipse->GetR2();//ellipse axis along y
	double theta= ellipse->GetTheta();//rotation angle (wrt x axis)
	double bmaj= std::max(D1,D2);
	double bmin= std::min(D1,D2);
	double bpa= MathUtils::Mod(theta,180.);
	if(D2>D1) bpa= MathUtils::Mod(theta-90.,180.);

	//Get beam ellipse axis pars
	double D1_beam= 2*beam->GetR1();//ellipse axis along x
	double D2_beam= 2*beam->GetR2();//ellipse axis along y
	double theta_beam= beam->GetTheta();//rotation angle (wrt x axis)
	double bmaj_beam= std::max(D1,D2);
	double bmin_beam= std::min(D1,D2);
	double bpa_beam= MathUtils::Mod(theta_beam,180.);
	if(D2_beam>D1_beam) bpa_beam= MathUtils::Mod(theta_beam-90.,180.);

	double bmaj_deconv= 0;
	double bmin_deconv= 0;
	double bpa_deconv= 0;
	int status= GetBeamDeconvolvedEllipsePars(
		bmaj_deconv,bmin_deconv,bpa_deconv,
		bmaj,bmin,bpa,
		bmaj_beam,bmin_beam,bpa_beam
	);

	if(status<0){
		WARN_LOG("Failed to compute deconvolved ellipse pars (check given ellipse/beam pars)!");
		return nullptr;
	}

	double R1_deconv= bmaj_deconv/2.;
	double R2_deconv= bmin_deconv/2.;
	double theta_deconv= bpa_deconv;
	TEllipse* ellipse_deconv= new TEllipse(X0,Y0,R1_deconv,R2_deconv,theta_deconv);

	return ellipse_deconv;

}//close GetBeamDeconvolvedEllipse()


double AstroUtils::GetWCSPointDist_Haversine(double ra1,double dec1,double ra2,double dec2)
{
	//Convert to radians
	ra1*= TMath::DegToRad();
	ra2*= TMath::DegToRad();
	dec1*= TMath::DegToRad();
	dec2*= TMath::DegToRad();

	//Compute angular distance using Haversine formula (e.g. see https://en.wikipedia.org/wiki/Great-circle_distance)
	double dlon = ra2 - ra1;
  double dlat = dec2 - dec1;
	double a= pow(sin(dlat/2.),2);
	double b= cos(dec1)*cos(dec2)*pow(sin(dlon)/2.,2);
	double w= sqrt(a+b);
	double ang_sep= 2 * asin(std::min(1.,w));

	//Convert to degrees
	ang_sep*= TMath::RadToDeg();

	return ang_sep;

}//close GetWCSPointDist_Haversine()


double AstroUtils::GetWCSPointDist_Vincenty(double ra1,double dec1,double ra2,double dec2)
{
	//Convert to radians
	ra1*= TMath::DegToRad();
	ra2*= TMath::DegToRad();
	dec1*= TMath::DegToRad();
	dec2*= TMath::DegToRad();

	//Compute angular distance using Vincenty formula (e.g. see https://en.wikipedia.org/wiki/Great-circle_distance)
	double dlon = ra2 - ra1;
  double dlat = dec2 - dec1;
	
	double a= pow(cos(dec2)*sin(dlon),2);
	double b= pow(cos(dec1)*sin(dec2)-sin(dec1)*cos(dec2)*cos(dlon),2);
	double numer= sqrt(a+b); 
	double denom= sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(dlon);
	double ang_sep= atan2(numer,denom);

	//Convert to degrees
	ang_sep*= TMath::RadToDeg();

	return ang_sep;

}//close GetWCSPointDist_Haversine()


double AstroUtils::GetWCSPointBearing(double ra1,double dec1,double ra2,double dec2)
{
	//Calculate the bearing of point 2 from point 1 along a great circle.
  //The bearing is East of North and is in [0, 360), whereas position angle is also East of North but (-180,180]
	//Convert in radians

	DEBUG_LOG("(ra1,dec1)=("<<std::setprecision(9)<<ra1<<","<<dec1<<"), (ra2,dec2)=("<<ra2<<","<<dec2<<")");

	ra1*= TMath::DegToRad();
	ra2*= TMath::DegToRad();
	dec1*= TMath::DegToRad();
	dec2*= TMath::DegToRad();

	double dlon= ra2-ra1;
	double y= sin(dlon)*cos(dec2);
	double x= cos(dec1)*sin(dec2) - sin(dec1)*cos(dec2)*cos(dlon);
	
	double bear= atan2(y,x);
	DEBUG_LOG("Bear dlon="<<std::setprecision(9)<<dlon<<", (y,x)=("<<y<<","<<x<<"), bear(rad)="<<bear<<", bear(deg)="<<bear*TMath::RadToDeg());

	bear*= TMath::RadToDeg();//Convert to degrees
	
	return bear;

}//close GetWCSPointBearing()

std::string AstroUtils::EllipseToDS9Region(TEllipse* ellipse,std::string text,std::string color,std::vector<std::string> tags,bool useImageCoords)
{
	//Check input data
	if(!ellipse){
		ERROR_LOG("Null ptr to input ellipse given!");
		return std::string("");
	}
	
	//Get ellipse pars
	double x0= ellipse->GetX1();
	double y0= ellipse->GetY1();
	if(useImageCoords){//DS9 starts pix numbering from 1
		x0++;
		y0++;
	}
	double R1= ellipse->GetR1()*2;//DS9 wants axis (not semi-axis)
	double R2= ellipse->GetR2()*2;//DS9 wants axis (not semi-axis)
	double theta= ellipse->GetTheta();
	//theta-= 90;//DS9 format??

	//Encode to string
	std::stringstream sstream;
	sstream<<"ellipse "<<x0<<" "<<y0<<" "<<R1<<" "<<R2<<" "<<theta<<" ";
	sstream<<"# ";
	sstream<<"text={"<<text<<"} ";
	sstream<<"color="<<color<<" ";
	for(size_t k=0;k<tags.size();k++){
		sstream<<"tag={"<<tags[k]<<"} ";
	}

	return sstream.str();

}//close EllipseToDS9Region()


std::string AstroUtils::ContourToDS9Region(Contour* contour,std::string text,std::string color,std::vector<std::string> tags,bool useImageCoords)
{
	//Check input data
	if(!contour){
		ERROR_LOG("Null ptr to input contour given!");
		return std::string("");
	}

	//Check size of contour
	//NB: DS9 crashes miserably when given a polygon region with one point 
	int nPoints= contour->GetN();
	if(nPoints<=1){
		WARN_LOG("Too few points (<=1) present in contour, will return empty string!");
		return std::string("");
	}
	
	//Loop over contour points and encode to string
	std::stringstream sstream;
	sstream<<"polygon ";
	for(int j=0;j<nPoints;j++){
		TVector2* contPnt= contour->GetPoint(j);
		if(!contPnt) {
			WARN_LOG("Null ptr to contour point, skip to next!");
			continue;
		}
		if(useImageCoords){
			long int x= static_cast<long int>(contPnt->X()) + 1;//DS9 starts pix numbering from 1
			long int y= static_cast<long int>(contPnt->Y()) + 1;//DS9 starts pix numbering from 1
			sstream<<x<<" "<<y<<" ";
		}
		else{
			double x= contPnt->X();
			double y= contPnt->Y();
			sstream<<std::fixed<<std::setprecision(4)<<x<<" "<<y<<" ";
		}
		
	}//end loop contour points
	sstream<<"# ";

	sstream<<"text={"<<text<<"} ";
	sstream<<"color="<<color<<" ";
	for(size_t k=0;k<tags.size();k++){
		sstream<<"tag={"<<tags[k]<<"} ";
	}
	
	return sstream.str();

}//close ContourToDS9Region()


}//close namespace



