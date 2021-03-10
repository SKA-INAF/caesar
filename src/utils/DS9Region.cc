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
* @file DS9Region.cc
* @class DS9Region
* @brief DS9 region class
*
* DS9 region class
* @author S. Riggi
* @date 27/02/2019
*/

#include <DS9Region.h>
#include <MathUtils.h>
#include <CodeUtils.h>
#include <Contour.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//OpenCV headers
#include <opencv2/core/core.hpp>

//ROOT headers
#include <TMath.h>
#include <TEllipse.h>

//OpenCV headers
#include <opencv2/core/core.hpp>

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

ClassImp(Caesar::DS9RegionMetaData)
ClassImp(Caesar::DS9Region)
ClassImp(Caesar::DS9PolygonRegion)
ClassImp(Caesar::DS9BoxRegion)
ClassImp(Caesar::DS9CircleRegion)
ClassImp(Caesar::DS9EllipseRegion)

namespace Caesar {

//===========================
//==    DS9REGION CLASS
//===========================
DS9Region::DS9Region()
{
	shapeType= eUNKNOWN_SHAPE;
	csType= eUNKNOWN_CS;
	hasMetaDataSet= false;
	x0= 0;
	y0= 0;
	x1= 0;
	x2= 0;
	y1= 0;
	y2= 0;

}//close constructor

DS9Region::DS9Region(int shape,int cs)
	: shapeType(shape), csType(cs)
{
	hasMetaDataSet= false;
	x0= 0;
	y0= 0;
	x1= 0;
	x2= 0;
	y1= 0;
	y2= 0;

}//close constrcutor

DS9Region::~DS9Region()
{

}//close destructor

//================================
//==    DS9 POLYGON REGION CLASS
//================================
DS9PolygonRegion::DS9PolygonRegion(int cs)
	: DS9Region(ePOLYGON_SHAPE,cs)
{

}//close constructor

DS9PolygonRegion::~DS9PolygonRegion()
{

}//close destructor

void DS9PolygonRegion::Print()
{
	cout<<std::setprecision(12)<<"POLYGON REGION: {";
	for(size_t i=0;i<points.size()-1;i++) cout<<"("<<points[i].X()<<","<<points[i].Y()<<"), ";
	cout<<"("<<points[points.size()-1].X()<<","<<points[points.size()-1].Y()<<")}"<<endl;

}//close Print()


bool DS9PolygonRegion::IsPointInsideRegion(double x,double y)
{
	return MathUtils::IsPointInsidePolygon(x,y,points);

}//close IsPointInsideRegion()


bool DS9PolygonRegion::IsPointInsideRegion(TVector2 p)
{
	return MathUtils::IsPointInsidePolygon(p,points);

}//close IsPointInsideRegion()


TEllipse* DS9PolygonRegion::GetEllipse()
{
	//Check if polygon has points
	if(points.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No points found in polygon, returning nullptr!");
		#endif
		return nullptr;
	}

	//Create OpenCV point collection
	std::vector<cv::Point2f> points_cv;
	for(size_t i=0;i<points.size();i++){
		double x= points[i].X();
		double y= points[i].Y();
		points_cv.push_back( cv::Point2f(x,y) );	
	}

	//Compute rotated bounding box
	cv::RotatedRect MinBoundingRect= cv::minAreaRect(points_cv);

	double MinBoundingRect_height= MinBoundingRect.size.height;
	double MinBoundingRect_width= MinBoundingRect.size.width;
	TVector2 BoundingBoxCenter= TVector2(MinBoundingRect.center.x,MinBoundingRect.center.y);
	double BoundingBoxMaj= std::max(MinBoundingRect_height,MinBoundingRect_width);
	double BoundingBoxMin= std::min(MinBoundingRect_height,MinBoundingRect_width);
	double BoundingBoxAngle= MinBoundingRect.angle;//counterclockwise in degree 
	if(MinBoundingRect.size.width < MinBoundingRect.size.height){
    BoundingBoxAngle+= 180;
  }
	else{
		BoundingBoxAngle+= 90;  
	}

	//If coord system is WCS modify theta to 180-theta as x-axis (RA) increases in the opposite direction
	if(csType!=eUNKNOWN_CS && csType!=eIMG_CS){
		BoundingBoxAngle= 180.-BoundingBoxAngle;
	}

	//Return the ellipse
	TEllipse* ellipse= new TEllipse;
	ellipse->SetX1(BoundingBoxCenter.X());
	ellipse->SetY1(BoundingBoxCenter.Y());
	ellipse->SetTheta(BoundingBoxAngle);
	ellipse->SetR1(BoundingBoxMaj/2.);
	ellipse->SetR2(BoundingBoxMin/2.);

	return ellipse;

}//close GetEllipse()

Contour* DS9PolygonRegion::GetContour(bool computePars)
{
	//Check if polygon has points
	if(points.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No points found in polygon, returning nullptr!");
		#endif
		return nullptr;
	}
	
	//Create and fill contour
	Contour* contour= new Contour;
	for(size_t i=0;i<points.size();i++){
		contour->AddPoint(points[i]);
	}

	//Compute contour parameters
	if(computePars && contour->ComputeParameters()<0){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/more failures occurred while computing contour parameters!");
		#endif
	}

	return contour;

}//close GetContour()

int DS9PolygonRegion::ComputeBoundingBox()
{
	//Check if polygon has points
	if(points.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No points found in polygon!");
		#endif
		return -1;
	}

	//Create OpenCV point collection
	std::vector<cv::Point2f> points_cv;
	
	for(size_t i=0;i<points.size();i++){
		double x= points[i].X();
		double y= points[i].Y();
		points_cv.push_back( cv::Point2f(x,y) );
		//cout<<"P="<<cv::Point2f(x,y)<<endl;	
	}

	//Compute moments
	cv::Moments moments = cv::moments(points_cv);
	double m00= moments.m00;
	double m10= moments.m10;
	double m01= moments.m01;
	
	//Compute centroid
	x0= m10/m00;
	y0= m01/m00;
	
	//Compute bounding box
	cv::Rect bbox= cv::boundingRect(points_cv);
	x1= bbox.x;
	y1= bbox.y;
	x2= x1 + bbox.width;
	y2= y1 + bbox.height;

	return 0;

}//close ComputeBoundingBox()


//================================
//==    DS9 BOX REGION CLASS
//================================
DS9BoxRegion::DS9BoxRegion(int cs)
	: DS9Region(eBOX_SHAPE,cs)
{
	cx= 0;
	cy= 0;
	height= 0;
	width= 0;
	theta= 0;
	xmin= 0;
	xmax= 0;
	ymin= 0;
	ymax= 0;

}//close constructor

DS9BoxRegion::~DS9BoxRegion()
{

}//close destructor

void DS9BoxRegion::Print()
{
	cout<<std::setprecision(12)<<"BOX REGION: "<<" x["<<xmin<<","<<xmax<<"], y["<<ymin<<","<<ymax<<"]"<<endl;

}//close Print()

//void DS9BoxRegion::ComputeBoxCoords()
int DS9BoxRegion::ComputeBoundingBox()
{
	//Set centroid & unrotated bbox
	double xmin_unrot= std::min(cx - width/2.,cx + width/2.); 
	double xmax_unrot= std::max(cx - width/2.,cx + width/2.); 
	double ymin_unrot= std::min(cy - height/2.,cy + height/2.);
	double ymax_unrot= std::max(cy - height/2.,cy + height/2.);

	std::vector<TVector2> vertex 
	{
		TVector2(xmin_unrot,ymin_unrot),//bottom left
		TVector2(xmax_unrot,ymin_unrot),//bottom right
		TVector2(xmax_unrot,ymax_unrot),//top right
		TVector2(xmin_unrot,ymax_unrot)//top left
	};

	x0= cx;
	y0= cy;
	x1= xmin_unrot;
	x2= xmax_unrot;
	y1= ymin_unrot;
	y2= ymax_unrot;
	
	//Compute vertex points rotated
	points.clear();
	points= vertex;
	for(size_t i=0;i<vertex.size();i++){
		double xrot= 0;
		double yrot= 0;
		MathUtils::ComputeRotatedCoords(xrot,yrot,vertex[i].X(),vertex[i].Y(),cx,cy,theta);
		points[i].SetX(xrot);
		points[i].SetY(yrot);
	}

	//Compute min/max coordinates
	xmin= 1.e+99;
	xmax= -1.e+99;
	ymin= 1.e+99;
	ymax= -1.e+99;

	for(size_t i=0;i<points.size();i++){
		double x= points[i].X();
		double y= points[i].Y();
		if(x<xmin) xmin= x;
		if(x>xmax) xmax= x;
		if(y<ymin) ymin= y;
		if(y>ymax) ymax= y;
	}

	return 0;
	
}//close ComputeBoundingBox()

bool DS9BoxRegion::IsPointInsideRegion(double x,double y)
{
	return MathUtils::IsPointInsidePolygon(x,y,points);

}//close IsPointInsideRegion()


bool DS9BoxRegion::IsPointInsideRegion(TVector2 p)
{
	return MathUtils::IsPointInsidePolygon(p,points);

}//close IsPointInsideRegion()

TEllipse* DS9BoxRegion::GetEllipse()
{
	//Compute ellipse bmaj/bmin
	double bmaj= std::max(height,width);
	double bmin= std::min(height,width);
	double angle= theta;
	if(width < height){
    angle+= 180;
  }
	else{
		angle+= 90;  
	}

	//If coord system is WCS modify theta to 180-theta as x-axis (RA) increases in the opposite direction
	if(csType!=eUNKNOWN_CS && csType!=eIMG_CS){
		angle= 180.-angle;
	}

	//Return the ellipse
	TEllipse* ellipse= new TEllipse;
	ellipse->SetX1(cx);
	ellipse->SetY1(cy);
	ellipse->SetTheta(angle);
	ellipse->SetR1(bmaj/2.);
	ellipse->SetR2(bmin/2.);

	return ellipse;

}//close GetEllipse()

Contour* DS9BoxRegion::GetContour(bool computePars)
{
	//Check if polygon has points
	if(points.empty()){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No points found in box (check if computed), returning nullptr!");
		#endif
		return nullptr;
	}
	
	//Create and fill contour
	Contour* contour= new Contour;
	for(size_t i=0;i<points.size();i++){
		contour->AddPoint(points[i]);
	}

	//Compute contour parameters
	if(computePars && contour->ComputeParameters()<0){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/more failures occurred while computing contour parameters!");
		#endif
	}

	return contour;
	
}//close GetContour()

//================================
//==    DS9 CIRCLE REGION CLASS
//================================
DS9CircleRegion::DS9CircleRegion(int cs)
	: DS9Region(eBOX_SHAPE,cs)
{
	cx= 0;
	cy= 0;
	r= 0;
	
}//close constructor

DS9CircleRegion::~DS9CircleRegion()
{

}//close destructor

void DS9CircleRegion::Print()
{
	cout<<std::setprecision(12)<<"CIRCLE REGION: C("<<cx<<","<<cy<<"), R="<<r<<endl;

}//close Print()


bool DS9CircleRegion::IsPointInsideRegion(double x,double y)
{
	double d= sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy));
	bool isInside= (d<r);
	
	return isInside;

}//close IsPointInsideRegion()


bool DS9CircleRegion::IsPointInsideRegion(TVector2 p)
{
	return this->IsPointInsideRegion(p.X(),p.Y());

}//close IsPointInsideRegion()

TEllipse* DS9CircleRegion::GetEllipse()
{
	//Compute ellipse bmaj/bmin
	double bmaj= 2*r;
	double bmin= 2*r;
	double angle= 0;
	
	//Return the ellipse
	TEllipse* ellipse= new TEllipse;
	ellipse->SetX1(cx);
	ellipse->SetY1(cy);
	ellipse->SetTheta(angle);
	ellipse->SetR1(bmaj/2.);
	ellipse->SetR2(bmin/2.);

	return ellipse;

}//close GetEllipse()

Contour* DS9CircleRegion::GetContour(bool computePars)
{
	//Create contour
	Contour* contour= new Contour;

	//Fill contour
	double t= 0;
	double t_step= 2;
	double t_max= 360;
	while(true){
		if(t>=t_max) break;	
		double t_rad= t*TMath::DegToRad();
		double x= cx + r*cos(t_rad);
		double y= cy + r*sin(t_rad);
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

int DS9CircleRegion::ComputeBoundingBox()
{
	x0= cx;
	y0= cy;
	x1= x0-r;
	x2= x0+r;
	y1= y0-r;
	y2= y0+r;

	return 0;

}//close ComputeBoundingBox()

//================================
//==    DS9 ELLIPSE REGION CLASS
//================================
DS9EllipseRegion::DS9EllipseRegion(int cs)
	: DS9Region(eELLIPSE_SHAPE,cs)
{
	cx= 0;
	cy= 0;
	a= 0;
	b= 0;
	theta= 0;

}//close constructor

DS9EllipseRegion::~DS9EllipseRegion()
{

}//close destructor

void DS9EllipseRegion::Print()
{
	cout<<std::setprecision(12)<<"ELLIPSE REGION: C("<<cx<<","<<cy<<"), a="<<a<<", b="<<b<<", theta="<<theta<<endl;

}//close Print()


bool DS9EllipseRegion::IsPointInsideRegion(double x,double y)
{
	double xdiff= x-cx;
	double ydiff= y-cy;
	double theta_rad= theta*TMath::DegToRad();
	double dx= cos(theta_rad)*xdiff + sin(theta_rad)*ydiff;
	double dy= sin(theta_rad)*xdiff - cos(theta_rad)*ydiff;
	double dx2= dx*dx;
	double dy2= dy*dy;
	double a2= a*a;
	double b2= b*b;
	
	double d= dx2/a2 + dy2/b2;
	bool isInside= (d<1);
	
	return isInside;

}//close IsPointInsideRegion()


bool DS9EllipseRegion::IsPointInsideRegion(TVector2 p)
{
	return this->IsPointInsideRegion(p.X(),p.Y());

}//close IsPointInsideRegion()

TEllipse* DS9EllipseRegion::GetEllipse()
{
	//Return the ellipse
	TEllipse* ellipse= new TEllipse;
	ellipse->SetX1(cx);
	ellipse->SetY1(cy);
	ellipse->SetTheta(theta);
	ellipse->SetR1(a);
	ellipse->SetR2(b);

	return ellipse;

}//close GetEllipse()

Contour* DS9EllipseRegion::GetContour(bool computePars)
{
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
		double x= cx + a*cos(t_rad)*cos(theta_rad) - b*sin(t_rad)*sin(theta_rad);
		double y= cy + a*cos(t_rad)*sin(theta_rad) + b*sin(t_rad)*cos(theta_rad);
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

int DS9EllipseRegion::ComputeBoundingBox()
{
	//Set center
	x0= cx;
	y0= cy;
	double width= a;
	double height= b;

	//Comput bbox (see http://stackoverflow.com/a/88020)
	double theta_rad= theta*TMath::DegToRad();
	double cos_angle = cos(theta_rad);
  double sin_angle = sin(theta_rad);
  double tan_angle = tan(theta_rad);

  double t1 = atan(-height * tan_angle / width);
  double t2 = t1 + TMath::Pi();

  double dx1 = 0.5 * width * cos_angle * cos(t1) - 0.5 * height * sin_angle * sin(t1);
  double dx2 = 0.5 * width * cos_angle * cos(t2) - 0.5 * height * sin_angle * sin(t2);

  if(dx1 > dx2){
		double dx1_tmp= dx1;
		double dx2_tmp= dx2;
		dx1= dx2_tmp;
		dx2= dx1_tmp;
	}

  t1 = atan(height / tan_angle / width);
  t2 = t1 + TMath::Pi();

  double dy1 = 0.5 * height * cos_angle * sin(t1) + 0.5 * width * sin_angle * cos(t1);
  double dy2 = 0.5 * height * cos_angle * sin(t2) + 0.5 * width * sin_angle * cos(t2);

  if(dy1 > dy2){
		double dy1_tmp= dy1;
		double dy2_tmp= dy2;
		dy1= dy2_tmp;
		dy2= dy1_tmp;
	}

  x1= x0 + dx1;
  x2= x0 + dx2;
  y1= y0 + dy1;
  y2= y0 + dy2;

	return 0;

}//close ComputeBoundingBox()

}//close namespace
