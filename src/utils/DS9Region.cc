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

ClassImp(Caesar::DS9Region)
ClassImp(Caesar::DS9PolygonRegion)
ClassImp(Caesar::DS9BoxRegion)
ClassImp(Caesar::DS9CircleRegion)

namespace Caesar {

//===========================
//==    DS9REGION CLASS
//===========================
DS9Region::DS9Region()
{
	shapeType= eUNKNOWN_SHAPE;
	csType= eUNKNOWN_CS;

}//close constructor

DS9Region::DS9Region(int shape,int cs)
	: shapeType(shape), csType(cs)
{
			
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
	cv::RotatedRect MinBoundingRect= cv::minAreaRect(points);

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

void DS9BoxRegion::ComputeBoxCoords()
{
	double xmin_unrot= std::min(cx - width/2.,cx + width/2.); 
	double xmax_unrot= std::max(cx - width/2.,cx + width/2.); 
	double ymin_unrot= std::min(cy - height/2.,cy + height/2.);
	double ymax_unrot= std::max(cy - height/2.,cy + height/2.);

	//VERTEX NOT ROTATED
	std::vector<TVector2> vertex 
	{
		TVector2(xmin_unrot,ymin_unrot),//bottom left
		TVector2(xmax_unrot,ymin_unrot),//bottom right
		TVector2(xmax_unrot,ymax_unrot),//top right
		TVector2(xmin_unrot,ymax_unrot)//top left
	};
	
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

	/*	
	MathUtils::ComputeRotatedCoords(xmin_rot,ymin_rot,xmin_unrot,ymin_unrot,cx,cy,theta);
	MathUtils::ComputeRotatedCoords(xmax_rot,ymax_rot,xmax_unrot,ymax_unrot,cx,cy,theta);

	xmin= std::min(xmin_rot,xmax_rot);
	xmax= std::max(xmin_rot,xmax_rot);
	ymin= std::min(ymin_rot,ymax_rot);
	ymax= std::max(ymin_rot,ymax_rot);

	points.push_back(TVector2(xmin,ymin));
	points.push_back(TVector2(xmax,ymin));
	points.push_back(TVector2(xmax,ymax));
	points.push_back(TVector2(xmin,ymax));
	*/
	
}//close ComputeBoxCoords()

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
	double theta= 0;
	double theta_step= 2;
	double theta_max= 360;
	while(true){
		if(theta>=theta_max) break;	
		double theta_rad= theta*TMath::DegToRad();
		double x= cx + r*cos(theta_rad);
		double y= cy + r*sin(theta_rad);
		contour->AddPoint(TVector2(x,y));
		theta+= theta_step;		
	}//end loop 
	
	//Compute contour parameters
	if(computePars && contour->ComputeParameters()<0){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("One/more failures occurred while computing contour parameters!");
		#endif
	}

	return contour;

}//close GetContour()

}//close namespace
