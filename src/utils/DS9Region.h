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
* @file DS9Region.h
* @class DS9Region
* @brief DS9 region class
*
* DS9 region class
* @author S. Riggi
* @date 27/02/2019
*/


#ifndef _DS9_REGION_h
#define _DS9_REGION_h 1

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

//ROOT headers
#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>

//C++ headers
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
#include <time.h>
#include <ctime>

namespace Caesar {

//===========================
//==    DS9REGION CLASS
//===========================
class DS9Region : public TObject
{
	public:
		/**
		* \brief Standard constructor
		*/
		DS9Region();
		/**
		* \brief Parametric constructor
		*/
		DS9Region(int shape=eUNKNOWN_SHAPE,int cs=eUNKNOWN_CS);
		/**
		* \brief Destructor
		*/
		virtual ~DS9Region();

	public:
		/**
		* \brief Region shape type
		*/
		enum ShapeType {
			eUNKNOWN_SHAPE=0,
			eCIRCLE_SHAPE=1,
			eBOX_SHAPE=2,
			ePOLYGON_SHAPE=3
		};

		/**
		* \brief Region coordinate system type
		*/
		enum CSType {
			eUNKNOWN_CS=0,
			eIMG_CS=1,
			eFK5_CS=2,
			eFK4_CS=3,
			eGAL_CS=4
		};

	public:
		/**
		* \brief Print region info
		*/
		virtual void Print(){};

		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(double x,double y){return false;}
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(TVector2 p){return false;}
		

	public:
		int shapeType;
		int csType;

	ClassDef(DS9Region,1)
		
};//close class DS9Region()

//=================================
//==    DS9 POLYGON REGION CLASS
//=================================
class DS9PolygonRegion : public DS9Region
{
	public:
		/**
		* \brief Parametric constructor
		*/
		DS9PolygonRegion(int cs=eUNKNOWN_CS);
		/**
		* \brief Destructor
		*/
		virtual ~DS9PolygonRegion();

	public:
		/**
		* \brief Print region info
		*/
		virtual void Print();
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(double x,double y);
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(TVector2 p);

	public:

		std::vector<TVector2> points;

	ClassDef(DS9PolygonRegion,1)

};//close class DS9PolygonRegion


//=================================
//==    DS9 BOX REGION CLASS
//=================================
class DS9BoxRegion : public DS9Region
{
	public:
		/**
		* \brief Parametric constructor
		*/
		DS9BoxRegion(int cs=eUNKNOWN_CS);
		/**
		* \brief Destructor
		*/
		virtual ~DS9BoxRegion();

	public:
		/**
		* \brief Print region info
		*/
		virtual void Print();
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(double x,double y);
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(TVector2 p);
		/**
		* \brief Compute box coordinates
		*/
		void ComputeBoxCoords();

	public:
	
		double cx;
		double cy;	
		double width;
		double height;
		double theta;
		double xmin;
		double xmax;
		double ymin;
		double ymax;
		std::vector<TVector2> points;

	ClassDef(DS9BoxRegion,1)

};//close class DS9BoxRegion


//=================================
//==    DS9 CIRCLE REGION CLASS
//=================================
class DS9CircleRegion : public DS9Region
{
	public:
		/**
		* \brief Parametric constructor
		*/
		DS9CircleRegion(int cs=eUNKNOWN_CS);
		/**
		* \brief Destructor
		*/
		virtual ~DS9CircleRegion();

	public:
		/**
		* \brief Print region info
		*/
		virtual void Print();
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(double x,double y);
		/**
		* \brief Is point inside region
		*/
		virtual bool IsPointInsideRegion(TVector2 p);

	public:

		double cx;
		double cy;
		double r;
		
	ClassDef(DS9CircleRegion,1)

};//close class DS9CircleRegion


#ifdef __MAKECINT__
#pragma link C++ class DS9Region+;
#pragma link C++ class vector<DS9Region>+;
#pragma link C++ class vector<DS9Region*>+;
#pragma link C++ class DS9PolygonRegion+;
#pragma link C++ class vector<DS9PolygonRegion>+;
#pragma link C++ class vector<DS9PolygonRegion*>+;
#pragma link C++ class DS9BoxRegion+;
#pragma link C++ class vector<DS9BoxRegion>+;
#pragma link C++ class vector<DS9BoxRegion*>+;
#pragma link C++ class DS9CircleRegion+;
#pragma link C++ class vector<DS9CircleRegion>+;
#pragma link C++ class vector<DS9CircleRegion*>+;
#endif

}//close namespace

#endif
