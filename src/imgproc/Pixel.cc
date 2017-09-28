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
* @file Pixel.cc
* @class Pixel
* @brief Pixel data class
*
* Pixel class
* @author S. Riggi
* @date 20/01/2015
*/

#include <Pixel.h>
#include <Logger.h>

#include <TObject.h>

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

using namespace std;

ClassImp(Caesar::Pixel)

namespace Caesar {

Pixel::Pixel() : TObject() {
	
	Init();
	
}//close costructor

Pixel::Pixel(long int _gbin,long int _ix,long int _iy,double _x,double _y,double _S)
	: id(_gbin), S(_S), x(_x), y(_y), ix(_ix), iy(_iy)
{
	type= eNormal;
	S_curv= 0; 
	S_edge= 0;
	//isOnEdge= false; 
	//distanceToEdge= std::numeric_limits<double>::infinity();
	bkgLevel= 0;
	noiseLevel= 0;

}//close constructor


Pixel::~Pixel(){
	

}//close destructor


Pixel::Pixel(const Pixel& pixel) : TObject(pixel) {
  // Contour copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((Pixel&)pixel).Copy(*this);
}


void Pixel::Copy(TObject &obj) const {

	//Copy object
	TObject::Copy((Pixel&)obj);

	//Copy pixel members
  ((Pixel&)obj).id = id;
	((Pixel&)obj).type = type;
	((Pixel&)obj).S = S;
	((Pixel&)obj).x = x;
	((Pixel&)obj).y = y;
	((Pixel&)obj).ix = ix;
	((Pixel&)obj).iy = iy;
	//((Pixel&)obj).isOnEdge = isOnEdge;
	//((Pixel&)obj).distanceToEdge = distanceToEdge;
	((Pixel&)obj).S_curv = S_curv;
	((Pixel&)obj).S_edge = S_edge;
	((Pixel&)obj).bkgLevel = bkgLevel;
	((Pixel&)obj).noiseLevel = noiseLevel;

}//close Copy()

Pixel& Pixel::operator=(const Pixel& pixel) { 
	// Operator =
  if (this != &pixel)  ((Pixel&)pixel).Copy(*this);
  return *this;
}

void Pixel::Init(){

	//Init values
	id= -1;
	type= eNormal;
	S= 0; 
	S_curv= 0; 
	S_edge= 0;
	x= -1; 
	y=-1; 
	ix= -1; 
	iy= -1; 
	//isOnEdge= false; 
	//distanceToEdge= std::numeric_limits<double>::infinity();
	bkgLevel= 0;
	noiseLevel= 0;

}//close Init()

}//close namespace

