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
* @file SourceCube.cc
* @class SourceCube
* @brief Source cube class
*
* Class representing an image source cube
* @author S. Riggi
* @date 26/03/2018
*/

#include <SourceCube.h>
#include <Source.h>
#include <CodeUtils.h>

#include <Image.h>
#include <Contour.h>

#include <TObject.h>
#include <TMatrixD.h>

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

ClassImp(Caesar::SourceCube)

namespace Caesar {


SourceCube::SourceCube() 
	: TNamed()
{
	//Initialize 
	Init();

}//close costructor

SourceCube::SourceCube(std::string name)
	: TNamed(name.c_str(),name.c_str())
{
	//Init pars
	Init();

}//close constructor

SourceCube::~SourceCube()
{
	//Delete source collection
	DEBUG_LOG("Deleting source collection added in cube...");
	CodeUtils::DeletePtrCollection<Source>(m_sources);
	DEBUG_LOG("done!");
	
}//close destructor


SourceCube::SourceCube(const SourceCube& sourceCube) 
{
  // Copy constructor
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((SourceCube&)sourceCube).Copy(*this);
}

void SourceCube::Copy(TObject &obj) const 
{
	DEBUG_LOG("Copying parent TNamed...");
	TNamed::Copy((SourceCube&)obj);

	// Copy this source to source obj	
  //((SourceCube&)obj).Type = Type;
	
	//Copy source collection
	//Delete first any existing collection
	for(size_t i=0;i<(((SourceCube&)obj).m_sources).size();i++){
		if( (((SourceCube&)obj).m_sources)[i] ){
			delete (((SourceCube&)obj).m_sources)[i];
			(((SourceCube&)obj).m_sources)[i]= 0;
		}
	}
	(((SourceCube&)obj).m_sources).clear();

	Source* aSource= 0;
	for(size_t i=0;i<m_sources.size();i++){
		aSource= new Source;
		*aSource= *(m_sources[i]);
		(((SourceCube&)obj).m_sources).push_back(aSource);
	}

}//close Copy()

SourceCube& SourceCube::operator=(const SourceCube& sourceCube) 
{ 
	// Operator =
  if (this != &sourceCube) ((SourceCube&)sourceCube).Copy(*this);
  return *this;
}


void SourceCube::Init()
{
	//Init source cube data
	m_sources.clear();

}//close Init()

}//close namespace 


