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
* @file WorkerData.cc
* @class WorkerData
* @brief WorkerData class
*
* Worker data class
* @author S. Riggi
* @date 20/01/2015
*/

#include <WorkerData.h>
#include <WorkerTask.h>
#include <Logger.h>
#include <Source.h>

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <exception>

#include <chrono>

using namespace std;

//ClassImp(Caesar::WorkerDataHeader)
ClassImp(Caesar::WorkerData)

namespace Caesar {

//Standard Constructor
WorkerData::WorkerData() : TObject() {	
	data_type= eUNKNOWN_DATA;
	source= 0;
	sources.clear();
	edge_source= 0;
	edge_sources.clear();
}

//Copy constructor
WorkerData::WorkerData(const WorkerData& info) : TObject() {
	((WorkerData&)info).Copy(*this);
}

//Destructor
WorkerData::~WorkerData(){
	
	for(unsigned int i=0;i<sources.size();i++){
		if(sources[i]){
			delete sources[i];
			sources[i]= 0;
		}
	}	
	sources.clear();

	for(unsigned int i=0;i<edge_sources.size();i++){
		if(edge_sources[i]){
			delete edge_sources[i];
			edge_sources[i]= 0;
		}
	}	
	edge_sources.clear();

}//close destructor

// Operator =
WorkerData& WorkerData::operator=(const WorkerData& data) { 
	if (this != &data) ((WorkerData&)data).Copy(*this);
  return *this;
}

//Copy
void WorkerData::Copy(TObject &obj) const {
			
	TObject::Copy((WorkerData&)obj);
	((WorkerData&)obj).info = info;
	((WorkerData&)obj).data_type = data_type;		
	
	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((WorkerData&)obj).sources).size();i++){
		if( (((WorkerData&)obj).sources)[i] ){
			delete (((WorkerData&)obj).sources)[i];
			(((WorkerData&)obj).sources)[i]= 0;
		}
	}
	(((WorkerData&)obj).sources).clear();

	//Copy sources
	((WorkerData&)obj).source= 0;
	for(unsigned int i=0;i<sources.size();i++){
		((WorkerData&)obj).source= new Source;
		*(((WorkerData&)obj).source)= *(sources[i]);
		(((WorkerData&)obj).sources).push_back( ((WorkerData&)obj).source );
	}

	//Delete first a previously existing edge source vector
	for(unsigned int i=0;i<(((WorkerData&)obj).edge_sources).size();i++){
		if( (((WorkerData&)obj).edge_sources)[i] ){
			delete (((WorkerData&)obj).edge_sources)[i];
			(((WorkerData&)obj).edge_sources)[i]= 0;
		}
	}
	(((WorkerData&)obj).edge_sources).clear();

	//Copy sources
	((WorkerData&)obj).edge_source= 0;
	for(unsigned int i=0;i<edge_sources.size();i++){
		((WorkerData&)obj).edge_source= new Source;
		*(((WorkerData&)obj).edge_source)= *(edge_sources[i]);
		(((WorkerData&)obj).edge_sources).push_back( ((WorkerData&)obj).edge_source );
	}

}//close Copy()



}//close namespace 
