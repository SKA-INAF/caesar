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
* @file TaskData.cc
* @class TaskData
* @brief TaskData class
*
* Task data class
* @author S. Riggi
* @date 20/01/2015
*/

#include <TaskData.h>
#include <Logger.h>
#include <Source.h>
#include <CodeUtils.h>

#include <json/json.h>

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

ClassImp(Caesar::TaskData)

namespace Caesar {

//Standard Constructor
TaskData::TaskData() : TObject() {	
	
	//Init task info
	workerId= -1;
	ix_min= -1;
	ix_max= -1;
	iy_min= -1;
	iy_max= -1;

	//Init sources
	source= 0;
	sources.clear();
	ext_sources.clear();
	sources_edge.clear();
	ext_sources_edge.clear();
	
	neighborTaskId.clear();
	neighborWorkerId.clear();
}

//Copy constructor
TaskData::TaskData(const TaskData& info) : TObject() {
	((TaskData&)info).Copy(*this);
}

//Destructor
TaskData::~TaskData(){
	
	//Clear sources
	ClearSources();

}//close destructor

// Operator =
TaskData& TaskData::operator=(const TaskData& data) { 
	if (this != &data) ((TaskData&)data).Copy(*this);
  return *this;
}

//Copy
void TaskData::Copy(TObject &obj) const {
			
	TObject::Copy((TaskData&)obj);
	((TaskData&)obj).workerId = workerId;
	((TaskData&)obj).ix_min = ix_min;	
	((TaskData&)obj).ix_max = ix_max;		
	((TaskData&)obj).iy_min = iy_min;	
	((TaskData&)obj).iy_max = iy_max;

	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((TaskData&)obj).sources).size();i++){
		if( (((TaskData&)obj).sources)[i] ){
			delete (((TaskData&)obj).sources)[i];
			(((TaskData&)obj).sources)[i]= 0;
		}
	}
	(((TaskData&)obj).sources).clear();

	//Copy sources
	((TaskData&)obj).source= 0;
	for(unsigned int i=0;i<sources.size();i++){
		((TaskData&)obj).source= new Source;
		*(((TaskData&)obj).source)= *(sources[i]);
		(((TaskData&)obj).sources).push_back( ((TaskData&)obj).source );
	}


	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((TaskData&)obj).sources_edge).size();i++){
		if( (((TaskData&)obj).sources_edge)[i] ){
			delete (((TaskData&)obj).sources_edge)[i];
			(((TaskData&)obj).sources_edge)[i]= 0;
		}
	}
	(((TaskData&)obj).sources_edge).clear();

	//Copy sources
	for(unsigned int i=0;i<sources_edge.size();i++){
		((TaskData&)obj).source= new Source;
		*(((TaskData&)obj).source)= *(sources_edge[i]);
		(((TaskData&)obj).sources_edge).push_back( ((TaskData&)obj).source );
	}

	
	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((TaskData&)obj).ext_sources).size();i++){
		if( (((TaskData&)obj).ext_sources)[i] ){
			delete (((TaskData&)obj).ext_sources)[i];
			(((TaskData&)obj).ext_sources)[i]= 0;
		}
	}
	(((TaskData&)obj).ext_sources).clear();

	//Copy sources
	for(unsigned int i=0;i<ext_sources.size();i++){
		((TaskData&)obj).source= new Source;
		*(((TaskData&)obj).source)= *(ext_sources[i]);
		(((TaskData&)obj).ext_sources).push_back( ((TaskData&)obj).source );
	}


	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((TaskData&)obj).ext_sources_edge).size();i++){
		if( (((TaskData&)obj).ext_sources_edge)[i] ){
			delete (((TaskData&)obj).ext_sources_edge)[i];
			(((TaskData&)obj).ext_sources_edge)[i]= 0;
		}
	}
	(((TaskData&)obj).ext_sources_edge).clear();

	//Copy sources
	for(unsigned int i=0;i<ext_sources_edge.size();i++){
		((TaskData&)obj).source= new Source;
		*(((TaskData&)obj).source)= *(ext_sources_edge[i]);
		(((TaskData&)obj).ext_sources_edge).push_back( ((TaskData&)obj).source );
	}
	
	//Clear first a previously existing vector
	((TaskData&)obj).neighborWorkerId.clear();
	((TaskData&)obj).neighborTaskId.clear();

	//Copy vector
	for(unsigned int i=0;i<neighborTaskId.size();i++){
		(((TaskData&)obj).neighborTaskId).push_back( neighborTaskId[i] );
	}
	for(unsigned int i=0;i<neighborWorkerId.size();i++){
		(((TaskData&)obj).neighborWorkerId).push_back( neighborWorkerId[i] );
	}

}//close Copy()

void TaskData::ClearSources(){

	//Clear
	for(size_t i=0;i<sources.size();i++){
		if(sources[i]){
			delete sources[i];
			sources[i]= 0;
		}
	}	
	sources.clear();

	for(size_t i=0;i<ext_sources.size();i++){
		if(ext_sources[i]){
			delete ext_sources[i];
			ext_sources[i]= 0;
		}
	}	
	ext_sources.clear();

	for(size_t i=0;i<sources_edge.size();i++){
		if(sources_edge[i]){
			delete sources_edge[i];
			sources_edge[i]= 0;
		}
	}	
	sources_edge.clear();

	for(size_t i=0;i<ext_sources_edge.size();i++){
		if(ext_sources_edge[i]){
			delete ext_sources_edge[i];
			ext_sources_edge[i]= 0;
		}
	}	
	ext_sources_edge.clear();

}//close ClearSources()



}//close namespace 

