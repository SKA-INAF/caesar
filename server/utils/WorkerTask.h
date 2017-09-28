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
* @file WorkerTask.h
* @class WorkerTask
* @brief WorkerTask class
*
* Class for serializing task
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef WorkerTask_H
#define WorkerTask_H

#include <Source.pb.h>
#include <json/json.h>

#include <TObject.h>

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <time.h>
#include <sys/time.h>

namespace Caesar {

class WorkerTask : public TObject {

	public:
		//Standard constructor
		WorkerTask();
				
		//Copy constructor
		WorkerTask(const WorkerTask& task);

		//Destructor
		virtual ~WorkerTask();

	public:
		// Operator =
		WorkerTask& operator=(const WorkerTask& data);

		//Copy method
		void Copy(TObject &obj) const;

		//Serialize methods
		int SerializeToJson(Json::Value& jsonObj);
		int SerializeToJsonString(std::string& jsonString,bool isMinified=true);
		int SerializeToProtobuf(SourcePB::WorkerTask& taskPB);
		int EncodeFromJson(const Json::Value& jsonObj);
		int EncodeFromJsonString(std::string& jsonString);
		int EncodeFromProtobuf(const SourcePB::WorkerTask& taskPB);

	public:
		std::string filename; 
		std::string jobId;
		std::string worker_name;
		std::string broker_name;
		long int IdX;
		long int IdY;
		long int ix_min;
		long int ix_max;
		long int iy_min;
		long int iy_max;	
		timeval timestamp;	

	ClassDef(WorkerTask,1)
};

#ifdef __MAKECINT__
#pragma link C++ class WorkerTask+;
#endif

}//close namespace 

#endif


