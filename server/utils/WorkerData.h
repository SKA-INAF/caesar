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
* @file WorkerData.h
* @class WorkerData
* @brief WorkerData class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef WorkerData_H
#define WorkerData_H

#include <WorkerTask.h>
#include <Source.h>

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

/*
class WorkerDataHeader : public TObject {
	public:
		//Standard constructor
		WorkerDataHeader() : TObject() {
			jobId= "";
			IdX= -1;
			IdY= -1;
			gettimeofday(&timestamp, NULL);
		};
		
		//Param constructor
		WorkerDataHeader(std::string id,long int ix,long int iy,timeval* tv=0) : TObject(){
			jobId= id;
			IdX= ix;
			IdY= iy;
			if(tv) {
				timestamp= *tv;	
			}
			else{
				gettimeofday(&timestamp, NULL);	
			}
		}
		
		//Copy constructor
		WorkerDataHeader(const WorkerDataHeader& info) : TObject() {
 			((WorkerDataHeader&)info).Copy(*this);
		}
		
		virtual ~WorkerDataHeader(){};

	public:	
		// Operator =
		WorkerDataHeader& operator=(const WorkerDataHeader& info) { 
			if (this != &info) ((WorkerDataHeader&)info).Copy(*this);
  		return *this;
		}

		//Copy method
		void Copy(TObject &obj) const {
			TObject::Copy((WorkerDataHeader&)obj);
  		((WorkerDataHeader&)obj).jobId = jobId;
			((WorkerDataHeader&)obj).IdX = IdX;	
			((WorkerDataHeader&)obj).IdY = IdY;	
		}

	public: 
		std::string jobId;
		long int IdX;
		long int IdY;
		timeval timestamp;

	ClassDef(WorkerDataHeader,1)
};
*/


class WorkerData : public TObject {
	public:
		//Standard constructor
		WorkerData();
				
		//Copy constructor
		WorkerData(const WorkerData& info);

		//Destructor
		virtual ~WorkerData();

		enum SourceDataType {
			eUNKNOWN_DATA= 0,
			eCOMPACT_SOURCE_DATA= 1,
			eEXTENDED_SOURCE_DATA= 2
		};

	public:
		// Operator =
		WorkerData& operator=(const WorkerData& data);

		//Copy method
		void Copy(TObject &obj) const;


	public: 
		
		WorkerTask info;
		int data_type;
		Source* source;
		std::vector<Source*> sources;
		Source* edge_source;
		std::vector<Source*> edge_sources;

	ClassDef(WorkerData,1)
};


#ifdef __MAKECINT__
//#pragma link C++ class WorkerDataHeader+;
#pragma link C++ class WorkerData+;
#endif

}//close namespace

#endif



