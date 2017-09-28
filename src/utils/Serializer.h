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
* @file Serializer.h
* @class Serializer
* @brief Serializer class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef _SERIALIZER_H
#define _SERIALIZER_H

#include <Source.h>
#include <SourceFitter.h>
#include <TaskData.h>
#include <Contour.h>

#include <Source.pb.h>
#include <TaskData.pb.h>


#include <TVector2.h> 

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>

namespace Caesar {

class SBuffer : public TObject {
	public:
		SBuffer(){};
		virtual ~SBuffer(){};
	public: 
		std::string data;
		long int size;

	ClassDef(SBuffer,1)
};

class Serializer : public TObject {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Serializer();
		/**
		* \brief Class destructor: free allocated memory
		*/
		~Serializer();

	public: 
		
			//## SOURCE SERIALIZATION
			//Source --> Buffer
			static int EncodePointToProtobuf(CaesarPB::Point& point_pb,TVector2& point);
			static int EncodeContourToProtobuf(CaesarPB::Contour& contour_pb,Contour* contour);
			static int EncodePixelToProtobuf(CaesarPB::Pixel& pixel_pb,Pixel* pixel);
			static int EncodeBlobToProtobuf(CaesarPB::Blob& blob_pb,Source* source);
			static int EncodeSourceToProtobuf(CaesarPB::Source& source_pb,Source* source);		
			static int SourceToBuffer(SBuffer& buffer,Source* source);
			static int EncodeSourceComponentParsToProtobuf(CaesarPB::SourceComponentPars* sourceCompPars_pb,SourceComponentPars& sourceCompPars);
			static int EncodeSourceFitParsToProtobuf(CaesarPB::SourceFitPars& sourceFitPars_pb,SourceFitPars& sourceFitPars);

			//Buffer --> Source
			static int EncodeProtobufToPoint(TVector2& point,const CaesarPB::Point& point_pb);
			static int EncodeProtobufToContour(Contour& contour,const CaesarPB::Contour& contour_pb);
			static int EncodeProtobufToPixel(Pixel& pixel,const CaesarPB::Pixel& pixel_pb);	
			static int EncodeProtobufToBlob(Source& source,const CaesarPB::Blob& blob_pb);
			static int EncodeProtobufToSource(Source& source,const CaesarPB::Source& source_pb);			
			static int BufferToSource(Source& source,SBuffer& buffer);

			static int EncodeProtobufToSourceComponentPars(SourceComponentPars& sourceComponentPars,const CaesarPB::SourceComponentPars& sourceComponentPars_pb);
			static int EncodeProtobufToSourceFitPars(SourceFitPars& sourceFitPars,CaesarPB::SourceFitPars& sourceFitPars_pb);

			//## WORKER DATA SERIALIZATION ###
			//TaskData --> Buffer
			static int EncodeTaskDataToProtobuf(CaesarPB::TaskData& taskData_pb,TaskData* taskData);
			static int TaskDataToBuffer(SBuffer& buffer,TaskData* taskData);
			static char* TaskDataToCharArray(long int& buffer_size,TaskData* taskData);

			static int EncodeTaskDataCollectionToProtobuf(CaesarPB::TaskDataCollection& taskDataCollection_pb,std::vector<TaskData*> taskDataCollection);
			static int TaskDataCollectionToBuffer(SBuffer& buffer,std::vector<TaskData*> taskDataCollection);
			static char* TaskDataCollectionToCharArray(long int& buffer_size,std::vector<TaskData*> taskDataCollection);

			//Buffer --> TaskData
			static int EncodeProtobufToTaskData(TaskData& taskData,const CaesarPB::TaskData& taskData_pb);
			static int BufferToTaskData(TaskData& taskData,SBuffer& buffer);
			static int CharArrayToTaskData(TaskData& taskData,char* buffer,long int buffer_size);

			static int EncodeProtobufToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,const CaesarPB::TaskDataCollection& taskDataCollection_pb,bool isTaskCollectionPreAllocated=false);
			static int BufferToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,SBuffer& buffer,bool isTaskCollectionPreAllocated=false);
			static int CharArrayToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,char* buffer,long int buffer_size,bool isTaskCollectionPreAllocated=false);

	public:
    
		//ClassDef(Serializer,1)

};

/*
#ifdef __MAKECINT__
#pragma link C++ class Serializer+;
#endif
*/

}//close namespace

#endif



