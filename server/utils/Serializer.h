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
#ifndef Serializer_H
#define Serializer_H

#include <Source.pb.h>
#include <Source.h>
#include <Contour.h>
#include <WorkerData.h>
#include <WorkerTask.h>

//Tango header
#include <tango.h>

//Msgpack
#include <msgpack.hpp>


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
			static int EncodePointToProtobuf(SourcePB::Point& point_pb,TVector2& point);
			static int EncodeContourToProtobuf(SourcePB::Contour& contour_pb,Contour* contour);
			static int EncodePixelToProtobuf(SourcePB::Pixel& pixel_pb,Pixel* pixel);
			static int EncodeBlobToProtobuf(SourcePB::Blob& blob_pb,Source* source);
			static int EncodeSourceToProtobuf(SourcePB::Source& source_pb,Source* source);		
			static int SourceToBuffer(SBuffer& buffer,Source* source);
			
			//Buffer --> Source
			static int EncodeProtobufToPoint(TVector2& point,const SourcePB::Point& point_pb);
			static int EncodeProtobufToContour(Contour& contour,const SourcePB::Contour& contour_pb);
			static int EncodeProtobufToPixel(Pixel& pixel,const SourcePB::Pixel& pixel_pb);	
			static int EncodeProtobufToBlob(Source& source,const SourcePB::Blob& blob_pb);
			static int EncodeProtobufToSource(Source& source,const SourcePB::Source& source_pb);			
			static int BufferToSource(Source& source,SBuffer& buffer);

			//Source --> Tango pipe
			static int SourceToTangoPipe(Tango::DevicePipeBlob& pipe_blob,Source* source);
			static int SourceCollectionToTangoPipe(Tango::DevicePipeBlob& pipe_blob,std::vector<Source*>& sources);

			//Tango pipe --> Source
			static int TangoPipeToSource(Source& source,Tango::DevicePipeBlob& pipe_blob);
			static int TangoPipeToSourceCollection(std::vector<Source*>& source_list,Tango::DevicePipeBlob& pipe_blob);
			static int MakeSourceDataPipe(Tango::DevicePipeBlob& pipe_blob,std::string jobId,long int IdX,long int IdY,std::vector<Source*>& sources,std::vector<Source*>& edge_sources);

			//## WORKER TASK SERIALIZATION #####
			//WorkerTask collection --> Json
			static int WorkerTasksToJson(Json::Value& root,std::vector<WorkerTask*>& tasks);
			static int WorkerTasksToJsonString(std::string& jsonString,std::vector<WorkerTask*>& tasks,bool isMinified=true);

			//Json --> WorkerTask collection
			static int JsonToWorkerTasks(std::vector<WorkerTask*>& tasks,Json::Value& root);
			static int JsonStringToWorkerTasks(std::vector<WorkerTask*>& tasks,std::string& jsonString);

			//## WORKER DATA SERIALIZATION ###
			//WorkerData --> Buffer
			static int EncodeWorkerDataToProtobuf(SourcePB::WorkerData& workerData_pb,WorkerData* workerData);
			static int WorkerDataToBuffer(SBuffer& buffer,WorkerData* workerData);
			static int WorkerDataToCharArray(unsigned char* buffer,long int& buffer_size,WorkerData* workerData);

			//Buffer --> WorkerData
			static int EncodeProtobufToWorkerData(WorkerData& workerData,const SourcePB::WorkerData& workerData_pb);
			static int BufferToWorkerData(WorkerData& workerData,SBuffer& buffer);


			

			//Msgpack
			static int SourceToBuffer(Source* source,msgpack::sbuffer& buffer);
			static int SourceToString(Source* source,std::string& msg);
			static int SourceToDevString(Source* source,Tango::DevString& msg);
			
			static int SourceCollectionToBuffer(std::vector<Source*>& sources,msgpack::sbuffer& buffer);
			static int SourceCollectionToString(std::vector<Source*>& sources,std::string& msg);
			static int SourceCollectionToDevString(std::vector<Source*>& sources,Tango::DevString& msg);
		

		
	public:
    
		ClassDef(Serializer,1)


};

#ifdef __MAKECINT__
#pragma link C++ class Serializer+;
#endif

}//close namespace

#endif



