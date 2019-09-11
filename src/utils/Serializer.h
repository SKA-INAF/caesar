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
#include <AstroObject.h>

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

class ImgMetaData;

class SBuffer : public TObject 
{
	public:
		/** 
		\brief SBuffer constructor
 		*/
		SBuffer(){};
		
		/** 
		\brief SBuffer destructor
 		*/
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
		/**
		* \brief Serialize image metadata to protobuf
		*/
		static int EncodeMetaDataToProtobuf(CaesarPB::ImgMetaData& metadata_pb,ImgMetaData* metadata);
		/**
		* \brief Serialize TVector2 point object to protobuf
		*/
		static int EncodePointToProtobuf(CaesarPB::Point& point_pb,TVector2& point);
		/**
		* \brief Serialize contour object to protobuf
		*/			
		static int EncodeContourToProtobuf(CaesarPB::Contour& contour_pb,Contour* contour);
		/**
		* \brief Serialize image pixel object to protobuf
		*/
		static int EncodePixelToProtobuf(CaesarPB::Pixel& pixel_pb,Pixel* pixel);
		/**
		* \brief Serialize blob object to protobuf
		*/
		static int EncodeBlobToProtobuf(CaesarPB::Blob& blob_pb,Source* source);
		/**
		* \brief Serialize source object to protobuf
		*/
		static int EncodeSourceToProtobuf(CaesarPB::Source& source_pb,Source* source);
		/**
		* \brief Serialize source collection to protobuf
		*/
		static int EncodeSourceCollectionToProtobuf(CaesarPB::SourceCollection& sources_pb,const std::vector<Source*>& sources);
		/**
		* \brief Serialize source object to SBuffer object
		*/
		static int SourceToBuffer(SBuffer& buffer,Source* source);
		/**
		* \brief Serialize source fit component pars to protobuf
		*/
		static int EncodeSourceComponentParsToProtobuf(CaesarPB::SourceComponentPars* sourceCompPars_pb,SourceComponentPars& sourceCompPars);
		/**
		* \brief Serialize source fit pars to protobuf
		*/	
		static int EncodeSourceFitParsToProtobuf(CaesarPB::SourceFitPars& sourceFitPars_pb,SourceFitPars& sourceFitPars);
		/**
		* \brief Serialize spectral index data to protobuf
		*/
		static int EncodeSpectralIndexDataToProtobuf(CaesarPB::SpectralIndexData& spectralIndexData_pb,const SpectralIndexData& spectralIndexData);
		/**
		* \brief Serialize spectral index collection to protobuf
		*/
		static int EncodeSpectralIndexDataCollectionToProtobuf(CaesarPB::SpectralIndexDataCollection& spectralIndexDataCollection_pb,const std::vector<SpectralIndexData>& spectralIndexDataCollection);
		/**
		* \brief Serialize astronomical object to protobuf
		*/
		static int EncodeAstroObjectToProtobuf(CaesarPB::AstroObject& astroObject_pb,const AstroObject& astroObject);
		/**
		* \brief Serialize astro object collection to protobuf
		*/
		static int EncodeAstroObjectCollectionToProtobuf(CaesarPB::AstroObjectCollection& astroObjectCollection_pb,const std::vector<AstroObject>& astroObjectCollection);

		/**
		* \brief Deserialize protobuf metadata object
		*/
		static int EncodeProtobufToMetaData(ImgMetaData& metadata,const CaesarPB::ImgMetaData& metadata_pb);
		/**
		* \brief Deserialize protobuf point object
		*/
		static int EncodeProtobufToPoint(TVector2& point,const CaesarPB::Point& point_pb);
		/**
		* \brief Deserialize protobuf contour object
		*/
		static int EncodeProtobufToContour(Contour& contour,const CaesarPB::Contour& contour_pb);
		/**
		* \brief Deserialize protobuf pixel object
		*/
		static int EncodeProtobufToPixel(Pixel& pixel,const CaesarPB::Pixel& pixel_pb);	
		/**
		* \brief Deserialize protobuf blob object
		*/
		static int EncodeProtobufToBlob(Source& source,const CaesarPB::Blob& blob_pb);
		/**
		* \brief Deserialize protobuf source object
		*/
		static int EncodeProtobufToSource(Source& source,const CaesarPB::Source& source_pb);			
		/**
		* \brief Convert SBuffer to Source object
		*/
		static int BufferToSource(Source& source,SBuffer& buffer);
		/**
		* \brief Deserialize protobuf source fit component pars object
		*/
		static int EncodeProtobufToSourceComponentPars(SourceComponentPars& sourceComponentPars,const CaesarPB::SourceComponentPars& sourceComponentPars_pb);
		/**
		* \brief Deserialize protobuf source fit pars object
		*/
		static int EncodeProtobufToSourceFitPars(SourceFitPars& sourceFitPars,const CaesarPB::SourceFitPars& sourceFitPars_pb);
		/**
		* \brief Deserialize protobuf spectral index data object
		*/
		static int EncodeProtobufToSpectralIndexData(SpectralIndexData& spectralIndexData,const CaesarPB::SpectralIndexData& spectralIndexData_pb);
		/**
		* \brief Deserialize protobuf spectral index data collection object
		*/
		static int EncodeProtobufToSpectralIndexDataCollection(std::vector<SpectralIndexData>& spectralIndexDataCollection,const CaesarPB::SpectralIndexDataCollection& spectralIndexDataCollection_pb);
		/**
		* \brief Deserialize protobuf astro object object
		*/
		static int EncodeProtobufToAstroObject(AstroObject& astroObject,const CaesarPB::AstroObject& astroObject_pb);
		/**
		* \brief Deserialize protobuf astro object collection object
		*/
		static int EncodeProtobufToAstroObjectCollection(std::vector<AstroObject>& astroObjectCollection,const CaesarPB::AstroObjectCollection& astroObjectCollection_pb);

		//## WORKER DATA SERIALIZATION ###
		/**
		* \brief Serialize task data object to protobuf
		*/
		static int EncodeTaskDataToProtobuf(CaesarPB::TaskData& taskData_pb,TaskData* taskData);
		/**
		* \brief Serialize task data to SBuffer
		*/
		static int TaskDataToBuffer(SBuffer& buffer,TaskData* taskData);
		/**
		* \brief Serialize task data object to char array
		*/
		static char* TaskDataToCharArray(long int& buffer_size,TaskData* taskData);
		/**
		* \brief Serialize task data collection to protobuf
		*/
		static int EncodeTaskDataCollectionToProtobuf(CaesarPB::TaskDataCollection& taskDataCollection_pb,const std::vector<TaskData*>& taskDataCollection);
		/**
		* \brief Convert task data collection to SBuffer
		*/
		static int TaskDataCollectionToBuffer(SBuffer& buffer,const std::vector<TaskData*>& taskDataCollection);
		/**
		* \brief Serialize task data collection to char array
		*/
		static char* TaskDataCollectionToCharArray(long int& buffer_size,const std::vector<TaskData*>& taskDataCollection);

			
	
		/**
		* \brief Deserialize task data protobuf object
		*/
		static int EncodeProtobufToTaskData(TaskData& taskData,const CaesarPB::TaskData& taskData_pb);
		/**
		* \brief Convert SBuffer to task data
		*/
		static int BufferToTaskData(TaskData& taskData,SBuffer& buffer);
		/**
		* \brief Convert char array to task data
		*/
		static int CharArrayToTaskData(TaskData& taskData,char* buffer,long int buffer_size);
		/**
		* \brief Deserialize task data collection protobuf object
		*/
		static int EncodeProtobufToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,const CaesarPB::TaskDataCollection& taskDataCollection_pb,bool isTaskCollectionPreAllocated=false);
		/**
		* \brief Convert SBuffer to task data collection
		*/
		static int BufferToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,SBuffer& buffer,bool isTaskCollectionPreAllocated=false);
		/**
		* \brief Convert char array to task data collection
		*/
		static int CharArrayToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,char* buffer,long int buffer_size,bool isTaskCollectionPreAllocated=false);

		/**
		* \brief Convert source collection to char array
		*/
		static char* SourceCollectionToCharArray(long int& buffer_size,const std::vector<Source*>& sourceCollection);

		/**
		* \brief Deserialize source collection protobuf object
		*/
		static int EncodeProtobufToSourceCollection(std::vector<Source*>& sources,const CaesarPB::SourceCollection& sources_pb,bool isSourceCollectionPreAllocated=false);
		/**
		* \brief Convert char array to source collection
		*/
		static int CharArrayToSourceCollection(std::vector<Source*>& sources,char* buffer,long int buffer_size,bool isSourceCollectionPreAllocated=false);


};//close Serializer class


}//close namespace

#endif


