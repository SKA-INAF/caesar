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
* @file WorkerTask.cc
* @class WorkerTask
* @brief WorkerTask class
*
* Worker task class
* @author S. Riggi
* @date 20/01/2015
*/

#include <WorkerTask.h>
#include <Logger.h>
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
#include <sys/time.h>
#include <chrono>

using namespace std;

ClassImp(Caesar::WorkerTask)

namespace Caesar {

//Standard Constructor
WorkerTask::WorkerTask() : TObject() {	
	filename= "";
	jobId= "";
	worker_name= "";
	broker_name= "";
	IdX= -1;
	IdY= -1;
	ix_min= -1;
	ix_max= -1;
	iy_min= -1;
	iy_max= -1;
	gettimeofday(&timestamp, NULL);
}

//Copy constructor
WorkerTask::WorkerTask(const WorkerTask& info) : TObject() {
	((WorkerTask&)info).Copy(*this);
}

//Destructor
WorkerTask::~WorkerTask(){
	

}//close destructor

// Operator =
WorkerTask& WorkerTask::operator=(const WorkerTask& data) { 
	if (this != &data) ((WorkerTask&)data).Copy(*this);
  return *this;
}

//Copy
void WorkerTask::Copy(TObject &obj) const {
			
	TObject::Copy((WorkerTask&)obj);
	((WorkerTask&)obj).filename = filename;
	((WorkerTask&)obj).jobId = jobId;
	((WorkerTask&)obj).worker_name = worker_name;
	((WorkerTask&)obj).broker_name = broker_name;
	((WorkerTask&)obj).IdX = IdX;	
	((WorkerTask&)obj).IdY = IdY;	
	((WorkerTask&)obj).ix_min = ix_min;	
	((WorkerTask&)obj).ix_max = ix_max;		
	((WorkerTask&)obj).iy_min = iy_min;	
	((WorkerTask&)obj).iy_max = iy_max;
	((WorkerTask&)obj).timestamp = timestamp;		

}//close Copy()

int WorkerTask::SerializeToJson(Json::Value& root){

	try {
		root["filename"]= filename;
		root["jobId"]= jobId;	
		root["worker_name"]= worker_name;
		root["broker_name"]= broker_name;
		root["IdX"]= Json::Value( static_cast<Json::Int64>(IdX) );
		root["IdY"]= Json::Value( static_cast<Json::Int64>(IdY) );
		root["ix_min"]= Json::Value( static_cast<Json::Int64>(ix_min) );
		root["ix_max"]= Json::Value( static_cast<Json::Int64>(ix_max) );	
		root["iy_min"]= Json::Value( static_cast<Json::Int64>(iy_min) );
		root["iy_max"]= Json::Value( static_cast<Json::Int64>(iy_max) );
		root["t_sec"]= Json::Value( static_cast<Json::Int64>(timestamp.tv_sec) );
		root["t_nsec"]= Json::Value( static_cast<Json::Int64>(timestamp.tv_usec*1000) );
	}
	catch(const std::exception& e){
		ERROR_LOG("Failed to encode to json object (err="<<e.what()<<")");
		return -1;
	}

	return 0;

}//close SerializeToJson()

int WorkerTask::SerializeToJsonString(std::string& jsonString,bool isMinified){

	//Encode first to json object
	Json::Value jsonObj;
	if(SerializeToJson(jsonObj)<0) {
		ERROR_LOG("Failed to encode to json object!");
		return -1;
	}

	//Encode to string
	if(CodeUtils::JsonToString(jsonString,jsonObj,isMinified)<0){
		ERROR_LOG("Failed to encode json object to string!");
		return -1;
	}

	return 0;

}//close SerializeToJsonString()


int WorkerTask::SerializeToProtobuf(SourcePB::WorkerTask& taskPB){

	SourcePB::Timestamp* timestampPB= new SourcePB::Timestamp;
	timestampPB->set_seconds( timestamp.tv_sec );
	timestampPB->set_nanos( timestamp.tv_usec*1000 );

	taskPB.set_filename(filename);
	taskPB.set_worker_name(worker_name);	
	taskPB.set_broker_name(broker_name);
	taskPB.set_jobid(jobId);
	taskPB.set_idx(IdX);
	taskPB.set_idy(IdY);
	taskPB.set_ix_min(ix_min);
	taskPB.set_ix_max(ix_max);
	taskPB.set_iy_min(iy_min);
	taskPB.set_iy_max(iy_max);
	taskPB.set_allocated_timestamp(timestampPB);
	
	return 0;

}//close SerializeToProtobuf()


int WorkerTask::EncodeFromJson(const Json::Value& root){

	//Check input json object
	if(root.isNull() || root.empty()){
		ERROR_LOG("Null/empty json, failed to encode!");
		return -1;
	}

	//Get fields
	//--> filename
	Json::Value filenameObj= root["filename"];
	if(filenameObj.isNull() || filenameObj.isNull()){
		ERROR_LOG("Null/empty filename field, failed to encode!");
		return -1;
	}
	filename= filenameObj.asString();
	if(filename==""){
		ERROR_LOG("Empty string filename field, failed to encode!");
		return -1;
	}

	//--> jobId
	Json::Value jobIdObj= root["jobId"];
	if(jobIdObj.isNull() || jobIdObj.isNull()){
		ERROR_LOG("Null/empty jobId field, failed to encode!");
		return -1;
	}
	jobId= jobIdObj.asString();
	if(jobId==""){
		ERROR_LOG("Empty string jobId field, failed to encode!");
		return -1;
	}

	//--> worker_name
	Json::Value workerNameObj= root["worker_name"];
	if(workerNameObj.isNull() || workerNameObj.isNull()){
		ERROR_LOG("Null/empty worker_name field, failed to encode!");
		return -1;
	}
	worker_name= workerNameObj.asString();
	if(worker_name==""){
		ERROR_LOG("Empty string worker_name field, failed to encode!");
		return -1;
	}
		
	//--> broker_name
	Json::Value brokerNameObj= root["broker_name"];
	if(brokerNameObj.isNull() || brokerNameObj.isNull()){
		ERROR_LOG("Null/empty broker_name field, failed to encode!");
		return -1;
	}
	broker_name= brokerNameObj.asString();
	if(broker_name==""){
		ERROR_LOG("Empty string broker_name field, failed to encode!");
		return -1;
	}

	//--> IdX
	Json::Value IdXObj= root["IdX"];
	if(IdXObj.isNull() || IdXObj.isNull()){
		ERROR_LOG("Null/empty broker_name field, failed to encode!");
		return -1;
	}
	IdX= IdXObj.asInt64();
	if(IdX<0){
		ERROR_LOG("Invalid IdX field, failed to encode!");
		return -1;
	}
	
	//--> IdY
	Json::Value IdYObj= root["IdY"];
	if(IdYObj.isNull() || IdYObj.isNull()){
		ERROR_LOG("Null/empty broker_name field, failed to encode!");
		return -1;
	}
	IdY= IdYObj.asInt64();
	if(IdY<0){
		ERROR_LOG("Invalid IdY field, failed to encode!");
		return -1;
	}

	//--> ix_min
	Json::Value ix_minObj= root["ix_min"];
	if(ix_minObj.isNull() || ix_minObj.isNull()){
		ERROR_LOG("Null/empty ix_min field, failed to encode!");
		return -1;
	}
	ix_min= ix_minObj.asInt64();
	if(ix_min<0 && ix_min!=-1){
		ERROR_LOG("Invalid ix_min field, failed to encode!");
		return -1;
	}
	

	//--> ix_max
	Json::Value ix_maxObj= root["ix_max"];
	if(ix_maxObj.isNull() || ix_maxObj.isNull()){
		ERROR_LOG("Null/empty ix_max field, failed to encode!");
		return -1;
	}
	ix_max= ix_maxObj.asInt64();
	if(ix_max<0 && ix_max!=-1){
		ERROR_LOG("Invalid ix_max field, failed to encode!");
		return -1;
	}

	//--> iy_min
	Json::Value iy_minObj= root["iy_min"];
	if(iy_minObj.isNull() || iy_minObj.isNull()){
		ERROR_LOG("Null/empty iy_min field, failed to encode!");
		return -1;
	}
	iy_min= iy_minObj.asInt64();
	if(iy_min<0 && iy_min!=-1){
		ERROR_LOG("Invalid iy_min field, failed to encode!");
		return -1;
	}
	

	//--> iy_max
	Json::Value iy_maxObj= root["iy_max"];
	if(iy_maxObj.isNull() || iy_maxObj.isNull()){
		ERROR_LOG("Null/empty iy_max field, failed to encode!");
		return -1;
	}
	iy_max= iy_maxObj.asInt64();
	if(iy_max<0 && iy_max!=-1){
		ERROR_LOG("Invalid iy_max field, failed to encode!");
		return -1;
	}
	

	//--> timesec
	Json::Value TimeSecObj= root["t_sec"];
	if(TimeSecObj.isNull() || TimeSecObj.isNull()){
		ERROR_LOG("Null/empty t_sec field, failed to encode!");
		return -1;
	}
	timestamp.tv_sec= TimeSecObj.asInt64();
	if(timestamp.tv_sec<0){
		ERROR_LOG("Invalid timestamp.tv_sec field, failed to encode!");
		return -1;
	}

	Json::Value TimeNSecObj= root["t_nsec"];
	if(TimeNSecObj.isNull() || TimeNSecObj.isNull()){
		ERROR_LOG("Null/empty t_nsec field, failed to encode!");
		return -1;
	}
	timestamp.tv_usec= TimeNSecObj.asInt64()*1.e-3;
	if(timestamp.tv_usec<0){
		ERROR_LOG("Invalid timestamp.tv_usec field, failed to encode!");
		return -1;
	}
		

	return 0;

}//close EncodeFromJson()

int WorkerTask::EncodeFromJsonString(std::string& jsonString){

	//Parse string to json object
	Json::Value root;
	if(CodeUtils::StringToJson(root,jsonString)<0){
		ERROR_LOG("Failed to parse input string to json!");
		return -1;
	}
	
	//Encode from json object
	if(EncodeFromJson(root)<0){
		ERROR_LOG("Failed to encode from json!");
		return -1;
	}

	return 0;

}//close EncodeFromJsonString()


int WorkerTask::EncodeFromProtobuf(const SourcePB::WorkerTask& taskPB){

	//Fill fields
	//--> filename
	if(!taskPB.has_filename()){
		ERROR_LOG("Missing filename field, failed to encode!");
		return -1;
	}
	filename= taskPB.filename();
	if(filename==""){
		ERROR_LOG("Empty string filename field, failed to encode!");
		return -1;
	}

	//--> jobId
	if(!taskPB.has_jobid()){
		ERROR_LOG("Missing filename field, failed to encode!");
		return -1;
	}
	jobId= taskPB.jobid();
	if(jobId==""){
		ERROR_LOG("Empty string jobId field, failed to encode!");
		return -1;
	}

	//--> worker_name
	if(!taskPB.has_worker_name()){
		ERROR_LOG("Missing worker_name field, failed to encode!");
		return -1;
	}
	worker_name= taskPB.worker_name();
	if(worker_name==""){
		ERROR_LOG("Empty string worker_name field, failed to encode!");
		return -1;
	}

	//--> broker_name
	if(!taskPB.has_broker_name()){
		ERROR_LOG("Missing broker_name field, failed to encode!");
		return -1;
	}
	broker_name= taskPB.broker_name();
	if(broker_name==""){
		ERROR_LOG("Empty string broker_name field, failed to encode!");
		return -1;
	}

	//--> IdX
	if(!taskPB.has_idx()){
		ERROR_LOG("Missing IdX field, failed to encode!");
		return -1;
	}
	IdX= taskPB.idx();
	if(IdX<0){
		ERROR_LOG("Invalid IdX field, failed to encode!");
		return -1;
	}

	//--> IdY
	if(!taskPB.has_idy()){
		ERROR_LOG("Missing IdY field, failed to encode!");
		return -1;
	}
	IdY= taskPB.idy();
	if(IdY<0){
		ERROR_LOG("Invalid IdY field, failed to encode!");
		return -1;
	}

	//--> ix_min
	if(!taskPB.has_ix_min()){
		ERROR_LOG("Missing ix_min field, failed to encode!");
		return -1;
	}
	ix_min= taskPB.ix_min();
	if(ix_min<0 && ix_min!=-1){
		ERROR_LOG("Invalid ix_min field, failed to encode!");
		return -1;
	}

	//--> ix_max
	if(!taskPB.has_ix_max()){
		ERROR_LOG("Missing ix_max field, failed to encode!");
		return -1;
	}
	ix_max= taskPB.ix_max();
	if(ix_max<0 && ix_max!=-1){
		ERROR_LOG("Invalid ix_max field, failed to encode!");
		return -1;
	}

	//--> iy_min
	if(!taskPB.has_iy_min()){
		ERROR_LOG("Missing iy_min field, failed to encode!");
		return -1;
	}
	iy_min= taskPB.iy_min();
	if(iy_min<0 && iy_min!=-1){
		ERROR_LOG("Invalid iy_min field, failed to encode!");
		return -1;
	}

	//--> iy_max
	if(!taskPB.has_iy_max()){
		ERROR_LOG("Missing iy_max field, failed to encode!");
		return -1;
	}
	iy_max= taskPB.iy_max();
	if(iy_max<0 && iy_max!=-1){
		ERROR_LOG("Invalid iy_max field, failed to encode!");
		return -1;
	}


	//--> timestamp
	if(!taskPB.has_timestamp()){
		ERROR_LOG("Missing timestamp field, failed to encode!");
		return -1;
	}
	const SourcePB::Timestamp& timestampPB= taskPB.timestamp();
	timestamp.tv_sec= timestampPB.seconds();
	timestamp.tv_usec= timestampPB.nanos()*1.e-3;

	if(timestamp.tv_sec<0){
		ERROR_LOG("Invalid timestamp.tv_sec field, failed to encode!");
		return -1;
	}
	if(timestamp.tv_usec<0){
		ERROR_LOG("Invalid timestamp.tv_usec field, failed to encode!");
		return -1;
	}

	return 0;

}//close EncodeFromProtobuf()

}//close namespace 
