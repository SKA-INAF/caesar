
#include <WorkerSourceDataCallBack.h>
#include <SFinderBroker.h>
#include <WorkerManager.h>

#include <Logger.h>
#include <WorkerData.h>
#include <Serializer.h>

#include <tango.h>

//## Standard headers
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <signal.h>
#include <ctime>
#include <stdexcept>
#include <stdlib.h>

#include <ratio>
#include <chrono>
#include <map>
#include <vector>

using namespace Caesar;
using namespace std;


namespace SFinderBroker_ns {


SFinderBroker* WorkerSourceDataCallBack::device;


WorkerSourceDataCallBack::WorkerSourceDataCallBack (SFinderBroker* dev) {
   
	device = dev;
	
}//close constructor
	

WorkerSourceDataCallBack::~WorkerSourceDataCallBack(){

}//destructor


void WorkerSourceDataCallBack::push_event(Tango::EventData* event) { 
    
	//## Get access to event info
	/*
	[device] : The DeviceProxy object on which the call was executed (Tango::DeviceProxy *)
	[attr_name] : The attribute name (std::string &)
	[event] : The event name (std::string &)
	[attr_value] : The attribute data (DeviceAttribute *)
	[err] : A boolean flag set to true if the request failed. False otherwise (bool)
	[errors] : The error stack (DevErrorList &)
	*/

	if(!event){
		INFO_LOG("Received nullptr EventData!"); 
		return;
	}
	
	
  try {   	
    if (!event->err) { 
			//Get event data
			SBuffer buffer;	
			std::vector<unsigned char> event_value;

			//std::string event_value;
			
			Tango::DeviceAttribute attr_value= *(event->attr_value);
    	attr_value >> event_value; 
			//attr_value >> buffer.data;

			//buffer.data= std::string(event_value,422417);
			//buffer.size= (buffer.data).size(); 
			//buffer.size= 422417;//DEBUG
			
			buffer.data= std::string(event_value.begin(),event_value.end());
			buffer.size= event_value.size();

			int status= event->err;
			std::string dev_name= (event->device)->name();
			std::string attr_name= event->attr_name;
			std::string event_name= event->event;
			
			INFO_LOG("Worker sourceData event received: name="<< event->attr_name << " event " << event->event << " value="<<buffer.data<<" (size="<<buffer.size<<") (err= " << event->err << ")");

			//Parse sourceData
			WorkerData sourceData;
			if(Serializer::BufferToWorkerData(sourceData,buffer)<0){
				throw std::runtime_error("Failed to parse from buffer to WorkerData");
			}
			std::string worker_name= (sourceData.info).worker_name;
			std::string jobId= (sourceData.info).jobId;
			long int IdX= (sourceData.info).IdX;
			long int IdY= (sourceData.info).IdY;
			int dataType= sourceData.data_type;
			INFO_LOG("Job "<<jobId<<", worker="<<worker_name<<", task (IdX,IdY)=("<<IdX<<","<<IdY<<"), data_type="<<dataType<<", #sources="<<(sourceData.sources).size()<<", #edge_sources="<<(sourceData.edge_sources).size());

			//Update worker state in collection
			//if(device && (device->m_workerManager) && (device->m_workerManager)->UpdateWorkerState(dev_name,event_value)<0 ){
			//	ERROR_LOG("Failed to update worker state!");
			//}

    }//close if success event 
  }//close try 	
	catch(Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Failed to access to event data!");
		return;
	}	
	catch (const std::exception& e) { 
  	ERROR_LOG("C++ exception occurred while accessing to event data (err="<<e.what()<<")");
		return;
  }
  catch (...) { 
  	ERROR_LOG("Unknown exception occurred while accessing to event data");
		return;
  }
 
	return;
	
}//close EventCallBack::push_event()



}//close namespace


