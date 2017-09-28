
#include <WorkerStateCallBack.h>
#include <SFinderBroker.h>
#include <WorkerManager.h>

#include <Logger.h>

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


SFinderBroker* WorkerStateCallBack::device;


WorkerStateCallBack::WorkerStateCallBack (SFinderBroker* dev) {
   
	device = dev;
	
}//close constructor
	

WorkerStateCallBack::~WorkerStateCallBack(){

}//destructor


void WorkerStateCallBack::push_event(Tango::EventData* event) { 
    
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
			Tango::DevState event_value;
			Tango::DeviceAttribute attr_value= *(event->attr_value);
    	attr_value >> event_value; 
			
			int status= event->err;
			std::string dev_name= (event->device)->name();
			std::string attr_name= event->attr_name;
			std::string event_name= event->event;
			
			INFO_LOG("Worker State event received: name="<< event->attr_name << " event " << event->event << " value="<<event_value<<" (err= " << event->err << ")");

			//Update worker state in collection
			if(device && (device->m_workerManager) && (device->m_workerManager)->UpdateWorkerState(dev_name,event_value)<0 ){
				ERROR_LOG("Failed to update worker state!");
			}

    }//close if success event 
  }//close try 	
	catch(Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Failed to access to event data!");
		return;
	}	
  catch (...) { 
  	ERROR_LOG("Unknown exception occurred while accessing to event data");
		return;
  }
 
	return;
	
}//close EventCallBack::push_event()



}//close namespace


