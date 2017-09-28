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
* @file WorkerManager.cc
* @class WorkerManager
* @brief WorkerManager class
*
* WorkerManager class
* @author S. Riggi
* @date 20/01/2015
*/

#include <WorkerManager.h>
#include <Logger.h>

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


namespace Caesar {

//Standard Constructor
WorkerManager::WorkerManager() {	
	Init();
}


//Destructor
WorkerManager::~WorkerManager(){
	
	if(m_freeWorkers){
		delete m_freeWorkers;
		m_freeWorkers= 0;
	}
	if(m_busyWorkers){
		delete m_busyWorkers;
		m_busyWorkers= 0;
	}
	if(m_workers){
		delete m_workers;
		m_workers= 0;
	}

	for(unsigned int i=0;i<m_subscriptions.size();i++){
		if(m_subscriptions[i]){
			delete m_subscriptions[i];
			m_subscriptions[i]= 0;
		}
	}
	m_subscriptions.clear();

}//close destructor

void WorkerManager::Init(){

	m_workers= 0;
	m_workers= new Tango::Group("Workers");

	m_freeWorkers= 0;
	m_freeWorkers= new Tango::Group("FreeWorkers");

	m_busyWorkers= 0;
	m_busyWorkers= new Tango::Group("BusyWorkers");
	
	//m_workers->add(m_freeWorkers);	
	//m_workers->add(m_busyWorkers);	

	m_subscriptionData= 0;
	m_subscriptions.clear();

}//close Init()


bool WorkerManager::HasWorker(const std::string& device_name){
	
	bool hasWorker= false;
	try{
		hasWorker= m_workers->contains(device_name);
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Exception occurred while checking if worker ("<<device_name<<") exists in Group returning false!");
		hasWorker= false;
	}
	return hasWorker;

}//close HasWorker()

int WorkerManager::AddWorker(const std::string& device_name){

	std::lock_guard<std::mutex> lock(m_mutex);

	//Check given device name
	if(device_name=="") return -1;
	
	//Check if device already exists
	DEBUG_LOG("Check if device "<<device_name<<" already exist in the list...");
	if(HasWorker(device_name)){
		WARN_LOG("Worker ("<<device_name<<") already added!");
		return 0;
	}

	//Add worker in Group (add in busy group initially)
	DEBUG_LOG("Adding device "<<device_name<<" to worker group...");
	try{
		m_workers->add(device_name,20000);
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Exception occurred while adding worker "<<device_name<<" to Group...");
		return -1;
	}
	
	return 0;

}//close AddWorker()


int WorkerManager::FlushWorkers(){

	std::lock_guard<std::mutex> lock(m_mutex);

	try{
		m_workers->remove_all();
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Exception occurred while removing all workers in Group...");
		return -1;
	}

	return 0;

}//close FlushWorkers()


int WorkerManager::GetFreeWorkers(Tango::Group& group,int nMaxWorkers,bool requireExactly){
	
	std::lock_guard<std::mutex> lock(m_mutex);

	//Check for free workers
	//int nFreeWorkers= static_cast<int>(m_workers->get_group("FreeWorkers")->get_size());
	int nFreeWorkers= static_cast<int>(m_freeWorkers->get_size());
	if(nFreeWorkers<=0) {
		INFO_LOG("No free workers available!");
		return -1;
	}
	
	//Return all free worker group in this case
	//std::vector<std::string> freeWorkersNames= m_workers->get_group("FreeWorkers")->get_device_list();
	std::vector<std::string> freeWorkersNames= m_freeWorkers->get_device_list();
	if(freeWorkersNames.empty()){
		WARN_LOG("No free workers available!");
		return -1;
	}

	if(nMaxWorkers==-1){
		try {		
			group.add(freeWorkersNames);
		}
		catch(Tango::DevFailed& e){
			ERROR_LOG("Failed to add all workers in group!");
			return -1;
		}
		catch(const std::exception& e){
			ERROR_LOG("C++ exception when adding all workers in group (err="<<e.what()<<")");
			return -1;
		}
		return 0;
	}//close if
	
	
	if(nMaxWorkers>nFreeWorkers && requireExactly){
		WARN_LOG("Number of requested workers exceed the available free workers!");
		return -1;
	}

	//Return group
	try {	
		for(int i=0;i<freeWorkersNames.size();i++){
			group.add(freeWorkersNames[i]);
		}
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Failed to add selected workers in group!");
		return -1;
	}
	catch(const std::exception& e){
		ERROR_LOG("C++ exception when adding selected workers in group (err="<<e.what()<<")");
		return -1;
	}

	return 0;

}//close GetFreeWorkers()

//Get worker names
void WorkerManager::GetWorkerNames(std::vector<std::string>& worker_names){
			
	std::lock_guard<std::mutex> lock(m_mutex);
	worker_names.clear();
	worker_names= m_workers->get_device_list();

}//close GetWOrkerNames()

void WorkerManager::GetFreeWorkerNames(std::vector<std::string>& worker_names){
			
	std::lock_guard<std::mutex> lock(m_mutex);
	worker_names.clear();
	//worker_names= m_workers->get_group("FreeWorkers")->get_device_list();
	worker_names= m_freeWorkers->get_device_list();

}//close GetFreeWOrkerNames()

void WorkerManager::GetBusyWorkerNames(std::vector<std::string>& worker_names){
			
	std::lock_guard<std::mutex> lock(m_mutex);
	worker_names.clear();
	//worker_names= m_workers->get_group("BusyWorkers")->get_device_list();
	worker_names= m_busyWorkers->get_device_list();

}//close GetFreeWOrkerNames()

int WorkerManager::SubscribeGroupToEvent(std::string attr_name,Tango::EventType event_type,Tango::CallBack* event_cb){

	//std::lock_guard<std::mutex> lock(m_mutex);

	//Check callback
	if(!event_cb){
		ERROR_LOG("Null ptr to given callback!");
		return -1;
	}

	//Check group size
	if(m_workers->get_size()<=0){
		WARN_LOG("No workers present in group...");
		return -1;
	}
	

	//Subscribe to event and store subscriptions made
	bool hasFailures= true;
	std::vector<std::string> device_names;
	try {
		device_names= m_workers->get_device_list();
	}
	catch(const Tango::DevFailed& e){
		ERROR_LOG("Cannot get device names, failed to subscribe to events!");
		return -1;
	}

	for(unsigned int i=0;i<device_names.size();i++){
		DEBUG_LOG("Retrieving worker "<<device_names[i]<<" from Group...");
		
		Tango::DeviceProxy* device= 0;
		std::string device_name= device_names[i];
		try {
			device= m_workers->get_device(device_name);
			if(!device) throw std::runtime_error("Null ptr to retrieved device proxy!");
			device_name= device->name();
		}
		catch(const Tango::DevFailed& e){
			ERROR_LOG("Cannot get device name, failed to subscribe worker no. "<<i+1<<" in group to requested event!");
			continue;
		}
		catch(const std::exception& e){
			ERROR_LOG("Exception ("<<e.what()<<"), failed to subscribe worker no. "<<i+1<<" in group to requested event!");
			continue;
		}
		catch(...){
			ERROR_LOG("Unkwown exception, failed to subscribe worker no. "<<i+1<<" in group to requested event!");
			continue;
		}

		//Check if worker is already subscribed
		int index= -1;
		bool hasSubscription= false;
		DEBUG_LOG("Check if worker "<<device_name<<" is already subscribed...");
		EvtSubscriptionData* subscriptionData= FindSubscription(index,device_name,attr_name,event_type);

		if(subscriptionData){
			hasSubscription= true;
			DEBUG_LOG("Subscription to event "<<attr_name<<" ("<<event_type<<") for device "<<device_name<<" already present, checking if already subscribed...");
			if(subscriptionData->is_subscribed){
				DEBUG_LOG("Subscription to event "<<attr_name<<" ("<<event_type<<") for device "<<device_name<<" already done, go to next device...");
				continue;
			}
		}
		else{
			DEBUG_LOG("Subscription to event "<<attr_name<<" ("<<event_type<<") for device "<<device_name<<" not found, adding to the list...");
		}

		//Subscribe to event
		int subscriptionId= -1;
		try {
			DEBUG_LOG("Subscribing to event...");
			subscriptionId= device->subscribe_event(attr_name,event_type,event_cb);
			DEBUG_LOG("done!");
		}
		catch(const Tango::DevFailed& e){
			ERROR_LOG("Failed to subscribe worker no. "<<i+1<<" in group to requested event!");
			continue;
		}

		//Add subscription to the list
		m_mutex.lock();
		if(hasSubscription){	
			DEBUG_LOG("Update subscription data...");
			m_subscriptions[index]->subscription_id = subscriptionId;
			m_subscriptions[index]->is_subscribed= true;		
		}
		else{
			DEBUG_LOG("Add subscription to the list...");
			m_subscriptionData= new EvtSubscriptionData(device_name,attr_name,event_type);
			m_subscriptionData->subscription_id = subscriptionId;
			m_subscriptionData->is_subscribed= true;
			m_subscriptions.push_back(m_subscriptionData);
		}
		m_mutex.unlock();

		DEBUG_LOG("Go to next worker in group...");
	}//end loop devices in group

	return 0;

}//close SubscribeGroupToEvent()


EvtSubscriptionData* WorkerManager::FindSubscription(int& index,std::string& device_name,std::string& attr_name,Tango::EventType& event_type){

	//Init		
	index= -1;

	//Check for empty collection
	int nSubscriptions= static_cast<int>(m_subscriptions.size());
	if(nSubscriptions<=0) return 0;
			
	//Find item
	EvtSubscriptionData subscriptionData(device_name,attr_name,event_type);
	std::vector<EvtSubscriptionData*>::iterator it = std::find_if(m_subscriptions.begin(),m_subscriptions.end(), MatchSubscription(subscriptionData));
	if (it==m_subscriptions.end()) return 0;//not found in collection
			
	size_t pos = it-m_subscriptions.begin();
	index= pos;			
	return m_subscriptions[index];

}//close FindSubscription()


int WorkerManager::UpdateWorkerState(std::string& device_name,Tango::DevState& state){

	//std::lock_guard<std::mutex> lock(m_mutex);

	//Check
	if(device_name=="") return -1;
	if(state!=Tango::RUNNING && state!=Tango::INIT && state!=Tango::ON){
		WARN_LOG("Invalid state ("<<state<<") given!");
		return -1;	
	}

	//Find device in group
	DEBUG_LOG("Check if device "<<device_name<<" is in the Group...");
	if(!m_workers->contains(device_name)){
		WARN_LOG("Given device ("<<device_name<<" does not exists in the worker list!");
		return -1;
	}

	//Check if a change of group is required
	DEBUG_LOG("Check in which group the device is found...");
	//bool isInBusyGroup= m_workers->get_group("BusyWorkers")->contains(device_name);
	//bool isInFreeGroup= m_workers->get_group("FreeWorkers")->contains(device_name);
	bool isInBusyGroup= m_busyWorkers->contains(device_name);
	bool isInFreeGroup= m_freeWorkers->contains(device_name);
	/*
	if(!isInBusyGroup && !isInFreeGroup){
		WARN_LOG("Given device ("<<device_name<<" is present but not added to busy/free worker list!");
		return -1;
	} 
	*/

	DEBUG_LOG("Check if a change of group is required...");
	if( (isInBusyGroup && state==Tango::RUNNING) ||  
			(isInFreeGroup && state==Tango::ON) ||
			(isInBusyGroup && state==Tango::INIT)
	){
		INFO_LOG("No group change needed for worker device ("<<device_name<<")!");
		return 0;
	}

	//Change group
	if(state==Tango::ON){
		
		//Remove from busy group	
		if(isInBusyGroup){
			DEBUG_LOG("Removing device "<<device_name<<" from busy group...");
			m_mutex.lock();
			try{
				//m_workers->get_group("BusyWorkers")->remove(device_name);
				m_busyWorkers->remove(device_name);
			}
			catch(const Tango::DevFailed& e){
				m_mutex.unlock();
				ERROR_LOG("Failed to remove worker device ("<<device_name<<") from busy group!");
				return -1;
			}
			m_mutex.unlock();
		}//close if

		//Add to free group
		if(!isInFreeGroup){
			DEBUG_LOG("Adding device "<<device_name<<" too free group...");
			m_mutex.lock();
			try{
				//m_workers->get_group("FreeWorkers")->add(device_name,20000);
				m_freeWorkers->add(device_name,20000);
			}
			catch(const Tango::DevFailed& e){
				m_mutex.unlock();
				ERROR_LOG("Failed to add worker device ("<<device_name<<") to free group!");
				return -1;
			}
			m_mutex.unlock();
		}//close if

	}//close if
	
	if(state==Tango::RUNNING || state==Tango::INIT){
		
		//Remove from free group
		if(isInFreeGroup){
			DEBUG_LOG("Removing device "<<device_name<<" from free group...");
			m_mutex.lock();
			try{
				//m_workers->get_group("FreeWorkers")->remove(device_name);
				m_freeWorkers->remove(device_name);
			}
			catch(const Tango::DevFailed& e){
				m_mutex.unlock();
				ERROR_LOG("Failed to remove worker device ("<<device_name<<") from free group!");
				return -1;
			}
			m_mutex.unlock();
		}//close if

		//Add to busy group
		if(!isInBusyGroup){
			DEBUG_LOG("Adding device "<<device_name<<" to busy group...");
			m_mutex.lock();
			try{
				//m_workers->get_group("BusyWorkers")->add(device_name,20000);
				m_busyWorkers->add(device_name,20000);
			}
			catch(const Tango::DevFailed& e){
				m_mutex.unlock();
				ERROR_LOG("Failed to add worker device ("<<device_name<<") to busy group!");
				return -1;
			}			
			m_mutex.unlock();
		}//close if
	}//close if
	

	return 0;

}//close UpdateWorkerState()


}//close namespace
