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
* @file WorkerManager.h
* @class WorkerManager
* @brief WorkerManager
*
* WorkerManager
* @author S. Riggi
* @date 16/06/2016
*/
#ifndef WorkerManager_H
#define WorkerManager_H

#include <Logger.h>

#include <tango.h>


#include <chrono>
#include <thread>
#include <mutex>

namespace Caesar {

class EvtSubscriptionData {

	public:
		EvtSubscriptionData(std::string devicename,std::string attrname,Tango::EventType eventtype) 
			: device_name(devicename), attr_name(attrname),event_type(eventtype) 
		{
			is_subscribed= false;
			subscription_id= -1;
		};
		~EvtSubscriptionData(){};

	public:
		std::string device_name;
		bool is_subscribed;
		int subscription_id;
		std::string attr_name;
		Tango::EventType event_type; 

};//close EvtSubscriptionData()

struct MatchSubscription {
	MatchSubscription(const EvtSubscriptionData& subscriptionData) 
		: m_subscriptionData(subscriptionData) {}
 	bool operator()(const EvtSubscriptionData* obj) const {
		bool areEqual= ( obj->device_name==m_subscriptionData.device_name && 
										 obj->attr_name==m_subscriptionData.attr_name &&
										 obj->event_type==m_subscriptionData.event_type);
  	return areEqual;
 	}
 	private:
  	const EvtSubscriptionData& m_subscriptionData;
};


typedef std::vector<EvtSubscriptionData*> Subscriptions;


class WorkerManager {

	public :
		//Constructor
		WorkerManager();
		
		//Destructor
		~WorkerManager();

	public:
		//Has worker
		bool HasWorker(const std::string& device_name);

		//Add worker to the list
		int AddWorker(const std::string& device_name);
		
		//Flush workers
		int FlushWorkers();

		//Get free worker group
		int GetFreeWorkers(Tango::Group& group,int nMaxWorkers=-1,bool requireExactly=false);		

		//Get worker names
		void GetWorkerNames(std::vector<std::string>& worker_names);
		void GetFreeWorkerNames(std::vector<std::string>& worker_names);
		void GetBusyWorkerNames(std::vector<std::string>& worker_names);

		//Subscribe to event
		int SubscribeGroupToEvent(std::string attr_name,Tango::EventType event_type,Tango::CallBack* event_cb);

		//Find subscription
		EvtSubscriptionData* FindSubscription(int& index,std::string& device_name,std::string& attr_name,Tango::EventType& event_type);
	
		//Update worker status
		int UpdateWorkerState(std::string& device_name,Tango::DevState& state);

	private:
		void Init();

	private:
		mutable std::mutex m_mutex;
		Tango::Group* m_workers;
		Tango::Group* m_freeWorkers;
		Tango::Group* m_busyWorkers;
		EvtSubscriptionData* m_subscriptionData;
		Subscriptions m_subscriptions;

};//close WorkerManager

}//close namespace

#endif
