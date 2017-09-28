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
* @file SFinderTask.h
* @class SFinderTask
* @brief SFinderTask class
*
* Class for source finder tasks
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef SFinderTask_H
#define SFinderTask_H

#include <SFinder.h>

//Tango headers
#include <tango.h>
#include <yat4tango/DeviceTask.h>
#include <yat/threading/Task.h>
#include <yat/threading/Mutex.h>


//Caesar headers
#include <Logger.h>
#include <Img.h>
#include <BkgData.h>
#include <Source.h>
#include <WorkerData.h>
#include <WorkerTask.h>
namespace Caesar {
	class Img;
	class BkgData;
	class Source;
	class WorkerData;
	class WorkerTask;
	class Logger;
	class ScopedLogger;
}

//C++ headers
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

using namespace Caesar;

namespace SFinder_ns {

class SFinder; 

class SFinderTask : public yat4tango::DeviceTask {
  
  public:
		//- Constructor
    SFinderTask (const yat::Task::Config& cfg,SFinder* dev);

    //- Destructor
    virtual ~SFinderTask ();

	public:
		//- Specialization of the exit behaviour
    virtual void exit () throw (Tango::DevFailed);

		//- start source finder (timeout_ms=0 means asynchronous exec)
    void Start(const std::vector<Caesar::WorkerTask*>& tasks,size_t timeout_ms = 0) throw (Tango::DevFailed);

		//- Stop source finder
		void Stop(size_t timeout_ms = 0) throw (Tango::DevFailed);
			

	protected:
    //- process_message (implements yat4tango::DeviceTask pure virtual method)
    virtual void process_message (yat::Message& msg) throw (Tango::DevFailed);
	
	private:
		/** 
		\brief Main thread function 
 		*/
		void Run(const std::vector<Caesar::WorkerTask*>& tasks);

		/** 
		\brief Run source finder task for a single image
 		*/
		int RunTask(Caesar::WorkerTask& task);

		/** 
		\brief Read image
 		*/
    Caesar::Img* ReadImage(const std::string& filename,long int tileMinX=-1,long int tileMaxX=-1,long int tileMinY=-1,long int tileMaxY=-1);

 		/** 
		\brief Compute stats & bkg
 		*/
		Caesar::BkgData* ComputeStatsAndBkg(Caesar::Img* img);

		/** 
		\brief Find compact sources
 		*/
		int FindCompactSources(Caesar::WorkerData& workerData,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);

		/** 
		\brief Find sources
 		*/
		int FindSources(std::vector<Caesar::Source*>& sources,Caesar::Img* inputImg,bool computeStatsAndBkg=true,Caesar::BkgData* inputBkgData=0);
		
		/** 
		\brief Select sources
 		*/
		int SelectSources(std::vector<Caesar::Source*>& sources);
		bool IsGoodSource(Caesar::Source* aSource);
		bool IsPointLikeSource(Caesar::Source* aSource);

		/** 
		\brief ThrowProgressEvent
 		*/
		int PushWorkerProgressEvent();

		/** 
		\brief PushWorkerDataEvent
 		*/
		int PushWorkerDataEvent(Caesar::WorkerData* workerData);

	private:
		SFinder* m_device;	

		//- the task's configuration
    yat::Task::Config m_cfg;

		//Bool flags
		std::atomic<bool> m_stopped;

		//- Mutex
		yat::Mutex m_mutex;
	
};//close class


}//close namespace


#endif
