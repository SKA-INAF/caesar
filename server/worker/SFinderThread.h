#ifndef SFinderThread_H
#define SFinderThread_H

#include <SFinder.h>


#include <tango.h>

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


#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>


using namespace Caesar;

namespace SFinder_ns {
  
class SFinder;

class SFinderThread : public Tango::LogAdapter {
  
	private:
  	static SFinder* device;	
			
  public :
    SFinderThread(SFinder* dev);
    ~SFinderThread();

	public:
		/** 
		\brief Start the thread
 		*/		
		//void Start(std::vector<Caesar::WorkerTask*>& tasks){
		void Start(const std::vector<Caesar::WorkerTask*> tasks){
			m_stopThread = false;
    	//m_thread = std::thread(&SFinderThread::Run,this,std::ref(tasks));
			m_thread = std::thread(&SFinderThread::Run,this,tasks);
    }
		
		/** 
		\brief Stop
 		*/
		void Stop(){
			DEBUG_LOG("Called thread Stop()...");
			std::lock_guard<std::mutex> lock( m_mutex );
			m_stopThread = true;
			if(m_thread.joinable()) m_thread.join();
			DEBUG_LOG("Called thread Stop(): done");
		}
     		
	private:
		/** 
		\brief Main thread function 
 		*/
		//void Run(std::vector<Caesar::WorkerTask*>& tasks);
		void Run(const std::vector<Caesar::WorkerTask*> tasks);

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
		//int PushWorkerDataEvent(Caesar::WorkerData& workerData);
		int PushWorkerDataEvent(Caesar::WorkerData* workerData);
		

	private:	
  	static log4tango::Logger* m_logger;
    static SFinder* m_device;	
		
		mutable std::mutex m_mutex;	
		std::atomic<bool> m_stopThread;
		std::thread m_thread;

	friend class SFinder;

};//close class


}//close namespace

#endif
