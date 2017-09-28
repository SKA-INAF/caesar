#ifndef WorkerSourceDataCallBack_H
#define WorkerSourceDataCallBack_H

#include <tango.h>
#include <SFinderBroker.h>


namespace SFinderBroker_ns {
  
class SFinderBroker;

class WorkerSourceDataCallBack : public Tango::CallBack { 
		
	public:
		WorkerSourceDataCallBack(SFinderBroker* dev);
		~WorkerSourceDataCallBack();
		   
	public:
		void push_event(Tango::EventData*);
			
	private:
		static SFinderBroker* device;	

	friend class SFinderBroker;

};

}//close namespace

#endif
