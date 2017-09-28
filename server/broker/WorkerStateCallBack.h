#ifndef WorkerStateCallBack_H
#define WorkerStateCallBack_H

#include <tango.h>
#include <SFinderBroker.h>


namespace SFinderBroker_ns {
  
class SFinderBroker;

class WorkerStateCallBack : public Tango::CallBack { 
		
	public:
		WorkerStateCallBack(SFinderBroker* dev);
		~WorkerStateCallBack();
		   
	public:
		void push_event(Tango::EventData*);
			
	private:
		static SFinderBroker* device;	

	friend class SFinderBroker;

};

}//close namespace

#endif
