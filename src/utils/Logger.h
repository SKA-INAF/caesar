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
* @file Logger.h
* @class Logger
* @brief Logger
*
* Logger class for managing logging msg
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef Logger_h
#define Logger_h 1

//#include <CodeUtils.h>

#include <TObject.h>
#include <TMath.h>

#ifdef USE_TANGO
#include <tango.h>
#include <log4tango.h>
#endif

//== LOG4CXX HEADERS ==
#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/simplelayout.h>
#include <log4cxx/patternlayout.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/rollingfileappender.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/net/syslogappender.h>
#include <log4cxx/helpers/exception.h>

//Boost
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

//#include <syslog.h>


#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <stdlib.h>

#include <vector>
#include <algorithm>
#include <map>
#include <string>

using namespace std;


namespace Caesar {


class Logger : public TObject {

	public:
		/**
		* \brief Constructor
		*/
		Logger(std::string level="OFF",std::string tag="logger")
			: m_level(level), m_tag(tag)		
		{

		};

		/**
		* \brief Destructor
		*/
		virtual ~Logger(){};

	public:

		/**
		* \brief Returns log4cxx log level from log4tango log level
		*/
		#ifdef USE_TANGO
			static log4cxx::LevelPtr GetMappedLogLevel(int level) {	
				if(level==log4tango::Level::DEBUG) return log4cxx::Level::getDebug();
				else if(level==log4tango::Level::INFO) return log4cxx::Level::getInfo();
				else if(level==log4tango::Level::WARN) return log4cxx::Level::getWarn();
				else if(level==log4tango::Level::ERROR) return log4cxx::Level::getError();
				else if(level==log4tango::Level::FATAL) return log4cxx::Level::getFatal();
				else return log4cxx::Level::getOff();
				return log4cxx::Level::getOff();
			}//close GetMappedLogLevel()
			
			/**
			* \brief Returns log4tango log level from string 
			*/
			static log4tango::Level::Value GetTangoLogLevelFromString(std::string sLevel) {
  			if (sLevel == "DEBUG") return log4tango::Level::DEBUG;
				else if (sLevel == "INFO") return log4tango::Level::INFO; 
				else if (sLevel == "WARN") return log4tango::Level::WARN;
				else if (sLevel == "ERROR") return log4tango::Level::ERROR;
				else if (sLevel == "FATAL") return log4tango::Level::FATAL;
				else if (sLevel == "OFF") return log4tango::Level::OFF;
				else return log4tango::Level::OFF;
				return log4tango::Level::OFF;
			}	
		#endif

	public:
		/**
		* \brief Pure abstract method to initialize logger
		*/
		virtual int Init() = 0;
		
		/**
		* \brief Get machine host name
		*/
		virtual std::string GetHost() const {
			char hostname[HOST_NAME_MAX];
			gethostname(hostname, HOST_NAME_MAX);
			std::string syslog_host(hostname); 
			return syslog_host;
		}

		/**
		* \brief Log a message (optionally prepended by a prefix) at given log level in string format
		*/
		virtual void Log(const std::string sLevel, const std::string msg, const std::string msg_prefix=""){
			log4cxx::LevelPtr logLevel= log4cxx::Level::toLevel(sLevel,log4cxx::Level::getOff());
			if(logger && logger->isEnabledFor(logLevel) ){
				log4cxx::helpers::MessageBuffer oss_; \
      	logger->forcedLog(logLevel, oss_.str(oss_ << msg_prefix << msg), log4cxx::spi::LocationInfo(m_tag.c_str(),__LOG4CXX_FUNC__,__LINE__) );
			}
		}	

		/**
		* \brief Set logger log level threshold
		*/
		virtual void SetLogLevel(const std::string sLevel){
			log4cxx::LevelPtr thisLogLevel= logger->getLevel();//retrieve the current log level
			if(!thisLogLevel) thisLogLevel= log4cxx::Level::getOff();
			logger->setLevel(log4cxx::Level::toLevel(sLevel,thisLogLevel));// In case of wrong string the Log level is set to the old or to OFF
		}//close SetLogLevel()

		/**
		* \brief Log a message (optionally prepended by a prefix) at given Tango log level
		*/
		#ifdef USE_TANGO
			virtual void Log(int level, const std::string msg, const std::string msg_prefix=""){
				log4cxx::LevelPtr logLevel= GetMappedLogLevel(level);					
				if(logger && logger->isEnabledFor(logLevel) ){
					log4cxx::helpers::MessageBuffer oss_; \
      		logger->forcedLog(logLevel, oss_.str(oss_ << msg_prefix << msg), log4cxx::spi::LocationInfo(m_tag.c_str(),__LOG4CXX_FUNC__,__LINE__) );
				}
			}//close Log()
		#endif

		/**
		* \brief Set logger log level threshold from Tango log level
		*/
		#ifdef USE_TANGO	
			virtual void SetLogLevel(int level){
				log4cxx::LevelPtr logLevel= GetMappedLogLevel(level);
				logger->setLevel(logLevel);
			}//close SetLogLevel()
		#endif


	protected:
		std::string m_level;
		std::string m_tag;		
		log4cxx::LoggerPtr logger;
		log4cxx::LayoutPtr layout;
		log4cxx::AppenderPtr appender;		
		
	ClassDef(Logger,1)		

};//close Logger class


#ifdef __MAKECINT__
#pragma link C++ class Logger+;
#endif


class SysLogger : public Logger {

	public:
		/**
		* \brief Syslog logger constructor
		*/
		SysLogger(const std::string level,const std::string tag,const std::string facility)
			: Logger(level,tag), m_facility(facility)
		{
			
		};

		/**
		* \brief Destructor
		*/
		virtual ~SysLogger(){};

	private:
		

	public:
		/**
		* \brief Initialize syslog logger
		*/
		virtual int Init(){
			//Create logger
			logger= log4cxx::LoggerPtr(log4cxx::Logger::getLogger(m_tag.c_str()));
			if(!logger) return -1;
			return 0;

			//Define log layout
			layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%c %m%n") );
			if(!layout) return -1;

			//Get host & facility code	
			std::string syslog_host= GetHost();
			int syslog_facility_code= GetFacilityCode(m_facility);

			//Create and add appender
			appender= log4cxx::AppenderPtr( new log4cxx::net::SyslogAppender(layout, syslog_host, syslog_facility_code) );
			appender->setOption("SYSLOGHOST",syslog_host);
			logger->addAppender(appender);
	
			//Set logging level (OFF is set if a wrong name is given!)
			log4cxx::LevelPtr level= log4cxx::Level::toLevel(m_level,log4cxx::Level::getOff());
			logger->setLevel(level);  
			
			return 0;
		}//close Init()

		/**
		* \brief Returns syslog facility code
		*/
		virtual int GetFacilityCode(const std::string syslog_facility) const {
			std::string default_syslog_facility= "LOCAL6";
			int syslog_facility_code= log4cxx::net::SyslogAppender::getFacility(syslog_facility);
			if(syslog_facility_code==-1) syslog_facility_code= log4cxx::net::SyslogAppender::getFacility(default_syslog_facility);
			
			return syslog_facility_code;
		}//close GetFacilityCode()
		
	protected:
		std::string m_facility;
		
	ClassDef(SysLogger,1)		

};//close SysLogger class

#ifdef __MAKECINT__
#pragma link C++ class SysLogger+;
#endif


class FileLogger : public Logger {

	public:
		/**
		* \brief File logger constructor
		*/
		FileLogger(const std::string level,const std::string tag,const std::string filename,bool appendFlag,const std::string maxFileSize,int maxBackupFiles)
			: Logger(level,tag), m_filename(filename), m_appendFlag(appendFlag), m_maxFileSize(maxFileSize), m_maxBackupFiles(maxBackupFiles)
		{
			
		};

		/**
		* \brief Destructor
		*/
		virtual ~FileLogger(){};

	public:

		/**
		* \brief Initialize file logger
		*/
		virtual int Init(){
			//Create logger
			logger= log4cxx::LoggerPtr(log4cxx::Logger::getLogger(m_tag));
			if(!logger) return -1;
			
			//Define log layout
			//layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] [%l] %m%n") );
			layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] %m%n") );
			if(!layout) return -1;

			//Create and add appender
			//appender= log4cxx::AppenderPtr( new log4cxx::FileAppender(layout,m_filename,m_appendFlag) );	
			appender= log4cxx::AppenderPtr( new log4cxx::RollingFileAppender(layout,m_filename,m_appendFlag) );	
			//appender->setMaxBackupIndex(m_maxBackupFiles);
			//appender->setMaxFileSize(m_maxFileSize);
			std::stringstream maxBackupFiles_stream;
			maxBackupFiles_stream << m_maxBackupFiles;
			appender->setOption("MAXBACKUPINDEX",maxBackupFiles_stream.str());//obsolete??
			appender->setOption("MAXFILESIZE",m_maxFileSize);
			appender->setOption("FILENAME",m_filename);
			logger->addAppender(appender);
	
			//Set logging level (OFF is set if a wrong name is given!)
			log4cxx::LevelPtr level= log4cxx::Level::toLevel(m_level,log4cxx::Level::getInfo());
			logger->setLevel(level);  

			return 0;
		}//close Init()

		
	protected:
		std::string m_filename;
		bool m_appendFlag;
		std::string m_maxFileSize;
		int m_maxBackupFiles; 

	ClassDef(FileLogger,1)		

};//close FileLogger class

#ifdef __MAKECINT__
#pragma link C++ class FileLogger+;
#endif


class ConsoleLogger : public Logger {

	public:
		/**
		* \brief Console logger constructor
		*/
		ConsoleLogger(const std::string level,const std::string tag,const std::string target)
			: Logger(level,tag), m_target(target)
		{
			
		};

		/**
		* \brief Destructor
		*/
		virtual ~ConsoleLogger(){};

	public:
		/**
		* \brief Initialize console logger
		*/
		virtual int Init(){
			//Create logger
			logger= log4cxx::LoggerPtr(log4cxx::Logger::getLogger(m_tag));
			if(!logger) return -1;
			
			//Define log layout
			//layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] [%l] %m%n") );
			layout= log4cxx::LayoutPtr( new log4cxx::PatternLayout("%d %-5p [%c] %m%n") );
			if(!layout) return -1;

			//Create and add appender
			appender= log4cxx::AppenderPtr( new log4cxx::ConsoleAppender(layout,m_target) );	
			appender->setOption("TARGET",m_target);
			logger->addAppender(appender);
	
			//Set logging level (OFF is set if a wrong name is given!)
			log4cxx::LevelPtr level= log4cxx::Level::toLevel(m_level,log4cxx::Level::getInfo());
			logger->setLevel(level);  

			return 0;
		}//close Init()

	
	protected:
		std::string m_target;
		

	ClassDef(ConsoleLogger,1)		

};//close ConsoleLogger class

#ifdef __MAKECINT__
#pragma link C++ class ConsoleLogger+;
#endif

/**
* \brief Logger target enumerations
*/
enum LoggerTarget {
	eCONSOLE_TARGET= 1,
	eFILE_TARGET= 2,
	eSYSLOG_TARGET= 3
};

class LoggerManager : public TObject {
	
	public:
		/**
		* \brief Returns unique instance of logger
		*/
		static LoggerManager& Instance() {
    	// Since it's a static variable, if the class has already been created,
      // It won't be created again.
      // And it is thread-safe in C++11.
      static LoggerManager myInstance;
 
      // Return a reference to our instance.
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    LoggerManager(LoggerManager const&) = delete;             // Copy construct
    LoggerManager(LoggerManager&&) = delete;                  // Move construct
    LoggerManager& operator=(LoggerManager const&) = delete;  // Copy assign
    LoggerManager& operator=(LoggerManager &&) = delete;      // Move assign

	public:	
		/**
		* \brief Define log target enumerations
		*/
		enum LogTarget{eSysLog=1,eConsole=2,eFile=3};

		/**
		* \brief Create a syslog logger
		*/
  	static int CreateSysLogger(const std::string level,const std::string tag="syslogger",const std::string facility="local6") 
		{
			//Skip if already created
			if(m_logger) {	
				cerr<<"CreateSysLogger(): WARN: Already created...skip!"<<endl;
				return -1;
			}

			//Create syslogger
			m_target= eSysLog;
			m_logger= new SysLogger(level,tag,facility);
			m_logger->Init();

			return 0;
		}//close CreateSysLogger()

		static int CreateFileLogger(const std::string level,const std::string tag="filelogger",const std::string filename="out.log",bool appendFlag=false,std::string maxFileSize="10MB",int maxBackupFiles=2) {
				
			//Skip if already created
			if(m_logger) {	
				return -1;
			}

			//Create syslogger
			m_target= eFile;
			m_logger= new FileLogger(level,tag,filename,appendFlag,maxFileSize,maxBackupFiles);
			m_logger->Init();
			
			return 0;
		}//close CreateFileLogger()

		/**
		* \brief Create a console logger
		*/
		static int CreateConsoleLogger(const std::string level,const std::string tag="consolelogger",const std::string target="System.out") 
		{				
			//Skip if already created
			if(m_logger) {	
				return -1;
			}

			//Create syslogger
			m_target= eFile;
			m_logger= new ConsoleLogger(level,tag,target);
			m_logger->Init();
			
			return 0;
		}//close CreateConsoleLogger()

		Logger* GetLogger(){return m_logger;}
		int GetLoggerTarget(){return m_target;}

		

	protected:
		/**
		* \brief Logger manager protected constructor
		*/
		LoggerManager(){
			
		};
		/**
		* \brief Destructor
		*/
		virtual ~LoggerManager(){
			
		};

	private:
		static int m_target;
		static Logger* m_logger;

	ClassDef(LoggerManager,1)	

};//close LoggerManager()

#ifdef __MAKECINT__
#pragma link C++ class LoggerManager+;
#endif

class ScopedLogger {
	public:
		/**
		* \brief Scoped logger constructor
		*/
		ScopedLogger(std::string level,std::string prefix="",std::string device_name="")
  		: m_level(level), m_msgprefix(prefix), m_deviceName(device_name) 
		{}

			
		/**
		* \brief Destructor
		*/
  	~ScopedLogger(){ 
			//logging command
			Logger* logger= LoggerManager::Instance().GetLogger();
			if(logger) logger->Log(m_level, m_sstream.str(), m_msgprefix);			

			//Tango logging
			#ifdef USE_TANGO	
				try {
					Tango::DeviceImpl* ds_impl= Tango::Util::instance(false)->get_device_by_name(m_deviceName);
					log4tango::Level::Value tango_level= Logger::GetTangoLogLevelFromString(m_level);
					if(ds_impl && ds_impl->get_logger()->is_level_enabled(tango_level) ) {
						ds_impl->get_logger()->log(tango_level, m_msgprefix+m_sstream.str());
					}
				}
				catch(Tango::DevFailed& e){
					//cerr<<"ScopedLogger::~ScopedLogger(): ERROR: Failed to log with Tango logger!"<<endl;
					//Tango::Except::print_exception(e);
				}
			#endif
		}//close destructor

	public:
		/**
		* \brief Returns message string stream to be logged
		*/
		std::stringstream& stream(){ 	
			return m_sstream; 
		}	
	
	private:
		std::stringstream m_sstream;
  	std::string m_level;
		std::string m_msgprefix;
		std::string m_deviceName;

};//close ScopedLogger

/**
* \brief Returns Tango device name in which log is emitted
*/
#ifdef USE_TANGO	
inline std::string GetDeviceName(){	
	std::string dev_name= "";
	try {
		//dev_name= Tango::Util::instance(false)->get_ds_name();
		std::vector<Tango::DeviceImpl*> dev_list= Tango::Util::instance(false)->get_device_list("*");	
		if(dev_list.empty()) dev_name= "";
		else dev_name= dev_list[0]->get_name();
	}
	catch(Tango::DevFailed& e){
		dev_name= "";
	}
	return dev_name;
}

inline bool areInTangoServer(){
	try {
		std::vector<Tango::DeviceImpl*> dev_list= Tango::Util::instance(false)->get_device_list("*");	
		if(!dev_list.empty()) return true;
	}
	catch(Tango::DevFailed& e){
		return false;
	}
	return false;
}
#endif

/**
* \brief Returns class name in which log is emitted
*/
inline std::string getClassName(std::string fullFuncName,std::string funcName){

	//Init pattern to be searched
	std::string result= "";
	//std::string pattern("::(.*)::");//FIX ME!!!
	std::string pattern("([-A-Za-z0-9_]+)::");
	pattern+= funcName;

	//Create regex
	boost::regex expression;
	try {
  	expression = pattern;
  }
  catch (boost::regex_error& e) {
  	return result;
  }

	//Find match
	boost::smatch matches;
	if (boost::regex_search(fullFuncName, matches, expression) && matches.size()>1) {
		result= std::string(matches[1].first, matches[1].second);
		//result= std::string(matches[matches.size()-1].first, matches[matches.size()-1].second);
		//for(int i=0;i<matches.size();i++) cout<<"match no. "<<i<<"="<<matches[i]<<endl;
  }//close if
	
	return result;

}//close function

/**
* \brief Returns class name prefix in which log is emitted
*/
inline std::string getClassNamePrefix(std::string fullFuncName,std::string funcName){
	std::string className= getClassName(fullFuncName,funcName);
	std::string sprefix= "::";
	if(className=="") return className;
	return className+sprefix;
}

/**
* \brief Shortcut macro to get actual class name
*/
#define __CLASS__ getClassName(__PRETTY_FUNCTION__,__FUNCTION__)

/**
* \brief Shortcut macro to get actual class name prefix
*/
#define __CLASS_PREFIX__ getClassNamePrefix(__PRETTY_FUNCTION__,__FUNCTION__)

/**
* \brief Shortcut macro to get actual Tango device class name
* \param deviceInstance - Tango device server class instance
*/
#define __DEVICE_CLASS(deviceInstance) deviceInstance->get_device_class()->get_name()

/**
* \brief Shortcut macro to get actual Tango device name
* \param deviceInstance - Tango device server class instance
*/
#define __DEVICE_NAME(deviceInstance) deviceInstance->get_name()

//== LOG MACROS ===
/**
* \brief Shortcut macro to set a prefix for log message
*/
#define LOG_PREFIX \
	__CLASS_PREFIX__ + __FUNCTION__ + std::string("() - ")

#ifdef USE_TANGO
	#define DEVICE_NAME GetDeviceName()
	#define __INSIDE_TANGO_SERVER areInTangoServer()
#else 
	#define	DEVICE_NAME std::string("")
	#define __INSIDE_TANGO_SERVER std::string("")
#endif

/**
* \brief Shortcut macro to emit a log message with desired log level inside a given device
* \param DeviceName - Tango device server name
* \param Level - Message log level
* \param What - Message stream to be sent
*/
#define CAESAR_LOG(DeviceName, Level, What) \
	ScopedLogger(Level,LOG_PREFIX,DeviceName).stream() << What

/**
* \brief Shortcut macro to emit a log message with desired log level
* \param Level - Message log level
* \param What - Message stream to be sent
*/
#define LOG(Level, What) CAESAR_LOG(DEVICE_NAME,Level,What)

/**
* \brief Shortcut macro to emit a log message with INFO log level
* \param What - Message stream to be sent
*/
#define INFO_LOG(What) LOG("INFO",What)

/**
* \brief Shortcut macro to emit a log message with WARN log level
* \param What - Message stream to be sent
*/
#define WARN_LOG(What) LOG("WARN",What)

/**
* \brief Shortcut macro to emit a log message with DEBUG log level
* \param What - Message stream to be sent
*/
#define DEBUG_LOG(What) LOG("DEBUG",What)

/**
* \brief Shortcut macro to emit a log message with ERROR log level
* \param What - Message stream to be sent
*/
#define ERROR_LOG(What) LOG("ERROR",What)

/**
* \brief Shortcut macro to emit a log message with FATAL log level
* \param What - Message stream to be sent
*/
#define FATAL_LOG(What) LOG("FATAL",What)



//--> DEVICE LOG MACROS ===
/**
* \brief Shortcut macro to emit a log message with desired log level inside a Tango device server instance
* \param What - Message stream to be sent
* \param Level - Message log level
* \param deviceInstance - Tango device server class instance
*/
#define DEV_LOG(deviceInstance,Level, What) \
	ScopedLogger(Level,LOG_PREFIX,__DEVICE_NAME(deviceInstance)).stream() << What

#define _INFO_LOG(What) DEV_LOG(this,"INFO",What)
#define _WARN_LOG(What) DEV_LOG(this,"WARN",What)
#define _DEBUG_LOG(What) DEV_LOG(this,"DEBUG",What)
#define _ERROR_LOG(What) DEV_LOG(this,"ERROR",What)
#define _FATAL_LOG(What) DEV_LOG(this,"FATAL",What)

#define __INFO_LOG(deviceInstance,What) DEV_LOG(deviceInstance,"INFO",What)
#define __WARN_LOG(deviceInstance,What) DEV_LOG(deviceInstance,"WARN",What)
#define __DEBUG_LOG(deviceInstance,What) DEV_LOG(deviceInstance,"DEBUG",What)
#define __ERROR_LOG(deviceInstance,What) DEV_LOG(deviceInstance,"ERROR",What)
#define __FATAL_LOG(deviceInstance,What) DEV_LOG(deviceInstance,"FATAL",What)

}//close namespace 

#endif

