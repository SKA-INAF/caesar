#include <Logger.h>

#define BOOST_LOG_DYN_LINK 1

ClassImp(Caesar::Logger)
ClassImp(Caesar::SysLogger)
ClassImp(Caesar::FileLogger)
ClassImp(Caesar::ConsoleLogger)
ClassImp(Caesar::LoggerManager)

namespace Caesar {

int LoggerManager::m_target;
Logger* LoggerManager::m_logger= 0;

}//close namespace
