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
* @file SysUtils.cc
* @class SysUtils
* @brief Utility functions for system tasks
*
* Utility functions for system tasks
* @author S. Riggi
* @date 23/08/2010
*/


#include <SysUtils.h>
#include <CodeUtils.h>
#include <Logger.h>

//CFITSIO headers
#include <fitsio.h>

//Boost headers
#include <boost/filesystem.hpp>

//ROOT headers
#include <TObject.h>

//OpenMP headers
#ifdef OPENMP_ENABLED
  #include <omp.h>
#endif

//MPI headers
#ifdef MPI_ENABLED
	#include <mpi.h>
#endif

//C++ headers
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>
#include <thread>

using namespace std;

ClassImp(Caesar::SysUtils)
ClassImp(Caesar::FileInfo)

namespace Caesar {


SysUtils::SysUtils(){

}

SysUtils::~SysUtils(){

}


bool SysUtils::CheckDir(std::string path)
{
	//Check input file path
	if(path==""){	
		WARN_LOG("Empty filename given!");
		return false;
	}

	try {
		//Check if file exists on filesystem
		boost::filesystem::path file_path(path.c_str());
		if (!boost::filesystem::exists(file_path)){
			ERROR_LOG("File "<<path<<" not found in local filesystem!");
			return false;
		}

		//Check if directory
		if (!boost::filesystem::is_directory(file_path)){
    	ERROR_LOG("File ("<<file_path<<") is not a directory!");
			return false;
    }
	
	}//close try
	catch (const boost::filesystem::filesystem_error& ex) {
    ERROR_LOG("Exception detected while checking file (err: "<<ex.what()<<")!");
		return false;
  }

	return true;

}//close CheckDir()


bool SysUtils::CheckFile(std::string path,Caesar::FileInfo& info,bool match_extension,std::string extension)
{
	//Check input file path
	if(path==""){	
		WARN_LOG("Empty filename given!");
		return false;
	}
	if(!(&info)){
		ERROR_LOG("Null ptr to file info struct given!");
		return false;
	}

	//Check if file actually exists on filesystem
  try {
		//Check if file exists on filesystem
		boost::filesystem::path file_path(path.c_str());
		if (!boost::filesystem::exists(file_path)){
			ERROR_LOG("File "<<path<<" not found in local filesystem!");
			return false;
		}
		if (!boost::filesystem::is_regular_file(file_path)){
			ERROR_LOG("File "<<path<<" is not a regular file!");
			return false;
		}
		if (boost::filesystem::is_directory(file_path)){
    	ERROR_LOG("File ("<<file_path<<") is a directory!");
			return false;
    }
	
		//Get filename and extension
		if(!file_path.has_filename()){
			ERROR_LOG("File ("<<file_path<<") does not have a filename!");
			return false;
		}

		//Set file info
		info.filename= file_path.filename().string();
		info.filename_wext= file_path.stem().string();
		info.size= boost::filesystem::file_size(file_path);
        	
		//Check extension
		if(!file_path.has_extension()){
			ERROR_LOG("Given file without extension!");
			return false;
		}

		std::string file_extension= file_path.extension().string();	
		info.extension= file_extension;
						
		if(match_extension && file_extension!=extension){
			ERROR_LOG("Invalid file extension detected ("<<file_extension<<"!="<<extension<<")...");
			return false;
		}
	
		//Dump file info
		//info.Print();
		std::string info_printable= info.GetPrintable();
		DEBUG_LOG(info_printable);

  }//close try block

  catch (const boost::filesystem::filesystem_error& ex) {
    ERROR_LOG("Exception detected while checking file (err: "<<ex.what()<<")!");
		return false;
  }

	return true;

}//close CheckFile()

int SysUtils::GetFITSImageSize(const std::string& filename,long int& Nx,long int& Ny){

	//Init	
	Nx= 0;
	Ny= 0;

	//Open file
  int status = 0;
	fitsfile *fptr;//pointer to the FITS file, defined in fitsio.h
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if(status){
		ERROR_LOG("Failed to open given fits file "<<filename<<"!");	
		return -1;
	}

  //Read the NAXIS1 and NAXIS2 keyword to get image size
	int nfound;
	long naxes[2];
	fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  if (status){
		ERROR_LOG("Failed to get NAXIS keyword from given fits file "<<filename<<"!");	
		return -1;
	}
        
	Nx= naxes[0];
	Ny= naxes[1];
	
	return 0;

}//close GetFITSImageSize()


int SysUtils::GetNCores(){

	return std::thread::hardware_concurrency();

}//close GetNCores()

void SysUtils::SetOMPThreads(int nthreads){

	#ifdef OPENMP_ENABLED
		omp_set_dynamic(false);
  	omp_set_num_threads(nthreads);
	#endif

}//close SetOMPThreads()

int SysUtils::GetOMPThreads(){
	
	int nthreads= 0;

	#ifdef OPENMP_ENABLED
		nthreads= omp_get_num_threads();
	#endif

	return nthreads;

}//close GetOMPThreads()

int SysUtils::GetOMPMaxThreads(){
	
	int nthreads_max= 1;

	#ifdef OPENMP_ENABLED
		nthreads_max= omp_get_max_threads();
	#endif

	return nthreads_max;

}//close GetOMPThreads()

int SysUtils::GetOMPCores(){

	#ifdef OPENMP_ENABLED
		return omp_get_num_procs();
	#endif

	return 0;

}//close GetOMPCores()

int SysUtils::GetOMPThreadId(){

	#ifdef OPENMP_ENABLED
		return omp_get_thread_num();
	#endif

	return 0;

}//close GetOMPThreadId()

bool SysUtils::IsMPIInitialized(){

	bool isInitialized= false;
	#ifdef MPI_ENABLED
		int mpi_init_flag= 0;
		MPI_Initialized(&mpi_init_flag);
		if(mpi_init_flag==1) {
			//DEBUG_LOG("MPI is initialized for this run...");
			isInitialized= true;	
		}
		else {
			//WARN_LOG("MPI was not initialized for this run (hint: call MPI_Init), will run on single processor...");
			isInitialized= false;
		}
	#endif
	
	return isInitialized;

}//close IsMPIInitialized()

int SysUtils::GetProcMemoryInfo(ProcMemInfo& info)
{
	//Stores each word in status file
	char buffer[1024] = "";
	unsigned long currRealMem= 0;
	unsigned long peakRealMem= 0;
	unsigned long currVirtMem= 0;
	unsigned long peakVirtMem= 0;
	
	//Open linux file contains this process info
	FILE* file= fopen("/proc/self/status", "r");
	if (file == NULL) {
		ERROR_LOG("File /proc/self/status not found!");
		return -1;
	}

	//Read the entire file, recording mems in kB
	while (fscanf(file, " %1023s", buffer) == 1) 	
	{		
		int read_status= 0;
		if (strcmp(buffer, "VmRSS:") == 0) {
			read_status= fscanf(file, " %lu", &currRealMem);
		}
		else if (strcmp(buffer, "VmHWM:") == 0) {
			read_status= fscanf(file, " %lu", &peakRealMem);
		}
		else if (strcmp(buffer, "VmSize:") == 0) {
			read_status= fscanf(file, " %lu", &currVirtMem);
		}
		else if (strcmp(buffer, "VmPeak:") == 0) {
			read_status= fscanf(file, " %lu", &peakVirtMem);
		}
		else{//not interested in this field
			continue;
		}

		//Check read was succesful
		if(read_status!=1){
			ERROR_LOG("Failed to read memory field from file /proc/self/status!");
			return -1;
		}
    
	}//end loop

	//Close file
	fclose(file);

	//Set proc mem field
	info.realMem= currRealMem;
	info.realPeakMem= peakRealMem;	
	info.virtMem= currVirtMem;
	info.virtPeakMem= peakVirtMem;

	return 0;

}//close GetProcMemoryInfo()


std::string SysUtils::GetHost() 
{
	char hostname[HOST_NAME_MAX];
	gethostname(hostname, HOST_NAME_MAX);
	return std::string(hostname);
		
}

int SysUtils::GetProcId() 
{
	int procid= 0;			
	#ifdef MPI_ENABLED
		bool isMPIInitialized= IsMPIInitialized();
		if(isMPIInitialized){
			MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		}
		else{
			procid= 0;
		}
	#endif

	return procid;
		
}


std::string SysUtils::GetAsciiLogo()
{
	std::stringstream ss;	
	ss<<R"(   ______                              _____ _______           __         )"<<endl;
  ss<<R"(  / ____/___ ____  _________ ______   / ___// ____(_)___  ____/ /__  _____)"<<endl;
 	ss<<R"( / /   / __ `/ _ \/ ___/ __ `/ ___/   \__ \/ /_  / / __ \/ __  / _ \/ ___/)"<<endl;
	ss<<R"(/ /___/ /_/ /  __(__  ) /_/ / /      ___/ / __/ / / / / / /_/ /  __/ /    )"<<endl;
	ss<<R"(\____/\__,_/\___/____/\__,_/_/      /____/_/   /_/_/ /_/\__,_/\___/_/     )"<<endl;
                                                                          
	return ss.str();
}


void SysUtils::PrintAsciiLogo()
{
	std::string logo_str= GetAsciiLogo();
	cout<<logo_str<<endl;
}


}//close namespace



