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
* @file SysUtils.h
* @class SysUtils
* @brief Utility functions for system tasks
*
* Utility functions for system tasks
* @author S. Riggi
* @date 23/08/2010
*/


#ifndef SysUtils_h
#define SysUtils_h 1

#include <TObject.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <time.h>
#include <ctime>


using namespace std;


namespace Caesar {

class FileInfo : public TObject {

	public:
		FileInfo() {
			path= "";
			filename= "";
			filename_wext= "";
			extension= "";
			size= 0;
		}
		virtual ~FileInfo(){};

	public:
		std::string GetPrintable() {
			std::stringstream ss;
			ss<<"FileInfo: ";
			ss<<"path: "<<path<<", ";
			ss<<"name: "<<filename<<", ";
			ss<<"namewext: "<<filename_wext<<", ";
			ss<<"ext: "<<extension<<", ";
			ss<<"size (kB): "<<size/1024.;
			return ss.str();
		}
		void Print(){
			cout <<"*** FILE INFO ***"<<endl;
			cout << "File Path: "<<path<<endl;
			cout << "File Name: "<<filename<<endl;
			cout << "File WExt: "<<filename_wext<<endl;
			cout << "File Extension: "<< extension<<endl;
			cout << "File size (kB): "<<size/1024.<<endl;
			cout <<"***********************"<<endl;
		}

		
	public: 
		std::string path;
		std::string filename;
		std::string filename_wext;
		std::string extension;
		int size;

	ClassDef(FileInfo,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class FileInfo+;
#endif	
	

class SysUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SysUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~SysUtils();

		
	public:

		/**
		* \brief Check if a file exists in filesystem
		*/
		static bool CheckFile(std::string path,FileInfo& info,bool match_extension=false,std::string extension="");
	
		/**
		* \brief Compute the difference between two timestamps
		*/
		static timespec TimeDiff(timespec start, timespec end) {
			timespec temp;
			if ((end.tv_nsec-start.tv_nsec)<0) {
				temp.tv_sec = end.tv_sec-start.tv_sec-1;
				temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
			} 
			else {
				temp.tv_sec = end.tv_sec-start.tv_sec;
				temp.tv_nsec = end.tv_nsec-start.tv_nsec;
			}
			return temp;
		}

		/**
		* \brief Compute the sum of two timestamps
		*/
		static timespec TimeSum (timespec time1, timespec time2) {
			timespec  result ;
			result.tv_sec = time1.tv_sec + time2.tv_sec ;
    	result.tv_nsec = time1.tv_nsec + time2.tv_nsec ;
    	if (result.tv_nsec >= 1000000000L) {
      	result.tv_sec++ ;  result.tv_nsec = result.tv_nsec - 1000000000L ;
    	}
    	return (result) ;
		}
		
		/**
		* \brief Convert a timestamp to seconds
		*/
		static double TimeToSec (timespec time){
    	return ((double)time.tv_sec + (time.tv_nsec/1.e+09)) ;
		}

		/**
		* \brief Convert a timestamp to nanoseconds
		*/
		static double TimeToNSec (timespec time){
			return (time.tv_sec*1.e+09 + (double)time.tv_nsec) ;
		}

		/**
		* \brief Get the size of a FITS image
		*/
		static int GetFITSImageSize(const std::string& filename,long int& Nx,long int& Ny);
		
		/**
		* \brief Get the number of cores in system
		*/
		static int GetNCores();
	
		/**
		* \brief Set the number of threads to be used by OpenMP (if enabled at build)
		*/
		static void SetOMPThreads(int nthreads);
		/**
		* \brief Get the number of threads currently used by OpenMP (return 0 if OMP is disabled at build)
		*/
		static int GetOMPThreads();
		/**
		* \brief Get the maximum number of threads that can be used by OpenMP (return 0 if OMP is disabled at build)
		*/
		static int GetOMPMaxThreads();

		/**
		* \brief Get the number of cores currently used by OpenMP (return 0 if OMP is disabled at build)
		*/
		static int GetOMPCores();	
		/**
		* \brief Get thread id (return 0 if OMP is disabled at build)
		*/
		static int GetOMPThreadId();

		/**
		* \brief Is MPI run initialized (return 0 if MPI is disabled at build)
		*/
		static bool IsMPIInitialized();

	private:
	
		ClassDef(SysUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class SysUtils+;
#endif	

}//close namespace


#endif 
