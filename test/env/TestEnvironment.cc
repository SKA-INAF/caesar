// **************************************************************************
// * License and Disclaimer                                                 *
// *                                                                        *
// * Copyright 2016 Simone Riggi																			      *
// *																																	      *
// * This file is part of SKA DSH.LMC 																		  *
// * SKA DSH.LMC is free software: you can redistribute it and/or modify it *
// * under the terms of the GNU General Public License as published by      *
// * the Free Software Foundation, either * version 3 of the License,       * 
// * or (at your option) any later version.                                 * 
// * Caesar is distributed in the hope that it will be useful, but 			    *
// * WITHOUT ANY WARRANTY; without even the implied warranty of             * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
// * See the GNU General Public License for more details. You should        * 
// * have received a copy of the GNU General Public License along with      * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                      *
// **************************************************************************

/**
* @file TestEnvironment.cc
* @class TestEnvironment
* @brief Class to setup the test enviroment 
*
* Class to setup the test enviroment (performed at beginning and end of all tests)
* @author S. Riggi
* @date 15/01/2016
*/

#include <TestEnvironment.h>

#include "gtest/gtest.h"


// Check if crash reporting is used.
#if defined(ENABLE_CRASH_REPORT)
#  include <crashreporting/crash_report.h>
#else
#  define DECLARE_CRASH_HANDLER
#  define INSTALL_CRASH_HANDLER
#endif

DECLARE_CRASH_HANDLER;

//Init static members
//...

TestEnvironment::TestEnvironment()
{


}

TestEnvironment::~TestEnvironment()
{


}


void TestEnvironment::SetUp(){

	//Setup 
	//...

}//close SetUp()


void TestEnvironment::TearDown(){


}//close TearDown()



