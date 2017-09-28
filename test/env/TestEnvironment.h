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
* @file TestEnvironment.h
* @class TestEnvironment
* @brief Class to setup the test enviroment 
*
* Class to setup the test enviroment (performed at beginning and end of all tests)
* @author S. Riggi
* @date 15/01/2016
*/

#ifndef TEST_ENVIRONMENT_H
#define TEST_ENVIRONMENT_H

#include <gtest/gtest.h>
#include <gtest/gtest_prod.h>


class TestEnvironment : public ::testing::Environment { 

	public: 
		/** 
		\brief Standard constructor
 		*/
  	TestEnvironment();
		/** 
		\brief Destructor
 		*/
		virtual ~TestEnvironment();

	public:  
		/** 
		\brief Method to set up the environment (executed before all tests)
 		*/
		virtual void SetUp(); 
		/** 
		\brief Method to destroy the environment (executed after all tests)
 		*/
  	virtual void TearDown();

	private:
		
	public:
		
  	
};//close class 


#endif

