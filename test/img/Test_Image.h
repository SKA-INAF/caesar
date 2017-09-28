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
* @file Test_Image.h
* @class Test_Image
* @brief Class for testing Image class
*
* Test class for Image class
* @author S. Riggi
* @date 08/08/2017
*/

#ifndef IMAGE_TEST_H
#define IMAGE_TEST_H

#include <gtest/gtest.h>
#include <gtest/gtest_prod.h>
#include <TestEnvironment.h>

#include <Image.h>

namespace Caesar {


class ImageTest : public ::testing::Test {

	protected:

		ImageTest() {
    	// You can do set-up work for each test here.
			SetUpData();

		}//close constructor

  	virtual ~ImageTest() {
    	// You can do clean-up work that doesn't throw exceptions here.
			if(m_img) {
				delete m_img;
				m_img= 0;
			}

  	}//close destructor

	protected:
		// If the constructor and destructor are not enough for setting up
  	// and cleaning up each test, you can define the following methods:
  	virtual void SetUp() {
			// Code here will be called immediately after the constructor (right
    	// before each test).
			

  	}

  	virtual void TearDown() {
			// Code here will be called immediately after each test (right
    	// before the destructor).

			
 	 	}

		// Per-test-case set-up.
  	// Called before the first test in this test case.
  	// Can be omitted if not needed.
		static void SetUpTestCase() {
			
			
		}//close SetUpTestCase()
	
		// Per-test-case tear-down.
  	// Called after the last test in this test case.
  	// Can be omitted if not needed.
  	static void TearDownTestCase() {
    	
  	}//close TearDownTestCase()

		static void SetUpData(){
		
			//Initialize image
			m_img= new Caesar::Image();

			//Fill image
			
		}//close SetUpData()

	// Objects declared here can be used by all tests in the test case for Foo.
	//...
	//...
	static Caesar::Image* m_img;
	static double m_numeric_tol;	

};//close class 


}//close namespace

#endif
