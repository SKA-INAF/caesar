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
* @file Test_Image.cc
* @class Test_Image
* @brief Class for testing Image class
*
* Test class for Image class
* @author S. Riggi
* @date 08/08/2017
*/


#include <Test_Image.h>
#include <Image.h>

#include "gtest/gtest.h"
#include "gmock/gmock-matchers.h"

#include <iostream>


namespace Caesar {

//Set static vars
Image* ImageTest::m_img= 0;
double ImageTest::m_numeric_tol= 1.e-6;

//===========================================
//==           DEFINE TESTS              ====
//===========================================

TEST_F(ImageTest, TestComputeStats)
{
	
	/*
	//Expected value
	//NB: Using data generated in astropy example (ref. http://docs.astropy.org/en/stable/stats/robust.html#sigma-clipping)
	double median_exp= 0.032658645;//computed using python numpy
	//double median_exp= 0.03265865;//computed using R

	//Compute median
	bool useParallelVersion= false;
	double median= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);
	double diff= fabs(median-median_exp);

	//Compute median with parallel algorithms
	useParallelVersion= true;
	double median_parallel= Caesar::StatsUtils::GetMedianFast(m_data,useParallelVersion);
	double diff_parallel= fabs(median_parallel-median_exp);

	//Print text
	std::stringstream ss;
	ss<<"StatsUtilsTest::TestMedianFast(): INFO: median_exp="<<median_exp<<", median="<<median<<" (diff="<<diff<<"), median_parallel="<<median_parallel<<" (diff="<<diff_parallel<<")";
	std::cerr<<ss.str()<<std::endl;

	//Check if median computation is correct within numerical tolerance
	ASSERT_LT(diff,m_numeric_tol);
	ASSERT_LT(diff_parallel,m_numeric_tol);
	*/

}//close TestComputeStats()



}//close namespace


