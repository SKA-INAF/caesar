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
* @file Test_StatsUtils.h
* @class Test_StatsUtils
* @brief Class for testing StatsUtils class
*
* Test class for StatsUtils class
* @author S. Riggi
* @date 08/08/2017
*/


#ifndef STATS_UTILS_TEST_H
#define STATS_UTILS_TEST_H

#include <gtest/gtest.h>
#include <gtest/gtest_prod.h>


#include <TestEnvironment.h>


namespace Caesar {


class StatsUtilsTest : public ::testing::Test {

	protected:

		StatsUtilsTest() {
    	// You can do set-up work for each test here.
			SetUpData();

		}//close constructor

  	virtual ~StatsUtilsTest() {
    	// You can do clean-up work that doesn't throw exceptions here.
			m_data.clear();
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
	
			//Clear data
			m_data.clear();

			//Fill test data
			m_data=
			{
				2.25327184e-01,   3.94143678e-01,  -2.29493730e-01,
        -8.75640089e-02,  -9.96064901e-02,   3.85906411e-01,
         1.89884161e-01,   4.71733517e+00,   2.67650234e+00,
         1.68872595e-01,   1.84668946e+00,  -3.08954219e-01,
         2.37605958e-01,  -2.57666813e+00,   1.84171765e-01,
         6.37455306e-02,   1.71366122e-01,   1.62250198e+00,
         3.04062634e+00,   5.60550279e+00,   6.05589039e+00,
        -4.99102700e+00,  -9.11065007e-02,   7.40439039e+00,
        -7.07987823e-02,  -2.74990259e-01,  -1.28723681e-01,
        -2.55721485e+00,   1.25046290e-01,  -3.20411531e-01,
        -2.20876668e-01,  -4.34874202e+00,  -1.47912599e-01,
         3.08602919e-01,  -2.58571382e-01,   5.34101739e-02,
        -7.85656365e-03,  -2.33618700e-01,   1.64546852e+00,
        -8.41770325e-01,   1.54358110e-01,   1.64700831e-01,
         1.07655408e+01,   2.67305590e-01,  -6.71060522e-02,
        -2.62736077e+00,   2.19931919e-01,   1.31052746e-01,
         1.28026305e-01,  -3.23391209e-01,  -4.86522488e-03,
        -1.47606182e-01,   2.13071656e+00,  -1.96300779e-02,
         1.82035782e-01,   6.34436430e-02,   1.61169337e+00,
        -9.32838193e-02,  -1.88889251e-01,  -8.20099386e-02,
        -3.40408277e-03,   7.58303471e-02,   1.10518361e+01,
        -8.45143033e-03,  -1.91189000e-01,  -6.91963551e-02,
        -7.56573455e+00,   9.62962948e-02,   2.04130583e+00,
         1.26523988e-02,   7.74653174e+00,   4.64362072e-02,
         3.39765189e-02,  -4.75843459e-02,   3.80960088e-01,
        -9.86639767e-02,  -1.08572295e-01,   8.32100093e-02,
        -2.31236486e-01,   1.56239620e-01,   2.98896909e-01,
        -4.13997005e-01,   8.52517462e-02,   3.52668248e+00,
        -1.27487405e-01,  -7.94543629e-02,  -2.65761155e-02,
        -5.95581759e-02,  -6.18025938e-02,   2.08926190e+00,
         2.30466313e-01,  -1.60713635e+00,  -1.62672852e-01,
         2.15476864e+00,   1.04212975e-01,  -1.15157594e-01,
         2.83906327e-02,  -6.38656834e-02,   4.83282838e+00,
         1.38949829e-01,   5.10559584e-01,  -2.76672791e-01,
         3.67991331e+00,   2.86905817e+00,  -2.37771852e-01,
        -1.01363271e-01,  -1.19262808e-01,  -1.05134593e-02,
        -3.87255961e-01,   1.51000241e+01,   1.04778205e-01,
        -9.47902396e-01,  -8.50527743e+00,   1.94800333e-02,
        -7.00222387e+00,  -5.54518551e-01,   2.01282979e+00,
         7.80186645e-02,   1.15638817e+01,   7.89378128e+00,
         9.69443274e+00,  -2.32207878e-02,   8.23683343e+00,
         4.12898572e-01,  -2.21081314e-02,   2.04034542e-01,
        -1.38409970e-01,   3.07275411e-01,   7.37024879e+00,
         1.21768767e-01,  -2.09050673e-01,   2.42229058e-01,
         1.37963633e-01,   2.60369246e-01,  -1.25617512e-01,
        -9.62054237e-02,   4.60783340e-01,   9.27810670e+00,
         1.64539303e+01,   2.27378273e-01,  -2.73219840e-01,
         1.16590736e-01,  -7.98898059e-02,   2.68473740e+00,
         4.65235678e+00,   3.16041473e+00,  -2.36328090e-02,
         1.69288535e+00,   1.39602335e+00,   1.51365822e-03,
        -2.66851694e-01,  -5.10818212e+00,   1.38754631e-01,
        -3.19146876e-02,  -2.67403119e-02,   1.00903691e+01,
        -3.88642343e+00,  -1.46135551e-01,  -7.69759618e-02,
         1.88703179e-02,   5.48501708e+00,  -5.73774385e-02,
         4.40907307e+00,   1.16918828e+01,   1.74305072e+00,
        -1.72799469e+00,   5.49032715e-02,  -1.78183017e-01,
         3.97911875e+00,  -6.24584502e-02,  -3.15334032e-02,
         4.51344699e-01,  -1.40940055e-01,   1.88652145e-01,
        -8.04858263e+00,  -2.37788991e-01,   1.54650595e-01,
         3.08750440e+00,  -3.93028031e+00,   1.21263905e-01,
        -3.51178117e-01,   9.01868924e-02,  -1.36802180e-01,
         3.31910159e-01,   4.28110212e+00,  -9.06771608e-02,
        -1.37567522e-01,  -2.42815481e-01,   9.88588358e-01,
        -5.60710990e-02,  -7.29387089e-02,   3.13407711e-02,
         1.15704300e-01,  -8.95878986e-02,   3.66681408e+00,
        -2.87558295e-01,   2.72906370e-01,  -1.37889837e-01,
        -1.30458720e-01,  -1.04237862e-01
			};

		}//close SetUpData()

	// Objects declared here can be used by all tests in the test case for Foo.
	//...
	//...
	static std::vector<double> m_data;
	static double m_numeric_tol;	

};//close class 


}//close namespace

#endif
