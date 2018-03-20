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
* @file EllipseUtils.h
* @class EllipseUtils
* @brief EllipseUtils
*
* EllipseUtils
* @author S. Riggi
* @date 05/08/2015
*/

#ifndef _ELLIPSE_UTILS_h
#define _ELLIPSE_UTILS_h 1

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TF1.h>
#include <TF12.h>
#include <TF2.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TGraph2D.h>
#include <TEllipse.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>


#include <deque>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>


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
#include <string>

namespace Caesar {

class EllipseUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    EllipseUtils();
		
		/**
		* \brief Class destructor: free allocated memory
		*/
    virtual ~EllipseUtils();

	public:
	
		//==========================================================================
		//== DEFINE PROGRAM CONSTANTS ==============================================
		//==========================================================================
		static const int NORMAL_TERMINATION= 0;
		static const int NO_INTERSECTION_POINTS= 100;
		static const int ONE_INTERSECTION_POINT= 101;
		static const int LINE_TANGENT_TO_ELLIPSE= 102;
		static const int DISJOINT_ELLIPSES= 103;
		static const int ELLIPSE2_OUTSIDETANGENT_ELLIPSE1= 104;
		static const int ELLIPSE2_INSIDETANGENT_ELLIPSE1= 105;
		static const int ELLIPSES_INTERSECT= 106;
		static const int TWO_INTERSECTION_POINTS= 107;
		static const int THREE_INTERSECTION_POINTS= 108;
		static const int FOUR_INTERSECTION_POINTS= 109;
		static const int ELLIPSE1_INSIDE_ELLIPSE2= 110;
		static const int ELLIPSE2_INSIDE_ELLIPSE1= 111;
		static const int ELLIPSES_ARE_IDENTICAL= 112;
		static const int INTERSECTION_POINT= 113;
		static const int TANGENT_POINT= 114;
		static const int ERROR_ELLIPSE_PARAMETERS= -100;
		static const int ERROR_DEGENERATE_ELLIPSE= -101;
		static const int ERROR_POINTS_NOT_ON_ELLIPSE= -102;
		static const int ERROR_INVERSE_TRIG= -103;
		static const int ERROR_LINE_POINTS= -104;
		static const int ERROR_QUARTIC_CASE= -105;
		static const int ERROR_POLYNOMIAL_DEGREE= -107;
		static const int ERROR_POLYNOMIAL_ROOTS= -108;
		static const int ERROR_INTERSECTION_PTS= -109;
		static const int ERROR_CALCULATIONS= -112;
	
		typedef boost::geometry::model::d2::point_xy<double,boost::geometry::cs::cartesian> point_2d;
		typedef boost::geometry::model::d2::point_xy<double> point_xy;
		typedef boost::geometry::model::polygon<point_2d> polygon_2d;

	public:

		/**
		* \brief Compute ellipse overlap
		*/
		//static double ComputeEllipseOverlap(TEllipse* ellipse1, TEllipse* ellipse2,TGraph& overlapAreaGraph,double& err,bool returnRelArea=true,int algoChoice=1);
		
		/**
		* \brief Compute overlap area polygon
		*/
		//static float getOverlapingAreaPoly(polygon_2d poly, polygon_2d poly2);

	
		static double ellipse_ellipse_overlap (double PHI_1, double A1, double B1, 
                                double H1, double K1, double PHI_2, 
                                double A2, double B2, double H2, double K2,
                                double X[4], double Y[4], int * NROOTS,
                                int *rtnCode, int choice);	
	
	private:
		
		/* Solve for real or complex roots of the quartic equation
 		* x^4 + a x^3 + b x^2 + c x + d = 0,
 		* returning the number of such roots.
 		* Roots are returned ordered.
 		*/
		static int gsl_poly_solve_quartic (double a, double b, double c, double d,double * x0, double * x1,double * x2, double * x3);
		static int gsl_poly_complex_solve_quartic(double a, double b, double c, double d,gsl_complex * z0, gsl_complex * z1,gsl_complex * z2, gsl_complex * z3);

		static int double_cmp(const void *a, const void *b);
		static void print_double_array(const double *array, size_t len);

		
		static double ellipse2tr (double x, double y, double AA, double BB, double CC, double DD, double EE, double FF);
		static int istanpt (double x, double y, double A1, double B1, double AA, double BB,double CC, double DD, double EE, double FF);
	
		static double nointpts (double A1, double B1, double A2, double B2, double H1, 
                 double K1, double H2, double K2, double PHI_1, double PHI_2,
                 double H2_TR, double K2_TR, double AA, double BB, 
                 double CC, double DD, double EE, double FF, int *rtnCode);

		static double twointpts (double x[], double y[], double A1, double B1, double PHI_1, 
                  double A2, double B2, double H2_TR, double K2_TR, 
                  double PHI_2, double AA, double BB, double CC, double DD, 
                  double EE, double FF, int *rtnCode);

		static double threeintpts (double xint[], double yint[], double A1, double B1, 
                    double PHI_1, double A2, double B2, double H2_TR, 
                    double K2_TR, double PHI_2, double AA, double BB, 
                    double CC, double DD, double EE, double FF,
                    int *rtnCode);

		static double fourintpts (double xint[], double yint[], double A1, double B1, 
                   double PHI_1, double A2, double B2, double H2_TR, 
                   double K2_TR, double PHI_2, double AA, double BB, 
                   double CC, double DD, double EE, double FF, int *rtnCode);

		static void QUADROOTS (double p[], double r[][5]);
		static void CUBICROOTS(double p[], double r[][5]);
		static void BIQUADROOTS(double p[],double r[][5]);

		//static int ellipse2poly(float PHI_1, float A1, float B1, float H1, float K1, polygon_2d * po,int n=20);

	private:
		
		//Some private methods
		static void SWAPD(double a,double b){
			do { 
				double tmp= b; 
				b= a; 
				a= tmp; 
			} 
			while(0);
		}//close SWAPD()

		static void SWAP(gsl_complex a,gsl_complex b){
			do { 
				gsl_complex tmp= b; 
				b= a; 
				a= tmp; 
			} 
			while(0);
		}//close SWAP()

		static double REAL(double* z,int i){
			return ((z)[2*(i)]);
		}

		static double IMAG(double* z,int i){
			return ((z)[2*(i)+1]);
		}

		//#define SWAPD(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)
		//#define SWAP(a,b) do { gsl_complex tmp = b ; b = a ; a = tmp ; } while(0)
		//#define REAL(z,i) ((z)[2*(i)])
		//#define IMAG(z,i) ((z)[2*(i)+1])
		//#define DEBUG 1
		//#define MIN(a,b) ((a<b)?a:b)
		//#define MAX(a,b) ((a<b)?b:a)

	private:

		//Some private consts
		static double EPS;
		static double pi;
		static double twopi;
		static int DEBUG;

	ClassDef(EllipseUtils,1)

};//close EllipseUtils class

#ifdef __MAKECINT__
#pragma link C++ class EllipseUtils+;
#endif

}//close namespace

#endif
