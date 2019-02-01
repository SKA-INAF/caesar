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
* @file EllipseUtils.cc
* @class EllipseUtils
* @brief EllipseUtils
*
* EllipseUtils class
* @author S. Riggi
* @date 20/01/2015
*/

#include <EllipseUtils.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>


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
#include <boost/scoped_array.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

using namespace std;
using namespace boost::geometry;


ClassImp(Caesar::EllipseUtils)


namespace Caesar {

double EllipseUtils::EPS= +1.0E-05;
double EllipseUtils::pi= (2.0*asin (1.0));	//-- a maximum-precision value of pi
double EllipseUtils::twopi= (2.0*pi);//-- a maximum-precision value of 2*pi
int EllipseUtils::DEBUG= 1;

EllipseUtils::EllipseUtils(){

}//close costructor


EllipseUtils::~EllipseUtils(){

}//close destructor

/*
double EllipseUtils::ComputeEllipseOverlap(TEllipse* ellipse1, TEllipse* ellipse2,TGraph& overlapAreaGraph,double& err,bool returnRelArea,int algoChoice){

	if(!ellipse1 || !ellipse2) return -1;
	
	//Get ellipse params
	double Theta_1= ellipse1->GetTheta();
	double R1_1= ellipse1->GetR1();
	double R2_1= ellipse1->GetR2(); 
	double Cx_1= ellipse1->GetX1();
	double Cy_1= ellipse1->GetY1();
	double Area_1= TMath::Pi()*R1_1*R2_1;

	double Theta_2= ellipse2->GetTheta();
	double R1_2= ellipse2->GetR1();
	double R2_2= ellipse2->GetR2(); 
	double Cx_2= ellipse2->GetX1();
	double Cy_2= ellipse2->GetY1(); 
	double Area_2= TMath::Pi()*R1_2*R2_2;

	double x[4], y[4];
  int nroots= 0;
	int rtn;
	
	double area= ellipse_ellipse_overlap (
			Theta_1, R1_1, R2_1, Cx_1, Cy_1, 
			Theta_2, R1_2, R2_2, Cx_2, Cy_2,
			x,y,&nroots, &rtn, algoChoice
	);

	//Convert ellipse in polygons and find polygon overlapping areas
	polygon_2d poly, poly2;
	ellipse2poly(Theta_1, R1_1, R2_1, Cx_1, Cy_1, &poly);
  ellipse2poly(Theta_2, R1_2, R2_2, Cx_2, Cy_2, &poly2);
	float areaPoly = getOverlapingAreaPoly(poly, poly2);

	//Estimate errors
	float eps = 0.0001;
	err = (fabs(areaPoly)<eps)?0.0:fabs(areaPoly-area)/areaPoly;

	//Fill area graph
	overlapAreaGraph.Set(0);
	for(int i=0; i<nroots; i++) {
		double theta= -Theta_1;	
		double Xroot= x[i]*cos(theta) + y[i]*sin(theta) + Cx_1;
		double Yroot= -x[i]*sin(theta) + y[i]*cos(theta) + Cy_1;
		overlapAreaGraph.SetPoint(i,Xroot,Yroot);
	}
  overlapAreaGraph.SetMarkerColor(kBlue-6);
	overlapAreaGraph.SetMarkerStyle(24);
	overlapAreaGraph.SetMarkerSize(1.3);
	overlapAreaGraph.SetFillColor(kBlue-6);
	overlapAreaGraph.SetLineColor(kBlue-6);
	
	double A= area;
	if(returnRelArea){
		A= area/Area_1;
	}

	return A;

}//close EllipseUtils::ComputeEllipseOverlap()
*/

//===========================================================================
//== ELLIPSE-ELLIPSE OVERLAP ================================================
//===========================================================================
double EllipseUtils::ellipse_ellipse_overlap (double PHI_1, double A1, double B1, 
                                double H1, double K1, double PHI_2, 
                                double A2, double B2, double H2, double K2,
                                double X[4], double Y[4], int * NROOTS,
                                int *rtnCode, int choice)
{

	//Convert theta to radians
	PHI_1*= TMath::DegToRad();
	PHI_2*= TMath::DegToRad();	

  //=======================================================================
  //== DEFINE LOCAL VARIABLES =============================================
  //=======================================================================
  int i, nroots, nychk, nintpts, fnRtnCode;
  double AA, BB, CC, DD, EE, FF, H2_TR, K2_TR, A22, B22, PHI_2R;
  double cosphi, cosphi2, sinphi, sinphi2, cosphisinphi;
  double tmp0, tmp1, tmp2, tmp3;
  double cy[5] = {0.0, 0.0, 0.0, 0.0, 0.0},  py[5] = {0.0, 0.0, 0.0, 0.0, 0.0}, r[3][5] ={{ 0.0, 0.0, 0.0, 0.0, 0.0}};
  double x1, x2;
  double ychk[4] = {0.0, 0.0, 0.0, 0.0};
  double xint[4], yint[4];
  double OverlapArea= -1;

  //==================================
  //====        DATA CHECK        ====
  //==================================
  // Each of the ellipse axis lengths must be positive
  if ( (!(A1 > EPS) || !(B1 > EPS)) || (!(A2 > EPS) || !(B2 > EPS)) ) {
  	(*rtnCode) = ERROR_ELLIPSE_PARAMETERS;
    return -1.0;
  }

  // The rotation angles should be between -2pi and 2pi (?)
  if ( (fabs (PHI_1) > (pi)) ) {
  	PHI_1 = fmod (PHI_1, pi);		
  }
  if ( (fabs (PHI_2) > (pi)) ) {	
    PHI_2 = fmod (PHI_2, pi);	
  }

  //=======================================================================
  //== DETERMINE THE TWO ELLIPSE EQUATIONS FROM INPUT PARAMETERS ==========
  //=======================================================================
  //-- Finding the points of intersection between two general ellipses
  //-- requires solving a quartic equation.  Before attempting to solve the
  //-- quartic, several quick tests can be used to eliminate some cases
  //-- where the ellipses do not intersect.  Optionally, can whittle away
  //-- at the problem, by addressing the easiest cases first.

  //-- Working with the translated+rotated ellipses simplifies the
  //-- calculations.  The ellipses are translated then rotated so that the
  //-- first ellipse is centered at the origin and oriented with the 
  //-- coordinate axes.  Then, the first ellipse will have the implicit 
  //-- (polynomial) form of
  //--   x^2/A1^2 + y+2/B1^2 = 1
	
  //-- For the second ellipse, the center is first translated by the amount
  //-- required to put the first ellipse at the origin, e.g., by (-H1, -K1)  
  //-- Then, the center of the second ellipse is rotated by the amount
  //-- required to orient the first ellipse with the coordinate axes, e.g.,
  //-- through the angle -PHI_1.
  //-- The translated and rotated center point coordinates for the second
  //-- ellipse are found with the rotation matrix, derivations are 
  //-- described in the reference.
  cosphi = cos (PHI_1);
  sinphi = sin (PHI_1);
  H2_TR = (H2 - H1)*cosphi + (K2 - K1)*sinphi;
  K2_TR = (H1 - H2)*sinphi + (K2 - K1)*cosphi;
  PHI_2R = PHI_2 - PHI_1;
  if ( (fabs (PHI_2R) > (twopi)) )
  	PHI_2R = fmod (PHI_2R, twopi);

  //-- Calculate implicit (Polynomial) coefficients for the second ellipse
  //-- in its translated-by (-H1, -H2) and rotated-by -PHI_1 postion
  //--       AA*x^2 + BB*x*y + CC*y^2 + DD*x + EE*y + FF = 0
  //-- Formulas derived in the reference
  //-- To speed things up, store multiply-used expressions first 
  cosphi = cos (PHI_2R);
  cosphi2 = cosphi*cosphi;
  sinphi = sin (PHI_2R);
  sinphi2 = sinphi*sinphi;
  cosphisinphi = 2.0*cosphi*sinphi;
  A22 = A2*A2;
  B22 = B2*B2;
  tmp0 = (cosphi*H2_TR + sinphi*K2_TR)/A22;
  tmp1 = (sinphi*H2_TR - cosphi*K2_TR)/B22;
  tmp2 = cosphi*H2_TR + sinphi*K2_TR;
  tmp3 = sinphi*H2_TR - cosphi*K2_TR;

  //-- implicit polynomial coefficients for the second ellipse
  AA = cosphi2/A22 + sinphi2/B22;
  BB = cosphisinphi/A22 - cosphisinphi/B22;
  CC = sinphi2/A22 + cosphi2/B22;
  DD = -2.0*cosphi*tmp0 - 2.0*sinphi*tmp1;
  EE = -2.0*sinphi*tmp0 + 2.0*cosphi*tmp1;
  FF = tmp2*tmp2/A22 + tmp3*tmp3/B22 - 1.0;  
  
	//=======================================================================
  //== CREATE AND SOLVE THE QUARTIC EQUATION TO FIND INTERSECTION POINTS ==
  //=======================================================================
  //-- If execution arrives here, the ellipses are at least 'close' to
  //-- intersecting.
  //-- Coefficients for the Quartic Polynomial in y are calculated from
  //-- the two implicit equations.
  //-- Formulas for these coefficients are derived in the reference.
  cy[4] = pow (A1, 4.0)*AA*AA + B1*B1*(A1*A1*(BB*BB - 2.0*AA*CC) + B1*B1*CC*CC);
  cy[3] = 2.0*B1*(B1*B1*CC*EE + A1*A1*(BB*DD - AA*EE));
 	cy[2] = A1*A1*((B1*B1*(2.0*AA*CC - BB*BB) + DD*DD  - 2.0*AA*FF)  - 2.0*A1*A1*AA*AA) + B1*B1*(2.0*CC*FF + EE*EE);
  cy[1] = 2.0*B1*(A1*A1*(AA*EE - BB*DD) + EE*FF);
  cy[0] = (A1*(A1*AA - DD) + FF)*(A1*(A1*AA + DD) + FF);

  //-- Once the coefficients for the Quartic Equation in y are known, the
  //-- roots of the quartic polynomial will represent y-values of the 
  //-- intersection points of the two ellipse curves.
  //-- The quartic sometimes degenerates into a polynomial of lesser 
  //-- degree, so handle all possible cases.

	
  if (fabs (cy[4]) > EPS){
  	//== QUARTIC COEFFICIENT NONZERO, USE QUARTIC FORMULA ===============
    for (i = 0; i <= 3; i++) py[4-i] = cy[i]/cy[4];
    py[0] = 1.0;
            
		
		double z[10]; //ret = s GSL_SUCCESS if all the roots are found and GSL_EFAILED if the QR reduction does not converge.
    gsl_complex  zz[4];
    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (5);//5 coeff
    if (choice==1){
    	gsl_poly_complex_solve (cy, 5, w, z);
      nroots=4;
    }
    else if(choice==2){
    	// zsolve_quartic.c - finds the complex roots of x^4 + a x^3 + b x^2 + c x + d = 0
      nroots = gsl_poly_complex_solve_quartic (py[1], py[2], py[3], py[4], &zz[0], &zz[1], &zz[2], &zz[3]);
    }

    for (i = 0; i < 4; i++){
    	if(choice==1){
      	r[1][i+1] = REAL(z,i); 
        r[2][i+1] = IMAG(z,i);
    	}
    	else if(choice==2) {
    		r[1][i+1] = GSL_REAL(zz[i]); 
        r[2][i+1] = GSL_IMAG(zz[i]);
      }
		}//end loop i<4

    if (choice==1) gsl_poly_complex_workspace_free(w); //free all the memory associated with the workspace w
	}//close if fabs

  else if (fabs (cy[3]) > EPS){ 
  	//== QUARTIC DEGENERATES TO CUBIC, USE CUBIC FORMULA ================
    for (i = 0; i <= 2; i++) py[3-i] = cy[i]/cy[3];
    py[0] = 1.0;
		
    gsl_complex Z0, Z1, Z2;
    nroots =  gsl_poly_complex_solve_cubic(py[1], py[2], py[3], &Z0, &Z1, &Z2);
    r[1][1] = GSL_REAL(Z0);
    r[2][1] = GSL_IMAG(Z0);
    r[1][2] = GSL_REAL(Z1);
    r[2][2] = GSL_IMAG(Z1);
    r[1][3] = GSL_REAL(Z2);
    r[2][3] = GSL_IMAG(Z2);
	}//close else if
     
	else if (fabs (cy[2]) > EPS) {
		//== QUARTIC DEGENERATES TO QUADRATIC, USE QUADRATIC FORMULA ========
    for (i = 0; i <= 1; i++) py[2-i] = cy[i]/cy[2];
    py[0] = 1.0;
		 
		gsl_complex z0, z1;
    nroots = gsl_poly_complex_solve_quadratic (py[0], py[1], py[2], &z0, &z1);
    r[1][1] = GSL_REAL(z0);
    r[2][1] = GSL_IMAG(z0);
    r[1][2] = GSL_REAL(z1);
    r[2][2] = GSL_IMAG(z1);
	}//close else if 

  else if (fabs (cy[1]) > EPS) {
  	//== QUARTIC DEGENERATES TO LINEAR: SOLVE DIRECTLY ==================
    //-- cy[1]*Y + cy[0] = 0
    r[1][1] = (-cy[0]/cy[1]);
    r[2][1] = 0.0;
   	nroots = 1;
	}//close else if
     
	else {
  	//== COMPLETELY DEGENERATE QUARTIC: ELLIPSES IDENTICAL??? ===========
    //-- a completely degenerate quartic, which would seem to
    //-- indicate that the ellipses are identical.  However, some
    //-- configurations lead to a degenerate quartic with no
    //-- points of intersection.
    nroots = 0;
	}

	//=======================================================================
	//== CHECK ROOTS OF THE QUARTIC: ARE THEY POINTS OF INTERSECTION? =======
	//=======================================================================
	//-- determine which roots are real, discard any complex roots
	nychk = 0;

	// GSL returns roots sorted in ascending order. We need descending sorted roots
  for (i = nroots; i >= 1; i--) {
 		if (fabs (r[2][i]) < EPS) {
    	nychk++;
      ychk[ nychk - 1 ] = r[1][i]*B1;
    }//close if fabs() 
	}//end loop for i

	//-- sort the real roots by straight insertion

	//-- determine whether polynomial roots are points of intersection for the two ellipses
  nintpts = 0;
  for (i = 0; i < nychk; i++) {
		//-- check for multiple roots
    if ( ( i < nychk -1 ) && fabs (ychk[i] - ychk[i+1]) < EPS/2.0 ) {
    	continue;
    }
          
		//-- check intersection points for ychk[i]
    if (fabs (ychk[i]) > B1) x1 = 0.0;
    else x1 = A1*sqrt (1.0 - (ychk[i]*ychk[i])/(B1*B1));
    x2 = -x1;
          
		if (fabs(ellipse2tr(x1, ychk[i], AA, BB, CC, DD, EE, FF)) < EPS) {
    	nintpts++;
      if (nintpts > 4) {
      	(*rtnCode) = ERROR_INTERSECTION_PTS;
        return -1.0;
      }
      
     	xint[nintpts-1] = x1;
      yint[nintpts-1] = ychk[i];
		}//close if fabs

    if ((fabs(ellipse2tr(x2, ychk[i], AA, BB, CC, DD, EE, FF)) < EPS) && (fabs (x2 - x1) > EPS)) {				
    	nintpts++;
      if (nintpts > 4) {
      	(*rtnCode) = ERROR_INTERSECTION_PTS;
        return -1.0;
      }
      
      xint[nintpts-1] = x2;
      yint[nintpts-1] = ychk[i];
		}//close if fabs
	}//end loop i<nychk
               
	//write intersection points
  for (i=0; i<nintpts;i++) {
  	X[i] = xint[i];
    Y[i] = yint[i];
  }
  *NROOTS = nintpts;

	//=======================================================================
	//== HANDLE ALL CASES FOR THE NUMBER OF INTERSCTION POINTS ==============
	//=========================================xs==============================
  switch (nintpts) {
  	case 0:
    case 1:
    	OverlapArea = nointpts (A1, B1, A2, B2, H1, K1, H2, K2, PHI_1, PHI_2, H2_TR, K2_TR, AA,BB, CC, DD, EE, FF, rtnCode);
      return OverlapArea;
			
    case 2:
    	//-- when there are two intersection points, it is possible for
      //-- them to both be tangents, in which case one of the ellipses
      //-- is fully contained within the other.  Check the points for
      //-- tangents; if one of the points is a tangent, then the other
      //-- must be as well, otherwise there would be more than 2 
      //-- intersection points.
      fnRtnCode = istanpt (xint[0], yint[0], A1, B1, AA, BB, CC, DD, EE, FF);
			if (fnRtnCode == TANGENT_POINT) {
      	OverlapArea = nointpts (A1, B1, A2, B2, H1, K1, H2, K2,  PHI_1, PHI_2, H2_TR, K2_TR, AA,BB, CC,  DD, EE, FF, rtnCode);
			}
      else {
      	OverlapArea = twointpts (xint, yint, A1, B1, PHI_1, A2, B2, H2_TR, K2_TR, PHI_2, AA, BB, CC, DD, EE, FF, rtnCode);
      }
      return OverlapArea;
			
		case 3:
    	//-- when there are three intersection points, one and only one
      //-- of the points must be a tangent point.
      OverlapArea = threeintpts (xint, yint,  A1, B1, PHI_1, A2, B2, H2_TR, K2_TR, PHI_2, AA, BB, CC, DD, EE, FF, rtnCode);
      return OverlapArea;
			
		case 4:
    	//-- four intersections points has only one case.
      OverlapArea = fourintpts (xint, yint,  A1, B1, PHI_1, A2, B2, H2_TR, K2_TR, PHI_2, AA, BB, CC, DD, EE, FF, rtnCode);
      return OverlapArea;
		
    default:
    	//-- should never get here (but get compiler warning for missing return value if this line is omitted)
      (*rtnCode) = ERROR_INTERSECTION_PTS;
			return -1.0;
	} //switch

}//close function


/*
//convert ellipse to polygon
int EllipseUtils::ellipse2poly(float PHI_1, float A1, float B1, float H1, float K1, polygon_2d * po,int n)
{
    //using namespace boost::geometry;
    polygon_2d  poly;
    
    float xc = H1;
    float yc = K1;
    float a = A1;
    float b = B1;
    float w = PHI_1;
    //	----
    if(!n)
    {
        std::cout << "error ellipse(): nshould be >0\n" <<std::endl;
        return 0;
    }
    float t = 0;
    int i = 0;
    double coor[n*2+1][2];
		//std::vector<point_xy> points; 
		const std::size_t points_size = n*2+1;
		boost::scoped_array<point_xy> points(new point_xy[points_size]); 

    double x, y;
    float step = pi/n;
    float sinphi = sin(w);
    float cosphi = cos(w);
    // std::cout << "step: " << step << "  sin=  " << sinphi << "  cos= " << cosphi << std::endl;
    for(i=0; i<2*n+1; i++)
    {   
        x = xc + a*cos(t)*cosphi - b*sin(t)*sinphi;
        y = yc + a*cos(t)*sinphi + b*sin(t)*cosphi;
        if(fabs(x) < 1e-4) x = 0;
        if(fabs(y) < 1e-4) y = 0;

        coor[i][0] = x;
        coor[i][1] = y;
				//points+= point_xy(x,y);
				points[i]= point_xy(x,y);
        t += step;
        // std::cout << "x=  " << x << " | y= " << y << std::endl;
    }

    //assign_points(poly, coor);
		assign_points(poly,std::make_pair(&points[0], &points[0] + points_size));

    correct(poly);
    *po = poly;
    return 1;
}
*/

double EllipseUtils::ellipse2tr (double x, double y, double AA, double BB, double CC, double DD, double EE, double FF){
	return (AA*x*x + BB*x*y + CC*y*y + DD*x + EE*y + FF);
}


double EllipseUtils::nointpts (double A1, double B1, double A2, double B2, double H1, 
                 double K1, double H2, double K2, double PHI_1, double PHI_2,
                 double H2_TR, double K2_TR, double AA, double BB, 
                 double CC, double DD, double EE, double FF, int *rtnCode)
{
     //some tmp-variables to avoid doing things several times.
     double A1B1 = A1*B1;
     double A2B2 = A2*B2;
     double Area_1 = pi*A1B1;
     double Area_2 = pi*A2B2;
     //-- The relative size of the two ellipses can be found from the axis
     //-- lengths 
     double relsize = A1B1 - A2B2; 
     if (relsize > 0.0)
     {
          //-- First Ellipse is larger than second ellipse.
          //-- If second ellipse center (H2_TR, K2_TR) is inside
          //-- first ellipse, then ellipse 2 is completely inside 
          //-- ellipse 1. Otherwise, the ellipses are disjoint.
          if ( ((H2_TR*H2_TR) / (A1*A1) 
                + (K2_TR*K2_TR) / (B1*B1)) < 1.0 )
          {
               (*rtnCode) = ELLIPSE2_INSIDE_ELLIPSE1;
               return Area_2;
          }
          else
          {
               (*rtnCode) = DISJOINT_ELLIPSES;
               return 0.0;
          }
     }
     else if (relsize < 0.0)
     {
          //-- Second Ellipse is larger than first ellipse
          //-- If first ellipse center (0, 0) is inside the
          //-- second ellipse, then ellipse 1 is completely inside
          //-- ellipse 2. Otherwise, the ellipses are disjoint
          //--   AA*x^2 + BB*x*y + CC*y^2 + DD*x + EE*y + FF = 0
          if (FF < 0.0)
          {
               (*rtnCode) = ELLIPSE1_INSIDE_ELLIPSE2;
               return Area_1;
          }
          else
          {
               (*rtnCode) = DISJOINT_ELLIPSES;
               return 0.0;
          }
     }
     else
     {
          //-- If execution arrives here, the relative sizes are identical.
          //-- Are the ellipses the same?  Check the parameters to see.
          //MC. Ellipses are the same if: H1=H2 And K1==K2 And Area_1 == Area_2
          //if ((((fabs (H1 - H2)) < EPS) && (fabs (K1 - K2) < EPS))
          //	  && (fabs (PHI_1 - PHI_2) < EPS))
          if( (fabs (H1 - H2) < EPS) && (fabs (K1 - K2) < EPS) && (fabs (Area_1 - Area_2) < EPS))
          {
               (*rtnCode) = ELLIPSES_ARE_IDENTICAL;
               return Area_1; 
          }
          else
          {
               //-- ellipses must be disjoint
               (*rtnCode) = DISJOINT_ELLIPSES;
               return 0.0;
          }
     }//-- end if (relsize > 0.0)
}

//-- two distinct intersection points (x1, y1) and (x2, y2) find overlap area
double EllipseUtils::twointpts (double x[], double y[], double A1, double B1, double PHI_1, 
                  double A2, double B2, double H2_TR, double K2_TR, 
                  double PHI_2, double AA, double BB, double CC, double DD, 
                  double EE, double FF, int *rtnCode)
{
     double area1, area2;
     double xmid, ymid, xmid_rt, ymid_rt;
     double theta1, theta2;
     double tmp, trsign;
     double x1_tr, y1_tr, x2_tr, y2_tr;
     //double discr;
     double cosphi, sinphi;

     //-- if execution arrives here, the intersection points are not
     //-- tangents.
	
     //-- determine which direction to integrate in the ellipse_segment
     //-- routine for each ellipse.

     //-- find the parametric angles for each point on ellipse 1
     if (fabs (x[0]) > A1)
          x[0] = (x[0] < 0) ? -A1 : A1;
     if (y[0] < 0.0) 	 //-- Quadrant III or IV
          theta1 = twopi - acos (x[0] / A1);
     else             //-- Quadrant I or II      
          theta1 = acos (x[0] / A1);
		
     if (fabs (x[1]) > A1)
          x[1] = (x[1] < 0) ? -A1 : A1;
     if (y[1] < 0.0) 	 //-- Quadrant III or IV
          theta2 = twopi - acos (x[1] / A1);
     else             //-- Quadrant I or II      
          theta2 = acos (x[1] / A1);

     //-- logic is for proceeding counterclockwise from theta1 to theta2
     if (theta1 > theta2)
     {
          tmp = theta1;
          theta1 = theta2;
          theta2 = tmp;
     }

     //-- find a point on the first ellipse that is different than the two
     //-- intersection points.
     xmid = A1*cos ((theta1 + theta2)/2.0);	
     ymid = B1*sin ((theta1 + theta2)/2.0);	
	
     //-- the point (xmid, ymid) is on the first ellipse 'between' the two
     //-- intersection points (x[1], y[1]) and (x[2], y[2]) when travelling 
     //-- counter- clockwise from (x[1], y[1]) to (x[2], y[2]).  If the point
     //-- (xmid, ymid) is inside the second ellipse, then the desired segment
     //-- of ellipse 1 contains the point (xmid, ymid), so integrate 
     //-- counterclockwise from (x[1], y[1]) to (x[2], y[2]).  Otherwise, 
     //-- integrate counterclockwise from (x[2], y[2]) to (x[1], y[1])
     if (ellipse2tr (xmid, ymid, AA, BB, CC, DD, EE, FF) > 0.0)
     {
          tmp = theta1;
          theta1 = theta2;
          theta2 = tmp;
     }

     //-- here is the ellipse segment routine for the first ellipse
     if (theta1 > theta2)
          theta1 -= twopi;
     if ((theta2 - theta1) > pi)
          trsign = 1.0;
     else
          trsign = -1.0;
     /* area1 = 0.5*(A1*B1*(theta2 - theta1)  */
     /* 			 + trsign*fabs (x[1]*y[2] - x[2]*y[1]));	 */

     area1 = 0.5*(A1*B1*(theta2 - theta1) 
                  + trsign*fabs (x[0]*y[1] - x[1]*y[0]));	
	
     if (area1 < 0)
     {
          printf("TWO area1=%f\n",area1);
          area1 += A1*B1;
          getc(stdin);
     }
     //-- find ellipse 2 segment area.  The ellipse segment routine
     //-- needs an ellipse that is centered at the origin and oriented
     //-- with the coordinate axes.  The intersection points (x[1], y[1]) and
     //-- (x[2], y[2]) are found with both ellipses translated and rotated by
     //-- (-H1, -K1) and -PHI_1.  Further translate and rotate the points
     //-- to put the second ellipse at the origin and oriented with the
     //-- coordinate axes.  The translation is (-H2_TR, -K2_TR), and the
     //-- rotation is -(PHI_2 - PHI_1) = PHI_1 - PHI_2
     cosphi = cos (PHI_1 - PHI_2);
     sinphi = sin (PHI_1 - PHI_2);
     x1_tr = (x[0] - H2_TR)*cosphi + (y[0] - K2_TR)*-sinphi;
     y1_tr = (x[0] - H2_TR)*sinphi + (y[0] - K2_TR)*cosphi;
     x2_tr = (x[1] - H2_TR)*cosphi + (y[1] - K2_TR)*-sinphi;
     y2_tr = (x[1] - H2_TR)*sinphi + (y[1] - K2_TR)*cosphi;
	
     //-- determine which branch of the ellipse to integrate by finding a
     //-- point on the second ellipse, and asking whether it is inside the
     //-- first ellipse (in their once-translated+rotated positions)
     //-- find the parametric angles for each point on ellipse 1
     if (fabs (x1_tr) > A2)
          x1_tr = (x1_tr < 0) ? -A2 : A2;
     if (y1_tr < 0.0) 	 //-- Quadrant III or IV
          theta1 = twopi - acos (x1_tr/A2);
     else             //-- Quadrant I or II      
          theta1 = acos (x1_tr/A2);
		
     if (fabs (x2_tr) > A2)
          x2_tr = (x2_tr < 0) ? -A2 : A2;
     if (y2_tr < 0.0) 	 //-- Quadrant III or IV
          theta2 = twopi - acos (x2_tr/A2);
     else             //-- Quadrant I or II      
          theta2 = acos (x2_tr/A2);

     //-- logic is for proceeding counterclockwise from theta1 to theta2
     if (theta1 > theta2)
     {
          tmp = theta1;
          theta1 = theta2;
          theta2 = tmp;
     }

     //-- find a point on the second ellipse that is different than the two
     //-- intersection points.
     xmid = A2*cos ((theta1 + theta2)/2.0);	
     ymid = B2*sin ((theta1 + theta2)/2.0);
	
     //-- translate the point back to the second ellipse in its once-
     //-- translated+rotated position
     cosphi = cos (PHI_2 - PHI_1);
     sinphi = sin (PHI_2 - PHI_1);
     xmid_rt = xmid*cosphi + ymid*-sinphi + H2_TR;
     ymid_rt = xmid*sinphi + ymid*cosphi + K2_TR;

     //-- the point (xmid_rt, ymid_rt) is on the second ellipse 'between' the
     //-- intersection points (x[1], y[1]) and (x[2], y[2]) when travelling
     //-- counterclockwise from (x[1], y[1]) to (x[2], y[2]).  If the point
     //-- (xmid_rt, ymid_rt) is inside the first ellipse, then the desired 
     //-- segment of ellipse 2 contains the point (xmid_rt, ymid_rt), so 
     //-- integrate counterclockwise from (x[1], y[1]) to (x[2], y[2]).  
     //-- Otherwise, integrate counterclockwise from (x[2], y[2]) to 
     //-- (x[1], y[1])
     if (((xmid_rt*xmid_rt)/(A1*A1) + (ymid_rt*ymid_rt)/(B1*B1)) > 1.0)
     {
          tmp = theta1;
          theta1 = theta2;
          theta2 = tmp;
     }

     //-- here is the ellipse segment routine for the second ellipse
     if (theta1 > theta2)
          theta1 -= twopi;
     if ((theta2 - theta1) > pi)
          trsign = 1.0;
     else
          trsign = -1.0;
     area2 = 0.5*(A2*B2*(theta2 - theta1) 
                  + trsign*fabs (x1_tr*y2_tr - x2_tr*y1_tr));	
	
     if (area2 < 0)
     {
#if DEBUG
          printf("TWO area2=%f\n",area2);
#endif
          area2 += A2*B2;
          
     }
	
     (*rtnCode) = TWO_INTERSECTION_POINTS;
#if DEBUG
     printf("Twointpts: \t area1 =%f,  area2=%f\n",area1, area2);
#endif
     return area1 + area2;
}//close function


//-- three distinct intersection points, must have two intersections
//-- and one tangent, which is the only possibility
double EllipseUtils::threeintpts (double xint[], double yint[], double A1, double B1, 
                    double PHI_1, double A2, double B2, double H2_TR, 
                    double K2_TR, double PHI_2, double AA, double BB, 
                    double CC, double DD, double EE, double FF,
                    int *rtnCode)
{
     int i, tanpts, tanindex, fnRtn;
     double OverlapArea;

     //-- need to determine which point is a tangent, and which two points
     //-- are intersections
     tanpts = 0;
     for (i = 0; i < 3; i++)
     {
          fnRtn = istanpt (xint[i], yint[i], A1, B1, AA, BB, CC, DD, EE, FF);

          if (fnRtn == TANGENT_POINT)
          {
               tanpts++;
               tanindex = i;
          }
     }
#if DEBUG
     printf("tanindex=%d\n",tanindex);
#endif
     //-- there MUST be 2 intersection points and only one tangent
     if (tanpts != 1)
     {
          //-- should never get here unless there is a problem discerning
          //-- whether or not a point is a tangent or intersection
          (*rtnCode) = ERROR_INTERSECTION_PTS;
          return -1.0;
     }
	
     //-- store the two interesection points into (x[1], y[1]) and 
     //-- (x[2], y[2])
     switch (tanindex)
     {
     case 0:
          xint[0] = xint[2];
          yint[0] = yint[2];
			
          break;

     case 1:
          xint[1] = xint[2];
          yint[1] = yint[2];
			
          break;

     case 2:
          //-- intersection points are already in the right places
          break;
     }

     OverlapArea = twointpts (xint, yint, A1, B1, PHI_1, A2, B2, H2_TR, K2_TR,
                              PHI_2, AA, BB, CC, DD, EE, FF, rtnCode);
     (*rtnCode) = THREE_INTERSECTION_POINTS;
     return OverlapArea;
}//close function

//-- four intersection points
double EllipseUtils::fourintpts (double xint[], double yint[], double A1, double B1, 
                   double PHI_1, double A2, double B2, double H2_TR, 
                   double K2_TR, double PHI_2, double AA, double BB, 
                   double CC, double DD, double EE, double FF, int *rtnCode)
{
     int i, j, k;
     double xmid, ymid, xint_tr[4], yint_tr[4], OverlapArea;
     double theta[4], theta_tr[4], cosphi, sinphi, tmp0, tmp1, tmp2;
     double area1=0, area2=0, area3=0, area4=0, area5=0;


     //some tmp-variables to avoid calculating the same thing several times.
     double A1B1 = A1*B1;
     double A2B2 = A2*B2;
     double Area_1 = pi*A1B1;
     double Area_2 = pi*A2B2;

	
     //-- only one case, which involves two segments from each ellipse, plus
     //-- two triangles.
     //-- get the parametric angles along the first ellipse for each of the
     //-- intersection points
//		for (i = 1; i <= 4; i++)
     for (i = 0; i <= 3; i++)
     {
          if (fabs (xint[i]) > A1)
               xint[i] = (xint[i] < 0) ? -A1 : A1;
          if (yint[i] < 0.0) 	 //-- Quadrant III or IV
               theta[i] = twopi - acos (xint[i] / A1);
          else             //-- Quadrant I or II      
               theta[i] = acos (xint[i] / A1);
     }
		
     //-- sort the angles by straight insertion, and put the points in 
     //-- counter-clockwise order
#if DEBUG
     for (int k=0; k<=3; k++)
     {
          printf("k=%d:  Theta = %f, xint=%f, yint=%f\n",k,theta[k], xint[k], yint[k]);
     }
#endif
     for (j = 1; j <= 3; j++)
     {
          tmp0 = theta[j];
          tmp1 = xint[j];
          tmp2 = yint[j];
		
          for (k = j - 1; k >= 0; k--)
          {
               if (theta[k] <= tmp0)
                    break;
			
               theta[k+1] = theta[k];
               xint[k+1] = xint[k];
               yint[k+1] = yint[k];
          }
		
          theta[k+1] = tmp0;
          xint[k+1] = tmp1;
          yint[k+1] = tmp2;
     }
#if DEBUG
     printf("AFTER sorting\n");
     for (int k=0; k<=3; k++)
     {
          printf("k=%d:  Theta = %f, xint=%f, yint=%f\n",k,theta[k], xint[k], yint[k]);
     }    
#endif
	
     //-- find the area of the interior quadrilateral
     /* area1 = 0.5*fabs ((xint[3] - xint[1])*(yint[4] - yint[2]) */
     /* 				  - (xint[4] - xint[2])*(yint[3] - yint[1])); */
     area1 = 0.5*fabs ((xint[2] - xint[0])*(yint[3] - yint[1])
                       - (xint[3] - xint[1])*(yint[2] - yint[0]));
	
     //-- the intersection points lie on the second ellipse in its once
     //-- translated+rotated position.  The segment algorithm is implemented
     //-- for an ellipse that is centered at the origin, and oriented with
     //-- the coordinate axes; so, in order to use the segment algorithm
     //-- with the second ellipse, the intersection points must be further
     //-- translated+rotated by amounts that put the second ellipse centered
     //-- at the origin and oriented with the coordinate axes.
     cosphi = cos (PHI_1 - PHI_2);
     sinphi = sin (PHI_1 - PHI_2);
     for (i = 0; i <= 3; i++)
     {
          xint_tr[i] = (xint[i] - H2_TR)*cosphi + (yint[i] - K2_TR)*-sinphi;
          yint_tr[i] = (xint[i] - H2_TR)*sinphi + (yint[i] - K2_TR)*cosphi;
		
          if (fabs (xint_tr[i]) > A2)
               xint_tr[i] = (xint_tr[i] < 0) ? -A2 : A2;
          if (yint_tr[i] < 0.0) 	 //-- Quadrant III or IV
               theta_tr[i] = twopi - acos (xint_tr[i]/A2);
          else             //-- Quadrant I or II      
               theta_tr[i] = acos (xint_tr[i]/A2);
     }

     //-- get the area of the two segments on ellipse 1
     xmid = A1*cos ((theta[0] + theta[1])/2.0);	
     ymid = B1*sin ((theta[0] + theta[1])/2.0);
     //-- the point (xmid, ymid) is on the first ellipse 'between' the two
     //-- sorted intersection points (xint[1], yint[1]) and (xint[2], yint[2])
     //-- when travelling counter- clockwise from (xint[1], yint[1]) to 
     //-- (xint[2], yint[2]).  If the point (xmid, ymid) is inside the second 
     //-- ellipse, then one desired segment of ellipse 1 contains the point 
     //-- (xmid, ymid), so integrate counterclockwise from (xint[1], yint[1])
     //-- to (xint[2], yint[2]) for the first segment, and from 
     //-- (xint[3], yint[3] to (xint[4], yint[4]) for the second segment.
     if (ellipse2tr (xmid, ymid, AA, BB, CC, DD, EE, FF) < 0.0)
     {
          area2 = 0.5*(A1B1*(theta[1] - theta[0])
                       - fabs (xint[0]*yint[1] - xint[1]*yint[0]));
	    
          area2 = 0.5*(A1B1*(theta[1] - theta[0])
                       - fabs (xint[0]*yint[1] - xint[1]*yint[0]));

          area4 = 0.5*(A2B2*(theta_tr[2] - theta_tr[1]) - fabs (xint_tr[1]*yint_tr[2] - xint_tr[2]*yint_tr[1]) );
			
          if (theta_tr[3] > theta_tr[0])
               area5 = 0.5*(A2B2*(theta_tr[0] - (theta_tr[3] - twopi))
                            - fabs (xint_tr[3]*yint_tr[0] - xint_tr[0]*yint_tr[3]));
          else
               area5 = 0.5*(A2B2*(theta_tr[0] - theta_tr[3])
                            - fabs (xint_tr[3]*yint_tr[0] - xint_tr[0]*yint_tr[3]));
     }
     else
     {
          area2 = 0.5*(A1B1*(theta[2] - theta[1])
                       - fabs (xint[1]*yint[2] - xint[2]*yint[1]));
          area3 = 0.5*(A1B1*(theta[0] - (theta[3] - twopi))
                       - fabs (xint[3]*yint[0] - xint[0]*yint[3]));
          area4 = 0.5*(A2B2*(theta_tr[1] - theta_tr[0])
                       - fabs (xint_tr[0]*yint_tr[1] - xint_tr[1]*yint_tr[0]));
          area5 = 0.5*(A2B2*(theta_tr[3] - theta_tr[2])
                       - fabs (xint_tr[2]*yint_tr[3] - xint_tr[3]*yint_tr[2]));
     }
     if(area5<0)
     {
#if DEBUG
          printf("\n\t\t-------------> area5 is negativ (%f). Add: pi*A2*B2=%f <------------\n",area5, Area_2);
#endif
          area5 += Area_2;
     }
     if(area4<0)
     {
#if DEBUG
          printf("\n\t\t-------------> area4 is negativ (%f). Add: pi*A2*B2=%f <------------\n",area4, Area_2);
#endif
          area4 += Area_2;
     }
     if(area3<0)
     {
#if DEBUG
          printf("\n\t\t-------------> area3 is negativ (%f). Add: pi*A2*B2=%f <------------\n",area3, Area_1);
#endif
          area3 += Area_1;
     }
     if(area2<0)
     {
#if DEBUG
          printf("\n\t\t-------------> area2 is negativ (%f). Add: pi*A2*B2=%f <------------\n",area2, Area_1);
#endif
          area2 += Area_1;
     }
		
#if DEBUG
     printf("\narea1=%f, area2=%f area3=%f, area4=%f, area5=%f\n\n",area1, area2, area3, area4, area5);
#endif
     OverlapArea = area1 + area2 + area3 + area4 + area5;
	
     (*rtnCode) = FOUR_INTERSECTION_POINTS;
     return OverlapArea;
}//close function

//-- check whether an intersection point is a tangent or a cross-point
int EllipseUtils::istanpt (double x, double y, double A1, double B1, double AA, double BB,
             double CC, double DD, double EE, double FF)
{
     double x1, y1, x2, y2, theta, test1, test2, eps_radian;

     //-- Avoid inverse trig calculation errors: there could be an error 
     //-- if |x1/A| > 1.0 when calling acos().  If execution arrives here, 
     //-- then the point is on the ellipse within EPS.
     if (fabs (x) > A1)
          x = (x < 0) ? -A1 : A1;

     //-- Calculate the parametric angle on the ellipse for (x, y)
     //-- The parametric angles depend on the quadrant where each point
     //-- is located.  See Table 1 in the reference.
     if (y < 0.0) 	 //-- Quadrant III or IV
          theta = twopi - acos (x / A1);
     else             //-- Quadrant I or II      
          theta = acos (x / A1);

     //-- determine the distance from the origin to the point (x, y)
     /* branch = sqrt (x*x + y*y); */

     /* //-- use the distance to find a small angle, such that the distance */
     /* //-- along ellipse 1 is approximately 2*EPS */
     /* if (branch < 100.0*EPS) */
     /* 	eps_radian = 2.0*EPS; */
     /* else */
     /* 	eps_radian = asin (2.0*EPS/branch); */
	
     //fix 24.11.12
     eps_radian = 0.1; //arbitrary value
	
     //-- determine two points that are on each side of (x, y) and lie on
     //-- the first ellipse
     x1 = A1*cos (theta + eps_radian);
     y1 = B1*sin (theta + eps_radian);
     x2 = A1*cos (theta - eps_radian);
     y2 = B1*sin (theta - eps_radian);
	
     //-- evaluate the two adjacent points in the second ellipse equation
     test1 = ellipse2tr (x1, y1, AA, BB, CC, DD, EE, FF);
     test2 = ellipse2tr (x2, y2, AA, BB, CC, DD, EE, FF);

     //-- if the ellipses are tangent at the intersection point, then
     //-- points on both sides will either both be inside ellipse 1, or
     //-- they will both be outside ellipse 1
#if DEBUG
     printf("\t\t--- debug istanpt with (x,y)=(%f, %f), A1=%f, B1=%f\n", x, y, A1, B1);
     printf("theta=%f\n", theta);
     printf("eps_Radian=%f\n", eps_radian);
     printf("(x1, y1)=(%f, %f)\n", x1, y1);
     printf("(x2, y2)=(%f, %f)\n", x2, y2);
     printf("test1=%f\n", test1);
     printf("test2=%f\n", test2);
#endif

     if ((test1*test2) > 0.0)
          return TANGENT_POINT;
     else
          return INTERSECTION_POINT;
}//close function



//===========================================================================
//-- CACM Algorithm 326: Roots of low order polynomials.
//-- Nonweiler, Terence R.F., CACM Algorithm 326: Roots of low order 
//-- polynomials, Communications of the ACM, vol. 11 no. 4, pages 
//-- 269-270 (1968). Translated into c and programmed by M. Dow, ANUSF,
//-- Australian National University, Canberra, Australia.
//-- Accessed at http://www.netlib.org/toms/326.
//-- Modified to void functions, integers replaced with floating point
//-- where appropriate, some other slight modifications for readability
//-- and debugging ease.
//===========================================================================
void EllipseUtils::QUADROOTS (double p[], double r[][5])
{
     /*
      Array r[3][5]  p[5]
      Roots of poly p[0]*x^2 + p[1]*x + p[2]=0
      x=r[1][k] + i r[2][k]  k=1,2
    */
     double b,c,d;
     b=-p[1]/(2.0*p[0]);
     c=p[2]/p[0];
     d=b*b-c;
     if(d>=0.0)
     {
          if(b>0.0)
               b=(r[1][2]=(sqrt(d)+b));
          else
               b=(r[1][2]=(-sqrt(d)+b));
          r[1][1]=c/b;
          r[2][1]=(r[2][2]=0.0);
     }
     else
     {
          d=(r[2][1]=sqrt(-d));
          r[2][2]=-d;
          r[1][1]=(r[1][2]=b);
     }
     return;
}//close function

void EllipseUtils::CUBICROOTS(double p[], double r[][5])
{
     /*
      Array r[3][5]  p[5]
      Roots of poly p[0]*x^3 + p[1]*x^2 + p[2]*x + p[3] = 0
      x=r[1][k] + i r[2][k]  k=1,...,3
      Assumes 0<arctan(x)<pi/2 for x>0
    */
     double s,t,b,c,d;
     int k;
     if(p[0]!=1.0)
     {
          for(k=1;k<4;k++)
               p[k]=p[k]/p[0];
          p[0]=1.0;
     }
     s=p[1]/3.0;
     t=s*p[1];
     b=0.5*(s*(t/1.5-p[2])+p[3]);
     t=(t-p[2])/3.0;
     c=t*t*t;
     d=b*b-c;
     if(d>=0.0)
     {
          d=pow((sqrt(d)+fabs(b)),1.0/3.0);
          if(d!=0.0)
          {
               if(b>0.0)
                    b=-d;
               else
                    b=d;
               c=t/b;
          }
          d=r[2][2]=sqrt(0.75)*(b-c);
          b=b+c;
          c=r[1][2]=-0.5*b-s;
          if((b>0.0 && s<=0.0) || (b<0.0 && s>0.0))
          {
               r[1][1]=c;
               r[2][1]=-d;
               r[1][3]=b-s;
               r[2][3]=0.0;
          }
          else
          {
               r[1][1]=b-s;
               r[2][1]=0.0;
               r[1][3]=c;
               r[2][3]=-d;
          }
     }  /* end 2 equal or complex roots */
     else
     {
          if(b==0.0)
               d=atan(1.0)/1.5;
          else
               d=atan(sqrt(-d)/fabs(b))/3.0;
          if(b<0.0)
               b=2.0*sqrt(t);
          else
               b=-2.0*sqrt(t);
          c=cos(d)*b;
          t=-sqrt(0.75)*sin(d)*b-0.5*c;
          d=-t-c-s;
          c=c-s;
          t=t-s;
      
          if(fabs(c)>fabs(t))
          {
               r[1][3]=c;
          }
          else
          {
               r[1][3]=t;
               t=c;
          }
          if(fabs(d)>fabs(t))
          {
               r[1][2]=d;
          }
          else
          {
               r[1][2]=t;
               t=d;
          }
          r[1][1]=t;
          for(k=1;k<4;k++)
               r[2][k]=0.0;
     }
     return;
}//close function

void EllipseUtils::BIQUADROOTS(double p[],double r[][5])
{
     /*
      Array r[3][5]  p[5]
      Roots of poly p[0]*x^4 + p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4] = 0
      x=r[1][k] + i r[2][k]  k=1,...,4
    */
     double a,b,c,d,e;
     int k,j;
     if(p[0] != 1.0)
     {
          for(k=1;k<5;k++)
               p[k]=p[k]/p[0];
          p[0]=1.0;
     }
     e=0.25*p[1];
     b=2.0*e;
     c=b*b;
     d=0.75*c;
     b=p[3]+b*(c-p[2]);
     a=p[2]-d;
     c=p[4]+e*(e*a-p[3]);
     a=a-d;
     p[1]=0.5*a;
     p[2]=(p[1]*p[1]-c)*0.25;
     p[3]=b*b/(-64.0);
     if(p[3]<0.0)
     {
          CUBICROOTS(p,r);
          for(k=1;k<4;k++)
          {
               if(r[2][k]==0.0 && r[1][k]>0.0)
               {
                    d=r[1][k]*4.0;
                    a=a+d;
                    if(a>=0.0 && b>=0.0)
                         p[1]=sqrt(d);
                    else if(a<=0.0 && b<=0.0)
                         p[1]=sqrt(d);
                    else
                         p[1]=-sqrt(d);
                    b=0.5*(a+b/p[1]);
                    goto QUAD;
               }
          }
     }
     if(p[2]<0.0)
     {
          b=sqrt(c);
          d=b+b-a;
          p[1]=0.0;
          if(d>0.0)
               p[1]=sqrt(d);
     }
     else
     {
          if(p[1]>0.0)
               b=sqrt(p[2])*2.0+p[1];
          else
               b=-sqrt(p[2])*2.0+p[1];
          if(b!=0.0)
          {
               p[1]=0.0;
          }
          else
          {
               for(k=1;k<5;k++)
               {
                    r[1][k]=-e;
                    r[2][k]=0.0;
               }
               goto END;
          }
     }
QUAD:

     p[2]=c/b;
     QUADROOTS(p,r);
     for(k=1;k<3;k++)
          for(j=1;j<3;j++)
               r[j][k+2]=r[j][k];
     p[1]=-p[1];
     p[2]=b;
     QUADROOTS(p,r);
     for(k=1;k<5;k++)
     {
          r[1][k]=r[1][k]-e;
     }
END:
     for(k=1;k<5;k++)
          return;
}//close function


int EllipseUtils::double_cmp(const void *a, const void *b) { 
     const double *ia = (const double *)a; // casting pointer types 
     const double *ib = (const double *)b;
     return *ia  - *ib; 
     /* double comparison: returns negative if b > a 
       and positive if a > b */ 
} 

void EllipseUtils::print_double_array(const double *array, size_t len) 
{ 
     size_t i;
 
     for(i=0; i<len; i++) 
          printf("%f | ", array[i]);
 
     putchar('\n');
} 



int EllipseUtils::gsl_poly_solve_quartic (double a, double b, double c, double d,
                        double *x0, double *x1, double *x2, double *x3)
{
  /* 
   * This code is based on a simplification of
   * the algorithm from zsolve_quartic.c for real roots
   */
  double u[3];
  double aa, pp, qq, rr, rc, sc, tc, mt;
  double w1r, w1i, w2r, w2i, w3r;
  double v[3], v1, v2, arg, theta;
  double disc, h;
  int k1= 0; 
	int k2= 0;
  double zarr[4];

  /* Deal easily with the cases where the quartic is degenerate. The
   * ordering of solutions is done explicitly. */
  if (0 == b && 0 == c)
    {
      if (0 == d)
        {
          if (a > 0)
            {
              *x0 = -a;
              *x1 = 0.0;
              *x2 = 0.0;
              *x3 = 0.0;
            }
          else
            {
              *x0 = 0.0;
              *x1 = 0.0;
              *x2 = 0.0;
              *x3 = -a;
            }
          return 4;
        }
      else if (0 == a)
        {
          if (d > 0)
            {
              return 0;
            }
          else
            {
              *x1 = sqrt (sqrt (-d));
              *x0 = -(*x1);
              return 2;
            }
        }
    }

  if (0.0 == c && 0.0 == d)
    {
      *x0=0.0;
      *x1=0.0;
      if (gsl_poly_solve_quadratic(1.0,a,b,x2,x3)==0) {
	mt=3;
      } else {
	mt=1;
      }
    }
  else 
    {
      /* For non-degenerate solutions, proceed by constructing and
       * solving the resolvent cubic */
      aa = a * a;
      pp = b - (3.0/8.0) * aa;
      qq = c - (1.0/2.0) * a * (b - (1.0/4.0) * aa);
      rr = d - (1.0/4.0) * (a * c - (1.0/4.0) * aa * (b - (3.0/16.0) * aa));
      rc = (1.0/2.0) * pp;
      sc = (1.0/4.0) * ((1.0/4.0) * pp * pp - rr);
      tc = -((1.0/8.0) * qq * (1.0/8.0) * qq);

      /* This code solves the resolvent cubic in a convenient fashion
       * for this implementation of the quartic. If there are three real
       * roots, then they are placed directly into u[].  If two are
       * complex, then the real root is put into u[0] and the real
       * and imaginary part of the complex roots are placed into
       * u[1] and u[2], respectively. Additionally, this
       * calculates the discriminant of the cubic and puts it into the
       * variable disc. */
      {
	double qcub = (rc * rc - 3 * sc);
	double rcub = (2 * rc * rc * rc - 9 * rc * sc + 27 * tc);

	double Q = qcub / 9;
	double R = rcub / 54;

	double Q3 = Q * Q * Q;
	double R2 = R * R;

	double CR2 = 729 * rcub * rcub;
	double CQ3 = 2916 * qcub * qcub * qcub;

	disc = (CR2 - CQ3) / 2125764.0;

	if (0 == R && 0 == Q)
	  {
	    u[0] = -rc / 3;
	    u[1] = -rc / 3;
	    u[2] = -rc / 3;
	  }
	else if (CR2 == CQ3)
	  {
	    double sqrtQ = sqrt (Q);
	    if (R > 0)
	      {
		u[0] = -2 * sqrtQ - rc / 3;
		u[1] = sqrtQ - rc / 3;
		u[2] = sqrtQ - rc / 3;
	      }
	    else
	      {
		u[0] = -sqrtQ - rc / 3;
		u[1] = -sqrtQ - rc / 3;
		u[2] = 2 * sqrtQ - rc / 3;
	      }
	  }
	else if (CR2 < CQ3)
	  {
	    double sqrtQ = sqrt (Q);
	    double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
	    double theta = acos (R / sqrtQ3);
	    if (R / sqrtQ3 >= 1.0) theta = 0.0;
	    {
	      double norm = -2 * sqrtQ;
	  
	      u[0] = norm * cos (theta / 3) - rc / 3;
	      u[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - rc / 3;
	      u[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - rc / 3;
	    }
	  }
	else
	  {
	    double sgnR = (R >= 0 ? 1 : -1);
	    double modR = fabs (R);
	    double sqrt_disc = sqrt (R2 - Q3);
	    double A = -sgnR * pow (modR + sqrt_disc, 1.0 / 3.0);
	    double B = Q / A;
	    double mod_diffAB = fabs (A - B);

	    u[0] = A + B - rc / 3;
	    u[1] = -0.5 * (A + B) - rc / 3;
	    u[2] = -(sqrt (3.0) / 2.0) * mod_diffAB;
	  }
      }
      /* End of solution to resolvent cubic */

      /* Combine the square roots of the roots of the cubic 
       * resolvent appropriately. Also, calculate 'mt' which 
       * designates the nature of the roots:
       * mt=1 : 4 real roots (disc == 0)
       * mt=2 : 0 real roots (disc < 0)
       * mt=3 : 2 real roots (disc > 0)
       */

      if (0.0 == disc) 
	u[2] = u[1];

      if (0 >= disc)
	{
	  mt = 2; 

	  /* One would think that we could return 0 here and exit,
	   * since mt=2. However, this assignment is temporary and
	   * changes to mt=1 under certain conditions below.  
	   */
	  
	  v[0] = fabs (u[0]);
	  v[1] = fabs (u[1]);
	  v[2] = fabs (u[2]);

	  v1 = GSL_MAX (GSL_MAX (v[0], v[1]), v[2]);
	  /* Work out which two roots have the largest moduli */
	  k1 = 0, k2 = 0;
	  if (v1 == v[0])
	    {
	      k1 = 0;
	      v2 = GSL_MAX (v[1], v[2]);
	    }
	  else if (v1 == v[1])
	    {
	      k1 = 1;
	      v2 = GSL_MAX (v[0], v[2]);
	    }
	  else
	    {
	      k1 = 2;
	      v2 = GSL_MAX (v[0], v[1]);
	    }

	  if (v2 == v[0])
	    {
	      k2 = 0;
	    }
	  else if (v2 == v[1])
	    {
	      k2 = 1;
	    }
	  else
	    {
	      k2 = 2;
	    }
	  
	  if (0.0 <= u[k1]) 
	    {
	      w1r=sqrt(u[k1]);
	      w1i=0.0;
	    } 
	  else 
	    {
	      w1r=0.0;
	      w1i=sqrt(-u[k1]);
	    }
	  if (0.0 <= u[k2]) 
	    {
	      w2r=sqrt(u[k2]);
	      w2i=0.0;
	    } 
	  else 
	    {
	      w2r=0.0;
	      w2i=sqrt(-u[k2]);
	    }
	}
      else
	{
	  mt = 3;

	  if (0.0 == u[1] && 0.0 == u[2]) 
	    {
	      arg = 0.0;
	    } 
	  else 
	    {
	      arg = sqrt(sqrt(u[1] * u[1] + u[2] * u[2]));
	    }
	  theta = atan2(u[2], u[1]);
	  
	  w1r = arg * cos(theta / 2.0);
	  w1i = arg * sin(theta / 2.0);
	  w2r = w1r;
	  w2i = -w1i;
	}
  
      /* Solve the quadratic to obtain the roots to the quartic */
      w3r = qq / 8.0 * (w1i * w2i - w1r * w2r) / 
	(w1i * w1i + w1r * w1r) / (w2i * w2i + w2r * w2r);
      h = a / 4.0;

      zarr[0] = w1r + w2r + w3r - h;
      zarr[1] = -w1r - w2r + w3r - h;
      zarr[2] = -w1r + w2r - w3r - h;
      zarr[3] = w1r - w2r - w3r - h;
      
      /* Arrange the roots into the variables z0, z1, z2, z3 */
      if (2 == mt)
        {
          if (u[k1] >= 0 && u[k2] >= 0)
            {
              mt = 1;
	      *x0 = zarr[0];
	      *x1 = zarr[1];
	      *x2 = zarr[2];
	      *x3 = zarr[3];
            }
	  else
	    {
	      return 0;
	    }
	}
      else 
        {
	  *x0 = zarr[0];
	  *x1 = zarr[1];
        }
    }
  
  /* Sort the roots as usual */
  if (1 == mt)
    {
      /* Roots are all real, sort them by the real part */
      if (*x0 > *x1)
        SWAPD (*x0, *x1);
      if (*x0 > *x2)
        SWAPD (*x0, *x2);
      if (*x0 > *x3)
        SWAPD (*x0, *x3);

      if (*x1 > *x2)
        SWAPD (*x1, *x2);
      if (*x2 > *x3)
        {
          SWAPD (*x2, *x3);
          if (*x1 > *x2)
            SWAPD (*x1, *x2);
        }
      return 4;
    }
  else
    {
      /* 2 real roots */
      if (*x0 > *x1)
        SWAPD (*x0, *x1);
    }

  return 2;
}//close function


int EllipseUtils::gsl_poly_complex_solve_quartic (double a, double b, double c, double d,
                                gsl_complex * z0, gsl_complex * z1,
                                gsl_complex * z2, gsl_complex * z3)
{
  gsl_complex i, zarr[4], w1, w2, w3;
  double r4 = 1.0 / 4.0;
  double q2 = 1.0 / 2.0, q4 = 1.0 / 4.0, q8 = 1.0 / 8.0;
  double q1 = 3.0 / 8.0, q3 = 3.0 / 16.0;
  double u[3], v[3], v1, v2, disc;
  double aa, pp, qq, rr, rc, sc, tc, h;
  int k1 = 0, k2 = 0, mt;

  GSL_SET_COMPLEX (&i, 0.0, 1.0);
  GSL_SET_COMPLEX (&zarr[0], 0.0, 0.0);
  GSL_SET_COMPLEX (&zarr[1], 0.0, 0.0);
  GSL_SET_COMPLEX (&zarr[2], 0.0, 0.0);
  GSL_SET_COMPLEX (&zarr[3], 0.0, 0.0);
  GSL_SET_COMPLEX (&w1, 0.0, 0.0);
  GSL_SET_COMPLEX (&w2, 0.0, 0.0);
  GSL_SET_COMPLEX (&w3, 0.0, 0.0);

  /* Deal easily with the cases where the quartic is degenerate. The
   * ordering of solutions is done explicitly. */
  if (0 == b && 0 == c)
    {
      if (0 == d)
        {
          if (a > 0)
            {
              GSL_SET_COMPLEX (z0, -a, 0.0);
              GSL_SET_COMPLEX (z1, 0.0, 0.0);
              GSL_SET_COMPLEX (z2, 0.0, 0.0);
              GSL_SET_COMPLEX (z3, 0.0, 0.0);
            }
          else
            {
              GSL_SET_COMPLEX (z0, 0.0, 0.0);
              GSL_SET_COMPLEX (z1, 0.0, 0.0);
              GSL_SET_COMPLEX (z2, 0.0, 0.0);
              GSL_SET_COMPLEX (z3, -a, 0.0);
            }
          return 4;
        }
      else if (0 == a)
        {
          if (d > 0)
            {
              double sqrt_d = sqrt (d);
              gsl_complex i_sqrt_d = gsl_complex_mul_real (i, sqrt_d);
              gsl_complex minus_i = gsl_complex_conjugate (i);
              *z3 = gsl_complex_sqrt (i_sqrt_d);
              *z2 = gsl_complex_mul (minus_i, *z3);
              *z1 = gsl_complex_negative (*z2);
              *z0 = gsl_complex_negative (*z3);
            }
          else
            {
              double sqrt_abs_d = sqrt (-d);
              *z3 = gsl_complex_sqrt_real (sqrt_abs_d);
              *z2 = gsl_complex_mul (i, *z3);
              *z1 = gsl_complex_negative (*z2);
              *z0 = gsl_complex_negative (*z3);
            }
          return 4;
        }
    }

  if (0.0 == c && 0.0 == d)
    {
      disc = (a * a - 4.0 * b);
      if (disc < 0.0)
        {
	  /* 2 complex roots and 2 real roots for the quartic
	   */
          mt = 3;
        }
      else
        {
	  /* 4 real roots for the quartic
	   */
          mt = 1;
        }
      *z0 = zarr[0];
      *z1 = zarr[0];
      gsl_poly_complex_solve_quadratic (1.0, a, b, z2, z3);
    }
  else
    {
      /* For non-degenerate solutions, proceed by constructing and
       * solving the resolvent cubic */
      aa = a * a;
      pp = b - q1 * aa;
      qq = c - q2 * a * (b - q4 * aa);
      rr = d - q4 * (a * c - q4 * aa * (b - q3 * aa));
      rc = q2 * pp;
      sc = q4 * (q4 * pp * pp - rr);
      tc = -(q8 * qq * q8 * qq);

      /* This code solves the resolvent cubic in a convenient fashion
       * for this implementation of the quartic. If there are three real
       * roots, then they are placed directly into u[].  If two are
       * complex, then the real root is put into u[0] and the real
       * and imaginary part of the complex roots are placed into
       * u[1] and u[2], respectively. Additionally, this
       * calculates the discriminant of the cubic and puts it into the
       * variable disc. */
      {
        double qcub = (rc * rc - 3 * sc);
        double rcub = (2 * rc * rc * rc - 9 * rc * sc + 27 * tc);

        double Q = qcub / 9;
        double R = rcub / 54;

        double Q3 = Q * Q * Q;
        double R2 = R * R;

        double CR2 = 729 * rcub * rcub;
        double CQ3 = 2916 * qcub * qcub * qcub;

        disc = (CR2 - CQ3) / 2125764.0;
	
        if (0 == R && 0 == Q)
          {
            u[0] = -rc / 3;
            u[1] = -rc / 3;
            u[2] = -rc / 3;
          }
        else if (CR2 == CQ3)
          {
            double sqrtQ = sqrt (Q);
            if (R > 0)
              {
                u[0] = -2 * sqrtQ - rc / 3;
                u[1] = sqrtQ - rc / 3;
                u[2] = sqrtQ - rc / 3;
              }
            else
              {
                u[0] = -sqrtQ - rc / 3;
                u[1] = -sqrtQ - rc / 3;
                u[2] = 2 * sqrtQ - rc / 3;
              }
          }
        else if (CR2 < CQ3)
          {
            double sqrtQ = sqrt (Q);
            double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
            double theta = acos (R / sqrtQ3);
            if (R / sqrtQ3 >= 1.0)
              theta = 0.0;
	    {
	      double norm = -2 * sqrtQ;
	      
	      u[0] = norm * cos (theta / 3) - rc / 3;
	      u[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - rc / 3;
	      u[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - rc / 3;
	    }
          }
        else
          {
            double sgnR = (R >= 0 ? 1 : -1);
            double modR = fabs (R);
            double sqrt_disc = sqrt (R2 - Q3);
            double A = -sgnR * pow (modR + sqrt_disc, 1.0 / 3.0);
            double B = Q / A;

            double mod_diffAB = fabs (A - B);
            u[0] = A + B - rc / 3;
            u[1] = -0.5 * (A + B) - rc / 3;
            u[2] = -(sqrt (3.0) / 2.0) * mod_diffAB;
          }
      }
      /* End of solution to resolvent cubic */

      /* Combine the square roots of the roots of the cubic 
       * resolvent appropriately. Also, calculate 'mt' which 
       * designates the nature of the roots:
       * mt=1 : 4 real roots 
       * mt=2 : 0 real roots 
       * mt=3 : 2 real roots 
       */

      if (0 == disc)
        {
          u[2] = u[1];
        }
      if (0 >= disc)
        {
          mt = 2;
          v[0] = fabs (u[0]);
          v[1] = fabs (u[1]);
          v[2] = fabs (u[2]);

          v1 = GSL_MAX (GSL_MAX (v[0], v[1]), v[2]);
          if (v1 == v[0])
            {
              k1 = 0;
              v2 = GSL_MAX (v[1], v[2]);
            }
          else if (v1 == v[1])
            {
              k1 = 1;
              v2 = GSL_MAX (v[0], v[2]);
            }
          else
            {
              k1 = 2;
              v2 = GSL_MAX (v[0], v[1]);
            }

          if (v2 == v[0])
            {
              k2 = 0;
            }
          else if (v2 == v[1])
            {
              k2 = 1;
            }
          else
            {
              k2 = 2;
            }
          w1 = gsl_complex_sqrt_real (u[k1]);
          w2 = gsl_complex_sqrt_real (u[k2]);
        }
      else
        {
          mt = 3;
          GSL_SET_COMPLEX (&w1, u[1], u[2]);
          GSL_SET_COMPLEX (&w2, u[1], -u[2]);
          w1 = gsl_complex_sqrt (w1);
          w2 = gsl_complex_sqrt (w2);
        }

      /* Solve the quadratic in order to obtain the roots 
       * to the quartic */
      if (0.0 != gsl_complex_abs (gsl_complex_mul (w1, w2))) {
	w3 = gsl_complex_mul_real (gsl_complex_inverse 
				   (gsl_complex_mul (w1, w2)), -qq / 8.0);
      }
      h = r4 * a;
      zarr[0] = gsl_complex_add_real (gsl_complex_add 
				      (gsl_complex_add (w1, w2), w3), -h);
      zarr[1] = gsl_complex_add_real (gsl_complex_add 
				      (gsl_complex_negative 
				       (gsl_complex_add (w1, w2)), w3), -h);
      zarr[2] = gsl_complex_add_real (gsl_complex_sub 
				      (gsl_complex_sub (w2, w1), w3), -h);
      zarr[3] = gsl_complex_add_real (gsl_complex_sub 
				      (gsl_complex_sub (w1, w2), w3), -h);

      /* Arrange the roots into the variables z0, z1, z2, z3 */
      if (2 == mt)
        {
          if (u[k1] >= 0 && u[k2] >= 0)
            {
              mt = 1;
              GSL_SET_COMPLEX (z0, GSL_REAL (zarr[0]), 0.0);
              GSL_SET_COMPLEX (z1, GSL_REAL (zarr[1]), 0.0);
              GSL_SET_COMPLEX (z2, GSL_REAL (zarr[2]), 0.0);
              GSL_SET_COMPLEX (z3, GSL_REAL (zarr[3]), 0.0);
            }
          else if (u[k1] >= 0 && u[k2] < 0)
            {
              *z0 = zarr[0];
              *z1 = zarr[3];
              *z2 = zarr[2];
              *z3 = zarr[1];
            }
          else if (u[k1] < 0 && u[k2] >= 0)
            {
              *z0 = zarr[0];
              *z1 = zarr[2];
              *z2 = zarr[3];
              *z3 = zarr[1];
            }
          else if (u[k1] < 0 && u[k2] < 0)
            {
              *z0 = zarr[0];
              *z1 = zarr[1];
              *z2 = zarr[3];
              *z3 = zarr[2];
            }
        }
      else if (3 == mt)
        {
          GSL_SET_COMPLEX (z0, GSL_REAL (zarr[0]), 0.0);
          GSL_SET_COMPLEX (z1, GSL_REAL (zarr[1]), 0.0);
          *z2 = zarr[3];
          *z3 = zarr[2];
        }
    }

  /* 
   * Sort the roots as usual.
   * This code is most likely not optimal. 
   */

  if (1 == mt)
    {
      /* Roots are all real, sort them by the real part */

      if (GSL_REAL (*z0) > GSL_REAL (*z1)) SWAP (*z0, *z1);
      if (GSL_REAL (*z0) > GSL_REAL (*z2)) SWAP (*z0, *z2);
      if (GSL_REAL (*z0) > GSL_REAL (*z3)) SWAP (*z0, *z3);

      if (GSL_REAL (*z1) > GSL_REAL (*z2)) SWAP (*z1, *z2);
      if (GSL_REAL (*z2) > GSL_REAL (*z3))
        {
          SWAP (*z2, *z3);
          if (GSL_REAL (*z1) > GSL_REAL (*z2)) SWAP (*z1, *z2);
        }
    }
  else if (2 == mt)
    {
      /* Roots are all complex. z0 and z1 are conjugates
       * and z2 and z3 are conjugates. */

      /* If all of the real parts are equal, just sort
	 by the imaginary parts */
      if (GSL_REAL (*z0) == GSL_REAL (*z2)) 
	{
	
	  /* Ensure that the pairs are ordered so that the
	     root with negative imaginary part is first 
	  */
	  if (GSL_IMAG (*z2) > GSL_IMAG (*z3)) SWAP (*z2, *z3);
	  if (GSL_IMAG (*z0) > GSL_IMAG (*z1)) SWAP (*z0, *z1);
	
	  if (GSL_IMAG (*z0) < GSL_IMAG (*z2)) 
	    {
	      SWAP (*z1, *z2);
	      SWAP (*z2, *z3);
	    } 
	  else 
	    {
	      SWAP (*z0, *z2);
	      SWAP (*z0, *z1);
	    }
	
	} 
      else 
	{
	  /* Otherwise, sort the real parts first */
	  if (GSL_REAL (*z0) > GSL_REAL (*z2))
	    {
	      SWAP (*z0, *z2);
	      SWAP (*z1, *z3);
	    }
	  /* Then sort by the imaginary parts */
	  if (GSL_IMAG (*z0) > GSL_IMAG (*z1)) SWAP (*z0, *z1);
	  if (GSL_IMAG (*z2) > GSL_IMAG (*z3)) SWAP (*z2, *z3);
	}
    }
  else
    {
      /* 2 real roots. z2 and z3 are conjugates. */

      /* Swap complex roots, if necessary. */
      if (GSL_IMAG (*z2) > GSL_IMAG (*z3)) SWAP (*z2, *z3);

      /* Sort real parts */
      if (GSL_REAL (*z0) == GSL_REAL (*z2)) 
	{
	  if (GSL_REAL (*z0) < GSL_REAL(*z1)) 
	    {
	      SWAP (*z1, *z3);
	      SWAP (*z1, *z2);
	      SWAP (*z0, *z1);
	    }
	  else if (GSL_REAL (*z0) == GSL_REAL (*z1)) 
	    {
	      SWAP (*z0,*z2);
	    }
	  else 
	    {
	      SWAP (*z0, *z1);
	      SWAP (*z1, *z2);
	    }
	} 
      else if (GSL_REAL (*z1) == GSL_REAL (*z2)) 
	{
	  if (GSL_REAL (*z0) < GSL_REAL(*z1)) 
	    {
	      SWAP (*z1, *z2);
	    }
	  else 
	    {
	      SWAP (*z0, *z3);
	      SWAP (*z0, *z2);
	    }
	} 
      else 
	{
	  if (GSL_REAL (*z0) > GSL_REAL (*z1)) SWAP (*z0, *z1);
	  if (GSL_REAL (*z1) > GSL_REAL (*z2))
	    {
	      if (GSL_REAL (*z0) > GSL_REAL (*z2))
		{
		  SWAP (*z0, *z2);
		  SWAP (*z1, *z3);
		}
	      else
		{
		  SWAP (*z1, *z2);
		  SWAP (*z2, *z3);
		}
	    }
	}
    }

  return 4;
}//close gsl_poly_complex_solve_quartic()


}//close namespace
