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
* @file MathUtils.h
* @class MathUtils
* @brief Utility functions for math tasks
*
* Utility functions for math tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef _MATH_UTILS_h
#define _MATH_UTILS_h 1

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

//OpenCV
/*
#if OPENCV_MAJOR_VERSION < 4
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
#endif
*/
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//ROOT headers
#include <TObject.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TEllipse.h>
#include <TVector2.h>

//Boost headers
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

//C++ headers
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

#include <complex>

using namespace std;

namespace Caesar {

class Contour;

class MathUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    MathUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~MathUtils();

		
		typedef boost::geometry::model::d2::point_xy<double,boost::geometry::cs::cartesian> point_2d;
		typedef boost::geometry::model::d2::point_xy<double> point_xy;
		typedef boost::geometry::model::polygon<point_2d> polygon_2d;	

	public:

		/**
		* \brief Return 2D grid partition given Nx x Ny pixels and box sizes
		*/
		static int Compute2DGrid(
			std::vector<long int>& ix_min,std::vector<long int>& ix_max,
			std::vector<long int>& iy_min,std::vector<long int>& iy_max,
			long int Nx,long int Ny,long int boxSizeX,long int boxSizeY,float gridStepSizeX,float gridStepSizeY
		);

		/**
		* \brief Return 2D grid partition given x & y range and step sizes
		*/
		static int Compute2DFloatGrid(
			std::vector<float>& ix_min,std::vector<float>& ix_max,
			std::vector<float>& iy_min,std::vector<float>& iy_max,
			float xmin,float xmax,float xstep,float ymin,float ymax,float ystep
		);

		/**
		* \brief Return 2D grid axis bin given grid parameters
		*/
		static int FindGrid2DAxisBin(float x,long int nx,float xmin,float xmax,float xstep);

		/**
		* \brief Return 2D grid global bin given grid parameters
		*/
		static long int FindGrid2DBin(float x,float y,long int nx,float xmin,float xmax,float xstep,long int ny,float ymin,float ymax,float ystep);

		/**
		* \brief Perform bilinear interpolation on regular grid
		*/
		static int BiLinearInterpolation(
			std::vector<double>const& sampled_gridX,std::vector<double>const& sampled_gridY,
			std::vector<double>const& sampledZ,
			std::vector<double>const& interp_gridX,std::vector<double>const& interp_gridY,
			std::vector<double>& interpZ
		);

		/**
		* \brief Get opencv mat convolutions
		*/
		static cv::Mat GetConvolution(cv::Mat I, cv::Mat kernel);

		/**
		* \brief Get convolution (2nd version)
		*/
		static cv::Mat GetConvolution2(cv::Mat I, cv::Mat kernel);

		/**
		* \brief Get atrous convolution
		*/
		static cv::Mat GetATrousConvolution(cv::Mat I, cv::Mat kernel,int scale);

		/**
		* \brief Get mirror index
		*/
		static int GetMirrorIndex(int index,int N);

		/**
		* \brief Compute matrix trace
		*/
		static double GetMatrixTrace(TMatrixD* T){			
			double trace= 0;
			for(int i=0;i<T->GetNrows();i++){
				trace+= (*T)(i,i);		
			}
  		return trace;
		}

		/**
		* \brief Compute DFT shifted
		*/
		static std::vector< std::complex<double> > DFTShifted(std::vector< std::complex<double> > data, int n);

		/**
		* \brief Compute DFT
		*/
		static std::vector< std::complex<double> > DFT(std::vector< std::complex<double> > data, int n);

		/**
		* \brief Compute IDFT
		*/
		static std::vector< std::complex<double> > IDFT(std::vector< std::complex<double> > data, int n);

		/**
		* \brief Eta function definition
		*/
		static int EtaAuxiliaryFcn(int s,int N){
			int thr= -floor(N/2.) + N-1;
			int fval= 0;
			if(s<=thr) fval= s;
			else fval= N-s;
			return fval;
		}

		/**
		* \brief Compute contout curvature
		*/
		static std::vector<double> GetContourCurvature(std::vector< std::complex<double> > data,double sigma);
		
		/**
		* \brief Compute factorial
		*/
		static long int factorial(int n) {
  		return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
		}

		/**
		* \brief Compute 2d gaussian integral
		*/
		static double Compute2DGausIntegral(double A,double sigmaX,double sigmaY){
			double fluxDensity= 2*TMath::Pi()*A*sigmaX*sigmaY;
			return fluxDensity;
		}

		/**
		* \brief Compute 2d gaussian ellipse integral
		*/
		static double Compute2DGausEllipseIntegral(double A,double Bmaj,double Bmin){
			double fluxDensity= TMath::Pi()*A*Bmaj*Bmin/(4.*log(2));
			return fluxDensity;
		}

		/**
		* \brief Compute ellipse area from ROOT TEllipse
		*/
		static double ComputeEllipseArea(TEllipse* ellipse)
		{
			if(!ellipse) return -1;
			double r1= ellipse->GetR1();
			double r2= ellipse->GetR2();
			double A= TMath::Pi()*r1*r2;
			return A;
		}

		/**
		* \brief Compute ellipse area
		*/
		static double ComputeEllipseArea(double bmaj,double bmin)
		{
			double a= bmaj/2.;//semi-major axis
			double b= bmin/2.;//semi-minor axis
			double A= TMath::Pi()*a*b;
			return A;
		}
		
		/**
		* \brief Compute overlap area between two ellipses
		*/
		static int ComputeEllipseOverlapArea(double& overlapArea,double& err,int& rtn,TEllipse* ellipse1, TEllipse* ellipse2,int method=1,Contour* overlapContour=0);

		/**
		* \brief Compute polygon area
		*/
		static double ComputePolygonArea(polygon_2d& poly);

		/**
		* \brief Compute contour area
		*/
		static int ComputeContourArea(double& area,Contour* contour);

		/**
		* \brief Compute contour overlap area
		*/
		static int ComputeContourOverlapArea(double& overlapArea,int& overlapFlag,Contour* contour,Contour* contour2,Contour* overlapContour=0);

		/**
		* \brief Compute overlap area between two polygons (in BOOST format)
		*/
		static int ComputePolygonOverlapArea(double& overlapArea,polygon_2d& overlap_poly,polygon_2d& poly, polygon_2d& poly2);

		/**
		* \brief Convert contour to polygon (in BOOST format)
		*/
		static int Contour2Polygon(polygon_2d& poly,Contour* contour);

		/**
		* \brief Convert an ellipse to polygon (in BOOST format)
		*/
		static int Ellipse2Polygon(polygon_2d& poly,double Cx, double Cy, double R1, double R2, double theta, int n=20);

		/**
		* \brief Convert a TEllipse to polygon (in BOOST format)
		*/
		static int Ellipse2Polygon(polygon_2d& poly,TEllipse* ellipse, int n);

		/**
		* \brief Compute ellipse eccentricity from ROOT TEllipse
		*/
		static double ComputeEllipseEccentricity(TEllipse* ellipse);
		/**
		* \brief Compute ellipse eccentricity from (bmaj, bmin)
		*/
		static double ComputeEllipseEccentricity(double bmaj,double bmin);

		/**
		* \brief Modulus operator
		*/
		template<typename T>
		static T Mod(T x, T y)
		{
			//Check if given type is not a float/double
			if(std::numeric_limits<T>::is_exact){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Mod function requires float/double arguments, will return 0!");
				#endif
				return 0;
			}

    	if (0 == y) return x;

    	double m= x - y * std::floor(x/y);

    	// handle boundary cases resulting from floating-point limited accuracy:
    	if (y > 0)              // modulo range: [0..y)
    	{
      	if (m>=y)           // Mod(-1e-16             , 360.    ): m= 360.
        	return 0;

        if (m<0 )
        {
        	if (y+m == y)
          	return 0  ; // just in case...
          else
            return y+m; // Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14 
        }
    	}//close if y>0
    	else                    // modulo range: (y..0]
    	{
      	if (m<=y)           // Mod(1e-16              , -360.   ): m= -360.
        	return 0;

        if (m>0 )
        {
        	if (y+m == y)
          	return 0  ; // just in case...
          else
            return y+m; // Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14 
        }
    	}//close else

    	return m;
 		}//close Mod function

		/**
		* \brief Limit angle in range [-theta_limit,theta_limit]
		*/
		static double GetAngleInRange(double theta,double theta_limit)
		{
			while (theta <= -theta_limit) theta+= 2*theta_limit;
  		while (theta > theta_limit) theta-= 2*theta_limit;
  		return theta;
		}

		/**
		* \brief Evaluate 2d elliptical gaussian at given (x,y) (NB: theta in radians)
		*/
		static double EvalGaus2D(double X,double Y,double A,double X0,double Y0,double sigmaX,double sigmaY,double theta)
		{
			double cost2= cos(theta)*cos(theta);
			double sint2= sin(theta)*sin(theta);
			double sin2t = sin(2. * theta);
			double xstd2= sigmaX*sigmaX;
			double ystd2= sigmaY*sigmaY;
			double xdiff = X - X0;
  		double ydiff = Y - Y0;
			double a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2));
  		double b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2));
  		double c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2));
			double fcn= A * exp(-((a * xdiff*xdiff) + (b * xdiff * ydiff) + (c * ydiff*ydiff)));
	
			return fcn;
		}//close EvalGaus2D()

		/**
		* \brief Check if 2D point with coordinates (x,y) is inside a polygon
		*/
		static bool IsPointInsidePolygon(double x,double y,const std::vector<TVector2>& polygon);
		
		/**
		* \brief Check if 2D point with coordinates (x,y) is inside a polygon
		*/
		static bool IsPointInsidePolygon(const TVector2& point,const std::vector<TVector2>& polygon)
		{
			return IsPointInsidePolygon(point.X(),point.Y(),polygon);
		}

		/**
		* \brief Check if 2D point is lying on segment a-b including vertex
		*/
		static bool IsPointInsideSegment(TVector2 point, TVector2 a, TVector2 b, double epsilon=1.e-12);

		/**
		* \brief Check if 2D point is lying on polygon boundary, including vertices
		*/
		static bool IsPointOnPolygonBoundary(const TVector2& point, const std::vector<TVector2>& polygon);

		/**
		* \brief Check if 2D point is lying on polygon boundary, including vertices
		*/
		static bool IsPointOnPolygonBoundary(double x, double y, const std::vector<TVector2>& polygon)
		{
			TVector2 point(x, y);
			return IsPointOnPolygonBoundary(point, polygon);
		}

		/**
		* \brief Check if 2D point with coordinates (x,y) is inside a polygon (version 2)
		*/
		/*
		static bool IsPointInsidePolygonV2(const TVector2& point,const std::vector<TVector2>& polygon)
		{
			return IsPointInsidePolygonV2(point.X(),point.Y(),polygon);
		} 
		*/
		
		/**
		* \brief Check if 2D point with coordinates (x,y) is inside a polygon (version 2)
		*/
		//static bool IsPointInsidePolygonV2(double x,double y,const std::vector<TVector2>& polygon, bool include_borders=false);

		
		/**
 		* @brief Check if a point lies inside, on or outside any polygon.
 		*
 		* Winding number algorithm can be used to check if any point lies inside a
 		* polygon. A more detailed explanation can be found in the blog post. The link
 		* is attached at the top of the file.
		*
 		*
 		* @param query_point Point to check.
 		* @param vertices Vertices making up the polygon in anticlockwise direction.
 		* @return  = 1: query_point lies inside the polygon.
 		*          = 0: query_point lies on the polygon.
 		*          =-1: query_point lies outside the polygon.
 		*/
		//static int query_point_inside_polygon(const TVector2& query_point, const std::vector<TVector2>& vertices);


		/**
 		* @brief The result can be used to test if the query point lies on the left or
 		*        right side of the line formed by pt1 and pt2 when viewed in
 		*        anticlockwise  direction.
 		*
 		* @param pt1: First point to form equation of line.
 		* @param pt2: Second point to form equation of line.
 		* @param query_point: Query point
 		* @return: > 0: Query point lies on left of the line.
 		*          = 0: Query point lies on the line.
 		*          < 0: Query point lies on right of the line.
 		*/	
		/*
		static double substitute_point_in_line(const TVector2& pt1, const TVector2& pt2, const TVector2& query_point) 
		{
    	return ((query_point.Y() - pt1.Y()) * (pt2.X() - pt1.X())) -
           	 ((query_point.X() - pt1.X()) * (pt2.Y() - pt1.Y()));
		};
		*/

		/**
		* \brief Compute rotated coordinates (x,y)
		*/
		static void ComputeRotatedCoords(double& xrot,double& yrot,double x,double y,double cx,double cy,double theta);
		
		/**
		* \brief Compute Euclidean distance between two points
		*/
		static double GetEuclideanDist(double x1,double y1,double x2,double y2)
		{
			double d= sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
			return d;
		}

		/**
		* \brief SynchrotronSelfAbsSED model
		*/
		static double SynchrotronSelfAbsSED(double* x,double* pars);
		/**
		* \brief SynchrotronExtFreeFreeAbsSED model
		*/
		static double SynchrotronExtFreeFreeAbsSED(double* x,double* pars);
		/**
		* \brief SynchrotronIntFreeFreeAbsSED model
		*/
		static double SynchrotronIntFreeFreeAbsSED(double* x,double* pars);
		/**
		* \brief FreeFreeSED model
		*/
		static double FreeFreeSED(double* x,double* pars);

	private:
	
		ClassDef(MathUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class MathUtils+;
#endif	

}//close namespace


#endif 
