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

//OpenCV
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//ROOT headers
#include <TObject.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TGraph.h>
#include <TEllipse.h>

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
		* \brief Compute ellipse area
		*/
		static double ComputeEllipseArea(TEllipse* ellipse){
			if(!ellipse) return -1;
			double r1= ellipse->GetR1();
			double r2= ellipse->GetR2();
			double A= TMath::Pi()*r1*r2;
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
		static int ComputeContourOverlapArea(double& overlapArea,Contour* contour,Contour* contour2,Contour* overlapContour=0);

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

		
	private:
	
		ClassDef(MathUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class MathUtils+;
#endif	

}//close namespace


#endif 
