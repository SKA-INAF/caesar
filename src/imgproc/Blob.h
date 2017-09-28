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
* @file Blob.h
* @class Blob
* @brief Blob class
*
* Class representing an image blob
* @author S. Riggi
* @date 20/01/2015
*/

#ifndef _BLOB_h
#define _BLOB_h 1

#include <Pixel.h>
#include <Image.h>
#include <Consts.h>


#include <TObject.h>
#include <TMatrixD.h>
#include <TColor.h>

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
#include <iostream>
#include <time.h>
#include <ctime>

using namespace std;

namespace Caesar {

class Contour;


struct MatchPixelType {
	MatchPixelType(const int& type) : m_type(type) {}
 	bool operator()(const Pixel* obj) const {
  	return obj->type == m_type;
 	}
 	private:
  	const int& m_type;
};


class Blob : public TNamed {

	public:
		/** 
		\brief Class constructor: initialize structures.
 		*/
		Blob();
		/** 
		\brief Parametric constructor
 		*/
		Blob(std::string name);			
		
		/** 
		\brief Parametric constructor
 		*/
		Blob(std::vector<Pixel*>const& pixels,std::string name="");			

		/**
		* \brief Copy constructor
		*/
		Blob(const Blob& blob);
		/**
		* \brief Class destructor: free allocated memory
		*/
		virtual ~Blob();
		/**
		* \brief Assignment Operator
		*/
		Blob& operator=(const Blob &blob);
		/**
		* \brief Copy method
		*/
		void Copy(TObject& blob) const;


	public:
	
		//================================================
		//==         MAIN PARAMETERS 
		//================================================
		/**
		* \brief Set blob id
		*/
		void SetId(int id){Id=id;}
		
		/**
		* \brief Set blob name
		*/
		//void SetName(std::string name){Name=name;}
		void SetName(std::string name){TNamed::SetName(name.c_str());}

		/**
		* \brief Set image range
		*/
		//void SetImageRange(Caesar::ImgRange range){m_ImgRange= range;}
		/**
		* \brief Set image range
		*/
		/*
		void SetImageRange(float xmin,float xmax,float ymin,float ymax){
			m_ImgRange= Caesar::ImgRange(xmin,xmax,ymin,ymax);
			//m_ImageMinX= xmin;
			//m_ImageMaxX= xmax;
			//m_ImageMinY= ymin;
			//m_ImageMaxY= ymax;
		}
		*/

		/**
		* \brief Get image range
		*/
		/*
		void GetImageRange(float& xmin,float& xmax,float& ymin,float& ymax){
			m_ImgRange.GetRange(xmin,xmax,ymin,ymax);
			//xmin= m_ImageMinX;
			//xmax= m_ImageMaxX;
			//ymin= m_ImageMinY;
			//ymax= m_ImageMaxY;
		}
		*/

		/**
		* \brief Get image range
		*/	
		//const ImgRange& GetImageRange() const {return m_ImgRange;}

		/**
		* \brief Is blob at image edge
		*/
		bool IsAtEdge(){return HasPixelsAtEdge;}	

		/**
		* \brief Set edge flag 
		*/
		void SetEdgeFlag(bool choice){HasPixelsAtEdge=choice;}

		

		//================================================
		//==         PIXELS 
		//================================================
		/**
		* \brief Get number of pixels
		*/
		int GetNPixels(){return static_cast<int>(m_Pixels.size());}

		/**
		* \brief Get pixel collection
		*/	
		const PixelCollection& GetPixels() const {return m_Pixels;}

		/**
		* \brief Set pixel collection
		*/
		void SetPixels(PixelCollection& pixels){m_Pixels= pixels;}		

		/**
		* \brief Add a pixel to collection
		*/
		int AddPixel(Pixel* pixel,bool makeCopy=false);

		/**
		* \brief Has pixels stored?
		*/
		bool HasPixels(){return (!m_Pixels.empty());}

		/**
		* \brief Get access to pixel by index
		*/
		Pixel* GetPixel(int index) {
			if(m_Pixels.empty() || index<0 || index>=(signed)(m_Pixels.size()) ) return nullptr;
			return m_Pixels[index];
		}

		/**
		* \brief Get access to seed pixel indexes
		*/
		std::vector<int> GetSeedPixelIndexes(){
			std::vector<int> seedPixelIndexes;
			if(m_Pixels.empty()) return seedPixelIndexes;
			std::vector<Pixel*>::iterator it = m_Pixels.begin();
			while ((it = std::find_if(it, m_Pixels.end(), MatchPixelType(Pixel::eSeed))) != m_Pixels.end()) {
				int index = std::distance(m_Pixels.begin(), it);
				seedPixelIndexes.push_back(index);
    		it++;
			}//end loop
			return seedPixelIndexes;
		}//close GetSeedPixelIndexes()

		
		//================================================
		//==         BLOB PARAMETERS 
		//================================================
		/**
		* \brief Compute stats
		*/
		int ComputeStats(bool computeRobustStats=true,bool forceRecomputing=false);

		/**
		* \brief Has stats?
		*/
		bool HasStats(){return m_HasStats;}

		/**
		* \brief Has stats?
		*/
		void SetHasStats(bool value){m_HasStats=value;}

		/**
		* \brief Compute morphology parameters
		*/
		int ComputeMorphologyParams();

		/**
		* \brief Compute Zernike moments
		*/
		int ComputeZernikeMoments(int order=6);

		/**
		* \brief Has morphology parameters computed
		*/
		bool HasParameters(){return m_HasParameters;}
		/**
		* \brief Set has morphology pars
		*/
		void SetHasParameters(bool value){m_HasParameters= value;}
	
		/**
		* \brief Get 1st pixel moment (= mean)
		*/
		double GetM1(){return m_M1;}
		/**
		* \brief Set 1st pixel moment (= mean)
		*/
		void SetM1(double value){m_M1= value;}
		/**
		* \brief Get 2nd pixel moment 
		*/
		double GetM2(){return m_M2;}
		/**
		* \brief Set 2nd pixel moment 
		*/
		void SetM2(double value){m_M2= value;}

		/**
		* \brief Get 3rd pixel moment 
		*/
		double GetM3(){return m_M3;}
		/**
		* \brief Set 3rd pixel moment 
		*/
		void SetM3(double value){m_M3= value;}
		/**
		* \brief Get 4th pixel moment 
		*/
		double GetM4(){return m_M4;}

		/**
		* \brief Set 4th pixel moment 
		*/
		void SetM4(double value){m_M4=value;}
		/**
		* \brief Get 1st pixel curvature moment (= mean)
		*/
		double GetM1Curv(){return m_M1_curv;}	
		/**
		* \brief Set 1st pixel curvature moment (= mean)
		*/
		void SetM1Curv(double value){m_M1_curv=value;}
		/**
		* \brief Get 2nd pixel curvature moment 
		*/
		double GetM2Curv(){return m_M2_curv;}	
		/**
		* \brief Set 2nd pixel curvature moment 
		*/
		void SetM2Curv(double value){m_M2_curv=value;}

		/**
		* \brief Get pixel flux sum
		*/
		double GetS(){return m_S;}
		/**
		* \brief Set pixel flux sum
		*/
		void SetS(double value){m_S= value;}
		/**
		* \brief Get pixel flux max
		*/
		double GetSmax(){return m_Smax;}
		/**
		* \brief Set pixel flux max
		*/
		void SetSmax(double value){m_Smax= value;}
		/**
		* \brief Get pixel flux min
		*/
		double GetSmin(){return m_Smin;}
		/**
		* \brief Set pixel flux min
		*/
		void SetSmin(double value){m_Smin= value;}
		/**
		* \brief Get pixel flux XX correlation
		*/
		double GetSxx(){return m_Sxx;}
		/**
		* \brief Set pixel flux XX correlation
		*/
		void SetSxx(double value){m_Sxx=value;}
		/**
		* \brief Get pixel flux YY correlation
		*/
		double GetSyy(){return m_Syy;}
		/**
		* \brief Get pixel flux YY correlation
		*/
		void SetSyy(double value){m_Syy=value;}
		/**
		* \brief Get pixel flux XY correlation
		*/
		double GetSxy(){return m_Sxy;}
		/**
		* \brief Set pixel flux XY correlation
		*/
		void SetSxy(double value){m_Sxy=value;}

		/**
		* \brief Get pixel flux sum weighted by position x
		*/
		double GetSx(){return m_Sx;}
		/**
		* \brief Set pixel flux sum weighted by position x
		*/
		void SetSx(double value){m_Sx= value;}
		/**
		* \brief Get pixel flux sum weighted by position y
		*/
		double GetSy(){return m_Sy;}
		/**
		* \brief Set pixel flux sum weighted by position y
		*/
		void SetSy(double value){m_Sy= value;}
		/**
		* \brief Get id of maximum flux pixel
		*/
		long int GetSmaxPixId() {return m_PixIdmax;}
		/**
		* \brief Set id of maximum flux pixel
		*/
		void SetSmaxPixId(long int value){m_PixIdmax= value;}
		/**
		* \brief Get id of minimum flux pixel
		*/
		long int GetSminPixId() {return m_PixIdmin;}
		/**
		* \brief Set id of minimum flux pixel
		*/
		void SetSminPixId(long int value){m_PixIdmin= value;}
		/**
		* \brief Get pixel sum of curvature
		*/
		double GetScurv(){return m_S_curv;}
		/**
		* \brief Set pixel sum of curvature
		*/
		void SetScurv(double value){m_S_curv= value;}
		/**
		* \brief Get pixel sum of edgeness
		*/
		double GetSedge(){return m_S_edge;}	
		/**
		* \brief Set pixel sum of edgeness
		*/	
		void SetSedge(double value){m_S_edge=value;}
		
		/**
		* \brief Get source x-y range
		*/
		void GetSourceRange(float& xmin,float& xmax,float& ymin,float& ymax){
			xmin= m_Xmin; 
			xmax= m_Xmax;
			ymin= m_Ymin;
			ymax= m_Ymax;
		}
		/**
		* \brief Set source x-y range
		*/
		void SetSourceRange(float xmin,float xmax,float ymin,float ymax){
			m_Xmin= xmin; 
			m_Xmax= xmax;
			m_Ymin= ymin;
			m_Ymax= ymax;
		}

		/**
		* \brief Get source pixel coordinate range
		*/
		void GetSourcePixelRange(long int& ixmin,long int& ixmax,long int& iymin,long int& iymax){
			ixmin= m_Ix_min; 
			ixmax= m_Ix_max;
			iymin= m_Iy_min;
			iymax= m_Iy_max;
		}

		/**
		* \brief Set source pixel coordinate range
		*/
		void SetSourcePixelRange(long int ixmin,long int ixmax,long int iymin,long int iymax){
			m_Ix_min= ixmin; 
			m_Ix_max= ixmax;
			m_Iy_min= iymin;
			m_Iy_max= iymax;
		}

		

		/**
		* \brief Dump blob info
		*/
		void Print(){
			cout<<"*** BLOB NO. "<<Id<<" ***"<<endl;
			cout<<"N= "<<NPix<<" Smin/Smax="<<m_Smin<<"/"<<m_Smax<<" Xmin/Xmax="<<m_Xmin<<"/"<<m_Xmax<<", Ymin/Ymax="<<m_Ymin<<"/"<<m_Ymax<<endl;
			cout<<"X0="<<X0<<" Y0="<<Y0<<" Mean="<<Mean<<" RMS="<<RMS<<" Median="<<Median<<" MedianRMS="<<MedianRMS<<endl;
			cout<<"****************************"<<endl;
		}
		/**
		* \brief Generate an image from source pixel
		*/
		Image* GetImage(ImgType mode,int pixMargin=1);

		//================================================
		//==         CONTOURS
		//================================================
		/**
		* \brief Return contours
		*/
		std::vector<Contour*> GetContours(){return m_Contours;}
		/**
		* \brief Return contour with index
		*/
		Contour* GetContour(int index) {
			if(index<0 || index>=(int)m_Contours.size() ) return 0;
			return m_Contours[index];
		}
		/**
		* \brief Add contour
		*/
		void AddContour(Contour* aContour){
			m_Contours.push_back(aContour);
		}

		

	private:

		/**
		* \brief Initialize data
		*/
		void Init();

		/**
		* \brief Update moments with given pixel
		*/
		void UpdateMoments(Pixel* pixel);
		/**
		* \brief Reset moments
		*/
		void ResetMoments();
		/**
		* \brief Reset stats
		*/
		void ResetStats();
		/**
		* \brief Clear pixels
		*/
		void ClearPixels();
		/**
		* \brief Clear contours
		*/
		void ClearContours();

	public:
		
		bool HasPixelsAtEdge;
		long int MinDistanceToEdgeX;
		long int MinDistanceToEdgeY;

		//Main params
		long int Id;//Blob id
		//std::string Name;//Blob name
			
		//Stats params
		long int NPix;//Number of pixels in blob
		double Mean;//mean = M1/N
		double RMS;
		double Skewness;
		double Median;
		double MedianRMS;
		double X0;//X position average
		double Y0;//Y position average
			
		//Curvature Moments
		double Mean_curv;
		double RMS_curv;
		double Median_curv;
		double MedianRMS_curv;

		//2D morphological pars
		std::vector<double> Moments;
		std::vector<double> HuMoments;
		std::vector<double> ZMMoments;	

	protected:


		bool m_HasStats;	
		bool m_HasParameters;	

		//Pixel intensity moments
		double m_M1;//1st moment
		double m_M2;//2nd moment
		double m_M3;//3rd moment
		double m_M4;//4th moment
				
		//Pixel curvature moments
		double m_M1_curv;//1st moment
		double m_M2_curv;//2nd moment

		//Moments accumulator
		double m_S;//sum of pixel signals
		double m_Smax;//max of pixel signals
		double m_Smin;//min of pixel signals
		double m_Sxx;
		double m_Syy;
		double m_Sxy;
		double m_Sx;//Signal-weighted X position average
		double m_Sy;//Signal-weighted Y position average
		long int m_PixIdmax;//id of pixel with max signal
		long int m_PixIdmin;//id of pixel with min signal
		
		double m_S_curv;//sum of pixel curvature
		double m_S_edge;//sum of edge estimator


		float m_Xmin;
		float m_Xmax;
		float m_Ymin;
		float m_Ymax;
		long int m_Ix_min;
		long int m_Ix_max;
		long int m_Iy_min;
		long int m_Iy_max;

		//Pixel collection
		PixelCollection m_Pixels;
			
		//Contour collection
		std::vector<Contour*> m_Contours;
	
	ClassDef(Blob,1)

};//close Blob()



#ifdef __MAKECINT__
#pragma link C++ class Blob+;
#pragma link C++ class vector<Blob>+;
#pragma link C++ class vector<Blob*>+;
#endif

}//close namespace

#endif

