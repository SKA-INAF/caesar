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

class ImgMetaData;
class Contour;
class WCS;

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
		void SetName(std::string name){TNamed::SetName(name.c_str());}

		/**
		* \brief Is blob at image edge
		*/
		bool IsAtEdge(){return HasPixelsAtEdge;}	

		/**
		* \brief Set edge flag 
		*/
		void SetEdgeFlag(bool choice){HasPixelsAtEdge=choice;}

		/**
		* \brief Set image metadata
		*/
		int SetImageMetaData(ImgMetaData* data,double crpix1Offset=0,double crpix2Offset=0)
		{
			//################################################################################
			//## With tile image the crpix keyword is shifted in image coordinates (ix,iy). 
			//## Source mostly uses original phys coordinates (x,y) referred to original image, 
			//## so we allow for shifting the crpix val, otherwise skycoordinates will be shifted
			//################################################################################
			if(!data) return -1;
			CodeUtils::DeletePtr<ImgMetaData>(m_imgMetaData);//delete existing
			m_imgMetaData= new ImgMetaData;
			*m_imgMetaData = *data;
			m_imgMetaData->Cx+= crpix1Offset;
			m_imgMetaData->Cy+= crpix2Offset;
			return 0;
		}

		/**
		* \brief get image metadata
		*/
		ImgMetaData* GetImageMetaData(){return m_imgMetaData;}

		/**
		* \brief Set image metadata
		*/
		bool HasImageMetaData(){
			if(!m_imgMetaData) return false; 
			return true;
		}

		//================================================
		//==         PIXELS 
		//================================================
		/**
		* \brief Get number of pixels
		*/
		long int GetNPixels(){return static_cast<long int>(m_Pixels.size());}

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

		/**
		* \brief Find and return largest-significance seed pixel 
		*/
		Pixel* GetSeedPixel()
		{
			//Loop over pixels and return the largest
			long int seedPixelIndex= -1;
			double Zmax= -1.e+99;
			for(size_t i=0;i<m_Pixels.size();i++){
				int pixType= m_Pixels[i]->type;
				if(pixType!=Pixel::eSeed) continue;
				std::pair<double,double> bkgInfo= m_Pixels[i]->GetBkg();
				double bkgLevel= bkgInfo.first;
				double noiseLevel= bkgInfo.second;
				double S= m_Pixels[i]->S;
				double Z= (S-bkgLevel)/noiseLevel;
				if(fabs(Z)>Zmax) {
					Zmax= Z;
					seedPixelIndex= i;
				}
			}//end loop pixels

			//Check if seed exists and was found
			if(seedPixelIndex<0) return nullptr;

			//Return seed pixel
			return m_Pixels[seedPixelIndex];

		}//close GetSeedPixel()

		/**
		* \brief Mark pixels with significance within given thresholds as halo pixels
		*/
		void MarkHaloPixels(double ZThr){
			for(size_t i=0;i<m_Pixels.size();i++){
				double bkg= m_Pixels[i]->bkgLevel;
				double rms= m_Pixels[i]->noiseLevel;
				double S= m_Pixels[i]->S;
				if(bkg!=0 && rms!=0 && fabs((S-bkg)/rms)<ZThr ) {
					m_Pixels[i]->type= Pixel::eHalo;
				}
			}//end pixel loop
		}//close MarkHaloPixels()
		
		//================================================
		//==         BLOB PARAMETERS 
		//================================================
		/**
		* \brief Compute stats
		*/
		int ComputeStats(bool computeRobustStats=true,bool forceRecomputing=false,bool useParallelMedian=false);

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
		double GetSmax() const {return m_Smax;}
		/**
		* \brief Set pixel flux max
		*/
		void SetSmax(double value){m_Smax= value;}
		/**
		* \brief Get pixel flux min
		*/
		double GetSmin() const {return m_Smin;}
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
		* \brief Get source pixel x_min
		*/
		long int GetXmin(){return m_Xmin;}

		/**
		* \brief Get source pixel x_max
		*/
		long int GetXmax(){return m_Xmax;}

		/**
		* \brief Get source pixel y_min
		*/
		long int GetYmin(){return m_Ymin;}

		/**
		* \brief Get source pixel y_max
		*/
		long int GetYmax(){return m_Ymax;}

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
		* \brief Get source pixel ix_min
		*/
		long int GetIxMin(){return m_Ix_min;}

		/**
		* \brief Get source pixel ix_max
		*/
		long int GetIxMax(){return m_Ix_max;}

		/**
		* \brief Get source pixel iy_min
		*/
		long int GetIyMin(){return m_Iy_min;}

		/**
		* \brief Get source pixel iy_max
		*/
		long int GetIyMax(){return m_Iy_max;}


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
		* \brief Get sample source standard deviations
		*/
		int GetSampleStdDev(double& sigmaX,double& sigmaY,double& covXY);

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
		Image* GetImage(ImgType mode,int pixMargin=1,bool includeHaloPixels=false);	

		/**
		* \brief Generate an image from source pixel
		*/
		TH2D* GetWCSHisto(ImgType mode,int pixMargin=1,int coordSyst=-1);

		//================================================
		//==         CONTOURS
		//================================================
		/**
		* \brief Has contours
		*/
		bool HasContours(){return !m_Contours.empty();}

		/**
		* \brief Return contours. NB: Do not delete pointers.
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
		* \brief Return contours converted in WCS (TO BE DEPRECATED)
		*/
		//std::vector<Contour*> GetWCSContours(WorldCoor* wcs=0,int coordSystem=-1,int pixOffset=0,bool computePars=false);
		/**
		* \brief Return contours converted in WCS 
		*/
		std::vector<Contour*> GetWCSContours(WCS* wcs=0,int coordSystem=-1,int pixOffset=0,bool computePars=false);
		
		/**
		* \brief Return contour with index and convert to WCS (TO BE DEPRECATED)
		*/
		//Contour* GetWCSContour(int index,WorldCoor* wcs=0,int coordSystem=-1,int pixOffset=0,bool computePars=false);
		/**
		* \brief Return contour with index and convert to WCS
		*/
		Contour* GetWCSContour(int index,WCS* wcs=0,int coordSystem=-1,int pixOffset=0,bool computePars=false);



		/**
		* \brief Add contour
		*/
		void AddContour(Contour* aContour){
			m_Contours.push_back(aContour);
		}

		/**
		* \brief Is point on contour?
		*/
		bool IsPointOnContour(double x,double y,double tol=1);

		//================================================
		//==         BKG INFO
		//================================================
		
		/**
		* \brief Get bkg sum
		*/
		double GetBkgSum(){return m_bkgSum;}
		/**
		* \brief Get bkg rms sum
		*/
		double GetBkgRMSSum(){return m_bkgRMSSum;}
		/**
		* \brief Set bkg sum (USED BY SERIALIZER)
		*/
		void SetBkgSum(double bkg){m_bkgSum=bkg;}
		/**
		* \brief Set bkg rms sum (USED BY SERIALIZER)
		*/
		void SetBkgRMSSum(double rms){m_bkgRMSSum=rms;}
		/**
		* \brief Check if bkg info in a box around the blob has ben set
		*/
		bool HasBoxBkgInfo(){return m_hasBoxBkgInfo;}
		/**
		* \brief Set has box bkg info flag (USED BY SERIALIZER)
		*/
		void SetHasBoxBkgInfo(bool flag){m_hasBoxBkgInfo=flag;}

		/**
		* \brief Get bkg estimate in a box around the blob
		*/
		double GetBoxBkg(){return m_boxBkg;}
		/**
		* \brief Get bkg rms estimate in a box around the blob
		*/
		double GetBoxBkgRMS(){return m_boxBkgRMS;}
		/**
		* \brief Set bkg and rms estimate in a box around the blob
		*/
		void SetBoxBkgInfo(double bkg,double rms)
		{
			m_boxBkg= bkg;
			m_boxBkgRMS= rms;
			m_hasBoxBkgInfo= true;			
		}

		//================================================
		//==         WCS
		//================================================
		/**
		* \brief Get WCS from stored metadata (TO BE DEPRECATED)
		*/
		/*
		WorldCoor* GetWCS(int coordSystem=-1){
			if(!m_imgMetaData) return nullptr;
			return m_imgMetaData->GetWorldCoord(coordSystem);
		}
		*/
		/**
		* \brief Get WCS from stored metadata
		*/
		WCS* GetWCS(int coordSystem=-1){
			if(!m_imgMetaData) return nullptr;
			return m_imgMetaData->GetWCS(coordSystem);
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

		//Bkg/noise sum	
		double m_bkgSum;
		double m_bkgRMSSum;
		bool m_hasBoxBkgInfo;
		double m_boxBkg;
		double m_boxBkgRMS;

		//Image range
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

		//Image metadata
		ImgMetaData* m_imgMetaData;
	
	ClassDef(Blob,3)

};//close Blob()





#ifdef __MAKECINT__
#pragma link C++ class Blob+;
#pragma link C++ class vector<Blob>+;
#pragma link C++ class vector<Blob*>+;
#endif

}//close namespace

#endif

