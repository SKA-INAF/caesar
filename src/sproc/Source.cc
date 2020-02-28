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
* @file Source.cc
* @class Source
* @brief Source class
*
* Class representing an image source
* @author S. Riggi
* @date 20/01/2015
*/


#include <Source.h>

#include <Blob.h>
#include <Image.h>
#include <Contour.h>
#include <GraphicsUtils.h>
#include <WCSUtils.h>
#include <Consts.h>
#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <SpectralIndexData.h>
#include <AstroObject.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TPad.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TExec.h>

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

ClassImp(Caesar::Source)

namespace Caesar {


Source::Source() 
	: Blob() 
{
	//Initialize 
	Init();

}//close costructor

Source::Source(std::string name)
	: Blob(name)
{

}//close parametric constructor


Source::Source(std::vector<Pixel*>const& pixels,std::string name)
	: Blob(pixels,name)
{

}//close parametric constructor


Source::~Source()
{
	
}//close destructor


Source::Source(const Source& source) 
//	: Blob() 
{
  // Contour copy constructor
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy constuctor called...");
	#endif
  Init();
  ((Source&)source).Copy(*this);
}

void Source::Copy(TObject &obj) const 
{
	//Copy mother blob 
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy blob...");
	#endif
	Blob::Copy((Source&)obj);

	// Copy this source to source obj	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy source variables...");
	#endif
  ((Source&)obj).Type = Type;
	((Source&)obj).Flag = Flag;	
	((Source&)obj).SimType= SimType;
	((Source&)obj).SimMaxScale= SimMaxScale;
	((Source&)obj).m_BeamFluxIntegral = m_BeamFluxIntegral;
	((Source&)obj).m_IsGoodSource = m_IsGoodSource;	
	((Source&)obj).m_DepthLevel = m_DepthLevel;
	((Source&)obj).m_HasNestedSources = m_HasNestedSources;

	((Source&)obj).m_S_true= m_S_true;
	((Source&)obj).m_X0_true= m_X0_true;
	((Source&)obj).m_Y0_true= m_Y0_true;
	((Source&)obj).m_HasTrueInfo= m_HasTrueInfo;

	((Source&)obj).m_HasFitInfo= m_HasFitInfo;
	((Source&)obj).m_fitPars= m_fitPars;
	((Source&)obj).m_fitStatus= m_fitStatus;

	//Spectral index data
	((Source&)obj).m_hasSpectralIndexData= m_hasSpectralIndexData;
	((Source&)obj).m_spectralIndexData= m_spectralIndexData;
	((Source&)obj).m_hasComponentSpectralIndexData= m_hasComponentSpectralIndexData;
	((Source&)obj).m_componentSpectralIndexData= m_componentSpectralIndexData;

	//Astro object data
	((Source&)obj).m_hasAstroObjectData= m_hasAstroObjectData;
	((Source&)obj).m_astroObjects= m_astroObjects;
	((Source&)obj).m_hasComponentAstroObjectData= m_hasComponentAstroObjectData;
	((Source&)obj).m_componentAstroObjects= m_componentAstroObjects;

	//Set object class ids
	((Source&)obj).ObjClassId= ObjClassId;
	((Source&)obj).ObjClassSubId= ObjClassSubId;
	((Source&)obj).ObjLocationId= ObjLocationId;
	((Source&)obj).ObjConfirmed= ObjConfirmed;
	
	//Set object component class ids
	((Source&)obj).componentObjLocationIds= componentObjLocationIds;
	((Source&)obj).componentObjClassIds= componentObjClassIds;
	((Source&)obj).componentObjClassSubIds= componentObjClassSubIds;
	((Source&)obj).componentObjConfirmed= componentObjConfirmed;
	

	//Delete first a previously existing vector
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Delete existing nested source list ...");
	#endif
	for(size_t i=0;i<(((Source&)obj).m_NestedSources).size();i++){
		if( (((Source&)obj).m_NestedSources)[i] ){
			delete (((Source&)obj).m_NestedSources)[i];
			(((Source&)obj).m_NestedSources)[i]= 0;
		}
	}
	(((Source&)obj).m_NestedSources).clear();

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Copy nested source list ...");
	#endif
	((Source&)obj).m_NestedSource= 0;
	for(unsigned int i=0;i<m_NestedSources.size();i++){
		((Source&)obj).m_NestedSource= new Source;
		*(((Source&)obj).m_NestedSource)= *(m_NestedSources[i]);
		(((Source&)obj).m_NestedSources).push_back( ((Source&)obj).m_NestedSource );
	}
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("End copy source ...");
	#endif
	
}//close Copy()

Source& Source::operator=(const Source& source) { 
	// Operator =
  if (this != &source)  ((Source&)source).Copy(*this);
  return *this;
}


void Source::Init(){

	//Init source flags
	//Type= eUnknown;
	Type= eUnknownType;
	Flag= eCandidate;
	SimType= eUnknownSimClass;
	SimMaxScale= 0;
	m_BeamFluxIntegral= 0;
	
	m_IsGoodSource= true;

	//Init nested source info
	m_DepthLevel= 0;
	m_HasNestedSources= false;
	m_NestedSource= 0;
	m_NestedSources.clear();

	//Init source true info
	m_HasTrueInfo= false;
	m_S_true= -1;
	m_X0_true= -1;
	m_Y0_true= -1;

	//Init fit info
	m_HasFitInfo= false;
	m_fitStatus= eFitUnknownStatus;

	//Init spectral index data
	m_hasSpectralIndexData= false;
	m_hasComponentSpectralIndexData= false;
	m_componentSpectralIndexData.clear();

	//Init astro objects
	m_hasAstroObjectData= false;
	m_astroObjects.clear();
	m_hasComponentAstroObjectData= false;
	m_componentAstroObjects.clear();

	//Init object class ids
	ObjLocationId= eUNKNOWN_OBJECT_LOCATION;
	ObjClassId= eUNKNOWN_OBJECT;
	ObjClassSubId= eUNKNOWN_OBJECT;
	ObjConfirmed= false;

	//Init component object class ids
	componentObjLocationIds.clear();
	componentObjClassIds.clear();
	componentObjClassSubIds.clear();
	componentObjConfirmed.clear();

}//close Init()


int Source::Draw(int pixMargin,ImgType imgType,bool drawImg,bool drawContours,bool drawNested,bool drawFitComponents,int lineColor,int lineStyle,bool useWCS,int coordSyst)
{
	//Draw image
	if(drawImg){
		//Get source image
		Image* simg= GetImage(imgType,pixMargin);
		if(!simg){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to get image from source!");
			#endif
			return -1;
		}

		//Get metadata info
		std::string units= "";
		TString freq= "";
		if(m_imgMetaData){
			units= m_imgMetaData->BUnit;
			freq= Form("%1.2f %s",m_imgMetaData->Freq,(m_imgMetaData->FreqUnit).c_str());
		}

		//Get histo from image
		TH2D* htemp= simg->GetHisto2D(std::string(this->GetName()));
		htemp->SetStats(0);
		htemp->GetXaxis()->SetTitle("X");
		htemp->GetYaxis()->SetTitle("Y");
		htemp->GetZaxis()->SetTitle(units.c_str());
		htemp->GetZaxis()->SetTitleSize(0.05);
		htemp->GetZaxis()->SetTitleOffset(0.9);
		htemp->Draw();
		gPad->Update();
		if(useWCS) htemp->Draw("COLAZ");
		else htemp->Draw("COLZ");

		//Append Caesar image to current pad
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Appending source "<<this->GetName()<<" image to current pad ...");
		#endif
		simg->AppendPad();

		//Draw main source info
		TPaveText* sourceMainInfoText = new TPaveText(0.4,0.15,0.8,0.3,"NDC");
		sourceMainInfoText->AddText(Form("Name: %s, Type: %d, Freq: %s",this->GetName(),Type,freq.Data()));
		sourceMainInfoText->AddText(Form("NPix: %ld, Xmin/Xmax: %1.2f/%1.2f, Ymin/Ymax: %1.2f/%1.2f",NPix,m_Xmin,m_Xmax,m_Ymin,m_Ymax));
		sourceMainInfoText->AddText(Form("C(%1.2f,%1.2f), Wc(%1.2f,%1.2f)",X0,Y0,m_Sx,m_Sy));
		sourceMainInfoText->AddText(Form("Stot(mJy): %1.1f, Smin/Smax(mJy): %1.2f/%1.2f",m_S*1.e+3,m_Smin*1.e+3,m_Smax*1.e+3));
		sourceMainInfoText->SetTextAlign(12);
		sourceMainInfoText->SetTextSize(0.02);
		sourceMainInfoText->SetTextFont(52);
		sourceMainInfoText->SetFillColor(0);
		sourceMainInfoText->SetBorderSize(1);	
		sourceMainInfoText->Draw("same");
		
		//Set WCS axis	
		if(useWCS){
			gPad->Update();
			TGaxis* xaxis_wcs= new TGaxis;
			TGaxis* yaxis_wcs= new TGaxis;
			bool useImageCoords= false;
			int status= GraphicsUtils::SetWCSAxis(simg,*xaxis_wcs,*yaxis_wcs,coordSyst,useImageCoords);
			if(status>=0){
				TExec* ex = new TExec("ex","GraphicsUtils::PadUpdater_PhysCoords()");
   			htemp->GetListOfFunctions()->Add(ex);
				xaxis_wcs->Draw("same");
				yaxis_wcs->Draw("same");
			}
			else{
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to set gAxis!");
				#endif
			}	
		}//close if useWCS
	}//close if draw image

	//Drawing contours?
	if(drawContours){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Drawing #"<<m_Contours.size()<<" contours present for source "<<this->GetName()<<"...");
		#endif
		for(size_t i=0;i<m_Contours.size();i++){		
			TGraph* thisContourGraph= m_Contours[i]->GetGraph();
			if(thisContourGraph) {
				thisContourGraph->SetMarkerSize(8);
				thisContourGraph->SetMarkerSize(0.3);
				thisContourGraph->SetMarkerColor(lineColor);
				thisContourGraph->SetLineColor(lineColor);
				thisContourGraph->SetLineStyle(lineStyle);
				thisContourGraph->SetLineWidth(2);
				thisContourGraph->Draw("Lsame");
			}//close if 
		}//end loop contours
	}//close if draw contours

	//Draw fitted components (if any)
	if(drawFitComponents && m_HasFitInfo){
		std::vector<TEllipse*> ellipses= m_fitPars.GetFittedEllipses();
		#ifdef LOGGING_ENABLED
			INFO_LOG("Drawing #"<<ellipses.size()<<" fitted components present for source "<<this->GetName()<<"...");
		#endif
		for(size_t i=0;i<ellipses.size();i++){
			ellipses[i]->SetLineColor(lineColor);
			ellipses[i]->SetLineStyle(kDotted);
			ellipses[i]->SetLineWidth(2);
			ellipses[i]->SetFillColor(0);
			ellipses[i]->SetFillStyle(0);
			ellipses[i]->Draw("lsame");
		}//end loop fit ellipses
	}//close if


	//Drawing nested sources recursively?
	if(drawNested){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Draw source nested ...");
		#endif
		bool drawNestedImage= false;
		for(size_t i=0;i<m_NestedSources.size();i++){
			if(m_NestedSources[i]) {
				m_NestedSources[i]->Draw(pixMargin,imgType,drawNestedImage,drawContours,drawNested,drawFitComponents,lineColor+2,lineStyle,useWCS,coordSyst);
			}
		}
	}//close if nested
	
	return 0;

}//close Draw()


void Source::Draw(bool drawBoundingBox,bool drawEllipse,bool drawNested,int lineColor,int lineStyle)
{
	//Drawing contours?
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("#"<<m_Contours.size()<<" contours present for source "<<Id<<"...");
	#endif
	for(size_t i=0;i<m_Contours.size();i++){		
		TGraph* thisContourGraph= m_Contours[i]->GetGraph();
		if(thisContourGraph) {
			thisContourGraph->SetMarkerSize(8);
			thisContourGraph->SetMarkerSize(0.3);
			thisContourGraph->SetMarkerColor(lineColor);
			thisContourGraph->SetLineColor(lineColor);
			thisContourGraph->SetLineStyle(lineStyle);
			thisContourGraph->SetLineWidth(2);
			thisContourGraph->Draw("Lsame");
		}//close if 
		
		//Get bounding box
		if(drawBoundingBox){
			TPolyLine* thisBoundingBox= m_Contours[i]->GetBoundingBoxLine();
			if(thisBoundingBox){
				thisBoundingBox->SetLineColor(lineColor);
				thisBoundingBox->SetLineStyle(kDashed);
				thisBoundingBox->Draw("lsame");
			}
		}//close if

		//Get fitted ellipse
		if(drawEllipse){
			TEllipse* thisFittedEllipse= m_Contours[i]->GetFittedEllipse();
			if(thisFittedEllipse){
				thisFittedEllipse->SetLineColor(lineColor);
				thisFittedEllipse->SetLineStyle(kDotted);
				thisFittedEllipse->SetLineWidth(2);
				thisFittedEllipse->SetFillColor(0);
				thisFittedEllipse->SetFillStyle(0);
				thisFittedEllipse->Draw("lsame");
			}
		}//close if draw ellipse
	}//end loop contours

	if(drawNested){
		for(unsigned int i=0;i<m_NestedSources.size();i++){
			if(m_NestedSources[i]) m_NestedSources[i]->Draw(drawBoundingBox,drawEllipse,drawNested,lineColor+2,lineStyle);
		}
	}//close if nested

}//close Draw()


const std::string Source::GetDS9Region(bool dumpNestedSourceInfo,bool convertToWCS,WCS* wcs,int coordSystem)
{
	//Check if has pixels
	//NB: DS9 crashes miserably when given a polygon region with one point 
	if(NPix<=1) return std::string("");

	//Convert contours to WCS?
	std::stringstream sstream;
	
	//Define source tags
	//- source name
	std::string regionText= this->GetName();

	//- source color
	std::string regionColor= this->GetDS9RegionColor();

	//- source type
	std::string sourceTypeStr= GetSourceTypeStr(this->Type);

	//- source flag
	std::string sourceFlagStr= GetSourceFlagStr(this->Flag);

	//std::vector<std::string> regionTags {this->GetDS9RegionTag()};
	std::vector<std::string> regionTags {sourceTypeStr,sourceFlagStr};


	std::string region= "";
	bool useImageCoords= true;
	if(convertToWCS) useImageCoords= false;
	
	if(convertToWCS){
		int pixOffset= 1;
		std::vector<Contour*> contours_wcs= GetWCSContours(wcs,coordSystem,pixOffset);		
		if(contours_wcs.empty()){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to convert contours in WCS, region will be empty!");
			#endif
		}
		else{
			//Loop over WCS contours
			for(size_t i=0; i<contours_wcs.size(); i++){ 
				region= AstroUtils::ContourToDS9Region(contours_wcs[i],regionText,regionColor,regionTags,useImageCoords);
			}

			//Delete contours at the end
			CodeUtils::DeletePtrCollection<Contour>(contours_wcs);
		}

	}//close if
	else{
		for(size_t i=0; i<m_Contours.size(); i++){ 
			region= AstroUtils::ContourToDS9Region(m_Contours[i],regionText,regionColor,regionTags,useImageCoords);
		}
	}//close else

	sstream<<region;

	//###### FILL NESTED SOURCE REGIONS ###########
	//Fill nested source regions
	if(dumpNestedSourceInfo && m_HasNestedSources){			
		sstream<<endl;
		for(size_t k=0;k<m_NestedSources.size();k++){

			std::string regionText_nested= m_NestedSources[k]->GetName();
			std::string regionColor_nested= m_NestedSources[k]->GetDS9RegionColor();
			std::string sourceTypeStr_nested= GetSourceTypeStr(m_NestedSources[k]->Type);	
			std::string sourceFlagStr_nested= GetSourceFlagStr(m_NestedSources[k]->Flag);
			//std::vector<std::string> regionTags_nested {m_NestedSources[k]->GetDS9RegionTag()};
			std::vector<std::string> regionTags_nested {sourceTypeStr_nested,sourceFlagStr_nested};

			std::string region_nested= "";
			if(convertToWCS){
				int pixOffset= 1;
				std::vector<Contour*> contours_wcs= m_NestedSources[k]->GetWCSContours(wcs,coordSystem,pixOffset);
				if(contours_wcs.empty()){
					#ifdef LOGGING_ENABLED
						WARN_LOG("Failed to convert contours in WCS, region will be empty!");
					#endif
				}
				else{
					//Loop over WCS contours
					for(size_t i=0;i<contours_wcs.size(); i++){ 
						region_nested= AstroUtils::ContourToDS9Region(contours_wcs[i],regionText_nested,regionColor_nested,regionTags_nested,useImageCoords);
					}

					//Delete contours at the end
					CodeUtils::DeletePtrCollection<Contour>(contours_wcs);
				}
			}//close if
			else{
				std::vector<Contour*> nestedContours= m_NestedSources[k]->m_Contours;
				for(size_t i=0;i<nestedContours.size(); i++){ 
					region_nested= AstroUtils::ContourToDS9Region(nestedContours[i],regionText_nested,regionColor_nested,regionTags_nested,useImageCoords);
				}
			}//close else

			sstream<<region_nested;
			if(k!=m_NestedSources.size()-1) sstream<<endl;

		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close GetDS9Region()


const std::string Source::GetDS9FittedEllipseRegion(bool useFWHM,bool dumpNestedSourceInfo,bool convertToWCS,WCS* wcs,int coordSystem,bool useWCSSimpleConversion)
{
	//Check WCS & metadata
	if(convertToWCS && !wcs && !m_imgMetaData){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Requested to convert ellipse to WCS but no wcs was provided and no metadata to compute it are available!");
		#endif
		return std::string("");
	}

	bool useImageCoords= true;
	int pixOffset= 0;
	if(convertToWCS) {
		useImageCoords= false;
		pixOffset= 1;
	}

	//Check if source has fit info
	std::stringstream sstream;
	
	if(m_HasFitInfo){
		//Get fit ellipses
		std::vector<TEllipse*> ellipses;
		
		if(GetFitEllipses(ellipses,useFWHM,convertToWCS,wcs,coordSystem,pixOffset,useWCSSimpleConversion)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to get WorldCoord system from metadata!");
			#endif
			return std::string("");
		}

		//Get fit info
		//- fit quality flag
		int fitQuality= m_fitPars.GetFitQuality();
		std::string fitQualityFlagStr= GetSourceFitQualityStr(fitQuality);
	
		//Loop over fit ellipses and convert to DS9 regions
		for(size_t i=0;i<ellipses.size();i++){
			if(!ellipses[i]) continue;

			//Get fit component flag
			int fitComponentFlag= -1;
			m_fitPars.GetComponentFlag(fitComponentFlag,i);
			std::string fitComponentFlagStr= GetSourceFlagStr(fitComponentFlag);

			//Get fit component type
			int fitComponentType= -1;
			m_fitPars.GetComponentType(fitComponentType,i);
			std::string fitComponentTypeStr= GetSourceTypeStr(fitComponentType);

			//Get encoded string region
			std::string regionText(Form("%s_fitcomp%d",this->GetName(),(int)(i+1)));
			std::string regionColor= "red";
			std::vector<std::string> regionTags {fitComponentTypeStr,"fit-component",fitComponentFlagStr,fitQualityFlagStr};
			std::string region= AstroUtils::EllipseToDS9Region(ellipses[i],regionText,regionColor,regionTags,useImageCoords);
			sstream<<region;

			if(i!=ellipses.size()-1) sstream<<endl;
		}//end loop ellipses

		//Delete ellipses
		CodeUtils::DeletePtrCollection<TEllipse>(ellipses);

	}//close if has fit info

	//Loop over nested components and get fit ellipse regions
	if(dumpNestedSourceInfo && m_HasNestedSources){			
		for(size_t k=0;k<m_NestedSources.size();k++){	
			std::string nestedRegionStr= m_NestedSources[k]->GetDS9FittedEllipseRegion(useFWHM,false,convertToWCS,wcs,coordSystem);
			if(nestedRegionStr!="") sstream<<nestedRegionStr<<endl;
		}//end loop nested sources
	}//close if

	return sstream.str();
	
}//close GetDS9FittedEllipseRegion()



const std::string Source::GetDS9EllipseRegion(bool dumpNestedSourceInfo)
{			
	//ellipse x y radius radius angle
	std::stringstream sstream;
	sstream<<"ellipse ";
	for(unsigned int i=0; i<m_Contours.size(); i++){ 
		if(!m_Contours[i]->HasEllipseFit) continue;
		double EllX= m_Contours[i]->EllipseCenter.X();
		double EllY= m_Contours[i]->EllipseCenter.Y();
		double EllMajAxis= m_Contours[i]->EllipseMajAxis;
		double EllMinAxis= m_Contours[i]->EllipseMinAxis;
		double EllRotAxis= m_Contours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
		//TVector2 BBoxCenter= m_Contours[i]->BoundingBoxCenter;
		//double BBoxMinAxis=  m_Contours[i]->BoundingBoxMin;	
		//double BBoxMajAxis= m_Contours[i]->BoundingBoxMaj;
		//double BBoxAngle= m_Contours[i]->BoundingBoxAngle;	
		sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
	}
	sstream<<"# text={S"<<Id<<"}";

	if(dumpNestedSourceInfo && m_HasNestedSources){			
		sstream<<endl;
		for(unsigned int k=0;k<m_NestedSources.size();k++){
			sstream<<"ellipse ";
			std::vector<Contour*> nestedContours= m_NestedSources[k]->m_Contours;
			for(unsigned int i=0; i<nestedContours.size(); i++){ 
				if(!nestedContours[i]->HasEllipseFit) continue;
				double EllX= nestedContours[i]->EllipseCenter.X();
				double EllY= nestedContours[i]->EllipseCenter.Y();
				double EllMajAxis= nestedContours[i]->EllipseMajAxis;
				double EllMinAxis= nestedContours[i]->EllipseMinAxis;
				double EllRotAxis= nestedContours[i]->EllipseRotAngle-90;//-90 comes from DS9 strange format!
				//TVector2 BBoxCenter= nestedContours[i]->BoundingBoxCenter;
				//double BBoxMinAxis=  nestedContours[i]->BoundingBoxMin;	
				//double BBoxMajAxis= nestedContours[i]->BoundingBoxMaj;
				//double BBoxAngle= nestedContours[i]->BoundingBoxAngle;	
				sstream<<EllX+1<<" "<<EllY+1<<" "<<(EllMajAxis/2)<<" "<<(EllMinAxis/2)<<" "<<(EllRotAxis)<<" ";
			}//end loop contours
			//sstream<<"# text={S"<<Id<<"_Nest"<<k<<"}";
			sstream<<"# text={"<<this->GetName()<<"_Nest"<<k<<"}";
			if(k!=m_NestedSources.size()-1) sstream<<endl;
		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close GetDS9EllipseRegion()

bool Source::IsAdjacentSource(Source* aSource)
{		
	//Check input sources
	if(!aSource){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input Source!");
		#endif
		return false;
	}		

	//Check if pixel collections are empty
	if(GetNPixels()<=0 || aSource->GetNPixels()<=0){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("This or given source have no pixels, return not adjacent.");
		#endif
		return false;
	}	

	//Check if bouding boxes are overlapping
	//NB: If not skip the adjacency check
	bool areBoundingBoxesOverlapping= CheckBoxOverlapping(aSource);
	if(!areBoundingBoxesOverlapping){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Sources not adjacent as their bounding boxes (S1(x["<<m_Xmin<<","<<m_Xmax<<"), y["<<m_Ymin<<","<<m_Ymax<<"]), S2(x["<<aSource->m_Xmin<<","<<aSource->m_Xmax<<"), y["<<aSource->m_Ymin<<","<<aSource->m_Ymax<<"]) do not overlap...");
		#endif
		return false;		
	}

	//Find if there are adjacent pixels
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Finding if this source (pos=("<<X0<<","<<Y0<<"), #"<<m_Pixels.size()<<" pix) is adjacent to source (pos("<<aSource->X0<<","<<aSource->Y0<<"), #"<<(aSource->m_Pixels).size()<<" pix)");
	#endif
	auto it = std::find_first_of(
		(aSource->m_Pixels).begin(), (aSource->m_Pixels).end(), 
		m_Pixels.begin(), m_Pixels.end(),
		PixelMatcher::AreAdjacent
	);
	
	bool isAdjacent= false;
	if (it == (aSource->m_Pixels).end() ) {
		isAdjacent= false;
  } 
	else {
		isAdjacent= true;
		#ifdef LOGGING_ENABLED
  		DEBUG_LOG("Sources ("<<this->Id<<","<<aSource->Id<<") are adjacent (found a match at " << std::distance((aSource->m_Pixels).begin(), it)<<")");
		#endif
  }

	return isAdjacent;

}//close IsAdjacentSource()

long int Source::GetNMatchingPixels(std::vector<Pixel*>& matching_pixels,Source* aSource,bool sorted)
{
	//Initialize empty list of matching pixels
	matching_pixels.clear();

	//Check input sources
	if(!aSource){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input Source!");
		#endif
		return 0;
	}	

	//Check if pixel collections are empty
	if(GetNPixels()<=0 || aSource->GetNPixels()<=0){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("This or given source have no pixels, return 0.");
		#endif
		return 0;
	}	
	
	//Check if bounding boxes are overlapping
	//NB: If not skip the adjacency check
	bool areBoundingBoxesOverlapping= CheckBoxOverlapping(aSource);
	if(!areBoundingBoxesOverlapping){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Sources not overlapping as their bounding boxes (S1(x["<<m_Xmin<<","<<m_Xmax<<"), y["<<m_Ymin<<","<<m_Ymax<<"]), S2(x["<<aSource->m_Xmin<<","<<aSource->m_Xmax<<"), y["<<aSource->m_Ymin<<","<<aSource->m_Ymax<<"]) do not overlap...");
		#endif
		return 0;		
	}

	//Find intersection in pixel collections
	//NB: Sort pixel collections is mandatory for set_intersection()
	if(!sorted){
		std::sort(m_Pixels.begin(), m_Pixels.end(),PixelMatcher());
		std::sort((aSource->m_Pixels).begin(), (aSource->m_Pixels).end(),PixelMatcher());
	}

	std::set_intersection (
		(aSource->m_Pixels).begin(), (aSource->m_Pixels).end(), 
		m_Pixels.begin(), m_Pixels.end(),
		std::back_inserter(matching_pixels),
		PixelMatcher()
	);
  	
	long int nMatchingPixels= static_cast<long int>(matching_pixels.size());

	return nMatchingPixels;

}//close GetNMatchingPixels()


bool Source::FindSourceMatchByOverlapArea(SourceOverlapMatchPars& pars, const std::vector<Source*>& sources, float matchOverlapThr)
{
	//Reset pars
	pars.ResetPars();

	//Check given threshold
	if(matchOverlapThr<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid threshold ("<<matchOverlapThr<<") given (hint: should be >0)!");
		#endif
		return false;
	}
 
	//Return if given source collection to be searched is empty
	if(sources.empty()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Given source collection is empty, no match can be searched!");
		#endif
		return false;
	}

	//Return if this source has no pixels
	if(m_Pixels.empty() || NPix<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("This source has no pixels stored, cannot search for a match with a source catalog!");
		#endif
		return false;
	}	

	//Loop over sources to find matches
	//NB: Two pixel fractions are compared (M=overlap pixels, N=pixels of true source, N_rec=pixels of rec source)
	// 1) f= M/N>thr (to ensure rec source is not a tiny fraction of true source)  
	// 2) f_rec= M/N_rec>thr (to ensure rec source is not encompassing by a large amount the true source) 
	std::vector<SourceOverlapMatchPars> tmpMatchPars;
	std::vector<long int> overlappingSourceIndexes;

	for(size_t j=0;j<sources.size();j++){
		long int NPix_rec= sources[j]->GetNPixels();
		double S_rec= sources[j]->GetS();	
		double Sx_rec= sources[j]->GetSx(); 
		double Sy_rec= sources[j]->GetSy();		

		std::vector<Pixel*> matching_pixels;
		long int NMatchingPixels= this->GetNMatchingPixels(matching_pixels,sources[j]);
		double f= (double)(NMatchingPixels)/(double)(NPix);
		double f_rec= (double)(NMatchingPixels)/(double)(NPix_rec);		

		//Skip if no match is found
		if(NMatchingPixels<=0) continue;

		//Store overlapping source index
		overlappingSourceIndexes.push_back(j);
			
		//Check if overlap is above required threshold.
		//If so compute and store overlap info
		if( f>=matchOverlapThr && f_rec>=matchOverlapThr ){
			
			//Compute flux sum of overlapping pixels
			double S_match= 0.;
			for (auto item : matching_pixels) S_match+= item->S;	

			//Compute flux ratio overlapping/tot
			double Sratio= S_match/m_S;	
			double Sratio_rec= S_match/S_rec;	
	
			//Compute difference in position
			double dX= Sx_rec - m_Sx;
			double dY= Sy_rec - m_Sy;

			//Add match source pars to tmp list
			SourceOverlapMatchPars overlapPars(j,f,f_rec);
			overlapPars.Sratio= Sratio;
			overlapPars.Sratio_rec= Sratio_rec;
			overlapPars.dX= dX;
			overlapPars.dY= dY;		
			tmpMatchPars.push_back(overlapPars);
	
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source "<<this->GetName()<<" (X0="<<X0<<", Y0="<<Y0<<", N="<<NPix<<"): found match with source "<<sources[j]->GetName()<<" (X0="<<sources[j]->X0<<", Y0="<<sources[j]->Y0<<", N="<<sources[j]->NPix<<"), NMatchingPixels="<<NMatchingPixels<<", f="<<f<<", f_rec="<<f_rec<<" (t="<<matchOverlapThr<<"), Sratio="<<Sratio<<", Sratio_rec="<<Sratio_rec<<", dX="<<dX<<", dY="<<dY);
			#endif
		}
		
	}//end loop sources	

	//Check if match is found
	if(tmpMatchPars.empty()) return false;

	//Search for best overlap in case of multiple matches
	//NB: Consider 2 rec match sources: 1st) has larger f, 2nd) has larger f_rec. Which is the best one?
	//    It is assumed here that the best is those that covers true source with a larger fraction, e.g. the 1st)
	if(tmpMatchPars.size()>1) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<tmpMatchPars.size()<<" source matches found, searching for the best one ...");
		#endif

		double overlap_best= -1.e+99;
		int best_index= 0;

		for(size_t j=0;j<tmpMatchPars.size();j++){
			double overlap= tmpMatchPars[j].overlapFraction;
			if(overlap>=overlap_best){
				best_index= j;
				overlap_best= overlap;
			}
		}//end loop matched sources

		pars= tmpMatchPars[best_index];

	}//close if
	else{
		pars= tmpMatchPars[0];
	}

	//Set overlapping source indexes
	pars.overlappingSourceIndexes= overlappingSourceIndexes;

	return true;

}//close FindSourceMatchByOverlapArea()


bool Source::FindSourceMatchByPos(std::vector<SourcePosMatchPars>& matchPars, const std::vector<Source*>& sources, float matchPosThr)
{
	//Reset pars
	matchPars.clear();
	//pars.ResetPars();

	//Check given threshold
	if(matchPosThr<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid threshold ("<<matchPosThr<<") given (hint: should be >0)!");
		#endif
		return false;
	}
 
	//Return if given source collection to be searched is empty
	if(sources.empty()) {
		#ifdef LOGGING_ENABLED
			WARN_LOG("Given source collection is empty, no match can be searched!");
		#endif
		return false;
	}

	//Return if this source has no pixels
	if(m_Pixels.empty() || NPix<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("This source has no pixels stored, cannot search for a match with a source catalog!");
		#endif
		return false;
	}	


	//Loop over all reconstructed sources and match them in position
	//std::vector<SourcePosMatchPars> tmpMatchPars;	
	
	for(size_t j=0;j<sources.size();j++){
	
		//Find match for mother source
		std::vector<SourcePosMatchPars> matchInfo;
		bool foundMatch= FindSourceMatchByPos(matchInfo,j,-1,sources[j],matchPosThr);
		if(foundMatch){
			matchPars.insert(matchPars.end(),matchInfo.begin(),matchInfo.end());		
		}

		//Find match for nested source (if any present)
		std::vector<Source*> nestedSources= sources[j]->GetNestedSources();
		for(size_t k=0;k<nestedSources.size();k++){
			std::vector<SourcePosMatchPars> matchInfo_nested;
			bool foundMatch_nested= FindSourceMatchByPos(matchInfo_nested,j,k,nestedSources[k],matchPosThr);
			if(foundMatch_nested){
				matchPars.insert(matchPars.end(),matchInfo_nested.begin(),matchInfo_nested.end());		
			}
		}//end loop nested sources

	}//end loop sources

	//Check if match is found
	if(matchPars.empty()) return false;

	//Search for best overlap in case of multiple matches
	if(matchPars.size()>1) {
		#ifdef LOGGING_ENABLED
			INFO_LOG("#"<<matchPars.size()<<" source matches found, sorting by position difference ...");
		#endif
		std::sort(matchPars.begin(),matchPars.end());		
	}

	return true;

}//close FindSourceMatchByPos()

bool Source::FindSourceMatchByPos(std::vector<SourcePosMatchPars>& matchPars, long int source_index, long int nested_source_index, Source* source, float matchPosThr)
{
	//Set reference position to be compared with collection
	double X0_ref= X0;
	double Y0_ref= Y0;
	//double X0_ref= m_Sx;
	//double Y0_ref= m_Sy;
	if(m_HasTrueInfo){
		X0_ref= m_X0_true;
		Y0_ref= m_Y0_true;
	}	

	//If source has fit info loop over fitted components to find best match
	bool hasFitInfo= source->HasFitInfo();
	bool foundMatch= false;
	if(hasFitInfo){
		SourceFitPars fitPars= source->GetFitPars();
		//double dist_best= 1.e+99;
		//long int best_component= -1;
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("# fit components="<<fitPars.GetNComponents());
		#endif
		for(int k=0;k<fitPars.GetNComponents();k++){
			double X0_fitcomp= fitPars.GetParValue(k,"x0");	
			double Y0_fitcomp= fitPars.GetParValue(k,"y0");
			double dx= fabs(X0_fitcomp-X0_ref);
			double dy= fabs(Y0_fitcomp-Y0_ref);
			double dist= sqrt(dx*dx+dy*dy);
			
			if(dx<=matchPosThr && dy<=matchPosThr){//match found
				foundMatch= true;
				SourcePosMatchPars matchInfo(source_index,dist,k,nested_source_index);
				matchPars.push_back(matchInfo);
				//best_component= k;
				//dist_best= dist;
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Match info for source "<<this->GetName()<<" (pos("<<X0_ref<<","<<Y0_ref<<") against source "<<source->GetName()<<" (pos("<<X0_fitcomp<<","<<Y0_fitcomp<<"), dx="<<dx<<", dy="<<dy<<" dist="<<dist<<", matchPosThr="<<matchPosThr);	
				#endif
			}
			
		}//end loop fitted components	
			
		/*
		if(best_component!=-1){
			foundMatch= true;
			pars.posDiff= dist_best;
			pars.fitComponentIndex= best_component;
		}
		#ifdef LOGGING_ENABLED	
			INFO_LOG("Match info for source "<<this->GetName()<<" (pos("<<X0_ref<<","<<Y0_ref<<") against source "<<source->GetName()<<" (pos("<<source->X0<<","<<source->Y0<<"), found? "<<foundMatch<<", hasFitInfo? "<<hasFitInfo<<", dist_best="<<dist_best<<" best_component="<<best_component<<", matchPosThr="<<matchPosThr);
		#endif
		*/

	}//close if has fit info
	else{
		//No fit info available so compute offset using source barycenter 
		double dx= fabs(source->X0 - X0_ref);
		double dy= fabs(source->Y0 - Y0_ref);
		//double dx= fabs(source->GetSx() - X0_ref);
		//double dy= fabs(source->GetSy() - Y0_ref);
		double dist= sqrt(dx*dx+dy*dy);
		
		if(dx<=matchPosThr && dy<=matchPosThr){//match found
			foundMatch= true;
			SourcePosMatchPars matchInfo(source_index,dist,-1,nested_source_index);
			matchPars.push_back(matchInfo);
			//pars.posDiff= dist;
			//pars.fitComponentIndex= -1;
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("Match info for source "<<this->GetName()<<" (pos("<<X0_ref<<","<<Y0_ref<<") against source "<<source->GetName()<<" (pos("<<source->X0<<","<<source->Y0<<"), found? "<<foundMatch<<", dx="<<dx<<", dy="<<dy<<" matchPosThr="<<matchPosThr); 
			#endif
		}

	}//close else !has fit info

	
	return foundMatch;

}//close FindSourceMatchByPos()



int Source::MergeSource(Source* aSource,bool copyPixels,bool checkIfAdjacent,bool computeStatPars,bool computeMorphPars,bool sumMatchingPixels)
{
	//Check input sources
	if(!aSource){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given input Source!");
		#endif
		return -1;
	}	

	//If adjacency check is enabled check if sources are mergeable (e.g. if there is at least
	//a pixel adjacent to each other)
	if(checkIfAdjacent) {
		bool areAdjacent= IsAdjacentSource(aSource);
		if(!areAdjacent){
			#ifdef LOGGING_ENABLED	
				WARN_LOG("Sources are not adjacent nor overlapping, no merging will be performed!");
			#endif
			return -1;
		}
	}

	//Find differences in pixel collections
	//Pixel found in both collections are not merged to this source
	//First sort pixel collections (MANDATORY FOR set_difference() routine)
	std::sort(m_Pixels.begin(), m_Pixels.end(),PixelMatcher());
	std::sort((aSource->m_Pixels).begin(), (aSource->m_Pixels).end(),PixelMatcher());
	
	std::vector<Pixel*> pixelsToBeMerged;
	std::set_difference (
		(aSource->m_Pixels).begin(), (aSource->m_Pixels).end(), 
		m_Pixels.begin(), m_Pixels.end(),
		std::back_inserter(pixelsToBeMerged),
		PixelMatcher()
	);

	//Find and sum common pixels (if option enabled)
	//NB: Do this before adding pixels
	if(sumMatchingPixels){	
		//Find common pixels
		typedef std::vector< std::pair<long int,long int> > IndexPairs;
		IndexPairs intersect_indexes= CodeUtils::FindIntersectionIndexes ( 
			m_Pixels.begin(), m_Pixels.end(),
			(aSource->m_Pixels).begin(), (aSource->m_Pixels).end(), 
			PixelMatcher(),
			true
		);
		
		//Sum common pixels fluxes
		for(size_t i=0;i<intersect_indexes.size();i++){
			long int index1= intersect_indexes[i].first;
			long int index2= intersect_indexes[i].second;
			m_Pixels[index1]->AddPixelFlux( (aSource->m_Pixels)[index2] );
		}

		//Sum true fluxes
		if(m_HasTrueInfo && aSource->HasTrueInfo()){
			m_S_true+= aSource->GetTrueFlux();
		}

		
	}//close if sumMatchingPixels
  	
	//If no pixels are to be merged (e.g. all pixels overlapping) return?
	int nMergedPixels= (int)pixelsToBeMerged.size();
	if(nMergedPixels<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No pixels to be merged (perfectly overlapping sources?)!");
		#endif
		//return -1;
	}
	else{
		#ifdef LOGGING_ENABLED
			INFO_LOG("# "<<nMergedPixels<<"/"<<aSource->GetNPixels()<<" pixels to be merged to this source (Npix="<<this->GetNPixels()<<")...");
		#endif
	}

	//Now merge the pixels (if any to be merged)
	//Source moments will be updated in AddPixel()
	for(int i=0;i<nMergedPixels;i++){
		if(this->AddPixel(pixelsToBeMerged[i],copyPixels)<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to add pixel no. "<<i<<" to this source, skip to next pixel in list...");
			#endif
			continue;			
		}
	}//end loop pixels to be merged

	
	
	//Set new source type
	int mergedSourceType= aSource->Type;
	if( this->Type==eCompact || this->Type==ePointLike ){
		if(mergedSourceType==eCompact || mergedSourceType==ePointLike) this->Type= eCompact; 
		else if(mergedSourceType==eExtended || mergedSourceType==eCompactPlusExtended) this->Type= eCompactPlusExtended; 
		//else this->Type= eUnknown;
		else this->Type= eUnknownType; 
	}
	else if( this->Type==eExtended ){
		if(mergedSourceType==eCompact || mergedSourceType==ePointLike) this->Type= eCompactPlusExtended; 
		else if(mergedSourceType==eExtended) this->Type= eExtended;  
		else if(mergedSourceType==eCompactPlusExtended) this->Type= eCompactPlusExtended; 
		//else this->Type= eUnknown;
		else this->Type= eUnknownType; 
	}
	else if( this->Type==eCompactPlusExtended ){
		//if(mergedSourceType==eUnknown) this->Type= eUnknown; 
		if(mergedSourceType==eUnknownType) this->Type= eUnknownType; 
		else this->Type= eCompactPlusExtended;  
	}
	else{
		//this->Type= eUnknown;
		this->Type= eUnknownType; 
	}

	//Set sim max scale to max of the two sources
	SimMaxScale= max(SimMaxScale,aSource->SimMaxScale);

	//Set sim type if different to combination of both types
	if(SimType!=aSource->SimType){
		//Search if merged source sim type is already present
		//NB: This is to prevent adding the same types (e.g. when multiple sources with same type are added)
		std::string simtype1_str= std::to_string(SimType);
		std::string simtype2_str= std::to_string(aSource->SimType);
		std::size_t found = simtype1_str.find(simtype2_str);
		if (found==std::string::npos){//not found, add it
			//Add sim type string and sort
			std::string simtype12_str= simtype1_str + simtype2_str;
			std::sort(simtype12_str.begin(), simtype12_str.end());
			
			//Set new sim type (catch for errors)
			int simtype_merged= SimType;
			try{
				simtype_merged= std::stoi(simtype12_str);
				#ifdef LOGGING_ENABLED
					INFO_LOG("Changing simtype from (simtype1="<<SimType<<", simtype2="<<aSource->SimType<<") to simtype="<<simtype_merged);
				#endif
				SimType= simtype_merged;
			}
			catch(...){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("C++ exception occurred while converting merged stringified simtype "<<simtype12_str<<" to int code (will not add merged soure to simtype!");
				#endif
			}
		}//close if found
	}//close if sim type

	//Invalidate stats, pars if some pixels have been merged
	if(nMergedPixels>0){
		//At this stage stats (mean/median/etc...) are invalid and need to be recomputed if desired
		this->SetHasStats(false);//set stats to false to remember that current stats are not valid anymore and need to be recomputed
		if(computeStatPars){
			bool computeRobustStats= true;
			bool forceRecomputing= false;//no need to re-compute moments (already updated in AddPixel())
			if(sumMatchingPixels) forceRecomputing= true;
			this->ComputeStats(computeRobustStats,forceRecomputing);
		}

		//Contour and other parameters are also invalid
		this->SetHasParameters(false);
		if(computeMorphPars){
			this->ComputeMorphologyParams();
		}

		//At this stage fitting information (if present) is invalid (there are new pixels that have not been fitted)
		this->m_HasFitInfo= false;
		this->m_fitStatus= eFitUnknownStatus;
		this->m_fitPars.Reset();

	}//close if

	return 0;
	
}//close MergeSource()

bool Source::CheckBoxOverlapping(Source* aSource)
{
	//Check input blob
	if(!aSource){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Null ptr to input source given!");
		#endif
		return false;
	}
	
	//Get bounding box pars
	float xmin= aSource->m_Xmin;
	float xmax= aSource->m_Xmax;
	float ymin= aSource->m_Ymin;
	float ymax= aSource->m_Ymax;

	return HasBoxOverlap(xmin,xmax,ymin,ymax);

}//close CheckBoxOverlapping()


float Source::GetCentroidDistance(Source* aSource)
{
	//Check given source
	if(!aSource){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input source given, returning inf!");
		#endif
		return std::numeric_limits<float>::infinity();
	}

	//If one or both sources has no pars computed return inf
	if(!this->HasStats() || aSource->HasStats()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("One or both sources has no stats computed (you must compute pars before getting distance, returning inf)");
		#endif
		return std::numeric_limits<float>::infinity();
	}

	//Compute distance
	float dist2_x= pow(this->X0 - aSource->X0,2);
	float dist2_y= pow(this->Y0 - aSource->Y0,2);
	float dist= sqrt(dist2_x + dist2_y); 

	return dist;

}//close GetCentroidDistance()


int Source::Fit(SourceFitOptions& fitOptions)
{
	//Create source fitter
	SourceFitter fitter;
	if(fitter.FitSource(this,fitOptions)<0){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to fit source "<<this->GetName()<<" ...");
		#endif
		m_HasFitInfo= false;
		m_fitStatus= fitter.GetFitStatus();
		return -1;
	}
	
	//Get fit results
	SourceFitPars fitPars= fitter.GetFitPars();
	fitPars.Print();

	int fitStatus= fitter.GetFitStatus();
	if(fitStatus==eFitConverged || fitStatus==eFitConvergedWithWarns){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Fit of source "<<this->GetName()<<" converged (status="<<fitStatus<<"), storing fit parameters...");	
		#endif
		m_fitPars= fitPars;
		m_HasFitInfo= true;
	}
	
	//Store latest fit status
	m_fitStatus= fitStatus;

	return 0;

}//close Fit()


int Source::Fit(SourceFitOptions& fitOptions,SourceFitPars& initfitPars)
{
	//Get start fit pars vector
	std::vector<std::vector<double>> fitPars_start;
	initfitPars.GetFitParVec(fitPars_start);

	//Create source fitter
	SourceFitter fitter;
	if(fitter.FitSource(this,fitOptions,fitPars_start)<0){	
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to fit source "<<this->GetName()<<" ...");
		#endif
		m_HasFitInfo= false;
		m_fitStatus= fitter.GetFitStatus();
		return -1;
	}
	
	//Get fit results
	SourceFitPars fitPars= fitter.GetFitPars();
	fitPars.Print();

	int fitStatus= fitter.GetFitStatus();
	if(fitStatus==eFitConverged || fitStatus==eFitConvergedWithWarns){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Fit of source "<<this->GetName()<<" converged (status="<<fitStatus<<"), storing fit parameters...");	
		#endif
		m_fitPars= fitPars;
		m_HasFitInfo= true;
	}
	
	//Store latest fit status
	m_fitStatus= fitStatus;

	return 0;

}//close Fit()


int Source::GetFitEllipses(std::vector<TEllipse*>& fitEllipses,bool useFWHM,bool convertToWCS,WCS* wcs,int coordSystem,int pixOffset,bool useWCSSimpleConversion)
{
	//Init data 
	fitEllipses.clear();

	//Return empty vector if no fit info are available
	if(!m_HasFitInfo){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source "<<this->GetName()<<" has no fit info available, returning empty vector...");
		#endif
		return 0;
	}

	//Check if metadata are not available and convertToWCS is requested
	if(!m_imgMetaData && convertToWCS){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Requested to convert to WCS but not metadata have been set on this source to build the WCS!");
		#endif
		return -1;
	}

	//Get ellipses from fit parameters in pixel coordinates
	fitEllipses= m_fitPars.GetFittedEllipses(useFWHM);
	if(fitEllipses.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No fitted ellipse returned (there should be at least one component fitted, check!)...");
		#endif
		return -1;
	}

	//Convert to sky coordinates?
	if(convertToWCS){
		//Build the WCS with metadata if not given
		bool deleteWCS= false;
		if(!wcs){
			wcs= m_imgMetaData->GetWCS(coordSystem);
			if(!wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to get WorldCoord system from metadata!");
				#endif
				CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
				return -1;
			}
			deleteWCS= true;
		}
		
		//Convert ellipses to WCS
		std::vector<TEllipse*> fitEllipses_wcs;
		for(size_t i=0;i<fitEllipses.size();i++){
			TEllipse* fitEllipse_wcs= 0;
			if(useWCSSimpleConversion) fitEllipse_wcs= AstroUtils::PixelToWCSEllipseSimple(fitEllipses[i],wcs,pixOffset);
			else fitEllipse_wcs= AstroUtils::PixelToWCSEllipse(fitEllipses[i],wcs,pixOffset);

			if(!fitEllipse_wcs){
				#ifdef LOGGING_ENABLED
					ERROR_LOG("Failed to convert fit ellipse no. "<<i+1<<" to WCS!");
				#endif
				CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
				if(deleteWCS) WCSUtils::DeleteWCS(&wcs);
				return -1;
			}
			fitEllipses_wcs.push_back(fitEllipse_wcs);
		}//end loop ellipses

		//Delete fitEllipses and replace with wcs ones
		CodeUtils::DeletePtrCollection<TEllipse>(fitEllipses);
		fitEllipses.insert(fitEllipses.end(),fitEllipses_wcs.begin(),fitEllipses_wcs.end());
	
		//Delete wcs (if allocated)
		if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	}//close if convert to WCS

	
	return 0;

}//close GetFitEllipses()
	


int Source::FindComponentPeaks(std::vector<ImgPeak>& peaks,double peakZThr,int maxPeaks,int peakShiftTolerance,std::vector<int> kernels,int peakKernelMultiplicityThr,bool invertSearch)
{
	//Init
	peaks.clear();

	//Get source image
	int pixMargin= 0;
	Image* peakSearchMap= this->GetImage(eFluxMap,pixMargin);
	if(!peakSearchMap){
		#ifdef LOGGING_ENABLED	
			ERROR_LOG("Failed to get source image!");	
		#endif
		return -1;
	}

	//Invert search (e.g. to search valleys rather than peaks)
	if(invertSearch){
		peakSearchMap->Scale(-1);
	}

	//Get source pars (median, average bkg & noise)
	double Smedian= this->Median;
	double Smad= this->MedianRMS;
	double bkgMean= m_bkgSum/(double)(m_Pixels.size());
	double rmsMean= m_bkgRMSSum/(double)(m_Pixels.size());

	//Finding peaks
	#ifdef LOGGING_ENABLED
		INFO_LOG("Finding peaks in source (id="<<Id<<", name="<<this->GetName()<<")");	
	#endif
	//std::vector<TVector2> peakPoints;
	std::vector<ImgPeak> peakPoints;
	bool skipBorders= true;
	if(peakSearchMap->FindPeaks(peakPoints,kernels,peakShiftTolerance,skipBorders,peakKernelMultiplicityThr)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to find peaks in source (id="<<Id<<", name="<<this->GetName()<<") image!");
		#endif
		CodeUtils::DeletePtr<Image>(peakSearchMap);
		return -1;
	}

	//Select peaks (skip peaks at boundary or faint peaks)
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<peakPoints.size()<<" peaks found in source (id="<<Id<<", name="<<this->GetName()<<"), apply selection...");	
	#endif
	//std::vector<TVector2> peakPoints_selected;
	std::vector<ImgPeak> peakPoints_selected;
	std::vector<double> peakFluxes_selected;
	for(size_t i=0;i<peakPoints.size();i++){
		//double x= peakPoints[i].X();
		//double y= peakPoints[i].Y();
		double x= peakPoints[i].x;
		double y= peakPoints[i].y;
		long int gbin= peakSearchMap->FindBin(x,y);
		if(gbin<0){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to find gbin of peak("<<x<<","<<y<<"), this should not occur!");
			#endif
			CodeUtils::DeletePtr<Image>(peakSearchMap);
			return -1;
		}
		double Speak= peakSearchMap->GetBinContent(gbin);

		//Remove faint peaks in case more than one peak is found
		if(peakPoints.size()>1){			
			double Zpeak_imgbkg= 0;	
			double Zpeak_sourcebkg= 0;
			if(rmsMean!=0) Zpeak_imgbkg= (Speak-bkgMean)/rmsMean;
			if(Smad!=0) Zpeak_sourcebkg= (Speak-Smedian)/Smad;
			if(Zpeak_imgbkg<peakZThr) {
			//if(Zpeak_sourcebkg<peakZThr) {	
				#ifdef LOGGING_ENABLED
					INFO_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak_imgbkg="<<Zpeak_imgbkg<<", Zpeak_sourcebkg="<<Zpeak_sourcebkg<<"<"<<peakZThr<<")");
				#endif
				continue;
			}
		}//close if

		//Remove peaks lying on the source contour
		if(this->HasContours() && this->IsPointOnContour(x,y,0.5)) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Removing peak ("<<x<<","<<y<<") from the list as lying on source contour");
			#endif
			continue;
		}

		//Add peak to selected peak
		peakPoints_selected.push_back(peakPoints[i]);
		peakFluxes_selected.push_back(Speak);
	}//end loop peaks

	//Sort peaks by flux
	std::vector<size_t> sort_index;//sorting index
	std::vector<double> peakFluxes_sorted;
	CodeUtils::sort_descending(peakFluxes_selected,peakFluxes_sorted,sort_index);

	//Select peaks if more than max allowed
	int nPeaks_selected= static_cast<int>(peakPoints_selected.size());
	if(nPeaks_selected<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No components left in source (id="<<Id<<", name="<<this->GetName()<<") after selection!");	
		#endif
		CodeUtils::DeletePtr<Image>(peakSearchMap);
		return 0;
	}
	
	int nComponents= nPeaks_selected;
	if(maxPeaks>0) nComponents= std::min(nPeaks_selected,maxPeaks);
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<nComponents<<" components found in source (id="<<Id<<", name="<<this->GetName()<<"), max peaks="<<maxPeaks<<") ...");
	#endif

	peaks.clear();
	for(int i=0;i<nComponents;i++){
		size_t index= sort_index[i];
		peaks.push_back(peakPoints_selected[index]);
	}//end loop peaks


	//Delete peak map
	CodeUtils::DeletePtr<Image>(peakSearchMap);

	return 0;

}//close FindComponentPeaks()


int Source::FindBlendedComponents(std::vector<Source*>& deblendedComponents,std::vector<ImgPeak>& deblendedPeaks,double peakZThr,int maxPeaks,double sigmaMin,double sigmaMax,double sigmaStep,int minBlobSize,double thrFactor,int kernelFactor,int pixMargin)
{
	//Init
	deblendedComponents.clear();
	deblendedPeaks.clear();
	
	//Get source image
	//int pixMargin= 10;//NB: If set to 0 some peaks could be lost
	Image* sourceImg= this->GetImage(eFluxMap,pixMargin);
	if(!sourceImg){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get source image!");	
		#endif
		return -1;
	}

	//Get source pars (median, average bkg & noise)
	double Smedian= this->Median;
	double Smad= this->MedianRMS;
	double bkgMean= m_bkgSum/(double)(m_Pixels.size());
	double rmsMean= m_bkgRMSSum/(double)(m_Pixels.size());

	//Finding blended components
	std::vector<Source*> sources;
	std::vector<ImgPeak> peakPoints;
	int status= sourceImg->FindBlendedSources(
		sources,peakPoints,
		sigmaMin,sigmaMax,sigmaStep,
		minBlobSize,
		thrFactor,kernelFactor
	);
	if(status<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to find blended components in source "<<this->GetName()<<"!");
		#endif
		return -1;
	}

	if(peakPoints.empty()){
		#ifdef LOGGING_ENABLED
			INFO_LOG("No components found in source (name="<<this->GetName()<<")");	
		#endif
		CodeUtils::DeletePtr<Image>(sourceImg);
		CodeUtils::DeletePtrCollection<Source>(sources);
		return 0;
	}

	//Select components by peak flux (skip peaks at boundary or faint peaks)
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<peakPoints.size()<<" peaks found in source (name="<<this->GetName()<<"), apply selection...");	
	#endif
	std::vector<ImgPeak> peakPoints_selected;
	std::vector<double> peakFluxes_selected;
	std::vector<Source*> sources_selected;
	for(size_t i=0;i<peakPoints.size();i++){
		double x= peakPoints[i].x;
		double y= peakPoints[i].y;
		double Speak= peakPoints[i].S;

		//Remove faint peaks in case more than one peak is found
		if(peakPoints.size()>1){			
			double Zpeak_imgbkg= 0;	
			double Zpeak_sourcebkg= 0;
			if(rmsMean!=0) Zpeak_imgbkg= (Speak-bkgMean)/rmsMean;
			if(Smad!=0) Zpeak_sourcebkg= (Speak-Smedian)/Smad;
			if(Zpeak_imgbkg<peakZThr) {
			//if(Zpeak_sourcebkg<peakZThr) {
				#ifdef LOGGING_ENABLED
					INFO_LOG("Removing peak ("<<x<<","<<y<<") from the list as below peak significance thr (Zpeak_imgbkg="<<Zpeak_imgbkg<<", Zpeak_sourcebkg="<<Zpeak_sourcebkg<<"<"<<peakZThr<<")");
				#endif
				continue;
			}
		}//close if

		//Remove peaks lying on the source contour
		if(this->HasContours() && this->IsPointOnContour(x,y,0.5)) {
			#ifdef LOGGING_ENABLED
				INFO_LOG("Removing peak ("<<x<<","<<y<<") from the list as lying on source contour...");
			#endif
			continue;
		}

		//Add peak to selected peak
		peakPoints_selected.push_back(peakPoints[i]);
		peakFluxes_selected.push_back(Speak);
		sources_selected.push_back(sources[i]);
	}//end loop peaks

	//Sort peaks by flux
	std::vector<size_t> sort_index;//sorting index
	std::vector<double> peakFluxes_sorted;
	CodeUtils::sort_descending(peakFluxes_selected,peakFluxes_sorted,sort_index);

	//Select peaks if more than max allowed
	int nPeaks_selected= static_cast<int>(peakPoints_selected.size());
	if(nPeaks_selected<=0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No components left in source (id="<<Id<<", name="<<this->GetName()<<") after selection!");	
		#endif
		CodeUtils::DeletePtr<Image>(sourceImg);
		CodeUtils::DeletePtrCollection<Source>(sources);
		return 0;
	}
	
	int nComponents= nPeaks_selected;
	if(maxPeaks>0) nComponents= std::min(nPeaks_selected,maxPeaks);
	#ifdef LOGGING_ENABLED
		INFO_LOG("#"<<nComponents<<" components found in source (id="<<Id<<", name="<<this->GetName()<<"), max peaks="<<maxPeaks<<") ...");
	#endif

	deblendedPeaks.clear();
	deblendedComponents.clear();
	Source* aBlendedSource= 0;
	for(int i=0;i<nComponents;i++){
		size_t index= sort_index[i];
		deblendedPeaks.push_back(peakPoints_selected[index]);
	
		aBlendedSource= new Source;
		*aBlendedSource= *(sources_selected[index]);
		deblendedComponents.push_back(aBlendedSource);
	}//end loop peaks


	//Delete peak map
	CodeUtils::DeletePtr<Image>(sourceImg);
	CodeUtils::DeletePtrCollection<Source>(sources);

	return 0;

}//close FindBlendedComponents()


//std::string Source::GetIAUName(bool useWeightedPos,WorldCoor* wcs,int coordSystem)
std::string Source::GetIAUName(bool useWeightedPos,WCS* wcs,int coordSystem)
{
	//Init name
	std::string iau= "";

	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(!m_imgMetaData){
			#ifdef LOGGING_ENABLED
				WARN_LOG("No metadata are available to retrieve WCS!");
			#endif
			return iau;
		}
		//wcs= m_imgMetaData->GetWorldCoord(coordSystem);	
		wcs= m_imgMetaData->GetWCS(coordSystem);
		if(!wcs){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to get WCS from metadata!");
			#endif
			return iau;
		}
		deleteWCS= true;
	}//close if

	//Compute WCS string coords	
	std::string wcspos_str= "";
	int status= 0;
	if(useWeightedPos) status= AstroUtils::PixelToWCSStrCoords(wcspos_str,wcs,m_Sx,m_Sy);
	else status= AstroUtils::PixelToWCSStrCoords(wcspos_str,wcs,X0,Y0);

	if(status<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute WCS source pos in string format!");
		#endif

		//if(deleteWCS) CodeUtils::DeletePtr<WCS>(wcs);
		if(deleteWCS) WCSUtils::DeleteWCS(&wcs);
		return iau;
	}

	//Compute IAU name
	if(AstroUtils::GetIAUCoords(iau,wcspos_str)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to compute IAU name from WCS string coords!");
		#endif

		iau= "";
		//if(deleteWCS) CodeUtils::DeletePtr<WCS>(wcs);
		if(deleteWCS) WCSUtils::DeleteWCS(&wcs);
		return iau;
	}

	//Delete WCS
	//if(deleteWCS) CodeUtils::DeletePtr<WCS>(wcs);
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	return iau;

}//close GetIAUName()

int Source::GetWCSCoords(double& xwcs,double& ywcs,double x,double y,WCS* wcs,int coordSystem)
{
	//Init pos
	xwcs= -999;
	ywcs= -999;

	//If wcs is not given, retrieve it from metadata
	bool deleteWCS= false;
	if(!wcs){
		if(!m_imgMetaData){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Requested to get WCS centroid coords but no wcs was provided and no metadata are available to retrieve it!");
			#endif
			return -1;
		}
		
		wcs= m_imgMetaData->GetWCS(coordSystem);
		if(!wcs){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to get WorldCoord system from metadata!");
			#endif
			return -1;
		}
		deleteWCS= true;
	}//close if

	//Get WCS centroid coords
	int status= AstroUtils::PixelToWCSCoords(xwcs,ywcs,wcs,x,y);
	if(status<0){
		#ifdef LOGGING_ENABLED	
			WARN_LOG("Failed to get WCS coordinate corresponding to ("<<x<<","<<y<<")!");
		#endif
		return -1;
	}

	//Delete WCS
	//if(deleteWCS) CodeUtils::DeletePtr<WCS>(wcs);
	if(deleteWCS) WCSUtils::DeleteWCS(&wcs);

	return 0;

}//close GetWCSCoords()

int Source::GetSpectralAxisInfo(double& val,double& dval,std::string& units)
{
	//Init values
	val= -999;
	dval= -999;
	units= "";

	//Check if metadata are available
	if(!m_imgMetaData){
		#ifdef LOGGING_ENABLED
			WARN_LOG("No metadata available to get spectral axis info!");
		#endif
		return -1;
	}

	//Get info
	val= m_imgMetaData->Freq;	
	dval= m_imgMetaData->dFreq;
	units= m_imgMetaData->FreqUnit;

	return 0;

}//close GetSpectralAxisInfo()


int Source::ComputeObjClassId()
{
	return ComputeObjClassId(ObjClassId,ObjClassSubId,ObjConfirmed,m_astroObjects);
}


int Source::ComputeComponentObjClassId()
{
	//Check if has fit component
	if(!m_HasFitInfo || m_componentAstroObjects.empty()){
		return 0;
	}

	//Check for mismatch between number of fit components & astro object collection size
	int nComponents= m_fitPars.GetNComponents();
	int collectionSize= static_cast<int>(m_componentAstroObjects.size());
	if(collectionSize!=nComponents){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Mismatch between number of fit components stored (="<<nComponents<<") and component object collection size (="<<collectionSize<<")!");
		#endif
		return -1;
	}

	//Loop over components and compute 
	componentObjClassIds.clear();
	componentObjClassSubIds.clear();
	componentObjConfirmed.clear();

	for(size_t i=0;i<m_componentAstroObjects.size();i++)
	{
		int classId= eUNKNOWN_OBJECT;
		int classSubId= eUNKNOWN_OBJECT; 
		bool confirmed= false;
		ComputeObjClassId(classId,classSubId,confirmed,m_componentAstroObjects[i]);
		componentObjClassIds.push_back(classId);
		componentObjClassSubIds.push_back(classSubId);
		componentObjConfirmed.push_back(confirmed);
	}


	return 0;

}//close ComputeComponentObjClassId()


int Source::ComputeObjClassId(int& classId,int& classSubId,bool& confirmed,std::vector<AstroObject>& astroObjects)
{
	//Set to unknown if no astro object data present
	if(astroObjects.empty()){
		classId= eUNKNOWN_OBJECT; 
		classSubId= eUNKNOWN_OBJECT; 
		confirmed= false;
		return 0;
	}

	//Set to obj type if only one object data present
	//otherwise set to the major type.
	int nObjs= static_cast<int>(astroObjects.size());
	if(nObjs==1){
		classId= astroObjects[0].id;
		classSubId= astroObjects[0].subid;
		confirmed= astroObjects[0].confirmed;
	}
	else{
		std::vector<int> objIds;
		std::vector<int> objSubIds;
		std::vector<bool> objConfirmeds;		
		for(int i=0;i<nObjs;i++){
			objIds.push_back(astroObjects[i].id);
			objSubIds.push_back(astroObjects[i].subid);	
			objConfirmeds.push_back(astroObjects[i].confirmed);
		}
		int nmodes_id= 0;
		int mode_id= CodeUtils::FindVectorMode(objIds.begin(),objIds.end(),nmodes_id);
		int nmodes_subid= 0;
		int mode_subid= CodeUtils::FindVectorMode(objSubIds.begin(),objSubIds.end(),nmodes_subid);
		if(nmodes_id==1) classId= mode_id;
		else classId= eMULTI_CLASS_OBJECT;
		if(nmodes_subid==1) classSubId= mode_subid;
		else classSubId= eMULTI_CLASS_OBJECT;
		int nmodes_confirmed= 0;
		bool mode_confirmed= CodeUtils::FindVectorMode(objConfirmeds.begin(),objConfirmeds.end(),nmodes_confirmed);
		if(nmodes_confirmed==1) confirmed= mode_confirmed;
		else confirmed= false;
	}

	return 0;

}//close ComputeObjClassId()




int Source::AddAstroObject(AstroObject& astroObject)
{
	//Add if no objects present and compute obj class id
	if(m_astroObjects.empty()){
		m_astroObjects.push_back(astroObject);
		m_hasAstroObjectData= true;
		ComputeObjClassId();
		return 0;
	}

	//Search if object is already present 
	auto it= std::find_if(
		m_astroObjects.begin(),
		m_astroObjects.end(), 
		[&astroObject](const AstroObject& obj) { 
    	return obj.name == astroObject.name; 
		}
	);

	//Add to collection if not present and re-compute class id
	if(it==m_astroObjects.end()){
		m_astroObjects.push_back(astroObject);
		m_hasAstroObjectData= true;
		ComputeObjClassId();
	}

	return 0;

}//close AddAstroObject()

int Source::AddComponentAstroObject(int componentIndex,AstroObject& astroObject)
{
	//Check if has fit components
	if(!m_HasFitInfo){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source "<<this->GetName()<<" has no fit components, no astro object will be added!");
		#endif
		return -1;
	}

	//Check if component index is in range
	int nComponents= m_fitPars.GetNComponents();
	if(componentIndex<0 || componentIndex>=nComponents){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Invalid component index ("<<componentIndex<<") given when nComponents="<<nComponents<<", no astro object will be added!");
		#endif
		return -1;
	}

	//If collection size is empty allocate space for all components
	if(m_componentAstroObjects.empty()){
		for(int i=0;i<nComponents;i++){
			m_componentAstroObjects.push_back( std::vector<AstroObject>() );
		}
	}

	//Add if no objects present and compute obj class id
	if(m_componentAstroObjects[componentIndex].empty()){
		m_componentAstroObjects[componentIndex].push_back(astroObject);
		m_hasComponentAstroObjectData= true;
		ComputeComponentObjClassId();
		return 0;
	}

	//Search if object is already present 
	auto it= std::find_if(
		m_componentAstroObjects[componentIndex].begin(),
		m_componentAstroObjects[componentIndex].end(), 
		[&astroObject](const AstroObject& obj) { 
    	return obj.name == astroObject.name; 
		}
	);

	//Add to collection if not present and re-compute class id
	if(it==m_componentAstroObjects[componentIndex].end()){
		m_componentAstroObjects[componentIndex].push_back(astroObject);
		m_hasComponentAstroObjectData= true;
		ComputeComponentObjClassId();
	}
	
	return 0;

}//close AddComponentAstroObject()


}//close namespace


