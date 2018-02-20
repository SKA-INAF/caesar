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

#include <TObject.h>
#include <TMatrixD.h>

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
	DEBUG_LOG("Copy constuctor called...");
  Init();
  ((Source&)source).Copy(*this);
}

void Source::Copy(TObject &obj) const {

	//Copy mother blob 
	Blob::Copy((Source&)obj);

	// Copy this source to source obj	
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
	
	//Delete first a previously existing vector
	for(unsigned int i=0;i<(((Source&)obj).m_NestedSources).size();i++){
		if( (((Source&)obj).m_NestedSources)[i] ){
			delete (((Source&)obj).m_NestedSources)[i];
			(((Source&)obj).m_NestedSources)[i]= 0;
		}
	}
	(((Source&)obj).m_NestedSources).clear();

	((Source&)obj).m_NestedSource= 0;
	for(unsigned int i=0;i<m_NestedSources.size();i++){
		((Source&)obj).m_NestedSource= new Source;
		*(((Source&)obj).m_NestedSource)= *(m_NestedSources[i]);
		(((Source&)obj).m_NestedSources).push_back( ((Source&)obj).m_NestedSource );
	}


}//close Copy()

Source& Source::operator=(const Source& source) { 
	// Operator =
  if (this != &source)  ((Source&)source).Copy(*this);
  return *this;
}


void Source::Init(){

	//Init source flags
	Type= eUnknown;
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

}//close Init()


void Source::Draw(bool drawBoundingBox,bool drawEllipse,bool drawNested,int lineColor,int lineStyle){

	//Drawing contours?
	DEBUG_LOG("#"<<m_Contours.size()<<" contours present for source "<<Id<<"...");
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


const std::string Source::GetDS9Region(bool dumpNestedSourceInfo){

	//Check if has pixels
	//NB: DS9 crashes miserably when given a polygon region with one point 
	if(NPix<=1) return std::string("");

	//global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 image
	std::stringstream sstream;
	sstream<<"polygon ";
	for(unsigned int i=0; i<m_Contours.size(); i++){ 
		int nPoints= m_Contours[i]->GetN();
		for(int j=0;j<nPoints;j++){
			TVector2* contPnt= m_Contours[i]->GetPoint(j);
			if(!contPnt) continue;
			sstream<<(int)contPnt->X()+1<<" "<<(int)contPnt->Y()+1<<" ";
			//sstream<<(int)contPnt->X()<<" "<<(int)contPnt->Y()<<" ";
		}
	}
	//sstream<<"# text={S"<<Id<<"}";
	sstream<<"# text={"<<this->GetName()<<"}";

	if(dumpNestedSourceInfo && m_HasNestedSources){			
		sstream<<endl;
		for(unsigned int k=0;k<m_NestedSources.size();k++){
			sstream<<"polygon ";
			std::vector<Contour*> nestedContours= m_NestedSources[k]->m_Contours;
			for(unsigned int i=0; i<nestedContours.size(); i++){ 
				int nPoints= nestedContours[i]->GetN();
				for(int j=0;j<nPoints;j++){
					TVector2* contPnt= nestedContours[i]->GetPoint(j);
					if(!contPnt) continue;
					sstream<<(int)contPnt->X()+1<<" "<<(int)contPnt->Y()+1<<" ";
					//sstream<<(int)contPnt->X()<<" "<<(int)contPnt->Y()<<" ";
				}
			}//end loop contours
			//sstream<<"# text={S"<<Id<<"_Nest"<<k<<"}";
			sstream<<"# text={"<<this->GetName()<<"_Nest"<<k<<"}";
			if(k!=m_NestedSources.size()-1) sstream<<endl;
		}//end loop nested sources
	}//close dumpNestedSourceInfo

	const std::string dsregions= sstream.str();
	return dsregions;

}//close GetDS9Region()

const std::string Source::GetDS9FittedEllipseRegion(bool useFWHM)
{
	//Check if source has fit info
	const std::string dsregions= "";
	if(!m_HasFitInfo) return dsregions;

	//Loop over fitted components and get their ellipses
	std::vector<TEllipse*> ellipses= m_fitPars.GetFittedEllipses();

	//ellipse x y radius radius angle
	std::stringstream sstream;
	
	for(size_t i=0;i<ellipses.size();i++){
		if(!ellipses[i]) continue;

		//Get ellipse pars
		double x0= ellipses[i]->GetX1();
		double y0= ellipses[i]->GetY1();
		double R1= ellipses[i]->GetR1();
		double R2= ellipses[i]->GetR2();
		double theta= ellipses[i]->GetTheta();
		//theta-= 90;//DS9 format
		//sstream<<"ellipse "<<x0+1<<" "<<y0+1<<" "<<R1<<" "<<R2<<" "<<theta<<" # text={S"<<Id<<"_"<<i+1<<"}";
		sstream<<"ellipse "<<x0+1<<" "<<y0+1<<" "<<R1<<" "<<R2<<" "<<theta<<" # text={"<<this->GetName()<<"_"<<i+1<<"}";
		if(i!=ellipses.size()-1) sstream<<endl;
	}//end loop ellipses
	
	return sstream.str();

}//close GetDS9FittedEllipseRegion()


const std::string Source::GetDS9EllipseRegion(bool dumpNestedSourceInfo){
			
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

bool Source::IsAdjacentSource(Source* aSource){
		
	//Check input sources
	if(!aSource){
		ERROR_LOG("Null ptr to given input Source!");
		return false;
	}		

	//Check if pixel collections are empty
	if(GetNPixels()<=0 || aSource->GetNPixels()<=0){
		DEBUG_LOG("This or given source have no pixels, return not adjacent.");
		return false;
	}	

	//Check if bouding boxes are overlapping
	//NB: If not skip the adjacency check
	bool areBoundingBoxesOverlapping= CheckBoxOverlapping(aSource);
	if(!areBoundingBoxesOverlapping){
		DEBUG_LOG("Sources not adjacent as their bounding boxes (S1(x["<<m_Xmin<<","<<m_Xmax<<"), y["<<m_Ymin<<","<<m_Ymax<<"]), S2(x["<<aSource->m_Xmin<<","<<aSource->m_Xmax<<"), y["<<aSource->m_Ymin<<","<<aSource->m_Ymax<<"]) do not overlap...");
		return false;		
	}

	//Find if there are adjacent pixels
	DEBUG_LOG("Finding if this source (pos=("<<X0<<","<<Y0<<"), #"<<m_Pixels.size()<<" pix) is adjacent to source (pos("<<aSource->X0<<","<<aSource->Y0<<"), #"<<(aSource->m_Pixels).size()<<" pix)");
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
  	DEBUG_LOG("Sources ("<<this->Id<<","<<aSource->Id<<") are adjacent (found a match at " << std::distance((aSource->m_Pixels).begin(), it)<<")");
  }

	return isAdjacent;

}//close IsAdjacentSource()

long int Source::GetNMatchingPixels(std::vector<Pixel*>& matching_pixels,Source* aSource,bool sorted)
{
	//Initialize empty list of matching pixels
	matching_pixels.clear();

	//Check input sources
	if(!aSource){
		ERROR_LOG("Null ptr to given input Source!");
		return 0;
	}	

	//Check if pixel collections are empty
	if(GetNPixels()<=0 || aSource->GetNPixels()<=0){
		DEBUG_LOG("This or given source have no pixels, return 0.");
		return 0;
	}	
	
	//Check if bounding boxes are overlapping
	//NB: If not skip the adjacency check
	bool areBoundingBoxesOverlapping= CheckBoxOverlapping(aSource);
	if(!areBoundingBoxesOverlapping){
		DEBUG_LOG("Sources not overlapping as their bounding boxes (S1(x["<<m_Xmin<<","<<m_Xmax<<"), y["<<m_Ymin<<","<<m_Ymax<<"]), S2(x["<<aSource->m_Xmin<<","<<aSource->m_Xmax<<"), y["<<aSource->m_Ymin<<","<<aSource->m_Ymax<<"]) do not overlap...");
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
		WARN_LOG("Invalid threshold ("<<matchOverlapThr<<") given (hint: should be >0)!");
		return false;
	}
 
	//Return if given source collection to be searched is empty
	if(sources.empty()) {
		WARN_LOG("Given source collection is empty, no match can be searched!");
		return false;
	}

	//Return if this source has no pixels
	if(m_Pixels.empty() || NPix<=0){
		WARN_LOG("This source has no pixels stored, cannot search for a match with a source catalog!");
		return false;
	}	

	//Loop over sources to find matches
	//NB: Two pixel fractions are compared (M=overlap pixels, N=pixels of true source, N_rec=pixels of rec source)
	// 1) f= M/N<thr (to ensure rec source is not a tiny fraction of true source)  
	// 2) f_rec= M/N_rec<thr (to ensure rec source is not encompassing by a large amount the true source) 
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
	
			INFO_LOG("Source "<<this->GetName()<<" (X0="<<X0<<", Y0="<<Y0<<", N="<<NPix<<"): found match with source "<<sources[j]->GetName()<<" (X0="<<sources[j]->X0<<", Y0="<<sources[j]->Y0<<", N="<<sources[j]->NPix<<"), NMatchingPixels="<<NMatchingPixels<<", f="<<f<<", f_rec="<<f_rec<<" (t="<<matchOverlapThr<<"), Sratio="<<Sratio<<", Sratio_rec="<<Sratio_rec<<", dX="<<dX<<", dY="<<dY);
		}
		
	}//end loop sources	

	//Check if match is found
	if(tmpMatchPars.empty()) return false;

	//Search for best overlap in case of multiple matches
	//NB: Consider 2 rec match sources: 1st) has larger f, 2nd) has larger f_rec. Which is the best one?
	//    It is assumed here that the best is those that covers true source with a larger fraction, e.g. the 1st)
	if(tmpMatchPars.size()>1) {
		INFO_LOG("#"<<tmpMatchPars.size()<<" source matches found, searching for the best one ...");
		
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


bool Source::FindSourceMatchByPos(SourcePosMatchPars& pars, const std::vector<Source*>& sources, float matchPosThr)
{
	//Reset pars
	pars.ResetPars();

	//Check given threshold
	if(matchPosThr<=0){
		WARN_LOG("Invalid threshold ("<<matchPosThr<<") given (hint: should be >0)!");
		return false;
	}
 
	//Return if given source collection to be searched is empty
	if(sources.empty()) {
		WARN_LOG("Given source collection is empty, no match can be searched!");
		return false;
	}

	//Return if this source has no pixels
	if(m_Pixels.empty() || NPix<=0){
		WARN_LOG("This source has no pixels stored, cannot search for a match with a source catalog!");
		return false;
	}	

	//Set reference position to be compared with collection
	double X0_ref= X0;
	double Y0_ref= Y0;
	if(m_HasTrueInfo){
		X0_ref= m_X0_true;
		Y0_ref= m_Y0_true;
	}	

	//Loop over all reconstructed sources and match them in position
	std::vector<SourcePosMatchPars> tmpMatchPars;	
	
	for(size_t j=0;j<sources.size();j++){
		//If source has fit info loop over fitted components to find best match
		bool hasFitInfo= sources[j]->HasFitInfo();
		if(hasFitInfo){
			SourceFitPars fitPars= sources[j]->GetFitPars();
			double dist_best= 1.e+99;
			long int best_component= -1;
			for(int k=0;k<fitPars.GetNComponents();k++){
				double X0_fitcomp= fitPars.GetParValue(k,"x0");	
				double Y0_fitcomp= fitPars.GetParValue(k,"y0");
				double dx= fabs(X0_fitcomp-X0_ref);
				double dy= fabs(Y0_fitcomp-Y0_ref);
				double dist= sqrt(dx*dx+dy*dy);
				if(dx<=matchPosThr && dy<=matchPosThr && dist<dist_best){//match found
					best_component= k;
					dist_best= dist;
				}
			}//end loop fitted components	
			
			if(best_component!=-1){
				tmpMatchPars.push_back(SourcePosMatchPars(j,dist_best,best_component));	
			}
	
		}//close if has fit info
		else{
			//No fit info available so compute offset using source barycenter 
			double dx= fabs(sources[j]->X0 - X0_ref);
			double dy= fabs(sources[j]->Y0 - Y0_ref);
			double dist= sqrt(dx*dx+dy*dy);
			if(dx<=matchPosThr && dy<=matchPosThr){//match found
				tmpMatchPars.push_back(SourcePosMatchPars(j,dist,-1));
			}
		}//close else !has fit info
	}//end loop sources

	//Check if match is found
	if(tmpMatchPars.empty()) return false;

	//Search for best overlap in case of multiple matches
	if(tmpMatchPars.size()>1) {
		INFO_LOG("#"<<tmpMatchPars.size()<<" source matches found, searching for the best one ...");
		
		double dist_best= 1.e+99;
		int best_index= -1;

		for(size_t j=0;j<tmpMatchPars.size();j++){
			double dist= tmpMatchPars[j].posDiff;	
			if(dist<=dist_best){
				best_index= j;
				dist_best= dist;
			}
		}//end loop matched sources

		pars= tmpMatchPars[best_index];

	}//close if
	else{
		pars= tmpMatchPars[0];
	}

	return true;

}//close FindSourceMatchByPos()


int Source::MergeSource(Source* aSource,bool copyPixels,bool checkIfAdjacent,bool computeStatPars,bool computeMorphPars,bool sumMatchingPixels){

	//Check input sources
	if(!aSource){
		ERROR_LOG("Null ptr to given input Source!");
		return -1;
	}	

	//If adjacency check is enabled check if sources are mergeable (e.g. if there is at least
	//a pixel adjacent to each other)
	if(checkIfAdjacent) {
		bool areAdjacent= IsAdjacentSource(aSource);
		if(!areAdjacent){
			WARN_LOG("Sources are not adjacent nor overlapping, no merging will be performed!");
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
		WARN_LOG("No pixels to be merged (perfectly overlapping sources?)!");
		//return -1;
	}
	else{
		INFO_LOG("# "<<nMergedPixels<<"/"<<aSource->GetNPixels()<<" pixels to be merged to this source (Npix="<<this->GetNPixels()<<")...");
	}

	//Now merge the pixels (if any to be merged)
	//Source moments will be updated in AddPixel()
	for(int i=0;i<nMergedPixels;i++){
		if(this->AddPixel(pixelsToBeMerged[i],copyPixels)<0){
			WARN_LOG("Failed to add pixel no. "<<i<<" to this source, skip to next pixel in list...");
			continue;			
		}
	}//end loop pixels to be merged

	
	
	//Set new source type
	int mergedSourceType= aSource->Type;
	if( this->Type==eCompact || this->Type==ePointLike ){
		if(mergedSourceType==eCompact || mergedSourceType==ePointLike) this->Type= eCompact; 
		else if(mergedSourceType==eExtended || mergedSourceType==eCompactPlusExtended) this->Type= eCompactPlusExtended; 
		else this->Type= eUnknown; 
	}
	else if( this->Type==eExtended ){
		if(mergedSourceType==eCompact || mergedSourceType==ePointLike) this->Type= eCompactPlusExtended; 
		else if(mergedSourceType==eExtended) this->Type= eExtended;  
		else if(mergedSourceType==eCompactPlusExtended) this->Type= eCompactPlusExtended; 
		else this->Type= eUnknown; 
	}
	else if( this->Type==eCompactPlusExtended ){
		if(mergedSourceType==eUnknown) this->Type= eUnknown; 
		else this->Type= eCompactPlusExtended;  
	}
	else{
		this->Type= eUnknown; 
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
				INFO_LOG("Changing simtype from (simtype1="<<SimType<<", simtype2="<<aSource->SimType<<") to simtype="<<simtype_merged);
				SimType= simtype_merged;
			}
			catch(...){
				ERROR_LOG("C++ exception occurred while converting merged stringified simtype "<<simtype12_str<<" to int code (will not add merged soure to simtype!");
			}
		}//close if found
	}//close if sim type

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

	return 0;
	
}//close MergeSource()

bool Source::CheckBoxOverlapping(Source* aSource)
{
	//Check input blob
	if(!aSource){
		ERROR_LOG("Null ptr to input source given!");
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
		ERROR_LOG("Null ptr to input source given, returning inf!");
		return std::numeric_limits<float>::infinity();
	}

	//If one or both sources has no pars computed return inf
	if(!this->HasStats() || aSource->HasStats()){
		WARN_LOG("One or both sources has no stats computed (you must compute pars before getting distance, returning inf)");
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
	//if(fitter.FitSource(this,blobPars,nMaxComponents)<0){
	if(fitter.FitSource(this,fitOptions)<0){	
		WARN_LOG("Failed to fit source "<<this->GetName()<<" ...");
		return -1;
	}
	
	//Get fit results
	SourceFitPars fitPars= fitter.GetFitPars();
	fitPars.Print();

	int fitStatus= fitter.GetFitStatus();
	if(fitStatus==SourceFitter::eFitConverged || fitStatus==SourceFitter::eFitConvergedWithWarns){
		INFO_LOG("Fit of source "<<this->GetName()<<" converged (status="<<fitStatus<<"), storing fit parameters...");	
		m_fitPars= fitPars;
		m_HasFitInfo= true;
	}

	return 0;

}//close Fit()


}//close namespace


