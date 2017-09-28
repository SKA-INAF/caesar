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
* @file Serializer.cc
* @class Serializer
* @brief Serializer class
*
* Class for serializing objects
* @author S. Riggi
* @date 20/01/2015
*/

#include <Serializer.h>
#include <Source.pb.h>
#include <TaskData.pb.h>

#include <Logger.h>
#include <CodeUtils.h>

#include <TVector2.h>

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <exception>

#include <chrono>

using namespace std;

ClassImp(Caesar::SBuffer)
//ClassImp(Caesar::Serializer)

namespace Caesar {

Serializer::Serializer() 
{	
	
}//close constructor


Serializer::~Serializer()
{

}//close destructor

int Serializer::EncodeSourceComponentParsToProtobuf(CaesarPB::SourceComponentPars* sourceCompPars_pb,SourceComponentPars& sourceCompPars)
{
	//Check pointer
	if(!sourceCompPars_pb){
		ERROR_LOG("Null ptr given to source fit pars protobuf object!");
		return -1;
	}

	//Loop over fit parameters and add to protobuf object
	try {
		std::map<std::string,double> FitPars= sourceCompPars.GetFitPars();
		std::map<std::string,double> FitParsErr= sourceCompPars.GetFitParErrors();

		std::map<std::string,double>::iterator it = FitPars.begin();
		for (it=FitPars.begin(); it!=FitPars.end(); ++it){
			std::string parName= it->first;
			double parValue= it->second;
			CaesarPB::Dict* fitPar_pb = sourceCompPars_pb->add_fitpars();
			fitPar_pb->set_key(parName);
			fitPar_pb->set_val(parValue);
		}//end loop pars

		for (it=FitParsErr.begin(); it!=FitParsErr.end(); ++it){
			std::string parName= it->first;
			double parErr= it->second;
			CaesarPB::Dict* fitParErr_pb = sourceCompPars_pb->add_fitparserr();
			fitParErr_pb->set_key(parName);
			fitParErr_pb->set_val(parErr);
		}//end loop par errors

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Source component pars encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeSourceComponentParsToProtobuf()

int Serializer::EncodeSourceFitParsToProtobuf(CaesarPB::SourceFitPars& sourceFitPars_pb,SourceFitPars& sourceFitPars)
{
	//Add main pars
	sourceFitPars_pb.set_ncomponents(sourceFitPars.GetNComponents());
	sourceFitPars_pb.set_chi2(sourceFitPars.GetChi2());
	sourceFitPars_pb.set_ndof(sourceFitPars.GetNDF());
	sourceFitPars_pb.set_npars_free(sourceFitPars.GetNFreePars());
	sourceFitPars_pb.set_nfit_points(sourceFitPars.GetNFitPoints());
	sourceFitPars_pb.set_status(sourceFitPars.GetStatus());
	sourceFitPars_pb.set_minimizer_status(sourceFitPars.GetMinimizerStatus());
	sourceFitPars_pb.set_offset(sourceFitPars.GetOffsetPar());
	sourceFitPars_pb.set_offset_err(sourceFitPars.GetOffsetParErr());

	//Add component fit pars
	std::vector<SourceComponentPars> pars= sourceFitPars.GetPars();
	for(size_t i=0;i<pars.size();i++){
		CaesarPB::SourceComponentPars* compFitPars_pb= sourceFitPars_pb.add_pars();
		if(EncodeSourceComponentParsToProtobuf(compFitPars_pb,pars[i])<0){
			ERROR_LOG("Failed to encode fit pars for component "<<i+1<<" in protobuf!");
			return -1;
		}
	}//end loop component fit pars

	return 0;

}//close EncodeSourceFitParsToProtobuf()


int Serializer::EncodePointToProtobuf(CaesarPB::Point& point_pb,TVector2& point){

	try {
		point_pb.set_x(point.X());
		point_pb.set_y(point.Y());
	}
	catch(std::exception const & e) {
		ERROR_LOG("Point encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodePointToProtobuf()

int Serializer::EncodeContourToProtobuf(CaesarPB::Contour& contour_pb,Contour* contour){
	
	if(!contour) return -1;

	try {
		contour_pb.set_hasparameters(contour->HasParameters);
		contour_pb.set_area(contour->Area);
		contour_pb.set_perymeter(contour->Perymeter);
		contour_pb.set_isconvexcontour(contour->IsConvexContour);
		contour_pb.set_circularityratio(contour->CircularityRatio);
			
		CaesarPB::Point* BoundingBoxCenter= new CaesarPB::Point;
		if(EncodePointToProtobuf(*BoundingBoxCenter,contour->BoundingBoxCenter)<0){
			throw std::runtime_error("Failed to encode bounding box center field");
		}
		contour_pb.set_allocated_boundingboxcenter(BoundingBoxCenter);

		contour_pb.set_boundingboxmaj(contour->BoundingBoxMaj);
		contour_pb.set_boundingboxmin(contour->BoundingBoxMin);
		contour_pb.set_boundingboxangle(contour->BoundingBoxAngle);
		contour_pb.set_elongation(contour->Elongation);
		contour_pb.set_rectangularity(contour->Rectangularity);
		contour_pb.set_roundness(contour->Roundness);
		contour_pb.set_eccentricity(contour->Eccentricity);
		contour_pb.set_tiltangle(contour->TiltAngle);
		contour_pb.set_hasellipsefit(contour->HasEllipseFit);

		CaesarPB::Point* EllipseCenter= new CaesarPB::Point;
		if(EncodePointToProtobuf(*EllipseCenter,contour->EllipseCenter)<0){
			throw std::runtime_error("Failed to encode ellipse center field");
		}
		contour_pb.set_allocated_ellipsecenter(EllipseCenter);

		contour_pb.set_ellipsemajaxis(contour->EllipseMajAxis);
		contour_pb.set_ellipseminaxis(contour->EllipseMinAxis);
		contour_pb.set_ellipserotangle(contour->EllipseRotAngle);
		contour_pb.set_ellipsefitredchi2(contour->EllipseFitRedChi2);
		contour_pb.set_ellipsearearatio(contour->EllipseAreaRatio);

		contour_pb.set_m00(contour->m00);
		contour_pb.set_m10(contour->m10);
		contour_pb.set_m01(contour->m01);
		contour_pb.set_m20(contour->m20);
		contour_pb.set_m11(contour->m11);
		contour_pb.set_m02(contour->m02);
		contour_pb.set_m30(contour->m30);
		contour_pb.set_m21(contour->m21);
		contour_pb.set_m12(contour->m12);
		contour_pb.set_m03(contour->m03);

		contour_pb.set_mu20(contour->mu20);
		contour_pb.set_mu11(contour->mu11);
		contour_pb.set_mu02(contour->mu02);
		contour_pb.set_mu30(contour->mu30);
		contour_pb.set_mu21(contour->mu21);
		contour_pb.set_mu12(contour->mu12);
		contour_pb.set_mu03(contour->mu03);
		
		contour_pb.set_nu20(contour->nu20);
		contour_pb.set_nu11(contour->nu11);
		contour_pb.set_nu02(contour->nu02);
		contour_pb.set_nu30(contour->nu30);
		contour_pb.set_nu21(contour->nu21);
		contour_pb.set_nu12(contour->nu12);
		contour_pb.set_nu03(contour->nu03);

		for(unsigned int i=0;i<contour->HuMoments.size();i++){		
 			contour_pb.add_humoments(contour->HuMoments[i]);
		}
	
		for(unsigned int i=0;i<contour->BoundingBoxVertex.size();i++){		
			CaesarPB::Point* thisBoundingBoxVertex = contour_pb.add_boundingboxvertex();
			if(EncodePointToProtobuf(*thisBoundingBoxVertex,(contour->BoundingBoxVertex)[i])<0){
				std::stringstream errMsg;
				errMsg<<"BoundingBoxVertex no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	
		CaesarPB::Point* Centroid= new CaesarPB::Point;
		//CaesarPB::Point Centroid;
		if(EncodePointToProtobuf(*Centroid,contour->Centroid)<0){
		//if(EncodePointToProtobuf(Centroid,contour->Centroid)<0){
			throw std::runtime_error("Failed to encode ellipse center field");
		}
		contour_pb.set_allocated_centroid(Centroid);

		contour_pb.set_hasfdpars(contour->HasFDPars);
		for(size_t i=0;i<contour->RealFDs.size();i++){		
 			contour_pb.add_realfds(contour->RealFDs[i]);
		}
		for(size_t i=0;i<contour->ImagFDs.size();i++){		
 			contour_pb.add_imagfds(contour->ImagFDs[i]);
		}
		for(size_t i=0;i<contour->ModFDs.size();i++){		
 			contour_pb.add_modfds(contour->ModFDs[i]);
		}

		contour_pb.set_hasbepars(contour->HasBEPars);
		for(size_t i=0;i<contour->BendingEnergies.size();i++){		
 			contour_pb.add_bendingenergies(contour->BendingEnergies[i]);
		}	

		contour_pb.set_hascentroiddistancefdpars(contour->HasCentroidDistanceFDPars);
		for(size_t i=0;i<contour->CentroidDistanceModFDs.size();i++){		
 			contour_pb.add_centroiddistancemodfds(contour->CentroidDistanceModFDs[i]);
		}

		for(int i=0;i<contour->GetN();i++){		
			CaesarPB::Point* thisPoint = contour_pb.add_m_points();
			if(EncodePointToProtobuf(*thisPoint,*contour->GetPoint(i))<0){
				std::stringstream errMsg;
				errMsg<<"Point no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Point encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeContourToProtobuf()


int Serializer::EncodePixelToProtobuf(CaesarPB::Pixel& pixel_pb,Pixel* pixel){
	
	if(!pixel) return -1;

	try {
		//Set public source fields
		pixel_pb.set_id(pixel->id);
		if(pixel->type==Pixel::eNormal){
			pixel_pb.set_type(CaesarPB::Pixel::eNormal);
		}
		else if(pixel->type==Pixel::eSeed){
			pixel_pb.set_type(CaesarPB::Pixel::eSeed);
		}
		else if(pixel->type==Pixel::eHalo){
			pixel_pb.set_type(CaesarPB::Pixel::eHalo);
		}
		else{
			std::stringstream errMsg;
			errMsg<<"Invalid pixel type ("<<pixel->type<<"), failed to encode!";
			throw std::runtime_error(errMsg.str().c_str());
		}

		pixel_pb.set_s(pixel->S);
		pixel_pb.set_x(pixel->x);
		pixel_pb.set_y(pixel->y);
		pixel_pb.set_ix(pixel->ix);
		pixel_pb.set_iy(pixel->iy);
		//pixel_pb.set_isonedge(pixel->isOnEdge);
		//pixel_pb.set_distancetoedge(pixel->distanceToEdge);

		pixel_pb.set_s_curv(pixel->GetCurv());
		pixel_pb.set_s_edge(pixel->GetEdge());
		
		std::pair<double,double> bkginfo= pixel->GetBkg();
		pixel_pb.set_bkglevel(bkginfo.first);
		pixel_pb.set_noiselevel(bkginfo.second);

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Pixel encoding to protobuf failed with status "<<e.what());
		return -1;
	}


	return 0;
}

int Serializer::EncodeBlobToProtobuf(CaesarPB::Blob& blob_pb,Source* source){
	
	if(!source) return -1;

	try{
		blob_pb.set_haspixelsatedge(source->HasPixelsAtEdge);
		blob_pb.set_id(source->Id);
		blob_pb.set_name(source->GetName());
		blob_pb.set_npix(source->NPix);
		blob_pb.set_mean(source->Mean);
		blob_pb.set_rms(source->RMS);
		blob_pb.set_skewness(source->Skewness);	
		blob_pb.set_median(source->Median);
		blob_pb.set_medianrms(source->MedianRMS);
		blob_pb.set_x0(source->X0);
		blob_pb.set_y0(source->Y0);
		blob_pb.set_mean_curv(source->Mean_curv);
		blob_pb.set_rms_curv(source->RMS_curv);
		blob_pb.set_median_curv(source->Median_curv);
		blob_pb.set_medianrms_curv(source->MedianRMS_curv);
		
		for(unsigned int i=0;i<source->Moments.size();i++) blob_pb.add_moments((source->Moments)[i]);
		for(unsigned int i=0;i<source->HuMoments.size();i++) blob_pb.add_humoments((source->HuMoments)[i]);
		for(unsigned int i=0;i<source->ZMMoments.size();i++) blob_pb.add_zmmoments((source->ZMMoments)[i]);
		blob_pb.set_m_hasstats(source->HasStats());	
		blob_pb.set_m_hasparameters(source->HasParameters());
		blob_pb.set_m_m1(source->GetM1());
		blob_pb.set_m_m2(source->GetM2());
		blob_pb.set_m_m3(source->GetM3());
		blob_pb.set_m_m4(source->GetM4());
		blob_pb.set_m_m1_curv(source->GetM1Curv());
		blob_pb.set_m_m2_curv(source->GetM2Curv());

		blob_pb.set_m_s(source->GetS());
		blob_pb.set_m_smax(source->GetSmax());
		blob_pb.set_m_smin(source->GetSmin());
		blob_pb.set_m_sxx(source->GetSxx());
		blob_pb.set_m_syy(source->GetSyy());
		blob_pb.set_m_sxy(source->GetSxy());
		blob_pb.set_m_sx(source->GetSx());
		blob_pb.set_m_sy(source->GetSy());
		blob_pb.set_m_pixidmax(source->GetSmaxPixId());
		blob_pb.set_m_pixidmin(source->GetSminPixId());
		blob_pb.set_m_s_curv(source->GetScurv());
		blob_pb.set_m_s_edge(source->GetSedge());
	
	
		
		/*
		//===== MARKED FOR REMOVAL ====
		float xmin, xmax, ymin, ymax;
		source->GetImageRange(xmin,xmax,ymin,ymax);
		blob_pb.set_m_imageminx(xmin);
		blob_pb.set_m_imagemaxx(xmax);
		blob_pb.set_m_imageminy(ymin);
		blob_pb.set_m_imagemaxy(ymax);

		long int imgsizex, imgsizey;
		source->GetImageSize(imgsizex,imgsizey);
		blob_pb.set_m_imagesizex(imgsizex);
		blob_pb.set_m_imagesizey(imgsizey);

		double imgSmin, imgSmax;
		source->GetImageSRange(imgSmin, imgSmax);
		blob_pb.set_m_imagemins(imgSmin);
		blob_pb.set_m_imagemaxs(imgSmax);

		double imgSmin_curv, imgSmax_curv;
		source->GetImageScurvRange(imgSmin_curv, imgSmax_curv);
		blob_pb.set_m_imageminscurv(imgSmin_curv);
		blob_pb.set_m_imagemaxscurv(imgSmax_curv);

		double imgSmin_edge, imgSmax_edge;
		source->GetImageSedgeRange(imgSmin_edge, imgSmax_edge);
		blob_pb.set_m_imageminsedge(imgSmin_edge);
		blob_pb.set_m_imagemaxsedge(imgSmax_edge);

		blob_pb.set_m_imagerms(source->GetImageRMS());
		//=================================
		*/

		float sxmin, sxmax, symin, symax;
		source->GetSourceRange(sxmin,sxmax,symin,symax);
		blob_pb.set_m_xmin(sxmin);
		blob_pb.set_m_xmax(sxmax);
		blob_pb.set_m_ymin(symin);
		blob_pb.set_m_ymax(symax);

		long int ixmin, ixmax, iymin, iymax;
		source->GetSourcePixelRange(ixmin, ixmax, iymin, iymax);
		blob_pb.set_m_ix_min(ixmin);
		blob_pb.set_m_ix_max(ixmax);
		blob_pb.set_m_iy_min(iymin);
		blob_pb.set_m_iy_max(iymax);
		
		//Add pixel collection to blob
		for(int k=0;k<source->GetNPixels();k++){
			CaesarPB::Pixel* thisPixel = blob_pb.add_m_pixels();
			if(EncodePixelToProtobuf(*thisPixel,source->GetPixel(k))<0){
				std::stringstream errMsg;
				errMsg<<"Pixel no. "<<k+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		//Add contout to blob
		for(unsigned int k=0;k<source->GetContours().size();k++){
			CaesarPB::Contour* thisContour = blob_pb.add_m_contours();
			if(EncodeContourToProtobuf(*thisContour,source->GetContour(k))<0){
				std::stringstream errMsg;
				errMsg<<"Contour no. "<<k+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}	
		}//end loop contours

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Blob encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeBlobToProtobuf()


int Serializer::EncodeSourceToProtobuf(CaesarPB::Source& source_pb,Source* source){

	if(!source) return -1;

	try {
		
		source_pb.set_type(source->Type);
		source_pb.set_flag(source->Flag);
	
		//Set private source fields
		source_pb.set_m_beamfluxintegral(source->GetBeamFluxIntegral());
		source_pb.set_m_isgoodsource(source->IsGoodSource());
		source_pb.set_m_depthlevel(source->GetDepthLevel());
		source_pb.set_m_hasnestedsources(source->HasNestedSources());

		//Set true info	
		source_pb.set_m_hastrueinfo(source->HasTrueInfo());
		source_pb.set_m_s_true(source->GetTrueFlux());
		double x0_true, y0_true;
		source->GetTruePos(x0_true,y0_true);	
		source_pb.set_m_x0_true(x0_true);
		source_pb.set_m_y0_true(y0_true);

		//Set fit pars
		SourceFitPars fitPars= source->GetFitPars();
		CaesarPB::SourceFitPars* fitPars_pb= new CaesarPB::SourceFitPars;
		if(EncodeSourceFitParsToProtobuf(*fitPars_pb,fitPars)<0){
			throw std::runtime_error("Failed to encode source fit pars!");
		}
		
		
		//Set blob field
		CaesarPB::Blob* blob= new CaesarPB::Blob;
		if(EncodeBlobToProtobuf(*blob,source)<0){
			throw std::runtime_error("Failed to encode blob!");
		}
		source_pb.set_allocated_blob(blob);

		//Add nested sources
		for(int i=0;i<source->GetNestedSourceNumber();i++){
			CaesarPB::Source* thisNestedSource = source_pb.add_m_nestedsources();
			if(EncodeSourceToProtobuf(*thisNestedSource,source->GetNestedSource(i))<0){
				std::stringstream errMsg;
				errMsg<<"Nested source no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeSourceToProtobuf()


int Serializer::EncodeTaskDataToProtobuf(CaesarPB::TaskData& taskData_pb,TaskData* taskData){

	//## Check input data
	if(!taskData) return -1;

	try {
		//Fill task info
		//taskData_pb.set_filename(taskData->filename);
		//taskData_pb.set_jobid(taskData->jobId);
		taskData_pb.set_workerid(taskData->workerId);	
		//taskData_pb.set_taskid(taskData->taskId);
		taskData_pb.set_ix_min(taskData->ix_min);
		taskData_pb.set_ix_max(taskData->ix_max);
		taskData_pb.set_iy_min(taskData->iy_min);
		taskData_pb.set_iy_max(taskData->iy_max);
		//taskData_pb.set_x_min(taskData->x_min);
		//taskData_pb.set_x_max(taskData->x_max);
		//taskData_pb.set_y_min(taskData->y_min);
		//taskData_pb.set_y_max(taskData->y_max);

		//Fill neighbor list
		for(unsigned int i=0;i<taskData->neighborTaskId.size();i++) {
			taskData_pb.add_neighbortaskid((taskData->neighborTaskId)[i]);
		}
		for(unsigned int i=0;i<taskData->neighborWorkerId.size();i++) {
			taskData_pb.add_neighborworkerid((taskData->neighborWorkerId)[i]);
		}

		//Set source collections
		for(unsigned int i=0;i<(taskData->sources).size();i++){
			CaesarPB::Source* thisCaesarPB = taskData_pb.add_sources();
			if(EncodeSourceToProtobuf(*thisCaesarPB,(taskData->sources)[i])<0){
				std::stringstream errMsg;
				errMsg<<"Encoding of source no. "<<i+1<<" in collection to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		for(unsigned int i=0;i<(taskData->sources_edge).size();i++){
			CaesarPB::Source* thisCaesarPB = taskData_pb.add_sources_edge();
			if(EncodeSourceToProtobuf(*thisCaesarPB,(taskData->sources_edge)[i])<0){
				std::stringstream errMsg;
				errMsg<<"Encoding of edge source no. "<<i+1<<" in collection to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		for(unsigned int i=0;i<(taskData->ext_sources).size();i++){
			CaesarPB::Source* thisCaesarPB = taskData_pb.add_ext_sources();
			if(EncodeSourceToProtobuf(*thisCaesarPB,(taskData->ext_sources)[i])<0){
				std::stringstream errMsg;
				errMsg<<"Encoding of extended source no. "<<i+1<<" in collection to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		for(unsigned int i=0;i<(taskData->ext_sources_edge).size();i++){
			CaesarPB::Source* thisCaesarPB = taskData_pb.add_ext_sources_edge();
			if(EncodeSourceToProtobuf(*thisCaesarPB,(taskData->ext_sources_edge)[i])<0){
				std::stringstream errMsg;
				errMsg<<"Encoding of edge extended source no. "<<i+1<<" in collection to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Task data encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeTaskDataToProtobuf()

int Serializer::EncodeTaskDataCollectionToProtobuf(CaesarPB::TaskDataCollection& taskDataCollection_pb,std::vector<TaskData*> taskDataCollection){

	//Fill task collections
	for(size_t i=0;i<taskDataCollection.size();i++){
		CaesarPB::TaskData* thisTaskPB = taskDataCollection_pb.add_tasks();
		if(EncodeTaskDataToProtobuf(*thisTaskPB,taskDataCollection[i])<0){
			std::stringstream errMsg;
			errMsg<<"Encoding of task data no. "<<i+1<<" in collection to protobuf failed!";
			throw std::runtime_error(errMsg.str().c_str());
		}
	}

	return 0;

}//close EncodeTaskDataCollectionToProtobuf()


int Serializer::SourceToBuffer(SBuffer& buffer,Source* source){

	//## Check input source
	if(!source) return -1;

	try {
		//## Create google protobuf source message
		CaesarPB::Source source_pb;
		if(EncodeSourceToProtobuf(source_pb,source)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer.size = source_pb.ByteSize();
		buffer.data= source_pb.SerializeAsString();		
		
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close SourceToBuffer()


int Serializer::TaskDataToBuffer(SBuffer& buffer,TaskData* taskData){

	//## Check input data
	if(!taskData) {
		return -1;
	}

	try {
		//## Create google protobuf source message
		CaesarPB::TaskData taskData_pb;
		if(EncodeTaskDataToProtobuf(taskData_pb,taskData)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer.size = taskData_pb.ByteSize();
		buffer.data= taskData_pb.SerializeAsString();		
		
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Task data encoding failed with status "<<e.what());
		return -1;
	}
	return 0;

}//close TaskDataToBuffer()

int Serializer::TaskDataCollectionToBuffer(SBuffer& buffer,std::vector<TaskData*> taskDataCollection){

	try {
		//## Create google protobuf source message
		CaesarPB::TaskDataCollection taskDataCollection_pb;
		if(EncodeTaskDataCollectionToProtobuf(taskDataCollection_pb,taskDataCollection)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer.size = taskDataCollection_pb.ByteSize();
		buffer.data= taskDataCollection_pb.SerializeAsString();		
		
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Task data collection encoding failed with status "<<e.what());
		return -1;
	}
	return 0;

}//close TaskDataCollectionToBuffer()


char* Serializer::TaskDataToCharArray(long int& buffer_size,TaskData* taskData){

	//## Check input data
	if(!taskData) {
		return 0;
	}

	char* buffer= 0;

	try {
		//## Create google protobuf source message
		CaesarPB::TaskData taskData_pb;
		if(EncodeTaskDataToProtobuf(taskData_pb,taskData)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer_size = taskData_pb.ByteSize();
		if(!buffer) buffer = (char*)malloc(buffer_size);
		taskData_pb.SerializeToArray(buffer, buffer_size);
		
	}//close try blocks
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed with status "<<e.what());
		return 0;
	}

	return buffer;

}//close TaskDataToCharArray()


char* Serializer::TaskDataCollectionToCharArray(long int& buffer_size,std::vector<TaskData*> taskDataCollection){

	char* buffer= 0;

	try {
		//## Create google protobuf source message
		CaesarPB::TaskDataCollection taskDataCollection_pb;
		if(EncodeTaskDataCollectionToProtobuf(taskDataCollection_pb,taskDataCollection)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer_size = taskDataCollection_pb.ByteSize();
		if(!buffer) buffer = (char*)malloc(buffer_size);
		taskDataCollection_pb.SerializeToArray(buffer, buffer_size);
		
	}//close try blocks
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed with status "<<e.what());
		return 0;
	}

	return buffer;

}//close TaskDataCollectionToCharArray()


int Serializer::EncodeProtobufToSourceComponentPars(SourceComponentPars& sourceComponentPars,const CaesarPB::SourceComponentPars& sourceComponentPars_pb)
{

	try {
		for(int i=0;i<sourceComponentPars_pb.fitpars_size();i++){
			const CaesarPB::Dict& par_pb = sourceComponentPars_pb.fitpars(i);
			const CaesarPB::Dict& parErr_pb = sourceComponentPars_pb.fitparserr(i);

			std::string parName= par_pb.key();
			float parVal= par_pb.val();
			float parErr= parErr_pb.val();
			sourceComponentPars.SetParValueAndError(parName,parVal,parErr);
		}//end loop fit pars

	}//close try blocl
	catch(std::exception const & e) {
		ERROR_LOG("Encoding protobuf to source component pars failed with status "<<e.what());
		return -1;
	}
	
	return 0;

}//close EncodeProtobufToSourceComponentPars()

int Serializer::EncodeProtobufToSourceFitPars(SourceFitPars& sourceFitPars,CaesarPB::SourceFitPars& sourceFitPars_pb)
{
	try {	
		if(sourceFitPars_pb.has_ncomponents()) sourceFitPars.SetNComponents(sourceFitPars_pb.ncomponents());	
		if(sourceFitPars_pb.has_chi2()) sourceFitPars.SetChi2(sourceFitPars_pb.chi2());	
		if(sourceFitPars_pb.has_ndof()) sourceFitPars.SetNDF(sourceFitPars_pb.ndof());	
		if(sourceFitPars_pb.has_npars_free()) sourceFitPars.SetNFreePars(sourceFitPars_pb.npars_free());	
		if(sourceFitPars_pb.has_nfit_points()) sourceFitPars.SetNFitPoints(sourceFitPars_pb.nfit_points());	
		if(sourceFitPars_pb.has_status()) sourceFitPars.SetStatus(sourceFitPars_pb.status());	
		if(sourceFitPars_pb.has_minimizer_status()) sourceFitPars.SetMinimizerStatus(sourceFitPars_pb.minimizer_status());	
		if(sourceFitPars_pb.has_offset()) sourceFitPars.SetOffsetPar(sourceFitPars_pb.offset());	
		if(sourceFitPars_pb.has_offset_err()) sourceFitPars.SetOffsetParErr(sourceFitPars_pb.offset_err());	
		
		//Set fit component pars	
		for(int i=0;i<sourceFitPars_pb.pars_size();i++){
			const CaesarPB::SourceComponentPars& componentPars_pb = sourceFitPars_pb.pars(i);
			SourceComponentPars componentPars;
			if(EncodeProtobufToSourceComponentPars(componentPars,componentPars_pb)<0){
				throw std::runtime_error("Failed to encode component fit pars!");
			}
			sourceFitPars.SetComponentPars(i,componentPars);

		}//end loop fit component pars
	
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Encoding protobuf to source fit pars failed with status "<<e.what());
		return -1;
	}
	return 0;

}//close EncodeProtobufToSourceFitPars()



int Serializer::EncodeProtobufToSource(Source& source,const CaesarPB::Source& source_pb){

	try {
		//Set public fields
		if(source_pb.has_type()) source.Type= source_pb.type();
		if(source_pb.has_flag()) source.Flag= source_pb.flag();
	
		//Set private source fields
		if(source_pb.has_m_beamfluxintegral()) source.SetBeamFluxIntegral(source_pb.m_beamfluxintegral());
		if(source_pb.has_m_isgoodsource()) source.SetGoodSourceFlag(source_pb.m_isgoodsource());
		if(source_pb.has_m_depthlevel()) source.SetDepthLevel(source_pb.m_depthlevel());
		if(source_pb.has_m_hasnestedsources()) source.SetHasNestedSources(source_pb.m_hasnestedsources());

		if(source_pb.has_m_hastrueinfo()) {
			source.SetTrueInfo(source_pb.m_s_true(),source_pb.m_x0_true(),source_pb.m_y0_true());
		}
		
		
		//Set blob fields
		if(source_pb.has_blob()){
			const CaesarPB::Blob& blob_pb= source_pb.blob();
			if(EncodeProtobufToBlob(source,blob_pb)<0){
				ERROR_LOG("Blob encoding failed!");
				return -1;
			}
		}//close if has blob
	

		//Set nested sources (to be done after setting blob fields)
		Source* aNestedSource= 0;
		std::vector<Source*> nested_sources;
		for(int i=0;i<source_pb.m_nestedsources_size();i++){
			const CaesarPB::Source& thisNestedCaesarPB = source_pb.m_nestedsources(i);
			aNestedSource= new Source;
			if(EncodeProtobufToSource(*aNestedSource,thisNestedCaesarPB)<0){
				//Clear sources				
				for(unsigned int k=0;k<nested_sources.size();k++){
					delete nested_sources[k];
					nested_sources[k]= 0;	
				}
				nested_sources.clear();
	
				std::stringstream errMsg;
				errMsg<<"Nested source no. "<<i+1<<" encoding from protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
			nested_sources.push_back(aNestedSource);
			source.AddNestedSource(aNestedSource);
		}//end loop nested source


	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Encoding protobuf to source failed with status "<<e.what());
		return -1;
	}


	return 0;

}//close EncodeProtobufToSource()



int Serializer::EncodeProtobufToBlob(Source& source,const CaesarPB::Blob& blob_pb){
		
	try{

		//## Set public fields
		//Main params
		if(blob_pb.has_haspixelsatedge()) source.HasPixelsAtEdge= blob_pb.haspixelsatedge();
		if(blob_pb.has_id()) source.Id= blob_pb.id();
		if(blob_pb.has_name()) source.SetName(blob_pb.name());

		if(blob_pb.has_npix()) source.NPix= blob_pb.npix();		
		if(blob_pb.has_mean()) source.Mean= blob_pb.mean();
		if(blob_pb.has_rms()) source.RMS= blob_pb.rms();
		if(blob_pb.has_skewness()) source.Skewness= blob_pb.skewness();
		if(blob_pb.has_median()) source.Median= blob_pb.median();
		if(blob_pb.has_medianrms()) source.MedianRMS= blob_pb.medianrms();
		if(blob_pb.has_x0()) source.X0= blob_pb.x0();
		if(blob_pb.has_y0()) source.Y0= blob_pb.y0();
		
		//Curvature Moments
		if(blob_pb.has_mean_curv()) source.Mean_curv= blob_pb.mean_curv();
		if(blob_pb.has_rms_curv()) source.RMS_curv= blob_pb.rms_curv();
		if(blob_pb.has_median_curv()) source.Median_curv= blob_pb.median_curv();
		if(blob_pb.has_medianrms_curv()) source.MedianRMS_curv= blob_pb.medianrms_curv();
	
		//2D morphological pars
		for(int i=0;i<blob_pb.moments_size();i++){
			(source.Moments).push_back(blob_pb.moments(i));
		}
		for(int i=0;i<blob_pb.humoments_size();i++){
			(source.HuMoments).push_back(blob_pb.humoments(i));
		}
		for(int i=0;i<blob_pb.zmmoments_size();i++){
			(source.ZMMoments).push_back(blob_pb.zmmoments(i));
		}

		//## Set private fields
		if(blob_pb.has_m_hasstats()) source.SetHasStats(blob_pb.m_hasstats());
		if(blob_pb.has_m_hasparameters()) source.SetHasParameters(blob_pb.m_hasparameters());

		if(blob_pb.has_m_m1()) source.SetM1(blob_pb.m_m1());
		if(blob_pb.has_m_m2()) source.SetM2(blob_pb.m_m2());
		if(blob_pb.has_m_m3()) source.SetM3(blob_pb.m_m3());
		if(blob_pb.has_m_m4()) source.SetM4(blob_pb.m_m4());
		if(blob_pb.has_m_m1_curv()) source.SetM1Curv(blob_pb.m_m1_curv());
		if(blob_pb.has_m_m2_curv()) source.SetM2Curv(blob_pb.m_m2_curv());
		if(blob_pb.has_m_s()) source.SetS(blob_pb.m_s());
		if(blob_pb.has_m_smax()) source.SetSmax(blob_pb.m_smax());
		if(blob_pb.has_m_smin()) source.SetSmin(blob_pb.m_smin());
		if(blob_pb.has_m_sxx()) source.SetSxx(blob_pb.m_sxx());
		if(blob_pb.has_m_syy()) source.SetSyy(blob_pb.m_syy());
		if(blob_pb.has_m_sxy()) source.SetSxy(blob_pb.m_sxy());
		if(blob_pb.has_m_sx()) source.SetSx(blob_pb.m_sx());
		if(blob_pb.has_m_sy()) source.SetSy(blob_pb.m_sy());
		if(blob_pb.has_m_pixidmax()) source.SetSmaxPixId(blob_pb.m_pixidmax());
		if(blob_pb.has_m_pixidmin()) source.SetSminPixId(blob_pb.m_pixidmin());
		if(blob_pb.has_m_s_curv()) source.SetScurv(blob_pb.m_s_curv());
		if(blob_pb.has_m_s_edge()) source.SetSedge(blob_pb.m_s_edge());

		/*
		//============= MARKED FOR REMOVAL =================
		if(blob_pb.has_m_imagesizex() && blob_pb.has_m_imagesizey()) source.SetImageSize(blob_pb.m_imagesizex(),blob_pb.m_imagesizey());
		if(blob_pb.has_m_imageminx() && blob_pb.has_m_imagemaxx() &&
			 blob_pb.has_m_imageminy() && blob_pb.has_m_imagemaxy()
		){
			source.SetImageRange(blob_pb.m_imageminx(),blob_pb.m_imagemaxx(),blob_pb.m_imageminy(),blob_pb.m_imagemaxy());
		}
		if(blob_pb.has_m_imagemins() && blob_pb.has_m_imagemaxs()) {
			source.SetImageSRange(blob_pb.m_imagemins(),blob_pb.m_imagemaxs());
		}

		if(blob_pb.has_m_imageminscurv() && blob_pb.has_m_imagemaxscurv()) {
			source.SetImageScurvRange(blob_pb.m_imageminscurv(),blob_pb.m_imagemaxscurv());
		}

		if(blob_pb.has_m_imageminsedge() && blob_pb.has_m_imagemaxsedge()) {
			source.SetImageSedgeRange(blob_pb.m_imageminsedge(),blob_pb.m_imagemaxsedge());
		}
	
		if(blob_pb.has_m_imagerms()) source.SetImageRMS(blob_pb.m_imagerms());
		//=======================================================
		*/

		
		if( blob_pb.has_m_xmin() && blob_pb.has_m_xmax() &&
				blob_pb.has_m_ymin() && blob_pb.has_m_ymax()
		) {
			source.SetSourceRange(blob_pb.m_xmin(),blob_pb.m_xmax(),blob_pb.m_ymin(),blob_pb.m_ymax());
		}

		if( blob_pb.has_m_ix_min() && blob_pb.has_m_ix_max() &&
				blob_pb.has_m_iy_min() && blob_pb.has_m_iy_max()
		) {
			source.SetSourcePixelRange(blob_pb.m_ix_min(),blob_pb.m_ix_max(),blob_pb.m_iy_min(),blob_pb.m_iy_max());
		}

		//Add pixel collection to blob
		Pixel* aPixel= 0;
		std::vector<Pixel*> pixel_list;
		for(int i=0;i<blob_pb.m_pixels_size();i++){
			const CaesarPB::Pixel& thisPixelPB= blob_pb.m_pixels(i);
			aPixel= new Pixel;
			if(EncodeProtobufToPixel(*aPixel,thisPixelPB)<0){
				//Clear pixel list				
				for(unsigned int k=0;k<pixel_list.size();k++){
					delete pixel_list[k];
					pixel_list[k]= 0;	
				}
				pixel_list.clear();
	
				std::stringstream errMsg;
				errMsg<<"Pixel no. "<<i+1<<" encoding from protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
			pixel_list.push_back(aPixel);
			//source.AddPixel(aPixel);//NB: Do not use this otherwise NPix and all other moments will be incremented
		}	
		source.SetPixels(pixel_list);

		//Add contour collection to blob
		Contour* aContour= 0;
		std::vector<Contour*> contour_list;
		for(int i=0;i<blob_pb.m_contours_size();i++){
			const CaesarPB::Contour& thisContourPB= blob_pb.m_contours(i);
			aContour= new Contour;
			if(EncodeProtobufToContour(*aContour,thisContourPB)<0){
				//Clear contour list				
				for(unsigned int k=0;k<contour_list.size();k++){
					delete contour_list[k];
					contour_list[k]= 0;	
				}
				contour_list.clear();
	
				std::stringstream errMsg;
				errMsg<<"Contour no. "<<i+1<<" encoding from protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
			contour_list.push_back(aContour);
			source.AddContour(aContour);
		}	

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Blob encoding from protobuf failed with status "<<e.what());
		return -1;
	}


	return 0;

}//close EncodeProtobufToBlob()


int Serializer::EncodeProtobufToPixel(Pixel& pixel,const CaesarPB::Pixel& pixel_pb){

	try {
		//Set public source fields
		if(pixel_pb.has_id()) pixel.id= pixel_pb.id();

		if(pixel_pb.has_type()){
			if(pixel_pb.type()==CaesarPB::Pixel::eNormal) pixel.type= Pixel::eNormal;
			else if(pixel_pb.type()==CaesarPB::Pixel::eSeed) pixel.type= Pixel::eSeed;
			else if(pixel_pb.type()==CaesarPB::Pixel::eHalo) pixel.type= Pixel::eHalo;
			else{
				throw std::runtime_error("Invalid pixel type, failed to encode!");
			}
		}

		if(pixel_pb.has_s()) pixel.S= pixel_pb.s();
		if(pixel_pb.has_x()) pixel.x= pixel_pb.x();
		if(pixel_pb.has_y()) pixel.y= pixel_pb.y();
		if(pixel_pb.has_ix()) pixel.ix= pixel_pb.ix();
		if(pixel_pb.has_iy()) pixel.iy= pixel_pb.iy();
		//if(pixel_pb.has_isonedge()) pixel.isOnEdge= pixel_pb.isonedge();
		//if(pixel_pb.has_distancetoedge()) pixel.distanceToEdge= pixel_pb.distancetoedge();
		
		//Set private fields
		if(pixel_pb.has_s_curv()) pixel.SetCurv(pixel_pb.s_curv());
		if(pixel_pb.has_s_edge()) pixel.SetEdge(pixel_pb.s_edge());
		if(pixel_pb.has_bkglevel() && pixel_pb.has_noiselevel()){
			pixel.SetBkg(pixel_pb.bkglevel(),pixel_pb.noiselevel());
		}
	
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Pixel encoding to protobuf failed with status "<<e.what());
		return -1;
	}


	return 0;

}//close EncodeProtobufToPixel()


int Serializer::EncodeProtobufToContour(Contour& contour,const CaesarPB::Contour& contour_pb){

	try {		
		if(contour_pb.has_hasparameters()) contour.HasParameters= contour_pb.hasparameters();
		if(contour_pb.has_area()) contour.Area= contour_pb.area();
		if(contour_pb.has_perymeter()) contour.Perymeter= contour_pb.perymeter();
		if(contour_pb.has_isconvexcontour()) contour.IsConvexContour= contour_pb.isconvexcontour();
		if(contour_pb.has_circularityratio()) contour.CircularityRatio= contour_pb.circularityratio();

		
		if(contour_pb.has_boundingboxcenter()) {
			const CaesarPB::Point& BoundingBoxCenterPB= contour_pb.boundingboxcenter();
			if(EncodeProtobufToPoint(contour.BoundingBoxCenter,BoundingBoxCenterPB)<0){
				throw std::runtime_error("Failed to encode bounding box center field");	
			}
		}

		if(contour_pb.has_boundingboxmaj()) contour.BoundingBoxMaj= contour_pb.boundingboxmaj();
		if(contour_pb.has_boundingboxmin()) contour.BoundingBoxMin= contour_pb.boundingboxmin();
		if(contour_pb.has_boundingboxangle()) contour.BoundingBoxAngle= contour_pb.boundingboxangle();
		if(contour_pb.has_elongation()) contour.Elongation= contour_pb.elongation();
		if(contour_pb.has_rectangularity()) contour.Rectangularity= contour_pb.rectangularity();
		if(contour_pb.has_roundness()) contour.Roundness= contour_pb.roundness();
		if(contour_pb.has_eccentricity()) contour.Eccentricity= contour_pb.eccentricity();
		if(contour_pb.has_tiltangle()) contour.TiltAngle= contour_pb.tiltangle();
		if(contour_pb.has_hasellipsefit()) contour.HasEllipseFit= contour_pb.hasellipsefit();

		if(contour_pb.has_ellipsecenter()) {
			const CaesarPB::Point& EllipseCenterPB= contour_pb.ellipsecenter();
			if(EncodeProtobufToPoint(contour.EllipseCenter,EllipseCenterPB)<0){
				throw std::runtime_error("Failed to encode ellipse center field");	
			}
		}

		if(contour_pb.has_ellipsemajaxis()) contour.EllipseMajAxis= contour_pb.ellipsemajaxis();
		if(contour_pb.has_ellipseminaxis()) contour.EllipseMinAxis= contour_pb.ellipseminaxis();
		if(contour_pb.has_ellipserotangle()) contour.EllipseRotAngle= contour_pb.ellipserotangle();
		if(contour_pb.has_ellipsefitredchi2()) contour.EllipseFitRedChi2= contour_pb.ellipsefitredchi2();
		if(contour_pb.has_ellipsearearatio()) contour.EllipseAreaRatio= contour_pb.ellipsearearatio();
		
		if(contour_pb.has_m00()) contour.m00= contour_pb.m00();
		if(contour_pb.has_m10()) contour.m10= contour_pb.m10();
		if(contour_pb.has_m01()) contour.m01= contour_pb.m01();
		if(contour_pb.has_m20()) contour.m20= contour_pb.m20();
		if(contour_pb.has_m11()) contour.m11= contour_pb.m11();
		if(contour_pb.has_m02()) contour.m02= contour_pb.m02();
		if(contour_pb.has_m30()) contour.m30= contour_pb.m30();
		if(contour_pb.has_m21()) contour.m21= contour_pb.m21();
		if(contour_pb.has_m12()) contour.m12= contour_pb.m12();
		if(contour_pb.has_m03()) contour.m03= contour_pb.m03();
		
		if(contour_pb.has_mu20()) contour.m00= contour_pb.mu20();
		if(contour_pb.has_mu11()) contour.m00= contour_pb.mu11();
		if(contour_pb.has_mu02()) contour.m00= contour_pb.mu02();
		if(contour_pb.has_mu30()) contour.m00= contour_pb.mu30();
		if(contour_pb.has_mu21()) contour.m00= contour_pb.mu21();
		if(contour_pb.has_mu12()) contour.m00= contour_pb.mu12();
		if(contour_pb.has_mu03()) contour.m00= contour_pb.mu03();

		if(contour_pb.has_nu20()) contour.m00= contour_pb.nu20();
		if(contour_pb.has_nu11()) contour.m00= contour_pb.nu11();
		if(contour_pb.has_nu02()) contour.m00= contour_pb.nu02();
		if(contour_pb.has_nu30()) contour.m00= contour_pb.nu30();
		if(contour_pb.has_nu21()) contour.m00= contour_pb.nu21();
		if(contour_pb.has_nu12()) contour.m00= contour_pb.nu12();
		if(contour_pb.has_nu03()) contour.m00= contour_pb.nu03();
		
	
		for(int i=0;i<contour_pb.humoments_size();i++){		
			(contour.HuMoments)[i]= contour_pb.humoments(i);
		}
		
		for(int i=0;i<contour_pb.boundingboxvertex_size();i++){		
			const CaesarPB::Point& thisBoundingBoxVertexPB= contour_pb.boundingboxvertex(i);
			if(EncodeProtobufToPoint((contour.BoundingBoxVertex)[i],thisBoundingBoxVertexPB)<0){
				throw std::runtime_error("Failed to encode bounding box vertex field");	
			}
		}

		if(contour_pb.has_centroid()) {
			const CaesarPB::Point& CentroidPB= contour_pb.centroid();
			if(EncodeProtobufToPoint(contour.Centroid,CentroidPB)<0){
				throw std::runtime_error("Failed to encode centroid field");	
			}
		}
	
		if(contour_pb.has_hasfdpars()) contour.HasFDPars= contour_pb.hasfdpars(); 
		for(int i=0;i<contour_pb.realfds_size();i++){		
			(contour.RealFDs).push_back(contour_pb.realfds(i));
		}
		for(int i=0;i<contour_pb.imagfds_size();i++){		
			(contour.ImagFDs).push_back(contour_pb.imagfds(i));
		}
		for(int i=0;i<contour_pb.modfds_size();i++){		
			(contour.ModFDs).push_back(contour_pb.modfds(i));
		}
	
		if(contour_pb.has_hasbepars()) contour.HasBEPars= contour_pb.hasbepars(); 
		for(int i=0;i<contour_pb.bendingenergies_size();i++){		
			(contour.BendingEnergies).push_back(contour_pb.bendingenergies(i));
		}

		if(contour_pb.has_hascentroiddistancefdpars()) contour.HasCentroidDistanceFDPars= contour_pb.hascentroiddistancefdpars(); 
		for(int i=0;i<contour_pb.centroiddistancemodfds_size();i++){		
			(contour.CentroidDistanceModFDs).push_back(contour_pb.centroiddistancemodfds(i));
		}

		//Add points
		for(int i=0;i<contour_pb.m_points_size();i++){		
			const CaesarPB::Point& thisContourPointPB= contour_pb.m_points(i);
			TVector2 thisContourPoint;
			if(EncodeProtobufToPoint(thisContourPoint,thisContourPointPB)<0){
				throw std::runtime_error("Failed to encode centroid field");	
			}
			contour.AddPoint(thisContourPoint);
		}	

	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Contour encoding from protobuf failed with status "<<e.what());
		return -1;
	}
	
	return 0;

}//close EncodeProtobufToContour()

int Serializer::EncodeProtobufToPoint(TVector2& point,const CaesarPB::Point& point_pb){

	try {
		if(point_pb.has_x()) point.SetX(point_pb.x());
		if(point_pb.has_y()) point.SetY(point_pb.y());
	}
	catch(std::exception const & e) {
		ERROR_LOG("Point encoding from protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeProtobufToPoint()

int Serializer::EncodeProtobufToTaskData(TaskData& taskData,const CaesarPB::TaskData& taskData_pb){

	try {		
		//Fill task info
		/*
		//--> filename
		if(!taskData_pb.has_filename()){
			ERROR_LOG("Missing filename field, failed to encode!");
			return -1;
		}
		taskData.filename= taskData_pb.filename();
		if(taskData.filename==""){
			ERROR_LOG("Empty string filename field, failed to encode!");
			return -1;
		}

		//--> jobId
		if(!taskData_pb.has_jobid()){
			ERROR_LOG("Missing filename field, failed to encode!");
			return -1;
		}
		taskData.jobId= taskData_pb.jobid();
		if(taskData.jobId==""){
			ERROR_LOG("Empty string jobId field, failed to encode!");
			return -1;
		}
		*/

		//--> workerId
		if(!taskData_pb.has_workerid()){
			ERROR_LOG("Missing workerId field, failed to encode!");
			return -1;
		}
		taskData.workerId= taskData_pb.workerid();

		//--> ix_min
		if(!taskData_pb.has_ix_min()){
			ERROR_LOG("Missing ix_min field, failed to encode!");
			return -1;
		}
		taskData.ix_min= taskData_pb.ix_min();
		if(taskData.ix_min<0 && taskData.ix_min!=-1){
			ERROR_LOG("Invalid ix_min field, failed to encode!");
			return -1;
		}

		//--> ix_max
		if(!taskData_pb.has_ix_max()){
			ERROR_LOG("Missing ix_max field, failed to encode!");
			return -1;
		}
		taskData.ix_max= taskData_pb.ix_max();
		if(taskData.ix_max<0 && taskData.ix_max!=-1){
			ERROR_LOG("Invalid ix_max field, failed to encode!");
			return -1;
		}

		//--> iy_min
		if(!taskData_pb.has_iy_min()){
			ERROR_LOG("Missing iy_min field, failed to encode!");
			return -1;
		}
		taskData.iy_min= taskData_pb.iy_min();
		if(taskData.iy_min<0 && taskData.iy_min!=-1){
			ERROR_LOG("Invalid iy_min field, failed to encode!");
			return -1;
		}

		//--> iy_max
		if(!taskData_pb.has_iy_max()){
			ERROR_LOG("Missing iy_max field, failed to encode!");
			return -1;
		}
		taskData.iy_max= taskData_pb.iy_max();
		if(taskData.iy_max<0 && taskData.iy_max!=-1){
			ERROR_LOG("Invalid iy_max field, failed to encode!");
			return -1;
		}

		/*
		//--> x_min
		if(!taskData_pb.has_x_min()){
			ERROR_LOG("Missing x_min field, failed to encode!");
			return -1;
		}
		taskData.x_min= taskData_pb.x_min();
		
		//--> x_max
		if(!taskData_pb.has_x_max()){
			ERROR_LOG("Missing x_max field, failed to encode!");
			return -1;
		}
		taskData.x_max= taskData_pb.x_max();
		
		//--> y_min
		if(!taskData_pb.has_y_min()){
			ERROR_LOG("Missing y_min field, failed to encode!");
			return -1;
		}
		taskData.y_min= taskData_pb.y_min();

		//--> y_max
		if(!taskData_pb.has_y_max()){
			ERROR_LOG("Missing y_max field, failed to encode!");
			return -1;
		}
		taskData.y_max= taskData_pb.y_max();
		*/

	
		//Encode neighbour list
		for(int i=0;i<taskData_pb.neighbortaskid_size();i++){	
			(taskData.neighborTaskId).push_back(taskData_pb.neighbortaskid(i));
		}//end loop

		for(int i=0;i<taskData_pb.neighborworkerid_size();i++){	
			(taskData.neighborTaskId).push_back(taskData_pb.neighborworkerid(i));
		}//end loop

		//Parse source collection
		Source* aSource= 0;
		for(int i=0;i<taskData_pb.sources_size();i++){		
			const CaesarPB::Source& thisCaesarPB= taskData_pb.sources(i);
			aSource= new Source;
			if(EncodeProtobufToSource(*aSource,thisCaesarPB)<0){
				delete aSource;
				aSource= 0;
				throw std::runtime_error("Failed to encode source in collection from protobuf!");	
			}
			(taskData.sources).push_back(aSource);
		}	


		Source* aSource_edge= 0;
		for(int i=0;i<taskData_pb.sources_edge_size();i++){		
			const CaesarPB::Source& thisCaesarPB= taskData_pb.sources_edge(i);
			aSource_edge= new Source;
			if(EncodeProtobufToSource(*aSource_edge,thisCaesarPB)<0){
				delete aSource_edge;
				aSource_edge= 0;
				throw std::runtime_error("Failed to encode edge source in collection from protobuf!");	
			}
			(taskData.sources_edge).push_back(aSource_edge);
		}	


		Source* aExtSource= 0;
		for(int i=0;i<taskData_pb.ext_sources_size();i++){		
			const CaesarPB::Source& thisCaesarPB= taskData_pb.ext_sources(i);
			aExtSource= new Source;
			if(EncodeProtobufToSource(*aExtSource,thisCaesarPB)<0){
				delete aExtSource;
				aExtSource= 0;
				throw std::runtime_error("Failed to encode source in collection from protobuf!");	
			}
			(taskData.ext_sources).push_back(aExtSource);
		}	

		Source* aExtSource_edge= 0;
		for(int i=0;i<taskData_pb.ext_sources_edge_size();i++){		
			const CaesarPB::Source& thisCaesarPB= taskData_pb.ext_sources_edge(i);
			aExtSource_edge= new Source;
			if(EncodeProtobufToSource(*aExtSource_edge,thisCaesarPB)<0){
				delete aExtSource_edge;
				aExtSource_edge= 0;
				throw std::runtime_error("Failed to encode source in collection from protobuf!");	
			}
			(taskData.ext_sources_edge).push_back(aExtSource_edge);
		}	
		
	}//close try block
	catch(std::exception const & e) {
		//Clear allocated sources
		for(unsigned int i=0;i<(taskData.sources).size();i++){
			delete (taskData.sources)[i];
			(taskData.sources)[i]= 0;
		}
		taskData.sources.clear();

		for(unsigned int i=0;i<(taskData.sources_edge).size();i++){
			delete (taskData.sources_edge)[i];
			(taskData.sources_edge)[i]= 0;
		}
		taskData.sources_edge.clear();

		for(unsigned int i=0;i<(taskData.ext_sources).size();i++){
			delete (taskData.ext_sources)[i];
			(taskData.ext_sources)[i]= 0;
		}
		taskData.ext_sources.clear();

		for(unsigned int i=0;i<(taskData.ext_sources_edge).size();i++){
			delete (taskData.ext_sources_edge)[i];
			(taskData.ext_sources_edge)[i]= 0;
		}
		taskData.ext_sources_edge.clear();

		ERROR_LOG("Task data encoding from protobuf failed with status "<<e.what());
		return -1;
	}
	
	return 0;

}//close EncodeProtobufToWorkerData()

int Serializer::EncodeProtobufToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,const CaesarPB::TaskDataCollection& taskDataCollection_pb,bool isTaskCollectionPreAllocated){

	
	if(!isTaskCollectionPreAllocated){
		try {		
			//First clear existing vector
			for(unsigned int i=0;i<taskDataCollection.size();i++){
				if(taskDataCollection[i]){
					delete taskDataCollection[i];
					taskDataCollection[i]= 0;
				}
			}//end loop tasks
			taskDataCollection.clear();

			//Now fill vector with new tasks
			TaskData* aTaskData= 0;
			for(int i=0;i<taskDataCollection_pb.tasks_size();i++){		
				const CaesarPB::TaskData& thisTaskPB= taskDataCollection_pb.tasks(i);
				aTaskData= new TaskData;
				if(EncodeProtobufToTaskData(*aTaskData,thisTaskPB)<0){
					delete aTaskData;
					aTaskData= 0;
					throw std::runtime_error("Failed to encode task data in collection from protobuf!");	
				}
				taskDataCollection.push_back(aTaskData);
			}//end loop tasks
			
		}//close try block

		catch(std::exception const & e) {
			//Clear allocated tasks
			for(unsigned int i=0;i<taskDataCollection.size();i++){
				if(taskDataCollection[i]){
					delete taskDataCollection[i];
					taskDataCollection[i]= 0;
				}
			}//end loop tasks
			taskDataCollection.clear();
			ERROR_LOG("Task data collection encoding from protobuf failed with status "<<e.what());
			return -1;
		}
	}//close if
	else {	
		try{
			//Check first if number of pre-allocated tasks and tasks to be set is equal
			int nTasks= taskDataCollection_pb.tasks_size();
			int nAllocatedTasks= (int)taskDataCollection.size();
			if(nTasks!=nAllocatedTasks){
				throw std::runtime_error("Pre-Allocated tasks in vector is different from the number of tasks to be set!");	
			}

			//Now fill vector with new tasks
			for(int i=0;i<taskDataCollection_pb.tasks_size();i++){		
				const CaesarPB::TaskData& thisTaskPB= taskDataCollection_pb.tasks(i);
				if(!taskDataCollection[i]){
					throw std::runtime_error("Null pointer to task item in collection!");
				}

				if(EncodeProtobufToTaskData(*taskDataCollection[i],thisTaskPB)<0){//this update the collection
					throw std::runtime_error("Failed to encode task data in collection from protobuf!");	
				}
				
			}//end loop tasks	
		}//close try block
		catch(std::exception const & e) {
			//Do not clear allocated tasks in this case (because they have been previously allocated)
			//NB: Some of the tasks (before the crash) have been updated in the vector
			ERROR_LOG("Task data collection encoding from protobuf failed with status "<<e.what());
			return -1;
		}
	}//close else

	return 0;

}//close EncodeProtobufToTaskDataCollection()


int Serializer::BufferToSource(Source& source,SBuffer& buffer){

	//## Check for empty data
	if((buffer.data).empty() || buffer.size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		CaesarPB::Source source_pb;
  	if( !source_pb.ParseFromString(buffer.data) ) {
			throw std::runtime_error("Parsing to protobuf failed!");
		}

		//## Convert protobuf to Source
		if( EncodeProtobufToSource(source,source_pb)<0) {
			throw std::runtime_error("Encoding from protobuf to source failed!");
		}
		
	}//close try
	catch(std::exception const & e) {
		ERROR_LOG("Parsing source from buffer failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close BufferToSource()



int Serializer::BufferToTaskData(TaskData& taskData,SBuffer& buffer){

	//## Check for empty data
	if((buffer.data).empty() || buffer.size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		CaesarPB::TaskData taskData_pb;
  	if( !taskData_pb.ParseFromArray(buffer.data.c_str(),buffer.size) ) {
			throw std::runtime_error("Parsing of buffer to TaskData protobuf failed!");
		}

		//## Convert protobuf to TaskData
		if( EncodeProtobufToTaskData(taskData,taskData_pb)<0) {
			throw std::runtime_error("Encoding from protobuf to taskData failed!");
		}
		
	}//close try
	catch(std::exception const & e) {
		ERROR_LOG("Parsing taskData from buffer failed (err="<<e.what()<<")");
		return -1;
	}

	return 0;

}//close BufferToTaskData()


int Serializer::CharArrayToTaskData(TaskData& taskData,char* buffer,long int buffer_size){

	//## Check for bad buffer
	if(!buffer || buffer_size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		CaesarPB::TaskData taskData_pb;
  	if( !taskData_pb.ParseFromArray(buffer,buffer_size) ) {
			throw std::runtime_error("Parsing of char array to TaskData protobuf failed!");
		}

		//## Convert protobuf to TaskData
		if( EncodeProtobufToTaskData(taskData,taskData_pb)<0) {
			throw std::runtime_error("Encoding from protobuf to taskData failed!");
		}
		
	}//close try
	catch(std::exception const & e) {
		ERROR_LOG("Parsing taskData from buffer failed (err="<<e.what()<<")");
		return -1;
	}

	return 0;

}//close CharArrayToTaskData()


int Serializer::BufferToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,SBuffer& buffer,bool isTaskCollectionPreAllocated){

	//## Check for empty data
	if((buffer.data).empty() || buffer.size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		CaesarPB::TaskDataCollection taskDataCollection_pb;
  	if( !taskDataCollection_pb.ParseFromArray(buffer.data.c_str(),buffer.size) ) {
			throw std::runtime_error("Parsing of buffer to TaskDataCollection protobuf failed!");
		}

		//## Convert protobuf to TaskDataCollection
		if( EncodeProtobufToTaskDataCollection(taskDataCollection,taskDataCollection_pb,isTaskCollectionPreAllocated)<0) {
			throw std::runtime_error("Encoding from protobuf to taskDataCollection failed!");
		}
		
	}//close try
	catch(std::exception const & e) {
		ERROR_LOG("Parsing taskDataCollection from buffer failed (err="<<e.what()<<")");
		return -1;
	}

	return 0;

}//close BufferToTaskDataCollection()

int Serializer::CharArrayToTaskDataCollection(std::vector<TaskData*>& taskDataCollection,char* buffer,long int buffer_size,bool isTaskCollectionPreAllocated){

	//## Check for empty data
	if(!buffer || buffer_size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		CaesarPB::TaskDataCollection taskDataCollection_pb;
  	if( !taskDataCollection_pb.ParseFromArray(buffer,buffer_size) ) {
			throw std::runtime_error("Parsing of char array to TaskDataCollection protobuf failed!");
		}

		//## Convert protobuf to TaskDataCollection
		if( EncodeProtobufToTaskDataCollection(taskDataCollection,taskDataCollection_pb,isTaskCollectionPreAllocated)<0) {
			throw std::runtime_error("Encoding from protobuf to taskDataCollection failed!");
		}
		
	}//close try
	catch(std::exception const & e) {
		ERROR_LOG("Parsing taskDataCollection from char array failed (err="<<e.what()<<")");
		return -1;
	}
	return 0;

}//close CharArrayToTaskDataCollection()



}//close namespace
