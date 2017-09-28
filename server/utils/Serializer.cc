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
#include <Logger.h>
#include <CodeUtils.h>


#include <TVector2.h>

//# MSG PACK
//#ifdef BUILD_CAESAR_SERVER
	#include <tango.h>
	#include <msgpack.hpp>
//#endif

#include <Source.pb.h>

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
ClassImp(Caesar::Serializer)

namespace Caesar {

Serializer::Serializer() {	
	
}   


Serializer::~Serializer(){

}

//#ifdef BUILD_CAESAR_SERVER
int Serializer::SourceToTangoPipe(Tango::DevicePipeBlob& pipe_blob,Source* source){

	//Check source
	if(!source) return -1;
	
	//Encode
	try {
		//## Create google protobuf source message
		SourcePB::Source source_pb;
		if(EncodeSourceToProtobuf(source_pb,source)<0){
			throw std::runtime_error("Encoding to protobuf failed!");
		}
		Tango::DevLong buffer_size= static_cast<Tango::DevLong>(source_pb.ByteSize());

		//## Encode source buffer to pipe
		pipe_blob.set_name(source->Name.c_str());
		
		std::vector<std::string> field_names {"size","data"};
		pipe_blob.set_data_elt_names(field_names);
		//pipe_blob["data"] << source_pb.SerializeAsString();
		//pipe_blob["size"] << buffer_size;

		pipe_blob<< buffer_size << source_pb.SerializeAsString();

	}//close try block
	catch(Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Exception occurred while creating and filling source pipe!");
		return -1;
	}
	catch(std::exception &e){
		ERROR_LOG("Runtime exception occurred (e="<<e.what()<<") while creating and filling source pipe!");
		return -1;
	}
	catch(...){
		ERROR_LOG("Unknown exception occurred while creating and filling source pipe (name: "<<source->Name<<"!");
		return -1;
	}

	return 0;

}//close SourceToTangoPipe()


int Serializer::SourceCollectionToTangoPipe(Tango::DevicePipeBlob& pipe_blob,std::vector<Source*>& source_list){
	
	//Check source list
	if(source_list.empty()) return 0;

	//Encode
	try {
		//Make source collection blob
		INFO_LOG("Init source blob collection with "<<source_list.size()<<" elements ...");
		pipe_blob.set_name("Sources");
		std::vector<std::string> field_names;
		for(unsigned int i=0;i<source_list.size();i++){
			std::string source_name= source_list[i]->Name;
			field_names.push_back(source_name);
		}
		pipe_blob.set_data_elt_names(field_names);
		
		//Filling source collection blob
		INFO_LOG("Filling source blob collection with "<<source_list.size()<<" elements ...");
		
		for(unsigned int i=0;i<source_list.size();i++){		
			std::string source_name= source_list[i]->Name;
			std::shared_ptr<Tango::DevicePipeBlob> aSourceBlobItem= std::make_shared<Tango::DevicePipeBlob>();
			if(SourceToTangoPipe(*aSourceBlobItem,source_list[i])<0){
				std::stringstream errMsg;
				errMsg<<"Failed to encode source pipe blob no. "<<i<<" in collection!";
				throw std::runtime_error(errMsg.str().c_str());
			}
			pipe_blob << *aSourceBlobItem;
		}//end loop sources

	}//close try block
	catch(const Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Exception occurred while creating and filling source collection pipe!");
		return -1;
	}
	catch(std::exception &e){
		ERROR_LOG("Runtime exception occurred (e="<<e.what()<<") while creating and filling source collection pipe!");
		return -1;
	}
	catch(...){
		ERROR_LOG("Unknown exception occurred while creating and filling source collection pipe!");
		return -1;
	}

	return 0;

}//close SourceCollectionToTangoPipe()


int Serializer::MakeSourceDataPipe(Tango::DevicePipeBlob& pipe_blob,std::string runId,long int IdX,long int IdY,std::vector<Source*>& sources,std::vector<Source*>& edge_sources){

	//Encode source collection to pipe blob
	std::shared_ptr<Tango::DevicePipeBlob> source_blob= std::make_shared<Tango::DevicePipeBlob>();
	if(SourceCollectionToTangoPipe(*source_blob,sources)<0){
		ERROR_LOG("Failed to encode source collection to blob!");
		return -1;
	}
	
	//Encode source collection to pipe blob
	std::shared_ptr<Tango::DevicePipeBlob> edge_source_blob= std::make_shared<Tango::DevicePipeBlob>();
	if(SourceCollectionToTangoPipe(*edge_source_blob,edge_sources)<0){
		ERROR_LOG("Failed to encode edge source collection to blob!");
		return -1;
	}

	//Make final blob
	pipe_blob.set_name("SourceData");
	std::vector<std::string> field_names {"runId","IdX","IdY","sources","edgeSources"};
	pipe_blob.set_data_elt_names(field_names);
	try {
		pipe_blob << runId << IdX << IdY << *source_blob << *edge_source_blob;
	}
	catch(const Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Exception occurred while inserting pipe blob fields!");
		return -1;
	}
	catch(...){
		ERROR_LOG("Unknown exception occurred while inserting pipe blob fields!");
		return -1;
	}

	return 0;

}//close MakeSourceDataPipe()


int Serializer::TangoPipeToSource(Source& source,Tango::DevicePipeBlob& pipe_blob){

	//Extract pipe in buffer
	SBuffer buffer;
	try {
		//pipe_blob >> buffer.size >> buffer.data;
		pipe_blob["size"] >> buffer.size;
		pipe_blob["data"] >> buffer.data;
	}
	catch(const Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Exception occurred while extracting pipe blob fields!");
		return -1;
	}
	catch(std::exception &e){
		ERROR_LOG("Runtime exception occurred (e="<<e.what()<<") while extracting pipe blob fields!");
		return -1;
	}
	catch(...){
		ERROR_LOG("Unknown exception occurred while extracting pipe blob fields!");
		return -1;
	}

	//Encode buffer to source
	if(BufferToSource(source,buffer)<0){
		ERROR_LOG("Failed to encode buffer to source!");
		return -1;
	}

	return 0;

}//close TangoPipeToSource()

int Serializer::TangoPipeToSourceCollection(std::vector<Source*>& source_list,Tango::DevicePipeBlob& pipe_blob){
	
	//Init
	source_list.clear();

	Source* aSource= 0;
	try{
		//Get blob size
		size_t nb = pipe_blob.get_data_elt_nb();
	
		//Extract blob items into Source
		for(unsigned int i=0;i<nb;i++){
			Tango::DevicePipeBlob thisBlob;
			pipe_blob >> thisBlob;

			aSource= new Source;
			if(TangoPipeToSource(*aSource,thisBlob)<0){
				delete aSource;
				aSource= 0;
				std::stringstream errMsg;
				errMsg<<"Failed to encode pipe blob no. "<<i<<" to source in collection!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}//end loop items in collection

		//Add source to list
		source_list.push_back(aSource);
		
	}//close try block
	catch(const Tango::DevFailed& e){
		Tango::Except::print_exception(e);
		ERROR_LOG("Exception occurred while creating and filling source collection pipe!");
		for(unsigned int i=0;i<source_list.size();i++){
			delete source_list[i];
			source_list[i]= 0;
		}
		return -1;
	}
	catch(std::exception &e){
		ERROR_LOG("Runtime exception occurred (e="<<e.what()<<") while creating and filling source collection pipe!");
		for(unsigned int i=0;i<source_list.size();i++){
			delete source_list[i];
			source_list[i]= 0;
		}
		return -1;
	}
	catch(...){
		ERROR_LOG("Unknown exception occurred while creating and filling source collection pipe!");
		for(unsigned int i=0;i<source_list.size();i++){
			delete source_list[i];
			source_list[i]= 0;
		}
		return -1;
	}

	return 0;

}//close TangoPipeToSourceCollection()

int Serializer::EncodePointToProtobuf(SourcePB::Point& point_pb,TVector2& point){

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

int Serializer::EncodeContourToProtobuf(SourcePB::Contour& contour_pb,Contour* contour){
	
	if(!contour) return -1;

	try {
		contour_pb.set_hasparameters(contour->HasParameters);
		contour_pb.set_area(contour->Area);
		contour_pb.set_perymeter(contour->Perymeter);
		contour_pb.set_isconvexcontour(contour->IsConvexContour);
		contour_pb.set_circularityratio(contour->CircularityRatio);
			
		SourcePB::Point* BoundingBoxCenter= new SourcePB::Point;
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

		SourcePB::Point* EllipseCenter= new SourcePB::Point;
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
			SourcePB::Point* thisBoundingBoxVertex = contour_pb.add_boundingboxvertex();
			if(EncodePointToProtobuf(*thisBoundingBoxVertex,(contour->BoundingBoxVertex)[i])<0){
				std::stringstream errMsg;
				errMsg<<"BoundingBoxVertex no. "<<i+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	
		SourcePB::Point* Centroid= new SourcePB::Point;
		//SourcePB::Point Centroid;
		if(EncodePointToProtobuf(*Centroid,contour->Centroid)<0){
		//if(EncodePointToProtobuf(Centroid,contour->Centroid)<0){
			throw std::runtime_error("Failed to encode ellipse center field");
		}
		contour_pb.set_allocated_centroid(Centroid);

		for(unsigned int i=0;i<contour->RealFDs.size();i++){		
 			contour_pb.add_realfds(contour->RealFDs[i]);
		}
		for(unsigned int i=0;i<contour->ImagFDs.size();i++){		
 			contour_pb.add_imagfds(contour->ImagFDs[i]);
		}
		for(unsigned int i=0;i<contour->ModFDs.size();i++){		
 			contour_pb.add_modfds(contour->ModFDs[i]);
		}
		for(unsigned int i=0;i<contour->BendingEnergies.size();i++){		
 			contour_pb.add_bendingenergies(contour->BendingEnergies[i]);
		}	
		for(unsigned int i=0;i<contour->CentroidDistanceModFDs.size();i++){		
 			contour_pb.add_centroiddistancemodfds(contour->CentroidDistanceModFDs[i]);
		}

		for(int i=0;i<contour->GetN();i++){		
			SourcePB::Point* thisPoint = contour_pb.add_m_points();
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


int Serializer::EncodePixelToProtobuf(SourcePB::Pixel& pixel_pb,Pixel* pixel){
	
	if(!pixel) return -1;

	try {
		//Set public source fields
		pixel_pb.set_id(pixel->id);
		if(pixel->type==Pixel::eNormal){
			pixel_pb.set_type(SourcePB::Pixel::eNormal);
		}
		else if(pixel->type==Pixel::eSeed){
			pixel_pb.set_type(SourcePB::Pixel::eSeed);
		}
		else if(pixel->type==Pixel::eHalo){
			pixel_pb.set_type(SourcePB::Pixel::eHalo);
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
		pixel_pb.set_isonedge(pixel->isOnEdge);
		pixel_pb.set_distancetoedge(pixel->distanceToEdge);

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

int Serializer::EncodeBlobToProtobuf(SourcePB::Blob& blob_pb,Source* source){
	
	if(!source) return -1;

	try{
		blob_pb.set_haspixelsatedge(source->HasPixelsAtEdge);
		blob_pb.set_id(source->Id);
		blob_pb.set_name(source->Name);
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
	
		long int imgsizex, imgsizey;
		source->GetImageSize(imgsizex,imgsizey);
		blob_pb.set_m_imagesizex(imgsizex);
		blob_pb.set_m_imagesizey(imgsizey);

		double xmin, xmax, ymin, ymax;
		source->GetImageRange(xmin,xmax,ymin,ymax);
		blob_pb.set_m_imageminx(xmin);
		blob_pb.set_m_imagemaxx(xmax);
		blob_pb.set_m_imageminy(ymin);
		blob_pb.set_m_imagemaxy(ymax);

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

		double sxmin, sxmax, symin, symax;
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
			SourcePB::Pixel* thisPixel = blob_pb.add_m_pixels();
			if(EncodePixelToProtobuf(*thisPixel,source->GetPixel(k))<0){
				std::stringstream errMsg;
				errMsg<<"Pixel no. "<<k+1<<" encoding to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		//Add contout to blob
		for(unsigned int k=0;k<source->GetContours().size();k++){
			SourcePB::Contour* thisContour = blob_pb.add_m_contours();
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


int Serializer::EncodeSourceToProtobuf(SourcePB::Source& source_pb,Source* source){

	if(!source) return -1;

	try {
		
		source_pb.set_type(source->Type);
		source_pb.set_flag(source->Flag);
	
		//Set private source fields
		source_pb.set_m_beamfluxintegral(source->GetBeamFluxIntegral());
		source_pb.set_m_isgoodsource(source->IsGoodSource());
		source_pb.set_m_depthlevel(source->GetDepthLevel());
		source_pb.set_m_hasnestedsources(source->HasNestedSources());
		
		//Set blob field
		SourcePB::Blob* blob= new SourcePB::Blob;
		if(EncodeBlobToProtobuf(*blob,source)<0){
			throw std::runtime_error("Failed to encode blob!");
		}
		source_pb.set_allocated_blob(blob);

		//Add nested sources
		for(int i=0;i<source->GetNestedSourceNumber();i++){
			SourcePB::Source* thisNestedSource = source_pb.add_m_nestedsources();
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


int Serializer::EncodeWorkerDataToProtobuf(SourcePB::WorkerData& workerData_pb,WorkerData* workerData){

	//## Check input data
	if(!workerData) return -1;

	try {
		//Set timestamp
		SourcePB::Timestamp* timestamp= new SourcePB::Timestamp;
		timestamp->set_seconds( ((workerData->info).timestamp).tv_sec );
		timestamp->set_nanos( ((workerData->info).timestamp).tv_usec*1000 );

		//Set task
		SourcePB::WorkerTask* taskPB= new SourcePB::WorkerTask;
		if( (workerData->info).SerializeToProtobuf(*taskPB)<0 ){
			throw std::runtime_error("Failed to serialize task to protobuf!");
		}
		workerData_pb.set_allocated_task(taskPB);

		//Set source data type
		if(workerData->data_type==WorkerData::eUNKNOWN_DATA) workerData_pb.set_data_type(SourcePB::WorkerData::eUNKNOWN_DATA);
		else if(workerData->data_type==WorkerData::eCOMPACT_SOURCE_DATA) workerData_pb.set_data_type(SourcePB::WorkerData::eCOMPACT_SOURCE_DATA);
		else if(workerData->data_type==WorkerData::eEXTENDED_SOURCE_DATA) workerData_pb.set_data_type(SourcePB::WorkerData::eEXTENDED_SOURCE_DATA);
		else{
			throw std::runtime_error("Invalid source data type detected!");	
		}

		//Set source collections
		for(unsigned int i=0;i<(workerData->sources).size();i++){
			SourcePB::Source* thisSourcePB = workerData_pb.add_sources();
			if(EncodeSourceToProtobuf(*thisSourcePB,(workerData->sources)[i])<0){
				std::stringstream errMsg;
				errMsg<<"Encoding of source no. "<<i+1<<" in collection to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}

		for(unsigned int i=0;i<(workerData->edge_sources).size();i++){
			SourcePB::Source* thisSourcePB = workerData_pb.add_edge_sources();
			if(EncodeSourceToProtobuf(*thisSourcePB,(workerData->edge_sources)[i])<0){
				std::stringstream errMsg;
				errMsg<<"Encoding of edge source no. "<<i+1<<" in collection to protobuf failed!";
				throw std::runtime_error(errMsg.str().c_str());
			}
		}
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Worker data encoding to protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeWorkerDataToProtobuf()


int Serializer::SourceToBuffer(SBuffer& buffer,Source* source){

	//## Check input source
	if(!source) return -1;

	try {
		//## Create google protobuf source message
		SourcePB::Source source_pb;
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


int Serializer::WorkerDataToBuffer(SBuffer& buffer,WorkerData* workerData){

	//## Check input data
	if(!workerData) {
		return -1;
	}

	try {
		//## Create google protobuf source message
		SourcePB::WorkerData workerData_pb;
		if(EncodeWorkerDataToProtobuf(workerData_pb,workerData)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer.size = workerData_pb.ByteSize();
		buffer.data= workerData_pb.SerializeAsString();		
		
	}//close try block
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed with status "<<e.what());
		return -1;
	}
	return 0;

}//close WorkerDataToBuffer()


int Serializer::WorkerDataToCharArray(unsigned char* buffer,long int& buffer_size,WorkerData* workerData){

	//## Check input data
	if(!workerData) {
		return -1;
	}

	try {
		//## Create google protobuf source message
		SourcePB::WorkerData workerData_pb;
		if(EncodeWorkerDataToProtobuf(workerData_pb,workerData)<0){
			throw std::runtime_error("Encoding failed!");
		}
		
		//## Fill buffer 
		buffer_size = workerData_pb.ByteSize();
		if(!buffer) buffer = (unsigned char*)malloc(buffer_size);
		workerData_pb.SerializeToArray(buffer, buffer_size);
		
	}//close try blocks
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close WorkerDataToCharArray()

int Serializer::EncodeProtobufToSource(Source& source,const SourcePB::Source& source_pb){

	try {
		//Set public fields
		if(source_pb.has_type()) source.Type= source_pb.type();
		if(source_pb.has_flag()) source.Flag= source_pb.flag();
	
		//Set private source fields
		if(source_pb.has_m_beamfluxintegral()) source.SetBeamFluxIntegral(source_pb.m_beamfluxintegral());
		if(source_pb.has_m_isgoodsource()) source.SetGoodSourceFlag(source_pb.m_isgoodsource());
		if(source_pb.has_m_depthlevel()) source.SetDepthLevel(source_pb.m_depthlevel());
		if(source_pb.has_m_hasnestedsources()) source.SetHasNestedSources(source_pb.m_hasnestedsources());
		
		//Set blob fields
		if(source_pb.has_blob()){
			const SourcePB::Blob& blob_pb= source_pb.blob();
			if(EncodeProtobufToBlob(source,blob_pb)<0){
				ERROR_LOG("Blob encoding failed!");
				return -1;
			}
		}//close if has blob
	

		//Set nested sources (to be done after setting blob fields)
		Source* aNestedSource= 0;
		std::vector<Source*> nested_sources;
		for(int i=0;i<source_pb.m_nestedsources_size();i++){
			const SourcePB::Source& thisNestedSourcePB = source_pb.m_nestedsources(i);
			aNestedSource= new Source;
			if(EncodeProtobufToSource(*aNestedSource,thisNestedSourcePB)<0){
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
		return false;
	}


	return 0;

}//close EncodeProtobufToSource()



int Serializer::EncodeProtobufToBlob(Source& source,const SourcePB::Blob& blob_pb){
		
	try{

		//## Set public fields
		//Main params
		if(blob_pb.has_haspixelsatedge()) source.HasPixelsAtEdge= blob_pb.haspixelsatedge();
		if(blob_pb.has_id()) source.Id= blob_pb.id();
		if(blob_pb.has_name()) source.Name= blob_pb.name();

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
			const SourcePB::Pixel& thisPixelPB= blob_pb.m_pixels(i);
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
			source.AddPixel(aPixel);
		}	

		//Add contour collection to blob
		Contour* aContour= 0;
		std::vector<Contour*> contour_list;
		for(int i=0;i<blob_pb.m_contours_size();i++){
			const SourcePB::Contour& thisContourPB= blob_pb.m_contours(i);
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


int Serializer::EncodeProtobufToPixel(Pixel& pixel,const SourcePB::Pixel& pixel_pb){

	try {
		//Set public source fields
		if(pixel_pb.has_id()) pixel.id= pixel_pb.id();

		if(pixel_pb.has_type()){
			if(pixel_pb.type()==SourcePB::Pixel::eNormal) pixel.type= Pixel::eNormal;
			else if(pixel_pb.type()==SourcePB::Pixel::eSeed) pixel.type= Pixel::eSeed;
			else if(pixel_pb.type()==SourcePB::Pixel::eHalo) pixel.type= Pixel::eHalo;
			else{
				throw std::runtime_error("Invalid pixel type, failed to encode!");
			}
		}

		if(pixel_pb.has_s()) pixel.S= pixel_pb.s();
		if(pixel_pb.has_x()) pixel.x= pixel_pb.x();
		if(pixel_pb.has_y()) pixel.y= pixel_pb.y();
		if(pixel_pb.has_ix()) pixel.ix= pixel_pb.ix();
		if(pixel_pb.has_iy()) pixel.iy= pixel_pb.iy();
		if(pixel_pb.has_isonedge()) pixel.isOnEdge= pixel_pb.isonedge();
		if(pixel_pb.has_distancetoedge()) pixel.distanceToEdge= pixel_pb.distancetoedge();
		if(pixel_pb.has_s()) pixel.S= pixel_pb.s();
		if(pixel_pb.has_s()) pixel.S= pixel_pb.s();
		if(pixel_pb.has_s()) pixel.S= pixel_pb.s();

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


int Serializer::EncodeProtobufToContour(Contour& contour,const SourcePB::Contour& contour_pb){

	
	try {
		
		if(contour_pb.has_hasparameters()) contour.HasParameters= contour_pb.hasparameters();
		if(contour_pb.has_area()) contour.Area= contour_pb.area();
		if(contour_pb.has_perymeter()) contour.Perymeter= contour_pb.perymeter();
		if(contour_pb.has_isconvexcontour()) contour.IsConvexContour= contour_pb.isconvexcontour();
		if(contour_pb.has_circularityratio()) contour.CircularityRatio= contour_pb.circularityratio();

		
		if(contour_pb.has_boundingboxcenter()) {
			const SourcePB::Point& BoundingBoxCenterPB= contour_pb.boundingboxcenter();
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
			const SourcePB::Point& EllipseCenterPB= contour_pb.ellipsecenter();
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
			const SourcePB::Point& thisBoundingBoxVertexPB= contour_pb.boundingboxvertex(i);
			if(EncodeProtobufToPoint((contour.BoundingBoxVertex)[i],thisBoundingBoxVertexPB)<0){
				throw std::runtime_error("Failed to encode bounding box vertex field");	
			}
		}

		if(contour_pb.has_centroid()) {
			const SourcePB::Point& CentroidPB= contour_pb.centroid();
			if(EncodeProtobufToPoint(contour.Centroid,CentroidPB)<0){
				throw std::runtime_error("Failed to encode centroid field");	
			}
		}
	
		for(int i=0;i<contour_pb.realfds_size();i++){		
			(contour.RealFDs).push_back(contour_pb.realfds(i));
		}
		for(int i=0;i<contour_pb.imagfds_size();i++){		
			(contour.ImagFDs).push_back(contour_pb.imagfds(i));
		}
		for(int i=0;i<contour_pb.modfds_size();i++){		
			(contour.ModFDs).push_back(contour_pb.modfds(i));
		}
		for(int i=0;i<contour_pb.bendingenergies_size();i++){		
			(contour.BendingEnergies).push_back(contour_pb.bendingenergies(i));
		}
		for(int i=0;i<contour_pb.centroiddistancemodfds_size();i++){		
			(contour.CentroidDistanceModFDs).push_back(contour_pb.centroiddistancemodfds(i));
		}

		//Add points
		for(int i=0;i<contour_pb.m_points_size();i++){		
			const SourcePB::Point& thisContourPointPB= contour_pb.m_points(i);
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

int Serializer::EncodeProtobufToPoint(TVector2& point,const SourcePB::Point& point_pb){

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

int Serializer::EncodeProtobufToWorkerData(WorkerData& workerData,const SourcePB::WorkerData& workerData_pb){

	try {
		
		//Check task
		if(!workerData_pb.has_task()) {
			throw std::runtime_error("No task info present in the data!");	
		}

		//Fill task info
		const SourcePB::WorkerTask& taskPB= workerData_pb.task();
		if( (workerData.info).EncodeFromProtobuf(taskPB)<0 ){
			throw std::runtime_error("Failed to encode task info from protobuf!");	
		}

		//Set source data type
		if(workerData_pb.has_data_type()) {
			SourcePB::WorkerData_SourceDataType sourceDataTypePB= workerData_pb.data_type();
			if(sourceDataTypePB==SourcePB::WorkerData::eUNKNOWN_DATA) workerData.data_type= WorkerData::eUNKNOWN_DATA;
			else if(sourceDataTypePB==SourcePB::WorkerData::eCOMPACT_SOURCE_DATA) workerData.data_type= WorkerData::eCOMPACT_SOURCE_DATA;
			else if(sourceDataTypePB==SourcePB::WorkerData::eEXTENDED_SOURCE_DATA) workerData.data_type= WorkerData::eEXTENDED_SOURCE_DATA;
			else{
				throw std::runtime_error("Invalid source data type detected!");	
			}
		}//close if

		//Parse source collection
		Source* aSource= 0;
		for(int i=0;i<workerData_pb.sources_size();i++){		
			const SourcePB::Source& thisSourcePB= workerData_pb.sources(i);
			aSource= new Source;
			if(EncodeProtobufToSource(*aSource,thisSourcePB)<0){
				delete aSource;
				aSource= 0;
				throw std::runtime_error("Failed to encode source in collection from protobuf!");	
			}
			(workerData.sources).push_back(aSource);
		}	

		Source* aEdgeSource= 0;
		for(int i=0;i<workerData_pb.edge_sources_size();i++){		
			const SourcePB::Source& thisSourcePB= workerData_pb.edge_sources(i);
			aEdgeSource= new Source;
			if(EncodeProtobufToSource(*aEdgeSource,thisSourcePB)<0){
				delete aEdgeSource;
				aEdgeSource= 0;
				throw std::runtime_error("Failed to encode edge source in collection from protobuf!");	
			}
			(workerData.edge_sources).push_back(aEdgeSource);
		}	
	
		
	}//close try block
	catch(std::exception const & e) {
		//Clear allocated sources
		for(unsigned int i=0;i<(workerData.sources).size();i++){
			delete (workerData.sources)[i];
			(workerData.sources)[i]= 0;
		}
		workerData.sources.clear();

		for(unsigned int i=0;i<(workerData.edge_sources).size();i++){
			delete (workerData.edge_sources)[i];
			(workerData.edge_sources)[i]= 0;
		}
		workerData.edge_sources.clear();

		ERROR_LOG("Worker data encoding from protobuf failed with status "<<e.what());
		return -1;
	}

	return 0;

}//close EncodeProtobufToWorkerData()


int Serializer::BufferToSource(Source& source,SBuffer& buffer){

	//## Check for empty data
	if((buffer.data).empty() || buffer.size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		SourcePB::Source source_pb;
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



int Serializer::BufferToWorkerData(WorkerData& workerData,SBuffer& buffer){

	//## Check for empty data
	if((buffer.data).empty() || buffer.size<=0) {
		return -1;
	}

	try {
		//## Parse input and encode to protobuf message
		SourcePB::WorkerData workerData_pb;
  	//if( !workerData_pb.ParseFromString(buffer.data) ) {
		if( !workerData_pb.ParseFromArray(buffer.data.c_str(),buffer.size) ) {
			throw std::runtime_error("Parsing of buffer to WorkerData protobuf failed!");
		}

		//## Convert protobuf to WorkerData
		if( EncodeProtobufToWorkerData(workerData,workerData_pb)<0) {
			throw std::runtime_error("Encoding from protobuf to workerData failed!");
		}
		
	}//close try
	catch(std::exception const & e) {
		ERROR_LOG("Parsing workerData from buffer failed (err="<<e.what()<<")");
		return -1;
	}

	return 0;

}//close BufferToWorkerData()




int Serializer::WorkerTasksToJson(Json::Value& root,std::vector<WorkerTask*>& tasks){
			
	Json::Value taskArray;
	bool hasFailed= false;
	for (unsigned int i=0;i<tasks.size();i++){
		Json::Value jsonObj;
		if(tasks[i]->SerializeToJson(jsonObj)<0){
			hasFailed= true;
			break;
		}
		taskArray.append(jsonObj);
	}//end loop

	if(hasFailed){
		ERROR_LOG("Failed to serialize worker tasks to json object!");
		return -1;
	}

	root["tasks"]= taskArray;
			
	return 0;

}//close WorkerTasksToJson()

int Serializer::WorkerTasksToJsonString(std::string& jsonString,std::vector<WorkerTask*>& tasks,bool isMinified){

	//Encode first to json object
	Json::Value jsonObj;
	if(WorkerTasksToJson(jsonObj,tasks)<0) {
		ERROR_LOG("Failed to encode worker task collection to json object!");
		return -1;
	}

	//Encode to string
	if(CodeUtils::JsonToString(jsonString,jsonObj,isMinified)<0){
		ERROR_LOG("Failed to encode json object to string!");
		return -1;
	}

	return 0;
}//close WorkerTasksToJsonString()


int Serializer::JsonToWorkerTasks(std::vector<WorkerTask*>& tasks,Json::Value& root){

	//## Check object
	if(root.isNull() || root.empty()){
		ERROR_LOG("Invalid JSON (null or empty)!");
		return -1;
	}

	//## Check array
	const Json::Value& taskArray= root["tasks"];
	if(taskArray.isNull() || taskArray.empty() || !taskArray.isArray()){
		ERROR_LOG("Invalid json task object (null/empty or not an array)!");
		return -1;
	}

	//Loop array and fill values
	WorkerTask* aWorkerTask= 0;
	bool hasFailed= false;
	for (Json::ValueConstIterator it = taskArray.begin(); it != taskArray.end(); ++it) {
  	const Json::Value& taskItem = *it;

		aWorkerTask= new WorkerTask;
		if(aWorkerTask->EncodeFromJson(taskItem)<0){
			hasFailed= true;
			delete aWorkerTask;
			aWorkerTask= 0;
			break;
		}

		tasks.push_back(aWorkerTask);
	}//end loop array entries

	//Check if failed
	if(hasFailed){
		ERROR_LOG("Failed to encode from Json object!");
		for(unsigned int i=0;i<tasks.size();i++){
			if(tasks[i]){
				delete tasks[i];
				tasks[i]= 0;
			}
		}
		tasks.clear();
		return -1;
	}

	return 0;

}//close JsonToWorkerTasks()

int Serializer::JsonStringToWorkerTasks(std::vector<WorkerTask*>& tasks,std::string& jsonString){

	//Convert input string to json object
	Json::Reader reader;
	Json::Value root;
  if(!reader.parse(jsonString, root)) {
		ERROR_LOG("Failed to encode string to json ("<<reader.getFormattedErrorMessages()<<")");
		return -1;
	}

	//Convert from json object to WorkerTasks
	if(JsonToWorkerTasks(tasks,root)<0){
		ERROR_LOG("Failed to convert json object to worker tasks");
		return -1;
	}
	
	return 0;
}//close JsonStringToWorkerTasks()


int Serializer::SourceToBuffer(Source* source,msgpack::sbuffer& buffer){

	/*
	try{	
		msgpack::pack(&buffer,*source);
	}
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed (err: "<<e.what()<<")");
		return -1;
	}
	*/
	return 0;

}//close SourceToBuffer()

int Serializer::SourceToString(Source* source,std::string& msg){
	/*
	msgpack::sbuffer buffer;
	if(SourceToBuffer(source,buffer)<0) return -1;

	msg= std::string(buffer.data(),buffer.size());
	*/
	return 0;

}//close SourceToString()


int Serializer::SourceToDevString(Source* source,Tango::DevString& msg){

	/*
	msgpack::sbuffer buffer;
	if(SourceToBuffer(source,buffer)<0)	return -1;

	msg= CORBA::string_alloc(buffer.size());
	memcpy (msg, buffer.data(), buffer.size());	
	msg[buffer.size()]= 0;
	*/
	return 0;

}//close SourceToDevString()

int Serializer::SourceCollectionToBuffer(std::vector<Source*>& sources,msgpack::sbuffer& buffer){

	/*
	try{	
		msgpack::pack(&buffer,sources);
	}
	catch(std::exception const & e) {
		ERROR_LOG("Source encoding failed (err: "<<e.what()<<")");
		return -1;
	}
	*/
	return 0;

}//close SourceCollectionToBuffer()

int Serializer::SourceCollectionToString(std::vector<Source*>& sources,std::string& msg){

	/*
	msgpack::sbuffer buffer;
	if(SourceCollectionToBuffer(sources,buffer)<0) return -1;

	msg= std::string(buffer.data(),buffer.size());
	*/
	return 0;

}//close SourceCollectionToString()

int Serializer::SourceCollectionToDevString(std::vector<Source*>& sources,Tango::DevString& msg){

	/*
	msgpack::sbuffer buffer;
	if(SourceCollectionToBuffer(sources,buffer)<0)	return -1;

	msg= CORBA::string_alloc(buffer.size());
	memcpy (msg, buffer.data(), buffer.size());	
	msg[buffer.size()]= 0;
	*/
	return 0;
	
}//close SourceCollectionToDevString()

//#endif



}//close namespace
