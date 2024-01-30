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
* @file SourceImporter.cc
* @class SourceImporter
* @brief SourceImporter class
*
* Class to import image sources from different formats
* @author S. Riggi
* @date 17/05/2022
*/


#include <SourceImporter.h>
#include <Source.h>
#include <Contour.h>
#include <CodeUtils.h>
#include <AstroUtils.h>
#include <WCSUtils.h>
#include <Consts.h>
#include <ImgMetaData.h>

#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <TObject.h>
#include <TEllipse.h>

//WCSTOOLS (TO BE DEPRECATED)
//#include <wcs.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::SourceImporter)

namespace Caesar {

SourceImporter::SourceImporter()
{

}

SourceImporter::~SourceImporter()
{

}

//=================================================
//==        JSON IMPORTER
//=================================================
int SourceImporter::ImportFromJson(std::string filename, std::vector<Source*>& sources, Image* img)
{
	// - Init input collection
	sources.clear();

	// - Read json file
	Json::Value root;
	
	#ifdef LOGGING_ENABLED
		INFO_LOG("Reading and parsing json file " << filename <<" ...");
	#endif
	try{
		std::ifstream fin(filename);
		fin >> root;
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read and parse file "<< filename << "!");
		#endif
		return -1;
	}

	// - Read image
	if(!img){
		#ifdef LOGGING_ENABLED
			INFO_LOG("Reading image listed in json file " << filename <<" ...");
		#endif
		img= ReadImageFromJson(root);
	}

	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to input image or failed to read image from json file!");
		#endif
		return -1;
	}

	// - Read sources
	auto root_slist= root["sources"];
	if(root_slist.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source list is empty in parsed json file, returning empty source collection!");
		#endif
		return 0;
	}

	
	for (Json::Value::ArrayIndex i=0; i!=root_slist.size(); i++)
	{
		Source* source= ReadSourceFromJson(root_slist[i], img);
		if(!source){
			#ifdef LOGGING_ENABLED
				WARN_LOG("Failed to read source no. "<<i+1<<" from json, skip to next...");
			#endif
			continue;
		}

		sources.push_back(source);

	}//end loop sources

	return 0;

}//close ImportFromJson()


Image* SourceImporter::ReadImageFromJson(Json::Value& root)
{
	//- Get image path from json
	std::string filename= "";
	try {
		filename= root["metadata"]["inputdata"]["filename"].asString();
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to get image filename from json!");
		#endif
		return nullptr;
	}	

	//- Read image
	Image* img= new Image;
	if(img->ReadFITS(filename)<0)
	{
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read image from file "<< filename << "!");
		#endif
		delete img;
		img= 0;
		return nullptr;
	}

	//- Set beam info
	double bmaj= root["metadata"]["inputdata"]["bmaj"].asDouble();
	double bmin= root["metadata"]["inputdata"]["bmin"].asDouble();
	double pa= root["metadata"]["inputdata"]["pa"].asDouble();
  ImgMetaData* metadata= img->GetMetaData();
	if(metadata && !metadata->HasBeamInfo()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Input image has no beam information, setting user-supplied beam set in json ("<<bmaj<<","<<bmin<<","<<pa<<") ...");
		#endif
		metadata->SetBeamInfo(bmaj/3600., bmin/3600., pa);//convert bmaj/bmin from arcsec to degree
	}


	return img;

}//close ReadImageFromJson()

Source* SourceImporter::ReadSourceFromJson(Json::Value& json, Image* img)
{
	// - Retrieve some fields
	std::string sname= json["name"].asString();
	//int sindex= json["index"].asInt();

	// - Check number of islands
	unsigned int island_size= json["islands"].size();
	unsigned int nislands= json["nislands"].asUInt();
	if(island_size==0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("No islands for source "<<sname<<", this cannot occur, returning nullptr!");
		#endif
		return nullptr;
	}
	if(island_size!=nislands){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Mismatch between nislands field and island array size for source "<<sname<<", this cannot occur, returning nullptr!");
		#endif
		return nullptr;
	}

	// - Fill source data
	//   NB: If only one island is present and no parent, create standard source, otherwise create one source with nested 
	//Source* source= new Source(sname);
	Source* source= nullptr;
	auto root_ilist= json["islands"];

	if(nislands==1)
	{
		source= ReadSourceIslandFromJson(root_ilist[0], img);
	}
	else
	{
		//Create source
		source= new Source(sname);

		//Loop over islands
		std::vector<Source*> nested_sources;
		for (Json::Value::ArrayIndex i=0; i!=root_ilist.size(); i++)
		{
			Source* source_nested= ReadSourceIslandFromJson(root_ilist[i], img);
			if(!source_nested){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Failed to read island no. "<<i+1<<" of source "<<sname<<" from json, skip to next...");
				#endif
				continue;
			}
			nested_sources.push_back(source_nested);

			//Add nested source
			source->AddNestedSource(source_nested);

			//Merge nested source into mother
			bool copyPixels= true;
			bool checkIfAdjacent= false;
			bool computeStatPars= false;
			bool computeMorphPars= false;
			bool sumMatchingPixels= false;
			source->MergeSource(
				source_nested, 
				copyPixels, checkIfAdjacent,
				computeStatPars, computeMorphPars, 
				sumMatchingPixels
			);

		}//end loop islands

		//Set composite source flag
		source->SetAreNestedComponentsOfCompositeSource(true);

		//Recompute merged source stats
		bool computeRobustStats= true;
		bool forceRecomputing= true;
		source->ComputeStats(computeRobustStats,forceRecomputing);

		//Recompute merged source pars
		source->ComputeMorphologyParams();
			
	}//close else

	// - Set nest level
	source->SetDepthLevel(json["nest_level"].asInt());

	
	return source;

}//close ReadSourceFromJson()

Source* SourceImporter::ReadSourceIslandFromJson(Json::Value& json, Image* img)
{
	// - Return if no image is given
	if(!img){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Input image is null, cannot retrieve pixels and create source, returning nullptr ...");
		#endif
		return nullptr;
	}

	// - Retrieve image metadata if available
	ImgMetaData* metadata= img->GetMetaData();
	float Xmin= img->GetXmin();
	float Ymin= img->GetYmin();
	float Xmax= img->GetXmax();
	float Ymax= img->GetYmax();
	double beamArea= 0;
	double beam_bmaj= 0;
	double beam_bmin= 0;
	double beam_pa= 0;
	double pixSizeX= 0;
	double pixSizeY= 0;
	double pixSize= 0;

	if(metadata->HasBeamInfo()) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Beam info available, retrieving beam ellipse pars ...");
		#endif
		beamArea= metadata->GetBeamFluxIntegral();
		beam_bmaj= metadata->Bmaj*3600;//in arcsec
		beam_bmin= metadata->Bmin*3600;//in arcsec
		beam_pa= metadata->Bpa;//defined from north
		pixSizeX= metadata->dX*3600;//in arcsec 
		pixSizeY= metadata->dY*3600;//in arcsec 
		pixSize= std::min(fabs(pixSizeX),fabs(pixSizeY));
	}

	// - Check vertices and/or pixel list field are present
	std::string sname= json["name"].asString();
	auto root_vertices= json["vertices"];
	auto root_pixels= json["pixels"];
	bool has_pixels= !root_pixels.empty();
	bool has_contour= !root_vertices.empty();
	if(!has_pixels && !has_contour){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Source island "<<sname<<" has empty vertices & pixel fields, cannot create source, returning nullptr ...");
		#endif
		return nullptr;
	}

	// - Init source
	Source* source= new Source(sname);

	// - Retrieve pixels
	std::vector<Pixel*> pixels;
	Pixel* pixel= nullptr;
	bool hasBorderPixels= false;
	if(has_pixels)
	{
		for (Json::Value::ArrayIndex i=0; i!=root_pixels.size(); i++)
		{
			long int ix= root_pixels[i][0].asInt64();
			long int iy= root_pixels[i][1].asInt64();
			bool hasBin= img->HasBin(ix,iy);
			if(!hasBin){
				#ifdef LOGGING_ENABLED
					WARN_LOG("Pixel "<<i+1<<" of source island "<<sname<<" does not exist in image, skipping it ...");
				#endif
				continue;
			}
			double x= img->GetX(ix);
			double y= img->GetX(iy);
			long int gbin= img->GetBin(ix,iy);
			double S= img->GetPixelValue(ix,iy);
			if( img->IsEdgeBin(ix,iy) ) {
				source->SetEdgeFlag(true);
				hasBorderPixels= true;
			}

			pixel= new Pixel;
			pixel->S= S;
			pixel->id= gbin;
			pixel->SetPhysCoords(x,y);
			pixel->SetCoords(ix,iy);

			source->AddPixel(pixel);
		}
	}
	else
	{//retrieve pixels from vertices
		
		// - Create contour from vertex points
		Contour* contour= new Contour;
		for (Json::Value::ArrayIndex i=0; i!=root_vertices.size(); i++)
		{
			double x= root_vertices[i][0].asDouble();
			double y= root_vertices[i][1].asDouble();
			contour->AddPoint(TVector2(x,y));

		}//end loop vertex

		// - Compute contour pars
		if(contour->ComputeParameters()<0){	
			#ifdef LOGGING_ENABLED
				WARN_LOG("One/more failures occurred while computing contour parameter for source island "<<sname<<", returning nullptr!");
			#endif
			delete source;
			source= 0;
			return nullptr; 
		}
	
		// - Get contour bbox vertex points
		std::vector<TVector2> cvert= contour->BoundingBoxVertex;
		std::vector<TVector2> cvert_norot= contour->BoundingBoxVertex_noRot;
		#ifdef LOGGING_ENABLED
			INFO_LOG("Source island "<<sname<<": cvert=[("<<cvert[0].X()<<","<<cvert[0].Y()<<"),("<<cvert[1].X()<<","<<cvert[1].Y()<<"),("<<cvert[2].X()<<","<<cvert[2].Y()<<"),("<<cvert[3].X()<<","<<cvert[3].Y()<<")]");
			INFO_LOG("Source island "<<sname<<": cvert_norot=[("<<cvert_norot[0].X()<<","<<cvert_norot[0].Y()<<"),("<<cvert_norot[1].X()<<","<<cvert_norot[1].Y()<<"),("<<cvert_norot[2].X()<<","<<cvert_norot[2].Y()<<"),("<<cvert_norot[3].X()<<","<<cvert_norot[3].Y()<<")]");
		#endif

		std::vector<double> cvert_x;
		std::vector<double> cvert_y;
		for(size_t k=0;k<cvert_norot.size();k++){
			cvert_x.push_back(cvert_norot[k].X());
			cvert_y.push_back(cvert_norot[k].Y());
		}
		
		double xmin= *min_element(cvert_x.begin(), cvert_x.end());
		double xmax= *max_element(cvert_x.begin(), cvert_x.end());
		double ymin= *min_element(cvert_y.begin(), cvert_y.end());
		double ymax= *max_element(cvert_y.begin(), cvert_y.end());

		long int gbin_min= img->FindBin(xmin,ymin);
		long int gbin_max= img->FindBin(xmax,ymax);
		long int ix_min= img->GetBinX(gbin_min);
		long int ix_max= img->GetBinX(gbin_max);
		long int iy_min= img->GetBinY(gbin_min);
		long int iy_max= img->GetBinY(gbin_max);

		
		#ifdef LOGGING_ENABLED
			INFO_LOG("Source island "<<sname<<": cvert xmin/xmax="<<xmin<<"/"<<xmax<<", ymin/ymax="<<ymin<<"/"<<ymax<<", ix_min/ix_max="<<ix_min<<"/"<<ix_max<<", iy_min/iy_max="<<iy_min<<"/"<<iy_max);
		#endif

		// - Sort contour counter clock wise
		contour->SortPointsCounterClockWise();

		std::vector<TVector2> points= contour->GetPoints();
		for(size_t k=0;k<points.size();k++){
			#ifdef LOGGING_ENABLED
				INFO_LOG("Source island "<<sname<<": contour point("<<points[k].X()<<","<<points[k].Y()<<") ...");
			#endif
		}
		

		// - Loop over contour bbox pixels and take those inside contour
		bool includeBorders= true;
		bool sortContour= false;

		for(long int ix=ix_min;ix<=ix_max;ix++)
		{		
			double x= img->GetX(ix);

			for(long int iy=iy_min;iy<=iy_max;iy++)
			{
				double y= img->GetY(iy);
				//bool isInside= contour->IsPointInsideContourV2(x,y,includeBorders,sortContour);
				bool isInside= contour->IsPointInsideContour(x,y,includeBorders);
				if(!isInside) continue;
				#ifdef LOGGING_ENABLED
					INFO_LOG("Source island "<<sname<<": pixel("<<x<<","<<y<<") is inside contour ...");
				#endif

				long int gbin= img->GetBin(ix,iy);
				double S= img->GetPixelValue(ix,iy);
				if( img->IsEdgeBin(ix,iy) ) {
					source->SetEdgeFlag(true);
					hasBorderPixels= true;
				}

				pixel= new Pixel;
				pixel->S= S;
				pixel->id= gbin;
				pixel->SetPhysCoords(x,y);
				pixel->SetCoords(ix,iy);

				source->AddPixel(pixel);
			}//end loop bbox pixels x
		}//end loop bbox pixels y

	}//close else

	// - Set beam area
	source->SetBeamFluxIntegral(beamArea);
	
	// - Compute source stats
	source->ComputeStats();

	// - Compute source morph pars
	source->ComputeMorphologyParams();

	// - Adding image metadata to image (needed for WCS)
	source->SetImageMetaData(metadata, Xmin, Ymin);

	//- Set Bkg/noise estimators
	//  NB: pixel bkg & rms are not stored in json, so we need to reset them
	double nPixels= static_cast<double>(source->NPix);
	double bkgLevel= json["bkg"].asDouble();
	double bkgRMS= json["rms"].asDouble();
	double bkgSum= nPixels*bkgLevel;
	double bkgRMSSum= nPixels*bkgRMS;
	source->SetBkgSum(bkgSum);
	source->SetBkgRMSSum(bkgRMSSum);
	
	//- Source flags
	std::string morph_label= json["morph_label"].asString();
	int sourceMorphId= GetSourceMorphId(morph_label);
	source->MorphId= sourceMorphId;

	std::string sourceness_label= json["sourceness_label"].asString();
	int sourcenessId= GetSourcenessId(sourceness_label);
	source->SourcenessId= sourcenessId;

	//bool at_edge= json["border"].asBool();
	bool at_edge= static_cast<bool>(json["border"].asInt());
	source->SetEdgeFlag(at_edge);

	float sourcenessScore= json["sourceness_score"].asFloat();
	source->SourcenessScore= sourcenessScore;
	
	std::vector<std::string> user_tags;
	auto tags_list= json["tags"]; 
	for (Json::Value::ArrayIndex i=0; i!=tags_list.size(); i++) user_tags.push_back(tags_list[i].asString());
	source->SetTags(user_tags);

	// - Set classification info
	std::string objLabel= json["class_label"].asString();
	source->ObjClassId= GetAstroObjectType(objLabel);
	source->ObjClassStrId= objLabel;
	source->ObjClassScore= json["class_score"].asFloat();

	// - Set cross-match info
	//...
	//...

	// - Set spectral info
	auto json_spectralinfo= json["spectral_info"];
	if(!json_spectralinfo.empty())
	{
		SpectralIndexData sid;
		sid.spectralIndex= json_spectralinfo["alpha"].asDouble();
		sid.spectralIndexErr= json_spectralinfo["alpha_err"].asDouble();
		sid.isSpectralIndexFit= static_cast<bool>(json_spectralinfo["fit"].asInt());
		sid.spectralFitChi2= json_spectralinfo["chi2"].asDouble();
		sid.spectralFitNDF= json_spectralinfo["ndf"].asInt();

		source->SetSpectralIndexData(sid);
	}
	
	// - Set fit info
	auto json_fitinfo= json["fit_info"];
	source->SetHasFitInfo(false);

	if(!json_fitinfo.empty())
	{
		//Set fit status
		std::string fitQuality_str= json_fitinfo["fit_quality"].asString();
		int fitQuality= GetSourceFitQuality(fitQuality_str);
		int fitStatus= GetSourceFitStatusFromFitQuality(fitQuality);
		
		source->SetFitStatus(fitStatus);

		//Set fit pars
		SourceFitPars fitPars;
		int nPars= json_fitinfo["npars"].asInt();
		fitPars.SetNComponents(json_fitinfo["ncomponents"].asInt());
		fitPars.SetChi2(json_fitinfo["chi2"].asDouble());
		fitPars.SetNDF(json_fitinfo["ndf"].asInt());
		fitPars.SetStatus(fitStatus);
		fitPars.SetNPars(nPars);
		fitPars.SetNFreePars(json_fitinfo["npars_free"].asInt());
		fitPars.SetNFitPoints(json_fitinfo["ndata"].asInt());
		fitPars.SetFitQuality(fitQuality);

		double flux= json_fitinfo["flux"].asDouble();
		double flux_err= json_fitinfo["flux_err"].asDouble();
		double fluxDensity= flux;
		double fluxDensityErr= flux_err;
		if(beamArea>0){
			fluxDensity= flux*beamArea;
			fluxDensityErr= flux_err*beamArea;
		}

		fitPars.SetFluxDensity(fluxDensity);
		fitPars.SetFluxDensityErr(fluxDensityErr);

		if(metadata->HasBeamInfo()) {
			fitPars.SetComponentBeamEllipsePars(beam_bmaj,beam_bmin,beam_pa);
			fitPars.SetComponentImagePixSize(pixSize);
		}
	
		//Set components
		auto json_fitcomponents= json_fitinfo["components"];
		std::vector<SpectralIndexData> sids;

		if(!json_fitcomponents.empty())
		{
			int parCounter= 0;
			std::vector<std::string> parNames= {"A","x0","y0","sigmaX","sigmaY","theta"};

			for (Json::Value::ArrayIndex k=0; k!=json_fitcomponents.size(); k++) 
			{
				//Set par values
				std::vector<double> parVals= {
					json_fitcomponents[k]["Speak"].asDouble(),
					json_fitcomponents[k]["x"].asDouble(),
					json_fitcomponents[k]["y"].asDouble(),
					json_fitcomponents[k]["sx"].asDouble(),
					json_fitcomponents[k]["sy"].asDouble(),
					json_fitcomponents[k]["theta"].asDouble()
				};
				std::vector<double> parErrs= {
					json_fitcomponents[k]["Speak_err"].asDouble(),
					json_fitcomponents[k]["x_err"].asDouble(),
					json_fitcomponents[k]["y_err"].asDouble(),
					json_fitcomponents[k]["sx_err"].asDouble(),
					json_fitcomponents[k]["sy_err"].asDouble(),
					json_fitcomponents[k]["theta_err"].asDouble()
				};
	
				for(size_t i=0;i<parNames.size();i++) {
					fitPars.SetParValueAndError(k, parNames[i], parVals[i], parErrs[i]);
					parCounter++;
				}
			}//end loop components

			//Compute derivative matrix
			fitPars.ComputeFluxDensityDerivMatrix();

			//Get covariance matrix
			auto json_covmatrix= json_fitinfo["cov_matrix"];
			std::vector<double> covmatrix_triang_data;
			for (Json::Value::ArrayIndex i=0; i!=json_covmatrix.size(); i++){
				covmatrix_triang_data.push_back( json_covmatrix[i].asDouble() ); 
			} 
	
			TMatrixD CovMatrix(nPars,nPars);
			for(int i=0;i<nPars;i++){
				for(int j=i;j<nPars;j++){
					int index= CodeUtils::GetTriuIndexFrom2DIndex(i,j,nPars);
					#ifdef LOGGING_ENABLED
						DEBUG_LOG("Retrieving cov matrix array at index "<<index<<" ("<<i<<","<<j<<"), size="<<covmatrix_triang_data.size()<<" ...");
					#endif
					double c_ij= covmatrix_triang_data[index];
					CovMatrix(i,j)= c_ij;
					CovMatrix(j,i)= c_ij;
				}
			}

			fitPars.SetCovarianceMatrix(CovMatrix);

			//Compute other pars
			for (Json::Value::ArrayIndex k=0; k!=json_fitcomponents.size(); k++) 
			{
				//Compute and set ellipse pars
				fitPars.ComputeComponentEllipsePars(k);

				//Set WCS ellipse pars
				double x0_wcs= json_fitcomponents[k]["ra"].asDouble();
				double y0_wcs= json_fitcomponents[k]["dec"].asDouble();
				double bmaj_wcs= json_fitcomponents[k]["bmaj"].asDouble();
				double bmin_wcs= json_fitcomponents[k]["bmin"].asDouble();
				double pa_wcs= json_fitcomponents[k]["pa"].asDouble();
				fitPars.SetComponentFitWCSEllipsePars(k,x0_wcs,y0_wcs,bmaj_wcs,bmin_wcs,pa_wcs);

				//Set WCS ellipse par errors
				double x0_wcs_err= json_fitcomponents[k]["ra_err"].asDouble();
				double y0_wcs_err= json_fitcomponents[k]["dec_err"].asDouble();
				double bmaj_wcs_err= json_fitcomponents[k]["bmaj_err"].asDouble();
				double bmin_wcs_err= json_fitcomponents[k]["bmin_err"].asDouble();
				double pa_wcs_err= json_fitcomponents[k]["pa_err"].asDouble();
				fitPars.SetComponentFitWCSEllipseParErrors(k,x0_wcs_err,y0_wcs_err,bmaj_wcs_err,bmin_wcs_err,pa_wcs_err);

				//Set WCS deconv pars
				double bmaj_wcs_deconv= json_fitcomponents[k]["bmaj_deconv"].asDouble();
				double bmin_wcs_deconv= json_fitcomponents[k]["bmin_deconv"].asDouble();
				double pa_wcs_deconv= json_fitcomponents[k]["pa_deconv"].asDouble();
				#ifdef LOGGING_ENABLED
					DEBUG_LOG("Source island "<<sname<<" deconv pars ("<<bmaj_wcs_deconv<<","<<bmin_wcs_deconv<<","<<pa_wcs_deconv<<") ...");
				#endif
				fitPars.SetComponentFitWCSDeconvolvedEllipsePars(k, bmaj_wcs_deconv, bmin_wcs_deconv, pa_wcs_deconv);

				//Set sourceness
				std::string compSourcenessLabel= json_fitcomponents[k]["sourceness_label"].asString();
				double compSourceness= GetSourcenessId(compSourcenessLabel);
				double compSourcenessScore= json_fitcomponents[k]["sourceness_score"].asDouble();
				fitPars.SetComponentSourcenessId(k, compSourceness);
				fitPars.SetComponentSourcenessScore(k, compSourcenessScore);

				//Set morph id
				std::string compMorphLabel= json_fitcomponents[k]["morph_label"].asString();
				int compMorphId= GetSourceMorphId(compMorphLabel);
				fitPars.SetComponentMorphId(k, compMorphId);

				//Fill spectral info
				auto json_compspectralinfo= json_fitcomponents[k]["spectral_info"];
				if(!json_spectralinfo.empty())
				{
					SpectralIndexData sid;
					sid.spectralIndex= json_compspectralinfo["alpha"].asDouble();
					sid.spectralIndexErr= json_compspectralinfo["alpha_err"].asDouble();
					sid.isSpectralIndexFit= static_cast<bool>(json_compspectralinfo["fit"].asInt());
					sid.spectralFitChi2= json_compspectralinfo["chi2"].asDouble();
					sid.spectralFitNDF= json_compspectralinfo["ndf"].asInt();
					sids.push_back(sid);
				}

			}//end fit components
		}//close if fit components

		//Set fit pars
		source->SetFitPars(fitPars); 		

		//Set component spectral index data
		source->SetComponentSpectralIndexData(sids);

	}//close fit info


	return source;

}//close ReadSourceIslandFromJson()

}//close namespace
