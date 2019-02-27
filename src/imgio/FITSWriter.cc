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
* @file FITSWriter.cc
* @class FITSWriter
* @brief FITSWriter
*
* FITS Image Writer class
* @author S. Riggi
* @date 20/01/2015
*/

#include <FITSWriter.h>
#include <Image.h>

//CFITSIO headers
#include <fitsio.h>

//Python interface
//#include <TPython.h>
//#include <Python.h>

#include <TObject.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
using namespace std;

ClassImp(Caesar::FITSWriter)

namespace Caesar {


FITSWriter::FITSWriter() {

}//close costructor


FITSWriter::~FITSWriter() {

}//close destructor

/*
int FITSWriter::Init(){
	
	try {
		//Initialize
		//Py_Initialize();

		//## Import python modules
		//== Numpy ==
		if(!TPython::Exec("import numpy as np")){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to import numpy python module!");
			#endif
			return -1;
		}
		//== pyfits ==
		if(!TPython::Exec("import pyfits")){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to import pyfits python module!");
			#endif
			return -1;
		}
	}
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to initialize python interface!");
		#endif
		return -1;
	}

	return 0;

}//close Init()
*/



int FITSWriter::WriteFITS(Image* img,std::string outfilename,bool recreate)
{
	//## Check image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return -1;
	}
	if(outfilename==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output filename!");
		#endif
		return -1;
	}
	if(recreate){
		outfilename= std::string("!") + outfilename;
	}

	//## Init vars
	fitsfile* fp= 0;
  int status = 0;

	//## Open file with filter
	#ifdef LOGGING_ENABLED	
		DEBUG_LOG("Creating FITS file "<<outfilename<<" to be written ...");
	#endif
  fits_create_file(&fp, outfilename.c_str(), &status);
  if (status) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to create FITS file "<<outfilename<<" for writing image data!");
		#endif
		HandleError(status,fp);
		return -1;
	}

	//## Create image
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Making image in FITS file "<<outfilename<<" ot be written ...");
	#endif
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	long naxes[2]= {Nx,Ny};
	fits_create_img(fp,FLOAT_IMG,2, naxes, &status);
	if (status) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to create image in FITS file "<<outfilename<<"!");
		#endif
		HandleError(status,fp);
		return -1;
	}
	
	//## Fill image data
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Filling image data in FITS file "<<outfilename<<" ot be written ...");
	#endif
	long fpixel[2]= {1,1};
	long nelements= Nx*Ny;
	float* imgdata= (float*)((img->GetPixels()).data());
	fits_write_pix(fp, TFLOAT, fpixel, nelements, imgdata, &status);
	if (status) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to fill image data in FITS file "<<outfilename<<"!");
		#endif
		HandleError(status,fp);
		return -1;
	}

	//## Writing additional keywords in header
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Filling HDU keywords in FITS file "<<outfilename<<" ot be written ...");
	#endif
	if(WriteFITSKeywords(img,fp)<0){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Failed to write keywords in FITS file!");
		#endif
	}

	//## Close image
	if (fp) {
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Closing FITS file "<<outfilename<<"...");
		#endif
		fits_close_file(fp, &status);
	}

	return 0;

}//close WriteFITS()


int FITSWriter::WriteFITSKeywords(Image* img,fitsfile* fp)
{
	//Check if image has metadata
	if(!img->HasMetaData()){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Input image has no metadata, nothing to be written to FITS keywords...");
		#endif
		return 0;
	}
	
	//Fill metadata
	int status= 0;
	ImgMetaData* metadata= img->GetMetaData();

	//Reference pixel id
	double Cx= static_cast<double>(metadata->Cx);
	double Cy= static_cast<double>(metadata->Cy);
	fits_write_key(fp,TDOUBLE,"CRPIX1", &Cx,"Pixel reference x image coord",&status);
	fits_write_key(fp,TDOUBLE,"CRPIX2", &Cy,"Pixel reference y image coord",&status);

	//Reference pixel coords
	double Xc= static_cast<double>(metadata->Xc);
	double Yc= static_cast<double>(metadata->Yc);
	fits_write_key(fp,TDOUBLE,"CRVAL1", &Xc,"Pixel reference x phys coord",&status);
	fits_write_key(fp,TDOUBLE,"CRVAL2", &Yc,"Pixel reference y phys coord",&status);
	
	//Pixel size
	double dX= static_cast<double>(metadata->dX);
	double dY= static_cast<double>(metadata->dY);
	fits_write_key(fp,TDOUBLE,"CDELT1", &dX,"Pixel size x (arcsec)",&status);
	fits_write_key(fp,TDOUBLE,"CDELT2", &dY,"Pixel size y (arcsec)",&status);
	
	//System rotation info
	double RotX= static_cast<double>(metadata->RotX);
	double RotY= static_cast<double>(metadata->RotY);
	fits_write_key(fp,TDOUBLE,"CROTA1", &RotX,"Coord system rotation x angle (deg)",&status);
	fits_write_key(fp,TDOUBLE,"CROTA2", &RotY,"Coord system rotation y angle (deg)",&status);
	
	//Type of astro coords
	std::string CoordTypeX= metadata->CoordTypeX;
	std::string CoordTypeY= metadata->CoordTypeY;
	fits_write_key(fp,TSTRING,"CTYPE1", (char*)(CoordTypeX.c_str()),"Coord system x type",&status);
	fits_write_key(fp,TSTRING,"CTYPE2", (char*)(CoordTypeY.c_str()),"Coord system y type",&status);

	//Units
	std::string BUnit= metadata->BUnit;
	fits_write_key(fp,TSTRING,"CTYPE1", (char*)(BUnit.c_str()),"Physical units of the image data",&status);
	
	//Beam info
	double Bmaj= metadata->Bmaj;
	double Bmin= metadata->Bmin;
	double Bpa= metadata->Bpa;
	fits_write_key(fp,TDOUBLE,"BMAJ", &Bmaj,"Beam ellipse FWHM major axis (deg)",&status);
	fits_write_key(fp,TDOUBLE,"BMIN", &Bmin,"Beam ellipse FWHM minor axis (deg)",&status);
	fits_write_key(fp,TDOUBLE,"BPA", &Bpa,"Beam ellipse rotation angle of the major axis measuring from North (up) counterclockwise (deg)",&status);
	
	//Obs Epoch (DEPRECATED in FITS STANDARD (use EQUINOX?)
	//double Epoch;

	return 0;

}//close WriteFITSKeywords()


void FITSWriter::HandleError(int& status,fitsfile* fp)
{
	char errdescr[FLEN_STATUS+1];
	fits_get_errstatus(status, errdescr);
	#ifdef LOGGING_ENABLED
 		ERROR_LOG("Error while processing FITS file (reason="<<errdescr<<")");
	#endif
  status = 0;
  if (fp) fits_close_file(fp, &status);

}//close HandleError()

/*
int FITSWriter::WriteFITS(Image* img,std::string outfilename){

	//## Check image
	if(!img){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Null ptr to given image!");
		#endif
		return -1;
	}
	if(outfilename==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty output filename!");
		#endif
		return -1;
	}

	//## Initialize and load needed python modules
	if(Init()<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Initialization failed!");
		#endif
		return -1;	
	}

	//## Load data to python
	//outfilename
	PyObject* pystr = TPython::ObjectProxy_FromVoidPtr(&outfilename, "std::string");

	//image data
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	std::vector< std::vector<float> > v;
	for(long int j=0;j<Ny;j++){
		v.push_back( std::vector<float>() );
		for(long int i=0;i<Nx;i++){
			double w= img->GetPixelValue(i,j);
			v[j].push_back(w);
		}
	}
	PyObject* pyvec2d= TPython::ObjectProxy_FromVoidPtr(&v,"std::vector< std::vector<float> >");

	PyObject* pymain = PyImport_ImportModule("__main__");
	PyModule_AddObject(pymain, "outfile", pystr);	
	PyModule_AddObject(pymain, "mat", pyvec2d);
	Py_DECREF(pymain);
	
	//## Print the imported data
	if(!TPython::Exec("print 'outfile:',outfile")){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to import outfilename in python!");
		#endif
		return -1;
	}
	if(!TPython::Exec("print 'mat: ',mat")){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to import mat in python!");
		#endif
		return -1;
	}
	
  	
	//## Bind ROOT img obj
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Binding image to python...");
	#endif
	if(!TPython::Bind(img,"img")){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to bind image on the ROOT-Python interface!");
		#endif
		return -1;	
	}
	bool hasMetaData= img->HasMetaData();
	Caesar::ImgMetaData* metadata= 0;
	if(hasMetaData){
		#ifdef LOGGING_ENABLED
			DEBUG_LOG("Binding image meta-data to python...");
		#endif
		metadata= img->GetMetaData();
		if(!TPython::Bind(metadata,"metadata")){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to bind image metadata on the ROOT-Python interface!");
			#endif
			hasMetaData= false;
		}
	}

	
	//## Create a HDU object 	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Creating a HDU from the input image...");
	#endif
	if(!TPython::Exec("hdu = pyfits.PrimaryHDU(mat)")){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to create hdu from bin array!");
		#endif
		return -1;
	}

	//## Fill header info	
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Filling hdu header...");
	#endif
	try {
		TPython::Exec("hdu.header[\'NAXIS1\']= metadata.Nx"); 
		TPython::Exec("hdu.header[\'NAXIS2\']= metadata.Ny"); 
		TPython::Exec("hdu.header[\'CTYPE1\']= metadata.CoordTypeX"); 
		TPython::Exec("hdu.header[\'CTYPE2\']= metadata.CoordTypeY"); 
		TPython::Exec("hdu.header[\'CRPIX1\']= metadata.Cx"); 
		TPython::Exec("hdu.header[\'CRPIX2\']= metadata.Cy"); 
		TPython::Exec("hdu.header[\'CRVAL1\']= metadata.Xc"); 
		TPython::Exec("hdu.header[\'CRVAL2\']= metadata.Yc"); 
		TPython::Exec("hdu.header[\'CDELT1\']= metadata.dX"); 
		TPython::Exec("hdu.header[\'CDELT2\']= metadata.dY"); 
		TPython::Exec("hdu.header[\'CROTA1\']= metadata.RotX"); 
		TPython::Exec("hdu.header[\'CROTA2\']= metadata.RotY"); 
		TPython::Exec("hdu.header[\'BUNIT\']= metadata.BUnit"); 
		TPython::Exec("hdu.header[\'BMAJ\']= metadata.Bmaj"); 
		TPython::Exec("hdu.header[\'BMIN\']= metadata.Bmin"); 
		TPython::Exec("hdu.header[\'BPA\']= metadata.Bpa"); 
		TPython::Exec("hdu.header[\'EPOCH\']= metadata.Epoch"); 
	}//close try
	catch(...){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to fill hdu header vars!");
		#endif
		return -1;
	}

	//## Write to FITS file
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Writing to fits...");
	#endif
	if(!TPython::Exec("hdu.writeto(str(outfile))")){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to write FITS file!");
		#endif
		return -1;
	}

	//## Finalize
	//cout<<"FITSWriter::WriteFITS(): INFO: Finalizing..."<<endl;
	//Py_DECREF(pymain);
	//Py_Finalize();
	
	return 0;

}//close WriteFITS()
*/


}//close namespace
