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

//Python interface
#include <TPython.h>
#include <Python.h>
//#include <numpy/ndarrayobject.h>
//#include <numpy/ndarraytypes.h>

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

int FITSWriter::Init(){
	
	try {
		//Initialize
		//Py_Initialize();

		//## Import python modules
		//== Numpy ==
		if(!TPython::Exec("import numpy as np")){
			ERROR_LOG("Failed to import numpy python module!");
			return -1;
		}
		//== pyfits ==
		if(!TPython::Exec("import pyfits")){
			ERROR_LOG("Failed to import pyfits python module!");
			return -1;
		}
	}
	catch(...){
		ERROR_LOG("Failed to initialize python interface!");
		return -1;
	}

	return 0;

}//close Init()


int FITSWriter::WriteFITS(Image* img,std::string outfilename){

	//## Check image
	if(!img){
		ERROR_LOG("Null ptr to given image!");
		return -1;
	}
	if(outfilename==""){
		ERROR_LOG("Empty output filename!");
		return -1;
	}

	//## Initialize and load needed python modules
	if(Init()<0){
		ERROR_LOG("Initialization failed!");
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
		ERROR_LOG("Failed to import outfilename in python!");
		return -1;
	}
	if(!TPython::Exec("print 'mat: ',mat")){
		ERROR_LOG("Failed to import mat in python!");
		return -1;
	}
	
  	
	//## Bind ROOT img obj
	INFO_LOG("Binding image to python...");
	if(!TPython::Bind(img,"img")){
		ERROR_LOG("Failed to bind image on the ROOT-Python interface!");
		return -1;	
	}
	bool hasMetaData= img->HasMetaData();
	Caesar::ImgMetaData* metadata= 0;
	if(hasMetaData){
		INFO_LOG("Binding image meta-data to python...");
		metadata= img->GetMetaData();
		if(!TPython::Bind(metadata,"metadata")){
			ERROR_LOG("Failed to bind image metadata on the ROOT-Python interface!");
			hasMetaData= false;
		}
	}

	
	//## Create a HDU object 	
	INFO_LOG("Creating a HDU from the input image...");
	if(!TPython::Exec("hdu = pyfits.PrimaryHDU(mat)")){
		ERROR_LOG("Failed to create hdu from bin array!");
		return -1;
	}

	//## Fill header info	
	INFO_LOG("Filling hdu header...");
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
		ERROR_LOG("Failed to fill hdu header vars!");
		return -1;
	}

	//## Write to FITS file
	INFO_LOG("Writing to fits...");
	if(!TPython::Exec("hdu.writeto(str(outfile))")){
		ERROR_LOG("Failed to write FITS file!");
		return -1;
	}

	//## Finalize
	//cout<<"FITSWriter::WriteFITS(): INFO: Finalizing..."<<endl;
	//Py_DECREF(pymain);
	//Py_Finalize();
	
	return 0;

}//close WriteFITS()



}//close namespace
