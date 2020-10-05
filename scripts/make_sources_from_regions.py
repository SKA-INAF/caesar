#!/usr/bin/env python

from __future__ import print_function

## STANDARD MODULES
import sys
import numpy as np
import os
import re
import json
from collections import defaultdict
import operator as op

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## ASTROPY MODULES
from astropy.io import fits
import regions
from radio_beam import Beam
#from regions import read_ds9

## GRAPHICS MODULES
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches


## LOGGER
import logging
import logging.config
logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)-15s %(levelname)s - %(message)s",datefmt='%Y-%m-%d %H:%M:%S')
logger= logging.getLogger(__name__)
logger.setLevel(logging.INFO)

## ROOT + CAESAR
try:
	## ROOT
	import ROOT
	from ROOT import gSystem, gROOT, AddressOf

	## CAESAR
	gSystem.Load('libCaesar.so')
	from ROOT import Caesar
	from ROOT.Caesar import Image

except:
	logger.error("Cannot load ROOT & Caesar modules!")
	exit(1)


def make_caesar_source(source_data,source_name,source_id,source_type,offsetx=0,offsety=0,beamArea=0):
	""" Create Caesar source from source data array """

	# - Create Caesar source
	source= Caesar.Source()

	# - Get source indexes and fill pixels in Caesar source
	source_indexes= np.column_stack(np.where(source_data!=0))
	nRows= (source_data.shape)[0]	
	nCols= (source_data.shape)[1]	

	for index in source_indexes:
		rowId= index[0]
		colId= index[1]
		S= source_data[rowId,colId]
		ix= colId
		iy= rowId
		#iy= nRows-1-rowId
		gbin= ix + iy*nCols
		pixel= Caesar.Pixel(gbin,ix,iy,ix+offsetx,iy+offsety,S)
		source.AddPixel(pixel)

		# Is at edge
		if (ix==0) or (ix==nCols-1) or (iy==0) or (iy==nRows-1):
			source.SetEdgeFlag(True)

	# - Retun None if npixels is too small	
	nPix= source.GetNPixels()
	npixels_min= 3
	if nPix<npixels_min:
		logger.warn("Too few pixels (%d) for this source, return None!" % nPix)
		return None
			
	# Set some flags
	source.SetName(source_name)
	source.SetId(source_id)
	source.SetType(source_type)
	source.SetFlag(Caesar.eCandidate)
	
	# Set flux correction factor
	source.SetBeamFluxIntegral(beamArea)

	# Compute stats & morph pars
	source.ComputeStats();
	source.ComputeMorphologyParams();

	return source


#### GET SCRIPT ARGS ####
def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')


def find_duplicates(seq):
	""" Return dict with duplicated item in list"""
	tally = defaultdict(list)
	for i,item in enumerate(seq):
		tally[item].append(i)

  #return ({key:locs} for key,locs in tally.items() if len(locs)>0)
	return (locs for key,locs in tally.items() if len(locs)>0)


###########################
##     ARGS
###########################
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")

	# - Input options
	parser.add_argument('-img','--img', dest='img', required=True, type=str, help='Input image filename (.fits)') 
	parser.add_argument('-region','--region', dest='region', required=True, type=str, help='DS9 region filename (.reg)') 
	
	# - Output options
	parser.add_argument('-outfile','--outfile', dest='outfile', required=False, type=str, default='sources.root', help='Output Caesar ROOT filename (.root)') 
	
	args = parser.parse_args()	

	return args



##############
##   MAIN   ##
##############
def main():
	"""Main function"""
	
	#===========================
	#==   PARSE ARGS
	#===========================
	logger.info("Get script args ...")
	try:
		args= get_args()
	except Exception as ex:
		logger.error("Failed to get and parse options (err=%s)",str(ex))
		return 1

	# - Input filelist
	inputfile= args.img
	regionfile= args.region

	# - Output file
	outfilename= args.outfile
	
	#===========================
	#==   INIT
	#===========================
	
	# - Initialize output tree & file
	outfile= ROOT.TFile(outfilename,'RECREATE')
	outtree= ROOT.TTree('SourceInfo','SourceInfo')
	cs = Caesar.Source()
	outtree.Branch('Source',cs)		


	#===========================
	#==   READ IMAGE
	#===========================
	logger.info("Read image ...")
	hdu= fits.open(inputfile)
	data= hdu[0].data
	header= hdu[0].header
	shape= data.shape

	# - Compute beam area
	dX= header['CDELT1']
	dY= header['CDELT2']
	bmaj= header['BMAJ']
	bmin= header['BMIN']
	pixelArea= np.abs(dX*dY);
	A= np.pi*bmaj*bmin/(4*np.log(2))
	beamArea= A/pixelArea
	
	#===========================
	#==   READ REGIONS
	#===========================
	logger.info("Read regions ...")
	region_list= regions.read_ds9(regionfile)

	logger.info("#%d regions found ..." % len(region_list))
	
	#=========================================
	#==   EXTRACT SOURCES FROM REGIONS
	#=========================================
	# - Get region names and tags
	counter= 0
	caesar_sources= ROOT.std.vector("Caesar::Source*")()
	
	for region in region_list:
		print(type(region))
		sname= region.meta['text']
	
		# - Extract source masks
		#mask= region.to_mask(mode='center')
		mask= region.to_mask(mode='subpixels')
		#mask= region.to_mask(mode='exact')
		source_data= mask.to_image(shape)
				
		# - Create Caesar source object
		source_name= sname
		source_id= counter+1
		source_type= Caesar.eCompact
		caesar_source= make_caesar_source(source_data,source_name,source_id,source_type,0,0,beamArea)
		caesar_sources.push_back(caesar_source)

		# - Fill branch
		caesar_source.Copy(cs)
		cs.Print()
		outtree.Fill()

		counter+= 1

	# - Write to file
	logger.info("Writing source Tree to file ...")
	outfile.cd()
	outtree.Write()
	outfile.Close()

	# - Write to DS9
	logger.info("Writing to DS9 ...")
	Caesar.SourceExporter.WriteToDS9("ds9.reg",caesar_sources)
	
	return 0


###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

