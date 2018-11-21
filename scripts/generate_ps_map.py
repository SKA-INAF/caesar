#!/usr/bin/env python

##################################################
###          MODULE IMPORT
##################################################
## STANDARD MODULES
import os
import sys
import subprocess
import string
import time
import signal
from threading import Thread
import datetime
import numpy as np
import random
import math
##from ctypes import *

## ASTRO
from scipy import ndimage
##import pyfits
from astropy.io import fits
from astropy.units import Quantity
from astropy.modeling.parameters import Parameter
from astropy.modeling.core import Fittable2DModel
from astropy.modeling.models import Box2D, Gaussian2D, Ring2D, Ellipse2D, TrapezoidDisk2D, Disk2D, AiryDisk2D, Sersic2D
from astropy import wcs

## ROOT
import ROOT
from ROOT import gSystem, TFile, TTree, gROOT, AddressOf

## CAESAR
gSystem.Load('libCaesar')
from ROOT import Caesar

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## Graphics modules
import matplotlib.pyplot as plt
import pylab
##################################################

#### GET SCRIPT ARGS ####
def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

###########################
##     ARGS
###########################
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")

	## Needed to trick CASA when running in batch	
	parser.add_argument('-c', dest='scriptname', required=False, type=str, default='',action='store',help='Script name')

	# - RUN
	parser.add_argument('--no-convolve', dest='do_convolve', action='store_false')	
	parser.set_defaults(do_convolve=True)

	# - INPUTS
	parser.add_argument('-mosaic', '--mosaic', dest='mosaic', required=True, type=str, default='', action='store',help='Mosaic Image in units Jy/beam (.fits)')
	parser.add_argument('-residual', '--residual', dest='residual', required=True, type=str, default='', action='store',help='Residual Image in units Jy/pixel (.fits)')
	parser.add_argument('-model', '--model', dest='model', required=False, type=str, default='', action='store',help='Sky Model Image (if given override model generation) in units Jy/pixel (.fits)')
	parser.add_argument('-img', '--img', dest='img', required=True, type=str, default='', action='store',help='Restored Image (.fits)')
		

	# - POINT SOURCE GENERATION OPTIONS
	parser.add_argument('-marginx', '--marginx', dest='marginx', required=False, type=int, default=0,action='store',help='Image x margin in pixels')
	parser.add_argument('-marginy', '--marginy', dest='marginy', required=False, type=int, default=0,action='store',help='Image y margin in pixels')
	parser.add_argument('-nsources', '--nsources', dest='nsources', required=False, type=int, default=0, action='store',help='Number of point-source to be generated (if >0 overrides the density generation) (default=0)')
	parser.add_argument('-source_density', '--source_density', dest='source_density', required=False, type=float, default=1000, action='store',help='Compact source density (default=1000)')
	parser.add_argument('-Smin', '--Smin', dest='Smin', required=False, type=float, default=1.e-6, action='store',help='Minimum source flux in Jy (default=1.e-6)')
	parser.add_argument('-Smax', '--Smax', dest='Smax', required=False, type=float, default=1, action='store',help='Maximum source flux in Jy (default=1)')
	
	# - OUTPUT FILE OPTIONS
	parser.add_argument('-outfile_img', '--outfile_img', dest='outfile_img', required=False, type=str, default='simmap.fits',action='store',help='Output filename')
	parser.add_argument('-outfile_model', '--outfile_model', dest='outfile_model', required=False, type=str, default='skymodel.fits', action='store',help='Model filename')
	parser.add_argument('-outfile_modelconv', '--outfile_modelconv', dest='outfile_modelconv', required=False, type=str, default='skymodel_conv.fits', action='store',help='Model convolved filename')
	parser.add_argument('-outfile_sources', '--outfile_sources', dest='outfile_sources', required=False, type=str, default='sources.dat',action='store',help='Ascii file with list of generated sources')
	parser.add_argument('-outfile_ds9region', '--outfile_ds9region', dest='outfile_ds9region', required=False, type=str, default='dsregion.reg',action='store',help='DS9 source region filename')
	
	args = parser.parse_args()	

	return args




###########################
##     SIMULATOR CLASS
###########################

class SkyMapSimulator(object):

	""" Sky map simulator class

			Attributes:
				mosaic_file: mosaic image (.fits)
				ny: image height in pixels
				pixsize: pixel size in arcsec (default=1)
	"""

	def __init__(self, mosaic_file):
		""" Return a SkyMapGenerator object """

		## Input maps
		self.mosaic_file= mosaic_file
		self.mosaic_im= None
		self.mosaic_cs= None

		self.model_data= None
		self.model_im= None

		## Image parameters
		self.nx= 0
		self.ny= 0
		self.marginx= 0 # in pixels (no margin)
		self.marginy= 0 # in pixels (no margin)
		self.pixsize= 0 # in arcsec
		
		self.beam_bmaj= 6.5 # in arcsec
		self.beam_bmin= 6.5 # in arcsec
		self.beam_bpa= 0 # in deg
		
		
		## Compact source parameters
		self.generate_skymodel= True
		self.ps_list= []
		self.nsources= 0 # default is density generator
		self.source_density= 2000. # in sources/deg^2
		self.Smin= 1.e-6 # in Jy 
		self.Smax= 1 # in Jy
		
		## Map output file
		self.model_outfile= 'skymodel.fits'
		
		## DS9 & ascii output file
		self.ps_outfile= 'sources.dat'
		self.ds9_outfile= 'ds9region.reg'

	
	def set_margins(self,marginx,marginy):
		""" Set margin in X & Y """
		#if (marginx<0 or marginy<0 or marginx>=self.nx/2 or marginy>=self.ny/2) :
		#	raise ValueError('Invalid margin specified (<0 or larger than image half size!')
		self.marginx= marginx
		self.marginy= marginy

	def set_ds9region_outfile(self,filename):
		""" Set the output DS9 region filename """
		self.ds9_outfile= filename

	def set_source_outfile(self,filename):
		""" Set the output source list filename """
		self.ps_outfile= filename


	def set_model_outfile(self,filename):
		""" Set the output model filename """
		self.model_outfile= filename

	
	def set_source_flux_range(self,Smin,Smax):
		""" Set source flux range """
		self.Smin= Smin
		self.Smax= Smax

	def set_nsources(self,n):
		""" Set number of sources to be generated """
		if n<0:	
			raise ValueError('Invalid number of sources specified (shall be >=0)')
		self.nsources= n

	def set_source_density(self,density):
		""" Set compact source density in deg^-2 """
		self.source_density= density


	def write_source_list(self):
		""" Write source list to file """

		## Open file	
		fout = open(self.ps_outfile, 'wb')
	
		## Write data
		for i in range(len(self.ps_list)):
			name= self.ps_list[i][0]
			x= self.ps_list[i][1]
			y= self.ps_list[i][2]
			S= self.ps_list[i][3]
			data= (("%s %s %s %s") % (x,y,S,name) )

			fout.write(data)
			fout.write('\n')			

		fout.close();


	def write_ds9_regions(self):
		""" Write DS9 region """
	
		## Open file	
		fout = open(self.ds9_outfile, 'wb')
	
		## Write file header
		fout.write('global color=red font=\"helvetica 8 normal\" edit=1 move=1 delete=1 include=1\n')
		fout.write('image\n')

		for i in range(len(self.ps_list)):
			name= self.ps_list[i][0]
			x= self.ps_list[i][1]+1 # NB: DS9 starts index from 1
			y= self.ps_list[i][2]+1
			S= self.ps_list[i][3]
			# point=[circle|box|diamond|cross|x|arrow|boxcircle] [size]
			#region= (("circle point %s %s # %s") % (x,y,name) )
			bmaj_pix= self.beam_bmaj/self.pixsize # in pixels
			bmin_pix= self.beam_bmin/self.pixsize # in pixels
			bpa= self.beam_bpa
			region= (("ellipse %s %s %s %s %s # %s") % (x,y,bmaj_pix,bmin_pix,bpa,name) )

			fout.write(region)
			fout.write('\n')			

		fout.close();

	
	def write_data_to_fits(self,im,outputfile):
		""" Write casa image to FITS """
		im.tofits(outputfile,overwrite=True)
		

	def init(self):
		""" Initialize data """
		
		## Check if residual map is empty
		if not self.mosaic_file:
			raise ValueError('Missing mosaic map filename!')

		## Read residual map from file
		file_ext= os.path.splitext(self.mosaic_file)[1]
		if file_ext=='.fits':
			print 'INFO: Importing mosaic FITS file %s as CASA image ...' % self.mosaic_file
			self.mosaic_im= ia.newimagefromfits(infile=self.mosaic_file)
		else:
			print 'INFO: Importing mosaic CASA image %s ' % self.mosaic_file
			self.mosaic_im= ia.newimagefromfile(infile=self.mosaic_file)

		## Check image read
		if not self.mosaic_im:
			errmsg= 'ERROR: Failed to import given mosaic input file!'
			print(errmsg)
			raise ValueError(errmsg)
			
		## Read and store mosaic image info 
		self.mosaic_cs= self.mosaic_im.coordsys() 	
		self.mosaic_data= self.mosaic_im.getregion()
		data_shape= self.mosaic_im.shape()
		self.nx= data_shape[0]
		self.ny= data_shape[1]
		dx= np.abs(np.rad2deg(self.mosaic_cs.increment(type='direction')['numeric'][0]))*3600 # in arcsec
		dy= np.abs(np.rad2deg(self.mosaic_cs.increment(type='direction')['numeric'][1]))*3600 # in arcsec
		self.pixsize= min(dx,dy)
		beam= self.mosaic_im.restoringbeam()
		self.beam_bmaj= beam['major']['value'] # in arcsec
		self.beam_bmin= beam['minor']['value'] # in arcsec
		self.beam_bpa= beam['positionangle']['value'] # in deg	
		print 'INFO: Mosaic image shape=', data_shape
		print ('INFO: Restored image beam= (%s,%s,%s)') % (self.beam_bmaj,self.beam_bmin,self.beam_bpa)
					
		## Check margins
		if (self.marginx<0 or self.marginy<0 or self.marginx>=self.nx/2 or self.marginy>=self.ny/2) :
			raise ValueError('Invalid margin specified (<0 or larger than image half size!')
					

	def generate_compact_sources(self):
		""" Generate list of compact sources in the map.
					- Uniform spatial distribution
					- Uniform flux distribution

				Arguments:	
					density: source density in #sources/deg^2 (e.g. 2000)
		"""
		
		# Get mosaic pixels with non-nan values
		mask= ~np.isnan(self.mosaic_data)
		npix_good= np.count_nonzero(mask)

		# Compute number of sources to be generated given map area in pixels
		area= ((self.nx-2*self.marginx)*(self.ny-2*self.marginy))*self.pixsize/(3600.*3600.) # in deg^2
		 

		if self.nsources>0:
			nsources= self.nsources
		else: # density generator
			nsources= int(round(self.source_density*area))
		S_min= self.Smin
		S_max= self.Smax
		lgS_min= np.log(S_min)
		lgS_max= np.log(S_max)
		randomize_flux= False
		if self.Smin<self.Smax:
			randomize_flux= True

		print 'INFO: Generating #',nsources,' compact sources in map...'

		
		## Start generation loop
		self.model_data= ia.makearray(0,[self.nx,self.ny,1,1])
		##self.model_data= Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.nx, y_width=2*self.ny)(self.gridx, self.gridy)
		self.ps_list= []
		index= 0

		#for index in range(0,nsources):	
		while index < nsources:
			if index%100==0 :
				print ("INFO: Generating compact source no. %s/%s" % (index+1,nsources))
			
			## Generate random coordinates
			x0= np.random.uniform(self.marginx,self.nx-self.marginx-1)
			y0= np.random.uniform(self.marginy,self.ny-self.marginy-1)

			## Compute amplitude given significance level and bkg
			## Generate flux uniform in log
			if randomize_flux:
				lgS= np.random.uniform(lgS_min,lgS_max)
				S= np.exp(lgS)
			else:
				S= S_min
				
			## Fill sky model map
			ix= int(np.round(x0))
			iy= int(np.round(y0))
			##if self.model_data[iy,ix]>0:
			if self.model_data[ix][iy][0][0]>0 or not mask[ix][iy][0][0]:
				print ("INFO: Skip source generated @ pixel (%s,%s) as already filled or corresponding to nan value in mosaic ..." % (ix,iy) ) 
				continue

			index+= 1
			source_name= 'S' + str(index+1)
			#self.model_data[iy,ix]+= S
			self.model_data[ix][iy][0][0]+= S
			self.ps_list.append([source_name,x0,y0,S])
			print ('INFO: Source %s: Pos(%s,%s), ix=%s, iy=%s, S=%s' % (source_name,str(x0),str(y0),str(ix),str(iy),str(S)))

		## Create skymodel image
		self.model_im= ia.newimagefromarray(pixels=self.model_data, csys=self.mosaic_cs.torecord())  
		self.model_im.setbrightnessunit('Jy/pixel')
		print ('INFO: model units=%s') % self.model_im.brightnessunit()	
		print 'INFO: model shape=', self.model_im.shape()


		print ('INFO: End source generation')

		

	#####################################
	###         RUN           ##
	#####################################
	def run(self):
		""" Generate sky map and convolve with beam """
		
		## == INITIALIZE DATA ==
		print ('INFO: Initializing simulator data...')
		self.init()

		## == GENERATE SKYMODEL WITH POINT SOURCES ==
		print ('INFO: Generating compact sources...')
		self.generate_compact_sources()
			
		## == WRITE MAPS TO FITS FILES ==
		print ('INFO: Writing model to FITS...')
		self.write_data_to_fits(self.model_im,self.model_outfile)
		
		## == WRITE DS9 REGION FILE ==
		print ('INFO: Writing DS9 regions ...')
		self.write_ds9_regions()

		print ('INFO: Writing source list ...')
		self.write_source_list()

		## Close open CASA images
		self.mosaic_im.done()
		self.model_im.done()











###########################
##     IMAGER CLASS
###########################

class SkyMapImager(object):

	""" Sky map imager class

			Attributes:
				residual: beam residual image (.fits)
				model: skymodel image (.fits)
				img_file: beam restored image (.fits)
	"""

	def __init__(self, residual_file, model_file, img_file):
		""" Return a SkyMapImager object """

		## Input maps
		self.residual_file= residual_file
		self.residual_im= None
		self.residual_data= None
		
		self.model_file= model_file
		self.model_im= None
		self.model_data= None

		self.img_file= img_file
		self.img_im= None
		
		self.model_im_conv= None
		self.model_conv_data= None
		
		self.img_conv_im= None

		
		## Initialize restoring beam
		self.beam_bmaj= 6.5 # in arcsec
		self.beam_bmin= 6.5 # in arcsec
		self.beam_bpa= 0 # in deg
		
		## Map output file
		self.img_outfile= 'simmap.fits'
		self.model_conv_outfile= 'skymodel_conv.fits'
		
		
	def set_map_outfile(self,filename):
		""" Set the output map filename """
		self.img_outfile= filename

	
	def set_modelconv_outfile(self,filename):
		""" Set the output model colvolved filename """
		self.model_conv_outfile= filename

	
	def write_data_to_fits(self,im,outputfile):
		""" Write casa image to FITS """
		im.tofits(outputfile,overwrite=True)
		

	def init(self):
		""" Initialize imager """
		
		## Read residual map from file
		if self.residual_file:
			print 'INFO: Importing residual FITS file %s as CASA image ...' % self.residual_file
			self.residual_im= ia.newimagefromfits(infile=self.residual_file)
			self.residual_data= self.residual_im.getregion()
			print 'INFO: Residual image shape=', self.residual_im.shape()	

			## Get bounding box (used to extract submap from model)
			bb= self.residual_im.boundingbox()
			bl_coords= bb['blcf'].split(',')
			tr_coords= bb['trcf'].split(',')
			x0_wcs= bl_coords[0]
			x1_wcs= tr_coords[0]
			y0_wcs= bl_coords[1]
			y1_wcs= tr_coords[1]
			cs= self.residual_im.coordsys()
			xref_wcs= np.rad2deg(cs.referencevalue(type='direction')['numeric'][0]) # in deg
			yref_wcs= np.rad2deg(cs.referencevalue(type='direction')['numeric'][1]) # in deg
			nx= self.residual_im.shape()[0]
			ny= self.residual_im.shape()[1]

			print('INFO: (nx,ny)=(%s,%s), pixref_wcs(%s,%s)') % (nx,ny,xref_wcs,yref_wcs)

		else:
			raise ValueError('Missing residual map filename!')
			

		## Read restored image from file
		if self.img_file:
			print 'INFO: Importing restored image FITS file %s as CASA image ...' % self.img_file
			self.img_im= ia.newimagefromfits(infile=self.img_file)

			beam= self.img_im.restoringbeam()
			self.beam_bmaj= beam['major']['value'] # in arcsec
			self.beam_bmin= beam['minor']['value'] # in arcsec
			self.beam_bpa= beam['positionangle']['value'] # in deg	
			print ('INFO: Restored image beam= (%s,%s,%s)') % (self.beam_bmaj,self.beam_bmin,self.beam_bpa)
			print 'INFO: Restored image shape=', self.img_im.shape()	
		else:
			raise ValueError('Missing restored image input filename!')
			
			
		## Import model map from file
		if self.model_file:
			## Read model image
			print 'INFO: Importing model FITS map from file %s as CASA image ' % self.model_file
			ia.open(self.model_file)

			## Create region from bounding box of beam residual map
			#print ('INFO: Creating region from bounding box of beam residual map [(x0,y0),(x1,y1)]=[(%s,%s)(%s,%s)]') % (x0_wcs,y0_wcs,x1_wcs,y1_wcs)
			#x0= ia.topixel([x0_wcs,y0_wcs])['numeric'][0]
			#y0= ia.topixel([x0_wcs,y0_wcs])['numeric'][1]
			#x1= ia.topixel([x1_wcs,y1_wcs])['numeric'][0]
			#y1= ia.topixel([x1_wcs,y1_wcs])['numeric'][1]
			xref= ia.topixel([str(xref_wcs)+'deg',str(yref_wcs)+'deg'])['numeric'][0]
			yref= ia.topixel([str(xref_wcs)+'deg',str(yref_wcs)+'deg'])['numeric'][1]
			x0= xref-nx/2
			x1= xref+nx/2-1
			y0= yref-ny/2
			y1= yref+ny/2-1
			print ('INFO: Creating region from bounding box of beam residual map pixref_wcs(%s,%s), pixref(%s,%s), [(x0,y0),(x1,y1)]=[(%s,%s)(%s,%s)]') % (xref_wcs,yref_wcs,xref,yref,x0,y0,x1,y1)
			
			reg= rg.box(blc=[x0,y0],trc=[x1,y1])
			
			## Extract model subimage
			print ('INFO: Extracting model subimage using region [(x0,y0),(x1,y1)]=[(%s,%s)(%s,%s)]') % (x0,y0,x1,y1)
			self.model_im= ia.subimage(region=reg)
			self.model_data= self.model_im.getregion()
			print ('INFO: model units=%s') % self.model_im.brightnessunit()	
			print 'INFO: model shape=', self.model_im.shape()	

			## Close model image
			ia.done

		else:
			raise ValueError('Missing skymodel image input filename!')

		
		
	def convolve_skymodel(self):
		""" Convolve sky model with restored beam """
				
		## Convolving skymodel with restoring beam
		## NB: If input map is in Jy/pixel units, the convolve task automatically scale results to Jy/beam (no further scaling needed)
		bmaj= str(self.beam_bmaj) + 'arcsec'
		bmin= str(self.beam_bmin) + 'arcsec'
		bpa= str(self.beam_bpa) + 'deg'
		print ('INFO: Convolving skymodel with restored beam (%s,%s,%s)') % (bmaj,bmin,bpa)
		self.model_im_conv= self.model_im.convolve2d(type='gauss',major=bmaj,minor=bmin,pa=bpa)

		## Get pixel data from convolved image
		self.model_conv_data= self.model_im_conv.getregion()


	#####################################
	###         RUN           ##
	#####################################
	def run(self):
		""" Convolve skymodel with beam """
		
		## == INITIALIZE DATA ==
		print ('INFO: Init and read input maps ...')
		self.init()
	
		## == CONVOLVE SKY MODEL WITH BEAM ==
		print ('INFO: Convolving skymodel with beam ...')
		self.convolve_skymodel()		

		
		## == MAKE FINAL MAP ==
		print ('INFO: Creating final map I=res + model_convolved ...')
		print 'INFO: residual data size=',self.residual_data.shape
		print 'INFO: model conv data size=',self.model_conv_data.shape	
		data= self.residual_data + self.model_conv_data
		cs = self.img_im.coordsys()
		self.img_conv_im= ia.newimagefromarray(pixels=data, csys=cs.torecord())  		

		## == WRITE MAPS TO FITS FILES ==
		print ('INFO: Writing images to FITS...')
		self.write_data_to_fits(self.img_conv_im,self.img_outfile)
		self.write_data_to_fits(self.model_im_conv,self.model_conv_outfile)

		## Close open CASA images
		self.model_im.done()
		self.img_im.done()
		self.residual_im.done()



##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   Get script args
	#===========================
	print('Get script args')
	try:
		args= get_args()
	except Exception as ex:
		print("Failed to get and parse options (err=%s)",str(ex))
		return 1
	
	# - Run options
	do_convolve= args.do_convolve	

	# - Input maps
	mosaic_file= args.mosaic
	residual_file= args.residual
	img_file= args.img
	model_file= args.model
	
	# - Compact source args
	marginX= args.marginx
	marginY= args.marginy
	Smin= args.Smin
	Smax= args.Smax
	source_density= args.source_density
	nsources= args.nsources

	# - Output args
	outfile_img= args.outfile_img
	outfile_sources= args.outfile_sources
	outfile_ds9region= args.outfile_ds9region
	outfile_model= args.outfile_model
	outfile_modelconv= args.outfile_modelconv

	print("*** ARGS ***")
	print("Margin X: %s" % marginX)
	print("Margin Y: %s" % marginY)
	print("Source flux range: (%s,%s)" % (Smin, Smax))
	print("Source density (deg^-2): %s" % source_density)	
	print("Image output filename : %s " % outfile_img)
	print("************")

	## Check if generate sky model
	if model_file:
		generate_model= False
		model_file_input= model_file
	else:
		generate_model= True
		model_file_input= outfile_model


	## Generate simulated sky map
	if generate_model:
		print ('INFO: Generate simulated skymodel...')
		simulator= SkyMapSimulator(mosaic_file)
		simulator.set_margins(marginX,marginY)
		simulator.set_model_outfile(outfile_model)
		simulator.set_ds9region_outfile(outfile_ds9region)
		simulator.set_source_outfile(outfile_sources)
		simulator.set_nsources(nsources)
		simulator.set_source_flux_range(Smin,Smax)
		simulator.set_source_density(source_density)
		simulator.run()

	## Generate restored image by convolving model with restored beam
	if do_convolve:
		print ('INFO: Convolve skymodel with beam ...')
		imager= SkyMapImager(residual_file,model_file_input,img_file)
		imager.set_map_outfile(outfile_img)
		imager.set_modelconv_outfile(outfile_modelconv)
		imager.run()
	
	print 'INFO: End run'
	
###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

