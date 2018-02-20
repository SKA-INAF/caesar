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
import errno

## ASTRO MODULES
from astropy.io import fits

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

#### GET SCRIPT ARGS ####
def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")
	
	# - MANDATORY OPTIONS
	parser.add_argument('-vis', '--vis', dest='vis', required=True, type=str,action='store',help='Input visibility CASA table name produced by simobserve')
	parser.add_argument('-mapsize', '--mapsize', dest='mapsize', required=True, type=int,action='store',help='Map size in pixels')
	
	# OPTIONAL OPTIONS	
	parser.add_argument('-outimage', '--outimage', dest='outimage', required=False, type=str, default='mosaic',action='store',help='Output mosaic image')
	parser.add_argument('-mask', '--mask', dest='mask', required=False, type=str, default='',action='store',help='Mask file (default=none)')
	parser.add_argument('-outproject', '--outproject', dest='outproject', required=False, type=str, default='rec',action='store',help='Output project name (default=rec)')
	parser.add_argument('-pixsize', '--pixsize', dest='pixsize', required=False, type=str, default='1arcsec',action='store',help='Pixel size (default=1 arcsec)')
	parser.add_argument('-phasecenter', '--phasecenter', dest='phasecenter', required=False, type=str, default='J2000 254.85067091580214deg -41.47631111052697deg',action='store',help='Map phase center')
	parser.add_argument('-deconvolver', '--deconvolver', dest='deconvolver', required=False, type=str, default='clark',action='store',help='Deconvolver (default=clark)')
	parser.add_argument('-gridder', '--gridder', dest='gridder', required=False, type=str, default='standard',action='store',help='Gridder (default=standard)')
	parser.add_argument('-weighting', '--weighting', dest='weighting', required=False, type=str, default='briggs',action='store',help='weighting (default=briggs)')
	parser.add_argument('-projection', '--projection', dest='projection', required=False, type=str, default='SIN',action='store',help='projection (default=SIN)')
	
	parser.add_argument('-niter', '--niter', dest='niter', required=False, type=int, default=1000,action='store',help='Clean tot number of iterations (default=1000)')
	parser.add_argument('-cycleniter', '--cycleniter', dest='cycleniter', required=False, type=int, default=-1,action='store',help='Max cycle niter (default=-1)')
	parser.add_argument('-threshold', '--threshold', dest='threshold', required=False, type=float, default=0,action='store',help='Stopping threshold in Jy (default=0)')
	##parser.add_argument('--mosaic', dest='enable_mosaic', action='store_true')		
	parser.add_argument('--no-mosaic', dest='enable_mosaic', action='store_false')
	parser.set_defaults(enable_mosaic=True)
	##parser.add_argument('--fitsout', dest='enable_fitsout', action='store_true')		
	parser.add_argument('--no-fitsout', dest='enable_fitsout', action='store_false')	
	parser.set_defaults(enable_fitsout=True)
	parser.add_argument('-fitsout', '--fitsout', dest='fitsout', required=False, type=str, default='output.fits',action='store',help='Output FITS file (default=output.fits)')

	parser.add_argument('-scales', '--scales', dest='scales', required=False, type=str, default='0',action='store',help='List of scales (in pixels) for multiscale deconvolver (comma separated) (default=0)')
	parser.add_argument('-restoringbeam', '--restoringbeam', dest='restoringbeam', required=False, type=str, default='common',action='store',help='Restoring beam to be used (default=common)')
	
	parser.add_argument('--interactive', dest='enable_interactive', action='store_true')	
	parser.set_defaults(enable_interactive=False)
	
	parser.add_argument('-c', dest='scriptname', required=False, type=str, default='',action='store',help='Script name')

	args = parser.parse_args()	

	return args


def clean_observation(vis,image_size,niter=1000,cycleniter=-1,threshold=0,mask='',imagename='',field='',cell_size='1arcsec',phase_center='',weighting='briggs',projection='SIN',deconvolver='clark',gridder='standard',scales=[0],restoringbeam='common',savefits=True,fitsout='output.fits',interactive=False):
	""" Clean observation """
	
	## Set image name if empty
	if not imagename:
		project_name= 'rec'
		recimg= 'recmap'
		imagename=project_name + '/' + recimg	

	## Cleaning
	print ('INFO: Cleaning simulated data and produce final image...')
	tclean(
		imagename=imagename,
		vis=vis,
		field=field, 
		mask=mask,
		imsize=image_size,
		phasecenter=phase_center,
		projection=projection,
		cell=cell_size,	
		niter=niter,
		cycleniter=cycleniter,
		threshold=threshold,
		deconvolver=deconvolver,
		gridder=gridder,
		weighting=weighting,
		scales=scales,
		interactive=interactive
	)

	## Exporting to FITS
	if savefits:
		exported_map= imagename + '.image' 
		print ('INFO: Exporting CASA map %s to FITS...' % exported_map)
		exportfits(imagename=exported_map, fitsimage=fitsout, history=False, overwrite=True)


##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	## Get script args
	print 'INFO: #%d arguments found: %s' % (len(sys.argv),str(sys.argv))
	try:
		args= get_args()
	except Exception as ex:
		print("ERROR: Failed to get and parse options (err=%s)",str(ex))
		return 1

	vis= args.vis
	niter= args.niter
	mask= args.mask
	mapsize= args.mapsize	
	outproject= args.outproject
	pixsize= args.pixsize
	phasecenter= args.phasecenter
	deconvolver= args.deconvolver
	gridder= args.gridder
	weighting= args.weighting
	projection= args.projection
	outimage= args.outimage
	cycleniter= args.cycleniter
	threshold= args.threshold

	enable_mosaic= args.enable_mosaic
	##mosaicimage= outproject + '_' + outimage
	mosaicimage= outproject + '/' + outimage
	
	enable_fitsout= args.enable_fitsout	
	fitsout= args.fitsout
	
	enable_interactive= args.enable_interactive
	scales_str = [str(item) for item in args.scales.split(',')]
	scales= []
	for item in scales_str:
		scale_int= int(item)
		scales.append(scale_int)

	restoringbeam= args.restoringbeam

	print("*** ARGS ***")
	print 'vis: ', vis
	print 'niter: ', niter
	print 'cycleniter: ', cycleniter
	print 'threshold: ', threshold
	print 'mask: ', mask
	print 'mapsize: ', mapsize
	print 'outproject: ', outproject
	print 'pixsize: ', pixsize
	print 'phasecenter: ', phasecenter
	print 'deconvolver: ', deconvolver
	print 'gridder: ', gridder
	print 'weighting: ', weighting
	print 'scales: ', scales
	print 'projection: ', projection
	print 'restoringbeam: ', restoringbeam
	print 'outimage: ', outimage
	print 'mosaicimage: ', mosaicimage
	print 'fitsout: ', fitsout
	print("************")
	
	# Retrieve number of pointings in simulation
	#tb= casac.table()	
	tb.open(vis + '/FIELD')
	npointings= tb.nrows()
	tb.close()
	print 'INFO: #' + str(npointings) + ' pointings present in input data...'
	
	## Image each pointing
	print ("INFO: Starting imaging pointing-by-pointing...")
	t_start = time.time()

	weightmap_list= []
	cleanmap_list= []
	bmaj_list= []
	bmin_list= []
	bpa_list= []
	bunit_list= []
	for field in range(0,npointings):
		
		imgname= outproject + '_field' + str(field) + '/recmap'
		##imgname= outproject + '/field' + str(field) + '/recmap'
		weightmap= imgname + '.pb'
		cleanmap= imgname + '.image'
		weightmap_list.append(weightmap) 
		cleanmap_list.append(cleanmap) 
		print 'INFO: Imaging field no. ', field
		clean_observation(
			vis=vis,
			niter=niter,
			cycleniter=cycleniter,
			threshold=threshold,
			mask=mask,
			imagename=imgname,
			field=str(field),
			image_size=mapsize,
			cell_size=pixsize,
			phase_center=phasecenter,
			weighting=weighting,
			projection=projection,	
			deconvolver=deconvolver,
			gridder=gridder,
			restoringbeam=restoringbeam,
			scales=scales,
			savefits=False,
			fitsout='tmpfield.fits',
			interactive=enable_interactive
		)
		
		## Retrieve beam information per each cleaned map
		print 'INFO: Get CASA image header...'
		imghead= imhead(imagename=cleanmap,mode='list')
		bunit= imghead['bunit']
		bmaj= imghead['beammajor']['value']
		bmin= imghead['beamminor']['value']
		bpa= imghead['beampa']['value']
		bmaj_list.append(bmaj)
		bmin_list.append(bmin)
		bpa_list.append(bpa)
		bunit_list.append(bunit)
	
	print 'weightmap_list=',weightmap_list
	print 'cleanmap_list=',cleanmap_list
	print 'bmaj_list=',bmaj_list
	print 'bmin_list=',bmin_list
	print 'bpa_list=',bpa_list
	print 'bunit_list=',bunit_list
	
	## Make linear mosaic with all fields
	if enable_mosaic:
		lm= casac.linearmosaic()
		lm.defineoutputimage(nx=mapsize,ny=mapsize,cellx=pixsize,celly=pixsize, imagecenter=phasecenter,outputimage=mosaicimage)
		lm.setlinmostype('optimal')
		##lm.makemosaic(images=['rec_field0/recmap.image','rec_field1/recmap.image'],weightimages=['rec_field0/recmap.pb','rec_field1/recmap.pb'])
		lm.makemosaic(images=cleanmap_list,weightimages=weightmap_list)
		#lm.saultweightimage('test_sault.linmos') 

		## Set missing field in mosaic image header
		## NB: beam info are missing. Each field has its own restoring beam info. Which one to use? Using the largest bmaj. Alternatively one can set the restoring beam in tclean (FIX ME)
		largest_bmaj_index= bmaj_list.index(max(bmaj_list))	
		cleanmap_head= imhead(imagename=cleanmap_list[largest_bmaj_index],mode='list')
		imhead(mosaicimage, mode="put", hdkey='beammajor', hdvalue=str(cleanmap_head['beammajor']['value']) + cleanmap_head['beammajor']['unit'])
		imhead(mosaicimage, mode="put", hdkey='beamminor', hdvalue=str(cleanmap_head['beamminor']['value']) + cleanmap_head['beamminor']['unit'])
		imhead(mosaicimage, mode="put", hdkey='beampa', hdvalue=str(cleanmap_head['beampa']['value']) + cleanmap_head['beampa']['unit'])
		imhead(mosaicimage, mode="put", hdkey='bunit', hdvalue=cleanmap_head['bunit'])

		## Exporting to FITS
		if enable_fitsout:
			exported_map= mosaicimage 
			print ('INFO: Exporting mosaic map %s to FITS...' % exported_map)
			exportfits(imagename=exported_map, fitsimage=fitsout, history=False, overwrite=True)

	
	t_stop = time.time()
	t_elapsed = t_stop - t_start

	## Print performance info
	print ('INFO: Imaging completed after dt(s)=%s' % str(t_elapsed))
	
###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())



