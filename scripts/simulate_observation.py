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
	
	# OPTIONAL OPTIONS
	parser.add_argument('-vis', '--vis', dest='vis', required=False, type=str, default='vis.ms', action='store',help='Output visibility CASA table name (default=vis.ms)')
	parser.add_argument('-skymodel', '--skymodel', dest='skymodel', required=False, type=str, default='skymodel.fits',action='store',help='Input FITS skymodel file (default=skymodel.fits')
	parser.add_argument('-outproject', '--outproject', dest='outproject', required=False, type=str, default='sim',action='store',help='Output simulation project name (default=sim)')
	
	parser.add_argument('-total_time', '--total_time', dest='total_time', required=False, type=int, default=43200,action='store',help='Total observation time in seconds (default=43200)')
	#parser.add_argument('-telconfigs', '--telconfigs', dest='telconfigs', required=False, type=list, default=['atca_all.cfg'],action='store',help='Telescope configuration list (default=atca_all.cfg)')
	#parser.add_argument('-telconfigs', '--telconfigs', dest='telconfigs', required=False, nargs='+', type=str, default=['atca_all.cfg'],action='store',help='Telescope configuration list (default=atca_all.cfg)')
	parser.add_argument('-telconfigs', '--telconfigs', dest='telconfigs', required=False, type=str, default='atca_all.cfg',action='store',help='Telescope configuration list (comma separated) (default=atca_all.cfg)')
	parser.add_argument('-mapsize', '--mapsize', dest='mapsize', required=False, type=str, default='',action='store',help='Map size (default=equal to skymodel')
	parser.add_argument('-pixsize', '--pixsize', dest='pixsize', required=False, type=str, default='',action='store',help='Pixel size (default=equal to skymodel')
	parser.add_argument('-addnoise','--addnoise', dest='addnoise', action='store_true')
	parser.set_defaults(addnoise=False)
	parser.add_argument('-frequency_center', '--frequency_center', dest='frequency_center', required=False, type=str, default='2.1GHz',action='store',help='Frequency center (default=2.1GHz)')
	parser.add_argument('-frequency_bandwidth', '--frequency_bandwidth', dest='frequency_bandwidth', required=False, type=str, default='10MHz',action='store',help='Frequency bandwidth (default=10MHz)')
	parser.add_argument('-maptype', '--maptype', dest='maptype', required=False, type=str, default='square',action='store',help='Map type {square|hexagonal} (default=square)')
	parser.add_argument('-obsmode', '--obsmode', dest='obsmode', required=False, type=str, default='int',action='store',help='Observation mode (default=int)')
	parser.add_argument('-graphics','--graphics', dest='enable_graphics', action='store_true')
	parser.add_argument('-no-graphics','--no-graphics', dest='enable_graphics', action='store_false')
	parser.set_defaults(enable_graphics=True)
	
	parser.add_argument('-c', dest='scriptname', required=False, type=str, default='',action='store',help='Script name')

	args = parser.parse_args()	

	return args

def mkdir_p(path):
	""" mkdir -p functionality """
	try:
		os.makedirs(path)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise


#### SIMULATE OBSERVATION ####
def simulate_observation(concat_vis='vis.ms',skymodel='skymodel.fits',exec_simobs_step=True,project_name='sim',total_time='43200s',telconfigs=['atca_6a.cfg','atca_6b.cfg','atca_ew352.cfg','atca_ew367.cfg'],obsmode='int',maptype='square',mapsize='',direction='',indirection='',incell='',incenter= '2.1GHz',inwidth='10MHz',integration='10s',imgaxes=['254.851041667','-41.4765888889','2.1GHz','I'],use_noise_vis=True,add_thermal_noise=False,graphics='both'):
	"""Simulate an observation from a sky model"""

	# Create output dir
	mkdir_p(project_name)

	## Strip filename and get basename
	skymodel_file_base, skymodel_file_ext= os.path.splitext(os.path.basename(skymodel))

	## Import sky model FITS file?
	##skymodel_img= project_name + '/' + str(skymodel_file_base) + '-casa'
	skymodel_img= project_name + '/skymodel'
	##skymodel_img= project_name + '.skymodel'
	print ('INFO: Importing sky model FITS file %s in CASA as image %s...' % (skymodel,skymodel_img))
	importfits(fitsimage=skymodel,imagename=skymodel_img,overwrite=True)
		
	## Generate simulated sky map visibilities
	print ('INFO: Generate simulated observation from sky model for given telescope configurations (n=%d configs)...' % len(telconfigs))
	vis_list= []	
	index = 0 
	for config in telconfigs:
		config_base, config_ext= os.path.splitext(os.path.basename(config)) 
		#current_obsmode= obsmode[index]

		if use_noise_vis:
			vis= project_name + '/' + project_name + '.' + config_base + '.noisy.ms'
			#vis= project_name + '.' + config_base + '.noisy.ms'
		else:
			vis= project_name + '/' + project_name + '.' + config_base + '.ms'
			#vis= project_name + '.' + config_base + '.ms'
		vis_list.append(vis) 

		# Execute the simulation step?
		if exec_simobs_step:

			# Add thermal noise?
			if add_thermal_noise:
				thermalnoise= 'tsys-atm'
			else:
				thermalnoise= ''

			# Simulate observation
			print ('INFO: Simulating observation for tel config %s (vis=%s)' % (str(config),str(vis)))
			simobserve(
				project=project_name, 
				skymodel=skymodel_img, 
				incenter=incenter, 
				incell=incell,
				inwidth=inwidth,
				indirection=indirection,
				direction=direction,	
				antennalist=config,
				sdantlist=config,
				totaltime=total_time,
				integration=integration,
				obsmode=obsmode,
				maptype=maptype,
				mapsize=mapsize,
				thermalnoise=thermalnoise,
				graphics=graphics,
				overwrite=True
			)
	
		## Update index
		index += 1

	## Concatenate visibilities
	print 'INFO: Concatenating visibilities: ', vis_list
	concat(vis=vis_list, concatvis=concat_vis)	


##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	## Get script args
	print 'INFO: #%d arguments found: %s' % (len(sys.argv),str(sys.argv))

	#if "-c" in sys.argv:
	#	index = sys.argv.index("-c")
	#	print 'index=',index

	try:
		args= get_args()
	except Exception as ex:
		print("ERROR: Failed to get and parse options (err=%s)",str(ex))
		return 1

	outproject= args.outproject
	vis= args.vis
	visout= outproject + '/' + vis
	
	skymodel= args.skymodel
	total_time= args.total_time
	t_tot= str(total_time) + 's'
	#telconfigs= args.telconfigs
	telconfigs = [str(item) for item in args.telconfigs.split(',')]

	mapsize= args.mapsize
	maptype= args.maptype
	obsmode= args.obsmode
	##obsmode= [str(item) for item in args.obsmode.split(',')]
	pixsize= args.pixsize
	frequency_center= args.frequency_center
	frequency_bandwidth= args.frequency_bandwidth
	enable_graphics= args.enable_graphics
	graphics= 'both'
	if not enable_graphics:
		graphics= 'none'
	addnoise= args.addnoise
	use_noise_vis= False
	add_thermal_noise= False
	if addnoise:
		use_noise_vis= True
		add_thermal_noise= True

	print("*** ARGS ***")
	print 'visout: ', visout
	print 'outproject: ', outproject
	print 'skymodel: ', skymodel
	print 'total_time: ', t_tot
	print 'obsmode: ', obsmode
	print 'maptype: ', maptype
	print 'mapsize: ', mapsize
	print 'pixsize: ', pixsize
	print 'telconfigs: ', telconfigs
	print 'addnoise? ', addnoise
	print 'use_noise_vis? ', use_noise_vis
	print 'add_thermal_noise? ', add_thermal_noise
	print 'frequency_center=',frequency_center
	print 'frequency_bandwidth=',frequency_bandwidth
	print 'graphics=',graphics
	print("************")
		
	## Simulate observation
	print("INFO: Starting simulation...")
	t_start = time.time()
	simulate_observation(
		project_name=outproject,
		concat_vis=visout,
		skymodel=skymodel,
		telconfigs=telconfigs,
		mapsize=mapsize,
		incell=pixsize,
		obsmode=obsmode,
		maptype=maptype,
		incenter=frequency_center,
		inwidth=frequency_bandwidth,
		use_noise_vis=use_noise_vis,
		add_thermal_noise=add_thermal_noise,
		total_time=t_tot,	
		graphics=graphics
	)
	t_stop = time.time()
	t_elapsed = t_stop - t_start

	## Print performance info
	print ('INFO: Simulation completed after dt(s)=%s' % str(t_elapsed))
	
###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())


