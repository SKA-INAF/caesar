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
import io

## ASTRO
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from astropy.io import ascii 
from astropy.coordinates import ICRS, Galactic, FK4, FK5

## COMMAND-LINE ARG MODULES
import getopt
import argparse
##################################################

#### GET SCRIPT ARGS ####
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")
	
	# - GENERAL IMAGE OPTIONS
	parser.add_argument('-inputfile', '--inputfile', dest='inputfile', required=True, type=str, action='store',help='Input Wise HII ascii catalog file')
	parser.add_argument('-outfile', '--outfile', dest='outfile', required=False, type=str, default='catalog.dat',action='store',help='Output filename')
	
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
	print('Get script args')
	try:
		args= get_args()
	except Exception as ex:
		print("Failed to get and parse options (err=%s)",str(ex))
		return 1

	input_file= args.inputfile
	output_file= args.outfile

	#===========================
	#==   READ TABLE
	#===========================
	row_start= 0
	table= ascii.read(input_file,data_start=row_start, delimiter='|')
	
	print(table.colnames)

	#===========================
	#==   MODIFY TABLE
	#===========================
	#table.remove_columns(['col1','col5','col6','col16','col17','col18','col19'])
	print(table)

	rowIndex= 0
	for data in table:
		ra_gal= data['Galactic Longitude']	
		dec_gal= data['Galactic Latitude']
		c= SkyCoord(ra_gal,dec_gal,frame='galactic', unit=(u.deg, u.deg))
		ra= c.fk5.ra.deg
		dec= c.fk5.dec.deg
		table[rowIndex]['Galactic Longitude']= str(ra)
		table[rowIndex]['Galactic Latitude']= str(dec)
		#print("Gal(%s,%s), J2000(%s,%s)" % (str(ra_gal),str(dec_gal),str(ra),str(dec)))
		rowIndex+= 1
		print("row=%d" % rowIndex)

	print(table)

	
	#===========================
	#==   SAVE TABLE
	#===========================	
	#buf= io.StringIO()
	buf= io.BytesIO()
	#ascii.write(table=table,output=output_file,quotechar=' ',format='no_header',delimiter='|',overwrite=True)
	ascii.write(table,buf,format='no_header',delimiter='|',overwrite=True)

	#- Remove astropy quotes around strings
	print('INFO: Saving table to file %s ' % output_file)
	with open(output_file,"w") as f:
		f.write(buf.getvalue().replace('"',""))

###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

