#!/usr/bin/env python

from __future__ import print_function

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

## ASTRO
from astropy.io import fits
from astropy.modeling.models import Box2D, Gaussian2D


## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## KERAS MODULES
import keras
from keras import layers
from keras import models
from keras import optimizers
from keras.utils import plot_model
from keras import backend as K
from keras.models import Model
from keras.layers.normalization import BatchNormalization
from keras.layers.convolutional import Conv2D
from keras.layers.convolutional import MaxPooling2D
from keras.layers.core import Activation
from keras.layers.core import Dropout
from keras.layers.core import Lambda
from keras.layers.core import Dense
from keras.layers import Flatten
from keras.layers import Input
import tensorflow as tf

## ADDON ML MODULES
from sklearn.model_selection import train_test_split

## GRAPHICS MODULES
import matplotlib.pyplot as plt



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

	## INPUT DATA
	parser.add_argument('-inputimg', '--inputimg', dest='inputimg', required=False, type=str,action='store',help='Mosaic residual image from which to extract train data')
	parser.add_argument('-filelist_bkg', '--filelist_bkg', dest='filelist_bkg', required=False, type=str,action='store',help='List of files with bkg train images')
	parser.add_argument('-filelist_source', '--filelist_source', dest='filelist_source', required=False, type=str,action='store',help='List of files with source train images')	
	parser.add_argument('-nnarcfile', '--nnarcfile', dest='nnarcfile', required=False, type=str,action='store',help='Name of file with NN architecture')	
	parser.add_argument('--generate_train_data', dest='generate_train_data', action='store_true',help='Generate train data randomly instead of reading train data from disk')	
	parser.set_defaults(generate_train_data=False)	
	
	parser.add_argument('--no-training', dest='no_training', action='store_true',help='Disable NN build and training')	
	parser.set_defaults(no_training=False)	

	parser.add_argument('-filelist_sourcepars', '--filelist_sourcepars', dest='filelist_sourcepars', required=False, type=str,action='store',help='List of files with source target pars')
	parser.add_argument('-marginx', '--marginx', dest='marginx', required=False, type=int, default=0,action='store',help='Input image x margin in pixels used in source generation')
	parser.add_argument('-marginy', '--marginy', dest='marginy', required=False, type=int, default=0,action='store',help='Input image y margin in pixels used in source generation')
	parser.add_argument('-marginx_source', '--marginx_source', dest='marginx_source', required=False, type=int, default=2,action='store',help='Train image x margin in pixels used in source generation')
	parser.add_argument('-marginy_source', '--marginy_source', dest='marginy_source', required=False, type=int, default=2,action='store',help='Train image y margin in pixels used in source generation')
	parser.add_argument('-Smin', '--Smin', dest='Smin', required=False, type=float, default=1.e-6, action='store',help='Minimum source flux in Jy (default=1.e-6)')
	parser.add_argument('-Smax', '--Smax', dest='Smax', required=False, type=float, default=1, action='store',help='Maximum source flux in Jy (default=1)')
	parser.add_argument('-Smodel', '--Smodel', dest='Smodel', required=False, type=str, default='uniform', action='store',help='Source flux generation model (default=uniform)')
	parser.add_argument('-Sslope', '--Sslope', dest='Sslope', required=False, type=float, default=1.6, action='store',help='Slope par in expo source flux generation model (default=1.6)')
	parser.add_argument('-bmaj_min', '--bmaj_min', dest='bmaj_min', required=False, type=float, default=4, action='store',help='Gaussian components min bmaj in arcsec (default=4)')
	parser.add_argument('-bmaj_max', '--bmaj_max', dest='bmaj_max', required=False, type=float, default=10, action='store',help='Gaussian components max bmaj in arcsec (default=10)')
	parser.add_argument('-bmin_min', '--bmin_min', dest='bmin_min', required=False, type=float, default=4, action='store',help='Gaussian components  min bmin in arcsec (default=4)')
	parser.add_argument('-bmin_max', '--bmin_max', dest='bmin_max', required=False, type=float, default=10, action='store',help='Gaussian components  max bmin in arcsec (default=10)')
	parser.add_argument('-pa_min', '--pa_min', dest='pa_min', required=False, type=float, default=-90, action='store',help='Gaussian components  min position angle in deg (default=0)')
	parser.add_argument('-pa_max', '--pa_max', dest='pa_max', required=False, type=float, default=90, action='store',help='Gaussian components  max position angle in deg (default=180)')
	parser.add_argument('-nsources_max', '--nsources_max', dest='nsources_max', required=False, type=int, default=5, action='store',help='Maximum number of sources generated per crop image (default=5)')
	
	parser.add_argument('--generate_bkg', dest='generate_bkg', action='store_true',help='Generate bkg instead of reading it from input image')	
	parser.set_defaults(generate_bkg=False)	
	parser.add_argument('-bkg_rms', '--bkg_rms', dest='bkg_rms', required=False, type=float, default=300.e-6, action='store',help='Generated bkg rms (default=300 muJy/beam)')
	parser.add_argument('-bkg_mean', '--bkg_mean', dest='bkg_mean', required=False, type=float, default=0, action='store',help='Generated bkg average (default=0 muJy/beam)')


	parser.add_argument('-normdatamin', '--normdatamin', dest='normdatamin', required=False, type=float, default=-0.0100, action='store',help='Normalization min used to scale data in (0,1) range (default=-100 mJy/beam)')	
	parser.add_argument('-normdatamax', '--normdatamax', dest='normdatamax', required=False, type=float, default=10, action='store',help='Normalization max used to scale data in (0,1) range (default=10 Jy/beam)')

	parser.add_argument('--normalize_targets', dest='normalize_targets', action='store_true',help='Normalize target data before training')	
	parser.set_defaults(normalize_targets=False)
	
	parser.add_argument('-nx', '--nx', dest='nx', required=False, type=int, default=101, action='store',help='Image width in pixels (default=101)')
	parser.add_argument('-ny', '--ny', dest='ny', required=False, type=int, default=101, action='store',help='Image height in pixels (default=101)')	
	parser.add_argument('-nsamples_bkg', '--nsamples_bkg', dest='nsamples_bkg', required=False, type=int, default=10, action='store',help='Number of train images for bkg extracted from input maps (default=10)')
	parser.add_argument('-nsamples_source', '--nsamples_source', dest='nsamples_source', required=False, type=int, default=-1, action='store',help='Number of train images extracted around sources from input maps (default=-1)')	
	parser.add_argument('-nmaxobjects', '--nmaxobjects', dest='nmaxobjects', required=False, type=int, default=5, action='store',help='Max number of predicted objects in target (default=5)')
	parser.add_argument('-ntargetpars', '--ntargetpars', dest='ntargetpars', required=False, type=int, default=6, action='store',help='Nmber of pars per objects in target (default=6)')	

	parser.add_argument('-conv_kern_size_min', '--conv_kern_size_min', dest='conv_kern_size_min', required=False, type=int, default=3, action='store',help='Min size of conv kernel (default=3)')	
	parser.add_argument('-conv_kern_size_max', '--conv_kern_size_max', dest='conv_kern_size_max', required=False, type=int, default=7, action='store',help='Max size of conv kernel (default=3)')	

	parser.add_argument('-conv_nfilt_min', '--conv_nfilt_min', dest='conv_nfilt_min', required=False, type=int, default=16, action='store',help='Number of min convolution filters used (default=16)')	
	parser.add_argument('-conv_nfilt_max', '--conv_nfilt_max', dest='conv_nfilt_max', required=False, type=int, default=32, action='store',help='Number of max convolution filters used (default=32)')
	parser.add_argument('-conv_pool_kernel_size', '--conv_pool_kernel_size', dest='conv_pool_kernel_size', required=False, type=int, default=2, action='store',help='Size of max pool kernel in conv layers (default=2)')

	parser.add_argument('-dense_size_min', '--dense_size_min', dest='dense_size_min', required=False, type=int, default=16, action='store',help='Number of min neurons used in dense layer(default=16)')
	parser.add_argument('-dense_size_max', '--dense_size_max', dest='dense_size_max', required=False, type=int, default=32, action='store',help='Number of max neurons used in dense layer(default=32)')

	parser.add_argument('-conv_activation', '--conv_activation', dest='conv_activation', required=False, type=str, default='relu', action='store',help='Activation used in convolution layers (default=relu)')
	parser.add_argument('-dense_activation', '--dense_activation', dest='dense_activation', required=False, type=str, default='relu', action='store',help='Activation used in dense layers (default=relu)')

	parser.add_argument('-optimizer', '--optimizer', dest='optimizer', required=False, type=str, default='rmsprop', action='store',help='Optimizer used (default=rmsprop)')
	parser.add_argument('-learning_rate', '--learning_rate', dest='learning_rate', required=False, type=float, default=1.e-4, action='store',help='Learning rate (default=1.e-4)')

	parser.add_argument('--use_standard_nn', dest='use_standard_nn', action='store_true')	
	parser.set_defaults(use_standard_nn=False)	

	parser.add_argument('--no-conv', dest='no_conv', action='store_true')	
	parser.set_defaults(no_conv=False)	
	parser.add_argument('--no-dropout', dest='no_dropout', action='store_true')	
	parser.set_defaults(no_dropout=False)	
	parser.add_argument('--no-maxpool', dest='no_maxpool', action='store_true')	
	parser.set_defaults(no_maxpool=False)
	parser.add_argument('--no-batchnorm', dest='no_batchnorm', action='store_true')	
	parser.set_defaults(no_batchnorm=False)

	parser.add_argument('-test_size', '--test_size', dest='test_size', required=False, type=float, default=0.2, action='store',help='Fraction of input data used for testing the network (default=0.2)')
	parser.add_argument('-spars_loss_weight', '--spars_loss_weight', dest='spars_loss_weight', required=False, type=float, default=1, action='store',help='Loss weight to be given to source pars learning (default=1)')
	parser.add_argument('-labels_loss_weight', '--labels_loss_weight', dest='labels_loss_weight', required=False, type=float, default=1, action='store',help='Loss weight to be given to source labels learning (default=1)')
	parser.add_argument('-nepochs', '--nepochs', dest='nepochs', required=False, type=int, default=100, action='store',help='Number of epochs used in network training (default=100)')
	parser.add_argument('--saveimg', dest='saveimg', action='store_true')	
	parser.set_defaults(saveimg=False)

	parser.add_argument('--flip_train', dest='flip_train', action='store_true')	
	parser.set_defaults(flip_train=False)
	parser.add_argument('--flip_test', dest='flip_test', action='store_true')	
	parser.set_defaults(flip_test=False)

	parser.add_argument('-outfile_loss', '--outfile_loss', dest='outfile_loss', required=False, type=str, default='nn_loss.png', action='store',help='Name of NN loss plot file (default=nn_loss.png)')
	parser.add_argument('-outfile_accuracy', '--outfile_accuracy', dest='outfile_accuracy', required=False, type=str, default='nn_accuracy.png', action='store',help='Name of NN accuracy plot file (default=nn_accuracy.png)')
	parser.add_argument('-outfile_model', '--outfile_model', dest='outfile_model', required=False, type=str, default='nn_model.png', action='store',help='Name of NN model plot file (default=nn_model.png)')
	parser.add_argument('-outfile_posaccuracy', '--outfile_posaccuracy', dest='outfile_posaccuracy', required=False, type=str, default='nn_posaccuracy.png', action='store',help='Name of NN source position accuracy plot file (default=nn_posaccuracy.png)')
	parser.add_argument('-outfile_fluxaccuracy', '--outfile_fluxaccuracy', dest='outfile_fluxaccuracy', required=False, type=str, default='nn_fluxaccuracy.png', action='store',help='Name of NN source flux accuracy plot file (default=nn_fluxaccuracy.png)')

	parser.add_argument('-outfile_nnout_train', '--outfile_nnout_train', dest='outfile_nnout_train', required=False, type=str, default='train_nnout.dat', action='store',help='Name of output file with NN output for train data (default=train_nnout.dat)')
	parser.add_argument('-outfile_nnout_test', '--outfile_nnout_test', dest='outfile_nnout_test', required=False, type=str, default='test_nnout.dat', action='store',help='Name of output file with NN output for test data (default=test_nnout.dat)')
	parser.add_argument('-outfile_nnout_metrics', '--outfile_nnout_metrics', dest='outfile_nnout_metrics', required=False, type=str, default='nnout_metrics.dat', action='store',help='Name of output file with NN train metrics (default=nnout_metrics.dat)')

	args = parser.parse_args()	

	return args


	
###########################
##     NN PRINT CLASS
###########################
class NNPrinter(keras.callbacks.Callback):

	def on_train_begin(self, logs={}):
		self.losses = []

	def on_batch_end(self, batch, logs={}):
		#y_true= self.model.outputs
		#print("type(y_true)=",type(y_true))
		#print("len(y_true)=",len(y_true))
		#y_true_shape= K.int_shape(y_true)
		#print("y_true shape=",y_true_shape)
		
		loss= logs.get('loss')
		#print("LOSS=", loss)
		self.losses.append(loss)

###########################
##     NN TRAINER CLASS
###########################
SIGMA_TO_FWHM= np.sqrt(8*np.log(2))

class CNNTrainer(object):

	""" Network trainer class

			Attributes:
				residual: beam residual image (.fits)
				model: skymodel image (.fits)
				img_file: beam restored image (.fits)
	"""

	def __init__(self):
		""" Return a CNNTrainer object """

		# - Input file
		self.img_file= None 
		self.img_data= None
		self.img_sizex= 0
		self.img_sizey= 0
		self.pixsize= 4 # in arcsec
		self.img_bkg_filelist= ''
		self.img_source_filelist= ''
		self.sourcepars_filelist= ''
		self.nnarc_file= ''

		# - Bkg generation
		self.gen_bkg_from_img= True
		self.bkg_rms= 300.e-6
		self.bkg_mean= 0
		
		# - Source generation
		self.gridx= None
		self.gridy= None
		self.nsources_max= 5
		self.source_gen_marginx= 2
		self.source_gen_marginy= 2
		self.Smin= 1.e-6 # in Jy 
		self.Smax= 1 # in Jy
		self.Smodel= 'uniform'
		self.Sslope= 1.6
		self.truncate_models= False
		self.trunc_thr= 0.01 # 1% flux truncation at maximum

		# - Train data
		self.npars= 6 # Number of parameters per object
		self.nobjects= 5 # Max number of objects
		self.nsamples_bkg= 10
		self.nsamples_source= 10
		self.train_img_sizex= 101
		self.train_img_sizey= 101
		self.normalize_inputs= True
		self.normmin= 0.001
		self.normmax= 10
		self.train_data_gen_enabled= True
		self.inputs_bkg= None		
		self.outputs_bkg= None
		self.outputs_labels_bkg= None
		self.inputs_source= None		
		self.outputs_source= None
		self.outputs_labels_source= None
		self.normalize_targets= False
		self.theta_min= -90
		self.theta_max= 90
		self.sigma_min= 0
		self.sigma_max= 20
		self.normmin_pars= np.array([0,0,self.normmin,self.sigma_min,self.sigma_min,np.radians(self.theta_min)])
		self.normmax_pars= np.array([self.train_img_sizex,self.train_img_sizey,self.normmax,self.sigma_max,self.sigma_max,np.radians(self.theta_max)])
		self.test_size= 0.2
		self.inputs_train= None
		self.inputs_test= None 
		self.outputs_train= None
		self.outputs_test= None
		self.outputs_labels_train= None 
		self.outputs_labels_test= None

		self.flipped_outputs_labels_train= None
		self.flipped_outputs_train= None
		self.flipped_outputs_labels_test= None
		self.flipped_outputs_test= None
		
		# - Network architecture & train options
		self.do_training= True
		self.model= None
		self.fitsout= None
		self.standard_nn_build_enabled= False
		self.nepochs= 10
		self.dropout= 0.25
		self.dropout_dense= 0.5
		self.padding= "same"
		#self.padding= "valid"	
		self.conv_nfilt_min= 16
		self.conv_nfilt_max= 32
		self.conv_kern_size_min= 3
		self.conv_kern_size_max= 3
		self.conv_act= 'relu'
		self.pool_size= 2
		self.dense_size_min= 16
		self.dense_size_max= 32
		self.dense_act= 'relu'
		self.flip_test_data= True
		self.flip_train_data= True
		self.spars_loss_weight= 1
		self.labels_loss_weight= 1
		self.conv_enabled= True
		self.dropout_enabled= True
		self.maxpool_enabled= True
		self.batchnorm_enabled= True
		self.train_loss_vs_epoch= None
		self.test_loss_vs_epoch= None
		self.train_accuracy_vs_epoch= None
		self.test_accuracy_vs_epoch= None	
		self.optimizer= 'rmsprop'
		self.learning_rate= 1.e-4

		# - NN callback
		self.nnprinter_cb= NNPrinter()

		# - Output file options
		self.writeimg= False
		self.outfile_loss= 'nn_loss.png'
		self.outfile_accuracy= 'nn_accuracy.png'
		self.outfile_model= 'nn_model.png'
		self.outfile_posaccuracy= 'nn_posaccuracy.png'
		self.outfile_fluxaccuracy= 'nn_fluxaccuracy.png'
		self.outfile_nnout_train= 'train_nnout.dat'
		self.outfile_nnout_test= 'test_nnout.dat'
		self.outfile_nnout_metrics= 'nnout_metrics.dat'

	#################################
	##     SETTER/GETTER METHODS
	#################################
	def enable_bkg_generation_from_img(self,choice):
		""" Turn on/off bkg generation from input image. """
		self.gen_bkg_from_img= choice

	def set_gen_bkg_rms(self,rms):
		""" Set generated bkg rms """
		self.bkg_rms= rms

	def set_gen_bkg_mean(self,mean):
		""" Set generated bkg mean"""
		self.bkg_mean= mean

	def set_img_filename(self,filename):
		""" Set the input residual image used to generate train data """
		self.img_file= filename

	def set_img_bkg_filelist(self,filename):
		""" Set the name of filelist with bkg train images """
		self.img_bkg_filelist= filename

	def set_img_source_filelist(self,filename):
		""" Set the name of filelist with source train images """
		self.img_source_filelist= filename
	
	def set_sourcepars_filelist(self,filename):
		""" Set the name of filelist with source target pars """
		self.sourcepars_filelist= filename

	def set_nnarc_filename(self,filename):
		""" Set NN architecture filename """
		self.nnarc_file= filename

	def set_outfile_loss(self,filename):
		""" Set output file name for loss plot """
		self.outfile_loss= filename

	def set_outfile_accuracy(self,filename):
		""" Set output file name for accuracy plot """
		self.outfile_accuracy= filename

	def set_outfile_model(self,filename):
		""" Set output file name for model plot """
		self.outfile_model= filename

	def set_outfile_posaccuracy(self,filename):
		""" Set output file name for pos accuracy plot """
		self.outfile_posaccuracy= filename

	def set_outfile_fluxaccuracy(self,filename):
		""" Set output file name for flux accuracy plot """
		self.outfile_fluxaccuracy= filename

	def set_outfile_nnout_train(self,filename):	
		""" Set output file name where to store NN output for train data"""
		self.outfile_nnout_train= filename

	def set_outfile_nnout_test(self,filename):	
		""" Set output file name where to store NN output for test data"""	
		self.outfile_nnout_test= filename

	def set_outfile_nnout_metrics(self,filename):	
		""" Set output file name where to store NN output metrics"""	
		self.outfile_nnout_metrics= filename

	def enable_training(self,choice):
		""" Enable/disable training step"""
		self.do_training= choice

	def set_margins(self,marginx,marginy):
		""" Set margin in X & Y """
		self.marginx= marginx
		self.marginy= marginy

	def set_source_margins(self,marginx,marginy):
		""" Set margin in X & Y for source generation """
		self.source_gen_marginx= marginx
		self.source_gen_marginy= marginy

	def set_input_data_norm_range(self,datamin,datamax):
		""" Set input data normalization range """
		self.normmin= datamin
		self.normmax= datamax
		self.normmin_pars= np.array([0,0,self.normmin,self.sigma_min,self.sigma_min,np.radians(self.theta_min)])
		self.normmax_pars= np.array([self.train_img_sizex,self.train_img_sizey,self.normmax,self.sigma_max,self.sigma_max,np.radians(self.theta_max)])

	def set_nobjects(self,n):
		""" Set maximum number of detected object in image """
		self.nobjects= n

	def set_npars(self,n):
		""" Set number of source parameters to be fitted in model """
		self.npars= n

	def set_bkg_sample_size(self,n):
		""" Set number of images for bkg used in training """
		self.nsamples_bkg= n

	def set_test_sample_size(self,f):
		""" Set test sample proportion """
		self.test_size= f

	def set_source_sample_size(self,n):
		""" Set number of images for sources used in training """
		self.nsamples_source= n

	def set_train_img_size(self,nx,ny):
		""" Set size of input image given to the network for training """
		self.train_img_sizex= nx
		self.train_img_sizey= ny
		self.normmin_pars= np.array([0,0,self.normmin,self.sigma_min,self.sigma_min,np.radians(self.theta_min)])
		self.normmax_pars= np.array([self.train_img_sizex,self.train_img_sizey,self.normmax,self.sigma_max,self.sigma_max,np.radians(self.theta_max)])

	def use_standard_nn(self,choice):	
		""" Use standard nn architecture instead of reading from file """
		self.standard_nn_build_enabled= choice

	def set_dense_activation(self,act):
		""" Set dense layer activation fucntion """
		self.dense_act= act

	def set_conv_activation(self,act):
		""" Set convolution layer activation fucntion """
		self.conv_act= act

	def set_conv_pool_kernel_size(self,n):
		""" Set size of max pool kernel in conv layer """
		self.pool_size= n	

	def set_conv_kern_size_min(self,n):
		""" Set min size of conv kernel in conv layer """
		self.conv_kern_size_min= n		
	
	def set_conv_kern_size_max(self,n):
		""" Set max size of conv kernel in conv layer """
		self.conv_kern_size_max= n	

	def set_conv_nfilt_min(self,n):
		""" Set min number of filters in conv layer """
		self.conv_nfilt_min= n		
		
	def set_conv_nfilt_max(self,n):
		""" Set max number of filters in conv layer """
		self.conv_nfilt_max= n

	def set_dense_size_min(self,n):
		""" Set min number of neurons in dense layer """
		self.dense_size_min= n

	def set_dense_size_max(self,n):
		""" Set max number of neurons in dense layer """
		self.dense_size_max= n

	def enable_conv(self,choice):
		""" Turn on/off convolution layers """
		self.conv_enabled= choice		

	def enable_dropout(self,choice):
		""" Turn on/off dropout layers """
		self.dropout_enabled= choice	

	def enable_maxpool(self,choice):
		""" Turn on/off maxpool layers """
		self.maxpool_enabled= choice		

	def enable_batchnorm(self,choice):
		""" Turn on/off batch normalization layers """
		self.batchnorm_enabled= choice		

	def enable_target_normalization(self,choice):
		""" Turn on/off target normalization """
		self.normalize_targets= choice

	def set_spars_loss_weight(self,w):
		""" Set source par loss weight """
		self.spars_loss_weight= w

	def set_labels_loss_weight(self,w):
		""" Set source labels loss weight """
		self.labels_loss_weight= w

	def set_nepochs(self,w):
		""" Set number of train epochs """
		self.nepochs= w

	def set_optimizer(self,opt):
		""" Set optimizer """
		self.optimizer= opt

	def set_learning_rate(self,lr):
		""" Set learning rate """
		self.learning_rate= lr

	def save_train_img_to_disk(self,choice):
		""" Turn on/off writing to disk of train images for bkg and source """
		self.writeimg= choice

	def enable_train_data_generation(self,choice):
		""" Turn on/off generation of train data """
		self.train_data_gen_enabled= choice

	def enable_train_data_flip(self,choice):
		""" Turn on/off flipping of train data during training """
		self.flip_train_data= choice

	def enable_test_data_flip(self,choice):
		""" Turn on/off flipping of test data during training """
		self.flip_test_data= choice

	def set_nsources_max(self,n):
		""" Set the maximum number of sources to be generated in train image """
		self.nsources_max= n

	def set_source_flux_rand_model(self,model):
		""" Set the source flux random model """
		self.Smodel= model

	def set_source_flux_rand_exp_slope(self,slope):
		""" Set the source flux expo model slope par """
		self.Sslope= slope	
	
	def set_source_flux_range(self,Smin,Smax):
		""" Set source flux range """
		self.Smin= Smin
		self.Smax= Smax	

	def set_beam_bmaj_range(self,bmaj_min,bmaj_max):
		""" Set beam bmaj range """
		self.beam_bmaj_min= bmaj_min
		self.beam_bmaj_max= bmaj_max	

	def set_beam_bmin_range(self,bmin_min,bmin_max):
		""" Set beam bmin range """
		self.beam_bmin_min= bmin_min
		self.beam_bmin_max= bmin_max	

	def set_beam_pa_range(self,pa_min,pa_max):
		""" Set beam pa range """
		self.beam_bpa_min= pa_min
		self.beam_bpa_max= pa_max	
	
	
	#################################
	##     HELPER METHODS
	#################################
	def has_patterns_in_string(self,s,patterns):
		""" Return true if patterns are found in string """
		if not patterns:		
			return False

		found= False
		for pattern in patterns:
			found= pattern in s
			if found:
				break

		return found

	#################################
	##     WRITE DATA TO ASCII FILE
	#################################
	def write_ascii(self,data,filename,header=''):
		""" Write data to ascii file """
		# - Skip if data is empty
		if data.size<=0:
			return

		# - Open file and write header
		fout = open(filename, 'wt')
		if header:
			fout.write(header)
			fout.write('\n')	
			fout.flush()	
		
		# - Write t
		nrows= data.shape[0]
		ncols= data.shape[1]
		for i in range(nrows):
			fields= '  '.join(map(str, data[i,:]))
			fout.write(fields)
			fout.write('\n')	
			fout.flush()	

		fout.close();

	#################################
	##     READ ASCII DATA
	#################################
	def read_ascii(self,filename,skip_patterns=[]):
		""" Read an ascii file line by line """
	
		try:
			f = open(filename, 'r')
		except IOError:
			errmsg= 'Could not read file: ' + filename
			print ("ERROR: " + errmsg)
			raise IOError(errmsg)

		fields= []
		for line in f:
			line = line.strip()
			line_fields = line.split()
			if not line_fields:
				continue

			# Skip pattern
			skipline= self.has_patterns_in_string(line_fields[0],skip_patterns)
			if skipline:
				continue 		

			fields.append(line_fields)

		f.close()	

		return fields

	#################################
	##     WRITE IMAGE DATA TO FITS
	#################################
	def write_fits(self,data,filename):
		""" Read data to FITS image """
		hdu= fits.PrimaryHDU(data)
		hdul= fits.HDUList([hdu])
		hdul.writeto(filename,overwrite=True)

	
	#################################
	##     READ FITS IMAGE
	#################################
	def read_fits(self,filename):
		""" Read FITS image and return data """

		try:
			hdu= fits.open(filename,memmap=False)
		except Exception as ex:
			errmsg= 'Cannot read image file: ' + filename
			print ("ERROR: " + errmsg)
			raise IOError(errmsg)

		data= hdu[0].data
		data_size= np.shape(data)
		nchan= len(data.shape)
		if nchan==4:
			output_data= data[0,0,:,:]
		elif nchan==2:
			output_data= data	
		else:
			errmsg= 'Invalid/unsupported number of channels found in file ' + filename + ' (nchan=' + str(nchan) + ')!'
			print ("ERROR: " + errmsg)
			hdu.close()
			raise IOError(errmsg)

		hdu.close()

		return output_data


	#################################
	##     READ INPUT IMAGE
	#################################
	def read_img(self):
		""" Read FITS image and set image data """

		# - Read FITS image
		try:
			hdu= fits.open(self.img_file,memmap=False)
		except Exception as ex:
			errmsg= 'Cannot read image file: ' + self.img_file
			print ("ERROR: " + errmsg)
			return -1

		data= hdu[0].data
		data_size= np.shape(data)
		nchan= len(data.shape)
		if nchan==4:
			self.img_data= data[0,0,:,:]
		elif nchan==2:
			self.img_data= data	
		else:
			errmsg= 'Invalid/unsupported number of channels found in file ' + self.img_file + ' (nchan=' + str(nchan) + ')!'
			print ("ERROR: " + errmsg)
			hdu.close()
			return -1

		imgsize= np.shape(self.img_data)
		self.img_sizex= imgsize[1]
		self.img_sizey= imgsize[0]

		# - Read fits metadata
		header= hdu[0].header
		dx= np.abs(header['CDELT1']*3600.) # in arcsec
		dy= np.abs(header['CDELT2']*3600.) # in arcsec
		self.pixsize= min(dx,dy)
			
		hdu.close()

		return 0

	#################################
	##     CROP IMAGE
	#################################
	def crop_img(self,x0,y0,dx,dy):
		""" Extract sub image of size (dx,dy) around pixel (x0,y0) """

		#- Extract crop data
		xmin= int(x0-dx/2)
		xmax= int(x0+dx/2)
		ymin= int(y0-dy/2)
		ymax= int(y0+dy/2)		
		crop_data= self.img_data[ymin:ymax,xmin:xmax]
	
		#- Replace NAN with zeros and inf with large numbers
		np.nan_to_num(crop_data,False)

		return crop_data


	#################################
	##     READ BKG TRAIN DATA
	#################################
	def read_bkg_train_data(self,filelist):
		""" Read background train data """
		
		# - Init data
		input_data= []	
		output_size= self.nobjects*self.npars
		output_data= []  
		output_label_size= self.nobjects
		output_label_data= []
		nchannels= 1
		filelist_data= []

		# - Read list with files		
		try:
			filelist_data= self.read_ascii(filelist,['#'])
		except IOError:
			errmsg= 'Cannot read file: ' + filelist
			print ("ERROR: " + errmsg)
			return -1

		# - Read image files in list	
		imgcounter= 0
		for item in filelist_data:
			imgcounter+= 1
			filename= item[0]
			print("INFO: Reading file %s ..." % filename) 

			data= None
			try:
				data= self.read_fits(filename)
			except Exception as ex:
				errmsg= 'Failed to read bkg image data (err=' + str(ex) + ')'
				print ("ERROR: " + errmsg)
				return -1
	
			imgsize= np.shape(data)
			nx= imgsize[1]
			ny= imgsize[0]	
			#print("INFO: Bkg image no. %d has size (%d,%d)" % (imgcounter,nx,ny) )	

			# - Check bkg image size is equal to desired train image
			if nx!=self.train_img_sizex or ny!=self.train_img_sizey:
				errmsg= 'Bkg image no. ' + str(imgcounter) + ' has size different from desired train image!'
				print ("ERROR: " + errmsg)
				return -1

			# - Set train data as a tensor of size [Nsamples,Nx,Ny,Nchan] Nchan=1
			data= data.reshape(imgsize[0],imgsize[1],nchannels)
			input_data.append(data)

			# - Set train target & labels
			output_data.append( np.zeros((1,output_size)) )
			output_label_data.append( np.zeros((1,output_label_size)) )

		#- Convert list to array
		self.inputs_bkg= np.array(input_data)
		self.inputs_bkg= self.inputs_bkg.astype('float32')
	
		self.outputs_bkg= np.array(output_data)
		self.outputs_bkg= self.outputs_bkg.astype('float32')

		outputs_shape= self.outputs_bkg.shape
		N= outputs_shape[0]
		self.outputs_bkg= self.outputs_bkg.reshape((N,output_size))

		self.outputs_labels_bkg= np.array(output_label_data)
		self.outputs_labels_bkg= self.outputs_labels_bkg.astype('float32')
		self.outputs_labels_bkg= self.outputs_labels_bkg.reshape((N,output_label_size))

		# - Normalize to [0,1]
		if self.normalize_inputs:
			print("DEBUG: inputs_bkg (BEFORE NORMALIZATION): min/max=%s/%s" % (str(np.min(self.inputs_bkg)),str(np.max(self.inputs_bkg))))
			self.inputs_bkg= (self.inputs_bkg - self.normmin)/(self.normmax-self.normmin)
			print("DEBUG: inputs_bkg (AFTER NORMALIZATION): min/max=%s/%s" % (str(np.min(self.inputs_bkg)),str(np.max(self.inputs_bkg))))
		

		# - Normalize targets to [0,1]
		#if self.normalize_targets:	
		#	targets_normmin= np.zeros(self.nobjects*self.npars)
		#	targets_normmax= np.zeros(self.nobjects*self.npars)
		#	par_counter= 0
		#	for k in range(self.nobjects):
		#		for l in range(self.npars):
		#			targets_normmin[par_counter]= self.normmin_pars[l]
		#			targets_normmax[par_counter]= self.normmax_pars[l]
		#			par_counter+= 1

		#	print("DEBUG: targets_normmin=", targets_normmin)
		#	print("DEBUG: targets_normmax=", targets_normmax)
		#	print("DEBUG: outputs_bkg (BEFORE NORMALIZATION): min/max=%s/%s" % (str(np.min(self.outputs_bkg)),str(np.max(self.outputs_bkg))))
		#	self.outputs_bkg= (self.outputs_bkg - targets_normmin)/(targets_normmax-targets_normmin)
		#	print("DEBUG: outputs_bkg (AFTER NORMALIZATION): min/max=%s/%s" % (str(np.min(self.outputs_bkg)),str(np.max(self.outputs_bkg))))

		print("DEBUG: outputs_bkg: min/max=%s/%s" % (str(np.min(self.outputs_bkg)),str(np.max(self.outputs_bkg))))
		
		print("DEBUG: inputs_bkg size=", np.shape(self.inputs_bkg))
		print("DEBUG: outputs_bkg size=", np.shape(self.outputs_bkg))
		print("DEBUG: outputs_labels_bkg size=", np.shape(self.outputs_labels_bkg))
		print("DEBUG: outputs_bkg=",self.outputs_bkg)
		print("DEBUG: outputs_labels_bkg=",self.outputs_labels_bkg)

		return 0


	#################################
	##     READ SOURCE TRAIN DATA
	#################################
	def read_source_train_data(self,filelist,filelist_pars):
		""" Read source train data """
				
		# - Init data
		input_data= []
		output_size= self.nobjects*self.npars
		output_data= []  
		output_label_size= self.nobjects
		output_label_data= []
		nchannels= 1
		filelist_data= []
		filelist_pars_data= []

		# - Read list with image files		
		try:
			filelist_data= self.read_ascii(filelist,['#'])
		except IOError:
			errmsg= 'Cannot read file: ' + filelist
			print ("ERROR: " + errmsg)
			return -1

		# - Read list with source pars files		
		try:
			filelist_pars_data= self.read_ascii(filelist_pars,['#'])
		except IOError:
			errmsg= 'Cannot read file: ' + filelist_pars
			print ("ERROR: " + errmsg)
			return -1

		# - Check lists have the same number of entries
		if len(filelist_data)!=len(filelist_pars_data):
			print("ERROR: Source img and pars filelist have different number of entries (%s!=%s)" % (len(filelist_data),len(filelist_pars_data)))
			return -1

		# - Read source images & pars
		imgcounter= 0

		for item, item_pars in zip(filelist_data,filelist_pars_data):
			imgcounter+= 1
			filename= item[0]
			filename_pars= item_pars[0]
			print("INFO: Reading files: %s, %s ..." % (filename,filename_pars) )

			# - Read source img
			try:
				data= self.read_fits(filename)
			except Exception as ex:
				errmsg= 'Failed to read source image data (err=' + str(ex) + ')'
				print ("ERROR: " + errmsg)
				return -1
	
			imgsize= np.shape(data)
			nx= imgsize[1]
			ny= imgsize[0]
			#print("INFO: Source image no. %d has size (%d,%d)" % (imgcounter,nx,ny) )

			# - Check source image size is equal to desired train image
			if nx!=self.train_img_sizex or ny!=self.train_img_sizey:
				errmsg= 'Source image no. ' + str(imgcounter) + ' has size different from desired train image!'
				print ("ERROR: " + errmsg)
				return -1

			# - Set train data as a tensor of size [Nsamples,Nx,Ny,Nchan] Nchan=1
			data= data.reshape(imgsize[0],imgsize[1],nchannels)
			input_data.append(data)

			# - Read source pars
			source_pars= []
			skip_patterns= ['#']
			try:
				source_pars= self.read_ascii(filename_pars,skip_patterns)
			except IOError:
				errmsg= 'Cannot read file: ' + filename_spar
				print ("ERROR: " + errmsg)
				return -1

			source_pars_size= np.shape(source_pars)
			nsources= source_pars_size[0]
			npars= source_pars_size[1]

			# - Check source pars number is >= desired pars
			if npars<self.npars:
				print("ERROR: Source pars read from file no. %s (%d) smaller than desired number of source pars (%d)" % (str(imgcounter),npars,self.npars) )
				return -1

			# - Set train targets
			targets= np.zeros((1,output_size))
			target_labels= np.zeros((1,output_label_size))
			par_counter= 0

			for k in range(nsources):
				target_labels[0,k]= 1
				for l in range(self.npars):			
					targets[0,par_counter+l]= source_pars[k][1+l]	
				par_counter+= self.npars

			output_data.append(targets)
			output_label_data.append(target_labels)


		#- Convert list to array
		self.inputs_source= np.array(input_data)
		self.inputs_source= self.inputs_source.astype('float32')

		self.outputs_source= np.array(output_data)
		self.outputs_source= self.outputs_source.astype('float32')

		outputs_shape= self.outputs_source.shape
		N= outputs_shape[0]
		self.outputs_source= self.outputs_source.reshape((N,output_size))

		self.outputs_labels_source= np.array(output_label_data)
		self.outputs_labels_source= self.outputs_labels_source.astype('float32')
		self.outputs_labels_source= self.outputs_labels_source.reshape((N,output_label_size))

		# - Normalize to [0,1]
		if self.normalize_inputs:
			print("DEBUG: inputs_source (BEFORE NORMALIZATION): min/max=%s/%s" % (str(np.min(self.inputs_source)),str(np.max(self.inputs_source))))
			self.inputs_source= (self.inputs_source - self.normmin)/(self.normmax-self.normmin)
			print("DEBUG: inputs_source (AFTER NORMALIZATION): min/max=%s/%s" % (str(np.min(self.inputs_source)),str(np.max(self.inputs_source))))
			
		# - Normalize targets to [0,1]
		if self.normalize_targets:
			
			targets_normmin= np.zeros(self.nobjects*self.npars)
			targets_normmax= np.zeros(self.nobjects*self.npars)
			par_counter= 0
			for k in range(self.nobjects):
				for l in range(self.npars):
					targets_normmin[par_counter]= self.normmin_pars[l]
					targets_normmax[par_counter]= self.normmax_pars[l]
					par_counter+= 1

			print("DEBUG: targets_normmin=", targets_normmin)
			print("DEBUG: targets_normmax=", targets_normmax)
			print("DEBUG: outputs_source (BEFORE NORMALIZATION): min/max=%s/%s" % (str(np.min(self.outputs_source)),str(np.max(self.outputs_source))))
			self.outputs_source= (self.outputs_source - targets_normmin)/(targets_normmax-targets_normmin)
			print("DEBUG: outputs_source (AFTER NORMALIZATION): min/max=%s/%s" % (str(np.min(self.outputs_source)),str(np.max(self.outputs_source))))

		print("DEBUG: inputs_source size=", np.shape(self.inputs_source))
		print("DEBUG: outputs_source size=", np.shape(self.outputs_source))
		print("DEBUG: outputs_labels_source size=", np.shape(self.outputs_labels_source))
		print("DEBUG: outputs_source=",self.outputs_source)
		print("DEBUG: outputs_labels_source=",self.outputs_labels_source)


		return 0

	#################################
	##     MAKE BKG TRAIN DATA
	#################################
	def make_bkg_train_data(self,writeimg=False):
		""" Prepare bkg train data """

		# - Init data
		input_data= []	
		output_size= self.nobjects*self.npars
		output_data= []  
		output_label_size= self.nobjects
		output_label_data= []
	
		nchannels= 1
		nx= self.img_sizex
		ny= self.img_sizey
		marginx= self.train_img_sizex/2
		marginy= self.train_img_sizey/2
		print("INFO: Input image size (%s,%s), margins(%s,%s), crop image size(%s,%s)" % (nx,ny,marginx,marginy,self.train_img_sizex,self.train_img_sizey))
		
		# - Extract nsamples img
		index= 0

		while index < self.nsamples_bkg:
			if index%100==0 :
				print("INFO: Generating bkg train image no. %s/%s ..." % (index+1,self.nsamples_bkg))
	
			if self.gen_bkg_from_img:
				# - Generate crop img center randomly
				x0= int(np.random.uniform(marginx,nx-marginx-1))
				y0= int(np.random.uniform(marginy,ny-marginy-1))
				print("INFO: Extract crop image around pos(%s,%s)" % (x0,y0))
			
				# - Extract crop img data
				data_crop= self.crop_img(x0,y0,self.train_img_sizex,self.train_img_sizey)

			else:
				# - Generate random bkg data
				data_crop= self.generate_noise(self.train_img_sizex,self.train_img_sizey,self.bkg_rms,self.bkg_mean)			

			imgcropsize= np.shape(data_crop)
			print("INFO: img crop shape=",imgcropsize)
			
			# - Check data integrity (skip if all zeros or nan/inf)
			n_nonzero= np.count_nonzero(data_crop)
			n_finite= (np.isfinite(data_crop)).sum()
			if n_nonzero<=0 or n_finite<=0:
				print ("WARN: Skip sample image (all pixels NaN/inf/zero)...")
				continue

			# - Save crop img to file?
			outfilename= 'train_bkg-RUN' + str(index+1) + '.fits'
			if writeimg:
				self.write_fits(data_crop,outfilename)

			# - Set train data as a tensor of size [Nsamples,Nx,Ny,Nchan] Nchan=1
			data_crop= data_crop.reshape(imgcropsize[0],imgcropsize[1],nchannels)
			input_data.append(data_crop)

			# - Set train target & labels
			output_data.append( np.zeros((1,output_size)) )
			output_label_data.append( np.zeros((1,output_label_size)) )

			index+= 1

		#- Convert list to array
		self.inputs_bkg= np.array(input_data)
		self.inputs_bkg= self.inputs_bkg.astype('float32')
	
		self.outputs_bkg= np.array(output_data)
		self.outputs_bkg= self.outputs_bkg.astype('float32')

		outputs_shape= self.outputs_bkg.shape
		N= outputs_shape[0]
		self.outputs_bkg= self.outputs_bkg.reshape((N,output_size))

		self.outputs_labels_bkg= np.array(output_label_data)
		self.outputs_labels_bkg= self.outputs_labels_bkg.astype('float32')
		self.outputs_labels_bkg= self.outputs_labels_bkg.reshape((N,output_label_size))

		# - Normalize to [0,1]
		if self.normalize_inputs:
			self.inputs_bkg= (self.inputs_bkg - self.normmin)/(self.normmax-self.normmin)

		print("DEBUG: inputs_bkg size=", np.shape(self.inputs_bkg))
		print("DEBUG: outputs_bkg size=", np.shape(self.outputs_bkg))
		print("DEBUG: outputs_labels_bkg size=", np.shape(self.outputs_labels_bkg))
		print("DEBUG: outputs_bkg=",self.outputs_bkg)
		print("DEBUG: outputs_labels_bkg=",self.outputs_labels_bkg)

		return 0

	#################################
	##     GENERATE NOISE IMAGE
	#################################
	def generate_noise(self,nx,ny,sigma,mean=0):
		""" Generate image data from random noise """
		data= np.random.normal(mean,sigma,(ny,nx))
		return data

	#################################
	##     GENERATE BLOB
	#################################
	def generate_blob(self,ampl,x0,y0,sigmax,sigmay,theta,trunc_thr=0.01):
		""" Generate a blob 
				Arguments: 
					ampl: peak flux in Jy
					x0, y0: gaussian means in pixels
					sigmax, sigmay: gaussian sigmas in pixels
					theta: rotation in degrees
					trunc_thr: truncation significance threshold
		"""
		data= Gaussian2D(ampl,x0,y0,sigmax,sigmay,theta=math.radians(theta))(self.gridx, self.gridy)
		
		## Truncate data such that sum(data)_trunc/sum(data)<f
		f= trunc_thr 
		if self.truncate_models:
			totFlux= (float)(np.sum(data,axis=None))
			
			data_vect_sorted= np.ravel(data)
			data_csum= np.cumsum(data_vect_sorted)/totFlux
			fluxThr= data_vect_sorted[np.argmin(data_csum<f)]
			data[data<fluxThr] = 0		

		return data

	#################################
	##     MAKE SOURCE TRAIN DATA
	#################################
	def make_source_train_data(self,writeimg=False):
		""" Prepare source train data """

		# - Init data
		input_data= []
		output_size= self.nobjects*self.npars
		output_data= []  
		output_label_size= self.nobjects
		output_label_data= []
	
		nchannels= 1
		nx= self.img_sizex
		ny= self.img_sizey
		marginx= self.train_img_sizex/2
		marginy= self.train_img_sizey/2
		marginx_source= self.source_gen_marginx
		marginy_source= self.source_gen_marginy

		# - Initialize grid for source generation	
		print("INFO: Generating grid for source generation ...")
		self.gridy, self.gridx = np.mgrid[0:self.train_img_sizey, 0:self.train_img_sizex]


		# - Set source randomization pars
		S_min= self.Smin 
		S_max= self.Smax 
		lgS_min= np.log10(S_min)
		lgS_max= np.log10(S_max)
		randomize_flux= False
		if self.Smin<self.Smax:
			randomize_flux= True

		randomize_gaus= False
		Bmaj_min= self.beam_bmaj_min
		Bmaj_max= self.beam_bmaj_max
		Bmin_min= self.beam_bmin_min
		Bmin_max= self.beam_bmin_max
		Pa_min= self.beam_bpa_min
		Pa_max= self.beam_bpa_max	
		if self.beam_bmaj_min<self.beam_bmaj_max:
			randomize_gaus= True
		if self.beam_bmin_min<self.beam_bmin_max:
			randomize_gaus= True
		if self.beam_bpa_min<self.beam_bpa_max:
			randomize_gaus= True
		
		# - Extract nsamples img
		index= 0

		while index < self.nsamples_source:
			if index%100==0 :
				print("INFO: Generating source train image no. %s/%s ..." % (index+1,self.nsamples_source))
	
			if self.gen_bkg_from_img:
				# - Generate crop img center randomly
				x0= int(np.random.uniform(marginx,nx-marginx-1))
				y0= int(np.random.uniform(marginy,ny-marginy-1))
				ix= int(np.round(x0))
				iy= int(np.round(y0))
				if self.img_data[iy,ix]==0 or np.isnan(self.img_data[iy,ix]):
					print("WARN: Skip sample image crop centered on (%s,%s) (pixel is zero or nan) ..." % (x0,y0))
					continue
			
				# - Extract crop img data
				data_crop= self.crop_img(x0,y0,self.train_img_sizex,self.train_img_sizey)

			else:
				# - Generate random bkg data
				data_crop= self.generate_noise(self.train_img_sizex,self.train_img_sizey,self.bkg_rms,self.bkg_mean)	


			imgcropsize= np.shape(data_crop)
			
			# - Check data integrity (skip if all zeros or nan/inf)
			n_nonzero= np.count_nonzero(data_crop)
			n_finite= (np.isfinite(data_crop)).sum()
			if n_nonzero<=0 or n_finite<=0:
				print ("WARN: Skip sample image crop centered on (%s,%s) (all pixels NaN/inf/zero) ..." % (x0,y0))
				continue

			# - Generate and add sources to cropped image
			nsources_max= int(round(np.random.uniform(1,self.nsources_max)))
			sources_data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.train_img_sizex, y_width=2*self.train_img_sizey)(self.gridx, self.gridy)
			mask_data = Box2D(amplitude=0,x_0=0,y_0=0,x_width=2*self.train_img_sizex, y_width=2*self.train_img_sizey)(self.gridx, self.gridy)
			source_pars= []
			nsources= 0

			print("INFO: Generating #%d sources in image ..." % (nsources_max))

			while nsources < nsources_max:
				# Generate source position
				x0_source= np.random.uniform(marginx_source,self.train_img_sizex-marginx_source-1)
				y0_source= np.random.uniform(marginy_source,self.train_img_sizey-marginy_source-1)
				ix= int(np.round(x0_source))
				iy= int(np.round(y0_source))
				
				# Skip if pixel already filled by a source or if crop data is nan
				if mask_data[iy,ix]!=0:
					print("WARN: Generated source position (%s,%s) is already taken (mask=%s) in image crop centered on (%s,%s), skip generation..." % (x0_source,y0_source,mask_data[iy,ix],x0,y0))
					continue
				if data_crop[iy,ix]==0 or np.isnan(data_crop[iy,ix]):
					print("WARN: Generated source position (%s,%s) on zero or nan pixel (data=%s) in image crop centered on (%s,%s), skip generation..." % (x0_source,y0_source,data_crop[iy,ix],x0,y0))
					continue
				

				## Generate flux uniform or expo in log
				S= S_min
				if randomize_flux:
					if self.Smodel=='uniform':
						lgS= np.random.uniform(lgS_min,lgS_max)
					elif self.Smodel=='exp':
						x= np.random.exponential(scale=1./self.Sslope)
						lgS= x + lgS_min
						if lgS>lgS_max:
							continue
					else:
						lgS= np.random.uniform(lgS_min,lgS_max)
					S= np.power(10,lgS)
				
				## Generate gaus pars
				if randomize_gaus:
					bmin= random.uniform(Bmin_min,Bmin_max)
					bmaj= random.uniform(bmin,Bmaj_max)
					pa= random.uniform(Pa_min,Pa_max)
				else:
					bmin= self.beam_bmin_min
					bmaj= self.beam_bmaj_min
					pa= self.beam_bpa_min

				sigmax= bmaj/(self.pixsize * SIGMA_TO_FWHM)
				sigmay= bmaj/(self.pixsize * SIGMA_TO_FWHM)
				theta = 90 + pa
				theta_rad= np.radians(theta)

				## Generate gaus 2D data
				print("DEBUG: Generating source no. %d: (x0,y0,S,sigmax,sigmay,theta)=(%s,%s,%s,%s,%s,%s)" % (nsources,x0_source,y0_source,S,sigmax,sigmay,theta))
				blob_data= self.generate_blob(ampl=S,x0=x0_source,y0=y0_source,sigmax=sigmax,sigmay=sigmay,theta=theta,trunc_thr=self.trunc_thr)
				if blob_data is None:
					print("WARN: Failed to generate Gaus2D (hint: too large trunc threshold), skip and regenerate...")
					continue

				#print("blob_data=",blob_data)

				## Update source mask & counts				
				sources_data+= blob_data
				mask_data[iy,ix]+= S
				nsources+= 1
				sname= 'S' + str(nsources)
				source_pars.append([sname,x0_source,y0_source,S,sigmax,sigmay,theta_rad])
				
			
			# - Add generated sources to train image
			data_crop+= sources_data

			# - Save crop img and source pars to file?
			outfilename= 'train_source-RUN' + str(index+1) + '.fits'
			outfilename_pars= 'train_source_pars-RUN' + str(index+1) + '.dat'
			if writeimg:
				self.write_fits(data_crop,outfilename)
				self.write_ascii(np.array(source_pars),outfilename_pars,'# sname x0(pix) y0(pix) S(Jy/beam) sigmaX(pix) sigmaY(pix) theta(rad)')
	
			
			# - Set train data as a tensor of size [Nsamples,Nx,Ny,Nchan] Nchan=1
			data_crop= data_crop.reshape(imgcropsize[0],imgcropsize[1],nchannels)
			input_data.append(data_crop)

			# - Set train targets
			targets= np.zeros((1,output_size))
			target_labels= np.zeros((1,output_label_size))
			
			nobjs= min(nsources_max,self.nobjects)
			par_counter= 0
			for k in range(nobjs):
				target_labels[0,k]= 1
				#x0= source_pars[k][1]
				#y0= source_pars[k][2]
				#S= source_pars[k][3]
				#sigmaX= source_pars[k][4]
				#sigmaY= source_pars[k][5]
				#theta= np.radians(source_pars[k][6])	
				#targets[0,par_counter+0]= x0
				#targets[0,par_counter+1]= y0
				#targets[0,par_counter+2]= S
				#targets[0,par_counter+3]= sigmaX
				#targets[0,par_counter+4]= sigmaY
				#targets[0,par_counter+5]= theta
				#par_counter+= 6

				for l in range(self.npars):			
					targets[0,par_counter+l]= source_pars[k][1+l]
					
				par_counter+= self.npars

			output_data.append(targets)
			output_label_data.append(target_labels)
			
			# Update sample counter
			index+= 1

		
		#- Convert list to array
		self.inputs_source= np.array(input_data)
		self.inputs_source= self.inputs_source.astype('float32')

		self.outputs_source= np.array(output_data)
		self.outputs_source= self.outputs_source.astype('float32')

		outputs_shape= self.outputs_source.shape
		N= outputs_shape[0]
		self.outputs_source= self.outputs_source.reshape((N,output_size))

		self.outputs_labels_source= np.array(output_label_data)
		self.outputs_labels_source= self.outputs_labels_source.astype('float32')
		self.outputs_labels_source= self.outputs_labels_source.reshape((N,output_label_size))

		# - Normalize to [0,1]
		if self.normalize_inputs:
			self.inputs_source= (self.inputs_source - self.normmin)/(self.normmax-self.normmin)

		print("DEBUG: inputs_source size=", np.shape(self.inputs_source))
		print("DEBUG: outputs_source size=", np.shape(self.outputs_source))
		print("DEBUG: outputs_labels_source size=", np.shape(self.outputs_labels_source))
		print("DEBUG: outputs_source=",self.outputs_source)
		print("DEBUG: outputs_labels_source=",self.outputs_labels_source)

		return 0

	#############################
	##     READ TRAIN DATA
	#############################
	def read_train_data(self):	
		""" Read train data from disk using input filelists """
		
		# - Read train data for bkg
		print("INFO: Reading train data for bkg ...")
		status= self.read_bkg_train_data(self.img_bkg_filelist)
		if status<0:
			return -1

		# - Read train data for source
		print("INFO: Reading train data for source ...")
		status= self.read_source_train_data(self.img_source_filelist,self.sourcepars_filelist)
		if status<0:
			return -1

		# - Merge data for bkg & sources
		print("INFO: Merging train data for bkg & sources ...")
		inputs= np.concatenate((self.inputs_bkg,self.inputs_source))
		outputs= np.concatenate((self.outputs_bkg,self.outputs_source))
		outputs_labels= np.concatenate((self.outputs_labels_bkg,self.outputs_labels_source))

		# - Shuffle data before splitting in test & validation sample
		print("INFO: Shuffling train data ...")
		indices= np.arange(inputs.shape[0])
		np.random.shuffle(indices)
		inputs= inputs[indices]
		outputs= outputs[indices]
		outputs_labels= outputs_labels[indices]
	
		print("DEBUG: inputs size=", np.shape(inputs))
		print("DEBUG: outputs size=", np.shape(outputs))
		print("DEBUG: outputs_labels size=", np.shape(outputs_labels))

		# - Partition the data into training and cross-validation splits
		print("INFO: Splitting data into train & test samples ...")
		split= train_test_split(
			inputs,outputs,outputs_labels, 
			test_size=self.test_size, 
			random_state=None
		)
		(self.inputs_train, self.inputs_test, self.outputs_train, self.outputs_test, self.outputs_labels_train, self.outputs_labels_test) = split

		print("DEBUG: inputs_train size=", np.shape(self.inputs_train))
		print("DEBUG: inputs_test size=", np.shape(self.inputs_test))
		print("DEBUG: outputs_train size=", np.shape(self.outputs_train))
		print("DEBUG: outputs_test size=", np.shape(self.outputs_test))
		print("DEBUG: outputs_labels_train size=", np.shape(self.outputs_labels_train))
		print("DEBUG: outputs_labels_test size=", np.shape(self.outputs_labels_test))

		return 0


	#############################
	##     GENERATE TRAIN DATA
	#############################
	def generate_train_data(self):
		""" Generate training data """
	
		# - Read input image
		if self.gen_bkg_from_img:
			print("INFO: Reading input image %s ..." % self.img_file)
			status= self.read_img()
			if status<0:
				return -1	

		# - Generate train data for bkg
		print("INFO: Generating train data for bkg ...")
		status= self.make_bkg_train_data(self.writeimg)
		if status<0:
			return -1	

		# - Generate train data for sources
		print("INFO: Generating train data for sources ...")
		status= self.make_source_train_data(self.writeimg)
		if status<0:
			return -1	

		# - Merge data for bkg & sources
		print("INFO: Merging train data for bkg & sources ...")
		inputs= np.concatenate((self.inputs_bkg,self.inputs_source))
		outputs= np.concatenate((self.outputs_bkg,self.outputs_source))
		outputs_labels= np.concatenate((self.outputs_labels_bkg,self.outputs_labels_source))

		# - Shuffle data before splitting in test & validation sample
		print("INFO: Shuffling train data ...")
		indices= np.arange(inputs.shape[0])
		np.random.shuffle(indices)
		inputs= inputs[indices]
		outputs= outputs[indices]
		outputs_labels= outputs_labels[indices]
	
		print("DEBUG: inputs size=", np.shape(inputs))
		print("DEBUG: outputs size=", np.shape(outputs))
		print("DEBUG: outputs_labels size=", np.shape(outputs_labels))

		# - Partition the data into training and cross-validation splits
		print("INFO: Splitting data into train & test samples ...")
		split= train_test_split(
			inputs,outputs,outputs_labels, 
			test_size=self.test_size, 
			random_state=None
		)
		(self.inputs_train, self.inputs_test, self.outputs_train, self.outputs_test, self.outputs_labels_train, self.outputs_labels_test) = split

		print("DEBUG: inputs_train size=", np.shape(self.inputs_train))
		print("DEBUG: inputs_test size=", np.shape(self.inputs_test))
		print("DEBUG: outputs_train size=", np.shape(self.outputs_train))
		print("DEBUG: outputs_test size=", np.shape(self.outputs_test))
		print("DEBUG: outputs_labels_train size=", np.shape(self.outputs_labels_train))
		print("DEBUG: outputs_labels_test size=", np.shape(self.outputs_labels_test))

		return 0

	#####################################
	##     NN LOSS FUNCTION FOR PARS
	#####################################
	def tot_loss(self,y_true,y_pred):

		####y_true = K.print_tensor(y_true, message='y_true = ')	
		####y_pred = K.print_tensor(y_pred, message='y_pred = ')
		####tf.print('y_true = ',y_true,output_stream=sys.stdout)
		####tf.print('y_pred = ',y_pred,output_stream=sys.stdout)
		#y_true = tf.Print(input_=y_true,data=[y_true], message='y_true = ',summarize=100)
		#y_pred = tf.Print(input_=y_pred,data=[y_pred], message='y_pred = ',summarize=100)
		
		# - Extract tensors relative to pars
		y_true_type= y_true[:,:self.nobjects]
		####y_true_type = K.print_tensor(y_true_type, message='y_true(type) = ')	
		####tf.print('y_true_type = ',y_true_type,output_stream=sys.stdout)
		#y_true_type= tf.Print(input_=y_true_type, data=[y_true_type], message='y_true_type = ',summarize=100)
		
		y_pred_type= y_pred[:,:self.nobjects]
		####y_pred_type = K.print_tensor(y_pred_type, message='y_pred(type) = ')		
		####tf.print('y_pred_type = ',y_pred_type)
		#y_pred_type= tf.Print(input_=y_pred_type, data=[y_pred_type], message='y_pred_type = ',summarize=100)
		
		# - Extract tensors relative to parameters
		y_true_pars= y_true[:,self.nobjects:]
		###y_true_pars = K.print_tensor(y_true_pars, message='y_true(pars) = ')
		###tf.print('y_true_pars = ',y_true_pars)
		#y_true_pars= tf.Print(input_=y_true_pars, data=[y_true_pars], message='y_true_pars = ',summarize=100)
		
		y_pred_pars= y_pred[:,self.nobjects:]
		####y_pred_pars = K.print_tensor(y_pred_pars, message='y_pred(pars) = ')	
		####tf.print('y_pred_pars = ',y_pred_pars)
		#y_pred_pars= tf.Print(input_=y_pred_pars, data=[y_pred_pars], message='y_pred_pars = ',summarize=100)
		

		# - Replicate true labels to have a tensor of size (N,nobjects*npars)
		#   This is multiplied by pars so that objects with label=0 are not entering in pars MSE computation (they sum zero)
		w= K.repeat_elements(y_true_type,self.npars,axis=1)
		#w= tf.Print(input_=w, data=[w], message='w= ',summarize=100)

		#y_diff= w*(y_pred_pars - y_true_pars)
		#y_diff= tf.Print(input_=y_diff, data=[y_diff], message='y_diff= ',summarize=100)
		#y_diff_sqr= K.square(y_diff)
		#y_diff_sqr= tf.Print(input_=y_diff_sqr, data=[y_diff_sqr], message='y_diff_sqr= ',summarize=100)

		# - Count number of objects
		N= tf.count_nonzero(y_true_type,dtype=tf.dtypes.float32)
		#N= tf.Print(input_=N, data=[N], message='N= ',summarize=100)

		# - Compute MSE for pars relative to target objects (not background)
		#mse= keras.metrics.mean_squared_error(y_pred,y_true)
		#mse= K.mean(K.square( w*(y_pred_pars - y_true_pars) ), axis=-1)
		#mse= K.mean(K.square( y_diff ), axis=-1)
		#mse= K.sum(y_diff_sqr, axis=-1)/N
		mse= K.sum(K.square( w*(y_pred_pars - y_true_pars) ), axis=-1)/N
		#tf.Print(input_=mse, data=[mse], message='mse= ',summarize=100)

		# - Compute binary crossentropy for labels
		binaryCE = keras.metrics.binary_crossentropy(y_true_type,y_pred_type)

		# - Compute total loss as weighted sum of MSE + binaryCE		
		tot= self.spars_loss_weight*mse + self.labels_loss_weight*binaryCE
		#tot= mse

	
		return tot

	#####################################
	##     BUILD NETWORK
	#####################################
	def build_network(self):

		if self.standard_nn_build_enabled:
			#- Create a standard network
			print("INFO: Building standard network architecture ...")
			status= self.build_standard_network()
			if status<0:
				return -1
		else:
			#- Create the network from file
			print("INFO: Building network architecture from file %s ..." % self.nnarc_file)
			status= self.build_network_from_file(self.nnarc_file)
			if status<0:
				return -1

		return 0

	#####################################
	##     BUILD NETWORK FROM SPEC FILE
	#####################################
	def build_network_from_file(self,filename):
		""" Building deep network from file """

		# - Read NN architecture file
		nn_data= []
		skip_patterns= ['#']
		try:
			nn_data= self.read_ascii(filename,skip_patterns)
		except IOError:
			print("ERROR: Failed to read nn arc file %d!" % filename)
			return -1

		nlayers= np.shape(nn_data)[0]
		
		# - Input layer
		nchan= 1	
		inputShape = (self.train_img_sizey, self.train_img_sizex, nchan)
		inputs= Input(shape=inputShape,dtype='float', name='input')
		x= inputs

		# - Parse NN architecture file and create intermediate layers
		for index in range(nlayers):
			layer_info= nn_data[index]
			print("Layer no. %d: %s" % (index,layer_info))

			layer_type= layer_info[0]

			# - Add Conv2D layer?
			if layer_type=='Conv2D':
				nfields= len(layer_info)
				if nfields!=5:
					print("ERROR: Invalid number of fields (n=%d) given in Conv2D layer specification (5 expected)" % nfields)
				nfilters= int(layer_info[1])
				kernSize= int(layer_info[2])
				activation= str(layer_info[3])
				padding= str(layer_info[4])
				x = layers.Conv2D(filters=nfilters, kernel_size=(kernSize,kernSize), activation=activation, padding=padding)(x)			
	
			# - Add MaxPooling2D layer?
			elif layer_type=='MaxPooling2D':
				nfields= len(layer_info)
				if nfields!=3:
					print("ERROR: Invalid number of fields (n=%d) given in MaxPooling2D layer specification (3 expected)" % nfields)
				poolSize= int(layer_info[1])
				padding= str(layer_info[2])
				x = layers.MaxPooling2D(pool_size=(poolSize,poolSize),strides=None,padding=padding)(x)

			# - Add Dropout layer?
			elif layer_type=='Dropout':
				nfields= len(layer_info)
				if nfields!=2:
					print("ERROR: Invalid number of fields (n=%d) given in Dropout layer specification (2 expected)" % nfields)
				dropout= float(layer_info[1])
				x = layers.Dropout(dropout)(x)

			# - Add BatchNormalization layer?
			elif layer_type=='BatchNormalization':
				x = layers.BatchNormalization()(x)
	
			# - Add Flatten layer?
			elif layer_type=='Flatten':
				x = layers.Flatten()(x)

			# - Add Dense layer?
			elif layer_type=='Dense':
				nfields= len(layer_info)
				if nfields!=3:
					print("ERROR: Invalid number of fields (n=%d) given in Dense layer specification (3 expected)" % nfields)
				nNeurons= int(layer_info[1])
				activation= str(layer_info[2])
				x = layers.Dense(nNeurons, activation=activation)(x)

			else:
				print("ERROR: Invalid/unknown layer type parsed (%s)!" % layer_type)
				return -1
			
		# - Output layers
		type_prediction = layers.Dense(self.nobjects, activation='sigmoid', name='type')(x)
		pars_prediction = layers.Dense(self.nobjects*self.npars, activation='linear', name='pars')(x)

		# - Concatenate output layers
		nn_prediction= layers.concatenate([type_prediction,pars_prediction],name='nnout')

		# - Create NN model
		self.model = Model(
				inputs=inputs,
				#outputs=[type_prediction, pars_prediction],
				outputs=nn_prediction,
				name="SourceNet"
		)

		# - Print network architecture
		self.model.summary()

		#- Set optimizer & loss function per each output
		print("INFO: Compiling network ...")
		if self.optimizer=='rmsprop':
			opt= optimizers.RMSprop(lr=self.learning_rate)
		elif self.optimizer=='sgd':
			#opt= optimizers.SGD(lr=self.learning_rate, decay=1e-6, momentum=0.9, nesterov=True)
			opt= optimizers.SGD(lr=self.learning_rate, nesterov=False)
		elif self.optimizer=='sgdn':
			opt= optimizers.SGD(lr=self.learning_rate, nesterov=True)
		else:
			opt= optimizers.RMSprop(lr=self.learning_rate)
		
		#opt= Adam(lr=INIT_LR, decay=INIT_LR / nepochs)
	
		losses = {
			"type": "binary_crossentropy",
			"pars": "mse"
		}
		lossWeights = {
			"type": self.labels_loss_weight,
			"pars": self.spars_loss_weight
		}

		#self.model.compile(optimizer=opt,loss=losses, loss_weights=lossWeights, metrics=['accuracy'])
		self.model.compile(optimizer=opt,loss=self.tot_loss, metrics=['accuracy'])
		
		return 0

	###########################
	##     BUILD NETWORK
	###########################
	def build_standard_network(self):
		""" Building a standard deep network """
		
		# - Input layer
		nchan= 1	
		inputShape = (self.train_img_sizey, self.train_img_sizex, nchan)
		inputs= Input(shape=inputShape,dtype='float', name='input')

		# - Convolutional layers
		x= inputs

		if self.conv_enabled:
			# - CONV LAYER 1
			x = layers.Conv2D(filters=self.conv_nfilt_min, kernel_size=(self.conv_kern_size_min,self.conv_kern_size_min), activation=self.conv_act, padding=self.padding)(inputs)
			if self.maxpool_enabled:
				x = layers.MaxPooling2D(pool_size=(self.pool_size,self.pool_size),strides=None,padding='valid')(x)
			if self.dropout_enabled:
				x = layers.Dropout(self.dropout)(x)

			# - CONV LAYER 2
			x = layers.Conv2D(filters=self.conv_nfilt_max, kernel_size=(self.conv_kern_size_min,self.conv_kern_size_min), activation=self.conv_act, padding=self.padding)(x)
			if self.maxpool_enabled:
				x = layers.MaxPooling2D(pool_size=(self.pool_size,self.pool_size),strides=None,padding='valid')(x)
			if self.dropout_enabled:
				x = layers.Dropout(self.dropout)(x)

			# - CONV LAYER 3
			x = layers.Conv2D(filters=self.conv_nfilt_max, kernel_size=(self.conv_kern_size_max,self.conv_kern_size_max), activation=self.conv_act, padding=self.padding)(x)
			if self.maxpool_enabled:	
				x = layers.MaxPooling2D(pool_size=(self.pool_size,self.pool_size),strides=None,padding='valid')(x)
			if self.dropout_enabled:
				x = layers.Dropout(self.dropout)(x)

			# - Batch normalization
			if self.batchnorm_enabled:
				x = layers.BatchNormalization()(x)
		

		# - Fully connected layers
		x = layers.Flatten()(x)
		x = layers.Dense(self.dense_size_max, activation=self.dense_act)(x)
		if self.dropout_enabled:
			x = layers.Dropout(self.dropout_dense)(x)
		x = layers.Dense(self.dense_size_min, activation=self.dense_act)(x)
		if self.dropout_enabled:
			x = layers.Dropout(self.dropout_dense)(x)

		

		# - Batch normalization
		if self.batchnorm_enabled:
			x = layers.BatchNormalization()(x)

		# - Output layers
		type_prediction = layers.Dense(self.nobjects, activation='sigmoid', name='type')(x)
		pars_prediction = layers.Dense(self.nobjects*self.npars, activation='linear', name='pars')(x)

		# - Create NN model
		self.model = Model(
				inputs=inputs,
				outputs=[type_prediction, pars_prediction],
				name="SourceNet"
		)

		# - Print network architecture
		self.model.summary()

		#- Set optimizer & loss function per each output
		print("INFO: Compiling network ...")
		opt= optimizers.RMSprop(lr=1.e-4)
		#opt= Adam(lr=INIT_LR, decay=INIT_LR / nepochs)
	
		losses = {
			"type": "binary_crossentropy",
			"pars": "mse"
		}
		lossWeights = {
			"type": self.labels_loss_weight,
			"pars": self.spars_loss_weight
		}

		##self.model.compile(optimizer=opt,loss=losses, loss_weights=lossWeights, metrics=['accuracy',precision,recall,f1_score])
		#self.model.compile(optimizer=opt,loss=losses, loss_weights=lossWeights, metrics=['accuracy',my_mse_loss])
		self.model.compile(optimizer=opt,loss=self.tot_loss, metrics=['accuracy'])
		

		return 0


	###########################
	##     TRAIN NETWORK
	###########################
	def train_network(self):
		""" Training deep network """

		#self.fitout= self.model.fit(
			#	x=self.inputs_train, 
			#	y={"type": self.outputs_labels_train,"pars": self.outputs_train},
			#	validation_data=(self.inputs_test,{"type": self.outputs_labels_test, "pars": self.outputs_test}),
			#	##batch_size=64
			#	epochs=self.nepochs,
			#	verbose=1
		#)

		self.flipped_outputs_labels_train= self.outputs_labels_train
		self.flipped_outputs_train= self.outputs_train
		self.flipped_outputs_labels_test= self.outputs_labels_test
		self.flipped_outputs_test= self.outputs_test

		#self.train_loss_vs_epoch= np.zeros((3,self.nepochs))	
		#self.test_loss_vs_epoch= np.zeros((3,self.nepochs))
		#self.train_accuracy_vs_epoch= np.zeros((2,self.nepochs))
		#self.test_accuracy_vs_epoch= np.zeros((2,self.nepochs))	
		self.train_loss_vs_epoch= np.zeros((1,self.nepochs))	
		self.test_loss_vs_epoch= np.zeros((1,self.nepochs))
		self.train_accuracy_vs_epoch= np.zeros((1,self.nepochs))
		self.test_accuracy_vs_epoch= np.zeros((1,self.nepochs))	

		for epoch in range(self.nepochs):
		
			# - Train for 1 epoch
			self.fitout= self.model.fit(
				x=self.inputs_train, 
				#y={"type": self.flipped_outputs_labels_train,"pars": self.flipped_outputs_train},
				y={"nnout": np.concatenate((self.flipped_outputs_labels_train,self.flipped_outputs_train),axis=1) },
				#validation_data=(self.inputs_test,{"type": self.outputs_labels_test,"pars": self.outputs_test}),
				validation_data=(self.inputs_test,{"nnout": np.concatenate((self.outputs_labels_test,self.outputs_test),axis=1) }),
				#batch_size=64
				epochs=1,
				verbose=2,
				#callbacks=[self.nnprinter_cb]
			)

			# - Save epoch loss
			print ('== EPOCH %d ==' % epoch)
			print (self.fitout.history)
			self.train_loss_vs_epoch[0,epoch]= self.fitout.history['loss'][0]
			#self.train_loss_vs_epoch[1,epoch]= self.fitout.history['type_loss'][0]
			#self.train_loss_vs_epoch[2,epoch]= self.fitout.history['pars_loss'][0]
			self.test_loss_vs_epoch[0,epoch]= self.fitout.history['val_loss'][0]
			#self.test_loss_vs_epoch[1,epoch]= self.fitout.history['val_type_loss'][0]
			#self.test_loss_vs_epoch[2,epoch]= self.fitout.history['val_pars_loss'][0]
			self.train_accuracy_vs_epoch[0,epoch]= self.fitout.history['acc'][0]
			#self.train_accuracy_vs_epoch[0,epoch]= self.fitout.history['type_acc'][0]
			#self.train_accuracy_vs_epoch[1,epoch]= self.fitout.history['pars_acc'][0]
			self.test_accuracy_vs_epoch[0,epoch]= self.fitout.history['val_acc'][0]
			#self.test_accuracy_vs_epoch[0,epoch]= self.fitout.history['val_type_acc'][0]
			#self.test_accuracy_vs_epoch[1,epoch]= self.fitout.history['val_pars_acc'][0]

			# - Get predictions for train data and flip targets according to smallest MSE match
			npars_flip= 2 # only (x,y) pars used for flipping
			#npars_flip= self.npars

			if self.flip_train_data:
				#(outputs_labels_pred, outputs_pred)= self.model.predict(self.inputs_train)
				nnout_pred= self.model.predict(self.inputs_train)
				outputs_labels_pred= nnout_pred[:,:self.nobjects]
				outputs_pred= nnout_pred[:,self.nobjects:]

				nsamples= outputs_pred.shape[0]
				
				for sample in range(nsamples):
			
					mses= np.zeros((self.nobjects,self.nobjects))
					predout_mse= np.zeros((self.nobjects,npars_flip))
					expout_mse= np.zeros((self.nobjects,npars_flip))
					predout= np.zeros((self.nobjects,self.npars))
					expout= np.zeros((self.nobjects,self.npars))
					expout_labels= np.zeros((self.nobjects,1))
			
					for i in range(self.nobjects):
						expout_labels[i,0]= self.flipped_outputs_labels_train[sample,i]
						for j in range(self.npars):
							predout[i,j]= outputs_pred[sample,j+i*self.npars]
							expout[i,j]= self.flipped_outputs_train[sample,j+i*self.npars]
						for j in range(npars_flip):
							predout_mse[i,j]= outputs_pred[sample,j+i*self.npars]
							expout_mse[i,j]= self.flipped_outputs_train[sample,j+i*self.npars]
					
					for i in range(self.nobjects):
						for j in range(self.nobjects):
							mse= np.mean(np.square(expout_mse[i,:]-predout_mse[j,:]))	
							mses[i,j]= mse
	
			
					# - Find new ordering according to smallest MSE
					mses_copy= mses
					reorder_indexes= np.zeros(self.nobjects,dtype=int)
					for i in range(self.nobjects):
						ind_exp, ind_pred= np.unravel_index(mses.argmin(),mses.shape) # Find index of smallest mse
						mses[ind_exp]= np.Inf # Set mse to largest value so that it is not re-assigned anymore
						mses[:,ind_pred]= np.Inf 
						reorder_indexes[ind_pred]= ind_exp	
				
					#- Save before flipping
					#target= ','.join(map(str, self.flipped_outputs_train[sample,:]))
					#target_labels= ','.join(map(str, self.flipped_outputs_labels_train[sample,:]))

					#- Flip target
					self.flipped_outputs_train[sample]= expout[reorder_indexes].flatten() 
					self.flipped_outputs_labels_train[sample]= expout_labels[reorder_indexes].flatten() 
			
					#- Print
					#flipped_target= ','.join(map(str, flipped_outputs_train[sample,:]))
					#flipped_target_labels= ','.join(map(str, flipped_outputs_labels_train[sample,:]))
					#pred= ','.join(map(str, outputs_pred[sample,:]))
					#pred_labels= ','.join(map(str, outputs_labels_pred[sample,:]))
					#mse= ','.join(map(str, mses_copy[:,:]))
					#print("DEBUG: Entry no. %d: reorder_indexes=[%s], target_labels=[%s], flipped_target_labels=[%s]" % (sample+1,reorder_indexes,target_labels,flipped_target_labels) )
					#print("DEBUG: Entry no. %d: pred_labels=[%s], target_labels=[%s], flipped_target_labels=[%s], pred=[%s], target=[%s], flipped_target=[%s], mse=[%s], reorder_indexes=[%s]" % (sample+1,pred_labels,target_labels,flipped_target_labels,pred,target,flipped_target,mse,reorder_indexes) )
		
			# - Get predictions for test sample and flip according to smallest MSE match
			if self.flip_test_data:
				#(outputs_labels_pred, outputs_pred)= self.model.predict(self.inputs_test)
				nnout_pred= self.model.predict(self.inputs_test)
				outputs_labels_pred= nnout_pred[:,:self.nobjects]
				outputs_pred= nnout_pred[:,self.nobjects:]
				nsamples= outputs_pred.shape[0]
		
				for sample in range(nsamples):
			
					mses= np.zeros((self.nobjects,self.nobjects))
					predout_mse= np.zeros((self.nobjects,npars_flip))
					expout_mse= np.zeros((self.nobjects,npars_flip))
					predout= np.zeros((self.nobjects,self.npars))
					expout= np.zeros((self.nobjects,self.npars))
					expout_labels= np.zeros((self.nobjects,1))
			
					for i in range(self.nobjects):
						expout_labels[i,0]= self.flipped_outputs_labels_test[sample,i]
						for j in range(self.npars):
							predout[i,j]= outputs_pred[sample,j+i*self.npars]
							expout[i,j]= self.flipped_outputs_test[sample,j+i*self.npars]
						for j in range(npars_flip):
							predout_mse[i,j]= outputs_pred[sample,j+i*self.npars]
							expout_mse[i,j]= self.flipped_outputs_test[sample,j+i*self.npars]
					
					for i in range(self.nobjects):
						for j in range(self.nobjects):
							mse= np.mean(np.square(expout_mse[i,:]-predout_mse[j,:]))	
							mses[i,j]= mse
	
					# - Find new ordering according to smallest MSE
					mses_copy= mses
					reorder_indexes= np.zeros(self.nobjects,dtype=int)
					for i in range(self.nobjects):
						ind_exp, ind_pred= np.unravel_index(mses.argmin(),mses.shape) # Find index of smallest mse
						mses[ind_exp]= np.Inf # Set mse to largest value so that it is not re-assigned anymore
						mses[:,ind_pred]= np.Inf 
						reorder_indexes[ind_pred]= ind_exp	
				
					#- Save before flipping
					#target= ','.join(map(str, flipped_outputs_test[sample,:]))
					#target_labels= ','.join(map(str, flipped_outputs_labels_test[sample,:]))

					#- Flip target
					self.flipped_outputs_test[sample]= expout[reorder_indexes].flatten() 
					self.flipped_outputs_labels_test[sample]= expout_labels[reorder_indexes].flatten() 
			
					#- Print
					#flipped_target= ','.join(map(str, flipped_outputs_test[sample,:]))
					#flipped_target_labels= ','.join(map(str, flipped_outputs_labels_test[sample,:]))
					#pred= ','.join(map(str, outputs_pred[sample,:]))
					#pred_labels= ','.join(map(str, outputs_labels_pred[sample,:]))
					#mse= ','.join(map(str, mses_copy[:,:]))

		#===========================
		#==   SAVE NN
		#===========================
		#- Save the model weights
		print("INFO: Saving model weights ...")
		self.model.save_weights('model_weights.h5')

		# -Save the model architecture in json format
		with open('model_architecture.json', 'w') as f:
			f.write(self.model.to_json())
		
		#- Save the model
		print("INFO: Saving model ...")
		self.model.save('model.h5')

		return 0

	###########################
	##     EVALUATE NETWORK
	###########################
	def evaluate_network(self):
		""" Evaluating network """

		#================================
		#==   SAVE TRAIN METRICS
		#================================
		N= self.train_loss_vs_epoch.shape[1]
		epoch_ids= np.array(range(N))
		epoch_ids+= 1
		epoch_ids= epoch_ids.reshape(N,1)

		#print("DEBUG: train_loss_vs_epoch shape=",self.train_loss_vs_epoch.reshape(N,3).shape)
		print("DEBUG: epoch_ids shape=",epoch_ids.shape)

		metrics_data= np.concatenate(
			#(epoch_ids,self.train_loss_vs_epoch.reshape(N,3),self.test_loss_vs_epoch.reshape(N,3),self.train_accuracy_vs_epoch.reshape(N,2),self.test_accuracy_vs_epoch.reshape(N,2)),
			(epoch_ids,self.train_loss_vs_epoch.reshape(N,1),self.test_loss_vs_epoch.reshape(N,1),self.train_accuracy_vs_epoch.reshape(N,1),self.test_accuracy_vs_epoch.reshape(N,1)),
			axis=1
		)
			
		self.write_ascii(metrics_data,self.outfile_nnout_metrics,'# epoch - tot loss - type loss - pars loss - tot loss (test) - type loss (test) - pars loss (test) - type acc - pars acc - type acc (test) - pars acc (test)')	

		
		#================================
		#==   EVALUATE NN ON TRAIN DATA
		#================================
		print("INFO: Classifying train data ...")
		#(predictions_labels_train, predictions_train)= self.model.predict(self.inputs_train)
		nnout_train= self.model.predict(self.inputs_train)
		predictions_labels_train= nnout_train[:,:self.nobjects]
		predictions_train= nnout_train[:,self.nobjects:]

		for i in range(predictions_labels_train.shape[0]):
			#target= ','.join(map(str, self.outputs_labels_train[i,:]))
			target= ','.join(map(str, self.flipped_outputs_labels_train[i,:]))	
			pred= ','.join(map(str, predictions_labels_train[i,:]))
			print("DEBUG: Train labels entry no. %d: target=[%s], pred=[%s]" % (i+1,target,pred) )

		for i in range(predictions_train.shape[0]):
			#target= ','.join(map(str, self.outputs_train[i,:]))
			target= ','.join(map(str, self.flipped_outputs_train[i,:]))
			pred= ','.join(map(str, predictions_train[i,:]))
			print("DEBUG: Train spars entry no. %d: target=[%s], pred=[%s]" % (i+1,target,pred) )


		# - Computing true & false detections
		nsamples_train= self.outputs_labels_train.shape[0]
		detThr= 0.5
		nobjs_tot= 0
		nobjs_true= 0
		nobjs_rec= 0
		nobjs_rec_true= 0
		nobjs_rec_false= 0
		s_list= []
		xpull_list= []
		ypull_list= []
		spull_list= []

		nnout_train= []
	
		for i in range(nsamples_train):
			#target= self.outputs_labels_train[i,:]
			target= self.flipped_outputs_labels_train[i,:]
			pred= predictions_labels_train[i,:]	
			#target_pars= self.outputs_train[i,:]
			target_pars= self.flipped_outputs_train[i,:]
			pred_pars= predictions_train[i,:]

			true_obj_indexes= np.argwhere(target==1).flatten()
			rec_obj_indexes= np.argwhere(pred>detThr).flatten()
			n= len(true_obj_indexes)
			nrec= len(rec_obj_indexes)
			ntrue= 0
			nrec_true= 0
			nrec_false= 0

			
			for index in range(self.nobjects):
				obj_data= []
				obj_data.append(target[index])
				obj_data.append(pred[index])
				for k in range(self.npars):
					obj_data.append(target_pars[k+index*self.npars])
					obj_data.append(pred_pars[k+index*self.npars])
				
				nnout_train.append(obj_data)
			
	
			for index in true_obj_indexes:
				x0_true= target_pars[0 + index*self.npars]
				y0_true= target_pars[1 + index*self.npars]
				S_true= target_pars[2 + index*self.npars]
				if pred[index]>detThr:
					ntrue+= 1
					x0_rec= pred_pars[0 + index*self.npars]
					y0_rec= pred_pars[1 + index*self.npars]
					S_rec= pred_pars[2 + index*self.npars]
					s_list.append(np.log10(S_true))
					spull_list.append(S_rec/S_true-1)
					xpull_list.append(x0_rec-x0_true)
					ypull_list.append(y0_rec-y0_true)

			for index in rec_obj_indexes:
				if target[index]==1:
					nrec_true+= 1
				else:
					nrec_false+= 1
	
			nobjs_tot+= n
			nobjs_rec+= nrec
			nobjs_true+= ntrue
			nobjs_rec_true+= nrec_true 
			nobjs_rec_false+= nrec_false

		completeness_train= float(nobjs_true)/float(nobjs_tot)
		reliability_train= float(nobjs_rec_true)/float(nobjs_rec)

		print("INFO: NN Train Results: Completeness(det/tot=%d/%d)=%s, Reliability(true/rec=%d/%d)=%s" % (nobjs_true,nobjs_tot,str(completeness_train),nobjs_rec_true,nobjs_rec,str(reliability_train)))

		# - Write ascii file with results
		self.write_ascii(np.array(nnout_train),self.outfile_nnout_train,'# target label - predicted label - target pars - predicted pars')	

		#================================
		#==   EVALUATE NN ON TEST DATA
		#================================
		print("INFO: Classifying test data ...")
		#(predictions_labels_test, predictions_test)= self.model.predict(self.inputs_test)
		nnout_test= self.model.predict(self.inputs_test)
		predictions_labels_test= nnout_test[:,:self.nobjects]
		predictions_test= nnout_test[:,self.nobjects:]

		nsamples_test= self.outputs_labels_test.shape[0]
		detThr= 0.5
		nobjs_tot= 0
		nobjs_true= 0
		nobjs_rec= 0
		nobjs_rec_true= 0
		nobjs_rec_false= 0
		s_list_test= []
		xpull_list_test= []
		ypull_list_test= []
		spull_list_test= []
		nnout_test= []

		for i in range(nsamples_test):
			#target= self.outputs_labels_test[i,:]
			target= self.flipped_outputs_labels_test[i,:]
			pred= predictions_labels_test[i,:]	
			#target_pars= self.outputs_test[i,:]
			target_pars= self.flipped_outputs_test[i,:]
			pred_pars= predictions_test[i,:]

			true_obj_indexes= np.argwhere(target==1).flatten()
			rec_obj_indexes= np.argwhere(pred>detThr).flatten()
			n= len(true_obj_indexes)
			nrec= len(rec_obj_indexes)
			ntrue= 0
			nrec_true= 0
			nrec_false= 0
		
			for index in range(self.nobjects):
				obj_data= []
				obj_data.append(target[index])
				obj_data.append(pred[index])
				for k in range(self.npars):
					obj_data.append(target_pars[k+index*self.npars])
					obj_data.append(pred_pars[k+index*self.npars])
				
				nnout_test.append(obj_data)

			for index in true_obj_indexes:
				x0_true= target_pars[0 + index*self.npars]
				y0_true= target_pars[1 + index*self.npars]
				S_true= target_pars[2 + index*self.npars]
				if pred[index]>detThr:
					ntrue+= 1
					x0_rec= pred_pars[0 + index*self.npars]
					y0_rec= pred_pars[1 + index*self.npars]
					S_rec= pred_pars[2 + index*self.npars]
					s_list_test.append(np.log10(S_true))
					spull_list_test.append(S_rec/S_true-1)
					xpull_list_test.append(x0_rec-x0_true)
					ypull_list_test.append(y0_rec-y0_true)

			for index in rec_obj_indexes:
				if target[index]==1:
					nrec_true+= 1
				else:
					nrec_false+= 1

			nobjs_tot+= n
			nobjs_rec+= nrec
			nobjs_true+= ntrue
			nobjs_rec_true+= nrec_true 
			nobjs_rec_false+= nrec_false

		completeness_test= float(nobjs_true)/float(nobjs_tot)
		reliability_test= float(nobjs_rec_true)/float(nobjs_rec)

		print("INFO: NN Test Results: Completeness(det/tot=%d/%d)=%s, Reliability(true/rec=%d/%d)=%s" % (nobjs_true,nobjs_tot,str(completeness_test),nobjs_rec_true,nobjs_rec,str(reliability_test)))

		#for i in range(predictions_labels_test.shape[0]):
		#	target= ','.join(map(str, self.outputs_labels_test[i,:]))
		#	pred= ','.join(map(str, predictions_labels_test[i,:]))
		#	print("INFO: Test labels entry no. %d: target=[%s], pred=[%s]" % (i+1,target,pred) )

		#for i in range(predictions_test.shape[0]):
		#	target= ','.join(map(str, self.outputs_test[i,:]))
		#	pred= ','.join(map(str, predictions_test[i,:]))
		#	print("INFO: Test spars entry no. %d: target=[%s], pred=[%s]" % (i+1,target,pred) )

		# - Write ascii file with results
		self.write_ascii(np.array(nnout_test),self.outfile_nnout_test,'# target label - predicted label - target pars - predicted pars')	


		#===========================
		#==   PLOT NN RESULTS
		#===========================
		# - Plot the network	
		print("INFO: Printing network model architecture to file ...")
		plot_model(self.model, to_file=self.outfile_model)

		# - Plot the total loss, type loss, spars loss
		#lossNames = ["loss", "type_loss", "pars_loss"]
		lossNames = ["loss"]
		plt.style.use("ggplot")
		#(fig, ax) = plt.subplots(3, 1, figsize=(13, 13))
		(fig, ax) = plt.subplots(1, 1, figsize=(13, 13),squeeze=False)
		
		for (i, lossName) in enumerate(lossNames):
			# Plot the loss for both the training and validation data
			#title = "Loss for {}".format(lossName) if lossName != "loss" else "Total loss"
			#ax[i].set_title(title)
			#ax[i].set_xlabel("Epoch #")
			#ax[i].set_ylabel("Loss")
			#ax[i].plot(np.arange(0, self.nepochs), self.train_loss_vs_epoch[i], label="TRAIN SAMPLE - " + lossName)
			#ax[i].plot(np.arange(0, self.nepochs), self.test_loss_vs_epoch[i], label="TEST SAMPLE - " + lossName)
			#ax[i].legend()	

			title = "Loss for {}".format(lossName) if lossName != "loss" else "Total loss"
			ax[i,0].set_title(title)
			ax[i,0].set_xlabel("Epoch #")
			ax[i,0].set_ylabel("Loss")
			ax[i,0].plot(np.arange(0, self.nepochs), self.train_loss_vs_epoch[i], label="TRAIN SAMPLE - " + lossName)
			ax[i,0].plot(np.arange(0, self.nepochs), self.test_loss_vs_epoch[i], label="TEST SAMPLE - " + lossName)
			ax[i,0].legend()

		plt.tight_layout()
		plt.savefig(self.outfile_loss)
		plt.close()

		# - Plot the accuracy
		print("INFO: Plot the network accuracy metric to file ...")
		#accuracyNames = ["type_acc", "pars_acc"]
		accuracyNames = ["acc"]
		plt.style.use("ggplot")
		#(fig, ax) = plt.subplots(2, 1, figsize=(8, 8))
		(fig, ax) = plt.subplots(1, 1, figsize=(8, 8),squeeze=False)

		for (i, accuracyName) in enumerate(accuracyNames):
			# Plot the loss for both the training and validation data
			#ax[i].set_title("Accuracy for {}".format(accuracyName))
			#ax[i].set_xlabel("Epoch #")
			#ax[i].set_ylabel("Accuracy")
			#ax[i].plot(np.arange(0, self.nepochs), self.train_accuracy_vs_epoch[i], label="TRAIN SAMPLE - " + accuracyName)
			#ax[i].plot(np.arange(0, self.nepochs), self.test_accuracy_vs_epoch[i], label="TEST SAMPLE - " + accuracyName)
			#ax[i].legend()

			ax[i,0].set_title("Accuracy for {}".format(accuracyName))
			ax[i,0].set_xlabel("Epoch #")
			ax[i,0].set_ylabel("Accuracy")
			ax[i,0].plot(np.arange(0, self.nepochs), self.train_accuracy_vs_epoch[i], label="TRAIN SAMPLE - " + accuracyName)
			ax[i,0].plot(np.arange(0, self.nepochs), self.test_accuracy_vs_epoch[i], label="TEST SAMPLE - " + accuracyName)
			ax[i,0].legend()

		plt.tight_layout()
		plt.savefig(self.outfile_accuracy)
		plt.close()

		# - Plot x, y position reco accuracy for detected sources
		print("INFO: Plot the source (x, y) position accuracy ...")
		plt.style.use("ggplot")
		(fig, ax) = plt.subplots(2, 2, figsize=(8, 8))

		ax[0,0].set_title("x Position Accuracy")
		ax[0,0].set_xlabel("logS (Jy/beam)")
		ax[0,0].set_ylabel("dx")
		ax[0,0].scatter(np.array(s_list),np.array(xpull_list),label="TRAIN SAMPLE")
		ax[0,0].legend()

		ax[0,1].set_title("y Position Accuracy")
		ax[0,1].set_xlabel("logS (Jy/beam)")
		ax[0,1].set_ylabel("dy")
		ax[0,1].scatter(np.array(s_list),np.array(ypull_list),label="TRAIN SAMPLE")
		ax[0,1].legend()

		ax[1,0].set_title("x Position Accuracy")
		ax[1,0].set_xlabel("logS (Jy/beam)")
		ax[1,0].set_ylabel("dx")
		ax[1,0].scatter(np.array(s_list_test),np.array(xpull_list_test),label="TEST SAMPLE")
		ax[1,0].legend()

		ax[1,1].set_title("y Position Accuracy")
		ax[1,1].set_xlabel("logS (Jy/beam)")
		ax[1,1].set_ylabel("dy")
		ax[1,1].scatter(np.array(s_list_test),np.array(ypull_list_test),label="TEST SAMPLE")
		ax[1,1].legend()

		plt.tight_layout()
		plt.savefig(self.outfile_posaccuracy)
		plt.close()

		# - Plot flux reco accuracy for detected sources
		print("INFO: Plot the source flux accuracy ...")
		plt.style.use("ggplot")
		(fig, ax) = plt.subplots(2, 1, figsize=(8, 8))

		ax[0].set_title("Flux Accuracy")
		ax[0].set_xlabel("logS (Jy/beam)")
		ax[0].set_ylabel("dS")
		ax[0].scatter(np.array(s_list),np.array(spull_list),label="TRAIN SAMPLE")
		ax[0].legend()

		ax[1].set_title("Flux Accuracy")
		ax[1].set_xlabel("logS (Jy/beam)")
		ax[1].set_ylabel("dS")
		ax[1].scatter(np.array(s_list_test),np.array(spull_list_test),label="TEST SAMPLE")
		ax[1].legend()

		plt.tight_layout()
		plt.savefig(self.outfile_fluxaccuracy)
		plt.close()
	

		return 0

	###########################
	##     RUN NN 
	###########################
	def run(self):
		""" Run CNN network training """
	
		#===========================
		#==   SET TRAINING DATA
		#===========================
		# - Set training data (generate or read from disk)
		if self.train_data_gen_enabled:
			print("INFO: Generating training data ...")
			status= self.generate_train_data()
			if status<0:
				return -1		
		else:
			print("INFO: Reading training data from disk ...")
			status= self.read_train_data()
			if status<0:
				return -1
			
		
		#===========================
		#==   RUN NN
		#===========================
		if self.do_training:
			#===========================
			#==   BUILD NN
			#===========================
			#- Create the network
			print("INFO: Building network architecture ...")
			status= self.build_network()
			if status<0:
				return -1

			#===========================
			#==   TRAIN NN
			#===========================
			print("INFO: Training network ...")
			status= self.train_network()
			if status<0:
				return -1

			#===========================
			#==   EVALUATE NN
			#===========================
			print("INFO: Evaluating network results ...")
			status= self.evaluate_network()
			if status<0:
				return -1


		return 0
	
##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	#===========================
	#==   PARSE ARGS
	#===========================
	print("INFO: Get script args ...")
	try:
		args= get_args()
	except Exception as ex:
		print("ERROR: Failed to get and parse options (err=%s)",str(ex))
		return 1

	# - Input file
	inputimg= args.inputimg
	filelist_bkg= args.filelist_bkg
	filelist_source= args.filelist_source
	filelist_sourcepars= args.filelist_sourcepars
	nnarcfile= args.nnarcfile

	enable_train_data_generation= False
	if args.generate_train_data:
		enable_train_data_generation= True

	# - Source generation
	marginX= args.marginx
	marginY= args.marginy
	marginX_source= args.marginx_source
	marginY_source= args.marginy_source
	Smin= args.Smin
	Smax= args.Smax
	Smodel= args.Smodel
	Sslope= args.Sslope
	nsources_max= args.nsources_max
	bmaj_min= args.bmaj_min
	bmaj_max= args.bmaj_max
	bmin_min= args.bmin_min
	bmin_max= args.bmin_max
	pa_min= args.pa_min
	pa_max= args.pa_max	
	generate_bkg_from_img= True
	if args.generate_bkg:
		generate_bkg_from_img= False
	bkg_rms= args.bkg_rms
	bkg_mean= args.bkg_mean

	# - Train image generation
	nx= args.nx
	ny= args.ny
	nsamples_bkg= args.nsamples_bkg
	nsamples_source= args.nsamples_source
	nmaxobjects= args.nmaxobjects
	ntargetpars= args.ntargetpars
	test_size= args.test_size	
	flip_train= args.flip_train
	flip_test= args.flip_test

	# - Network architecture
	do_training= True
	if args.no_training:
		do_training= False
	normdatamin= args.normdatamin
	normdatamax= args.normdatamax
	conv_kern_size_min= args.conv_kern_size_min
	conv_kern_size_max= args.conv_kern_size_max
	conv_nfilt_min= args.conv_nfilt_min
	conv_nfilt_max= args.conv_nfilt_max
	conv_pool_kernel_size= args.conv_pool_kernel_size
	dense_size_min= args.dense_size_min
	dense_size_max= args.dense_size_max
	spars_loss_weight= args.spars_loss_weight
	labels_loss_weight= args.labels_loss_weight
	nepochs= args.nepochs
	
	normalize_targets= args.normalize_targets
	optimizer= args.optimizer
	learning_rate= args.learning_rate

	use_standard_nn= args.use_standard_nn
	

	dropout_enabled= True
	if args.no_dropout:
		dropout_enabled= False

	maxpool_enabled= True
	if args.no_maxpool:
		maxpool_enabled= False

	batchnorm_enabled= True
	if args.no_batchnorm:
		batchnorm_enabled= False

	conv_enabled= True
	if args.no_conv:
		conv_enabled= False

	conv_activation= args.conv_activation
	dense_activation= args.dense_activation
	
	# - Output file
	saveimg= args.saveimg
	outfile_loss= args.outfile_loss
	outfile_accuracy= args.outfile_accuracy
	outfile_model= args.outfile_model
	outfile_posaccuracy= args.outfile_posaccuracy
	outfile_fluxaccuracy= args.outfile_fluxaccuracy
	outfile_nnout_train= args.outfile_nnout_train
	outfile_nnout_test= args.outfile_nnout_test
	outfile_nnout_metrics= args.outfile_nnout_metrics

	#===========================
	#==   CHECK ARGS
	#===========================
	# - Check if input file is needed and not given
	if enable_train_data_generation:
		if generate_bkg_from_img and not inputimg:
			print("ERROR: Missing input file argument (needed to generate train data)!")
			return 1
	else:
		if not filelist_bkg:
			print("ERROR: Missing input filelist_bkg argument (needed to read bkg train data)!")
			return 1			
		if not filelist_source:
			print("ERROR: Missing input filelist_source argument (needed to read source train data)!")
			return 1	
		if not filelist_sourcepars:
			print("ERROR: Missing input filelist_sourcepars argument (needed to read source target pars data)!")
			return 1

	# - Check if nnarc file is needed and not given
	if do_training and not use_standard_nn and not nnarcfile:
		print("ERROR: No nn architecture file given in input (needed to build the net for training)!")
		return 1

	#===========================
	#==   RUN NN TRAINING
	#===========================
	print("INFO: Creating and running CNN network ...")
	cnn= CNNTrainer()

	# - Set train data options
	cnn.set_img_filename(inputimg)
	cnn.set_img_bkg_filelist(filelist_bkg)
	cnn.set_img_source_filelist(filelist_source)
	cnn.set_sourcepars_filelist(filelist_sourcepars)
	cnn.set_nnarc_filename(nnarcfile)

	cnn.enable_train_data_generation(enable_train_data_generation)
	cnn.set_margins(marginX,marginY)
	cnn.set_source_margins(marginX_source,marginY_source)
	cnn.set_nsources_max(nsources_max)
	cnn.set_source_flux_range(Smin,Smax)
	cnn.set_source_flux_rand_model(Smodel)
	cnn.set_source_flux_rand_exp_slope(Sslope)
	cnn.set_beam_bmaj_range(bmaj_min,bmaj_max)
	cnn.set_beam_bmin_range(bmin_min,bmin_max)
	cnn.set_beam_pa_range(pa_min,pa_max)
	cnn.enable_bkg_generation_from_img(generate_bkg_from_img)
	cnn.set_gen_bkg_rms(bkg_rms)
	cnn.set_gen_bkg_mean(bkg_mean)

	# - Set NN architecture & train options
	cnn.enable_training(do_training)
	cnn.use_standard_nn(use_standard_nn)
	
	cnn.set_input_data_norm_range(normdatamin,normdatamax)
	cnn.enable_target_normalization(normalize_targets)
	cnn.set_nobjects(nmaxobjects)
	cnn.set_npars(ntargetpars)
	cnn.set_bkg_sample_size(nsamples_bkg)
	cnn.set_source_sample_size(nsamples_source)
	cnn.set_train_img_size(nx,ny)
	cnn.set_test_sample_size(test_size)
	cnn.enable_train_data_flip(flip_train)
	cnn.enable_test_data_flip(flip_test)
	
	cnn.set_conv_kern_size_min(conv_kern_size_min)
	cnn.set_conv_kern_size_max(conv_kern_size_max)
	cnn.set_conv_nfilt_min(conv_nfilt_min)
	cnn.set_conv_nfilt_max(conv_nfilt_max)
	cnn.set_conv_pool_kernel_size(conv_pool_kernel_size)
	cnn.set_conv_activation(conv_activation)
	cnn.set_dense_size_min(dense_size_min)
	cnn.set_dense_size_max(dense_size_max)
	cnn.set_dense_activation(dense_activation)
	cnn.set_spars_loss_weight(spars_loss_weight)
	cnn.set_labels_loss_weight(labels_loss_weight)
	cnn.set_nepochs(nepochs)

	cnn.set_optimizer(optimizer)
	cnn.set_learning_rate(learning_rate)
	cnn.enable_conv(conv_enabled)
	cnn.enable_dropout(dropout_enabled)
	cnn.enable_maxpool(maxpool_enabled)
	cnn.enable_batchnorm(batchnorm_enabled)
	
	# - Set output options
	cnn.save_train_img_to_disk(saveimg)
	cnn.set_outfile_loss(outfile_loss)
	cnn.set_outfile_accuracy(outfile_accuracy)
	cnn.set_outfile_model(outfile_model)
	cnn.set_outfile_posaccuracy(outfile_posaccuracy)
	cnn.set_outfile_fluxaccuracy(outfile_fluxaccuracy)
	cnn.set_outfile_nnout_train(outfile_nnout_train)
	cnn.set_outfile_nnout_test(outfile_nnout_test)
	cnn.set_outfile_nnout_metrics(outfile_nnout_metrics)
	
	# - Run network train
	cnn.run()
	

	

###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	sys.exit(main())

