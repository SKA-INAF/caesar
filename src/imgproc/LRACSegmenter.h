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
* @file LRACSegmenter.h
* @class LRACSegmenter
* @brief Class implementing Local Region-based Active Contour segmentation algorithm
*
* @author S. Lankton, S. Riggi
* @date 15/06/2015
*/

#ifndef _LRAC_SEGMENTER_H
#define _LRAC_SEGMENTER_H 1

#include <Contour.h>

#include <TObject.h>
#include <TMatrixD.h>


#include <stdio.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <float.h>


namespace Caesar{

class Image;

class LRACSegmenter : public TObject {

	public:
  	/**
		\brief Class constructor
		*/
    LRACSegmenter();
		/**
		\brief Class destructor
		*/
    ~LRACSegmenter();
		
		
	public:
		
		/**
		* \brief Find the ChanVese segmentation of input image
		*/
		static Image* FindSegmentation(Image* img,Image* initSegmImg,int iterations=1000,double lambda=0.1,double radius=1,double eps=1.e-2);

	private:

		//pt struct
		struct pt{
    	long x;
    	long y;
    	long z;
    	long idx;
    	struct pt *prev;
    	struct pt *next;
		};
		typedef struct pt PT;

		//List struct
		struct ll{
  		PT *head;
  		PT *curr;
  		long length;
		};
		typedef struct ll LL;	

	private:

		/**
		* \brief Set level set to checker board model
		*/
		static Image* GetCheckerBoardLevelSet(Image* inputImg,double square_size=10);

		/**
		* \brief Set level set to circle model
		*/
		static Image* GetCircleLevelSet(Image* inputImg,double radius_to_image_ratio=0.2);

		//Check convergence	
		

		// functions to minimize lrbac (chanvese) energy
		static void lrbac_chanvese(double *img, double *phi, double *label, long *dims,
                    LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                    int iter, double rad, double lambda,double eps);

		static  double *en_lrbac_compute(LL *Lz,double *phi, double *img, long *dims, double *scale, double lam, double rad);
		//static  void en_lrbac_init(LL *Lz, double *img,double *phi, long *dims, double rad);
		static  void en_lrbac_init(long *dims, double rad);
		static  void en_lrbac_init_point(double* img, double* phi, int idx, int x, int y, int z, long *dims, double rad);
		static  double *en_lrbac_gball(double rad);
		static  void en_lrbac_destroy();
		static  void en_lrbac_update(double* img, long *dims, LL *Lin2out, LL *Lout2in, double rad);
	
		static  void ls_mask2phi3c(double* mask, double* phi, double* label, long* dims, 
                   LL* Lz, LL* Ln1, LL* Ln2, LL* Lp1, LL* Lp2);

		// returns curvature (kappa) at the list point p (x,y,z)
		static double en_kappa_pt(PT* p, double *phi, long *dims);

		// returns curvature (kappa) at p (x,y,z) and sets dx and dy values for the norm
		static double en_kappa_norm_pt(PT* p, double *phi, long *dims, double *dx, double *dy, double *dz);


		//List methods
		static  void ls_iteration(double *F, double *phi, double* label, long* dims, 
                  LL* Lz, LL* Ln1, LL* Lp1, LL *Ln2, LL *Lp2, 
                  LL *Lin2out, LL* Lout2in);
		static void ll_init(LL *list);
		static void ll_step(LL *list);
		static void ll_destroy(LL *list);
		static LL* ll_create();
		static void ll_push(LL *list, PT *add);
		static void ll_pushnew(LL *list, long x, long y, long z, long idx);
		static PT *ll_pop(LL *list);
		static void ll_pop_free(LL *list);
		static void ll_remcurr_free(LL *list);
		static PT *ll_remcurr(LL* list);
		static PT *pt_create(long x, long y, long z, long idx);

		static double ls_min_hood_onlevel(int idx, long x, long y, long z, long *dims, 
                           double *phi,double *label, double level);

		static double ls_max_hood_onlevel(int idx, long x, long y, long z, long *dims, 
                           double *phi,double *label, double level);

	private:

		static double* gball;
		static double* Ain;//local means
		static double* Aout;
		static double* Sin;
		static double* Sout; 
		
		static double uin;// means
		static double uout;
		static double sumin;
		static double sumout;
		static double ain; 
		static double aout;

		static double Fsum;
		static double Q;
		static double dQ;
		
	ClassDef(LRACSegmenter,1)

};//close class
		 

#ifdef __MAKECINT__
#pragma link C++ class LRACSegmenter+;
#endif

}//close namespace 
                    
#endif


