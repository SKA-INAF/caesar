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
* @file LRACSegmenter.cc
* @class LRACSegmenter
* @brief Class implementing Local Region-based Active Contour segmentation algorithm
*
* @author S. Lankton, S. Riggi
* @date 15/06/2015
*/

#include <LRACSegmenter.h>
#include <Image.h>
#include <Source.h>
#include <CodeUtils.h>

#include <TObject.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <iostream>

using namespace std;

ClassImp(Caesar::LRACSegmenter)

namespace Caesar {

//Static vars
double* LRACSegmenter::gball= 0;
double* LRACSegmenter::Ain= 0;
double* LRACSegmenter::Aout= 0;
double* LRACSegmenter::Sin= 0;
double* LRACSegmenter::Sout= 0;

double LRACSegmenter::uin= 0;// means
double LRACSegmenter::uout= 0;
double LRACSegmenter::sumin= 0;
double LRACSegmenter::sumout= 0;
double LRACSegmenter::ain= 0; 
double LRACSegmenter::aout= 0;

double LRACSegmenter::Fsum= 0;
double LRACSegmenter::Q= 0;
double LRACSegmenter::dQ= 0;

LRACSegmenter::LRACSegmenter() {

}//close constructor


LRACSegmenter::~LRACSegmenter() {

}//close destructor

Image* LRACSegmenter::GetCheckerBoardLevelSet(Image* inputImg,double square_size)
{
	//Init level set
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	Image* segmImg= inputImg->GetCloned("");
	segmImg->Reset();

	//Fill level set
	for(long int i=0;i<Nx;i++){
    for(long int j=0;j<Ny;j++){
			double xx= sin(TMath::Pi()/square_size*i);
			double yy= sin(TMath::Pi()/square_size*j);
			double w= xx*yy;
			if(w>0) segmImg->SetPixelValue(i,j,1);
			else segmImg->SetPixelValue(i,j,0);
    }//end loop biny
  }//end loop binx

	return segmImg;

}//close GetCheckerBoardLevelSet()

Image* LRACSegmenter::GetCircleLevelSet(Image* inputImg,double radius_to_image_ratio)
{
	//Init level set
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	Image* segmImg= inputImg->GetCloned("");
	segmImg->Reset();

	//Fill level set
	int centerX= Nx/2;
	int centerY= Ny/2;
	double R= std::min(Nx,Ny) * radius_to_image_ratio;

	for(long int i=0;i<Nx;i++){
    for(long int j=0;j<Ny;j++){
			double x= i - centerX;
			double y= j - centerY;
			double w = R - sqrt(x*x + y*y) ;
			if(w>0) segmImg->SetPixelValue(i,j,1);
			else segmImg->SetPixelValue(i,j,0);
    }//end loop biny
  }//end loop binx

	return segmImg;

}//close GetCircleLevelSet()


Image* LRACSegmenter::FindSegmentation(Image* inputImg,Image* inputSegmMap,int niters,double lambda,double radius,double eps){

	//Check inputs
	if(!inputImg){
		ERROR_LOG("Null ptr to input/mask images!");
		return nullptr;
	}

	//Check if image has data
	if(!inputImg->HasPixelData()){
		WARN_LOG("Input image has no pixel data, returning nullptr!");
		return nullptr;
	}
	long int Nx= inputImg->GetNx();
	long int Ny= inputImg->GetNy();
	long int npixels= inputImg->GetNPixels();

	//Check img & mask size	
	if(inputSegmMap && !inputSegmMap->HasSameBinning(inputImg)){
		ERROR_LOG("Initial segmentation map has different binning wrt to input image!");
		return nullptr;
	}

	//Normalize image in [0,255]
	Image* inputImg_norm= inputImg->GetNormalizedImage("LINEAR",0,255);
	if(!inputImg_norm){
		ERROR_LOG("Failed to normalize input image in range [0,255]!");
		return nullptr;
	}	

	//## Init algorithm data
	long dims[5];
	dims[0] = Ny;
	dims[1] = Nx;
	dims[2] = 1; 
	dims[3] = dims[0]*dims[1]; 
	dims[4] = dims[0]*dims[1]*dims[2];
	
	//Init arrays
	double* phi= new double[npixels];	
	double* label= new double[npixels];	
	double* img= new double[npixels];	
	double* mask= new double[npixels];	

	//Create linked lists
  LL* Lz = ll_create();
  LL* Ln1 = ll_create();
  LL* Ln2 = ll_create();
  LL* Lp1 = ll_create();
  LL* Lp2 = ll_create();
  LL* Lin2out = ll_create();
  LL* Lout2in = ll_create();

	//Create level set if not given
	if(!inputSegmMap){
		inputSegmMap= GetCheckerBoardLevelSet(inputImg_norm);
		//inputSegmMap= GetCircleLevelSet(inputImg_norm);
	}

	//Fill arrays
	long int index= 0;

	//if(inputSegmMap){
		for(long int i=0;i<Nx;i++){
			long int ix= i;
			for(long int j=0;j<Ny;j++){
				long int iy= Ny-1-j;
				double w= inputImg_norm->GetPixelValue(ix,iy);
				double segm= inputSegmMap->GetPixelValue(ix,iy);
				img[index]= w;
				mask[index]= segm;
				phi[index]= segm;
				label[index]= segm; 
				//INFO_LOG("img("<<index+1<<")="<<img[index]);
				index++; 
			}//end loop bins y
		}//end loop bins x
	//}//close if
	
	/*
	else{
		
		//## Set up initial circular contour for a 256x256 image
		double rowCenter= Nx/2.0;
		double colCenter= Ny/2.0;
		double initContourRadius= 0.5;//1
		INFO_LOG("Initializing level set from dummy gaussian (center("<<rowCenter<<","<<colCenter<<", radius="<<initContourRadius<<")");
	
		for(long int i=0;i<Nx;i++){
			long int ix= i;
			for(long int j=0;j<Ny;j++){
				long int iy= Ny-1-j;
				double w= inputImg_norm->GetPixelValue(ix,iy);

				double x= double(i) - rowCenter;
				double y= double(j) - colCenter;
				double segm= 900.0/(900.0 + x*x + y*y ) - initContourRadius;
				if(segm>0) segm= 1;	
				else segm= 0;
				
				img[index]= w;
				mask[index]= segm;
				phi[index]= segm;
				label[index]= segm; 
				//INFO_LOG("img("<<index+1<<")="<<img[index]);
				index++; 

			}//end loop bins y
		}//end loop bins x
	}//close else
	*/

  //Initialize lists, phi, and labels
	INFO_LOG("Initialize lists, phi and labels...");
  ls_mask2phi3c(mask,phi,label,dims,Lz,Ln1,Ln2,Lp1,Lp2);
	
	//## Run segmentation
	INFO_LOG("Run LRAC segmentation...");
	//lrbac_chanvese( img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,niters,radius,lambda,plhs,display);
	lrbac_chanvese( img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,niters,radius,lambda,eps);

	//## Get segmentation map
	INFO_LOG("Retrieve segmentation map...");
	Image* segmMap= inputImg->GetCloned("",true,true);
	segmMap->Reset();

	index= 0;
	for(long int i=0;i<Nx;i++){
		long int ix= i;
		for(long int j=0;j<Ny;j++){
			long int iy= Ny-1-j;
			double phiVal= phi[index];	
			if(phiVal<=0){
				segmMap->FillPixel(ix,iy,1.);
			}
			else{
				segmMap->FillPixel(ix,iy,0.);
			}
			index++;
		}//end loop bins y
	}//end loop bins x
	

	//Free data
	DEBUG_LOG("Freeing algorithm data...");
	if(inputImg_norm){
		delete inputImg_norm;
		inputImg_norm;
	}
	if(img){
		delete[] img;
	}
	if(mask){
		delete[] mask;
	}
	if(phi){
		delete[] phi;
	}
	if(label){
		delete[] label;
	}

	//Destroy linked lists
  ll_destroy(Lz);
  ll_destroy(Ln1);
  ll_destroy(Ln2);
  ll_destroy(Lp1);
  ll_destroy(Lp2);
  ll_destroy(Lin2out);
  ll_destroy(Lout2in);

	return segmMap;

}//close FindSegmentation()



void LRACSegmenter::lrbac_chanvese(
	double *img, double *phi, double *label, long *dims,
  LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
  int iter, double rad, double lambda,double eps
)
{
  double *F;
  double scale[1]; scale[0]=0;
  int countdown;

  //initialize datastructures and statistics
  //en_lrbac_init(Lz,img,phi,dims,rad);
	en_lrbac_init(dims,rad);

	int niters_performed= 0;

  for(int i=0;i<iter;i++){
		if(i%100==0) INFO_LOG("Running LRAC algorithm [iter "<<i+1<<"/"<<iter<<", Q="<<Q<<", dQ="<<dQ<<"] ...");

    //compute force
    F = en_lrbac_compute(Lz,phi,img,dims, scale,lambda,rad);
    
		//perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    
		//update statistics
    en_lrbac_update(img, dims, Lin2out, Lout2in, rad);

    //display stuff (maybe)
    //if(display) show_countdown(iter,i,&countdown,plhs);	
		
		niters_performed++;

		//Check convergence
		if(dQ<eps){
			INFO_LOG("Convergence reached (Q="<<Q<<", dQ="<<dQ<<"<"<<eps<<"), stopping algorithm and return segmentation...");
			break;
		}

  }//end iter loop

	//Check if convergence has been reached
	if(niters_performed>=iter){
		WARN_LOG("Convergence was not reached, maximum number of iterations ("<<iter<<") reached...");
	}

  //destroy old datastructures
  en_lrbac_destroy();

}//close lrbac_chanvese()


//void LRACSegmenter::en_lrbac_init(LL *Lz,double *img,double *phi, long *dims, double rad)
void LRACSegmenter::en_lrbac_init(long *dims, double rad)
{
	//Define vars
	long NUMEL= dims[4];

  //create ball
  gball = en_lrbac_gball(rad);
  
  //allocate memory for lookups
  Ain  = (double*)malloc(NUMEL*sizeof(double)); if(Ain==NULL) return;
  Sin  = (double*)malloc(NUMEL*sizeof(double)); if(Sin==NULL) return;
  Aout = (double*)malloc(NUMEL*sizeof(double)); if(Aout==NULL) return;
  Sout = (double*)malloc(NUMEL*sizeof(double)); if(Sout==NULL) return;
  
  //poision "uninitialized" points
  for(int i=0;i<NUMEL;i++){
    Ain[i] = -1; Aout[i] = -1;
  }

}//close en_lrbac_init()


double* LRACSegmenter::en_lrbac_compute(LL *Lz,double *phi, double *img, long *dims, double *scale, double lam, double rad)
{
	//Define vars  
	int x,y,z,idx,n;
  double *F, *kappa;
  double a,Fmax,I;

	double u= 0;
	double v= 0;//NB: not initialized in original code (CHECK!!!)

  // allocate space for F
  F = (double*)malloc(Lz->length*sizeof(double));     if(F==NULL) return NULL;
  kappa = (double*)malloc(Lz->length*sizeof(double)); if(kappa==NULL) return NULL;

	//begining of list
  ll_init(Lz); n=0; Fmax = 0.00001; 

	//Loop through list
  while(Lz->curr != NULL){          
    x = Lz->curr->x; y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    I = img[idx];

    if(Ain[idx] <0){
      en_lrbac_init_point(img,phi,idx,x,y,z,dims,rad);
    }
    if(Ain[idx] >0) u = Sin[idx] /Ain[idx];
    if(Aout[idx]>0) v = Sout[idx]/Aout[idx];
    a = (I-u)*(I-u)-(I-v)*(I-v);
    if(fabs(a)>Fmax) Fmax = fabs(a);
    F[n] = a;
    kappa[n] = en_kappa_pt(Lz->curr, phi, dims); //compute kappa
    ll_step(Lz); n++;       //next point
  }
  if(scale[0]==0) scale[0] = Fmax;

	double Fsum_update= 0;
  for(int j=0;j<Lz->length;j++){
    F[j] = F[j]/scale[0] + lam*kappa[j];
		Fsum_update+= F[j];
  }
	double FsumChange= Fsum_update/Fsum-1;

	//INFO_LOG("Fsum="<<Fsum_update<<" (previous="<<Fsum<<", change="<<FsumChange<<")");
	Fsum= Fsum_update;
	
	//Free data
  free(kappa);

  return F;

}//close en_lrbac_compute()

void LRACSegmenter::en_lrbac_update(double* img, long *dims, LL *Lin2out, LL *Lout2in, double rad)
{
	//Define vars
  int x,y,z,idx;
  int i,j,k,irad,idia,ridx,bidx;

	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];

  irad = (int)(floor(rad));
  idia = irad*2+1;

	//Init list
  ll_init(Lin2out);
  
	while(Lin2out->curr != NULL){
    x = Lin2out->curr->x; y = Lin2out->curr->y; z = Lin2out->curr->z; idx = Lin2out->curr->idx;

    for(i=-irad;i<=irad;i++){
      if((x+i)<0 || (x+i)>=DIMX) continue;
      for(j=-irad;j<=irad;j++){ 
        if((y+j)<0 || (y+j)>=DIMY) continue;
        for(k=-irad;k<=irad;k++){
          if((z+k)<0 || (z+k)>=DIMZ) continue;
          ridx = idx+(i*OFFX)+(j*OFFY)+(k*OFFZ);
          bidx = (j+irad)+((i+irad)*idia)+((k+irad)*idia*idia);

          if(Ain[ridx]>=0)
          {
            Sin[ridx]  -= img[idx]*gball[bidx];
            Ain[ridx]  -= gball[bidx];
            Sout[ridx] += img[idx]*gball[bidx];
            Aout[ridx] += gball[bidx];
          }
        }//end loop z
      }//end loop y
    }//end loop x
    ll_remcurr_free(Lin2out);
  }//end while loop

	//Init list
  ll_init(Lout2in);

  while(Lout2in->curr != NULL){
    x = Lout2in->curr->x; y = Lout2in->curr->y; z = Lout2in->curr->z; idx = Lout2in->curr->idx;
    
    for(i=-irad;i<=irad;i++){
      if((x+i)<0 || (x+i)>=DIMX) continue;
      for(j=-irad;j<=irad;j++){ 
        if((y+j)<0 || (y+j)>=DIMY) continue;
        for(k=-irad;k<=irad;k++){
          if((z+k)<0 || (z+k)>=DIMZ) continue;
          ridx = idx+(i*OFFX)+(j*OFFY)+(k*OFFZ);
          bidx = (j+irad)+((i+irad)*idia)+((k+irad)*idia*idia);

          if(Ain[ridx]>=0)
          {
            Sin[ridx]  += img[idx]*gball[bidx];
            Ain[ridx]  += gball[bidx];
            Sout[ridx] -= img[idx]*gball[bidx];
            Aout[ridx] -= gball[bidx];
          }
        }//end loop z
      }//end loop y
    }//end loop x
    ll_remcurr_free(Lout2in);

  }//end while loop

  if(uin>0)  uin  = sumin/ain;
  if(uout>0) uout = sumout/aout;

}//close en_lrbac_update()

void LRACSegmenter::en_lrbac_init_point(double* img, double* phi, int idx, int x, int y, int z, long *dims, double rad)
{
	//Define vars  
	double usum,vsum,au,av;
  int i,j,k,irad,idia,ridx,bidx;

	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];

  usum= vsum= au= av= 0;
  irad = (int)(floor(rad));
  idia = irad*2+1;

  for(i=-irad;i<=irad;i++){
    if((x+i)<0 || (x+i)>=DIMX) continue;
    for(j=-irad;j<=irad;j++){ 
      if((y+j)<0 || (y+j)>=DIMY) continue;
      for(k=-irad;k<=irad;k++){
        if((z+k)<0 || (z+k)>=DIMZ) continue;
        ridx = idx+(i*OFFX)+(j*OFFY)+(k*OFFZ);
        bidx = (j+irad)+((i+irad)*idia)+((k+irad)*idia*idia);
        
        if(phi[ridx]<=0){
          usum += img[ridx]*gball[bidx];
          au   += gball[bidx];
        }
        else{
          vsum += img[ridx]*gball[bidx];
          av   += gball[bidx];
        }
      }
    }
  }
  Ain[idx] = au;   Aout[idx] = av;
  Sin[idx] = usum; Sout[idx] = vsum;

}//close en_lrbac_init_point()



void LRACSegmenter::ls_mask2phi3c(double* mask, double* phi, double* label, long* dims, LL *Lz, LL *Ln1, LL *Ln2, LL *Lp1, LL *Lp2)
{
	//Define vars
  int x,y,z,idx;
  int i,j,k;
  int u,d,r,l,f,b;
  int  flag=0;

	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long DIMXY= dims[3];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];
	
  //find 'interface' and mark as 0, create Lz
  for(x=0;x<DIMX;x++) for(y=0;y<DIMY;y++) for(z=0;z<DIMZ;z++){
    idx = (int)(z*DIMXY+x*DIMY+y); 

    //mark the inside and outside of label and phi
    if(mask[idx]==0){ label[idx] =  3; phi[idx]= 3; }
    else            { label[idx] = -3; phi[idx]=-3; }

    if(mask[idx] == 1){
      flag = 0;
      //if any neighbors are 1;
      if(((y+1)<DIMY) && mask[idx+OFFY]==0){flag = 1;}//up
      if(((y-1)>=0)   && mask[idx-OFFY]==0){flag = 1;}//down
      if(((x+1)<DIMX) && mask[idx+OFFX]==0){flag = 1;}//right
      if(((x-1)>=0)   && mask[idx-OFFX]==0){flag = 1;}//left
      if(((z+1)<DIMZ) && mask[idx+OFFZ]==0){flag = 1;}//front
      if(((z-1)>=0)   && mask[idx-OFFZ]==0){flag = 1;}//back
      if(flag){
        ll_pushnew(Lz,x,y,z,idx);
        label[idx] = 0;
        phi[idx] = 0;
      }
    }
  }

  //scan Lz to create Ln1 and Lp1
  ll_init(Lz);
  while(Lz->curr != NULL){
    x = Lz->curr->x; y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    
    if(((y+1)<DIMY) && abs((int)label[idx+OFFY])==3){//up
      if(phi[idx+OFFY]<0){//in
        label[idx+OFFY]=-1; phi[idx+OFFY]=-1;
        ll_pushnew(Ln1,x,y+1,z,idx+OFFY);
      }
      else{                //out
        label[idx+OFFY]=1; phi[idx+OFFY]=1;
        ll_pushnew(Lp1,x,y+1,z,idx+OFFY);
      }
    }
    if(((y-1)>=0)   && abs((int)label[idx-OFFY])==3){//down
      if(phi[idx-OFFY]<0){ //in
        label[idx-OFFY]=-1; phi[idx-OFFY]=-1;
        ll_pushnew(Ln1,x,y-1,z,idx-OFFY);
      }
      else{                //out
        label[idx-OFFY]=1; phi[idx-OFFY]=1;
        ll_pushnew(Lp1,x,y-1,z,idx-OFFY);
      }
    }
    if(((x+1)<DIMX) && abs((int)label[idx+OFFX])==3){//right
      if(phi[idx+OFFX]<0){//in
        label[idx+OFFX]=-1; phi[idx+OFFX]=-1;
        ll_pushnew(Ln1,x+1,y,z,idx+OFFX);
      }
      else{                //out
        label[idx+OFFX]=1; phi[idx+OFFX]=1;
        ll_pushnew(Lp1,x+1,y,z,idx+OFFX);
      }
    }
    if(((x-1)>=0)   && abs((int)label[idx-OFFX])==3){//left
      if(phi[idx-OFFX]<0){ //in
        label[idx-OFFX]=-1; phi[idx-OFFX]=-1;
        ll_pushnew(Ln1,x-1,y,z,idx-OFFX);
      }
      else{                //out
        label[idx-OFFX]=1; phi[idx-OFFX]=1;
        ll_pushnew(Lp1,x-1,y,z,idx-OFFX);
      }
    }
    if(((z+1)<DIMZ) && abs((int)label[idx+OFFZ])==3){//front
      if(phi[idx+OFFZ]<0){//in
        label[idx+OFFZ]=-1; phi[idx+OFFZ]=-1;
        ll_pushnew(Ln1,x,y,z+1,idx+OFFZ);
      }
      else{                //out
        label[idx+OFFZ]=1; phi[idx+OFFZ]=1;
        ll_pushnew(Lp1,x,y,z+1,idx+OFFZ);
      }
    }
    if(((z-1)>=0) && abs((int)label[idx-OFFZ])==3 ){//back
      if(phi[idx-OFFZ]<0){ //in
        label[idx-OFFZ]=-1; phi[idx-OFFZ]=-1;
        ll_pushnew(Ln1,x,y,z-1,idx-OFFZ);
      }
      else{                //out
        label[idx-OFFZ]=1; phi[idx-OFFZ]=1;
        ll_pushnew(Lp1,x,y,z-1,idx-OFFZ);
      }
    }

    ll_step(Lz);
  }


  //scan Ln1 to create Ln2
  ll_init(Ln1);
  while(Ln1->curr != NULL){
    x = Ln1->curr->x; y = Ln1->curr->y; z = Ln1->curr->z; idx = Ln1->curr->idx;
    
    if(((y+1)<DIMY) && label[idx+OFFY]==-3){//up
        label[idx+OFFY]=-2; phi[idx+OFFY]=-2;
        ll_pushnew(Ln2,x,y+1,z,idx+OFFY);
    }
    if(((y-1)>=0) && label[idx-OFFY]==-3){//down
        label[idx-OFFY]=-2; phi[idx-OFFY]=-2;
        ll_pushnew(Ln2,x,y-1,z,idx-OFFY);
    }
    if(((x+1)<DIMX) && label[idx+OFFX]==-3){//right
        label[idx+OFFX]=-2; phi[idx+OFFX]=-2;
        ll_pushnew(Ln2,x+1,y,z,idx+OFFX);
    }
    if(((x-1)>=0) && label[idx-OFFX]==-3){//left
        label[idx-OFFX]=-2; phi[idx-OFFX]=-2;
        ll_pushnew(Ln2,x-1,y,z,idx-OFFX);
    }
    if(((z+1)<DIMZ) && label[idx+OFFZ]==-3){//front
        label[idx+OFFZ]=-2; phi[idx+OFFZ]=-2;
        ll_pushnew(Ln2,x,y,z+1,idx+OFFZ);
    }
    if(((z-1)>=0) && label[idx-OFFZ]==-3){//back
        label[idx-OFFZ]=-2; phi[idx-OFFZ]=-2;
        ll_pushnew(Ln2,x,y,z-1,idx-OFFZ);
    }
    ll_step(Ln1);
  }

  //scan Lp1 to create Lp2
  ll_init(Lp1);
  while(Lp1->curr != NULL){
    x = Lp1->curr->x; y = Lp1->curr->y; z = Lp1->curr->z; idx = Lp1->curr->idx;
    
    if(((y+1)<DIMY) && label[idx+OFFY]==3){//up
        label[idx+OFFY]=2; phi[idx+OFFY]=2;
        ll_pushnew(Lp2,x,y+1,z,idx+OFFY);
    }
    if(((y-1)>=0) && label[idx-OFFY]==3){//down
        label[idx-OFFY]=2; phi[idx-OFFY]=2;
        ll_pushnew(Lp2,x,y-1,z,idx-OFFY);
    }
    if(((x+1)<DIMX) && label[idx+OFFX]==3){//right
        label[idx+OFFX]=2; phi[idx+OFFX]=2;
        ll_pushnew(Lp2,x+1,y,z,idx+OFFX);
    }
    if(((x-1)>=0) && label[idx-OFFX]==3){//left
        label[idx-OFFX]=2; phi[idx-OFFX]=2;
        ll_pushnew(Lp2,x-1,y,z,idx-OFFX);
    }
    if(((z+1)<DIMZ) && label[idx+OFFZ]==3){//front
        label[idx+OFFZ]=2; phi[idx+OFFZ]=2;
        ll_pushnew(Lp2,x,y,z+1,idx+OFFZ);
    }
    if(((z-1)>=0) && label[idx-OFFZ]==3){//back
        label[idx-OFFZ]=2; phi[idx-OFFZ]=2;
        ll_pushnew(Lp2,x,y,z-1,idx-OFFZ);
    }
    ll_step(Lp1);
  }

}//close ls_mask2phi3c()


// allocates and populates memory for a 3D Gaussian, 
// size (floor(rad)*2+1)^3 centered in the middle with sigma = rad/2.
double* LRACSegmenter::en_lrbac_gball(double rad)
{
  //Define vars
	double* gball;
  int dia,dia2,i,j,k,idx;
  double cen,x2,y2,z2,sig2;
  double gsum;
  dia = (int)(floor(rad)*2+1);
  dia2 = dia*dia;
  cen = (int)(floor(rad));
  sig2 = (rad/2)*(rad/2);

  gball = (double*)malloc(sizeof(double)*dia*dia*dia);
  if(gball == NULL) return NULL;

  gsum = 0;
  for(i=0;i<dia;i++){
    for(j=0;j<dia;j++){
      for(k=0;k<dia;k++){
        idx = i+j*dia+k*dia2;
        x2 = ((double)i-cen)*((double)i-cen);
        y2 = ((double)j-cen)*((double)j-cen);
        z2 = ((double)k-cen)*((double)k-cen);
        gball[idx] = exp(-(x2+y2+z2)/(2*sig2));
        gsum += gball[idx];
      }
    }
  }

  for(i=0;i<(dia2*dia);i++){
    gball[i] = gball[i]/gsum;
  }

  return gball;

}//close en_lrbac_gball()

void LRACSegmenter::en_lrbac_destroy(){
  if(gball!=NULL) free(gball);
  if(Ain!=NULL) free(Ain); if(Aout!=NULL) free(Aout);
  if(Sin!=NULL) free(Sin); if(Sout!=NULL) free(Sout);
}


double LRACSegmenter::en_kappa_pt(PT* p, double *phi, long *dims){
  double dx,dy,dz;
  return en_kappa_norm_pt(p, phi, dims, &dx, &dy, &dz);
}

double LRACSegmenter::en_kappa_norm_pt(PT* p, double *phi, long *dims, double *pdx, double *pdy, double *pdz)
{
	//Define vars
  double kappa;
  double dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz,dx2,dy2,dz2;
  int idx,x,y,z,n;
  int xok,yok,zok;

	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];

  x = p->x; y = p->y; z = p->z; idx = p->idx;

  dx=dy=dz=dxx=dyy=dzz=dxy=dxz=dyz=dx2=dy2=dz2=0;
  xok = yok = zok = 0;

  if((x+1)<DIMX && (x-1)>=0) xok = 1;
  if((y+1)<DIMY && (y-1)>=0) yok = 1;
  if((z+1)<DIMZ && (z-1)>=0) zok = 1;

  if(xok){
    dx  = (phi[idx-OFFX]-phi[idx+OFFX])/2;// (l-r)/2
    dxx = (phi[idx-OFFX]-2*phi[idx]+phi[idx+OFFX]); // l-2c+r
    dx2 = dx*dx;
  }
  if(yok){
    dy  = (phi[idx-OFFY]-phi[idx+OFFY])/2;// (u-d)/2
    dyy = (phi[idx-OFFY]-2*phi[idx]+phi[idx+OFFY]);// u-2c+d
    dy2 = dy*dy;
  }
  if(zok){
    dz  = (phi[idx-OFFZ]-phi[idx+OFFZ])/2;// (b-f)/2
    dzz = (phi[idx-OFFZ]-2*phi[idx]+phi[idx+OFFZ]);// b-2c+f
    dz2 = dz*dz;
  }
  if(xok && yok){// (ul+dr-ur-dl)/4
    dxy = (phi[idx-OFFY-OFFX]+phi[idx+OFFY+OFFX]-phi[idx-OFFY+OFFX]-phi[idx+OFFY-OFFX])/4;
  }
  if(xok && zok){// (lf+rb-rf-lb)/4
    dxz = (phi[idx+OFFZ-OFFX]+phi[idx-OFFZ+OFFX]-phi[idx+OFFZ+OFFX]-phi[idx-OFFZ-OFFX])/4;
  }
  if(zok && yok){// (uf+db-df-ub)/4
    dyz = (phi[idx-OFFY+OFFZ]+phi[idx+OFFY-OFFZ]-phi[idx+OFFY+OFFZ]-phi[idx-OFFY-OFFZ])/4;
  }
  
  kappa = (dxx*(dy2+dz2)+dyy*(dx2+dz2)+dzz*(dx2+dy2)-
           2*dx*dy*dxy-2*dx*dz*dxz-2*dy*dz*dyz)/
           (dx2+dy2+dz2+.00000001);

  pdx[0] = dx;
  pdy[0] = dy;
  pdz[0] = dz;

  return kappa;

}//close en_kappa_norm_pt()

//=======================================
//==    LIST FUNCTIONS
//=======================================
void LRACSegmenter::ll_init(LL *list){
  if(list == NULL) return;
  list->curr = list->head;
}

void LRACSegmenter::ll_step(LL *list){
  if(list == NULL) return;
  if(list->curr != NULL){
    list->curr = list->curr->next;
  }
}

LRACSegmenter::LL* LRACSegmenter::ll_create(){
  LL *newll = (LL*)malloc(sizeof(LL));
  if(newll == NULL) return NULL;
  newll->head = NULL;
  newll->curr = NULL;
  newll->length = 0;
  return newll;
}

void LRACSegmenter::ll_destroy(LL *list){
  if(list==NULL) return;
  while(list->head != NULL){
    ll_pop_free(list);
  }
  free(list);
}

void LRACSegmenter::ll_push(LL *list, PT *add){
  if(add == NULL || list == NULL) return;
  add->next = list->head;
  add->prev = NULL;
  if(add->next != NULL){
    add->next->prev = add;
  }
  list->head = add;
  list->length++;
}

void LRACSegmenter::ll_pushnew(LL *list, long x, long y, long z, long idx){
  if(list == NULL) return;
  PT* add = pt_create(x,y,z,idx);
  if(add == NULL) return;
  ll_push(list,add);
}

LRACSegmenter::PT* LRACSegmenter::ll_pop(LL *list){
  if(list == NULL) return NULL;
  PT *out = list->head;
  if(out != NULL){
    list->head = out->next;
    if(list->curr == out) list->curr = list->head;
    if(list->head != NULL) list->head->prev = NULL;
    list->length--;
  }
  return out;
}

void LRACSegmenter::ll_pop_free(LL *list){
  PT* p = ll_pop(list);
  if(p != NULL) free(p);
}

LRACSegmenter::PT* LRACSegmenter::pt_create(long x, long y, long z, long idx){
  PT* newpt = (PT*)malloc(sizeof(PT));
  if(newpt == NULL) return NULL;

  newpt->x = x;
  newpt->y = y;
  newpt->z = z;
  newpt->idx = idx;
  newpt->prev = NULL;
  newpt->next = NULL;
  return newpt;
}

void LRACSegmenter::ll_remcurr_free(LL *list){
  PT* p = ll_remcurr(list);
  if(p != NULL) free(p);
}

LRACSegmenter::PT* LRACSegmenter::ll_remcurr(LL *list){
  if(list == NULL) return NULL;
  PT* out = list->curr;
  if(out == list->head){
    return ll_pop(list);
  }
  else
  {
    if(out != NULL)
    {
      if(out->next != NULL){
        out->next->prev = out->prev;
      }
      if(out->prev != NULL){
        out->prev->next = out->next;
      }
      list->curr = out->next;
      list->length--;
    }
    return out;
  }
}//close ll_remcurr


void LRACSegmenter::ls_iteration(double *F, double *phi, double* label, long* dims, 
                  LL* Lz, LL* Ln1, LL* Lp1, LL *Ln2, LL *Lp2, 
                  LL *Lin2out, LL* Lout2in)
{
	//Define vars
  int x,y,z,i,idx;
  int u,d,r,l,f,b;
  double p, phi_old;
  LL *Sz, *Sn1, *Sp1, *Sn2, *Sp2;

	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];

	//Copy current phi to compute Q
	long NPIXELS= DIMX*DIMY;
	double* phi_prev= new double[NPIXELS];
	for(long t=0;t<NPIXELS;t++) phi_prev[t]= phi[t];

  // create 'changing status' lists
  Sz  = ll_create();
  Sn1 = ll_create();
  Sp1 = ll_create();
  Sn2 = ll_create();
  Sp2 = ll_create();

  // #1) Normalize F
  double Fmax = .001;
  for(i=0;i<Lz->length;i++){
    if(fabs(F[i])>Fmax) Fmax = fabs(F[i]);
  }

  for(i=0;i<Lz->length;i++){ 
    F[i] = F[i]/Fmax*0.4;
  }

  // #2) add F to phi(Lz), create Sn1 & Sp1 
  //                                             ========
  //     (a) scan Lz values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Lz); i = 0;
  while(Lz->curr != NULL){
    x= Lz->curr->x; y= Lz->curr->y; z= Lz->curr->z; idx= Lz->curr->idx;
    phi_old = phi[idx];

    phi[idx] = phi[idx]+F[i];

		

    //check to see if point crossed interface
    if(phi_old<=0 && phi[idx]>0 ){
      ll_pushnew(Lin2out,x,y,z,idx);
    }
    if(phi_old>0  && phi[idx]<=0){
      ll_pushnew(Lout2in,x,y,z,idx);
    }
    
    if(phi[idx] > .5){
      ll_push(Sp1, ll_remcurr(Lz)); 
    }
    else if(phi[idx] < -.5){
      ll_push(Sn1, ll_remcurr(Lz));
    }
    else{
      ll_step(Lz);
    }
    i++; //increment index into F
  }
  if(F!= NULL) free(F); // free F (no longer needed);

  // #3) update Ln1,Ln2,Lp1,Lp2
  //                                    ==========
  //     (c) scan Ln1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Ln1);
  while(Ln1->curr != NULL){
    x= Ln1->curr->x; y= Ln1->curr->y; z= Ln1->curr->z; idx= Ln1->curr->idx;
    p = ls_max_hood_onlevel(idx, x, y, z, dims, phi,label,0);
    if(p>=-0.5){        // found something
      phi[idx] = p-1;

      if(phi[idx]>=-0.5){
        ll_push(Sz,ll_remcurr(Ln1));
      }
      else if(phi[idx]<-1.5){
        ll_push(Sn2,ll_remcurr(Ln1));
      }
      else ll_step(Ln1);
    }
    else{
      ll_push(Sn2,ll_remcurr(Ln1));
    }
  }

  //                                                      ========
  //     (c) scan Lp1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Lp1);
  while(Lp1->curr != NULL){
    x= Lp1->curr->x; y= Lp1->curr->y; z= Lp1->curr->z; idx= Lp1->curr->idx;
    p = ls_min_hood_onlevel(idx, x, y, z, dims, phi,label,0);
    if(p<=0.5){         // found something
      phi_old = phi[idx];
      phi[idx] = p+1;

      if(phi[idx]<=0.5){
        ll_push(Sz,ll_remcurr(Lp1));
      }
      else if(phi[idx]>1.5){
        ll_push(Sp2,ll_remcurr(Lp1));
      }
      else ll_step(Lp1);
    }
    else{
      ll_push(Sp2,ll_remcurr(Lp1));
    }
  }

  //                         ===========
  //     (c) scan Ln2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Ln2);
  while(Ln2->curr != NULL){
    x = Ln2->curr->x; y = Ln2->curr->y; z = Ln2->curr->z; idx= Ln2->curr->idx;
    p = ls_max_hood_onlevel(idx, x, y, z, dims, phi,label,-1);
    if(p>=-1.5){         // found something
      phi[idx] = p-1;
      if(phi[idx]>=-1.5){
        ll_push(Sn1,ll_remcurr(Ln2));
      }
      else if(phi[idx]<-2.5){
        ll_remcurr_free(Ln2);
        phi[idx] = -3; label[idx] = -3;
      }
      else ll_step(Ln2);
    }
    else{
      ll_remcurr_free(Ln2);
      phi[idx] = -3; label[idx] = -3;
    }
  }

  //                                                              =========
  //     (d) scan Lp2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Lp2);
  while(Lp2->curr != NULL){
    x = Lp2->curr->x; y = Lp2->curr->y; z = Lp2->curr->z; idx= Lp2->curr->idx;
    p = ls_min_hood_onlevel(idx, x, y, z, dims, phi,label,1);
    if(p<=1.5){         // found something
      phi[idx] = p+1;
      if(phi[idx]<=1.5){
        ll_push(Sp1,ll_remcurr(Lp2));
      }
      else if(phi[idx]>2.5){
        ll_remcurr_free(Lp2);
        phi[idx] = 3; label[idx] = 3;
      }
      else ll_step(Lp2);
    }
    else{
      ll_remcurr_free(Lp2);
      phi[idx] = 3; label[idx] = 3;
    }
  }

  // #4) Deal with S-lists Sz,Sn1,Sp1,Sn2,Sp2
  //     (a) Scan Sz
  ll_init(Sz);
  while(Sz->curr != NULL){
    idx= Sz->curr->idx;
    ll_push(Lz,ll_remcurr(Sz));
    label[idx] = 0;
  }

  //     (b) Scan Sn1
  ll_init(Sn1);
  while(Sn1->curr != NULL){
    x = Sn1->curr->x; y = Sn1->curr->y; z = Sn1->curr->z; idx = Sn1->curr->idx;
    ll_push(Ln1,ll_remcurr(Sn1));
    label[idx] = -1;
    if(((y+1)<DIMY) && phi[idx+OFFY]==-3){
      ll_pushnew(Sn2,x,y+1,z,idx+OFFY); 
      phi[idx+OFFY] = phi[idx] - 1; 
    }//up
    if(((y-1)>=0)   && phi[idx-OFFY]==-3){
      ll_pushnew(Sn2,x,y-1,z,idx-OFFY);
      phi[idx-OFFY] = phi[idx] - 1; 
    }//down
    if(((x+1)<DIMX) && phi[idx+OFFX]==-3){
      ll_pushnew(Sn2,x+1,y,z,idx+OFFX);
      phi[idx+OFFX] = phi[idx] - 1; 
    }//right
    if(((x-1)>=0)   && phi[idx-OFFX]==-3){
      ll_pushnew(Sn2,x-1,y,z,idx-OFFX);
      phi[idx-OFFX] = phi[idx] - 1; 
    }//left
    if(((z+1)<DIMZ) && phi[idx+OFFZ]==-3){
      ll_pushnew(Sn2,x,y,z+1,idx+OFFZ);
      phi[idx+OFFZ] = phi[idx] - 1; 
    }//front
    if(((z-1)>=0)   && phi[idx-OFFZ]==-3){
      ll_pushnew(Sn2,x,y,z-1,idx-OFFZ);
      phi[idx-OFFZ] = phi[idx] - 1; 
    }//back
  }

  //     (c) Scan Sp1
  ll_init(Sp1);
  while(Sp1->curr != NULL){
    x = Sp1->curr->x; y = Sp1->curr->y; z = Sp1->curr->z; idx=Sp1->curr->idx;
    ll_push(Lp1,ll_remcurr(Sp1));
    label[idx] = 1;
    if(((y+1)<DIMY) && phi[idx+OFFY]==3){
      ll_pushnew(Sp2,x,y+1,z,idx+OFFY); 
      phi[idx+OFFY] = phi[idx] + 1;
    }//up
    if(((y-1)>=0)   && phi[idx-OFFY]==3){
      ll_pushnew(Sp2,x,y-1,z,idx-OFFY); 
      phi[idx-OFFY] = phi[idx] + 1;
    }//down
    if(((x+1)<DIMX) && phi[idx+OFFX]==3){
      ll_pushnew(Sp2,x+1,y,z,idx+OFFX); 
      phi[idx+OFFX] = phi[idx] + 1;
    }//right
    if(((x-1)>=0)   && phi[idx-OFFX]==3){
      ll_pushnew(Sp2,x-1,y,z,idx-OFFX); 
      phi[idx-OFFX] = phi[idx] + 1;
    }//left
    if(((z+1)<DIMZ) && phi[idx+OFFZ]==3){
      ll_pushnew(Sp2,x,y,z+1,idx+OFFZ); 
      phi[idx+OFFZ] = phi[idx] + 1;
    }//front
    if(((z-1)>=0)   && phi[idx-OFFZ]==3){
      ll_pushnew(Sp2,x,y,z-1,idx-OFFZ); 
      phi[idx-OFFZ] = phi[idx] + 1;
    }//back
  }

  //     (d) Scan Sn2
  ll_init(Sn2);
  while(Sn2->curr != NULL){
    idx = Sn2->curr->idx;
    ll_push(Ln2,ll_remcurr(Sn2));
    label[idx] = -2;
  }

  //     (e) Scan Sp2
  ll_init(Sp2);
  while(Sp2->curr != NULL){
    idx = Sp2->curr->idx;
    ll_push(Lp2,ll_remcurr(Sp2));
    label[idx] = 2;
  }


	//Compute Q
	double Q_old= Q;
	for(long t=0;t<NPIXELS;t++){
		double w= 0;
		double w_old= 0;
		if(phi[t]<0) w= 1;		
		if(phi_prev[t]<0) w_old= 1;

		double segm_diff= fabs(w-w_old);
		Q+= segm_diff;
		
		//double phi_diff= fabs(phi[t]-phi_prev[t]);
		//Q+= phi_diff;
	}
	dQ= fabs(Q/Q_old-1.);

	//Clear up
  ll_destroy(Sz);
  ll_destroy(Sn1);
  ll_destroy(Sp1);
  ll_destroy(Sn2);
  ll_destroy(Sp2);

	if(phi_prev){
		delete[] phi_prev;
	}

}//close ls_iteration()


double LRACSegmenter::ls_max_hood_onlevel(int idx, long x, long y, long z, long *dims, double *phi, double *label, double level)
{
	//Define vars
  double pmax = -3;
	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];

  if(((y+1)<DIMY) && (label[idx+OFFY]>=level) && (phi[idx+OFFY]>pmax)) pmax = phi[idx+OFFY]; 
  if(((y-1)>=0  ) && (label[idx-OFFY]>=level) && (phi[idx-OFFY]>pmax)) pmax = phi[idx-OFFY]; 
  if(((x+1)<DIMX) && (label[idx+OFFX]>=level) && (phi[idx+OFFX]>pmax)) pmax = phi[idx+OFFX]; 
  if(((x-1)>=0  ) && (label[idx-OFFX]>=level) && (phi[idx-OFFX]>pmax)) pmax = phi[idx-OFFX]; 
  if(((z+1)<DIMZ) && (label[idx+OFFZ]>=level) && (phi[idx+OFFZ]>pmax)) pmax = phi[idx+OFFZ]; 
  if(((z-1)>=0  ) && (label[idx-OFFZ]>=level) && (phi[idx-OFFZ]>pmax)) pmax = phi[idx-OFFZ]; 

  return pmax;

}//close ls_max_hood_onlevel()

double LRACSegmenter::ls_min_hood_onlevel(int idx, long x, long y, long z, long *dims, double *phi, double *label, double level)
{
	//Define vars
  double pmin = 3;
	long DIMX= dims[1];
	long DIMY= dims[0];
	long DIMZ= dims[2];
	long OFFX= dims[0];
	long OFFY= 1;
	long OFFZ= dims[3];

  if(((y+1)<DIMY) && (label[idx+OFFY]<=level) && (phi[idx+OFFY]<pmin)) pmin = phi[idx+OFFY]; 
  if(((y-1)>=0  ) && (label[idx-OFFY]<=level) && (phi[idx-OFFY]<pmin)) pmin = phi[idx-OFFY]; 
  if(((x+1)<DIMX) && (label[idx+OFFX]<=level) && (phi[idx+OFFX]<pmin)) pmin = phi[idx+OFFX]; 
  if(((x-1)>=0  ) && (label[idx-OFFX]<=level) && (phi[idx-OFFX]<pmin)) pmin = phi[idx-OFFX]; 
  if(((z+1)<DIMZ) && (label[idx+OFFZ]<=level) && (phi[idx+OFFZ]<pmin)) pmin = phi[idx+OFFZ]; 
  if(((z-1)>=0  ) && (label[idx-OFFZ]<=level) && (phi[idx-OFFZ]<pmin)) pmin = phi[idx-OFFZ]; 

  return pmin;

}//close ls_min_hood_onlevel()

}//close namespace

