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
* @file ZernikeMoments.cc
* @class ZernikeMoments
* @brief ZernikeMoments
*
* ZernikeMoments class
* @author S. Riggi
* @date 20/01/2015
*/

#include <ZernikeMoments.h>
#include <Image.h>
#include <Logger.h>

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TCut.h>
#include <TEventList.h>
#include <TMath.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TColor.h>

#include <Math/WrappedTF1.h>
#include <Math/GSLIntegrator.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/WrappedFunction.h>
#include <Math/WrappedParamFunction.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>


#include <deque>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/foreach.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/scoped_array.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <time.h>
#include <ctime>
#include <queue>

#include <complex>
#include <cmath>
#include <cfloat> // Has definition of DBL_EPSILON
#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_specfunc.h>

//#include "../../cmatrix.h"


namespace Caesar {



#define MAX_L 32

// This sets the maximum D parameter (15)
// The D parameter has to match MAX_D. See mb_Znl() below.
#define MAX_D 4 //15
// This is based on the maximum D parameter being 15, which sets the number of returned values.
#define MAX_Z 72
// This is also based on the maximum D parameter - contains pre-computed factorials
#define MAX_LUT 240

//ClassImp(ZernikeMoments)

ZernikeMoments::ZernikeMoments(){

}//close costructor


ZernikeMoments::~ZernikeMoments(){

}//close destructor


std::vector<double> ZernikeMoments::GetZernike2D_Direct(Image* img, double order, double radius) {


	//Init radius
	// N is the smaller of (Nx,Ny)
	long int Nx= img->GetNx();
	long int Ny= img->GetNy();
	int N = Nx < Ny ? Nx : Ny;

	double D= MAX_L;//order
	if (order > 0) D = order;
		
	double R= radius;//radius
	//if (radius<0.0) R = N;
	//if (radius<0.0) R= sqrt(Nx*Nx+Ny*Ny)/2.;
	if (radius<0.0) R= sqrt(Nx*Nx+Ny*Ny)/2.;
	//R= N*sqrt(2.);

	//Compute image moments
	double moment10 = 0.0;
	double moment00 = 0.0;
	double moment01 = 0.0;
	double wsum= 0;
	double wmin= 1.e+99;
	double wmax= -1.e+99;
	std::vector<double> xList;
	std::vector<double> yList;
	std::vector<double> wList;

	for (long int i=0;i<Nx;i++) {
		double x= img->GetX(i);
		for (long int j=0;j<Ny;j++) {
			double y= img->GetY(j);
			double w = img->GetPixelValue(i,j);
			if(w==0) continue;
			xList.push_back(x);
			yList.push_back(y);
			wList.push_back(w);
			wsum += w;
			moment10 += x*w;
			moment00 += w;
			moment01 += y*w;
			if(w>wmax) wmax= w;
			if(w<wmin) wmin= w;
		}
	}
	double m00= moment00;
	double m10_m00 = moment10/moment00;
	double m01_m00 = moment01/moment00;
	double Rmax= sqrt(m00)*sqrt(Nx/Ny+Ny/Nx)/2.;

	//Build factorials
	std::vector<double> Fnm;
	std::vector<int> n_s;
	std::vector<int> l_s;
	std::vector< std::complex<double> > zsum;

	for (int n = 0; n <= D; n++) {
		for (int l = 0; l <= n; l++) {
			if ( (n-l)%2==0 ) {
		
				std::stringstream ss;
				ss<<"(n,l)=("<<n<<","<<l<<"), F(";
				for (int m = 0; m <= (n-l)/2; m++) {
					long double a1= pow((double)-1.0,(double)m);
					long double a2= ((long double)gsl_sf_fact(n-m)/((long double)gsl_sf_fact(m)*(long double)gsl_sf_fact((n-2*m+l)/2)*(long double)gsl_sf_fact((n-2*m-l)/2)) );
					long double F= a1*a2;
					Fnm.push_back(F);
					ss<<F<<",";
				}//end loop m
				ss<<")";
				DEBUG_LOG(ss.str());
					
				n_s.push_back(n);
				l_s.push_back(l);
				zsum.push_back( std::complex<double>(0.0,0.0) );
			}//close if
		}//end loop l
	}//end loop n

	//Compute Rnm polinomials
	std::complex<double> Vnl(0.0,0.0);
	double npixInCircle= 0;
	double wSumInCircle= 0;
	double normmin= 0;
	double normmax= 1;

	for (int i=0;i<Nx;i++) {
		// In the paper, the center of the unit circle was the center of the image
		double ix= i;
			
		for (int j=0;j<Ny;j++) {
			double iy= j;
			double p= img->GetBinContent(i+1,j+1);
			double p_norm= normmin + (normmax-normmin)*(p-wmin)/(wmax-wmin);
			if(p==0) continue;			
			//p= p_norm;

			/*
			//Method 1 (Image Center)
			double x= 2.*ix + 1 - Nx;
			double y= 2.*iy + 1 - Ny;
			double r= sqrt(x*x+y*y)/R;
			double theta = atan2((Ny-1-2.*iy),(2.*ix-Nx+1));
			*/
			
			//Method 2 (Image Center of mass)
			double x= img->GetX(i) - m10_m00;
			double y= img->GetY(j) - m01_m00;
			double r= sqrt(x*x+y*y)/Rmax;
			double theta = atan2(y,x);
			
			//Skip pixels outside unit-radius circle
			if (r>1.) continue;

			//DEBUG_LOG("(Nx,Ny)="<<Nx<<","<<Ny<<") M("<<m10_m00<<","<<m01_m00<<") (ix,iy)=("<<ix<<","<<iy<<"), (x,y)=("<<x<<","<<y<<") (x0,y0)="<<x0<<","<<y0<<") r="<<r<<", r0="<<r0<<" theta="<<theta<<" theta0="<<theta0<<" p="<<p<<" R="<<R);

			//Count number of pixels in unit-radius circle
			if(r>0){
				npixInCircle++;
				wSumInCircle+= p;
			}

			int theLUT = 0;
			for (size_t theZ=0;theZ<zsum.size();theZ++) {
				int n = n_s[theZ];
				int l = l_s[theZ];
				Vnl = std::complex<double>(0.0,0.0);
				for(int m=0;m<=(n-l)/2; m++) {	
					double rpow= pow(r,(double)(n-2*m));
					Vnl+= ( std::polar(1.0,l*theta) * Fnm[theLUT] * rpow );
					theLUT++;
				}
				zsum[theZ] += (std::conj(Vnl) * p);
			}//end loop z entries
		}//end loop bins Y
	}//end loop bins X

	npixInCircle+= 1;
	DEBUG_LOG("npixInCircle="<<npixInCircle<<" wSumInCircle="<<wSumInCircle);

	double Re, Im;
	std::vector<double> moments;
	for (size_t theZ=0;theZ<zsum.size();theZ++) {
		int n= n_s[theZ];
		int m= l_s[theZ];
		zsum[theZ]*= (n+1)/TMath::Pi();
		//zsum[theZ]*= (n+1)/npixInCircle;
		double Re = std::real(zsum[theZ]);
		double Im = std::imag (zsum[theZ]);
		double phase= atan2(Im,Re)*TMath::RadToDeg();
		double ampl= fabs(sqrt(Re*Re+Im*Im));
		moments.push_back(ampl);
		DEBUG_LOG("(n,m)=("<<n<<","<<m<<"), Z="<<Re<<"+"<<Im<<"i, A="<<ampl<<" phase="<<phase);
	}

	return moments;

}//close GetZernike2D_Direct()


/* mb_Znl
  Zernike moment generating function.  The moment of degree n and
  angular dependence l for the pixels defined by coordinate vectors
  X and Y and intensity vector P.  X, Y, and P must have the same
  length
*/
std::vector<double> ZernikeMoments::mb_Znl(double *X, double *Y, double *P, int size, double D, double m10_m00, double m01_m00, double R, double psum) {

	static double LUT[MAX_LUT];
	static int n_s[MAX_Z], l_s[MAX_Z];
	static char init_lut=0;

	double x, y, p ;   //individual values of X, Y, P
	int i, theZ, theLUT, numZ=0;
	int n= 0;
	int m= 0;
	int l= 0;
	
	std::complex<double> sum [MAX_Z];
	std::complex<double> Vnl;

	// The LUT indexes don't work unless D == MAX_D
	// To make it more flexible, store the LUT by [m][n][l].  Needs [(D+1)/2][D+1][D+1] of storage.
	// Other hard-coded D values should just need changing MAX_D, MAX_Z and MAX_LUT above.
	assert (D == MAX_D);

	if (!init_lut) {
		theZ=0;
		theLUT=0;
		for (n = 0; n <= MAX_D; n++) {
			for (l = 0; l <= n; l++) {
				if ( (n-l) % 2 == 0 ) {
					for (m = 0; m <= (n-l)/2; m++) {
						LUT[theLUT] = pow((double)-1.0,(double)m) * ( (long double) gsl_sf_fact(n-m) / ( (long double)gsl_sf_fact(m) * (long double)gsl_sf_fact((n - 2*m + l) / 2) *
							(long double)gsl_sf_fact((n - 2*m - l) / 2) ) );
						theLUT++;
					}
					n_s[theZ] = n;
					l_s[theZ] = l;
					theZ++;
				}
			}
		}
		init_lut = 1;
	}

	// Get the number of Z values, and clear the sums.
	for (n = 0; n <= D; n++) {
		for (l = 0; l <= n; l++) {
			if ( (n-l) % 2 == 0 ) {
				sum [numZ] = std::complex<double>(0.0,0.0);
				numZ++;
			}
		}
	}

	// int nfoo=0;
	for(i = 0 ; i < size ; i++) {
		x = (X[i] - m10_m00)/R;
		y = (Y[i] - m01_m00)/R;
		double sqr_x2y2 = sqrt (x*x + y*y);
		if (sqr_x2y2 > 1.0) continue;

		p = P[i] / psum;

		double atan2yx = atan2(y,x);
		theLUT = 0;
		for (theZ = 0; theZ < numZ; theZ++) {
			n = n_s[theZ];
			l = l_s[theZ];
			Vnl = std::complex<double>(0.0,0.0);
			for( m = 0; m <= (n-l)/2; m++ ) {
				Vnl += ( std::polar (1.0, l*atan2yx) * LUT[theLUT] * pow( sqr_x2y2, (double)(n - 2*m)) );
				theLUT++;
			}
			sum [theZ] += (conj(Vnl) * p);
		}
	}

	double preal, pimag;
	std::vector<double> zvalues;
	for (theZ = 0; theZ < numZ; theZ++) {
		sum [theZ] *= ((n_s[theZ]+1)/TMath::Pi());
		preal = std::real ( sum [theZ] );
		pimag = std::imag ( sum [theZ] );
		double phase= atan2(pimag,preal)*TMath::RadToDeg();
		double ampl= fabs(sqrt(preal*preal+pimag*pimag));
		zvalues.push_back(ampl);
		DEBUG_LOG("(n,m)=("<<n<<","<<m<<"), Z="<<preal<<"+"<<pimag<<"i, A="<<ampl<<" phase="<<phase);
	}

	return zvalues;

}//close mb_Znl()


std::vector<double> ZernikeMoments::GetZernike2DOld(Image* img, double D, double R) {
	
	double *Y,*X,*P,psum;
	double intensity;
	int x,y,size;

	long int rows = img->GetNy();
	long int cols = img->GetNx();
	if (D<=0) D=15;
	if (R<=0) R=rows/2;

	Y= new double[rows*cols];
	X= new double[rows*cols];
	P= new double[rows*cols];

	// Find all non-zero pixel coordinates and values 
	size=0;
	psum=0;
	double moment10 = 0.0, moment00 = 0.0, moment01 = 0.0;
	for (y=0;y<rows;y++){
		double binY= img->GetY(y);
		for (x=0;x<cols;x++) {
			double binX= img->GetX(x);

			intensity = img->GetPixelValue(x,y);
			if (intensity != 0) {
				Y[size] = binY;
				X[size] = binX;
				P[size] = intensity;
				psum   += intensity;
				size++;
			}
		
			// moments
			moment10 += binX * intensity;
			moment00 += intensity;
			moment01 += binY * intensity;
		}//end loop bins X
	}//end loop bins Y

  // Normalize the coordinates to the center of mass and normalize
  //pixel distances using the maximum radius argument (R)

	double m10_m00 = moment10/moment00;
	double m01_m00 = moment01/moment00;
	std::vector<double> zvalues= mb_Znl (X,Y,P,size,D,m10_m00,m01_m00,R,psum);

	delete Y;
	delete X;
	delete P;

	return zvalues;

}//close GetZernike2DOld()



std::vector<double> ZernikeMoments::GetZernike2D (Image* image, double order, double rad) {
	
	//zvalue -array of double- a pre-allocated array of double of a suficient size
  //                          (the actual size is returned by "output_size))
  //output_size -* long- the number of enteries in the array "zvalues" (normally 72)

	//Normalize image to range [0,1]
	Image* img= image->GetNormalizedImage(0,1);

	long int Nx= img->GetNx();
	long int Ny= img->GetNy();

	// N is the smaller of I->width and I->height
	int N = Nx < Ny ? Nx : Ny;

	int L= 15;
	if (order > 0) L = (int)order;
	assert (L < MAX_L);

	if (rad > 0.0) rad = N;
	//int D = (int)(rad * 2);
	int D= N;
	//int D= N*sqrt(2);
	double area = TMath::Pi()*rad*rad;
	//double area = TMath::Pi()*D*D/4;

	static double H1[MAX_L][MAX_L];
	static double H2[MAX_L][MAX_L];
	static double H3[MAX_L][MAX_L];
	static char init= 1;

	double COST[MAX_L], SINT[MAX_L], R[MAX_L];
	double Rn, Rnm, Rnm2, Rnnm2, Rnmp2, Rnmp4;

	double a,b, x, y, r, r2, f, const_t;
	double AR[MAX_L][MAX_L], AI[MAX_L][MAX_L];
	
	double sum = 0;
	int cols = Nx;
	int rows = Ny;
	
	// Compute x/0, y/0 and 0/0 moments to center the unit circle on the centroid
	double moment10 = 0.0, moment00 = 0.0, moment01 = 0.0;
	double intensity;
	std::vector<double> xList;
	std::vector<double> yList;
	std::vector<double> wList;
	double wsum= 0;

	for (long int i=0;i<Nx;i++) {
		double x= img->GetX(i);
		for (long int j=0;j<Ny;j++) {
			double y= img->GetY(j);
			double w = img->GetPixelValue(i,j);
			if(w==0) continue;
			xList.push_back(x);
			yList.push_back(y);
			wList.push_back(w);
			wsum += w;
			moment10 += x*w;
			moment00 += w;
			moment01 += y*w;
		}
	}
	double m00= moment00;
	double m10_m00 = moment10/moment00;
	double m01_m00 = moment01/moment00;
	double Rmax= sqrt(m00)*sqrt(Nx/Ny+Ny/Nx)/2.;
	
		
	// Pre-initialization of statics
	if (init) {
		for (int n = 0; n < MAX_L; n++) {
			for (int m = 0; m <= n; m++) {
				if (n != m) {
					H3[n][m] = -(double)(4.0 * (m+2.0) * (m + 1.0) ) / (double)( (n+m+2.0) * (n - m) ) ;
					H2[n][m] = ( (double)(H3[n][m] * (n+m+4.0)*(n-m-2.0)) / (double)(4.0 * (m+3.0)) ) + (m+2.0);
					H1[n][m] = ( (double)((m+4.0)*(m+3.0))/2.0) - ( (m+4.0)*H2[n][m] ) + ( (double)(H3[n][m]*(n+m+6.0)*(n-m-4.0)) / 8.0 );
				}
			}
		}
		init = 0;
	}

	// Zero-out the Zernike moment accumulators
	for (int n = 0; n <= L; n++) {
		for (int m = 0; m <= n; m++) {
			AR[n][m] = AI[n][m] = 0.0;
		}
	}

	
	double npixInCircle= 0;
	double wSumInCircle= 0;


	for (long int i=0;i<Nx;i++) {
		// In the paper, the center of the unit circle was the center of the image
		//x = (double)(2*i+1-N)/(double)D;
		//x = (img->GetXaxis()->GetBinCenter(i+1)-m10_m00)/(double)D;
		x= (img->GetX(i)-m10_m00)/Rmax;

		for (long int j=0;j<Ny;j++) {
			// In the paper, the center of the unit circle was the center of the image
			//y = (double)(2*j+1-N)/(double)D;
			//y= (img->GetYaxis()->GetBinCenter(j+1)-m01_m00)/(double)D;
			y= (img->GetY(j)-m01_m00)/Rmax;
				
			double binContent= img->GetPixelValue(i,j);
			if(binContent==0) continue;
				
			r2 = x*x + y*y;
			r = sqrt (r2); 

			if ( r<DBL_EPSILON || r > 1.0) continue;
				
			//compute all powers of r and save in a table
			R[0] = 1;
			for (int n=1; n <= L; n++) R[n] = r*R[n-1];
				// compute COST SINT and save in tables 
				a = COST[0] = x/r;
				b = SINT[0] = y/r;
				for (int m = 1; m <= L; m++) {
					COST[m] = a * COST[m-1] - b * SINT[m-1];
					SINT[m] = a * SINT[m-1] + b * COST[m-1];
				}

				// compute contribution to Zernike moments for all 
				// orders and repetitions by the pixel at (i,j)
				// In the paper, the intensity was the raw image intensity
				f = binContent/wsum;

				Rnmp2 = Rnm2 = 0;
				for (int n = 0; n <= L; n++) {
					// In the paper, this was divided by the area in pixels
					// seemed that pi was supposed to be the area of a unit circle.
					const_t = (n+1) * f/TMath::Pi();
					Rn = R[n];
					if (n >= 2) Rnm2 = R[n-2];
					for (int m = n; m >= 0; m -= 2) {
						if (m == n) {
							Rnm = Rn;
							Rnmp4 = Rn;
						} 
						else if (m == n-2) {
							Rnnm2 = n*Rn - (n-1)*Rnm2;
							Rnm = Rnnm2;
							Rnmp2 = Rnnm2;
						}	 
						else {
							Rnm = H1[n][m] * Rnmp4 + ( H2[n][m] + (H3[n][m]/r2) ) * Rnmp2;
							Rnmp4 = Rnmp2;
							Rnmp2 = Rnm;
						}
						AR[n][m] += const_t * Rnm * COST[m];
						AI[n][m] -= const_t * Rnm * SINT[m];
					}//end loop m
				}//end loop n
		}//end loop y
	}//end loop x

	int numZ=0;
	std::vector<double> zvalues;
	for (int n = 0; n <= L; n++) {
		for (int m = 0; m <= n; m++) {
			if ( (n-m) % 2 == 0 ) {
				double Re= AR[n][m];
				double Im= AI[n][m];
				double phase= atan2(Im,Re)*TMath::RadToDeg();

				//AR[n][m] *= AR[n][m];
				//AI[n][m] *= AI[n][m];
				double ampl= fabs (sqrt ( Re*Re + Im*Im ));			
				zvalues.push_back(ampl);
				DEBUG_LOG("(n,m)=("<<n<<","<<m<<"), Z="<<Re<<"+"<<Im<<"i, A="<<ampl<<" phase="<<phase);
				numZ++;
			}
		}
	}
	//*output_size = numZ;

	img->Delete();

	return zvalues;

}//close GetZernike2D()

}//close namespace
