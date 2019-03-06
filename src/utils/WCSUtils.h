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
* @file WCSUtils.h
* @class WCSUtils
* @brief Utility functions for world coordinate system tasks
*
* Utility functions for world coordinate system tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef _WCS_UTILS_h
#define _WCS_UTILS_h 1


#include <TObject.h>
#include <TMath.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

#include <vector>
#include <algorithm>
#include <map>
#include <time.h>
#include <ctime>

using namespace std;


namespace Caesar {

class ImgMetaData;
class WCS;

//========================================
//==   WCS UTILS
//========================================
class WCSUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    WCSUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~WCSUtils();

		
	public:

		/**
		* \brief Return WCS from image metadata 
		*/
		static WCS* ComputeWCSFromImgMetaData(ImgMetaData* metadata,int coordSystem=-1);

		/**
		* \brief Get WCS coordinates corresponding to image coordinates
		*/
		static int PixelToWCSCoords(double& xpos, double& ypos,WCS* wcs,double ix,double iy);
		/**
		* \brief Get WCS coordinates in string format corresponding to image coordinates
		*/
		static int PixelToWCSStrCoords(std::string& wcs_str,WCS* wcs,double ix,double iy,int max_str_length=4096);

		/**
		* \brief Delete WCS class
		*/
		static void DeleteWCS(WCS** wcs);

	public:
		/**
		* \brief Initialize WCS output coordinate system for use by PIX2WCS
		*/
		static void wcsoutinit(WCS* wcs,char* coorsys);

		/**
		* \brief Return current value of WCS output coordinate system set by -wcsout
		*/
		static char* getwcsout(WCS* wcs);

		/**
		* \brief Return current value of WCS output coordinate system set by -wcsout
		*/
		static std::string GetWCSTypeStr(WCS* wcs);

		/**
		* \brief Convert pixel coordinates to World Coordinate string 
		*/
		static int pix2wcst(WCS* wcs,double xpix,double ypix,char *wcstring,int lstr);

		/**
		* \brief Convert pixel coordinates to World Coordinates
		*/
		static void pix2wcs(WCS* wcs,double xpix,double ypix,double *xpos,double *ypos);	

		/**
		* \brief Convert World Coordinates to pixel coordinates
		*/
		static void wcs2pix(WCS* wcs,double xpos,double ypos,double *xpix,double *ypix,int *offscl);

		/**
		* \brief Set up WCS structure from keyword values
		*/
		static WCS* wcskinit( 
			int naxis1,	/* Number of pixels along x-axis */
			int naxis2,	/* Number of pixels along y-axis */
			char *ctype1,	/* FITS WCS projection for axis 1 */
			char *ctype2,	/* FITS WCS projection for axis 2 */
			double crpix1,	/* Reference pixel coordinates */
			double crpix2,	/* Reference pixel coordinates */
			double crval1,	/* Coordinate at reference pixel in degrees */
			double crval2,	/* Coordinate at reference pixel in degrees */
			double *cd,	/* Rotation matrix, used if not NULL */
			double cdelt1,	/* scale in degrees/pixel, if cd is NULL */
			double cdelt2,	/* scale in degrees/pixel, if cd is NULL */
			double crota,	/* Rotation angle in degrees, if cd is NULL */
			int equinox,	/* Equinox of coordinates, 1950 and 2000 supported */
			double epoch /* Epoch of coordinates, for FK4/FK5 conversion */
		);	

		/**
		* \brief Initialize WCS input coordinate system for use by WCS2PIX
		*/
		static void wcsininit(WCS* wcs,char* coorsys);


	protected:

		//========================================
		//==   WCS enums
		//========================================
		enum WCSType {
			eWCS_J2000= 1,//J2000(FK5) right ascension and declination
			eWCS_B1950= 2,//B1950(FK4) right ascension and declination
			eWCS_GALACTIC= 3,//Galactic longitude and latitude
			eWCS_ECLIPTIC= 4,//Ecliptic longitude and latitude
			eWCS_ALTAZ= 5,//Azimuth and altitude/elevation
			eWCS_LINEAR= 6,//Linear with optional units
			eWCS_NPOLE= 7,//Longitude and north polar angle
			eWCS_SPA= 8,//Longitude and south polar angle
			eWCS_PLANET= 9,//Longitude and latitude on planet
			eWCS_XY= 10,//X-Y Cartesian coordinates 
			eWCS_ICRS= 11//ICRS right ascension and declination
		};

		//Projections (1-26 are WCSLIB) (values for wcs->prjcode)
		enum WCSProjType {
			eWCS_PIX= -1,//Pixel WCS 
			eWCS_LIN= 0,//Linear projection
			eWCS_AZP= 1,//Zenithal/Azimuthal Perspective
			eWCS_SZP= 2,//Zenithal/Azimuthal Perspective
			eWCS_TAN= 3,//Gnomonic = Tangent Plane
			eWCS_SIN= 4,//Orthographic/synthesis
			eWCS_STG= 5,//Stereographic
			eWCS_ARC= 6,//Zenithal/azimuthal equidistant
			eWCS_ZPN= 7,//Zenithal/azimuthal PolyNomial
			eWCS_ZEA= 8,//Zenithal/azimuthal Equal Area
			eWCS_AIR= 9,//Airy
			eWCS_CYP= 10,//CYlindrical Perspective
			eWCS_CAR= 11,//Cartesian
			eWCS_MER= 12,//Mercator
			eWCS_CEA= 13,//Cylindrical Equal Area
			eWCS_COP= 14,//Conic PerSpective (COP)
			eWCS_COD= 15,//COnic equiDistant
			eWCS_COE= 16,//COnic Equal area
			eWCS_COO= 17,//COnic Orthomorphic
			eWCS_BON= 18,//Bonne
			eWCS_PCO= 19,//Polyconic
			eWCS_SFL= 20,//Sanson-Flamsteed (GLobal Sinusoidal)
			eWCS_PAR= 21,//Parabolic
			eWCS_AIT= 22,//Hammer-Aitoff
			eWCS_MOL= 23,//Mollweide
			eWCS_CSC= 24,//COBE quadrilateralized Spherical Cube
			eWCS_QSC= 25,//Quadrilateralized Spherical Cube
			eWCS_TSC= 26,//Tangential Spherical Cube
			eWCS_NCP= 27,//Special case of SIN from AIPS
			eWCS_GLS= 28,//Same as SFL from AIPS
			eWCS_DSS= 29,//Digitized Sky Survey plate solution
			eWCS_PLT= 30,//Plate fit polynomials (SAO)
			eWCS_TNX= 31,//Tangent Plane (NOAO corrections)
			eWCS_ZPX= 32,//Zenithal Azimuthal Polynomial (NOAO corrections)
			eWCS_TPV= 33,//Tangent Plane (SCAMP corrections)
			eNWCSTYPE= 34//Number of WCS types (-1 really means no WCS)
		};

		enum WCSProjCode {
			AZP = 101,
			SZP = 102,
			TAN = 103,
			STG = 104,
			SIN = 105,
			ARC = 106,
			ZPN = 107,
			ZEA = 108,
			AIR = 109,
			CYP = 201,
			CEA = 202,
			CAR = 203,
			MER = 204,
			SFL = 301,	
			PAR = 302,
			MOL = 303,
			AIT = 401,
			COP = 501,
			COE = 502,
			COD = 503,
			COO = 504,
			BON = 601,
			PCO = 602,
			TSC = 701,
			CSC = 702,
			QSC = 703
		};

		//Method to use
		enum WCSMethod {
			eWCS_BEST= 0,// Use best WCS projections
			eWCS_ALT= 1,//Use not best WCS projections
			eWCS_OLD= 2,//Use AIPS WCS projections
			eWCS_NEW= 3//Use WCSLIB 2.5 WCS projections
		};
		

		// Distortion codes (values for wcs->distcode)
		enum DistorsionCode {
			eDISTORT_NONE= 0,//No distortion coefficients
			eDISTORT_SIRTF= 1//SIRTF distortion matrix
		};

		// TNX/ZPX permitted types of surfaces 
		enum TnxSurfaceType {
			eTNX_CHEBYSHEV= 1,
			eTNX_LEGENDRE= 2,
			eTNX_POLYNOMIAL= 3
		};

		// TNX/ZPX cross-terms flags
		enum TnxCrossTermsFlag {
			eTNX_XNONE= 0,//no x-terms (old no)
			eTNX_XFULL= 1,//full x-terms (new yes)
			eTNX_XHALF= 2//half x-terms (new)
		};

		//========================================
		//==   WCS structs
		//========================================
		typedef struct poly
		{
			double        *basis;         /* Current values of the basis functions */
  		double        *coeff;         /* Polynom coefficients */
  		int           ncoeff;         /* Number of coefficients */
  		int           *group;         /* Groups */
  		int           ndim;           /* dimensionality of the polynom */
  		int           *degree;        /* Degree in each group */
  		int           ngroup;         /* Number of different groups */
		} polystruct;

		struct wcsprm {
   		int flag;
   		char pcode[4];
   		char lngtyp[5], lattyp[5];
   		int lng, lat;
   		int cubeface;
		};

		struct linprm {
   		int flag;
   		int naxis;
  	 	double *crpix;
   		double *pc;
   		double *cdelt;

   		// Intermediates
   		double *piximg;
   		double *imgpix;
		};

		struct celprm {
   		int flag;
   		double ref[4];
   		double euler[5];
		};

		struct prjprm 
		{
			//Defines
			static const int _MAXPV= 100;

			//Struct members
			char code[4];
  		int flag;
  		double phi0, theta0;
  		double r0;
  		double p[10];
  		double w[20];
  		int    n;
  		int npv;
  		double ppv[2*_MAXPV];
  		poly* inv_x;
  		poly* inv_y;

			int (*prjfwd)(const double, const double,prjprm*,double*, double*);
  		int (*prjrev)(const double, const double,prjprm*,double*, double*);
		};

		
		// TNX/ZPX surface fitting structure and flags
		struct IRAFsurface {
  		double xrange;	/* 2. / (xmax - xmin), polynomials */
  		double xmaxmin;	/* - (xmax + xmin) / 2., polynomials */
  		double yrange;	/* 2. / (ymax - ymin), polynomials */
  		double ymaxmin;	/* - (ymax + ymin) / 2., polynomials */
  		int	 type;		/* type of curve to be fitted */
  		int    xorder;	/* order of the fit in x */
  		int    yorder;	/* order of the fit in y */
  		int    xterms;	/* cross terms for polynomials */
  		int    ncoeff;	/* total number of coefficients */
  		double *coeff;	/* pointer to coefficient vector */
  		double *xbasis;	/* pointer to basis functions (all x) */
  		double *ybasis;	/* pointer to basis functions (all y) */
		};

		/* SIRTF distortion matrix coefficients */
		//#define DISTMAX 10
		struct Distort 
		{
			static const int _DISTMAX= 10;
  		int    a_order;                /* max power for the 1st dimension */
  		double a[_DISTMAX][_DISTMAX];  /* coefficient array of 1st dimension */
  		int    b_order;                /* max power for 1st dimension */
  		double b[_DISTMAX][_DISTMAX];  /* coefficient array of 2nd dimension */
  		int    ap_order;               /* max power for the 1st dimension */
  		double ap[_DISTMAX][_DISTMAX]; /* coefficient array of 1st dimension */
  		int    bp_order;               /* max power for 1st dimension */
  		double bp[_DISTMAX][_DISTMAX]; /* coefficient array of 2nd dimension */
		};

		

		//========================================
		//==   WCS functions
		//========================================
		/**
		* \brief Convert World Coordinates to pixel coordinates
		*/
		static void wcsc2pix(WCS*wcs,double xpos,double ypos,char *coorsys,double *xpix,double *ypix,int *offscl);
		/**
		* \brief Convert World Coordinates to pixel coordinates
		*/
		static int wcspix(double xpos,double ypos,WCS* wcs,double* xpix,double* ypix);

		/**
		* \brief Convert from RA,Dec to pixel location
		*/
		static int dsspix(double xpos,double ypos,WCS* wcs,double *xpix,double *ypix);
		/**
		* \brief Convert from RA,Dec to pixel location
		*/
		static int platepix(double xpos,double ypos,WCS* wcs,double *xpix,double *ypix);

		/**
		* \brief Convert focal plane to pixel coordinates
		*/
		static void foc2pix(WCS* wcs,double u,double v,double *x,double *y);

		/**
		* \brief Inverse transform (world to physical) gnomonic projection
		*/
		static int tnxpix(double xpos,double ypos,WCS* wcs,double *xpix,double *ypix);

		
		

		/**
		* \brief Set projection in WCS structure from FITS keyword values
		*/
		static int wcstype(WCS* wcs,char* ctype1,char* ctype2);

		/**
		* \brief Set scale and rotation in WCS structure from axis scale and rotation
		*/
		static void wcsdeltset(WCS* wcs,double cdelt1,double cdelt2,double crota);

		/**
		* \brief Set error message
		*/
		static void setwcserr(const char* errmsg){ strcpy (g_wcserrmsg, errmsg); return; }

		/**
		* \brief Set WCS command
		*/
		static void setwcscom(WCS* wcs);
		/**
		* \brief Compute image rotation 
		*/
		static void wcsrotset(WCS* wcs);
		
		/**
		* \brief Initialize catalog search command set by -wcscom
		*/
		static void wcscominit(WCS* wcs,int i,const char* command);

		/**
		* \brief Set coordinate system from string 
		*/
		static int wcscsys(char* wcstring);

		/**
		* \brief Get WCS position
		*/
		static int wcspos (double xpix, double ypix, WCS* wcs,double* xpos,double* ypos);

		/**
		* \brief Convert pixel to focal plane coordinates
		*/
		static void pix2foc(WCS* wcs,double x,double y,double *u,double *v);


		/**
		* \brief Set scale and rotation in WCS structure
		*/
		static void wcscdset(WCS* wcs,double* cd);

		/**
		* \brief Set up rotation matrix for WCSLIB projection subroutines
		*/
		static void wcslibrot(WCS* wcs);

		/**
		* \brief Delete WCS struct
		*/
		static void wcsfree(WCS* wcs);

		/**
		* \brief Return 0 if WCS structure is filled, else 1
		*/
		static int nowcs(WCS* wcs);

		
		/**
		* \brief Return 1 if WCS structure is filled, else 0
		*/
		static int iswcs(WCS* wcs);

		/**
		* \brief Delete wcs command
		*/
		static void freewcscom (WCS* wcs);


		/**
		* \brief Set WCS distortion code from CTYPEi in FITS header
		*/
		static void setdistcode(WCS* wcs,char* ctype);

		/**
		* \brief Free a polynom structure and everything it contains
		*/
		static void poly_end(polystruct *poly);

		/**
		* \brief Evaluate a multidimensional polynom
		*/
		static double poly_func(polystruct *poly, double *pos);

		/**
		* \brief Invert matrix
		*/
		static int matinv(const int, const double [], double []);
		
		/**
		* \brief Routine to determine accurate position for pixel coordinates 
		* returns 0 if successful otherwise 1 = angle too large for projection
		* based on amdpos() from getimage 
		*/
		static int dsspos(double xpix,double ypix,WCS* wcs,double* xpos,double* ypos);
		/**
		* \brief Routine to determine accurate position for pixel coordinates
		* returns 0 if successful otherwise 1 = angle too large for projection;
		* based on amdpos() from getimage
		*/
		static int platepos (double xpix,double ypix,WCS* wcs,double* xpos,double* ypos);
		/**
		* \brief Forward transform (physical to world) gnomonic projection
		*/
		static int tnxpos(double xpix,double ypix,WCS* wcs,double* xpos,double* ypos);

		
		/**
		* \brief Procedure to evaluate the fitted surface at a single point
		*/
		static double wf_gseval(IRAFsurface* sf,double x,double y);
	
		/**
		* \brief Procedure to evaluate all the non-zero chebyshev function coefficients for a given x and order
		*/
		static void wf_gsb1cheb (double x,int order,double k1,double k2,double* basis);

		/**
		* \brief Procedure to evaluate all the non-zero legendre functions for a single point and given order.
		*/
		static void wf_gsb1leg (double x,int order,double k1,double k2,double* basis);
		
		/**
		* \brief Procedure to evaluate all the non-zero polynomial functions for a single point and given order
		*/
		static void wf_gsb1pol (double x,int order,double* basis);

		/**
		* \brief Procedure to calculate a new surface which is a derivative of the input surface
		*/
		static double wf_gsder(IRAFsurface* sf1,double x,double y,int nxd,int nyd);

		/**
		* \brief Procedure to fetch the number and magnitude of the coefficients
		*/
		static int wf_gscoeff (IRAFsurface* sf,double* coeff);

		/**
		* \brief Procedure to free the surface descriptor
		*/
		static void wf_gsclose (IRAFsurface* sf);

		/**
		* \brief Forward transform (physical to world) gnomonic projection
		*/
		static int zpxpos (double xpix,double ypix,WCS* wcs,double* xpos,double* ypos);

		
		/**
		* \brief Inverse transform (world to physical) for the zenithal azimuthal polynomial projection
		*/
		static int zpxpix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix);

		/**
		* \brief Routine to determine accurate position for pixel coordinates
		* returns 0 if successful otherwise 1 = angle too large for projection;
		* does: -SIN, -TAN, -ARC, -NCP, -GLS or -SFL, -MER, -AIT projections
		* anything else is linear
		*/
		static int worldpos(double xpix,double ypix,WCS* wcs,double* xpos,double* ypos);

		/**
		* \brief Routine to determine accurate pixel coordinates for an RA and Dec
		*/
		static int worldpix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix);

		/**
		* \brief Convert from coordinate system sys1 to coordinate system sys2
		*/
		static void wcscon(int sys1,int sys2,double eq1,double eq2,double* dtheta,double* dphi,double epoch);
		/**
		* \brief Precession -  FK5 (Fricke, post-IAU1976)
		* This routine will not correctly convert between FK5 and FK4
		* For output in FK4, precess to 2000.0 and use fk524() on result
		* Based on slaPreces(), P.T.Wallace   Starlink   22 December 1993
		*/
		static void fk5prec(double ep0,double ep1,double* ra,double* dec);

		/**
		* \brief The following routines are modified from Patrick Wallace's SLALIB
		* Precess coordinates between epochs in FK4
		* Precession - FK4 (Bessel-Newcomb, pre-IAU1976)
		* This routine will not correctly convert between FK4 and FK5
		* For output in FK5, precess to 1950.0 and use fk425() on result
		* Based on slaPreces(), P.T.Wallace   Starlink   22 December 1993
		*/
		static void fk4prec(double ep0,double ep1,double* ra,double* dec);

		/**
		* \brief Conversion function 
		*/
		static void fk524 (double* ra,double* dec);

		/**
		* \brief Conversion function 
		*/
		static void fk524e(double* ra,double* dec,double epoch);
		/**
		* \brief Conversion function 
		*/
		static void fk524m(double* ra,double* dec,double* rapm,double* decpm);
		/**
		* \brief This routine converts stars from the IAU 1976 FK5 Fricke
    * system, to the old Bessel-Newcomb FK4 system, using Yallop's
    * implementation
		*/
		static void fk524pv(double* ra,double* dec,double* rapm,double* decpm,double* parallax,double* rv);

		/**
		* \brief Transform IAU 1958 galactic coordinates to B1950.0 'FK4' equatorial coordinates
		*/
		static void gal2fk4 (double* dtheta,double* dphi);

		
		/**
		* \brief Convert B1950 right ascension and declination to ecliptic coordinates
		*/
		static void fk42ecl (double *dtheta,double *dphi,double epoch);
		/**
		* \brief Conversion function 
		*/
		static void fk425m (double* ra,double* dec,double* rapm,double* decpm);
	
		/**
		* \brief Convert J2000 right ascension and declination to ecliptic coordinates
		*/
		static void fk52ecl (double *dtheta,double *dphi,double epoch);


		/**
		* \brief Transform J2000 equatorial coordinates to IAU 1958 galactic coordinates
		*/
		static void fk52gal (double *dtheta,double *dphi);

		/**
		* \brief Transform B1950.0 FK4 equatorial coordinates to IAU 1958 galactic coordinates
		*/
		static void fk42gal (double *dtheta,double *dphi);

		/**
		* \brief Transform IAU 1958 galactic coordinates to J2000 equatorial coordinates
		*/
		static void gal2fk5 (double *dtheta,double *dphi);
	
		/**
		* \brief Convert ecliptic coordinates to B1950 right ascension and declination
		*/
		static void ecl2fk4 (double *dtheta,double *dphi,double epoch);

		/**
		* \brief Convert ecliptic coordinates to J2000 right ascension and declination
		*/
		static void ecl2fk5 (double *dtheta,double *dphi,double epoch);

		/**
		* \brief Conversion function 
		*/
		static void fk425 (double	*ra,double	*dec);

		/**
		* \brief Conversion function 
		*/
		static void fk425e (double* ra,double* dec,double epoch);

		/**
		* \brief This routine converts stars from the old Bessel-Newcomb FK4 system to the IAU 1976 FK5 Fricke system, using Yallop's implementation
		*/
		static void fk425pv (double* ra,double* dec,double* rapm,double* decpm,double* parallax,double* rv);

		/**
		* \brief Form the matrix of precession between two epochs (IAU 1976, FK5)
		* Notes:
		*	  1) The epochs are TDB (loosely ET) Julian epochs.
		*   2) The matrix is in the sense   v(ep1)  =  rmatp * v(ep0)
		* References:
		*   - Lieske,J.H., 1979. Astron. Astrophys.,73,282 (equations (6) & (7), p283)
		*   - Kaplan,G.H., 1981. USNO circular no. 163, pa2
		* Based on slaPrec(), P.T.Wallace   Starlink   31 October 1993
		*/
		static void mprecfk5(double ep0,double ep1,double rmatp[9]);

		/**
		* \brief Generate the matrix of precession between two epochs,
		*  using the old, pre-IAU1976, Bessel-Newcomb model, using
		*  Kinoshita's formulation (double precision)
		*
		*  The matrix is in the sense   v(bep1)  =  rmatp * v(bep0)	
		*
		*  Reference:
		*     Kinoshita, H. (1975) 'Formulas for precession', SAO Special
		*     Report No. 364, Smithsonian Institution Astrophysical
		*     Observatory, Cambridge, Massachusetts.
		*
		*  Based on slaPrebn() by P.T.Wallace   Starlink   30 October 1993
		*/
		static void mprecfk4(double bep0,double bep1,double rmatp[9]);

		/**
		* \brief Make 3-D rotation matrix from up to three rotations
		*/
		static void rotmat(int axes,double rot1,double rot2,double rot3,double* matrix);

		/**
		* \brief Set equinox from string (return 0.0 if not obvious) 
		*/
		static double wcsceq(char* wcstring);

		/**
		* \brief Return string with right ascension in hours and declination in degrees
		*/
		static char* eqstrn (double dra,double ddec);

		
		/**
		* \brief Convert right ascension, declination, and distance to geocentric equatorial rectangular coordinates 
		*/
		static void s2v3 (double rra,double rdec,double r,double pos[3]);

		
		/**
		* \brief Convert geocentric equatorial rectangular coordinates to right ascension, declination, and distance
		*/
		static void v2s3 (double pos[3],double *rra,double *rdec,double *r);

		static int wcsrev(const char[][16],
              wcsprm*,
              const double[], 
              linprm*,
              double[], 
             	prjprm*, 
              double *,
              double *, 
              const double[], 
              celprm*, 
              double[]);

		static int celrev(const char *,
              const double, const double,
              prjprm *,
              double *, double *,
              celprm *,
              double *, double *);
		
		static int celset(const char *, celprm *, prjprm *);
		static int prjset(const char [], prjprm *);

		static int azpset(prjprm* prj);
		static int szpset(prjprm* prj);
		static int tanset(prjprm* prj);
		static int stgset(prjprm * prj);
		static int sinset(prjprm * prj);
		static int arcset(prjprm *prj);
		static int zpnset(prjprm* prj);
		static int zeaset(prjprm* prj);
		static int airset(prjprm * prj);
		static int cypset(prjprm * prj);
		static int ceaset(prjprm * prj);
		static int carset(prjprm * prj);
		static int merset(prjprm * prj);
		static int sflset(prjprm * prj);
		static int parset(prjprm * prj);
		static int molset(prjprm * prj);
		static int aitset(prjprm *prj);
		static int copset(prjprm * prj);
		static int coeset(prjprm * prj);
		static int codset(prjprm * prj);
		static int cooset(prjprm * prj);
		static int bonset(prjprm * prj);
		static int pcoset(prjprm * prj);
		static int tscset(prjprm * prj);
		static int cscset(prjprm* prj);
		static int qscset(prjprm * prj);
	
		static int parfwd(const double, const double, prjprm *, double *, double *);
		static int parrev(const double, const double, prjprm *, double *, double *);
		static int molfwd(const double, const double, prjprm *, double *, double *);
		static int molrev(const double, const double, prjprm *, double *, double *);
		static int aitfwd(const double, const double, prjprm *, double *, double *);
  	static int aitrev(const double, const double, prjprm *, double *, double *);
		static int copfwd(const double, const double, prjprm *, double *, double *);
   	static int coprev(const double, const double, prjprm *, double *, double *);
		static int coefwd(const double, const double, prjprm *, double *, double *);
   	static int coerev(const double, const double, prjprm *, double *, double *);
		static int codfwd(const double, const double, prjprm *, double *, double *);
   	static int codrev(const double, const double, prjprm *, double *, double *);
		static int coofwd(const double, const double, prjprm *, double *, double *);
   	static int coorev(const double, const double, prjprm *, double *, double *);
		static int bonfwd(const double, const double, prjprm *, double *, double *);
   	static int bonrev(const double, const double, prjprm *, double *, double *);
		static int sflfwd(const double, const double, prjprm *, double *, double *);
   	static int sflrev(const double, const double, prjprm *, double *, double *);
		static int pcofwd(const double, const double, prjprm *, double *, double *);
   	static int pcorev(const double, const double, prjprm *, double *, double *);
		static int tscfwd(const double, const double, prjprm *, double *, double *);
   	static int tscrev(const double, const double, prjprm *, double *, double *);
		static int cscfwd(const double, const double, prjprm *, double *, double *);
   	static int cscrev(const double, const double, prjprm *, double *, double *);
		static int qscfwd(const double, const double, prjprm *, double *, double *);
   	static int qscrev(const double, const double, prjprm *, double *, double *);
		static int azpfwd(const double, const double, prjprm *, double *, double *);
   	static int azprev(const double, const double, prjprm *, double *, double *);
		static int szpfwd(const double, const double, prjprm *, double *, double *);
   	static int szprev(const double, const double, prjprm *, double *, double *);
		static int tanfwd(const double, const double, prjprm *, double *, double *);
   	static int tanrev(const double, const double, prjprm *, double *, double *);
		static int stgfwd(const double, const double, prjprm *, double *, double *);
   	static int stgrev(const double, const double, prjprm *, double *, double *);
		static int raw_to_pv(prjprm *prj, double x, double y, double *xo, double *yo);
		static int sinfwd(const double, const double, prjprm *, double *, double *);
   	static int sinrev(const double, const double, prjprm *, double *, double *);
		static int arcfwd(const double, const double, prjprm *, double *, double *);
   	static int arcrev(const double, const double, prjprm *, double *, double *);
		static int zeafwd(const double, const double, prjprm *, double *, double *);
   	static int zearev(const double, const double, prjprm *, double *, double *);
		static int airfwd(const double, const double, prjprm *, double *, double *);
   	static int airrev(const double, const double, prjprm *, double *, double *);
		static int ceafwd(const double, const double, prjprm *, double *, double *);
   	static int cearev(const double, const double, prjprm *, double *, double *);
		static int carfwd(const double, const double, prjprm *, double *, double *);
   	static int carrev(const double, const double, prjprm *, double *, double *);
		static int merfwd(const double, const double, prjprm *, double *, double *);
   	static int merrev(const double, const double, prjprm *, double *, double *);
		static int cypfwd(const double, const double, prjprm *, double *, double *);
   	static int cyprev(const double, const double, prjprm *, double *, double *);
		static int zpnfwd(const double, const double, prjprm *, double *, double *);
   	static int zpnrev(const double, const double, prjprm *, double *, double *);

		static int wcsset(const int,const char[][16],wcsprm *);
		static int wcsfwd(const char[][16],wcsprm *,const double[],const double[],celprm *,double *,double *, prjprm *, double[], linprm *,double[]);
		static int celfwd(const char *,
              const double, const double,
              celprm *,
              double *, double *,
              prjprm *,
              double *, double *);
		static int linfwd(const double[], linprm *, double[]);
		static int sphfwd(const double, const double,
              const double [],
              double *, double *);

		static int sphrev(const double, const double,const double [],double *, double *);


		static int prjfwd(const double, const double, prjprm*, double*, double*);
		static int prjrev(const double, const double, prjprm*, double*, double*);
		static int linrev(const double[], linprm *, double[]);
		static int linset(linprm*);

	private:
		/**
		* \brief Find string s2 within null-terminated string s1
		*/
		static char* strsrch (const char* s1,const char* s2);
		/**
		* \brief Find string s2 within string s1
		*/
		static char* strnsrch (const char* s1,const char* s2, const int ls1);

		/**
		* \brief Return 1 if string is an integer number,
		* 2 if floating point,
		* 3 if sexigesimal, with or without decimal point
		* 4 if yyyy-mm-dd date
		* else 0
		*/
		static int isnum (const char* string);

		/**
		* \brief Convert cos to deg
		*/
		static double cosdeg(const double angle);
		/**
		* \brief Convert sin to deg
		*/
		static double sindeg(const double angle);
		/**
		* \brief Compute arctan in deg
		*/
		static double atan2deg(const double y,const double x);
		/**
		* \brief Compute tan in deg
		*/
		static double tandeg(const double angle);
		/**
		* \brief Compute arctan in deg
		*/
		static double atandeg (const double v);
		/**
		* \brief Compute arccos in deg
		*/
		static double acosdeg(const double v);
		/**
		* \brief Compute arcsin in deg
		*/
		static double asindeg (const double v);

		/**
		* \brief Format angle into decimal degrees string
		*/
		static void deg2str(char* str,int lstr,const double deg,const int ndec);
		/**
		* \brief Convert degrees to hh:mm:ss.ss
		*/
		static void ra2str(char* str,int lstr,const double ra,const int ndec);
		/**
		* \brief Convert degrees to dd:mm:ss.ss
		*/
		static void dec2str(char* str,int lstr,const double dec,const int ndec);
		/**
		* \brief Format number into string
		*/
		static void num2str(char *str,const double  num,const int field,const int ndec);

		/**
		* \brief Convert deg to rad
		*/
		template<typename T>
		static T degrad_i(T x) 
		{
			return (x)*kPI/180.;
		}

		/**
		* \brief Convert rad to deg
		*/
		template<typename T>
		static T raddeg_i(T x) 
		{
			return (x)*180./kPI;
		}

		/**
		* \brief Convert seconds to rad
		*/
		template<typename T>
		static T secrad_i(T x) 
		{
			return (x)*kAS2R;
		}

		template<typename T,typename K>
		static T copysgn_i(T X,K Y){
			return ((Y) < 0.0 ? -fabs(X) : fabs(X));
		}

		template<typename T,typename K>
		static T copysgni_i(T X,K Y){
			return ((Y) < 0 ? -abs(X) : abs(X));
		}



		
	protected:

		static double g_zpix;
		static int g_izpix;
		static char* g_wcscom0[10];
		static int g_wcsproj0;
		static char g_wcserrmsg[80];

		static const int kWCSSET= 137;
		static const int kLINSET= 137;
		static const int kCELSET= 137;
		static constexpr double kPI= 3.1415926535898;
		static constexpr double d2pi = 6.283185307179586476925287;//two PI
		static constexpr double kAS2R= 4.8481368110953e-6;//arcsec to rad
		//static constexpr double d2r= kPI / 180.0;
		//static constexpr double r2d= 180.0 / kPI;
		static constexpr double kD2R= kPI/180.0;
		static constexpr double kR2D= 180.0/kPI;
		static constexpr double kSQRT2= 1.4142135623730950488;
		static constexpr double kSQRT2INV= 1.0/kSQRT2;
		static constexpr double kWCSTRIG_TOL= 1e-10;//Domain tolerance for asin and acos functions
		static constexpr double kBADCVAL= 0.0;
		static constexpr double kSPHTOL= 0.00001;
		//static constexpr double TOL= 1.0e-11;
		static const int kDISTMAX= 10;
		static const int kMAXPV= 100;
		static const int kPOLY_MAXDIM= 4;//Max dimensionality of polynom
		static const int kPOLY_MAXDEGREE= 10;//Max degree of the polynom
		//static const int max_niter= 500;
		static const int kMAX_NITER= 500;
		

		static constexpr double g_a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6};
		static constexpr double g_ad[3] = {1.245e-3,  -1.580e-3,  -0.659e-3};
		static constexpr double g_tiny= 1.e-30;//small number to avoid arithmetic problems

		static const int g_npcode = 26;
		static constexpr char g_pcodes[26][4] =
      {"AZP", "SZP", "TAN", "STG", "SIN", "ARC", "ZPN", "ZEA", "AIR", "CYP",
       "CEA", "CAR", "MER", "COP", "COE", "COD", "COO", "SFL", "PAR", "MOL",
       "AIT", "BON", "PCO", "TSC", "CSC", "QSC"};

		
		static constexpr double g_jgal[3][3] = {
			{-0.054875539726,-0.873437108010,-0.483834985808},
			{0.494109453312,-0.444829589425, 0.746982251810},
			{-0.867666135858,-0.198076386122, 0.455983795705}
		};

		static constexpr double g_bgal[3][3] = {
			{-0.066988739415,-0.872755765852,-0.483538914632},
			{0.492728466075,-0.450346958020, 0.744584633283},
			{-0.867600811151,-0.188374601723, 0.460199784784}
		};

		static const int g_idg= 0;


		/* 
		Convert B1950.0 FK4 star data to J2000.0 FK5 
		*/
		static constexpr double g_em[6][6] = {
    	{	0.9999256782,		/* em[0][0] */
				-0.0111820611,		/* em[0][1] */
				-0.0048579477,		/* em[0][2] */
	 			0.00000242395018,	/* em[0][3] */
				-0.00000002710663,	/* em[0][4] */
				-0.00000001177656 
			},	/* em[0][5] */
    	{	0.0111820610,		/* em[1][0] */
	 			0.9999374784,		/* em[1][1] */
				-0.0000271765,		/* em[1][2] */
		 		0.00000002710663,	/* em[1][3] */
	 			0.00000242397878,	/* em[1][4] */
				-0.00000000006587 
			},	/* em[1][5] */
    	{	0.0048579479,		/* em[2][0] */
				-0.0000271474,		/* em[2][1] */
	 			0.9999881997,		/* em[2][2] */
	 			0.00000001177656,	/* em[2][3] */
				-0.00000000006582,	/* em[2][4] */
	 			0.00000242410173 
			},	/* em[2][5] */
    	{ -0.000551,		/* em[3][0] */
				-0.238565,		/* em[3][1] */
	 			0.435739,		/* em[3][2] */
	 			0.99994704,		/* em[3][3] */
				-0.01118251,		/* em[3][4] */
				-0.00485767 
			},		/* em[3][5] */
    	{	0.238514,		/* em[4][0] */
				-0.002667,		/* em[4][1] */
				-0.008541,		/* em[4][2] */
	 			0.01118251,		/* em[4][3] */
	 			0.99995883,		/* em[4][4] */
				-0.00002718 
			},		/* em[4][5] */
    	{	-0.435623,		/* em[5][0] */
	 			0.012254,		/* em[5][1] */
	 			0.002117,		/* em[5][2] */
	 			0.00485767,		/* em[5][3] */
				-0.00002714,		/* em[5][4] */
	 			1.00000956 
			}		/* em[5][5] */
    };


		/* 
		FK524  convert J2000 FK5 star data to B1950 FK4
   	based on Starlink sla_fk524 by P.T.Wallace 27 October 1987 
		*/
		static constexpr double g_emi[6][6] = {
    	{	0.9999256795,		/* emi[0][0] */
	 			0.0111814828,		/* emi[0][1] */
	 			0.0048590039,		/* emi[0][2] */
				-0.00000242389840,	/* emi[0][3] */
				-0.00000002710544,	/* emi[0][4] */
				-0.00000001177742 
			},	/* emi[0][5] */ 
    	{	-0.0111814828,		/* emi[1][0] */
	 			0.9999374849,		/* emi[1][1] */
				-0.0000271771,		/* emi[1][2] */
	 			0.00000002710544,	/* emi[1][3] */
				-0.00000242392702,	/* emi[1][4] */
	 			0.00000000006585 		
			},	/* emi[1][5] */
    	{	-0.0048590040,		/* emi[2][0] */
				-0.0000271557,		/* emi[2][1] */
	 			0.9999881946,		/* emi[2][2] */
	 			0.00000001177742,	/* emi[2][3] */
	 			0.00000000006585,	/* emi[2][4] */
				-0.00000242404995 
			},	/* emi[2][5] */ 
    	{	-0.000551,		/* emi[3][0] */
	 			0.238509,		/* emi[3][1] */
				-0.435614,		/* emi[3][2] */
	 			0.99990432,		/* emi[3][3] */
	 			0.01118145,		/* emi[3][4] */
	 			0.00485852 
			},		/* emi[3][5] */ 
   	 	{	-0.238560,		/* emi[4][0] */
				-0.002667,		/* emi[4][1] */
	 			0.012254,		/* emi[4][2] */
				-0.01118145,		/* emi[4][3] */
	 			0.99991613,		/* emi[4][4] */
				-0.00002717 
			},		/* emi[4][5] */ 
    	{	 0.435730,		/* emi[5][0] */
				-0.008541,		/* emi[5][1] */
	 			0.002117,		/* emi[5][2] */
				-0.00485852,		/* emi[5][3] */
				-0.00002716,		/* emi[5][4] */
	 			0.99996684 
			}		/* emi[5][5] */
    };


	friend struct WCS;

	ClassDef(WCSUtils,1)

};//close class


//========================================
//==   WorldCoor
//========================================
class WCS : public TObject
{
	public:
		WCS(){};
		virtual ~WCS()
		{
			
		};//close destructor


	private:
		//Defines
		static const int kMAXPV= 100;
	
	public:
	
		//Struct members
  	double	xref;		/* X reference coordinate value (deg) */
  	double	yref;		/* Y reference coordinate value (deg) */
  	double	xrefpix;	/* X reference pixel */
  	double	yrefpix;	/* Y reference pixel */
  	double	xinc;		/* X coordinate increment (deg) */
  	double	yinc;		/* Y coordinate increment (deg) */
  	double	rot;		/* rotation around axis (deg) (N through E) */
  	double	cd[4];		/* rotation matrix */
  	double	dc[4];		/* inverse rotation matrix */
  	double	equinox;	/* Equinox of coordinates default to 1950.0 */
  	double	epoch;		/* Epoch of coordinates default to equinox */
  	double	nxpix;		/* Number of pixels in X-dimension of image */
  	double	nypix;		/* Number of pixels in Y-dimension of image */
  	double	plate_ra;	/* Right ascension of plate center */
  	double	plate_dec;	/* Declination of plate center */
  	double	plate_scale;	/* Plate scale in arcsec/mm */
  	double	x_pixel_offset;	/* X pixel offset of image lower right */
  	double	y_pixel_offset;	/* Y pixel offset of image lower right */
  	double	x_pixel_size;	/* X pixel_size */
  	double	y_pixel_size;	/* Y pixel_size */
  	double	ppo_coeff[6];	/* pixel to plate coefficients for DSS */
 	 	double	x_coeff[20];	/* X coefficients for plate model */
  	double	y_coeff[20];	/* Y coefficients for plate model */
  	double	xpix;		/* X (RA) coordinate (pixels) */
  	double	ypix;		/* Y (dec) coordinate (pixels) */
  	double	zpix;		/* Z (face) coordinate (pixels) */
  	double	xpos;		/* X (RA) coordinate (deg) */
  	double	ypos;		/* Y (dec) coordinate (deg) */
  	double	crpix[9];	/* Values of CRPIXn keywords */
  	double	crval[9];	/* Values of CRVALn keywords */
  	double	cdelt[9];	/* Values of CDELTn keywords */
  	double	pc[81];		/* Values of PCiiijjj keywords */
  	double	projp[10];	/* Constants for various projections */
  	int		pvfail;		/* If non-zero, significant inaccuracy likely to occur in projection */
  	double	projppv[2*kMAXPV]; /* SCAMP constants for the PV coordinates */
  	WCSUtils::poly* inv_x;		/* SCAMP projection correction polynom in x */
  	WCSUtils::poly* inv_y;		/* SCAMP projection correction polynom in y */
  	double	longpole;	/* Longitude of North Pole in degrees */
  	double	latpole;	/* Latitude of North Pole in degrees */
  	double	rodeg;		/* Radius of the projection generating sphere */
  	double	imrot;		/* Rotation angle of north pole */
  	double	pa_north;	/* Position angle of north (0=horizontal) */
  	double	pa_east;	/* Position angle of east (0=horizontal) */
  	double	radvel;		/* Radial velocity (km/sec away from observer)*/
  	double	zvel;		/* Radial velocity (v/c away from observer)*/
  	double	zpzd;		/* Colat of FIP (degs) */
  	double	zpr;		/* Radius of FIP (degs) */
  	int		imflip;		/* If not 0, image is reflected around axis */
  	int		prjcode;	/* projection code (-1-32) */
  	int		latbase;	/* Latitude base 90 (NPA), 0 (LAT), -90 (SPA) */
  	int		ncoeff1;	/* Number of x-axis plate fit coefficients */
  	int		ncoeff2;	/* Number of y-axis plate fit coefficients */
 	 	int		zpnp;		 /* ZP polynomial order (0-9) */
  	int		changesys;	/* 1 for FK4->FK5, 2 for FK5->FK4 */
  				/* 3 for FK4->galactic, 4 for FK5->galactic */
  	int		printsys;	/* 1 to print coordinate system, else 0 */
  	int		ndec;		/* Number of decimal places in PIX2WCST */
  	int		degout;		/* 1 to always print degrees in PIX2WCST */
  	int		tabsys;		/* 1 to put tab between RA & Dec, else 0 */
  	int		rotmat;		/* 0 if CDELT, CROTA; 1 if CD */
  	int		coorflip;	/* 0 if x=RA, y=Dec; 1 if x=Dec, y=RA */
  	int		offscl;		/* 0 if OK, 1 if offscale */
  	int		wcson;		/* 1 if WCS is set, else 0 */
  	int		naxis;		/* Number of axes in image (for WCSLIB 3.0) */
  	int		naxes;		/* Number of axes in image */
  	int		wcsproj;	/* eWCS_OLD: AIPS worldpos() and worldpix()
				   eWCS_NEW: Mark Calabretta's WCSLIB subroutines
				   eWCS_BEST: WCSLIB for all but CAR,COE,NCP
				   eWCS_ALT:  AIPS for all but CAR,COE,NCP */
  	int		linmode;	/* 0=system only, 1=units, 2=system+units */
  	int		detector;	/* Instrument detector number */
  	char		instrument[32];	/* Instrument name */
  	char		ctype[9][16];	/* Values of CTYPEn keywords */
  	char		c1type[8];	/*  1st coordinate type code:
					RA--, GLON, ELON */
  	char		c2type[8];	/*  2nd coordinate type code:
					DEC-, GLAT, ELAT */
  	char		ptype[8];	/*  projection type code:
				    SIN, TAN, ARC, NCP, GLS, MER, AIT, etc */
  	char		units[9][32];	/* Units if LINEAR */
  	char		radecsys[32];	/* Reference frame: FK4, FK4-NO-E, FK5, GAPPT*/
  	char		radecout[32];	/* Output reference frame: FK4,FK5,GAL,ECL */
  	char		radecin[32];	/* Input reference frame: FK4,FK5,GAL,ECL */
  	double	eqin;		/* Input equinox (match sysin if 0.0) */
  	double	eqout;		/* Output equinox (match sysout if 0.0) */
  	int		sysin;		/* Input coordinate system code */
  	int		syswcs;		/* WCS coordinate system code */
  	int		sysout;		/* Output coordinate system code */
				/* eWCS_B1950, eWCS_J2000, eWCS_ICRS, eWCS_GALACTIC,
				 * eWCS_ECLIPTIC, eWCS_LINEAR, eWCS_ALTAZ  */
  	char		center[32];	/* Center coordinates (with frame) */
  	//struct wcsprm wcsl;		/* WCSLIB main projection parameters */
  	//struct linprm lin;		/* WCSLIB image/pixel conversion parameters */
  	//struct celprm cel;		/* WCSLIB projection type */
  	//struct prjprm prj;		/* WCSLIB projection parameters */
  	//struct IRAFsurface *lngcor;	/* RA/longitude correction structure */
  	//struct IRAFsurface *latcor;	/* Dec/latitude correction structure */
		WCSUtils::wcsprm wcsl;		/* WCSLIB main projection parameters */
  	WCSUtils::linprm lin;		/* WCSLIB image/pixel conversion parameters */
  	WCSUtils::celprm cel;		/* WCSLIB projection type */
  	WCSUtils::prjprm prj;		/* WCSLIB projection parameters */
  	WCSUtils::IRAFsurface *lngcor;	/* RA/longitude correction structure */
  	WCSUtils::IRAFsurface *latcor;	/* Dec/latitude correction structure */
  	int		distcode;	/* Distortion code 0=none 1=SIRTF */
  	//struct WCSUtils::Distort distort;	/* SIRTF distortion coefficients */
		WCSUtils::Distort distort;	/* SIRTF distortion coefficients */
  	char *command_format[10];	/* WCS command formats */
				/* where %s is replaced by WCS coordinates */
				/* where %f is replaced by the image filename */
				/* where %x is replaced by image coordinates */
  	double	ltm[4];		/* Image rotation matrix */
  	double	ltv[2];		/* Image offset */
  	int		idpix[2];	/* First pixel to use in image (x, y) */
  	int		ndpix[2];	/* Number of pixels to use in image (x, y) */
  	//struct WCS *wcs;	/* WCS upon which this WCS depends */
  	//struct WCS *wcsdep;	/* WCS depending on this WCS */
		WCS* wcs;	/* WCS upon which this WCS depends */
  	WCS* wcsdep;	/* WCS depending on this WCS */
  	char* wcsname;	/* WCS name (defaults to NULL pointer) */
  	char	wcschar;	/* WCS character (A-Z, null, space) */
  	int		logwcs;		/* 1 if DC-FLAG is set for log wavelength */

	ClassDef(WCS,1)

};//close class WCS

#ifdef __MAKECINT__
#pragma link C++ class WCSUtils+;
#pragma link C++ class WCS+;
#endif



}//close namespace

#endif
