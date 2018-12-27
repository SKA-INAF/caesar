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
* @file WCSUtils.cc
* @class WCSUtils
* @brief Utility functions for world coordinate system tasks
*
* Utility functions for world coordinate system tasks
* @author S. Riggi
* @date 15/01/2016
*/


#include <WCSUtils.h>
#include <Consts.h>
#include <Logger.h>
#include <ImgMetaData.h>

#include <TObject.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <ctime>

using namespace std;

ClassImp(Caesar::WCSUtils)
ClassImp(Caesar::WCS)

namespace Caesar {


double WCSUtils::g_zpix= 0.0;
int WCSUtils::g_izpix= 0;
int WCSUtils::g_wcsproj0= 0;
char WCSUtils::g_wcserrmsg[80];
char* WCSUtils::g_wcscom0[10];

constexpr double WCSUtils::g_em[6][6];
constexpr double WCSUtils::g_emi[6][6];
constexpr char WCSUtils::g_pcodes[26][4];

WCSUtils::WCSUtils()
{

}

WCSUtils::~WCSUtils()
{

}


WCS* WCSUtils::ComputeWCSFromImgMetaData(ImgMetaData* metadata,int coordSystem)
{
	//Check metadata
	if(!metadata){
		ERROR_LOG("Null ptr to input image metadata given!");
		return nullptr;
	}

	//## Compute the wcs from vars
	//## NB: FITS keywords CRPIX are defined from [1,N] not 0-based
	//##     Since we use 0-based convention we set Cx->Cx-1 Cy->Cy-1
	WCS* wcs= wcskinit(
		metadata->Nx, metadata->Ny,
		(char*)metadata->CoordTypeX.c_str(),(char*)metadata->CoordTypeY.c_str(),
		metadata->Cx-1, metadata->Cy-1,
		metadata->Xc, metadata->Yc,
		NULL,
		metadata->dX,metadata->dY,
		metadata->RotY,(int)(metadata->Epoch),metadata->Epoch
	);
	//std::string wcsType= std::string(getwcsout(wcs));
	std::string wcsType= GetWCSTypeStr(wcs);
	std::string wcsType_default= metadata->GetWCSType();
	
	//Convert wcs to desired type
	char* flag = (char*)("");
	if(coordSystem==eGALACTIC)
		flag = (char*)("GALACTIC");	
	else if(coordSystem==eJ2000)
		flag = (char*)("FK5");
	else if(coordSystem==eB1950)
		flag = (char*)("FK4");
	else if(coordSystem==-1 && wcsType_default!="")					
		flag = (char*)(wcsType_default.c_str());
	else{
		ERROR_LOG("Invalid coord system type ("<<coordSystem<<") specified, will not build WCS!");
		return nullptr;
	}

	if(strcmp(flag,"")!=0) {
		wcsoutinit (wcs,flag);
	}
			
	return wcs;

}//close ComputeWCSFromImgMetaData()


int WCSUtils::PixelToWCSCoords(double& xpos, double& ypos,WCS* wcs,double ix,double iy) 
{
	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}

	//Convert coords
	pix2wcs (wcs,ix,iy,&xpos, &ypos);
	
	return 0;

}//close PixelToWCSCoords()


int WCSUtils::PixelToWCSStrCoords(std::string& wcs_str,WCS* wcs,double ix,double iy,int max_str_length) 
{
	//Init str
	wcs_str= "";

	//Check WCS
	if(!wcs){
		ERROR_LOG("Null ptr to given WCS!");
		return -1;
	}
	
	if(max_str_length<=0){
		ERROR_LOG("Invalid max wcs string size given (must be >0)!");
		return -1;
	}

	//Convert coords
	char data[max_str_length];
	int status= pix2wcst (wcs,ix,iy,data,max_str_length);
	if(status==0){
		WARN_LOG("Failed to convert pixel coords ("<<ix<<","<<iy<<") to WCS string!");
		return -1;
	}	
	wcs_str= std::string(data);

	return 0;

}//close PixelToWCSStrCoords()


char* WCSUtils::getwcsout(WCS* wcs)
{
  if (nowcs (wcs))
		return (NULL);
  else
		return(wcs->radecout);

}//close getwcsout()


std::string WCSUtils::GetWCSTypeStr(WCS* wcs)
{
	char* wcsType= getwcsout(wcs);
	std::string wcsTypeStr= "";	
	if(wcsType) wcsTypeStr= std::string(wcsType);

	return wcsTypeStr;

}//close GetWCSTypeStr()


/* Convert pixel coordinates to World Coordinate string */

int WCSUtils::pix2wcst(WCS* wcs,double xpix,double ypix,char* wcstring,int lstr)
{
	double	xpos,ypos;
	char	rastr[32], decstr[32];
	int	minlength, lunits, lstring;

	if (nowcs (wcs)) {
	    if (lstr > 0)
		wcstring[0] = 0;
	    return(0);
	    }

	pix2wcs (wcs,xpix,ypix,&xpos,&ypos);

	/* If point is off scale, set string accordingly */
	if (wcs->offscl) {
	    (void)sprintf (wcstring,"Off map");
	    return (1);
	    }

	/* Print coordinates in degrees */
	else if (wcs->degout == 1) {
	    minlength = 9 + (2 * wcs->ndec);
	    if (lstr > minlength) {
		deg2str (rastr, 32, xpos, wcs->ndec);
		deg2str (decstr, 32, ypos, wcs->ndec);
		if (wcs->tabsys)
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
		lstr = lstr - minlength;
		}
	    else {
		if (wcs->tabsys)
		    strncpy (wcstring,"*********	**********",lstr);
		else
		    strncpy (wcstring,"*******************",lstr);
		lstr = 0;
		}
	    }

	/* print coordinates in sexagesimal notation */
	else if (wcs->degout == 0) {
	    minlength = 18 + (2 * wcs->ndec);
	    if (lstr > minlength) {
		if (wcs->sysout == eWCS_J2000 || wcs->sysout == eWCS_B1950) {
		    ra2str (rastr, 32, xpos, wcs->ndec);
		    dec2str (decstr, 32, ypos, wcs->ndec-1);
		    }
		else {
		    dec2str (rastr, 32, xpos, wcs->ndec);
		    dec2str (decstr, 32, ypos, wcs->ndec);
		    }
		if (wcs->tabsys) {
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		    }
		else {
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
		    }
	        lstr = lstr - minlength;
		}
	    else {
		if (wcs->tabsys) {
		    strncpy (wcstring,"*************	*************",lstr);
		    }
		else {
		    strncpy (wcstring,"**************************",lstr);
		    }
		lstr = 0;
		}
	    }

	/* Label galactic coordinates */
	if (wcs->sysout == eWCS_GALACTIC) {
	    if (lstr > 9 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	galactic");
		else
		    strcat (wcstring," galactic");
		}
	    }

	/* Label ecliptic coordinates */
	else if (wcs->sysout == eWCS_ECLIPTIC) {
	    if (lstr > 9 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	ecliptic");
		else
		    strcat (wcstring," ecliptic");
		}
	    }

	/* Label planet coordinates */
	else if (wcs->sysout == eWCS_PLANET) {
	    if (lstr > 9 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	planet");
		else
		    strcat (wcstring," planet");
		}
	    }

	/* Label alt-az coordinates */
	else if (wcs->sysout == eWCS_ALTAZ) {
	    if (lstr > 7 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	alt-az");
		else
		    strcat (wcstring," alt-az");
		}
	    }

	/* Label north pole angle coordinates */
	else if (wcs->sysout == eWCS_NPOLE) {
	    if (lstr > 7 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	long-npa");
		else
		    strcat (wcstring," long-npa");
		}
	    }

	/* Label south pole angle coordinates */
	else if (wcs->sysout == eWCS_SPA) {
	    if (lstr > 7 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	long-spa");
		else
		    strcat (wcstring," long-spa");
		}
	    }

	/* Label equatorial coordinates */
	else if (wcs->sysout==eWCS_B1950 || wcs->sysout==eWCS_J2000) {
	    if (lstr > (int) strlen(wcs->radecout)+1 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	");
		else
		    strcat (wcstring," ");
		strcat (wcstring, wcs->radecout);
		}
	    }

	/* Output linear coordinates */
	else {
	    num2str (rastr, xpos, 0, wcs->ndec);
	    num2str (decstr, ypos, 0, wcs->ndec);
	    lstring = strlen (rastr) + strlen (decstr) + 1;
	    lunits = strlen (wcs->units[0]) + strlen (wcs->units[1]) + 2;
	    if (wcs->syswcs == eWCS_LINEAR && wcs->linmode == 1) {
		if (lstr > lstring + lunits) {
		    if (strlen (wcs->units[0]) > 0) {
			strcat (rastr, " ");
			strcat (rastr, wcs->units[0]);
			}
		    if (strlen (wcs->units[1]) > 0) {
			strcat (decstr, " ");
			strcat (decstr, wcs->units[1]);
			}
		    lstring = lstring + lunits;
		    }
		}
	    if (lstr > lstring) {
		if (wcs->tabsys)
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
		}
	    else {
		if (wcs->tabsys)
		    strncpy (wcstring,"**********	*********",lstr);
		else
		    strncpy (wcstring,"*******************",lstr);
		}
	    if (wcs->syswcs == eWCS_LINEAR && wcs->linmode != 1 &&
		lstr > lstring + 7)
		strcat (wcstring, " linear");
	    if (wcs->syswcs == eWCS_LINEAR && wcs->linmode == 2 &&
		lstr > lstring + lunits + 7) {
		if (strlen (wcs->units[0]) > 0) {
		    strcat (wcstring, " ");
		    strcat (wcstring, wcs->units[0]);
		    }
		if (strlen (wcs->units[1]) > 0) {
		    strcat (wcstring, " ");
		    strcat (wcstring, wcs->units[1]);
		    }
		    
		}
	    }

	return (1);

}//close pix2wcst()


void WCSUtils::pix2wcs(WCS* wcs,double xpix,double ypix,double* xpos,double* ypos)
{
	double	xpi, ypi, xp, yp;
  double	eqin, eqout;
  //int wcspos();

  if (nowcs (wcs))
		return;
    
	wcs->xpix = xpix;
  wcs->ypix = ypix;
  wcs->zpix = g_zpix;
  wcs->offscl = 0;

  // If this WCS is converted from another WCS rather than pixels, convert now
  if (wcs->wcs != NULL) {
		pix2wcs (wcs->wcs, xpix, ypix, &xpi, &ypi);
	}
  else {
		pix2foc (wcs, xpix, ypix, &xpi, &ypi);
	}

  // Convert image coordinates to sky coordinates

  // Use Digitized Sky Survey plate fit
  if (wcs->prjcode == eWCS_DSS) {
		if (dsspos (xpi, ypi, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

  // Use SAO plate fit
  else if (wcs->prjcode == eWCS_PLT) {
		if (platepos (xpi, ypi, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

  // Use NOAO IRAF corrected plane tangent projection 
  else if (wcs->prjcode == eWCS_TNX) {
		if (tnxpos (xpi, ypi, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

  // Use NOAO IRAF corrected zenithal projection
  else if (wcs->prjcode == eWCS_ZPX) {
		if (zpxpos (xpi, ypi, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

  // Use Classic AIPS projections
  else if (wcs->wcsproj == eWCS_OLD || wcs->prjcode <= 0) {
		if (worldpos (xpi, ypi, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

  //Use Mark Calabretta's WCSLIB projections
  else if (wcspos (xpi, ypi, wcs, &xp, &yp))
		wcs->offscl = 1;
	    	

  // Do not change coordinates if offscale
  if (wcs->offscl) {
		*xpos = 0.0;
		*ypos = 0.0;
	}
  else {

		//Convert coordinates to output system, if not LINEAR
    if (wcs->prjcode > 0) {

	    // Convert coordinates to desired output system
	    eqin = wcs->equinox;
	    eqout = wcs->eqout;
	    wcscon (wcs->syswcs,wcs->sysout,eqin,eqout,&xp,&yp,wcs->epoch);
	  }
		if (wcs->latbase == 90)
	  	yp = 90.0 - yp;
		else if (wcs->latbase == -90)
	    yp = yp - 90.0;
		wcs->xpos = xp;
		wcs->ypos = yp;
		*xpos = xp;
		*ypos = yp;
	}

  // Keep RA/longitude within range if spherical coordinate output (Not LINEAR or XY)
  if (wcs->sysout > 0 && wcs->sysout != 6 && wcs->sysout != 10) {
		if (*xpos < 0.0)
	    *xpos = *xpos + 360.0;
		else if (*xpos > 360.0)
	    *xpos = *xpos - 360.0;
	}
  
	return;

}//close pix2wcs()


void WCSUtils::wcs2pix (WCS* wcs,double xpos,double ypos,double* xpix,double* ypix,int* offscl)
{
  wcsc2pix (wcs, xpos, ypos, wcs->radecin, xpix, ypix, offscl);  
	return;

}//close wcs2pix()


void WCSUtils::wcsc2pix (WCS* wcs,double xpos,double ypos,char* coorsys,double* xpix,double* ypix,int* offscl)
{
	double xp, yp, xpi, ypi;
  double eqin, eqout;
  int sysin;
  //int wcspix();

  if (nowcs (wcs))
		return;

  *offscl = 0;
  xp = xpos;
  yp = ypos;
  if (wcs->latbase == 90)
		yp = 90.0 - yp;
  else if (wcs->latbase == -90)
		yp = yp - 90.0;
  if (coorsys == NULL) {
		sysin = wcs->syswcs;
		eqin = wcs->equinox;
	}
  else {
		sysin = wcscsys (coorsys);
		eqin = wcsceq (coorsys);
	}
  wcs->zpix = 1.0;

    /* Convert coordinates to same system as image */
    if (sysin > 0 && sysin != 6 && sysin != 10) {
	eqout = wcs->equinox;
	wcscon (sysin, wcs->syswcs, eqin, eqout, &xp, &yp, wcs->epoch);
	}

    /* Convert sky coordinates to image coordinates */

    /* Use Digitized Sky Survey plate fit */
    if (wcs->prjcode == eWCS_DSS) {
	if (dsspix (xp, yp, wcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use SAO polynomial plate fit */
    else if (wcs->prjcode == eWCS_PLT) {
	if (platepix (xp, yp, wcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use NOAO IRAF corrected plane tangent projection */
    else if (wcs->prjcode == eWCS_TNX) {
	if (tnxpix (xp, yp, wcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use NOAO IRAF corrected zenithal projection */
    else if (wcs->prjcode == eWCS_ZPX) {
	if (zpxpix (xp, yp, wcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use Classic AIPS projections */
    else if (wcs->wcsproj == eWCS_OLD || wcs->prjcode <= 0) {
	if (worldpix (xp, yp, wcs, &xpi, &ypi))
	    *offscl = 1;
	}

    /* Use Mark Calabretta's WCSLIB projections */
    else if (wcspix (xp, yp, wcs, &xpi, &ypi)) {
	*offscl = 1;
	}

    /* If this WCS is converted from another WCS rather than pixels, convert now */
    if (wcs->wcs != NULL) {
	wcsc2pix (wcs->wcs, xpi, ypi, NULL, xpix, ypix, offscl);
	}
    else {
	foc2pix (wcs, xpi, ypi, xpix, ypix);

	/* Set off-scale flag to 2 if off image but within bounds of projection */
	if (!*offscl) {
	    if (*xpix < 0.5 || *ypix < 0.5)
		*offscl = 2;
	    else if (*xpix > wcs->nxpix + 0.5 || *ypix > wcs->nypix + 0.5)
		*offscl = 2;
	    }
	}

  wcs->offscl = *offscl;
  wcs->xpos = xpos;
  wcs->ypos = ypos;
  wcs->xpix = *xpix;
  wcs->ypix = *ypix;

	return;

}//close wcsc2pix()


int WCSUtils::wcspix(double xpos,double ypos,WCS* wcs,double* xpix,double* ypix)
{
	int offscl;
  //int wcsfwd();
  double wcscrd[4], imgcrd[4], pixcrd[4];
  double phi, theta;

  *xpix = 0.0;
  *ypix = 0.0;
  if (wcs->wcsl.flag != kWCSSET) {
		//if (wcsset (wcs->lin.naxis, (void *)&wcs->ctype, &wcs->wcsl) )
		if (wcsset (wcs->lin.naxis, wcs->ctype, &wcs->wcsl) )
	  	return (1);
	}

    /* Set input for WCSLIB subroutines */
    wcscrd[0] = 0.0;
    wcscrd[1] = 0.0;
    wcscrd[2] = 0.0;
    wcscrd[3] = 0.0;
    wcscrd[wcs->wcsl.lng] = xpos;
    wcscrd[wcs->wcsl.lat] = ypos;

    /* Initialize output for WCSLIB subroutines */
    pixcrd[0] = 0.0;
    pixcrd[1] = 0.0;
    pixcrd[2] = 1.0;
    pixcrd[3] = 1.0;
    imgcrd[0] = 0.0;
    imgcrd[1] = 0.0;
    imgcrd[2] = 1.0;
    imgcrd[3] = 1.0;

    /* Invoke WCSLIB subroutines for coordinate conversion */
    //offscl = wcsfwd ((void *)&wcs->ctype, &wcs->wcsl, wcscrd, wcs->crval, &wcs->cel, &phi, &theta, &wcs->prj, imgcrd, &wcs->lin, pixcrd);
		offscl = wcsfwd( wcs->ctype, &wcs->wcsl, wcscrd, wcs->crval, &wcs->cel, &phi, &theta, &wcs->prj, imgcrd, &wcs->lin, pixcrd);

    if (!offscl) {
	*xpix = pixcrd[0];
	*ypix = pixcrd[1];
	if (wcs->prjcode == eWCS_CSC || wcs->prjcode == eWCS_QSC ||
	    wcs->prjcode == eWCS_TSC)
	    wcs->zpix = pixcrd[2] - 1.0;
	else
	    wcs->zpix = pixcrd[2];
	}

  return (offscl);

}//close wcspix()

int WCSUtils::wcsfwd(
	const char ctype[][16],
	wcsprm* wcs,	
	const double world[], 
	const double crval[], 
	celprm* cel, 
	double* phi,double* theta, 
	prjprm* prj, 
	double imgcrd[], 
	linprm* lin, 
	double pixcrd[]
)
{
   int    err, j;
   double offset;

   /* Initialize if required. */
   if (wcs->flag != kWCSSET) {
      if (wcsset(lin->naxis, ctype, wcs)) return 1;
   }

   /* Convert to relative physical coordinates. */
   for (j = 0; j < lin->naxis; j++) {
      if (j == wcs->lng) continue;
      if (j == wcs->lat) continue;
      imgcrd[j] = world[j] - crval[j];
   }

   if (wcs->flag != 999) {
      /* Compute projected coordinates. */
      if (strcmp(wcs->pcode, "NCP") == 0) {
         /* Convert NCP to SIN. */
         if (cel->ref[1] == 0.0) {
            return 2;
         }

         strcpy(wcs->pcode, "SIN");
         prj->p[1] = 0.0;
         prj->p[2] = cosdeg (cel->ref[1])/sindeg (cel->ref[1]);
         prj->flag = (prj->flag < 0) ? -1 : 0;
      }

      if ((err = celfwd(wcs->pcode, world[wcs->lng], world[wcs->lat], cel,
                   phi, theta, prj, &imgcrd[wcs->lng], &imgcrd[wcs->lat]))) {
         return err;
      }

      /* Do we have a CUBEFACE axis? */
      if (wcs->cubeface != -1) {
         /* Separation between faces. */
         if (prj->r0 == 0.0) {
            offset = 90.0;
         } else {
            offset = prj->r0*kPI/2.0;
         }

         /* Stack faces in a cube. */
         if (imgcrd[wcs->lat] < -0.5*offset) {
            imgcrd[wcs->lat] += offset;
            imgcrd[wcs->cubeface] = 5.0;
         } else if (imgcrd[wcs->lat] > 0.5*offset) {
            imgcrd[wcs->lat] -= offset;
            imgcrd[wcs->cubeface] = 0.0;
         } else if (imgcrd[wcs->lng] > 2.5*offset) {
            imgcrd[wcs->lng] -= 3.0*offset;
            imgcrd[wcs->cubeface] = 4.0;
         } else if (imgcrd[wcs->lng] > 1.5*offset) {
            imgcrd[wcs->lng] -= 2.0*offset;
            imgcrd[wcs->cubeface] = 3.0;
         } else if (imgcrd[wcs->lng] > 0.5*offset) {
            imgcrd[wcs->lng] -= offset;
            imgcrd[wcs->cubeface] = 2.0;
         } else {
            imgcrd[wcs->cubeface] = 1.0;
         }
      }
   }

   /* Apply forward linear transformation. */
   if (linfwd(imgcrd, lin, pixcrd)) {
      return 4;
   }

	return 0;

}//close wcsfwd()

int WCSUtils::celfwd(const char pcode[4],const double lng,const double lat,celprm* cel,double* phi,double* theta,prjprm* prj,double* x,double* y)
{
   int    err;

   if (cel->flag != kCELSET) {
      if (celset(pcode, cel, prj)) return 1;
   }

   /* Compute native coordinates. */
   sphfwd(lng, lat, cel->euler, phi, theta);

   /* Apply forward projection. */
   if ((err = prj->prjfwd(*phi, *theta, prj, x, y))) {
      return err == 1 ? 2 : 3;
   }

	return 0;

}//close celfwd()

int WCSUtils::linfwd(const double imgcrd[],linprm* lin, double pixcrd[])
{
   int i, ij, j, n;

   n = lin->naxis;

   if (lin->flag != kLINSET) {
      if (linset(lin)) return 1;
   }

   for (i = 0, ij = 0; i < n; i++) {
      pixcrd[i] = 0.0;
      for (j = 0; j < n; j++, ij++) {
         pixcrd[i] += lin->imgpix[ij] * imgcrd[j];
      }
   }

   for (j = 0; j < n; j++) {
      pixcrd[j] += lin->crpix[j];
   }
   
	return 0;

}//close linfwd()


int WCSUtils::sphfwd (const double lng, const double lat,const double eul[5],double* phi,double* theta)
{	
	const double tol = 1.0e-5;
  double coslat, coslng, dlng, dphi, sinlat, sinlng, x, y, z;

   coslat = cosdeg (lat);
   sinlat = sindeg (lat);

   dlng = lng - eul[0];
   coslng = cosdeg (dlng);
   sinlng = sindeg (dlng);

   /* Compute the native longitude. */
   x = sinlat*eul[4] - coslat*eul[3]*coslng;
   if (fabs(x) < tol) {
      /* Rearrange formula to reduce roundoff errors. */
      x = -cosdeg (lat+eul[1]) + coslat*eul[3]*(1.0 - coslng);
   }
   y = -coslat*sinlng;
   if (x != 0.0 || y != 0.0) {
      dphi = atan2deg (y, x);
   } else {
      /* Change of origin of longitude. */
      dphi = dlng - 180.0;
   }
   *phi = eul[2] + dphi;

   /* Normalize the native longitude. */
   if (*phi > 180.0) {
      *phi -= 360.0;
   } else if (*phi < -180.0) {
      *phi += 360.0;
   }

   /* Compute the native latitude. */
   if (fmod(dlng,180.0) == 0.0) {
      *theta = lat + coslng*eul[1];
      if (*theta >  90.0) *theta =  180.0 - *theta;
      if (*theta < -90.0) *theta = -180.0 - *theta;
   } else {
      z = sinlat*eul[3] + coslat*eul[4]*coslng;
      /* Use an alternative formula for greater numerical accuracy. */
      if (fabs(z) > 0.99) {
	if (z < 0)
           *theta = -acosdeg (sqrt(x*x+y*y));
	else
           *theta =  acosdeg (sqrt(x*x+y*y));
      } else {
         *theta = asindeg (z);
      }
   }

   return 0;

}//close sphfwd()


WCS* WCSUtils::wcskinit (
	int naxis1, int naxis2, 
	char* ctype1, char* ctype2, 
	double crpix1, double crpix2, 
	double crval1, double crval2,
	double* cd, 
	double cdelt1, double cdelt2, 
	double crota, 
	int equinox, 
	double epoch
)
{

	//Allocate WCS struct
	WCS* wcs = (WCS*)calloc (1, sizeof(WCS));
	//WCS* wcs = new WCS;

  // Set WCSLIB flags so that structures will be reinitialized 
  wcs->cel.flag = 0;
  wcs->lin.flag = 0;
  wcs->wcsl.flag = 0;

  //Image dimensions
  wcs->naxis = 2;
  wcs->naxes = 2;
  wcs->lin.naxis = 2;
  wcs->nxpix = naxis1;
  wcs->nypix = naxis2;

  wcs->wcsproj = g_wcsproj0;

  wcs->crpix[0] = crpix1;
  wcs->crpix[1] = crpix2;
  wcs->xrefpix = wcs->crpix[0];
  wcs->yrefpix = wcs->crpix[1];
  wcs->lin.crpix = wcs->crpix;

  if (wcstype (wcs, ctype1, ctype2)) {
		wcsfree (wcs);
		return (NULL);
	}
    
	if (wcs->latbase == 90)
		crval2 = 90.0 - crval2;
  else if (wcs->latbase == -90)
		crval2 = crval2 - 90.0;

  wcs->crval[0] = crval1;
  wcs->crval[1] = crval2;
  wcs->xref = wcs->crval[0];
  wcs->yref = wcs->crval[1];
  wcs->cel.ref[0] = wcs->crval[0];
  wcs->cel.ref[1] = wcs->crval[1];
  wcs->cel.ref[2] = 999.0;

  if (cd != NULL)
		wcscdset (wcs, cd);
	else if (cdelt1 != 0.0)
		wcsdeltset (wcs, cdelt1, cdelt2, crota);
	else {
		wcsdeltset (wcs, 1.0, 1.0, crota);
		setwcserr ("WCSRESET: setting CDELT to 1");
	}
    
	wcs->lin.cdelt = wcs->cdelt;
  wcs->lin.pc = wcs->pc;

  //Coordinate reference frame and equinox
  wcs->equinox =  (double) equinox;
  if (equinox > 1980)
		strcpy (wcs->radecsys,"FK5");
  else
		strcpy (wcs->radecsys,"FK4");
    
	if (epoch > 0)
		wcs->epoch = epoch;
  else
		wcs->epoch = 0.0;
    
	wcs->wcson = 1;

  strcpy (wcs->radecout, wcs->radecsys);
  wcs->syswcs = wcscsys (wcs->radecsys);
    
	wcsoutinit (wcs, wcs->radecsys);
  wcsininit (wcs, wcs->radecsys);
    
	wcs->eqout = 0.0;
  wcs->printsys = 1;
  wcs->tabsys = 0;

  // Initialize special WCS commands
  setwcscom (wcs);

  return (wcs);

}//close wcskinit()



void WCSUtils::wcsoutinit(WCS* wcs,char* coorsys)
{
	int sysout, i;

  if (nowcs (wcs))
		return;

  /* If argument is null, set to image system and equinox */
  if (coorsys == NULL || strlen (coorsys) < 1 ||
	!strcmp(coorsys,"IMSYS") || !strcmp(coorsys,"imsys")) {
	sysout = wcs->syswcs;
	strcpy (wcs->radecout, wcs->radecsys);
	wcs->eqout = wcs->equinox;
	if (sysout == eWCS_B1950) {
	    if (wcs->eqout != 1950.0) {
		wcs->radecout[0] = 'B';
		sprintf (wcs->radecout+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		}
	    else
		strcpy (wcs->radecout, "B1950");
	    }
	else if (sysout == eWCS_J2000) {
	    if (wcs->eqout != 2000.0) {
		wcs->radecout[0] = 'J';
		sprintf (wcs->radecout+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		}
	    else
		strcpy (wcs->radecout, "J2000");
	    }
	}

    /* Ignore new coordinate system if it is not supported */
    else {
	if ((sysout = wcscsys (coorsys)) < 0)
	return;

	/* Do not try to convert linear or alt-az coordinates */
	if (sysout != wcs->syswcs &&
	    (wcs->syswcs == eWCS_LINEAR || wcs->syswcs == eWCS_ALTAZ))
	    return;

	strcpy (wcs->radecout, coorsys);
	wcs->eqout = wcsceq (coorsys);
	}

    wcs->sysout = sysout;
    if (wcs->wcson) {

	/* Set output in degrees flag and number of decimal places */
	if (wcs->sysout == eWCS_GALACTIC || wcs->sysout == eWCS_ECLIPTIC ||
	    wcs->sysout == eWCS_PLANET) {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else if (wcs->sysout == eWCS_ALTAZ) {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else if (wcs->sysout == eWCS_NPOLE || wcs->sysout == eWCS_SPA) {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else {
	    wcs->degout = 0;
	    wcs->ndec = 3;
	    }
	}
    
	return;

}//close wcsoutinit()




void WCSUtils::wcsininit(WCS* wcs,char* coorsys)
{
	int sysin, i;

  if (nowcs (wcs))
		return;

  /* If argument is null, set to image system and equinox */
  if (coorsys == NULL || strlen (coorsys) < 1) {
	wcs->sysin = wcs->syswcs;
	strcpy (wcs->radecin, wcs->radecsys);
	wcs->eqin = wcs->equinox;
	if (wcs->sysin == eWCS_B1950) {
	    if (wcs->eqin != 1950.0) {
		wcs->radecin[0] = 'B';
		sprintf (wcs->radecin+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		}
	    else
		strcpy (wcs->radecin, "B1950");
	    }
	else if (wcs->sysin == eWCS_J2000) {
	    if (wcs->eqin != 2000.0) {
		wcs->radecin[0] = 'J';
		sprintf (wcs->radecin+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		}
	    else
		strcpy (wcs->radecin, "J2000");
	    }
	}

    /* Ignore new coordinate system if it is not supported */
    if ((sysin = wcscsys (coorsys)) < 0)
	return;

    wcs->sysin = sysin;
    wcs->eqin = wcsceq (coorsys);
    strcpy (wcs->radecin, coorsys);
    
	return;

}//close wcsininit()



void WCSUtils::setwcscom(WCS* wcs)
{
	char envar[16];
  int i;
  char *str;
  if (nowcs(wcs))
		return;
  
	for (i = 0; i < 10; i++) {
		if (i == 0)
	    strcpy (envar, "WCS_COMMAND");
		else
	    sprintf (envar, "WCS_COMMAND%d", i);
	if (g_wcscom0[i] != NULL)
	    wcscominit (wcs, i, g_wcscom0[i]);
	else if ((str = getenv (envar)) != NULL)
	    wcscominit (wcs, i, str);
	else if (i == 1)
	    wcscominit (wcs, i, "sua2 -ah %s");	/* F1= Search USNO-A2.0 Catalog */
	else if (i == 2)
	    wcscominit (wcs, i, "sgsc -ah %s");	/* F2= Search HST GSC */
	else if (i == 3)
	    wcscominit (wcs, i, "sty2 -ah %s"); /* F3= Search Tycho-2 Catalog */
	else if (i == 4)
	    wcscominit (wcs, i, "sppm -ah %s");	/* F4= Search PPM Catalog */
	else if (i == 5)
	    wcscominit (wcs, i, "ssao -ah %s");	/* F5= Search SAO Catalog */
	else
	    wcs->command_format[i] = NULL;
	}
    
	return;

}//close setwcscom()



void WCSUtils::wcscominit(WCS* wcs,int i,const char* command)
{
	int lcom,icom;

  if (iswcs(wcs)) {
	lcom = strlen (command);
	if (lcom > 0) {
	    if (wcs->command_format[i] != NULL)
		free (wcs->command_format[i]);
	    wcs->command_format[i] = (char *) calloc (lcom+2, 1);
	    if (wcs->command_format[i] == NULL)
		return;
	    for (icom = 0; icom < lcom; icom++) {
		if (command[icom] == '_')
		    wcs->command_format[i][icom] = ' ';
		else
		    wcs->command_format[i][icom] = command[icom];
		}
	    wcs->command_format[i][lcom] = 0;
	    }
	}
    
	return;

}//close wcscominit()

int WCSUtils::wcstype(WCS* wcs,char* ctype1,char* ctype2)
{
    int i, iproj;
    int nctype = eNWCSTYPE;
    char ctypes[eNWCSTYPE][4];
    char dtypes[10][4];

    /* Initialize projection types */
    strcpy (ctypes[0], "LIN");
    strcpy (ctypes[1], "AZP");
    strcpy (ctypes[2], "SZP");
    strcpy (ctypes[3], "TAN");
    strcpy (ctypes[4], "SIN");
    strcpy (ctypes[5], "STG");
    strcpy (ctypes[6], "ARC");
    strcpy (ctypes[7], "ZPN");
    strcpy (ctypes[8], "ZEA");
    strcpy (ctypes[9], "AIR");
    strcpy (ctypes[10], "CYP");
    strcpy (ctypes[11], "CAR");
    strcpy (ctypes[12], "MER");
    strcpy (ctypes[13], "CEA");
    strcpy (ctypes[14], "COP");
    strcpy (ctypes[15], "COD");
    strcpy (ctypes[16], "COE");
    strcpy (ctypes[17], "COO");
    strcpy (ctypes[18], "BON");
    strcpy (ctypes[19], "PCO");
    strcpy (ctypes[20], "SFL");
    strcpy (ctypes[21], "PAR");
    strcpy (ctypes[22], "AIT");
    strcpy (ctypes[23], "MOL");
    strcpy (ctypes[24], "CSC");
    strcpy (ctypes[25], "QSC");
    strcpy (ctypes[26], "TSC");
    strcpy (ctypes[27], "NCP");
    strcpy (ctypes[28], "GLS");
    strcpy (ctypes[29], "DSS");
    strcpy (ctypes[30], "PLT");
    strcpy (ctypes[31], "TNX");
    strcpy (ctypes[32], "ZPX");
    strcpy (ctypes[33], "TPV");

    /* Initialize distortion types */
    strcpy (dtypes[1], "SIP");

    if (!strncmp (ctype1, "LONG",4))
	strncpy (ctype1, "XLON",4);

    strcpy (wcs->ctype[0], ctype1);

    /* This is only to catch special non-standard projections */
    strncpy (wcs->ptype, ctype1, 3);
    wcs->ptype[3] = 0;

    /* Linear coordinates */
    if (!strncmp (ctype1,"LINEAR",6)) {
	wcs->prjcode = eWCS_LIN;
	strcpy (wcs->c1type, "LIN");
	strcpy (wcs->ptype, "LIN");
	}

    /* Pixel coordinates */
    else if (!strncmp (ctype1,"PIXEL",6)) {
	wcs->prjcode = eWCS_PIX;
	strcpy (wcs->c1type, "PIX");
	strcpy (wcs->ptype, "PIX");
	}

    /*Detector pixel coordinates */
    else if (strsrch (ctype1,"DET")) {
	wcs->prjcode = eWCS_PIX;
	strcpy (wcs->c1type, "PIX");
	strcpy (wcs->ptype, "PIX");
	}

    /* Set up right ascension, declination, latitude, or longitude */
    else if (ctype1[0] == 'R' || ctype1[0] == 'D' ||
	     ctype1[0] == 'A' || ctype1[1] == 'L') {
	wcs->c1type[0] = ctype1[0];
	wcs->c1type[1] = ctype1[1];
	if (ctype1[2] == '-') {
	    wcs->c1type[2] = 0;
	    iproj = 3;
	    }
	else {
	    wcs->c1type[2] = ctype1[2];
	    iproj = 4;
	    if (ctype1[3] == '-') {
		wcs->c1type[3] = 0;
		}
	    else {
		wcs->c1type[3] = ctype1[3];
		wcs->c1type[4] = 0;
		}
	    }
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	if (ctype1[iproj] == '-') iproj = iproj + 1;
	wcs->ptype[0] = ctype1[iproj];
	wcs->ptype[1] = ctype1[iproj+1];
	wcs->ptype[2] = ctype1[iproj+2];
	wcs->ptype[3] = 0;
	sprintf (wcs->ctype[0],"%-4s %3s",wcs->c1type,wcs->ptype);
	for (i = 0; i < 8; i++)
	    if (wcs->ctype[0][i] == ' ') wcs->ctype[0][i] = '-';

	/*  Find projection type  */
	wcs->prjcode = 0;  /* default type is linear */
	for (i = 1; i < nctype; i++) {
	    if (!strncmp(wcs->ptype, ctypes[i], 3))
		wcs->prjcode = i;
	    }

	/* Handle "obsolete" NCP projection (now WCSLIB should be OK)
	if (wcs->prjcode == eWCS_NCP) {
	    if (wcs->wcsproj == eWCS_BEST)
		wcs->wcsproj = eWCS_OLD;
	    else if (wcs->wcsproj == eWCS_ALT)
		wcs->wcsproj = eWCS_NEW;
	    } */

	/* Work around bug in WCSLIB handling of CAR projection
	else if (wcs->prjcode == eWCS_CAR) {
	    if (wcs->wcsproj == eWCS_BEST)
		wcs->wcsproj = eWCS_OLD;
	    else if (wcs->wcsproj == eWCS_ALT)
		wcs->wcsproj = eWCS_NEW;
	    } */

	/* Work around bug in WCSLIB handling of COE projection
	else if (wcs->prjcode == eWCS_COE) {
	    if (wcs->wcsproj == eWCS_BEST)
		wcs->wcsproj = eWCS_OLD;
	    else if (wcs->wcsproj == eWCS_ALT)
		wcs->wcsproj = eWCS_NEW;
	    }

	else if (wcs->wcsproj == eWCS_BEST) */
	if (wcs->wcsproj == eWCS_BEST)
	    wcs->wcsproj = eWCS_NEW;

	else if (wcs->wcsproj == eWCS_ALT)
	    wcs->wcsproj = eWCS_OLD;

	/* if (wcs->wcsproj == eWCS_OLD && (
	    wcs->prjcode != eWCS_STG && wcs->prjcode != eWCS_AIT &&
	    wcs->prjcode != eWCS_MER && wcs->prjcode != eWCS_GLS &&
	    wcs->prjcode != eWCS_ARC && wcs->prjcode != eWCS_TAN &&
	    wcs->prjcode != eWCS_TNX && wcs->prjcode != eWCS_SIN &&
	    wcs->prjcode != eWCS_PIX && wcs->prjcode != eWCS_LIN &&
	    wcs->prjcode != eWCS_CAR && wcs->prjcode != eWCS_COE &&
	    wcs->prjcode != eWCS_NCP && wcs->prjcode != eWCS_ZPX))
	    wcs->wcsproj = eWCS_NEW; */

	/* Handle NOAO corrected TNX as uncorrected TAN if oldwcs is set */
	if (wcs->wcsproj == eWCS_OLD && wcs->prjcode == eWCS_TNX) {
	    wcs->ctype[0][6] = 'A';
	    wcs->ctype[0][7] = 'N';
	    wcs->prjcode = eWCS_TAN;
	    }

	/* Handle NOAO corrected ZPX as uncorrected ZPN if oldwcs is set */
	if (wcs->wcsproj == eWCS_OLD && wcs->prjcode == eWCS_ZPX) {
	    wcs->ctype[0][6] = 'P';
	    wcs->ctype[0][7] = 'N';
	    wcs->prjcode = eWCS_ZPN;
	    }
	}

    /* If not sky coordinates, assume linear */
    else {
	wcs->prjcode = eWCS_LIN;
	strcpy (wcs->c1type, "LIN");
	strcpy (wcs->ptype, "LIN");
	return (0);
	}

    /* Second coordinate type */
    if (!strncmp (ctype2, "NPOL",4)) {
	ctype2[0] = ctype1[0];
	strncpy (ctype2+1, "LAT",3);
	wcs->latbase = 90;
	strcpy (wcs->radecsys,"NPOLE");
	wcs->syswcs = eWCS_NPOLE;
	}
    else if (!strncmp (ctype2, "SPA-",4)) {
	ctype2[0] = ctype1[0];
	strncpy (ctype2+1, "LAT",3);
	wcs->latbase = -90;
	strcpy (wcs->radecsys,"SPA");
	wcs->syswcs = eWCS_SPA;
	}
    else
	wcs->latbase = 0;
    strcpy (wcs->ctype[1], ctype2);

    /* Linear coordinates */
    if (!strncmp (ctype2,"LINEAR",6)) {
	wcs->prjcode = eWCS_LIN;
	strcpy (wcs->c2type, "LIN");
	}

    /* Pixel coordinates */
    else if (!strncmp (ctype2,"PIXEL",6)) {
	wcs->prjcode = eWCS_PIX;
	strcpy (wcs->c2type, "PIX");
	}

    /* Detector coordinates */
    else if (!strncmp (ctype2,"DET",3)) {
	wcs->prjcode = eWCS_PIX;
	strcpy (wcs->c2type, "PIX");
	}

    /* Set up right ascension, declination, latitude, or longitude */
    else if (ctype2[0] == 'R' || ctype2[0] == 'D' ||
	     ctype2[0] == 'A' || ctype2[1] == 'L') {
	wcs->c2type[0] = ctype2[0];
	wcs->c2type[1] = ctype2[1];
	if (ctype2[2] == '-') {
	    wcs->c2type[2] = 0;
	    iproj = 3;
	    }
	else {
	    wcs->c2type[2] = ctype2[2];
	    iproj = 4;
	    if (ctype2[3] == '-') {
		wcs->c2type[3] = 0;
		}
	    else {
		wcs->c2type[3] = ctype2[3];
		wcs->c2type[4] = 0;
		}
	    }
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	if (ctype2[iproj] == '-') iproj = iproj + 1;
	if (ctype2[iproj] == '-') iproj = iproj + 1;

	if (!strncmp (ctype1, "DEC", 3) ||
	    !strncmp (ctype1+1, "LAT", 3))
	    wcs->coorflip = 1;
	else
	    wcs->coorflip = 0;
	if (ctype2[1] == 'L' || ctype2[0] == 'A') {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else {
	    wcs->degout = 0;
	    wcs->ndec = 3;
	    }
	sprintf (wcs->ctype[1],"%-4s %3s",wcs->c2type,wcs->ptype);
	for (i = 0; i < 8; i++)
	    if (wcs->ctype[1][i] == ' ') wcs->ctype[1][i] = '-';
	}

    /* If not sky coordinates, assume linear */
    else {
	strcpy (wcs->c2type, "LIN");
	wcs->prjcode = eWCS_LIN;
	}

    /* Set distortion code from CTYPE1 extension */
    setdistcode (wcs, ctype1);
    
	return (0);

}//close wcstype()


int WCSUtils::wcscsys(char* wcstring)
{
	double equinox;

    if (wcstring[0] == 'J' || wcstring[0] == 'j' ||
	!strcmp (wcstring,"2000") || !strcmp (wcstring, "2000.0") ||
	!strcmp (wcstring,"ICRS") || !strcmp (wcstring, "icrs") ||
	!strncmp (wcstring,"FK5",3) || !strncmp (wcstring, "fk5",3))
	return eWCS_J2000;

    if (wcstring[0] == 'B' || wcstring[0] == 'b' ||
	!strcmp (wcstring,"1950") || !strcmp (wcstring, "1950.0") ||
	!strncmp (wcstring,"FK4",3) || !strncmp (wcstring, "fk4",3))
	return eWCS_B1950;

    else if (wcstring[0] == 'I' || wcstring[0] == 'i' )
	return eWCS_ICRS;

    else if (wcstring[0] == 'G' || wcstring[0] == 'g' )
	return eWCS_GALACTIC;

    else if (wcstring[0] == 'E' || wcstring[0] == 'e' )
	return eWCS_ECLIPTIC;

    else if (wcstring[0] == 'A' || wcstring[0] == 'a' )
	return eWCS_ALTAZ;

    else if (wcstring[0] == 'N' || wcstring[0] == 'n' )
	return eWCS_NPOLE;

    else if (wcstring[0] == 'L' || wcstring[0] == 'l' )
	return eWCS_LINEAR;

    else if (!strncasecmp (wcstring, "pixel", 5))
	return eWCS_XY;

    else if (wcstring[0] == 'P' || wcstring[0] == 'p' )
	return eWCS_PLANET;

    else if (isnum (wcstring) == 1 || isnum (wcstring) == 2) {
	equinox = atof (wcstring);
	if (equinox > 1980.0)
	    return eWCS_J2000;
	else if (equinox > 1900.0)
	    return eWCS_B1950;
	else
	    return -1;
	}
    else

	return -1;

}//close wcscsys()

void WCSUtils::wcsdeltset(WCS* wcs,double cdelt1,double cdelt2,double crota)
{
    double *pci;
    double crot, srot;
    int i, j, naxes;

    naxes = wcs->naxis;
    if (naxes > 2)
	naxes = 2;
    wcs->cdelt[0] = cdelt1;
    if (cdelt2 != 0.0)
	wcs->cdelt[1] = cdelt2;
    else
	wcs->cdelt[1] = cdelt1;
    wcs->xinc = wcs->cdelt[0];
    wcs->yinc = wcs->cdelt[1];
    pci = wcs->pc;
    for (i = 0; i < naxes; i++) {
	for (j = 0; j < naxes; j++) {
	    if (i ==j)
		*pci = 1.0;
	    else
		*pci = 0.0;
	    pci++;
	    }
	}
    wcs->rotmat = 0;

    /* If image is reversed, value of CROTA is flipped, too */
    wcs->rot = crota;
    if (wcs->rot < 0.0)
	wcs->rot = wcs->rot + 360.0;
    if (wcs->rot >= 360.0)
	wcs->rot = wcs->rot - 360.0;
    crot = cos (degrad_i(wcs->rot));
    if (cdelt1 * cdelt2 > 0)
	srot = sin (-degrad_i(wcs->rot));
    else
	srot = sin (degrad_i(wcs->rot));

    /* Set CD matrix */
    wcs->cd[0] = wcs->cdelt[0] * crot;
    if (wcs->cdelt[0] < 0)
	wcs->cd[1] = -fabs (wcs->cdelt[1]) * srot;
    else
	wcs->cd[1] = fabs (wcs->cdelt[1]) * srot;
    if (wcs->cdelt[1] < 0)
	wcs->cd[2] = fabs (wcs->cdelt[0]) * srot;
    else
	wcs->cd[2] = -fabs (wcs->cdelt[0]) * srot;
    wcs->cd[3] = wcs->cdelt[1] * crot;
    (void) matinv (2, wcs->cd, wcs->dc);

    /* Set rotation matrix */
    wcslibrot (wcs);

    /* Set image rotation and mirroring */
    if (wcs->coorflip) {
	if (wcs->cdelt[0] < 0 && wcs->cdelt[1] > 0) {
	    wcs->imflip = 1;
	    wcs->imrot = wcs->rot - 90.0;
	    if (wcs->imrot < -180.0) wcs->imrot = wcs->imrot + 360.0;
	    wcs->pa_north = wcs->rot;
	    wcs->pa_east = wcs->rot - 90.0;
	    if (wcs->pa_east < -180.0) wcs->pa_east = wcs->pa_east + 360.0;
	    }
	else if (wcs->cdelt[0] > 0 && wcs->cdelt[1] < 0) {
	    wcs->imflip = 1;
	    wcs->imrot = wcs->rot + 90.0;
	    if (wcs->imrot > 180.0) wcs->imrot = wcs->imrot - 360.0;
	    wcs->pa_north = wcs->rot;
	    wcs->pa_east = wcs->rot - 90.0;
	    if (wcs->pa_east < -180.0) wcs->pa_east = wcs->pa_east + 360.0;
	    }
	else if (wcs->cdelt[0] > 0 && wcs->cdelt[1] > 0) {
	    wcs->imflip = 0;
	    wcs->imrot = wcs->rot + 90.0;
	    if (wcs->imrot > 180.0) wcs->imrot = wcs->imrot - 360.0;
	    wcs->pa_north = wcs->imrot;
	    wcs->pa_east = wcs->rot + 90.0;
	    if (wcs->pa_east > 180.0) wcs->pa_east = wcs->pa_east - 360.0;
	    }
	else if (wcs->cdelt[0] < 0 && wcs->cdelt[1] < 0) {
	    wcs->imflip = 0;
	    wcs->imrot = wcs->rot - 90.0;
	    if (wcs->imrot < -180.0) wcs->imrot = wcs->imrot + 360.0;
	    wcs->pa_north = wcs->imrot;
	    wcs->pa_east = wcs->rot + 90.0;
	    if (wcs->pa_east > 180.0) wcs->pa_east = wcs->pa_east - 360.0;
	    }
	}
    else {
	if (wcs->cdelt[0] < 0 && wcs->cdelt[1] > 0) {
	    wcs->imflip = 0;
	    wcs->imrot = wcs->rot;
	    wcs->pa_north = wcs->rot + 90.0;
	    if (wcs->pa_north > 180.0) wcs->pa_north = wcs->pa_north - 360.0;
	    wcs->pa_east = wcs->rot + 180.0;
	    if (wcs->pa_east > 180.0) wcs->pa_east = wcs->pa_east - 360.0;
	    }
	else if (wcs->cdelt[0] > 0 && wcs->cdelt[1] < 0) {
	    wcs->imflip = 0;
	    wcs->imrot = wcs->rot + 180.0;
	    if (wcs->imrot > 180.0) wcs->imrot = wcs->imrot - 360.0;
	    wcs->pa_north = wcs->imrot + 90.0;
	    if (wcs->pa_north > 180.0) wcs->pa_north = wcs->pa_north - 360.0;
	    wcs->pa_east = wcs->imrot + 180.0;
	    if (wcs->pa_east > 180.0) wcs->pa_east = wcs->pa_east - 360.0;
	    }
	else if (wcs->cdelt[0] > 0 && wcs->cdelt[1] > 0) {
	    wcs->imflip = 1;
	    wcs->imrot = -wcs->rot;
	    wcs->pa_north = wcs->imrot + 90.0;
	    if (wcs->pa_north > 180.0) wcs->pa_north = wcs->pa_north - 360.0;
	    wcs->pa_east = wcs->rot;
	    }
	else if (wcs->cdelt[0] < 0 && wcs->cdelt[1] < 0) {
	    wcs->imflip = 1;
	    wcs->imrot = wcs->rot + 180.0;
	    if (wcs->imrot > 180.0) wcs->imrot = wcs->imrot - 360.0;
	    wcs->pa_north = wcs->imrot + 90.0;
	    if (wcs->pa_north > 180.0) wcs->pa_north = wcs->pa_north - 360.0;
	    wcs->pa_east = wcs->rot + 90.0;
	    if (wcs->pa_east > 180.0) wcs->pa_east = wcs->pa_east - 360.0;
	    }
	}

  return;

}//close wcsdeltset()


int WCSUtils::wcspos (double xpix, double ypix, WCS* wcs,double* xpos,double* ypos)
{
	int offscl;
  int i;
  //int wcsrev();
  double wcscrd[4], imgcrd[4], pixcrd[4];
  double phi, theta;
    
  *xpos = 0.0;
  *ypos = 0.0;

  pixcrd[0] = xpix;
  pixcrd[1] = ypix;
  if (wcs->prjcode == eWCS_CSC || wcs->prjcode == eWCS_QSC || wcs->prjcode == eWCS_TSC)
		pixcrd[2] = (double) (g_izpix + 1);
  else
		pixcrd[2] = g_zpix;
  pixcrd[3] = 1.0;
    
	for (i = 0; i < 4; i++)
		imgcrd[i] = 0.0;
  //offscl = wcsrev ((void *)&wcs->ctype, &wcs->wcsl, pixcrd, &wcs->lin, imgcrd,&wcs->prj, &phi, &theta, wcs->crval, &wcs->cel, wcscrd);
	offscl = wcsrev (wcs->ctype, &wcs->wcsl, pixcrd, &wcs->lin, imgcrd,&wcs->prj, &phi, &theta, wcs->crval, &wcs->cel, wcscrd);
  if (offscl == 0) {
		*xpos = wcscrd[wcs->wcsl.lng];
		*ypos = wcscrd[wcs->wcsl.lat];
	}

  return (offscl);

}//close wcspos()


void WCSUtils::pix2foc(WCS* wcs,double u,double v,double* x,double* y)
{
    int m, n, i, j, k;
    double s[kDISTMAX], sum;
    double temp_u, temp_v;

    /* Spitzer distortion */
    if (wcs->distcode == eDISTORT_SIRTF) {
	m = wcs->distort.a_order;
	n = wcs->distort.b_order;

	temp_u = u - wcs->xrefpix;
	temp_v = v - wcs->yrefpix;

	/* compute u */
	for (j = 0; j <= m; j++) {
	    s[j] = wcs->distort.a[m-j][j];
	    for (k = j-1; k >= 0; k--) {
		s[j] = (temp_v * s[j]) + wcs->distort.a[m-j][k];
		}
	    }
  
	sum = s[0];
	for (i=m; i>=1; i--){
	    sum = temp_u*sum + s[m-i+1];
	    }
	*x = sum;

	/* compute v*/
	for (j=0; j<=n; j++) {
	    s[j] = wcs->distort.b[n-j][j];
	    for (k=j-1; k>=0; k--) {
		s[j] =temp_v*s[j] + wcs->distort.b[n-j][k];
		}
	    }
   
	sum = s[0];
	for (i=n; i>=1; i--)
	    sum = temp_u*sum + s[n-i+1];

	*y = sum;
  
	*x = u + *x;
	*y = v + *y;

/*	*x = u + *x + coeff.crpix1; */
/*	*y = v + *y + coeff.crpix2; */
	}

    /* If no distortion, return pixel positions unchanged */
    else {
	*x = u;
	*y = v;
	}

  return;

}//close pix2foc()





void WCSUtils::wcscdset(WCS* wcs,double* cd)
{
	double tcd;
	if (cd == NULL)
		return;

	wcs->rotmat = 1;
  wcs->cd[0] = cd[0];
  wcs->cd[1] = cd[1];
  wcs->cd[2] = cd[2];
  wcs->cd[3] = cd[3];
  (void) matinv (2, wcs->cd, wcs->dc);

  // Compute scale
  wcs->xinc = sqrt (cd[0]*cd[0] + cd[2]*cd[2]);
  wcs->yinc = sqrt (cd[1]*cd[1] + cd[3]*cd[3]);

  // Deal with x=Dec/y=RA case 
  if (wcs->coorflip) {
		tcd = cd[1];
		cd[1] = -cd[2];
		cd[2] = -tcd;
	}
    
	wcslibrot (wcs);
  wcs->wcson = 1;

  //Compute image rotation
  wcsrotset (wcs);

  wcs->cdelt[0] = wcs->xinc;
  wcs->cdelt[1] = wcs->yinc;

	return;

}//close wcscdset()


void WCSUtils::wcslibrot (WCS* wcs)
{    
	int i, mem, naxes;

  naxes = wcs->naxis;
  if (naxes > 2)
		naxes = 2;
  if (naxes < 1 || naxes > 9) {
		naxes = wcs->naxes;
		wcs->naxis = naxes;
	}
  mem = naxes * naxes * sizeof(double);
  if (wcs->lin.piximg == NULL)
		wcs->lin.piximg = (double*)malloc(mem);
  if (wcs->lin.piximg != NULL) {
		if (wcs->lin.imgpix == NULL)
	  	wcs->lin.imgpix = (double*)malloc(mem);
		if (wcs->lin.imgpix != NULL) {
	    wcs->lin.flag = kLINSET;
	  	if (naxes == 2) {
				for (i = 0; i < 4; i++) {
			    wcs->lin.piximg[i] = wcs->cd[i];
			  }
			}
	  	else if (naxes == 3) {
				for (i = 0; i < 9; i++)
			    wcs->lin.piximg[i] = 0.0;
				wcs->lin.piximg[0] = wcs->cd[0];
				wcs->lin.piximg[1] = wcs->cd[1];
				wcs->lin.piximg[3] = wcs->cd[2];
				wcs->lin.piximg[4] = wcs->cd[3];
				wcs->lin.piximg[8] = 1.0;
			}
	  	else if (naxes == 4) {
				for (i = 0; i < 16; i++)
			    wcs->lin.piximg[i] = 0.0;
				wcs->lin.piximg[0] = wcs->cd[0];
				wcs->lin.piximg[1] = wcs->cd[1];
				wcs->lin.piximg[4] = wcs->cd[2];
				wcs->lin.piximg[5] = wcs->cd[3];
				wcs->lin.piximg[10] = 1.0;
				wcs->lin.piximg[15] = 1.0;
			}
	    
			(void) matinv (naxes, wcs->lin.piximg, wcs->lin.imgpix);
	  	wcs->lin.crpix = wcs->crpix;
	  	wcs->lin.cdelt = wcs->cdelt;
	  	wcs->lin.pc = wcs->pc;
	  	wcs->lin.flag = kLINSET;
		}
	}
    	
	return;

}//close wcslibrot()


void WCSUtils::wcsrotset(WCS* wcs)
{
	int off;
  double cra, cdec, xc, xn, xe, yc, yn, ye;

  // If image is one-dimensional, leave rotation angle alone
  if (wcs->nxpix < 1.5 || wcs->nypix < 1.5) {
		wcs->imrot = wcs->rot;
		wcs->pa_north = wcs->rot + 90.0;
		wcs->pa_east = wcs->rot + 180.0;
		return;
	}

	// Do not try anything if image is LINEAR (not Cartesian projection)
  if (wcs->syswcs == eWCS_LINEAR)
		return;

  wcs->xinc = fabs (wcs->xinc);
  wcs->yinc = fabs (wcs->yinc);

  // Compute position angles of North and East in image
  xc = wcs->xrefpix;
  yc = wcs->yrefpix;
   
	pix2wcs (wcs, xc, yc, &cra, &cdec);
    
	if (wcs->coorflip) {
		wcs2pix (wcs, cra+wcs->yinc, cdec, &xe, &ye, &off);
		wcs2pix (wcs, cra, cdec+wcs->xinc, &xn, &yn, &off);
	}
  else {
		wcs2pix (wcs, cra+wcs->xinc, cdec, &xe, &ye, &off);
		wcs2pix (wcs, cra, cdec+wcs->yinc, &xn, &yn, &off);
	}
    
	wcs->pa_north = raddeg_i (atan2 (yn-yc, xn-xc));
  if (wcs->pa_north < -90.0)
		wcs->pa_north = wcs->pa_north + 360.0;
  wcs->pa_east = raddeg_i (atan2 (ye-yc, xe-xc));
    
	if (wcs->pa_east < -90.0)
		wcs->pa_east = wcs->pa_east + 360.0;

  // Compute image rotation angle from North
  if (wcs->pa_north < -90.0)
		wcs->imrot = 270.0 + wcs->pa_north;
  else
		wcs->imrot = wcs->pa_north - 90.0;

  // Compute CROTA
  if (wcs->coorflip) {
		wcs->rot = wcs->imrot + 90.0;
		if (wcs->rot < 0.0)
	  	wcs->rot = wcs->rot + 360.0;
	}
  else
		wcs->rot = wcs->imrot;
  if (wcs->rot < 0.0)
		wcs->rot = wcs->rot + 360.0;
  if (wcs->rot >= 360.0)
		wcs->rot = wcs->rot - 360.0;

  // Set image mirror flag based on axis orientation
  wcs->imflip = 0;
  if (wcs->pa_east - wcs->pa_north < -80.0 &&
		wcs->pa_east - wcs->pa_north > -100.0)
	wcs->imflip = 1;
  if (wcs->pa_east - wcs->pa_north < 280.0 &&
		wcs->pa_east - wcs->pa_north > 260.0)
	wcs->imflip = 1;
  if (wcs->pa_north - wcs->pa_east > 80.0 &&
		wcs->pa_north - wcs->pa_east < 100.0)
	wcs->imflip = 1;
  if (wcs->coorflip) {
		if (wcs->imflip)
	    wcs->yinc = -wcs->yinc;
	}
  else {
		if (!wcs->imflip)
	  	wcs->xinc = -wcs->xinc;
	}

	return;

}//close wcsrotset()




void WCSUtils::wcsfree(WCS* wcs)
{
	//Free WCS structure if allocated but not filled
	if (nowcs (wcs)) {
		if (wcs) free (wcs);
		return;
	}

  // Free WCS on which this WCS depends
  if (wcs->wcs) {
		wcsfree (wcs->wcs);
		wcs->wcs = NULL;
	}

  freewcscom (wcs);
  if (wcs->wcsname != NULL)
		free (wcs->wcsname);
  if (wcs->lin.imgpix != NULL)
		free (wcs->lin.imgpix);
  if (wcs->lin.piximg != NULL)
		free (wcs->lin.piximg);
  if (wcs->inv_x != NULL)
		poly_end (wcs->inv_x);
  if (wcs->inv_y != NULL)
		poly_end (wcs->inv_y);
  free (wcs);
    
	return;

}//close wcsfree()


int WCSUtils::nowcs(WCS* wcs)
{
	if (wcs == NULL)
		return (1);
  else
		return (!wcs->wcson);

}//close nowcs()



int WCSUtils::iswcs(WCS* wcs)
{
  if (wcs == NULL)
		return (0);
  else
		return (wcs->wcson);

}//close iswcs()

void WCSUtils::freewcscom (WCS* wcs)
{
	int i;
  for (i = 0; i < 10; i++) {
		if (g_wcscom0[i] != NULL) {
	  	free (g_wcscom0[i]);
	    g_wcscom0[i] = NULL;
	  }
	}
    
	if (iswcs (wcs)) {
		for (i = 0; i < 10; i++) {
	  	if (wcs->command_format[i] != NULL) {
				free (wcs->command_format[i]);
			}
	  }
	}
    
	return;

}//close freewcscom()



char* WCSUtils::strsrch (const char* s1,const char* s2)
{
	int ls1;
  ls1 = strlen (s1);
  return (strnsrch (s1, s2, ls1));
}


char* WCSUtils::strnsrch (const char* s1,const char* s2, const int ls1)
{
    char *s,*s1e;
    char cfirst,clast;
    int i,ls2;

    /* Return null string if either pointer is NULL */
    if (s1 == NULL || s2 == NULL)
	return (NULL);

    /* A zero-length pattern is found in any string */
    ls2 = strlen (s2);
    if (ls2 ==0)
	return ((char *) s1);

    /* Only a zero-length string can be found in a zero-length string */
    if (ls1 ==0)
	return (NULL);

    cfirst = (char) s2[0];
    clast = (char) s2[ls2-1];
    s1e = (char *) s1 + (int) ls1 - ls2 + 1;
    s = (char *) s1;
    while (s < s1e) { 

	/* Search for first character in pattern string */
	if (*s == cfirst) {

	    /* If single character search, return */
	    if (ls2 == 1)
		return (s);

	    /* Search for last character in pattern string if first found */
	    if (s[ls2-1] == clast) {

		/* If two-character search, return */
		if (ls2 == 2)
		    return (s);

		/* If 3 or more characters, check for rest of search string */
		i = 1;
		while (i < ls2 && s[i] == s2[i])
		    i++;

		/* If entire string matches, return */
		if (i >= ls2)
		    return (s);
		}
	    }
	s++;
	}
    return (NULL);

}//close strnsrch()


int WCSUtils::isnum (const char* string)
{
    int lstr, i, nd, cl;
    char cstr, cstr1, cstr2;
    int fpcode;

    /* Return 0 if string is NULL */
    if (string == NULL)
	return (0);

    lstr = strlen (string);
    nd = 0;
    cl = 0;
    fpcode = 1;

    /* Return 0 if string starts with a D or E */
    cstr = string[0];
    if (cstr == 'D' || cstr == 'd' ||
	cstr == 'E' || cstr == 'e') {
	return (0);
	}

    /* Remove trailing spaces */
    while (string[lstr-1] == ' ')
	lstr--;

    /* Numeric strings contain 0123456789-+ and d or e for exponents */
    for (i = 0; i < lstr; i++) {
	cstr = string[i];
	if (cstr == '\n')
	    break;

	/* Ignore leading spaces */
	if (cstr == ' ' && nd == 0)
	    continue;

	if ((cstr < 48 || cstr > 57) &&
	    cstr != '+' && cstr != '-' &&
	    cstr != 'D' && cstr != 'd' &&
	    cstr != 'E' && cstr != 'e' &&
	    cstr != ':' && cstr != '.')
	    return (0);
	else if (cstr == '+' || cstr == '-') {
	    if (string[i+1] == '-' || string[i+1] == '+')
		return (0);
	    else if (i > 0) {
		cstr1 = string[i-1];
		cstr2 = string[i+1];
		if (cstr == '-' && cstr1 > 47 && cstr1 < 58 &&
		    cstr2 > 47 && cstr2 < 58)
		    return (4);
		else if (cstr1 != 'D' && cstr1 != 'd' &&
		    cstr1 != 'E' && cstr1 != 'e' &&
		    cstr1 != ':' && cstr1 != ' ')
		    return (0);
		}
	    }
	else if (cstr >= 47 && cstr <= 57)
	    nd++;

	/* Check for colon */
	else if (cstr == 58)
	    cl++;
	if (cstr=='.' || cstr=='d' || cstr=='e' || cstr=='d' || cstr=='e')
	    fpcode = 2;
	}
    if (nd > 0) {
	if (cl)
	    fpcode = 3;
	return (fpcode);
	}
    else
	
	return (0);

}//close isnum()

void WCSUtils::setdistcode(WCS* wcs,char* ctype)
{
    char *extension;
    int lctype;

    lctype = strlen (ctype);
    if (lctype < 9)
	wcs->distcode = eDISTORT_NONE;
    else {
	extension = ctype + 8;
	if (!strncmp (extension, "-SIP", 4))
	    wcs->distcode = eDISTORT_SIRTF;
	else
	    wcs->distcode = eDISTORT_NONE;
	}
    return;

}//close setdistcode()


int WCSUtils::dsspix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix)
{
  double div,xi,eta,x,y,xy,x2,y2,x2y,y2x,x3,y3,x4,y4,x2y2,cjunk,dx,dy;
  double sypos,cypos,syplate,cyplate,sxdiff,cxdiff;
  double f,fx,fy,g,gx,gy, xmm, ymm;
  double conr2s = 206264.8062470964;
  double tolerance = 0.0000005;
  int    max_iterations = 50;
  int    i;
  double xr, yr; 	/* position in radians */

  *xpix = 0.0;
  *ypix = 0.0;

/* Convert RA and Dec in radians to standard coordinates on a plate */
  xr = degrad_i (xpos);
  yr = degrad_i (ypos);
  sypos = sin (yr);
  cypos = cos (yr);
  if (wcs->plate_dec == 0.0)
    wcs->plate_dec = degrad_i (wcs->yref);
  syplate = sin (wcs->plate_dec);
  cyplate = cos (wcs->plate_dec);
  if (wcs->plate_ra == 0.0)
    wcs->plate_ra = degrad_i (wcs->yref);
  sxdiff = sin (xr - wcs->plate_ra);
  cxdiff = cos (xr - wcs->plate_ra);
  div = (sypos * syplate) + (cypos * cyplate * cxdiff);
  if (div == 0.0)
    return (1);
  xi = cypos * sxdiff * conr2s / div;
  eta = ((sypos * cyplate) - (cypos * syplate * cxdiff)) * conr2s / div;

/* Set initial value for x,y */
  if (wcs->plate_scale == 0.0)
    return (1);
  xmm = xi / wcs->plate_scale;
  ymm = eta / wcs->plate_scale;

/* Iterate by Newton's method */
  for (i = 0; i < max_iterations; i++) {

    /* X plate model */
    xy = xmm * ymm;
    x2 = xmm * xmm;
    y2 = ymm * ymm;
    x2y = x2 * ymm;
    y2x = y2 * xmm;
    x2y2 = x2 + y2;
    cjunk = x2y2 * x2y2;
    x3 = x2 * xmm;
    y3 = y2 * ymm;
    x4 = x2 * x2;
    y4 = y2 * y2;
    f = wcs->x_coeff[0]*xmm      + wcs->x_coeff[1]*ymm +
        wcs->x_coeff[2]          + wcs->x_coeff[3]*x2 +
        wcs->x_coeff[4]*xy       + wcs->x_coeff[5]*y2 +
        wcs->x_coeff[6]*x2y2     + wcs->x_coeff[7]*x3 +
        wcs->x_coeff[8]*x2y      + wcs->x_coeff[9]*y2x +
        wcs->x_coeff[10]*y3      + wcs->x_coeff[11]*xmm*x2y2 +
        wcs->x_coeff[12]*xmm*cjunk;
    /* magnitude and color terms ignored
      + wcs->x_coeff[13]*mag +
        wcs->x_coeff[14]*mag*mag   + wcs->x_coeff[15]*mag*mag*mag +
        wcs->x_coeff[16]*mag*xmm   + wcs->x_coeff[17]*mag*(x2+y2) +
        wcs->x_coeff[18]*mag*xmm*(x2+y2)  + wcs->x_coeff[19]*color;
    */

    /*  Derivative of X model wrt x */
    fx = wcs->x_coeff[0]           + wcs->x_coeff[3]*2.0*xmm +
         wcs->x_coeff[4]*ymm       + wcs->x_coeff[6]*2.0*xmm +
         wcs->x_coeff[7]*3.0*x2    + wcs->x_coeff[8]*2.0*xy +
         wcs->x_coeff[9]*y2        + wcs->x_coeff[11]*(3.0*x2+y2) +
         wcs->x_coeff[12]*(5.0*x4 +6.0*x2*y2+y4);
    /* magnitude and color terms ignored
         wcs->x_coeff[16]*mag      + wcs->x_coeff[17]*mag*2.0*xmm +
         wcs->x_coeff[18]*mag*(3.0*x2+y2);
    */

    /* Derivative of X model wrt y */
    fy = wcs->x_coeff[1]           + wcs->x_coeff[4]*xmm +
         wcs->x_coeff[5]*2.0*ymm   + wcs->x_coeff[6]*2.0*ymm +
         wcs->x_coeff[8]*x2        + wcs->x_coeff[9]*2.0*xy +
         wcs->x_coeff[10]*3.0*y2   + wcs->x_coeff[11]*2.0*xy +
         wcs->x_coeff[12]*4.0*xy*x2y2;
    /* magnitude and color terms ignored
         wcs->x_coeff[17]*mag*2.0*ymm +
         wcs->x_coeff[18]*mag*2.0*xy;
    */

    /* Y plate model */
    g = wcs->y_coeff[0]*ymm       + wcs->y_coeff[1]*xmm +
       wcs->y_coeff[2]            + wcs->y_coeff[3]*y2 +
       wcs->y_coeff[4]*xy         + wcs->y_coeff[5]*x2 +
       wcs->y_coeff[6]*x2y2       + wcs->y_coeff[7]*y3 +
       wcs->y_coeff[8]*y2x        + wcs->y_coeff[9]*x2y +
       wcs->y_coeff[10]*x3        + wcs->y_coeff[11]*ymm*x2y2 +
       wcs->y_coeff[12]*ymm*cjunk;
    /* magnitude and color terms ignored
       wcs->y_coeff[13]*mag        + wcs->y_coeff[14]*mag*mag +
       wcs->y_coeff[15]*mag*mag*mag + wcs->y_coeff[16]*mag*ymm +
       wcs->y_coeff[17]*mag*x2y2 +
       wcs->y_coeff[18]*mag*ymm*x2y2 + wcs->y_coeff[19]*color;
    */

    /* Derivative of Y model wrt x */
    gx = wcs->y_coeff[1]           + wcs->y_coeff[4]*ymm +
         wcs->y_coeff[5]*2.0*xmm   + wcs->y_coeff[6]*2.0*xmm +
         wcs->y_coeff[8]*y2       + wcs->y_coeff[9]*2.0*xy +
         wcs->y_coeff[10]*3.0*x2  + wcs->y_coeff[11]*2.0*xy +
         wcs->y_coeff[12]*4.0*xy*x2y2;
    /* magnitude and color terms ignored
         wcs->y_coeff[17]*mag*2.0*xmm +
         wcs->y_coeff[18]*mag*ymm*2.0*xmm;
    */

    /* Derivative of Y model wrt y */
    gy = wcs->y_coeff[0]            + wcs->y_coeff[3]*2.0*ymm +
         wcs->y_coeff[4]*xmm        + wcs->y_coeff[6]*2.0*ymm +
         wcs->y_coeff[7]*3.0*y2     + wcs->y_coeff[8]*2.0*xy +
         wcs->y_coeff[9]*x2         + wcs->y_coeff[11]*(x2+3.0*y2) +
         wcs->y_coeff[12]*(5.0*y4 + 6.0*x2*y2 + x4);
    /* magnitude and color terms ignored
         wcs->y_coeff[16]*mag       + wcs->y_coeff[17]*mag*2.0*ymm +
         wcs->y_coeff[18]*mag*(x2+3.0*y2);
    */

    f = f - xi;
    g = g - eta;
    dx = ((-f * gy) + (g * fy)) / ((fx * gy) - (fy * gx));
    dy = ((-g * fx) + (f * gx)) / ((fx * gy) - (fy * gx));
    xmm = xmm + dx;
    ymm = ymm + dy;
    if ((fabs(dx) < tolerance) && (fabs(dy) < tolerance)) break;
    }

/* Convert mm from plate center to plate pixels */
  if (wcs->x_pixel_size == 0.0 || wcs->y_pixel_size == 0.0)
    return (1);
  x = (wcs->ppo_coeff[2] - xmm*1000.0) / wcs->x_pixel_size;
  y = (wcs->ppo_coeff[5] + ymm*1000.0) / wcs->y_pixel_size;

/* Convert from plate pixels to image pixels */
  *xpix = x - wcs->x_pixel_offset + 1.0 - 0.5;
  *ypix = y - wcs->y_pixel_offset + 1.0 - 0.5;

/* If position is off of the image, return offscale code */
  if (*xpix < 0.5 || *xpix > wcs->nxpix+0.5)
    return -1;
  if (*ypix < 0.5 || *ypix > wcs->nypix+0.5)
    return -1;

  return 0;

}//close dsspix()

int WCSUtils::platepix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix)
{
    double xi,eta,x,y,xy,x2,y2,x2y,y2x,x3,y3,r2,dx,dy;
    double tdec,ctan,ccos,traoff, craoff, etar, xir;
    double f,fx,fy,g,gx,gy;
    double ra0, dec0, ra, dec;
    double tolerance = 0.0000005;
    int    max_iterations = 50;
    int    i;
    int	ncoeff1 = wcs->ncoeff1;
    int	ncoeff2 = wcs->ncoeff2;

    /* Convert RA and Dec in radians to standard coordinates on a plate */
    ra = degrad_i (xpos);
    dec = degrad_i (ypos);
    tdec = tan (dec);
    ra0 = degrad_i (wcs->crval[0]);
    dec0 = degrad_i (wcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    traoff = tan (ra - ra0);
    craoff = cos (ra - ra0);
    etar = (1.0 - ctan * craoff / tdec) / (ctan + (craoff / tdec));
    xir = traoff * ccos * (1.0 - (etar * ctan));
    xi = raddeg_i (xir);
    eta = raddeg_i (etar);

    /* Set initial value for x,y */
    x = xi * wcs->dc[0] + eta * wcs->dc[1];
    y = xi * wcs->dc[2] + eta * wcs->dc[3];

    /* if (wcs->x_coeff[1] == 0.0)
	x = xi - wcs->x_coeff[0];
    else
	x = (xi - wcs->x_coeff[0]) / wcs->x_coeff[1];
    if (wcs->y_coeff[2] == 0.0)
	y = eta - wcs->y_coeff[0];
    else
	y = (eta - wcs->y_coeff[0]) / wcs->y_coeff[2]; */

    /* Iterate by Newton's method */
    for (i = 0; i < max_iterations; i++) {

	/* X plate model */
	xy = x * y;
	x2 = x * x;
	y2 = y * y;
	x3 = x2 * x;
	y3 = y2 * y;
	x2y = x2 * y;
	y2x = y2 * x;
	r2 = x2 + y2;

	f = wcs->x_coeff[0]	+ wcs->x_coeff[1]*x +
	    wcs->x_coeff[2]*y	+ wcs->x_coeff[3]*x2 +
	    wcs->x_coeff[4]*y2	+ wcs->x_coeff[5]*xy;

	/*  Derivative of X model wrt x */
	fx = wcs->x_coeff[1]	+ wcs->x_coeff[3]*2.0*x +
	     wcs->x_coeff[5]*y;

	/* Derivative of X model wrt y */
	fy = wcs->x_coeff[2]	+ wcs->x_coeff[4]*2.0*y +
	     wcs->x_coeff[5]*x;

	if (ncoeff1 > 6) {
	    f = f + wcs->x_coeff[6]*x3	+ wcs->x_coeff[7]*y3;
	    fx = fx + wcs->x_coeff[6]*3.0*x2;
	    fy = fy + wcs->x_coeff[7]*3.0*y2;
	    }

	if (ncoeff1 > 8) {
	    f = f +
		wcs->x_coeff[8]*x2y	+ wcs->x_coeff[9]*y2x +
		wcs->x_coeff[10]*r2 + wcs->x_coeff[11]*x*r2 +
		wcs->x_coeff[12]*y*r2;

	    fx = fx +	wcs->x_coeff[8]*2.0*xy + 
			wcs->x_coeff[9]*y2 +
	 		wcs->x_coeff[10]*2.0*x +
			wcs->x_coeff[11]*(3.0*x2+y2) +
			wcs->x_coeff[12]*2.0*xy;

	    fy = fy +	wcs->x_coeff[8]*x2 +
			wcs->x_coeff[9]*2.0*xy +
			wcs->x_coeff[10]*2.0*y +
			wcs->x_coeff[11]*2.0*xy +
			wcs->x_coeff[12]*(3.0*y2+x2);
	    }

	/* Y plate model */
	g = wcs->y_coeff[0]	+ wcs->y_coeff[1]*x +
	    wcs->y_coeff[2]*y	+ wcs->y_coeff[3]*x2 +
	    wcs->y_coeff[4]*y2	+ wcs->y_coeff[5]*xy;

	/* Derivative of Y model wrt x */
	gx = wcs->y_coeff[1]	+ wcs->y_coeff[3]*2.0*x +
	     wcs->y_coeff[5]*y;

	/* Derivative of Y model wrt y */
	gy = wcs->y_coeff[2]	+ wcs->y_coeff[4]*2.0*y +
	     wcs->y_coeff[5]*x;

	if (ncoeff2 > 6) {
	    g = g + wcs->y_coeff[6]*x3	+ wcs->y_coeff[7]*y3;
	    gx = gx + wcs->y_coeff[6]*3.0*x2;
	    gy = gy + wcs->y_coeff[7]*3.0*y2;
	    }

	if (ncoeff2 > 8) {
	    g = g +
		wcs->y_coeff[8]*x2y	+ wcs->y_coeff[9]*y2x +
		wcs->y_coeff[10]*r2	+ wcs->y_coeff[11]*x*r2 +
		wcs->y_coeff[12]*y*r2;

	    gx = gx +	wcs->y_coeff[8]*2.0*xy + 
			wcs->y_coeff[9]*y2 +
	 		wcs->y_coeff[10]*2.0*x +
			wcs->y_coeff[11]*(3.0*x2+y2) +
			wcs->y_coeff[12]*2.0*xy;

	    gy = gy +	wcs->y_coeff[8]*x2 +
			wcs->y_coeff[9]*2.0*xy +
			wcs->y_coeff[10]*2.0*y +
			wcs->y_coeff[11]*2.0*xy +
			wcs->y_coeff[12]*(3.0*y2+x2);
	    }

	f = f - xi;
	g = g - eta;
	dx = ((-f * gy) + (g * fy)) / ((fx * gy) - (fy * gx));
	dy = ((-g * fx) + (f * gx)) / ((fx * gy) - (fy * gx));
	x = x + dx;
	y = y + dy;
	if ((fabs(dx) < tolerance) && (fabs(dy) < tolerance)) break;
	}

    /* Convert from plate pixels to image pixels */
    *xpix = x + wcs->crpix[0];
    *ypix = y + wcs->crpix[1];

    /* If position is off of the image, return offscale code */
    if (*xpix < 0.5 || *xpix > wcs->nxpix+0.5)
	return -1;
    if (*ypix < 0.5 || *ypix > wcs->nypix+0.5)
	return -1;

  return 0;

}//close platepix()


void WCSUtils::poly_end(polystruct *poly)
{
  
	if (poly) {
  	free(poly->coeff);
    free(poly->basis);
    free(poly->degree);
    free(poly->group);
    free(poly);
  }
  
}//close poly_end()


double WCSUtils::poly_func(polystruct *poly, double *pos)
{
   
	double xpol[kPOLY_MAXDIM+1];
  double *post, *xpolt, *basis, *coeff, xval;
  long double	val;
  int	expo[kPOLY_MAXDIM+1], gexpo[kPOLY_MAXDIM+1];
  int	*expot, *degree,*degreet, *group,*groupt, *gexpot,d,g,t, ndim;

	/* Prepare the vectors and counters */
  ndim = poly->ndim;
  basis = poly->basis;
  coeff = poly->coeff;
  group = poly->group;
  degree = poly->degree;
  if (ndim)
    {
    for (xpolt=xpol, expot=expo, post=pos, d=ndim; --d;)
      {
      *(++xpolt) = 1.0;
      *(++expot) = 0;
      }
    for (gexpot=gexpo, degreet=degree, g=poly->ngroup; g--;)
      *(gexpot++) = *(degreet++);
    if (gexpo[*group])
      gexpo[*group]--;
    }

/* The constant term is handled separately */
  val = *(coeff++);
  *(basis++) = 1.0;
  *expo = 1;
  *xpol = *pos;

/* Compute the rest of the polynom */
  for (t=poly->ncoeff; --t; )
    {
/*-- xpol[0] contains the current product of the x^n's */
    val += (*(basis++)=*xpol)**(coeff++);
/*-- A complex recursion between terms of the polynom speeds up computations */
/*-- Not too good for roundoff errors (prefer Horner's), but much easier for */
/*-- multivariate polynomials: this is why we use a long double accumulator */
    post = pos;
    groupt = group;
    expot = expo;
    xpolt = xpol;
    for (d=0; d<ndim; d++, groupt++)
      if (gexpo[*groupt]--)
        {
        ++*(expot++);
        xval = (*(xpolt--) *= *post);
        while (d--)
          *(xpolt--) = xval;
        break;
        }
      else
        {
        gexpo[*groupt] = *expot;
        *(expot++) = 0;
        *(xpolt++) = 1.0;
        post++;
        }
    }

  return (double)val;

}//close poly_func()


int WCSUtils::matinv(const int n, const double mat[], double inv[])
{
	register int i, ij, ik, j, k, kj, pj;
  int    itemp, mem, *mxl, *lxm, pivot;
  double colmax, *lu, *rowmax, dtemp;

  // Allocate memory for internal arrays
  mem = n * sizeof(int);
  if ((mxl = (int*)malloc(mem)) == (int*)0) return 1;
  if ((lxm = (int*)malloc(mem)) == (int*)0) {
  	free(mxl);
    return 1;
  }

  mem = n * sizeof(double);
  if ((rowmax = (double*)malloc(mem)) == (double*)0) {
  	free(mxl);
    free(lxm);
    return 1;
  }

  mem *= n;
  if ((lu = (double*)malloc(mem)) == (double*)0) {
  	free(mxl);
    free(lxm);
    free(rowmax);
    return 1;
  }


  // Initialize arrays
  for (i = 0, ij = 0; i < n; i++) {
  	//Vector which records row interchanges
    mxl[i] = i;

    rowmax[i] = 0.0;

    for (j = 0; j < n; j++, ij++) {
    	dtemp = fabs(mat[ij]);
      if (dtemp > rowmax[i]) rowmax[i] = dtemp;

      lu[ij] = mat[ij];
    }

    // A row of zeroes indicates a singular matrix
    if (rowmax[i] == 0.0) {
    	free(mxl);
      free(lxm);
      free(rowmax);
      free(lu);
      return 2;
    }
	}


  // Form the LU triangular factorization using scaled partial pivoting
  for (k = 0; k < n; k++) {
  	// Decide whether to pivot
    colmax = fabs(lu[k*n+k]) / rowmax[k];
    pivot = k;

    for (i = k+1; i < n; i++) {
    	ik = i*n + k;
      dtemp = fabs(lu[ik]) / rowmax[i];
      if (dtemp > colmax) {
      	colmax = dtemp;
        pivot = i;
      }
    }

    if (pivot > k) {
    	// We must pivot, interchange the rows of the design matrix
      for (j = 0, pj = pivot*n, kj = k*n; j < n; j++, pj++, kj++) {
      	dtemp = lu[pj];
        lu[pj] = lu[kj];
        lu[kj] = dtemp;
      }

      // Amend the vector of row maxima
      dtemp = rowmax[pivot];
      rowmax[pivot] = rowmax[k];
      rowmax[k] = dtemp;

      // Record the interchange for later use
      itemp = mxl[pivot];
      mxl[pivot] = mxl[k];
      mxl[k] = itemp;
 		}

    // Gaussian elimination
    for (i = k+1; i < n; i++) {
    	ik = i*n + k;

      // Nothing to do if lu[ik] is zero
      if (lu[ik] != 0.0) {
      	// Save the scaling factor
        lu[ik] /= lu[k*n+k];

        // Subtract rows
        for (j = k+1; j < n; j++) {
        	lu[i*n+j] -= lu[ik]*lu[k*n+j];
        }
     	}
   	}
	}


  // mxl[i] records which row of mat corresponds to row i of lu
  // lxm[i] records which row of lu  corresponds to row i of mat
  for (i = 0; i < n; i++) {
  	lxm[mxl[i]] = i;
  }


  // Determine the inverse matrix
  for (i = 0, ij = 0; i < n; i++) {
  	for (j = 0; j < n; j++, ij++) {
    	inv[ij] = 0.0;
    }
 	}

 	for (k = 0; k < n; k++) {
  	inv[lxm[k]*n+k] = 1.0;

    // Forward substitution
    for (i = lxm[k]+1; i < n; i++) {
    	for (j = lxm[k]; j < i; j++) {
      	inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
      }
    }

    // Backward substitution
    for (i = n-1; i >= 0; i--) {
    	for (j = i+1; j < n; j++) {
      	inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
      }
      inv[i*n+k] /= lu[i*n+i];
    }
	}

  free(mxl);
  free(lxm);
  free(rowmax);
  free(lu);

  return 0;

}//close matinv()


double WCSUtils::cosdeg(const double angle)
{
   double resid;

   resid = fabs(fmod(angle,360.0));
   if (resid == 0.0) {
      return 1.0;
   } else if (resid == 90.0) {
      return 0.0;
   } else if (resid == 180.0) {
      return -1.0;
   } else if (resid == 270.0) {
      return 0.0;
   }

   return cos(angle*kD2R);

}//close cosdeg()


double WCSUtils::sindeg(const double angle)
{
   double resid;

   resid = fmod(angle-90.0,360.0);
   if (resid == 0.0) {
      return 1.0;
   } else if (resid == 90.0) {
      return 0.0;
   } else if (resid == 180.0) {
      return -1.0;
   } else if (resid == 270.0) {
      return 0.0;
   }

   return sin(angle*kD2R);

}//close sindeg()


double WCSUtils::atan2deg(const double y,const double x)
{
   if (y == 0.0) {
      if (x >= 0.0) {
         return 0.0;
      } else if (x < 0.0) {
         return 180.0;
      }
   } else if (x == 0.0) {
      if (y > 0.0) {
         return 90.0;
      } else if (y < 0.0) {
         return -90.0;
      }
   }

   return atan2(y,x)*kR2D;

}//close atan2deg()


double WCSUtils::tandeg(const double angle)
{
   double resid;

   resid = fmod(angle,360.0);
   if (resid == 0.0 || fabs(resid) == 180.0) {
      return 0.0;
   } else if (resid == 45.0 || resid == 225.0) {
      return 1.0;
   } else if (resid == -135.0 || resid == -315.0) {
      return -1.0;
   }

   return tan(angle*kD2R);

}//close tandeg()


double WCSUtils::atandeg (const double v)
{
   if (v == -1.0) {
      return -45.0;
   } else if (v == 0.0) {
      return 0.0;
   } else if (v == 1.0) {
      return 45.0;
   }

   return atan(v)*kR2D;

}//close atandeg()



double WCSUtils::acosdeg(const double v)
{
   if (v >= 1.0) {
      if (v-1.0 <  kWCSTRIG_TOL) return 0.0;
   } else if (v == 0.0) {
      return 90.0;
   } else if (v <= -1.0) {
      if (v+1.0 > -kWCSTRIG_TOL) return 180.0;
   }

   return acos(v)*kR2D;

}//close acosdeg()

double WCSUtils::asindeg (const double v)
{
   if (v <= -1.0) {
      if (v+1.0 > -kWCSTRIG_TOL) return -90.0;
   } else if (v == 0.0) {
      return 0.0;
   } else if (v >= 1.0) {
      if (v-1.0 <  kWCSTRIG_TOL) return 90.0;
   }

	return asin(v)*kR2D;

}//close asindeg()

int WCSUtils::linrev(const double pixcrd[],linprm* lin,double imgcrd[])
{
   int i, ij, j, n;
   double temp;

   n = lin->naxis;

   if (lin->flag != kLINSET) {
      if (linset(lin)) return 1;
   }

   for (i = 0; i < n; i++) {
      imgcrd[i] = 0.0;
   }

   for (j = 0; j < n; j++) {
      temp = pixcrd[j] - lin->crpix[j];
      for (i = 0, ij = j; i < n; i++, ij+=n) {
         imgcrd[i] += lin->piximg[ij] * temp;
      }
   }

   return 0;

}//close linrev()


int WCSUtils::linset(linprm* lin)
{
   int i, ij, j, mem, n;

   n = lin->naxis;

   /* Allocate memory for internal arrays. */
   mem = n * n * sizeof(double);
   lin->piximg = (double*)malloc(mem);
   if (lin->piximg == (double*)0) return 1;

   lin->imgpix = (double*)malloc(mem);
   if (lin->imgpix == (double*)0) {
      free(lin->piximg);
      return 1;
   }

   /* Compute the pixel-to-image transformation matrix. */
   for (i = 0, ij = 0; i < n; i++) {
      for (j = 0; j < n; j++, ij++) {
         lin->piximg[ij] = lin->cdelt[i] * lin->pc[ij];
      }
   }

   /* Compute the image-to-pixel transformation matrix. */
   if (matinv(n, lin->piximg, lin->imgpix)) return 2;

   lin->flag = kLINSET;

   return 0;

}//close linset()

int WCSUtils::wcsrev(const char ctype[][16],wcsprm* wcs,const double pixcrd[],linprm* lin,double imgcrd[],prjprm* prj,double* phi,double* theta,const double crval[],celprm* cel,double world[])
{
	int    err, face, j;
  double offset;

  // Initialize if required
  if (wcs->flag != kWCSSET) {
  	if (wcsset(lin->naxis, ctype, wcs)) return 1;
  }

  // Apply reverse linear transformation
  if (linrev(pixcrd, lin, imgcrd)) {
  	return 4;
  }

  // Convert to world coordinates
 	for (j = 0; j < lin->naxis; j++) {
  	if (j == wcs->lng) continue;
    if (j == wcs->lat) continue;
    world[j] = imgcrd[j] + crval[j];
  }

  if (wcs->flag != 999) {
  	// Do we have a CUBEFACE axis?
    if (wcs->cubeface != -1) {
    	face = (int)(imgcrd[wcs->cubeface] + 0.5);
      if (fabs(imgcrd[wcs->cubeface]-face) > 1e-10) {
      	return 3;
      }

      // Separation between faces
      if (prj->r0 == 0.0) {
      	offset = 90.0;
      } 
			else {
      	offset = prj->r0*kPI/2.0;
      }

      // Lay out faces in a plane
      switch (face) {
      case 0:
      	imgcrd[wcs->lat] += offset;
        break;
         case 1:
            break;
         case 2:
            imgcrd[wcs->lng] += offset;
            break;
         case 3:
            imgcrd[wcs->lng] += offset*2;
            break;
         case 4:
            imgcrd[wcs->lng] += offset*3;
            break;
         case 5:
            imgcrd[wcs->lat] -= offset;
            break;
         default:
            return 3;
         }
      }

      //Compute celestial coordinates
      if (strcmp(wcs->pcode, "NCP") == 0) {
      	//Convert NCP to SIN
       	if (cel->ref[1] == 0.0) {
        	return 2;
        }

        strcpy(wcs->pcode, "SIN");
        prj->p[1] = 0.0;
        prj->p[2] = cosdeg (cel->ref[1])/sindeg (cel->ref[1]);
        prj->flag = (prj->flag < 0) ? -1 : 0;
      }

      if ((err = celrev(wcs->pcode, imgcrd[wcs->lng], imgcrd[wcs->lat], prj, phi, theta, cel, &world[wcs->lng], &world[wcs->lat]))) {
      	return err;
      }
   }

   return 0;

}//close wcsrev()


int WCSUtils::celrev(const char pcode[4],const double x, const double y,prjprm* prj,double* phi,double* theta,celprm* cel,double* lng,double* lat)
{
   int    err;

   if (cel->flag != kCELSET) {
      if(celset(pcode, cel, prj)) return 1;
   }

   /* Apply reverse projection. */
   if ((err = prj->prjrev(x, y, prj, phi, theta))) {
      return err == 1 ? 2 : 3;
   }

   /* Compute native coordinates. */
   sphrev(*phi, *theta, cel->euler, lng, lat);

   return 0;

}//close celrev()


int WCSUtils::prjset(const char pcode[4], prjprm* prj)
{
   /* Set pointers to the forward and reverse projection routines. */
   if (strcmp(pcode, "AZP") == 0) {
      azpset(prj);
   } else if (strcmp(pcode, "SZP") == 0) {
      szpset(prj);
   } else if (strcmp(pcode, "TAN") == 0) {
      tanset(prj);
   } else if (strcmp(pcode, "STG") == 0) {
      stgset(prj);
   } else if (strcmp(pcode, "SIN") == 0) {
      sinset(prj);
   } else if (strcmp(pcode, "ARC") == 0) {
      arcset(prj);
   } else if (strcmp(pcode, "ZPN") == 0) {
      zpnset(prj);
   } else if (strcmp(pcode, "ZEA") == 0) {
      zeaset(prj);
   } else if (strcmp(pcode, "AIR") == 0) {
      airset(prj);
   } else if (strcmp(pcode, "CYP") == 0) {
      cypset(prj);
   } else if (strcmp(pcode, "CEA") == 0) {
      ceaset(prj);
   } else if (strcmp(pcode, "CAR") == 0) {
      carset(prj);
   } else if (strcmp(pcode, "MER") == 0) {
      merset(prj);
   } else if (strcmp(pcode, "SFL") == 0) {
      sflset(prj);
   } else if (strcmp(pcode, "PAR") == 0) {
      parset(prj);
   } else if (strcmp(pcode, "MOL") == 0) {
      molset(prj);
   } else if (strcmp(pcode, "AIT") == 0) {
      aitset(prj);
   } else if (strcmp(pcode, "COP") == 0) {
      copset(prj);
   } else if (strcmp(pcode, "COE") == 0) {
      coeset(prj);
   } else if (strcmp(pcode, "COD") == 0) {
      codset(prj);
   } else if (strcmp(pcode, "COO") == 0) {
      cooset(prj);
   } else if (strcmp(pcode, "BON") == 0) {
      bonset(prj);
   } else if (strcmp(pcode, "PCO") == 0) {
      pcoset(prj);
   } else if (strcmp(pcode, "TSC") == 0) {
      tscset(prj);
   } else if (strcmp(pcode, "CSC") == 0) {
      cscset(prj);
   } else if (strcmp(pcode, "QSC") == 0) {
      qscset(prj);
   } else {
      /* Unrecognized projection code. */
      return 1;
   }

	return 0;

}//close prjset()


int WCSUtils::celset(const char pcode[4],celprm* cel,prjprm* prj)
{
   int dophip;
   const double tol = 1.0e-10;
   double clat0, cphip, cthe0, slat0, sphip, sthe0;
   double latp, latp1, latp2;
   double u, v, x, y, z;

   /* Initialize the projection driver routines. */
   if (prjset(pcode, prj)) {
      return 1;
   }

   /* Set default for native longitude of the celestial pole? */
   dophip = (cel->ref[2] == 999.0);

   /* Compute celestial coordinates of the native pole. */
   if (prj->theta0 == 90.0) {
      /* Reference point is at the native pole. */

      if (dophip) {
         /* Set default for longitude of the celestial pole. */
         cel->ref[2] = 180.0;
      }

      latp = cel->ref[1];
      cel->ref[3] = latp;

      cel->euler[0] = cel->ref[0];
      cel->euler[1] = 90.0 - latp;
   } else {
      /* Reference point away from the native pole. */

      /* Set default for longitude of the celestial pole. */
      if (dophip) {
         cel->ref[2] = (cel->ref[1] < prj->theta0) ? 180.0 : 0.0;
      }

      clat0 = cosdeg (cel->ref[1]);
      slat0 = sindeg (cel->ref[1]);
      cphip = cosdeg (cel->ref[2]);
      sphip = sindeg (cel->ref[2]);
      cthe0 = cosdeg (prj->theta0);
      sthe0 = sindeg (prj->theta0);

      x = cthe0*cphip;
      y = sthe0;
      z = sqrt(x*x + y*y);
      if (z == 0.0) {
         if (slat0 != 0.0) {
            return 1;
         }

         /* latp determined by LATPOLE in this case. */
         latp = cel->ref[3];
      } else {
         if (fabs(slat0/z) > 1.0) {
            return 1;
         }

         u = atan2deg (y,x);
         v = acosdeg (slat0/z);

         latp1 = u + v;
         if (latp1 > 180.0) {
            latp1 -= 360.0;
         } else if (latp1 < -180.0) {
            latp1 += 360.0;
         }

         latp2 = u - v;
         if (latp2 > 180.0) {
            latp2 -= 360.0;
         } else if (latp2 < -180.0) {
            latp2 += 360.0;
         }

         if (fabs(cel->ref[3]-latp1) < fabs(cel->ref[3]-latp2)) {
            if (fabs(latp1) < 90.0+tol) {
               latp = latp1;
            } else {
               latp = latp2;
            }
         } else {
            if (fabs(latp2) < 90.0+tol) {
               latp = latp2;
            } else {
               latp = latp1;
            }
         }

         cel->ref[3] = latp;
      }

      cel->euler[1] = 90.0 - latp;

      z = cosdeg (latp)*clat0;
      if (fabs(z) < tol) {
         if (fabs(clat0) < tol) {
            /* Celestial pole at the reference point. */
            cel->euler[0] = cel->ref[0];
            cel->euler[1] = 90.0 - prj->theta0;
         } else if (latp > 0.0) {
            /* Celestial pole at the native north pole.*/
            cel->euler[0] = cel->ref[0] + cel->ref[2] - 180.0;
            cel->euler[1] = 0.0;
         } else if (latp < 0.0) {
            /* Celestial pole at the native south pole. */
            cel->euler[0] = cel->ref[0] - cel->ref[2];
            cel->euler[1] = 180.0;
         }
      } else {
         x = (sthe0 - sindeg (latp)*slat0)/z;
         y =  sphip*cthe0/clat0;
         if (x == 0.0 && y == 0.0) {
            return 1;
         }
         cel->euler[0] = cel->ref[0] - atan2deg (y,x);
      }

      /* Make euler[0] the same sign as ref[0]. */
      if (cel->ref[0] >= 0.0) {
         if (cel->euler[0] < 0.0) cel->euler[0] += 360.0;
      } else {
         if (cel->euler[0] > 0.0) cel->euler[0] -= 360.0;
      }
   }

   cel->euler[2] = cel->ref[2];
   cel->euler[3] = cosdeg (cel->euler[1]);
   cel->euler[4] = sindeg (cel->euler[1]);
   cel->flag = kCELSET;

   /* Check for ill-conditioned parameters. */
   if (fabs(latp) > 90.0+tol) {
      return 2;
   }

   return 0;
}//close celset()




int WCSUtils::wcsset (const int naxis, const char ctype[][16],wcsprm* wcs)

/*
const int naxis;
const char ctype[][16];
struct wcsprm *wcs;
*/
{
   int  nalias = 2;
   char aliases [2][4] = {"NCP", "GLS"};

   int j, k;
   int *ndx = NULL;
   char requir[16];

   strcpy(wcs->pcode, "");
   strcpy(requir, "");
   wcs->lng = -1;
   wcs->lat = -1;
   wcs->cubeface = -1;

   for (j = 0; j < naxis; j++) {
      if (ctype[j][4] != '-') {
         if (strcmp(ctype[j], "CUBEFACE") == 0) {
            if (wcs->cubeface == -1) {
               wcs->cubeface = j;
            } else {
               /* Multiple CUBEFACE axes! */
               return 1;
            }
         }
         continue;
      }

      /* Got an axis qualifier, is it a recognized WCS projection? */
      for (k = 0; k < g_npcode; k++) {
         if (strncmp(&ctype[j][5], g_pcodes[k], 3) == 0) break;
      }

      if (k == g_npcode) {
         /* Maybe it's a projection alias. */
         for (k = 0; k < nalias; k++) {
            if (strncmp(&ctype[j][5], aliases[k], 3) == 0) break;
         }

         /* Not recognized. */
         if (k == nalias) {
            continue;
         }
      }

      /* Parse the celestial axis type. */
      if (strcmp(wcs->pcode, "") == 0) {
         sprintf(wcs->pcode, "%.3s", &ctype[j][5]);

         if (strncmp(ctype[j], "RA--", 4) == 0) {
            wcs->lng = j;
            strcpy(wcs->lngtyp, "RA");
            strcpy(wcs->lattyp, "DEC");
            ndx = &wcs->lat;
            sprintf(requir, "DEC--%s", wcs->pcode);
         } else if (strncmp(ctype[j], "DEC-", 4) == 0) {
            wcs->lat = j;
            strcpy(wcs->lngtyp, "RA");
            strcpy(wcs->lattyp, "DEC");
            ndx = &wcs->lng;
            sprintf(requir, "RA---%s", wcs->pcode);
         } else if (strncmp(&ctype[j][1], "LON", 3) == 0) {
            wcs->lng = j;
            sprintf(wcs->lngtyp, "%cLON", ctype[j][0]);
            sprintf(wcs->lattyp, "%cLAT", ctype[j][0]);
            ndx = &wcs->lat;
            sprintf(requir, "%s-%s", wcs->lattyp, wcs->pcode);
         } else if (strncmp(&ctype[j][1], "LAT", 3) == 0) {
            wcs->lat = j;
            sprintf(wcs->lngtyp, "%cLON", ctype[j][0]);
            sprintf(wcs->lattyp, "%cLAT", ctype[j][0]);
            ndx = &wcs->lng;
            sprintf(requir, "%s-%s", wcs->lngtyp, wcs->pcode);
         } else if (strncmp(&ctype[j][2], "LN", 2) == 0) {
            wcs->lng = j;
            sprintf(wcs->lngtyp, "%c%cLN", ctype[j][0], ctype[j][1]);
            sprintf(wcs->lattyp, "%c%cLT", ctype[j][0], ctype[j][1]);
            ndx = &wcs->lat;
            sprintf(requir, "%s-%s", wcs->lattyp, wcs->pcode);
         } else if (strncmp(&ctype[j][2], "LT", 2) == 0) {
            wcs->lat = j;
            sprintf(wcs->lngtyp, "%c%cLN", ctype[j][0], ctype[j][1]);
            sprintf(wcs->lattyp, "%c%cLT", ctype[j][0], ctype[j][1]);
            ndx = &wcs->lng;
            sprintf(requir, "%s-%s", wcs->lngtyp, wcs->pcode);
         } else {
            /* Unrecognized celestial type. */
            return 1;
         }
      } else {
         if (strncmp(ctype[j], requir, 8) != 0) {
            /* Inconsistent projection types. */
            return 1;
         }

	if (ndx == NULL)
	    return 1;
         *ndx = j;
         strcpy(requir, "");
      }
   }

   if (strcmp(requir, "")) {
      /* Unmatched celestial axis. */
      return 1;
   }

   /* Do simple alias translations. */
   if (strncmp(wcs->pcode, "GLS", 3) == 0) {
      strcpy(wcs->pcode, "SFL");
   }

   if (strcmp(wcs->pcode, "")) {
      wcs->flag = kWCSSET;
   } else {
      /* Signal for no celestial axis pair. */
      wcs->flag = 999;
   }

   return 0;

}//close wcsset()


int WCSUtils::parset(prjprm * prj)
{
   strcpy(prj->code, "PAR");
   prj->flag   = PAR;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
      prj->w[2] = 180.0;
      prj->w[3] = 1.0/prj->w[2];
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = 1.0/prj->w[0];
      prj->w[2] = kPI*prj->r0;
      prj->w[3] = 1.0/prj->w[2];
   }

   prj->prjfwd = parfwd;
   prj->prjrev = parrev;

   return 0;

}//close parset()

int WCSUtils::parfwd(const double phi,const double theta,prjprm * prj,double* x,double* y)
{
   double s;

   if (prj->flag != PAR) {
      if (parset(prj)) return 1;
   }

   s = sindeg (theta/3.0);
   *x = prj->w[0]*phi*(1.0 - 4.0*s*s);
   *y = prj->w[2]*s;

   return 0;

}//close parfwd()


int WCSUtils::parrev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double s, t;

   if (prj->flag != PAR) {
      if (parset(prj)) return 1;
   }

   s = y*prj->w[3];
   if (s > 1.0 || s < -1.0) {
      return 2;
   }

   t = 1.0 - 4.0*s*s;
   if (t == 0.0) {
      if (x == 0.0) {
         *phi = 0.0;
      } else {
         return 2;
      }
   } else {
      *phi = prj->w[1]*x/t;
   }

   *theta = 3.0*asindeg (s);

   return 0;

}//close parrev()

int WCSUtils::molset(prjprm * prj)
{
   strcpy(prj->code, "MOL");
   prj->flag   = MOL;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = kSQRT2*prj->r0;
   prj->w[1] = prj->w[0]/90.0;
   prj->w[2] = 1.0/prj->w[0];
   prj->w[3] = 90.0/prj->r0;
   prj->w[4] = 2.0/kPI;

   prj->prjfwd = molfwd;
   prj->prjrev = molrev;

   return 0;

}//close molset()


int WCSUtils::molfwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   int   j;
   double gamma, resid, u, v, v0, v1;
   const double tol = 1.0e-13;

   if (prj->flag != MOL) {
      if (molset(prj)) return 1;
   }

   if (fabs(theta) == 90.0) {
      *x = 0.0;
      *y = copysgn_i (prj->w[0],theta);
   } else if (theta == 0.0) {
      *x = prj->w[1]*phi;
      *y = 0.0;
   } else {
      u  = kPI*sindeg (theta);
      v0 = -kPI;
      v1 =  kPI;
      v  = u;
      for (j = 0; j < 100; j++) {
         resid = (v - u) + sin(v);
         if (resid < 0.0) {
            if (resid > -tol) break;
            v0 = v;
         } else {
            if (resid < tol) break;
            v1 = v;
         }
         v = (v0 + v1)/2.0;
      }

      gamma = v/2.0;
      *x = prj->w[1]*phi*cos(gamma);
      *y = prj->w[0]*sin(gamma);
   }

   return 0;

}//close molfwd()

int WCSUtils::molrev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double s, y0, z;
   const double tol = 1.0e-12;

   if (prj->flag != MOL) {
      if (molset(prj)) return 1;
   }

   y0 = y/prj->r0;
   s  = 2.0 - y0*y0;
   if (s <= tol) {
      if (s < -tol) {
         return 2;
      }
      s = 0.0;

      if (fabs(x) > tol) {
         return 2;
      }
      *phi = 0.0;
   } else {
      s = sqrt(s);
      *phi = prj->w[3]*x/s;
   }

   z = y*prj->w[2];
   if (fabs(z) > 1.0) {
      if (fabs(z) > 1.0+tol) {
         return 2;
      }
      z = copysgn_i (1.0,z) + y0*s/kPI;
   } else {
      z = asin(z)*prj->w[4] + y0*s/kPI;
   }

   if (fabs(z) > 1.0) {
      if (fabs(z) > 1.0+tol) {
         return 2;
      }
      z = copysgn_i (1.0,z);
   }

   *theta = asindeg (z);

   return 0;

}//close molrev()

int WCSUtils::aitset(prjprm *prj)
{
   strcpy(prj->code, "AIT");
   prj->flag   = AIT;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = 2.0*prj->r0*prj->r0;
   prj->w[1] = 1.0/(2.0*prj->w[0]);
   prj->w[2] = prj->w[1]/4.0;
   prj->w[3] = 1.0/(2.0*prj->r0);

   prj->prjfwd = aitfwd;
   prj->prjrev = aitrev;

   return 0;

}//close aitset()

int WCSUtils::aitfwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double cthe, w;

   if (prj->flag != AIT) {
      if (aitset(prj)) return 1;
   }

   cthe = cosdeg (theta);
   w = sqrt(prj->w[0]/(1.0 + cthe*cosdeg (phi/2.0)));
   *x = 2.0*w*cthe*sindeg (phi/2.0);
   *y = w*sindeg (theta);

   return 0;

}//close aitfwd()


int WCSUtils::aitrev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double s, u, xp, yp, z;
   const double tol = 1.0e-13;

   if (prj->flag != AIT) {
      if (aitset(prj)) return 1;
   }

   u = 1.0 - x*x*prj->w[2] - y*y*prj->w[1];
   if (u < 0.0) {
      if (u < -tol) {
         return 2;
      }

      u = 0.0;
   }

   z = sqrt(u);
   s = z*y/prj->r0;
   if (fabs(s) > 1.0) {
      if (fabs(s) > 1.0+tol) {
         return 2;
      }
      s = copysgn_i (1.0,s);
   }

   xp = 2.0*z*z - 1.0;
   yp = z*x*prj->w[3];
   if (xp == 0.0 && yp == 0.0) {
      *phi = 0.0;
   } else {
      *phi = 2.0*atan2deg (yp, xp);
   }
   *theta = asindeg (s);

   return 0;

}//close aitrev()




int WCSUtils::copset(prjprm * prj)
{
   strcpy(prj->code, "COP");
   prj->flag   = copysgni_i (static_cast<int>(COP), prj->flag);
   prj->phi0   = 0.0;
   prj->theta0 = prj->p[1];

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = sindeg (prj->p[1]);
   if (prj->w[0] == 0.0) {
      return 1;
   }

   prj->w[1] = 1.0/prj->w[0];

   prj->w[3] = prj->r0*cosdeg (prj->p[2]);
   if (prj->w[3] == 0.0) {
      return 1;
   }

   prj->w[4] = 1.0/prj->w[3];
   prj->w[5] = 1.0/tandeg (prj->p[1]);

   prj->w[2] = prj->w[3]*prj->w[5];

   prj->prjfwd = copfwd;
   prj->prjrev = coprev;

   return 0;

}//close copset()


int WCSUtils::copfwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double a, r, s, t;

   if (abs(prj->flag) != COP) {
      if (copset(prj)) return 1;
   }

   t = theta - prj->p[1];
   s = cosdeg (t);
   if (s == 0.0) {
      return 2;
   }

   a = prj->w[0]*phi;
   r = prj->w[2] - prj->w[3]*sindeg (t)/s;

   *x =             r*sindeg (a);
   *y = prj->w[2] - r*cosdeg (a);

   if (prj->flag > 0 && r*prj->w[0] < 0.0) {
      return 2;
   }

   return 0;

}//close copfwd()


int WCSUtils::coprev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double a, dy, r;

   if (abs(prj->flag) != COP) {
      if (copset(prj)) return 1;
   }

   dy = prj->w[2] - y;
   r  = sqrt(x*x + dy*dy);
   if (prj->p[1] < 0.0) r = -r;

   if (r == 0.0) {
      a = 0.0;
   } else {
      a = atan2deg (x/r, dy/r);
   }

   *phi   = a*prj->w[1];
   *theta = prj->p[1] + atandeg (prj->w[5] - r*prj->w[4]);

   return 0;

}//close coprev()



int WCSUtils::coeset(prjprm * prj)
{
   double theta1, theta2;

   strcpy(prj->code, "COE");
   prj->flag   = COE;
   prj->phi0   = 0.0;
   prj->theta0 = prj->p[1];

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   theta1 = prj->p[1] - prj->p[2];
   theta2 = prj->p[1] + prj->p[2];

   prj->w[0] = (sindeg (theta1) + sindeg (theta2))/2.0;
   if (prj->w[0] == 0.0) {
      return 1;
   }

   prj->w[1] = 1.0/prj->w[0];

   prj->w[3] = prj->r0/prj->w[0];
   prj->w[4] = 1.0 + sindeg (theta1)*sindeg (theta2);
   prj->w[5] = 2.0*prj->w[0];
   prj->w[6] = prj->w[3]*prj->w[3]*prj->w[4];
   prj->w[7] = 1.0/(2.0*prj->r0*prj->w[3]);
   prj->w[8] = prj->w[3]*sqrt(prj->w[4] + prj->w[5]);

   prj->w[2] = prj->w[3]*sqrt(prj->w[4] - prj->w[5]*sindeg (prj->p[1]));

   prj->prjfwd = coefwd;
   prj->prjrev = coerev;

   return 0;

}//close coeset()


int WCSUtils::coefwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double a, r;

   if (prj->flag != COE) {
      if (coeset(prj)) return 1;
   }

   a = phi*prj->w[0];
   if (theta == -90.0) {
      r = prj->w[8];
   } else {
      r = prj->w[3]*sqrt(prj->w[4] - prj->w[5]*sindeg (theta));
   }

   *x =             r*sindeg (a);
   *y = prj->w[2] - r*cosdeg (a);

   return 0;

}//close coefwd()


int WCSUtils::coerev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double a, dy, r, w;
   const double tol = 1.0e-12;

   if (prj->flag != COE) {
      if (coeset(prj)) return 1;
   }

   dy = prj->w[2] - y;
   r  = sqrt(x*x + dy*dy);
   if (prj->p[1] < 0.0) r = -r;

   if (r == 0.0) {
      a = 0.0;
   } else {
      a = atan2deg (x/r, dy/r);
   }

   *phi = a*prj->w[1];
   if (fabs(r - prj->w[8]) < tol) {
      *theta = -90.0;
   } else {
      w = (prj->w[6] - r*r)*prj->w[7];
      if (fabs(w) > 1.0) {
         if (fabs(w-1.0) < tol) {
            *theta = 90.0;
         } else if (fabs(w+1.0) < tol) {
            *theta = -90.0;
         } else {
            return 2;
         }
      } else {
         *theta = asindeg (w);
      }
   }

   return 0;

}//close coerev()


int WCSUtils::codset(prjprm * prj)
{
   strcpy(prj->code, "COD");
   prj->flag   = COD;
   prj->phi0   = 0.0;
   prj->theta0 = prj->p[1];

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   if (prj->p[2] == 0.0) {
      prj->w[0] = prj->r0*sindeg (prj->p[1])*kD2R;
   } else {
      prj->w[0] = prj->r0*sindeg (prj->p[1])*sindeg (prj->p[2])/prj->p[2];
   }

   if (prj->w[0] == 0.0) {
      return 1;
   }

   prj->w[1] = 1.0/prj->w[0];
   prj->w[2] = prj->r0*cosdeg (prj->p[2])*cosdeg (prj->p[1])/prj->w[0];
   prj->w[3] = prj->w[2] + prj->p[1];

   prj->prjfwd = codfwd;
   prj->prjrev = codrev;

   return 0;

}//close codset()


int WCSUtils::codfwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double a, r;

   if (prj->flag != COD) {
      if (codset(prj)) return 1;
   }

   a = prj->w[0]*phi;
   r = prj->w[3] - theta;

   *x =             r*sindeg (a);
   *y = prj->w[2] - r*cosdeg (a);

   return 0;

}//close codfwd()


int WCSUtils::codrev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double a, dy, r;

   if (prj->flag != COD) {
      if (codset(prj)) return 1;
   }

   dy = prj->w[2] - y;
   r  = sqrt(x*x + dy*dy);
   if (prj->p[1] < 0.0) r = -r;

   if (r == 0.0) {
      a = 0.0;
   } else {
      a = atan2deg (x/r, dy/r);
   }

   *phi   = a*prj->w[1];
   *theta = prj->w[3] - r;

   return 0;

}//close codrev()


int WCSUtils::cooset(prjprm * prj)
{
   double cos1, cos2, tan1, tan2, theta1, theta2;

   strcpy(prj->code, "COO");
   prj->flag   = COO;
   prj->phi0   = 0.0;
   prj->theta0 = prj->p[1];

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   theta1 = prj->p[1] - prj->p[2];
   theta2 = prj->p[1] + prj->p[2];

   tan1 = tandeg ((90.0 - theta1)/2.0);
   cos1 = cosdeg (theta1);

   if (theta1 == theta2) {
      prj->w[0] = sindeg (theta1);
   } else {
      tan2 = tandeg ((90.0 - theta2)/2.0);
      cos2 = cosdeg (theta2);
      prj->w[0] = log(cos2/cos1)/log(tan2/tan1);
   }
   if (prj->w[0] == 0.0) {
      return 1;
   }

   prj->w[1] = 1.0/prj->w[0];

   prj->w[3] = prj->r0*(cos1/prj->w[0])/pow(tan1,prj->w[0]);
   if (prj->w[3] == 0.0) {
      return 1;
   }
   prj->w[2] = prj->w[3]*pow(tandeg ((90.0 - prj->p[1])/2.0),prj->w[0]);
   prj->w[4] = 1.0/prj->w[3];

   prj->prjfwd = coofwd;
   prj->prjrev = coorev;

   return 0;

}//close cooset()


int WCSUtils::coofwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double a, r;

   if (prj->flag != COO) {
      if (cooset(prj)) return 1;
   }

   a = prj->w[0]*phi;
   if (theta == -90.0) {
      if (prj->w[0] < 0.0) {
         r = 0.0;
      } else {
         return 2;
      }
   } else {
      r = prj->w[3]*pow(tandeg ((90.0 - theta)/2.0),prj->w[0]);
   }

   *x =             r*sindeg (a);
   *y = prj->w[2] - r*cosdeg (a);

   return 0;

}//close coofwd()


int WCSUtils::coorev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double a, dy, r;

   if (prj->flag != COO) {
      if (cooset(prj)) return 1;
   }

   dy = prj->w[2] - y;
   r  = sqrt(x*x + dy*dy);
   if (prj->p[1] < 0.0) r = -r;

   if (r == 0.0) {
      a = 0.0;
   } else {
      a = atan2deg (x/r, dy/r);
   }

   *phi = a*prj->w[1];
   if (r == 0.0) {
      if (prj->w[0] < 0.0) {
         *theta = -90.0;
      } else {
         return 2;
      }
   } else {
      *theta = 90.0 - 2.0*atandeg (pow(r*prj->w[4],prj->w[1]));
   }

   return 0;
}//close coorev()

int WCSUtils::bonset(prjprm * prj)
{
   strcpy(prj->code, "BON");
   prj->flag   = BON;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[1] = 1.0;
      prj->w[2] = prj->r0*cosdeg (prj->p[1])/sindeg (prj->p[1]) + prj->p[1];
   } else {
      prj->w[1] = prj->r0*kD2R;
      prj->w[2] = prj->r0*(cosdeg (prj->p[1])/sindeg (prj->p[1]) + prj->p[1]*kD2R);
   }

   prj->prjfwd = bonfwd;
   prj->prjrev = bonrev;

   return 0;

}//close bonset()


int WCSUtils::bonfwd(const double phi,const double theta,prjprm *prj,double* x,double* y)
{
   double a, r;

   if (prj->p[1] == 0.0) {
      /* Sanson-Flamsteed. */
      return sflfwd(phi, theta, prj, x, y);
   }

   if (prj->flag != BON) {
      if (bonset(prj)) return 1;
   }

   r = prj->w[2] - theta*prj->w[1];
   a = prj->r0*phi*cosdeg (theta)/r;

   *x =             r*sindeg (a);
   *y = prj->w[2] - r*cosdeg (a);

   return 0;

}//close bonfwd()


int WCSUtils::bonrev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double a, cthe, dy, r;

   if (prj->p[1] == 0.0) {
      /* Sanson-Flamsteed. */
      return sflrev(x, y, prj, phi, theta);
   }

   if (prj->flag != BON) {
      if (bonset(prj)) return 1;
   }

   dy = prj->w[2] - y;
   r = sqrt(x*x + dy*dy);
   if (prj->p[1] < 0.0) r = -r;

   if (r == 0.0) {
      a = 0.0;
   } else {
      a = atan2deg (x/r, dy/r);
   }

   *theta = (prj->w[2] - r)/prj->w[1];
   cthe = cosdeg (*theta);
   if (cthe == 0.0) {
      *phi = 0.0;
   } else {
      *phi = a*(r/prj->r0)/cthe;
   }

   return 0;

}//close borev()

int WCSUtils::pcoset(prjprm * prj)
{
   strcpy(prj->code, "PCO");
   prj->flag   = PCO;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
      prj->w[2] = 360.0/kPI;
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = 1.0/prj->w[0];
      prj->w[2] = 2.0*prj->r0;
   }

   prj->prjfwd = pcofwd;
   prj->prjrev = pcorev;

   return 0;

}//close pcoset()


int WCSUtils::pcofwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double a, cthe, cotthe, sthe;

   if (prj->flag != PCO) {
      if (pcoset(prj)) return 1;
   }

   cthe = cosdeg (theta);
   sthe = sindeg (theta);
   a = phi*sthe;

   if (sthe == 0.0) {
      *x = prj->w[0]*phi;
      *y = 0.0;
   } else {
      cotthe = cthe/sthe;
      *x = prj->r0*cotthe*sindeg (a);
      *y = prj->r0*(cotthe*(1.0 - cosdeg (a)) + theta*kD2R);
   }

   return 0;

}//close pcofwd()


int WCSUtils::pcorev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   int   j;
   double f, fneg, fpos, lambda, tanthe, theneg, thepos, w, xp, xx, ymthe, yp;
   const double tol = 1.0e-12;

   if (prj->flag != PCO) {
      if (pcoset(prj)) return 1;
   }

   w = fabs(y*prj->w[1]);
   if (w < tol) {
      *phi = x*prj->w[1];
      *theta = 0.0;
   } else if (fabs(w-90.0) < tol) {
      *phi = 0.0;
      *theta = copysgni_i (90.0,y);
   } else {
      /* Iterative solution using weighted division of the interval. */
      if (y > 0.0) {
         thepos =  90.0;
      } else {
         thepos = -90.0;
      }
      theneg = 0.0;

      xx = x*x;
      ymthe = y - prj->w[0]*thepos;
      fpos = xx + ymthe*ymthe;
      fneg = -999.0;

      for (j = 0; j < 64; j++) {
         if (fneg < -100.0) {
            /* Equal division of the interval. */
            *theta = (thepos+theneg)/2.0;
         } else {
            /* Weighted division of the interval. */
            lambda = fpos/(fpos-fneg);
            if (lambda < 0.1) {
               lambda = 0.1;
            } else if (lambda > 0.9) {
               lambda = 0.9;
            }
            *theta = thepos - lambda*(thepos-theneg);
         }

         /* Compute the residue. */
         ymthe = y - prj->w[0]*(*theta);
         tanthe = tandeg (*theta);
         f = xx + ymthe*(ymthe - prj->w[2]/tanthe);

         /* Check for convergence. */
         if (fabs(f) < tol) break;
         if (fabs(thepos-theneg) < tol) break;

         /* Redefine the interval. */
         if (f > 0.0) {
            thepos = *theta;
            fpos = f;
         } else {
            theneg = *theta;
            fneg = f;
         }
      }

      xp = prj->r0 - ymthe*tanthe;
      yp = x*tanthe;
      if (xp == 0.0 && yp == 0.0) {
         *phi = 0.0;
      } else {
         *phi = atan2deg (yp, xp)/sindeg (*theta);
      }
   }

   return 0;

}//close pcorev()

int WCSUtils::tscset(prjprm * prj)
{
   strcpy(prj->code, "TSC");
   prj->flag   = TSC;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 45.0;
      prj->w[1] = 1.0/45.0;
   } else {
      prj->w[0] = prj->r0*kPI/4.0;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = tscfwd;
   prj->prjrev = tscrev;

   return 0;

}//close tscset()

int WCSUtils::tscfwd(const double phi,const double theta,prjprm *prj,double* x,double* y)
{
   int   face;
   double cthe, l, m, n, rho;
   double x0 = 0.0;
   double y0 = 0.0;
   double xf = 0.0;
   double yf = 0.0;
   const double tol = 1.0e-12;

   if (prj->flag != TSC) {
      if (tscset(prj)) return 1;
   }

   cthe = cosdeg (theta);
   l = cthe*cosdeg (phi);
   m = cthe*sindeg (phi);
   n = sindeg (theta);

   face = 0;
   rho  = n;
   if (l > rho) {
      face = 1;
      rho  = l;
   }
   if (m > rho) {
      face = 2;
      rho  = m;
   }
   if (-l > rho) {
      face = 3;
      rho  = -l;
   }
   if (-m > rho) {
      face = 4;
      rho  = -m;
   }
   if (-n > rho) {
      face = 5;
      rho  = -n;
   }

   if (face == 0) {
      xf =  m/rho;
      yf = -l/rho;
      x0 =  0.0;
      y0 =  2.0;
   } else if (face == 1) {
      xf =  m/rho;
      yf =  n/rho;
      x0 =  0.0;
      y0 =  0.0;
   } else if (face == 2) {
      xf = -l/rho;
      yf =  n/rho;
      x0 =  2.0;
      y0 =  0.0;
   } else if (face == 3) {
      xf = -m/rho;
      yf =  n/rho;
      x0 =  4.0;
      y0 =  0.0;
   } else if (face == 4) {
      xf =  l/rho;
      yf =  n/rho;
      x0 =  6.0;
      y0 =  0.0;
   } else if (face == 5) {
      xf =  m/rho;
      yf =  l/rho;
      x0 =  0.0;
      y0 = -2.0;
   }

   if (fabs(xf) > 1.0) {
      if (fabs(xf) > 1.0+tol) {
         return 2;
      }
      xf = copysgn_i (1.0,xf);
   }
   if (fabs(yf) > 1.0) {
      if (fabs(yf) > 1.0+tol) {
         return 2;
      }
      yf = copysgn_i (1.0,yf);
   }

   *x = prj->w[0]*(xf + x0);
   *y = prj->w[0]*(yf + y0);

   return 0;

}//close tscfwd()


int WCSUtils::tscrev(const double x,const double y,prjprm *prj,double* phi,double* theta)
{
   double l, m, n, xf, yf;

   if (prj->flag != TSC) {
      if (tscset(prj)) return 1;
   }

   xf = x*prj->w[1];
   yf = y*prj->w[1];

   /* Check bounds. */
   if (fabs(xf) <= 1.0) {
      if (fabs(yf) > 3.0) return 2;
   } else {
      if (fabs(xf) > 7.0) return 2;
      if (fabs(yf) > 1.0) return 2;
   }

   /* Map negative faces to the other side. */
   if (xf < -1.0) xf += 8.0;

   /* Determine the face. */
   if (xf > 5.0) {
      /* face = 4 */
      xf = xf - 6.0;
      m  = -1.0/sqrt(1.0 + xf*xf + yf*yf);
      l  = -m*xf;
      n  = -m*yf;
   } else if (xf > 3.0) {
      /* face = 3 */
      xf = xf - 4.0;
      l  = -1.0/sqrt(1.0 + xf*xf + yf*yf);
      m  =  l*xf;
      n  = -l*yf;
   } else if (xf > 1.0) {
      /* face = 2 */
      xf = xf - 2.0;
      m  =  1.0/sqrt(1.0 + xf*xf + yf*yf);
      l  = -m*xf;
      n  =  m*yf;
   } else if (yf > 1.0) {
      /* face = 0 */
      yf = yf - 2.0;
      n  = 1.0/sqrt(1.0 + xf*xf + yf*yf);
      l  = -n*yf;
      m  =  n*xf;
   } else if (yf < -1.0) {
      /* face = 5 */
      yf = yf + 2.0;
      n  = -1.0/sqrt(1.0 + xf*xf + yf*yf);
      l  = -n*yf;
      m  = -n*xf;
   } else {
      /* face = 1 */
      l  =  1.0/sqrt(1.0 + xf*xf + yf*yf);
      m  =  l*xf;
      n  =  l*yf;
   }

   if (l == 0.0 && m == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (m, l);
   }
   *theta = asindeg (n);

   return 0;

}//close tscrev()

int WCSUtils::cscset(prjprm* prj)
{
   strcpy(prj->code, "CSC");
   prj->flag   = CSC;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 45.0;
      prj->w[1] = 1.0/45.0;
   } else {
      prj->w[0] = prj->r0*kPI/4.0;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = cscfwd;
   prj->prjrev = cscrev;

   return 0;

}//close cscset()


int WCSUtils::cscfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   int   face;
   double cthe, eta, l, m, n, rho, xi;
   const float tol = 1.0e-7;

   float a, a2, a2b2, a4, ab, b, b2, b4, ca2, cb2;
   float x0 = 0.0;
   float y0 = 0.0;
   float xf = 0.0;
   float yf = 0.0;
   const float gstar  =  1.37484847732;
   const float mm     =  0.004869491981;
   const float gamma  = -0.13161671474;
   const float omega1 = -0.159596235474;
   const float d0  =  0.0759196200467;
   const float d1  = -0.0217762490699;
   const float c00 =  0.141189631152;
   const float c10 =  0.0809701286525;
   const float c01 = -0.281528535557;
   const float c11 =  0.15384112876;
   const float c20 = -0.178251207466;
   const float c02 =  0.106959469314;

   if (prj->flag != CSC) {
      if (cscset(prj)) return 1;
   }

   cthe = cosdeg (theta);
   l = cthe*cosdeg (phi);
   m = cthe*sindeg (phi);
   n = sindeg (theta);

   face = 0;
   rho  = n;
   if (l > rho) {
      face = 1;
      rho  = l;
   }
   if (m > rho) {
      face = 2;
      rho  = m;
   }
   if (-l > rho) {
      face = 3;
      rho  = -l;
   }
   if (-m > rho) {
      face = 4;
      rho  = -m;
   }
   if (-n > rho) {
      face = 5;
      rho  = -n;
   }

   if (face == 0) {
      xi  =  m;
      eta = -l;
      x0  =  0.0;
      y0  =  2.0;
   } else if (face == 1) {
      xi  =  m;
      eta =  n;
      x0  =  0.0;
      y0  =  0.0;
   } else if (face == 2) {
      xi  = -l;
      eta =  n;
      x0  =  2.0;
      y0  =  0.0;
   } else if (face == 3) {
      xi  = -m;
      eta =  n;
      x0  =  4.0;
      y0  =  0.0;
   } else if (face == 4) {
      xi  =  l;
      eta =  n;
      x0  =  6.0;
      y0  =  0.0;
   } else if (face == 5) {
      xi  =  m;
      eta =  l;
      x0  =  0.0;
      y0  = -2.0;
   }

   a =  xi/rho;
   b = eta/rho;

   a2 = a*a;
   b2 = b*b;
   ca2 = 1.0 - a2;
   cb2 = 1.0 - b2;

   /* Avoid floating underflows. */
   ab   = fabs(a*b);
   a4   = (a2 > 1.0e-16) ? a2*a2 : 0.0;
   b4   = (b2 > 1.0e-16) ? b2*b2 : 0.0;
   a2b2 = (ab > 1.0e-16) ? a2*b2 : 0.0;

   xf = a*(a2 + ca2*(gstar + b2*(gamma*ca2 + mm*a2 +
          cb2*(c00 + c10*a2 + c01*b2 + c11*a2b2 + c20*a4 + c02*b4)) +
          a2*(omega1 - ca2*(d0 + d1*a2))));
   yf = b*(b2 + cb2*(gstar + a2*(gamma*cb2 + mm*b2 +
          ca2*(c00 + c10*b2 + c01*a2 + c11*a2b2 + c20*b4 + c02*a4)) +
          b2*(omega1 - cb2*(d0 + d1*b2))));

   if (fabs(xf) > 1.0) {
      if (fabs(xf) > 1.0+tol) {
         return 2;
      }
      xf = copysgn_i (1.0,xf);
   }
   if (fabs(yf) > 1.0) {
      if (fabs(yf) > 1.0+tol) {
         return 2;
      }
      yf = copysgn_i (1.0,yf);
   }

   *x = prj->w[0]*(x0 + xf);
   *y = prj->w[0]*(y0 + yf);

   return 0;

}//close cscfwd()


int WCSUtils::cscrev(const double x,const double y,prjprm *prj,double* phi,double* theta)
{
   int   face;
   double l = 0.0;
   double m = 0.0;
   double n = 0.0;

   float     a, b, xf, xx, yf, yy, z0, z1, z2, z3, z4, z5, z6;
   const float p00 = -0.27292696;
   const float p10 = -0.07629969;
   const float p20 = -0.22797056;
   const float p30 =  0.54852384;
   const float p40 = -0.62930065;
   const float p50 =  0.25795794;
   const float p60 =  0.02584375;
   const float p01 = -0.02819452;
   const float p11 = -0.01471565;
   const float p21 =  0.48051509;
   const float p31 = -1.74114454;
   const float p41 =  1.71547508;
   const float p51 = -0.53022337;
   const float p02 =  0.27058160;
   const float p12 = -0.56800938;
   const float p22 =  0.30803317;
   const float p32 =  0.98938102;
   const float p42 = -0.83180469;
   const float p03 = -0.60441560;
   const float p13 =  1.50880086;
   const float p23 = -0.93678576;
   const float p33 =  0.08693841;
   const float p04 =  0.93412077;
   const float p14 = -1.41601920;
   const float p24 =  0.33887446;
   const float p05 = -0.63915306;
   const float p15 =  0.52032238;
   const float p06 =  0.14381585;

   if (prj->flag != CSC) {
      if (cscset(prj)) return 1;
   }

   xf = x*prj->w[1];
   yf = y*prj->w[1];

   /* Check bounds. */
   if (fabs(xf) <= 1.0) {
      if (fabs(yf) > 3.0) return 2;
   } else {
      if (fabs(xf) > 7.0) return 2;
      if (fabs(yf) > 1.0) return 2;
   }

   /* Map negative faces to the other side. */
   if (xf < -1.0) xf += 8.0;

   /* Determine the face. */
   if (xf > 5.0) {
      face = 4;
      xf = xf - 6.0;
   } else if (xf > 3.0) {
      face = 3;
      xf = xf - 4.0;
   } else if (xf > 1.0) {
      face = 2;
      xf = xf - 2.0;
   } else if (yf > 1.0) {
      face = 0;
      yf = yf - 2.0;
   } else if (yf < -1.0) {
      face = 5;
      yf = yf + 2.0;
   } else {
      face = 1;
   }

   xx  =  xf*xf;
   yy  =  yf*yf;

   z0 = p00 + xx*(p10 + xx*(p20 + xx*(p30 + xx*(p40 + xx*(p50 + xx*(p60))))));
   z1 = p01 + xx*(p11 + xx*(p21 + xx*(p31 + xx*(p41 + xx*(p51)))));
   z2 = p02 + xx*(p12 + xx*(p22 + xx*(p32 + xx*(p42))));
   z3 = p03 + xx*(p13 + xx*(p23 + xx*(p33)));
   z4 = p04 + xx*(p14 + xx*(p24));
   z5 = p05 + xx*(p15);
   z6 = p06;

   a = z0 + yy*(z1 + yy*(z2 + yy*(z3 + yy*(z4 + yy*(z5 + yy*z6)))));
   a = xf + xf*(1.0 - xx)*a;

   z0 = p00 + yy*(p10 + yy*(p20 + yy*(p30 + yy*(p40 + yy*(p50 + yy*(p60))))));
   z1 = p01 + yy*(p11 + yy*(p21 + yy*(p31 + yy*(p41 + yy*(p51)))));
   z2 = p02 + yy*(p12 + yy*(p22 + yy*(p32 + yy*(p42))));
   z3 = p03 + yy*(p13 + yy*(p23 + yy*(p33)));
   z4 = p04 + yy*(p14 + yy*(p24));
   z5 = p05 + yy*(p15);
   z6 = p06;

   b = z0 + xx*(z1 + xx*(z2 + xx*(z3 + xx*(z4 + xx*(z5 + xx*z6)))));
   b = yf + yf*(1.0 - yy)*b;

   if (face == 0) {
      n =  1.0/sqrt(a*a + b*b + 1.0);
      l = -b*n;
      m =  a*n;
   } else if (face == 1) {
      l =  1.0/sqrt(a*a + b*b + 1.0);
      m =  a*l;
      n =  b*l;
   } else if (face == 2) {
      m =  1.0/sqrt(a*a + b*b + 1.0);
      l = -a*m;
      n =  b*m;
   } else if (face == 3) {
      l = -1.0/sqrt(a*a + b*b + 1.0);
      m =  a*l;
      n = -b*l;
   } else if (face == 4) {
      m = -1.0/sqrt(a*a + b*b + 1.0);
      l = -a*m;
      n = -b*m;
   } else if (face == 5) {
      n = -1.0/sqrt(a*a + b*b + 1.0);
      l = -b*n;
      m = -a*n;
   }

   if (l == 0.0 && m == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (m, l);
   }
   *theta = asindeg (n);

   return 0;
}//close cscrev()


int WCSUtils::qscset(prjprm * prj)
{
   strcpy(prj->code, "QSC");
   prj->flag   = QSC;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 45.0;
      prj->w[1] = 1.0/45.0;
   } else {
      prj->w[0] = prj->r0*kPI/4.0;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = qscfwd;
   prj->prjrev = qscrev;

   return 0;

}//close qscset()


int WCSUtils::qscfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   int   face;
   double cthe, l, m, n, omega, p, rho, rhu, t, tau;
   double xi = 0.0;
   double eta = 0.0;
   double x0 = 0.0;
   double y0 = 0.0;
   double xf = 0.0;
   double yf = 0.0;
   const double tol = 1.0e-12;

   if (prj->flag != QSC) {
      if (qscset(prj)) return 1;
   }

   if (fabs(theta) == 90.0) {
      *x = 0.0;
      *y = copysgn_i (2.0*prj->w[0],theta);
      return 0;
   }

   cthe = cosdeg (theta);
   l = cthe*cosdeg (phi);
   m = cthe*sindeg (phi);
   n = sindeg (theta);

   face = 0;
   rho  = n;
   if (l > rho) {
      face = 1;
      rho  = l;
   }
   if (m > rho) {
      face = 2;
      rho  = m;
   }
   if (-l > rho) {
      face = 3;
      rho  = -l;
   }
   if (-m > rho) {
      face = 4;
      rho  = -m;
   }
   if (-n > rho) {
      face = 5;
      rho  = -n;
   }

   rhu = 1.0 - rho;

   if (face == 0) {
      xi  =  m;
      eta = -l;
      if (rhu < 1.0e-8) {
         /* Small angle formula. */
         t = (90.0 - theta)*kD2R;
         rhu = t*t/2.0;
      }
      x0  =  0.0;
      y0  =  2.0;
   } else if (face == 1) {
      xi  =  m;
      eta =  n;
      if (rhu < 1.0e-8) {
         /* Small angle formula. */
         t = theta*kD2R;
         p = fmod(phi,360.0);
         if (p < -180.0) p += 360.0;
         if (p >  180.0) p -= 360.0;
         p *= kD2R;
         rhu = (p*p + t*t)/2.0;
      }
      x0  =  0.0;
      y0  =  0.0;
   } else if (face == 2) {
      xi  = -l;
      eta =  n;
      if (rhu < 1.0e-8) {
         /* Small angle formula. */
         t = theta*kD2R;
         p = fmod(phi,360.0);
         if (p < -180.0) p += 360.0;
         p = (90.0 - p)*kD2R;
         rhu = (p*p + t*t)/2.0;
      }
      x0  =  2.0;
      y0  =  0.0;
   } else if (face == 3) {
      xi  = -m;
      eta =  n;
      if (rhu < 1.0e-8) {
         /* Small angle formula. */
         t = theta*kD2R;
         p = fmod(phi,360.0);
         if (p < 0.0) p += 360.0;
         p = (180.0 - p)*kD2R;
         rhu = (p*p + t*t)/2.0;
      }
      x0  =  4.0;
      y0  =  0.0;
   } else if (face == 4) {
      xi  =  l;
      eta =  n;
      if (rhu < 1.0e-8) {
         /* Small angle formula. */
         t = theta*kD2R;
         p = fmod(phi,360.0);
         if (p > 180.0) p -= 360.0;
         p *= (90.0 + p)*kD2R;
         rhu = (p*p + t*t)/2.0;
      }
      x0  =  6;
      y0  =  0.0;
   } else if (face == 5) {
      xi  =  m;
      eta =  l;
      if (rhu < 1.0e-8) {
         /* Small angle formula. */
         t = (90.0 + theta)*kD2R;
         rhu = t*t/2.0;
      }
      x0  =  0.0;
      y0  = -2;
   }

   if (xi == 0.0 && eta == 0.0) {
      xf  = 0.0;
      yf  = 0.0;
   } else if (-xi >= fabs(eta)) {
      omega = eta/xi;
      tau = 1.0 + omega*omega;
      xf  = -sqrt(rhu/(1.0-1.0/sqrt(1.0+tau)));
      yf  = (xf/15.0)*(atandeg (omega) - asindeg (omega/sqrt(tau+tau)));
   } else if (xi >= fabs(eta)) {
      omega = eta/xi;
      tau = 1.0 + omega*omega;
      xf  =  sqrt(rhu/(1.0-1.0/sqrt(1.0+tau)));
      yf  = (xf/15.0)*(atandeg (omega) - asindeg (omega/sqrt(tau+tau)));
   } else if (-eta > fabs(xi)) {
      omega = xi/eta;
      tau = 1.0 + omega*omega;
      yf  = -sqrt(rhu/(1.0-1.0/sqrt(1.0+tau)));
      xf  = (yf/15.0)*(atandeg (omega) - asindeg (omega/sqrt(tau+tau)));
   } else if (eta > fabs(xi)) {
      omega = xi/eta;
      tau = 1.0 + omega*omega;
      yf  =  sqrt(rhu/(1.0-1.0/sqrt(1.0+tau)));
      xf  = (yf/15.0)*(atandeg (omega) - asindeg (omega/sqrt(tau+tau)));
   }

   if (fabs(xf) > 1.0) {
      if (fabs(xf) > 1.0+tol) {
         return 2;
      }
      xf = copysgn_i (1.0,xf);
   }
   if (fabs(yf) > 1.0) {
      if (fabs(yf) > 1.0+tol) {
         return 2;
      }
      yf = copysgn_i (1.0,yf);
   }

   *x = prj->w[0]*(xf + x0);
   *y = prj->w[0]*(yf + y0);


   return 0;

}//close qscfwd()


int WCSUtils::qscrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   int   direct, face;
   double omega, rho, rhu, tau, xf, yf, w;
   double l = 0.0;
   double m = 0.0;
   double n = 0.0;
   const double tol = 1.0e-12;

   if (prj->flag != QSC) {
      if (qscset(prj)) return 1;
   }

   xf = x*prj->w[1];
   yf = y*prj->w[1];

   /* Check bounds. */
   if (fabs(xf) <= 1.0) {
      if (fabs(yf) > 3.0) return 2;
   } else {
      if (fabs(xf) > 7.0) return 2;
      if (fabs(yf) > 1.0) return 2;
   }

   /* Map negative faces to the other side. */
   if (xf < -1.0) xf += 8.0;

   /* Determine the face. */
   if (xf > 5.0) {
      face = 4;
      xf = xf - 6.0;
   } else if (xf > 3.0) {
      face = 3;
      xf = xf - 4.0;
   } else if (xf > 1.0) {
      face = 2;
      xf = xf - 2.0;
   } else if (yf > 1.0) {
      face = 0;
      yf = yf - 2.0;
   } else if (yf < -1.0) {
      face = 5;
      yf = yf + 2.0;
   } else {
      face = 1;
   }

   direct = (fabs(xf) > fabs(yf));
   if (direct) {
      if (xf == 0.0) {
         omega = 0.0;
         tau = 1.0;
         rho = 1.0;
         rhu = 0.0;
      } else {
         w = 15.0*yf/xf;
         omega = sindeg (w)/(cosdeg (w) - kSQRT2INV);
         tau = 1.0 + omega*omega;
         rhu = xf*xf*(1.0 - 1.0/sqrt(1.0 + tau));
         rho = 1.0 - rhu;
      }
   } else {
      if (yf == 0.0) {
         omega = 0.0;
         tau = 1.0;
         rho = 1.0;
         rhu = 0.0;
      } else {
         w = 15.0*xf/yf;
         omega = sindeg (w)/(cosdeg (w) - kSQRT2INV);
         tau = 1.0 + omega*omega;
         rhu = yf*yf*(1.0 - 1.0/sqrt(1.0 + tau));
         rho = 1.0 - rhu;
      }
   }

   if (rho < -1.0) {
      if (rho < -1.0-tol) {
         return 2;
      }

      rho = -1.0;
      rhu =  2.0;
      w   =  0.0;
   } else {
      w = sqrt(rhu*(2.0-rhu)/tau);
   }

   if (face == 0) {
      n = rho;
      if (direct) {
         m = w;
         if (xf < 0.0) m = -m;
         l = -m*omega;
      } else {
         l = w;
         if (yf > 0.0) l = -l;
         m = -l*omega;
      }
   } else if (face == 1) {
      l = rho;
      if (direct) {
         m = w;
         if (xf < 0.0) m = -m;
         n = m*omega;
      } else {
         n = w;
         if (yf < 0.0) n = -n;
         m = n*omega;
      }
   } else if (face == 2) {
      m = rho;
      if (direct) {
         l = w;
         if (xf > 0.0) l = -l;
         n = -l*omega;
      } else {
         n = w;
         if (yf < 0.0) n = -n;
         l = -n*omega;
      }
   } else if (face == 3) {
      l = -rho;
      if (direct) {
         m = w;
         if (xf > 0.0) m = -m;
         n = -m*omega;
      } else {
         n = w;
         if (yf < 0.0) n = -n;
         m = -n*omega;
      }
   } else if (face == 4) {
      m = -rho;
      if (direct) {
         l = w;
         if (xf < 0.0) l = -l;
         n = l*omega;
      } else {
         n = w;
         if (yf < 0.0) n = -n;
         l = n*omega;
      }
   } else if (face == 5) {
      n = -rho;
      if (direct) {
         m = w;
         if (xf < 0.0) m = -m;
         l = m*omega;
      } else {
         l = w;
         if (yf < 0.0) l = -l;
         m = l*omega;
      }
   }

   if (l == 0.0 && m == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (m, l);
   }
   *theta = asindeg (n);

   return 0;

}//close qscrev()


int WCSUtils::azpset(prjprm* prj)
{
   strcpy(prj->code, "AZP");
   prj->flag   = copysgni_i (static_cast<int>(AZP), prj->flag);
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = prj->r0*(prj->p[1] + 1.0);
   if (prj->w[0] == 0.0) {
      return 1;
   }

   prj->w[3] = cosdeg (prj->p[2]);
   if (prj->w[3] == 0.0) {
      return 1;
   }

   prj->w[2] = 1.0/prj->w[3];
   prj->w[4] = sindeg (prj->p[2]);
   prj->w[1] = prj->w[4] / prj->w[3];

   if (fabs(prj->p[1]) > 1.0) {
      prj->w[5] = asindeg (-1.0/prj->p[1]);
   } else {
      prj->w[5] = -90.0;
   }

   prj->w[6] = prj->p[1] * prj->w[3];
   prj->w[7] = (fabs(prj->w[6]) < 1.0) ? 1.0 : 0.0;

   prj->prjfwd = azpfwd;
   prj->prjrev = azprev;

	return 0;

}//close azpset()

int WCSUtils::azpfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   double a, b, cphi, cthe, r, s, t;

   if (abs(prj->flag) != AZP) {
      if (azpset(prj)) return 1;
   }

   cphi = cosdeg (phi);
   cthe = cosdeg (theta);

   s = prj->w[1]*cphi;
   t = (prj->p[1] + sindeg (theta)) + cthe*s;
   if (t == 0.0) {
      return 2;
   }

   r  =  prj->w[0]*cthe/t;
   *x =  r*sindeg (phi);
   *y = -r*cphi*prj->w[2];

   /* Bounds checking. */
   if (prj->flag > 0) {
      /* Overlap. */
      if (theta < prj->w[5]) {
         return 2;
      }

      /* Divergence. */
      if (prj->w[7] > 0.0) {
         t = prj->p[1] / sqrt(1.0 + s*s);

         if (fabs(t) <= 1.0) {
            s = atandeg (-s);
            t = asindeg (t);
            a = s - t;
            b = s + t + 180.0;

            if (a > 90.0) a -= 360.0;
            if (b > 90.0) b -= 360.0;

            if (theta < ((a > b) ? a : b)) {
               return 2;
            }
         }
      }
   }

   return 0;

}//close azpfwd()


int WCSUtils::azprev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double a, b, r, s, t, ycosg;
   const double tol = 1.0e-13;

   if (abs(prj->flag) != AZP) {
      if (azpset(prj)) return 1;
   }

   ycosg = y*prj->w[3];

   r = sqrt(x*x + ycosg*ycosg);
   if (r == 0.0) {
      *phi   =  0.0;
      *theta = 90.0;
   } else {
      *phi = atan2deg (x, -ycosg);

      s = r / (prj->w[0] + y*prj->w[4]);
      t = s*prj->p[1]/sqrt(s*s + 1.0);

      s = atan2deg (1.0, s);

      if (fabs(t) > 1.0) {
         t = copysgn_i (90.0,t);
         if (fabs(t) > 1.0+tol) {
            return 2;
         }
      } else {
         t = asindeg (t);
      }

      a = s - t;
      b = s + t + 180.0;

      if (a > 90.0) a -= 360.0;
      if (b > 90.0) b -= 360.0;

      *theta = (a > b) ? a : b;
   }

   return 0;

}//close azprev()


int WCSUtils::szpset(prjprm* prj)
{
   strcpy(prj->code, "SZP");
   prj->flag   = copysgni_i (static_cast<int>(SZP), prj->flag);
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = 1.0/prj->r0;

   prj->w[3] = prj->p[1] * sindeg (prj->p[3]) + 1.0;
   if (prj->w[3] == 0.0) {
      return 1;
   }

   prj->w[1] = -prj->p[1] * cosdeg (prj->p[3]) * sindeg (prj->p[2]);
   prj->w[2] =  prj->p[1] * cosdeg (prj->p[3]) * cosdeg (prj->p[2]);
   prj->w[4] =  prj->r0 * prj->w[1];
   prj->w[5] =  prj->r0 * prj->w[2];
   prj->w[6] =  prj->r0 * prj->w[3];
   prj->w[7] =  (prj->w[3] - 1.0) * prj->w[3] - 1.0;

   if (fabs(prj->w[3] - 1.0) < 1.0) {
      prj->w[8] = asindeg (1.0 - prj->w[3]);
   } else {
      prj->w[8] = -90.0;
   }

   prj->prjfwd = szpfwd;
   prj->prjrev = szprev;

	return 0;

}//close szpset()


int WCSUtils::szpfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   double a, b, cphi, cthe, s, sphi, t;

   if (abs(prj->flag) != SZP) {
      if (szpset(prj)) return 1;
   }

   cphi = cosdeg (phi);
   sphi = sindeg (phi);
   cthe = cosdeg (theta);
   s = 1.0 - sindeg (theta);

   t = prj->w[3] - s;
   if (t == 0.0) {
      return 2;
   }

   *x =  (prj->w[6]*cthe*sphi - prj->w[4]*s)/t;
   *y = -(prj->w[6]*cthe*cphi + prj->w[5]*s)/t;

   /* Bounds checking. */
   if (prj->flag > 0) {
      /* Divergence. */
      if (theta < prj->w[8]) {
         return 2;
      }

      /* Overlap. */
      if (fabs(prj->p[1]) > 1.0) {
         s = prj->w[1]*sphi - prj->w[2]*cphi;
         t = 1.0/sqrt(prj->w[7] + s*s);

         if (fabs(t) <= 1.0) {
            s = atan2deg (s, prj->w[3] - 1.0);
            t = asindeg (t);
            a = s - t;
            b = s + t + 180.0;

            if (a > 90.0) a -= 360.0;
            if (b > 90.0) b -= 360.0;

            if (theta < ((a > b) ? a : b)) {
               return 2;
            }
         }
      }
   }

   return 0;

}//close szpfwd()


int WCSUtils::szprev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double a, b, c, d, r2, sth1, sth2, sthe, sxy, t, x1, xp, y1, yp, z;
   const double tol = 1.0e-13;

   if (abs(prj->flag) != SZP) {
      if (szpset(prj)) return 1;
   }

   xp = x*prj->w[0];
   yp = y*prj->w[0];
   r2 = xp*xp + yp*yp;

   x1 = (xp - prj->w[1])/prj->w[3];
   y1 = (yp - prj->w[2])/prj->w[3];
   sxy = xp*x1 + yp*y1;

   if (r2 < 1.0e-10) {
      /* Use small angle formula. */
      z = r2/2.0;
      *theta = 90.0 - kR2D*sqrt(r2/(1.0 + sxy));

   } else {
      t = x1*x1 + y1*y1;
      a = t + 1.0;
      b = sxy - t;
      c = r2 - sxy - sxy + t - 1.0;
      d = b*b - a*c;

      /* Check for a solution. */
      if (d < 0.0) {
         return 2;
      }
      d = sqrt(d);

      /* Choose solution closest to pole. */
      sth1 = (-b + d)/a;
      sth2 = (-b - d)/a;
      sthe = (sth1 > sth2) ? sth1 : sth2;
      if (sthe > 1.0) {
         if (sthe-1.0 < tol) {
            sthe = 1.0;
         } else {
            sthe = (sth1 < sth2) ? sth1 : sth2;
         }
      }

      if (sthe < -1.0) {
         if (sthe+1.0 > -tol) {
            sthe = -1.0;
         }
      }

      if (sthe > 1.0 || sthe < -1.0) {
         return 2;
      }

      *theta = asindeg (sthe);

      z = 1.0 - sthe;
   }

   *phi = atan2deg (xp - x1*z, -(yp - y1*z));

   return 0;

}//close szprev()

int WCSUtils::tanset(prjprm* prj)
{
  int k;

   strcpy(prj->code, "TAN");
   prj->flag   = copysgni_i (static_cast<int>(TAN), prj->flag);
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->prjfwd = tanfwd;
   prj->prjrev = tanrev;

   for (k = (kMAXPV-1); k >= 0 && prj->ppv[k] == 0.0 && prj->ppv[k+kMAXPV] == 0.0; k--);
   if (k < 0)
      k = 0;
   prj->npv = k;

	return 0;

}//close tanset()


int WCSUtils::tanfwd(const double phi,const double theta,prjprm *prj,double* x,double* y)
{
   double r, s;
   double xp[2];

   if (abs(prj->flag) != TAN) {
      if(tanset(prj)) return 1;
   }

   s = sindeg (theta);
   if (s <= 0.0) {
      return 2;
   }

   r =  prj->r0*cosdeg (theta)/s;
   xp[0] =  r*sindeg (phi);
   xp[1] = -r*cosdeg (phi);
   *x = prj->inv_x? poly_func(prj->inv_x, xp) : xp[0];
   *y = prj->inv_y? poly_func(prj->inv_y, xp) : xp[1];

   if (prj->flag > 0 && s < 0.0) {
      return 2;
   }

   return 0;

}//close tanfwd()


int WCSUtils::tanrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double r;
   double xp;
   double yp;

   if (abs(prj->flag) != TAN) {
      if (tanset(prj)) return 1;
   }

   if (prj->npv) {
      raw_to_pv(prj, x,y, &xp, &yp);
   } else {
      xp = x;
      yp = y;
   }

   r = sqrt(xp*xp + yp*yp);
   if (r == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (xp, -yp);
   }

   *theta = atan2deg (prj->r0, r);

   return 0;

}//close tanrev()


int WCSUtils::raw_to_pv(struct prjprm *prj, double x, double y, double *xo, double *yo)
{
   int k;
   double	*a,*b, r,r3,r5,r7,xy,x2,x3,x4,x5,x6,x7,y2,y3,y4,y5,y6,y7,xp,yp;


   k=prj->npv;
   a = prj->ppv+kMAXPV;		/* Latitude comes first for compatibility */
   b = prj->ppv;			/* Longitude */
   xp = *(a++);
   xp += *(a++)*x;
   yp = *(b++);
   yp += *(b++)*y;
   if (!--k) goto poly_end;
   xp += *(a++)*y;
   yp += *(b++)*x;
   if (!--k) goto poly_end;
   r = sqrt(x*x + y*y);
   xp += *(a++)*r;
   yp += *(b++)*r;
   if (!--k) goto poly_end;
   xp += *(a++)*(x2=x*x);
   yp += *(b++)*(y2=y*y);
   if (!--k) goto poly_end;
   xp += *(a++)*(xy=x*y);
   yp += *(b++)*xy;
   if (!--k) goto poly_end;
   xp += *(a++)*y2;
   yp += *(b++)*x2;
   if (!--k) goto poly_end;
   xp += *(a++)*(x3=x*x2);
   yp += *(b++)*(y3=y*y2);
   if (!--k) goto poly_end;
   xp += *(a++)*x2*y;
   yp += *(b++)*y2*x;
   if (!--k) goto poly_end;
   xp += *(a++)*x*y2;
   yp += *(b++)*y*x2;
   if (!--k) goto poly_end;
   xp += *(a++)*y3;
   yp += *(b++)*x3;
   if (!--k) goto poly_end;
   xp += *(a++)*(r3=r*r*r);
   yp += *(b++)*r3;
   if (!--k) goto poly_end;
   xp += *(a++)*(x4=x2*x2);
   yp += *(b++)*(y4=y2*y2);
   if (!--k) goto poly_end;
   xp += *(a++)*x3*y;
   yp += *(b++)*y3*x;
   if (!--k) goto poly_end;
   xp += *(a++)*x2*y2;
   yp += *(b++)*x2*y2;
   if (!--k) goto poly_end;
   xp += *(a++)*x*y3;
   yp += *(b++)*y*x3;
   if (!--k) goto poly_end;
   xp += *(a++)*y4;
   yp += *(b++)*x4;
   if (!--k) goto poly_end;
   xp += *(a++)*(x5=x4*x);
   yp += *(b++)*(y5=y4*y);
   if (!--k) goto poly_end;
   xp += *(a++)*x4*y;
   yp += *(b++)*y4*x;
   if (!--k) goto poly_end;
   xp += *(a++)*x3*y2;
   yp += *(b++)*y3*x2;
   if (!--k) goto poly_end;
   xp += *(a++)*x2*y3;
   yp += *(b++)*y2*x3;
   if (!--k) goto poly_end;
   xp += *(a++)*x*y4;
   yp += *(b++)*y*x4;
   if (!--k) goto poly_end;
   xp += *(a++)*y5;
   yp += *(b++)*x5;
   if (!--k) goto poly_end;
   xp += *(a++)*(r5=r3*r*r);
   yp += *(b++)*r5;
   if (!--k) goto poly_end;
   xp += *(a++)*(x6=x5*x);
   yp += *(b++)*(y6=y5*y);
   if (!--k) goto poly_end;
   xp += *(a++)*x5*y;
   yp += *(b++)*y5*x;
   if (!--k) goto poly_end;
   xp += *(a++)*x4*y2;
   yp += *(b++)*y4*x2;
   if (!--k) goto poly_end;
   xp += *(a++)*x3*y3;
   yp += *(b++)*y3*x3;
   if (!--k) goto poly_end;
   xp += *(a++)*x2*y4;
   yp += *(b++)*y2*x4;
   if (!--k) goto poly_end;
   xp += *(a++)*x*y5;
   yp += *(b++)*y*x5;
   if (!--k) goto poly_end;
   xp += *(a++)*y6;
   yp += *(b++)*x6;
   if (!--k) goto poly_end;
   xp += *(a++)*(x7=x6*x);
   yp += *(b++)*(y7=y6*y);
   if (!--k) goto poly_end;
   xp += *(a++)*x6*y;
   yp += *(b++)*y6*x;
   if (!--k) goto poly_end;
   xp += *(a++)*x5*y2;
   yp += *(b++)*y5*x2;
   if (!--k) goto poly_end;
   xp += *(a++)*x4*y3;
   yp += *(b++)*y4*x3;
   if (!--k) goto poly_end;
   xp += *(a++)*x3*y4;
   yp += *(b++)*y3*x4;
   if (!--k) goto poly_end;
   xp += *(a++)*x2*y5;
   yp += *(b++)*y2*x5;
   if (!--k) goto poly_end;
   xp += *(a++)*x*y6;
   yp += *(b++)*y*x6;
   if (!--k) goto poly_end;
   xp += *(a++)*y7;
   yp += *(b++)*x7;
   if (!--k) goto poly_end;
   xp += *a*(r7=r5*r*r);
   yp += *b*r7;

poly_end:

  *xo = xp;
  *yo = yp;

   return 0;

}//close raw_to_pv()

int WCSUtils::stgset(prjprm * prj)
{
   strcpy(prj->code, "STG");
   prj->flag   =  STG;
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 360.0/kPI;
      prj->w[1] = kPI/360.0;
   } else {
      prj->w[0] = 2.0*prj->r0;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = stgfwd;
   prj->prjrev = stgrev;

	return 0;

}//close stgset()


int WCSUtils::stgfwd(const double phi,const double theta,prjprm *prj,double* x,double* y)
{
   double r, s;

   if (prj->flag != STG) {
      if (stgset(prj)) return 1;
   }

   s = 1.0 + sindeg (theta);
   if (s == 0.0) {
      return 2;
   }

   r =  prj->w[0]*cosdeg (theta)/s;
   *x =  r*sindeg (phi);
   *y = -r*cosdeg (phi);

   return 0;

}//close stgfwd()


int WCSUtils::stgrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double r;

   if (prj->flag != STG) {
      if (stgset(prj)) return 1;
   }

   r = sqrt(x*x + y*y);
   if (r == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (x, -y);
   }
   *theta = 90.0 - 2.0*atandeg (r*prj->w[1]);

   return 0;

}//close stgrev()



int WCSUtils::sinset(prjprm * prj)
{
   strcpy(prj->code, "SIN");
   prj->flag   = copysgni_i (static_cast<int>(SIN), prj->flag);
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = 1.0/prj->r0;
   prj->w[1] = prj->p[1]*prj->p[1] + prj->p[2]*prj->p[2];
   prj->w[2] = prj->w[1] + 1.0;
   prj->w[3] = prj->w[1] - 1.0;

   prj->prjfwd = sinfwd;
   prj->prjrev = sinrev;

	return 0;

}//close sinset()


int WCSUtils::sinfwd(const double phi,const double theta,prjprm *prj,double* x,double* y)
{
   double cphi, cthe, sphi, t, z;

   if (abs(prj->flag) != SIN) {
      if (sinset(prj)) return 1;
   }

   t = (90.0 - fabs(theta))*kD2R;
   if (t < 1.0e-5) {
      if (theta > 0.0) {
         z = t*t/2.0;
      } else {
         z = 2.0 - t*t/2.0;
      }
      cthe = t;
   } else {
      z =  1.0 - sindeg (theta);
      cthe = cosdeg (theta);
   }

   cphi = cosdeg (phi);
   sphi = sindeg (phi);
   *x =  prj->r0*(cthe*sphi + prj->p[1]*z);
   *y = -prj->r0*(cthe*cphi - prj->p[2]*z);

   /* Validate this solution. */
   if (prj->flag > 0) {
      if (prj->w[1] == 0.0) {
         /* Orthographic projection. */
         if (theta < 0.0) {
            return 2;
         }
      } else {
         /* "Synthesis" projection. */
         t = -atandeg (prj->p[1]*sphi - prj->p[2]*cphi);
         if (theta < t) {
            return 2;
         }
      }
   }

   return 0;

}//close sinfwd()


int WCSUtils::sinrev (const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   const double tol = 1.0e-13;
   double a, b, c, d, r2, sth1, sth2, sthe, sxy, x0, x1, xp, y0, y1, yp, z;

   if (abs(prj->flag) != SIN) {
      if (sinset(prj)) return 1;
   }

   /* Compute intermediaries. */
   x0 = x*prj->w[0];
   y0 = y*prj->w[0];
   r2 = x0*x0 + y0*y0;

   if (prj->w[1] == 0.0) {
      /* Orthographic projection. */
      if (r2 != 0.0) {
         *phi = atan2deg (x0, -y0);
      } else {
         *phi = 0.0;
      }

      if (r2 < 0.5) {
         *theta = acosdeg (sqrt(r2));
      } else if (r2 <= 1.0) {
         *theta = asindeg (sqrt(1.0 - r2));
      } else {
         return 2;
      }

   } else {
      /* "Synthesis" projection. */
      x1 = prj->p[1];
      y1 = prj->p[2];
      sxy = x0*x1 + y0*y1;

      if (r2 < 1.0e-10) {
         /* Use small angle formula. */
         z = r2/2.0;
         *theta = 90.0 - kR2D*sqrt(r2/(1.0 + sxy));

      } else {
         a = prj->w[2];
         b = sxy - prj->w[1];
         c = r2 - sxy - sxy + prj->w[3];
         d = b*b - a*c;

         /* Check for a solution. */
         if (d < 0.0) {
            return 2;
         }
         d = sqrt(d);

         /* Choose solution closest to pole. */
         sth1 = (-b + d)/a;
         sth2 = (-b - d)/a;
         sthe = (sth1 > sth2) ? sth1 : sth2;
         if (sthe > 1.0) {
            if (sthe-1.0 < tol) {
               sthe = 1.0;
            } else {
               sthe = (sth1 < sth2) ? sth1 : sth2;
            }
         }

         if (sthe < -1.0) {
            if (sthe+1.0 > -tol) {
               sthe = -1.0;
            }
         }

         if (sthe > 1.0 || sthe < -1.0) {
            return 2;
         }

         *theta = asindeg (sthe);
         z = 1.0 - sthe;
      }

      xp = -y0 + prj->p[2]*z;
      yp =  x0 - prj->p[1]*z;
      if (xp == 0.0 && yp == 0.0) {
         *phi = 0.0;
      } else {
         *phi = atan2deg (yp,xp);
      }
   }

   return 0;

}//close sinrev()


int WCSUtils::arcset(prjprm *prj)
{
   strcpy(prj->code, "ARC");
   prj->flag   =  ARC;
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = arcfwd;
   prj->prjrev = arcrev;

	return 0;

}//close arcset()

int WCSUtils::arcfwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   double r;

   if (prj->flag != ARC) {
      if (arcset(prj)) return 1;
   }

   r =  prj->w[0]*(90.0 - theta);
   *x =  r*sindeg (phi);
   *y = -r*cosdeg (phi);

   return 0;

}//close arcfwd()


int WCSUtils::arcrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double r;

   if (prj->flag != ARC) {
      if (arcset(prj)) return 1;
   }

   r = sqrt(x*x + y*y);
   if (r == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (x, -y);
   }
   *theta = 90.0 - r*prj->w[1];

   return 0;

}//close arcrev()

int WCSUtils::zeaset(prjprm* prj)
{
   strcpy(prj->code, "ZEA");
   prj->flag   =  ZEA;
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 360.0/kPI;
      prj->w[1] = kPI/360.0;
   } else {
      prj->w[0] = 2.0*prj->r0;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = zeafwd;
   prj->prjrev = zearev;

   return 0;

}//close zeaset()

int WCSUtils::zeafwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   double r;

   if (prj->flag != ZEA) {
      if (zeaset(prj)) return 1;
   }

   r =  prj->w[0]*sindeg ((90.0 - theta)/2.0);
   *x =  r*sindeg (phi);
   *y = -r*cosdeg (phi);

   return 0;

}//close zeafwd()


int WCSUtils::zearev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double r, s;
   const double tol = 1.0e-12;

   if (prj->flag != ZEA) {
      if (zeaset(prj)) return 1;
   }

   r = sqrt(x*x + y*y);
   if (r == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (x, -y);
   }

   s = r*prj->w[1];
   if (fabs(s) > 1.0) {
      if (fabs(r - prj->w[0]) < tol) {
         *theta = -90.0;
      } else {
         return 2;
      }
   } else {
      *theta = 90.0 - 2.0*asindeg (s);
   }

   return 0;

}//close zearev()

int WCSUtils::airset(prjprm * prj)
{
   const double tol = 1.0e-4;
   double cxi;

   strcpy(prj->code, "AIR");
   prj->flag   =  AIR;
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   prj->w[0] = 2.0*prj->r0;
   if (prj->p[1] == 90.0) {
      prj->w[1] = -0.5;
      prj->w[2] =  1.0;
   } else if (prj->p[1] > -90.0) {
      cxi = cosdeg ((90.0 - prj->p[1])/2.0);
      prj->w[1] = log(cxi)*(cxi*cxi)/(1.0-cxi*cxi);
      prj->w[2] = 0.5 - prj->w[1];
   } else {
      return 1;
   }

   prj->w[3] = prj->w[0] * prj->w[2];
   prj->w[4] = tol;
   prj->w[5] = prj->w[2]*tol;
   prj->w[6] = kR2D/prj->w[2];

   prj->prjfwd = airfwd;
   prj->prjrev = airrev;

   return 0;

}//close airset()


int WCSUtils::airfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   double cxi, r, txi, xi;

   if (prj->flag != AIR) {
      if (airset(prj)) return 1;
   }

   if (theta == 90.0) {
      r = 0.0;
   } else if (theta > -90.0) {
      xi = kD2R*(90.0 - theta)/2.0;
      if (xi < prj->w[4]) {
         r = xi*prj->w[3];
      } else {
         cxi = cosdeg ((90.0 - theta)/2.0);
         txi = sqrt(1.0-cxi*cxi)/cxi;
         r = -prj->w[0]*(log(cxi)/txi + prj->w[1]*txi);
      }
   } else {
      return 2;
   }

   *x =  r*sindeg (phi);
   *y = -r*cosdeg (phi);

   return 0;

}//close airfwd()


int WCSUtils::airrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   int   j;
   double cxi, lambda, r, r1, r2, rt, txi, x1, x2, xi;
   const double tol = 1.0e-12;

   if (prj->flag != AIR) {
      if (airset(prj)) return 1;
   }

   r = sqrt(x*x + y*y)/prj->w[0];

   if (r == 0.0) {
      xi = 0.0;
   } else if (r < prj->w[5]) {
      xi = r*prj->w[6];
   } else {
      /* Find a solution interval. */
      x1 = 1.0;
      r1 = 0.0;
      for (j = 0; j < 30; j++) {
         x2 = x1/2.0;
         txi = sqrt(1.0-x2*x2)/x2;
         r2 = -(log(x2)/txi + prj->w[1]*txi);

         if (r2 >= r) break;
         x1 = x2;
         r1 = r2;
      }
      if (j == 30) return 2;

      for (j = 0; j < 100; j++) {
         /* Weighted division of the interval. */
         lambda = (r2-r)/(r2-r1);
         if (lambda < 0.1) {
            lambda = 0.1;
         } else if (lambda > 0.9) {
            lambda = 0.9;
         }
         cxi = x2 - lambda*(x2-x1);

         txi = sqrt(1.0-cxi*cxi)/cxi;
         rt = -(log(cxi)/txi + prj->w[1]*txi);

         if (rt < r) {
             if (r-rt < tol) break;
             r1 = rt;
             x1 = cxi;
         } else {
             if (rt-r < tol) break;
             r2 = rt;
             x2 = cxi;
         }
      }
      if (j == 100) return 2;

      xi = acosdeg (cxi);
   }

   if (r == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (x, -y);
   }
   *theta = 90.0 - 2.0*xi;

   return 0;

}//close airrev()

int WCSUtils::ceaset(prjprm * prj)
{
   strcpy(prj->code, "CEA");
   prj->flag   = CEA;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
      if (prj->p[1] <= 0.0 || prj->p[1] > 1.0) {
         return 1;
      }
      prj->w[2] = prj->r0/prj->p[1];
      prj->w[3] = prj->p[1]/prj->r0;
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = kR2D/prj->r0;
      if (prj->p[1] <= 0.0 || prj->p[1] > 1.0) {
         return 1;
      }
      prj->w[2] = prj->r0/prj->p[1];
      prj->w[3] = prj->p[1]/prj->r0;
   }

   prj->prjfwd = ceafwd;
   prj->prjrev = cearev;

   return 0;

}//close ceaset()

int WCSUtils::ceafwd(const double phi,const double theta,prjprm *prj,double* x,double* y)
{
   if (prj->flag != CEA) {
      if (ceaset(prj)) return 1;
   }

   *x = prj->w[0]*phi;
   *y = prj->w[2]*sindeg (theta);

   return 0;

}//close ceafwd()


int WCSUtils::cearev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double s;
   const double tol = 1.0e-13;

   if (prj->flag != CEA) {
      if (ceaset(prj)) return 1;
   }

   s = y*prj->w[3];
   if (fabs(s) > 1.0) {
      if (fabs(s) > 1.0+tol) {
         return 2;
      }
      s = copysgn_i (1.0,s);
   }

   *phi   = x*prj->w[1];
   *theta = asindeg (s);

   return 0;

}//close cearev()


int WCSUtils::carset(prjprm * prj)
{
   strcpy(prj->code, "CAR");
   prj->flag   = CAR;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = carfwd;
   prj->prjrev = carrev;

   return 0;

}//close carset()

int WCSUtils::carfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   if (prj->flag != CAR) {
      if (carset(prj)) return 1;
   }

   *x = prj->w[0]*phi;
   *y = prj->w[0]*theta;

   return 0;

}//close carfwd()


int WCSUtils::carrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   if (prj->flag != CAR) {
      if (carset(prj)) return 1;
   }

   *phi   = prj->w[1]*x;
   *theta = prj->w[1]*y;

   return 0;

}//close carrev()

int WCSUtils::merset(prjprm * prj)
{
   strcpy(prj->code, "MER");
   prj->flag   = MER;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = merfwd;
   prj->prjrev = merrev;

   return 0;

}//close merset()

int WCSUtils::merfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   if (prj->flag != MER) {
      if (merset(prj)) return 1;
   }

   if (theta <= -90.0 || theta >= 90.0) {
      return 2;
   }

   *x = prj->w[0]*phi;
   *y = prj->r0*log(tandeg ((90.0+theta)/2.0));

   return 0;

}//close merfwd()


int WCSUtils::merrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   if (prj->flag != MER) {
      if (merset(prj)) return 1;
   }

   *phi   = x*prj->w[1];
   *theta = 2.0*atandeg (exp(y/prj->r0)) - 90.0;

   return 0;

}//close merrev()

int WCSUtils::sflset(prjprm * prj)
{
   strcpy(prj->code, "SFL");
   prj->flag   = SFL;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;
      prj->w[0] = 1.0;
      prj->w[1] = 1.0;
   } else {
      prj->w[0] = prj->r0*kD2R;
      prj->w[1] = 1.0/prj->w[0];
   }

   prj->prjfwd = sflfwd;
   prj->prjrev = sflrev;

   return 0;

}//close sflset()


int WCSUtils::sflfwd(const double phi,const double theta,prjprm* prj,double* x,double* y)
{
   if (prj->flag != SFL) {
      if (sflset(prj)) return 1;
   }

   *x = prj->w[0]*phi*cosdeg (theta);
   *y = prj->w[0]*theta;

   return 0;

}//close sflfwd()


int WCSUtils::sflrev(const double x,const double y,prjprm* prj,double* phi,double* theta)
{
   double w;

   if (prj->flag != SFL) {
      if (sflset(prj)) return 1;
   }

   w = cos(y/prj->r0);
   if (w == 0.0) {
      *phi = 0.0;
   } else {
      *phi = x*prj->w[1]/cos(y/prj->r0);
   }
   *theta = y*prj->w[1];

   return 0;

}//close sflrev()

int WCSUtils::cypset(prjprm * prj)
{
   strcpy(prj->code, "CYP");
   prj->flag   = CYP;
   prj->phi0   = 0.0;
   prj->theta0 = 0.0;

   if (prj->r0 == 0.0) {
      prj->r0 = kR2D;

      prj->w[0] = prj->p[2];
      if (prj->w[0] == 0.0) {
         return 1;
      }

      prj->w[1] = 1.0/prj->w[0];

      prj->w[2] = kR2D*(prj->p[1] + prj->p[2]);
      if (prj->w[2] == 0.0) {
         return 1;
      }

      prj->w[3] = 1.0/prj->w[2];
   } else {
      prj->w[0] = prj->r0*prj->p[2]*kD2R;
      if (prj->w[0] == 0.0) {
         return 1;
      }

      prj->w[1] = 1.0/prj->w[0];

      prj->w[2] = prj->r0*(prj->p[1] + prj->p[2]);
      if (prj->w[2] == 0.0) {
         return 1;
      }

      prj->w[3] = 1.0/prj->w[2];
   }

   prj->prjfwd = cypfwd;
   prj->prjrev = cyprev;

   return 0;

}//close cypset()


int WCSUtils::cypfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   double s;

   if (prj->flag != CYP) {
      if (cypset(prj)) return 1;
   }

   s = prj->p[1] + cosdeg (theta);
   if (s == 0.0) {
         return 2;
      }

   *x = prj->w[0]*phi;
   *y = prj->w[2]*sindeg (theta)/s;

   return 0;

}//close cypfwd()


int WCSUtils::cyprev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   double eta;

   if (prj->flag != CYP) {
      if (cypset(prj)) return 1;
   }

   *phi   = x*prj->w[1];
   eta    = y*prj->w[3];
   *theta = atan2deg (eta,1.0) + asindeg (eta*prj->p[1]/sqrt(eta*eta+1.0));

   return 0;

}//close cyprev()


int WCSUtils::zpnset(prjprm* prj)
{
   int   i, j, k;
   double d, d1, d2, r, zd, zd1, zd2;
   const double tol = 1.0e-13;

   strcpy(prj->code, "ZPN");
   prj->flag   = copysgni_i (static_cast<int>(ZPN), prj->flag);
   prj->phi0   =  0.0;
   prj->theta0 = 90.0;

   if (prj->r0 == 0.0) prj->r0 = kR2D;

   /* Find the highest non-zero coefficient. */
   for (k = 9; k >= 0 && prj->p[k] == 0.0; k--){
	i = 0; }
   /* if (k < 0) return 1; */

   /* if (k < 0 switch to ARC projection */
   if (k < 0) {
	return (arcset (prj));
	}

   prj->n = k;

   /* No negative derivative -> no point of inflection. */
   zd = kPI;

   /* Processing subroutines */
   prj->prjfwd = zpnfwd;
   prj->prjrev = zpnrev;

   if (k >= 3) {
      /* Find the point of inflection closest to the pole. */
      zd1 = 0.0;
      d1  = prj->p[1];
      if (d1 <= 0.0) {
         return 1;
      }

      /* Find the point where the derivative first goes negative. */
      for (i = 0; i < 180; i++) {
         zd2 = i*kD2R;
         d2  = 0.0;
         for (j = k; j > 0; j--) {
            d2 = d2*zd2 + j*prj->p[j];
         }

         if (d2 <= 0.0) break;
         zd1 = zd2;
         d1  = d2;
      }

      if (i == 180) {
         /* No negative derivative -> no point of inflection. */
         zd = kPI;
      } else {
         /* Find where the derivative is zero. */
         for (i = 1; i <= 10; i++) {
            zd = zd1 - d1*(zd2-zd1)/(d2-d1);

            d = 0.0;
            for (j = k; j > 0; j--) {
               d = d*zd + j*prj->p[j];
            }

            if (fabs(d) < tol) break;

            if (d < 0.0) {
               zd2 = zd;
               d2  = d;
            } else {
               zd1 = zd;
               d1  = d;
            }
         }
      }

      r = 0.0;
      for (j = k; j >= 0; j--) {
         r = r*zd + prj->p[j];
      }
      prj->w[0] = zd;
      prj->w[1] = r;
   }

	return 0;

}//close zpnset()

int WCSUtils::zpnfwd(const double phi,const double theta,prjprm *prj,double *x,double *y)
{
   int   j;
   double r, s;

   if (abs(prj->flag) != ZPN) {
      if (zpnset(prj)) return 1;
   }

   s = (90.0 - theta)*kD2R;

   r = 0.0;
   for (j = 9; j >= 0; j--) {
      r = r*s + prj->p[j];
   }
   r = prj->r0*r;

   *x =  r*sindeg (phi);
   *y = -r*cosdeg (phi);

   if (prj->flag > 0 && s > prj->w[0]) {
      return 2;
   }

   return 0;

}//close zpnfwd()


int WCSUtils::zpnrev(const double x,const double y,prjprm *prj,double *phi,double *theta)
{
   int   i, j, k;
   double a, b, c, d, lambda, r, r1, r2, rt, zd, zd1, zd2;
   const double tol = 1.0e-13;

   if (abs(prj->flag) != ZPN) {
      if (zpnset(prj)) return 1;
   }

   k = prj->n;

   r = sqrt(x*x + y*y)/prj->r0;

   if (k < 1) {
      /* Constant - no solution. */
      return 1;
   } else if (k == 1) {
      /* Linear. */
      zd = (r - prj->p[0])/prj->p[1];
   } else if (k == 2) {
      /* Quadratic. */
      a = prj->p[2];
      b = prj->p[1];
      c = prj->p[0] - r;

      d = b*b - 4.0*a*c;
      if (d < 0.0) {
         return 2;
      }
      d = sqrt(d);

      /* Choose solution closest to pole. */
      zd1 = (-b + d)/(2.0*a);
      zd2 = (-b - d)/(2.0*a);
      zd  = (zd1<zd2) ? zd1 : zd2;
      if (zd < -tol) zd = (zd1>zd2) ? zd1 : zd2;
      if (zd < 0.0) {
         if (zd < -tol) {
            return 2;
         }
         zd = 0.0;
      } else if (zd > kPI) {
         if (zd > kPI+tol) {
            return 2;
         }
         zd = kPI;
      }
   } else {
      /* Higher order - solve iteratively. */
      zd1 = 0.0;
      r1  = prj->p[0];
      zd2 = prj->w[0];
      r2  = prj->w[1];

      if (r < r1) {
         if (r < r1-tol) {
            return 2;
         }
         zd = zd1;
      } else if (r > r2) {
         if (r > r2+tol) {
            return 2;
         }
         zd = zd2;
      } else {
         /* Disect the interval. */
         for (j = 0; j < 100; j++) {
            lambda = (r2 - r)/(r2 - r1);
            if (lambda < 0.1) {
               lambda = 0.1;
            } else if (lambda > 0.9) {
               lambda = 0.9;
            }

            zd = zd2 - lambda*(zd2 - zd1);

            rt = 0.0;
            for (i = k; i >= 0; i--) {
                rt = (rt * zd) + prj->p[i];
            }

            if (rt < r) {
                if (r-rt < tol) break;
                r1 = rt;
                zd1 = zd;
            } else {
                if (rt-r < tol) break;
                r2 = rt;
                zd2 = zd;
            }

            if (fabs(zd2-zd1) < tol) break;
         }
      }
   }

   if (r == 0.0) {
      *phi = 0.0;
   } else {
      *phi = atan2deg (x, -y);
   }
   *theta = 90.0 - zd*kR2D;

   return 0;

}//close zpnrev()


int WCSUtils::dsspos(double xpix,double ypix,WCS* wcs,double* xpos,double* ypos)
{
  double x, y, xmm, ymm, xmm2, ymm2, xmm3, ymm3, x2y2;
  double xi, xir, eta, etar, raoff, ra, dec;
  double cond2r = 1.745329252e-2;
  double cons2r = 206264.8062470964;
  double twopi = 6.28318530717959;
  double ctan, ccos;

	/*  
	Ignore magnitude and color terms 
  double mag = 0.0;
  double color = 0.0; 
	*/

	// Convert from image pixels to plate pixels
  x = xpix + wcs->x_pixel_offset - 1.0 + 0.5;
  y = ypix + wcs->y_pixel_offset - 1.0 + 0.5;

	// Convert from pixels to millimeters
  xmm = (wcs->ppo_coeff[2] - x * wcs->x_pixel_size) / 1000.0;
  ymm = (y * wcs->y_pixel_size - wcs->ppo_coeff[5]) / 1000.0;
  xmm2 = xmm * xmm;
  ymm2 = ymm * ymm;
  xmm3 = xmm * xmm2;
  ymm3 = ymm * ymm2;
  x2y2 = xmm2 + ymm2;

	//Compute coordinates from x,y and plate model
  xi =  wcs->x_coeff[ 0]*xmm	+ wcs->x_coeff[ 1]*ymm +
	wcs->x_coeff[ 2]		+ wcs->x_coeff[ 3]*xmm2 +
	wcs->x_coeff[ 4]*xmm*ymm	+ wcs->x_coeff[ 5]*ymm2 +
	wcs->x_coeff[ 6]*(x2y2)	+ wcs->x_coeff[ 7]*xmm3 +
	wcs->x_coeff[ 8]*xmm2*ymm	+ wcs->x_coeff[ 9]*xmm*ymm2 +
	wcs->x_coeff[10]*ymm3	+ wcs->x_coeff[11]*xmm*(x2y2) +
	wcs->x_coeff[12]*xmm*x2y2*x2y2;

	/*  
	Ignore magnitude and color terms 
	+ wcs->x_coeff[13]*mag	+ wcs->x_coeff[14]*mag*mag +
	wcs->x_coeff[15]*mag*mag*mag + wcs->x_coeff[16]*mag*xmm +
	wcs->x_coeff[17]*mag*x2y2	+ wcs->x_coeff[18]*mag*xmm*x2y2 +
	wcs->x_coeff[19]*color; 
	*/

  eta =	wcs->y_coeff[ 0]*ymm	+ wcs->y_coeff[ 1]*xmm +
	wcs->y_coeff[ 2]		+ wcs->y_coeff[ 3]*ymm2 +
	wcs->y_coeff[ 4]*xmm*ymm	+ wcs->y_coeff[ 5]*xmm2 +
	wcs->y_coeff[ 6]*(x2y2)	+ wcs->y_coeff[ 7]*ymm3 +
	wcs->y_coeff[ 8]*ymm2*xmm	+ wcs->y_coeff[ 9]*ymm*xmm2 +
	wcs->y_coeff[10]*xmm3	+ wcs->y_coeff[11]*ymm*(x2y2) +
	wcs->y_coeff[12]*ymm*x2y2*x2y2;

	/*  
	Ignore magnitude and color terms 
	+ wcs->y_coeff[13]*mag	+ wcs->y_coeff[14]*mag*mag +
	wcs->y_coeff[15]*mag*mag*mag + wcs->y_coeff[16]*mag*ymm +
	wcs->y_coeff[17]*mag*x2y2)	+ wcs->y_coeff[18]*mag*ymm*x2y2 +
	wcs->y_coeff[19]*color; 
	*/

	// Convert to radians
  xir = xi / cons2r;
  etar = eta / cons2r;

	//Convert to RA and Dec
  ctan = tan (wcs->plate_dec);
  ccos = cos (wcs->plate_dec);
  raoff = atan2 (xir / ccos, 1.0 - etar * ctan);
  ra = raoff + wcs->plate_ra;
  if (ra < 0.0) ra = ra + twopi;
  *xpos = ra / cond2r;

  dec = atan (cos (raoff) * ((etar + ctan) / (1.0 - (etar * ctan))));
  *ypos = dec / cond2r;

  return 0;

}//close dsspos()


int WCSUtils::platepos (double xpix,double ypix,WCS* wcs,double* xpos,double* ypos)
{
    double x, y, x2, y2, x3, y3, r2;
    double xi, xir, eta, etar, raoff, ra, dec, ra0, dec0;
    double twopi = 6.28318530717959;
    double ctan, ccos;
    int ncoeff1 = wcs->ncoeff1;
    int ncoeff2 = wcs->ncoeff2;

    /*  
		//Ignore magnitude and color terms 
    double mag = 0.0;
    double color = 0.0; 
		*/

    //Convert from pixels to millimeters
    x = xpix - wcs->crpix[0];
    y = ypix - wcs->crpix[1];
    x2 = x * x;
    y2 = y * y;
    x3 = x * x2;
    y3 = y * y2;
    r2 = x2 + y2;

    //Compute xi,eta coordinates in degrees from x,y and plate model
    xi =  wcs->x_coeff[ 0]		+ wcs->x_coeff[ 1]*x +
	  wcs->x_coeff[ 2]*y	+ wcs->x_coeff[ 3]*x2 +
	  wcs->x_coeff[ 4]*y2	+ wcs->x_coeff[ 5]*x*y;

    if (ncoeff1 > 6)
	  xi = xi + wcs->x_coeff[ 6]*x3	+ wcs->x_coeff[ 7]*y3;

    if (ncoeff1 > 8) {
	xi = xi + wcs->x_coeff[ 8]*x2*y	+ wcs->x_coeff[ 9]*x*y2 +
		  wcs->x_coeff[10]*(r2)	+ wcs->x_coeff[11]*x*r2 +
		  wcs->x_coeff[12]*y*r2;
	}

    eta = wcs->y_coeff[ 0]		+ wcs->y_coeff[ 1]*x +
	  wcs->y_coeff[ 2]*y	+ wcs->y_coeff[ 3]*x2 +
	  wcs->y_coeff[ 4]*y2	+ wcs->y_coeff[ 5]*x*y;

    if (ncoeff2 > 6)
	eta = eta + wcs->y_coeff[ 6]*x3	+ wcs->y_coeff[ 7]*y3;

    if (ncoeff2 > 8) {
	eta = eta + wcs->y_coeff[ 8]*x2*y + wcs->y_coeff[ 9]*y2*x +
		    wcs->y_coeff[10]*r2   + wcs->y_coeff[11]*x*r2 +
		    wcs->y_coeff[12]*y*r2;
	}

    //Convert to radians
    xir = degrad_i (xi);
    etar = degrad_i (eta);

    //Convert to RA and Dec
    ra0 = degrad_i (wcs->crval[0]);
    dec0 = degrad_i (wcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    raoff = atan2 (xir / ccos, 1.0 - etar * ctan);
    ra = raoff + ra0;
    if (ra < 0.0) ra = ra + twopi;
    *xpos = raddeg_i (ra);

    dec = atan (cos (raoff) / ((1.0 - (etar * ctan)) / (etar + ctan)));
    *ypos = raddeg_i (dec);
   
 	return 0;

}//close platepos()


int WCSUtils::tnxpos(double xpix,double ypix,WCS* wcs,double* xpos,double* ypos)
{
    int	ira, idec;
    double x, y, r, phi, theta, costhe, sinthe, dphi, cosphi, sinphi, dlng, z;
    double colatp, coslatp, sinlatp, longp;
    double xs, ys, ra, dec, xp, yp;
    //double wf_gseval();

    /* Convert from pixels to image coordinates */
    xpix = xpix - wcs->crpix[0];
    ypix = ypix - wcs->crpix[1];

    /* Scale and rotate using CD matrix */
    if (wcs->rotmat) {
	x = xpix * wcs->cd[0] + ypix * wcs->cd[1];
	y = xpix * wcs->cd[2] + ypix * wcs->cd[3];
	}

    else {

	/* Check axis increments - bail out if either 0 */
	if (wcs->cdelt[0] == 0.0 || wcs->cdelt[1] == 0.0) {
	    *xpos = 0.0;
	    *ypos = 0.0;
	    return 2;
	    }

	/* Scale using CDELT */
	xs = xpix * wcs->cdelt[0];
	ys = ypix * wcs->cdelt[1];

	/* Take out rotation from CROTA */
	if (wcs->rot != 0.0) {
	    double cosr = cos (degrad_i (wcs->rot));
	    double sinr = sin (degrad_i (wcs->rot));
	    x = xs * cosr - ys * sinr;
	    y = xs * sinr + ys * cosr;
    	    }
	else {
	    x = xs;
	    y = ys;
	    }
	}

    /* get the axis numbers */
    if (wcs->coorflip) {
	ira = 1;
	idec = 0;
	}
    else {
	ira = 0;
	idec = 1;
	}
    colatp = degrad_i (90.0 - wcs->crval[idec]);
    coslatp = cos(colatp);
    sinlatp = sin(colatp);
    longp = degrad_i(wcs->longpole);

    /*  Compute native spherical coordinates phi and theta in degrees from the
	projected coordinates. this is the projection part of the computation */
    if (wcs->lngcor != NULL)
	xp = x + wf_gseval (wcs->lngcor, x, y);
    else
	xp = x;
    if (wcs->latcor != NULL)
	yp = y + wf_gseval (wcs->latcor, x, y);
    else
	yp = y;
    x = xp;
    y = yp;
    r = sqrt (x * x + y * y);

    /* Compute phi */
    if (r == 0.0)
	phi = 0.0;
    else
	phi = atan2 (x, -y);

    /* Compute theta */
    theta = atan2 (wcs->rodeg, r);

    /*  Compute the celestial coordinates ra and dec from the native
	coordinates phi and theta. this is the spherical geometry part
	of the computation */

    costhe = cos (theta);
    sinthe = sin (theta);
    dphi = phi - longp;
    cosphi = cos (dphi);
    sinphi = sin (dphi);

    /* Compute the ra */
    x = sinthe * sinlatp - costhe * coslatp * cosphi;
    if (fabs (x) < kSPHTOL)
	x = -cos (theta + colatp) + costhe * coslatp * (1.0 - cosphi);
    y = -costhe * sinphi;
    if (x != 0.0 || y != 0.0)
	dlng = atan2 (y, x);
    else
	dlng = dphi + kPI ;
    ra =  wcs->crval[ira] + raddeg_i(dlng);

    /* normalize ra */
    if (wcs->crval[ira] >= 0.0) {
	if (ra < 0.0)
	    ra = ra + 360.0;
	}
    else {
	if (ra > 0.0)
	    ra = ra - 360.0;
	}
    if (ra > 360.0)
	ra = ra - 360.0;
    else if (ra < -360.0)
	ra = ra + 360.0;

    /* compute the dec */
    if (fmod (dphi, kPI) == 0.0) {
	dec = raddeg_i(theta + cosphi * colatp);
	if (dec > 90.0)
	    dec = 180.0 - dec;
	if (dec < -90.0)
	    dec = -180.0 - dec;
	}
    else {
	z = sinthe * coslatp + costhe * sinlatp * cosphi;
	if (fabs(z) > 0.99) {
	    if (z >= 0.0)
		dec = raddeg_i(acos (sqrt(x * x + y * y)));
	    else
		dec = raddeg_i(-acos (sqrt(x * x + y * y)));
	    }
	else
		dec = raddeg_i(asin (z));
	}

    /* store the results */
    *xpos  = ra;
    *ypos = dec;
    
	return (0);

}//close tnxpos()


double WCSUtils::wf_gseval(IRAFsurface* sf,double x,double y)
{
    double sum, accum;
    int i, ii, k, maxorder, xorder;

    /* Calculate the basis functions */
    switch (sf->type) {
        case eTNX_CHEBYSHEV:
            wf_gsb1cheb (x, sf->xorder, sf->xmaxmin, sf->xrange, sf->xbasis);
            wf_gsb1cheb (y, sf->yorder, sf->ymaxmin, sf->yrange, sf->ybasis);
	    break;
        case eTNX_LEGENDRE:
            wf_gsb1leg (x, sf->xorder, sf->xmaxmin, sf->xrange, sf->xbasis);
            wf_gsb1leg (y, sf->yorder, sf->ymaxmin, sf->yrange, sf->ybasis);
	    break;
        case eTNX_POLYNOMIAL:
            wf_gsb1pol (x, sf->xorder, sf->xbasis);
            wf_gsb1pol (y, sf->yorder, sf->ybasis);
	    break;
        default:
            fprintf (stderr,"TNX_GSEVAL: unknown surface type\n");
	    return (0.0);
        }

    /* Initialize accumulator basis functions */
    sum = 0.0;

    /* Loop over y basis functions */
    if (sf->xorder > sf->yorder)
	maxorder = sf->xorder + 1;
    else
	maxorder = sf->yorder + 1;
    xorder = sf->xorder;
    ii = 0;

    for (i = 0; i < sf->yorder; i++) {

	/* Loop over the x basis functions */
	accum = 0.0;
	for (k = 0; k < xorder; k++) {
	    accum = accum + sf->coeff[ii] * sf->xbasis[k];
	    ii = ii + 1;
	    }
	accum = accum * sf->ybasis[i];
	sum = sum + accum;

        /* Elements of the coefficient vector where neither k = 1 or i = 1
           are not calculated if sf->xterms = no. */
        if (sf->xterms == eTNX_XNONE)
            xorder = 1;
        else if (sf->xterms == eTNX_XHALF) {
            if ((i + 1 + sf->xorder + 1) > maxorder)
                xorder = xorder - 1;
	    }
        }

    return (sum);

}//close wf_gseval()


void WCSUtils::wf_gsb1cheb (double x,int order,double k1,double k2,double* basis)
{
    int i;
    double xnorm;

    basis[0] = 1.0;
    if (order == 1)
	return;

    xnorm = (x + k1) * k2;
    basis[1] = xnorm;
    if (order == 2)
	return;

    for (i = 2; i < order; i++)
	basis[i] = 2. * xnorm * basis[i-1] - basis[i-2];

    return;

}//close wf_gsb1cheb()


void WCSUtils::wf_gsb1leg (double x,int order,double k1,double k2,double* basis)
{
    int i;
    double ri, xnorm;

    basis[0] = 1.0;
    if (order == 1)
	return;

    xnorm = (x + k1) * k2 ;
    basis[1] = xnorm;
    if (order == 2)
        return;

    for (i = 2; i < order; i++) {
	ri = i;
        basis[i] = ((2.0 * ri - 1.0) * xnorm * basis[i-1] -
                       (ri - 1.0) * basis[i-2]) / ri;
        }

    return;

}//close wf_gsb1leg()



void WCSUtils::wf_gsb1pol (double x,int order,double* basis)
{
    int     i;

    basis[0] = 1.0;
    if (order == 1)
	return;

    basis[1] = x;
    if (order == 2)
	return;

    for (i = 2; i < order; i++)
	basis[i] = x * basis[i-1];

    return;

}//close wf_gsb1pol()


int WCSUtils::zpxpos (double xpix,double ypix,WCS* wcs,double* xpos,double* ypos)
{
	//Define this here (not as global)
	double TOL= 1e-13;

    int	i, j, k, ira, idec;
    double x, y, r, phi, theta, costhe, sinthe, dphi, cosphi, sinphi, dlng, z;
    double colatp, coslatp, sinlatp, longp;
    double xs, ys, ra, dec, xp, yp;
    double a, b, c, d, zd, zd1, zd2, r1, r2, rt, lambda;
    //double wf_gseval();

    /* Convert from pixels to image coordinates */
    xpix = xpix - wcs->crpix[0];
    ypix = ypix - wcs->crpix[1];

    /* Scale and rotate using CD matrix */
    if (wcs->rotmat) {
	x = xpix * wcs->cd[0] + ypix * wcs->cd[1];
	y = xpix * wcs->cd[2] + ypix * wcs->cd[3];
	}

    else {

	/* Check axis increments - bail out if either 0 */
	if (wcs->cdelt[0] == 0.0 || wcs->cdelt[1] == 0.0) {
	    *xpos = 0.0;
	    *ypos = 0.0;
	    return 2;
	    }

	/* Scale using CDELT */
	xs = xpix * wcs->cdelt[0];
	ys = ypix * wcs->cdelt[1];

	/* Take out rotation from CROTA */
	if (wcs->rot != 0.0) {
	    double cosr = cos (degrad_i (wcs->rot));
	    double sinr = sin (degrad_i (wcs->rot));
	    x = xs * cosr - ys * sinr;
	    y = xs * sinr + ys * cosr;
    	    }
	else {
	    x = xs;
	    y = ys;
	    }
	}

    /* Get the axis numbers */
    if (wcs->coorflip) {
	ira = 1;
	idec = 0;
	}
    else {
	ira = 0;
	idec = 1;
	}
    colatp = degrad_i (90.0 - wcs->crval[idec]);
    coslatp = cos(colatp);
    sinlatp = sin(colatp);
    longp = degrad_i(wcs->longpole);

    /*  Compute native spherical coordinates phi and theta in degrees from the
	projected coordinates. this is the projection part of the computation */
    k = wcs->zpnp;
    if (wcs->lngcor != NULL)
	xp = x + wf_gseval (wcs->lngcor, x, y);
    else
	xp = x;
    if (wcs->latcor != NULL)
	yp = y + wf_gseval (wcs->latcor, x, y);
    else
	yp = y;
    x = xp;
    y = yp;
    r = sqrt (x * x + y * y) / wcs->rodeg;

    /* Solve */

    /* Constant no solution */
    if (k < 1) {
        *xpos = kBADCVAL;
        *ypos = kBADCVAL;
	return (1);
	}

    /* Linear */
    else if (k == 1) {
        zd = (r - wcs->prj.p[0]) / wcs->prj.p[1];
	}

    /* Quadratic */
    else if (k == 2) {

        a = wcs->prj.p[2];
        b = wcs->prj.p[1];
        c = wcs->prj.p[0] - r;
	d = b * b - 4. * a * c;
	if (d < 0.) {
	    *xpos = kBADCVAL;
	    *ypos = kBADCVAL;
	    return (1);
	    }
	d = sqrt (d);

	/* Choose solution closest to the pole */
	zd1 = (-b + d) / (2. * a);
	zd2 = (-b - d) / (2. * a);
	if (zd1 < zd2)
	    zd = zd1;
	else
	    zd = zd2;
	if (zd < -TOL) {
	    if (zd1 > zd2)
		zd = zd1;
	    else
		zd = zd2;
	    }
	if (zd < 0.) {
	    if (zd < -TOL) {
		*xpos = kBADCVAL;
		*ypos = kBADCVAL;
		return (1);
		}
	    zd = 0.;
	    }
	else if (zd > kPI) {
	    if (zd > (kPI + TOL)) {
		*xpos = kBADCVAL;
		*ypos = kBADCVAL;
		return (1);
		}
	    zd = kPI;
	    }
	}

    /* Higher order solve iteratively */
    else {

        zd1 = 0.;
	r1 = wcs->prj.p[0];
	zd2 = wcs->zpzd;
	r2 = wcs->zpr;

	if (r < r1) {
	    if (r < (r1 - TOL)) {
		*xpos = kBADCVAL;
		*ypos = kBADCVAL;
		return (1);
		}
	    zd = zd1;
	    }
	else if (r > r2) {
	    if (r > (r2 + TOL)) {
		*xpos = kBADCVAL;
		*ypos = kBADCVAL;
		return (1);
		}
	    zd = zd2;
	    }
	else {
	    for (j=0; j<100; j++) {
	        lambda = (r2 - r) / (r2 - r1);
		if (lambda < 0.1)
		    lambda = 0.1;
		else if (lambda > 0.9)
		    lambda = 0.9;
		zd = zd2 - lambda * (zd2 - zd1);
		rt = 0.;
		for (i=k; i>=0; i--)
		    rt = (rt * zd) + wcs->prj.p[i];
		if (rt < r) {
		    if ((r - rt) < TOL)
		        break;
		    r1 = rt;
		    zd1 = zd;
		    }
		else {
		    if ((rt - r) < TOL)
		        break;
		    r2 = rt;
		    zd2 = zd;
		    }
		lambda = zd2 - zd1;
		lambda = fabs (zd2 - zd1);
		if (fabs (zd2 - zd1) < TOL)
		    break;
		}
	    }
	}

    /* Compute phi */
    if (r == 0.0)
	phi = 0.0;
    else
	phi = atan2 (x, -y);

    /* Compute theta */
    theta = kPI / 2 - zd;

    /*  Compute the celestial coordinates ra and dec from the native
	coordinates phi and theta. this is the spherical geometry part
	of the computation */

    costhe = cos (theta);
    sinthe = sin (theta);
    dphi = phi - longp;
    cosphi = cos (dphi);
    sinphi = sin (dphi);

    /* Compute the ra */
    x = sinthe * sinlatp - costhe * coslatp * cosphi;
    if (fabs (x) < kSPHTOL)
	x = -cos (theta + colatp) + costhe * coslatp * (1.0 - cosphi);
    y = -costhe * sinphi;
    if (x != 0.0 || y != 0.0)
	dlng = atan2 (y, x);
    else
	dlng = dphi + kPI ;
    ra =  wcs->crval[ira] + raddeg_i(dlng);

    /* normalize ra */
    if (wcs->crval[ira] >= 0.0) {
	if (ra < 0.0)
	    ra = ra + 360.0;
	}
    else {
	if (ra > 0.0)
	    ra = ra - 360.0;
	}
    if (ra > 360.0)
	ra = ra - 360.0;
    else if (ra < -360.0)
	ra = ra + 360.0;

    /* compute the dec */
    if (fmod (dphi, kPI) == 0.0) {
	dec = raddeg_i(theta + cosphi * colatp);
	if (dec > 90.0)
	    dec = 180.0 - dec;
	if (dec < -90.0)
	    dec = -180.0 - dec;
	}
    else {
	z = sinthe * coslatp + costhe * sinlatp * cosphi;
	if (fabs(z) > 0.99) {
	    if (z >= 0.0)
		dec = raddeg_i(acos (sqrt(x * x + y * y)));
	    else
		dec = raddeg_i(-acos (sqrt(x * x + y * y)));
	    }
	else
		dec = raddeg_i(asin (z));
	}

    /* store the results */
    *xpos  = ra;
    *ypos = dec;
    
	return (0);

}//close zpxpos


int WCSUtils::worldpos(double xpix,double ypix,WCS* wcs,double* xpos,double* ypos)
{
  double cosr, sinr, dx, dy, dz, tx;
  double sins, coss, dt, l, m, mg, da, dd, cos0, sin0;
  double rat = 0.0;
  double dect = 0.0;
  double mt, a, y0, td, r2;  /* allan: for COE */
  double dec0, ra0, decout, raout;
  double geo1, geo2, geo3;
  double cond2r=1.745329252e-2;
  double twopi = 6.28318530717959;
  double deps = 1.0e-5;

  /* Structure elements */
  double xref;		/* X reference coordinate value (deg) */
  double yref;		/* Y reference coordinate value (deg) */
  double xrefpix;	/* X reference pixel */
  double yrefpix;	/* Y reference pixel */
  double xinc;		/* X coordinate increment (deg) */
  double yinc;		/* Y coordinate increment (deg) */
  double rot;		/* Optical axis rotation (deg)  (N through E) */
  int itype = wcs->prjcode;

  /* Set local projection parameters */
  xref = wcs->xref;
  yref = wcs->yref;
  xrefpix = wcs->xrefpix;
  yrefpix = wcs->yrefpix;
  xinc = wcs->xinc;
  yinc = wcs->yinc;
  rot = degrad_i (wcs->rot);
  cosr = cos (rot);
  sinr = sin (rot);

  /* Offset from ref pixel */
  dx = xpix - xrefpix;
  dy = ypix - yrefpix;

  /* Scale and rotate using CD matrix */
  if (wcs->rotmat) {
    tx = dx * wcs->cd[0] + dy * wcs->cd[1];
    dy = dx * wcs->cd[2] + dy * wcs->cd[3];
    dx = tx;
    }

  /* Scale and rotate using CDELTn and CROTA2 */
  else {

    /* Check axis increments - bail out if either 0 */
    if ((xinc==0.0) || (yinc==0.0)) {
      *xpos=0.0;
      *ypos=0.0;
      return 2;
      }

    /* Scale using CDELT */
    dx = dx * xinc;
    dy = dy * yinc;

    /* Take out rotation from CROTA */
    if (rot != 0.0) {
      tx = dx * cosr - dy * sinr;
      dy = dx * sinr + dy * cosr;
      dx = tx;
      }
    }

  /* Flip coordinates if necessary */
  if (wcs->coorflip) {
    tx = dx;
    dx = dy;
    dy = tx;
    }

  /* Default, linear result for error or pixel return  */
  *xpos = xref + dx;
  *ypos = yref + dy;
  if (itype <= 0)
    return 0;

  /* Convert to radians  */
  if (wcs->coorflip) {
    dec0 = degrad_i (xref);
    ra0 = degrad_i (yref);
    }
  else {
    ra0 = degrad_i (xref);
    dec0 = degrad_i (yref);
    }
  l = degrad_i (dx);
  m = degrad_i (dy);
  sins = l*l + m*m;
  decout = 0.0;
  raout = 0.0;
  cos0 = cos (dec0);
  sin0 = sin (dec0);

  /* Process by case  */
  switch (itype) {

    case eWCS_CAR:   /* -CAR Cartesian (was eWCS_PIX pixel and eWCS_LIN linear) */
      rat =  ra0 + l;
      dect = dec0 + m;
      break;

    case eWCS_SIN: /* -SIN sin*/ 
      if (sins>1.0) return 1;
      coss = sqrt (1.0 - sins);
      dt = sin0 * coss + cos0 * m;
      if ((dt>1.0) || (dt<-1.0)) return 1;
      dect = asin (dt);
      rat = cos0 * coss - sin0 * m;
      if ((rat==0.0) && (l==0.0)) return 1;
      rat = atan2 (l, rat) + ra0;
      break;

    case eWCS_TAN:   /* -TAN tan */
    case eWCS_TNX:   /* -TNX tan with polynomial correction */
    case eWCS_TPV:   /* -TPV tan with polynomial correction */
    case eWCS_ZPX:   /* -ZPX zpn with polynomial correction */
      if (sins>1.0) return 1;
      dect = cos0 - m * sin0;
      if (dect==0.0) return 1;
      rat = ra0 + atan2 (l, dect);
      dect = atan (cos(rat-ra0) * (m * cos0 + sin0) / dect);
      break;

    case eWCS_ARC:   /* -ARC Arc*/
      if (sins>=twopi*twopi/4.0) return 1;
      sins = sqrt(sins);
      coss = cos (sins);
      if (sins!=0.0) sins = sin (sins) / sins;
      else
	sins = 1.0;
      dt = m * cos0 * sins + sin0 * coss;
      if ((dt>1.0) || (dt<-1.0)) return 1;
      dect = asin (dt);
      da = coss - dt * sin0;
      dt = l * sins * cos0;
      if ((da==0.0) && (dt==0.0)) return 1;
      rat = ra0 + atan2 (dt, da);
      break;

    case eWCS_NCP:   /* -NCP North celestial pole*/
      dect = cos0 - m * sin0;
      if (dect==0.0) return 1;
      rat = ra0 + atan2 (l, dect);
      dt = cos (rat-ra0);
      if (dt==0.0) return 1;
      dect = dect / dt;
      if ((dect>1.0) || (dect<-1.0)) return 1;
      dect = acos (dect);
      if (dec0<0.0) dect = -dect;
      break;

    case eWCS_GLS:   /* -GLS global sinusoid */
    case eWCS_SFL:   /* -SFL Samson-Flamsteed */
      dect = dec0 + m;
      if (fabs(dect)>twopi/4.0) return 1;
      coss = cos (dect);
      if (fabs(l)>twopi*coss/2.0) return 1;
      rat = ra0;
      if (coss>deps) rat = rat + l / coss;
      break;

    case eWCS_MER:   /* -MER mercator*/
      dt = yinc * cosr + xinc * sinr;
      if (dt==0.0) dt = 1.0;
      dy = degrad_i (yref/2.0 + 45.0);
      dx = dy + dt / 2.0 * cond2r;
      dy = log (tan (dy));
      dx = log (tan (dx));
      geo2 = degrad_i (dt) / (dx - dy);
      geo3 = geo2 * dy;
      geo1 = cos (degrad_i (yref));
      if (geo1<=0.0) geo1 = 1.0;
      rat = l / geo1 + ra0;
      if (fabs(rat - ra0) > twopi) return 1; /* added 10/13/94 DCW/EWG */
      dt = 0.0;
      if (geo2!=0.0) dt = (m + geo3) / geo2;
      dt = exp (dt);
      dect = 2.0 * atan (dt) - twopi / 4.0;
      break;

    case eWCS_AIT:   /* -AIT Aitoff*/
      dt = yinc*cosr + xinc*sinr;
      if (dt==0.0) dt = 1.0;
      dt = degrad_i (dt);
      dy = degrad_i (yref);
      dx = sin(dy+dt)/sqrt((1.0+cos(dy+dt))/2.0) -
	  sin(dy)/sqrt((1.0+cos(dy))/2.0);
      if (dx==0.0) dx = 1.0;
      geo2 = dt / dx;
      dt = xinc*cosr - yinc* sinr;
      if (dt==0.0) dt = 1.0;
      dt = degrad_i (dt);
      dx = 2.0 * cos(dy) * sin(dt/2.0);
      if (dx==0.0) dx = 1.0;
      geo1 = dt * sqrt((1.0+cos(dy)*cos(dt/2.0))/2.0) / dx;
      geo3 = geo2 * sin(dy) / sqrt((1.0+cos(dy))/2.0);
      rat = ra0;
      dect = dec0;
      if ((l==0.0) && (m==0.0)) break;
      dz = 4.0 - l*l/(4.0*geo1*geo1) - ((m+geo3)/geo2)*((m+geo3)/geo2) ;
      if ((dz>4.0) || (dz<2.0)) return 1;;
      dz = 0.5 * sqrt (dz);
      dd = (m+geo3) * dz / geo2;
      if (fabs(dd)>1.0) return 1;;
      dd = asin (dd);
      if (fabs(cos(dd))<deps) return 1;;
      da = l * dz / (2.0 * geo1 * cos(dd));
      if (fabs(da)>1.0) return 1;;
      da = asin (da);
      rat = ra0 + 2.0 * da;
      dect = dd;
      break;

    case eWCS_STG:   /* -STG Sterographic*/
      dz = (4.0 - sins) / (4.0 + sins);
      if (fabs(dz)>1.0) return 1;
      dect = dz * sin0 + m * cos0 * (1.0+dz) / 2.0;
      if (fabs(dect)>1.0) return 1;
      dect = asin (dect);
      rat = cos(dect);
      if (fabs(rat)<deps) return 1;
      rat = l * (1.0+dz) / (2.0 * rat);
      if (fabs(rat)>1.0) return 1;
      rat = asin (rat);
      mg = 1.0 + sin(dect) * sin0 + cos(dect) * cos0 * cos(rat);
      if (fabs(mg)<deps) return 1;
      mg = 2.0 * (sin(dect) * cos0 - cos(dect) * sin0 * cos(rat)) / mg;
      if (fabs(mg-m)>deps) rat = twopi/2.0 - rat;
      rat = ra0 + rat;
      break;

    case eWCS_COE:    /* COE projection code from Andreas Wicenic, ESO */
      td = tan (dec0);
      y0 = 1.0 / td;
      mt = y0 - m;
      if (dec0 < 0.)
	a = atan2 (l,-mt);
      else
	a = atan2 (l, mt);
      rat = ra0 - (a / sin0);
      r2 = (l * l) + (mt * mt);
      dect = asin (1.0 / (sin0 * 2.0) * (1.0 + sin0*sin0 * (1.0 - r2)));
      break;
  }

  /* Return RA in range  */
  raout = rat;
  decout = dect;
  if (raout-ra0>twopi/2.0) raout = raout - twopi;
  if (raout-ra0<-twopi/2.0) raout = raout + twopi;
  if (raout < 0.0) raout += twopi; /* added by DCW 10/12/94 */

  /* Convert units back to degrees  */
  *xpos = raddeg_i (raout);
  *ypos = raddeg_i (decout);

  return 0;

}//close worldpos()


void WCSUtils::wcscon(int sys1,int sys2,double eq1,double eq2,double* dtheta,double* dphi,double epoch)
{
	//void fk5prec(), fk4prec();

  // Set equinoxes if 0.0
    if (eq1 == 0.0) {
	if (sys1 == eWCS_B1950)
	    eq1 = 1950.0;
	else
	    eq1 = 2000.0;
	}
    if (eq2 == 0.0) {
	if (sys2 == eWCS_B1950)
	    eq2 = 1950.0;
	else
	    eq2 = 2000.0;
	}

    /* Set systems and equinoxes so that ICRS coordinates are not precessed */
    if (sys1 == eWCS_ICRS && sys2 == eWCS_ICRS)
	eq2 = eq1;

    if (sys1 == eWCS_J2000 && sys2 == eWCS_ICRS && eq1 == 2000.0) {
	eq2 = eq1;
	sys1 = sys2;
	}

    if (sys1 == eWCS_ICRS && sys2 == eWCS_J2000 && eq2 == 2000.0) {
	eq1 = eq2;
	sys1 = sys2;
	}

    /* If systems and equinoxes are the same, return */
    if (sys2 == sys1 && eq1 == eq2)
	return;

    /* Precess from input equinox, if necessary */
    if (eq1 != eq2) {
	if (sys1 == eWCS_B1950 && eq1 != 1950.0)
	   fk4prec (eq1, 1950.0, dtheta, dphi);
	if (sys1 == eWCS_J2000 && eq1 != 2000.0)
	   fk5prec (eq1, 2000.0, dtheta, dphi);
	}

    /* Convert to B1950 FK4 */
    if (sys2 == eWCS_B1950) {
	if (sys1 == eWCS_J2000) {
	    if (epoch > 0)
		fk524e (dtheta, dphi, epoch);
	    else
		fk524 (dtheta, dphi);
	    }
	else if (sys1 == eWCS_GALACTIC) 
	    gal2fk4 (dtheta, dphi);
	else if (sys1 == eWCS_ECLIPTIC) {
	    if (epoch > 0)
		ecl2fk4 (dtheta, dphi, epoch);
	    else
		ecl2fk4 (dtheta, dphi, 1950.0);
	    }
	}

    else if (sys2 == eWCS_J2000) {
        if (sys1 == eWCS_B1950) {
            if (epoch > 0)
                fk425e (dtheta, dphi, epoch);
            else
                fk425 (dtheta, dphi);
            }
        else if (sys1 == eWCS_GALACTIC)
            gal2fk5 (dtheta, dphi);
	else if (sys1 == eWCS_ECLIPTIC) {
	    if (epoch > 0)
		ecl2fk5 (dtheta, dphi, epoch);
	    else
		ecl2fk5 (dtheta, dphi, 2000.0);
	    }
	}

    else if (sys2 == eWCS_GALACTIC) {
        if (sys1 == eWCS_B1950)
	    fk42gal (dtheta, dphi);
        else if (sys1 == eWCS_J2000)
	    fk52gal (dtheta, dphi);
        else if (sys1 == eWCS_ECLIPTIC) {
	    if (epoch > 0)
		ecl2fk5 (dtheta, dphi, epoch);
	    else
		ecl2fk5 (dtheta, dphi, 2000.0);
	    fk52gal (dtheta, dphi);
	    }
	}

    else if (sys2 == eWCS_ECLIPTIC) {
        if (sys1 == eWCS_B1950) {
	    if (epoch > 0)
		fk42ecl (dtheta, dphi, epoch);
	    else
		fk42ecl (dtheta, dphi, 1950.0);
	    }
        else if (sys1 == eWCS_J2000) {
	    if (epoch > 0)
		fk52ecl (dtheta, dphi, epoch);
	    else
		fk52ecl (dtheta, dphi, 2000.0);
	    }
        else if (sys1 == eWCS_GALACTIC) {
	    gal2fk5 (dtheta, dphi);
	    if (epoch > 0)
		fk52ecl (dtheta, dphi, epoch);
	    else
		fk52ecl (dtheta, dphi, 2000.0);
	    }
	}

    /* Precess to desired equinox, if necessary */
    if (eq1 != eq2) {
	if (sys2 == eWCS_B1950 && eq2 != 1950.0)
	    fk4prec (1950.0, eq2, dtheta, dphi);
	if (sys2 == eWCS_J2000 && eq2 != 2000.0)
	    fk5prec (2000.0, eq2, dtheta, dphi);
	}

    /* Keep latitude/declination between +90 and -90 degrees */
    if (*dphi > 90.0) {
	*dphi = 180.0 - *dphi;
	*dtheta = *dtheta + 180.0;
	}
    else if (*dphi < -90.0) {
	*dphi = -180.0 - *dphi;
	*dtheta = *dtheta + 180.0;
	}

    /* Keep longitude/right ascension between 0 and 360 degrees */
    if (*dtheta > 360.0)
	*dtheta = *dtheta - 360.0;
    else if (*dtheta < 0.0)
	*dtheta = *dtheta + 360.0;

  return;

}//close wcscon()



void WCSUtils::gal2fk5 (double *dtheta,double *dphi)
{
    double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
    //void v2s3(),s2v3();
    int i;
    //char *eqcoor, *eqstrn();
		char* eqcoor;

    /*  Spherical to Cartesian */
    dl = *dtheta;
    db = *dphi;
    rl = degrad_i (dl);
    rb = degrad_i (db);
    r = 1.0;
    s2v3 (rl,rb,r,pos);

    /*  Rotate to equatorial coordinates */
    for (i = 0; i < 3; i++) {
	    pos1[i] = pos[0]*g_jgal[0][i] + pos[1]*g_jgal[1][i] + pos[2]*g_jgal[2][i];
	    }

    /*  Cartesian to Spherical */
    v2s3 (pos1,&rra,&rdec,&r);
    dra = raddeg_i (rra);
    ddec = raddeg_i (rdec);
    *dtheta = dra;
    *dphi = ddec;

    /*  Print result if in diagnostic mode */
    if (g_idg) {
	fprintf (stderr,"GAL2FK5: long = %.5f lat = %.5f\n",dl,db);
	eqcoor = eqstrn (dra,ddec);
	fprintf (stderr,"GAL2FK5: J2000 RA,Dec= %s\n",eqcoor);
	free (eqcoor);
	}

  return;

}//close gal2fk5()


char* WCSUtils::eqstrn (double dra,double ddec)
{
	char	*eqcoor;	/* ASCII character string of position (returned) */
	char	decp;
	int	rah,irm,decd,decm;
	double	xpos,ypos,xp,yp,ras,decs;

    /*  Right ascension to hours, minutes, and seconds */
    xpos = dra / 15.0;
    rah = (int) xpos;
    xp = (double) 60.0 * (xpos - (double) rah);
    irm = (int) xp;
    ras = (double) 60.0 * (xp - (double) irm);

    /* Declination to degrees, minutes, seconds */
    if (ddec < 0) {
	ypos = -ddec;
	decp = '-';
	}
    else {
	decp = '+';
	ypos = ddec;
	}
    decd = (int) ypos;
    yp = (double) 60.0 * (ypos - (double) decd);
    decm = (int) yp;
    decs = (double) 60.0 * (yp - (double) decm);

    eqcoor = (char*)malloc (32);
    (void)sprintf (eqcoor,"%02d:%02d:%06.3f %c%02d:%02d:%05.2f",
		   rah,irm,ras,decp,decd,decm,decs);
    if (eqcoor[6] == ' ')
	eqcoor[6] = '0';
    if (eqcoor[20] == ' ')
	eqcoor[20] = '0';

  return (eqcoor);

}//close eqstrn()



void WCSUtils::s2v3 (double rra,double rdec,double r,double pos[3])
{
	pos[0] = r * cos (rra) * cos (rdec);
  pos[1] = r * sin (rra) * cos (rdec);
  pos[2] = r * sin (rdec);

  return;

}//close s2v3()


/* Convert geocentric equatorial rectangular coordinates to
   right ascension, declination, and distance */

void WCSUtils::v2s3 (double pos[3],double *rra,double *rdec,double *r)
{
	double x,y,z,rxy,rxy2,z2;

    x = pos[0];
    y = pos[1];
    z = pos[2];

    *rra = atan2 (y, x);

    /* Keep RA within 0 to 2pi range */
    if (*rra < 0.0)
	*rra = *rra + (2.0 * kPI);
    if (*rra > 2.0 * kPI)
	*rra = *rra - (2.0 * kPI);

    rxy2 = x*x + y*y;
    rxy = sqrt (rxy2);
    *rdec = atan2 (z, rxy);

    z2 = z * z;
    *r = sqrt (rxy2 + z2);

  return;

}//close v2s3()


void WCSUtils::fk52gal (double *dtheta,double *dphi)
{
    double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
    //void v2s3(),s2v3();
    //char *eqcoor, *eqstrn();
		char *eqcoor;
    int i;

    /*  Spherical to cartesian */
    dra = *dtheta;
    ddec = *dphi;
    rra = degrad_i (dra);
    rdec = degrad_i (ddec);
    r = 1.0;
    (void)s2v3 (rra,rdec,r,pos);

    /*  Rotate to galactic */
    for (i = 0; i < 3; i++) {
	pos1[i] = pos[0]*g_jgal[i][0] + pos[1]*g_jgal[i][1] + pos[2]*g_jgal[i][2];
	}

    /*  Cartesian to spherical */
    v2s3 (pos1,&rl,&rb,&r);

    dl = raddeg_i (rl);
    db = raddeg_i (rb);
    *dtheta = dl;
    *dphi = db;

    /*  Print result if in diagnostic mode */
    if (g_idg) {
	eqcoor = eqstrn (dra,ddec);
	fprintf (stderr,"FK52GAL: J2000 RA,Dec= %s\n",eqcoor);
	fprintf (stderr,"FK52GAL: long = %.5f lat = %.5f\n",dl,db);
	free (eqcoor);
	}

  return;

}//close fk52gal()

void WCSUtils::fk425 (double	*ra,double	*dec)
{
	double	rapm;		/* Proper motion in right ascension */
	double	decpm;		/* Proper motion in declination  */
	/* In: rad/trop.yr.  Out:  rad/jul.yr. */

  rapm = (double) 0.0;
  decpm = (double) 0.0;
  fk425m (ra, dec, &rapm, &decpm);
    
	return;

}//close fk425()

void WCSUtils::fk425e (double* ra,double* dec,double epoch)
{
	double	rapm;		/* Proper motion in right ascension */
	double	decpm;		/* Proper motion in declination  */
	/* In: rad/trop.yr.  Out:  rad/jul.yr. */

  rapm = (double) 0.0;
  decpm = (double) 0.0;
  fk425m (ra, dec, &rapm, &decpm);
  *ra = *ra + (rapm * (epoch - 2000.0));
  *dec = *dec + (decpm * (epoch - 2000.0));
   
	return;

}//close fk425e()

void WCSUtils::fk5prec(double ep0,double ep1,double* ra,double* dec)
{
    int i, j;
    double pm[9], *pmi, v1[3], v2[3], rra, rdec, r;
    //void v2s3(),s2v3(), mprecfk5();

    rra = degrad_i (*ra);
    rdec = degrad_i (*dec);
    r = 1.0;
 
    /* Generate appropriate precession matrix */
    mprecfk5 (ep0, ep1, pm);
 
    /* Convert RA,Dec to x,y,z */
    s2v3 (rra, rdec, r, v1);
 
    /* Multiply position vector by precession matrix */
    pmi = pm;
    for (i = 0; i < 3; i++) {
	v2[i] = 0;
	for (j = 0; j < 3; j++)
	    v2[i] = v2[i] + ( v1[j] * *pmi++ );
	}
 
    /* Back to RA,Dec */
    v2s3 (v2, &rra, &rdec, &r);

    /* Convert from radians to degrees */
    *ra = raddeg_i (rra);
    *dec = raddeg_i (rdec);
    
	return;

}//close fk5prec()



void WCSUtils::fk4prec(double ep0,double ep1,double* ra,double* dec)
{
	int i, j;
  double pm[9], *pmi, v1[3], v2[3], rra, rdec, r;
  //void v2s3(),s2v3(), mprecfk4();

  rra = degrad_i (*ra);
  rdec = degrad_i (*dec);
  r = 1.0;
 
  // Generate appropriate precession matrix
  mprecfk4 ( ep0, ep1, pm );
 
    /* Convert RA,Dec to x,y,z */
    s2v3 (rra, rdec, r, v1);
 
    /* Multiply position vector by precession matrix */
    pmi = pm;
    for (i = 0; i < 3; i++) {
	v2[i] = 0;
	for (j = 0; j < 3; j++)
	    v2[i] = v2[i] + (*pmi++ * v1[j]);
	}
 
    /* Back to RA,Dec */
    v2s3 (v2, &rra, &rdec, &r);

  /* Convert from radians to degrees */
  *ra = raddeg_i (rra);
  *dec = raddeg_i (rdec);

}//close fk4prec()


void WCSUtils::fk524 (double* ra,double* dec)
{
    double	rapm;	/* Proper motion in right ascension */
    double	decpm;	/* Proper motion in declination  */
			/* In:  deg/jul.yr.  Out: deg/trop.yr.  */

    rapm = (double) 0.0;
    decpm = (double) 0.0;
    fk524m (ra, dec, &rapm, &decpm);
    
	return;

}//close fk524()


void WCSUtils::fk524e(double* ra,double* dec,double epoch)
{    
	double	rapm;	/* Proper motion in right ascension */
  double	decpm;	/* Proper motion in declination  */
	/* In:  deg/jul.yr.  Out: deg/trop.yr.  */

  rapm = (double) 0.0;
  decpm = (double) 0.0;
  fk524m (ra, dec, &rapm, &decpm);
  *ra = *ra + (rapm * (epoch - 1950.0));
  *dec = *dec + (decpm * (epoch - 1950.0));
    
	return;

}//close fk524e()


void WCSUtils::fk524m(double* ra,double* dec,double* rapm,double* decpm)
{
	double parallax = 0.0;
  double rv = 0.0;

  fk524pv (ra, dec, rapm, decpm, &parallax, &rv);
  return;

}//close fk524m()

/*  
		This routine converts stars from the IAU 1976 FK5 Fricke
    system, to the old Bessel-Newcomb FK4 system, using Yallop's
    implementation (see ref 2) of a matrix method due to Standish
    (see ref 3).  The numerical values of ref 2 are used canonically.

 		Conversion from other than Julian epoch 2000.0 to other than Besselian
    epoch 1950.0 will require use of the appropriate precession, proper
    motion, and e-terms routines before and/or after fk524 is called.
 
 		In the FK4 catalogue the proper motions of stars within 10 degrees
    of the poles do not embody the differential e-term effect and should,
    strictly speaking, be handled in a different manner from stars outside
    these regions.  however, given the general lack of homogeneity of the
    star data available for routine astrometry, the difficulties of handling
    positions that may have been determined from astrometric fields spanning
    the polar and non-polar regions, the likelihood that the differential
    e-terms effect was not taken into account when allowing for proper motion
    in past astrometry, and the undesirability of a discontinuity in the
    algorithm, the decision has been made in this routine to include the
    effect of differential e-terms on the proper motions for all stars,
    whether polar or not, at epoch 2000, and measuring on the sky rather
    than in terms of dra, the errors resulting from this simplification are
    less than 1 milliarcsecond in position and 1 milliarcsecond per century
    in proper motion.

    References:

      1  "Mean and apparent place computations in the new IAU System.
          I. The transformation of astrometric catalog systems to the
 	  equinox J2000.0." Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Yallop, B.D.; Hohenkerk, C.Y.
 	  Astronomical Journal vol. 97, Jan. 1989, p. 265-273.

      2  "Mean and apparent place computations in the new IAU System.
	  II. Transformation of mean star places from FK4 B1950.0 to
 	  FK5 J2000.0 using matrices in 6-space."  Yallop, B.D.;
	  Hohenkerk, C.Y.; Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Astronomical Journal vol. 97, Jan. 1989,
	  p. 274-279.
 
      3  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement to
         the Astronomical Almanac", ISBN 0-935702-68-7.

      4  "Conversion of positions and proper motions from B1950.0 to the
	  IAU system at J2000.0", Standish, E.M.  Astronomy and
	  Astrophysics, vol. 115, no. 1, Nov. 1982, p. 20-22.

   P.T.Wallace   Starlink   19 December 1993
   Doug Mink     Smithsonian Astrophysical Observatory 1 November 2000 
*/
void WCSUtils::fk524pv(double* ra,double* dec,double* rapm,double* decpm,double* parallax,double* rv)
{
    double r2000,d2000;		/* J2000.0 ra,dec (radians) */
    double r1950,d1950;		/* B1950.0 ra,dec (rad) */

    /* Miscellaneous */
    double ur,ud;
    double sr, cr, sd, cd, x, y, z, w, wd;
    double v1[6],v2[6];
    double xd,yd,zd;
    double rxyz, rxysq, rxy;
    double dra,ddec;
    int	i,j;
    int	diag = 0;

    /* Constants */
    double zero = (double) 0.0;
    double vf = 21.095;	/* Km per sec to AU per tropical century */
			/* = 86400 * 36524.2198782 / 149597870 */

    /* Convert J2000 RA and Dec from degrees to radians */
    r2000 = degrad_i (*ra);
    d2000 = degrad_i (*dec);

    /* Convert J2000 RA and Dec proper motion from degrees/year to arcsec/tc */
    ur = *rapm  * 360000.0;
    ud = *decpm * 360000.0;

    /* Spherical to Cartesian */
    sr = sin (r2000);
    cr = cos (r2000);
    sd = sin (d2000);
    cd = cos (d2000);

    x = cr * cd;
    y = sr * cd;
    z = sd;

    v1[0] = x;
    v1[1] = y;
    v1[2] = z;
 
    if (ur != zero || ud != zero) {
	v1[3] = -(ur*y) - (cr*sd*ud);
	v1[4] =  (ur*x) - (sr*sd*ud);
	v1[5] =          (cd*ud);
	}
    else {
	v1[3] = zero;
	v1[4] = zero;
	v1[5] = zero;
	}
 
    /* Convert position + velocity vector to bn system */
    for (i = 0; i < 6; i++) {
	w = zero;
	for (j = 0; j < 6; j++) {
	    w = w + g_emi[i][j] * v1[j];
	    }
	v2[i] = w;
	}
 
    /* Vector components */
    x = v2[0];
    y = v2[1];
    z = v2[2];

    /* Magnitude of position vector */
    rxyz = sqrt (x*x + y*y + z*z);
 
    /* Apply e-terms to position */
    w = (x * g_a[0]) + (y * g_a[1]) + (z * g_a[2]);
    x = x + (g_a[0] * rxyz) - (w * x);
    y = y + (g_a[1] * rxyz) - (w * y);
    z = z + (g_a[2] * rxyz) - (w * z);
 
    /* Recompute magnitude of position vector */
    rxyz = sqrt (x*x + y*y + z*z);

    /* Apply e-terms to position and velocity */
    x = v2[0];
    y = v2[1];
    z = v2[2];
    w = (x * g_a[0]) + (y * g_a[1]) + (z * g_a[2]);
    wd = (x * g_ad[0]) + (y * g_ad[1]) + (z * g_ad[2]);
    x = x + (g_a[0] * rxyz) - (w * x);
    y = y + (g_a[1] * rxyz) - (w * y);
    z = z + (g_a[2] * rxyz) - (w * z);
    xd = v2[3] + (g_ad[0] * rxyz) - (wd * x);
    yd = v2[4] + (g_ad[1] * rxyz) - (wd * y);
    zd = v2[5] + (g_ad[2] * rxyz) - (wd * z);

    /*  Convert to spherical  */
    rxysq = (x * x) + (y * y);
    rxy = sqrt (rxysq);

    /* Convert back to spherical coordinates */
    if (x == zero && y == zero)
	r1950 = zero;
    else {
	r1950 = atan2 (y,x);
	if (r1950 < zero)
	    r1950 = r1950 + d2pi;
	}
    d1950 = atan2 (z,rxy);

    if (rxy > g_tiny) {
	ur = (x*yd - y*xd) / rxysq;
	ud = (zd*rxysq - z * (x*xd + y*yd)) / ((rxysq + z*z) * rxy);
	}

    if (*parallax > g_tiny) {
	*rv = ((x * xd) + (y * yd) + (z * zd)) / (*parallax * vf * rxyz);
	*parallax = *parallax / rxyz;
	}

    /* Return results */
    *ra = raddeg_i (r1950);
    *dec = raddeg_i (d1950);
    *rapm  = ur / 360000.0;
    *decpm = ud / 360000.0;

    if (diag) {
	dra = 240.0 * raddeg_i (r1950 - r2000);
	ddec = 3600.0 * raddeg_i (d1950 - d2000);
	fprintf(stderr,"B1950-J2000: dra= %11.5f sec  ddec= %f11.5f arcsec\n",
		dra, ddec);
	}
 
	return;

}//close fk524pv()


/*--- Transform IAU 1958 galactic coordinates to B1950.0 'FK4'
 *    equatorial coordinates */

void WCSUtils::gal2fk4 (double* dtheta,double* dphi)
{
	double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
  //void v2s3(),s2v3();
  //char *eqcoor, *eqstrn();
	char *eqcoor;	
  int i;

    /*  spherical to cartesian */
    dl = *dtheta;
    db = *dphi;
    rl = degrad_i (dl);
    rb = degrad_i (db);
    r = 1.0;
    s2v3 (rl,rb,r,pos);

    /*  rotate to equatorial coordinates */
    for (i = 0; i < 3; i++) {
	pos1[i] = pos[0]*g_bgal[0][i] + pos[1]*g_bgal[1][i] + pos[2]*g_bgal[2][i];
	}

    /*  cartesian to spherical */
    v2s3 (pos1,&rra,&rdec,&r);

/*  introduce e-terms */
/*	jpabe (rra,rdec,-1,g_idg); */

    dra = raddeg_i (rra);
    ddec = raddeg_i (rdec);
    *dtheta = dra;
    *dphi = ddec;

    /*  print result if in diagnostic mode */
    if (g_idg) {
	fprintf (stderr,"GAL2FK4: long = %.5f lat = %.5f\n",dl,db);
	eqcoor = eqstrn (dra,ddec);
	fprintf (stderr,"GAL2FK4: B1950 RA,Dec= %s\n",eqcoor);
	free (eqcoor);
	}

    return;

}//close gal2fk4()



void WCSUtils::fk425pv (double* ra,double* dec,double* rapm,double* decpm,double* parallax,double* rv)
{
    double r1950,d1950;		/* B1950.0 ra,dec (rad) */
    double r2000,d2000;		/* J2000.0 ra,dec (rad) */

    /* Miscellaneous */
    double ur,ud,sr,cr,sd,cd,w,wd;
    double x,y,z,xd,yd,zd, dra,ddec;
    double rxyz, rxysq, rxy, rxyzsq, spxy, spxyz;
    int	i,j;
    int	diag = 0;

    double r0[3],rd0[3];	/* star position and velocity vectors */
    double v1[6],v2[6];		/* combined position and velocity vectors */

    /* Constants */
    double zero = (double) 0.0;
    double vf = 21.095;	/* Km per sec to AU per tropical century */
			/* = 86400 * 36524.2198782 / 149597870 */

    /* Convert B1950 RA and Dec from degrees to radians */
    r1950 = degrad_i (*ra);
    d1950 = degrad_i (*dec);

    /* Convert B1950 RA and Dec proper motion from degrees/year to arcsec/tc */
    ur = *rapm  * 360000.0;
    ud = *decpm * 360000.0;

    /* Convert direction to Cartesian */
    sr = sin (r1950);
    cr = cos (r1950);
    sd = sin (d1950);
    cd = cos (d1950);
    r0[0] = cr * cd;
    r0[1] = sr * cd;
    r0[2] = sd;

    /* Convert motion to Cartesian */
    w = vf * *rv * *parallax;
    if (ur != zero || ud != zero || (*rv != zero && *parallax != zero)) {
	rd0[0] = (-sr * cd * ur) - (cr * sd * ud) + (w * r0[0]);
	rd0[1] =  (cr * cd * ur) - (sr * sd * ud) + (w * r0[1]);
	rd0[2] = 	                (cd * ud) + (w * r0[2]);
	}
    else {
	rd0[0] = zero;
	rd0[1] = zero;
	rd0[2] = zero;
	}

    /* Remove e-terms from position and express as position+velocity 6-vector */
    w = (r0[0] * g_a[0]) + (r0[1] * g_a[1]) + (r0[2] * g_a[2]);
    for (i = 0; i < 3; i++)
	v1[i] = r0[i] - g_a[i] + (w * r0[i]);

    /* Remove e-terms from proper motion and express as 6-vector */
    wd = (r0[0] * g_ad[0]) + (r0[1] * g_ad[1]) + (r0[2] * g_ad[2]);
    for (i = 0; i < 3; i++)
	v1[i+3] = rd0[i] - g_ad[i] + (wd * r0[i]);

    /* Alternately: Put proper motion in 6-vector without adding e-terms
    for (i = 0; i < 3; i++)
	v1[i+3] = rd0[i]; */

    /* Convert position + velocity vector to FK5 system */
    for (i = 0; i < 6; i++) {
	w = zero;
	for (j = 0; j < 6; j++) {
	    w += g_em[i][j] * v1[j];
	    }
	v2[i] = w;
	}

    /* Vector components */
    x = v2[0];
    y = v2[1];
    z = v2[2];
    xd = v2[3];
    yd = v2[4];
    zd = v2[5];

    /* Magnitude of position vector */
    rxysq = x*x + y*y;
    rxy = sqrt (rxysq);
    rxyzsq = rxysq + z*z;
    rxyz = sqrt (rxyzsq);

    spxy = (x * xd) + (y * yd);
    spxyz = spxy + (z * zd);

    /* Convert back to spherical coordinates */
    if (x == zero && y == zero)
	r2000 = zero;
    else {
	r2000 = atan2 (y,x);
	if (r2000 < zero)
	    r2000 = r2000 + d2pi;
	}
    d2000 = atan2 (z,rxy);

    if (rxy > g_tiny) {
	ur = ((x * yd) - (y * xd)) / rxysq;
	ud = ((zd * rxysq) - (z * spxy)) / (rxyzsq * rxy);
	}

    if (*parallax > g_tiny) {
	*rv = spxyz / (*parallax * rxyz * vf);
	*parallax = *parallax / rxyz;
	}

    /* Return results */
    *ra = raddeg_i (r2000);
    *dec = raddeg_i (d2000);
    *rapm  = ur / 360000.0;
    *decpm = ud / 360000.0;

    if (diag) {
	dra = 240.0 * raddeg_i (r2000 - r1950);
	ddec = 3600.0 * raddeg_i (d2000 - d1950);
	fprintf(stderr,"J2000-B1950: dra= %11.5f sec  ddec= %f11.5f arcsec\n",
		dra, ddec);
	}
    
	return;

}//close fk425pv()


void WCSUtils::fk42gal (double *dtheta,double *dphi)
{
    double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
    //void v2s3(),s2v3();
    int i;
    //char *eqcoor, *eqstrn();
		char *eqcoor;

    dra = *dtheta;
    ddec = *dphi;
    rra = degrad_i (dra);
    rdec = degrad_i (ddec);

    /*  remove e-terms */
    /*	call jpabe (rra,rdec,-1,g_idg) */

    /*  Spherical to Cartesian */
    r = 1.;
    s2v3 (rra,rdec,r,pos);

    /*  rotate to galactic */
    for (i = 0; i<3; i++) {
	pos1[i] = pos[0]*g_bgal[i][0] + pos[1]*g_bgal[i][1] + pos[2]*g_bgal[i][2];
	}

    /*  Cartesian to spherical */
    v2s3 (pos1,&rl,&rb,&r);

    dl = raddeg_i (rl);
    db = raddeg_i (rb);
    *dtheta = dl;
    *dphi = db;

    /*  Print result if in diagnostic mode */
    if (g_idg) {
	eqcoor = eqstrn (dra,ddec);
	fprintf (stderr,"FK42GAL: B1950 RA,Dec= %s\n",eqcoor);
	fprintf (stderr,"FK42GAL: long = %.5f lat = %.5f\n",dl,db);
	free (eqcoor);
	}

  return;

}//close fk42gal()

void WCSUtils::ecl2fk4 (double *dtheta,double *dphi,double epoch)
{
	//void ecl2fk5(), fk524e();

  // Convert from ecliptic to J2000 coordinates
  ecl2fk5 (dtheta, dphi, epoch);

  // Convert from J2000 to B1950 coordinates
  fk524e (dtheta, dphi, epoch);

	return;

}//close ecl2fk4()


void WCSUtils::ecl2fk5 (double *dtheta,double *dphi,double epoch)
{
    int i, j;
    double rtheta, rphi, v1[3], v2[3];
    double t, eps0, r;
    double rmat[9];	/* Rotation matrix */
    //void v2s3(),s2v3(), fk5prec(), rotmat();

    rtheta = degrad_i (*dtheta);
    rphi = degrad_i (*dphi);

    /* Convert RA,Dec to x,y,z */
    r = 1.0;
    s2v3 (rtheta, rphi, r, v1);

    /* Interval between basic epoch J2000.0 and current epoch (JC) in centuries*/
    t = (epoch - 2000.0) * 0.01;
 
    /* Mean obliquity */
    eps0 = secrad_i ((84381.448 + (-46.8150 + (-0.00059 + 0.001813*t) * t) * t));
 
    /* Form the equatorial to ecliptic rotation matrix (IAU 1980 theory).
     *  References: Murray, C.A., Vectorial Astrometry, section 4.3.
     *    The matrix is in the sense   v[ecl]  =  rmat * v[equ];  the
     *    equator, equinox and ecliptic are mean of date. */
    rotmat (1, eps0, 0.0, 0.0, rmat);
 
    /* Multiply position vector by ecliptic to equatorial rotation matrix */
    for (i = 0; i < 3; i++) {
	v2[i] = 0;
	for (j = 0; j < 3; j++)
	    v2[i] = v2[i] + (rmat[3*j + i] * v1[j]);
	}

    /* Cartesian to spherical */
    v2s3 (v2, &rtheta, &rphi, &r);

    /* Convert from radians to degrees */
    *dtheta = raddeg_i (rtheta);
    *dphi = raddeg_i (rphi);

    if (epoch != 2000.0)
	fk5prec (epoch, 2000.0, dtheta, dphi);

}//close ecl2fk5()


void WCSUtils::fk42ecl (double *dtheta,double *dphi,double epoch)
{
	//void fk425e(), fk52ecl();

  // Convert from B1950 to J2000 coordinates
  fk425e (dtheta, dphi, epoch);

  // Convert from J2000 to ecliptic coordinates
  fk52ecl (dtheta, dphi, epoch);

  return;

}//close fk42ecl()


void WCSUtils::fk425m (double* ra,double* dec,double* rapm,double* decpm)
{
	double parallax = 0.0;
  double rv = 0.0;

  fk425pv (ra, dec, rapm, decpm, &parallax, &rv);
  	
	return;

}//close fk425m()


void WCSUtils::fk52ecl (double *dtheta,double *dphi,double epoch)
{
    int i, j;
    double t, eps0, rphi, rtheta;
    double v1[3], v2[3], r;
    double rmat[9], *rmati;	/* Rotation matrix  */

    //void rotmat(), v2s3(), s2v3(), fk5prec();

    /* Precess coordinates from J2000 to epoch */
    if (epoch != 2000.0)
	fk5prec (2000.0, epoch, dtheta, dphi);

    /* Convert from degrees to radians */
    rtheta = degrad_i (*dtheta);
    rphi = degrad_i (*dphi);

    /* Convert RA,Dec to x,y,z */
    r = 1.0;
    s2v3 (rtheta, rphi, r, v1);

    /* Interval between basic epoch J2000.0 and current epoch (JC) in centuries*/
    t = (epoch - 2000.0) * 0.01;
 
    /* Mean obliquity */
    eps0 = secrad_i ((84381.448 + (-46.8150 + (-0.00059 + 0.001813*t) * t) * t));
 
    /* Form the equatorial to ecliptic rotation matrix (IAU 1980 theory).
     *  References: Murray, C.A., Vectorial Astrometry, section 4.3.
     *    The matrix is in the sense   v[ecl]  =  rmat * v[equ];  the
     *    equator, equinox and ecliptic are mean of date. */
    rotmat (1, eps0, 0.0, 0.0, rmat);

    /* Multiply position vector by equatoria to eccliptic rotation matrix */
    rmati = rmat;
    for (i = 0; i < 3; i++) {
	v2[i] = 0;
	for (j = 0; j < 3; j++)
	    v2[i] = v2[i] + (*rmati++ * v1[j]);
	}

    /* Convert x,y,z to latitude, longitude */
    v2s3 (v2, &rtheta, &rphi, &r);

    /* Convert from radians to degrees */
    *dtheta = raddeg_i (rtheta);
    *dphi = raddeg_i (rphi);

}//close fk52ecl()


void WCSUtils::mprecfk5(double ep0,double ep1,double rmatp[9])
{
    double t0, t, tas2r, w, zeta, z, theta;
    //void rotmat();
 
    /* Interval between basic epoch J2000.0 and beginning epoch (JC) */
    t0 = ( ep0 - 2000.0 ) / 100.0;
 
    /* Interval over which precession required (JC) */
    t =  ( ep1 - ep0 ) / 100.0;
 
    /* Euler angles */
    tas2r = secrad_i (t);
    w = 2306.2181 + ( ( 1.39656 - ( 0.000139 * t0 ) ) * t0 );
    zeta = (w + ( ( 0.30188 - 0.000344 * t0 ) + 0.017998 * t ) * t ) * tas2r;
    z = (w + ( ( 1.09468 + 0.000066 * t0 ) + 0.018203 * t ) * t ) * tas2r;
    theta = ( ( 2004.3109 + ( - 0.85330 - 0.000217 * t0 ) * t0 )
	  + ( ( -0.42665 - 0.000217 * t0 ) - 0.041833 * t ) * t ) * tas2r;
 
    /* Rotation matrix */
    rotmat (323, -zeta, theta, -z, rmatp);
    
	return;

}//close mprecfk5()


void WCSUtils::mprecfk4(double bep0,double bep1,double rmatp[9])
{
    double bigt, t, tas2r, w, zeta, z, theta;
    //void rotmat();
 
    /* Interval between basic epoch B1850.0 and beginning epoch in TC */
    bigt  = ( bep0 - 1850.0 ) / 100.0;
 
    /* Interval over which precession required, in tropical centuries */
    t = ( bep1 - bep0 ) / 100.0;
 
    /* Euler angles */
    tas2r = secrad_i (t);
    w = 2303.5548 + ( 1.39720 + 0.000059 * bigt ) * bigt;
    zeta = (w + ( 0.30242 - 0.000269 * bigt + 0.017996 * t ) * t ) * tas2r;
    z = (w + ( 1.09478 + 0.000387 * bigt + 0.018324 * t ) * t ) * tas2r;
    theta = ( 2005.1125 + ( - 0.85294 - 0.000365* bigt ) * bigt +
	    ( - 0.42647 - 0.000365 * bigt - 0.041802 * t ) * t ) * tas2r;
 
    /* Rotation matrix */
    rotmat (323, -zeta, theta, -z, rmatp);
    
	return;

}//close mprecfk4()



void WCSUtils::rotmat(int axes,double rot1,double rot2,double rot3,double* matrix)
{
    int i, j, k, naxis, iaxes, iaxis;
    double rot[3], srot, crot, *mati, w, wm[9], *wmi, matn[9];
    int axis[3];

    /* Initial final rotation matrix */
    mati = matrix;
    for (i = 0; i < 3; i++) {
	for (j=0; j < 3; j++) {
	    if (i == j)
		*mati++ = 1.0;
	    else
		*mati++ = 0.0;
	    }
	}

    /* Separate digits of rotation axis string and count rotations */
    naxis = 0;
    iaxes = axes;
    axis[0] = iaxes / 100;
    if (axis[0] > 0) {
	naxis++;
	iaxes = iaxes - (100 * axis[0]);
	}
    axis[naxis] = iaxes / 10;
    if (axis[naxis] > 0) {
	iaxes = iaxes - (10 * axis[naxis]);
	naxis++;
	}
    axis[naxis] = iaxes;
    if (axis[naxis] > 0)
	naxis++;

    /* Set up rotation angles */
    rot[0] = rot1;
    rot[1] = rot2;
    rot[2] = rot3;

    /* For each digit of axis string, set up matrix */
    for (iaxis = 0; iaxis < naxis; iaxis++) {

	/* Initialize current rotation matrix */
	mati = matn;
	for (i = 0; i < 3; i++) {
	    for (j=0; j < 3; j++) {
		if (i == j)
		    *mati++ = 1.0;
		else
		    *mati++ = 0.0;
		}
	    }

	srot = sin (rot[iaxis]);
	crot = cos (rot[iaxis]);
	
	/* Matrix for rotation in X */
	if (axis[iaxis] == 1) {
	    matn[4] = crot;
	    matn[5] = srot;
	    matn[7] = -srot;
	    matn[8] = crot;
	    }
	
	/* Matrix for rotation in Y */
	else if (axis[iaxis] == 2) {
	    matn[0] = crot;
	    matn[2] = -srot;
	    matn[6] = srot;
	    matn[8] = crot;
	    }
	
	/* Matrix for rotation in Z */
	else {
	    matn[0] = crot;
	    matn[1] = srot;
	    matn[3] = -srot;
	    matn[4] = crot;
	    }

	/* Multiply existing rotation matrix by new rotation matrix */
	for (i = 0; i < 3; i++) {
	    for (j = 0; j < 3; j++) {
		w = 0.0;
		for (k = 0; k < 3; k++)
		    w+= matn[3*i + k] * matrix[3*k + j];
		wm[3*i + j] = w;
		}
	    }

	/* Update output matrix */
	mati = matrix;
	wmi = wm;
	for (i = 0; i < 9; i++) {
	    *mati++ = *wmi++;
	    }
	}

  return;

}//close rotmat()




double WCSUtils::wcsceq(char* wcstring)
{
 	if (wcstring[0] == 'J' || wcstring[0] == 'j' ||
	wcstring[0] == 'B' || wcstring[0] == 'b')
	return (atof (wcstring+1));
    else if (!strncmp (wcstring, "FK4",3) ||
	     !strncmp (wcstring, "fk4",3))
	return (1950.0);
    else if (!strncmp (wcstring, "FK5",3) ||
	     !strncmp (wcstring, "fk5",3))
	return (2000.0);
    else if (!strncmp (wcstring, "ICRS",4) ||
	     !strncmp (wcstring, "icrs",4))
	return (2000.0);
    else if (wcstring[0] == '1' || wcstring[0] == '2')
	return (atof (wcstring));
    else

	return (0.0);

}//close wcsceq()


int WCSUtils::sphrev (const double phi, const double theta, const double eul[5], double* lng, double* lat)
{
	const double tol = 1.0e-5;
  double cosphi, costhe, dlng, dphi, sinphi, sinthe, x, y, z;

   costhe = cosdeg (theta);
   sinthe = sindeg (theta);

   dphi = phi - eul[2];
   cosphi = cosdeg (dphi);
   sinphi = sindeg (dphi);

   /* Compute the celestial longitude. */
   x = sinthe*eul[4] - costhe*eul[3]*cosphi;
   if (fabs(x) < tol) {
      /* Rearrange formula to reduce roundoff errors. */
      x = -cosdeg (theta+eul[1]) + costhe*eul[3]*(1.0 - cosphi);
   }
   y = -costhe*sinphi;
   if (x != 0.0 || y != 0.0) {
      dlng = atan2deg (y, x);
   } else {
      /* Change of origin of longitude. */
      dlng = dphi + 180.0;
   }
   *lng = eul[0] + dlng;

   /* Normalize the celestial longitude. */
   if (eul[0] >= 0.0) {
      if (*lng < 0.0) *lng += 360.0;
   } else {
      if (*lng > 0.0) *lng -= 360.0;
   }

   if (*lng > 360.0) {
      *lng -= 360.0;
   } else if (*lng < -360.0) {
      *lng += 360.0;
   }

   /* Compute the celestial latitude. */
   if (fmod(dphi,180.0) == 0.0) {
      *lat = theta + cosphi*eul[1];
      if (*lat >  90.0) *lat =  180.0 - *lat;
      if (*lat < -90.0) *lat = -180.0 - *lat;
   } else {
      z = sinthe*eul[3] + costhe*eul[4]*cosphi;

      /* Use an alternative formula for greater numerical accuracy. */
      if (fabs(z) > 0.99) {
	 if (z < 0)
            *lat = -acosdeg (sqrt(x*x+y*y));
	 else
            *lat =  acosdeg (sqrt(x*x+y*y));
      } else {
         *lat = asindeg (z);
      }
   }

	return 0;

}//close sphrev()



int WCSUtils::tnxpix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix)
{
    int	ira, idec, niter;
    double ra, dec, cosdec, sindec, cosra, sinra, x, y, phi, theta;
    double s, r, dphi, z, dpi, dhalfpi, twopi, tx;
    double xm, ym, f, fx, fy, g, gx, gy, denom, dx, dy;
    double colatp, coslatp, sinlatp, longp, sphtol;
    //double wf_gseval(), wf_gsder();

    /* get the axis numbers */
    if (wcs->coorflip) {
	ira = 1;
	idec = 0;
	}
    else {
	ira = 0;
	idec = 1;
	}

    /*  Compute the transformation from celestial coordinates ra and
	dec to native coordinates phi and theta. this is the spherical
	geometry part of the transformation */

    ra  = degrad_i (xpos - wcs->crval[ira]);
    dec = degrad_i (ypos);
    cosra = cos (ra);
    sinra = sin (ra);
    cosdec = cos (dec);
    sindec = sin (dec);
    colatp = degrad_i (90.0 - wcs->crval[idec]);
    coslatp = cos (colatp);
    sinlatp = sin (colatp);
    if (wcs->longpole == 999.0)
	longp = degrad_i (180.0);
    else
	longp = degrad_i(wcs->longpole);
    dpi = kPI;
    dhalfpi = dpi * 0.5;
    twopi = kPI + kPI;
    sphtol = kSPHTOL;

    /* Compute phi */
    x = sindec * sinlatp - cosdec * coslatp * cosra;
    if (fabs(x) < sphtol)
	x = -cos (dec + colatp) + cosdec * coslatp * (1.0 - cosra);
    y = -cosdec * sinra;
    if (x != 0.0 || y != 0.0)
	dphi = atan2 (y, x);
    else
	dphi = ra - dpi;
    phi = longp + dphi;
    if (phi > dpi)
	phi = phi - twopi;
    else if (phi < -dpi)
	phi = phi + twopi;

    /* Compute theta */
    if (fmod (ra, dpi) == 0.0) {
	theta = dec + cosra * colatp;
	if (theta > dhalfpi)
	    theta = dpi - theta;
	if (theta < -dhalfpi)
	    theta = -dpi - theta;
	}
    else {
	z = sindec * coslatp + cosdec * sinlatp * cosra;
	if (fabs (z) > 0.99) {
	    if (z >= 0.0)
		theta = acos (sqrt(x * x + y * y));
	    else
		theta = -acos (sqrt(x * x + y * y));
	    }
	else
	    theta = asin (z);
	}

    /*  Compute the transformation from native coordinates phi and theta
	to projected coordinates x and y */

    s = sin (theta);
    if (s == 0.0) {
	x = kBADCVAL;
	y = kBADCVAL;
	}
    else {
	r = wcs->rodeg * cos (theta) / s;
	if (wcs->lngcor == NULL && wcs->latcor == NULL) {
	    if (wcs->coorflip) {
		y  = r * sin (phi);
		x = -r * cos (phi);
		}
	    else {
		x  = r * sin (phi);
		y = -r * cos (phi);
		}
	    }
	else {
	    xm  = r * sin (phi);
	    ym = -r * cos (phi);
	    x = xm;
	    y = ym;
	    niter = 0;
	    while (niter < kMAX_NITER) {
		if (wcs->lngcor != NULL) {
		    f = x + wf_gseval (wcs->lngcor, x, y) - xm;
		    fx = wf_gsder (wcs->lngcor, x, y, 1, 0);
		    fx = 1.0 + fx;
		    fy = wf_gsder (wcs->lngcor, x, y, 0, 1);
		    }
		else {
		    f = x - xm;
		    fx = 1.0 ;
		    fy = 0.0;
		    }
		if (wcs->latcor != NULL) {
		    g = y + wf_gseval (wcs->latcor, x, y) - ym;
		    gx = wf_gsder (wcs->latcor, x, y, 1, 0);
		    gy = wf_gsder (wcs->latcor, x, y, 0, 1);
		    gy = 1.0 + gy;
		    }
		else {
		    g = y - ym;
		    gx = 0.0 ;
		    gy = 1.0;
		    }

		denom = fx * gy - fy * gx;
		if (denom == 0.0)
		    break;
		dx = (-f * gy + g * fy) / denom;
		dy = (-g * fx + f * gx) / denom;
		x = x + dx;
		y = y + dy;
		//if (MAX(MAX(fabs(dx),fabs(dy)),MAX(fabs(f),fabs(g))) < 2.80e-8)
		if (std::max(std::max(fabs(dx),fabs(dy)),std::max(fabs(f),fabs(g))) < 2.80e-8)
		    break;

		niter = niter + 1;
		}

	    /* Reverse x and y if axes flipped */
	    if (wcs->coorflip) {
		tx = x;
		x = y;
		y = tx;
		}
	    }
	}

    /* Scale and rotate using CD matrix */
    if (wcs->rotmat) {
	*xpix = x * wcs->dc[0] + y * wcs->dc[1];
	*ypix = x * wcs->dc[2] + y * wcs->dc[3];
	}

    else {

	/* Correct for rotation */
	if (wcs->rot!=0.0) {
	    double cosr = cos (degrad_i (wcs->rot));
	    double sinr = sin (degrad_i (wcs->rot));
	    *xpix = x * cosr + y * sinr;
	    *ypix = y * cosr - x * sinr;
	    }
	else {
	    *xpix = x;
	    *ypix = y;
	    }

	/* Scale using CDELT */
	if (wcs->xinc != 0.)
	    *xpix = *xpix / wcs->xinc;
	if (wcs->yinc != 0.)
	    *ypix = *ypix / wcs->yinc;
	}

    /* Convert to pixels  */
    *xpix = *xpix + wcs->xrefpix;
    *ypix = *ypix + wcs->yrefpix;

	return (0);

}//close tnxpix()

/* wf_gsder -- procedure to calculate a new surface which is a derivative of
 * the input surface.
 */

double WCSUtils::wf_gsder (IRAFsurface* sf1,double x,double y,int nxd,int nyd)
{
    int nxder, nyder, i, j, k, nbytes;
    int order, maxorder1, maxorder2, nmove1, nmove2;
    IRAFsurface *sf2 = 0;
    double *ptr1, *ptr2;
    double zfit, norm;
    //double wf_gseval();

		double *coeff = NULL;
		int nbcoeff = 0;

    if (sf1 == NULL)
	return (0.0);

    if (nxd < 0 || nyd < 0) {
	fprintf (stderr, "TNX_GSDER: order of derivatives cannot be < 0\n");
	return (0.0);
	}

    if (nxd == 0 && nyd == 0) {
	zfit = wf_gseval (sf1, x, y);
	return (zfit);
	}

    /* Allocate space for new surface */
    sf2 = (IRAFsurface *) malloc (sizeof (IRAFsurface));

    /* Check the order of the derivatives */
    nxder = std::min (nxd, sf1->xorder - 1);
    nyder = std::min (nyd, sf1->yorder - 1);

    /* Set up new surface */
    sf2->type = sf1->type;

    /* Set the derivative surface parameters */
    if (sf2->type == eTNX_LEGENDRE ||
	sf2->type == eTNX_CHEBYSHEV ||
	sf2->type == eTNX_POLYNOMIAL) {

	sf2->xterms = sf1->xterms;

	/* Find the order of the new surface */
	switch (sf2->xterms) {
	    case eTNX_XNONE: 
		if (nxder > 0 && nyder > 0) {
		    sf2->xorder = 1;
		    sf2->yorder = 1;
		    sf2->ncoeff = 1;
		    }
		else if (nxder > 0) {
		    sf2->xorder = std::max (1, sf1->xorder - nxder);
		    sf2->yorder = 1;
		    sf2->ncoeff = sf2->xorder;
		    }
		else if (nyder > 0) {
		    sf2->xorder = 1;
		    sf2->yorder = std::max (1, sf1->yorder - nyder);
		    sf2->ncoeff = sf2->yorder;
		    }
		break;

	    case eTNX_XHALF:
		maxorder1 = std::max (sf1->xorder+1, sf1->yorder+1);
		order = std::max(1, std::min(maxorder1-1-nyder-nxder,sf1->xorder-nxder));
		sf2->xorder = order;
		order = std::max(1, std::min(maxorder1-1-nyder-nxder,sf1->yorder-nyder));
		sf2->yorder = order;
		order = std::min (sf2->xorder, sf2->yorder);
		sf2->ncoeff = sf2->xorder * sf2->yorder - (order*(order-1)/2);
		break;

	    default:
		sf2->xorder = std::max (1, sf1->xorder - nxder);
		sf2->yorder = std::max (1, sf1->yorder - nyder);
		sf2->ncoeff = sf2->xorder * sf2->yorder;
	    }

	/* define the data limits */
	sf2->xrange = sf1->xrange;
	sf2->xmaxmin = sf1->xmaxmin;
	sf2->yrange = sf1->yrange;
	sf2->ymaxmin = sf1->ymaxmin;
	}

    else {
	fprintf (stderr, "TNX_GSDER: unknown surface type %d\n", sf2->type);
	return (0.0);
	}

    /* Allocate space for coefficients and basis functions */
    nbytes = sf2->ncoeff * sizeof(double);
    sf2->coeff = (double *) malloc (nbytes);
    nbytes = sf2->xorder * sizeof(double);
    sf2->xbasis = (double *) malloc (nbytes);
    nbytes = sf2->yorder * sizeof(double);
    sf2->ybasis = (double *) malloc (nbytes);

    /* Get coefficients */
    nbytes = sf1->ncoeff * sizeof(double);
    if (nbytes > nbcoeff) {
	if (nbcoeff > 0)
	    coeff = (double *) realloc (coeff, nbytes);
	else
	    coeff = (double *) malloc (nbytes);
	nbcoeff = nbytes;
	}
    (void) wf_gscoeff (sf1, coeff);

    /* Compute the new coefficients */
    switch (sf2->xterms) {
	case eTNX_XFULL:
	    ptr2 = sf2->coeff + (sf2->yorder - 1) * sf2->xorder;
	    ptr1 = coeff + (sf1->yorder - 1) * sf1->xorder;
	    for (i = sf1->yorder - 1; i >= nyder; i--) {
		for (j = i; j >= i-nyder+1; j--) {
		    for (k = 0; k < sf2->xorder; k++)
			ptr1[nxder+k] = ptr1[nxder+k] * (double)(j);
		    }
		for (j = sf1->xorder; j >= nxder+1; j--) {
		    for (k = j; k >= j-nxder+1; k--)
			ptr1[j-1] = ptr1[j-1] * (double)(k - 1);
		    }
		for (j = 0; j < sf2->xorder; j++)
		    ptr2[j] = ptr1[nxder+j];
		ptr2 = ptr2 - sf2->xorder;
		ptr1 = ptr1 - sf1->xorder;
		}
	    break;

	case eTNX_XHALF:
	    maxorder1 = std::max (sf1->xorder + 1, sf1->yorder + 1);
	    maxorder2 = std::max (sf2->xorder + 1, sf2->yorder + 1);
	    ptr2 = sf2->coeff + sf2->ncoeff;
	    ptr1 = coeff + sf1->ncoeff;
	    for (i = sf1->yorder; i >= nyder+1; i--) {
		nmove1 = std::max  (0, std::min (maxorder1 - i, sf1->xorder));
		nmove2 = std::max  (0, std::min (maxorder2 - i + nyder, sf2->xorder));
		ptr1 = ptr1 - nmove1;
		ptr2 = ptr2 - nmove2;
		for (j = i; j > i - nyder + 1; j--) {
		    for (k = 0; k < nmove2; k++)
			ptr1[nxder+k] = ptr1[nxder+k] * (double)(j-1);
		    }
		for (j = nmove1; j >= nxder+1; j--) {
		    for (k = j;  k >= j-nxder+1; k--)
			ptr1[j-1] = ptr1[j-1] * (double)(k - 1);
		    }
		for (j = 0; j < nmove2; j++)
		    ptr2[j] = ptr1[nxder+j];
		}
	    break;

	default:
	    if (nxder > 0 && nyder > 0)
		sf2->coeff[0] = 0.0;

	    else if (nxder > 0) { 
		ptr1 = coeff;
		ptr2 = sf2->coeff + sf2->ncoeff - 1;
		for (j = sf1->xorder; j >= nxder+1; j--) {
		    for (k = j; k >= j - nxder + 1; k--)
			ptr1[j-1] = ptr1[j-1] * (double)(k - 1);
		    ptr2[0] = ptr1[j-1];
		    ptr2 = ptr2 - 1;
		    }
		}

	    else if (nyder > 0) {
		ptr1 = coeff + sf1->ncoeff - 1;
		ptr2 = sf2->coeff;
		for (i = sf1->yorder; i >= nyder + 1; i--) {
		    for (j = i; j >= i - nyder + 1; j--)
			*ptr1 = *ptr1 * (double)(j - 1);
		    ptr1 = ptr1 - 1;
		    }
		for (i = 0; i < sf2->ncoeff; i++)
		    ptr2[i] = ptr1[i+1];
		}
	}

    /* evaluate the derivatives */
    zfit = wf_gseval (sf2, x, y);

    /* normalize */
    if (sf2->type != eTNX_POLYNOMIAL) { 
	norm = pow (sf2->xrange, (double)nxder) *
	       pow (sf2->yrange, (double)nyder);
	zfit = norm * zfit;
	}

    /* free the space */
    wf_gsclose (sf2);

	return (zfit);

}//close wf_gsder()



int WCSUtils::wf_gscoeff (IRAFsurface* sf,double* coeff)
{
	int ncoeff;		/* the number of coefficients */
  int i;

  /* Exctract coefficients from data structure and calculate their number */
  ncoeff = sf->ncoeff;
  for (i = 0; i < ncoeff; i++)
		coeff[i] = sf->coeff[i];
    
	return (ncoeff);

}//close wf_gscoeff()


void WCSUtils::wf_gsclose (IRAFsurface* sf)
{
  if (sf != NULL) {
	if (sf->xbasis != NULL)
	    free (sf->xbasis);
	if (sf->ybasis != NULL)
	    free (sf->ybasis);
	if (sf->coeff != NULL)
	    free (sf->coeff);
	free (sf);
	}
    
	return;

}//close wf_gsclose()


/* zpxpix -- inverse transform (world to physical) for the zenithal
 * azimuthal polynomial projection.
 */

int WCSUtils::zpxpix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix)
{
    int	i, ira, idec, niter;
    double ra, dec, cosdec, sindec, cosra, sinra, x, y, phi, theta;
    double s, r, dphi, z, dpi, dhalfpi, twopi, tx;
    double xm, ym, f, fx, fy, g, gx, gy, denom, dx, dy;
    double colatp, coslatp, sinlatp, longp, sphtol;
    //double wf_gseval(), wf_gsder();

    /* get the axis numbers */
    if (wcs->coorflip) {
	ira = 1;
	idec = 0;
	}
    else {
	ira = 0;
	idec = 1;
	}

    /*  Compute the transformation from celestial coordinates ra and
	dec to native coordinates phi and theta. this is the spherical
	geometry part of the transformation */

    ra  = degrad_i (xpos - wcs->crval[ira]);
    dec = degrad_i (ypos);
    cosra = cos (ra);
    sinra = sin (ra);
    cosdec = cos (dec);
    sindec = sin (dec);
    colatp = degrad_i (90.0 - wcs->crval[idec]);
    coslatp = cos (colatp);
    sinlatp = sin (colatp);
    if (wcs->longpole == 999.0)
	longp = degrad_i (180.0);
    else
	longp = degrad_i(wcs->longpole);
    dpi = kPI;
    dhalfpi = dpi * 0.5;
    twopi = kPI + kPI;
    sphtol = kSPHTOL;

    /* Compute phi */
    x = sindec * sinlatp - cosdec * coslatp * cosra;
    if (fabs(x) < sphtol)
	x = -cos (dec + colatp) + cosdec * coslatp * (1.0 - cosra);
    y = -cosdec * sinra;
    if (x != 0.0 || y != 0.0)
	dphi = atan2 (y, x);
    else
	dphi = ra - dpi;
    phi = longp + dphi;
    if (phi > dpi)
	phi = phi - twopi;
    else if (phi < -dpi)
	phi = phi + twopi;

    /* Compute theta */
    if (fmod (ra, dpi) == 0.0) {
	theta = dec + cosra * colatp;
	if (theta > dhalfpi)
	    theta = dpi - theta;
	if (theta < -dhalfpi)
	    theta = -dpi - theta;
	}
    else {
	z = sindec * coslatp + cosdec * sinlatp * cosra;
	if (fabs (z) > 0.99) {
	    if (z >= 0.0)
		theta = acos (sqrt(x * x + y * y));
	    else
		theta = -acos (sqrt(x * x + y * y));
	    }
	else
	    theta = asin (z);
	}

    /*  Compute the transformation from native coordinates phi and theta
	to projected coordinates x and y */

    s = dhalfpi - theta;
    r = 0.;
    for (i=9; i>=0; i--)
        r = r * s + wcs->prj.p[i];
    r = wcs->rodeg * r;

    if (wcs->lngcor == NULL && wcs->latcor == NULL) {
	if (wcs->coorflip) {
	    y  = r * sin (phi);
	    x = -r * cos (phi);
	} else {
	    x  = r * sin (phi);
	    y = -r * cos (phi);
	}
    } else {
	xm  = r * sin (phi);
	ym = -r * cos (phi);
	x = xm;
	y = ym;
	niter = 0;
	while (niter < kMAX_NITER) {
	    if (wcs->lngcor != NULL) {
		f = x + wf_gseval (wcs->lngcor, x, y) - xm;
		fx = wf_gsder (wcs->lngcor, x, y, 1, 0);
		fx = 1.0 + fx;
		fy = wf_gsder (wcs->lngcor, x, y, 0, 1);
		}
	    else {
		f = x - xm;
		fx = 1.0 ;
		fy = 0.0;
		}
	    if (wcs->latcor != NULL) {
		g = y + wf_gseval (wcs->latcor, x, y) - ym;
		gx = wf_gsder (wcs->latcor, x, y, 1, 0);
		gy = wf_gsder (wcs->latcor, x, y, 0, 1);
		gy = 1.0 + gy;
		}
	    else {
		g = y - ym;
		gx = 0.0 ;
		gy = 1.0;
		}

	    denom = fx * gy - fy * gx;
	    if (denom == 0.0)
		break;
	    dx = (-f * gy + g * fy) / denom;
	    dy = (-g * fx + f * gx) / denom;
	    x = x + dx;
	    y = y + dy;
	    if (std::max(std::max(fabs(dx),fabs(dy)),std::max(fabs(f),fabs(g))) < 2.80e-8)
		break;

	    niter = niter + 1;
	}

	/* Reverse x and y if axes flipped */
	if (wcs->coorflip) {
	    tx = x;
	    x = y;
	    y = tx;
	}
    }

    /* Scale and rotate using CD matrix */
    if (wcs->rotmat) {
	*xpix = x * wcs->dc[0] + y * wcs->dc[1];
	*ypix = x * wcs->dc[2] + y * wcs->dc[3];
	}

    else {

	/* Correct for rotation */
	if (wcs->rot!=0.0) {
	    double cosr = cos (degrad_i (wcs->rot));
	    double sinr = sin (degrad_i (wcs->rot));
	    *xpix = x * cosr + y * sinr;
	    *ypix = y * cosr - x * sinr;
	    }
	else {
	    *xpix = x;
	    *ypix = y;
	    }

	/* Scale using CDELT */
	if (wcs->xinc != 0.)
	    *xpix = *xpix / wcs->xinc;
	if (wcs->yinc != 0.)
	    *ypix = *ypix / wcs->yinc;
	}

    /* Convert to pixels  */
    *xpix = *xpix + wcs->xrefpix;
    *ypix = *ypix + wcs->yrefpix;

    return (0);
}//close 


int WCSUtils::worldpix (double xpos,double ypos,WCS* wcs,double* xpix,double* ypix)
{
  double dx, dy, ra0, dec0, ra, dec, coss, sins, dt, da, dd, sint;
  double l, m, geo1, geo2, geo3, sinr, cosr, tx, x, a2, a3, a4;
  double rthea,gamby2,a,b,c,phi,an,rap,v,tthea,co1,co2,co3,co4,ansq; /* COE */
  double cond2r=1.745329252e-2, deps=1.0e-5, twopi=6.28318530717959;

/* Structure elements */
  double xref;		/* x reference coordinate value (deg) */
  double yref;		/* y reference coordinate value (deg) */
  double xrefpix;	/* x reference pixel */
  double yrefpix;	/* y reference pixel */
  double xinc;		/* x coordinate increment (deg) */
  double yinc;		/* y coordinate increment (deg) */
  double rot;		/* Optical axis rotation (deg)  (from N through E) */
  int itype;

  /* Set local projection parameters */
  xref = wcs->xref;
  yref = wcs->yref;
  xrefpix = wcs->xrefpix;
  yrefpix = wcs->yrefpix;
  xinc = wcs->xinc;
  yinc = wcs->yinc;
  rot = degrad_i (wcs->rot);
  cosr = cos (rot);
  sinr = sin (rot);

  /* Projection type */
  itype = wcs->prjcode;

  /* Nonlinear position */
  if (itype > 0) {
    if (wcs->coorflip) {
      dec0 = degrad_i (xref);
      ra0 = degrad_i (yref);
      dt = xpos - yref;
      }
    else {
      ra0 = degrad_i (xref);
      dec0 = degrad_i (yref);
      dt = xpos - xref;
      }

    /* 0h wrap-around tests added by D.Wells 10/12/1994: */
    /* Modified to exclude weird reference pixels by D.Mink 2/3/2004 */
    if (xrefpix*xinc > 180.0 || xrefpix*xinc < -180.0) {
	if (dt > 360.0) xpos -= 360.0;
	if (dt < 0.0) xpos += 360.0;
	}
    else {
	if (dt > 180.0) xpos -= 360.0;
	if (dt < -180.0) xpos += 360.0;
	}
    /* NOTE: changing input argument xpos is OK (call-by-value in C!) */

    ra = degrad_i (xpos);
    dec = degrad_i (ypos);

    /* Compute direction cosine */
    coss = cos (dec);
    sins = sin (dec);
    l = sin(ra-ra0) * coss;
    sint = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
    }
  else {
    l = 0.0;
    sint = 0.0;
    sins = 0.0;
    coss = 0.0;
    ra = 0.0;
    dec = 0.0;
    ra0 = 0.0;
    dec0 = 0.0;
    m = 0.0;
    }

  /* Process by case  */
  switch (itype) {

    case eWCS_CAR:   /* -CAR Cartesian */
      l = ra - ra0;
      m = dec - dec0;
      break;

    case eWCS_SIN:   /* -SIN sin*/ 
	if (sint<0.0) return 1;
	m = sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0);
	break;

    case eWCS_TNX:   /* -TNX tan with polynomial correction */
    case eWCS_TPV:   /* -TPV tan with polynomial correction */
    case eWCS_ZPX:   /* -ZPX zpn with polynomial correction */
    case eWCS_TAN:   /* -TAN tan */
	if (sint<=0.0) return 1;
 	m = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
	l = l / m;
	m = (sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0)) / m;
	break;

    case eWCS_ARC:   /* -ARC Arc*/
	m = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
	if (m<-1.0) m = -1.0;
	if (m>1.0) m = 1.0;
	m = acos (m);
	if (m!=0) 
	    m = m / sin(m);
	else
	    m = 1.0;
	l = l * m;
	m = (sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0)) * m;
	break;

    case eWCS_NCP:   /* -NCP North celestial pole*/
	if (dec0==0.0) 
	    return 1;  /* can't stand the equator */
	else
	    m = (cos(dec0) - coss * cos(ra-ra0)) / sin(dec0);
	break;

    case eWCS_GLS:   /* -GLS global sinusoid */
    case eWCS_SFL:   /* -SFL Samson-Flamsteed */
	dt = ra - ra0;
	if (fabs(dec)>twopi/4.0) return 1;
	if (fabs(dec0)>twopi/4.0) return 1;
	m = dec - dec0;
	l = dt * coss;
	break;

    case eWCS_MER:   /* -MER mercator*/
	dt = yinc * cosr + xinc * sinr;
	if (dt==0.0) dt = 1.0;
	dy = degrad_i (yref/2.0 + 45.0);
	dx = dy + dt / 2.0 * cond2r;
	dy = log (tan (dy));
	dx = log (tan (dx));
	geo2 = degrad_i (dt) / (dx - dy);
	geo3 = geo2 * dy;
	geo1 = cos (degrad_i (yref));
	if (geo1<=0.0) geo1 = 1.0;
	dt = ra - ra0;
	l = geo1 * dt;
	dt = dec / 2.0 + twopi / 8.0;
	dt = tan (dt);
	if (dt<deps) return 2;
	m = geo2 * log (dt) - geo3;
	break;

    case eWCS_AIT:   /* -AIT Aitoff*/
	l = 0.0;
	m = 0.0;
	da = (ra - ra0) / 2.0;
	if (fabs(da)>twopi/4.0) return 1;
	dt = yinc*cosr + xinc*sinr;
	if (dt==0.0) dt = 1.0;
	dt = degrad_i (dt);
	dy = degrad_i (yref);
	dx = sin(dy+dt)/sqrt((1.0+cos(dy+dt))/2.0) -
	     sin(dy)/sqrt((1.0+cos(dy))/2.0);
	if (dx==0.0) dx = 1.0;
	geo2 = dt / dx;
	dt = xinc*cosr - yinc* sinr;
	if (dt==0.0) dt = 1.0;
	dt = degrad_i (dt);
	dx = 2.0 * cos(dy) * sin(dt/2.0);
	if (dx==0.0) dx = 1.0;
	geo1 = dt * sqrt((1.0+cos(dy)*cos(dt/2.0))/2.0) / dx;
	geo3 = geo2 * sin(dy) / sqrt((1.0+cos(dy))/2.0);
	dt = sqrt ((1.0 + cos(dec) * cos(da))/2.0);
	if (fabs(dt)<deps) return 3;
	l = 2.0 * geo1 * cos(dec) * sin(da) / dt;
	m = geo2 * sin(dec) / dt - geo3;
	break;

    case eWCS_STG:   /* -STG Sterographic*/
	da = ra - ra0;
	if (fabs(dec)>twopi/4.0) return 1;
	dd = 1.0 + sins * sin(dec0) + coss * cos(dec0) * cos(da);
	if (fabs(dd)<deps) return 1;
	dd = 2.0 / dd;
	l = l * dd;
	m = dd * (sins * cos(dec0) - coss * sin(dec0) * cos(da));
	break;

    case eWCS_COE:    /* allan: -COE projection added, AW, ESO*/
	gamby2 = sin (dec0);
	tthea = tan (dec0);
	rthea = 1. / tthea;
	a = -2. * tthea;
	b = tthea * tthea;
	c = tthea / 3.;
	a2 = a * a;
	a3 = a2 * a;
	a4 = a2 * a2;
	co1 = a/2.;
	co2 = -0.125 * a2 + b/2.;
	co3 = -0.25 * a*b + 0.0625 * a3 + c/2.0;
	co4 = -0.125 * b*b - 0.25 * a*c + 0.1875 * b*a2 - (5.0/128.0)*a4;
	phi = ra0 - ra;
	an = phi * gamby2;
	v = dec - dec0;
	rap = rthea * (1.0 + v * (co1+v * (co2+v * (co3+v * co4))));
	ansq = an * an;
	if (wcs->rotmat)
	    l = rap * an * (1.0 - ansq/6.0) * (wcs->cd[0] / fabs(wcs->cd[0]));
	else
	    l = rap * an * (1.0 - ansq/6.0) * (xinc / fabs(xinc));
	m = rthea - (rap * (1.0 - ansq/2.0));
	break;

    }  /* end of itype switch */

  /* Convert back to degrees  */
  if (itype > 0) {
    dx = raddeg_i (l);
    dy = raddeg_i (m);
    }

  /* For linear or pixel projection */
  else {
    dx = xpos - xref;
    dy = ypos - yref;
    }

  if (wcs->coorflip) {
    tx = dx;
    dx = dy;
    dy = tx;
    }

  /* Scale and rotate using CD matrix */
  if (wcs->rotmat) {
    tx = dx * wcs->dc[0] + dy * wcs->dc[1];
    dy = dx * wcs->dc[2] + dy * wcs->dc[3];
    dx = tx;
    }

  /* Scale and rotate using CDELTn and CROTA2 */
  else {

    /* Correct for rotation */
    if (rot!=0.0) {
      tx = dx*cosr + dy*sinr;
      dy = dy*cosr - dx*sinr;
      dx = tx;
      }

    /* Scale using CDELT */
    if (xinc != 0.)
      dx = dx / xinc;
    if (yinc != 0.)
      dy = dy / yinc;
    }

  /* Convert to pixels  */
  *xpix = dx + xrefpix;
  if (itype == eWCS_CAR) {
    if (*xpix > wcs->nxpix) {
      x = *xpix - (360.0 / xinc);
      if (x > 0.0) *xpix = x;
      }
    else if (*xpix < 0) {
      x = *xpix + (360.0 / xinc);
      if (x <= wcs->nxpix) *xpix = x;
      }
    }
  *ypix = dy + yrefpix;

  return 0;

}//close worldpix()

void WCSUtils::foc2pix(WCS* wcs,double x,double y,double* u,double* v)
{
    int m, n, i, j, k;
    double s[kDISTMAX], sum;
    double temp_x, temp_y;

    /* Spitzer distortion */
    if (wcs->distcode == eDISTORT_SIRTF) {
	m = wcs->distort.ap_order;
	n = wcs->distort.bp_order;

	temp_x = x - wcs->xrefpix;
	temp_y = y - wcs->yrefpix;

	/* compute u */
	for (j = 0; j <= m; j++) {
	    s[j] = wcs->distort.ap[m-j][j];
	    for (k = j-1; k >= 0; k--) {
	   	s[j] = (temp_y * s[j]) + wcs->distort.ap[m-j][k];
		}
	    }
  
	sum = s[0];
	for (i=m; i>=1; i--){
	    sum = (temp_x * sum) + s[m-i+1];
	    }
	*u = sum;

	/* compute v*/
	for (j = 0; j <= n; j++) {
	    s[j] = wcs->distort.bp[n-j][j];
	    for (k = j-1; k >= 0; k--) {
		s[j] = temp_y*s[j] + wcs->distort.bp[n-j][k];
		}
	    }
   
	sum = s[0];
	for (i = n; i >= 1; i--)
	    sum = temp_x * sum + s[n-i+1];

	*v = sum;

	*u = x + *u;
	*v = y + *v;
	}

    /* If no distortion, return pixel positions unchanged */
    else {
	*u = x;
	*v = y;
	}

  return;

}//close fox2pix()

void WCSUtils::deg2str(char* str,int lstr,double deg,int ndec)
{
    char degform[8];
    int field, ltstr;
    char tstring[64];
    double deg1;
    double dsgn;

    /* Keep angle between -180 and 360 degrees */
    deg1 = deg;
    if (deg1 < 0.0 ) {
	deg1 = -deg1;
	dsgn = -1.0;
	}
    else
	dsgn = 1.0;
    deg1 = fmod(deg1, 360.0);
    deg1 *= dsgn;
    if (deg1 <= -180.0)
	deg1 = deg1 + 360.0;

    /* Write angle to string, adding 4 digits to number of decimal places */
    field = ndec + 4;
    if (ndec > 0) {
	sprintf (degform, "%%%d.%df", field, ndec);
	sprintf (tstring, degform, deg1);
	}
    else {
	sprintf (degform, "%%%4d", field);
	sprintf (tstring, degform, (int)deg1);
	}

    /* Move formatted string to returned string */
    ltstr = (int) strlen (tstring);
    if (ltstr < lstr-1)
	strcpy (str, tstring);
    else {
	strncpy (str, tstring, lstr-1);
	str[lstr-1] = 0;
	}
    
	return;

}//close deg2str()

/* Write the right ascension ra in sexagesimal format into string*/

void WCSUtils::ra2str(char* str,int lstr,double ra,int ndec)
{
    double a,b;
    double seconds;
    char tstring[64];
    int hours;
    int minutes;
    int isec, ltstr;
    double dsgn;

    /* Keep RA between 0 and 360 */
    if (ra < 0.0 ) {
	ra = -ra;
	dsgn = -1.0;
	}
    else
	dsgn = 1.0;
    ra = fmod(ra, 360.0);
    ra *= dsgn;
    if (ra < 0.0)
	ra = ra + 360.0;

    a = ra / 15.0;

    /* Convert to hours */
    hours = (int) a;

    /* Compute minutes */
    b =  (a - (double)hours) * 60.0;
    minutes = (int) b;

    /* Compute seconds */
    seconds = (b - (double)minutes) * 60.0;

    if (ndec > 5) {
	if (seconds > 59.999999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%09.6f",hours,minutes,seconds);
	}
    else if (ndec > 4) {
	if (seconds > 59.99999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%08.5f",hours,minutes,seconds);
	}
    else if (ndec > 3) {
	if (seconds > 59.9999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%07.4f",hours,minutes,seconds);
	}
    else if (ndec > 2) {
	if (seconds > 59.999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%06.3f",hours,minutes,seconds);
	}
    else if (ndec > 1) {
	if (seconds > 59.99) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%05.2f",hours,minutes,seconds);
	}
    else if (ndec > 0) {
	if (seconds > 59.9) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%04.1f",hours,minutes,seconds);
	}
    else {
	isec = (int)(seconds + 0.5);
	if (isec > 59) {
	    isec = 0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    hours = hours + 1;
	    }
	hours = hours % 24;
	(void) sprintf (tstring,"%02d:%02d:%02d",hours,minutes,isec);
	}

    /* Move formatted string to returned string */
    ltstr = (int) strlen (tstring);
    if (ltstr < lstr-1)
	strcpy (str, tstring);
    else {
	strncpy (str, tstring, lstr-1);
	str[lstr-1] = 0;
	}
    
	return;

}//close ra2str()


void WCSUtils::dec2str(char* str,int lstr,double dec,int ndec)
{
    double a, b, dsgn, deg1;
    double seconds;
    char sign;
    int degrees;
    int minutes;
    int isec, ltstr;
    char tstring[64];

    /* Keep angle between -180 and 360 degrees */
    deg1 = dec;
    if (deg1 < 0.0 ) {
	deg1 = -deg1;
	dsgn = -1.0;
	}
    else
	dsgn = 1.0;
    deg1 = fmod(deg1, 360.0);
    deg1 *= dsgn;
    if (deg1 <= -180.0)
	deg1 = deg1 + 360.0;

    a = deg1;

    /* Set sign and do all the rest with a positive */
    if (a < 0) {
	sign = '-';
	a = -a;
	}
    else
	sign = '+';

    /* Convert to degrees */
    degrees = (int) a;

    /* Compute minutes */
    b =  (a - (double)degrees) * 60.0;
    minutes = (int) b;

    /* Compute seconds */
    seconds = (b - (double)minutes) * 60.0;

    if (ndec > 5) {
	if (seconds > 59.999999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%09.6f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 4) {
	if (seconds > 59.99999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%08.5f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 3) {
	if (seconds > 59.9999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%07.4f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 2) {
	if (seconds > 59.999) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%06.3f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 1) {
	if (seconds > 59.99) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%05.2f",sign,degrees,minutes,seconds);
	}
    else if (ndec > 0) {
	if (seconds > 59.9) {
	    seconds = 0.0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%04.1f",sign,degrees,minutes,seconds);
	}
    else {
	isec = (int)(seconds + 0.5);
	if (isec > 59) {
	    isec = 0;
	    minutes = minutes + 1;
	    }
	if (minutes > 59) {
	    minutes = 0;
	    degrees = degrees + 1;
	    }
	(void) sprintf (tstring,"%c%02d:%02d:%02d",sign,degrees,minutes,isec);
	}

    /* Move formatted string to returned string */
    ltstr = (int) strlen (tstring);
    if (ltstr < lstr-1)
	strcpy (str, tstring);
    else {
	strncpy (str, tstring, lstr-1);
	str[lstr-1] = 0;
	}

	return;

}//close dec2str()


void WCSUtils::num2str(char* string,double num,int field,int ndec)
{
    char numform[8];

    if (field > 0) {
	if (ndec > 0) {
	    sprintf (numform, "%%%d.%df", field, ndec);
	    sprintf (string, numform, num);
	    }
	else {
	    sprintf (numform, "%%%dd", field);
	    sprintf (string, numform, (int)num);
	    }
	}
    else {
	if (ndec > 0) {
	    sprintf (numform, "%%.%df", ndec);
	    sprintf (string, numform, num);
	    }
	else {
	    sprintf (string, "%d", (int)num);
	    }
	}
    
	return;

}//close num2str()


}//close namespace
