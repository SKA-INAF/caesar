
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
* @file ImgUtils.h
* @class ImgUtils
* @brief ImgUtils
*
* Image utility class
* @author S. Riggi
* @date 12/07/2018
*/

#ifndef _IMG_UTILS_h
#define _IMG_UTILS_h 1

#include <TObject.h>

namespace Caesar{

//Forward declarations
class Image;
class Source;

class ImgUtils : public TObject {

	public:
		/** 
		\brief Constructor
 		*/
		ImgUtils(){};
		
		/** 
		\brief Destructor
 		*/
		virtual ~ImgUtils(){};

	public:
			
		/** 
		\brief Get image with single circle level set
 		*/	
		static Image* GetCircleLevelSetImage(long int nX,long int nY,double f=0.1);

		/** 
		\brief Get image with checker board level set
 		*/
		static Image* GetCheckerBoardLevelSetImage(long int nX,long int nY,double f=0.1);

		
	ClassDef(ImgUtils,1)
		
};//close ImgUtils()

#ifdef __MAKECINT__
#pragma link C++ class ImgUtils+;
#endif

}//close namespace

#endif
