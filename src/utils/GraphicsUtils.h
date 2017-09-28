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
* @file GraphicsUtils.h
* @class GraphicsUtils
* @brief Utility functions for graphics tasks
*
* Utility functions for graphics tasks
* @author S. Riggi
* @date 15/01/2016
*/


#ifndef _GRAPHICS_UTILS_h
#define _GRAPHICS_UTILS_h 1

#include <TObject.h>
#include <TGaxis.h>
#include <TPolyLine.h>

namespace Caesar {

class Image;


class GraphicsUtils : public TObject {

  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    GraphicsUtils();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~GraphicsUtils();

		
	public:

		/**
		* \brief Get palette code
		*/
		static void SetPalette(int paletteStyle,int ncolors=999);

		/**
		* \brief Set thermal palette
		*/
		static int SetThermalPalette(int ncolors=999);
		/**
		* \brief Set hot-to-cold palette
		*/
		static int SetHotColdPalette(int ncolors=999);
		/**
		* \brief Set cold-to-hot palette
		*/
		static int SetColdHotPalette(int ncolors=999);
		/**
		* \brief Setblack & white palette
		*/
		static int SetBWPalette(int ncolors=999);


		
		//=========================================
		//==  NEW IMAGE METHODS 
		//=========================================
		
		/**
		* \brief Set WCS axis
		*/
		static int SetWCSAxis(Image* img,TGaxis& xaxis,TGaxis& yaxis,int coordSystem=-1);

		/**
		* \brief Set WCS proj grid
		*/
		static int SetWCSProjGrid(Image* img,std::vector<TPolyLine>& gridx,std::vector<TPolyLine>& gridy,int coordSystem);
		/**
		* \brief Retrieve image from pad
		*/
		static Image* FindImageFromPad();
		/**
		* \brief Update pad
		*/
		static int PadUpdater();
		/**
		* \brief Update ROOT canvas gaxis
		*/
		static int UpdateGAxis();


	private:
	
		ClassDef(GraphicsUtils,1)
};

#ifdef __MAKECINT__
#pragma link C++ class GraphicsUtils+;
#endif	

}//close namespace


#endif 
