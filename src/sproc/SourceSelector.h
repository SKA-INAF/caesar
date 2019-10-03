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
* @file SourceSelector.h
* @class SourceSelector
* @brief SourceSelector class
*
* Class to select sources using quality cuts provided in a cut file
* @author S. Riggi
* @date 21/06/2019
*/

#ifndef _SOURCE_SELECTOR_h
#define _SOURCE_SELECTOR_h 1


#include <Consts.h>

#include <TObject.h>
#include <TMatrixD.h>
#include <TTree.h>

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
#include <string>
#include <iostream>
#include <time.h>
#include <ctime>


namespace Caesar {

class Source;
class WCS;
class Cut;


class SourceSelector : public TObject 
{
  public:
		
		/** 
		\brief Class constructor: initialize structures.
 		*/
    SourceSelector();
		/**
		* \brief Class destructor: free allocated memory
		*/
   	virtual ~SourceSelector();

		
		typedef std::function<bool(Source*,Cut*)> CutFcn;
		typedef std::map<std::string,CutFcn> CutFcnRegistry;

	public:
		/**
		* \brief Select sources by applying quality cuts provided in cut file
		*/
		static int SelectSources(std::vector<Source*>& sources_sel,const std::vector<Source*>& sources,std::string cutFile);

	protected:

		/**
		* \brief Apply cut to source
		*/
		static bool ApplyCut(std::string cutName,Source* source,Cut* cut);

		//=====================================
		//==   SOURCE CUT DEFINITIONS
		//=====================================
		/**
		* \brief Select source by fit quality
		*/
		static bool HasFitCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source by fit converged 
		*/
		static bool SourceFitStatusCut(Source* source,Cut* cut);
		/**
		* \brief Select source by fit quality flag 
		*/
		static bool SourceFitQualityFlagCut(Source* source,Cut* cut);
		/**
		* \brief Select source by fit quality
		*/
		static bool HasGoodFitChi2Cut(Source* aSource,Cut* cut);
		/**
		* \brief Select source by flux
		*/
		static bool SourceFluxCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source by flux over island flux ratio
		*/
		static bool SourceFluxToIslandRatioCut(Source* aSource,Cut* cut);

		/**
		* \brief Select source fitted components by flux
		*/
		static bool SourceComponentFluxCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by centroid 
		*/
		static bool SourceComponentCentroidInsideIslandCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by peak amplitude flux 
		*/
		static bool SourceComponentPeakFluxCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by peak amplitude flux 
		*/
		static bool SourceComponentPeakSignificanceCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by type
		*/
		static bool SourceComponentTypeCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by flag
		*/
		static bool SourceComponentFlagCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by centroid distance
		*/
		static bool SourceComponentCentroidDistanceCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by fit ellipse/beam eccentricity ratio
		*/
		static bool SourceComponentEccentricityRatioCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source fitted components by fit ellipse/beam area ratio
		*/
		static bool SourceComponentBeamAreaRatioCut(Source* aSource,Cut* cut);


		/**
		* \brief Select source by number of pixels
		*/
		static bool NPixelsCut(Source* aSource,Cut* cut);
		/**
		* \brief Select source by min bounding box
		*/
		static bool MinBoundingBoxCut(Source* source,Cut* cut);
		/**
		* \brief Select source by good source flag
		*/
		static bool GoodSourceFlagCut(Source* source,Cut* cut);
		/**
		* \brief Select source by contour circularity ratio
		*/
		static bool CircularityRatioCut(Source* source,Cut* cut);	
		/**
		* \brief Select source by contour elongation
		*/
		static bool ElongationCut(Source* source,Cut* cut);
		/**
		* \brief Select source by ellipse area ratio
		*/
		static bool EllipseAreaRatioCut(Source* source,Cut* cut);
		/**
		* \brief Select source by beam area ratio cut
		*/
		static bool BeamAreaRatioCut(Source* source,Cut* cut);
		/**
		* \brief Select source by type flag
		*/
		static bool SourceTypeCut(Source* source,Cut* cut);
		/**
		* \brief Select source by sourceness flag
		*/
		static bool SourceFlagCut(Source* source,Cut* cut);
		/**
		* \brief Select source by source sim type
		*/
		static bool SourceSimTypeCut(Source* source,Cut* cut);

		
	
	private:
		/**
		* \brief Selection cut registry (e.g. list of defined source cuts)
		*/
		static CutFcnRegistry m_cutFcnRegistry; 
	
	ClassDef(SourceSelector,1)

};//close class

#ifdef __MAKECINT__
#pragma link C++ class SourceSelector+;
#endif

}//close namespace 


#endif
