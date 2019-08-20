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
* @file AstroObjectParser.cc
* @class AstroObjectParser
* @brief AstroObjectParser class
*
* AstroObjectParser class
* @author S. Riggi
* @date 06/08/2019
*/
#include <AstroObjectParser.h>
#include <AstroObject.h>
#include <MathUtils.h>
#include <CodeUtils.h>


#ifdef LOGGING_ENABLED
	#include <Logger.h>
#endif

#include <Consts.h>

//ROOT headers
#include <TMath.h>

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
#include <deque>

using namespace std;


ClassImp(Caesar::AstroObjectParser)

namespace Caesar {

/**
* \brief Init SIMBAD object identifier code map
*/
std::map<std::string,int> InitSimbadObjIdMap() {

	//Fill map items
	std::map<std::string,int> tmp
	{
  	std::pair<std::string,int> ("?",eUNKNOWN_OBJECT),

		//Radio source object
    std::pair<std::string,int> ("Rad",eRADIO_OBJ),
		std::pair<std::string,int> ("mR",eRADIO_OBJ),
		std::pair<std::string,int> ("cm",eRADIO_OBJ),
		std::pair<std::string,int> ("mm",eRADIO_OBJ),
		std::pair<std::string,int> ("smm",eRADIO_OBJ),
		std::pair<std::string,int> ("HI",eRADIO_OBJ),
		std::pair<std::string,int> ("rB",eRADIO_OBJ),
		std::pair<std::string,int> ("Mas",eRADIO_OBJ),

		//UV source object
		std::pair<std::string,int> ("UV",eUV_OBJ),

		//Infrared source object
		std::pair<std::string,int> ("IR",eIR_OBJ),
    std::pair<std::string,int> ("FIR",eIR_OBJ),
    std::pair<std::string,int> ("NIR",eIR_OBJ),

		//X source object
    std::pair<std::string,int> ("X",eX_OBJ),
		std::pair<std::string,int> ("UX?",eX_OBJ),
		std::pair<std::string,int> ("ULX",eX_OBJ),

		//Gamma source object
    std::pair<std::string,int> ("gam",eGAMMA_OBJ),
   	std::pair<std::string,int> ("gB",eGAMMA_RAY_BURST),
		
		//Star objects
		std::pair<std::string,int> ("*",eSTAR),
		std::pair<std::string,int> ("**",eSTAR),
		std::pair<std::string,int> ("*iC",eSTAR),
		std::pair<std::string,int> ("*iN",eSTAR),
		std::pair<std::string,int> ("*iA",eSTAR),
		std::pair<std::string,int> ("*i*",eSTAR),
		std::pair<std::string,int> ("V*?",eSTAR),
		std::pair<std::string,int> ("Pe*",eSTAR),		
		std::pair<std::string,int> ("HB*",eSTAR),
		std::pair<std::string,int> ("Y*O",eSTAR),
		std::pair<std::string,int> ("Ae*",eSTAR),		
		std::pair<std::string,int> ("Em*",eSTAR),
		std::pair<std::string,int> ("Be*",eSTAR),
		std::pair<std::string,int> ("BS*",eSTAR),
		std::pair<std::string,int> ("RG*",eSTAR),
		std::pair<std::string,int> ("AB*",eSTAR),
		std::pair<std::string,int> ("C*",eSTAR),
		std::pair<std::string,int> ("S*",eSTAR),
		std::pair<std::string,int> ("sg*",eSTAR),
		std::pair<std::string,int> ("s*r",eSTAR),		
		std::pair<std::string,int> ("s*y",eSTAR),		
		std::pair<std::string,int> ("s*b",eSTAR),		
		std::pair<std::string,int> ("HS*",eSTAR),		
		std::pair<std::string,int> ("pA*",eSTAR),		
		std::pair<std::string,int> ("WD*",eSTAR),		
		std::pair<std::string,int> ("ZZ*",eSTAR),		
		std::pair<std::string,int> ("LM*",eSTAR),		
		std::pair<std::string,int> ("BD*",eSTAR),		
		std::pair<std::string,int> ("N*",eSTAR),		
		std::pair<std::string,int> ("OH*",eSTAR),		
		std::pair<std::string,int> ("CH*",eSTAR),		
		std::pair<std::string,int> ("pr*",eSTAR),		
		std::pair<std::string,int> ("TT*",eSTAR),		
		std::pair<std::string,int> ("WR*",eSTAR),		
		std::pair<std::string,int> ("PM*",eSTAR),		
		std::pair<std::string,int> ("HV*",eSTAR),		
		std::pair<std::string,int> ("V*",eSTAR),		
		std::pair<std::string,int> ("Ir*",eSTAR),		
		std::pair<std::string,int> ("Or*",eSTAR),		
		std::pair<std::string,int> ("RI*",eSTAR),		
		std::pair<std::string,int> ("Er*",eSTAR),		
		std::pair<std::string,int> ("Fl*",eSTAR),		
		std::pair<std::string,int> ("FU*",eSTAR),		
		std::pair<std::string,int> ("RC*",eSTAR),		
		std::pair<std::string,int> ("RC?",eSTAR),		
		std::pair<std::string,int> ("Ro*",eSTAR),		
		std::pair<std::string,int> ("a2*",eSTAR),		
		std::pair<std::string,int> ("Psr",eSTAR),		
		std::pair<std::string,int> ("BY*",eSTAR),		
		std::pair<std::string,int> ("RS*",eSTAR),		
		std::pair<std::string,int> ("Pu*",eSTAR),		
		std::pair<std::string,int> ("RR*",eSTAR),		
		std::pair<std::string,int> ("Ce*",eSTAR),	
		std::pair<std::string,int> ("dS*",eSTAR),		
		std::pair<std::string,int> ("RV*",eSTAR),			
		std::pair<std::string,int> ("WV*",eSTAR),		
		std::pair<std::string,int> ("bC*",eSTAR),		
		std::pair<std::string,int> ("cC*",eSTAR),		
		std::pair<std::string,int> ("gD*",eSTAR),		
		std::pair<std::string,int> ("SX*",eSTAR),		
		std::pair<std::string,int> ("LP*",eSTAR),		
		std::pair<std::string,int> ("Mi*",eSTAR),		
		std::pair<std::string,int> ("sr*",eSTAR),	
		std::pair<std::string,int> ("SN*",eSTAR),		
		std::pair<std::string,int> ("su*",eSTAR),		
		std::pair<std::string,int> ("s?r",eSTAR),
		std::pair<std::string,int> ("AB?",eSTAR),
		std::pair<std::string,int> ("LP?",eSTAR),
		std::pair<std::string,int> ("pA?",eSTAR),	
		std::pair<std::string,int> ("BD?",eSTAR),	

		//Galaxy objects
		std::pair<std::string,int> ("G",eGALAXY),
		std::pair<std::string,int> ("G?",eGALAXY),
		std::pair<std::string,int> ("PoG",eGALAXY),
		std::pair<std::string,int> ("GiC",eGALAXY),
		std::pair<std::string,int> ("BiC",eGALAXY),
		std::pair<std::string,int> ("GiG",eGALAXY),
		std::pair<std::string,int> ("GiP",eGALAXY),
		std::pair<std::string,int> ("HzG",eGALAXY),
		std::pair<std::string,int> ("rG",eGALAXY),
		std::pair<std::string,int> ("H2G",eGALAXY),
		std::pair<std::string,int> ("LSB",eGALAXY),
		std::pair<std::string,int> ("AG?",eGALAXY),
		std::pair<std::string,int> ("Q?",eGALAXY),
		std::pair<std::string,int> ("Bz?",eGALAXY),
		std::pair<std::string,int> ("BL?",eGALAXY),
		std::pair<std::string,int> ("EmG",eGALAXY),	
		std::pair<std::string,int> ("SBG",eGALAXY),
		std::pair<std::string,int> ("bCG",eGALAXY),
		std::pair<std::string,int> ("AGN",eGALAXY),
		std::pair<std::string,int> ("SyG",eGALAXY),
		std::pair<std::string,int> ("Sy1",eGALAXY),
		std::pair<std::string,int> ("Sy2",eGALAXY),
		std::pair<std::string,int> ("Bla",eGALAXY),
		std::pair<std::string,int> ("BLL",eGALAXY),
		std::pair<std::string,int> ("OVV",eGALAXY),
		std::pair<std::string,int> ("QSO",eGALAXY),
		
		//Other objects
		std::pair<std::string,int> ("HII",eHII),
		std::pair<std::string,int> ("PN",ePN),
		std::pair<std::string,int> ("PN?",ePN),
		std::pair<std::string,int> ("SR?",eSNR),
		std::pair<std::string,int> ("SNR",eSNR),
		std::pair<std::string,int> ("bub",eBUBBLE),
		std::pair<std::string,int> ("MoC",eMOLECULAR_CLOUD),
		std::pair<std::string,int> ("cor",eMOLECULAR_CLOUD),
		std::pair<std::string,int> ("Pl?",ePLANET),
		std::pair<std::string,int> ("Pl",ePLANET),

		std::pair<std::string,int> ("ClG",eGALAXY_CLUSTER),
		std::pair<std::string,int> ("GrG",eGALAXY_GROUP),
		std::pair<std::string,int> ("CGG",eGALAXY_GROUP),
		std::pair<std::string,int> ("PaG",eGALAXY_GROUP),
		std::pair<std::string,int> ("IG",eGALAXY_GROUP),
		std::pair<std::string,int> ("SCG",eGALAXY_SUPERCLUSTER),

		std::pair<std::string,int> ("Gl?",eGLOBULAR_CLUSTER),
		std::pair<std::string,int> ("GlC",eGLOBULAR_CLUSTER),
		std::pair<std::string,int> ("Cl*",eSTAR_CLUSTER),
		std::pair<std::string,int> ("OpC",eSTAR_CLUSTER),
		

		std::pair<std::string,int> ("XB*",eX_BINARY),
		std::pair<std::string,int> ("LXB",eX_BINARY),
		std::pair<std::string,int> ("HXB",eX_BINARY),
		std::pair<std::string,int> ("XB?",eX_BINARY),
		std::pair<std::string,int> ("LX?",eX_BINARY),
		std::pair<std::string,int> ("HX?",eX_BINARY),

		std::pair<std::string,int> ("GNe",eNEBULA),
		std::pair<std::string,int> ("BNe",eNEBULA),
		std::pair<std::string,int> ("DNe",eNEBULA),
		std::pair<std::string,int> ("RNe",eNEBULA),
		std::pair<std::string,int> ("SFR",eSTAR_FORMING_REGION),
		std::pair<std::string,int> ("Cld",eCLOUD),
    std::pair<std::string,int> ("HVC",eCLOUD),
		std::pair<std::string,int> ("PoC",eCLOUD),
		std::pair<std::string,int> ("Sy*",eSTAR_BINARY),
		std::pair<std::string,int> ("Sy?",eSTAR_BINARY),
		std::pair<std::string,int> ("sh",eHI_SHELL),
		std::pair<std::string,int> ("No",eNOVA),
		std::pair<std::string,int> ("No?",eNOVA),
		std::pair<std::string,int> ("ev",eTRANSIENT_EVENT),
		

	};
	return tmp;

}//close InitSimbadObjIdMap()
std::map<std::string,int> AstroObjectParser::m_simbadObjIdMap(InitSimbadObjIdMap()); 


/**
* \brief Init SIMBAD object minor identifier code map
*/
std::map<std::string,int> InitSimbadObjSubIdMap() {

	//Fill map items
	std::map<std::string,int> tmp
	{
  	std::pair<std::string,int> ("?",eUNKNOWN_OBJECT),

		//Radio source object
    std::pair<std::string,int> ("Rad",eRADIO_OBJ),
		std::pair<std::string,int> ("mR",eRADIO_OBJ),
		std::pair<std::string,int> ("cm",eRADIO_OBJ),
		std::pair<std::string,int> ("mm",eRADIO_OBJ),
		std::pair<std::string,int> ("smm",eRADIO_OBJ),
		std::pair<std::string,int> ("HI",eRADIO_OBJ),
		std::pair<std::string,int> ("rB",eRADIO_OBJ),
		std::pair<std::string,int> ("Mas",eRADIO_OBJ),

		//Infrared source object
		std::pair<std::string,int> ("IR",eIR_OBJ),
    std::pair<std::string,int> ("FIR",eIR_OBJ),
    std::pair<std::string,int> ("NIR",eIR_OBJ),

		//UV source object
		std::pair<std::string,int> ("UV",eUV_OBJ),

		//X source object
    std::pair<std::string,int> ("X",eX_OBJ),
		std::pair<std::string,int> ("UX?",eX_OBJ),
		std::pair<std::string,int> ("ULX",eX_OBJ),

		//Gamma source object
    std::pair<std::string,int> ("gam",eGAMMA_OBJ),
   	std::pair<std::string,int> ("gB",eGAMMA_RAY_BURST),
		
		//Star objects
		std::pair<std::string,int> ("*",eSTAR),
		std::pair<std::string,int> ("*iC",eSTAR_IN_CLUSTER),
		std::pair<std::string,int> ("*iN",eSTAR_IN_NEBULA),
		std::pair<std::string,int> ("*iA",eSTAR_IN_ASSOCIATION),
		std::pair<std::string,int> ("*i*",eSTAR_IN_DOUBLE_SYSTEM),
		std::pair<std::string,int> ("V*?",eSTAR_VARIABLE),
		std::pair<std::string,int> ("Pe*",eSTAR_PECULIAR),		
		std::pair<std::string,int> ("HB*",eSTAR_HB),
		std::pair<std::string,int> ("Y*O",eSTAR_YSO),
		std::pair<std::string,int> ("Ae*",eSTAR_HERBIG),		
		std::pair<std::string,int> ("Em*",eSTAR_EMISS_LINE),
		std::pair<std::string,int> ("Be*",eSTAR_BE),
		std::pair<std::string,int> ("BS*",eSTAR_BS),
		std::pair<std::string,int> ("RG*",eSTAR_RG),
		std::pair<std::string,int> ("AB*",eSTAR_AB),
		std::pair<std::string,int> ("C*",eSTAR_C),
		std::pair<std::string,int> ("S*",eSTAR_S),
		std::pair<std::string,int> ("sg*",eSTAR_ESG),
		std::pair<std::string,int> ("s*r",eSTAR_RSG),		
		std::pair<std::string,int> ("s*y",eSTAR_YSG),		
		std::pair<std::string,int> ("s*b",eSTAR_BSG),		
		std::pair<std::string,int> ("HS*",eSTAR_HSD),		
		std::pair<std::string,int> ("pA*",eSTAR_PAGB),		
		std::pair<std::string,int> ("WD*",eSTAR_WD),		
		std::pair<std::string,int> ("ZZ*",eSTAR_PWD),		
		std::pair<std::string,int> ("LM*",eSTAR_LM),		
		std::pair<std::string,int> ("BD*",eSTAR_BD),
		std::pair<std::string,int> ("BD?",eSTAR_BD),		
		std::pair<std::string,int> ("N*",eSTAR_NEUTRON),		
		std::pair<std::string,int> ("OH*",eSTAR_OH),		
		std::pair<std::string,int> ("CH*",eSTAR_CH),		
		std::pair<std::string,int> ("pr*",eSTAR_PMS),		
		std::pair<std::string,int> ("TT*",eSTAR_TTAU),		
		std::pair<std::string,int> ("WR*",eSTAR_WR),		
		std::pair<std::string,int> ("PM*",eSTAR_PM),		
		std::pair<std::string,int> ("HV*",eSTAR_HV),		
		std::pair<std::string,int> ("V*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Ir*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Or*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("RI*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Er*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Fl*",eSTAR_FLARE),		
		std::pair<std::string,int> ("FU*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("RC*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("RC?",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Ro*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("a2*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Psr",eSTAR_PULSAR),		
		std::pair<std::string,int> ("BY*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("RS*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Pu*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("RR*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Ce*",eSTAR_VARIABLE),	
		std::pair<std::string,int> ("dS*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("RV*",eSTAR_VARIABLE),			
		std::pair<std::string,int> ("WV*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("bC*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("cC*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("gD*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("SX*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("LP*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("Mi*",eSTAR_VARIABLE),		
		std::pair<std::string,int> ("sr*",eSTAR_VARIABLE),	
		std::pair<std::string,int> ("SN*",eSTAR_SUPERNOVA),		
		std::pair<std::string,int> ("su*",eSTAR_SUBSTELLAR),
		std::pair<std::string,int> ("s?r",eSTAR_RSG),
		std::pair<std::string,int> ("AB?",eSTAR_AB),	
		std::pair<std::string,int> ("LP?",eSTAR_VARIABLE),	
		std::pair<std::string,int> ("pA?",eSTAR_PAGB),
		

		//Galaxy objects
		std::pair<std::string,int> ("G",eGALAXY),
		std::pair<std::string,int> ("G?",eGALAXY),
		std::pair<std::string,int> ("PoG",eGALAXY),
		std::pair<std::string,int> ("GiC",eGALAXY),
		std::pair<std::string,int> ("BiC",eGALAXY),
		std::pair<std::string,int> ("GiG",eGALAXY),
		std::pair<std::string,int> ("GiP",eGALAXY),
		std::pair<std::string,int> ("HzG",eGALAXY),
		std::pair<std::string,int> ("rG",eGALAXY),
		std::pair<std::string,int> ("H2G",eGALAXY),
		std::pair<std::string,int> ("LSB",eGALAXY),
		std::pair<std::string,int> ("AG?",eGALAXY),
		std::pair<std::string,int> ("Q?",eGALAXY),
		std::pair<std::string,int> ("Bz?",eGALAXY),
		std::pair<std::string,int> ("BL?",eGALAXY),
		std::pair<std::string,int> ("EmG",eGALAXY),	
		std::pair<std::string,int> ("SBG",eGALAXY),
		std::pair<std::string,int> ("bCG",eGALAXY),
		std::pair<std::string,int> ("AGN",eGALAXY),
		std::pair<std::string,int> ("SyG",eGALAXY),
		std::pair<std::string,int> ("Sy1",eGALAXY),
		std::pair<std::string,int> ("Sy2",eGALAXY),
		std::pair<std::string,int> ("Bla",eGALAXY),
		std::pair<std::string,int> ("BLL",eGALAXY),
		std::pair<std::string,int> ("OVV",eGALAXY),
		std::pair<std::string,int> ("QSO",eGALAXY),
		
		//Other objects
		std::pair<std::string,int> ("HII",eHII),
		std::pair<std::string,int> ("PN",ePN),
		std::pair<std::string,int> ("PN?",ePN),
		std::pair<std::string,int> ("SR?",eSNR),
		std::pair<std::string,int> ("SNR",eSNR),
		std::pair<std::string,int> ("bub",eBUBBLE),
		std::pair<std::string,int> ("MoC",eMOLECULAR_CLOUD),
		std::pair<std::string,int> ("cor",eMOLECULAR_CLOUD),
		std::pair<std::string,int> ("Pl?",ePLANET),
		std::pair<std::string,int> ("Pl",ePLANET),

		std::pair<std::string,int> ("ClG",eGALAXY_CLUSTER),
		std::pair<std::string,int> ("GrG",eGALAXY_GROUP),
		std::pair<std::string,int> ("CGG",eGALAXY_GROUP),
		std::pair<std::string,int> ("PaG",eGALAXY_GROUP),
		std::pair<std::string,int> ("IG",eGALAXY_GROUP),
		std::pair<std::string,int> ("SCG",eGALAXY_SUPERCLUSTER),

		std::pair<std::string,int> ("Gl?",eGLOBULAR_CLUSTER),
		std::pair<std::string,int> ("GlC",eGLOBULAR_CLUSTER),
		std::pair<std::string,int> ("Cl*",eSTAR_CLUSTER),
		std::pair<std::string,int> ("OpC",eSTAR_CLUSTER),
		

		std::pair<std::string,int> ("XB*",eX_BINARY),
		std::pair<std::string,int> ("LXB",eX_BINARY),
		std::pair<std::string,int> ("HXB",eX_BINARY),
		std::pair<std::string,int> ("XB?",eX_BINARY),
		std::pair<std::string,int> ("LX?",eX_BINARY),
		std::pair<std::string,int> ("HX?",eX_BINARY),

		std::pair<std::string,int> ("GNe",eNEBULA),
		std::pair<std::string,int> ("BNe",eNEBULA),
		std::pair<std::string,int> ("DNe",eNEBULA),
		std::pair<std::string,int> ("RNe",eNEBULA),
		std::pair<std::string,int> ("SFR",eSTAR_FORMING_REGION),
		std::pair<std::string,int> ("Cld",eCLOUD),
    std::pair<std::string,int> ("HVC",eCLOUD),
		std::pair<std::string,int> ("PoC",eCLOUD),
		std::pair<std::string,int> ("Sy*",eSTAR_BINARY),
		std::pair<std::string,int> ("Sy?",eSTAR_BINARY),
		std::pair<std::string,int> ("sh",eHI_SHELL),
		std::pair<std::string,int> ("No",eNOVA),
		std::pair<std::string,int> ("No?",eNOVA),
		std::pair<std::string,int> ("ev",eTRANSIENT_EVENT),
		
	};
	return tmp;

}//close InitSimbadObjSubIdMap()

std::map<std::string,int> AstroObjectParser::m_simbadObjSubIdMap(InitSimbadObjSubIdMap()); 


//===========================
//==    AstroObject CLASS
//===========================
AstroObjectParser::AstroObjectParser()
{
	
}//close constructor


AstroObjectParser::~AstroObjectParser()
{

}//close destructor


int AstroObjectParser::ParseSimbadData(std::vector<AstroObject*>& astroObjects,std::string filename,char delimiter)
{
	// Init data
	astroObjects.clear();	

	// Read region file and get region text lines
	std::vector<std::string> raw_data;
	if(Read(raw_data,filename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read region data from file "<<filename<<"!");
		#endif
		return -1;
	}

	// Check lines
	if(raw_data.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty data list read from file "<<filename<<"!");
		#endif
		return -1;
	}
	

	// Extract astro objects from parsed lines
	bool parseErr= false;

	for(size_t i=0;i<raw_data.size();i++)
	{
		AstroObject* astroObject= new AstroObject;
		if(ParseSimbadObjectData(*astroObject,raw_data[i],delimiter)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to parse astro object from data line no. "<<i+1<<", stop parsing and return error!");
			#endif
			parseErr= true;
			break;
		}		

		//Append region to list
		astroObjects.push_back(astroObject);

	
	}//end loop data

	//Clear data in case of errors
	if(parseErr){
		CodeUtils::DeletePtrCollection<AstroObject>(astroObjects);
		return -1;
	}
	
	return 0;

}//close ParseSimbadData()

int AstroObjectParser::ParseSimbadObjectData(AstroObject& astroObject,std::string data,char delimiter)
{
	//============================
	//==    SIMBAD FORMAT
	//============================
	//- Field 0: Object number (not needed)
	//- Field 1: Object name
	//- Field 2: Object identifier
	//- Field 3: Object coordinates (RA & DEC)
	//=============================

	//Check input
	if(data==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty catalog line given, nothing to be parsed!");
		#endif
		return -1;
	}
	
	
	//Split line using delimiter
	std::vector<std::string> fields= CodeUtils::SplitStringOnPattern(data,delimiter);
	size_t nRequestedFields= 4;

	if(fields.empty() || fields.size()!=nRequestedFields){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid number of fields parsed ("<<fields.size()<<" found when "<<nRequestedFields<<" expected)!");
		#endif
		return -1;
	}

	#ifdef LOGGING_ENABLED
		DEBUG_LOG("data="<<data<<": name="<<fields[1]<<", id="<<fields[2]<<", coords="<<fields[3]);
	#endif

	//Parse index
	std::string index_str= fields[0];
	CodeUtils::RemovePatternInString(index_str,"\t");
	CodeUtils::StripBlankSpaces(index_str);
	astroObject.index= atol(index_str.c_str());

	//Parse object name
	std::string name_str= fields[1];
	CodeUtils::RemovePatternInString(name_str,"\t");
	astroObject.name= name_str;
	
	//Parse object identifier
	std::string id_str= fields[2];
	CodeUtils::StripBlankSpaces(id_str);
	
	std::map<std::string,int>::iterator it= m_simbadObjIdMap.find(id_str);
	int id= eUNKNOWN_OBJECT;
	if(!m_simbadObjIdMap.empty() && it!=m_simbadObjIdMap.end()){
		id= it->second;
	}
	
	//If unknown id, check for star (should contain *)
	if(id==eUNKNOWN_OBJECT){
		if(CodeUtils::HasPatternInString(id_str,"*")){
			id= eSTAR;
		}
	}

	astroObject.id_str= id_str;
	astroObject.id= id;

	//Parse object minor identifier
	std::string subid_str= fields[2];
	CodeUtils::StripBlankSpaces(subid_str);
	
	it= m_simbadObjSubIdMap.find(subid_str);
	int subid= eUNKNOWN_OBJECT;
	if(!m_simbadObjSubIdMap.empty() && it!=m_simbadObjSubIdMap.end()){
		subid= it->second;
	}
	
	//If unknown id, check for star (should contain *)
	if(subid==eUNKNOWN_OBJECT){
		if(CodeUtils::HasPatternInString(subid_str,"*")){
			subid= eSTAR;
		}
	}

	astroObject.subid= subid;


	//Check if candidate or not (candidate object should contain a ?)
	if(CodeUtils::HasPatternInString(id_str,"?")){
		astroObject.confirmed= false;
	}
	else{
		astroObject.confirmed= true;
	}

	//Parse coordinates
	std::string coords_str= fields[3];
	istringstream coords_ss(coords_str);
	std::string ra_str= "";
	std::string dec_str= "";
	coords_ss >> ra_str >> dec_str;
	CodeUtils::StripBlankSpaces(ra_str);
	CodeUtils::StripBlankSpaces(dec_str);
	
	astroObject.x= atof(ra_str.c_str());
	astroObject.y= atof(dec_str.c_str());
	
	
	return 0;

}//close ParseSimbadObjectData()



int AstroObjectParser::ParseMGPSData(std::vector<AstroObject*>& astroObjects,std::string filename,char delimiter)
{
	// Init data
	astroObjects.clear();	

	// Read region file and get region text lines
	std::vector<std::string> raw_data;
	if(Read(raw_data,filename)<0){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to read region data from file "<<filename<<"!");
		#endif
		return -1;
	}

	// Check lines
	if(raw_data.empty()){
		#ifdef LOGGING_ENABLED
			WARN_LOG("Empty data list read from file "<<filename<<"!");
		#endif
		return -1;
	}
	

	// Extract astro objects from parsed lines
	bool parseErr= false;

	for(size_t i=0;i<raw_data.size();i++)
	{
		AstroObject* astroObject= new AstroObject;
		if(ParseMGPSObjectData(*astroObject,raw_data[i],delimiter)<0){
			#ifdef LOGGING_ENABLED
				ERROR_LOG("Failed to parse astro object from data line no. "<<i+1<<", stop parsing and return error!");
			#endif
			parseErr= true;
			break;
		}		

		//Append region to list
		astroObjects.push_back(astroObject);

	}//end loop data

	//Clear data in case of errors
	if(parseErr){
		CodeUtils::DeletePtrCollection<AstroObject>(astroObjects);
		return -1;
	}
	
	return 0;

}//close ParseMGPSData()



int AstroObjectParser::ParseMGPSObjectData(AstroObject& astroObject,std::string data,char delimiter)
{
	//============================
	//==    SIMBAD FORMAT
	//============================
	//- Field 0: Object name
	//- Field 1: Ra
	//- Field 2: Dec
	//- Field 3: RaErr
	//- Field 4: DecErr
	//- Field 5: FluxDensity
	//- Field 6: FluxDensityErr
	//- Field 7: Flux
	//- Field 8: FluxErr
	//- Field 9-11: bmaj, bmin, pa
	//- Field 12-14: deconv bmaj, bmin, pa (NOT NEEDED)
	//=============================

	//Check input
	if(data==""){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty catalog line given, nothing to be parsed!");
		#endif
		return -1;
	}
	
	
	//Split line using delimiter
	std::vector<std::string> fields= CodeUtils::SplitStringOnPattern(data,delimiter);
	size_t nRequestedFields= 12;

	if(fields.empty() || fields.size()!=nRequestedFields){
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Empty or invalid number of fields parsed ("<<fields.size()<<" found when "<<nRequestedFields<<" expected)!");
		#endif
		return -1;
	}

	#ifdef LOGGING_ENABLED
		INFO_LOG("data="<<data<<": name="<<fields[0]<<", ra="<<fields[1]<<", dec="<<fields[2]);
	#endif

	
	//Parse object name
	std::string name_str= fields[0];
	CodeUtils::RemovePatternInString(name_str,"\t");
	astroObject.name= name_str;
	
	//Set object identifier to RADIO source 
	astroObject.id_str= "";
	astroObject.id= eRADIO_OBJ;
	astroObject.subid= eRADIO_OBJ;


	//Parse coordinates
	std::string ra_str= fields[1];
	CodeUtils::StripBlankSpaces(ra_str);
	astroObject.x= atof(ra_str.c_str());

	std::string dec_str= fields[2];
	CodeUtils::StripBlankSpaces(dec_str);
	astroObject.y= atof(dec_str.c_str());

	//Parse coordinate errors
	std::string raErr_str= fields[3];
	CodeUtils::StripBlankSpaces(raErr_str);
	astroObject.xerr= atof(raErr_str.c_str());

	std::string decErr_str= fields[4];
	CodeUtils::StripBlankSpaces(decErr_str);
	astroObject.yerr= atof(decErr_str.c_str());

	//Parse flux and errors
	std::string fluxDensity_str= fields[5];
	std::string fluxDensityErr_str= fields[6];
	std::string flux_str= fields[7];
	std::string fluxErr_str= fields[8];
	CodeUtils::StripBlankSpaces(fluxDensity_str);
	CodeUtils::StripBlankSpaces(fluxDensityErr_str);
	CodeUtils::StripBlankSpaces(flux_str);
	CodeUtils::StripBlankSpaces(fluxErr_str);
	astroObject.fluxDensity= atof(fluxDensity_str.c_str());
	astroObject.fluxDensityErr= atof(fluxDensityErr_str.c_str());
	astroObject.flux= atof(flux_str.c_str());
	astroObject.fluxErr= atof(fluxErr_str.c_str());
	astroObject.hasFluxInfo= true;

	//Convert flux from mJy to Jy
	astroObject.fluxDensity/= 1000.;
	astroObject.fluxDensityErr/= 1000.;
	astroObject.flux/= 1000.;
	astroObject.fluxErr/= 1000.;

	//Parse ellipse info
	std::string bmaj_str= fields[9];
	std::string bmin_str= fields[10];
	std::string pa_str= fields[11];
	CodeUtils::StripBlankSpaces(bmaj_str);
	CodeUtils::StripBlankSpaces(bmin_str);
	CodeUtils::StripBlankSpaces(pa_str);
	astroObject.bmaj= atof(bmaj_str.c_str());
	astroObject.bmin= atof(bmin_str.c_str());
	astroObject.pa= atof(pa_str.c_str());
	astroObject.hasEllipseInfo= true;

	//Set frequency to MGPS frequency
	astroObject.nu= 0.843;
	astroObject.dnu= 0;//unknown
	astroObject.hasFrequencyInfo= true;

	/*
	//Parse deconv ellipse info
	bmaj_str= fields[12];
	bmin_str= fields[13];
	pa_str= fields[14];
	CodeUtils::StripBlankSpaces(bmaj_str);
	CodeUtils::StripBlankSpaces(bmin_str);
	CodeUtils::StripBlankSpaces(pa_str);
	astroObject.bmaj_deconv= atof(bmaj_str.c_str());
	astroObject.bmin_deconv= atof(bmin_str.c_str());
	astroObject.pa_deconv= atof(pa_str.c_str());
	astroObject.hasDeconvEllipseInfo= true;
	*/

	return 0;

}//close ParseMGPSObjectData()


int AstroObjectParser::Read(std::vector<std::string>& data,std::string filename)
{
	//Open file for reading
	ifstream in;  
  in.open(filename.c_str());
  if(!in.good()) {
		#ifdef LOGGING_ENABLED
			ERROR_LOG("Failed to open file "<<filename<<" for reading!");
		#endif
		return -1;
  }

	//Parsing file
	#ifdef LOGGING_ENABLED
		DEBUG_LOG("Reading and parsing file: "<<filename);
	#endif

	std::string parsedline= "";	
	
	while(std::getline(in,parsedline)) 
	{
		//Check file
		if (!in.good()) break;

		//Check first character
		std::string field= "";
		istringstream parsedline_stream(parsedline);
		parsedline_stream>>field;
		char first_char= *(field.c_str());
		if(first_char=='#' || first_char=='*' || first_char=='/' || first_char=='\n' || first_char==' '){//skip line
			#ifdef LOGGING_ENABLED
				DEBUG_LOG("--> skip line (line="<<parsedline<<")");
			#endif
			continue;
		}
	
		//Append data
		data.push_back(parsedline);
				
		//Exit at the end
		if (!in.good()) break;
		
	}//close while

	//Close file
	in.close();

	return 0;

}//close Read()


}//close namespace
