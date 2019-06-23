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
* @file CutParser.h
* @class CutParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#ifndef _CUT_PARSER_h
#define _CUT_PARSER_h 1

#include <Cut.h>


#include <TObject.h>
#include <TBranch.h>
#include <TTree.h>

#include <vector>
#include <string>


namespace Caesar {


//=======================================
//==           CUT PARSER
//=======================================
class CutParser : public TObject 
{  
	public:

		/** 
		\brief Get cut parser instance
 		*/
		static CutParser& Instance() 
		{
    	static CutParser myInstance; 
			/*
			if(!m_hasRegisteredCuts){
				RegisterPredefinedCuts();
			}
			*/
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    CutParser(CutParser const&) = delete;             // Copy construct
    CutParser(CutParser&&) = delete;                  // Move construct
    CutParser& operator=(CutParser const&) = delete;  // Copy assign
    CutParser& operator=(CutParser &&) = delete;      // Move assign

	protected:
		/** 
		\brief Class constructor
 		*/
  	CutParser();
		/** 
		\brief Class destructor
 		*/
  	virtual ~CutParser();

	public:
		/** 
		\brief Read the cut file, parse and set cuts to be used by other classes
 		*/
		int Parse(std::string filename="");
		
		/** 
		\brief Print registered options
 		*/
		void PrintCuts(){
			return CutFactory::Instance().PrintCuts();
		}

		/** 
		\brief Retrieve cut
 		*/
		Cut* GetCut(std::string name)
		{
			return CutFactory::Instance().GetCut(name);
		}

		/** 
		\brief Retrieve cuts
 		*/
		std::map<std::string,Cut*> GetCuts()
		{
			return CutFactory::Instance().GetCuts();
		}
	
		/** 
		\brief Check registered cut
 		*/
		bool HasCut(std::string name)
		{
			return CutFactory::Instance().HasCut(name);
		}

		/** 
		\brief Check is cut passed
 		*/
		template<typename T>
		bool IsCutPassed(std::string name,T val)
		{
			return CutFactory::Instance().IsCutPassed<T>(name,val);
		}

		
	private:

		static std::string m_cutFile;
		//static bool m_hasRegisteredCuts;
		
	ClassDef(CutParser,1)
		
};


#ifdef __MAKECINT__
#pragma link C++ class CutParser+;
#endif

}//close namespace 

#endif
 
