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
* @file ConfigParser.h
* @class ConfigParser
* @brief Parse the configuration file containing program parameters
* 
* @author S. Riggi
* @date 25/04/2010
*/

#ifndef ConfigParser_h
#define ConfigParser_h 1

#include <Option.h>

#include <json/json.h>

#include <TObject.h>
#include <TBranch.h>
#include <TTree.h>

#include <vector>
#include <string>


namespace Caesar{


class ConfigParser : public TObject {
  
	public:
		static ConfigParser& Instance() {
    	// Since it's a static variable, if the class has already been created,
      // It won't be created again.
      // And it is thread-safe in C++11.
      static ConfigParser myInstance;
 
			if(!m_HasRegisteredOptions){
				RegisterPredefinedOptions();
			}

      // Return a reference to our instance.
      return myInstance;
    }
 
    // delete copy and move constructors and assign operators
    ConfigParser(ConfigParser const&) = delete;             // Copy construct
    ConfigParser(ConfigParser&&) = delete;                  // Move construct
    ConfigParser& operator=(ConfigParser const&) = delete;  // Copy assign
    ConfigParser& operator=(ConfigParser &&) = delete;      // Move assign

	protected:
		/** 
		\brief Class constructor
 		*/
  	ConfigParser();
		/** 
		\brief Class destructor
 		*/
  	virtual ~ConfigParser();

	public:
		/** 
		\brief Read the config file, parse and set info to be used by other classes
 		*/
		int Parse(std::string filename="");
		
		/** 
		\brief Print registered options
 		*/
		void PrintOptions(){
			return OptionFactory::Instance().PrintOptions();
		}

		/** 
		\brief Register options
 		*/
		template <typename T>
			std::shared_ptr<Option<T>> GetOption(std::string name){
				return OptionFactory::Instance().GetOption<T>(name);
			}//close GetOption()
	
		/** 
		\brief Check registered options
 		*/
		bool HasOption(std::string name){//Find option name
			return OptionFactory::Instance().HasOption(name);
		}
		/** 
		\brief Set option from string
 		*/
		int SetOptionFromString(std::string name,std::string stringified_value){
			return OptionFactory::Instance().SetOptionFromString(name,stringified_value);
		}

		/** 
		\brief Get option value
 		*/
		template <typename T>
			int GetOptionValue(std::string name,T& value){
				return OptionFactory::Instance().GetOptionValue(name,value);
			}//close GetOptionValue()

		/** 
		\brief Get config option tree
 		*/
		TTree* GetConfigTree(std::string treeName="CfgInfo"){
			return OptionFactory::Instance().MakeTTree(treeName);
		}
	
	private:
		static int RegisterPredefinedOptions();
		
	private:

		static std::string m_ConfigFile;
		static bool m_HasRegisteredOptions;
		
	ClassDef(ConfigParser,1)
		
};


#ifdef __MAKECINT__
#pragma link C++ class ConfigParser+;
#endif

}//close namespace 

#endif
 
