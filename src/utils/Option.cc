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
* @file Option.cc
* @class Option
* @brief Define configuration option
* 
* @author S. Riggi
* @date 5/6/2016
*/

#include <Option.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

ClassImp(Caesar::OptionBase)
ClassImp(Caesar::Option<int>)
ClassImp(Caesar::Option<long int>)
ClassImp(Caesar::Option<double>)
ClassImp(Caesar::Option<float>)
ClassImp(Caesar::Option<std::string>)
ClassImp(Caesar::Option<bool>)
ClassImp(Caesar::Option<char*>)
ClassImp(Caesar::Option<std::vector<int>>)
ClassImp(Caesar::Option<std::vector<long int>>)
ClassImp(Caesar::Option<std::vector<float>>)
ClassImp(Caesar::Option<std::vector<double>>)
ClassImp(Caesar::Option<std::vector<std::string>>)
ClassImp(Caesar::Option<std::vector<bool>>)
ClassImp(Caesar::Option<std::vector<char*>>)

ClassImp(Caesar::OptionHelper<int>)
ClassImp(Caesar::OptionHelper<long int>)
ClassImp(Caesar::OptionHelper<double>)
ClassImp(Caesar::OptionHelper<float>)
ClassImp(Caesar::OptionHelper<std::string>)
ClassImp(Caesar::OptionHelper<char*>)
ClassImp(Caesar::OptionHelper<std::vector<int>>)
ClassImp(Caesar::OptionHelper<std::vector<long int>>)
ClassImp(Caesar::OptionHelper<std::vector<float>>)
ClassImp(Caesar::OptionHelper<std::vector<double>>)
ClassImp(Caesar::OptionHelper<std::vector<std::string>>)
ClassImp(Caesar::OptionHelper<std::vector<bool>>)
ClassImp(Caesar::OptionHelper<std::vector<char*>>)

ClassImp(Caesar::OptionFactory)


namespace Caesar {

/*
//== OPTION SPECIALIZATION ==
template <> 
int Option<long int>::GetJson(Json::Value& jsonObj) {
 	jsonObj["name"]= m_name;
	jsonObj["value"]= Json::Value( static_cast<Json::Int64>(m_value) );
	jsonObj["default"]= Json::Value( static_cast<Json::Int64>(m_defaultValue) );
	jsonObj["min"]= Json::Value( static_cast<Json::Int64>(m_minValue) );
	jsonObj["max"]= Json::Value( static_cast<Json::Int64>(m_maxValue) );
	return 0;
}

template <typename T>
int Option<T>::GetJson(Json::Value& jsonObj){	
	jsonObj["name"]= m_name;
	jsonObj["value"]= m_value;	
	jsonObj["default"]= m_defaultValue;
	jsonObj["min"]= m_minValue;
	jsonObj["max"]= m_maxValue;
	return 0;
}
template int Option<int>::GetJson(Json::Value& jsonObj);
template int Option<std::string>::GetJson(Json::Value& jsonObj);
template int Option<char*>::GetJson(Json::Value& jsonObj);
template int Option<double>::GetJson(Json::Value& jsonObj);
template int Option<float>::GetJson(Json::Value& jsonObj);
template int Option<bool>::GetJson(Json::Value& jsonObj);
*/


}//close namespace

