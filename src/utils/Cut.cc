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
* @file Cut.cc
* @class Cut
* @brief Define cuts
* 
* @author S. Riggi
* @date 20/6/2019
*/

#include <Cut.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>

using namespace std;

ClassImp(Caesar::Cut)

ClassImp(Caesar::EqualityCut<int>)
ClassImp(Caesar::EqualityCut<long int>)
ClassImp(Caesar::EqualityCut<double>)
ClassImp(Caesar::EqualityCut<float>)
ClassImp(Caesar::EqualityCut<bool>)

ClassImp(Caesar::BoundCut<int>)
ClassImp(Caesar::BoundCut<long int>)
ClassImp(Caesar::BoundCut<double>)
ClassImp(Caesar::BoundCut<float>)
ClassImp(Caesar::BoundCut<bool>)

ClassImp(Caesar::SingleBoundCut<int>)
ClassImp(Caesar::SingleBoundCut<long int>)
ClassImp(Caesar::SingleBoundCut<double>)
ClassImp(Caesar::SingleBoundCut<float>)
ClassImp(Caesar::SingleBoundCut<bool>)

ClassImp(Caesar::CutFactory)


namespace Caesar {


}//close namespace
