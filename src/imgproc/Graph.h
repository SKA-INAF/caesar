
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
* @file Graph.h
* @class Graph
* @brief Graph data structures
*
* Graph data structure class & utils
* @author S. Riggi
* @date 07/08/2017
*/


#ifndef _GRAPH_h
#define _GRAPH_h 1

#include <Logger.h>

#include <TObject.h>
#include <TMatrixD.h>

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
#include <time.h>
#include <ctime>



namespace Caesar {


//===============================================
//==         UNDIRECTED GRAPH CLASS
//===============================================

class Graph : public TObject {

	public:
		/** 
		\brief Standard constructor
 		*/
		Graph();

		/** 
		\brief Constructor with nvertex initializer
 		*/
		Graph(size_t nvertex);

		/** 
		\brief Construct from adjacency matrix
 		*/
		Graph(TMatrixD const& A);

		/** 
		\brief Standard destructor
 		*/
		virtual ~Graph();


	public:
		/** 
		\brief Get number of vertex
 		*/
		int GetNVertexes(){return static_cast<int>(m_adj.size());}

		/** 
		\brief Add vertex
 		*/
		void AddVertex(){
			m_adj.push_back( std::list<int>(0) );
		}

		/** 
		\brief Add edge
 		*/
		int AddEdge(int v, int w);

		/** 
		\brief Get adjacency matrix
 		*/
		TMatrixD* GetAdjacencyMatrix();

		/** 
		\brief Get connected components in the aundirected graph
 		*/
		int GetConnectedComponents(std::vector< std::vector<int> >& connected_items);

		
	private:
		/** 
		\brief Get connected components recursion util
 		*/
		void DFSUtil(std::vector<int>& connected_items,int index, std::vector<bool>& visited);

				
	protected:

		//- Adjancency lists
		std::vector< std::list<int> > m_adj;

	ClassDef(Graph,1)

};//close class Graph

#ifdef __MAKECINT__
#pragma link C++ class Graph+;
#endif	

}//close namespace

#endif
