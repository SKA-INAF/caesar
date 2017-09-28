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
* @file Graph.cc
* @class Graph
* @brief Graph data structures
*
* Graph data structure class & utils
* @author S. Riggi
* @date 07/08/2017
*/

#include <Graph.h>
#include <CodeUtils.h>
#include <Logger.h>

#include <TObject.h>
#include <TMatrixD.h>

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

ClassImp(Caesar::Graph)

namespace Caesar {

Graph::Graph()
{

}

Graph::Graph(size_t nvertex)
{
	m_adj.assign(nvertex,std::list<int>(0));
}

Graph::Graph(TMatrixD const& A)
{
	//Fill adajcency list from matrix			
	for(int i=0;i<A.GetNrows();i++){
		m_adj.push_back( std::list<int>() );
		for(int j=0;j<A.GetNcols();j++){
			if(A(i,j)>0) m_adj[i].push_back(j);
		}//end loop links per vertex
	}//end loop vertex
	
}//close constructor

Graph::~Graph()
{

	m_adj.clear();

}//close destructor


int Graph::AddEdge(int v, int w){
	
	//Check if vertexes has been allocated
	if(m_adj.empty() || v<0 || v>=m_adj.size() || w<0 || w>=m_adj.size() ){
		ERROR_LOG("Invalid edge specified ("<<v<<","<<w<<") (hint: check if graph has vertex)!");
		return -1;
	}
			
	//Check if already inserted
	std::list<int>::iterator it = std::find(m_adj[v].begin(), m_adj[v].end(), w);
	if(it!=m_adj[v].end()){
		WARN_LOG("Edge ("<<v<<","<<w<<") already inserted in graph!");
		return -1;
	}

   m_adj[v].push_back(w);
   m_adj[w].push_back(v);
			
	return 0;

}//close AddEdge()

TMatrixD* Graph::GetAdjacencyMatrix(){

	//Check if graph has vertexes
	int nVertexes= GetNVertexes();
	if(nVertexes<=0) {
		WARN_LOG("No vertexes present in this graph!");
		return nullptr;
	}

	//Create and fill adjacency matrix
	TMatrixD* A= new TMatrixD(nVertexes,nVertexes);
	A->Zero();
	for(size_t i=0;i<m_adj.size();i++){	
  	for(std::list<int>::iterator it = m_adj[i].begin(); it!=m_adj[i].end(); ++it){
			int index= *it;
			(*A)(i,index)= 1;
		}
	}	
	
	return A;
		
}//close GetAdjacencyMatrix()


void Graph::DFSUtil(std::vector<int>& connected_items, int index, std::vector<bool>& visited) 
{
	// Mark the current node as visited and print it
  visited[index] = true;
	connected_items.push_back(index);
  
  // Recur for all the vertices adjacent to this vertex
  for(std::list<int>::iterator it = m_adj[index].begin(); it != m_adj[index].end(); ++it){
		int list_index= *it;
  	if(!visited[list_index]) DFSUtil(connected_items,list_index, visited);
	}

}//close DFSUtil()


int Graph::GetConnectedComponents(std::vector< std::vector<int> >& connected_items)
{
	//Check if graph has data
	if(m_adj.empty()){
		WARN_LOG("No data present in graph!");
		return -1;
	}
 
	// Mark all the vertices as not visited
	std::vector<bool> visited(m_adj.size(),false);

 	//Loop over vertexes
	int counter= 0;
	for(size_t i=0;i<m_adj.size();i++) {
		//Skip if already visited
 		if (visited[i]) continue;
        
    // Print all reachable vertices from v
		connected_items.push_back( std::vector<int>() );
    DFSUtil(connected_items[counter],i, visited);
		counter++;

  }//end loop vertexes

	return 0;

}//close GetConnectedComponents()


}//close namespace



