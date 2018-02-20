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
* @file TaskData.h
* @class TaskData
* @brief TaskData class
*
* TaskData class
* @author S. Riggi
* @date 20/01/2015
*/
#ifndef _TASK_DATA_H
#define _TASK_DATA_H

#include <Source.h>

#include <string>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <iomanip>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <map>
#include <time.h>
#include <sys/time.h>

namespace Caesar {

//=======================================
//==    TILE CLASS                    ===
//=======================================
class Tile {
	
	public:
		Tile(long int xmin,long int xmax,long int ymin,long int ymax)
			: ix_min(xmin), ix_max(xmax), iy_min(ymin), iy_max(ymax)
		{};
		virtual ~Tile(){};

	public:	
	
		/*
		bool IsTileAdjacentX(Tile& const aTile){
			long int ix_min_N= aTile.ix_min;
			long int ix_max_N= aTile.ix_max;
			bool isAdjacentInX= (ix_max==ix_min_N-1 || ix_min==ix_max_N+1 || (ix_min==ix_min_N && ix_max==ix_max_N));
			return isAdjacentInX;
		}
		bool IsTileAdjacentY(Tile& const aTile){
			long int iy_min_N= aTile.iy_min;
			long int iy_max_N= aTile.iy_max;
			bool isAdjacentInY= (iy_max==iy_min_N-1 || iy_min==iy_max_N+1 || (iy_min==iy_min_N && iy_max==iy_max_N));
			return isAdjacentInY;
		}
		*/
		bool IsTileAdjacent(Tile const& aTile){
			long int ix_min_N= aTile.ix_min;
			long int ix_max_N= aTile.ix_max;
			long int iy_min_N= aTile.iy_min;
			long int iy_max_N= aTile.iy_max;
			bool isAdjacentInX= (ix_max==ix_min_N-1 || ix_min==ix_max_N+1 || (ix_min==ix_min_N && ix_max==ix_max_N));
			bool isAdjacentInY= (iy_max==iy_min_N-1 || iy_min==iy_max_N+1 || (iy_min==iy_min_N && iy_max==iy_max_N));
			bool isAdjacent= isAdjacentInX && isAdjacentInY;
			return isAdjacent;
		}

		bool IsTileOverlapping(Tile const& aTile){
			long int ix_min_N= aTile.ix_min;
			long int ix_max_N= aTile.ix_max;
			long int iy_min_N= aTile.iy_min;
			long int iy_max_N= aTile.iy_max;
			if (ix_max < ix_min_N) return false; // a is left of b
    	if (ix_min > ix_max_N) return false; // a is right of b
    	if (iy_max < iy_min_N) return false; // a is above b
    	if (iy_min > iy_max_N) return false; // a is below b
			return true;
		}

		/**
		* \brief Check if point is inside tile
		*/
		bool IsInsideTile(double x,double y){
			bool isInside= ( (x>=ix_min && x<=ix_max) && (y>=iy_min && y<=iy_max) );
			return isInside;
		}
		
	public:
		long int ix_min;
		long int ix_max;
		long int iy_min;
		long int iy_max;
		
};//close class Tile

//=======================================
//==    TASK DATA CLASS               ===
//=======================================
class TaskData : public TObject {
	
	public:
		/**
		* \brief Standard constructor
		*/		
		TaskData();
				
		/**
		* \brief Copy constructor
		*/
		TaskData(const TaskData& data);

		/**
		* \brief Destructor
		*/
		virtual ~TaskData();

	public:
		/**
		* \brief Operator =
		*/
		TaskData& operator=(const TaskData& data);

		/**
		* \brief Copy method
		*/
		void Copy(TObject &obj) const;

		/**
		* \brief Clear sources
		*/
		void ClearSources();

	public:
			
		/**
		* \brief Set tile
		*/
		int SetTile(long int xmin,long int xmax,long int ymin,long int ymax){
			if(xmin<0 || xmax<0 || xmin>=xmax) return -1;
			if(ymin<0 || ymax<0 || ymin>=ymax) return -1;
			ix_min= xmin;
			ix_max= xmax;
			iy_min= ymin;
			iy_max= ymax;
			return 0;
		}

		/**
		* \brief Add neighbor info
		*/
		void AddNeighborInfo(long int tid,long int wid){
			neighborTaskId.push_back(tid);
			neighborWorkerId.push_back(wid);
		}

		/**
		* \brief Check if point is inside task tile
		*/
		bool IsInsideTaskTile(double x,double y){
			bool isInside= ( (x>=ix_min && x<=ix_max) && (y>=iy_min && y<=iy_max) );
			return isInside;
		}

		/**
		* \brief Check if this task tile is adjacent to another given task
		*/
		bool IsTaskTileAdjacent(TaskData* aTask){
			if(!aTask) return false;
			long int ix_min_N= aTask->ix_min;
			long int ix_max_N= aTask->ix_max;
			long int iy_min_N= aTask->iy_min;
			long int iy_max_N= aTask->iy_max;
			bool isAdjacentInX= (ix_max==ix_min_N-1 || ix_min==ix_max_N+1 || (ix_min==ix_min_N && ix_max==ix_max_N));
			bool isAdjacentInY= (iy_max==iy_min_N-1 || iy_min==iy_max_N+1 || (iy_min==iy_min_N && iy_max==iy_max_N));
			bool isAdjacent= isAdjacentInX && isAdjacentInY;
			return isAdjacent;
		}

		/**
		* \brief Check if this task tile is overlapping with another given task
		*/
		bool IsTaskTileOverlapping(TaskData* aTask){
			if(!aTask) return false;
			long int ix_min_N= aTask->ix_min;
			long int ix_max_N= aTask->ix_max;
			long int iy_min_N= aTask->iy_min;
			long int iy_max_N= aTask->iy_max;
			if (ix_max < ix_min_N) return false; // a is left of b
    	if (ix_min > ix_max_N) return false; // a is right of b
    	if (iy_max < iy_min_N) return false; // a is above b
    	if (iy_min > iy_max_N) return false; // a is below b
			return true;
		}
	
		/**
		* \brief Check if this task tile is neighbor (adjacent or overlapping) to another given task
		*/
		bool IsTaskTileNeighbor(TaskData* aTask){
			if(!aTask) return false;
			bool isOverlapping= IsTaskTileOverlapping(aTask);
			bool isAdjacent= IsTaskTileAdjacent(aTask);
			bool isNeighbor= (isAdjacent || isOverlapping);
			return isNeighbor;
		}

	public: 

		//Task info		
		long int workerId;
		std::vector<long int> neighborTaskId;
		std::vector<long int> neighborWorkerId;
	
		//- Image coordinates
		long int ix_min;
		long int ix_max;
		long int iy_min;
		long int iy_max;

		//- Detected sources
		Source* source;
		std::vector<Source*> sources;
		std::vector<Source*> ext_sources;
		std::vector<Source*> sources_edge;
		std::vector<Source*> ext_sources_edge;
		
	ClassDef(TaskData,1)
};


#ifdef __MAKECINT__
#pragma link C++ class TaskData+;
#endif


}//close namespace

#endif



