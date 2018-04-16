/*======================================
  Head file for triangle mesh in 2D
=========================================*/

#ifndef PHG_PHG2D_H


#define Dim 2
#define GRID_TYPE_TRI 3
#define PHG_GRID_TYPE GRID_TYPE_TRI


/* Triangle */
#if PHG_GRID_TYPE == GRID_TYPE_TRI
#define NVert		(Dim + 1)	/**< Number of vertices / element */
#define NEdge		3
#define NVertFace	3	     /**< Number of vertex on face */
#define NEdgeFace	3	     /**< Number of edge on face */

#endif




#define PHG_PHG2D_H
#endif



