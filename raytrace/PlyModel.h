/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: PlyModel.h 2408 2010-09-14 19:47:33Z mshafae $
 *
 * Reads a PLY format model file and returns an indexed triangle face list.
 *
 */

#ifndef _PLYMODEL_H_
#define _PLYMODEL_H_

#include "FaceList.h"

FaceList* readPlyModel( const char* filename );

#endif