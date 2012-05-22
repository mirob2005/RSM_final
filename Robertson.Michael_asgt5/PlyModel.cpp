/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: PlyModel.cpp 2399 2010-09-09 20:47:35Z mshafae $
 *
 * Implementation of a PLY format reader. Assumes the input file
 * is ASCII.
 */

#include "PlyModel.h"
#include <cstring>
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

FaceList* readPlyModel( const char* filename ){
  char buffer[255], type[128], c;
  ifstream inputfile;
  unsigned int i;
  float v[3];
  int f[3], k;
  unsigned int nv;
  unsigned int nf;
  FaceList *fl;
  assert( filename );
  inputfile.open( filename, ios::in );
  if( inputfile.fail( ) ){
    //cerr << "File \"" << filename << "\" not found." << endl;
    //exit( 1 );
	return(0);
  }
  // Parse the header
  if(inputfile.getline(buffer, sizeof(buffer), '\n') != NULL){
    if( strcmp(buffer, "ply") != 0){
      cerr << "Error: Input file is not of .ply type." << endl;
      exit(1);
    }
  }else{
    cerr << "End of input?" << endl;
    exit( 1 );
  }
  if(inputfile.getline(buffer, sizeof(buffer), '\n') != NULL){
    if( strncmp(buffer, "format ascii", 12) != 0){
      cerr << "Error: Input file is not in ASCII format." << endl;
      exit(1);
    }
  }else{
    cerr << "End of input?" << endl;
    exit( 1 );
  }
  if(inputfile.getline(buffer, sizeof(buffer), '\n') != NULL){
    while (strncmp(buffer, "comment", 7) == 0){
      inputfile.getline(buffer, sizeof(buffer), '\n');
  }
  }else{
    cerr << "End of input?" << endl;
    exit( 1 );
  }

  if (strncmp(buffer, "element vertex", 14) == 0){
    sscanf(buffer, "element vertex %u\n", &nv);
  }else{
    cerr << "Error: number of vertices expected." << endl;
    exit(1);
  }

  i = 0;
  inputfile.getline(buffer, sizeof(buffer), '\n');
  while (strncmp(buffer, "property", 8) == 0) {
    if (i < 3) {
      sscanf(buffer, "property %s %c\n", type, &c);
      switch (i) {
      case 0:
        if (c != 'x') {
          cerr << "Error: first coordinate is not x." << endl;
          exit(1);
        }
        break;
      case 1:
        if (c != 'y') {
          cerr << "Error: first coordinate is not y." << endl;
          exit(1);
        }
        break;
      case 2:
        if (c != 'z') {
          cerr << "Error: first coordinate is not z." << endl;
          exit(1);
        }
        break;
      default:
        break;
      }
      i++;
    }
    inputfile.getline(buffer, sizeof(buffer), '\n');
  }
  
  if (strncmp(buffer, "element face", 12) == 0)
    sscanf(buffer, "element face %u\n", &nf);
  else {
    cerr << "Error: number of faces expected." << endl;
    exit(1);
  }

  inputfile.getline(buffer, sizeof(buffer), '\n');
  if (strncmp(buffer, "property list", 13) != 0) {
    cerr << "Error: property list expected." << endl;
    exit(1);
  }
  
  inputfile.getline(buffer, sizeof(buffer), '\n');
  while (strncmp(buffer, "end_header", 10) != 0){
    inputfile.getline(buffer, sizeof(buffer), '\n');
  }
  
  // Allocate FaceList object
  if( !(fl = new FaceList( filename, nv, nf )) ){
    cerr << "Could not allocate a new face list for the model." << endl;
    exit(1);
  }

  // Process the body of the input file
  // read vertex data from PLY file
  for (i = 0; i < nv; i++) {
    inputfile.getline(buffer, sizeof(buffer), '\n');
    sscanf(buffer,"%f %f %f", &v[0], &v[1], &v[2]);
    for(int j = 0; j < 3; j++ ){
      fl->vertices[i][j] = v[j];
    }
    fl->vertices[i][3] = 1.0;
  }

  // read the face data from the PLY file
  for (i = 0; i < nf; i++) {
    inputfile.getline(buffer, sizeof(buffer), '\n');
    sscanf(buffer, "%d %d %d %d", &k, &f[0], &f[1], &f[2] );
    if (k != 3) {
      fprintf(stderr, "Error: not a triangular face.\n");
      exit(1);
    }
    for(int j = 0; j < 3; j++ ){
      fl->faces[i][j] = f[j];
    }
  }

  fl->computeFaceNormals( );
  fl->computeVertexNormals( );
  
  inputfile.close( );
  return( fl );
}

