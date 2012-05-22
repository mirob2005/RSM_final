/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: FaceList.h 2404 2010-09-10 23:16:42Z mshafae $
 *
 * Indexed triangle face list data structure. Uses a 2D array of
 * floats for storage.
 */


#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
using namespace std;

#ifndef _FACELIST_H_
#define _FACELIST_H_

#ifndef MIN
#define MIN( x, y ) ((x) <= (y) ? (x) :  (y))
#endif
#ifndef MAX
#define MAX( x, y ) ((x) >= (y) ? (x) :  (y))
#endif
#ifndef ABS
#define ABS( x )    ((x) >= 0.0 ? (x) : -(x))
#endif


class FaceList{

private:
  // Model's file name
  std::string _filename;
  unsigned int _vertexCount;
  unsigned int _faceCount;

public:
  /*
   * Adjust the types of the following pointers
   * according to the name of your Vector class.
   */
  // Array of vertices (decimal number type)
  float **vertices;
  // Array of faces (indices into the array of vertices)
  int **faces;
  // Array of face normals (decimal number type)
  float **faceNormals;
  // Array of vertex normals (decimal number type)
  float **vertexNormals;

  FaceList( const char *filename, unsigned int vertexCount, unsigned int  faceCount ) : _filename( filename ), _vertexCount( vertexCount ), _faceCount( faceCount )
  {
    unsigned int i;
    if( !(vertices = (float**)calloc( _vertexCount, sizeof(float*))) ){
      std::cerr << "Didn't allocated memory for vertex array." << std::endl;
    }
    for( i = 0; i < _vertexCount; i++ ){
      if( !(vertices[i] = (float*)calloc( 4, sizeof(float) )) ){
        std::cerr << "Couldn't allocate memory for the vertices at location " << i << std::endl;
      }
    }
    
    if( !(faces = (int**)calloc( _faceCount, sizeof(int*))) ){
      std::cerr << "Didn't allocated memory for face array." << std::endl;
    }
    for( i = 0; i < _faceCount; i++ ){
      if( !(faces[i] = (int*)calloc( 3, sizeof(int) )) ){
        std::cerr << "Couldn't allocate memory for the vertices at location " << i << std::endl;
      }
    }
    
    if( !(faceNormals = (float**)calloc( _faceCount, sizeof(float*))) ){
      std::cerr << "Didn't allocated memory for vertex array." << std::endl;
    }
    for( i = 0; i < _faceCount; i++ ){
      if( !(faceNormals[i] = (float*)calloc( 3, sizeof(float) )) ){
        std::cerr << "Couldn't allocate memory for the vertices at location " << i << std::endl;
      }
    }

    if( !(vertexNormals = (float**)calloc( _vertexCount, sizeof(float*))) ){
      std::cerr << "Didn't allocated memory for vertex array." << std::endl;
    }
    for( i = 0; i < _vertexCount; i++ ){
      if( !(vertexNormals[i] = (float*)calloc( 3, sizeof(float) )) ){
        std::cerr << "Couldn't allocate memory for the vertices at location " << i << std::endl;
      }
    }
    
  };
  
  ~FaceList( ){
   	for( unsigned int i = 0; i < _vertexCount; i++ ){
   		free( vertices[i] );
      free( vertexNormals[i] );
   	}
    for( unsigned int i = 0; i < _faceCount; i++ ){
      free( faces[i] );
      free( faceNormals[i] );
    }
   	free( vertices );
    free( faces );
    free( faceNormals );
    free( vertexNormals );
  };

  unsigned int vertexCount( ){
    return _vertexCount;
  }

  unsigned int faceCount( ){
    return _faceCount;
  }
  
  std::string filename( ){
    return _filename;
  }
  
  void computeFaceNormals( ){
	float a[3];
	float b[3];
	float c[3];

	float u[3];
	float v[3];

	float N[3];
	float N_mag;

	for( unsigned int i = 0; i < _faceCount; i++ ){

		a[0] = vertices[faces[i][1]][0];
		a[1] = vertices[faces[i][1]][1];
		a[2] = vertices[faces[i][1]][2];

		b[0] = vertices[faces[i][0]][0];
		b[1] = vertices[faces[i][0]][1];
		b[2] = vertices[faces[i][0]][2];

		c[0] = vertices[faces[i][2]][0];
		c[1] = vertices[faces[i][2]][1];
		c[2] = vertices[faces[i][2]][2];

		u[0] = b[0]-a[0];
		u[1] = b[1]-a[1];
		u[2] = b[2]-a[2];

		v[0] = c[0]-a[0];
		v[1] = c[1]-a[1];
		v[2] = c[2]-a[2];

		N[0] = u[1]*v[2] - u[2]*v[1];
		N[1] = u[2]*v[0] - u[0]*v[2];
		N[2] = u[0]*v[1] - u[1]*v[0];

		N_mag = sqrt((N[0]*N[0]) + (N[1]*N[1]) + (N[2]*N[2]));

		faceNormals[i][0] = N[0]/N_mag;
		faceNormals[i][1] = N[1]/N_mag;
		faceNormals[i][2] = N[2]/N_mag;
	 
//		cout << "Face Normal = " << faceNormals[i][0] << " , " << faceNormals[i][1] << " , " << faceNormals[i][2] <<endl;
	}
  }
  
  void computeVertexNormals( ){

	float vertexNormals_mag;

	for( unsigned int i = 0; i < _faceCount; i++ ){
		for( unsigned int j = 0; j < 3; j++ ){

			vertexNormals[faces[i][j]][0] += faceNormals[i][0];
			vertexNormals[faces[i][j]][1] += faceNormals[i][1];
			vertexNormals[faces[i][j]][2] += faceNormals[i][2];			
		}
	}

	//NORMALIZE EACH VERTEX NORMAL
	for( unsigned int j = 0; j < _vertexCount; j++ ){
		vertexNormals_mag = sqrt(vertexNormals[j][0]*vertexNormals[j][0]+vertexNormals[j][1]*vertexNormals[j][1]+
			vertexNormals[j][2]*vertexNormals[j][2]);
		
		vertexNormals[j][0] = vertexNormals[j][0]/vertexNormals_mag;
		vertexNormals[j][1] = vertexNormals[j][1]/vertexNormals_mag;
		vertexNormals[j][2] = vertexNormals[j][2]/vertexNormals_mag;

	}
  }
};

#endif