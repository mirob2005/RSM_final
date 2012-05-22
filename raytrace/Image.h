/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: Image.h 1961 2010-02-24 08:46:53Z mshafae $
 *
 *
 * The image class defines a trivial encoding for images known as PPM
 * format; it simply consists of an array or RGB triples, with one byte
 * per component, preceeded by a simple header giving the size of the
 * image.
 */

#include <cstdio>
#ifdef __APPLE__
/* Apple's weird location of their OpenGL & GLUT implementation */
#include <OpenGL/OpenGL.h>
#else
#include <GL/glut.h>
#endif

#ifndef _IMAGE_H_
#define _IMAGE_H_


typedef unsigned char channel;  

class Pixel {
public:
  Pixel( ) { r = 0; g = 0; b = 0; };
  Pixel( channel _r, channel _g, channel _b ) { r = _r; g = _g; b = _b; };
  void write( std::ostream &out ) const{
    out << "[" << r << ", " << g << ", " << b <<", " << "]";
  };
  channel r;
  channel g;
  channel b;
};

std::ostream& operator <<( std::ostream &out, const Pixel &p );

class Image {
public:
  inline Image( int x_res, int y_res, float *background );
  inline ~Image( ) { delete[] pixels; }
  inline GLubyte* read( const char *file_name );
  inline bool write( const char *file_name );
  inline Pixel &operator()( int i, int j ) { return *( pixels + ( i * width + j ) ); }  
  Pixel *pixels;
  int    width;
  int    height;
};

inline Image::Image( int x_res, int y_res, float *background ){
  width  = x_res;
  height = y_res;
  pixels = new Pixel[ width * height ];
  Pixel *p = pixels;
  //printf("%f,%f,%f",background[0]*255,background[1]*255,background[2]*255);
  for( int i = 0; i < width * height; i++ ) *p++ = Pixel(background[0]*255,background[1]*255,background[2]*255);
}

inline GLubyte* Image::read( const char *file_name ){
  char buffer[100];
  FILE *fp;
  int size_x, size_y, maxval;
  unsigned char c;
  int i;

  // open file
  if ((fp=fopen (file_name, "rb"))==NULL){
      fprintf (stderr, "unable to open file%c\n", 7);
      exit (1);
  }

  // read file identifier (magic number)
  fgets (buffer, sizeof (buffer), fp);
  if ((buffer[0] != 'P') || (buffer[1] != '6')){
      fprintf (stderr, "incorrect file type%c\n", 7);
      exit (1);
  }

  // read image size
  do
      fgets (buffer, sizeof (buffer), fp);
  while (buffer[0] == '#' || buffer[0] == ' ');

  sscanf (buffer, "%d %d", &size_x, &size_y);

  printf( "Image width: %d, Image height: %d\n", size_x, size_y );

  // read maximum pixel value (usually 255)
  do{
      fgets (buffer, sizeof (buffer), fp);
  }while (buffer[0] == '#');

  sscanf (buffer, "%d", &maxval);

  // allocate RGBA texture buffer
  GLubyte *texture = (GLubyte *)malloc(size_x*size_y*4*sizeof(GLubyte));

  // read RGB data and calculate alpha value
  for (i=0; i < size_x*size_y*4; i++){
      // insert alpha value (0 or 255) after each RGB
      if ((i%4) == 3){
          // alpha channel example: make all nearly black pixels transparent
          texture[i]=(GLubyte)
(((texture[i-3]+texture[i-2]+texture[i-1])<10)?0:255);
      }else{
          c=fgetc(fp);
          texture[i]=(GLubyte) c;
      }
  }

  // close input file
  fclose(fp);
  return(texture);
}

inline bool Image::write( const char *file_name ){
  Pixel *p = pixels;
  FILE  *fp = fopen( file_name, "w+b" );
  if( fp == NULL ){
    return false;
  }
  fprintf( fp, "P6\n%d %d\n255\n", width, height );
  for( int i = 0; i < width * height; i++ ){
    fprintf( fp, "%c%c%c", p->r, p->g, p->b );
    p++;
  }
  fclose( fp );
  return true;
}

#endif
