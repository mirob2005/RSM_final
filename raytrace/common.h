/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: common.h 1961 2010-02-24 08:46:53Z mshafae $
 *
 */
 
#include <iostream>

#ifndef _COMMON_H_
#define _COMMON_H_

#define SQR(x) ((x)*(x))

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

typedef unsigned int uint;
typedef unsigned char uchar;


#define DELETEIFNOTNULL(pointer){ \
  if (pointer){ delete pointer; pointer = 0;} \
}

#define fopen_error(fname) { \
 std::cerr<<"Failed to open file "<<fname<<endl; \
 exit(0); \
}

#define fatal_error(err_string) { \
 std::cerr<<err_string<<endl; \
 exit(0); \
}


/* For use in hash_map */
struct eqint
{
  bool operator()(int n1, int n2) const
  {
    return n1 == n2;
  }
};


// Returns a random number between min and max, with the indicated precission between limits
double rangerand(double min, double max, int steps=10000);


// Description of the log levels:
// 0 - No debugging at all. Not even critical messages. No messages should be leveled at 0.
// 1 - Critical messages
// 2 - INFORMATIONAL : High level Information, eq coord of a created plane
// 3 - Low level info, eq about an image/vector  created

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 2
#endif

#define LOG(n, format, arg...) if (DEBUG_LEVEL >= (n)) fprintf(stderr, format , ## arg)
#define LOG2(n, arg) if (DEBUG_LEVEL >= (n)) cerr << arg << endl

#endif //_COMMON_H_

