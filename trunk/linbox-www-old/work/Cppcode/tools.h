

/***********************************************************************

   This software is for research and educational purposes only.

************************************************************************/



#ifndef TOOLS__H
#define TOOLS__H

#include <stdlib.h>
#include <iostream.h>


// This looks very wierd, but it was the only thing I could
// find that would work well across several compilers.

struct INIT_SIZE_TYPE { INIT_SIZE_TYPE(int,int) { } };
#define INIT_SIZE INIT_SIZE_TYPE(0,0)

struct INIT_VAL_TYPE { INIT_VAL_TYPE(int,int) { } };
#define INIT_VAL INIT_VAL_TYPE(0,0)


struct INIT_ADD_TYPE { INIT_ADD_TYPE(int,int) { } };
#define INIT_ADD INIT_ADD_TYPE(0,0)

struct INIT_SUB_TYPE { INIT_SUB_TYPE(int,int) { } };
#define INIT_SUB INIT_SUB_TYPE(0,0)

struct INIT_MUL_TYPE { INIT_MUL_TYPE(int,int) { } };
#define INIT_MUL INIT_MUL_TYPE(0,0)

struct INIT_DIV_TYPE { INIT_DIV_TYPE(int,int) { } };
#define INIT_DIV INIT_DIV_TYPE(0,0)

struct INIT_REM_TYPE { INIT_REM_TYPE(int,int) { } };
#define INIT_REM INIT_REM_TYPE(0,0)

struct INIT_NEG_TYPE { INIT_NEG_TYPE(int,int) { } };
#define INIT_NEG INIT_NEG_TYPE(0,0)

struct INIT_TRANS_TYPE { INIT_TRANS_TYPE(int,int) { } };
#define INIT_TRANS INIT_TRANS_TYPE(0,0)

struct INIT_LOOP_HOLE_TYPE { INIT_LOOP_HOLE_TYPE(int,int) { } };
#define INIT_LOOP_HOLE INIT_LOOP_HOLE_TYPE(0,0)



#define NAME2(x,y) _NAME2_AUX(x,y)
#define _NAME2_AUX(x,y) x ## y

#define __AKA__(x) NAME2(_XXX_,x)


inline int min(int a, int b) { return (a < b) ?  a : b; } 
inline int max(int a, int b) { return (a < b) ? b : a; }

inline long min(long a, long b) { return (a < b) ?  a : b; } 
inline long max(long a, long b) { return (a < b) ? b : a; }

inline long min(int a, long b) { return (a < b) ?  long(a) : b; } 
inline long max(int a, long b) { return (a < b) ? b : long(a); }

inline long min(long a, int b) { return (a < b) ?  a : long(b); } 
inline long max(long a, int b) { return (a < b) ? long(b) : a; }


inline void swap(long& a, long& b)  {  long t;  t = a; a = b; b = t; }
inline void swap(int& a, int& b)  {  int t;  t = a; a = b; b = t; }


void Error(const char *s);

#if (!defined(CPLUSPLUS_ONLY)) 
extern "C"
#endif
double GetTime();


#endif

