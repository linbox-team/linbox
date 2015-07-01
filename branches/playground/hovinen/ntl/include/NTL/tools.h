
#ifndef NTL_tools__H
#define NTL_tools__H

#include <stdlib.h>
#include <iostream.h>
#include <math.h>

#include <NTL/config.h>
#include <NTL/mach_desc.h>
#include <NTL/IsFinite.h>



struct INIT_SIZE_STRUCT { };
const INIT_SIZE_STRUCT INIT_SIZE = INIT_SIZE_STRUCT();
typedef const INIT_SIZE_STRUCT& INIT_SIZE_TYPE;

struct INIT_VAL_STRUCT { };
const INIT_VAL_STRUCT INIT_VAL = INIT_VAL_STRUCT();
typedef const INIT_VAL_STRUCT& INIT_VAL_TYPE;

struct INIT_TRANS_STRUCT { };
const INIT_TRANS_STRUCT INIT_TRANS = INIT_TRANS_STRUCT();
typedef const INIT_TRANS_STRUCT& INIT_TRANS_TYPE;


struct INIT_LOOP_HOLE_STRUCT { };
const INIT_LOOP_HOLE_STRUCT INIT_LOOP_HOLE = INIT_LOOP_HOLE_STRUCT();
typedef const INIT_LOOP_HOLE_STRUCT& INIT_LOOP_HOLE_TYPE;

struct INIT_FFT_STRUCT { };
const INIT_FFT_STRUCT INIT_FFT = INIT_FFT_STRUCT();
typedef const INIT_FFT_STRUCT& INIT_FFT_TYPE;


#ifdef NTL_NO_INIT_TRANS
#define NTL_OPT_RETURN(t, x) return x
#else
#define NTL_OPT_RETURN(t, x) return t(x, INIT_TRANS)
#endif


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

inline void conv(int& x, int a) { x = a; }
inline void conv(int& x, long a) { x = int(a); }
inline void conv(int& x, float a) { x = int(floor(double(a))); }
inline void conv(int& x, double a) { x = int(floor(a)); }

inline int to_int(int a) { return a; }
inline int to_int(long a) { return int(a); }
inline int to_int(float a) { return int(floor(double(a))); }
inline int to_int(double a) { return int(floor(a)); }


inline void conv(long& x, int a) { x = a; }
inline void conv(long& x, long a) { x = a; }
inline void conv(long& x, float a) { x = long(floor(double(a))); }
inline void conv(long& x, double a) { x = long(floor(a)); }

inline long to_long(int a) { return a; }
inline long to_long(long a) { return a; }
inline long to_long(float a) { return long(floor(double(a))); }
inline long to_long(double a) { return long(floor(a)); }

inline void conv(float& x, int a) { x = float(a); }
inline void conv(float& x, long a) { x = float(a); }
inline void conv(float& x, float a) { x = a; }
inline void conv(float& x, double a) { x = float(a); }

inline float to_float(int a) { return float(a); }
inline float to_float(long a) { return float(a); }
inline float to_float(float a) { return a; }
inline float to_float(double a) { return float(a); }

inline void conv(double& x, int a) { x = double(a); }
inline void conv(double& x, long a) { x = double(a); }
inline void conv(double& x, float a) { x = double(a); }
inline void conv(double& x, double a) { x = a; }

inline double to_double(int a) { return double(a); }
inline double to_double(long a) { return double(a); }
inline double to_double(float a) { return double(a); }
inline double to_double(double a) { return a; }


long SkipWhiteSpace(istream& s);



void Error(const char *s);


#if (!defined(NTL_CPLUSPLUS_ONLY)) 
extern "C"
#endif
double GetTime();


void PrintTime(ostream& s, double t);


#endif

