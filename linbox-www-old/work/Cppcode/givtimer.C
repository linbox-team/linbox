// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id$
// ==========================================================================
// Description:
// - various timer objects
// - to be rewritten to be more efficient

#include <math.h>
// #include "givconfig.h"
#define GIVARO_SYS UNIX

#if ( GIVARO_SYS == _SYS_MACOS )
#include <time.h> // for use of clock() 
#define USER_TIME ((double) ( (double)clock()/ (double)CLOCKS_PER_SEC))
#else
extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
  int getrusage (int, struct rusage*) ;
}
#endif


   // Return a value to initialize random generator 
long BaseTimer::seed() 
{
#if (GIVARO_SYS == _SYS_MACOS)
  return clock() ;
#else
  struct timeval tp;
  gettimeofday(&tp, 0) ;
  return(tp.tv_usec);
#endif   
}

   // Output the value of the timer :
ostream& BaseTimer::print( ostream& o ) const 
{ return o << _t ; }

   // Some arithmetic operator :
BaseTimer& BaseTimer::operator = (const BaseTimer & T) 
{  
   _t = T._t ; 
   return *this ; 
}
      
   // Computes and returns interval of time
   // beteween *this and T
const BaseTimer BaseTimer::operator - (const BaseTimer & T) const
{
   BaseTimer Tmp ;
   Tmp._t = _t - T._t ; 
   return Tmp ;
}

const BaseTimer BaseTimer::operator - () 
{
   BaseTimer Tmp ;
   Tmp._t = -_t ; 
   return Tmp ;
}

const BaseTimer BaseTimer::operator + (const BaseTimer & T)  const
{
   BaseTimer Tmp ;
   Tmp._t = _t + T._t ; 
   return Tmp ;
}

   // Start timer
void RealTimer::start()
{  
#if (GIVARO_SYS == _SYS_MACOS)
  _t = USER_TIME ;
#else
   struct timeval tmp2 ; 
   gettimeofday (&tmp2, 0) ;

   // real time 
   _t = (double) tmp2.tv_sec + 
         ((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC ; 
#endif   
}


   // Stop timer 
void RealTimer::stop()
{ 
#if (GIVARO_SYS == _SYS_MACOS)
  _t = USER_TIME - _t ;
#else
   struct timeval tmp2 ;  
   gettimeofday (&tmp2, 0) ;

   // real time 
   _t = (double) tmp2.tv_sec + 
         ((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC - _t ; 
#endif   
}

   // Start timer
void UserTimer::start()
{
#if (GIVARO_SYS == _SYS_MACOS)
  _t = USER_TIME ;
#else
   struct rusage  tmp1 ;  // to getrusage (sys+user times)
   getrusage (RUSAGE_SELF, &tmp1) ;
   // user time
   _t = (double) tmp1.ru_utime.tv_sec +
         ((double) tmp1.ru_utime.tv_usec)/ (double)MSPSEC ;
#endif   
}


   // Stop timer
void UserTimer::stop()
{
#if (GIVARO_SYS == _SYS_MACOS)
  _t = USER_TIME - _t;
#else 
   struct rusage  tmp1 ;  // to getrusage (sys+user times)
   getrusage (RUSAGE_SELF, &tmp1) ;
   // user time
   _t = (double) tmp1.ru_utime.tv_sec +
         ((double) tmp1.ru_utime.tv_usec)/ (double)MSPSEC - _t ;
#endif   
}


   // Start timer
void SysTimer::start()
{
#if (GIVARO_SYS == _SYS_MACOS)
  _t = USER_TIME ;
#else
  struct rusage  tmp1 ;  // to getrusage (sys+user times)
  getrusage (RUSAGE_SELF, &tmp1) ;
  // user time
  _t = (double) tmp1.ru_stime.tv_sec + 
       ((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC ;
#endif   
}


   // Stop timer
void SysTimer::stop()
{
#if (GIVARO_SYS == _SYS_MACOS)
  _t = USER_TIME - _t ;
#else
   struct rusage  tmp1 ;  // to getrusage (sys+user times)
   getrusage (RUSAGE_SELF, &tmp1) ;
   // user time
   _t = (double) tmp1.ru_stime.tv_sec +
         ((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC - _t ;
#endif   
}



   // Clear timer :
void Timer::clear() 
{ rt.clear() ; ut.clear(); st.clear() ; }

   // Start timer
void Timer::start() 
{ rt.start() ; ut.start(); st.start() ; }

  // Stop timer
void Timer::stop() 
{ rt.stop() ; ut.stop(); st.stop() ; }


ostream& Timer::print( ostream& o ) const
{
   o << "user time: " << usertime() << '\n' ;
   o << "sys. time: " << systime() << '\n' ;
   return o << "real time: " << realtime() << endl ;
}

   // Some arithmetic operator :
Timer& Timer::operator = (const Timer & T)
{
  ut = T.ut ; 
  st = T.st ; 
  rt = T.rt ;
  return *this ;
}

   // Comput._tes and returns interval of time
   // beteween *this and T
const Timer Timer::operator - (const Timer & T)  const
{
  Timer Tmp ;
  Tmp.ut = ut - T.ut ;
  Tmp.st = st - T.st ;
  Tmp.rt = rt - T.rt ;
  return Tmp ;
}

const Timer Timer::operator - ()
{
   Timer Tmp ;
   Tmp.ut = -ut ;
   Tmp.st = -st ;
   Tmp.rt = -rt ;
   return Tmp ;
}

const Timer Timer::operator + (const Timer & T)  const
{
   Timer Tmp ;
   Tmp.ut = ut + T.ut ;
   Tmp.st = st + T.st ;
   Tmp.rt = rt + T.rt ;
   return Tmp ;
}
