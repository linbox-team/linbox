/* linbox/util/timer.C
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
 *
 * See COPYING for license information.
 */
#ifndef __LINBOX_timer_INL
#define __LINBOX_timer_INL
// Description:
// - various timer objects
// - to be rewritten to be more efficient

#include <cmath>

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
//  int getrusage (int, struct rusage*) ;
}

#include <iostream>

#include "linbox/util/timer.h"

namespace LinBox 
{

// Return a value to initialize random generator 
long BaseTimer::seed() 
{
	struct timeval tp;
	gettimeofday(&tp, 0) ;
	return(tp.tv_usec);
}

// Output the value of the timer :
std::ostream& BaseTimer::print( std::ostream& o ) const 
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
	struct timeval tmp2 ; 
	gettimeofday (&tmp2, 0) ;

	// real time 
	_start_t = (double) tmp2.tv_sec + 
		((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC ; 
}


// Stop timer 
void RealTimer::stop()
{ 
	struct timeval tmp2 ;  
	gettimeofday (&tmp2, 0) ;

	// real time 
	_t = (double) tmp2.tv_sec + 
		((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC - _start_t ; 
}

// Start timer
void UserTimer::start()
{
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_start_t = (double) tmp1.ru_utime.tv_sec +
		((double) tmp1.ru_utime.tv_usec)/ (double)MSPSEC ;
}


// Stop timer
void UserTimer::stop()
{
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_t = (double) tmp1.ru_utime.tv_sec +
		((double) tmp1.ru_utime.tv_usec)/ (double)MSPSEC - _start_t ;
}


// Start timer
void SysTimer::start()
{
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_start_t = (double) tmp1.ru_stime.tv_sec + 
		((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC ;
}


// Stop timer
void SysTimer::stop()
{
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_t = (double) tmp1.ru_stime.tv_sec +
		((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC - _start_t ;
}



// Clear timer :
void Timer::clear() 
{ rt.clear() ; ut.clear(); st.clear(); _count = 0; }

// Start timer
void Timer::start() 
{ rt.start() ; ut.start(); st.start(); _count = 0; }

// Stop timer
void Timer::stop() 
{ rt.stop() ; ut.stop(); st.stop(); _count = 1; }


std::ostream& Timer::print( std::ostream& o ) const
{
	o << "user time: " << usertime() << '\n' ;
	o << "sys. time: " << systime() << '\n' ;
	return o << "real time: " << realtime() << std::endl ;
}

// Some arithmetic operator :
Timer& Timer::operator = (const Timer & T)
{
	ut = T.ut ; 
	st = T.st ; 
	rt = T.rt ;
	_count = T._count;
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
	Tmp._count = _count - T._count;
	return Tmp ;
}

const Timer Timer::operator - ()
{
	Timer Tmp ;
	Tmp.ut = -ut ;
	Tmp.st = -st ;
	Tmp.rt = -rt ;
	Tmp._count = - _count;
	return Tmp ;
}

const Timer Timer::operator + (const Timer & T)  const
{
	Timer Tmp ;
	Tmp.ut = ut + T.ut ;
	Tmp.st = st + T.st ;
	Tmp.rt = rt + T.rt ;
	Tmp._count = _count + T._count;
	return Tmp ;
}
 
}
#endif //__LINBOX_timer_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
