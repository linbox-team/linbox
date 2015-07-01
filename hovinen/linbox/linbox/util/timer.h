/* -*- mode: c; style: linux -*- */

/* linbox/src/util/timer.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * This file implements the C++ interface to commentators (for 
 * providing runtime commentary to the user)
 */

#ifndef __TIMER_H
#define __TIMER_H

#include <iostream.h>

// class RealTimer; class SysTimer; class UserTimer;

class BaseTimer { 
    public:
	enum { 
		MSPSEC = 1000000  // microsecond per second
	};

	// -- Clear timer :
	inline void clear() { _t = 0; }

	// -- total amount of second spent 
	inline double time() const { return _t; }

	// -- Return a value to initialize random generator
	static long seed();

	// -- basic methods:
	ostream& print( ostream& ) const;
       
	// -- Some arithmetic operators to compute cumulative time :
	BaseTimer& operator = (const BaseTimer & T) ;
	const BaseTimer operator - (const BaseTimer & T)  const;
	const BaseTimer operator - () ;
	const BaseTimer operator +  (const BaseTimer & T)  const;
	BaseTimer& operator += (const BaseTimer & T) { return *this = *this + T; };
	BaseTimer& operator -= (const BaseTimer & T) { return *this = *this - T; };

    public:
	double _t;  // time  
};
inline ostream& operator<< (ostream& o, const BaseTimer& BT)
{ return BT.print(o);}


class RealTimer : public BaseTimer {
    public:
	inline RealTimer( const BaseTimer& BT ): BaseTimer(BT) {};
	inline RealTimer( ){};
	void start();
	void stop();
};


class UserTimer : public BaseTimer {
    public:
	inline UserTimer( const BaseTimer& BT ): BaseTimer(BT) {};
	inline UserTimer( ) {};
	void start();
	void stop();
};


class SysTimer : public BaseTimer {
    public:
	inline SysTimer( const BaseTimer& BT ): BaseTimer(BT) {};
	inline SysTimer( ) {};
	void start();
	void stop();
};


class Timer {
	public :

		// Clear timer :
		void clear(); 

	// Start timer
	void start();

	// Stop timer 
	void stop();

	// total amount of second spent in user mode
	double usertime() const { return ut.time(); }

	// total amount of second spent in system mode
	double systime () const { return st.time(); }

	// real total amount of second spent.  
	double realtime () const { return rt.time(); }

	// retourne une petite graine
	// long seed() const { return RealTimer::seed(); }

	// Some arithmetic operators to compute cumulative time :
	Timer& operator = (const Timer & T) ;
	const Timer operator - (const Timer & T)  const;
	const Timer operator - () ;
	const Timer operator + (const Timer & T)  const;
	Timer& operator += (const Timer & T) { return *this = *this + T; };
	Timer& operator -= (const Timer & T) { return *this = *this - T; };


	// -- methods :
	ostream& print( ostream& ) const;

    public:
	RealTimer rt;
	UserTimer ut;
	SysTimer  st;
};
// inline ostream& operator<<( ostream& o, const Timer& T)
// { return T.print(o);}

inline ostream& operator<<( ostream& o, const Timer& T)
{ return o << T.realtime() << "s (" << T.usertime() << " cpu)"; }

#include "linbox/util/timer.C"

#endif 
