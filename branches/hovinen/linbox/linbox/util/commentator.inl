/* -*- mode: c; style: linux -*- */

/* linbox/src/util/commentator.inl
 * Copyright (C) 1999 B. David Saunders,
 *                    Jean-Guillaume Dumas
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#include <math.h>
#include <iostream.h>

#include "linbox/util/commentator.h"

// -----------------------------------------------------
// Mathematical routines
// -----------------------------------------------------
#ifndef __GIVMATHS__
#define __GIVMATHS__

#define polynomialize(k3,t3,c,a) { (t3) = (c)*pow((a),(k3));}

#define ispolynomial(c,a,k1,t1,k2,t2) { \
    (a) = nroot((t2)/(t1), (k2)-(k1), _NROOT_PRECISION_); \
    (c) = (t2)/pow((a),(k2)); \
}

#define isexponential(c,a,k1,t1,k2,t2) { \
    (a) = log((t2)/(t1)) / log((double)(k2)/(double)(k1)); \
    (c) = (t2)/exp((a)*log((double)(k2))); \
}

#define exponentialize(k3,t3,c,a) { (t3) = (c)*exp((a)*log(k3)); } 

#endif // __GIVMATHS__


// -------------------------------
// Inline members
// -------------------------------
inline void Commentator::start(char* id, char* msg, long msglevel, long msgclass)
{  
	IDs.push(id);
	if (strlen(msg) > 0) report(msg, msglevel, msgclass);
	Timer activ; activ.clear(); activ.start();
	Timers.push( activ );
	Estimator est; 
	Estimates.push( est );
}

inline ostream& Commentator::start(char* id, long msglevel, long msgclass)
{  
	IDs.push(id);
	Estimator est; 
	Estimates.push( est );
	Timer activ; activ.clear(); activ.start();
	Timers.push( activ );
	return report(msglevel, msgclass);
}

inline void Commentator::stop(char* msg, long msglevel, long msgclass, long time_type = SHORT_TIMING) 
{
	(Timers.top()).stop();
	if (header(msglevel, TIMING_MEASURE)) {
		if (time_type == LONG_TIMING) {
			cerr << "cpu  time : " << (Timers.top()).usertime() 
			     << " seconds" << endl;
			if (IDs.size() == 1) {
				header(msglevel, TIMING_MEASURE);
				cerr << "sys  time : " << (Timers.top()).systime() 
				     << " seconds" << endl;
				header(msglevel, TIMING_MEASURE);
				cerr << "real time : " << (Timers.top()).realtime() 
				     << " seconds" << endl;
			}
		} else {   
			cerr << (Timers.top()).realtime() 
			     << "s (" << (Timers.top()).usertime() << " cpu)" << endl;
		}
	}

	if (strlen(msg) > 0) report(msg, msglevel, msgclass);
	IDs.pop();
	Timers.pop();
	Estimates.pop();
}

inline ostream& Commentator::stop(long msglevel, long msgclass, long time_type = SHORT_TIMING) 
{
	(Timers.top()).stop();
	if (header(msglevel, TIMING_MEASURE)) {
		if (time_type == LONG_TIMING) {
			cerr << "cpu  time : " << (Timers.top()).usertime() 
			     << " seconds" << endl;
			if (IDs.size() == 1) {
				header(msglevel, TIMING_MEASURE);
				cerr << "sys  time : " << (Timers.top()).systime() 
				     << " seconds" << endl;
				header(msglevel, TIMING_MEASURE);
				cerr << "real time : " << (Timers.top()).realtime() 
				     << " seconds" << endl;
			}
		} else {
			cerr << (Timers.top()).realtime() 
			     << "s (" << (Timers.top()).usertime() << " cpu)" << endl;
		}
	}
	long isp = header(msglevel, msgclass);
	IDs.pop();
	Timers.pop();
	Estimates.pop();
	if (isp) 
		return cerr;
	else
		return cnull;
   
}

inline void Commentator::printheader(char information) const 
{
        long indent = indentationStep*IDs.size();
        cerr << "#I:" << information;
        for (long i = 0; i < indent; i++) cerr << " ";
        if (IDs.size()) cerr << IDs.top() << ": ";
        else cerr << " ";
}

inline long Commentator::header(long msglevel, long msgclass) const 
{
	long isp; char information;
	if ( (isp = intern_printed(information, msglevel, msgclass)) )
		printheader(information);
	return isp;
}




inline void Commentator::report(char* msg, long msglevel, long msgclass) const
{
	if (header(msglevel, msgclass))
		cerr << msg << endl;
}

inline ostream& Commentator::report(long msglevel, long msgclass)
{
	if (header(msglevel, msgclass)) 
		return cerr;
	else
		return cnull;
}

inline void Commentator::progress(char* msg, long msglevel, long k, long n)
{
	if ( header( msglevel, TIMING_ESTIMATE) ) {
		Timer star = Timers.top();
		(Timers.top()).stop();
		double t3 = (Timers.top()).realtime();
		if (k) Estimates.top().push_back( StepsAndTime( k, t3) );
		if (Estimates.top().size() >= 3) {

			Estimator::iterator Etop = Estimates.top().begin();

			double t1 = (*Etop).gettime();
			long k1 = (*Etop).getsteps();

			Etop += Estimates.top().size()/2;

			double t2 = (*Etop).gettime();
			long k2 = (*Etop).getsteps();

			double cp, ap, ce, ae;
			double guess1, guess2, guess;

			switch(estimationMethod) {
			case POLY_ESTIMATE:
				ispolynomial(cp, ap, k1, t1, k2, t2);
				polynomialize(n, guess, cp, ap);
				break;
			case EXPO_ESTIMATE:
				isexponential(ce, ae, k1, t1, k2, t2);
				exponentialize(n, guess, ce, ae);
				break;
			default:
				ispolynomial(cp, ap, k1, t1, k2, t2);
				polynomialize(k, guess1, cp, ap);
				isexponential(ce, ae, k1, t1, k2, t2);
				exponentialize(k, guess2, ce, ae);
				guess1 = t3 - guess1; if (guess1 < 0) guess1 = - guess1;
				guess2 = t3 - guess2; if (guess2 < 0) guess2 = - guess2;
				if (guess2 > guess1) {
					polynomialize(n,guess,cp,ap);
				} else {
					exponentialize(n,guess,ce,ae); 
				}
			}
			cerr << k << " / " << n << " " << msg << "; Elaps: " << t3 << "; Expect:" << guess-t3 << "." << endl;
		} else {
			cerr << k << " out of " << n << " " << msg << " done." << endl;
		}
        

		if (Estimates.top().size() > estimatesDepth)
			Estimates.top().pop_front();
		Timers.top() = star;
	}
}




inline long Commentator::intern_printed(char& information, long msglevel, long msgclass) const 
{
	if (msglevel == -1) {
		information = 'U';
		return 1;
	}
	switch(msgclass) {
        case INTERNAL_ERROR:
		information = 'E';
		return 1;
        case INTERNAL_WARNING:
		information = 'W';
		return gapHomologyLevel > msglevel;
        case PARTIAL_RESULT:
		information = 'R';
		return gapHomologyLevel > msglevel + (long)IDs.size() - 1;
        case INTERNAL_DESCRIPTION:
		information = 'D';
		return gapHomologyLevel > msglevel + (long)IDs.size();
        case TIMING_MEASURE:
		information = 'T';
		return gapTimingLevel > msglevel + (long)IDs.size() - 1;
        case TIMING_ESTIMATE:
		information = 'S';
		return gapTimingLevel > msglevel + (long)IDs.size();
        default:
		information = 'U';
		return 0;
	}
    
}

inline long Commentator::printed(long msglevel, long msgclass) const 
{
	char c; return intern_printed(c, msglevel, msgclass);
}
