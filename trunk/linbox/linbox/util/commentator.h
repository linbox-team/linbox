/* -*- mode: c; style: linux -*- */

/* linbox/src/util/commentator.h
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

#ifndef __COMMENTATOR_H
#define __COMMENTATOR_H

#include <deque>
#include <stack>

#include "linbox/util/timer.h"

#ifndef MAX
#  define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

// =============================================
// 
// Message classes definitions
// 
// =============================================

// Min and Max
// initializing commentator with PRINT_EVERYTHING, everything is printed.
// initializing commentator with PRINT_NOTHING, nothing is printed.
#define PRINT_EVERYTHING 100000 
#define PRINT_NOTHING 0

enum EstimationMethod {
	BEST_ESTIMATE,
	POLY_ESTIMATE,
	EXPO_ESTIMATE
};

enum TimingDelay {
	SHORT_TIMING,
	LONG_TIMING
};

enum MessageClass {
	MSG_CLASS_INVALID,
	TIMING_MEASURE,
	TIMING_ESTIMATE,
	PARTIAL_RESULT,
	INTERNAL_WARNING,
	INTERNAL_ERROR,
	INTERNAL_DESCRIPTION
};

enum MessageLevel {
	LVL_ALWAYS = -1,
	LVL_IMP    =  0,
	LVL_NORMAL =  1,
	LVL_UNIMP  =  2,
	LVL_BLABLA =  10,
	LVL_NEVER  = (2*PRINT_EVERYTHING)
};

// n th root of a double 
#define _NROOT_PRECISION_ 0.0001
double nroot(double a, long r, double precision);
long isnpower(long& l, long a) ;

// =============================================
// 
// Commentator class
// 
// =============================================

class Commentator {
    public: 
	Commentator(long a = 0) : 
		gapTimingLevel(0), 
		gapHomologyLevel(0), 
		indentationStep(0), 
		estimatesDepth(0), 
		cnull(new nullstreambuf()) {};
    
        // When copying, Timers and Estimators are not needed.
	Commentator(const Commentator& C) : 
		gapTimingLevel(C.gapTimingLevel), 
		gapHomologyLevel(C.gapHomologyLevel), 
		indentationStep(C.indentationStep),  
		IDs(C.IDs), 
		estimatesDepth(C.estimatesDepth),
		Timers(stack<Timer>()), 
		Estimates(stack<Estimator>()), 
		cnull(new nullstreambuf())   {};

	Commentator& operator= (const Commentator& C) {
		if (this == &C) return *this;
		gapTimingLevel = C.gapTimingLevel;
		gapHomologyLevel = C.gapHomologyLevel;
		indentationStep = C.indentationStep;
		IDs = C.IDs;
		estimatesDepth = C.estimatesDepth;
		Timers = stack<Timer>();
		Estimates = stack<Estimator>();
		return *this;
	}       
   

        // Sets user printing levels
        // tab is the indentation unit
        // depth is the number of steps stored for time estimates 
        //       (should be greater than 2)
	Commentator(long gt, long gh, long tab=3, long depth=2, long estimation=0):
		gapTimingLevel(gt), 
		gapHomologyLevel(gh), 
		indentationStep(tab), 
		estimatesDepth(MAX(depth,2)), 
		estimationMethod(estimation), 
		cnull( new nullstreambuf() ) {}

        // Starting and Stoping an Activity
	void start(char* id, char* msg, long msglevel, long msgclass);
	void stop(char* msg, long msglevel, long msgclass, long time_type);

        // Reporting progress : 'k' out of 'n' 'msg' done.
	void progress(char* msg, long msglevel, long k, long n);

        // Reporting.
	void report(char* msg, long msglevel, long msgclass) const;

        // Testing effective printing with those arguments.
	long printed(long msglevel, long msgclass) const;

    private:
        // Internal sub routines deciding wether things are printed or not
	long intern_printed(char& information, long msglevel, long msgclass) const;
	long header(long msglevel, long msgclass) const ;
	void printheader(char information) const ;

        // Printing decision is made from those
	long gapTimingLevel;
	long gapHomologyLevel;

        // Indentation unit
	long indentationStep;
        // Activity id
	stack<char*> IDs;

        // Activity Estimator
	unsigned long estimatesDepth, estimationMethod;

        // Activity Timer
	stack<Timer> Timers;

	class StepsAndTime {
	public:
		StepsAndTime(long k, double t) : _time(t), _steps(k) {}
		long getsteps() { return _steps; }
		double gettime() { return _time; }
	private:
		double _time;
		long _steps;
	};
        
	typedef deque< StepsAndTime > Estimator;
	stack<Estimator> Estimates;

        // Null ostream prints nothing
	struct nullstreambuf : public streambuf {
		nullstreambuf() {};
		streampos seekoff(long long, ios::seek_dir, int) {return 0;}
		streampos seekpos(long long, int) {return 0;}
		streampos sys_seek(long long, ios::seek_dir) {return 0;}
		int showmanyc(void) {return 0;}
		void imbue(void *) {}
	};

	ostream cnull;

    public:
        // Starting, stoping and reporting Activity, using CommentatorStream.
	ostream& start(char* id, long msglevel, long msgclass);
	ostream& stop(long msglevel, long msgclass, long time_type);
	ostream& report(long msglevel, long msgclass);    
};

// =============================================
// 
// wrappers for use by Pascal (& C) code.
// 
// =============================================
extern "C" Commentator* initializeCommentator(long timing, long homology);

extern "C" void startActivity(Commentator& C, char* id, char* msg, long msglevel, long msgclass);

extern "C" void stopActivity(Commentator& C, char* msg, long msglevel, long msgclass);

extern "C" void activityReport(const Commentator& C, char* msg, long msglevel, long msgclass);

extern "C" void progressReport(Commentator& C, char* msg, long msglevel, long k, long n);

extern "C" long isPrinted(const Commentator& C, long msglevel, long msgclass);

#include "commentator.inl"
#include "commentator.C"

#endif // __commentary_H
