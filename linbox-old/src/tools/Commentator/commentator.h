/****************************************************************************
 ** Time-stamp: <07 Feb 02 21:34:32 Jean-Guillaume.Dumas@imag.fr> 
 **
 *W  commentator.h  
 *A  B. D. Saunders, J.-G. Dumas
 **
 *H  @(#)$Id: commentator.h,v 0.0 1999/08
 **
 *Y  Copyright (C)  1999,  Us guys
 **
 **  This file implements the C++ interface to commentators (for 
 **  providing runtime commentary to the user)
 */

#ifndef __commentary_H
#define __commentary_H
#include <deque.h>
#include <stack.h>
#include <streambuf.h>
#include <LinBox/givtimer.C>

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

// Estimation methods
#define BEST_ESTIMATE 0
#define POLY_ESTIMATE 1
#define EXPO_ESTIMATE 2


// Timing display
#define SHORT_TIMING 0
#define LONG_TIMING  1


// Msg Class
#define TIMING_MEASURE          1
#define TIMING_ESTIMATE         2
#define PARTIAL_RESULT          3
#define INTERNAL_WARNING        4
#define INTERNAL_ERROR          5
#define INTERNAL_DESCRIPTION    6

// Msg Level
#define LVL_ALWAYS      -1
#define LVL_IMP         0
#define LVL_NORMAL      1
#define LVL_UNIMP       2
#define LVL_BLABLA      10
#define LVL_NEVER       (2*PRINT_EVERYTHING)


// n th root of a double 
#define _NROOT_PRECISION_ 0.0001
double nroot(double a, long r, double precision);
long isnpower(long& l, long a) ;

// MAX macro
#ifndef GIVMAX
#define GIVMAX(a,b) ((a)<(b) ? b:a)
#endif

   
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
            estimatesDepth(GIVMAX(depth,2)), 
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

#endif __commentary_H


