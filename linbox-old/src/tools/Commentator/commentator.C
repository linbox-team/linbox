#ifndef __commentary_C_
#define __commentary_C_

#include <math.h>
#include "commentator.h"

// -----------------------------------------------------
// Mathematical routines
// -----------------------------------------------------
double nroot(double a, long r, double precision) {
    long rm = r - 1 ;
    double c1 = double(rm)/double(r), c2 = a/double(r);
    double g = 1, pgr = 1, err = a - 1;
    
    while( err > precision ) {
        g = g*c1 + c2/pgr;
        pgr=pow(g,rm);
        err = a - pgr*g;
        if (err < 0)
            err = -err;
    }
    return g;
}

long isnpower(long& l, long a) {
    long r = 2;
    double g;
    while( (g = nroot(a, r, 0.1)) >= 2) {
        l = (long)floor(g);
        if (g-double(l) > 0.1)
            ++l;
        if (pow(l,r) == a)
            return r;
        ++r;
    }
    return 0;
}

    
// wrapper for use by Pascal (& C) code.
extern "C" 
Commentator* initializeCommentator(long timing, long homology)
{ return new Commentator(timing, homology); }

extern "C"
void startActivity(Commentator& C, char* id, char* msg, long msglevel, long msgclass)
{ C.start(id, msg, msglevel, msgclass); }

extern "C"
void stopActivity(Commentator& C, char* msg, long msglevel, long msgclass)
{ C.stop(msg, msglevel, msgclass); }

extern "C"
void activityReport(const Commentator& C, char* msg, long msglevel, long msgclass)
{ C.report(msg, msglevel, msgclass); }

extern "C"
void progressReport(Commentator& C, char* msg, long msglevel, long k, long n)
{ C.progress(msg, msglevel, k, n); }

extern "C"
long isPrinted(const Commentator& C, long msglevel, long msgclass)
{ return C.printed(msglevel, msgclass); }
#endif
