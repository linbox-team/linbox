// ======================================================================= //
// (C) 1999 The Linbox group
// Density : - array with incr, decr ...
//           - initialization with 0
//	     - resize is inefficient
//
// Time-stamp: <17 Mar 00 18:11:08 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#ifndef _LIN_DENSITY_H_
#define _LIN_DENSITY_H_

#ifdef _GPDEBUG_
#include <iostream.h>
#define _TESTDSIZE(i) ((i >= 0) && (_dsize>0) && (i < _dsize))
#define _GPSEGFAULT(j,str) { cerr << "GP " << str << " Seg fault. Got " << j << " and size is " << _dsize << endl; return -1;}
#else
#define _TESTDSIZE(i) (1)
#define _GPSEGFAULT(j,str) {}
#endif


class Density {
    typedef int Internal;
public:
    Density( )  : _density(0), _dsize(0) {}
    Density( const Internal d ) : _density(new Internal[d]), _dsize(d) {
        initialize();
    }

    ~Density() {
        deallocate();
    }
         
    void initialize() {
        for(Internal i=_dsize;i--;)
            _density[i] = 0;
    }       
   
    void resize(const Internal i) {
        if (_dsize != i) {
            deallocate();
            reallocate(i);
            initialize();
        }
    }        

// Elements access
    Internal getvalue( const Internal i) const {
        if ( _TESTDSIZE(i) ) return _density[i]; else _GPSEGFAULT(i,"get") ;
    }

    Internal incr( const Internal i) {
        if ( _TESTDSIZE(i) ) return (_density[i] += 1); else _GPSEGFAULT(i,"incr");
    }
  
    Internal decr( const Internal i) {
        if ( _TESTDSIZE(i) ) return (_density[i] -= 1); else _GPSEGFAULT(i,"decr");
    }
    
    Internal init( const Internal i) {
        if ( _TESTDSIZE(i) ) return (_density[i] = 0); else _GPSEGFAULT(i,"init");
    }

private:
    void reallocate(const Internal i) {
        _dsize = i;
        _density = new Internal(i);
    }       
   
    void deallocate() {
        if (_dsize) delete [] _density;
    }       
   
    Internal * _density;
    Internal _dsize;
};
    

#endif 
