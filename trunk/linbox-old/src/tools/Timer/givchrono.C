
class RealChrono  {
   long _t, _ut;  // time  
public:
    enum { 
        MSPSEC = 1000000  // microsecond per second
    };

    inline void clear() { _t = 0; _ut = 0; }

    inline void start() {  
        struct timeval tmp2 ; 
        gettimeofday (&tmp2, 0) ;
        
            // real time 
        _t = tmp2.tv_sec;
        _ut = tmp2.tv_usec;
    }

    inline void stop() { 
        struct timeval tmp2 ;  
        gettimeofday (&tmp2, 0) ;
        
            // real time 
        _t = tmp2.tv_sec - _t;
        _ut = tmp2.tv_usec - _ut;
    }

   

    
    inline double now() {
            struct timeval tmp2 ;  
            gettimeofday (&tmp2, 0) ;
                // real time 
            return ( 
                ((double) tmp2.tv_sec - _t)*(double)MSPSEC
                + (double)tmp2.tv_usec - _ut ) ; 
    }
    
    
    friend inline ostream& operator<< (ostream& o, const RealChrono& BT) { 
        return o << (double)(BT._t + (double)BT._ut/(double)MSPSEC); }

};


