// ======================================================================= //
// (C) Linbox 2000
// Pair of long and T : struct { column index, value }
// Time-stamp: <07 Mar 00 18:12:28 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#ifndef _LIN_PAIR_H_
#define _LIN_PAIR_H_
#include <iostream.h>

// ---------------------------------------------------
//
/// Pair of long and T : struct { column index, value }
template<class T> class Pair {
public:
    typedef Pair<T> Self_t;
    typedef T                      Type_t;

//    ~Pair() {};
    Pair() {};

    Pair(const long jj, const T& val) :_j(jj),_value(val){};
    Pair(const Self_t& p) :_j(p._j),_value(p._value){};


    T getvalue() const { return _value; };

    long getindex() const { return _j; };
    long j() const { return _j; };
    
    T affect(const T& val) { return _value = val; };
    T change_value(const T& val) { return _value = val; };
   
    long change_j(const long jj) { return _j = jj; };      
    long change_index(const long jj) { return _j = jj; };      
            
            
    Self_t assign(const T& val) {
        _value = val;
        return *this;
    };      
            
    Self_t assign(const long jj, const T& val) {
        _value = val;
        _j = jj;
        return *this;
    };      
            
    long decr() { return --_j; };      
    long operator--() { return --_j; };      
    long operator--(int) { return _j--; };      
    long incr() { return ++_j; };      
    long operator++() { return ++_j; };      
    long operator++(int) { return _j++; };      
            
    friend inline istream& operator>> (istream& is, Pair<T>& a) {
        long jj;
        T val;
        is >> jj >> val;
        a._value=val; a._j=jj;
//         a = Pair<T>(jj,val);
        return is;
};
    
    friend inline ostream& operator<< (ostream& o, const Pair<T> a){
//         return o << a.j() << " " << a.getvalue()  ;
        return o << a._j << " " << a._value ;
};
    

private:
    long _j;
    T _value;
};



#endif // _LIN_PAIR_H_
