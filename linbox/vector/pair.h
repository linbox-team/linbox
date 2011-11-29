// ======================================================================= 
// (C) Linbox 2000
// Pair of I and T : struct { column index, value }
// Time-stamp: <11 Sep 08 14:21:05 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= 
#ifndef _LIN_PAIR_H_
#define _LIN_PAIR_H_
#include <iostream>

// ---------------------------------------------------
//
/// Pair of I and T : struct { column index, value }
template<class I, class T> class Pair {
public:
    typedef Pair<I, T> Self_t;
    typedef T                      Type_t;
    typedef I                      first_type;
    typedef T                      second_type;
    
//    ~Pair() {};
    Pair() {};

    Pair(const I jj, const T& val) :first(jj),second(val){};
    Pair(const Self_t& p) :first(p.first),second(p.second){};


    T getvalue() const { return second; };

    I getindex() const { return first; };
    I j() const { return first; };
    
    T affect(const T& val) { return second = val; };
    T change_value(const T& val) { return second = val; };
   
    I change_j(const I jj) { return first = jj; };      
    I change_index(const I jj) { return first = jj; };      
            
            
    Self_t assign(const T& val) {
        second = val;
        return *this;
    };      
            
    Self_t assign(const I jj, const T& val) {
        second = val;
        first = jj;
        return *this;
    };      
            
    I decr() { return --first; };      
    I operator--() { return --first; };      
    I operator--(int) { return first--; };      
    I incr() { return ++first; };      
    I operator++() { return ++first; };      
    I operator++(int) { return first++; };      
            
    friend inline std::istream& operator>> (std::istream& is, Pair<I, T>& a) {
        is >> a.first >> a.second;
        return is;
    }
    
    friend inline std::ostream& operator<< (std::ostream& o, const Pair<I, T> a){
        return o << a.first << " " << a.second ;
    }
    

public:
    I first;
    T second;
};



#endif // _LIN_PAIR_H_
