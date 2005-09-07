// ======================================================================= // (C) Linbox 2000
// Pair of I and T : struct { column index, value }
// Time-stamp: <19 Sep 03 11:00:25 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= 
#ifndef _LIN_PAIR_H_
#define _LIN_PAIR_H_
#include <iostream>

// ---------------------------------------------------
//
/// Pair of I and T : struct { column index, value }
template<class T, class I = unsigned long> class Pair {
public:
    typedef Pair<T, I> Self_t;
    typedef T                      Type_t;

//    ~Pair() {};
    Pair() {};

    Pair(const I jj, const T& val) :_j(jj),_value(val){};
    Pair(const Self_t& p) :_j(p._j),_value(p._value){};


    T getvalue() const { return _value; };

    I getindex() const { return _j; };
    I j() const { return _j; };
    
    T affect(const T& val) { return _value = val; };
    T change_value(const T& val) { return _value = val; };
   
    I change_j(const I jj) { return _j = jj; };      
    I change_index(const I jj) { return _j = jj; };      
            
            
    Self_t assign(const T& val) {
        _value = val;
        return *this;
    };      
            
    Self_t assign(const I jj, const T& val) {
        _value = val;
        _j = jj;
        return *this;
    };      
            
    I decr() { return --_j; };      
    I operator--() { return --_j; };      
    I operator--(int) { return _j--; };      
    I incr() { return ++_j; };      
    I operator++() { return ++_j; };      
    I operator++(int) { return _j++; };      
            
    friend inline std::istream& operator>> (std::istream& is, Pair<T, I>& a) {
        I jj;
        T val;
        is >> jj >> val;
        a._value=val; a._j=jj;
//         a = Pair<T, I>(jj,val);
        return is;
};
    
    friend inline std::ostream& operator<< (std::ostream& o, const Pair<T, I> a){
//         return o << a.j() << " " << a.getvalue()  ;
        return o << a._j << " " << a._value ;
};
    

private:
    I _j;
    T _value;
};



#endif // _LIN_PAIR_H_
