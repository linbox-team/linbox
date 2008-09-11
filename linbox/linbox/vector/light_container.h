// =================================================================== //
// LightContainer : std::vector like container
// Time-stamp: <11 Sep 08 14:59:12 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef __Light_Container__ 
#define __Light_Container__

#include <iostream>
#include <cstdlib>
#include "linbox/util/contracts.h"

template<typename Elem> struct LightContainer {
private:
    typedef LightContainer<Elem> Self_t;
    size_t allocated;
    Elem * _container;
    Elem * _finish;
public:
    typedef Elem value_type;
    typedef Elem* iterator;
    typedef const Elem* const_iterator;

    LightContainer() : allocated(2), _container(new Elem[allocated]), _finish(_container) { 
        ENSURE( (allocated == 2) && (size() == 0) );

    }
    LightContainer(size_t s) : allocated(s), _container(new Elem[s]), _finish(_container+s) {
        
        ENSURE( (allocated == s) && (size() == s) );
        ENSURE( allocated >= size() );
     
    }

    Self_t& operator=(const Self_t& v) {
        this->resize(v.size());
        const_iterator vit = v.begin();
        for(iterator tit = begin(); tit != end(); ++tit, ++vit)
            *tit = *vit;
        return *this;
    }

    LightContainer(const Self_t& v) : allocated(v.allocated) {
        _container = new Elem[allocated];
        _finish = _container +v.size();
        const_iterator vit = v.begin();
        for(iterator tit = begin(); tit != end(); ++tit, ++vit)
            *tit = *vit;
    }

    void reserve(size_t s) { 
        STATE( size_t oldsize = size() );
        reallocate(s, size() );
        ENSURE( (allocated >= s) && (size() == oldsize) && ( allocated >= size() ) );
    }

    size_t size() const { return size_t(_finish - _container); }

    Elem& operator[](size_t i) { 
        REQUIRE( (i >= 0) && (i < allocated) && (i < size()) );
        return _container[i]; }

    const Elem& operator[](size_t i) const { 
        REQUIRE( (i >= 0) && (i < allocated) && (i < size()) );
        return _container[i]; }

    void clear() { _finish = _container; 
            ENSURE( (size() == 0) );

    }
    void resize(size_t s) {
        if (s>allocated) reallocate( s+(s>>1), s );
        else _finish = _container + s;
        ENSURE( allocated >= size() );
    }
    iterator begin() { return _container; }
    iterator end() { return _finish; }
    Elem& front() { return *_container; }
    const_iterator begin() const { return const_iterator(_container); }
    const_iterator end() const { return const_iterator(_finish); }

    void push_back(const Elem& c) {
        STATE( size_t  oldsize = size() );
        if (size() == allocated) reserve(allocated+(allocated>>1));
        *(_finish) = c; ++_finish; 

        ENSURE( size() == (oldsize+1) );
        ENSURE( allocated >= size() );
    }

    ~LightContainer() { delete[] _container; }


    iterator insert(iterator pos, const Elem& c) {
        REQUIRE( (pos-begin()) <= (end()-begin()) );
        REQUIRE( (pos-begin()) >= 0 );
        STATE( size_t oldsize = size() );
        if (pos == _finish) {
            push_back(c);
        } else {
            insert(pos+1, *(pos));
            *(pos) = c;
        }
        ENSURE( size() == oldsize+1 );
        ENSURE( allocated >= size() );
        return pos;
    }
    
    

    friend std::ostream& operator<< (std::ostream& o, const Self_t& C) {
        o << '[';
        for(const_iterator refs =  C.begin();
            refs != C.end() ;
            ++refs )
            o << (*refs) << ' ';
        return o << ']';
    }


protected:
    void reallocate(size_t s, size_t endc) {
        REQUIRE( (s >= endc) );
        if (allocated < s) {
            Elem * futur = new Elem[s];
            for(size_t i=0; (i<s) && (i < allocated); ++i)
                futur[i] = _container[i];
            size_t olds = size();
            delete [] _container;
            _container = futur;
            _finish = _container + olds;
            allocated = s;
        }
        _finish = _container + endc;
        ENSURE( allocated >= size() );
    }

};  

#endif
