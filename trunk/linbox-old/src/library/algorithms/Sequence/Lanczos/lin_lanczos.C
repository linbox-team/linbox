// ======================================================================= //
// Linbox project 1999
// Lanczos
// Time-stamp: <26 May 00 17:52:40 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef _LIN_LANCZOS_C_
#define _LIN_LANCZOS_C_

template<class Sequence> class Lanczos {
public:
    typedef          Sequence                 Sequence_t;
    typedef typename Sequence::Domain_t       Domain_t;
    typedef typename Sequence::Type_t         Type_t;
    typedef          Lanczos< Sequence >    Self_t;
private:
    Sequence_t * _container;
    Domain_t _domain;

public:
        //-- Constructors
    Lanczos() 
            : _container(), 
              _domain() 
        {}

    Lanczos(const Self_t& M) 
            : _container(M._container), 
              _domain(M._domain) 
        {}

    Lanczos(Sequence_t * D) 
            : _container(D), 
              _domain(D->getdomain()) 
        {}
  
        //-- Principal method
    long& operator() (long& r) {
        return sym_lanczos(r);
    };
    
        //-- Domains access
    const Domain_t& getdomain() const { return _domain; }
    Sequence_t * getsequence() const { return _container; }
    

private:

    long& sym_lanczos (long& r) {
        const long END = _container->size() >> 1;
        r = -1;

        typename Sequence_t::const_iterator _iter( _container->begin() );

        for(;(! _domain.iszero(*_iter)) && (r <= END) ; ++r, ++_iter) { }

        if (_container->check() )
            cerr << "Lanczos : " << r << endl;
        else
            cerr << "Échec   : " << r << endl;
    }   

};


    
#endif // _LIN_LANCZOS_C_
