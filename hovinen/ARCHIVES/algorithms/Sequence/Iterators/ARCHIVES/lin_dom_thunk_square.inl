template<class Domain, class Vecteur>
class Launch<Domain, Vecteur, Square> {
public:
    Launch() {}
    
    Launch(const Domain& d, long n_row, long n_col) : _domain(d), even(1) {
        u = Vecteur(n_col);
        for(long i=n_col;i--;)
            _domain.random(u[i]);
//         w.copy(u);
        w = u;
        v.resize(n_col);
    }        

    Launch(const Domain& d, long n_row, long n_col, const Vecteur& uu) : _domain(d), even(1) {
        u = uu;
        w = uu;
        v.resize(n_col);
    }        

    template< class Type_t, class BB >
    void operator() (Type_t& _value, BB * bbdom) {
        cerr << " Even : " << even << endl;
        if (even) {
            bbdom->Apply( v, w);  // GV
            DOTPROD(_value,u,v);       // GV 
            even = 0;
        } else {
            bbdom->Apply( w, v);  // GV
            DOTPROD(_value,u,w);       // GV
            even = 1;
        }  
    }
    
private:
    Domain _domain;
// -----------------------------------------------
// dot product
// Preconditions : u.size() > 0, v.size() > 0 ....
// -----------------------------------------------
    template <class Type_t, class A1, class A2>
    void DOTPROD(Type_t& coeff, const A1& u, const A2& v) {
        _domain.mul(coeff, u[0], v[0]);
        for (Indice k=u.size()-1 ; k>0; --k)
            _domain.axpyin(coeff, u[k], v[k]);
    }

    Vecteur u, v, w;
    bool even;
};


