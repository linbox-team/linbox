// ================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// Time-stamp: <30 May 00 14:53:51 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


#ifndef __Symmetric_Lanczos_H__
#define __Symmetric_Lanczos_H__
#include <LinBox/lin_rand.h>
#include <LinBox/lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class Symmetric_Lanczos : public Base_BB_Container< BlackBoxDomain, Vecteur > {
public:
    Symmetric_Lanczos() {} 

    Symmetric_Lanczos(BlackBoxDomain_t * D, const Vecteur& u0) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) {    
        init(u0); self_init();
    }
    
    Symmetric_Lanczos(BlackBoxDomain_t * D, RandIter& g ) 
            : Base_BB_Container< BlackBoxDomain, Vecteur>(D) { 
        init(g); self_init();
    }

    bool check () {
        bool zerow = _domain.iszero(mu);

        if (even)
            for( typename Vecteur::const_iterator i = u.begin(); (i != u.end()) && (zerow); ++i )
                zerow = _domain.iszero( *i );
        else
            for( typename Vecteur::const_iterator i = w.begin(); (i != w.end()) && (zerow); ++i )
                zerow = _domain.iszero( *i );
    

        return zerow;    
    }
    
protected:

    void self_init() {
        _BB_domain->Apply(v, u);
        _domain.assign(lambda,_domain.one);
        _domain.assign(nu,_domain.zero);
        w.resize(u.size());
    }
    
    Vecteur w;
    Type_t alpha, beta, mu, nu, lambda;

    void _launch () {
        if (even) {
            even = 0;
            
            DOTPROD(alpha, u, v);
            _domain.divin( alpha, _value );
        
                // mu_i+1 <-- nu_i - mu_i * alpha_i+1
            _domain.addin( _domain.mul(mu,lambda,alpha) ,nu); 
//             _domain.write(cerr << "mu: " , mu) << endl;


            _domain.negin( alpha );
            AXPY(w,alpha,u,v);
            if (_domain.iszero(_value) ) return ;

//             cerr << "w:=[ "; for(long k=0; k<(w.size()-1); ++k)
//                 _domain.write(cerr, w[k]) << ", ";
//             _domain.write(cerr , w[w.size()-1]) << " ];" << endl;
            

            alpha = _value; // temporary

            DOTPROD(_value,w,w);
            _domain.div(beta, _value, alpha) ;
        
                // nu_i+1 <-- - mu_i * beta_i+1 
            _domain.mul(nu,lambda,beta);

            _BB_domain->Apply(v,w);

            _domain.negin( beta );
            AXPYIN(v,beta,u);

            _domain.neg(lambda,mu);

        } else {
            even = 1;

            DOTPROD(alpha, w, v);
            _domain.divin( alpha, _value );
        
                // mu_i+1 <-- nu_i - mu_i * alpha_i+1
            _domain.addin( _domain.mul(mu,lambda,alpha) ,nu); 
//             _domain.write(cerr << "mu: " , mu) << endl;


            _domain.negin( alpha );
            AXPY(u,alpha,w,v);
            if (_domain.iszero(_value) ) return ;

//             cerr << "w:=[ "; for(long k=0; k<(u.size()-1); ++k)
//                 _domain.write(cerr, u[k]) << ", ";
//             _domain.write(cerr , u[u.size()-1]) << " ];" << endl;
            

            alpha = _value; // temporary

            DOTPROD(_value,u,u);
            _domain.div(beta, _value, alpha) ;
        
                // nu_i+1 <-- - mu_i * beta_i+1 
            _domain.mul(nu,lambda,beta);

            _BB_domain->Apply(v,u);

            _domain.negin( beta );
            AXPYIN(v,beta,w);

            _domain.neg(lambda,mu);
        }
    }
    
    void _wait () {}

};

#endif // __Symmetric_Lanczos_H__
