// ======================================================================= //
// Time-stamp: <09 Mar 07 17:57:48 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LAZY_PRODUCT__
#define __LAZY_PRODUCT__


namespace LinBox {

        // Lazy computation of the product of the moduli
    struct LazyProduct : public std::vector< Integer > {
        typedef std::vector< Integer > Father_t;
    protected:
        bool                _tobecomputed;
    public:

        LazyProduct() : Father_t(), _tobecomputed(true) {}

        void initialize(const Integer& i) {
            _tobecomputed = false;
            this->resize(0);
            this->push_back(i);
        }
            
        bool mulin(const Integer& i) {
            if (this->size()) {
                if (i != this->back()) {
                    this->push_back( i );
                    return _tobecomputed = true;
                } else {
                    return _tobecomputed;
                }
            
            } else {
                this->push_back( i );
                return _tobecomputed = false;
            }
        }
      
        bool mulin(const LazyProduct& i) {
            this->insert(this->end(), i.begin(), i.end());
            return _tobecomputed = (this->size()>1);
        }
      
        Integer & operator() () {
            if (_tobecomputed) {
                Father_t::const_iterator iter = this->begin();
                Father_t::iterator       prod = this->begin();
                for(++iter; iter != this->end(); ++iter)
                    *prod *= *iter;
                this->resize(1);
                _tobecomputed = false;
            }
            return this->back();
        }

        bool noncoprime(const Integer& i) const {
            Integer g;
            for(Father_t::const_iterator iter = this->begin(); iter != this->end(); ++iter)
                if ( gcd(g,i,*iter) > 1) return true;
            return false;
        }   
       
        friend std::ostream& operator<< (std::ostream& o, const LazyProduct& C) {
            o << "{";
            for(Father_t::const_iterator refs = C.begin();
                refs != C.end() ;
                ++refs )
                o << (*refs) << " " ;
            return o << "}";
        }
        
    };
    
}


#endif
