/* linbox/field/ntl-z_pE.h
 * Copyright (C) 2004  Pascal Giorgi
 * Copyright (C) 2011 LinBox
 *
 * Written by  Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * Modified by W. J. Turner <wjturner@acm.org>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file field/ntl/ntl-lzz_pE.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_lzz_pe_H
#define __LINBOX_field_ntl_lzz_pe_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pE.h>
#include <time.h>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include <givaro/zring.h>
#include "linbox/field/field-traits.h"
#include "linbox/field/hom.h"


#include "linbox/integer.h"

namespace Givaro
{
	template<>
	NTL::zz_pE& Caster(NTL::zz_pE &x, const Integer &y)
	{
		x=NTL::to_zz_pE(static_cast<int64_t>(y));
		return x;
	}
	template<>
	NTL::zz_pE& Caster(NTL::zz_pE &x, const double &y)
	{
		x=NTL::to_zz_pE(static_cast<long>(y));
		return x;
	}

	template<>
	Integer& Caster (Integer& x, const NTL::zz_pE &y) {
		NTL::zz_pX poly = rep(y);
		Integer base = static_cast<int64_t>(NTL::zz_p::modulus());
		long i = deg(poly)+1;
		x = 0;
		for( ; i-- ; ) {
			x *= base;
			x +=  int64_t(NTL::to_long(rep(coeff(poly, i))));
		}
		return x;
	}
} // namespace Givaro



// Namespace in which all LinBox library code resides
namespace LinBox
{

        //! use ZZ_pEBak mechanism too ?
	class NTL_zz_pE_Initialiser {
		NTL::zz_pX _f;
	public :
		NTL_zz_pE_Initialiser( const Integer & p, const Integer & k) {
			NTL::zz_p::init( (int64_t) p);
			NTL::zz_pX irredPoly = NTL::BuildIrred_zz_pX ((int64_t) k);
			NTL::zz_pE::init(irredPoly);
			_f = irredPoly;
		}
		
		NTL_zz_pE_Initialiser( const Integer & p, const NTL::zz_pX & f) {
			NTL::zz_p::init((int64_t) p);
			NTL::zz_pE::init(f);
			_f = f;
		}
		
		NTL_zz_pE_Initialiser(const NTL::zz_pX &f) {
			NTL::zz_pE::init(f);
			_f = f;
		}
		
		const NTL::zz_pX& modulus() const {
			return _f;
		}

            // template <class ElementInt>
            // NTL_zz_pE_Initialiser(const ElementInt& d) {
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
            // }

            // NTL_zz_pE_Initialiser (const NTL::ZZ& d) {
			// NTL::ZZ_p::init(d);
            // }

	};





        /*! @brief zz_pE
         * Define a parameterized class to easily handle Givaro::ZRing<NTL::zz_pE> field
         */

        /// \brief for large cardinality, small prime.  \ingroup field
	class NTL_zz_pE : public NTL_zz_pE_Initialiser, public Givaro::UnparametricOperations<NTL::zz_pE> {
	public:
		typedef NTL::zz_pE Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mOne ;

		NTL_zz_pE (const integer &p, const integer &k) :
                NTL_zz_pE_Initialiser(p,k),Father_t ()
        {
        	init(const_cast<Element &>(zero), 0);
		 	init(const_cast<Element &>(one), 1);
		 	init(const_cast<Element &>(mOne), p-1);
		}
        
		NTL_zz_pE (const integer &p, const NTL::zz_pX &f) :
                NTL_zz_pE_Initialiser(p,f), Father_t()
        {
            init(const_cast<Element &>(zero), 0);
		 	init(const_cast<Element &>(one), 1);
		 	init(const_cast<Element &>(mOne), p-1);
        }
        
        NTL_zz_pE(const NTL_zz_pE &F) :
        	NTL_zz_pE_Initialiser(F.modulus()), Father_t(),
        	zero(NTL::to_zz_pE(0)), one(NTL::to_zz_pE(1)), mOne(-one)
		{
        }

		Element& random (Element& x) const
            {
                NTL::random(x);
                return x;
            }


		bool isZero (const Element& a) const
            {
                return NTL::IsZero(a);
            }


		bool isOne (const Element& a) const
            {
                return NTL::IsOne(a);
            }

		bool isUnit (const Element& x) const
            {
            	if (isZero(x)) {
            		return false;
            	} 
            	
            	NTL::zz_pX g, tmp;
            	tmp = NTL::conv<NTL::zz_pX>(x);
            	NTL::GCD(g, tmp, modulus());
            	
            	return g == 1;
            }

		bool isMOne (const Element& x) const
            {
                Element y ; neg(y,x);
                return isOne(y);
            }

		NTL::zz_pX & init(NTL::zz_pX & f, integer n, int e) const
            {	integer base = characteristic();
			NTL::zz_pX x; SetCoeff(x, 1);
			if ( n != 0 and e > 0 )
				f = static_cast<int64_t>(n%base) + x*init(f, n/base, e-1);
			return f;
            }

		Element & init(Element & x, integer n = 0) const
            {   // assumes n >= 0.
                int e = exponent();
                n %= cardinality();
                init(const_cast<NTL::zz_pX &>(rep(x)), n, e);
                    //write(std::cout << "init-ed ", x) << std::endl;
                return x;
            }
        
        // documentation of NTL::conv:
        // http://www.shoup.net/ntl/doc/conversions.txt
        // XXX = long, ZZ, ZZ_p, ZZ_pE, ZZ_pX
        template<class XXX>
        Element &init(Element &x, const XXX &y) const {
        	x = NTL::conv<NTL::zz_pE>(y);
        	return x;
        }

        integer & convert(integer & x, const Element & y) const
            {
                NTL::zz_pX poly = rep(y);
                integer base = characteristic();
                long i = deg(poly)+1;
                x = 0;
                for( ; i-- ; ) {
                    x *= base;
                    x +=  static_cast<int64_t>(NTL::to_long(rep(coeff(poly, i))));
                }
                return x;
            }

		integer& characteristic (integer &c) const
            {	return c = characteristic(); }
		integer characteristic () const
            {
                return static_cast<int64_t>(NTL::zz_p::modulus());
            }


		integer& cardinality(integer& c) const
            {	return c = cardinality(); }

		int exponent() const
            {
                return (int) NTL::zz_pE::degree();
            }

		integer cardinality() const
            {
                int32_t e = exponent();
                int64_t p = characteristic();
                integer c(1);
                for(int32_t i = 0; i < e; ++i)
                    c *= p;
                return c;
            }

		static inline integer maxCardinality()
            { return NTL_zz_p::maxCardinality(); }


		Element& inv(Element& x, const Element& y) const
            {
                x=one/y;
                return x;
            }

		Element& invin(Element& x) const
            {
                x=one/x;
                return x;
            }
        
        Element& div(Element &x, const Element &y, const Element &z) const {
        	NTL::zz_pX g, zx;
        	conv(zx, z);
        	NTL::GCD(g, zx, modulus());
        	
        	NTL::zz_pE zg;
        	conv(zg, zx / g);
        	x = NTL::conv<NTL::zz_pE>(NTL::conv<NTL::zz_pX>(y) / g);
        	x /= zg;
        	
        	return x;
        }
        
        Element& divin(Element &x, const Element &y) const {
        	Element r;
        	div(r, x, y);
        	return x = r;
        }

        bool isDivisor(const Element &x, const Element &y) const {
        	if (isZero(y)) {
        		return false;
        	}
        	
        	if (isUnit(y)) {
        		return true;
        	}
        	
        	NTL::zz_pX a, b;
        	conv(a, x);
        	conv(b, y);
        	
        	NTL::zz_pX ga, gb;
        	NTL::GCD(ga, a, modulus());
        	NTL::GCD(gb, b, modulus());
        	
        	NTL::zz_pX r;
        	r = ga % gb;
        	return isZero(NTL::conv<Element>(r));
        }
		
		Element& gcd(Element &g, const Element &a, const Element &b) const {
			NTL::zz_pX r1, r2, x, y;
			conv(x, a);
			conv(y, b);
			
			NTL::GCD(r1, x, modulus());
			NTL::GCD(r2, y, r1);
			
			return g = NTL::conv<Element>(r2);
		}
		
		Element& gcdin(Element &a, const Element &b) const {
			NTL::zz_pE g;
			gcd(g, a, b);
			return a = g;
		}
		
		Element& dxgcd(Element &g, Element &s, Element &t, Element &u, Element &v, const Element &a, const Element &b) const {
			NTL::zz_pX gx, sx, tx, ux, vx, ax, bx;
			conv(ax, a);
			conv(bx, b);
			
			NTL::XGCD(gx, sx, tx, ax, bx);
			
			ux = ax / gx;
			vx = bx / gx;
			
			conv(g, gx);
			conv(s, sx);
			conv(t, tx);
			conv(u, ux);
			conv(v, vx);
			
			return g;
		}

		std::istream& read(std::istream& is, Element& x) const
            {
                long tmp = 0;
                is>>tmp;
                x=NTL::to_zz_pE(tmp);
                return is;
            }
            
		std::ostream& write( std::ostream& os ) const
            {
                return os << "Polynomial quotient ring using NTL::zz_pE";
            }
            
		std::ostream& write( std::ostream& os, const Element& x) const {
			os << x;
			return os;
		}
	}; // end of class NTL_zz_pE






	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::zz_pE> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::zz_pE> {
	public:
		typedef NTL::zz_pE Element;
        typedef size_t Residu_t;

		UnparametricRandIter<NTL::zz_pE>(const NTL_zz_pE & F ,
                                         const uint64_t seed = 0,
                                         const Residu_t& size = 0
                                         ) :
                _size(size), _seed(seed), _ring(F)
            {
                if(_seed == 0)
                    NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(std::time(nullptr))));
                else
                    NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(_seed)));
            }

        const NTL_zz_pE& ring() const { return _ring; }
		UnparametricRandIter<NTL::zz_pE>(const UnparametricRandIter<NTL::zz_pE>& R) :
                _size(R._size), _seed(R._seed), _ring(R._ring)

            {
                if(_seed == 0)
                    NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(std::time(nullptr))));
                else
                    NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(_seed)));
            }

		Element& random (Element& x) const
            {
                NTL::random(x);
                return x;
            }

	protected:
		size_t _size;
		uint64_t _seed;
        const NTL_zz_pE& _ring;
	}; // class UnparametricRandIters
	
} // LinBox

#endif //__LINBOX_ntl_lzz_pe_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
