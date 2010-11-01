/* Copyright (C) 1999 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr> 
 * Modified by Z. Wan to fit in linbox
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


// Givaro / Athapascan-1
// Valence computation

#ifndef __LINBOX_valence_H
#define __LINBOX_valence_H

#include <vector>
#include <linbox/blackbox/transpose.h>

#include <linbox/solutions/minpoly.h>

namespace LinBox 
{
	
	/*- @brief Valence of a blackbox linear operator A.
         * This is the coefficient of the smallest degree
         * non zero monomial of the minimal polynomial of A.
	 * The resulting value is a Field Element.
	 */
	template < class Blackbox, class DomainCategory, class MyMethod>
	typename Blackbox::Field::Element &valence (typename Blackbox::Field::Element & V,
                                                    const Blackbox& A,
                                                    const DomainCategory& tag,
                                                    const MyMethod& M);


	/** \brief Compute the valence of A
	 *
	 * The valence of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param v Field element into which to store the result
	 * @param A Black box of which to compute the determinant
	 * @param M may is a Method.
         \ingroup solutions
        */
    template <class Blackbox, class MyMethod>
    typename Blackbox::Field::Element &valence (typename Blackbox::Field::Element         &v, 
                                                const Blackbox                              &A,
                                                const MyMethod                           &M) 
    {
        return valence(v, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
    }

	// The valence with default Method 
    template<class Blackbox>
    typename Blackbox::Field::Element &valence (typename Blackbox::Field::Element         &v, 
                                                const Blackbox                               &A)
    {
        return valence(v, A, Method::Hybrid());
    }

    template<class Blackbox, class MyMethod>
    typename Blackbox::Field::Element &valence (
	typename Blackbox::Field::Element         &v, 
        const Blackbox                            &A,
        const RingCategories::ModularTag          &tag,
	const MyMethod& M)
    {
        typedef typename Blackbox::Field::Element Elt_t;
        std::vector<Elt_t> minp;
        minpoly(minp, A, tag, M);
        typename std::vector<Elt_t>::const_iterator it = minp.begin();
        for( ; it != minp.end(); ++it)
            if (! A.field().isZero(*it)) break;
        if (it != minp.end())
            return v=*it;
        else
            return A.field().init(v,0UL);
    }

}

#include "linbox/field/modular-double.h"
#include "linbox/field/givaro-zpz.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-early-single.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include <typeinfo>

namespace LinBox {
   
    template <class Blackbox, class MyMethod>
    struct IntegerModularValence {       
        const Blackbox &A;
        const MyMethod &M;

        IntegerModularValence(const Blackbox& b, const MyMethod& n) 
                : A(b), M(n) {}
        
        
        template<typename Field>
	typename Field::Element& operator()(typename Field::Element& v, const Field& F) const {
        commentator.start ("Modular Valence", "Mvalence");
	std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        F.write(report) << std::endl;
            typedef typename Blackbox::template rebind<Field>::other FBlackbox;
            report << typeid(A).name() << ", A is: " << A.rowdim() << 'x' << A.coldim() << std::endl;
            
            FBlackbox Ap(A, F);

            report << typeid(Ap).name() << ", Ap is: " << Ap.rowdim() << 'x' << Ap.coldim() << std::endl;

            valence( v, Ap, M);
        F.write( F.write(report << "one valence: ", v) << " mod " ) << std::endl;;
        commentator.stop ("done", NULL, "Mvalence");
            return v;
        }            
    };

    template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &valence (typename Blackbox::Field::Element &V, 
                                                    const Blackbox                     &A,
                                                    const RingCategories::IntegerTag   &tag,
                                                    const MyMethod                     &M)
    {
        commentator.start ("Integer Valence", "Ivalence");
#if __LINBOX_SIZEOF_LONG == 8
        RandomPrimeIterator genprime( 31 ); 
        ChineseRemainder< EarlySingleCRA< GivaroZpz<Std64> > > cra(3UL);
#else
        RandomPrimeIterator genprime( 26 ); 
        ChineseRemainder< EarlySingleCRA< Modular<double> > > cra(3UL);
#endif
        IntegerModularValence<Blackbox,MyMethod> iteration(A, M);
        cra(V, iteration, genprime);
        commentator.stop ("done", NULL, "Ivalence");
        return V;
    }

    
} //End of LinBox


namespace LinBox {
	

class Valence {
	public:

	// compute the bound for eigenvalue of AAT via oval of cassini
	// works with both SparseMatrix and DenseMatrix
	template <class Blackbox>
	static integer& cassini (integer& r, const Blackbox& A) {
		//commentator.start ("Cassini bound", "cassini");
    	integer _aat_diag, _aat_radius, _aat_radius1;
    	typedef typename Blackbox::Field Ring;
		_aat_diag = 0; _aat_radius = 0, _aat_radius1 = 0;

        std::vector< integer > d(A. rowdim()),w(A. coldim());
        std::vector<integer>::iterator di, wi;
        for(wi = w.begin();wi!= w.end();++wi) 
            *wi = 0;
        for(di = d.begin();di!= d.end();++di) 
            *di = 0;
		//typename Blackbox::ConstRowIterator row_p;
		typename Blackbox::Element tmp_e;
		Ring R(A. field());
		integer tmp; size_t i, j;
	
		for (j = 0, di = d. begin(); j < A. rowdim(); ++ j, ++ di) {
			// not efficient, but I am not tired of doing case by case
			for ( i = 0; i < A. coldim(); ++ i) {
				R. assign(tmp_e, A.getEntry( j, i));
				R. convert (tmp, tmp_e);
				if (tmp != 0) {
					*di += tmp * tmp;
					w [(int) i] += abs (tmp);
				}

            _aat_diag = _aat_diag >= *di ? _aat_diag : *di;
            }
        }

		for (j = 0, di = d. begin(); j < A. rowdim(); ++ j, ++ di) {
           	integer local_radius = 0;
			for (i = 0; i < A. coldim(); ++ i) {
				R. assign (tmp_e, A. getEntry (j, i));
				R. convert (tmp, tmp_e);
				if (tmp != 0) 
					local_radius += abs (tmp) * w[(int)i];
			}
			local_radius -= *di;
			if ( local_radius > _aat_radius1) {
				if ( local_radius > _aat_radius) {
					_aat_radius1 = _aat_radius;
					_aat_radius = local_radius;
				} else
					_aat_radius1 = local_radius;
			}
		}

        r = _aat_diag + (integer)sqrt( _aat_radius * _aat_radius1 );
	commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	//std::cout << "Cassini bound (AAT) =: " << r << std::endl;
		//commentator. stop ("done", NULL, "cassini");
		return r;
    }   

	// compute one valence of AAT over a field
	template <class Blackbox>
	static void one_valence(typename Blackbox::Element& v, unsigned long& r, const Blackbox& A) {
		//commentator.start ("One valence", "one valence");
		typedef std::vector<typename Blackbox::Element> Poly; Poly poly;
		typename Blackbox::Field F(A. field());
		Transpose<Blackbox> AT (&A);
		Compose<Blackbox, Transpose<Blackbox> > AAT(&A, &AT);
		// compute the minpoly of AAT
		minpoly(poly, AAT, Method::Wiedemann());
		typename Poly::iterator p;
		F. init (v, 0);

		for (p = poly. begin(); p != poly. end(); ++ p)
			if (! F. isZero (*p)) {
				F. assign (v, *p);
				break;
			}
		
		r = poly. size() -1;
		std::ostream& report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
			
		//	std::ostream& report = std::cout;
		report << "one valence =: " << v << " over ";
		A. field(). write(report); report << std::endl;
		//commentator. stop ("done", NULL, "one valence");
		return;
	}

	// compute the valence of AAT over an integer ring
	template <class Blackbox>
	static void valence(Integer& val, const Blackbox& A) {
		commentator. start ("Valence (AAT)", "Valence");
		typedef Modular<int32> Field;
		typedef typename MatrixHomTrait<Blackbox, Field>::value_type FBlackbox;
		int n_bit = (int)(log((double)Field::getMaxModulus()) / M_LN2 - 2);
		unsigned long d; 
                RandomPrimeIterator g(n_bit); Field::Element v;
		++g; Field F(*g);
		FBlackbox Ap(A, F);
		one_valence(v, d, Ap);
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		//std::cout<<"degree of minpoly of AAT: " << d << std::endl;
		valence (val, d, A);
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Integer valence =: " << val << std::endl;
		commentator. stop ("done", NULL, "Valence");
		return;
	}

	// compute the valence of AAT over an integer ring
	// d, the degree of min_poly of AAT
	template <class Blackbox>
	static void valence(Integer& val, unsigned long d, const Blackbox& A) {

		typedef Modular<int32> Field;
		typedef typename MatrixHomTrait<Blackbox, Field>::value_type FBlackbox;
		int n_bit = (int)(log((double)Field::getMaxModulus()) / M_LN2 - 2);
                RandomPrimeIterator rg(n_bit);
		std::vector<integer> Lv, Lm;
		unsigned long d1; Field::Element v; integer im = 1;
		//compute an upper bound for val.
		integer bound; cassini (bound, A); bound = pow (bound, d); bound *= 2;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Bound for valence: " << bound << std::endl;

		do {
			++rg;
			Field F(*rg);
                        FBlackbox Ap(A, F);
			one_valence(v, d1, Ap);
			if (d1 == d) {
				im *= *rg;
				Lm. push_back ( *rg ); Lv. push_back (integer(v));
			}
		} while (im < bound);

		val = 0;
		std::vector<integer>::iterator Lv_p, Lm_p; integer tmp, a, b, g;
		for (Lv_p = Lv. begin(), Lm_p = Lm. begin(); Lv_p != Lv. end(); ++ Lv_p, ++ Lm_p) {
			tmp = im / *Lm_p;
			gcd (g, *Lm_p, tmp, a, b);
			val += *Lv_p * b * tmp;
			val %= im;
		}

		if (sign (val) < 0)
			val += im;
		tmp = val - im;
		if (abs(tmp) < abs(val)) 
			val = tmp;

		return;
	}
};
} //End of LinBox
#endif //__LINBOX_valence_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
