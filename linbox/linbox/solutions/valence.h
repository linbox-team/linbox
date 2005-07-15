/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Copyright (C)  1999, Linbox project
// Givaro / Athapascan-1
// Valence computation
// Time-stamp: <15 Jul 05 10:53:43 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
// Modified by Z. Wan to fit in linbox
#ifndef __LINBOX_VALENCE_H__
#define __LINBOX_VALENCE_H__

#include <vector>
#include <linbox/solutions/minpoly.h>

namespace LinBox 
{
	
	/*- @memo Valence of a blackbox linear operator A.
         * This is the coefficient of the smallest degree
         * non zero monomial of the minimal polynomial of A.
	 * @doc The resulting value is a Field Element.
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

#include "linbox/field/modular.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

namespace LinBox {
   
    template <class Blackbox, class MyMethod>
    struct IntegerModularValence {       
        const Blackbox &A;
        const MyMethod &M;

        IntegerModularValence(const Blackbox& b, const MyMethod& n) 
                : A(b), M(n) {}
        
        
        template<typename Field>
	typename Field::Element& operator()(typename Field::Element& v, const Field& F) const {
            typedef typename Blackbox::template rebind<Field>::other FBlackbox;
            FBlackbox * Ap;
            MatrixHom::map(Ap, A, F);
            valence( v, *Ap, M);
            delete Ap;
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
        RandomPrime genprime( 26 ); 
        ChineseRemainder< Modular<double> > cra(3UL);
        IntegerModularValence<Blackbox,MyMethod> iteration(A, M);
        cra(V, iteration, genprime);
        commentator.stop ("done", NULL, "Ivalence");
        return V;
    }

    
} //End of LinBox
#endif
