/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/rank.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __RANK_H
#define __RANK_H

//#include "linbox-config.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/algorithms/blackbox-container-symmetrize.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/matrix/blas-matrix.h"


#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"

#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	
	template <class Blackbox>
	unsigned long &rank (unsigned long                   &r,
			     const Blackbox                  &A)
	{  return rank(r, A, Method::Wiedemann()); }

	const int BlasBound = 1 << 26;
	/** Compute the rank of a linear operator A, represented as a black box
	\ingroup solutions
	 */

	/// 
	/** \returns reference to r.
	For now the default method is Wiedemann(), using diagonal preconditioning and 
        the minpoly.  For small or dense matrices BlasElimination will be faster.
	\param r output rank of A.
	\param A linear transform, member of any blackbox class.
	\param M may be a Method::Wiedemann (default) or a Method::BlasElimination.
	\ingroup solutions
	*/
	template <class Blackbox, class Method>
	unsigned long &rank (unsigned long                   &r,
			     const Blackbox                  &A,
			     const Method    &M){}  // should be error here. 

	template <class Blackbox>
	unsigned long &rank (unsigned long                   &r,
			     const Blackbox                  &A,
			     const Method::Wiedemann    &M) 
	{
	    typedef typename Blackbox::Field Field;
		const Field F = A.field();
		typename Field::RandIter iter (F);

		commentator.start ("Rank", "rank");


		std::vector<typename Field::Element> d1, d2;
		size_t i;

		VectorWrapper::ensureDim (d1, A.coldim ());
		VectorWrapper::ensureDim (d2, A.rowdim ());

		for (i = 0; i < A.coldim (); i++)
			do iter.random (d1[i]); while (F.isZero (d1[i]));

		for (i = 0; i < A.rowdim (); i++)
			do iter.random (d2[i]); while (F.isZero (d2[i]));

		Diagonal<Field> D1 (F, d1), D2 (F, d2);
		Transpose<Blackbox> AT (&A);

		Compose<Diagonal<Field>,Transpose<Blackbox> > B1 (&D1, &AT);
		Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> > B2 (&B1, &D2);
		Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox> B3 (&B2, &A);
		Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox>, Diagonal<Field> > B (&B3, &D1);
                    // JGD 22.03.03
// 		BlackboxContainer<Field, Vector> TF (&B, F, iter);
// 		MasseyDomain<Field, BlackboxContainer<Field, Vector> > WD (&TF, M.earlyTermThreshold ());

		typedef Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox>, Diagonal<Field> > Blackbox1;
		BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);
		MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > WD (&TF, M.earlyTermThreshold ());

                    // Here there is an extra diagonal computation
                    // The probability of success is also divided by two, as 
                    // D2^2 contains only squares and squares are half the total elements
// 		BlackboxContainerSymmetrize<Field, Vector> TF (&B2, F, iter);
// 		MasseyDomain<Field, BlackboxContainerSymmetrize<Field, Vector> > WD (&TF, M.earlyTermThreshold ());

		WD.pseudo_rank (r);

		commentator.stop ("done", NULL, "rank");

		return r;
	}
	template <class Blackbox>
	unsigned long &symmetricRank (unsigned long                   &res,
				      const Blackbox                  &A,
				      const Method::Wiedemann    &M) 
	{
	    typedef typename Blackbox::Field Field;
	    const Field F = A.field();
	    typename Field::RandIter iter (F);
	    
	    commentator.start ("Symmetric Rank", "rank");
		

	    std::vector<typename Field::Element> d1, d2;
	    size_t i;

	    VectorWrapper::ensureDim (d1, A.coldim ());

	    for (i = 0; i < A.coldim (); i++)
		    do iter.random (d1[i]); while (F.isZero (d1[i]));


	    Diagonal<Field> D1 (F, d1);

	    
	    Compose<Diagonal<Field>,Blackbox > B1 (&D1, &A);
	    typedef Compose<Compose<Diagonal<Field>,Blackbox >, Diagonal<Field> > BlackBox1;
	    BlackBox1 B (&B1, &D1);

	    BlackboxContainerSymmetric<Field, BlackBox1> TF (&B, F, iter);
	    MasseyDomain<Field, BlackboxContainerSymmetric<Field, BlackBox1> > WD (&TF, M.earlyTermThreshold ());

	    WD.pseudo_rank (res);

	    commentator.stop ("done", NULL, "rank");

	    return res;
	}

	template <class Matrix>
	unsigned long &rank (unsigned long                       &r,
			     const Matrix                          &A,
			     const Method::SparseElimination  &M) 
	{
	    typedef typename Matrix::Field Field;
		const Field F = A.field();
		commentator.start ("Rank", "rank");

		GaussDomain<Field> GD (F);

		Matrix A1 (A);   // We make a copy as these data will be destroyed

		GD.rankin (r, A1, M.strategy ());
                
		commentator.stop ("done", NULL, "rank");
                
		return r;
	}
    
	// using BlasElimination method
	template <class Blackbox>
	unsigned long &rank (unsigned long                     &r,
			     const Blackbox                      &A,
			     const Method::BlasElimination  &M) 
	{
	    typedef typename Blackbox::Field Field;
		const Field F = A.field();
		integer a, b; F.characteristic(a); F.cardinality(b);
		linbox_check( a == b );
		linbox_check( a < LinBox::BlasBound);
		BlasMatrix<typename Field::Element> B(A);
		BlasMatrixDomain<Field> D(F);
		return r = D.rank(B);
	}


	template <class Matrix>
	unsigned long &rankin (unsigned long                      &r,
			       Matrix                               &A,
			       const Method::SparseElimination &M) 
	{
	    typedef typename Matrix::Field Field;
		const Field F = A.field();
		commentator.start ("Rank", "rank");

		GaussDomain<Field> GD (F);

		GD.rankin( r, A, M.strategy ());

		commentator.stop ("done", NULL, "rank");

		return r;
	}

	template <class Field>
	unsigned long &rankin (unsigned long                     &r,
			       BlasBlackbox<Field>                 &A,
			       const Method::BlasElimination  &M) 
	{
		const Field F = A.field();
		BlasMatrixDomain<Field> D(F);
		return r = D.rankin(static_cast< BlasMatrix<typename Field::Element>& >(A));
	}

} // LinBox
#endif // __RANK_H
