/* linbox/solutions/det.h
 * Copyright (C) 2001, 2002 LinBox
 * Time-stamp: <11 Sep 25 17:18:22 Jean-Guillaume.Dumas@imag.fr>
 *
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

/*! @file solutions/det.h
 * @ingroup solutions
 * @brief NO DOC
 */

#ifndef __LINBOX_det_H
#define __LINBOX_det_H

#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/getentry.h"
#include "linbox/vector/blas-vector.h"

#include "linbox/matrix/dense-matrix.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/prime-stream.h"
#include "linbox/util/debug.h"
#include "linbox/util/mpicpp.h"


// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** \brief Compute the determinant of A.
	 *
	 * The determinant of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param d    Field element into which to store the result
	 * @param A    Black box of which to compute the determinant
	 * @param tag  optional tag.  Specifies Integer, Rational or modular ring/field
	 * @param M    optional method.  The default is Method::Auto(), Other options
	 include Blackbox, Elimination, Wiedemann, DenseElimination and SparseElimination.
	 Sometimes it helps to	 indicate properties of the matrix in the method object
	 (for instance symmetry). See class Method for details.
	 \ingroup solutions
	 */
	template< class Blackbox, class DetMethod, class DomainCategory>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element          &d,
						const Blackbox                             &A,
						const DomainCategory                     &tag,
						const DetMethod                           &Meth);

	// The det where A can be modified in place
	// Default is to use the generic det (might copy)
	template< class Blackbox, class DetMethod, class DomainCategory>
	typename Blackbox::Field::Element &detInPlace (typename Blackbox::Field::Element	&d,
						  Blackbox	&A,
						  const DomainCategory			&tag,
						  const DetMethod			&Meth)
	{
		return det(d, A, tag, Meth);
	}

	// The det with default Method
	template<class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element	&d,
						const Blackbox				&A)
	{
		return det(d, A, Method::Auto());
	}

	// The det where A can be modified in place
	template<class Blackbox>
	typename Blackbox::Field::Element &detInPlace (typename Blackbox::Field::Element	&d,
						  Blackbox				&A)
	{
		return detInPlace(d, A, Method::Auto());
	}

	// The det with category specializer
	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element	&d,
						const Blackbox				&A,
						const MyMethod				&Meth)
	{
		return det(d, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), Meth);
	}

	// The in place det with category specializer
	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &detInPlace (typename Blackbox::Field::Element     &d,
						  Blackbox                              &A,
						  const MyMethod                        &Meth)
	{
		return detInPlace(d, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), Meth);
	}

	// The det with Auto Method
	template<class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element	&d,
						const Blackbox				&A,
						const RingCategories::ModularTag	&tag,
						const Method::Auto			&Meth)
	{
		if (useBlackboxMethod(A))
			return det(d, A, tag, Method::Blackbox(Meth));
		else

			return det(d, A, tag, Method::Elimination(Meth));
	}
	template<class Blackbox>
	typename Blackbox::Field::Element &detInPlace (typename Blackbox::Field::Element	&d,
						  Blackbox				&A,
						  const RingCategories::ModularTag	&tag,
						  const Method::Auto			&Meth)
	{
		/*
		   if (useBlackboxMethod(A))
		   return det(d, A, tag, Method::Blackbox(Meth));
		   else
		   */
		return detInPlace(d, A, tag, Method::Elimination(Meth));
	}

	// The det with Auto Method on BlasMatrix
	template<class Field>
	typename Field::Element &det (typename Field::Element	&d,
				      const BlasMatrix<Field>		&A,
				      const RingCategories::ModularTag	&tag,
				      const Method::Auto		&Meth)
	{
		return det(d, A, tag, Method::Elimination(Meth));
	}

	template<class Field>
	typename Field::Element &detInPlace (typename Field::Element	&d,
					BlasMatrix<Field>			&A,
					const RingCategories::ModularTag	&tag,
					const Method::Auto			&Meth)
	{
		return detInPlace(d, A, tag, Method::Elimination(Meth));
	}

	// Forward declaration saves us from including blackbox/toeplitz.h
	template<class A, class B> class Toeplitz;

	// Toeplitz determinant
	template<class CField, class PField >
	typename CField::Element& det(typename CField::Element		& res,
				      const Toeplitz<CField,PField>	& A )
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		return A.det(res);
	}
	template<class CField, class PField >
	typename CField::Element& detInPlace(typename CField::Element	& res,
					Toeplitz<CField,PField>		& A )
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		return A.det(res);
	}
	// The det with BlackBox Method
	template<class Blackbox>
	typename Blackbox::Field::Element &det (
						typename Blackbox::Field::Element       &d,
						const Blackbox                          &A,
						const RingCategories::ModularTag        &tag,
						const Method::Blackbox			&Meth)
	{
		return det(d, A, tag, Method::Wiedemann(Meth));
	}


	// The det with Wiedemann, finite field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element	&d,
						const Blackbox				&A,
						const RingCategories::ModularTag	&tag,
						const Method::Wiedemann			&Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		typedef typename Blackbox::Field Field;
		Field F = A.field();
		typedef BlasVector<Field> Polynomial;

		if(Meth.shapeFlags == Shape::Symmetric) {
            commentator().start ("Symmetric Wiedemann Determinant", "sdet");
			linbox_check (A.coldim () == A.rowdim ());
			Polynomial               phi(F);
			size_t            deg;
			typename Field::RandIter iter (F);

			// Precondition here to separate the eigenvalues, so that
			// minpoly (B) = charpoly (B) with high probability
			// Here there is an extra diagonal computation
			// The probability of success is also divided by two, as
			// diag^2 contains only squares and squares are half the total elements
			BlasVector<Field> diag (A.field(),A.coldim ());

			typename Field::Element pi;
			size_t i;
#if 0
			size_t iternum = 1;
#endif
			do {
				F.assign(pi, F.one);
				for (i = 0; i < A.coldim (); i++) {
					do iter.random (diag[i]); while (F.isZero (diag[i]));
					F.mulin (pi, diag[i]);
				}

				Diagonal<Field> D (diag);
				Compose<Blackbox,Diagonal<Field> > B_0 (&A, &D);
				typedef Compose<Diagonal<Field>,Compose<Blackbox,Diagonal<Field> > > Blackbox1;
				Blackbox1 B(&D, &B_0);

				BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);

				MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > WD (&TF, Meth.earlyTerminationThreshold);

				WD.minpoly (phi, deg);
#if 0
				std::cout << "\tdet: iteration # " << iternum << "\tMinpoly deg= "
				<< phi.size() << "\n" ;
				std::cout << "[" ;
				for(typename Polynomial::const_iterator refs =  phi.begin();
				    refs != phi.end() ;
				    ++refs )
					std::cout << (*refs) << " " ;
				std::cout << "]" << std::endl;
				++iternum;
#endif
			} while ( (phi.size () < A.coldim () + 1) && ( !F.isZero (phi[0]) ) );


				// Divided twice since multiplied twice by the diagonal matrix
            F.div (d, phi[0], pi);
            F.divin (d, pi);

            if ( (deg & 1) == 1)
                F.negin (d);

            commentator().stop ("done", NULL, "sdet");

            return d;
		}
		else {
            commentator().start ("Wiedemann Determinant", "wdet");
			linbox_check (A.coldim () == A.rowdim ());

			Polynomial               phi(F);
			size_t            deg;
			typename Field::RandIter iter (F);

			// Precondition here to separate the eigenvalues, so that
			// minpoly (B) = charpoly (B) with high probability
			BlasVector<Field> diag (F,A.coldim ());

			typename Field::Element pi;
			size_t i;
			do {
				F.assign(pi, F.one);
				for (i = 0; i < A.coldim (); i++) {
					do iter.random (diag[i]); while (F.isZero (diag[i]));
					F.mulin (pi, diag[i]);
				}

				Diagonal<Field> D (diag);

				Compose<Blackbox,Diagonal<Field> > B (&A, &D);

				typedef Compose<Blackbox,Diagonal<Field> > Blackbox1;

				BlackboxContainer<Field, Blackbox1> TF (&B, F, iter);

				MasseyDomain<Field, BlackboxContainer<Field, Blackbox1> > WD (&TF, Meth.earlyTerminationThreshold);

				WD.minpoly (phi, deg);
			} while ( (phi.size () < A.coldim () + 1) && ( !F.isZero (phi[0]) ) );

            F.div (d, phi[0], pi);

            if ( (deg & 1) == 1)
                F.negin (d);


            commentator().stop ("done", NULL, "wdet");

            return d;
		}
	}



	// the det with Blas, finite field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element       &d,
						const Blackbox                          &A,
						const RingCategories::ModularTag        &tag,
						const Method::DenseElimination           &Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		typedef typename Blackbox::Field Field;
		Field F = A.field();

		commentator().start ("Blas Determinant", "blasdet");

		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrix<Field> B(A);
		BlasMatrixDomain<Field> BMD(F);
		d= BMD.detInPlace(B);
		commentator().stop ("done", NULL, "blasdet");

		return d;
	}

	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element	&d,
						const Blackbox			&A,
						const RingCategories::ModularTag	&tag,
						const Method::SparseElimination		&Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		typedef typename Blackbox::Field Field;
		commentator().start ("Sparse Elimination Determinant", "SEDet");
		// We make a copy as these data will be destroyed
		SparseMatrix<Field, SparseMatrixFormat::SparseSeq> A1 (A.field(), A.rowdim(), A.coldim());
		typename Blackbox::Field::Element tmp;
		for(size_t i = 0; i < A.rowdim() ; ++i)
			for(size_t j = 0; j < A.coldim(); ++j)
				A1.setEntry(i,j,getEntry(tmp, A, i, j));
		GaussDomain<Field> GD ( A1.field() );
		GD.detInPlace (d, A1, Meth.pivotStrategy);
		commentator().stop ("done", NULL, "SEDet");
		return d;

	}


	template <class Field, class Vector>
	typename Field::Element &det (typename Field::Element	&d,
				      const SparseMatrix<Field, Vector>	&A,
				      const RingCategories::ModularTag	&tag,
				      const Method::SparseElimination		&Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		commentator().start ("Sparse Elimination Determinant", "SEDet");
		// We make a copy as these data will be destroyed
		SparseMatrix<Field, SparseMatrixFormat::SparseSeq> A1 (A);
		GaussDomain<Field> GD ( A.field() );
		GD.detInPlace (d, A1, Meth.pivotStrategy);
		commentator().stop ("done", NULL, "SEDet");
		return d;
	}

	template <class Field>
	typename Field::Element &detInPlace (typename Field::Element	&d,
					SparseMatrix<Field, SparseMatrixFormat::SparseSeq>  &A,
					const RingCategories::ModularTag	&tag,
					const Method::SparseElimination	&Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
        commentator().start ("Sparse Elimination Determinant in place", "SEDetin", A.rowdim() );
		GaussDomain<Field> GD ( A.field() );
		GD.detInPlace (d, A, Meth.pivotStrategy);
        commentator().stop ("done", NULL, "SEDetin");
		return d;
	}


	// The det with Elimination Method
	template<class Field, class Vector>
	typename Field::Element &det (typename Field::Element		&d,
				      const SparseMatrix<Field, Vector>	&A,
				      const RingCategories::ModularTag	&tag,
				      const Method::Elimination		&Meth)
	{
		return det(d, A, tag, Method::SparseElimination(Meth));
	}


	template <class Field>
	typename Field::Element &detInPlace (typename Field::Element	&d,
					SparseMatrix<Field, SparseMatrixFormat::SparseSeq>  &A,
					const RingCategories::ModularTag	&tag,
					const Method::Elimination		&Meth)
	{
		return detInPlace(d, A, tag, Method::SparseElimination(Meth));
	}

	template<class Field, class Vector>
	typename Field::Element &detInPlace (typename Field::Element			&d,
					SparseMatrix<Field, Vector>		&A,
					const RingCategories::ModularTag	&tag,
					const Method::Elimination		&Meth)
	{
		// Matrix is not of type SparseMatrix<..SparseSeq> otherwise previous specialization would occur
		// will copy A into SparseMatrix<..SparseSeq> or BlasMatrix
		const Field& F = A.field();
		integer c; F.characteristic(c);
		if ((c < LinBox::BlasBound) && ((A.rowdim() < 300) || (A.coldim() < 300) || (A.size() > (A.coldim()*A.rowdim()/100))))
			return det(d, A, tag, Method::DenseElimination(Meth));
		else
			return det(d, A, tag, Method::SparseElimination(Meth));
	}


	// The det with Elimination Method
	template<class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element	&d,
						const Blackbox				&A,
						const RingCategories::ModularTag	&tag,
						const Method::Elimination		&Meth)
	{
		// Matrix is not of type SparseMatrix otherwise previous specialization would occur
		// will copy A into BlasMatrix
		return det(d, A, tag, Method::DenseElimination(Meth));
	}


	template<class Blackbox>
	typename Blackbox::Field::Element &detInPlace (typename Blackbox::Field::Element	&d,
						  Blackbox	&A,
						  const RingCategories::ModularTag      &tag,
						  const Method::Elimination		&Meth)
	{
		// Matrix is not of type SparseMatrix not of type BlasMatrix
                // otherwise previous specialization would occur
		// will copy A into BlasMatrix
		return det(d, A, tag, Method::DenseElimination(Meth));
	}

	template<class Field>
	typename Field::Element &detInPlace (typename Field::Element			&d,
                                        BlasMatrix<Field>			&A,
                                        const RingCategories::ModularTag	&tag,
					const Method::Elimination		&Meth)
	{
		return detInPlace(d, A);
	}

	template<class Field>
	typename Field::Element &detInPlace (typename Field::Element			&d,
                                        BlasMatrix<Field>			&A,
                                        const RingCategories::ModularTag	&tag,
					const Method::DenseElimination		&Meth)
	{
		return detInPlace(d, A);
	}



	// This should work for a BlasMatrix too ?
	/** Rank of Blackbox \p A.
	  * \ingroup solutions
	  * A will be modified.
	  * \param[out]  d determinant of \p A.
	  * \param       A this BlasMatrix matrix will be modified in place in the process.
	  * \return \p d
	  */
	template <class Field>
	typename Field::Element &detInPlace (typename Field::Element             &d,
					BlasMatrix<Field>                  &A)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		Field F = A.field();

		commentator().start ("Determinant", "detInPlace");
		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrixDomain<Field> BMD(F);
		d= BMD.detInPlace(static_cast<BlasMatrix<Field>& > (A));
		commentator().stop ("done", NULL, "detInPlace");

		return d;
	}
} // end of LinBox namespace

#include "linbox/ring/modular.h"
//#include "linbox/field/givaro-zpz.h"

#ifdef __LINBOX_HAVE_MPI
#include "linbox/algorithms/cra-distributed.h"
#else
#ifdef __LINBOX_HAVE_KAAPI //use the kaapi version instead of the usual version if this macro is defined
#include "linbox/algorithms/cra-kaapi.h"
#else
#include "linbox/algorithms/cra-domain.h"
#endif
#endif

#include "linbox/algorithms/cra-builder-single.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

namespace LinBox
{

	template <class Blackbox, class MyMethod>
	struct IntegerModularDet {
		const Blackbox &A;
		const MyMethod &M;

		IntegerModularDet(const Blackbox& b, const MyMethod& n) :
			A(b), M(n)
		{}


		template<class Element, typename Field>
		IterationResult operator()(Element& d, const Field& F) const
		{
			typedef typename Blackbox::template rebind<Field>::other FBlackbox;
			FBlackbox Ap(A, F);
			detInPlace( d, Ap, RingCategories::ModularTag(), M);
			return IterationResult::CONTINUE;
		}
	};


	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &cra_det (typename Blackbox::Field::Element         &d,
						    const Blackbox                            &A,
						    const RingCategories::IntegerTag          &tag,
						    const MyMethod                            &Meth
#ifdef __LINBOX_HAVE_MPI
						    ,Communicator                             *C = NULL
#endif
						   )
	{
		//  if no parallelism or if this is the parent process
		//  begin the verbose output
#ifdef __LINBOX_HAVE_MPI
		if(!C || C->rank() == 0)
#endif
            commentator().start ("Integer Determinant", "idet");
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		IntegerModularDet<Blackbox, MyMethod> iteration(A, Meth);
                typedef Givaro::ModularBalanced<double> Field;
                PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<Field>::bestBitSize(A.coldim()));
		integer dd; // use of integer due to non genericity of cra. PG 2005-08-04

		//  will call regular cra if C=0
#ifdef __LINBOX_HAVE_MPI
		ChineseRemainderDistributed< CRABuilderEarlySingle< Field > > cra(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, C);
		cra(dd, iteration, genprime);
		if(!C || C->rank() == 0){
			A.field().init(d, dd); // convert the result from integer to original type
            commentator().stop ("done", NULL, "idet");
		}
#else
		ChineseRemainder< CRABuilderEarlySingle< Field > > cra(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);
		cra(dd, iteration, genprime);
		A.field().init(d, dd); // convert the result from integer to original type
        commentator().stop ("done", NULL, "idet");
#endif

		return d;
	}

} // end of LinBox namespace

//#if 0
#ifdef __LINBOX_HAVE_NTL
# include "linbox/algorithms/hybrid-det.h"
# define SOLUTION_CRA_DET lif_cra_det
#else
# define SOLUTION_CRA_DET cra_det
#endif

#include "linbox/algorithms/rational-cra-var-prec.h"
#include "linbox/algorithms/cra-builder-var-prec-early-single.h"
#include "linbox/algorithms/det-rational.h"
namespace LinBox
{

	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                            &A,
						const RingCategories::IntegerTag          &tag,
						const MyMethod                            &Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		return SOLUTION_CRA_DET(d, A, tag, Meth);
	}

	template< class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                            &A,
						const RingCategories::RationalTag       &tag,
						const MyMethod                          &Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		commentator().start ("Rational Determinant", "rdet");

		Integer num,den;

		IntegerModularDet<Blackbox, MyMethod> iteration(A, Meth);
                typedef Givaro::ModularBalanced<double> Field;
		PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<Field>::bestBitSize(A.coldim()));
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field > > rra(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);

		rra(num,den, iteration, genprime);

		A.field().init(d, num,den); // convert the result from integer to original type

		commentator().stop ("done", NULL, "rdet");
		return d;
	}

	template<class Field, class MyMethod>
	typename Field::Element &det (typename Field::Element                 &d,
				      const BlasMatrix<Field>                &A,
				      const RingCategories::RationalTag       &tag,
				      const MyMethod                          &Meth)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		commentator().start ("Dense Rational Determinant", "rdet");

		rational_det(d,A,Meth);

		commentator().stop ("done", NULL, "rdet");
		return d;
	}

} // end of LinBox namespace

#ifdef __LINBOX_HAVE_MPI
namespace LinBox
{

	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element       &d,
						const Blackbox                          &A,
						/*const*/ Communicator			&C)
	{
		return det(d, A, Method::Auto(C));
	}
}
#endif //__LINBOX_HAVE_MPI

#endif // __LINBOX_det_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
