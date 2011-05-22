/* linbox/solutions/rank.h
 * Copyright(C) LinBox
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __LINBOX_rank_H
#define __LINBOX_rank_H

//#include "linbox-config.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/algorithms/blackbox-container-symmetrize.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/whisart_trace.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/switch/cekstv.h"
#include "linbox/blackbox/butterfly.h"


#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/trace.h"
#include "linbox/solutions/methods.h"


#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox 
{


	/** 
            Compute the rank of a linear transform A over a field by selected method. 
            For very large and/or very sparse matrices the Wiedemann method will be faster
            (and is memory efficient).
            For some sparse matrices SparseElimination may outperform Wiedemann.
            For small or dense matrices BlasElimination will be faster.
	    \param r - output rank of A.
            \param - input A linear transform, member of any blackbox class.
            \param - input M may be a Method::Wiedemann (the default), a Method::BlasElimination, or a Method::SparseElimination..
            \returns  reference to r.
            \ingroup solutions
	*/
	
 	template <class Blackbox, class Method, class DomainCategory>
 	inline unsigned long &rank (unsigned long                   &r,
 			     const Blackbox                  &A,
			     const DomainCategory          &tag,
 			     const Method                   &M);

	// error hanlder for rational domain
	template <class Blackbox, class Method>
 	inline unsigned long &rank (unsigned long                           &r,
 			     const Blackbox                          &A,
			     const RingCategories::RationalTag     &tag,
 			     const Method                           &M){
            commentator.start ("Rational Rank", "Rrank");
                // Same mapping as the integer one
            rank(r, A, RingCategories::IntegerTag(), M);
            commentator.stop ("done", NULL, "Rrank");
            return r;
	}


	/** 
            Compute the rank of a linear transform A over a field. 

            The default method is Wiedemann(), using diagonal preconditioning and 
            the minpoly.  For small or dense matrices BlasElimination will be faster.
            \returns <em>r</em> rank of A.
            \param A linear transform, member of any blackbox class.
            \ingroup solutions
	*/
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A)
    {  return rank(r, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), Method::Hybrid()); 
    }

	/** A may be modified
	*/
    template <class Matrix>
    inline unsigned long &rankin (unsigned long                   &r,
                         Matrix                  &A)
    {  return rankin(r, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), Method::Elimination()); 
    }

    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const RingCategories::ModularTag     &tag,
                         const Method::Hybrid         &m)
    { // this should become a BB/Blas hybrid in the style of Duran/Saunders/Wan.  
        if (useBB(A)) return rank(r, A, tag, Method::Blackbox(m )); 
        else return rank(r, A, tag, Method::Elimination( m ));
    }

    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const RingCategories::ModularTag                  &tag,
                         const Method::Elimination    &m)
    {  
        typedef typename Blackbox::Field Field;
        const Field& F = A.field();
        integer a, b; F.characteristic(a); F.cardinality(b);
        if (a == b && a < LinBox::BlasBound)
		return rank(r, A, tag, Method::BlasElimination(m));
	else
		return rank(r, A, tag, Method::NonBlasElimination( m ));
    }


    template <class Field, class Vector>
    inline unsigned long &rank (unsigned long                   &r,
                         const SparseMatrix<Field, Vector>         &A,
                         const RingCategories::ModularTag                  &tag,
                         const Method::Elimination    &m)
    {  
        return rank(r, A, tag, Method::SparseElimination(m));
    }


	// specialization of NonBlas for SparseMatrix
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const   RingCategories::ModularTag           &tag,
                         const Method::NonBlasElimination& m)
    {
        return rank(r, A, tag, Method::SparseElimination(m));
    }


    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const  RingCategories::ModularTag                   &tag,
                         const Method::Blackbox& m);
    

	/** 
            Compute the rank of a linear transform A over a field. 

            The default method is Wiedemann(), using diagonal preconditioning and 
            the minpoly.  For small or dense matrices BlasElimination will be faster.
            \returns <em>r</em> rank of A.
            \param A linear transform, member of any blackbox class.
            \ingroup solutions
	*/
    template <class Blackbox, class Method>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const Method    &M)
    {  return rank(r, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M); 
    }

	/// M may be <code>Method::Wiedemann()</code>.
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &res,
                         const Blackbox                  &A,
                         const RingCategories::ModularTag          &tag,
                         const Method::Wiedemann    &M) 
// This is too much for solutions.  It belongs in algorithms
    {
            
        typedef typename Blackbox::Field Field;
        const Field F = A.field();
        typename Field::RandIter iter (F);
            
        if (M.symmetric()) {
            commentator.start ("Symmetric Rank", "srank");
		
                
            std::vector<typename Field::Element> d1;
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
            std::vector<typename Field::Element> phi;
            WD.pseudo_minpoly (phi, res);
            commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Pseudo Minpoly degree: " << res << std::endl;
               
            commentator.start ("Monte Carlo certification", "trace");
            typename Field::Element t, p2; F.init(p2, 0UL);
            trace(t, B);            
            if (phi.size() >= 2) F.neg(p2, phi[ phi.size()-2]);
                
            int nbperm = 0; unsigned long rk;
            int logn = 2*(unsigned long)floor( log( (double)A.rowdim() ) );
            bool tryagain = (! F.areEqual( t, p2 ));
            while( tryagain ) {
                commentator.stop ("fail", NULL, "trace");

//                 Permutation<Field> P(A.rowdim(), F);          
//                 for (i = 0; i < A.rowdim (); ++i)
//                     P.permute( rand() % A.rowdim() , rand() % A.rowdim() );
//                 for (i = 0; i < A.rowdim (); ++i)
//                     P.permute( rand() % A.rowdim() , rand() % A.rowdim() );

//                 Transpose< Permutation<Field> > TP(&P);
//                 typedef Compose< Permutation<Field>, Blackbox > BlackboxP;
//                 typedef Compose< Compose< Permutation<Field>, Blackbox >, Transpose< Permutation<Field> > > BlackboxPAP;
//                 BlackboxP PA(&P, &A);
//                 BlackboxPAP BP( &PA , &TP );

//                 for (i = 0; i < A.coldim (); i++)
//                     do iter.random (d1[i]); while (F.isZero (d1[i]));
//                 Diagonal<Field> D1 (F, d1);
//                 Compose<Diagonal<Field>,BlackboxPAP > B1 (&D1, &BP);
//                 typedef Compose<Compose<Diagonal<Field>,BlackboxPAP >, Diagonal<Field> > BlackBox2;
//                 BlackBox2 B (&B1, &D1);

                for (i = 0; i < A.coldim (); i++)
                    do iter.random (d1[i]); while (F.isZero (d1[i]));
                Diagonal<Field> D1 (F, d1);
                Compose<Diagonal<Field>,Blackbox > B1 (&D1, &A);
                BlackBox1 B (&B1, &D1);
                
                BlackboxContainerSymmetric<Field, BlackBox1> TF (&B, F, iter);
                MasseyDomain<Field, BlackboxContainerSymmetric<Field, BlackBox1> > WD (&TF, M.earlyTermThreshold ());
                
                WD.pseudo_minpoly (phi, rk);
                commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Permuted pseudo Minpoly degree: " << res << std::endl;
                commentator.start ("Monte Carlo certification", "trace");
                if (phi.size() >= 2) F.neg(p2, phi[ phi.size()-2]);
                
                trace(t, B);

                tryagain = (! F.areEqual( t, p2 ));
                if (res > rk) 
                    tryagain = true;
                else
                    res = rk;
                if( ++nbperm > logn) break;
            }
            commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "symm permutations : " << nbperm << std::endl;
            nbperm = 0;
            while(tryagain) {
                commentator.stop ("fail", NULL, "trace");
//             F.write( F.write( commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION)
//                               << "end trace: ", t) << ", p2: ", p2) << std::endl;
                typename Field::RandIter r (F);
                typename CekstvSwitch<Field>::Factory factory (r);
                typedef Butterfly<Field, CekstvSwitch<Field> > ButterflyP;
                ButterflyP P (F, A.rowdim(), factory);
                for (i = 0; i < A.coldim (); i++)
                    do iter.random (d1[i]); while (F.isZero (d1[i]));
                Diagonal<Field> D1 (F, d1);
                typedef Compose< ButterflyP, Diagonal<Field> > ButD;
                ButD PD(&P, &D1);

                Transpose< ButD > TP (&PD);
                
                Compose< ButD, Blackbox > B1( &PD, &A);

                typedef Compose< Compose< ButD, Blackbox > , Transpose< ButD > > BlackBoxBAB;
                BlackBoxBAB PAP(&B1, &TP);
                
                BlackboxContainerSymmetric<Field, BlackBoxBAB> TF (&PAP, F, iter);
                MasseyDomain<Field, BlackboxContainerSymmetric<Field, BlackBoxBAB> > WD (&TF, M.earlyTermThreshold ());
                
                WD.pseudo_minpoly (phi, rk);
                commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Butterfly pseudo Minpoly degree: " << res << std::endl;
                commentator.start ("Monte Carlo certification", "trace");
                if (phi.size() >= 2) F.neg(p2, phi[ phi.size()-2]);
                
                trace(t, PAP);

                tryagain = (! F.areEqual( t, p2 ));
                if (res > rk) 
                    tryagain = true;
                else
                    res = rk;
                ++nbperm;
            }
            
                // F.write( F.write( commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION)
                //              << "end trace: ", t) << ", p2: ", p2) << std::endl;
            
            commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "butterflies : " << nbperm << std::endl;

            commentator.stop ("success", NULL, "trace");
            commentator.stop ("done", NULL, "srank");
            
            return res;
        } else {
                
            commentator.start ("Rank", "wrank");

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
                // Here there is an extra diagonal computation
                // The probability of success is also divided by two, as 
                // D2^2 contains only squares and squares are half the total elements                
            typedef Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox>, Diagonal<Field> > Blackbox1;
            Blackbox1 B (&B3, &D1);

            BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);
            MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > WD (&TF, M.earlyTermThreshold ());
                
            std::vector<typename Field::Element> phi;
            WD.pseudo_minpoly (phi, res);
            commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Pseudo Minpoly degree: " << res << std::endl;
            commentator.start ("Monte Carlo certification", "trace");

            typename Field::Element t, p2; F.init(p2, 0UL);
//             trace(t, B);
            WhisartTraceTranspose(t, F, D1, A, D2);
            if (phi.size() >= 2) F.neg(p2, phi[ phi.size()-2]);
                
            int nbperm = 0; unsigned long rk;
            int logn = 2*(unsigned long)floor( log( (double)A.rowdim() ) );
            bool tryagain = (! F.areEqual( t, p2 ));
            while( tryagain ) {
                commentator.stop ("fail", NULL, "trace");
                Permutation<Field> P(A.rowdim(), F);          
                for (i = 0; i < A.rowdim (); ++i)
                    P.permute( rand() % A.rowdim() , rand() % A.rowdim() );
                for (i = 0; i < A.rowdim (); ++i)
                    P.permute( rand() % A.rowdim() , rand() % A.rowdim() );

                typedef Compose< Permutation<Field>, Blackbox > BlackboxP;
                BlackboxP BP(&P, &A);

                for (i = 0; i < A.coldim (); i++)
                    do iter.random (d1[i]); while (F.isZero (d1[i]));
                
                for (i = 0; i < A.rowdim (); i++)
                    do iter.random (d2[i]); while (F.isZero (d2[i]));
                
                Diagonal<Field> D1 (F, d1), D2 (F, d2);
                Transpose<BlackboxP> AT (&BP);
                
                Compose<Diagonal<Field>,Transpose<BlackboxP> > B1 (&D1, &AT);
                Compose<Compose<Diagonal<Field>,Transpose<BlackboxP> >, Diagonal<Field> > B2 (&B1, &D2);
                Compose<Compose<Compose<Diagonal<Field>,Transpose<BlackboxP> >, Diagonal<Field> >, BlackboxP> B3 (&B2, &BP);
                typedef Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<BlackboxP> >, Diagonal<Field> >, BlackboxP>, Diagonal<Field> > Blackbox1;
                Blackbox1 B (&B3, &D1);

                BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);
                MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > MD (&TF, M.earlyTermThreshold ());
                
                MD.pseudo_minpoly (phi, rk);
                commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Permuted pseudo Minpoly degree: " << rk << std::endl;
                commentator.start ("Monte Carlo certification", "trace");
                if (phi.size() >= 2) F.neg(p2, phi[ phi.size()-2]);
                
//                 trace(t, B);
                WhisartTraceTranspose(t, F, D1, BP, D2);
                tryagain = (! F.areEqual( t, p2 ));
                if (res > rk) 
                    tryagain = true;
                else
                    res = rk;
                if( ++nbperm > logn) break;
            }
            commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "permutations : " << nbperm << std::endl;
            nbperm = 0;
            while(tryagain) {
                commentator.stop ("fail", NULL, "trace");
                typename Field::RandIter r (F);
                typename CekstvSwitch<Field>::Factory factory (r);
                typedef Butterfly<Field, CekstvSwitch<Field> > ButterflyP;
                ButterflyP P (F, A.rowdim(), factory);
                
                typedef Compose< ButterflyP, Blackbox > BlackboxP;
                BlackboxP BP(&P, &A);

                for (i = 0; i < A.coldim (); i++)
                    do iter.random (d1[i]); while (F.isZero (d1[i]));
                
                for (i = 0; i < A.rowdim (); i++)
                    do iter.random (d2[i]); while (F.isZero (d2[i]));
                
                Diagonal<Field> D1 (F, d1), D2 (F, d2);
                Transpose<BlackboxP> AT (&BP);
                
                Compose<Diagonal<Field>,Transpose<BlackboxP> > B1 (&D1, &AT);
                Compose<Compose<Diagonal<Field>,Transpose<BlackboxP> >, Diagonal<Field> > B2 (&B1, &D2);
                Compose<Compose<Compose<Diagonal<Field>,Transpose<BlackboxP> >, Diagonal<Field> >, BlackboxP> B3 (&B2, &BP);
                typedef Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<BlackboxP> >, Diagonal<Field> >, BlackboxP>, Diagonal<Field> > Blackbox1;
                Blackbox1 B (&B3, &D1);

                BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);
                MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > MD (&TF, M.earlyTermThreshold ());
                
                MD.pseudo_minpoly (phi, rk);
                commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Butterfly pseudo Minpoly degree: " << rk << std::endl;
                commentator.start ("Monte Carlo certification", "trace");
                if (phi.size() >= 2) F.neg(p2, phi[ phi.size()-2]);
                
//                 trace(t, B);
                WhisartTraceTranspose(t, F, D1, BP, D2);
                tryagain = (! F.areEqual( t, p2 ));
                if (res > rk) 
                    tryagain = true;
                else
                    res = rk;
                ++nbperm;
            }
            commentator.report(Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "butterflies : " << nbperm << std::endl;
            commentator.stop ("success", NULL, "trace");
            commentator.stop ("done", NULL, "wrank");
            
            return res;
        }

    }

	/// M may be <code>Method::SparseElimination()</code>.
    template <class Field>
    inline unsigned long &rank (unsigned long                       &r,
                         const SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq>  &A,
                         const RingCategories::ModularTag    &tag,
                         const Method::SparseElimination     &M) 
    {
         // We make a copy as these data will be destroyed
        SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq> A1 (A);  
        return rankin(r, A1, tag, M);
    }
    
    template <class Field, class Method>
    inline unsigned long &rankin (unsigned long                   &r,
                           SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq>  &A,
                         const Method    &M)
    {  return rankin(r, A, typename FieldTraits<Field>::categoryTag(), M); 
    }


   template <class Blackbox, class Ring>
    inline unsigned long &rankin (unsigned long                     &r,
                           Blackbox &A,
                           const RingCategories::IntegerTag    &tag,
                           const Method::SparseElimination     &M)
    {
        commentator.start ("Integer Rank", "irank");
        typedef Modular<double> Field;
        integer mmodulus; 
        FieldTraits<Field>::maxModulus(mmodulus);
        RandomPrimeIterator genprime( (long) floor (log((double)mmodulus) ) );
        ++genprime;
        typedef typename Blackbox::template rebind< Field >::other FBlackbox;
        Field Fp(*genprime);
        FBlackbox Ap(A, Fp);
        commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Integer Rank is done modulo " << *genprime << std::endl;
        rankin(r, Ap, RingCategories::ModularTag(), M);
        commentator.stop ("done", NULL, "irank");
        return r;
    }
    
    template <class Field>
    inline unsigned long &rankin (unsigned long                       &r,
                           SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq>  &A,
                           const RingCategories::ModularTag    &tag,
                           const Method::SparseElimination     &M) 
    {
        commentator.start ("Sparse Elimination Rank", "serank");
        GaussDomain<Field> GD ( A.field() );
        GD.rankin (r, A, M.strategy ());
        commentator.stop ("done", NULL, "serank");
        return r;
    }

    inline unsigned long &rankin (unsigned long                       &r,
                           GaussDomain<GF2>::Matrix    &A,
                           const Method::SparseElimination     &)//M 
    {
        commentator.start ("Sparse Elimination Rank over GF2", "serankmod2");
        GaussDomain<GF2> GD ( A.field() );
        GD.rankin (r, A, Specifier::PIVOT_LINEAR);
        commentator.stop ("done", NULL, "serankmod2");
        return r;
    }

    inline unsigned long &rankin (unsigned long                       &r,
                           GaussDomain<GF2>::Matrix    &A,
                           const RingCategories::ModularTag    &,//tag
                           const Method::SparseElimination     &M) {
	return rankin(r, A, M);
    }

	// Change of representation to be able to call the sparse elimination
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                       &r,
                         const Blackbox  &A,
                         const RingCategories::ModularTag    &tag,
                         const Method::SparseElimination     &M) 
    {
        typedef typename Blackbox::Field Field;
        typedef SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq> SparseBB;
        SparseBB SpA(A.field(), A.rowdim(), A.coldim() );
	MatrixHom::map(SpA, A, A.field());
        return rankin(r, SpA, tag, M);
    }
    
	// M may be <code>Method::BlasElimination()</code>.
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                     &r,
                         const Blackbox                      &A,
                         const RingCategories::ModularTag          &tag,
                         const Method::BlasElimination  &M) 
    {
        
        commentator.start ("Blas Rank", "blasrank");
        typedef typename Blackbox::Field Field;
        const Field F = A.field();
        integer a, b; F.characteristic(a); F.cardinality(b);
        linbox_check( a == b );
        linbox_check( a < LinBox::BlasBound);
        BlasMatrix<typename Field::Element> B(A);
	BlasMatrixDomain<Field> D(F);
        r = D.rank(B);
	commentator.stop ("done", NULL, "blasrank");
        return r;
    }


// is this used?
	// A is modified.
    template <class Matrix>
    inline unsigned long &rankin (unsigned long                      &r,
                           Matrix                               &A,
                           const RingCategories::ModularTag          &tag,
                           const Method::SparseElimination &M) 
    {
        typedef typename Matrix::Field Field;
        const Field F = A.field();
        GaussDomain<Field> GD (F);
        GD.rankin( r, A, M.strategy ());
        return r;
    }

	/// A is modified.
    template <class Field>
    inline unsigned long &rankin (unsigned long                     &r,
                           BlasBlackbox<Field>                 &A,
                           const RingCategories::ModularTag          &tag,
                           const Method::BlasElimination  &M) 
    {

        commentator.start ("BlasBB Rank", "blasbbrank");
        const Field F = A.field();
        BlasMatrixDomain<Field> D(F);
        r = D.rankin(static_cast< BlasMatrix<typename Field::Element>& >(A));
        commentator.stop ("done", NULL, "blasbbrank");
        return r;
    }

    template <class Blackbox, class MyMethod>
    inline unsigned long &rank (unsigned long                     &r,
                         const Blackbox                      &A,
                         const RingCategories::IntegerTag          &tag,
                         const MyMethod                           &M)
    {
        commentator.start ("Integer Rank", "iirank");
        typedef Modular<double> Field;
        integer mmodulus; 
        FieldTraits<Field>::maxModulus(mmodulus);
        RandomPrimeIterator genprime( (long) floor (log((double)mmodulus) ) );
        ++genprime;
        typedef typename Blackbox::template rebind< Field >::other FBlackbox;
        Field Fp(*genprime);
        FBlackbox Ap(A, Fp );
        commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_DESCRIPTION) << "Integer Rank is done modulo " << *genprime << std::endl;
        
        rank(r, Ap, RingCategories::ModularTag(), M);
        commentator.stop ("done", NULL, "iirank");
        return r;
    }
} // LinBox


#ifdef __LINBOX_HAVE_GIVARO
#ifndef LINBOX_EXTENSION_DEGREE_MAX
#define LINBOX_EXTENSION_DEGREE_MAX 19
#endif

#include "linbox/field/givaro-extension.h"
namespace LinBox 
{
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const RingCategories::ModularTag                   &tag,
                         const Method::Blackbox& m)
    {  
        commentator.start ("BB Rank", "extend");
        if (m.certificate()) {
            typedef typename Blackbox::Field Field;
            const Field& F = A.field();
            integer a,c; F.cardinality(a); F.characteristic(c);
            if (a != c) {
                unsigned long extend = (unsigned long)FF_EXPONENT_MAX(a,(integer)LINBOX_EXTENSION_DEGREE_MAX);
                if (extend > 1) {
                    commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Extension of degree " << extend << std::endl;
                    GivaroExtension<Field> EF( F, extend);
                    typedef typename Blackbox::template rebind< GivaroExtension<Field>  >::other FBlackbox;
                    FBlackbox Ap(A, EF);
                    rank(r, Ap, tag, Method::Wiedemann(m));
                } else
                    rank(r, A, tag, Method::Wiedemann(m)); 
            } else {
                unsigned long extend = (unsigned long)FF_EXPONENT_MAX(c,(integer)LINBOX_EXTENSION_DEGREE_MAX);
                if (extend > 1) {
                    commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Word size extension : " << extend << std::endl;
                    GivaroGfq EF( (unsigned long)c, extend);                    
                    typedef typename Blackbox::template rebind< GivaroGfq >::other FBlackbox;
                    FBlackbox Ap(A, EF);
                    rank(r, Ap, tag, Method::Wiedemann(m));
                } else
                    rank(r, A, tag, Method::Wiedemann(m)); 
            }
        } else
            rank(r, A, tag, Method::Wiedemann(m)); 
        commentator.stop ("done", NULL, "extend");
        return r; 
    }
}
#else
namespace LinBox 
{
    template <class Blackbox>
    inline unsigned long &rank (unsigned long                   &r,
                         const Blackbox                  &A,
                         const  RingCategories::ModularTag                   &tag,
                         const Method::Blackbox& m)
    {  
            return rank(r, A, tag, Method::Wiedemann(m)); 
    }
}
#endif



#endif // __LINBOX_rank_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
