/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/wiedemann.h
 * Copyright (C) 2002 Zhendong Wan
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Move the Wiedemann stuff over to this file
 *
 * Create a singular and nonsingular version that is a bit intelligent about
 * which one to use in different circumstances
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __WIEDEMANN_H
#define __WIEDEMANN_H

#include <vector>
#include <algorithm>

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h" 
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/rank.h"

namespace LinBox 
{

template <class Field, class Vector>
Vector &solveWiedemann (const BlackboxArchetype<Vector> &A,
			Vector                          &x,		       
			const Vector                    &b,
			const Field                     &F,
			const SolverTraits              &traits)
{
	typedef std::vector<typename Field::Element> Polynomial;

	// FIXME: This is not generic wrt vector type. Do we care?
	linbox_check ((x.size () == A.coldim ()) &&
		      (b.size () == A.rowdim ()));
	linbox_check (traits.singular () != SolverTraits::NONSINGULAR || A.coldim () == A.rowdim ());

	if (traits.singular () == SolverTraits::SINGULAR || A.coldim () != A.rowdim ())
		solveWiedemannSingular (A, x, b, F, traits);

	commentator.start ("Solving system (Wiedemann)", "solveWiedemann");

	Polynomial                    P;
	Vector                        z;
	unsigned long                 deg;
	VectorDomain <Field>          VD (F);		
	typename Polynomial::iterator iter, iter_end;
	size_t                        count = 0;
	bool                          done = false;

	while (!done) {
		BlackboxContainer<Field, Vector> TF (&A, F, b);
		MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF);
		
		WD.minpoly (P, deg);

		// Nonsingular case
		if (!F.isZero (P.front ())) {
			if (traits.singular () == SolverTraits::UNKNOWN)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
					<< "System found to be nonsingular" << endl;
			
			{
				commentator.start ("Preparing polynomial for application");

				iter = P.begin ();

				while (++iter != P.end ()) {
					F.divin (*iter, P.front ());
					F.negin (*iter);
				}

				commentator.stop ("done");
			}
			
			{
				commentator.start ("Applying polynomial via Horner's rule", NULL, P.size () - 1);

				VD.mul (x, b, P.back ());

				iter_end = P.begin ();
				iter_end++;

				VectorWrapper::ensureDim (z, A.rowdim ());

				iter = P.end ();
				--iter;

				while (--iter > iter_end) {
					if (++count % 100 == 0)
						commentator.progress (count);

					A.apply (z, x);
					VD.axpy (x, *iter, z, b);
				}

				commentator.stop ("done");
			}

			if (traits.checkResult ()) {
				commentator.start ("Checking whether Ax=b");
				A.apply (z, x);

				if (VD.areEqual (z, b)) {
					commentator.stop ("yes");
					done = true;
				} else {
					commentator.stop ("no");

					if (P.size () <= A.rowdim ())
						commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
							<< "Minimal polynomial has low degree. Recomputing." << endl;
					else {
						commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
							<< "ERROR: Minimal polynomial has full degree. There's a bug here." << endl;
						return x;
					}
				}
			}
			else
				done = true;
		}

		// Singular case
		else {
			if (traits.singular () == SolverTraits::UNKNOWN)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
					<< "System found to be singular. Giving random solution." << endl;
			else if (traits.singular () == SolverTraits::NONSINGULAR)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
					<< "System found to be singular. Reverting to singular system solver." << endl;

			solveWiedemannSingular (A, x, b, F, traits);
			done = true;
		}
	}

	commentator.stop ("done", NULL, "solveWiedemann");

	return x;
}								       		

template <class Field, class Vector>
Vector &solveWiedemannSingular (const BlackboxArchetype<Vector> &A,
				Vector                          &x,		       
				const Vector                    &b,
				const Field                     &F,
				const SolverTraits              &traits)
{
	typedef std::vector<typename Field::Element> Polynomial;

	// FIXME: This is not generic wrt vector type. Do we care?
	linbox_check ((x.size () == A.coldim ()) &&
		      (b.size () == A.rowdim ()));

	commentator.start ("Solving singular system (Wiedemann)", "solveWiedemannSingular");

	if (!traits.checkResult ())
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
			<< "Not checking the result for a singular system! You will get incorrect results for singular systems." << endl;

	Polynomial                     P;
	Vector                         z;
	unsigned long                  deg;
	VectorDomain <Field>           VD (F);		
	long unsigned int              A_rank;
	typename Polynomial::iterator  iter, iter_end;
	BlackboxArchetype<Vector>     *Ap;
	bool                           done = false;
	bool                           first_iter = true;
	size_t                         count = 0;

	while (!done) {
		if (traits.rank () == SolverTraits::RANK_UNKNOWN || first_iter == false)
			rank (A_rank, A, F);
		else
			A_rank = traits.rank ();

		if (traits.precondition ()) {
			commentator.start ("Constructing butterfly preconditioner");

			// FIXME
			Ap = new Submatrix<Vector> (&A, 0, 0, A_rank, A_rank);

			commentator.stop ("done");
		}
		else 
			Ap = new Submatrix<Vector> (&A, 0, 0, A_rank, A_rank);
		
		{
			commentator.start ("Solving system with nonsingular leading principal minor");

			BlackboxContainer<Field, Vector> TF (Ap, F, b);
			MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF);
		
			WD.minpoly (P, deg);

			if (F.isZero (P.front ())) {
				first_iter = false;
				if (traits.rank () != SolverTraits::RANK_UNKNOWN && !traits.precondition ())
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
						<< "Rank specified in SolverTraits is incorrect. Recomputing.";
				else
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
						<< "Leading minor was singular. Assuming to be a bad preconditioner.";
			}

			commentator.start ("Preparing polynomial for application");

			iter = P.begin ();

			while (++iter != P.end ()) {
				F.divin (*iter, P.front ());
				F.negin (*iter);
			}

			commentator.stop ("done");
		}

		Vector Bb;

		if (traits.precondition ()) {
			commentator.start ("Preparing right hand sie");

			VectorWrapper::ensureDim (Bb, A.rowdim ());
			VD.copy (Bb, b);
		}
		else
			VD.copy (Bb, b);
		
		{
			commentator.start ("Applying polynomial via Horner's rule", NULL, P.size () - 1);

			// Make sure x is 0 first
			VD.subin (x, x);

			// FIXME: This code depends heavily on the vectors being dense. Do we really care?

			Vector y_wp (A_rank);
			Vector bp (A_rank);

			VD.copy (bp, b, 0, A_rank);
			VD.mul (y_wp, bp, P.back ());

			iter_end = P.begin ();
			iter_end++;

			VectorWrapper::ensureDim (z, A_rank);

			for (iter = P.end () - 2; iter != iter_end; iter--) {
				if (++count % 100 == 0)
					commentator.progress (count);

				Ap->apply (z, y_wp);
				VD.axpy (y_wp, *iter, z, bp);
			}

			VD.copy (x, y_wp, 0, A_rank);

			commentator.stop ("done");
		}

		delete Ap;

		{
			commentator.start ("Getting random solution");

			RandomDenseVectorFactory<Field> factory (F, A.coldim ());

			Vector v;

			VectorWrapper::ensureDim (v, A.coldim ());
			factory.next (v);

			VD.subin (x, v);

			commentator.stop ("done");
		}

		if (traits.checkResult ()) {
			commentator.start ("Checking whether Ax=b");
			A.apply (z, x);

			if (VD.areEqual (z, b)) {
				commentator.stop ("yes");
				done = true;
			} else {
				commentator.stop ("no");

				if (traits.certificate ()) {
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
						<< "System is likely inconsistent. Trying to produce certificate..." << endl;
					// FIXME - We need a certificate of inconsistency
				}
			}
		}
		else
			done = true;
	}

	commentator.stop ("done", NULL, "solveWiedemannSingular");

	return x;
}								       		

}

#endif // __WIEDEMANN_H
