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
#include "linbox/blackbox/butterfly.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h" 
#include "linbox/switch/cekstv.h"
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/rank.h"
#include "linbox/vector/stream.h"

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
	typename Polynomial::iterator iter;
	bool                          done = false;

	while (!done) {
		typename Field::RandIter r (F);

		BlackboxContainer<Field, Vector> TF (&A, F, r);
		MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF);
		
		WD.minpoly (P, deg);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Minimal polynomial coefficients: ";
		VD.write (report, P);
		report << endl;
			
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

				VectorWrapper::ensureDim (z, A.rowdim ());

				for (int i = P.size () - 1; --i > 0;) {
					if ((P.size () - i) & 0xff == 0)
						commentator.progress (P.size () - i);

					A.apply (z, x);
					VD.axpy (x, P[i], b, z);
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
						done = true;
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
	Vector                         z, v, Av, y_wp, bp, BAvpb;
	RandomDenseStream<Field>       stream (F, A.coldim ());
	unsigned long                  deg;
	VectorDomain <Field>           VD (F);
	long unsigned int              A_rank;
	typename Polynomial::iterator  iter;
	BlackboxArchetype<Vector>     *Ap, *B = NULL, *BA = NULL;
	CekstvSwitch<Field>           *s = NULL;
	typename Field::RandIter       r (F);
	bool                           done = false;
	bool                           first_iter = true;
	int                            tries = 0;

	VectorWrapper::ensureDim (v, A.rowdim ());
	VectorWrapper::ensureDim (Av, A.rowdim ());

	while (!done) {
		if (++tries > traits.maxTries ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Maximum tries exceeded with no resolution. Giving up." << endl;
			break;
		}

		if (traits.rank () == SolverTraits::RANK_UNKNOWN || first_iter == false)
			rank (A_rank, A, F);
		else
			A_rank = traits.rank ();

		if (traits.precondition ()) {
			commentator.start ("Constructing butterfly preconditioner");

			s = new CekstvSwitch<Field> (F, r);
			B = new Butterfly<Vector, CekstvSwitch<Field> > (A.rowdim (), *s);
			BA = new Compose<Vector> (B, &A);
			Ap = new Submatrix<Vector> (BA, 0, 0, A_rank, A_rank);

			commentator.stop ("done");
		} else {
			s = NULL;
			B = NULL;
			BA = NULL;
			Ap = new Submatrix<Vector> (&A, 0, 0, A_rank, A_rank);
		}

		commentator.start ("Solving system with nonsingular leading principal minor");

		{
			BlackboxContainer<Field, Vector> TF (Ap, F, v);
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

				delete Ap;
				if (BA != NULL) delete BA;
				if (B != NULL) delete B;
				if (s != NULL) delete s;
				continue;
			}

			commentator.start ("Preparing polynomial for application");

			iter = P.begin ();

			while (++iter != P.end ()) {
				F.divin (*iter, P.front ());
				F.negin (*iter);
			}

			commentator.stop ("done");
		}

		{
			commentator.start ("Preparing right hand side");

			stream >> v;

			A.apply (Av, v);
			VD.addin (Av, b);

			if (traits.precondition ()) {
				VectorWrapper::ensureDim (BAvpb, A.rowdim ());
				B->apply (BAvpb, Av);
			}

			commentator.stop ("done");
		}
		
		{
			commentator.start ("Applying polynomial via Horner's rule", NULL, P.size () - 1);

			// Make sure x is 0 first
			VD.subin (x, x);

			VectorWrapper::ensureDim (y_wp, A_rank);
			VectorWrapper::ensureDim (bp, A_rank);
			VectorWrapper::ensureDim (z, A_rank);

			if (traits.precondition ())
				VD.copy (bp, BAvpb, 0, A_rank);
			else
				VD.copy (bp, Av, 0, A_rank);

			VD.mul (y_wp, bp, P.back ());

			for (int i = P.size () - 1; --i > 0;) {
				if ((P.size () - i) & 0xff == 0)
					commentator.progress (P.size () - i);

				Ap->apply (z, y_wp);
				VD.axpy (y_wp, P[i], bp, z);
			}

			VD.copy (x, y_wp, 0, A_rank);
			VD.subin (x, v);

			commentator.stop ("done");
		}

		commentator.stop ("done");
		
		{
			commentator.start ("Deallocating space");

			delete Ap;
			if (BA != NULL) delete BA;
			if (B != NULL) delete B;
			if (s != NULL) delete s;

			commentator.stop ("done");
		}

		if (traits.checkResult ()) {
			commentator.start ("Checking whether Ax=b");
			A.apply (v, x);

			if (VD.areEqual (v, b)) {
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
