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
#include "linbox/blackbox/transpose.h"
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

/** Exception thrown when the system to be solved is inconsistent. Contains a
 * certificate of inconsistency.
 */

template <class Vector>
class InconsistentSystem 
{
    public:
	InconsistentSystem (Vector &u)
		: _cert (true), _u (u)
	{}

	InconsistentSystem ()
		: _cert (false)
	{}

	const Vector &u () const { return _u; }
	bool certified () const { return _cert; }

    private:

	bool _cert;
	Vector _u;
};

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

	Polynomial                     m_A;
	Vector                         z, v, Av, y_wp, bp, BAvpb, Qx, PTv;
	RandomDenseStream<Field>       stream (F, A.coldim ());
	unsigned long                  deg;
	VectorDomain <Field>           VD (F);
	long unsigned int              A_rank;
	typename Polynomial::iterator  iter;
	BlackboxArchetype<Vector>     *Ap = NULL, *P = NULL, *Q = NULL, *PAQ = NULL;
	typename Field::RandIter       r (F);
	bool                           done = false;
	bool                           first_iter = true;
	bool                           certificate = false;
	int                            tries = 0;

	typename LinBox::Vector<Field>::Dense switch_v (1);

	VectorWrapper::ensureDim (v, A.coldim ());
	VectorWrapper::ensureDim (Av, A.rowdim ());

	while (!done) {
		if (++tries > traits.maxTries ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
				<< "Maximum tries exceeded with no resolution. Giving up and assuming inconsistency." << endl;
			throw InconsistentSystem<Vector> ();
		}

		if (traits.rank () == SolverTraits::RANK_UNKNOWN || first_iter == false) {
			rank (A_rank, A, F);
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
				<< "Computed rank (A) = " << A_rank << endl;
		} else
			A_rank = traits.rank ();

		switch (traits.preconditioner ()) {
		    case SolverTraits::BUTTERFLY:
		    {
			commentator.start ("Constructing butterfly preconditioner");

			CekstvSwitchFactory<Field> factory (r);
			P = new Butterfly<Field, CekstvSwitch<Field> > (F, A.rowdim (), factory);
			PAQ = new Compose<Vector> (P, &A);
			Ap = new Submatrix<Field, Vector> (F, PAQ, 0, 0, A_rank, A_rank);

			commentator.stop ("done");
			break;
		    }

		    case SolverTraits::SPARSE:
		    {
			commentator.start ("Constructing sparse preconditioner");

			integer card;
			F.cardinality (card);

			const double LAMBDA = 3;

			NonzeroRandIter<Field> rp (F, r);
			double init_p = 1.0 - 1.0 / (double) card;

			SparseMatrix0<Field> *P_sparse = new SparseMatrix0<Field> (F, A.rowdim (), A.rowdim ());
			SparseMatrix0Base<typename Field::Element> Q_base (A.coldim (), A.coldim ());
			SparseMatrix0<Field> *Q_sparse = new SparseMatrix0<Field> (F, A.coldim (), A.coldim ());

			RandomSparseStream<Field, typename LinBox::Vector<Field>::Sparse, NonzeroRandIter<Field> >
				P_stream (F, rp, A.rowdim (), init_p, A.rowdim ());

			double log_m = LAMBDA * log ((double) A.rowdim ()) / M_LN2;

			for (unsigned int i = 0; i < A.rowdim (); ++i) {
				double new_p = log_m / (A.rowdim () - i + 1);
				if (init_p < new_p)
					P_stream.setP (init_p);
				else
					P_stream.setP (new_p);

				P_stream >> P_sparse->getRow (i);
			}

			RandomSparseStream<Field, typename LinBox::Vector<Field>::Sparse, NonzeroRandIter<Field> >
				Q_stream (F, rp, A.coldim (), init_p, A.rowdim ());

			log_m = LAMBDA * log ((double) A.coldim ()) / M_LN2;

			for (unsigned int i = 0; i < A.coldim (); ++i) {
				double new_p = log_m / (A.coldim () - i + 1);
				if (init_p < new_p)
					Q_stream.setP (init_p);
				else
					Q_stream.setP (new_p);

				Q_stream >> Q_base.getRow (i);
			}

			Q_base.transpose (*Q_sparse);

			Q = Q_sparse;
			P = P_sparse;

			Compose<Vector> AQ (&A, Q);
			PAQ = new Compose<Vector> (P, &AQ);
			Ap = new Submatrix<Field, Vector> (F, PAQ, 0, 0, A_rank, A_rank);

			commentator.stop ("done");
			break;
		    }

		    case SolverTraits::TOEPLITZ:
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Toeplitz preconditioner not implemented yet. Sorry." << endl;

		    case SolverTraits::NONE:
		    {
			P = Q = PAQ = NULL;
			Ap = new Submatrix<Field, Vector> (F, &A, 0, 0, A_rank, A_rank);
			break;
		    }
		}

		{
			commentator.start ("Preparing right hand side");

			stream >> v;

			A.apply (Av, v);
			VD.addin (Av, b);

			if (P != NULL) {
				VectorWrapper::ensureDim (BAvpb, P->rowdim ());
				P->apply (BAvpb, Av);
			}

			commentator.stop ("done");
		}

		commentator.start ("Solving system with nonsingular leading principal minor");

		{
			BlackboxContainer<Field, Vector> TF (Ap, F, Av);
			MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF);
		
			WD.minpoly (m_A, deg);

			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
			report << "Minimal polynomial coefficients: ";
			VD.write (report, m_A) << endl;

			if (F.isZero (m_A.front ())) {
				first_iter = false;
				if (traits.rank () != SolverTraits::RANK_UNKNOWN && traits.preconditioner () == SolverTraits::NONE)
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
						<< "Rank specified in SolverTraits is incorrect. Recomputing." << endl;
				else
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
						<< "Leading minor was singular. Assuming to be a bad preconditioner." << endl;

				commentator.stop ("done");

				delete Ap;
				if (P != NULL) delete P;
				if (Q != NULL) delete Q;
				if (PAQ != NULL) delete PAQ;
				continue;
			}

			commentator.start ("Preparing polynomial for application");

			iter = m_A.begin ();

			while (++iter != m_A.end ()) {
				F.divin (*iter, m_A.front ());
				F.negin (*iter);
			}

			commentator.stop ("done");
		}
		
		{
			commentator.start ("Applying polynomial via Horner's rule", NULL, m_A.size () - 1);

			// Make sure x is 0 first
			VD.subin (x, x);

			VectorWrapper::ensureDim (y_wp, A_rank);
			VectorWrapper::ensureDim (bp, A_rank);
			VectorWrapper::ensureDim (z, A_rank);

			if (traits.preconditioner () != SolverTraits::NONE)
				VD.copy (bp, BAvpb, 0, A_rank);
			else
				VD.copy (bp, Av, 0, A_rank);

			VD.mul (y_wp, bp, m_A.back ());

			for (int i = m_A.size () - 1; --i > 0;) {
				if ((m_A.size () - i) & 0xff == 0)
					commentator.progress (m_A.size () - i);

				Ap->apply (z, y_wp);
				VD.axpy (y_wp, m_A[i], bp, z);
			}

			VD.copy (x, y_wp, 0, A_rank);

			if (Q != NULL) {
				VectorWrapper::ensureDim (Qx, A.coldim ());
				Q->apply (Qx, x);
				VD.copy (x, Qx);
			}

			VD.subin (x, v);

			commentator.stop ("done");
		}

		commentator.stop ("done");

		if (traits.checkResult ()) {
			commentator.start ("Checking whether Ax=b");
			A.apply (Av, x);

			if (VD.areEqual (Av, b)) {
				commentator.stop ("yes");
				done = true;
			} else {
				commentator.stop ("no");

				if (traits.certificate ()) {
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
						<< "System is likely inconsistent. Trying to produce certificate..." << endl;

					commentator.start ("Obtaining certificate of inconsistency");

					// Get a random solution to (BA)^T v = 0; i.e. a random element of the right nullspace of (BA)^T
					VD.subin (Av, Av);

					SolverTraits cert_traits;

					cert_traits.method (SolverTraits::WIEDEMANN);
					cert_traits.preconditioner (SolverTraits::NONE);
					cert_traits.rank (A_rank);
					cert_traits.certificate (false);
					cert_traits.singular (SolverTraits::SINGULAR);
					cert_traits.maxTries (1);

					typename Field::Element vTb;

					ActivityState state = commentator.saveActivityState ();

					try {
						if (traits.preconditioner () != SolverTraits::NONE) {
							Transpose<Vector> PAQT (PAQ);
							solveWiedemannSingular (PAQT, v, Av, F, cert_traits);
						} else {
							Transpose<Vector> AT (&A);
							solveWiedemannSingular (A, v, Av, F, cert_traits);
						}

						VD.dot (vTb, v, b);

						if (!F.isZero (vTb)) {
							commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
								<< "Certificate found." << endl;
							certificate = true;

							if (P != NULL) {
								VectorWrapper::ensureDim (PTv, A.rowdim ());
								P->applyTranspose (PTv, v);
								VD.copy (v, PTv);
							}
						} else
							commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
								<< "Certificate not found. Continuing with next iteration..." << endl;
					}
					catch (InconsistentSystem<Vector> e) {
						commentator.restoreActivityState (state);
						commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
							<< "Could not find a vector in the right nullspace" << endl;
					}

					commentator.stop ("done");
				}
			}
		}
		else
			done = true;

		delete Ap;
		if (P != NULL) delete P;
		if (Q != NULL) delete Q;
		if (PAQ != NULL) delete PAQ;

		if (certificate)
			throw InconsistentSystem<Vector> (v);
	}

	commentator.stop ("done", NULL, "solveWiedemannSingular");

	return x;
}								       		

}

#endif // __WIEDEMANN_H
