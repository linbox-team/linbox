/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/wiedemann.inl
 * Copyright (C) 2002 Zhendong Wan
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-10-02  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Refactoring:
 * Put everything inside a WiedemannSolver class, with the following
 * interface:
 *    solve - Solve a general linear system
 *    solveNonsingular - Solve a nonsingular system
 *    solveSingular - Solve a general singular system
 *    findRandomSolution - Find a random solution to a singular preconditioned
 *                         problem
 *    findNullspaceElement - Find an element of the right nullspace
 *    certifyInconsistency - Find a certificate of inconsistency for a
 *                           linear system
 *    precondition - Form a preconditioner and return it
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

#ifndef __WIEDEMANN_INL
#define __WIEDEMANN_INL

#include <vector>
#include <algorithm>

#include "linbox/algorithms/wiedemann.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h" 
#include "linbox/switch/cekstv.h"
#include "linbox/solutions/rank.h"
#include "linbox/vector/stream.h"

namespace LinBox 
{

template <class Field, class Vector>
Vector &WiedemannSolver<Field, Vector>::solve (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b)
{
	linbox_check ((x.size () == A.coldim ()) &&
		      (b.size () == A.rowdim ()));
	linbox_check (_traits.singular () != SolverTraits::NONSINGULAR || A.coldim () == A.rowdim ());

	commentator.start ("Solving linear system (Wiedemann)", "WiedemannSolver::solve");

	SolverTraits::SingularState singular = _traits.singular ();
	bool done = false;

	unsigned int tries = _traits.maxTries ();

	while (!done && tries-- > 0) {
		switch (singular) {
		    case SolverTraits::UNKNOWN:
		    {
			ActivityState state = commentator.saveActivityState ();

			try {
				done = true;
				solveNonsingular (A, x, b);
			}
			catch (SolveFailed) {
				done = false;
				commentator.restoreActivityState (state);
			}
			catch (SingularSystem) {
				done = false;
				commentator.restoreActivityState (state);

				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
					<< "System found to be singular. Reverting to nonsingular solver." << endl;
				singular = SolverTraits::SINGULAR;
			}
			break;
		    }

		    case SolverTraits::NONSINGULAR:
		    {
			ActivityState state = commentator.saveActivityState ();

			try {
				done = true;
				solveNonsingular (A, x, b);
			}
			catch (SolveFailed) {
				done = false;
				commentator.restoreActivityState (state);
			}
			break;
		    }

		    case SolverTraits::SINGULAR:
		    {
			ActivityState state = commentator.saveActivityState ();

			try {
				done = true;
				solveSingular (A, x, b);
			}
			catch (SolveFailed) {
				done = false;
				commentator.restoreActivityState (state);
			}
			break;
		    }
		}
	}

	if (!done) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Maximum tries exceeded with no resolution. Giving up." << endl;
		throw SolveFailed ();
	}

	commentator.stop ("done", NULL, "WiedemannSolver::solve");

	return x;
}

template <class Field, class Vector>
Vector &WiedemannSolver<Field, Vector>::solveNonsingular (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b)
{
	typedef std::vector<typename Field::Element> Polynomial;
	typedef BlackboxContainer<Field, Vector>     BBContainer;
	typedef MasseyDomain<Field, BBContainer>     MDomain;
	typedef typename Field::RandIter             RandIter;
	typedef typename Polynomial::iterator        PolyIterator;

	commentator.start ("Solving nonsingular system (Wiedemann)", "WiedemannSolver::solveNonsingular");

	Polynomial m_A;
	Vector     z;
	bool       ret = false;

	{
		commentator.start ("Computing minimal polynomial");

		RandIter       r (_F);
		BBContainer    TF (&A, _F, r);
		MDomain        WD (&TF);
		unsigned long  deg;

		WD.minpoly (m_A, deg);

		commentator.stop ("done");
	}

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Minimal polynomial coefficients: ";
	_VD.write (report, m_A) << endl;

	if (_F.isZero (m_A.front ()))
		throw SingularSystem ();

	{
		commentator.start ("Preparing polynomial for application");

		PolyIterator iter = m_A.begin ();

		while (++iter != m_A.end ()) {
			_F.divin (*iter, m_A.front ());
			_F.negin (*iter);
		}

		commentator.stop ("done");
	}

	{
		commentator.start ("Applying polynomial via Horner's rule", NULL, m_A.size () - 1);

		_VD.mul (x, b, m_A.back ());

		VectorWrapper::ensureDim (z, A.rowdim ());

		for (int i = m_A.size () - 1; --i > 0;) {
			if ((m_A.size () - i) & 0xff == 0)
				commentator.progress (m_A.size () - i);

			A.apply (z, x);
			_VD.axpy (x, m_A[i], b, z);
		}

		commentator.stop ("done");
	}

	if (_traits.checkResult ()) {
		commentator.start ("Checking whether Ax=b");
		A.apply (z, x);

		if (_VD.areEqual (z, b))
			ret = true;
		else
			ret = false;

		commentator.stop (MSG_STATUS (ret));
	}

	commentator.stop (MSG_STATUS (ret), NULL, "WiedemannSolver::solveNonsingular");

	if (!ret)
		throw SolveFailed ();

	return x;
}

template <class Field, class Vector>
Vector &WiedemannSolver<Field, Vector>::solveSingular (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b)
{
	commentator.start ("Solving singular system (Wiedemann)", "WiedemannSolver::solveSingular");

	Vector u, Ax;
	unsigned long r;
	bool failed = false, cert = false;

	if (_traits.rank () != SolverTraits::UNKNOWN)
		r = _traits.rank ();
	else
		rank (r, A, _F);

	commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
		<< "Rank of A = " << r << endl;

	BlackboxArchetype<Vector> *P = NULL;
	BlackboxArchetype<Vector> *Q = NULL;
	BlackboxArchetype<Vector> *PAQ = NULL;
	const BlackboxArchetype<Vector> *B = precondition (A, PAQ, P, Q);

	try {
		findRandomSolution (*B, x, b, r, P, Q);
	}
	catch (BadPreconditioner) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Preconditioned matrix did not have generic rank profile" << endl;

		failed = true;
	}
	catch (SolveFailed) {
		if (_traits.certificate ()) {
			VectorWrapper::ensureDim (u, A.rowdim ());

			if (certifyInconsistency (u, *B, b, r, P))
				cert = true;
			else
				failed = true;
		} else
			failed = true;
	}

	if (!cert && !failed && _traits.checkResult ()) {
		commentator.start ("Checking system solution");

		VectorWrapper::ensureDim (Ax, A.rowdim ());

		A.apply (Ax, x);

		if (_VD.areEqual (Ax, b))
			commentator.stop ("passed");
		else {
			commentator.stop ("FAILED");

			if (_traits.certificate ()) {
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
					<< "Computed system solution is not correct. Attempting to find certificate of inconsistency." << endl;

				VectorWrapper::ensureDim (u, A.rowdim ());

				if (certifyInconsistency (u, *B, b, r, P))
					cert = true;
				else
					failed = true;
			} else
				failed = true;
		}
	}

	if (PAQ != NULL) delete PAQ;
	if (P != NULL) delete P;
	if (Q != NULL) delete Q;

	if (failed)
		throw SolveFailed ();
	if (cert)
		throw InconsistentSystem (u);

	commentator.stop ("done", NULL, "WiedemannSolver::solveSingular");

	return x;
}

template <class Field, class Vector>
Vector &WiedemannSolver<Field, Vector>::findRandomSolution (const BlackboxArchetype<Vector> &A,
							    Vector                          &x,
							    const Vector                    &b,
							    size_t                           r,
							    const BlackboxArchetype<Vector> *P,
							    const BlackboxArchetype<Vector> *Q)
{
	commentator.start ("Solving singular system with generic rank profile (Wiedemann)", "WiedemannSolver::findRandomSolution");

	Vector v, Avpb, PAvpb, bp, xp, Qinvx;

	RandomDenseStream<Field, Vector> stream (_F, A.coldim ());

	VectorWrapper::ensureDim (v, A.coldim ());
	VectorWrapper::ensureDim (Avpb, A.rowdim ());
	VectorWrapper::ensureDim (xp, r);
	VectorWrapper::ensureDim (bp, r);

	{
		commentator.start ("Preparing right hand side");

		stream >> v;

		A.apply (Avpb, v);
		_VD.addin (Avpb, b);

		if (P != NULL) {
			VectorWrapper::ensureDim (PAvpb, A.rowdim ());
			P->apply (PAvpb, Avpb);
			_VD.copy (bp, PAvpb, 0, r);
		} else {
			_VD.copy (bp, Avpb, 0, r);
		}

		commentator.stop ("done");
	}

	Submatrix<Field, Vector> Ap (_F, &A, 0, 0, r, r);

	try {
		solveNonsingular (Ap, xp, bp);
	}
	catch (SingularSystem) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Leading principal minor was found to be singular." << endl;
		throw BadPreconditioner ();
	}

	if (Q != NULL) {
		VectorWrapper::ensureDim (Qinvx, A.coldim ());
		_VD.copy (Qinvx, xp);
		Q->apply (x, Qinvx);
	} else {
		_VD.copy (x, xp);
	}

	_VD.subin (x, v);

	commentator.stop ("done", NULL, "WiedemannSolver::findRandomSolution");

	return x;
}

template <class Field, class Vector>
Vector &WiedemannSolver<Field, Vector>::findNullspaceElement (Vector                          &x,
							      const BlackboxArchetype<Vector> &A,
							      size_t                           r,
							      const BlackboxArchetype<Vector> *P,
							      const BlackboxArchetype<Vector> *Q)
{
	commentator.start ("Finding a nullspace element (Wiedemann)", "WiedemannSolver::findNullspaceElement");

	Vector v, Av, PAv, bp, xp, Qinvx;

	RandomDenseStream<Field, Vector> stream (_F, A.coldim ());

	VectorWrapper::ensureDim (v, A.coldim ());
	VectorWrapper::ensureDim (Av, A.rowdim ());
	VectorWrapper::ensureDim (bp, r);
	VectorWrapper::ensureDim (xp, r);

	{
		commentator.start ("Constructing right hand side");

		stream >> v;

		A.apply (Av, v);

		if (P != NULL) {
			VectorWrapper::ensureDim (PAv, A.rowdim ());
			P->apply (PAv, Av);
			_VD.copy (bp, PAv, 0, r);
		} else {
			_VD.copy (bp, Av, 0, r);
		}

		commentator.stop ("done");
	}

	Submatrix<Field, Vector> Ap (_F, &A, 0, 0, r, r);

	try {
		solveNonsingular (Ap, xp, bp);
	}
	catch (SingularSystem) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Leading principal minor was found to be singular." << endl;
		throw BadPreconditioner ();
	}

	if (Q != NULL) {
		VectorWrapper::ensureDim (Qinvx, A.coldim ());
		_VD.copy (Qinvx, xp);
		Q->apply (x, Qinvx);
	} else {
		_VD.copy (x, xp);
	}

	_VD.subin (x, v);

	commentator.stop ("done", NULL, "WiedemannSolver::findNullspaceElement");

	return x;
}

template <class Field, class Vector>
bool WiedemannSolver<Field, Vector>::certifyInconsistency (Vector                          &u,
							   const BlackboxArchetype<Vector> &A,
							   const Vector                    &b,
							   size_t                           r,
							   const BlackboxArchetype<Vector> *P)
{
	commentator.start ("Obtaining certificate of inconsistency (Wiedemann)", "WiedemannSolver::certifyInconsistency");

	Vector PTinvu;
	typename Field::Element uTb;

	SolverTraits cert_traits;

	bool ret = false;

	cert_traits.method (SolverTraits::WIEDEMANN);
	cert_traits.preconditioner (SolverTraits::NONE);
	cert_traits.rank (r);
	cert_traits.certificate (false);
	cert_traits.singular (SolverTraits::SINGULAR);
	cert_traits.maxTries (1);

	WiedemannSolver solver (_F, cert_traits);

	Transpose<Vector> AT (&A);

	if (P != NULL) {
		VectorWrapper::ensureDim (PTinvu, A.rowdim ());
		ActivityState state = commentator.saveActivityState ();

		bool failed = false;

		try {
			solver.findNullspaceElement (PTinvu, AT, r);
		}
		catch (BadPreconditioner) {
			commentator.restoreActivityState (state);
			failed = true;
			ret = false;
		}

		if (!failed) {
			_VD.dot (uTb, PTinvu, b);

			if (!_F.isZero (uTb)) {
				ret = true;
				VectorWrapper::ensureDim (u, A.coldim ());
				P->applyTranspose (u, PTinvu);
			}
		}
	} else {
		solver.findNullspaceElement (u, AT, r);
		_VD.dot (uTb, u, b);

		if (!_F.isZero (uTb))
			ret = true;
	}

	commentator.stop (MSG_STATUS (ret), NULL, "WiedemannSolver::certifyInconsistency");

	return ret;
}

template <class Field, class Vector>
const BlackboxArchetype<Vector> *WiedemannSolver<Field, Vector>::precondition (const BlackboxArchetype<Vector>  &A,
									       BlackboxArchetype<Vector>       *&PAQ,
									       BlackboxArchetype<Vector>       *&P,
									       BlackboxArchetype<Vector>       *&Q)
{
	switch (_traits.preconditioner ()) {
	    case SolverTraits::BUTTERFLY:
	    {
		    commentator.start ("Constructing butterfly preconditioner");

		    typename Field::RandIter r (_F);
		    CekstvSwitchFactory<Field> factory (r);
		    P = new Butterfly<Field, CekstvSwitch<Field> > (_F, A.rowdim (), factory);
		    Q = new Butterfly<Field, CekstvSwitch<Field> > (_F, A.coldim (), factory);
		    Compose<Vector> AQ (&A, Q);
		    PAQ = new Compose<Vector> (P, &AQ);

		    commentator.stop ("done");
		    break;
	    }

	    case SolverTraits::SPARSE:
	    {
		    commentator.start ("Constructing sparse preconditioner");

		    SparseMatrix0<Field> *QT;
		    SparseMatrix0<Field> *Q_sparse = new SparseMatrix0<Field> (_F, A.coldim (), A.coldim ());
		    P = makeLambdaSparseMatrix (A.rowdim ());
		    QT = makeLambdaSparseMatrix (A.coldim ());

		    QT->transpose (*Q_sparse);
		    Q = Q_sparse;

		    Compose<Vector> AQ (&A, Q);
		    PAQ = new Compose<Vector> (P, &AQ);

		    commentator.stop ("done");
		    break;
	    }

	    case SolverTraits::TOEPLITZ:
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Toeplitz preconditioner not implemented yet. Sorry." << endl;

	    case SolverTraits::NONE:
		return &A;
	}

	return PAQ;
}

template <class Field, class Vector>
SparseMatrix0<Field, Vector> *WiedemannSolver<Field, Vector>::makeLambdaSparseMatrix (size_t m)
{
	const double             LAMBDA = 3;
	integer                  card;

	_F.cardinality (card);

	typename Field::RandIter r (_F);
	NonzeroRandIter<Field>   rp (_F, r);
	double                   init_p = 1.0 - 1.0 / (double) card;
	double                   log_m = LAMBDA * log ((double) m) / M_LN2;
	double                   new_p;

	SparseMatrix0<Field>    *P = new SparseMatrix0<Field> (_F, m, m);

	RandomSparseStream<Field, typename LinBox::Vector<Field>::Sparse, NonzeroRandIter<Field> >
		stream (_F, rp, m, init_p, m);

	for (unsigned int i = 0; i < m; ++i) {
		new_p = log_m / (m - i + 1);

		if (init_p < new_p)
			stream.setP (init_p);
		else
			stream.setP (new_p);

		stream >> P->getRow (i);
	}

	return P;
}

}

#endif // __WIEDEMANN_INL
