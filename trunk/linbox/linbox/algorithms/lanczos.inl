/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/lanczos.inl
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LANCZOS_INL
#define __LANCZOS_INL

#include <vector>
#include <algorithm>

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

template <class Field, class Vector>
Vector &LanczosSolver<Field, Vector>::solve (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b) 
{
	linbox_check ((x.size () == A.coldim ()) &&
		      (b.size () == A.rowdim ()));
	linbox_check (!_traits.symmetric () || A.coldim () == A.rowdim ());

	commentator.start ("Solving linear system (Lanczos)", "LanczosSolver::solve");

	bool success;
	Vector d1, d2, b1, b2, bp, y;

	VectorWrapper::ensureDim (_w[0], A.coldim ());
	VectorWrapper::ensureDim (_w[1], A.coldim ());
	VectorWrapper::ensureDim (_Aw, A.coldim ());
	VectorWrapper::ensureDim (d1, A.coldim ());
	VectorWrapper::ensureDim (y, A.coldim ());
	VectorWrapper::ensureDim (bp, A.rowdim ());

	for (int i = 0; i < _traits.maxTries (); ++i) {
		RandomDenseStream<Field, Vector> stream (_F, _randiter, A.coldim ());

		if (_traits.symmetric ()) {
			stream >> d1;
			Diagonal<Field, Vector> D (_F, d1);
			Compose<Vector> B (&A, &D);

			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
			report << "Random D: ";
			_VD.write (report, d1) << endl;

			if ((success = iterate (B, y, b)))
				D.apply (x, y);
		} else {
			VectorWrapper::ensureDim (d2, A.coldim ());
			VectorWrapper::ensureDim (b1, A.rowdim ());
			VectorWrapper::ensureDim (b2, A.coldim ());

			stream >> d1 >> d2;
			Diagonal<Field, Vector> D1 (_F, d1);
			Diagonal<Field, Vector> D2 (_F, d2);
			Compose<Vector> B1 (&A, &D1);
			Compose<Vector> B2 (&D2, &B1);
			Transpose<Vector> AT (&A);
			Compose<Vector> B3 (&AT, &B2);
			Compose<Vector> B (&D1, &B3);

			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
			report << "Random D_1: ";
			_VD.write (report, d1) << endl;

			commentator.indent (report);
			report << "Random D_2: ";
			_VD.write (report, d2) << endl;

			D2.apply (b1, b);
			AT.apply (b2, b1);
			D1.apply (bp, b2);

			if ((success = iterate (B, y, bp)))
				D1.apply (x, y);
		}

		if (success) {
			if (_traits.checkResult ()) {
				commentator.start ("Checking whether Ax=b");

				A.apply (bp, x);

				if (_VD.areEqual (bp, b)) {
					commentator.stop ("passed");
					commentator.stop ("done", "Solve successful", "LanczosSolver::solve");
					return x;
				} else
					commentator.stop ("FAILED");
			} else {
				commentator.stop ("done", "Solve successful", "LanczosSolver::solve");
				return x;
			}
		}
	}

	commentator.stop ("done", "Solve failed", "LanczosSolver::solve");

	throw SolveFailed ();
}

template <class Field, class Vector>
bool LanczosSolver<Field, Vector>::iterate (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b) 
{
	commentator.start ("Lanczos iteration", "LanczosSolver::iterate", A.coldim ());

	// j is really a flip-flop: 0 means "even" and 1 means "odd". So "j" and
	// "j-2" are accessed with [j], while "j-1" and "j+1" are accessed via
	// [1-j]

	unsigned int j = 1, prods = 1;

	// N.B. For purposes of efficiency, I am purposefully making the
	// definitions of alpha and beta to be the *negatives* of what are given
	// in the Lambert thesis. This allows me to use stock vector AXPY
	// without any special modifications.

	typename Field::Element alpha, beta, delta[2], wb;

	// Zero out the vector _w[0]
	_VD.subin (_w[0], _w[0]);

	// Get a random vector _w[1]
	RandomDenseStream<Field, Vector> stream (_F, _randiter, A.coldim ());
	stream >> _w[1];

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Random w_1: ";
	_VD.write (report, _w[1]) << endl;

	A.apply (_Aw, _w[j]);                // Aw_j
	_VD.dot (delta[j], _w[j], _Aw);      // delta_j <- <w_j, Aw_j>

	if (_F.isZero (delta[j])) {
		commentator.stop ("FAILED", "<w_1, Aw_1> = 0", "LanczosSolver::iterate");
		return false;
	}

	_VD.dot (alpha, _Aw, _Aw);           //   alpha <- -<Aw_j, Aw_j> / delta_j
	_F.divin (alpha, delta[j]);
	_F.negin (alpha);

	_F.subin (beta, beta);               //    beta <- 0

	_VD.dot (wb, _w[j], b);              //       x <- <w_j, b> / delta_j w_j
	_F.divin (wb, delta[j]);
	_VD.mul (x, _w[j], wb);

	while (!_F.isZero (delta[j])) {
		commentator.progress ();

		commentator.indent (report);
		report << "Total matrix-vector products so far: " << prods << endl;

		commentator.indent (report);
		report << "alpha     = ";
		_F.write (report, alpha) << endl;

		commentator.indent (report);
		report << "beta      = ";
		_F.write (report, beta) << endl;

		commentator.indent (report);
		report << "w_j-1     = ";
		_VD.write (report, _w[1 - j]) << endl;

		commentator.indent (report);
		report << "w_j       = ";
		_VD.write (report, _w[j]) << endl;

		_VD.mulin (_w[1 - j], beta);    //   w_j+1 <- Aw_j + alpha w_j + beta w_j-1
		_VD.axpyin (_w[1 - j], alpha, _w[j]);
		_VD.addin (_w[1 - j], _Aw);

		commentator.indent (report);
		report << "w_j+1     = ";
		_VD.write (report, _w[1 - j]) << endl;

		commentator.indent (report);
		report << "Aw_j      = ";
		_VD.write (report, _Aw) << endl;

		j = 1 - j;                      //       j <- j + 1

		A.apply (_Aw, _w[j]);           // Aw_j

		_VD.dot (delta[j], _w[j], _Aw); // delta_j <- <w_j, Aw_j>

		commentator.indent (report);
		report << "delta_j-1 = ";
		_F.write (report, delta[1 - j]) << endl;

		commentator.indent (report);
		report << "delta_j   = ";
		_F.write (report, delta[j]) << endl;

		if (!_F.isZero (delta[j])) {
			_VD.dot (alpha, _Aw, _Aw);             // alpha <- -<Aw_j, Aw_j> / delta_j
			_F.divin (alpha, delta[j]);
			_F.negin (alpha);

			_F.div (beta, delta[j], delta[1 - j]); //  beta <- -delta_j / delta_j-1
			_F.negin (beta);

			_VD.dot (wb, _w[j], b);                //     x <- x + <w_j, b> / delta_j w_j
			_F.divin (wb, delta[j]);
			_VD.axpyin (x, wb, _w[j]);
		}

		++prods;
	}

	commentator.indent (report);
	report << "Total matrix-vector products: " << prods << endl;

	commentator.stop ("done", "delta_j = 0", "LanczosSolver::iterate");
	return true;
}
 
}  // namespace LinBox

#endif // __LANCZOS_INL
