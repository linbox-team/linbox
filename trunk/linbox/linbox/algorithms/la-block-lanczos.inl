
/* linbox/algorithms/la-block-lanczos.inl
 * Copyright 2002-2004 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.waterloo.ca>
 *
 * --------------------------------------------
 *
 * Licensed under the GNU Lesser General Public License. See COPYING for
 * details.
 *
 * Function definitions for block Lanczos iteration
 */

#ifndef __LINBOX_la_block_lanczos_INL
#define __LINBOX_la_block_lanczos_INL

#include "linbox/linbox-config.h"

#include <iostream>
#include <cassert>
#include <algorithm>
#include <cmath>

#include "linbox/util/debug.h"
#include "linbox/solutions/methods.h"
#include "linbox/matrix/dense-submatrix.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"
#include "linbox/algorithms/la-block-lanczos.h"

namespace LinBox 
{

#ifdef LABL_DETAILED_TRACE

template <class Field, class Matrix>
void LABLTraceReport (std::ostream &out, MatrixDomain<Field> &MD, const char *text, size_t iter, const Matrix &M)
{
	out << text << " [" << iter << "]:" << std::endl;
	MD.write (out, M);
}

template <class Field, class Vector>
void LABLTraceReport (std::ostream &out, VectorDomain<Field> &VD, const char *text, size_t iter, const Vector &v)
{
	out << text << " [" << iter << "]: ";
	VD.write (out, v) << std::endl;
}

void LABLReportPriorityIndices (std::ostream &out, std::list<size_t> &_priority_indices, const char *which)
{
	out << "Priority indices for " << which << ": ";
	std::list<std::size_t>::iterator pi_it = _priority_indices.begin ();
	while (pi_it != _priority_indices.end ()) {
		out << *pi_it;
		if (++pi_it != _priority_indices.end ())
			out << ", ";
		else
			break;
	}
	out << std::endl;
}

#else

template <class Domain, class Object>
inline void LABLTraceReport (std::ostream &out, Domain &D, const char *text, size_t iter, const Object &obj)
{}

void LABLReportPriorityIndices (std::ostream &, std::list<size_t> &, const char *)
{}

#endif

#ifdef DETAILED_PROFILE
#  define TIMER_DECLARE(part) UserTimer part##_timer; double part##_time = 0.0;
#  define TIMER_START(part) part##_timer.start ()
#  define TIMER_STOP(part) part##_timer.stop (); part##_time += part##_timer.time ()
#  define TIMER_REPORT(part) \
	commentator.report (Commentator::LEVEL_NORMAL, TIMING_MEASURE) \
		<< "Total " #part " time: " << part##_time << "s" << std::endl;
#else
#  define TIMER_DECLARE(part)
#  define TIMER_START(part)
#  define TIMER_STOP(part)
#  define TIMER_REPORT(part)
#endif



template <class Field, class Matrix>
LABlockLanczosSolver<Field, Matrix>::~LABlockLanczosSolver ()
{
	while (!_it_trashcan.empty ()) {
		delete _it_trashcan.top ();
		_it_trashcan.pop ();
	}

	while (!_ip_trashcan.empty ()) {
		delete _ip_trashcan.top ();
		_ip_trashcan.pop ();
	}
}

// N.B. This code was lifted from the Lanczos iteration in LinBox

template <class Field, class Matrix>
template <class Blackbox, class Vector>
bool LABlockLanczosSolver<Field, Matrix>::solve
	(const Blackbox &A, Vector &x, const Vector &b) 
{
	linbox_check (x.size () == A.coldim ());
	linbox_check (b.size () == A.rowdim ());
	linbox_check (A.rowdim () == A.coldim ());

	commentator.start ("Solving linear system (Biorthogonalising block Lanczos)",
			   "LABlockLanczosSolver::solve");

	std::ostream &reportU = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	Vector Ax;
	bool success = false;

	// Get the temporaries into the right sizes
	_b.resize (b.size (), 1);
	_x.resize (x.size (), 1);

	_v0.resize (A.coldim (), _traits.blockingFactor ());
	_ATu.resize (A.rowdim (), _traits.blockingFactor ());
	_Av.resize (A.coldim (), _traits.blockingFactor ());

	// Copy the right-hand side to the first column of _b
	_VD.copy (*(_b.colBegin ()), b);
	_VD.copy (*(_v0.colBegin ()), b);

	// Fill the remaining columns of _v0 with random data
	RandomDenseStream<Field, typename Matrix::Col> stream (_F, _randiter, A.coldim ());
	typename Matrix::ColIterator iter = _v0.colBegin ();

	for (++iter; iter != _v0.colEnd (); ++iter)
		stream >> *iter;

	// Run the iteration up to maximum number of tries
	for (unsigned int i = 0; !success && i < _traits.maxTries (); ++i) {
		iterate (A);
		_VD.copy (x, *(_x.colBegin ()));

		if (_traits.checkResult ()) {
			VectorWrapper::ensureDim (Ax, A.rowdim ());

			if (_traits.checkResult ()) {
				commentator.start ("Checking whether Ax=b");

				A.apply (Ax, x);

				if (_VD.areEqual (Ax, b)) {
					commentator.stop ("passed");
					success = true;
				} else {
					commentator.stop ("FAILED");
					success = false;

					reportU << "Ax = ";
					_VD.write (reportU, x) << std::endl;

					reportU << " b = ";
					_VD.write (reportU, b) << std::endl;
				}
			}
		}
	}

	// Clear out the collection of iterates
	while (!_it_trashcan.empty ()) {
		delete _it_trashcan.top ();
		_it_trashcan.pop ();
	}

	while (!_history.empty ()) {
		delete _history.front ();
		_history.pop_front ();
	}

	_uAv.reset ();

	commentator.stop ("done", (success ? "Solve successful" : "Solve failed"), "LABlockLanczosSolver::solve");

	return success;
}

template <class Field, class Matrix>
template <class Blackbox, class Matrix1>
unsigned int LABlockLanczosSolver<Field, Matrix>::sampleNullspace
	(const Blackbox &A, Matrix1 &x) 
{
	linbox_check (x.rowdim () == A.coldim ());
	linbox_check (x.coldim () == _traits.blockingFactor ());
	linbox_check (A.rowdim () == A.coldim ());

	commentator.start ("Sampling from nullspace (Lookahead-based block Lanczos)",
			   "LABlockLanczosSolver::sampleNullspace");

	// Get the temporaries into the right sizes
	_b.resize (x.rowdim (), _traits.blockingFactor ());
	_x.resize (x.rowdim (), _traits.blockingFactor ());
	_y.resize (x.rowdim (), _traits.blockingFactor ());

	_v0.resize (A.coldim (), _traits.blockingFactor ());
	_ATu.resize (A.rowdim (), _traits.blockingFactor ());
	_Av.resize (A.coldim (), _traits.blockingFactor ());

	// Fill y with random data
	RandomDenseStream<Field, typename Matrix::Col> stream (_F, _randiter, A.coldim ());
	typename Matrix::ColIterator iter;

	for (iter = _y.colBegin (); iter != _y.colEnd (); ++iter)
		stream >> *iter;

	// Fill v0 with random data
	stream.reset ();

	for (iter = _v0.colBegin (); iter != _v0.colEnd (); ++iter)
		stream >> *iter;

	// Find a right-hand side for the linear system
	_MD.blackboxMulLeft (_b, A, _y);

	// Run the iteration proper
	iterate (A);

	// Obtain the candidate nullspace vectors
	_MD.subin (_x, _y);

	// Copy vectors of _x that are true nullspace vectors into the solution
	_MD.blackboxMulLeft (_b, A, _x);

	typename Matrix::ColIterator bi, xi = x.colBegin (), xip;
	unsigned int number = 0;

	for (bi = _b.colBegin (), xip = _x.colBegin (); bi != _b.colEnd (); ++bi, ++xip) {
		if (_VD.isZero (*bi) && !_VD.isZero (*xip)) {
			_VD.copy (*xi, *xip);
			++number; ++xi;
		}
	}

	// Clear out the collection of iterates
	while (!_it_trashcan.empty ()) {
		delete _it_trashcan.top ();
		_it_trashcan.pop ();
	}

	while (!_history.empty ()) {
		delete _history.front ();
		_history.pop_front ();
	}

	_uAv.reset ();

	commentator.stop ("done", NULL, "LABlockLanczosSolver::sampleNullspace");
	return number;
}

template <class Field, class Matrix>
template <class Blackbox>
unsigned int LABlockLanczosSolver<Field, Matrix>::rank
	(const Blackbox &A) 
{
	linbox_check (A.rowdim () == A.coldim ());

	commentator.start ("Rank (Lookahead-based block Lanczos)",
			   "LABlockLanczosSolver::rank");

	// Get the temporaries into the right sizes
	_b.resize (A.rowdim (), 1);
	_x.resize (A.rowdim (), 1);
	_y.resize (A.rowdim (), _traits.blockingFactor ());

	_v0.resize (A.coldim (), _traits.blockingFactor ());
	_ATu.resize (A.rowdim (), _traits.blockingFactor ());
	_Av.resize (A.coldim (), _traits.blockingFactor ());

	// Fill v0 with random data
	RandomDenseStream<Field, typename Matrix::Col> stream (_F, _randiter, A.coldim ());
	typename Matrix::ColIterator iter;

	for (iter = _y.colBegin (); iter != _y.colEnd (); ++iter)
		stream >> *iter;

	_MD.blackboxMulLeft (_v0, A, _y);

	// Run the iteration proper
	iterate (A);

	// Clear out the collection of iterates
	while (!_it_trashcan.empty ()) {
		delete _it_trashcan.top ();
		_it_trashcan.pop ();
	}

	while (!_history.empty ()) {
		delete _history.front ();
		_history.pop_front ();
	}

	_uAv.reset ();

	commentator.stop ("done", NULL, "LABlockLanczosSolver::rank");
	return _rank;
}



template <class Field, class Matrix>
template <class Blackbox>
void LABlockLanczosSolver<Field, Matrix>::iterate (const Blackbox &A)
{
	linbox_check (_history.empty ());

	commentator.start ("Lookahead-based block Lanczos iteration",
			   "LABlockLanczosSolver::iterate", A.rowdim ());

	typename std::list<Iterate *>::iterator j;

	Iterate *iterate, *next_iterate;

	_iter = 0;
	_total_dim = 0;
	_rank = 0;

	const unsigned int N = _traits.blockingFactor ();

	unsigned int dead_iters = 0;
	Integer c;
	_F.characteristic (c);
	double logc = log (double (c));
	unsigned int max_dead_iters = (unsigned int) ceil (3 * log (double (A.rowdim ())) / (N * logc)) + 2;

	bool error = false;

	// FIXME: We need a mechanism to allow other procedures to access this
	TIMER_DECLARE(AV);
	TIMER_DECLARE(updateInnerProducts);
	TIMER_DECLARE(tailDecomp);
	TIMER_DECLARE(innerProducts);
	TIMER_DECLARE(projectionCoeff);
	TIMER_DECLARE(orthogonalization);
	TIMER_DECLARE(cleanup);
	TIMER_DECLARE(terminationCheck);

	_uAv.reset ();

	unsigned int history_total = 0, history_max = 0;

	// How many iterations between each progress update
	unsigned int progress_interval = A.rowdim () / _traits.blockingFactor () / 100;

	// Make sure there are a minimum of ten
	if (progress_interval == 0)
		progress_interval = 1;

	std::ostream &reportI = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &reportU = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

#ifdef LABL_DETAILED_TRACE
	std::ostream &reportN = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
#endif

	// Prepare the first iterate
	next_iterate = getNextIterate (0);

	// Get a random fat vectors _V[0] and _U[0]
	RandomDenseStream<Field, typename Matrix::Col> stream (_F, _randiter, A.coldim ());
	typename Matrix::ColIterator u_iter;

	for (u_iter = next_iterate->_u.colBegin (); u_iter != next_iterate->_u.colEnd (); ++u_iter)
		stream >> *u_iter;

	_MD.copy (next_iterate->_v, _v0);

#ifdef LABL_DETAILED_TRACE
	Matrix    u0 (A.rowdim (), _traits.blockingFactor ());
	Matrix    v0 (A.rowdim (), _traits.blockingFactor ());
#else  // Give those variables something just to satisfy the compiler
#  define AU0 iterate->_u
#  define AV0 iterate->_v
#endif

#ifdef LABL_DETAILED_TRACE
	_MD.copy (u0, next_iterate->_u);
	_MD.copy (v0, next_iterate->_v);
#endif

	_MD.subin (_x, _x);

	// Set up an initial inner product
	_uAv.extend ();
	_history.push_back (next_iterate);

	while (dead_iters < max_dead_iters) {
		TIMER_START(terminationCheck);
		if (_MD.isZero (next_iterate->_u))
			break;
		TIMER_STOP(terminationCheck);

		// Step 1: Obtain an iterate structure
		iterate = next_iterate;
		next_iterate = getNextIterate (_iter + 1);
		_uAv.extend ();

		LABLTraceReport (reportU, _MD, "u", _iter, iterate->_u);
		LABLTraceReport (reportU, _MD, "v", _iter, iterate->_v);

		LABLTraceReport (reportU, _MD, "x", _iter, _x);

		// Step 2: Compute A^T U, AV, and inner products
		TransposeMatrix<Matrix> uTA (_ATu);

		TIMER_START(AV);
		_MD.blackboxMulRight (uTA, transpose (iterate->_u), A);
		_MD.blackboxMulLeft (_Av, A, iterate->_v);
		TIMER_STOP(AV);

		TIMER_START(innerProducts);
		_MD.mul (*_uAv.get (_iter, _iter), transpose (iterate->_u), _Av);
		_MD.mul (*_uAv.get (_iter + 1, _iter), uTA, _Av);
		TIMER_STOP(innerProducts);

		if (_MD.isZero (*_uAv.get (_iter, _iter)) && _MD.isZero (*_uAv.get (_iter + 1, _iter)))
			++dead_iters;
		else
			dead_iters = 0;

		TIMER_START(updateInnerProducts);
		_MD.copy (*_uAv.get (_iter, _iter + 1), *_uAv.get (_iter + 1, _iter));

		for (j = _history.begin (); j != _history.end () && *j != iterate; ++j) {
			_MD.copy (*_uAv.get (_iter + 1, (*j)->_iter), *_uAv.get (_iter, (*j)->_iter + 1));
			_MD.copy (*_uAv.get ((*j)->_iter, _iter + 1), *_uAv.get ((*j)->_iter + 1, _iter));
			compute_alphaAvip1 (j, _iter);
			compute_uip1Abeta (j, _iter);
		}
		TIMER_STOP(updateInnerProducts);

		LABLTraceReport (reportU, _MD, "A^T u", _iter, _ATu);
		LABLTraceReport (reportU, _MD, "Av", _iter, _Av);

		_MD.copy (next_iterate->_u, _ATu);
		_MD.copy (next_iterate->_v, _Av);

		for (j = _history.begin (); j != _history.end (); ++j) {
			TIMER_START(tailDecomp);
			tailDecomp (j, iterate, A);
			TIMER_STOP(tailDecomp);

			// Step 6: Compute projection coefficients
			TIMER_START(projectionCoeff);
			DenseSubmatrix<Element> Cu (_Cu, 0, 0, N, (*j)->_rho_v);
			DenseSubmatrix<Element> Cv (_Cv, 0, 0, (*j)->_rho_u, N);

			DenseSubmatrix<Element> udotAvbarinv ((*j)->_udotAvbarinv, 0, 0, (*j)->_rho_v, (*j)->_rho_v);
			DenseSubmatrix<Element> ubarAvdotinv ((*j)->_ubarAvdotinv, 0, 0, (*j)->_rho_u, (*j)->_rho_u);

			DenseSubmatrix<Element> udot ((*j)->_udot, 0, 0, A.rowdim (), (*j)->_rho_v);
			DenseSubmatrix<Element> vdot ((*j)->_vdot, 0, 0, A.rowdim (), (*j)->_rho_u);

			_MD.copy (_T1, *_uAv.get (_iter + 1, (*j)->_iter));
			(*j)->_sigma_v.apply (_T1, false);

			DenseSubmatrix<Element> uip1Avbarj (_T1, 0, 0, N, (*j)->_rho_v);

			_MD.mul (Cu, uip1Avbarj, udotAvbarinv);
			_MD.negin (Cu);

#ifdef LABL_DETAILED_TRACE
			reportN << "Elimination step: C^u_(" << _iter + 1 << ", " << (*j)->_iter << "):" << std::endl;
			_MD.write (reportN, Cu);

			reportN << "Elimination step: (udot_" << (*j)->_iter << "^TAvbar_" << (*j)->_iter << ")^{-1}:" << std::endl;
			_MD.write (reportN, udotAvbarinv);
#endif

			_MD.copy (_T1, *_uAv.get ((*j)->_iter, _iter + 1));
			(*j)->_sigma_u.apply (_T1, true);

			DenseSubmatrix<Element> ubarjAvip1 (_T1, 0, 0, (*j)->_rho_u, N);

			_MD.mul (Cv, ubarAvdotinv, ubarjAvip1);
			_MD.negin (Cv);

#ifdef LABL_DETAILED_TRACE
			reportN << "Elimination step: C^v_(" << _iter + 1 << ", " << (*j)->_iter << "):" << std::endl;
			_MD.write (reportN, Cv);

			reportN << "Elimination step: (ubar_" << (*j)->_iter << "^TAvdot_" << (*j)->_iter << ")^{-1}:" << std::endl;
			_MD.write (reportN, ubarAvdotinv);
#endif

			TIMER_STOP(projectionCoeff);

			// Step 7: Eliminate blocks for iterate i + 1
			TIMER_START(orthogonalization);
			_MD.axpyin (next_iterate->_u, udot, transpose (Cu));
			_MD.axpyin (next_iterate->_v, vdot, Cv);
			TIMER_STOP(orthogonalization);

			// Step 8: Fix inner products
			TIMER_START(updateInnerProducts);
			fixInnerProducts (j, Cu, Cv, _iter);
			TIMER_STOP(updateInnerProducts);
		}

#ifdef LABL_DETAILED_TRACE
		unsigned int rho_u_0 = (_history.front ()->_iter == 0) ? _history.front ()->_rho_u : N;
		unsigned int rho_v_0 = (_history.front ()->_iter == 0) ? _history.front ()->_rho_v : N;

		checkAConjugacy (u0, next_iterate->_v, A, 0, _iter + 1, rho_u_0, N);
		checkAConjugacy (next_iterate->_u, v0, A, _iter + 1, 0, N, rho_v_0);
		checkAConjugacy (iterate->_u, next_iterate->_v, A, _iter, _iter + 1, iterate->_rho_u, N);
		checkAConjugacy (next_iterate->_u, iterate->_v, A, _iter + 1, _iter, N, iterate->_rho_v);
#endif

		// Step 9: Throw away unneeded iterates and update solution
		TIMER_START(cleanup);
		cleanup (false);
		TIMER_STOP(cleanup);

		// Step 10: Prepare for the next iteration
		history_total += _history.size ();
		history_max = max ((size_t) history_max, _history.size ());
		_total_dim += iterate->_rho_u;
		_history.push_back (next_iterate);
		++_iter;

		checkInnerProducts (A);

		if (!(_iter % progress_interval))
			commentator.progress (_total_dim);

		if (_total_dim > A.rowdim ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Maximum number of iterations passed without termination" << std::endl;
			error = true;
			break;
		}
	}

	// Now just take everything off the history and throw it in the trash
	TIMER_START(cleanup);
	cleanup (true);
	TIMER_STOP(cleanup);

	LABLTraceReport (reportU, _MD, "x", _iter, _x);

	TIMER_REPORT(AV);
	TIMER_REPORT(updateInnerProducts);
	TIMER_REPORT(tailDecomp);
	TIMER_REPORT(innerProducts);
	TIMER_REPORT(projectionCoeff);
	TIMER_REPORT(orthogonalization);
	TIMER_REPORT(cleanup);
	TIMER_REPORT(terminationCheck);

	reportI << "Maximum history length: " << history_max << std::endl
		<< "Average history length: " << (double) history_total / (double) _iter << std::endl;

	reportI << "Number of iterations: " << _iter << std::endl;

	if (error)
		commentator.stop ("ERROR", NULL, "LABlockLanczosSolver::iterate");
	else if (dead_iters >= max_dead_iters)
		commentator.stop ("breakdown", NULL, "LABlockLanczosSolver::iterate");
	else
		commentator.stop ("done", NULL, "LABlockLanczosSolver::iterate");
}

template <class Field, class Matrix>
template <class Matrix1>
void LABlockLanczosSolver<Field, Matrix>::fixInnerProducts
	(typename std::list<Iterate *>::iterator  l,
	 const Matrix1                           &Cu,
	 const Matrix1                           &Cv,
	 unsigned int                             iter)
{
	const unsigned int N = _traits.blockingFactor ();

	DenseSubmatrix<Element> udotAv ((*l)->_udotAv, 0, 0, Cu.coldim (), N);
	DenseSubmatrix<Element> uAvdot ((*l)->_uAvdot, 0, 0, N, Cv.rowdim ());
	_MD.axpyin (*_uAv.get (iter + 1, (*l)->_iter), Cu, udotAv);
	_MD.axpyin (*_uAv.get ((*l)->_iter, iter + 1), uAvdot, Cv);

#ifdef LABL_DETAILED_TRACE
	std::ostream &reportN = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	reportN << "fixInnerProducts: udot_" << (*l)->_iter << "^TAv_" << (*l)->_iter << ":" << std::endl;
	_MD.write (reportN, udotAv);

	reportN << "fixInnerProducts: u_" << (*l)->_iter << "^TAvdot_" << (*l)->_iter << ":" << std::endl;
	_MD.write (reportN, uAvdot);
#endif
}



template <class Field, class Matrix>
template <class Blackbox>
void LABlockLanczosSolver<Field, Matrix>::tailDecomp
	(typename std::list<Iterate *>::iterator l, Iterate *i, const Blackbox &A)
{
#ifdef LABL_DETAILED_TRACE
	commentator.start ("Tail decomposition", "LABlockLanczosSolver::tailDecomp");
	std::ostream &reportI = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::ostream &reportN = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	std::ostream &reportU = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
#endif

	typename std::list<Iterate *>::iterator j;

	const unsigned int N = _traits.blockingFactor ();

	unsigned int rho_u = 0, rho_v = 0;
	typename Field::Element d;

	if ((*l)->_rho_u < N) {
#ifdef LABL_DETAILED_TRACE
		reportI << "Performing tail decomposition using u_" << (*l)->_iter << "^TAv_" << i->_iter << std::endl;
#endif

		// Step 5.1: Do inversion step for the u side
		_MD.copy (_T1, *_uAv.get ((*l)->_iter, i->_iter));
		(*l)->_sigma_u.apply (_T1, true);
		i->_sigma_v.apply (_T1, false);

		DenseSubmatrix<Element> uhatAvhat (_T1, (*l)->_rho_u, i->_rho_v, N - (*l)->_rho_u, N - i->_rho_v);
		DenseSubmatrix<Element> uhatAvhatinv (_W, 0, 0, N - (*l)->_rho_u, N - (*l)->_rho_u);

#ifdef LABL_DETAILED_TRACE
		reportN << "ucheck_" << (*l)->_iter << "^TAvcheck_" << i->_iter << ":" << std::endl;
		_MD.write (reportN, _T1);

		reportN << "uhat_" << (*l)->_iter << "^TAvhat_" << i->_iter << ":" << std::endl;
		_MD.write (reportN, uhatAvhat);
#endif

		_eliminator.gaussJordan (uhatAvhatinv, _profile, _P, _T2, _Q, _T3, rho_u, d, uhatAvhat);

		DenseSubmatrix<Element> mu (uhatAvhatinv, 0, 0, rho_u, rho_u);
		DenseSubmatrix<Element> Tu (_T2, rho_u, 0, N - (*l)->_rho_u - rho_u, rho_u);
		DenseSubmatrix<Element> Tv (_T3, 0, rho_u, rho_u, N - i->_rho_v - rho_u);

		TransposeMatrix<DenseSubmatrix<Element> > TuT (Tu);

		(*l)->_sigma_u.append (_P, TuT, rho_u);
		(*l)->_sigma_u.applyLast ((*l)->_u, false);

		i->_sigma_v.append (_Q, Tv, rho_u);
		i->_sigma_v.applyLast (i->_v, false);

		DenseSubmatrix<Element> utildel ((*l)->_u, 0, (*l)->_rho_u, (*l)->_u.rowdim (), rho_u);

		for (j = l, ++j; j != _history.end (); ++j) {
#ifdef LABL_DETAILED_TRACE
			reportN << "u_" << (*j)->_iter << "^TAv_" << i->_iter << " before adjustment:" << std::endl;
			_MD.write (reportN, *_uAv.get ((*j)->_iter, i->_iter));
			reportN << "u_" << (*j)->_iter << "^TAv_" << i->_iter + 1 << " before adjustment:" << std::endl;
			_MD.write (reportN, *_uAv.get ((*j)->_iter, i->_iter + 1));
#endif

			DenseSubmatrix<Element> uhatj ((*j)->_u, 0, (*j)->_rho_u, (*j)->_u.rowdim (), N - (*j)->_rho_u);

			_MD.copy (_T2, *_uAv.get ((*j)->_iter, i->_iter));
			i->_sigma_v.apply (_T2, false);
			_MD.copy (_T5, _T2);
			(*j)->_sigma_u.apply (_T2, true);

			DenseSubmatrix<Element> ujhatAvieta (_T2, (*j)->_rho_u, i->_rho_v, N - (*j)->_rho_u, rho_u);

			if (!_MD.isZero (ujhatAvieta)) {
				DenseSubmatrix<Element> ujhatAvietamu (_T3, 0, 0, N - (*j)->_rho_u, rho_u);
				_MD.mul (ujhatAvietamu, ujhatAvieta, mu);
				_MD.negin (ujhatAvietamu);
				_MD.axpyin (uhatj, utildel, transpose (ujhatAvietamu));
			}

			_MD.copy (_T4, *_uAv.get ((*l)->_iter, i->_iter));
			(*l)->_sigma_u.apply (_T4, true);

			Matrix *ujAvietamu_block = newBlock ();
			(*j)->_steps.push_back (ElimStep ());
			(*j)->_steps.back ()._l = *l;
			(*j)->_steps.back ()._l_iter = (*l)->_iter;
			(*j)->_steps.back ()._ujAvkmu = ujAvietamu_block;
			(*j)->_steps.back ()._rho = (*l)->_rho_u;
			(*j)->_steps.back ()._rhop = rho_u;

			DenseSubmatrix<Element> ultildeAvi (_T4, (*l)->_rho_u, 0, rho_u, N);
			DenseSubmatrix<Element> ujAvieta (_T5, 0, i->_rho_v, N, rho_u);
			DenseSubmatrix<Element> ujAvietamu (*ujAvietamu_block, 0, 0, N, rho_u);
			_MD.mul (ujAvietamu, ujAvieta, mu);
			_MD.negin (ujAvietamu);
			_MD.axpyin (*_uAv.get ((*j)->_iter, i->_iter), ujAvietamu, ultildeAvi);

			adjust_alphaAvip1 (j, (*j)->_steps.back (), i->_iter);
		}

		if (*l != i) {
			DenseSubmatrix<Element> zeta ((*l)->_vdot, 0, (*l)->_rho_u, (*l)->_vdot.rowdim (), rho_u);
			DenseSubmatrix<Element> zeta_src (i->_v, 0, i->_rho_v, i->_v.rowdim (), rho_u);
			DenseSubmatrix<Element> mu_dest ((*l)->_ubarAvdotinv, (*l)->_rho_u, (*l)->_rho_u, rho_u, rho_u);

			_MD.copy (zeta, zeta_src);
			_MD.copy (mu_dest, mu);

			augmentuAvldot (*l, i, _profile, rho_u);
		}

		DenseSubmatrix<Element> utildei ((*l)->_u, 0, (*l)->_rho_u, (*l)->_u.rowdim (), rho_u);
		DenseSubmatrix<Element> utildei_dest (i->_udot, 0, i->_rho_v, i->_udot.rowdim (), rho_u);

		_MD.copy (utildei_dest, utildei);

		augmentuidotAv (i, *l, rho_u);

		DenseSubmatrix<Element> mu_dest (i->_udotAvbarinv, i->_rho_v, i->_rho_v, rho_u, rho_u);

		_MD.copy (mu_dest, mu);

#ifdef LABL_DETAILED_TRACE
		reportI << "u-rank of iterate " << (*l)->_iter << " before elimination: " << (*l)->_rho_u << std::endl;
		reportI << "v-rank of iterate " << i->_iter << " before elimination: " << i->_rho_v << std::endl;
		reportI << "Additional columns:                     " << rho_u << std::endl;

		reportU << "sigma^u_" << (*l)->_iter << ":" << std::endl;
		(*l)->_sigma_u.report (reportU);

		reportN << "Complete report on sigma^u_" << (*l)->_iter << ":" << std::endl;
		(*l)->_sigma_u.reportComplete (reportN);

		reportN << "sigma^v_" << i->_iter << ":" << std::endl;
		i->_sigma_v.report (reportN);

		reportU << "Complete report on sigma^v_" << i->_iter << ":" << std::endl;
		i->_sigma_v.reportComplete (reportU);
#endif
	}

	if ((*l)->_rho_v < N) {
#ifdef LABL_DETAILED_TRACE
		reportI << "Performing tail decomposition using u_" << i->_iter << "^TAv_" << (*l)->_iter << std::endl;
#endif

		// Step 5.2: Do inversion step for the v side
		_MD.copy (_T1, *_uAv.get (i->_iter, (*l)->_iter));
		i->_sigma_u.apply (_T1, true);
		(*l)->_sigma_v.apply (_T1, false);

		DenseSubmatrix<Element> uhatAvhat (_T1, i->_rho_u, (*l)->_rho_v, N - i->_rho_u, N - (*l)->_rho_v);

#ifdef LABL_DETAILED_TRACE
		reportN << "ucheck_" << i->_iter << "^TAvcheck_" << (*l)->_iter << ":" << std::endl;
		_MD.write (reportN, _T1);

		reportN << "uhat_" << i->_iter << "^TAvhat_" << (*l)->_iter << ":" << std::endl;
		_MD.write (reportN, uhatAvhat);
#endif

		DenseSubmatrix<Element> uhatAvhatinvT (_W, 0, 0, N - (*l)->_rho_v, N - (*l)->_rho_v);

		_eliminator.gaussJordan (uhatAvhatinvT, _profile, _P, _T2, _Q, _T3, rho_v, d, transpose (uhatAvhat));

		DenseSubmatrix<Element> nu (uhatAvhatinvT, 0, 0, rho_v, rho_v);
		DenseSubmatrix<Element> TuT (_T3, 0, rho_v, rho_v, N - i->_rho_u - rho_v);
		DenseSubmatrix<Element> TvT (_T2, rho_v, 0, N - (*l)->_rho_v - rho_v, rho_v);

		TransposeMatrix<DenseSubmatrix<Element> > Tv (TvT);

		(*l)->_sigma_v.append (_P, Tv, rho_v);
		(*l)->_sigma_v.applyLast ((*l)->_v, false);

		i->_sigma_u.append (_Q, TuT, rho_v);
		i->_sigma_u.applyLast (i->_u, false);

		DenseSubmatrix<Element> vtildel ((*l)->_v, 0, (*l)->_rho_v, (*l)->_v.rowdim (), rho_v);

		for (j = l, ++j; j != _history.end (); ++j) {
#ifdef LABL_DETAILED_TRACE
			reportN << "u_" << i->_iter << "^TAv_" << (*j)->_iter << " before adjustment:" << std::endl;
			_MD.write (reportN, *_uAv.get (i->_iter, (*j)->_iter));
			reportN << "u_" << i->_iter + 1 << "^TAv_" << (*j)->_iter << " before adjustment:" << std::endl;
			_MD.write (reportN, *_uAv.get (i->_iter + 1, (*j)->_iter));
#endif

			DenseSubmatrix<Element> vhatj ((*j)->_v, 0, (*j)->_rho_v, (*j)->_v.rowdim (), N - (*j)->_rho_v);

			_MD.copy (_T2, *_uAv.get (i->_iter, (*j)->_iter));
			i->_sigma_u.apply (_T2, true);
			_MD.copy (_T5, _T2);
			(*j)->_sigma_v.apply (_T2, false);

			DenseSubmatrix<Element> uizetaAjhat (_T2, i->_rho_u, (*j)->_rho_v, rho_v, N - (*j)->_rho_v);

			if (!_MD.isZero (uizetaAjhat)) {
				DenseSubmatrix<Element> nuTuizetaAjhat (_T3, 0, 0, rho_v, N - (*j)->_rho_v);
				_MD.mul (nuTuizetaAjhat, transpose (nu), uizetaAjhat);
				_MD.negin (nuTuizetaAjhat);
				_MD.axpyin (vhatj, vtildel, nuTuizetaAjhat);
			}

			_MD.copy (_T4, *_uAv.get (i->_iter, (*l)->_iter));
			(*l)->_sigma_v.apply (_T4, false);

			Matrix *nuukAvj_block = newBlock ();
			(*j)->_steps.push_back (ElimStep ());
			(*j)->_steps.back ()._l = *l;
			(*j)->_steps.back ()._l_iter = (*l)->_iter;
			(*j)->_steps.back ()._nuukAvj = nuukAvj_block;
			(*j)->_steps.back ()._rho = (*l)->_rho_v;
			(*j)->_steps.back ()._rhop = rho_v;

			DenseSubmatrix<Element> uiAvltilde (_T4, 0, (*l)->_rho_v, N, rho_v);
			DenseSubmatrix<Element> uizetaAvj (_T5, i->_rho_u, 0, rho_v, N);
			DenseSubmatrix<Element> nuTuizetaAvj (*nuukAvj_block, 0, 0, rho_v, N);
			_MD.mul (nuTuizetaAvj, transpose (nu), uizetaAvj);
			_MD.negin (nuTuizetaAvj);
			_MD.axpyin (*_uAv.get (i->_iter, (*j)->_iter), uiAvltilde, nuTuizetaAvj);

			adjust_uip1Abeta (j, (*j)->_steps.back (), i->_iter);
		}

		if (*l != i) {
			DenseSubmatrix<Element> eta ((*l)->_udot, 0, (*l)->_rho_v, (*l)->_udot.rowdim (), rho_v);
			DenseSubmatrix<Element> eta_src (i->_u, 0, i->_rho_u, i->_u.rowdim (), rho_v);
			DenseSubmatrix<Element> nu_dest ((*l)->_udotAvbarinv, (*l)->_rho_v, (*l)->_rho_v, rho_v, rho_v);

			_MD.copy (eta, eta_src);
			_MD.copy (nu_dest, transpose (nu));

			augmentuldotAv (*l, i, _profile, rho_v);
		}

		DenseSubmatrix<Element> vtildei ((*l)->_v, 0, (*l)->_rho_v, (*l)->_v.rowdim (), rho_v);
		DenseSubmatrix<Element> vtildei_dest (i->_vdot, 0, i->_rho_u, i->_vdot.rowdim (), rho_v);

		_MD.copy (vtildei_dest, vtildei);

		augmentuAvidot (i, *l, rho_v);

		DenseSubmatrix<Element> nu_dest (i->_ubarAvdotinv, i->_rho_u, i->_rho_u, rho_v, rho_v);

		_MD.copy (nu_dest, transpose (nu));

#ifdef LABL_DETAILED_TRACE
		reportI << "v-rank of iterate " << (*l)->_iter << " before elimination: " << (*l)->_rho_v << std::endl;
		reportI << "u-rank of iterate " << i->_iter << " before elimination: " << i->_rho_u << std::endl;
		reportI << "Additional columns:                     " << rho_v << std::endl;

		reportN << "sigma^v_" << (*l)->_iter << ":" << std::endl;
		(*l)->_sigma_v.report (reportN);

		reportU << "Complete report on sigma^v_" << (*l)->_iter << ":" << std::endl;
		(*l)->_sigma_v.reportComplete (reportU);

		reportN << "sigma^u_" << i->_iter << ":" << std::endl;
		i->_sigma_u.report (reportN);

		reportU << "Complete report on sigma^u_" << (*l)->_iter << ":" << std::endl;
		i->_sigma_u.reportComplete (reportU);
#endif
	}

	if (*l != i) {
		(*l)->_rho_u += rho_u;
		(*l)->_rho_v += rho_v;
	}

	i->_rho_v += rho_u;
	i->_rho_u += rho_v;

#ifdef LABL_DETAILED_TRACE
	commentator.stop ("done", NULL, "LABlockLanczosSolver::tailDecomp");
#endif
}



template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::cleanup (bool all)
{
	typename std::list<Iterate *>::iterator l = _history.begin (), second;

	const unsigned int N = _traits.blockingFactor ();

#ifdef LABL_DETAILED_TRACE
	int discard_count = 0;

	std::ostream &reportI = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
//	std::ostream &reportU = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
#endif

	if ((*l)->_rho_u < N || (*l)->_rho_v < N)
		return;

	second = l; ++second;

	while ((all && l != _history.end ()) ||
	       (second != _history.end () && (*second)->_rho_u == N && (*second)->_rho_v == N))
	{
		// Step 4: Update solution x
		DenseSubmatrix<Element> T1 (_T1, 0, 0, N, _b.coldim ());
		DenseSubmatrix<Element> T2 (_T2, 0, 0, N, _b.coldim ());

		_MD.mul (T1, transpose ((*l)->_u), _b);
		_MD.mul (T2, (*l)->_ubarAvdotinv, T1);
		_MD.axpyin (_x, (*l)->_vdot, T2);

		typename std::list<ElimStep>::iterator it;

		for (it = (*l)->_steps.begin (); it != (*l)->_steps.end (); ++it) {
			if (it->_ujAvkmu != NULL)
				_ip_trashcan.push (it->_ujAvkmu);
			if (it->_nuukAvj != NULL)
				_ip_trashcan.push (it->_nuukAvj);
		}

		(*l)->_steps.clear ();

		_uAv.contract ();

		_it_trashcan.push (*l);

		// Estimate of the rank
		_rank += (*l)->_rho_v;

#ifdef LABL_DETAILED_TRACE
		++discard_count;
#endif
		++l;
		if (!all) ++second;
	}

	_history.erase (_history.begin (), l);

#ifdef LABL_DETAILED_TRACE
	reportI << "Finished with cleanup: "
		<< discard_count << " iterate(s) discarded" << std::endl
		<< "History contents: ";

	for (l = _history.begin (); l != _history.end (); ++l)
		reportI << (*l)->_iter << " ";

	reportI << std::endl;
#endif	
}

template <class Field, class Matrix>
typename LABlockLanczosSolver<Field, Matrix>::Iterate *
LABlockLanczosSolver<Field, Matrix>::getNextIterate (unsigned int iter) 
{
	Iterate *ret;

	if (_it_trashcan.empty ()) {
#ifdef LABL_DETAILED_TRACE
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Allocating new iterate structure..." << std::endl;
#endif
		ret = new Iterate (*this, _x.rowdim (), _traits.blockingFactor (), iter);
		return ret;
	} else {
#ifdef LABL_DETAILED_TRACE
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Taking iterate structure from trash can..." << std::endl;
#endif
		ret = _it_trashcan.top ();
		_it_trashcan.pop ();
		ret->init (iter);
		_MD.subin (ret->_udotAvbarinv, ret->_udotAvbarinv);
		_MD.subin (ret->_ubarAvdotinv, ret->_ubarAvdotinv);
		return ret;
	}
}



template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::adjust_uip1Abeta
	(typename std::list<Iterate *>::iterator j,
	 ElimStep &step,
	 unsigned int iter)
{
	const unsigned int N = _traits.blockingFactor ();

	if (step._nuukAvj == NULL || step._l_iter < _history.front ()->_iter)
		return;

#ifdef LABL_DETAILED_TRACE
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Adjustment from iterate " << step._l->_iter << " for iterate " << (*j)->_iter << std::endl;
	report << "Block is: " << std::endl;
	_MD.write (report, *step._nuukAvj);
	report << "Start: " << step._rho << std::endl;
	report << "Length: " << step._rhop << std::endl;
#endif

	DenseSubmatrix<Element> nuukAvj (*step._nuukAvj, 0, 0, step._rhop, N);

	_MD.copy (_T1, *_uAv.get (iter + 1, step._l->_iter));
	step._l->_sigma_v.apply (_T1, false);

	DenseSubmatrix<Element> uip1Avltilde (_T1, 0, step._rho, N, step._rhop);

	_MD.axpyin (*_uAv.get (iter + 1, (*j)->_iter), uip1Avltilde, nuukAvj);
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::compute_uip1Abeta (typename std::list<Iterate *>::iterator j, unsigned int iter)
{
#ifdef LABL_DETAILED_TRACE
	commentator.start ("Applying history to u_{i+1}^TAv_j", "compute_uip1Abeta");
#endif

	typename std::list<ElimStep>::iterator l;

	for (l = (*j)->_steps.begin (); l != (*j)->_steps.end (); ++l)
		adjust_uip1Abeta (j, *l, iter);

#ifdef LABL_DETAILED_TRACE
	commentator.stop ("done", NULL, "compute_uip1Abeta");
#endif
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::adjust_alphaAvip1
	(typename std::list<Iterate *>::iterator j,
	 ElimStep &step,
	 unsigned int iter)
{
	const unsigned int N = _traits.blockingFactor ();

	if (step._ujAvkmu == NULL || step._l_iter < _history.front ()->_iter)
		return;

#ifdef LABL_DETAILED_TRACE
	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "Adjustment from iterate " << step._l->_iter << " for iterate " << (*j)->_iter << std::endl;
	report << "Block is: " << std::endl;
	_MD.write (report, *step._ujAvkmu);
	report << "Start: " << step._rho << std::endl;
	report << "Length: " << step._rhop << std::endl;
#endif

	DenseSubmatrix<Element> ujAvkmu (*step._ujAvkmu, 0, 0, N, step._rhop);

	_MD.copy (_T1, *_uAv.get (step._l->_iter, iter + 1));
	step._l->_sigma_u.apply (_T1, true);

	DenseSubmatrix<Element> ultildeAvip1 (_T1, step._rho, 0, step._rhop, N);

	_MD.axpyin (*_uAv.get ((*j)->_iter, iter + 1), ujAvkmu, ultildeAvip1);
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::compute_alphaAvip1 (typename std::list<Iterate *>::iterator j, unsigned int iter)
{
#ifdef LABL_DETAILED_TRACE
	commentator.start ("Applying history to u_j^TAv_{i+1}", "compute_alphaAvip1");
#endif

	typename std::list<ElimStep>::iterator l;

	for (l = (*j)->_steps.begin (); l != (*j)->_steps.end (); ++l)
		adjust_alphaAvip1 (j, *l, iter);

#ifdef LABL_DETAILED_TRACE
	commentator.stop ("done", NULL, "compute_alphaAvip1");
#endif
}



template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::augmentuidotAv
	(Iterate      *i,
	 Iterate      *l,
	 unsigned int  rho)
{
	const unsigned int N = _traits.blockingFactor ();

	_MD.copy (_T1, *_uAv.get (l->_iter, i->_iter));
	l->_sigma_u.apply (_T1, true);

	DenseSubmatrix<Element> utildeAv (_T1, l->_rho_u, 0, rho, N);
	DenseSubmatrix<Element> utildeAv_dest (i->_udotAv, i->_rho_v, 0, rho, N);

	_MD.copy (utildeAv_dest, utildeAv);
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::augmentuAvidot
	(Iterate      *i,
	 Iterate      *l,
	 unsigned int  rho)
{
	const unsigned int N = _traits.blockingFactor ();

	_MD.copy (_T1, *_uAv.get (i->_iter, l->_iter));
	l->_sigma_v.apply (_T1, false);

	DenseSubmatrix<Element> uAvtilde (_T1, 0, l->_rho_v, N, rho);
	DenseSubmatrix<Element> uAvtilde_dest (i->_uAvdot, 0, i->_rho_u, N, rho);

	_MD.copy (uAvtilde_dest, uAvtilde);
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::augmentuldotAv
	(Iterate                   *l,
	 Iterate                   *i,
	 std::vector<unsigned int> &profile,
	 unsigned int               rho)
{
	DenseSubmatrix<Element> zeta (l->_udotAv, l->_rho_v, 0, rho, l->_v.coldim ());
	_MD.copy (_T1, *_uAv.get (i->_iter, l->_iter));
	i->_sigma_u.apply (_T1, true);
	DenseSubmatrix<Element> zeta_src (_T1, i->_rho_u, 0, rho, _T1.coldim ());
	_MD.copy (zeta, zeta_src);
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::augmentuAvldot
	(Iterate                   *l,
	 Iterate                   *i,
	 std::vector<unsigned int> &profile,
	 unsigned int               rho)
{
	DenseSubmatrix<Element> zeta (l->_uAvdot, 0, l->_rho_u, l->_u.coldim (), rho);
	_MD.copy (_T1, *_uAv.get (l->_iter, i->_iter));
	i->_sigma_v.apply (_T1, false);
	DenseSubmatrix<Element> zeta_src (_T1, 0, i->_rho_v, _T1.coldim (), rho);
	_MD.copy (zeta, zeta_src);
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2>
void LABlockLanczosSolver<Field, Matrix>::extractMinor
	(Matrix1                   &M,
	 Matrix2                   &M1,
	 std::vector<unsigned int> &profile)
{
	std::vector<unsigned int>::const_iterator i;
	typename Matrix1::ColIterator ci = M.colBegin ();

	for (i = profile.begin (); i != profile.end (); ++i) {
		typename Matrix2::Column src = *(M1.colBegin () + *i);
		_VD.copy (*ci, src);
		++ci;
	}
}



template <class Field, class Matrix>
template <class Matrix1>
void LABlockLanczosSolver<Field, Matrix>::BasisTransformation::applyOne
	(Matrix1 &M, Permutation &P, Matrix *T, unsigned int rho, unsigned int s, bool left) 
{
	if (left) {
		DenseSubmatrix<Element> Mcheck (M, s, 0, _N - s, M.coldim ());

		_solver._MD.permuteRows (Mcheck, P.begin (), P.end ());

		DenseSubmatrix<Element> Mbar (M, s, 0, rho, M.coldim ());
		DenseSubmatrix<Element> Mhat (M, s + rho, 0, _N - s - rho, M.coldim ());

		DenseSubmatrix<Element> That (*T, _N - rho, s + rho, rho, _N - s - rho);
	
		_solver._MD.axpyin (Mhat, transpose (That), Mbar);
	} else {
		DenseSubmatrix<Element> Mcheck (M, 0, s, M.rowdim (), _N - s);

		_solver._MD.permuteColumns (Mcheck, P.begin (), P.end ());

		DenseSubmatrix<Element> Mbar (M, 0, s, M.rowdim (), rho);
		DenseSubmatrix<Element> Mhat (M, 0, s + rho, M.rowdim (), _N - s - rho);

		DenseSubmatrix<Element> That (*T, _N - rho, s + rho, rho, _N - s - rho);
	
		_solver._MD.axpyin (Mhat, Mbar, That);
	}
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &LABlockLanczosSolver<Field, Matrix>::BasisTransformation::apply (Matrix1 &M, bool left)
{
	linbox_check (M.coldim () == _N);

	typename std::vector<Permutation>::iterator Pi = _P.begin ();
	typename std::vector<Matrix *>::iterator Ti = _T.begin ();
	typename std::vector<unsigned int>::iterator rhoi = _rho.begin ();
	typename std::vector<unsigned int>::iterator si = _s.begin ();

	while (Ti != _T.end ()) {
		applyOne (M, *Pi, *Ti, *rhoi, *si, left);
		++Pi; ++Ti; ++rhoi; ++si;
	}

	return M;
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &LABlockLanczosSolver<Field, Matrix>::BasisTransformation::applyPermutation (Matrix1 &M, bool left)
{
	linbox_check (M.coldim () == _N);

	typename std::vector<unsigned int>::iterator si;
	typename std::vector<Permutation>::iterator Pi;

	if (left) {
		for (Pi = _P.begin (), si = _s.begin (); Pi != _P.end (); ++Pi, ++si) {
			DenseSubmatrix<Element> Mcheck (M, *si, 0, _N - *si, M.coldim ());
			_solver._MD.permuteRows (Mcheck, Pi->begin (), Pi->end ());
		}
	} else {
		for (Pi = _P.begin (), si = _s.begin (); Pi != _P.end (); ++Pi, ++si) {
			DenseSubmatrix<Element> Mcheck (M, 0, *si, M.rowdim (), _N - *si);
			_solver._MD.permuteColumns (Mcheck, Pi->begin (), Pi->end ());
		}
	}

	return M;
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &LABlockLanczosSolver<Field, Matrix>::BasisTransformation::applyLast (Matrix1 &M, bool left)
{
	linbox_check (M.coldim () == _N);

	applyOne (M, _P.back (), _T.back (), _rho.back (), _s.back (), left);
	return M;
}

template <class Field, class Matrix>
template <class Matrix1>
void LABlockLanczosSolver<Field, Matrix>::BasisTransformation::append
	(Permutation  &P,
	 Matrix1      &T,
	 unsigned int  rho) 
{
	linbox_check (T.rowdim () <= _N);
	linbox_check (T.coldim () <= _N);
	linbox_check (rho + T.coldim () <= _N);

	Matrix *Tnew = _solver.newBlock ();
	_solver._MD.subin (*Tnew, *Tnew);

	DenseSubmatrix<Element> Tnewhat (*Tnew, _N - T.rowdim (), _N - T.coldim (), T.rowdim (), T.coldim ());
	_solver._MD.copy (Tnewhat, T);

	_P.push_back (Permutation (P));
	_T.push_back (Tnew);
	_rho.push_back (rho);
	_s.push_back (_N - rho - T.coldim ());
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::BasisTransformation::reset ()
{
	for (typename std::vector<Matrix *>::iterator i = _T.begin (); i != _T.end (); ++i)
		_solver._ip_trashcan.push (*i);

	_P.clear ();
	_T.clear ();
	_rho.clear ();
	_s.clear ();
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::BasisTransformation::report (std::ostream &out)
{
	typename Matrix::RowIterator i;
	unsigned int idx;

	Matrix T (_N, _N);

	for (i = T.rowBegin (), idx = 0; i != T.rowEnd (); ++i, ++idx) {
		_solver._VD.subin (*i, *i);
		_solver._F.assign ((*i)[idx], _solver._one);
	}

	apply (T, false);

	_solver._MD.write (out, T);
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::BasisTransformation::reportComplete (std::ostream &out)
{
	typename std::vector<Permutation>::iterator Pi = _P.begin ();
	typename std::vector<Matrix *>::iterator Ti = _T.begin ();
	typename std::vector<unsigned int>::iterator rhoi = _rho.begin ();
	typename std::vector<unsigned int>::iterator si = _s.begin ();

	while (Pi != _P.end ()) {
		out << "Permutation: ";
		_solver._eliminator.writePermutation (out, *Pi) << std::endl;

		out << "Transform: " << std::endl;
		_solver._MD.write (out, **Ti);

		out << "rank:  " << *rhoi << std::endl;
		out << "start: " << *si << std::endl;

		++Pi; ++Ti; ++rhoi; ++si;
	}
}

template <class Field, class Matrix>
LABlockLanczosSolver<Field, Matrix>::BasisTransformation::~BasisTransformation () 
{
	typename std::vector<Matrix *>::iterator Ti;

	for (Ti = _T.begin (); Ti != _T.end (); ++Ti) {
		_solver._ip_trashcan.push (*Ti);
	}
}



template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::InnerProductArray::extend ()
{
	typename std::deque<std::deque<Matrix *> >::iterator i;

	for (i = _blocks.begin (); i != _blocks.end (); ++i)
		i->push_back (NULL);

	_blocks.push_back (std::deque<Matrix *> (_blocks.size () + 1, NULL));
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::InnerProductArray::contract ()
{
	typename std::deque<std::deque<Matrix *> >::iterator i;
	typename std::deque<Matrix *>::iterator j;

	++_base;

	for (j = _blocks.front ().begin (); j != _blocks.front ().end (); ++j)
		if (*j)
			_solver->_ip_trashcan.push (*j);

	_blocks.pop_front ();

	for (i = _blocks.begin (); i != _blocks.end (); ++i) {
		if (i->front ())
			_solver->_ip_trashcan.push (i->front ());

		i->pop_front ();
	}
}

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::InnerProductArray::reset ()
{
	_base = 0;

	while (!_blocks.empty ())
		contract ();
}

template <class Field, class Matrix>
Matrix *LABlockLanczosSolver<Field, Matrix>::InnerProductArray::get (int i, int j)
{
	linbox_check ((unsigned int) i >= _base);
	linbox_check ((unsigned int) j >= _base);
	linbox_check ((unsigned int) i < _base + _blocks.size ());
	linbox_check ((unsigned int) j < _base + _blocks.front ().size ());

	Matrix *ret = _blocks[i - _base][j - _base];

	if (ret == NULL)
		ret = _blocks[i - _base][j - _base] = _solver->newBlock ();

	linbox_check (ret != NULL);

	return ret;
}



template <class Field, class Matrix>
Matrix *LABlockLanczosSolver<Field, Matrix>::newBlock () 
{
	Matrix *ret;

	if (_ip_trashcan.empty ())
		ret = new Matrix (_traits.blockingFactor (), _traits.blockingFactor ());
	else {
		ret = _ip_trashcan.top ();
		_ip_trashcan.pop ();
	}

	linbox_check (ret != NULL);

#ifdef LABL_DETAILED_TRACE
	_MD.subin (*ret, *ret);
#endif

	return ret;
}



#ifdef LABL_DETAILED_TRACE

template <class Field, class Matrix>
template <class Matrix1, class Matrix2, class Blackbox>
void LABlockLanczosSolver<Field, Matrix>::checkAConjugacy
	(const Matrix1           &u,
	 const Matrix2           &v,
	 const Blackbox          &A,
	 size_t                   u_iter,
	 size_t                   v_iter,
	 size_t                   rho_u,
	 size_t                   rho_v)
{
	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Checking whether u_" << u_iter << " is A-conjugate to v_" << v_iter << "...";

	DenseSubmatrix<Element> T1p (_T1, 0, 0, rho_u, rho_v);

	Matrix Av (A.rowdim (), _traits.blockingFactor ());

	_MD.blackboxMulLeft (Av, A, v);
	_MD.mul (_T1, transpose (u), Av);

	if (_MD.isZero (T1p))
		report << "yes" << std::endl;
	else {
		report << "no" << std::endl;

		std::ostream &err_report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		err_report << "ERROR: u_" << u_iter << " is not A-conjugate to v_" << v_iter << std::endl;
		err_report << "Computed u_" << u_iter << "^T Av_" << v_iter << ":" << std::endl;
		_MD.write (report, _T1);
	}
}

template <class Field, class Matrix>
template <class Blackbox>
void LABlockLanczosSolver<Field, Matrix>::checkInnerProducts (const Blackbox &A)
{
	commentator.start ("Checking cached inner products", "LABlockLanczosSolver::checkInnerProducts");

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	typename std::list<Iterate *>::const_iterator i, j;

	Matrix Av (A.rowdim (), _traits.blockingFactor ());

	for (i = _history.begin (); i != _history.end (); ++i) {
		for (j = _history.begin (); j != _history.end (); ++j) {
			if (*i == _history.back () && *j == _history.back ())
				break;

			report << "Checking u_" << (*i)->_iter << "^TAv_" << (*j)->_iter << ": ";
			_MD.blackboxMulLeft (Av, A, (*j)->_v);
			_MD.mul (_T1, transpose ((*i)->_u), Av);
			_MD.copy (_T2, *_uAv.get ((*i)->_iter, (*j)->_iter));
			(*i)->_sigma_u.apply (_T2, true);
			(*j)->_sigma_v.apply (_T2, false);

			if (_MD.areEqual (_T1, _T2)) {
				report << "okay" << std::endl;
				report << "Inner product is" << std::endl;
				_MD.write (report, _T2);
			} else {
				report << "ERROR" << std::endl << "Computed: " << std::endl;
				_MD.write (report, _T1);
				report << "Cached: " << std::endl;
				_MD.write (report, _T2);
			}
		}
	}

	commentator.stop ("done", NULL, "LABlockLanczosSolver::checkInnerProducts");
}

#else // LABL_DETAILED_TRACE

template <class Field, class Matrix>
template <class Matrix1, class Matrix2, class Blackbox>
void LABlockLanczosSolver<Field, Matrix>::checkAConjugacy
	(const Matrix1           &u,
	 const Matrix2           &v,
	 const Blackbox          &A,
	 size_t                   u_iter,
	 size_t                   v_iter,
	 size_t                   rho_u,
	 size_t                   rho_v)
{
}

template <class Field, class Matrix>
template <class Blackbox>
void LABlockLanczosSolver<Field, Matrix>::checkInnerProducts (const Blackbox &A)
{
}

#endif // LABL_DETAILED_TRACE

template <class Field, class Matrix>
void LABlockLanczosSolver<Field, Matrix>::init_temps () 
{
	_T1.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_T2.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_T3.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_T4.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_T5.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_W.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_Cu.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_Cv.resize (_traits.blockingFactor (), _traits.blockingFactor ());
	_F.init (_one, 1);
}

} // namespace LinBox

#endif // __LINBOX_la_block_lanczos_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
