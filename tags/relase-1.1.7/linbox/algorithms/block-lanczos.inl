/* block-lanczos.inl
 * Copyright (C) 2002  LinBox
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

#ifndef __LINBOX_block_lanczos_INL
#define __LINBOX_block_lanczos_INL

#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/util/debug.h"
#include "linbox/solutions/methods.h"
#include "linbox/matrix/dense-submatrix.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"

#include "block-lanczos.h"

// I'm putting everything inside the LinBox namespace so that I can drop all of
// this in to LinBox easily at a later date, without any messy porting.

// Fix for Solaris wierdness
#undef _S
#undef _M
#undef _N

namespace LinBox 
{

#ifdef DETAILED_TRACE

std::ostream &operator << (std::ostream &out, const std::vector<bool> &S) 
{
	std::vector<bool>::const_iterator i;

	for (i = S.begin (); i != S.end (); ++i) {
		out << ((*i) ? "1" : "0");
		if (i != S.end () - 1)
			out << ", ";
	}

	return out;
}

template <class Field, class Matrix>
void BLTraceReport (std::ostream &out, MatrixDomain<Field> &MD, const char *text, size_t iter, const Matrix &M)
{
	out << text << " [" << iter << "]:" << std::endl;
	MD.write (out, M);
}

template <class Field, class Vector>
void BLTraceReport (std::ostream &out, VectorDomain<Field> &VD, const char *text, size_t iter, const Vector &v)
{
	out << text << " [" << iter << "]: ";
	VD.write (out, v) << std::endl;
}

void reportS (std::ostream &out, const std::vector<bool> &S, size_t iter) 
{
	out << "S_" << iter << ": [" << S << "]" << std::endl;
}

template <class Field, class Matrix>
void checkAConjugacy (const MatrixDomain<Field> &MD, const Matrix &AV, Matrix &V, Matrix &T,
		      size_t AV_iter, size_t V_iter) 
{
	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Checking whether V_" << V_iter << " is A-conjugate to V_" << AV_iter << "...";

	MD.mul (T, TransposeMatrix<Matrix> (V), AV);

	if (MD.isZero (T))
		report << "yes" << std::endl;
	else {
		report << "no" << std::endl;

		std::ostream &err_report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR);
		err_report << "ERROR: V_" << V_iter << " is not A-conjugate to V_" << AV_iter << std::endl;
		err_report << "Computed V_" << V_iter << "^T AV_" << AV_iter << ":" << std::endl;
		MD.write (report, T);
	}
}

#else

template <class Domain, class Object>
inline void BLTraceReport (std::ostream &out, Domain &D, const char *text, size_t iter, const Object &obj)
{}

void reportS (std::ostream &out, const std::vector<bool> &S, size_t iter) 
{}

template <class Field, class Matrix>
inline void checkAConjugacy (const MatrixDomain<Field> &MD, const Matrix &AV, const Matrix &V, Matrix &T,
			     size_t AV_iter, size_t V_iter) 
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

// N.B. This code was lifted from the Lanczos iteration in LinBox

template <class Field, class Matrix>
template <class Blackbox, class Vector>
Vector &BlockLanczosSolver<Field, Matrix>::solve (const Blackbox &A, Vector &x, const Vector &b) 
{
	linbox_check ((x.size () == A.coldim ()) &&
		      (b.size () == A.rowdim ()));

	commentator.start ("Solving linear system (Block Lanczos)", "BlockLanczosSolver::solve");

	bool success = false;
	Vector d1, d2, b1, b2, bp, y, Ax, ATAx, ATb;

	// Get the temporaries into the right sizes
	_V[0].resize (A.coldim (), _N);
	_V[1].resize (A.coldim (), _N);
	_V[2].resize (A.coldim (), _N);
	_AV.resize (A.coldim (), _N);

	NonzeroRandIter<Field> real_ri (_F, _randiter);
	RandomDenseStream<Field, Vector, NonzeroRandIter<Field> > stream (_F, real_ri, A.coldim ());

	for (unsigned int i = 0; !success && i < _traits.maxTries (); ++i) {
		std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		switch (_traits.preconditioner ()) {
		    case BlockLanczosTraits::NO_PRECONDITIONER:
			success = iterate (A, x, b);
			break;

		    case BlockLanczosTraits::SYMMETRIZE:
		    {
			VectorWrapper::ensureDim (bp, A.coldim ());

			Transpose<Blackbox> AT (&A);
			Compose< Transpose< Blackbox>, Blackbox> B (&AT, &A);

			AT.apply (bp, b);

			success = iterate (B, x, bp);

			break;
		    }

		    case BlockLanczosTraits::PARTIAL_DIAGONAL:
		    {
			VectorWrapper::ensureDim (d1, A.coldim ());
			VectorWrapper::ensureDim (y, A.coldim ());

			stream >> d1;
			Diagonal<Field> D (_F, d1);
			Compose<Blackbox, Diagonal<Field> > B (&A, &D);

			report << "Random D: ";
			_VD.write (report, d1) << std::endl;

			success = iterate (B, y, b);

			D.apply (x, y);

			break;
		    }

		    case BlockLanczosTraits::PARTIAL_DIAGONAL_SYMMETRIZE:
		    {
			VectorWrapper::ensureDim (d1, A.rowdim ());
			VectorWrapper::ensureDim (b1, A.rowdim ());
			VectorWrapper::ensureDim (bp, A.coldim ());

			typedef Diagonal<Field> PC1;
			typedef Transpose<Blackbox> PC2;
			typedef Compose<PC1, Blackbox> CO1;
			typedef Compose<PC2, CO1> CO2;

			stream >> d1;
			PC1 D (_F, d1);
			PC2 AT (&A);
			CO1 B1 (&D, &A);
			CO2 B (&AT, &B1);

			report << "Random D: ";
			_VD.write (report, d1) << std::endl;

			D.apply (b1, b);
			AT.apply (bp, b1);

			success = iterate (B, x, bp);

			break;
		    }

		    case BlockLanczosTraits::FULL_DIAGONAL:
		    {
			VectorWrapper::ensureDim (d1, A.coldim ());
			VectorWrapper::ensureDim (d2, A.rowdim ());
			VectorWrapper::ensureDim (b1, A.rowdim ());
			VectorWrapper::ensureDim (b2, A.coldim ());
			VectorWrapper::ensureDim (bp, A.coldim ());
			VectorWrapper::ensureDim (y, A.coldim ());

			typedef Diagonal<Field> PC1;
			typedef Transpose<Blackbox> PC2;
			typedef Compose<Blackbox, PC1> CO1;
			typedef Compose<PC1, CO1> CO2;
			typedef Compose<PC2, CO2> CO3;
			typedef Compose<PC1, CO3> CO4;

			stream >> d1 >> d2;
			PC1 D1 (_F, d1);
			PC1 D2 (_F, d2);
			PC2 AT (&A);
			CO1 B1 (&A, &D1);
			CO2 B2 (&D2, &B1);
 			CO3 B3 (&AT, &B2);
			CO4 B (&D1, &B3);

			report << "Random D_1: ";
			_VD.write (report, d1) << std::endl;

			report << "Random D_2: ";
			_VD.write (report, d2) << std::endl;

			D2.apply (b1, b);
			AT.apply (b2, b1);
			D1.apply (bp, b2);

			success = iterate (B, y, bp);
			D1.apply (x, y);

			break;
		    }

		    default:
			throw PreconditionFailed (__FUNCTION__, __LINE__,
						  "preconditioner is NO_PRECONDITIONER, SYMMETRIZE, PARTIAL_DIAGONAL_SYMMETRIZE, "
						  "PARTIAL_DIAGONAL, or FULL_DIAGONAL");
		}


                
                    // JGD  11.07.2005
                    // I DON'T KNOW WHY IT IS WORKING BUT IT DOES WORK
                    // Without this negin the results is minus the correct solution
                    // I have seen this sentence in mg-bla which comforts my choice:
                    // "Because we set Winv to -Winv, we have -x at the end of the
                    //  iteration. So negate the result and return it"
                    // I don't know when to negate it, so I check first,
                    // if the result is not correct I try to negate !!!
                VectorWrapper::ensureDim (Ax, A.rowdim ());
                        
                if ((_traits.preconditioner () == BlockLanczosTraits::SYMMETRIZE) ||
                    (_traits.preconditioner () == BlockLanczosTraits::PARTIAL_DIAGONAL_SYMMETRIZE) ||
                    (_traits.preconditioner () == BlockLanczosTraits::FULL_DIAGONAL))
			{
                            VectorWrapper::ensureDim (ATAx, A.coldim ());
                            VectorWrapper::ensureDim (ATb, A.coldim ());
                            
                            
                            A.apply (Ax, x);
                            A.applyTranspose (ATAx, Ax);
                            A.applyTranspose (ATb, b);
                            
                            if (_VD.areEqual (ATAx, ATb)) {
                                success = true;
                            } else {
                                success = false;
                            }
			}
                else {                    
                    A.apply (Ax, x);
                    
                    if (_VD.areEqual (Ax, b)) {
                        success = true;
                    } else {
                        success = false;
                    }
                }
                if(! success) _VD.negin(x);
                    // End of JGD INFAMOUS HACK


		if (_traits.checkResult ()) {
			VectorWrapper::ensureDim (Ax, A.rowdim ());

			if (_traits.checkResult () &&
			    ((_traits.preconditioner () == BlockLanczosTraits::SYMMETRIZE) ||
			     (_traits.preconditioner () == BlockLanczosTraits::PARTIAL_DIAGONAL_SYMMETRIZE) ||
			     (_traits.preconditioner () == BlockLanczosTraits::FULL_DIAGONAL)))
			{
				VectorWrapper::ensureDim (ATAx, A.coldim ());
				VectorWrapper::ensureDim (ATb, A.coldim ());

				commentator.start ("Checking whether A^T Ax = A^T b");

				A.apply (Ax, x);
				A.applyTranspose (ATAx, Ax);
				A.applyTranspose (ATb, b);

				if (_VD.areEqual (ATAx, ATb)) {
					commentator.stop ("passed");
					success = true;
				} else {
					commentator.stop ("FAILED");
					success = false;
				}
			}
			else if (_traits.checkResult ()) {
				commentator.start ("Checking whether Ax=b");

				A.apply (Ax, x);

				if (_VD.areEqual (Ax, b)) {
					commentator.stop ("passed");
					success = true;
				} else {
					commentator.stop ("FAILED");
					success = false;
				}
			}
		}
	}

	if (success) {
		commentator.stop ("done", "Solve successful", "BlockLanczosSolver::solve");
		return x;
	} else {
		commentator.stop ("done", "Solve failed", "BlockLanczosSolver::solve");
		throw SolveFailed ();
	}
}

template <class Field, class Matrix>
template <class Blackbox, class Vector>
bool BlockLanczosSolver<Field, Matrix>::iterate (const Blackbox &A, Vector &x, const Vector &b)
{
	linbox_check (_V[0].rowdim () == A.rowdim ());
	linbox_check (_V[1].rowdim () == A.rowdim ());
	linbox_check (_V[2].rowdim () == A.rowdim ());
	linbox_check (_V[0].coldim () == _V[1].coldim ());
	linbox_check (_V[0].coldim () == _V[2].coldim ());

	commentator.start ("Block Lanczos iteration", "BlockLanczosSolver::iterate", A.rowdim ());

	size_t    Ni;
	size_t    total_dim = 0;

	Vector    tmp, tmp1, tmp2;

	bool      ret = true;

	VectorWrapper::ensureDim (tmp, _traits.blockingFactor ());
	VectorWrapper::ensureDim (tmp1, _traits.blockingFactor ());
	VectorWrapper::ensureDim (tmp2, A.rowdim ());

	// How many iterations between each progress update
	unsigned int progress_interval = A.rowdim () / _traits.blockingFactor () / 100;

	// Make sure there are a minimum of ten
	if (progress_interval == 0)
		progress_interval = 1;

	// i is the index for temporaries where we need to go back to i - 1
	// j is the index for temporaries where we need to go back to j - 2
	int i = 0, j = 2, next_j, prev_j = 1, iter = 2;
	typename Matrix::ColIterator k;

	TIMER_DECLARE(AV);
	TIMER_DECLARE(Winv);
	TIMER_DECLARE(solution);
	TIMER_DECLARE(orthogonalization)

	// Get a random fat vector _V[0]
	RandomDenseStream<Field, typename Matrix::Col> stream (_F, _randiter, A.coldim ());

	for (k = _V[0].colBegin (); k != _V[0].colEnd (); ++k)
		stream >> *k;

	TIMER_START(AV);
	_MD.blackboxMulLeft (_AV, A, _V[0]);
	TIMER_STOP(AV);

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	// Initialize S_-1 to IN
	std::fill (_S.begin (), _S.end (), true);

	// Iteration 1
	TIMER_START(Winv);
	_MD.mul (_VTAV, transpose (_V[0]), _AV);
	Ni = compute_Winv_S (_Winv[0], _S, _VTAV);
	TIMER_STOP(Winv);

	// Check for catastrophic breakdown
	if (Ni == 0) {
		commentator.stop ("breakdown", NULL, "BlockLanczosSolver::iterate");
		return false;
	}

	total_dim += Ni;

#ifdef DETAILED_TRACE
	report << "Iteration " << iter << ": N_i = " << Ni << std::endl;
	report << "Iteration " << iter << ": Total dimension is " << total_dim << std::endl;
#endif

	TIMER_START(solution);
	vectorMulTranspose (tmp, _V[0], b, _S);
	_MD.vectorMul (tmp1, _Winv[0], tmp);
	vectorMul (x, _V[0], tmp1, _S);
	TIMER_STOP(solution);

	mul_SST (_V[1], _AV, _S);

	TIMER_START(orthogonalization);
	mul (_AVTAVSST_VTAV, transpose (_AV), _AV, _S);

	BLTraceReport (report, _MD, "V", 0, _V[0]);
	BLTraceReport (report, _MD, "AV", 0, _AV);
	BLTraceReport (report, _MD, "V^T A V", 0, _VTAV);
	BLTraceReport (report, _MD, "Winv", 0, _Winv[0]);
	reportS (report, _S, 0);
	BLTraceReport (report, _VD, "x", 0, x);
	BLTraceReport (report, _MD, "AVSS^T", 0, _V[1]);
	BLTraceReport (report, _MD, "V^T A^2 V", 0, _AVTAVSST_VTAV);

	_MD.addin (_AVTAVSST_VTAV, _VTAV);
	_MD.mul (_DEF, _Winv[0], _AVTAVSST_VTAV);
	addIN (_DEF);

	_MD.axpyin (_V[1], _V[0], _DEF);
	TIMER_START(orthogonalization);

	BLTraceReport (report, _MD, "D", 1, _DEF);

	BLTraceReport (report, _MD, "V", 1, _V[1]);
	checkAConjugacy (_MD, _AV, _V[1], _DEF, 0, 1);

	if (_MD.isZero (_V[1])) {
		commentator.stop ("done", NULL, "BlockLanczosSolver::iterate");
		return true;
	}

	// Iteration 2
	TIMER_START(AV);
	_MD.blackboxMulLeft (_AV, A, _V[1]);
	TIMER_STOP(AV);

#ifdef DETAILED_TRACE
 	// DEBUG: Save a copy of AV_1 for use later
	Matrix AV1_backup (_AV.rowdim (), _AV.coldim ());
	_MD.copy (AV1_backup, _AV);
#endif

	TIMER_START(Winv);
	_MD.mul (_VTAV, transpose (_V[1]), _AV);
	Ni = compute_Winv_S (_Winv[1], _S, _VTAV);
	TIMER_STOP(Winv);

	// Check for catastrophic breakdown
	if (Ni == 0) {
		commentator.stop ("breakdown", NULL, "BlockLanczosSolver::iterate");
		return false;
	}

	total_dim += Ni;

#ifdef DETAILED_TRACE
	report << "Iteration " << iter << ": N_i = " << Ni << std::endl;
	report << "Iteration " << iter << ": Total dimension is " << total_dim << std::endl;
#endif

	TIMER_START(solution);
	vectorMulTranspose (tmp, _V[1], b, _S);
	_MD.vectorMul (tmp1, _Winv[1], tmp);
	vectorMul (tmp2, _V[1], tmp1, _S);
	_VD.addin (x, tmp2);
	TIMER_STOP(solution);

	mul_SST (_V[2], _AV, _S);

	TIMER_START(orthogonalization);
	mul (_AVTAVSST_VTAV, transpose (_AV), _AV, _S);

	BLTraceReport (report, _MD, "AV", 1, _AV);
	BLTraceReport (report, _MD, "V^T A V", 1, _VTAV);
	BLTraceReport (report, _MD, "Winv", 1, _Winv[1]);
	reportS (report, _S, 1);
	BLTraceReport (report, _VD, "x", 1, x);
	BLTraceReport (report, _MD, "V^T A^2 V", 1, _AVTAVSST_VTAV);

	_MD.addin (_AVTAVSST_VTAV, _VTAV);
	_MD.mul (_DEF, _Winv[1], _AVTAVSST_VTAV);
	addIN (_DEF);
	_MD.axpyin (_V[2], _V[1], _DEF);

	BLTraceReport (report, _MD, "D", 2, _DEF);

	mul (_DEF, _Winv[0], _VTAV, _S);
	_MD.axpyin (_V[2], _V[0], _DEF);
	TIMER_STOP(orthogonalization);

	BLTraceReport (report, _MD, "E", 2, _DEF);
	BLTraceReport (report, _MD, "V", 2, _V[2]);

	checkAConjugacy (_MD, _AV, _V[2], _DEF, 1, 2);

	// Now we're ready to begin the real iteration
	while (!_MD.isZero (_V[j])) {
		next_j = j + 1;
		if (next_j > 2) next_j = 0;

		TIMER_START(AV);
		_MD.blackboxMulLeft (_AV, A, _V[j]);
		TIMER_STOP(AV);

		// First compute F_i+1, where we use Winv_i-2; then Winv_i and
		// Winv_i-2 can share storage, and we don't need the old _VTAV
		// and _AVTAVSST_VTAV any more. After this, F_i+1 is stored in
		// _DEF

		TIMER_START(orthogonalization);
		_MD.mul (_T, _VTAV, _Winv[1 - i]);
		addIN (_T);
		_MD.mul (_DEF, _Winv[i], _T);
		_MD.mulin (_DEF, _AVTAVSST_VTAV);
		TIMER_STOP(orthogonalization);

		// Now get the next VTAV, Winv, and S_i
		TIMER_START(Winv);
		_MD.mul (_VTAV, transpose (_V[j]), _AV);
		Ni = compute_Winv_S (_Winv[i], _S, _VTAV);
		TIMER_STOP(Winv);

		// Check for catastrophic breakdown
		if (Ni == 0) {
			ret = false;
			break;
		}

		total_dim += Ni;

#ifdef DETAILED_TRACE
		report << "Iteration " << iter << ": N_i = " << Ni << std::endl;
		report << "Iteration " << iter << ": Total dimension is " << total_dim << std::endl;
#endif

		BLTraceReport (report, _MD, "AV", iter, _AV);
		BLTraceReport (report, _MD, "F", iter + 1, _DEF);
		BLTraceReport (report, _MD, "V^T AV", iter, _VTAV);
		BLTraceReport (report, _MD, "Winv", iter, _Winv[i]);
		reportS (report, _S, iter);

		// Now that we have S_i, finish off with F_i+1
		TIMER_START(orthogonalization);
		mulin (_V[next_j], _DEF, _S);
		TIMER_STOP(orthogonalization);

		// Update x
		TIMER_START(solution);
		vectorMulTranspose (tmp, _V[j], b, _S);
		_MD.vectorMul (tmp1, _Winv[i], tmp);
		vectorMul (tmp2, _V[j], tmp1, _S);
		_VD.addin (x, tmp2);
		TIMER_STOP(solution);

		BLTraceReport (report, _VD, "x", iter, x);

		// Compute the next _AVTAVSST_VTAV
		TIMER_START(orthogonalization);
		mul (_AVTAVSST_VTAV, transpose (_AV), _AV, _S);

		BLTraceReport (report, _MD, "V^T A^2 V", iter, _AVTAVSST_VTAV);

		_MD.addin (_AVTAVSST_VTAV, _VTAV);

		// Compute D and update V_i+1
		_MD.mul (_DEF, _Winv[i], _AVTAVSST_VTAV);
		addIN (_DEF);
		_MD.axpyin (_V[next_j], _V[j], _DEF);

		BLTraceReport (report, _MD, "D", iter + 1, _DEF);

		// Compute E and update V_i+1
		mul (_DEF, _Winv[1 - i], _VTAV, _S);
		_MD.axpyin (_V[next_j], _V[prev_j], _DEF);

		BLTraceReport (report, _MD, "E", iter + 1, _DEF);

		// Add AV_i S_i S_i^T
		addin (_V[next_j], _AV, _S);
		TIMER_STOP(orthogonalization);

		BLTraceReport (report, _MD, "V", iter + 1, _V[next_j]);
		checkAConjugacy (_MD, _AV, _V[next_j], _DEF, iter, iter + 1);

#ifdef DETAILED_TRACE
		checkAConjugacy (_MD, AV1_backup, _V[next_j], _DEF, 1, iter + 1);
#endif

		i = 1 - i;
		prev_j = j;
		j = next_j;
		++iter;

		if (!(iter % progress_interval))
			commentator.progress (total_dim);

		if (total_dim > A.rowdim ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Maximum number of iterations passed without termination" << std::endl;
			commentator.stop ("ERROR", NULL, "BlockLanczosSolver::iterate");
			return false;
		}
	}

	// Because we set Winv to -Winv, we have -x at the end of the
	// iteration. So negate the result and return it
	_VD.negin (x);

	BLTraceReport (report, _VD, "x", iter, x);

	TIMER_REPORT(AV);
	TIMER_REPORT(Winv);
	TIMER_REPORT(solution);
	TIMER_REPORT(orthogonalization)

	commentator.stop (ret ? "done" : "breakdown", NULL, "BlockLanczosSolver::iterate");

	return ret;
}

template <class Field, class Matrix>
int BlockLanczosSolver<Field, Matrix>::compute_Winv_S
	(Matrix        &Winv,
	 std::vector<bool>                               &S,
	 const Matrix  &T)
{
	linbox_check (S.size () == Winv.rowdim ());
	linbox_check (S.size () == Winv.coldim ());
	linbox_check (S.size () == T.rowdim ());
	linbox_check (S.size () == T.coldim ());
	linbox_check (S.size () == _M.rowdim ());
	linbox_check (S.size () * 2 == _M.coldim ());

#ifdef DETAILED_TRACE
	commentator.start ("Computing Winv and S", "BlockLanczosSolver::compute_Winv_S", S.size ());

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input T:" << std::endl;
	_MD.write (report, T);
#endif

	DenseSubmatrix<Element> M1 (_M, 0, 0, T.rowdim (), T.coldim ());
	DenseSubmatrix<Element> M2 (_M, 0, T.coldim (), T.rowdim (), T.coldim ());

	_MD.copy (M1, T);
	setIN (M2);

	permute (_indices, S);

	typename Field::Element Mjj_inv;

	size_t row;
	size_t Ni = 0;

	for (row = 0; row < S.size (); ++row) {
#ifdef DETAILED_TRACE
		if (!(row & ((1 << 10) - 1)))
			commentator.progress (row);

		std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Iteration " << row << ": Matrix M = " << std::endl;
		_MD.write (report, _M);
#endif

		if (find_pivot_row (_M, row, 0, _indices)) {
#ifdef DETAILED_TRACE
			commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
				<< "Pivot found for column " << _indices[row] << std::endl;
#endif

			// Pivot element was found for (j, j)

			S[_indices[row]] = true;  // Use column j of V_i in W_i

			// Give the (j, j) entry unity
			_F.inv (Mjj_inv, _M.getEntry (_indices[row], _indices[row]));
			_VD.mulin (*(_M.rowBegin () + _indices[row]), Mjj_inv);

			// Zero the rest of the column j
			eliminate_col (_M, row, 0, _indices, Mjj_inv);

			++Ni;
		} else {
#ifdef DETAILED_TRACE
			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "No pivot found for column " << _indices[row] << std::endl;
#endif

			// No pivot element found

			S[_indices[row]] = false;  // Skip column j

			find_pivot_row (_M, row, _N, _indices);

			const typename Field::Element &Mjj = _M.refEntry (_indices[row], _indices[row] + _N);

			linbox_check (!_F.isZero (Mjj));

			// Zero the rest of the column j + N
			eliminate_col (_M, row, _N, _indices, _F.inv (Mjj_inv, Mjj));

			// Zero row j
			_VD.subin (*(_M.rowBegin () + _indices[row]), *(_M.rowBegin () + _indices[row]));
		}
	}

	_MD.neg (Winv, M2);

#ifdef DETAILED_TRACE
	report << "Computed Winv:" << std::endl;
	_MD.write (report, Winv);

	commentator.stop ("done", NULL, "BlockLanczosSolver::compute_Winv_S");
#endif

	return Ni;
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2>
Matrix1 &BlockLanczosSolver<Field, Matrix>::mul_SST
	(Matrix1                 &BSST,
	 const Matrix2           &B,
	 const std::vector<bool> &S) const
{
	linbox_check (B.rowdim () == BSST.rowdim ());
	linbox_check (B.coldim () == BSST.coldim ());
	linbox_check (B.coldim () == S.size ());

	typename Matrix2::ConstColIterator i;
	typename Matrix1::ColIterator j;
	std::vector<bool>::const_iterator k;

	for (i = B.colBegin (), j = BSST.colBegin (), k = S.begin ();
	     i != B.colEnd ();
	     ++i, ++j, ++k)
	{
		if (*k)
			_VD.copy (*j, *i);
		else
			_VD.subin (*j, *j);
	}

	return BSST;
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2, class Matrix3>
Matrix1 &BlockLanczosSolver<Field, Matrix>::mul
	(Matrix1                 &C,
	 const Matrix2           &A,
	 const Matrix3           &B,
	 const std::vector<bool> &S) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (A.rowdim () == C.rowdim ());
	linbox_check (B.coldim () == C.coldim ());

	typename Matrix2::ConstRowIterator i;
	typename Matrix3::ConstColIterator j;
	typename Matrix1::ColIterator k1;
	typename Matrix1::Col::iterator k2;
	std::vector<bool>::const_iterator l;

	for (j = B.colBegin (), l = S.begin (), k1 = C.colBegin ();
	     j != B.colEnd ();
	     ++j, ++l, ++k1)
	{
		if (*l) {
			for (i = A.rowBegin (), k2 = k1->begin (); i != A.rowEnd (); ++i, ++k2)
				_VD.dot (*k2, *i, *j);
		} else
			_VD.subin (*k1, *k1);
	}

	return C;
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2>
Matrix1 &BlockLanczosSolver<Field, Matrix>::mulin
	(Matrix1                 &A,
	 const Matrix2           &B,
	 const std::vector<bool> &S) const
{
	linbox_check (A.coldim () == B.rowdim ());
	linbox_check (B.rowdim () == B.coldim ());
	linbox_check (A.coldim () == S.size ());

	typename Matrix1::RowIterator i;
	typename Matrix2::ConstColIterator j;
	typename Vector<Field>::Dense::iterator k;
	std::vector<bool>::const_iterator l;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i) {
		for (j = B.colBegin (), k = _tmp.begin (), l = S.begin ();
		     j != B.colEnd ();
		     ++j, ++k, ++l)
		{
			if (*l)
				_VD.dot (*k, *i, *j);
			else
				_F.subin (*k, *k);
		}

		_VD.copy (*i, _tmp);
	}

	return A;
}

template <class Field, class Matrix>
template <class Vector1, class Matrix1, class Vector2>
Vector1 &BlockLanczosSolver<Field, Matrix>::vectorMul
	(Vector1                 &w,
	 const Matrix1           &A,
	 const Vector2           &v,
	 const std::vector<bool> &S) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix1::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = v.begin ();
	typename std::vector<bool>::const_iterator k = S.begin ();

	_VD.subin (w, w);

	for (; j != v.end (); ++j, ++i, ++k)
		if (*k)
			_VD.axpyin (w, *j, *i);

	return w;
}

template <class Field, class Matrix>
template <class Vector1, class Matrix1, class Vector2>
Vector1 &BlockLanczosSolver<Field, Matrix>::vectorMulTranspose
	(Vector1                 &w,
	 const Matrix1           &A,
	 const Vector2           &v,
	 const std::vector<bool> &S) const
{
	linbox_check (A.rowdim () == v.size ());
	linbox_check (A.coldim () == w.size ());

	typename Matrix1::ConstColIterator i = A.colBegin ();
	typename Vector1::iterator j = w.begin ();
	typename std::vector<bool>::const_iterator k = S.begin ();

	for (; j != w.end (); ++j, ++i, ++k)
		if (*k)
			_VD.dot (*j, *i, v);

	return w;
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &BlockLanczosSolver<Field, Matrix>::addIN (Matrix1 &A) const
{
	linbox_check (A.coldim () == A.rowdim ());

	typename Matrix1::RowIterator i;
	size_t idx = 0;

	for (i = A.rowBegin (); i != A.rowEnd (); ++i, ++idx)
		_F.addin ((*i)[idx], _one);

	return A;
}

template <class Field, class Matrix>
template <class Matrix1, class Matrix2>
Matrix1 &BlockLanczosSolver<Field, Matrix>::addin
	(Matrix1                 &A,
	 const Matrix2           &B,
	 const std::vector<bool> &S) const
{
	linbox_check (A.rowdim () == B.rowdim ());
	linbox_check (A.coldim () == B.coldim ());

	typename Matrix1::ColIterator i = A.colBegin ();
	typename Matrix2::ConstColIterator j = B.colBegin ();
	std::vector<bool>::const_iterator k = S.begin ();

	for (; i != A.colEnd (); ++i, ++j, ++k)
		if (*k) _VD.addin (*i, *j);

	return A;
}

template <class Field, class Matrix>
void BlockLanczosSolver<Field, Matrix>::permute (std::vector<size_t>     &indices,
						 const std::vector<bool> &S) const
{
	size_t idx;

	std::vector<size_t>::iterator i = indices.begin ();
	std::vector<bool>::const_iterator k;

	for (k = S.begin (), idx = 0; k != S.end (); ++k, ++idx) {
		if (!*k) {
			*i = idx;
			++i;
		}
	}

	for (k = S.begin (), idx = 0; k != S.end (); ++k, ++idx) {
		if (*k) {
			*i = idx;
			++i;
		}
	}
}

template <class Field, class Matrix>
template <class Matrix1>
Matrix1 &BlockLanczosSolver<Field, Matrix>::setIN (Matrix1 &A) const 
{
	linbox_check (A.coldim () == A.rowdim ());

	typename Matrix1::RowIterator i;
	size_t i_idx;

	for (i = A.rowBegin (), i_idx = 0; i != A.rowEnd (); ++i, ++i_idx) {
		_VD.subin (*i, *i);
		_F.assign ((*i)[i_idx], _one);
	}

	return A;
}

/* Find a row suitable for pivoting in column col and exchange that row with row
 * A - Matrix on which to operate
 * idx - Index of the row with which to exchange
 * Returns true if a pivot could be found; false otherwise
 */

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::find_pivot_row
	(Matrix                    &A,
	 size_t                     row,
	 int                        col_offset,
	 const std::vector<size_t> &indices)
{
	size_t idx;

	typename Matrix::Col col_vec;
	typename Matrix::Row row_vec;

	col_vec = *(_M.colBegin () + indices[row] + col_offset);
	row_vec = *(_M.rowBegin () + indices[row]);

	for (idx = row; idx < A.rowdim (); ++idx) {
		if (!_F.isZero (A.getEntry (indices[idx], indices[row] + col_offset))) {
			if (idx != row) {
				typename Matrix::Row row1 = *(A.rowBegin () + indices[idx]);
				std::swap_ranges (row_vec.begin (), row_vec.end (), row1.begin ());
			}

			return true;
		}
	}

	return false;
}

template <class Field, class Matrix>
void BlockLanczosSolver<Field, Matrix>::eliminate_col
	(Matrix                        &A,
	 size_t                         pivot,
	 int                            col_offset,
	 const std::vector<size_t>     &indices,
	 const typename Field::Element &Ajj_inv)
{
	// I'm assuming everything left of the column with the index of the pivot row is 0
	size_t row;

	typename DenseSubmatrix<Element>::Row pivot_row;
	typename Field::Element p;

	pivot_row = *(A.rowBegin () + indices[pivot]);

	for (row = 0; row < pivot; ++row) {
		const typename Field::Element &Aij = A.getEntry (indices[row], indices[pivot] + col_offset);

		if (!_F.isZero (Aij))
			_VD.axpyin (*(A.rowBegin () + indices[row]), _F.neg (p, Aij), pivot_row);
	}

	for (++row; row < A.rowdim (); ++row) {
		const typename Field::Element &Aij = A.getEntry (indices[row], indices[pivot] + col_offset);

		if (!_F.isZero (Aij))
			_VD.axpyin (*(A.rowBegin () + indices[row]), _F.neg (p, Aij), pivot_row);
	}
}

template <class Field, class Matrix>
void BlockLanczosSolver<Field, Matrix>::init_temps () 
{
	_VTAV.resize (_N, _N);
	_Winv[0].resize (_N, _N);
	_Winv[1].resize (_N, _N);
	_AVTAVSST_VTAV.resize (_N, _N);
	_T.resize (_N, _N);
	_DEF.resize (_N, _N);
	_S.resize (_N);
	_M.resize (_N, 2 * _N);
	_tmp.resize (_N);
	_indices.resize (_N);
}

// Check whether the given matrix is "almost" the identity, i.e. the identity
// with some zeros on the diagonal

template <class Field, class Matrix>
template <class Matrix1>
bool BlockLanczosSolver<Field, Matrix>::isAlmostIdentity (const Matrix1 &M) const 
{
	linbox_check (M.rowdim () == M.coldim ());

	typename Field::Element neg_one;

	_F.init (neg_one, -1);

	size_t i, j;

	for (i = 0; i < M.rowdim (); ++i) {
		for (j = 0; j < M.coldim (); ++j) {
			if (i != j && !_F.isZero (M.getEntry (i, j))) {
				if (!_F.isZero (M.getEntry (i, i))) {
					typename Matrix::ConstRowIterator row = M.rowBegin () + j;
					if (!_VD.isZero (*row))
						return false;
				}
				else if (!_F.isZero (M.getEntry (j, j))) {
					typename Matrix::ConstColIterator col = M.colBegin () + i;
					if (!_VD.isZero (*col))
						return false;
				} else
					return false;
			}
			else if (!_F.isZero (M.getEntry (i, j)) && !_F.areEqual (M.getEntry (i, j), neg_one))
				return false;
		}
	}

	return true;
}

// Test suite for BlockLanczosSolver
// All tests below return true on success and false on failure. They take a
// single argument: n for the row and column dimension of the matrices on which
// to operate

// Compute a random dense matrix and run compute_Winv_S on it. Check that the
// resulting S vector is all 'true' and then multiply the original matrix by the
// output. Add 1 to the result with addIN and check the the result is zero with
// isZero

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_compute_Winv_S_mul (int n) const
{
	commentator.start ("Testing compute_Winv_S, mul, addIN, and isZero", "test_compute_Winv_S_mul");

	Matrix A (n, n);
	Matrix AT (n, n);
	Matrix ATA (n, n);
	Matrix W (n, n);
	Matrix WA (n, n);
	std::vector<bool> S (n);

	bool ret = true;

	RandomDenseStream<Field, typename Matrix::Row> stream (_F, _randiter, n);
	typename Matrix::RowIterator i = A.rowBegin ();
	typename Matrix::ColIterator j = AT.colBegin ();

	// With very, very, very, very high probability, this will be
	// nonsingular
	for (; i != A.rowEnd (); ++i, ++j) {
		stream >> *i;
		_VD.copy (*j, *i);
	}

	_MD.mul (ATA, AT, A);

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Computed A^T A:" << std::endl;
	_MD.write (report, ATA);

	compute_Winv_S (W, S, ATA);

	report << "Computed W:" << std::endl;
	_MD.write (report, W);

	// Now W should be -A^-1
	_MD.mul (WA, W, ATA);

	report << "Computed WA^T A:" << std::endl;
	_MD.write (report, WA);

	if (!isAlmostIdentity (WA)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: WA^T A != I" << std::endl;
		ret = false;
	}

	// Now, just for kicks, do the same on the other side

	_MD.mul (WA, ATA, W);

	report << "Computed A^T A W:" << std::endl;
	_MD.write (report, WA);

	if (!isAlmostIdentity (WA)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A^T AW != I" << std::endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), NULL, "test_compute_Winv_S_mul");

	return ret;
}

// Same as above, but use mulin rather than mul

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_compute_Winv_S_mulin (int n) const
{
	commentator.start ("Testing compute_Winv_S, copy, mulin, addIN, and isZero", "test_compute_Winv_S_mulin");

	Matrix A (n, n);
	Matrix AT (n, n);
	Matrix ATA (n, n);
	Matrix W (n, n);
	Matrix WA (n, n);
	std::vector<bool> S (n);

	bool ret = true;

	RandomDenseStream<Field, typename Matrix::Row> stream (_F, _randiter, n);
	typename Matrix::RowIterator i = A.rowBegin ();
	typename Matrix::ColIterator j = AT.colBegin ();

	// With very, very, very, very high probability, this will be
	// nonsingular
	for (; i != A.rowEnd (); ++i, ++j) {
		stream >> *i;
		_VD.copy (*j, *i);
	}

	_MD.mul (ATA, AT, A);

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Computed A^T A:" << std::endl;
	_MD.write (report, ATA);

	compute_Winv_S (W, S, ATA);

	_MD.copy (WA, W);

	report << "Computed W:" << std::endl;
	_MD.write (report, W);

	// Now W should be -A^-1
	_MD.mulin (WA, ATA);

	report << "Computed WA^T A:" << std::endl;
	_MD.write (report, WA);

	if (!isAlmostIdentity (WA)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: WA^T A != I" << std::endl;
		ret = false;
	}

	// Now, just for kicks, do the same on the other side

	_MD.copy (WA, ATA);

	_MD.mulin (WA, W);

	report << "Computed A^T AW:" << std::endl;
	_MD.write (report, WA);

	if (!isAlmostIdentity (WA)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A^T AW != I" << std::endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), NULL, "test_compute_Winv_S_mulin");

	return ret;
}

// Compute a random nonsingular diagonal matrix and set an S vector so that
// every other entry is true and the rest are false. Call mul_SST and check that
// the entries on the resulting diagonal are correct.

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_mul_SST (int n) const
{
	commentator.start ("Testing addin", "test_mulTranspose");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), NULL, "test_mulTranspose");

	return ret;
}

// Same as test_compute_Winv_S_mul above, but zero out every other column using
// the method for test_mul_SST

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_mul_ABSST (int n) const
{
	commentator.start ("Testing addin", "test_mulTranspose");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), NULL, "test_mulTranspose");

	return ret;
}

// Compute a random dense matrix and two random vectors, and check that <A^T x,
// y> = <x, Ay>

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_mulTranspose (int m, int n) const
{
	commentator.start ("Testing mulTranspose, m-v mul", "test_mulTranspose");

	Matrix A (m, n);
	typename Vector<Field>::Dense x (m), y (n);
	typename Vector<Field>::Dense ATx (n), Ay (m);
	Element ATxy, xAy;

	bool ret = true;

	RandomDenseStream<Field, typename Matrix::Row> stream (_F, _randiter, n);
	typename Matrix::RowIterator i = A.rowBegin ();

	for (; i != A.rowEnd (); ++i)
		stream >> *i;

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Computed A:" << std::endl;
	_MD.write (report, A);

	RandomDenseStream<Field, Matrix> stream1 (_F, _randiter, m);
	stream1 >> x;

	RandomDenseStream<Field, Matrix> stream2 (_F, _randiter, n);
	stream1 >> y;

	report << "Computed     x: ";
	_VD.write (report, x) << std::endl;

	report << "Computed     y: ";
	_VD.write (report, y) << std::endl;

	_MD.vectorMul (ATx, transpose (A), x);

	report << "Computed A^T x: ";
	_VD.write (report, ATx) << std::endl;

	_MD.vectorMul (Ay, A, y);

	report << "Computed    Ay: ";
	_VD.write (report, Ay) << std::endl;

	_VD.dot (ATxy, ATx, y);

	report << "Computed  ATxy: ";
	_F.write (report, ATxy) << std::endl;

	_VD.dot (xAy, x, Ay);

	report << "Computed   xAy: ";
	_F.write (report, xAy) << std::endl;

	if (!_F.areEqual (ATxy, xAy)) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: <A^T x, y> != <x, Ay>" << std::endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), NULL, "test_mulTranspose");

	return ret;
}

// Same as test_mul_SST, but using mulTranspose

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_mulTranspose_ABSST (int n) const
{
	commentator.start ("Testing addin_ABSST", "test_mulTranspose_ABSST");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), NULL, "test_mulTranspose_ABSST");

	return ret;
}

// Same as test_mul_ABSST, but using mulin_ABSST

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_mulin_ABSST (int n) const
{
	commentator.start ("Testing addin_ABSST", "test_mulin_ABSST");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), NULL, "test_mulin_ABSST");

	return ret;
}

// Same as test_addin, but using test_addin_ABSST

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::test_addin_ABSST (int n) const
{
	commentator.start ("Testing addin_ABSST", "test_addin_ABSST");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), NULL, "test_addin_ABSST");

	return ret;
}

template <class Field, class Matrix>
bool BlockLanczosSolver<Field, Matrix>::runSelfCheck () const
{
	bool ret = true;

	commentator.start ("Running self check", "runSelfCheck", 10);

	if (!test_compute_Winv_S_mul (_N)) ret = false;
	if (!test_compute_Winv_S_mulin (_N)) ret = false;
	if (!test_mul_SST (_N)) ret = false;
	if (!test_mul_ABSST (_N)) ret = false;
	if (!test_mulTranspose (_N * 10, _N)) ret = false;
	if (!test_mulTranspose_ABSST (_N)) ret = false;
	if (!test_mulin_ABSST (_N)) ret = false;
	if (!test_addin_ABSST (_N)) ret = false;

	commentator.stop (MSG_STATUS (ret), NULL, "runSelfCheck");

	return ret;
}

} // namespace LinBox

#endif // __LINBOX_block_lanczos_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
