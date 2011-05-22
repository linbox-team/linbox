/* linbox/algorithms/la-block-lanczos.h
 * Copyright 2002-2004 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.waterloo.ca>
 *
 * --------------------------------------------
 *
 * Licensed under the GNU Lesser General Public License. See COPYING for
 * details.
 *
 * Class definitions for block Lanczos iteration
 */

#ifndef __LINBOX_la_block_lanczos_H
#define __LINBOX_la_block_lanczos_H

#include "linbox/linbox-config.h"

#include <vector>
#include <deque>

#include "linbox/field/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/dense.h"
#include "linbox/matrix/dense-submatrix.h"
#include "linbox/solutions/methods.h"
#include "linbox/algorithms/eliminator.h"

// Fix for Solaris wierdness
#undef _N
#undef _I
#undef _C
#undef _W
#undef _P
#undef _Q

namespace LinBox 
{

/** Biorthogonalising block Lanczos iteration
 *
 * This is a biorthogonalising variant of Montgomery's block Lanczos
 * iteration. The goal is to avoid having to symmetrise the input
 * matrix by constructing two sequences of block vectors that have
 * mutual orthogonality properties. This algorithm was proposed by
 * Bradford Hovinen.
 */
template <class Field, class Matrix = DenseMatrixBase<typename Field::Element> >
class LABlockLanczosSolver
{
    public:

	typedef typename Field::Element Element;

	/** Constructor
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 */
	LABlockLanczosSolver (const Field &F,
				    const BlockLanczosTraits &traits)
		: _traits (traits), _F (F), _VD (F), _MD (F), _randiter (F),
		  _uAv (this), _eliminator (F, _traits.blockingFactor ())
		{ init_temps (); }

	/** Constructor with a random iterator
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 * @param r Random iterator to use for randomization
	 */
	LABlockLanczosSolver (const Field &F,
			      const BlockLanczosTraits &traits,
			      typename Field::RandIter r)
		: _traits (traits), _F (F), _VD (F), _MD (F), _randiter (r),
		  _uAv (this), _eliminator (F, _traits.blockingFactor ())
		{ init_temps (); }

	/** Destructor
	 */
	~LABlockLanczosSolver ();

	/** Solve the linear system Ax = b.
	 *
	 * If the system is nonsingular, this method computes the unique
	 * solution to the system Ax = b. If the system is singular, it computes
	 * a random solution.
	 *
	 * @param A Black box for the matrix A
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @return True on success; false on failure
	 */
	template <class Blackbox, class Vector>
	bool solve (const Blackbox &A, Vector &x, const Vector &b);

	/** Sample uniformly from the (right) nullspace of A
	 *
	 * @param A Black box for the matrix A
	 * @param x Matrix into whose columns to store nullspace elements
	 * @return Number of nullspace vectors found
	 */
	template <class Blackbox, class Matrix1>
	unsigned int sampleNullspace (const Blackbox &A, Matrix1 &x);

	/** Estimate the rank of A
	 *
	 * @param A Black box for the matrix A
	 * @return Lower bound on the rank of A
	 */
	template <class Blackbox>
	unsigned int rank (const Blackbox &A);

    private:

	typedef typename MatrixDomain<Field>::Permutation Permutation;

	class BasisTransformation
	{
		LABlockLanczosSolver      &_solver;

		std::vector<Permutation>   _P;
		std::vector<Matrix *>      _T;
		std::vector<unsigned int>  _rho;
		std::vector<unsigned int>  _s;

		unsigned int _N;

		template <class Matrix1>
		void applyOne (Matrix1 &M, Permutation &P, Matrix *T, unsigned int rho, unsigned int s, bool left);

	    public:
		template <class Matrix1>
		Matrix1 &apply (Matrix1 &A, bool left);

		template <class Matrix1>
		Matrix1 &applyPermutation (Matrix1 &A, bool left);

		template <class Matrix1>
		Matrix1 &applyLast (Matrix1 &A, bool left);

		template <class Matrix1>
		void append (Permutation &P, Matrix1 &T, unsigned int rho);

		void reset ();

		void reportComplete (std::ostream &out);
		void report (std::ostream &out);

		BasisTransformation (LABlockLanczosSolver<Field, Matrix> &solver, unsigned int N)
			: _solver (solver), _N (N) { reset (); }

		~BasisTransformation ();
	};

	friend class BasisTransformation;

	struct Iterate;

	// Structure representing an elimination step
	struct ElimStep
	{
		Matrix *_ujAvkmu;
		Matrix *_nuukAvj;

		unsigned int _rho;
		unsigned int _rhop;

		Iterate *_l;
		int _l_iter;
	};

	// Structure representing an iterate
	struct Iterate 
	{
		// Record of the pseudoinverse
		Matrix _udotAvbarinv;       // N x N
		Matrix _ubarAvdotinv;       // N x N

		// Record of the iterate from this iteration
		Matrix _u;                  // N x n
		Matrix _v;                  // N x n

		// Record of the dot iterate from this iteration
		Matrix _udot;               // N x n
		Matrix _vdot;               // N x n

		// Record of the basis transformation
		BasisTransformation _sigma_u;
		BasisTransformation _sigma_v;

		// Record of udot_j^TAv_j and u_j^TAvdot_j
		Matrix _udotAv;             // N x N
		Matrix _uAvdot;             // N x N

		// Record of elimination steps used on _u and _v
		std::list<ElimStep> _steps;

		int _iter;
		unsigned int _rho_u, _rho_v;
		bool _done;

		Iterate (LABlockLanczosSolver &solver, size_t n, size_t N, unsigned int iter)
			: _udotAvbarinv (N, N), _ubarAvdotinv (N, N),
			  _u (n, N), _v (n, N), _udot (n, N), _vdot (n, N),
			  _sigma_u (solver, N), _sigma_v (solver, N),
			  _udotAv (N, N), _uAvdot (N, N)
			{ init (iter); }

		void init (unsigned int iter) 
		{
			_iter = iter;
			_rho_u = _rho_v = 0;
			_done = false;
			_sigma_u.reset ();
			_sigma_v.reset ();
			_steps.clear ();
		}
	};

	// Two-dimensional array of inner products

	class InnerProductArray 
	{
		LABlockLanczosSolver *_solver;
		std::deque<std::deque<Matrix *> > _blocks;
		unsigned int _base;

	    public:
		InnerProductArray (LABlockLanczosSolver *solver) : _solver (solver), _base (0) {}

		void extend ();
		void contract ();
		Matrix *get (int i, int j);
		void reset ();
	};

	// Run the block Lanczos iteration and return the result. Return false
	// if the method breaks down. Do not check that Ax = b in the end
	template <class Blackbox>
	void iterate (const Blackbox &A);

	template <class Matrix1>
	void fixInnerProducts (typename std::list<Iterate *>::iterator l, const Matrix1 &Cu, const Matrix1 &Cv, unsigned int iter);

	template <class Blackbox>
	void tailDecomp (typename std::list<Iterate *>::iterator l, Iterate *i, const Blackbox &A);

	// Clean up the queue of iterates and return an Iterate structure to
	// store the next iterate information
	void cleanup (bool all);

	void adjust_uip1Abeta (typename std::list<Iterate *>::iterator j,
			       ElimStep &step,
			       unsigned int iter);
	void compute_uip1Abeta (typename std::list<Iterate *>::iterator j, unsigned int iter);
	void adjust_alphaAvip1 (typename std::list<Iterate *>::iterator j,
				ElimStep &step,
				unsigned int iter);
	void compute_alphaAvip1 (typename std::list<Iterate *>::iterator j, unsigned int iter);

	void augmentuidotAv (Iterate *i, Iterate *l, unsigned int rho);
	void augmentuAvidot (Iterate *i, Iterate *l, unsigned int rho);

	// Augment the inner products in udotAv with information on a new profile
	void augmentuldotAv (Iterate *l, Iterate *i, std::vector<unsigned int> &profile, unsigned int rho);

	// Augment the inner products in uAvdot with information on a new profile
	void augmentuAvldot (Iterate *l, Iterate *i, std::vector<unsigned int> &profile, unsigned int rho);

	template <class Matrix1, class Matrix2>
	void extractMinor (Matrix1 &M, Matrix2 &M1, std::vector<unsigned int> &profile);

	// Retrieve a new iterate structure
	Iterate *getNextIterate (unsigned int iter);

	// Management of the grid of inner products
	Matrix *newBlock ();

	// Initialize the temporaries used in computation
	void init_temps ();

	template <class Blackbox>
	void checkInnerProducts (const Blackbox &A);

	template <class Matrix1, class Matrix2, class Blackbox>
	void checkAConjugacy (const Matrix1  &u,
			      const Matrix2  &v,
			      const Blackbox &A,
			      size_t          u_iter,
			      size_t          v_iter,
			      size_t          rho_u,
			      size_t          rho_v);

	template <class T>
	inline const T &max (const T &a, const T &b) const
	{
		return (a > b) ? a : b;
	}

	template <class T>
	inline const T &min (const T &a, const T &b) const
	{
		return (a < b) ? a : b;
	}

	// Private variables

	const BlockLanczosTraits  _traits;
	const Field               &_F;
	VectorDomain<Field>        _VD;
	MatrixDomain<Field>        _MD;
	typename Field::RandIter   _randiter;

	// Temporaries used in the computation

	mutable Matrix    _T1;           // N x N
	mutable Matrix    _T2;           // N x N
	mutable Matrix    _T3;           // N x N
	mutable Matrix    _T4;           // N x N
	mutable Matrix    _T5;           // N x N
	mutable Matrix    _W;            // N x N

	Matrix            _ATu;          // n x N
	Matrix            _Av;           // n x N

	Matrix            _Cu;           // N x N
	Matrix            _Cv;           // N x N

	Permutation       _P;
	Permutation       _Q;

	Matrix            _v0;           // n x N
	Matrix            _b;            // n x <=N
	Matrix            _x;            // n x <=N
	Matrix            _y;            // n x <=N

	Element           _one;

	unsigned int	  _iter;
	unsigned int	  _total_dim;
	unsigned int	  _rank;

	std::vector<unsigned int> _profile;

	// Records used in managing Ucirc and Vcirc
	std::map<unsigned int, unsigned int> _gamma;

	InnerProductArray _uAv;  // Array of inner products wrt. original bases

	std::list<Iterate *>  _history;
	std::stack<Iterate *> _it_trashcan; // Unused Iterate structures go here
	std::stack<Matrix *>  _ip_trashcan; // Unused inner product matrices
					    // go here

	Eliminator<Field, Matrix> _eliminator;

	// Construct a transpose matrix on the fly
	template <class Matrix1>
	static inline TransposeMatrix<Matrix1> transpose (Matrix1 &M)
		{ return TransposeMatrix<Matrix1> (M); }
};

} // namespace LinBox

#include "linbox/algorithms/la-block-lanczos.inl"

#endif // __LINBOX_la_block_lanczos_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
