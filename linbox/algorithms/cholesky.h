#ifndef LB_Cholesky_H
#define LB_Cholesky_H
// linbox/algorithms/cholesky.h
/* This file is part of LinBox. See Copying for license info.  */
#include "linbox/blackbox/permutation.h" // for P
#include "linbox/blackbox/triangular-fibb.h" // for L
#include "linbox/blackbox/local-snf.h" // for SNF
#include "linbox/blackbox/symmetric-block-diagonal.h" // for D
#include "linbox/blackbox/transpose-bb.h" // for P^T, L^T
//#include "linbox/matrix/dense-matrix.h" // for S

namespace LinBox
{
	/* Factor symmetric S as PLNDL^TP^T where 
	   F has 6 factors, P,L,SNF,D,L^T,P^T
	P is a permutation, 
	L is unit lower triangular, and 
	SNF is diagonal powers of p, the Smith normal form of S.
	D is diagonal of 1x1 and 2x2 blocks.  
	Rank of D is returned. L and D (most of it) are stored in S.
    */
	template<class SymmetricMatrix>
	size_t cholesky(FIBBProduct<typename SymmetricMatrix::Field> &F, 
			SymmetricMatrix & S,
			const typename SymmetricMatrix::Ring::Element & p)
	;
#if 0 
	// specializations of cholesky not needed, I hope!
	template<class Ring>
	size_t cholesky(FIBBProduct<Ring> &F, SymmetricEmbedded<Ring> & S,
			const typename Ring::Element & p)
	;
	template<class Ring>
	size_t cholesky(FIBBProduct<Ring> &F, SymmetricVSPacked<Ring> & S,
			const typename Ring::Element & p)
	;
	template<class Ring>
	size_t cholesky(FIBBProduct<Ring> &F, SymmetricCSPacked<Ring> & S,
			const typename Ring::Element & p)
	;
#endif
// implementations
	template<class SymmetricMatrix>
	size_t cholesky(FIBBProduct<typename SymmetricMatrix::Field> &F, 
			SymmetricMatrix & S,
			const typename SymmetricMatrix::Ring::Element & p)
	{
		size_t n = S.rowdim(); // = S.coldim() as well
		typedef typename SymmetricMatrix::Field Ring;
		const Ring & R = S.field();

		Permutation<Ring> * Pp = new Permutation<Ring>(R, n);

		TriangularFIBB<SymmetricMatrix> * Lp 
			= new TriangularFIBB<SymmetricMatrix>(S,Tag::UpLo::Lower, Tag::Diag::Unit);
		
		LocalSNF<Ring> * SNFp = new LocalSNF<Ring>(R, n);

		SymmetricBlockDiagonal<SymmetricMatrix> * Dp 
			= new SymmetricBlockDiagonal<SymmetricMatrix>(S,0);

		TransposeFIBB<Ring> * LTp = new TransposeFIBB<Ring>(*Lp);
		TransposeFIBB<Ring> * PTp = new TransposeFIBB<Ring>(*Pp);

		F.incorporate(Pp, Lp, SNFp, Dp, LTp, PTp);

		return choleskyBlock(*P,*L,*SNF,*D,S,p);
	}// cholesky

	/* returns the rank: rank(D) = rank(S) = rank(SNF).
	   (rank is immediate in D.)
	*/
	template<class SymmetricMatrix>
	size_t choleskyBlock(
		Permutation<typename SymmetricMatrix::Field> &	P, 
		TriangularFIBB<SymmetricMatrix> & 				L,
		LocalSNF<typename SymmetricMatrix::Field> &		SNF,
		SymmetricBlockDiagonal<SymmetricMatrix> & 		D,
		SymmetricMatrix & 								S,
		const typename SymmetricMatrix::Field::Element & p)
	{
		if (n < SymmetricMatrix::choleskyThreshold())
			return choleskyIncr(*P,*L,*SNF,*D,S,p);

		typedef typename SymmetricMatrix::Field Ring;
		const Ring & R = S.field();
		size_t n = S.rowdim(); // = also S.coldim()
		size_t m = n/2;
		SymmetricMatrix A1(S,0,0,m,m);
		Permutation P1(P,0,0,m,m);
		Triangular L1(L,0,0,m,m);
		SNF SNF1(SNF,0,0,m,m);
		SBD D1(D,0,0,m,m);
		size_t k = choleskyBlock(P1,L1,SNF1,D1,A1,p);
		Block X(L,k,0,m-k,k);
		Block Y(S,m,0,n-m,k);
		Block B(S,m,k,n-m, m-k);
		SymmetricMatrix C(S,m,m,n-m,n-m);
//		B becomes B - YDX^T
		Block Y1(R,n-m,k); // owns
		MA.mul(Y1,Y,D);
		MA.maxpyin(B,Y1,XT); // B -= Y1 X^T
		MA.maxpyin(BT,Y1T,XT); // B -= Y1 X^T
//		C becomes C - YDY^T
		MA.maxpyin(C,Y1,Y^T); // C -= Y1 X^T
		A2(S,k,k,n-k,n-k);
		P2(P,k,k,n-k,n-k);
		L2(L,k,k,n-k,n-k);
		SNF2(SNF,k,k,n-k,n-k);
		D2(D,k,k,n-k,n-k);
		k2 = choleskyBlock(P2,L2,SNF2,D2,A2,p);

	} // choleskyBlock

	template<class SymmetricMatrix>
	size_t choleskyIncr(
		Permutation<typename SymmetricMatrix::Field> &	P, 
		TriangularFIBB<SymmetricMatrix> & 				L,
		LocalSNF<typename SymmetricMatrix::Field> &		SNF,
		SymmetricBlockDiagonal<SymmetricMatrix> & 		D,
		SymmetricMatrix & 								S,
		const typename SymmetricMatrix::Field::Element & p)
	{
		typedef typename SymmetricMatrix::Field Ring;
		const Ring & R = S.field();

	/* In each step, for 1x1 or 2x2 invertible A, 
	(A,B^T|B,C) -> (I,0|X,I)(A,0|0,C-XAX^T)(I,X^T|0,I), 
	where X = BA^{-1} replaces B and C-XAX^T replaces C.
	*/
		BlasMatrixDomain<Ring> MA(R);
		size_t k = 0; 
		size_t n = S.rowdim(); // = S.coldim() as well
		while (k < n) 
		{
			typedef DenseMatrix<Ring> MotherMatrix;
			typedef DenseSubmatrix<Ring> Block;
		//cerr << "calling swaps, k is " << k << ", S in is " << endl;
		//P.write(cerr)<<endl;
		//S.write(cerr, Tag::FileFormat::Plain) << endl;
			int r = PLDswaps(P,S,k);
		//cerr << "after swaps, r is " << r << endl;
		//P.write(cerr)<<endl;
		//S.write(cerr, Tag::FileFormat::Plain) << endl;
		  if (r < 1 or r > 2) 
		  {
			SNF.push_back(k);
			if (r == 0 or R.isZero(p)) 
			  break; 
			else // r == 3 and p is real
			  for(size_t i = k; i < n; ++i)
			    for(size_t j = k; j < n; ++j)
			  	  R.divin(S.refEntry(i,j), p);
		  } else { // r is 1 or 2.
			if (r == 2) 
				D.twoBlockIndex.push_back(k); // values later
			SymmetricMatrix A(S,k,k,r,r);
			//Block A(S,k,k,r,r);
			Block B(S,k+r,k,n-k-r,r);
			Block BT(S,k,k+r,r,n-k-r);
			SymmetricMatrix C(S,k+r,k+r,n-k-r,n-k-r);
			//Block C(S,k+r,k+r,n-k-r,n-k-r);
		//cerr << "A original " << endl;
		//A.write(cerr, Tag::FileFormat::Plain);
		//cerr << "X as B " << endl;
		//X.write(cerr, Tag::FileFormat::Plain);
			MotherMatrix Bcopybase(S,k+r,k,n-k-r,r);
			Block Bcopy(Bcopybase,0,0,n-k-r,r);

			MA.invin(A);
		//cerr << "A inverted " << endl;
		//A.write(cerr, Tag::FileFormat::Plain);
			MA.mul(Bcopy,B,A); // X = BA^{-1}
		//cerr << "Xcopysub " << endl;
		//Xcopy.write(cerr, Tag::FileFormat::Plain);
			
			B.copy(Bcopy); 
		//cerr << "B" << endl;
		//B.write(cerr, Tag::FileFormat::Plain);
			
			MA.maxpyin(C,B,BT); // C <- C - (BA^{-1})B^T
			// restore BT to X^T = A^{-1}B^T
			// [would prefer simply MA.mul(BT,A,BT)]
			MotherMatrix BTcopybase(S,k,k+r,r,n-k-r);
			Block BTcopy(BTcopybase,0,0,r,n-k-r);
			MA.mul(BTcopy,A,BT); 
			BT.copy(BTcopy); 
			// restore A
			MA.invin(A); 
			k += r;
		//cerr << "k became " << k << ", S out is " << endl;
		//S.write(cerr, Tag::FileFormat::Plain);
		  } // r is 1 or 2
		} // while k < n
		// fixme: fix offdiags
		for (size_t c = 0; c < D.twoBlockIndex.size(); ++c)
		{
			size_t k = D.twoBlockIndex[c];

			//S.getEntry(D.twoBLockValue[c],k+1,k)
			typename Ring::Element x; R.init(x);
			R.assign(D.twoBlockValue[c],S.getEntry(x,k+1,k));

			S.setEntry(k+1,k,R.zero); // L's entry
		}
		return k; // At this point, k may be n
	} // choleskyIncr

	/* return values:
	   0, it is a zero block
	   1, leading entry a unit
	   2, leading 2x2 has unit det
	   3, it is a nonzero block but with no units (multiple of p)
	   */
	template<class SymmetricMatrix>
	int PLDswaps(Permutation<typename SymmetricMatrix::Field> & P, 
				SymmetricMatrix & S,
				size_t k)
	{
		size_t n = S.rowdim();
		bool allzero = true;
		typedef typename SymmetricMatrix::Field Ring;
		const Ring & R = S.field();
		MatrixDomain<Ring> MA(R);
		typename Ring::Element x; MA.init(x);
		if (R.isUnit(S.getEntry(x,k,k))) return 1;
		size_t i,j;
		// search for diagonal unit 
		for(i=k+1; i < n; ++i) 
		{
			if (R.isUnit(S.getEntry(x,i,i))) break;
			allzero = allzero and R.isZero(x);
		}
		if (i < n)
		{
			//S.symetricTransposition(i,k)
			swaprowcol(S,i,k);
			P.swapcols(i,k);
			return 1;
		} 
		if (k < n-1 and R.isUnit(S.getEntry(x,k+1,k))) return 2;
		// search for off diagonal unit 
		for(i=k+1; i < n; ++i) 
		{ for(j=k; j < i; ++j) 
			{
			if (R.isUnit(S.getEntry(x,i,j))) break; 
			allzero = allzero and R.isZero(x);
			}
			if (j < i) break;
		}
		if (i < n)
		{
			swaprowcol(S,j,k);
			P.swapcols(j,k);
			swaprowcol(S,i,k+1);
			P.swapcols(i,k+1);
			return 2;
		} 
		return allzero ? 0 : 3; // it is a zero block or mult of p.
	} // PLDswaps

	template<SymmetricMatrix>
	void swaprowcol(SymmetricMatrix & S, size_t i, size_t j)
	{ S.swaprowcol(i,j); }

	template <>
	void swaprowcol(DenseSubmatrix<Ring> & S, size_t i, size_t j)
	{
		// todo: needs block form, at least iterator form.
		for (size_t k = 0; k < S.rowdim(); ++k) 
			swap(S.refEntry(k,i), S.refEntry(k,j));
		for (size_t k = 0; k < S.coldim(); ++k) 
			swap(S.refEntry(i,k), S.refEntry(j,k));
	}

} // namespace LinBox


#endif // LB_Cholesky_H
