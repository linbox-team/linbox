// linbox/algorithms/PLDF.h
/*
	See Copying for license info
*/
#include "linbox/matrix/permutation-matrix.h"
#include "linbox/matrix/dense-matrix.h"

namespace LinBox
{
	void swapcols(MatrixPermutation<uint32_t> & P, size_t i, size_t j);

	/* Factor S as PLDL^TP^T where 
	P is a permutation, 
	L is unit lower triangular, and 
	D is diagonal of 1x1 and 2x2 blocks.  
	Rank of D is returned. L, L^T and D are stored in A.
	P is assumed to be I initially.
    */
	template<class Field>
	size_t PLD(MatrixPermutation<uint32_t> & P, DenseSubmatrix<Field> & S)// typename Field::Element & p)
	/* In each step, for 1x1 or 2x2 A, 
	(A,B^T|B,C) -> (I,0|X,I)(A,0|0,C-XAX^T)(I,X^T|0,I), 
	where X = BA^{-1} replaces B and C-XAX^T replaces C.
	*/
	{
		BlasMatrixDomain<Field> FC(S.field());
		size_t k = 0; 
		size_t n = S.rowdim(); 
		while (k < n-1) 
		{
			typedef DenseMatrix<Field> MotherMatrix;
			typedef DenseSubmatrix<Field> Block;
		//cerr << "k is " << k << ", S in is " << endl;
		//S.write(cerr, Tag::FileFormat::Plain);
			size_t r = PLDswaps(P,S,k);
		//cerr << "r is " << r << endl;
			if (r < 1 or r > 2) break;
			Block A(S,k,k,r,r);
			Block B(S,k+r,k,n-k-r,r);
			Block BT(S,k,k+r,r,n-k-r);
			Block C(S,k+r,k+r,n-k-r,n-k-r);
		//cerr << "A original " << endl;
		//A.write(cerr, Tag::FileFormat::Plain);
		//cerr << "X as B " << endl;
		//X.write(cerr, Tag::FileFormat::Plain);
			MotherMatrix Bcopybase(S,k+r,k,n-k-r,r);
			Block Bcopy(Bcopybase,0,0,n-k-r,r);

			FC.invin(A);
		//cerr << "A inverted " << endl;
		//A.write(cerr, Tag::FileFormat::Plain);
			FC.mul(Bcopy,B,A); // X = BA^{-1}
		//cerr << "Xcopysub " << endl;
		//Xcopy.write(cerr, Tag::FileFormat::Plain);
			
			B.copy(Bcopy); 
		//cerr << "B" << endl;
		//B.write(cerr, Tag::FileFormat::Plain);
			
			FC.maxpyin(C,B,BT); // C <- C - (BA^{-1})B^T
			// restore BT to X^T = A^{-1}B^T
			// [would prefer simply FC.mul(BT,A,BT)]
			MotherMatrix BTcopybase(S,k,k+r,r,n-k-r);
			Block BTcopy(BTcopybase,0,0,r,n-k-r);
			FC.mul(BTcopy,A,BT); 
			BT.copy(BTcopy); 
			// restore A
			FC.invin(A); 
			k += r;
		//cerr << "k became " << k << ", S out is " << endl;
		//S.write(cerr, Tag::FileFormat::Plain);
		}
		typename Field::Element x; S.field().init(x);
		if (k == n-1 and not S.field().isZero(S.getEntry(x,k,k)))
			return n;
		else
			return k;
	} // PLD

	template<class Field>
	size_t PLDswaps(MatrixPermutation<uint32_t> & P, 
				DenseSubmatrix<Field> & S, 
				size_t k)
	{
		size_t n = S.rowdim();
		const Field & F = S.field();
		MatrixDomain<Field> FC(S.field());
		typename Field::Element x; FC.init(x);
		if (F.isUnit(S.getEntry(x,k,k))) return 1;
		size_t i,j;
		// search for diagonal unit 
		for(i=k+1; i < n; ++i) 
			if (F.isUnit(S.getEntry(x,i,i))) break;
		if (i < n)
		{
			//S.symetricTransposition(i,k)
			swaprows(S,i,k);
			swapcols(S,i,k);
			swapcols(P,i,k);
			return 1;
		} 
		if (F.isUnit(S.getEntry(x,k+1,k))) return 2;
		// search for off diagonal unit 
		for(i=k+1; i < n; ++i) 
		{ for(j=k; j < i; ++j) if (F.isUnit(S.getEntry(x,i,j))) break; if (j < i) break;
		}
		if (i < n)
		{
			swaprows(S,i,k+1);
			swapcols(S,i,k+1);
			swapcols(P,i,k+1);
			swaprows(S,j,k);
			swapcols(S,j,k);
			swapcols(P,j,k);
			return 2;
		} 
		return -k; // it is a zero block
	} // PLDswaps

	template<class Ring>
	void swaprows(DenseSubmatrix<Ring> & S, size_t i, size_t j)
	{
		//typename DenseSubmatrix<Ring>::Row cofr(S);
		//typename DenseSubmatrix<Ring>::RowIterator r = cofr[i].rowBegin();
		//typename DenseSubmatrix<Ring>::RowIterator s = cofr[j].rowBegin();
		for (size_t k = 0; k < S.coldim(); ++k) 
			swap(S.refEntry(i,k), S.refEntry(j,k));
		//while (r != cofr[i].rowEnd()) swap(*(cofr[i] + k), *(cofr[j]+k));
			//swap(*(cofr[i] + k), *(cofr[j]+k));
	}

	template<class Ring> 
	void swapcols(DenseSubmatrix<Ring> & S, size_t i, size_t j)
	{
		for (size_t k = 0; k < S.rowdim(); ++k) 
			swap(S.refEntry(k,i), S.refEntry(k,j));
		//typename DenseSubmatrix<Ring>::Col rofc;
		//typename DenseSubmatrix<Ring>::ColIterator r = rofc[i].colBegin();
		//typename DenseSubmatrix<Ring>::ColIterator s = rofc[j].colBegin();
		//while (r != rofc[i].colEnd()) swap(*r++, *s++);
	}

	//template<class Ring> 
	void swapcols(MatrixPermutation<uint32_t> & P, size_t i, size_t j)
	{
		P.TransposeCols(i,j);
	}

} // namespace LinBox


