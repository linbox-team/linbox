/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/lifting-container-base.h
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef _LIFTING_CONTAINER_H
#define _LIFTING_CONTAINER_H
//#define DEBUG_LC
//#define DEBUG_CHUNK
//#define LC_CHECK_DIVISION

#include <vector>

#include <linbox-config.h>
#include <linbox/util/debug.h>
#include <linbox/blackbox/apply.h>
#include <linbox/algorithms/blackbox-container.h>
#include <linbox/algorithms/massey-domain.h>
#include <linbox/algorithms/blackbox-block-container.h>
#include <linbox/algorithms/block-massey-domain.h>
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/compose.h>
#include <linbox/algorithms/blas-domain.h>

namespace LinBox {

	/**
	 * BoundBlackbox: Sets
	 *      H_col_sqr <- H_col(A)^2,   short_col_sqr <- short_col(A)^2
	 * where H_col(A) is prod_j sqrt(sum_i a_ij^2)     ('Hadamard column bound')
	 *   short_col(A) is min_j  sqrt(sum_i a_ij^2)     ('shortest column')
	 *
	 * note: H_col is not actually a norm! but it is what we need for lifting bound computation
	 */
	template <class Ring, class ItMatrix>
	void SpecialBound(const Ring& R, typename Ring::Element& H_col_sqr, 
			  typename Ring::Element& short_col_sqr, const ItMatrix& A) {
		typedef typename Ring::Element Integer;
		Integer sqsum;
		size_t m, n, col=0;
		n=A.coldim();
		m=A.rowdim();
		R.init(H_col_sqr, 1);
		
		typename ItMatrix::ConstColIterator colIter;
		colIter = A.colBegin();

		for (; colIter != A.colEnd(); colIter++, col++) {
			typename ItMatrix::ConstCol::const_iterator elm;
			R.init(sqsum, 0);
			for (elm = colIter->begin(); elm != colIter->end(); elm++)
				R.axpyin(sqsum, *elm, *elm);
			R.mulin(H_col_sqr, sqsum);
			if (col == 0 || sqsum < short_col_sqr)
				short_col_sqr = sqsum;
		}
	}

	// in solveNonsingular, we may work with something that inherits from DenseMatrixBase
	template <class Ring>
	void BoundBlackbox(const Ring& R, typename Ring::Element& H_col_sqr, 
			   typename Ring::Element& short_col_sqr, const DenseMatrixBase<typename Ring::Element>& A) {
		SpecialBound(R, H_col_sqr, short_col_sqr, A);
	}

	// in other solvers we generally use BlasBlackbox which inherits from DenseSubmatrix
	template <class Ring>
	void BoundBlackbox(const Ring& R, typename Ring::Element& H_col_sqr, 
			   typename Ring::Element& short_col_sqr, const DenseSubmatrix<typename Ring::Element>& A) {
		SpecialBound(R, H_col_sqr, short_col_sqr, A);
	}


	template < class Ring>
	void BoundBlackbox (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const SparseMatrix<Ring>& A) {
		typedef typename Ring::Element Integer;
		Integer one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		R.init(H_col_sqr, 1);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.apply(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter)
				sqsum += (*iter)*(*iter);
			R.mulin(H_col_sqr, sqsum);
			if (i==0 || sqsum < short_col_sqr) 
				short_col_sqr = sqsum;
			e[i]=zero;
		}
	}

	template < class Ring, class Matrix1, class Matrix2>
	void BoundBlackbox (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const Compose<Matrix1,Matrix2> & A) {
		typedef typename Ring::Element Integer;
		Integer one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		R.init(H_col_sqr, 1);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.apply(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter)
				sqsum += (*iter)*(*iter);
			R.mulin(H_col_sqr, sqsum);
			if (i==0 || sqsum < short_col_sqr) 
				short_col_sqr = sqsum;
			e[i]=zero;
		}
	}
	
		template < class Ring, class Matrix>
		void BoundBlackbox (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const Transpose<Matrix> & A) {
		typedef typename Ring::Element Integer;
		Integer one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		R.init(H_col_sqr, 1);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.applyTranspose(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter)
				sqsum += (*iter)*(*iter);
			R.mulin(H_col_sqr, sqsum);
			if (i==0 || sqsum < short_col_sqr) 
				short_col_sqr = sqsum;
			e[i]=zero;
		}
	}


	/* 
	   This should work with blackboxes. However it is much slower if column iterators are available.
	   Furthermore the compiler always binds to this instead of the above faster version; so some 
	   trickier kind of specialization may have to be done when BoundBlackBox is to be used with true blackboxes.
	   (Or is the plural Blackboxen?)

	   template < class Ring, class IMatrix>
	   void BoundBlackbox 
	   (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const IMatrix& A) {
	   typedef typename Ring::Element Integer;
	   Integer one,zero,sqsum;
	   size_t m,n;
	   n=A.coldim();
	   m=A.rowdim();
	   R.init(one,1);
	   R.init(zero,0);
	   R.init(H_col_sqr, 1);
	   typename std::vector<Integer>::const_iterator iter;
	   std::vector<Integer> e(n,zero),tmp(m);
	   for (size_t i=0;i<n;i++){
	   e[i]=one;
	   A.apply(tmp,e);
	   sqsum=zero;
	   for (iter=tmp.begin();iter!=tmp.end();++iter)
	   sqsum += (*iter)*(*iter);
	   R.mulin(H_col_sqr, sqsum);
	   if (i==0 || sqsum < short_col_sqr) 
	   short_col_sqr = sqsum;
	   e[i]=zero;
	   }
	   }
	*/

	/**
	 * ApplyBound: computes
	 *         bound_A <- max_i(max(sum_{j|a_ij > 0} a_ij, sum_{j|a_ij < 0} |a_ij|))
	 * this is useful because for all u, v >= 0:
	 *     [b has all entries in -u..v] 
         *       => [each entry of A.b is at most (u+v)*bound_A in absolute value]
	 */
	template <class Ring, class ItMatrix> //iterable matrix
	void ApplyBound(const Ring& R, typename Ring::Element& A_bound, const ItMatrix& A) {
		typedef typename Ring::Element Integer;
		Integer possum, negsum, zero;
		R.init(bound_A, 0);
		R.init(zero, 0);
		
		typename ItMatrix::ConstRowIterator rowIter;
		rowIter = A.rowBegin();

		for (; rowIter != A.rowEnd(); rowIter++) {
			typename ItMatrix::ConstRow::const_iterator elm;
			R.init(possum, 0);
			R.init(negsum, 0);
			for (elm = rowIter->begin(); elm != rowIter->end(); elm++) 
				if (*elm > zero) 
					R.addin(possum, *elm);
				else
					R.subin(negsum, *elm);
			
			if (possum > bound_A)
				R.assign(bound_A, possum);
			if (negsum > bound_A)
				R.assign(bound_A, negsum);
		}
	}

	template< class _Ring>
	class LiftingContainer {
		
	public:
		typedef _Ring Ring;
		typedef typename Ring::Element Integer;
				
		// return the length of container
		virtual const size_t length() const =0;

		// return the size of the solution
		virtual const int size() const = 0;

		// return the ring
		virtual const Ring& ring() const = 0;
				
		// return the prime
		virtual const Integer& prime () const = 0; 		

		virtual ~LiftingContainer () {}
		
	};

	template< class _Ring, class _IMatrix>
	class LiftingContainerBase : public LiftingContainer< _Ring> {

	public:		
		typedef _IMatrix                  IMatrix;
		typedef _Ring                        Ring;
		typedef typename _Ring::Element   Integer;		
		typedef std::vector<Integer>      IVector;
#ifdef RSTIMING
		mutable Timer ttSetup, tRingApply, tRingOther, ttRingOther, ttRingApply;
#endif

	protected:
						
		const IMatrix&            _A;
		Ring                      _R;
		Integer                   _p;
		IVector                   _b;
		VectorDomain<Ring>      _VDR;
		size_t               _length;
		Integer            _numbound;
		Integer            _denbound;
		MatrixApplyDomain<Ring,IMatrix>    _MAD;
		//BlasApply<Ring>          _BA;

		// Pascal's chunk-and-blas optimization:
		// instead of doing the ring multiplication A.digit, we write
		//     A = A0 + A1*2^16 + A2*2^32 + ... 
		// where A0, A1, ... (the 'chunks') are matrices of double
		// Then, we compute A.digit by multiplying each Ai.digit, and shifting and adding the results
		
		// To extend the method to when A has a negative entry, we add a final chunk that 
		// corresponds to a big negative number, using 2s complementing to write the negative number
		// as  -(big power of 2) + (sum of positive terms)
		// for example,   -0x000123456789 --> -(1 << 48) + (0xFFFE) << 32 + (0xDCBA) << 16 + (0x9877) (4 chunks)
		bool                         use_chunks;
		bool                            use_neg;
		size_t                       chunk_size;
		size_t                       num_chunks;
		double *                         chunks; //size: n*n*(num_chunks+1)
		

	public:

		template <class Prime_Type, class Vector1>
		LiftingContainerBase (const Ring& R, const IMatrix& A, const Vector1& b, const Prime_Type& p):
			_A(A), _R(R), _VDR(R), _MAD(R,A) {
#ifdef RSTIMING
			ttSetup.start();
#endif
			linbox_check(A.rowdim() == b.size());
			int n,m;
			n=A.rowdim();
			m=A.coldim();

			assert(m == n); //logic may not work otherwise
			// initialise the prime as an Integer
			_R.init(_p,p);
			
			// initialize res = b
			_b.resize(b.size());
			typename Vector1::const_iterator         b_iter    = b.begin();
			typename std::vector<Integer>::iterator  res_iter  = _b.begin() ;
			for (; b_iter != b.end(); ++res_iter, ++b_iter) 
				_R.init(*res_iter, *b_iter);
						
			Integer had_sq, short_sq;
			BoundBlackbox(_R, had_sq, short_sq, A);
			
			typename std::vector<Integer>::const_iterator iterb = _b.begin();
			Integer normb_sq;
			_R.init(normb_sq, 0);
			for (;iterb!=_b.end();++iterb)
				normb_sq += (*iterb)*(*iterb);

			LinBox::integer had_sqi, short_sqi, normb_sqi, N, D, L, prime;
			_R.convert(had_sqi, had_sq);
			_R.convert(short_sqi, short_sq);
			_R.convert(normb_sqi, normb_sq);
			_R.convert(prime,_p);
			D = sqrt(had_sqi) + 1;
			N = sqrt(had_sqi * normb_sqi / short_sqi) + 1;
			L = N * D * 2; 
			_length = logp(L,prime) + 1;   // round up instead of down
#ifdef DEBUG_LC                                   
			cout<<" norms computed, p = "<<_p<<"\n";
			cout<<" N = "<<N<<", D = "<<D<<", length = "<<_length<<"\n";
#endif
			_R.init(_numbound,N);
			_R.init(_denbound,D);

			/*
			// (1 << maxChunkVal) * (p-1) * n <= 1 << 53
 			LinBox::integer maxChunkVal = 1;
 			maxChunkVal <<= 53;
 			maxChunkVal /= (prime-1) * n;
 			chunk_size = -1;
 			while (maxChunkVal > 0) {
 				maxChunkVal /= 2;
 				chunk_size++;
 			}

			// ideally we would use chunks with chunk_size bits in them to make best use of the 
			// double representation, but for now it is only implemented to use a chunk size of 16

 			use_chunks = (chunk_size >= 16); 

			if (use_chunks) {
				chunk_size = 16;
				
				LinBox::integer tmp=0;
				size_t maxBitSize = 0;				
				use_neg = false;
				typename IMatrix::ConstRawIterator it = A.rawBegin();
				for (int i=0; i<n*n; i++, ++it) {
					_R.convert(tmp, *it);
					maxBitSize = max(maxBitSize, tmp.bitsize());
					use_neg |= (tmp < 0);
				}
							
				num_chunks = (maxBitSize / chunk_size)+ (((maxBitSize % chunk_size) > 0)? 1:0);
				if (num_chunks ==1)
					use_neg= false;

				if (use_neg) 
					num_chunks++; //the leading chunk will be negative
				//cerr<<"max bit size    :"<<maxBitSize<<endl;
				//cerr<<"total of chunks :"<<num_chunks<<endl;

				int n2 = n*n;
				chunks = new double[n2*num_chunks];
 				memset(chunks, 0, sizeof(double)*n*n*num_chunks);
				it = A.rawBegin();

				if (num_chunks ==1)
					for (int i=0; i<n2; i++, ++it) {
						_R.convert(*(chunks+i), *it);
					}
				else
					for (int i=0; i<n2; i++, ++it) {
						integer tmp;
						double* pdbl = chunks + i;
						_R.convert(tmp, *it);
						if (tmp >= 0) {
							size_t tmpsize    = tmp.size();
							size_t tmpbitsize = tmp.bitsize();
							
							for (size_t j=0; j<tmpsize-1; j++) {
								*pdbl = tmp[j] & 0xFFFF;
								*(pdbl+n2) = tmp[j] >> 16;
								pdbl += 2*n2;
							}
							if ((tmpbitsize % 32) > 16 ) {
								*pdbl = tmp[tmpsize-1] & 0xFFFF;
								*(pdbl+n2) = tmp[tmpsize-1] >> 16;						
							}
							else {
								*pdbl = tmp[tmpsize-1] & 0xFFFF;
							}
							
						}
						else {
							++tmp;
							// 						tmp *= -1;
							size_t tmpsize    = tmp.size();
							size_t tmpbitsize = tmp.bitsize();
							size_t j;
							
							for (j=0; j<tmpsize-1; j++) {
								*pdbl = 0xFFFF ^ (tmp[j] & 0xFFFF);
								*(pdbl+n2) = 0xFFFF ^ (tmp[j] >> 16);
								pdbl += 2*n2;							
							}
							if ((tmpbitsize % 32) > 16){
								*pdbl = 0xFFFF ^ (tmp[tmpsize-1] & 0xFFFF);
								*(pdbl+n2) = 0xFFFF ^ (tmp[tmpsize-1] >> 16);
								pdbl += 2*n2;
								j=tmpsize<<1;
							}
							else {
								*pdbl = 0xFFFF ^ (tmp[tmpsize-1] & 0xFFFF);
								pdbl += n2;
								j = (tmpsize<<1) -1;
							}
							
							//j+=tmpbitsize ; //convert from a word count to a 16-bit count
							for (; j<num_chunks-1; j++, pdbl += n2) 
								*pdbl = 0xFFFF;
							*pdbl = 1; //set the leading negative chunk for this entry
						}
					}
#ifdef DEBUG_CHUNK
				cout << num_chunks << " chunks of "<< chunk_size << " bits each" << endl;
				if (!use_neg) cout << "not ";
				cout << "using negative leading chunk" << endl;
				cout << "Contents of chunks: " << endl;
				for (size_t i=0; i<num_chunks; i++) {
					cout << "chunk " << i << endl;
					for (int j=0; j<n*n; j++) {
						cout << static_cast<long long>(chunks[i*n*n+j]);
						if ((j+1)%n) cout << ' '; else cout << endl;
					}
				}
#endif			       
				use_neg = !(!use_neg);
			}
			*/
			_MAD.setup();
			
#ifdef DEBUG_LC		
			cout<<"lifting container initialized\n";			
#endif
#ifdef RSTIMING
			ttSetup.stop();
			ttRingOther.clear();
			ttRingApply.clear();
#endif
		}
		
		virtual IVector& nextdigit (IVector& , const IVector&) const = 0;

		class const_iterator {
		private:
			std::vector<Integer>          _res;
			const LiftingContainerBase    &_lc;
			size_t                   _position;
		public:
			const_iterator(const LiftingContainerBase& lc,size_t end=0)
				: _res(lc._b), _lc(lc), _position(end) {}					
			
			/**
			 * @returns False if the next digit cannot be computed
			 * (probably indicates modulus is bad)
			 */
			bool next (IVector& digit)  {
				
				linbox_check (digit.size() == _lc._A.rowdim());	  
				// compute next p-adic digit
				_lc.nextdigit(digit,_res);
#ifdef RSTIMING
				_lc.tRingApply.start();
#endif

				// prepare for computing next digit				
				// update tmp_integerv = _A * digit				
				IVector v2 (_lc._A.coldim());

				// v2 = _A.digit
				_lc._MAD.applyV(v2,digit);
				
				/*
				if (!_lc.use_chunks)
					_lc._BA.applyV (v2, _lc._A, digit);
				else {
					int n = _lc._A.rowdim();
					int chunksize = _lc.chunk_size;
					double* ddigit = new double[n];
					for (int i=0; i<n; i++) {
						_lc._R.convert(ddigit[i], digit[i]);
					}
#ifdef DEBUG_CHUNK
 					cout << "digits: ";
					for (int i=0; i<n; i++) 
 						cout << digit[i] << ' ';
 					cout << endl;
#endif


					if (_lc.num_chunks == 1) {
						double *ctd = new double[n];
						cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n,
							    1, _lc.chunks, n, ddigit, 1, 0, ctd, 1);
						for (int i=0;i<n;++i)
							_lc._R.init(v2[i],ctd[i]);
					}
					else {
						//rc: number of vectors to recombine
						//(the idea is that to compute a polynomial in the base 2^chunksize
						// with <= 53 bits in each coefficient, we can instead OR nonoverlapping blocks
						// of bits and then add them at the end, like this:
						//      AAAACCCCEEEEGGGG   instead  AAAA << 12 + BBBB << 10 + CCCC << 8 + ...
						//    +   BBBBDDDDFFFF00      of     
						// also note that we need separate blocks for positive and negative entries)

						int rc = (52 / chunksize) + 1; //constant at 4 for now

						//rclen: number of bytes in each of these OR-ed vectors
						// needs room to hold (max long long) << (num_chunks * chunksize) 

						int rclen = _lc.num_chunks*2 + 5;

						// 					cout << "rc= " << rc << ", rclen = " << rclen << endl;

						unsigned char* combined = new unsigned char[rc*n*rclen];
						memset(combined, 0, rc*n*rclen);

						//order from major index to minor: combining index, component of sol'n, byte
					
						//compute a product (chunk times digit) for each chunk
						double* ctd = new double[n];

						for (size_t i=0; i<_lc.num_chunks; i++) {
							//ctd <- A[i] . digit
							cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n,
								    1, _lc.chunks + (n*n*i), n, ddigit, 1, 0, ctd, 1);
							//cout << "chunk " << i << " times digit : ";
							//for (int j=0; j<n; j++) cout << (long long)ctd[j] << ' ';
							//cout << endl;
						
							if (!_lc.use_neg || i<_lc.num_chunks-1)
								for (int j=0; j<n; j++) {
									// up to 53 bits will be ored-in, to be summed later
									unsigned char* bitDest = combined;
									bitDest += rclen*((i % rc)*n+j);
									//{
									//cout << "rc[" << (i%rc) << ","<<
									//j<<"]:";
									//for (int i=0; i<rclen; i++) 
									//cout << (int)bitDest[i] << ' ';
									//cout << endl;
									//}
									//cout << "ctd[j]: " << (long long)ctd[j] << endl;
									long long mask = static_cast<long long>(ctd[j]);
									bitDest += 2*i;
									//mask <<= (i*chunksize) % 8; //useless when chunksize=16
									*((long long*) bitDest) |= mask; 
									//bitDest -= 2*i;
									//{
									//cout << "rc[" << (i%rc) << ","<<
									//j<<"]:";
									//for (int i=0; i<rclen; i++) 
									//cout << (int)bitDest[i] << ' ';
									//cout << endl;
									//}
								}
						}
						for (int i=0; i<n; i++) {
							LinBox::integer result, tmp;
							if (_lc.use_neg) {
								result = -ctd[i];
								result <<= (_lc.num_chunks-1)*16;
#ifdef DEBUG_CHUNK
								cout << "rcneg: " << result << endl;
#endif
							}
							else
								result = 0;

							for (int j=0; j<rc; j++) {
								unsigned char* thispos = combined + rclen*(j*n+i);
								importWords(tmp, rclen, -1, 1, 0, 0, thispos);
								result += tmp;
#ifdef DEBUG_CHUNK
								cout << "rc[" << j << "," << i << "]:" << tmp << endl;
#endif
							}
#ifdef DEBUG_CHUNK
							cout << "v2[" << i << "]:" << result  << endl;
#endif
							_lc._R.init(v2[i], result);
						}
						delete[] combined;
						delete[] ctd;
					}
				}
				*/				
#ifdef RSTIMING
				_lc.tRingApply.stop();
				_lc.ttRingApply += _lc.tRingApply;
				_lc.tRingOther.start();
#endif

				// update _res -= v2 
				_lc._VDR.subin (_res, v2);
				typename std::vector<Integer>::iterator p0;
				// update _res = _res / p
				int index=0;
				for ( p0 = _res.begin(); p0 != _res.end(); ++ p0, ++index){ 
#ifdef LC_CHECK_DIVISION
 					if (! _lc._R.isDivisor(*p0,_lc._p)) {
 						cerr<<"residue "<<*p0<<" not divisible by modulus "<<_lc._p<<endl;
 						cout<<"residue "<<*p0<<" not divisible by modulus "<<_lc._p<<endl;
 						return false;
 					}
#endif
					_lc._R.divin(*p0, _lc._p);
				}

				_position++;
#ifdef RSTIMING				
				_lc.tRingOther.stop();
				_lc.ttRingOther += _lc.tRingOther;
#endif
				return true;
			}

			bool operator != (const const_iterator& iterator) const {
				if ( &_lc != &iterator._lc)
					assert("try to compare differents LiftingContainer::const_iterator , abort\n");
				return _position != iterator._position;					
			}
	
			bool operator == (const const_iterator& iterator) const {
				if ( &_lc != &iterator._lc)
					assert("try to compare differents LiftingContainer::const_iterator , abort\n");
				return _position == iterator._position;					
			}

		};

		/** Bit manipulation function for possible use in optimization
		 // efficiently pulls out continuous blocks of bits, from lsb to msb inclusive
		 // least significant bits start at index 0, so msb >= lsb
		 // if any bits with index >= 8*numBytes are asked for they will be zeroes
		 static long long bytesToBits(unsigned char * byteArray, size_t numBytes, size_t lsb, size_t msb) { 
		 linbox_check(msb >= lsb);
		 size_t lsbi = lsb >> 3;
		 size_t msbi = msb >> 3;
		 if (msbi == lsbi) 
		 if (msbi >= numBytes) 
		 return 0;
		 else
		 return (byteArray[lsbi] >> (lsb & 7)) & ((1 << (msb - lsb + 1)) - 1);

		 long long result = (msbi < numBytes) ? (byteArray[msbi] & ((1 << ((msb & 7)+1)) - 1)) : 0;
		 for (size_t i=msbi-1; i>lsbi; i--) {
		 result <<= 8;
		 result |= (i < numBytes) ? byteArray[i] : 0;
		 }
		 result <<= 8 - (lsb & 7);
		 result |= (lsbi < numBytes) ? (byteArray[lsbi] >> (lsb & 7)) : 0;

		 return result;
		 }
		*/

		const_iterator begin() const {
			return const_iterator(*this);
		}		

		const_iterator end() const {
			return const_iterator (*this,_length);
		}
	
		virtual const size_t length() const
		{ return _length;}
		
		// return the size of the solution
		virtual const int size() const 
		{ return _A.coldim();}

		// return the ring
		virtual const Ring& ring() const 
		{ return _R; }
		
		// return the prime
		virtual const Integer& prime () const 
		{ return _p; }
		
		// return the bound for the numerator
		const Integer numbound() const
		{return _numbound;}
		
		// return the bound for the denominator
		const Integer denbound() const
		{return _denbound;}
		
	};

	template <class _Ring, class _Field, class _IMatrix, class _FMatrix>
	class DixonLiftingContainer : public LiftingContainerBase< _Ring, _IMatrix> {

	public:
		typedef _Field                               Field;
		typedef _Ring                                 Ring;
		typedef _IMatrix                           IMatrix;		
		typedef _FMatrix                           FMatrix;	
		typedef typename Field::Element            Element;
		typedef typename IMatrix::Element          Integer;
		typedef std::vector<Integer>               IVector;
		typedef std::vector<Element>               FVector;

	protected:
		
		const FMatrix&                  _Ap;	
		Field                            _F;
		const VectorDomain<Field>      _VDF;
		mutable FVector              _res_p;
		mutable FVector            _digit_p;		
		BlasApply<Field>                _BA;
		
	public:
#ifdef RSTIMING
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert;
#endif

		template <class Prime_Type, class VectorIn>
		DixonLiftingContainer (const Ring&       R, 
				       const Field&      F, 
				       const IMatrix&    A, 
				       const FMatrix&   Ap, 
				       const VectorIn&   b, 
				       const Prime_Type& p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap), _F(F), _VDF(F), 
			  _res_p(b.size()), _digit_p(A.coldim()), _BA(F) 
		{

			//Ap.write(cout,F);
#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
#endif
		}
		
		
		virtual ~DixonLiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit, const IVector& residu) const {		
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			LinBox::integer tmp;

			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, _R.convert(tmp,*iter));
			}
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
			tGetDigit.start();
#endif			
			// compute the solution by applying the inverse of A mod p
			_BA.applyV(_digit_p,_Ap,_res_p);
#ifdef RSTIMING
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					_R.init(*iter, _F.convert(tmp,*iter_p));
			}
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif
			return digit;
		}
		
	}; // end of class DixonLiftingContainerBase



	template <class _Ring, class _Field, class _IMatrix, class _FMatrix, class _FPolynomial>
	class WiedemannLiftingContainer : public LiftingContainerBase<_Ring, _IMatrix> {

	public:
		typedef _Field                                     Field;
		typedef _Ring                                       Ring;
		typedef _IMatrix                                 IMatrix;		
		typedef _FMatrix                                 FMatrix;
		typedef typename Field::Element                  Element;
		typedef typename Ring::Element                   Integer;
		typedef std::vector<Integer>                     IVector;
		typedef std::vector<Element>                     FVector;
		typedef _FPolynomial                         FPolynomial;
		typedef typename FPolynomial::iterator     FPolyIterator;

	protected:

		const FMatrix                  &_Ap;
		mutable FPolynomial        _MinPoly;
		Field                            _F;
		const VectorDomain<Field>      _VDF;
		mutable FVector              _res_p;
		mutable FVector            _digit_p;
		typename Field::RandIter      _rand;		
#ifdef RSTIMING
	public:
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert;
#endif			
	public:

		template <class Prime_Type, class VectorIn>
		WiedemannLiftingContainer (const Ring& R,
					   const Field& F, 
					   const IMatrix& A, 
					   const FMatrix& Ap, 
					   const FPolynomial& MinPoly,
					   const VectorIn& b, 
					   const Prime_Type& p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap), _MinPoly(MinPoly), _F(F), _VDF(F), _res_p(b.size()), _digit_p(A.coldim()), _rand(F) {

			// Normalize the minimal polynomial as f(x)=1- a1/a0 x - a2/a0 x^2 - ...
			FPolyIterator iter=_MinPoly.begin();
			while(++iter != _MinPoly.end ()){
				_F.divin (*iter, _MinPoly.front ());
				_F.negin (*iter);
			}
#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
#endif
		}

		virtual ~WiedemannLiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit,const IVector& residu) const {
			
			LinBox::integer tmp;
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, _R.convert(tmp,*iter));
			}
#ifdef RSTIMING			
			tGetDigitConvert.stop();
			ttGetDigitConvert+=tGetDigitConvert;
			tGetDigit.start();
#endif		
			// compute the solution of system by Minimal polynomial application
			_VDF.mul (_digit_p, _res_p, _MinPoly.back ());      
			FVector z(_Ap.rowdim ());      
			
			
			for (size_t i = _MinPoly.size () - 1; --i > 0;) {												
				_Ap.apply (z, _digit_p);
				_VDF.axpy (_digit_p, _MinPoly[i], _res_p, z);
			}      

			      
			// check results
			FVector error(_Ap.coldim());
			_Ap.apply(error,_digit_p);

			bool nosolution = false;
			int nosolution_threshold=5;
			int nst=0;
			size_t minpoly_degree;
			// until the digit is incorrect update the minpoly and recompute the digit
			while (!_VDF.areEqual(error,_res_p) && !nosolution ){
				minpoly_degree=_MinPoly.size();
				FPolynomial Poly;
				unsigned long deg;
				unsigned long size= (_Ap.rowdim() - _MinPoly.size())<<1 ;
				BlackboxContainer<Field, FMatrix > Sequence(&_Ap,_F,error,size);
				MasseyDomain<Field,BlackboxContainer<Field, FMatrix > > MD(&Sequence);
				MD.minpoly(Poly,deg);
				if (_F.isZero(Poly.front())) {
					// here we should stop the execution but not yet implemented
					cerr<<" the prime was not good \n, result will be wrong";
					break;
				}
				
				// denormalize the minimal polynomial
				FPolyIterator iter=_MinPoly.begin();
				while (++iter != _MinPoly.end()) {
					_F.mulin (*iter, _MinPoly.front());
					_F.negin (*iter);
				}	
		       	
				// update the minimal polynomial 
				Element zero;
				_F.init(zero,0);
				FPolynomial newMinPoly(_MinPoly.size()+Poly.size()-1,zero);
				for (size_t i=0; i < _MinPoly.size(); i++)
					for (size_t j=0 ; j < Poly.size(); j++)
						_F.axpyin(newMinPoly[i+j],_MinPoly[i],Poly[j]);	
				_MinPoly.clear();
				Poly.clear();
				_MinPoly=newMinPoly;

				// normalize the new minimal polynomial
				iter=_MinPoly.begin ();
				while (++iter != _MinPoly.end ()) {
					_F.divin (*iter, _MinPoly.front ());
					_F.negin (*iter);
				}
				
				_VDF.mul (_digit_p, _res_p, _MinPoly.back ());      
				FVector z(_Ap.rowdim ());      
				for (size_t i = _MinPoly.size () - 1; --i > 0;) {		
					_Ap.apply (z, _digit_p);
					_VDF.axpy (_digit_p, _MinPoly[i], _res_p, z);
				}    

				_Ap.apply(error,_digit_p);	 
				if (_MinPoly.size() > minpoly_degree){
					minpoly_degree = _MinPoly.size();
					nst=0;cerr<<"updating minpoly\n";
				}
				else {
					if (nst < nosolution_threshold) nst++;
					else{
						nosolution=true;
						throw PreconditionFailed (__FUNCTION__, __LINE__, "system is inconsistent or the choosen prime leads to inconsistent resolution");
					}
				}
			}
#ifdef RSTIMING
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					_R.init(*iter, _F.convert(tmp,*iter_p));
			}
			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif
			return digit;
		}


	}; // end of class WiedemannLiftingContainerBase


	// Block Wiedemann LiftingContianer
	template <class _Ring, class _Field, class _IMatrix, class _FMatrix>
	class BlockWiedemannLiftingContainer : public LiftingContainerBase<_Ring, _IMatrix> {

	public:
		typedef _Field                                	            Field;
		typedef _Ring                                 	             Ring;
		typedef _IMatrix                              	          IMatrix;       
		typedef _FMatrix                              	          FMatrix;
		typedef typename Field::Element               	          Element;
		typedef typename Ring::Element                            Integer;
		typedef std::vector<Integer>                              IVector;
		typedef std::vector<Element>                              FVector;
		typedef BlasMatrix<Element>                           Coefficient;
		typedef BlasMatrix<Element>                                 Block;
		typedef std::vector<Coefficient>                 FBlockPolynomial;
		typedef BlackboxBlockContainerRecord<Field, FMatrix>     Sequence;


	protected:

		const FMatrix                       &_Ap;
		Field                                 _F;
		const VectorDomain<Field>           _VDF;
		mutable FVector                   _res_p;
		mutable FVector                 _digit_p;
		typename Field::RandIter           _rand;
		size_t                              _row;
		size_t                              _col;
		size_t                                _m;
		size_t                                _n;
		Block                                 _U;
		BlasMatrixDomain<Field>             _BMD;
		Sequence                           *_Seq;
		BlockMasseyDomain<Field,Sequence>  *_Dom;
#ifdef RSTIMING
	public:
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert, tMinPoly, ttMinPoly;
#endif			
	public:

		template <class Prime_Type, class VectorIn>
		BlockWiedemannLiftingContainer (const Ring                         &R,
						const Field                        &F, 
						const IMatrix                      &A, 
						const FMatrix                     &Ap, 
						const VectorIn                     &b, 
						const Prime_Type                   &p,
						const size_t                        m,
						const size_t                        n)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap),
			  _F(F), 
			  _VDF(F), 
			  _res_p(b.size()), 
			  _digit_p(A.coldim()), 
			  _rand(F), 
			  _row(Ap.rowdim()), 
			  _col(Ap.coldim()), 
			  _m(m), 
			  _n(n), 
			  _U(m-1,Ap.rowdim()), 
			  _BMD(F) {
					
			for (size_t i=0;i<_m-1;++i)
				for (size_t j=0;j< _row;++j)
					_rand.random(_U.refEntry(i,j));
			
			Coefficient V(_col,n);
			for (size_t i=0;i< _col;++i)
				for (size_t j=0; j< n;++j)
					_rand.random(V.refEntry(i,j));
			
			
			Block UAp(_m, _row);

			typename Block::ConstRowIterator    iter_U   = _U.rowBegin();
			typename Block::RowIterator         iter_UAp = UAp.rowBegin();
			++iter_UAp;
			for (; iter_U != _U.rowEnd(); ++iter_UAp, ++iter_U) 
				Ap.applyTranspose( *iter_UAp , *iter_U );
				
			for (size_t i=0;i<m;++i)
				_rand.random(UAp.refEntry(0,i));			

			_Seq = new Sequence (&Ap, _F, UAp,V);
			_Dom = new BlockMasseyDomain<Field,Sequence> (_Seq);					      					      

#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
			ttMinPoly.clear();
#endif
		}

		virtual ~BlockWiedemannLiftingContainer() {
#ifdef _BM_TIMING
			_Dom->printTimer();
#endif
#ifdef _BBC_TIMING
			_Seq->printTimer();
#endif
			delete _Seq;
			delete _Dom;
		}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit,const IVector& residu) const {
			
			LinBox::integer tmp;
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, _R.convert(tmp,*iter));
			}
#ifdef RSTIMING			
			tGetDigitConvert.stop();
			ttGetDigitConvert+=tGetDigitConvert;
			tGetDigit.start();
#endif		
			// compute the Minimal polynomial of the modified Sequence
			_Seq->setU(_res_p,0);
			_Seq->recompute();
			FBlockPolynomial minpoly;
			std::vector<size_t> degree(_m);

#ifdef RSTIMING
			tMinPoly.start();
#endif
			_Dom->left_minpoly(minpoly,degree);
#ifdef RSTIMING
			tMinPoly.stop();
			ttMinPoly+=tMinPoly;
#endif		
	
			size_t idx=0;
			if ( _F.isZero(minpoly[0].getEntry(0,0))) {
				size_t i=1;
				while ( _F.isZero(minpoly[0].getEntry(i,0)))
					++i;
				if (i == _m)
					throw LinboxError(" block minpoly: matrix seems to be singular - abort");
				else 
					idx=i	;			
			}

			size_t deg = degree[idx];		
			BlasMatrix<Element> idx_poly(deg+1,_m-1);
			for (size_t i=0;i<deg+1;++i) 
				for (size_t j=0;j<_m-1;++j)
					idx_poly.setEntry(i,j,minpoly[i].getEntry(idx,j+1));

			BlasMatrix<Element> Combi(deg+1,_row);
			_BMD.mul(Combi,idx_poly,_U);
					

			FVector lhs(_col),row(_row);	
			for (size_t i=0;i<_row;++i)
				row[i]= Combi.getEntry(deg,i);
					
			_Ap.applyTranspose(lhs,row);					
			FVector lhsbis(lhs);
			for (int i = deg-1 ; i >= 0;--i) {
				for (size_t j=0;j<_row;++j)
					row[j]= Combi.getEntry(i,j);
				_VDF.add (lhs,row,lhsbis);
				_Ap.applyTranspose (lhsbis, lhs);			
			}   
					
			FVector accu (lhs);
			_Ap.applyTranspose(lhs,_res_p);
			_VDF.mulin(lhs,minpoly[deg].getEntry(idx,0));
			lhsbis=lhs;
			for (size_t i = deg-1 ; i > 0;--i) {
				_VDF.axpy (lhs,minpoly[i].getEntry(idx,0) , _res_p, lhsbis);
				_Ap.applyTranspose (lhsbis, lhs);			
			}  
					
			_VDF.addin(accu,lhs);
			Element scaling;
			_F.init(scaling);
			_F.neg(scaling,minpoly[0].getEntry(idx,0));
			_F.invin(scaling);
			_VDF.mul(_digit_p,accu,scaling);
					
			// check results
			FVector error(_Ap.coldim());
			_Ap.apply(error,_digit_p);		
			if (!_VDF.areEqual(error,_res_p)){
				cout<<"BlockMinpoly error\n";
				throw LinboxError("BlockMinpoly error\n");
			}
			
			
#ifdef RSTIMING			
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					_R.init(*iter, _F.convert(tmp,*iter_p));
			}
			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif
			return digit;
		}


	}; // end of class WiedemannLiftingContainerBase




} // end of namespace LinBox
#endif
