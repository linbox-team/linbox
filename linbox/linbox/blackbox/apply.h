/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* File: apply.h
 *  Author: Zhendong Wan
 *  Modified by Pascal Giorgi
 */

/* Reserve for possible optimal.
 */

#ifndef __LINBOX_APPLY_H__
#define __LINBOX_APPLY_H__

#include <linbox-config.h>
#include <linbox/integer.h>
#include <linbox/util/debug.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/field/ntl-ZZ.h>
#include <linbox/algorithms/lifting-container.h>
#include <vector>

#ifdef __LINBOX_BLAS_AVAILABLE
#include <linbox/fflas/fflas.h>
#endif

//#define DEBUG_CHUNK

namespace LinBox {
	
	// general case, y = A x
	template<class OutV, class Matrix, class InV>
	inline OutV& apply (OutV& y, const Matrix& A, const InV& x) {
		
		return A. apply (y, x);

	}

	template<class OutV, class Matrix, class InV>
	inline OutV& applyTranspose (OutV& y, const Matrix& A, const InV& x) {

		return A. applyTranspose (y, x);
		
	}

	
	template<class Domain>
	class BlasApply {
	  	 
	public:
		typedef typename Domain::Element    Element;
		typedef std::vector<Element>         Vector;

		BlasApply(const Domain& D) : _D(D), _MD(D) { 
			_D.characteristic(_prime); 
			_D.init(_one,1UL); 
			_D.init(_zero,0UL);
		}
	  	  		

		//#ifdef __LINBOX_BLAS_AVAILABLE
		inline Vector& applyV(Vector                        &y,
				      const BlasMatrix<Element>     &A, 
				      const Vector                  &x) const {
	    
			if (( _prime > 0) && ( _prime <  67108863)) {

				FFLAS::fgemv( _D, FFLAS::FflasNoTrans, 
					      A.rowdim(), A.coldim(),
					      _one,
					      A.getPointer(), A.getStride(),
					      &x[0],1,
					      _zero,
					      &y[0],1);  	      
			}
			else {
				_MD.vectorMul (y, A, x);	      
			}
			return y;
		}

		inline Vector& applyVTrans(Vector                        &y,
					   BlasMatrix<Element>           &A,
					   const Vector                  &x) const {
	    
			if (( _prime > 0) && ( _prime <  67108863)) {

				FFLAS::fgemv( _D, FFLAS::FflasTrans, 
					      A.rowdim(), A.coldim(),
					      _one,
					      A.getPointer(), A.getStride(),
					      &x[0],1,
					      _zero,
					      &y[0],1);  	      
			}
			else {
				TransposeMatrix<const BlasMatrix<Element> > B(A); 
				_MD.vectorMul (y, B, x);
			}
			return y;
		}		
		//#endif
	  
	private:
		Domain         _D;
		integer    _prime;
		Element _one,_zero;
		MatrixDomain<Domain> _MD;
	  

	};


	template <class Domain, class IMatrix>
	class MatrixApplyDomain {

	public:
		typedef typename Domain::Element    Element;
		typedef std::vector<Element>         Vector;
		
		MatrixApplyDomain(const Domain& D, const IMatrix &M) : _D(D), _M(M) {}

		void setup(LinBox::integer prime){}
				
		Vector& applyV(Vector& y, Vector& x) const { return _M.apply(y,x);}

		Vector& applyVTrans(Vector& y, Vector& x) const {return _M.applyTranspose(y,x);}

	private:
		Domain          _D;
		const IMatrix  &_M;
	};


	template <class Domain, class IMatrix>
	class BlasMatrixApplyDomain {
		
	public:
		typedef typename Domain::Element    Element;
		typedef std::vector<Element>         Vector;
		typedef IMatrix                       Matrix;
		
		
	       	BlasMatrixApplyDomain(const Domain& D, const IMatrix &M) : _D(D), _M(M), _MD(D), _m(M.rowdim()), _n(M.coldim()) {}
			

		~BlasMatrixApplyDomain () {
			if (use_chunks)
				delete[] chunks;
		}

		void setup(LinBox::integer prime){//setup
			// Compute the maximum size of chunks
			LinBox::integer maxChunkVal = 1;
 			maxChunkVal <<= 53;
 			maxChunkVal /= (prime-1) * _n;
 			chunk_size = -1;
 			while (maxChunkVal > 0) {
 				maxChunkVal /= 2;
 				chunk_size++;
 			}
					
			use_chunks = (chunk_size >= 16); 
			
			if (use_chunks){//usechunk
				// set maximum size of chunk to 16			
				chunk_size = 16;
				
				// compute the magnitude in bit of the matrix
				// check if we need negative representation
				LinBox::integer tmp=0;
				size_t maxBitSize = 0;				
				use_neg = false;
				typename Matrix::ConstRawIterator it = _M.rawBegin();
				for (size_t i=0; i<_m*_n; i++, ++it) {
					_D.convert(tmp, *it);
					maxBitSize = max(maxBitSize, tmp.bitsize());
					use_neg |= (tmp < 0);
				}
				
				// compute the number of chunk
				num_chunks = (maxBitSize / chunk_size)+ (((maxBitSize % chunk_size) > 0)? 1:0);
				if (num_chunks ==1)
					use_neg= false;

				if (use_neg) 
					num_chunks++; //the leading chunk will be negative
			
				int n2 = _m*_n;
				chunks = new double[n2*num_chunks];
				memset(chunks, 0, sizeof(double)*_n*_n*num_chunks);
				it = _M.rawBegin();

				if (num_chunks ==1)
					for (int i=0; i<n2; i++, ++it) {
						_D.convert(*(chunks+i), *it);
					}
				else
					for (int i=0; i<n2; i++, ++it) {
						integer tmp;
						double* pdbl = chunks + i;
						_D.convert(tmp, *it);
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
					for (size_t j=0; j<_m*_n; j++) {
						cout << static_cast<long long>(chunks[i*_m*_n+j]);
						if ((j+1)%n) cout << ' '; else cout << endl;
					}
				}
#endif			       
				use_neg = !(!use_neg);
			}
		}
	

			
		Vector& applyV(Vector& y, Vector& x) const {//applyV 		
		
			linbox_check( _n == x.size());
			linbox_check( _m == y.size());

			if (!use_chunks){
				_MD.vectorMul (y, _M, x);			
			}
			else{

				double* dx = new double[_n];
				for (size_t i=0; i<_n; i++) {
					_D.convert(dx[i], x[i]);
				}
#ifdef DEBUG_CHUNK
				cout << "x: ";
				for (size_t i=0; i<_n; i++) 
					cout << x[i] << ' ';
				cout << endl;
#endif

	
				if (num_chunks == 1) {
					double *ctd = new double[_n];
					cblas_dgemv(CblasRowMajor, CblasNoTrans, _m, _n,
						    1, chunks, _n, dx, 1, 0, ctd, 1);
					for (size_t i=0;i<_n;++i)
						_D.init(y[i],ctd[i]);
					delete[] ctd;
				}
				else {
					//rc: number of vectors to recombine
					//(the idea is that to compute a polynomial in the base 2^chunksize
					// with <= 53 bits in each coefficient, we can instead OR nonoverlapping blocks
					// of bits and then add them at the end, like this:
					//      AAAACCCCEEEEGGGG   instead  AAAA << 12 + BBBB << 10 + CCCC << 8 + ...
					//    +   BBBBDDDDFFFF00      of     
					// also note that we need separate blocks for positive and negative entries)
		
					int rc = (52 / chunk_size) + 1; //constant at 4 for now
		
					//rclen: number of bytes in each of these OR-ed vectors
					// needs room to hold (max long long) << (num_chunks * chunksize) 
		
					int rclen = num_chunks*2 + 5;
		
					// 					cout << "rc= " << rc << ", rclen = " << rclen << endl;
		
					unsigned char* combined = new unsigned char[rc*n*rclen];
					memset(combined, 0, rc*n*rclen);

					//order from major index to minor: combining index, component of sol'n, byte
		
					//compute a product (chunk times x) for each chunk
					double* ctd = new double[n];
		
					for (size_t i=0; i<num_chunks; i++) {
						cblas_dgemv(CblasRowMajor, CblasNoTrans, _m, _n, 1, chunks + (_m*_n*i), _n, dx, 1, 0, ctd, 1);
			
						if (!use_neg || i<num_chunks-1)
							for (size_t j=0; j<_n; j++) {
								// up to 53 bits will be ored-in, to be summed later
								unsigned char* bitDest = combined;
								bitDest += rclen*((i % rc)*_n+j);
								long long mask = static_cast<long long>(ctd[j]);
								bitDest += 2*i;
								*((long long*) bitDest) |= mask; 
							}
					}
		
					delete[] dx;
					
					for (size_t i=0; i<_n; i++) {
						LinBox::integer result, tmp;
						if (use_neg) {
							result = -ctd[i];
							result <<= (num_chunks-1)*16;
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
						_D.init(y[i], result);
					}
					delete[] combined;
					delete[] ctd;
				}
			}
			return y;
		}


		Vector& applyVTrans(Vector& y, Vector& x) const {
			TransposeMatrix<IMatrix> B(_M); 
			return _MD.vectorMul (y, B, x);
		}

	protected:
		Domain                             _D;
		const IMatrix                     &_M;
		MatrixDomain<Domain>              _MD;
		size_t                             _m;
		size_t                             _n;

		bool           use_chunks;
		bool              use_neg;
		size_t         chunk_size;
		size_t         num_chunks;
		double *           chunks; 


	};
	
	template<>
	template <class Domain>
	class MatrixApplyDomain<Domain, BlasMatrix<typename Domain::Element> > : public BlasMatrixApplyDomain<Domain, BlasMatrix<typename Domain::Element> > {

	public:
		MatrixApplyDomain (const Domain &D, const  BlasMatrix<typename Domain::Element> &M)
			: BlasMatrixApplyDomain<Domain, BlasMatrix<typename Domain::Element> > (D,M) {}
		
	};

	template<>
	template <class Domain>
	class MatrixApplyDomain<Domain, DenseMatrix<Domain> > : public BlasMatrixApplyDomain<Domain, DenseMatrix<Domain> > {
		
	public:
		MatrixApplyDomain (const Domain &D, const DenseMatrix<Domain> &M)
			: BlasMatrixApplyDomain<Domain, DenseMatrix<Domain> > (D,M) {}
	};
	


	/*
	template<>
	template <class Domain>
	class MatrixApplyDomain<Domain, BlasMatrix<typename Domain::Element> > {//class

	public:
		typedef typename Domain::Element    Element;
		typedef std::vector<Element>         Vector;
		typedef BlasMatrix<Element>          Matrix;
		
		MatrixApplyDomain(const Domain& D, const BlasMatrix<Element> &M) : _D(D), _M(M), _MD(D), _m(M.rowdim()), _n(M.coldim()) {}
			
		void setup(LinBox::integer prime){//setup
			// Compute the maximum size of chunks
			LinBox::integer maxChunkVal = 1;
 			maxChunkVal <<= 53;
 			maxChunkVal /= (prime-1) * _n;
 			chunk_size = -1;
 			while (maxChunkVal > 0) {
 				maxChunkVal /= 2;
 				chunk_size++;
 			}

			use_chunks = (chunk_size >= 16); 
			
			if (use_chunks){//usechunk
				// set maximum size of chunk to 16			
				chunk_size = 16;
				
				// compute the magnitude in bit of the matrix
				// check if we need negative representation
				LinBox::integer tmp=0;
				size_t maxBitSize = 0;				
				use_neg = false;
				typename Matrix::ConstRawIterator it = _M.rawBegin();
				for (int i=0; i<n*n; i++, ++it) {
					_D.convert(tmp, *it);
					maxBitSize = max(maxBitSize, tmp.bitsize());
					use_neg |= (tmp < 0);
				}
					
				// compute the number of chunk
				num_chunks = (maxBitSize / chunk_size)+ (((maxBitSize % chunk_size) > 0)? 1:0);
				if (num_chunks ==1)
					use_neg= false;

				if (use_neg) 
					num_chunks++; //the leading chunk will be negative
			

				int n2 = _n*_n;
				chunks = new double[n2*num_chunks];
				memset(chunks, 0, sizeof(double)*_n*_n*num_chunks);
				it = _M.rawBegin();

				if (num_chunks ==1)
					for (int i=0; i<n2; i++, ++it) {
						_D.convert(*(chunks+i), *it);
					}
				else
					for (int i=0; i<n2; i++, ++it) {
						integer tmp;
						double* pdbl = chunks + i;
						_D.convert(tmp, *it);
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
		}
	

			
		Vector& applyV(Vector& y, Vector& x) const {//applyV 		
		
			linbox-check( _n == x.size);
			linbox_check( _m == y.size);

			if (!use_chunks){
				_MD.vectorMul (y, _M, x);			
			}
			else{

				double* dx = new double[_n];
				for (int i=0; i<_n; i++) {
					_D.convert(dx[i], x[i]);
				}
#ifdef DEBUG_CHUNK
				cout << "x: ";
				for (int i=0; i<_n; i++) 
					cout << x[i] << ' ';
				cout << endl;
#endif

	
				if (num_chunks == 1) {
					double *ctd = new double[_n];
					cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n,
						    1, chunks, _n, dx, 1, 0, ctd, 1);
					for (int i=0;i<_n;++i)
						_D.init(y[i],ctd[i]);
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
		
					//compute a product (chunk times x) for each chunk
					double* ctd = new double[n];
		
					for (size_t i=0; i<num_chunks; i++) {
						cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1, chunks + (_n*_n*i), _n, dx, 1, 0, ctd, 1);
			
						if (!use_neg || i<num_chunks-1)
							for (int j=0; j<_n; j++) {
								// up to 53 bits will be ored-in, to be summed later
								unsigned char* bitDest = combined;
								bitDest += rclen*((i % rc)*_n+j);
								long long mask = static_cast<long long>(ctd[j]);
								bitDest += 2*i;
								*((long long*) bitDest) |= mask; 
							}
					}
		
					for (int i=0; i<_n; i++) {
						LinBox::integer result, tmp;
						if (use_neg) {
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
						_D.init(y[i], result);
					}
					delete[] combined;
					delete[] ctd;
				}
			}
			return y;
		}


		Vector& applyVTrans(Vector& y, Vector& x) const {
			TransposeMatrix<BlasMatrix<Element> > B(_M); 
			return _MD.vectorMul (y, B, x);
		}

	private:
		Domain                             _D;
		const BlasMatrix<Element>         &_M;
		MatrixDomain<Domain>              _MD;
		size_t                             _m;
		size_t                             _n;

		bool           use_chunks;
		bool              use_neg;
		size_t         chunk_size;
		size_t         num_chunks;
		double *           chunks; 
	};

	*/

} // end of namespace LinBox		
#endif
