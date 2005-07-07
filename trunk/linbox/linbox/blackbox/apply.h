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
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/field/ntl-ZZ.h>
#include <linbox/algorithms/lifting-container.h>
#include <vector>

#ifdef __LINBOX_BLAS_AVAILABLE
#include <linbox/fflas/fflas.h>
#endif

//#define DEBUG_CHUNK_SETUP
//#define DEBUG_CHUNK_APPLYV
//#define DEBUG_CHUNK_APPLYM

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

		inline Vector& applyVspecial (Vector                        &y,
					      BlasMatrix<Element>           &A,
					      const Vector                  &x) const {//toto
			
			size_t m,n;
			m = A.rowdim();
			n = A.coldim();
			linbox_check( x.size() == n);
			
			double * At_dbl = new double[m*n];
			
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					_D.convert(*(At_dbl+i+j*m), A.refEntry(i,j));
			
			integer tmp;
			bool use_neg=false;
			size_t maxword=0;
			for (size_t i=0;i<n;++i){
				_D.convert(tmp,x[i]);
				if (tmp <0)
					use_neg = true;
				if ( maxword < tmp.size())
					maxword= tmp.size();				
			}
			
			if (use_neg)
				maxword++;
			double *xdbl= new double[n*maxword];
			memset(xdbl, 0, sizeof(double)*n*maxword);
			
			for (size_t i=0;i<n;++i){
				_D.convert(tmp,x[i]);
				double * ptr= xdbl+i;
				if (tmp == 0)
					*ptr=0;
				else {
					if (tmp > 0) {
						for (size_t j=0;j<tmp.size();++j){
							*ptr= (double)x[i][j];
							ptr+= n;
						}
					}
					else {
						size_t j=0;
						for (;j<tmp.size();++j){
							*ptr= 0xFFFFFFFF^x[i][j];
							ptr+= n;
						}
						for (;j<maxword-1;++j){
							*ptr= 0xFFFFFFFF;
							ptr+= n;
						}
						*ptr=1;
					}
				}					
			}


			double *ydbl= new double[maxword*m];
					
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
				    maxword,m,n, 1,
				    xdbl,n, At_dbl, m, 0, ydbl, m);

			delete At_dbl;
			delete xdbl;

			size_t rclen=maxword*4+5;
			unsigned char* combined1 = new unsigned char[m*rclen];
			unsigned char* combined2 = new unsigned char[m*rclen];
			memset(combined1, 0, m*rclen);
			memset(combined2, 0, m*rclen);

			for (size_t i=0;i<m;++i){
				unsigned char *ptr= combined1+i*rclen;
				for (size_t j=0;j< maxword;j=j+2){
					if (!use_neg || j< maxword-1){						
						long long mask = static_cast<long long>(ydbl[j*m+i]);
						*((long long*) ptr) |= mask; 
						ptr+=4;
					}
				}
				ptr= combined2+4+i*rclen;
				for (size_t j=1;j< maxword;j=j+2){				
					if (!use_neg || j< maxword-1){						
						long long mask = static_cast<long long>(ydbl[j*m+i]);
						*((long long*) ptr) |= mask; 
						ptr+=4;
					}
				}
			}
			
			for (size_t i=0; i<m; i++) {
				LinBox::integer result, tmp;
				if (use_neg) {
					result = -ydbl[(maxword-1)*m+i];
					result <<= (maxword-1)*32;					
				}
				else
					result = 0;
				
				importWords(tmp, rclen, -1, 1, 0, 0, combined1+i*rclen);
				result += tmp;
				importWords(tmp, rclen, -1, 1, 0, 0, combined2+i*rclen);
				result += tmp;

				_D.init(y[i], result);
			}
			delete[] ydbl;
			delete[] combined1;
			delete[] combined2;

			return y;
		}// end of applyVspecial

	  
	private:
		Domain         _D;
		integer    _prime;
		Element _one,_zero;
		MatrixDomain<Domain> _MD;
	  

	};


	/**\brief blackbox apply optimizations
	\ingroup blackbox

	BlasApply and these MatrixApplyDomains are in blackbox/apply.h
	*/
	
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




	
	// special function to split an integer matrix in p-adic matrix store in an array of double
	template <class Domain, class IMatrix>
	void create_padic_chunk (const Domain           &D,
				 const IMatrix          &M,
				 double           *chunks, 
				 size_t         num_chunks);

	
	// Pascal's chunk-and-blas optimization:
	// instead of doing the ring multiplication A.digit, we write
	//     A = A0 + A1*2^16 + A2*2^32 + ... 
	// where A0, A1, ... (the 'chunks') are matrices of double
	// Then, we compute A.digit by multiplying each Ai.digit, and shifting and adding the results
	
	// To extend the method to when A has a negative entry, we add a final chunk that 
	// corresponds to a big negative number, using 2s complementing to write the negative number
	// as  -(big power of 2) + (sum of positive terms)
	// for example,   -0x000123456789 --> -(1 << 48) + (0xFFFE) << 32 + (0xDCBA) << 16 + (0x9877) (4 chunks)
	
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

			// ideally we would use chunks with chunk_size bits in them to make best use of the 
			// double representation, but for now it is only implemented to use a chunk size of 16	
			use_chunks = (chunk_size >= 16); 
			
			if (use_chunks){//usechunk
				// set maximum size of chunk to 16			
				chunk_size = 16;
				
				// compute the magnitude in bit of the matrix
				// check if we need negative representation
				LinBox::integer tmp=0, maxValue=0;
				size_t maxBitSize = 0;				
				use_neg = false;
				typename Matrix::ConstRawIterator it = _M.rawBegin();
				for (size_t i=0; i<_m*_n; i++, ++it) {
					_D.convert(tmp, *it);
					maxBitSize = max(maxBitSize, tmp.bitsize());
					maxValue   = (maxValue>tmp)? maxValue : tmp;
					use_neg |= (tmp < 0);
				}
				
				// compute the number of chunk
				if (maxValue*prime*_M.coldim() < integer("9007199254740992")){
					num_chunks=1;
					use_neg=false;
					//std::cout<<"using double apply\n";
				}
				else{
					num_chunks = (maxBitSize / chunk_size)+ (((maxBitSize % chunk_size) > 0)? 1:0);
					//std::cout<<"using padic double apply\n";
				}
					
				if (num_chunks ==1)
					use_neg= false;				
				
				if (use_neg) 
					num_chunks++; //the leading chunk will be negative
			
				int n2 = _m*_n;
				chunks = new double[n2*num_chunks];
				memset(chunks, 0, sizeof(double)*_m*_n*num_chunks);
			
				create_padic_chunk(_D, _M, chunks, num_chunks);

#ifdef DEBUG_CHUNK_SETUP			
				cout<<endl;
				cout << num_chunks << " chunks of "<< chunk_size << " bits each" << endl;
				if (!use_neg) cout << "not ";
				cout << "using negative leading chunk" << endl;
				cout << "Contents of chunks: " << endl;
				for (size_t i=0; i<num_chunks; i++) {
					cout << "chunk " << i << endl;
					for (size_t j=0; j<_m*_n; j++) {
						cout << static_cast<long long>(chunks[i*_m*_n+j]);
						if ((j+1)%_n) cout << ' '; else cout << endl;
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
#ifdef DEBUG_CHUNK_APPLYV
				cout << "x: ";
				for (size_t i=0; i<_n; i++) 
					cout << x[i] << ' ';
				cout << endl;
#endif

	
				if (num_chunks == 1) {
					double *ctd = new double[_m];
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
		
					unsigned char* combined = new unsigned char[rc*_n*rclen];
					memset(combined, 0, rc*_n*rclen);

					//order from major index to minor: combining index, component of sol'n, byte
		
					//compute a product (chunk times x) for each chunk

					double* ctd = new double[_n];
		
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
#ifdef DEBUG_CHUNK_APPLYV
							cout << "rcneg: " << result << endl;
#endif
						}
						else
							result = 0;
			
						for (int j=0; j<rc; j++) {
							unsigned char* thispos = combined + rclen*(j*_n+i);
							importWords(tmp, rclen, -1, 1, 0, 0, thispos);
							result += tmp;
#ifdef DEBUG_CHUNK_APPLYV
							cout << "rc[" << j << "," << i << "]:" << tmp << endl;
#endif
						}
#ifdef DEBUG_CHUNK_APPLYV
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


		IMatrix& applyM (IMatrix &Y, const IMatrix &X) const {

			linbox_check( _n == X.rowdim());
			linbox_check( _m == Y.rowdim());
			linbox_check( Y.coldim() == X.coldim()); 

			if (!use_chunks){
				_MD.mul (Y, _M, X);			
			}
			else{
				size_t _k= X.coldim();
				double* dX = new double[_n*_k];
				for (size_t i=0; i<_n; i++) 
					for(size_t j=0;j<_k;++j)
						_D.convert(dX[i*_k+j], X.getEntry(i,j));
					    
#ifdef DEBUG_CHUNK_APPLYM
				cout << "X: ";
				X.write(cout,_D);
				cout << endl;
#endif

	
				if (num_chunks == 1) {
					double *ctd = new double[_m*_k];
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
						    _m,_k,_n, 1,
						    chunks,_n, dX, _k, 0, ctd, _k);
					
					for (size_t i=0;i<_m;++i)
						for (size_t j=0;j<_k;++j)							
							_D.init(Y.refEntry(i,j),ctd[i*_k+j]);
					delete[] ctd;
					delete[] dX;
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
		
					unsigned char* combined = new unsigned char[rc*_m*_k*rclen];
					memset(combined, 0, rc*_m*_k*rclen);

					//order from major index to minor: combining index, component of sol'n, byte
		
					//compute a product (chunk times x) for each chunk
					double* ctd = new double[_m*_k];
		
					for (size_t i=0; i<num_chunks; i++) {
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							    _m,_k,_n, 1,
							    chunks+(_m*_n*i),_n, dX, _k, 0, ctd, _k);

						if (!use_neg || i<num_chunks-1)
							for (size_t j=0; j<_m*_k; j++) {
								// up to 53 bits will be ored-in, to be summed later
								unsigned char* bitDest = combined;
								bitDest += rclen*((i % rc)*_m*_k+j);
								long long mask = static_cast<long long>(ctd[j]);
								bitDest += 2*i;
								*((long long*) bitDest) |= mask; 
							}
					}
		
					delete[] dX;
					
					for (size_t i=0; i<_m*_k; i++) {
						LinBox::integer result, tmp;
						if (use_neg) {
							result = -ctd[i];
							result <<= (num_chunks-1)*16;
#ifdef DEBUG_CHUNK_APPLYM
							cout << "rcneg: " << result << endl;
#endif
						}
						else
							result = 0;
			
						for (int j=0; j<rc; j++) {
							unsigned char* thispos = combined + rclen*(j*_m*_k+i);
							importWords(tmp, rclen, -1, 1, 0, 0, thispos);
							result += tmp;
#ifdef DEBUG_CHUNK_APPLYM
							cout << "rc[" << j << "," << i << "]:" << tmp << endl;
#endif
						}
#ifdef DEBUG_CHUNK_APPLYM
						cout << "v2[" << i << "]:" << result  << endl;
#endif
						
						_D.init(*(Y.getWritePointer()+i), result);
					}
					delete[] combined;
					delete[] ctd;
				}
			}
			return Y;
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
	
	
	template<>
	template <class Domain>
	class MatrixApplyDomain<Domain, BlasBlackbox<Domain> > : 
		public BlasMatrixApplyDomain<Domain, BlasBlackbox<Domain> > {

	public:
		MatrixApplyDomain (const Domain &D, const  BlasBlackbox<Domain> &M)
			: BlasMatrixApplyDomain<Domain, BlasBlackbox<Domain> > (D,M) {}
		
	};



	/** \brief split an integer matrix into a padic chunk representation
	 *
	 */
	
	template <class Domain, class IMatrix>
	void create_padic_chunk (const Domain           &D,
				 const IMatrix          &M,
				 double            *chunks, 
				 size_t         num_chunks) {


		typename IMatrix::ConstRawIterator it= M.rawBegin();

		size_t m,n,mn;
		m  = M.rowdim();
		n  = M.coldim();
		mn = m*n;

		size_t tmpsize, tmpbitsize, j;
  

		if (num_chunks ==1)
			for (size_t i=0; i<mn; ++i, ++it) 
				D.convert(*(chunks+i), *it);  
		else
			for (size_t i=0; i<mn; ++i, ++it) {
				integer tmp;
				double* pdbl = chunks + i;
				D.convert(tmp, *it);
				tmpsize    = tmp.size();
				tmpbitsize = tmp.bitsize();
      
				if (tmp ==0)
					*pdbl=0;
				else  
					if (tmp > 0) {
					
		
#if LINBOX_SIZE_OF_LONG == 8
						// specialization for 64bits integer limbs
						for (j=0; j<tmpsize-1; j++) {
							*pdbl        =  tmp[j]        & 0xFFFF;
							*(pdbl+mn)   = (tmp[j] >> 16) & 0xFFFF;	
							*(pdbl+2*mn) = (tmp[j] >> 32) & 0xFFFF;
							*(pdbl+3*mn) = (tmp[j] >> 48) & 0xFFFF;
							pdbl      += 4*mn;
						}
						if ((tmpbitsize - j*64) > 0 ) {
							*pdbl = tmp[tmpsize-1]&0xFFFF; 
							pdbl+=mn;
						}
						if ((tmpbitsize - j*64) > 16 ) { 
							*pdbl = (tmp[tmpsize-1] >> 16)& 0xFFFF;
							pdbl+=mn;
						}
						if ((tmpbitsize - j*64) > 32 ) {
							*pdbl = (tmp[tmpsize-1] >> 32)& 0xFFFF;
							pdbl+=mn;
						}
						if ((tmpbitsize - j*64) > 48 ) 
							*pdbl = (tmp[tmpsize-1] >> 48)& 0xFFFF;
#else	     	    						
						// specialization for 32bits integer limbs	   	    
						for (j=0; j<tmpsize-1; j++) {
							*pdbl      = tmp[j] &  0xFFFF;
							*(pdbl+mn) = tmp[j] >> 16;									
							pdbl      += 2*mn;
						}
						if ((tmpbitsize - j*32) > 16 ) {
							*pdbl      = tmp[tmpsize-1] &  0xFFFF;
							*(pdbl+mn) = tmp[tmpsize-1] >> 16;									
						}
						else {
							*pdbl      = tmp[tmpsize-1] & 0xFFFF;								
						}						
#endif						
					}
					else {
						++tmp;
#if LINBOX_SIZE_OF_LONG == 8
						// specialization for 64bits integer limbs
						for (j=0; j<tmpsize-1; j++) {
							*pdbl        = 0xFFFF ^ ( tmp[j]        & 0xFFFF);
							*(pdbl+mn)   = 0xFFFF ^ ((tmp[j] >> 16) & 0xFFFF);
							*(pdbl+2*mn) = 0xFFFF ^ ((tmp[j] >> 32) & 0xFFFF);
							*(pdbl+3*mn) = 0xFFFF ^ ((tmp[j] >> 48) & 0xFFFF);
							pdbl        += 4*mn;
						}
						
						j=j<<2;
						if ((tmpbitsize - j*64) > 0 ) {
							*pdbl = 0xFFFF ^ (tmp[tmpsize-1]&0xFFFF); 
							pdbl+=mn;
							++j;
						}
						if ((tmpbitsize - j*64) > 16 ) { 
							*pdbl = 0xFFFF ^ ((tmp[tmpsize-1] >> 16)& 0xFFFF);
							pdbl+=mn;
							++j;
						}
						if ((tmpbitsize - j*64) > 32 ) {
							*pdbl = 0xFFFF ^ ((tmp[tmpsize-1] >> 32)& 0xFFFF);
								pdbl+=mn;
								++j;
						}
						if ((tmpbitsize - j*64) > 48 ) {
							*pdbl = 0xFFFF ^ ((tmp[tmpsize-1] >> 48)& 0xFFFF);
							++j;
						}
						
						for (; j<num_chunks-1; j++, pdbl += mn) 
							*pdbl      = 0xFFFF;
						*pdbl = 1; 
#else
						// specialization for 32bits integer limb
						for (j=0; j<tmpsize-1; j++) {
							*pdbl      = 0xFFFF ^ (tmp[j] & 0xFFFF);
							*(pdbl+mn) = 0xFFFF ^ (tmp[j] >> 16);
							pdbl      += 2*mn;							
						}
						j=j<<1;
						if ((tmpbitsize -j*32) > 16) {
							*pdbl      = 0xFFFF ^ (tmp[tmpsize-1] & 0xFFFF);
							*(pdbl+mn) = 0xFFFF ^ (tmp[tmpsize-1] >> 16);
							pdbl      += 2*mn;
							j+=2;
						}
						else {
							*pdbl      = 0xFFFF ^ (tmp[tmpsize-1] & 0xFFFF);
							pdbl      += mn;
							j+=1;
						}
						
						for (; j<num_chunks-1; j++, pdbl += mn) 
							*pdbl      = 0xFFFF;
						*pdbl = 1; 
#endif
					}
			}
		
	}



} // end of namespace LinBox		
#endif
