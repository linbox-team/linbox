
	/***************************************************
	 **** Polynomial Matrix Multiplication over Z[x] ***
	 ***************************************************/
	template<>
	class PolynomialMatrixFFTMulDomain<PID_integer> {
	public:
		typedef PID_integer       IntField;
		typedef Modular<uint32_t> ModField;
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,ModField> MatrixP_F; // Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,IntField> MatrixP_I; // Polynomial matrix stored as a polynomial of matrix

	private: 
		const IntField     *_field;
	public:
		inline const IntField & field() const { return *_field; }


		PolynomialMatrixFFTMulDomain (const IntField &F) :
			_field(&F) {}
		
		template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
		void mul (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b) {
			//compute a bound on the entry of the input matrix a and b
			FFT_PROFILE_START;
			integer maxA,maxB,maxC,minA,minB,mA,mB;
			int signA,signB,signC;
			maxA=minA=a.get(0,0);
			maxB=minB=b.get(0,0);
			for (size_t i=0;i<a.size();i++)
				for (size_t j=0;j<a.rowdim()*a.coldim();j++){
					maxA=std::max(maxA,a.get(j,i));
					minA=std::min(minA,a.get(j,i));
				}
			for (size_t i=0;i<b.size();i++)
				for (size_t j=0;j<b.rowdim()*b.coldim();j++){
					maxB=std::max(maxB,b.get(j,i));
					minB=std::min(minB,b.get(j,i));
				}			
			signA= -1*(minA<0)+ (maxA>0); // 1 -> only >0, -1 -> only <0, 0 otherwise 
			signB= -1*(minB<0)+ (maxB>0); // 1 -> only >0, -1 -> only <0, 0 otherwise 
			signC= signA*signB;
			mA=max(abs(maxA),abs(minA));
			mB=max(abs(maxB),abs(minB));
							
			FFT_PROFILING(2,"max norm computation");
#ifdef CRTNAIVE
			mul_crtnaive(c,a,b,mA,mB,signC);
#else
			mul_crtla(c,a,b,mA,mB,signC);
#endif
			
		}

		template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
		void midproduct (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
				 bool smallLeft=true, size_t n0=0, size_t n1=0) {
			//compute a bound on the entry of the input matrix a and b
			FFT_PROFILE_START;
			integer maxA,maxB,maxC,minA,minB,mA,mB;
			int signA,signB,signC;
			maxA=minA=a.get(0,0);
			maxB=minB=b.get(0,0);
			for (size_t i=0;i<a.size();i++)
				for (size_t j=0;j<a.rowdim()*a.coldim();j++){
					maxA=std::max(maxA,a.get(j,i));
					minA=std::min(minA,a.get(j,i));
				}
			for (size_t i=0;i<b.size();i++)
				for (size_t j=0;j<b.rowdim()*b.coldim();j++){
					maxB=std::max(maxB,b.get(j,i));
					minB=std::min(minB,b.get(j,i));
				}			
			signA= -1*(minA<0)+ (maxA>0); // 1 -> only >0, -1 -> only <0, 0 otherwise 
			signB= -1*(minB<0)+ (maxB>0); // 1 -> only >0, -1 -> only <0, 0 otherwise 
			signC= signA*signB;
			mA=max(abs(maxA),abs(minA));
			mB=max(abs(maxB),abs(minB));

			FFT_PROFILING(2,"max norm computation");

			midproduct_crtla(c,a,b,mA,mB,signC,smallLeft, n0,n1);
		} 

		template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
		void mul_crtnaive(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b, const integer& maxA, const integer& maxB, const int signC) {
			linbox_check(a.coldim() == b.rowdim());
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t s= a.size()+b.size()-1;
			size_t lpts=0;
			size_t pts  = 1; while (pts < s) { pts= pts<<1; ++lpts; }
			size_t _k=k,lk=0;
			// compute bit size of feasible prime for FFLAS
			while ( _k ) {_k>>=1; ++lk;}    
			size_t prime_bitsize= (53-lk)>>1;
			integer primes_prod;			
			vector<long> primes;
			integer bound=maxA*maxB*a.coldim()*(std::max(a.size(),b.size())-1);
			if (signC==0)
				bound*=2; // matrix entries are positive and negative 
			primes=getFFTModuli(bound, prime_bitsize,primes_prod);
#ifdef FFT_PROFILE
			std::cout << "num of primes " << primes.size() << endl;
#endif
			// allocate fftprime fields
			vector<ModField> f_i (primes.size(), ModField(3));

			vector<typename PMatrix1::template rebind<ModField>::Other_t > c_i (primes.size());

			double tDecomp=0.,tMul=0.;
			for (size_t i = 0; i < primes.size(); ++i){
				FFT_PROFILE_START;
				f_i[i] = ModField(primes[i]);
				c_i[i] = typename PMatrix1::template rebind<ModField>::Other_t (f_i[i],m, n,s); 
				typename PMatrix2::template rebind<ModField>::Other_t a_i (f_i[i], m, k, a.size());
				typename PMatrix3::template rebind<ModField>::Other_t b_i (f_i[i], k, n, b.size());
				typename PMatrix2::template rebind<ModField>()(a_i,a); 
				typename PMatrix3::template rebind<ModField>()(b_i,b); 
				PolynomialMatrixFFTPrimeMulDomain fftdomain (f_i[i]);
				FFT_PROFILE_GET(tDecomp);
				fftdomain.mul(c_i[i], a_i, b_i);
				FFT_PROFILE_GET(tMul);
			}
			FFT_PROFILE(2,"k prime decomposition",tDecomp);
			FFT_PROFILE(2,"FFTprime mult",tMul);
			FFT_PROFILE_START;
			// reconstruct the solution modulo the original prime
			if (primes.size() < 2) {				
				typename PMatrix1::template rebind<ModField>::Other_t::template rebind<IntField>()(c,c_i[0]);
				integer hM=(primes_prod-1)>>1;
				for (size_t j = 0; j < m * n; j++)
					for (size_t i = 0; i < s; i++) {
						// get the correct result according to the sign of C
						if ((signC==0 && c.ref(j,i)>hM) || signC<0) 
							c.ref(j,i)-=primes_prod;
					}			
			} else {
				vector<integer> crt (primes.size());
				vector<ModField::Element> crt_inv (primes.size());
				ModField::Element tmp;
				for (size_t i=0;i<primes.size(); ++i){
					crt[i]=primes_prod/primes[i];
					f_i[i].init(tmp,crt[i]);
					f_i[i].inv(crt_inv[i], tmp);
				}

				// Naive integer reconstruction algorithm
				// No binary tree since nbrprimes is low
				integer res,acc;
				integer hM=(primes_prod-1)>>1;
				for (size_t j = 0; j < m * n; j++)
					for (size_t i = 0; i < s; i++) {
						acc= integer(0);
						for (size_t l=0;l<primes.size(); ++l){
							f_i[l].mul(tmp, c_i[l].get(j,i) ,crt_inv[l]);
							acc+= integer(tmp)*crt[l]; 
							if (acc > primes_prod)
								acc-= primes_prod;
						}
						
						c.ref(j,i)=acc;
						// get the correct result according to the sign of C
						if ((signC==0 && c.ref(j,i)>hM) || signC<0) 
							c.ref(j,i)-=primes_prod;
					}			
			}
			FFT_PROFILING(2,"k prime reconstruction");
		}

		template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
		void mul_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b, const integer& maxA, const integer& maxB, const int signC) {
			// (convert to MatrixP representation)
			MatrixP_I a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_I b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_I c2(field(),c.rowdim(),c.coldim(),c.size());
			mul_crtla(c2,a2,b2,maxA,maxB,signC);
			c.copy(c2,0,c.size()-1);			
		}

		
		void mul_crtla(MatrixP_I &c, const MatrixP_I &a, const MatrixP_I &b, const integer& maxA, const integer& maxB, const int signC) {
			FFT_PROFILE_START;
			linbox_check(a.coldim() == b.rowdim());
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t s= a.size()+b.size()-1;
			size_t lpts=0;
			size_t pts  = 1; while (pts < s) { pts= pts<<1; ++lpts; }
			size_t _k=k,lk=0;
			// compute bit size of feasible prime for FFLAS
			while ( _k ) {_k>>=1; ++lk;}    
			size_t prime_bitsize= (53-lk)>>1;
			integer primes_prod;			
			vector<long> primes;
			integer bound=maxA*maxB*a.coldim()*(std::max(a.size(),b.size())-1);
			if (signC==0)
				bound*=2; // matrix entries are positive and negative 
			primes=getFFTModuli(bound, prime_bitsize,primes_prod);
			size_t num_primes = primes.size();
			size_t log_crt    = primes_prod.bitsize();
			size_t log_beta   = 16;
#ifdef FFT_PROFILER	
			double tMul=0.,tCopy=0;;		
			if (FFT_PROF_LEVEL<3){	
			cout << "number of FFT primes :" << primes.size() << endl;
			cout << "feasible prime bitsize : "<<prime_bitsize<<endl;
			cout << "bitsize of the output: "<<bound.bitsize()
			     <<"( "<< primes_prod.bitsize()<<" )"<<endl;
			cout <<" +++++++++++++++++++++++++++++++"<<endl;					
		}
#endif		
			FFT_PROFILING(2,"init of CRT approach");
			// reduce t_a and t_b modulo each FFT primes
			size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
			double* t_a_mod= new double[n_ta*num_primes];
			double* t_b_mod= new double[n_tb*num_primes];
			reduce(primes, t_a_mod, a.getPointer(), n_ta, log_beta, maxA.bitsize()/log_beta + (maxA.bitsize()%log_beta?1:0));
			reduce(primes, t_b_mod, b.getPointer(), n_tb, log_beta, maxB.bitsize()/log_beta + (maxB.bitsize()%log_beta?1:0));
			FFT_PROFILING(2,"reduction mod pi of input matrices");

			vector<MatrixP_F*> c_i (num_primes);
			vector<ModField> f(num_primes);
			for (size_t l=0;l<num_primes;l++)
				f[l]=ModField(primes[l]);
#ifdef HAVE_OPENMP
#pragma omp parallel for shared(c_i,f, t_a_mod,t_b_mod) schedule(dynamic)
#endif
			for (size_t l=0;l<num_primes;l++){ 
				FFT_PROFILE_START;
				MatrixP_F a_i (f[l], m, k, pts); 
				MatrixP_F b_i (f[l], k, n, pts); 
				c_i[l] = new MatrixP_F(f[l], m, n, pts); 
				// copy reduced data
				for (size_t i=0;i<m*k;i++)
					for (size_t j=0;j<a.size();j++)
						a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];
				for (size_t i=0;i<k*n;i++)
					for (size_t j=0;j<b.size();j++)
						b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];	
				FFT_PROFILE_GET(tCopy);
				PolynomialMatrixFFTPrimeMulDomain fftdomain (f[l]);
				fftdomain.mul_fft(lpts, *c_i[l], a_i, b_i);
				FFT_PROFILE_GET(tMul);
			}
			delete[] t_a_mod;
			delete[] t_b_mod;
			FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
			FFT_PROFILE(2,"FFTprime multiplication",tMul);
			FFT_PROFILE_START;
			
			if (primes.size() < 2) {
				FFT_PROFILE_START;
				c.copy(*c_i[0],0,s-1);
			} else {		
				FFT_PROFILE_START;
				// construct contiguous storage for c_i
				double *t_c_mod;
				size_t n_tc=m*n*s;
				t_c_mod = new double[n_tc*num_primes];
				for (size_t l=0;l<num_primes;l++)
					for (size_t i=0;i<m*n;i++)
						for (size_t j=0;j<s;j++)
							t_c_mod[l*n_tc + (j+i*s)]= c_i[l]->get(i,j);
				FFT_PROFILING(2,"linearization of results mod pi");

				// reconstruct the result in C 
				reconstruct(c.getWritePointer(),t_c_mod,n_tc,primes, primes_prod, log_beta,log_crt/log_beta+(log_crt/log_beta?1:0));
				delete[] t_c_mod;
			
			}	
				// get the correct result according to the sign of C
				if (signC==0){ 
					integer hM=(primes_prod-1)>>1;
					for (size_t i=0;i<m*n;i++)
						for (size_t j=0;j<s;j++)
							if (c.ref(i,j)>hM)
								c.ref(i,j)-=primes_prod;
				}
				if (signC<0)
					for (size_t i=0;i<m*n;i++)
						for (size_t j=0;j<s;j++)
							c.ref(i,j)-=primes_prod;
				       		
				for (size_t i=0;i<num_primes;i++)
					delete c_i[i];

				FFT_PROFILING(2,"k prime reconstruction");
		}
		 

		template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
		void midproduct_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
				      const integer& maxA, const integer& maxB, const int signC, 
				      bool smallLeft=true, size_t n0=0, size_t n1=0) {
			// (convert to MatrixP representation)
			MatrixP_I a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_I b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_I c2(field(),c.rowdim(),c.coldim(),c.size());
			midproduct_crtla(c2,a2,b2,maxA,maxB,smallLeft,signC,n0,n1);
			c.copy(c2,0,c2.size()-1);			
		}

		void midproduct_crtla(MatrixP_I &c, const MatrixP_I &a, const MatrixP_I &b,
				      const integer& maxA, const integer& maxB, const int signC,
				      bool smallLeft=true, size_t n0=0, size_t n1=0) {
			FFT_PROFILE_START;
			linbox_check(a.coldim() == b.rowdim());
			size_t m = a.rowdim();
			size_t k = a.coldim();
			size_t n = b.coldim();
			size_t hdeg = (n0==0?c.size():n0);
			size_t deg  = (n1==0?2*hdeg:n1);
			linbox_check(c.size()>=deg-hdeg);

			if (smallLeft){
				linbox_check(b.size()<hdeg+deg);
			}
			else
				linbox_check(a.size()<hdeg+deg);

			//linbox_check(2*c.size()-1 == b.size());
			//size_t deg= b.size()+1;
			//size_t hdeg= deg/2;
			size_t lpts=0;
			size_t pts  = 1; while (pts < deg) { pts= pts<<1; ++lpts; }
			size_t _k=k,lk=0;
			// compute bit size of feasible prime for FFLAS
			while ( _k ) {_k>>=1; ++lk;}    
			size_t prime_bitsize= (53-lk)>>1;
			integer primes_prod;			
			vector<long> primes;
			integer bound;
			if (smallLeft)
				bound=maxA*maxB*k* a.size();
			else
				bound=maxA*maxB*k* b.size();
			if (signC==0)
				bound*=2;
			primes=getFFTModuli(bound, prime_bitsize,primes_prod);
			size_t num_primes = primes.size();
			size_t log_crt    = primes_prod.bitsize();
			size_t log_beta   = 16;
#ifdef FFT_PROFILER	
			double tMul=0.,tCopy=0;;
			if (FFT_PROF_LEVEL<3){	
				cout << "number of FFT primes :" << primes.size() << endl;
				cout << "feasible prime bitsize : "<<prime_bitsize<<endl;
				cout << "bitsize of the output: "<<bound.bitsize()
				     <<"( "<< primes_prod.bitsize()<<" )"<<endl;
				cout <<" +++++++++++++++++++++++++++++++"<<endl;		
				}
#endif		
			FFT_PROFILING(2,"init of CRT approach");
			// reduce t_a and t_b modulo each FFT primes
			size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
			double* t_a_mod= new double[n_ta*num_primes];
			double* t_b_mod= new double[n_tb*num_primes];

			reduce(primes, t_a_mod, a.getPointer(), n_ta, log_beta, maxA.bitsize()/log_beta + (maxA.bitsize()%log_beta?1:0));
			reduce(primes, t_b_mod, b.getPointer(), n_tb, log_beta, maxB.bitsize()/log_beta + (maxB.bitsize()%log_beta?1:0));
			FFT_PROFILING(2,"reduction mod pi of input matrices");

			vector<MatrixP_F> c_i (primes.size());

			for (size_t l=0;l<num_primes;l++){ 
				FFT_PROFILE_START;
				ModField f(primes[l]);
				MatrixP_F a_i (f, m, k, pts); 
				MatrixP_F b_i (f, k, n, pts); 
				c_i[l] = MatrixP_F(f, m, n, pts); 
				// copy reduced data and reversed when necessary according to midproduct algo
 				for (size_t i=0;i<m*k;i++)
					for (size_t j=0;j<a.size();j++)
						if (smallLeft)  
							a_i.ref(i,hdeg-1-j)=t_a_mod[l*n_ta+j+i*a.size()];
						else
							a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];					
				for (size_t i=0;i<k*n;i++)
					for (size_t j=0;j<b.size();j++)	
						if (smallLeft)
							b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];	
						else 
							b_i.ref(i,hdeg-1-j)=t_b_mod[l*n_tb+j+i*b.size()];
				FFT_PROFILE_GET(tCopy);
				PolynomialMatrixFFTPrimeMulDomain fftdomain (f);
				fftdomain.midproduct_fft(lpts, c_i[l], a_i, b_i, smallLeft);
				FFT_PROFILE_GET(tMul);
			}
			delete[] t_a_mod;
			delete[] t_b_mod;
			FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
			FFT_PROFILE(2,"FFTprime multiplication",tMul);
			FFT_PROFILE_START;

			if (primes.size() < 2) {
				FFT_PROFILE_START;
				c.copy(c_i[0],0,c.size()-1);
			} else {		
				FFT_PROFILE_START;
				// construct contiguous storage for c_i
				double *t_c_mod;
				size_t n_tc=m*n*c.size()-1;
				t_c_mod = new double[n_tc*num_primes];
				for (size_t l=0;l<num_primes;l++)
					for (size_t i=0;i<m*n;i++)
						for (size_t j=0;j<c.size();j++)
							t_c_mod[l*n_tc + (j+i*c.size())]= c_i[l].get(i,j);
				FFT_PROFILING(2,"linearization of results mod pi");

				// reconstruct the result in C 
				reconstruct(c.getWritePointer(),t_c_mod,n_tc,primes, primes_prod, log_beta,log_crt/log_beta+(log_crt/log_beta?1:0));
				
				delete[] t_c_mod;
			}

			// get the correct result according to the sign of C			
			if (signC==0){ 
				integer hM=(primes_prod-1)>>1;
				for (size_t i=0;i<m*n;i++)
					for (size_t j=0;j<c.size();j++)
						if (c.ref(i,j)>hM)
							c.ref(i,j)-=primes_prod;
			}
			
			if (signC<0)
				for (size_t i=0;i<m*n;i++)
					for (size_t j=0;j<c.size();j++)
						c.ref(i,j)-=primes_prod; 
			
			FFT_PROFILING(2,"k prime reconstruction");
		} 
	};


	/***************************************************************************
	 **** Polynomial Matrix Multiplication over Fp[x], with p multiprecision ***
 	 ***************************************************************************/
	template <>
	class PolynomialMatrixFFTMulDomain<Modular<integer> > {
	public:
		typedef Modular<integer>              Field;
		typedef typename Field::Element     Element;
		typedef PID_integer                IntField;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP_F;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,IntField> MatrixP_I;
 
	private:
		const Field            *_field;  // Read only
		integer                     _p;
		
	public:
		inline const Field & field() const { return *_field; }
		
		PolynomialMatrixFFTMulDomain(const Field &F) : _field(&F) {
			field().cardinality(_p);
		}
		
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b) {
			MatrixP_F a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_F b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_F c2(field(),c.rowdim(),c.coldim(),c.size());
			mul(c2,a2,b2);
			c.copy(c2,0,c.size()-1);			
		}
		
		void mul (MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b) {
			PID_integer Z;
			PolynomialMatrixFFTMulDomain<PID_integer> Zmul(Z);
			const MatrixP_I* a2 = reinterpret_cast<const MatrixP_I*>(&a);
			const MatrixP_I* b2 = reinterpret_cast<const MatrixP_I*>(&b);
			MatrixP_I* c2       = reinterpret_cast<MatrixP_I*>(&c);
			Zmul.mul_crtla(*c2,*a2,*b2,_p-1,_p-1,1);
			// reduce the result mod p
			for (size_t i=0;i<c2->rowdim()*c2->coldim();i++)
				for (size_t j=0;j<c2->size();j++)
					c2->ref(i,j)%=_p;			

		}

		template<typename Matrix1, typename Matrix2, typename Matrix3>
		void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b, 
				 bool smallLeft=true, size_t n0=0, size_t n1=0) {
			
			MatrixP_F a2(field(),a.rowdim(),a.coldim(),a.size());
			MatrixP_F b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			MatrixP_F c2(field(),c.rowdim(),c.coldim(),c.size());
			midproduct(c2,a2,b2,smallLeft,n0,n1);
			c.copy(c2,0,c.size()-1);			
		}
		
		void midproduct (MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b,
				 bool smallLeft=true, size_t n0=0, size_t n1=0) {
			PID_integer Z;
			PolynomialMatrixFFTMulDomain<PID_integer> Zmul(Z);
			const MatrixP_I* a2 = reinterpret_cast<const MatrixP_I*>(&a);
			const MatrixP_I* b2 = reinterpret_cast<const MatrixP_I*>(&b);
			MatrixP_I* c2       = reinterpret_cast<MatrixP_I*>(&c);
			Zmul.midproduct_crtla(*c2,*a2,*b2,_p-1,_p-1,1,smallLeft,n0,n1);
			// reduce the result mod p
			for (size_t i=0;i<c2->rowdim()*c2->coldim();i++)
				for (size_t j=0;j<c2->size();j++)
					c2->ref(i,j)%=_p;			
		}
	};
