//----------------------------------------------------------------------
// LUdivine : LUP factorisation (Ibara & Moran) over a generic Field 
// rank version
//----------------------------------------------------------------------

template <class Field>
inline int LUdivine( Field& F,
	      	int M,  int N,		
		typename Field::element * A, int lda, int*jpiv)
{
  

  int MN = Min(M,N);
  int Nup,Ndown,R1,R2,ip,i,i0,j,k;
  
  if (MN > 1){ 

    Nup = MN>>1;
    Ndown =  M - Nup;
    typename Field::element * Ar, * Ac, * An;
    // Recursive call on NW block
    R1 = LUdivine(F, Nup, N, A, lda, jpiv);
    Ar = A + Nup*lda; // SW
    Ac = A + R1;      // NE
    An = Ar + R1;     // SE
    if (R1){
    
      flaswp(F,Ndown,Ar,lda,0,R1,jpiv,1);

#if DEBUG==2
      cerr<<"After first LUdivine recursive call and laswp"<<endl;
      cerr<<"R1="<<R1<<endl;
      write_field(F,cerr,A,M,N,lda);
#endif
      // Triangular block inversion and apply to SW
      Ftrsm(F, Ndown, R1, A, lda, Ar, lda);
      
      // Update of SE
      // An <- An - Ar*Ac
      FFFMMBLAS()( F, Ndown, N-R1, R1, FffmmblasNone, Ar, lda, Ac, lda, 
		   FffmmblasOne, An, lda, 0);
#if DEBUG==2
      cerr<<"After le FFFMMBLAS"<<endl;
      write_field(F,cerr,A,M,N,lda);
#endif
    }

    // Recursive call on SE

#if DEBUG==2
    cerr<<" Ndown="<<Ndown
    	<<" N-Nup="<<N-Nup<<" lda="<<lda<<endl;
#endif

    R2=LUdivine(F, Ndown, N-R1, An, lda, jpiv + R1);
    for ( i=R1;i!=MN;i++) jpiv[i] += R1;
    flaswp(F, Nup, A, lda, Nup, MN, jpiv, 1);
 

    // Non zero rows permutations
    if (R1<Nup){
      for (i=R1,i0=Nup; i0<(Nup+R2);i++,i0++){
	for (j=0; j<N; j++){
	  A[i*lda+j] = A[i0*lda+j];
	}
      }//FOR
    }//if (R1<Nup)
    
    return R1+R2;
  }    

   
  else if (MN==1){ // Last recursion level
    // Look for a non zero pivot
    ip=0;
    while (ip<N && iszero(*(A+ip))){ip++;}
    if (ip==N) // No non zero pivot found => row is zero
      return 0;
    *jpiv=ip;
    if (ip!=0){

#if DEBUG==2
      cerr<<"Permutation ip="<<ip<<endl;
#endif

      // Swap
      *A = *(A+ip);
      *(A+ip) = F.zero;
    }
    // Normalisation of the row
    for ( k=1; k<N; k++)
      F.divin(*(A+k), *A);
    return 1;
  }
  else //error
    return -1;

}
