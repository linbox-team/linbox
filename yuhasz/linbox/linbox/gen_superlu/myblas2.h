

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/*
 * File name:		myblas2.c
 * Purpose:
 *     Level 2 BLAS operations: solves and matvec, written in C.
 * Note:
 *     This is only used when the system lacks an efficient BLAS library.
 */

/*
 * Solves a dense UNIT lower triangular system. The unit lower 
 * triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
 * The solution will be returned in the rhs vector.
 */


template <class Field>
void lsolve ( int ldm, int ncol, typename Field::Element *M, typename Field::Element *rhs, Field& F )
{
    int k;
    typename Field::Element x0, x1, x2, x3, x4, x5, x6, x7;
    typename Field::Element temp_mul; // A.Duran 8/13/2002
    typename Field::Element temp_sub; // A.Duran 8/13/2002
    typename Field::Element *M0;
    register typename Field::Element *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
    register int firstcol = 0;

    M0 = &M[0];

#ifdef DEBUG
      cout << "At lsolve visited \n"; 
#endif

    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;
      Mki2 = Mki1 + ldm + 1;
      Mki3 = Mki2 + ldm + 1;
      Mki4 = Mki3 + ldm + 1;
      Mki5 = Mki4 + ldm + 1;
      Mki6 = Mki5 + ldm + 1;
      Mki7 = Mki6 + ldm + 1;

      x0 = rhs[firstcol];
      // x1 = rhs[firstcol+1] - x0 * *Mki0++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(x1, rhs[firstcol+1], temp_mul); // A.Duran 8/13/2002

      // x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+2], temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x1, *Mki1++); // A.Duran 8/13/2002
      F.sub(x2, temp_sub, temp_mul); // A.Duran 8/13/2002

      // x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+3], temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x1, *Mki1++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x2, *Mki2++); // A.Duran 8/13/2002
      F.sub(x3, temp_sub, temp_mul); // A.Duran 8/13/2002

      // x4 = rhs[firstcol+4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
      //	                   - x3 * *Mki3++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+4], temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x1, *Mki1++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x2, *Mki2++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x3, *Mki3++); // A.Duran 8/13/2002
      F.sub(x4, temp_sub, temp_mul); // A.Duran 8/13/2002

      // x5 = rhs[firstcol+5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
      //                  - x3 * *Mki3++ - x4 * *Mki4++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+5], temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x1, *Mki1++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x2, *Mki2++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x3, *Mki3++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x4, *Mki4++); // A.Duran 8/13/2002
      F.sub(x5, temp_sub, temp_mul); // A.Duran 8/13/2002

      // x6 = rhs[firstcol+6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
      //                  - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+6], temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x1, *Mki1++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x2, *Mki2++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x3, *Mki3++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x4, *Mki4++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x5, *Mki5++); // A.Duran 8/13/2002
      F.sub(x6, temp_sub, temp_mul); // A.Duran 8/13/2002

      // x7 = rhs[firstcol+7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
      //                  - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
      //		   - x6 * *Mki6++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+7], temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x1, *Mki1++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x2, *Mki2++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x3, *Mki3++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x4, *Mki4++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x5, *Mki5++); // A.Duran 8/13/2002
      F.subin(temp_sub, temp_mul); // A.Duran 8/13/2002
      F.mul(temp_mul, x6, *Mki6++); // A.Duran 8/13/2002
      F.sub(x7, temp_sub, temp_mul); // A.Duran 8/13/2002

      rhs[++firstcol] = x1;
      rhs[++firstcol] = x2;
      rhs[++firstcol] = x3;
      rhs[++firstcol] = x4;
      rhs[++firstcol] = x5;
      rhs[++firstcol] = x6;
      rhs[++firstcol] = x7;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++) {
	// rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
	// - x2 * *Mki2++ - x3 * *Mki3++
	// - x4 * *Mki4++ - x5 * *Mki5++
	// - x6 * *Mki6++ - x7 * *Mki7++;
	F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
	F.sub(temp_sub, rhs[k], temp_mul);
	F.mul(temp_mul, x1, *Mki1++);
	F.subin(temp_sub, temp_mul);
	F.mul(temp_mul, x2, *Mki2++); 
	F.subin(temp_sub, temp_mul); 
	F.mul(temp_mul, x3, *Mki3++); 
	F.subin(temp_sub, temp_mul); 
	F.mul(temp_mul, x4, *Mki4++);
	F.subin(temp_sub, temp_mul); 
	F.mul(temp_mul, x5, *Mki5++);
	F.subin(temp_sub, temp_mul); 
	F.mul(temp_mul, x6, *Mki6++);
	F.subin(temp_sub, temp_mul); 
	F.mul(temp_mul, x7, *Mki7++);
	F.sub(rhs[k], temp_sub, temp_mul);
#ifdef DEBUG
      cout << "At lsolve rhs[k] "<< rhs[k] <<"\n"; 
#endif 
      }
      M0 += 8 * ldm + 8;
    }

    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;
      Mki2 = Mki1 + ldm + 1;
      Mki3 = Mki2 + ldm + 1;

      x0 = rhs[firstcol];
      // x1 = rhs[firstcol+1] - x0 * *Mki0++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(x1, rhs[firstcol+1], temp_mul);

      // x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+2], temp_mul);
      F.mul(temp_mul, x1, *Mki1++);
      F.sub(x2, temp_sub, temp_mul);

      // x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(temp_sub, rhs[firstcol+3], temp_mul);
      F.mul(temp_mul, x1, *Mki1++);
      F.subin(temp_sub, temp_mul);
      F.mul(temp_mul, x2, *Mki2++);
      F.sub(x3, temp_sub, temp_mul);

      rhs[++firstcol] = x1;
      rhs[++firstcol] = x2;
      rhs[++firstcol] = x3;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++) {
	// rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
	//             - x2 * *Mki2++ - x3 * *Mki3++;
	F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
	F.sub(temp_sub, rhs[k], temp_mul);
	F.mul(temp_mul, x1, *Mki1++);
	F.subin(temp_sub, temp_mul);
	F.mul(temp_mul, x2, *Mki2++); 
	F.subin(temp_sub, temp_mul); 
	F.mul(temp_mul, x3, *Mki3++); 
	F.sub(rhs[k], temp_sub, temp_mul);
#ifdef DEBUG
      cout << "At lsolve rhs[k] "<< rhs[k] <<"\n"; 
#endif 
      }
      M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
      Mki0 = M0 + 1;
      Mki1 = Mki0 + ldm + 1;

      x0 = rhs[firstcol];
      // x1 = rhs[firstcol+1] - x0 * *Mki0++;
      F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
      F.sub(x1, rhs[firstcol+1], temp_mul);

      rhs[++firstcol] = x1;
      ++firstcol;
    
      for (k = firstcol; k < ncol; k++) {
	// rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
	F.mul(temp_mul, x0, *Mki0++); // A.Duran 8/13/2002
	F.sub(temp_sub, rhs[k], temp_mul);
	F.mul(temp_mul, x1, *Mki1++);
	F.sub(rhs[k], temp_sub, temp_mul);
#ifdef DEBUG
      cout << "At lsolve rhs[k] "<< rhs[k] <<"\n"; 
#endif 
     }
 
    }
    
}

/*
 * Solves a dense upper triangular system. The upper triangular matrix is
 * stored in a 2-dim array M(1:ldm,1:ncol). The solution will be returned
 * in the rhs vector.
 */


template <class Field>
void
usolve ( int ldm, int ncol, typename Field::Element *M, typename Field::Element *rhs, Field& F )
     /* int ldm;	 in */
     /* int ncol;	 in */
     /* double *M;	 in */
     /* double *rhs;	 modified */
{
  typename Field::Element xj;
  typename Field::Element temp_mul; // A.Duran 8/14/2002
  int jcol, j, irow;
  
#ifdef DEBUG
      cout << "At usolve visited \n"; 
#endif

  jcol = ncol - 1;
  
  for (j = 0; j < ncol; j++) {
    
    // xj = rhs[jcol] / M[jcol + jcol*ldm]; 		/* M(jcol, jcol) */
    F.div(xj, rhs[jcol], M[jcol + jcol*ldm]); // A.Duran 8/14/2002

#ifdef DEBUG
    cout << "At usolve xj" << xj <<"\n";
#endif

    rhs[jcol] = xj;
    
    for (irow = 0; irow < jcol; irow++) {
      // rhs[irow] -= xj * M[irow + jcol*ldm];	/* M(irow, jcol) */
      F.mul(temp_mul, xj, M[irow + jcol*ldm]); // A.Duran 8/14/2002
      F.subin(rhs[irow], temp_mul);
    }
    jcol--;
    
  }
}


/*
 * Performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.
 * The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
 */
template <class Field>
void matvec ( int ldm, int nrow, int ncol, typename Field::Element *M, typename Field::Element *vec, typename Field::Element *Mxvec, Field& F )
     
     /* int ldm;	 in -- leading dimension of M */
     /* int nrow;	 in */ 
     /* int ncol;	 in */
     /* double *M;	 in */
     /* double *vec;	 in */
     /* double *Mxvec;	 in/out */
     
{
  typename Field::Element vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
  typename Field::Element *M0;
  typename Field::Element temp_add; // A.Duran 8/14/2002
  typename Field::Element temp_mul; // A.Duran 8/14/2002
  
  register typename Field::Element *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
  register int firstcol = 0;
  int k;
 
#ifdef DEBUG
      cout << "At matvec visited \n"; 
#endif
 
  M0 = &M[0];
  while ( firstcol < ncol - 7 ) {	/* Do 8 columns */
    
    Mki0 = M0;
    Mki1 = Mki0 + ldm;
    Mki2 = Mki1 + ldm;
    Mki3 = Mki2 + ldm;
    Mki4 = Mki3 + ldm;
    Mki5 = Mki4 + ldm;
    Mki6 = Mki5 + ldm;
    Mki7 = Mki6 + ldm;
    
    vi0 = vec[firstcol++];
    vi1 = vec[firstcol++];
    vi2 = vec[firstcol++];
    vi3 = vec[firstcol++];	
    vi4 = vec[firstcol++];
    vi5 = vec[firstcol++];
    vi6 = vec[firstcol++];
    vi7 = vec[firstcol++];	
    
    for (k = 0; k < nrow; k++) {
      // Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
      // + vi2 * *Mki2++ + vi3 * *Mki3++ 
      // + vi4 * *Mki4++ + vi5 * *Mki5++
      // + vi6 * *Mki6++ + vi7 * *Mki7++;
      F.mul(temp_mul, vi0, *Mki0++); // A.Duran 8/14/2002
      F.add(temp_add, Mxvec[k], temp_mul);
      F.mul(temp_mul, vi1, *Mki1++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi2, *Mki2++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi3, *Mki3++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi4, *Mki4++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi5, *Mki5++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi6, *Mki6++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi7, *Mki7++);
      F.add(Mxvec[k], temp_add, temp_mul);
    }
    M0 += 8 * ldm;
  }
  
  while ( firstcol < ncol - 3 ) {	/* Do 4 columns */
    
    Mki0 = M0;
    Mki1 = Mki0 + ldm;
    Mki2 = Mki1 + ldm;
    Mki3 = Mki2 + ldm;
    
    vi0 = vec[firstcol++];
    vi1 = vec[firstcol++];
    vi2 = vec[firstcol++];
    vi3 = vec[firstcol++];	
    for (k = 0; k < nrow; k++) {
      // Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
      // + vi2 * *Mki2++ + vi3 * *Mki3++ ;
      F.mul(temp_mul, vi0, *Mki0++); // A.Duran 8/14/2002
      F.add(temp_add, Mxvec[k], temp_mul);
      F.mul(temp_mul, vi1, *Mki1++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi2, *Mki2++);
      F.addin(temp_add, temp_mul);
      F.mul(temp_mul, vi3, *Mki3++);
      F.add(Mxvec[k], temp_add, temp_mul);
    }
    M0 += 4 * ldm;
  }
  
  while ( firstcol < ncol ) {		/* Do 1 column */
    
    Mki0 = M0;
    vi0 = vec[firstcol++];
    for (k = 0; k < nrow; k++) {
      // Mxvec[k] += vi0 * *Mki0++;
      F.mul(temp_mul, vi0, *Mki0++); // A.Duran 8/14/2002
      F.addin(Mxvec[k], temp_mul);
    }
    M0 += ldm;
  }
  
}

