#include <stdlib.h>
// determinant finds the determinant of matrix A.
// The inputs L, vec_perm_r and vec_perm_c are obtained by the LU factorization
// Consider P_r * (A * P_c') = L * U
// So, det(A) = det(P_r) * det(L) * det(P_c')
// Implemented by A. Duran and D. Saunders
// University of Delaware

int swapNumber(int n, int *vec)
{
  int total_swap;
  int cycle_len;
  int *visited = (int *)malloc(n * sizeof(int)) ;
  int done = 0;
  int val_pre;
  int val_cur;
  int i, j;

  for (i = 0; i < n; i++)
    visited[i] = 0;
      
  // Finds the number of swap
  total_swap = 0;
  cycle_len = 0;
  i = 0;
          
  while (! done)
    {
      j = i;
      do
        {
          visited[j] = 1;
	  val_cur = vec[j];
	  j = val_cur;
	  cycle_len++;
	}
      while (i != val_cur);

      total_swap = total_swap + cycle_len - 1;
      do
        {
          i++;
          if (i >= n)
            done = 1;
        }
      while (visited[i] == 1);
      cout << "cycle_len :"<< cycle_len << endl;
      cycle_len = 0;
    }   
    return total_swap;
}


template <class Field>
void
determinant(char *what, SuperMatrix<Field> *A, int n, int *vec_perm_r, int *vec_perm_c, int rank, Field& F)
{
  typename Field::Element det;
  int lsb_det_perm_r;
  int lsb_det_perm_c;
  SCformat<Field>     *Astore;
  register int i, j, k, c, d, nsup;
  typename Field::Element *dp;
  int *sup_to_col, *rowind, *rowind_colptr;

  int nswap1;
  int nswap2;
  
  if (rank < n)
    F.assign(det, 0);
  else
    {
      det = 1;
      Astore = (SCformat<Field> *) A->Store;
      dp = (typename Field::Element *) Astore->nzval;
      sup_to_col = Astore->sup_to_col;
      rowind_colptr = Astore->rowind_colptr;
      rowind = Astore->rowind;
      
      for (k = 0; k <= Astore->nsuper+1; ++k) {
	c = sup_to_col[k];
	nsup = sup_to_col[k+1] - c;
	for (j = c; j < c + nsup; ++j) {
	  d = Astore->nzval_colptr[j];
	  for (i = rowind_colptr[c]; i < rowind_colptr[c+1]; ++i) {
	    if (rowind[i] == j)
	      F.mulin(det, dp[d++]);
	    else 
	      d++;
	  }
	}
      }
      
      nswap1 = swapNumber(n, vec_perm_r);
      lsb_det_perm_r = nswap1 & 1; // lsb of nswap1
      cout << "nswap1 :"<< nswap1 << endl;
      cout << "lsb_det_perm_r :"<< lsb_det_perm_r << endl << endl;
      
      nswap2 = swapNumber(n, vec_perm_c);
      cout << "nswap2 :"<< nswap2 << endl;
      lsb_det_perm_c = nswap2 & 1; // lsb of nswap2
      cout << "lsb_det_perm_c :"<< lsb_det_perm_c << endl;  
      
      if (lsb_det_perm_r != lsb_det_perm_c)
	F.negin(det);
    }    
  cout <<"determinant of the matrix " << what << " is : " << det << endl;
}
  
