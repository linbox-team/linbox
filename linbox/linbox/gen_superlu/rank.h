// rank finds the rank of matrix A.
// The inputs L and U are obtained by the LU factorization
// Consider P_r * (A * P_c') = L * U
// rank(A) = rank(U)
// Implemented by A. Duran and D. Saunders
// University of Delaware

// #define EXTR

template <class Field>
int
rank(char *what, SuperMatrix<Field> *A1, SuperMatrix<Field> *A2, Field& F)
{
  int rnk;
  int previous;
  int current; 
  SCformat<Field> *A1store;
  NCformat<Field> *A2store;
  register int i, j, k, c, d, n, nsup;
  int nullCol;
  int nullRow;
  int rowNum;
  typename Field::Element *dp;
  typename Field::Element value;
  int *sup_to_col, *rowind, *rowind_colptr;


  A2store = (NCformat<Field> *) A2->Store;
  n = A2->ncol;
  int *visitedRow = (int *)malloc(n * sizeof(int)) ;
  int *visitedCol = (int *)malloc(n * sizeof(int)) ;
  for (i = 0; i < n; i++)
    visitedRow[i] = 0;

  for (i = 0; i < A2store->colptr[n]; ++i)
     {
     rowNum = A2store->rowind[i];
     if (visitedRow[rowNum] != 1)
        visitedRow[rowNum] = 1;
     }

  for (i = 0; i < n; i++)
    visitedCol[i] = 1;

  previous = 0;
  for (i = 1; i <= n; ++i)
     {
     current = A2store->colptr[i]; 
     if (current == previous)
        visitedCol[i-1] = 0;
     previous = current;
     }

#ifdef EXTR
  cout <<"visitedCol \n";
  for (i = 0; i < n; ++i)
     {
     cout << visitedCol[i] <<" ";
     }
     cout << endl;
  cout <<" n "<< n << "\n";
#endif

  A1store = (SCformat<Field> *) A1->Store;
  dp = (typename Field::Element *) A1store->nzval;
  sup_to_col = A1store->sup_to_col;
  rowind_colptr = A1store->rowind_colptr;
  rowind = A1store->rowind;

  for (k = 0; k <= A1store->nsuper+1; ++k) {
    c = sup_to_col[k];
    nsup = sup_to_col[k+1] - c;
    for (j = c; (j < c + nsup) && (j < n); ++j) {
      d = A1store->nzval_colptr[j];
      for (i = rowind_colptr[c]; i < rowind_colptr[c+1]; ++i) {
	if (rowind[i] <= j)
          {
          F.assign(value, dp[d++]);
          if (!F.isZero(value)) // modify after comparing n and nnz
             {
             visitedCol[j] = 1;
             visitedRow[rowind[i]] = 1;
#ifdef EXTR
             cout << " value "<< value <<"\n";
#endif
             }
          }
	else 
	  d++;
      }
    }
  }

  nullCol = 0;
  nullRow = 0;
  for (i = 0; i < n; ++i)
     {
     if (visitedCol[i] == 0) 
        nullCol++ ;
     if (visitedRow[i] == 0)
        nullRow++;
     }

  if (nullRow > nullCol)
     rnk = n - nullRow;
  else
     rnk = n - nullCol;

  cout <<"Rank of the matrix "<< what <<" is : " << rnk << endl;
}
