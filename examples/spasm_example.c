/* indent -nfbs -nip -npsl -di0 rank_hybrid.c */
#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include <sys/time.h>
#include <stdbool.h>
#include <cstdint>
#include <stdlib.h>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif


#include <givaro/modular.h>
typedef Givaro::Modular<int32_t> Field;

#include <linbox/matrix/sparse-matrix.h>
typedef LinBox::SparseMatrix<Field, LinBox::SparseMatrixFormat::CSR > spasm;
 

/* --- primary SpaSM routines and data structures --- */

typedef int spasm_GFp;

// typedef struct {                /* matrix in compressed-sparse row format */
//   int nzmax;                    /* maximum number of entries */
//   int n;                        /* number of rows */
//   int m;                        /* number of columns */
//   int *p;                       /* row pointers (size n+1) */
//   int *j;                       /* column indices, size nzmax */
//   spasm_GFp *x;                 /* numerical values, size nzmax (optional) */
//   int prime;
// }      spasm;

typedef struct {                /* matrix in triplet form */
  int nzmax;                    /* maximum number of entries */
  int nz;                       /* # entries */
  int n;                        /* number of rows */
  int m;                        /* number of columns */
  int *i;                       /* row indices, size nzmax */
  int *j;                       /* column indices (size nzmax) */
  spasm_GFp *x;                 /* numerical values, size nzmax (optional) */
  int prime;
}      spasm_triplet;



typedef struct {                /* a PLUQ factorisation */
  spasm *L;
  spasm *U;
  int *qinv;                    /* the inverse Q is stored */
  int *p;
}      spasm_lu;

typedef struct {                /* a dense LU factorization */
  int n;                        /* number of rows */
  int m;                        /* number of columns */
  int prime; 
  int *p;                       /* positions of pivots in allocated rows */
  spasm_GFp **x;                /* pointers to the rows */
}      spasm_dense_lu;


#define SPASM_IDENTITY_PERMUTATION NULL
#define SPASM_IGNORE NULL
#define SPASM_IGNORE_VALUES 0
#define SPASM_WITH_NUMERICAL_VALUES 1
#define SPASM_KEEP_L 1
#define SPASM_DISCARD_L 0
#define SPASM_SUCCESS 0
#define SPASM_NO_SOLUTION 1



/* utilities */
static inline int spasm_max(int a, int b) {
  return (a > b) ? a : b;
}

static inline int spasm_min(int a, int b) {
  return (a < b) ? a : b;
}

static inline void spasm_lbswap(spasm *A, int i, int j) {
  int x = A->getColid(i);
  A->setColid(i, A->getColid(j));
  A->setColid(j, x);
  if(A->getData().size() > 0){
    x = A->getData(i);
    A->setData(i, A->getData(j));
    A->setData(j, x);
  }
}

static inline void spasm_swap(int *a, int i, int j) {
  int x = a[i];
  a[i] = a[j];
  a[j] = x;
}

//static inline int spasm_row_weight(const spasm * A, int i) {
static inline int spasm_row_weight(spasm * A, int i) {
//   int *Ap;
//   Ap = A->p;
//   return Ap[i + 1] - Ap[i];
    return A->rowLength(i);
}

double spasm_wtime() {
	struct timeval ts;

	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}


int spasm_nnz(const spasm * A) {
	assert(A != NULL);

	return A->getStart(A->rowdim());
//     return A->p[A->n];
}

/* return a string representing n in 4 bytes */
void spasm_human_format(int64_t n, char *target) {
	if (n < 1000) {
		sprintf(target, "%lld", n);
		return;
	}
	if (n < 1000000) {
		sprintf(target, "%.1fk", n / 1e3);
		return;
	}
	if (n < 1000000000) {
		sprintf(target, "%.1fm", n / 1e6);
		return;
	}
	if (n < 1000000000000ll) {
		sprintf(target, "%.1fg", n / 1e9);
		return;
	}
	if (n < 1000000000000000ll) {
		sprintf(target, "%.1ft", n / 1e12);
		return;
	}
}

void *spasm_malloc(size_t size) {
  void *x = std::malloc(size);
  if (x == NULL) {
    perror("malloc failed");
    exit(1);
  }
	
  return x;
}

void *spasm_calloc(size_t count, size_t size) {
  void *x = std::calloc(count, size);
  if (x == NULL) {
    perror("calloc failed");
    exit(1);
  }
  return x;
}

void *spasm_realloc(void *ptr, size_t size) {
  void *x = std::realloc(ptr, size);
  if (ptr != NULL && x == NULL && size != 0) {
    perror("realloc failed");
    exit(1);
  }
  return x;
}


/* allocate a sparse matrix (compressed-row form) */ 

spasm *spasm_csr_alloc(int n, int m, int nzmax, int prime, int with_values) {
  if (prime > 46337) {
    prime = 46337;
    fprintf(stderr, "WARNING: modulus has been set to 46337.\n");
  }
  
  // 	A = spasm_malloc(sizeof(spasm));	/* allocate the cs struct */
  // 	A->m = m;		/* define dimensions and nzmax */
  // 	A->n = n;
  // 	A->nzmax = nzmax;
  // 	A->prime = prime;
  // 	A->p = spasm_malloc((n + 1) * sizeof(int));
  // 	A->j = spasm_malloc(nzmax * sizeof(int));
  // 	A->x = (with_values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL);
  
  Field Fp(prime);
  spasm *A = new spasm(Fp,n,m,nzmax);
  return A;

}


/* allocate a sparse matrix (triplet form) */
spasm_triplet *spasm_triplet_alloc(int n, int m, int nzmax, int prime, int with_values) {
	spasm_triplet *A;

	A = (spasm_triplet*)spasm_malloc(sizeof(spasm_triplet));
	A->m = m;
	A->n = n;
	A->nzmax = nzmax;
	A->prime = prime;
	A->nz = 0;
	A->i = (int*)spasm_malloc(nzmax * sizeof(int));
	A->j = (int*)spasm_malloc(nzmax * sizeof(int));
	A->x = (spasm_GFp*)(with_values ? spasm_malloc(nzmax * sizeof(spasm_GFp)) : NULL);
	return A;
}



/*
 * change the max # of entries in a sparse matrix. If nzmax < 0, then the
 * matrix is trimmed to its current nnz.
 */
void spasm_csr_realloc(spasm * A, int nzmax) {
	assert(A != NULL);

	if (nzmax < 0) {
		nzmax = spasm_nnz(A);
	}
//     A->refColid().resize(nzmax);
//     A->refData().resize(nzmax);
//     A->setSize( static_cast<size_t>(nzmax) );
    A->resize(nzmax);

// 	A->j = spasm_realloc(A->j, nzmax * sizeof(int));
// 	if (A->x != NULL) {
// 		A->x = spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
// 	}
// 	A->nzmax = nzmax;
}


/*
 * change the max # of entries in a sparse matrix. If nzmax < 0, then the
 * matrix is trimmed to its current nnz.
 */
void spasm_triplet_realloc(spasm_triplet * A, int nzmax) {
	assert(A != NULL);

	if (nzmax < 0) {
		nzmax = A->nz;
	}
	A->i = (int*)spasm_realloc(A->i, nzmax * sizeof(int));
	A->j = (int*)spasm_realloc(A->j, nzmax * sizeof(int));
	if (A->x != NULL) {
	  A->x = (int*)spasm_realloc(A->x, nzmax * sizeof(spasm_GFp));
	}
	A->nzmax = nzmax;
}


/* free a sparse matrix */
void spasm_csr_free(spasm * A) {
	if (A == NULL) {
		return;
	}
// 	free(A->p);
// 	free(A->j);
// 	free(A->x);		/* trick : free does nothing on NULL pointer */
// 	free(A);
    delete A;
}

void spasm_triplet_free(spasm_triplet * A) {
	assert(A != NULL);
	free(A->i);
	free(A->j);
	free(A->x);		/* trick : free does nothing on NULL pointer */
	free(A);
}

void spasm_csr_lbresize(spasm * A, int n, int m, int nzmax) {
  int i;
	if (nzmax < 0) {
		nzmax = spasm_nnz(A);
	}
    size_t oldn = A->rowdim();
    A->resize(n,m,nzmax);
    if (oldn < n) {
		for (i = oldn; i < n + 1; i++) {
			A->setStart(i, A->getStart(oldn));
		}
    }
}

// void spasm_csr_resize(spasm * A, int n, int m) {
// 	int i, *Ap;
// 	assert(A != NULL);

// 	A->m = m;
// 	/* in case of a column shrink, check that no entries are left outside */
// 	A->p = spasm_realloc(A->p, (n + 1) * sizeof(int));

// 	if (A->n < n) {
// 		Ap = A->p;
// 		for (i = A->n; i < n + 1; i++) {
// 			Ap[i] = Ap[A->n];
// 		}
// 	}
// 	A->n = n;
// }

void spasm_vector_zero(spasm_GFp * x, int n) {
	int i;

	for (i = 0; i < n; i++) {
		x[i] = 0;
	}
}

void spasm_vector_set(spasm_GFp * x, int a, int b, spasm_GFp alpha) {
	int i;

	for (i = a; i < b; i++) {
		x[i] = alpha;
	}
}


/* add an entry to a triplet matrix; enlarge it if necessary */
void spasm_add_entry(spasm_triplet * T, int i, int j, spasm_GFp x) {
	int prime;
	spasm_GFp x_p;

	assert(T != NULL);
	assert((i >= 0) && (j >= 0));

	prime = T->prime;

	if (T->nz == T->nzmax) {
		spasm_triplet_realloc(T, 2 * T->nzmax);
	}
	if (T->x != NULL) {
		x_p = ((x % prime) + prime) % prime;
		if (x_p == 0) {
			return;
		}
		T->x[T->nz] = x_p;
	}
	T->i[T->nz] = i;
	T->j[T->nz] = j;
	T->nz += 1;
	T->n = spasm_max(T->n, i + 1);
	T->m = spasm_max(T->m, j + 1);
}

void spasm_triplet_transpose(spasm_triplet * T) {
	int i, j, k, nz, *Ti, *Tj;

	assert(T != NULL);
	nz = T->nz;
	Ti = T->i;
	Tj = T->j;

	for (k = 0; k < nz; k++) {
		i = Ti[k];
		j = Tj[k];
		Tj[k] = i;
		Ti[k] = j;
	}
	i = T->m;
	T->m = T->n;
	T->n = i;
}


/* C = compressed-row form of a triplet matrix T */

//spasm *spasm_compress(const spasm_triplet * T) {
spasm *spasm_compress(const spasm_triplet * T, const Field Fp) {
  int m, n, nz, sum, p, k, *w, *Ti, *Tj;
	spasm_GFp *Tx;
	//spasm *C;
	double start;

	m = T->m;
	n = T->n;
	Ti = T->i;
	Tj = T->j;
	Tx = T->x;
	nz = T->nz;
	//prime = T->prime;

	start = spasm_wtime();
	fprintf(stderr, "[CSR] Compressing... ");
	fflush(stderr);

	/* allocate result */
	spasm *C = new spasm(Fp, n, m, nz);
	//C = spasm_csr_alloc(n, m, nz, prime, Tx != NULL);	

	
	/* get workspace */
	w = (int*)spasm_calloc(n, sizeof(int));
// 	Cp = C->p;
// 	Cj = C->j;
// 	Cx = C->x;

	/* compute row counts */
	for (k = 0; k < nz; k++) {
		w[Ti[k]]++;
	}

	/* compute row pointers (in both Cp and w) */
	sum = 0;
	for (k = 0; k < n; k++) {
		C->setStart(k,sum);
// 		Cp[k] = sum;
		sum += w[k];
// 		w[k] = Cp[k];
		w[k] = C->getStart(k);
	}
// 	Cp[n] = sum;
    C->setStart(n,sum);

	/* dispatch entries */
	for (k = 0; k < nz; k++) {
		p = w[Ti[k]]++;	/* A(i,j) is the p-th entry in C */
// 		Cj[p] = Tj[k];
		C->setColid(p,Tj[k]);
// 		if (Cx != NULL) {
// 			Cx[p] = Tx[k];
// 		}
        if (C->getData().size() > 0) {
            C->setData(p,Tx[k]);
        }
	}

	/* success; free w and return C */
	char mem[16];
// 	int size = sizeof(int) * (n + nz) + sizeof(spasm_GFp) * ((Cx != NULL) ? nz : 0);
	int size = sizeof(int) * (n + nz) + sizeof(spasm_GFp) * ((C->getData().size() > 0 ) ? nz : 0);
	spasm_human_format(size, mem);
	fprintf(stderr, "Mem usage = %sbyte [%.2fs]\n", mem, spasm_wtime() - start);
	free(w);
	return C;
}


/*
 * load a matrix in SMS format from f. set prime == -1 to avoid loading
 * values.
 */
spasm_triplet *spasm_load_sms(FILE * f, int prime) {
	int i, j;
	spasm_GFp x;
	spasm_triplet *T;
	char type;
	double start;
	assert(f != NULL);

	start = spasm_wtime();
	if (fscanf(f, "%d %d %c\n", &i, &j, &type) != 3) {
		fprintf(stderr, "[spasm_load_sms] bad SMS file (header)\n");
		exit(1);
	}
	if (prime != -1 && type != 'M') {
		fprintf(stderr, "[spasm_load_sms] only ``Modular'' type supported\n");
		exit(1);
	}
	fprintf(stderr, "[IO] loading %d x %d matrix modulo %d... ", i, j, prime);
	fflush(stderr);

	/* allocate result */
	T = spasm_triplet_alloc(i, j, 1, prime, prime != -1);

	while (fscanf(f, "%d %d %d\n", &i, &j, &x) == 3) {
		if (i == 0 && j == 0 && x == 0) {
			break;
		}
		assert(i != 0);
		assert(j != 0);
		spasm_add_entry(T, i - 1, j - 1, x);
	}

	char nnz[16];
	spasm_human_format(T->nz, nnz);
	fprintf(stderr, "%s NNZ [%.1fs]\n", nnz, spasm_wtime() - start);
	return T;
}


/*
 * C = P.A.Q^-1 where P and Q^-1 are permutations of 0..n-1 and 0..m-1
 * respectively.
 * 
 */
spasm *spasm_permute(const spasm * A, const int *p, const int *qinv, int values) {
// 	int t, j, i, nz, m, n, *Ap, *Aj, *Cp, *Cj;
// 	spasm_GFp *Cx, *Ax;
	int t, j, i, nz, m, n ;
	//spasm *C;
	Field Fp = A->field();

	/* check inputs */
	assert(A != NULL);

// 	n = A->n;
// 	m = A->m;
// 	Ap = A->p;
// 	Aj = A->j;
// 	Ax = A->x;
	n = A->rowdim();
	m = A->coldim();

	/* alloc result */
// 	C = spasm_csr_alloc(n, m, A->nzmax, A->prime, values && (Ax != NULL));
// 	Cp = C->p;
// 	Cj = C->j;
// 	Cx = C->x;
	//C = spasm_csr_alloc(n, m, A->size(), A->field().characteristic(), values && (A->getData().size()>0));
	spasm *C = new spasm(Fp, n, m, A->size());
    nz = 0;

// 	for (i = 0; i < n; i++) {
// 		/* row i of C is row p[i] of A (denoted by j) */
// 		Cp[i] = nz;
// 		j = (p != NULL) ? p[i] : i;
// 		for (t = Ap[j]; t < Ap[j + 1]; t++) {
// 			/* col j of A is col qinv[j] of C */
// 			Cj[nz] = (qinv != NULL) ? qinv[Aj[t]] : Aj[t];
// 			if (Cx != NULL) {
// 				Cx[nz] = Ax[t];
// 			}
// 			nz++;
// 		}
// 	}
// 	/* finalize the last row of C */
// 	Cp[n] = nz;
	for (i = 0; i < n; i++) {
		/* row i of C is row p[i] of A (denoted by j) */
		C->setStart(i,nz);
		j = (p != NULL) ? p[i] : i;
		for (t = A->getStart(j); t < A->getEnd(j); t++) {
			/* col j of A is col qinv[j] of C */
            C->setColid(nz, (qinv != NULL) ? qinv[A->getColid(t)] : A->getColid(t));
			if (C->getData().size()>0) {
                C->setData(nz, A->getData(t));
			}
			nz++;
		}
	}
	/* finalize the last row of C */
    C->setStart(n,nz);
// 	Cp[n] = nz;
	return C;
}

/* 
   spasm_GFp x; F.inv(x,a);
*/

spasm_GFp spasm_GFp_inverse(spasm_GFp a, int prime) {
	int b0 = prime, t, q;
	int x0 = 0, x1 = 1;

	assert(prime > 1);
	while (a > 1) {
		q = a / prime;
		t = prime, prime = a % prime, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0)
		x1 += b0;
	return x1;
}





/*
 * x = x + beta * A[j], where x is a dense vector and A[j] is sparse
 * 
 * low-level operation for maximum flexibility;
 * 
 * This is where all the heavy lifting should take place.
 */
// void spasm_scatter(const int *Aj, const spasm_GFp * Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime) {
// 	int j, p;

// 	for (p = from; p < to; p++) {
// 		j = Aj[p];
// 		/* axpy-inplace */
// 		x[j] = (x[j] + ((beta * Ax[p]))) % prime /* ultra-naive */ ;
// 	}

// }
void spasm_scatter(const spasm::svector_t& Aj, const spasm::row_t& Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, int prime) {
	int j, p;

	for (p = from; p < to; p++) {
		j = Aj[p];
		/* axpy-inplace */
		x[j] = (x[j] + ((beta * Ax[p]))) % prime /* ultra-naive */ ;
		
	}

}

void spasm_lbscatter(const spasm::svector_t& Aj, const spasm::row_t& Ax, int from, int to, spasm_GFp beta, spasm_GFp * x, const Field F) {
	int j, p;

	for (p = from; p < to; p++) {
		j = Aj[p];
		/* axpy-inplace */
		//x[j] = (x[j] + ((beta * Ax[p]))) % prime /* ultra-naive */ ;
		F.axpyin(x[j], beta, Ax[p]);
	}

}

		/* axpy-inplace */
/* A->field().axpyin(x[j], beta, Ax[p]); */ 




int spasm_dfs(int j, const spasm * A, int top, int *xj, int *pstack, int *marks, const int *qinv) {
// 	int px, p2, inew, head, *Ap, *Aj;
	int px, p2, inew, head;

	/* check inputs */
	assert(A != NULL);
	assert(xj != NULL);
	assert(pstack != NULL);
	assert(marks != NULL);

// 	Ap = A->p;
// 	Aj = A->j;
	/*
	 * initialize the recursion stack (rows waiting to be traversed). The
	 * stack is held at the begining of xj, and has head elements.
	 */
	head = 0;
	xj[head] = j;

	/* stack empty ? */
	while (head >= 0) {
		/* get j from the top of the recursion stack */
		j = xj[head];
		inew = (qinv != NULL) ? qinv[j] : j;

		/*
		 * has row i been seen before ? adjacent columns are Gj[
		 * Gp[jnew]     : Gp[jnew + 1] ] UNSEEN columns are   Gj[
		 * pstack[head] : Gp[jnew + 1] ]
		 */

		if (!marks[j]) {
			/* mark node i as seen and initialize pstack. This is done only once. */
			marks[j] = 1;
// 			pstack[head] = (inew < 0) ? 0 : Ap[inew];
			pstack[head] = (inew < 0) ? 0 : A->getStart(inew);
		}
		/* index of last entry of row inew */
// 		p2 = (inew < 0) ? 0 : Ap[inew + 1];
		p2 = (inew < 0) ? 0 : A->getStart(inew + 1);


		/* examine all yet-unseen entries of row i */
		for (px = pstack[head]; px < p2; px++) {
// 			j = Aj[px];
			j = A->getColid(px);
			if (marks[j])
				continue;
			/*
			 * interrupt the enumeration of entries of row inew,
			 * and deal with column j instead. Save status of row
			 * inew in the stack.
			 */
			pstack[head] = px + 1;

			/*
			 * push column j onto the recursion stack. This will
			 * start a DFS from j
			 */
			xj[++head] = j;

			/* node i is not done, exit the loop */
			break;
		}

		/* depth-first search at node i done ? */
		if (px == p2)
			/* push initial column in the output stack and pop it from the recursion stack*/
			xj[--top] = xj[head--];
	}
	return top;
}


/*
 * Compute the set of nodes from G reachable from any node in B[k] (used to
 * determine the pattern of a sparse triangular solve)
 * 
 * G : graph to search
 * 
 * B : RHS (starting point of the search)
 * 
 * k : k-th row of B is used.
 * 
 * l : upper-bound on both dimensions of the matrix
 * 
 * xj: size 3l. Used as workspace. Output in xj[top:l]
 * 
 * pinv: mapping of rows to columns of G.
 * 
 * return value : top
 * 
 * xj [top...l-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * 
 * xj [l...3l-1] used as workspace
 */
int spasm_reach(const spasm * A, const spasm * B, int k, int l, int *xj, const int *qinv) {
// 	int  top, *Bp, *Bj, *pstack, *marks;
	int  top, *pstack, *marks;

	/* check inputs */
	assert(A != NULL);
	assert(B != NULL);
	assert(xj != NULL);

// 	Bp = B->p;
// 	Bj = B->j;
	top = l;
	pstack = xj + l;
	marks = pstack + l;

	/*
	 * iterates over the k-th row of B.  For each column index j present
	 * in B[k], check if j is in the pattern (i.e. if it is marked). If
	 * not, start a DFS from j and add to the pattern all columns
	 * reachable from j.
	 */
// 	for (int px = Bp[k]; px < Bp[k + 1]; px++)
// 		if (!marks[Bj[px]])
// 			top = spasm_dfs(Bj[px], A, top, xj, pstack, marks, qinv);
	for (int px = B->getStart(k); px < B->getEnd(k); px++)
		if (!marks[B->getColid(px)])
			top = spasm_dfs(B->getColid(px), A, top, xj, pstack, marks, qinv);

	/* unmark all marked nodes. */
	/*
	 * TODO : possible optimization : if stuff is marked "k", and
	 * initialized with -1, then this is not necessary
	 */
	for (int px = top; px < l; px++)
		marks[xj[px]] = 0;

	return top;
}


/*
 * (dense vector) * (sparse) Matrix y <--- y + x*A
 */
void spasm_gaxpy(const spasm * A, const spasm_GFp * x, spasm_GFp * y) {
	int i, n;
// 	int *Ap, *Aj, *Ax;

	/* check inputs */
	assert(x != NULL);
	assert(y != NULL);
	assert(A != NULL);

// 	n = A->n;
// 	Ap = A->p;
// 	Aj = A->j;
// 	Ax = A->x;
// 	prime = A->prime;
	n = A->rowdim();
	//prime = A->field().characteristic();

	for (i = 0; i < n; i++) {
// 		spasm_scatter(Aj, Ax, Ap[i], Ap[i + 1], x[i], y, prime);
	  spasm_lbscatter(A->getColid(), A->getData(), A->getStart(i), A->getEnd(i), x[i], y, A->field());
	}
}


/*************** Triangular solving with sparse RHS
 *
 * solve x * U = B[k], where U is (permuted) upper triangular.
 *
 * x has size m (number of columns of U, paradoxically).
 *
 * when this function returns, the solution is scattered in x, and its pattern
 * is given in xi[top : m].
 *
 * top is the return value.
 *
 */
int spasm_sparse_forward_solve(const spasm * U, const spasm * B, int k, int *xj, spasm_GFp * x, const int *qinv) {
// 	int top, m, prime, *Up, *Uj, *Bp, *Bj;
// 	spasm_GFp *Ux, *Bx;
  int top, m, prime = U->field().characteristic();

// 	m = U->m;
// 	Up = U->p;
// 	Uj = U->j;
// 	Ux = U->x;
// 	prime = U->prime;

// 	Bp = B->p;
// 	Bj = B->j;
// 	Bx = B->x;
	m = U->coldim();
	//prime = U->field().characteristic();

	/* xj[top : n] = Reach(U, B[k]) */
	top = spasm_reach(U, B, k, m, xj, qinv);


	/* clear x */
	for (int px = top; px < m; px++)
		x[xj[px]] = 0;

	/* scatter B[k] into x */
// 	for (int px = Bp[k]; px < Bp[k + 1]; px++)
// 		x[Bj[px]] = Bx[px];
	for (int px = B->getStart(k); px < B->getEnd(k); px++)
		x[B->getColid(px)] = B->getData(px);


	/* iterate over the (precomputed) pattern of x (= the solution) */
	for (int px = top; px < m; px++) {
		/* x[j] is nonzero */
		int j = xj[px];

		/* locate corresponding pivot if there is any */
		int i = (qinv != NULL) ? (qinv[j]) : j;
		if (i < 0)
			continue;

		/*
		 * the pivot entry on row i is 1, so we just have to multiply
		 * by -x[j]
		 */
// 		assert(Ux[Up[i]] == 1);
// 		spasm_scatter(Uj, Ux, Up[i] + 1, Up[i + 1], prime - x[j], x, prime);
		assert(U->getData(U->getStart(i)) == 1);
		spasm_lbscatter(U->getColid(), U->getData(), U->getStart(i) + 1, U->getEnd(i), prime - x[j], x, U->field());
	}

	return top;
}


/* eliminate everything in the (dense) vector x using the pivots found in A */
void spasm_eliminate_sparse_pivots(const spasm * A, const int npiv, const int *p, spasm_GFp * x) {
// 	int i, inew, j, prime, *Aj, *Ap;
// 	spasm_GFp *Ax;
  int i, inew, j, prime = A->field().characteristic();


// 	Aj = A->j;
// 	Ap = A->p;
// 	Ax = A->x;
// 	prime = A->prime;
	//prime = A->field().characteristic();

	for (i = 0; i < npiv; i++) {
		inew = (p != NULL) ? p[i] : i;
// 		j = Aj[Ap[inew]];
		j = A->getColid(A->getStart(inew));

		if (x[j] == 0) {
			continue;
		}
// 		spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], prime - x[j], x, prime);
		spasm_lbscatter(A->getColid(), A->getData(), A->getStart(inew), A->getEnd(inew), prime - x[j], x, A->field());
	}
}



/* make pivotal rows of A unitary */
void spasm_make_pivots_unitary(spasm * A, const int *p, const int npiv) {
// 	int prime = A->prime;
  //int prime = A->field().characteristic();
// 	int *Ap = A->p;
// 	spasm_GFp *Ax = A->x;

#pragma omp parallel for
	for (int i = 0; i < npiv; i++) {
		int inew;
		spasm_GFp diag, alpha;

		inew = p[i];
// 		diag = Ax[Ap[inew]];
		diag = A->getData(A->getStart(inew));
		if (diag == 1)
			continue;

		//alpha = spasm_GFp_inverse(diag, prime);
		A->field().inv(alpha,diag);
// 		for (int px = Ap[inew]; px < Ap[inew + 1]; px++)
// 			Ax[px] = (alpha * Ax[px]) % prime;
		for (int px = A->getStart(inew); px < A->getEnd(inew); px++)
		  //A->setData(px,(alpha * A->getData(px)) % prime);
		  
		  A->setData(px,A->field().mul(alpha, alpha, A->getData(px)));
	}
}


/**
 *   Computes a random linear combination of A[k:].
 *   returns TRUE iff it belongs to the row-space of U.
 *   This means that with proba >= 1-1/p, all pivots have been found.
 */
int spasm_early_abort(const spasm * A, const int *p, int k, const spasm * U, int nu) {
// 	int *Aj, *Ap;
	int i, j, inew, n, m, ok;
// 	spasm_GFp prime, *y, *Ax;
	spasm_GFp *y;

// 	n = A->n;
// 	m = A->m;
// 	prime = A->prime;
// 	Aj = A->j;
// 	Ap = A->p;
// 	Ax = A->x;
	n = A->rowdim();
	m = A->coldim();
	//prime = A->field().characteristic();

	y = (spasm_GFp*)spasm_malloc(m * sizeof(spasm_GFp));

	/*
	 * build a random (dense) linear combinations of the row of A, store
	 * in y
	 */
	for (j = 0; j < m; j++) {
		y[j] = 0;
	}
	for (i = k; i < n; i++) {
		inew = (p != NULL) ? p[i] : i;
// 		spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, y, prime);
		spasm_lbscatter(A->getColid(), A->getData(), A->getStart(inew), A->getEnd(inew), rand() % A->field().characteristic(), y, A->field());
	}

	spasm_eliminate_sparse_pivots(U, nu, SPASM_IDENTITY_PERMUTATION, y);

	/* if y != 0, then y does not belong to the row space of U */
	ok = 1;
	for (j = 0; j < m; j++) {
		if (y[j] != 0) {
			ok = 0;
			break;
		}
	}
	free(y);
	return ok;
}


/*
 * compute a (somewhat) LU decomposition using the GPLU algorithm.
 * 
 * r = min(n, m) is an upper-bound on the rank of A
 * 
 * L n * r U is r * m
 * 
L*U == row_permutation*A
 * 
 * qinv[j] = i if the pivot on column j is on row i. -1 if no pivot (yet) found
 * on column j.
 * 
 */
spasm_lu *spasm_LU(const spasm * A, const int *row_permutation, int keep_L) {
  //spasm *L, *U;
  spasm_lu *N;
  //spasm_GFp *Lx, *Ux, *x;
  spasm_GFp *x;
  Field Fp = A->field();
  
  //int *Lp, *Lj, *Up, *Uj, *p, *qinv, *xj;
  int *p, *qinv, *xj;
  int n, m, r, jpiv, i, inew, j, top, px, lnz, unz, old_unz, prime,
    defficiency, verbose_step;
  int rows_since_last_pivot, early_abort_done;
  
  
  
  /* check inputs */
  assert(A != NULL);
  
  n = A->rowdim();
  m = A->coldim();
  //n = A->n;
  //m = A->m;
  r = spasm_min(n, m);
  prime = A->field().characteristic();
  //prime = A->prime;
  defficiency = 0;
  verbose_step = spasm_max(1, n / 1000);
  
  /* educated guess of the size of L,U */
  lnz = 4 * A->size() + n;
  unz = 4 * A->size() + n;
  
  //lnz = 4 * spasm_nnz(A) + n;
  //unz = 4 * spasm_nnz(A) + n;

	/* get GFp workspace */
  x = (spasm_GFp*)spasm_malloc(m * sizeof(spasm_GFp));
  
  /* get int workspace */
  xj = (spasm_GFp*)spasm_malloc(3 * m * sizeof(int));
  spasm_vector_zero(xj, 3 * m);
  
  /* allocate result */
  N = (spasm_lu*)spasm_malloc(sizeof(spasm_lu));
  //N->L = L = (keep_L) ? spasm_csr_alloc(n, r, lnz, prime, true) : NULL;
  //N->U = U = spasm_csr_alloc(r, m, unz, prime, true);
  spasm *L = (keep_L) ? new spasm (Fp, n, r, lnz) : NULL;
  N->L = L;
  spasm *U = new spasm (Fp, r, m, unz);
  N->U = U;
  N->qinv = qinv = (int*)spasm_malloc(m * sizeof(int));
  N->p = p = (int*)spasm_malloc(n * sizeof(int));
  
  //Lp = (keep_L) ? L->p : NULL;
  //Up = U->p;
  
  for (i = 0; i < m; i++) {
    /* clear workspace */
    x[i] = 0;
    /* no rows pivotal yet */
    qinv[i] = -1;
  }
  
  for (i = 0; i < n; i++) {
    /* no rows exchange yet */
    p[i] = i;
  }
  
  /* no rows of U yet */
  for (i = 0; i <= r; i++) {
    U->setStart(i, 0);
  }
  old_unz = lnz = unz = 0;
  
  /* initialize early abort */
  rows_since_last_pivot = 0;
  early_abort_done = 0;
  
  /* --- Main loop : compute L[i] and U[i] ------------------- */
  for (i = 0; i < n; i++) {
    if (!keep_L && i - defficiency == r) {
	    fprintf(stderr, "\n[LU] full rank reached ; early abort\n");
	    break;
    }
    if (!keep_L && !early_abort_done && rows_since_last_pivot > 10 && (rows_since_last_pivot > (n / 100))) {
	    fprintf(stderr, "\n[LU] testing for early abort...");
	    fflush(stderr);
	    if (spasm_early_abort(A, row_permutation, i + 1, U, i - defficiency)) {
	      fprintf(stderr, "SUCCESS\n");
	      break;
	    } else {
	      fprintf(stderr, "FAILED\n");
				fflush(stderr);
	    }
	    early_abort_done = 1;
    }
    /*
     * --- Triangular solve: x * U = A[i]
     * ----------------------------------------
     */
    // if (keep_L) {
    // 	Lp[i] = lnz;	/* L[i] starts here */
    // }
    if (keep_L) {
      L->setStart(i,lnz);	/* L[i] starts here */
    }
    
    //Up[i - defficiency] = unz;	/* U[i] starts here */
    U->setStart(i - defficiency, unz);	/* U[i] starts here */
    
    /* not enough room in L/U ? realloc twice the size */
    // if (keep_L && lnz + m > L->nzmax) {
    // 	spasm_csr_realloc(L, 2 * L->nzmax + m);
    // }
    // if (unz + m > U->nzmax) {
    // 	spasm_csr_realloc(U, 2 * U->nzmax + m);
    // }
    if (keep_L && lnz + m > L->size()) {
      spasm_csr_realloc(L, 2 * L->size() + m);
    }
    if (unz + m > U->size()) {
      spasm_csr_realloc(U, 2 * U->size() + m);
    }
    //Lj = (keep_L) ? L->j : NULL;
    //Lx = (keep_L) ? L->x : NULL;
    //Uj = U->j;
    //Ux = U->x;
    
    inew = (row_permutation != NULL) ? row_permutation[i] : i;
    top = spasm_sparse_forward_solve(U, A, inew, xj, x, qinv);
    
    
    /*
     * --- Find pivot and dispatch coeffs into L and U
     * --------------------------
     */
    
    jpiv = -1;	/* column index of best pivot so far. */
    
    
    for (px = top; px < m; px++) {
      /* x[j] is (generically) nonzero */
      j = xj[px];
      
      /*
       * if x[j] == 0 (numerical cancelation), we just
       * ignore it
       */
      if (x[j] == 0) {
	continue;
      }
      if (qinv[j] < 0) {
	/* column j is not yet pivotal ? */
	
	/* have found the pivot on row i yet ? */
	if (jpiv == -1 || j < jpiv) {
	  jpiv = j;
	}
      } else if (keep_L) {
	/* column j is pivotal */
	/* x[j] is the entry L[i, qinv[j] ] */
	// Lj[lnz] = qinv[j];
	// Lx[lnz] = x[j];
	L->setColid(lnz, qinv[j]);
	L->setData(lnz, x[j]);
	lnz++;
      }
    }
    
    /* pivot found ? */
    if (jpiv != -1) {
      old_unz = unz;
      
      /* L[i,i] <--- x[jpiv]. Last entry of the row ! */
      if (keep_L) {
	// Lj[lnz] = i - defficiency;
	// Lx[lnz] = x[jpiv];
	L->setColid(lnz, i - defficiency);
	L->setColid(lnz, x[jpiv]);
	lnz++;
      }
      qinv[jpiv] = i - defficiency;
      p[i - defficiency] = i;
      
      /* pivot must be the first entry in U[i] */
      //Uj[unz] = jpiv;
      //Ux[unz] = 1;
      U->setColid(unz, jpiv);
      U->setData(unz, 1);
      unz++;
      
      /* send remaining non-pivot coefficients into U */
      //spasm_GFp beta = spasm_GFp_inverse(x[jpiv], prime);
      spasm_GFp beta;
      U->field().inv(beta, x[jpiv]);
      for (px = top; px < m; px++) {
	j = xj[px];
	
	if (qinv[j] < 0) {
	  // Uj[unz] = j;
	  // Ux[unz] = (x[j] * beta) % prime;
	  U->setColid(unz, j);
	  //U->setData(unz, (x[j] * beta) % prime);
	  U->field().mul(beta, beta, x[j]);
	  U->setData(unz, beta);
	  unz++;
	}
      }
      
      /* reset early abort */
      rows_since_last_pivot = 0;
      early_abort_done = 0;
    } else {
      defficiency++;
      p[n - defficiency] = i;
      rows_since_last_pivot++;
    }
    
    
    if ((i % verbose_step) == 0) {
      fprintf(stderr, "\rLU : %d / %d [|L| = %d / |U| = %d] -- current density= (%.3f vs %.3f) --- rank >= %d", i, n, lnz, unz, 1.0 * (m - top) / (m), 1.0 * (unz - old_unz) / m, i - defficiency);
      fflush(stderr);
    }
  }
  
  /*
   * --- Finalize L and U
   * -------------------------------------------------
   */
  fprintf(stderr, "\n");
  
  /* remove extra space from L and U */
  //Up[i - defficiency] = unz;
  U->setStart(i-defficiency,unz);
  // 	spasm_csr_resize(U, i - defficiency, m);
  // 	spasm_csr_realloc(U, -1);
  spasm_csr_lbresize(U,i - defficiency, m, -1);
  
  if (keep_L) {
    //Lp[n] = lnz;
    L->setStart(n, lnz);
    // 		spasm_csr_resize(L, n, n - defficiency);
    // 		spasm_csr_realloc(L, -1);
    spasm_csr_lbresize(L, n, n - defficiency, -1);
    
  }
  free(x);
  free(xj);
  return N;
}


void spasm_free_LU(spasm_lu * X) {
	assert(X != NULL);
	spasm_csr_free(X->L);
	spasm_csr_free(X->U);
	free(X->qinv);
	free(X->p);
	free(X);
}



spasm_dense_lu *spasm_dense_LU_alloc(int m, int prime) {
	spasm_dense_lu *R;

	R = (spasm_dense_lu*)spasm_malloc(sizeof(spasm_dense_lu));
	R->m = m;
	R->n = 0;
	R->prime = prime;
	R->x = (spasm_GFp**)spasm_malloc(m * sizeof(spasm_GFp *));
	R->p = (int*)spasm_malloc(m * sizeof(int));
	return R;
}

void spasm_dense_LU_free(spasm_dense_lu * A) {
	for (int i = 0; i < A->n; i++)
		free(A->x[i]);
	free(A->x);
	free(A->p);
	free(A);
}

int spasm_dense_LU_grow(spasm_dense_lu * A, const spasm_GFp * y, int k, int processed) {
	int n, m, status;
	spasm_GFp **Ax;

#pragma omp critical(dense_LU)
	{
#pragma omp atomic read
	  n = A->n;
	  status = (n == processed);
	  if (status) {
	    m = A->m;
	    A->p[n] = k;
	    Ax = A->x;
	    Ax[n] = (spasm_GFp*)spasm_malloc(m * sizeof(spasm_GFp));
	    for (int j = 0; j < m; j++) {
	      Ax[n][j] = y[j];
			}
#pragma omp atomic update
			A->n++;
		}
	}
	return status;
}


/*
 * if y belongs to the linear span of U, return 0. Else update U and return
 * 1. This function is THREAD-SAFE.
 */
int spasm_dense_LU_process(spasm_dense_lu * A, spasm_GFp * y) {
	int processed, k, n, m, prime, *p;
	spasm_GFp beta, **Ax;

	m = A->m;
	prime = A->prime;
	p = A->p;
	Ax = A->x;
	processed = 0;

	while (1) {
#pragma omp atomic read
		n = A->n;

		for (int i = processed; i < n; i++) {
			beta = prime - y[p[i]];
			for (int j = 0; j < m; j++)
				y[j] = (y[j] + beta * Ax[i][j]) % prime;
		}
		processed = n;

		for (k = 0; k < m; k++)
			if (y[k])
				break;
		if (k == m)
			return 0;

		/* make pivot unitary */
		beta = spasm_GFp_inverse(y[k], prime);
		
		for (int j = 0; j < m; j++)
			y[j] = (y[j] * beta) % prime;

		if (spasm_dense_LU_grow(A, y, k, processed))
			return 1;
	}
}

/*
 * Computes the Schur complement, by eliminating the pivots located on rows
 * p[0] ... p[n_pivots-1] of input matrix A. The pivots must be the entries
 * on the lines. This returns a sparse representation of S. The pivots must
 * be unitary.
 */
spasm *spasm_schur(spasm * A, const int *p, const int *qinv, const int npiv) {
  //spasm *S;
	//int k, *Sp, *Sj, Sn, Sm, m, n, snz, *xj, top, *q, verbose_step;
         size_t k, Sn, Sm, m, n, snz;
	 int *xj, top, *q, verbose_step;
	//spasm_GFp *Sx, *x;
	spasm_GFp *x;
	Field Fp = A->field();

	/* check inputs */
	assert(A != NULL);
	//n = A->n;
	//m = A->m;
	n = A->rowdim();
	m = A->coldim();

	assert(n >= npiv);
	assert(m >= npiv);

	

	/* Get Workspace */
	Sn = n - npiv;
	Sm = m - npiv;
	snz = 4 * (Sn + Sm);	/* educated guess */
	//S = spasm_csr_alloc(Sn, Sm, snz, A->field().characteristic(), SPASM_WITH_NUMERICAL_VALUES);
	spasm *S = new spasm(Fp, Sn, Sm, snz);

	x = (spasm_GFp*)spasm_malloc(m * sizeof(spasm_GFp));
	xj = (int*)spasm_malloc(3 * m * sizeof(size_t));
	spasm_vector_zero(xj, 3 * m);

	verbose_step = spasm_max(1, n / 1000);
	//Sp = S->p;
	//Sj = S->j;
	//Sx = S->x;

	

	/*
	 * q sends the non-pivotal columns of A to the columns of S. It is
	 * not the inverse of qinv...
	 */
	
       	
	q = (int*)spasm_malloc(m * sizeof(int));
	k = 0;
	for (size_t j = 0; j < m; j++){
	  
		q[j] = (qinv[j] < 0) ? k++ : -1;
		
	}

	snz = 0;		/* non-zero in S */
	Sn = 0;			/* rows in S */

	

	fprintf(stderr, "Starting Schur complement computation...\n");
	for (int i = npiv; i < n; i++) {
		int inew = p[i];

		/* triangular solve */
		top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);

		/* not enough room in S ? realloc twice the size */
		//if (snz + Sm > S->nzmax) {
		//spasm_csr_realloc(S, 2 * S->nzmax + Sm);
			//Sj = S->j;
			//Sx = S->x;
		//	}

		
		if (snz + Sm > S->size()) {
		  spasm_csr_realloc(S, 2 * S->size() + Sm);
		  printf("\n %zu + %zu vs %zu\n", snz, Sm, S->size());
		}
		//Sp[Sn] = snz;	/* S[i] starts here */
		
		S->setStart(Sn, snz);	/* S[i] starts here */

		
		for (int px = top; px < m; px++) {
		  int j = xj[px];
		  
		  if (x[j] == 0)
		    continue;
		  
		  /* save non-zero, non-pivot coefficients in S */
		  if (q[j] >= 0) {
		    //Sj[snz] = q[j];
		    //Sx[snz] = x[j];

		    S->setColid(snz, q[j]);
		    S->setData(snz, x[j]);
		    snz++;
		    
		  }
		  
		  
		  
		}
		
		Sn++;
		

		if ((i % verbose_step) == 0) {
			fprintf(stderr, "\rSchur : %d / %zu [S=%zu * %zu, %zu NNZ] -- current density= (%.3f)", i, n, Sn, Sm, snz, 1.0 * snz / (1.0 * Sm * Sn));
			fflush(stderr);
		}
	}
	

	//for(int i = 0; i< Sn; i++){
	// for(int px = 120; px < 150; px++){
	//     printf("%d: S[i, %zu] = %d\n", px, S->getColid(px), S->getData(px));
	//   }
	  //}
	
        /*special case: S empty*/
    if(snz == 0){
        fprintf(stderr, "Empty Schur\n");
        return NULL;
    }

	/* finalize S */
	fprintf(stderr, "\n");
	//Sp[S->n] = snz;
	S->setStart(n, S->rowdim());
	spasm_csr_realloc(S, -1);

	/* free extra workspace */
	free(q);
	free(x);
	free(xj);

	return S;
}


/** Samples R rows at random in the schur complement of A w.r.t. the pivots in p[0:n_pivots],
* and return the number that are non-zero (these rows of A are linearly independent from the pivots).
* The pivots must be unitary.
*/
double spasm_schur_probe_density(spasm * A, const int *p, const int *qinv, const int npiv, const int R) {
	int m, n, *xj, top, nnz;
	spasm_GFp *x;

	/* check inputs */
	//m = A->m;
	//n = A->n;
	m = A->coldim();
	n = A->rowdim();
	
	/* Get Workspace */
	x = (spasm_GFp*)spasm_malloc(m * sizeof(spasm_GFp));
	xj = (int*)spasm_malloc(3 * m * sizeof(int));
	spasm_vector_zero(xj, 3 * m);

	nnz = 0;
	for (int i = 0; i < R; i++) {
		/* pick a random row in S, check if non-zero */
		int inew = p[npiv + (rand() % (n - npiv))];
		top = spasm_sparse_forward_solve(A, A, inew, xj, x, qinv);
		for (int px = top; px < m; px++) {
			int j = xj[px];
			if (qinv[j] < 0 && x[j] != 0)
				nnz++;
		}
	}

	/* free extra workspace */
	free(x);
	free(xj);

	return ((double)nnz) / (m - npiv) / R;
}

/*
 * computes the rank of the schur complement, but not the schur complement
 * itself. The pivots must be unitary.
 */
int spasm_schur_rank(spasm * A, const int *p, const int *qinv, const int npiv) {
	int Sm, m, n, k, r, prime, step, threads, searched, prev_r;
	//int *q, *Ap, *Aj;
	int *q;
	double start;
	//spasm_GFp *Ax;

	n = A->rowdim();
	m = A->coldim();
	//Ap = A->p;
	//Aj = A->j;
	//Ax = A->x;
	//prime = A->prime;
	prime = A->field().characteristic();

	/* Get Workspace */
	Sm = m - npiv;
	q = (int*)spasm_malloc(Sm * sizeof(int));

	/* q sends columns of S to non-pivotal columns of A */
	k = 0;
	for (int j = 0; j < m; j++)
		if (qinv[j] < 0)
			q[k++] = j;

	spasm_dense_lu *U = spasm_dense_LU_alloc(Sm, prime);

	/* ---- compute Schur complement ----- */
	fprintf(stderr, "rank of dense schur complement...\n");

	start = spasm_wtime();
	r = 0;
	step = 1;
	k = 0;
	searched = 0;
	threads = 1;
	prev_r = 0;
#ifdef USE_OPENMP
	threads = omp_get_num_threads();
#endif

#pragma omp parallel
	{
	  spasm_GFp *x = (spasm_GFp*)spasm_malloc(m * sizeof(spasm_GFp));
	  spasm_GFp *y = (spasm_GFp*)spasm_malloc(Sm * sizeof(spasm_GFp));
		int gain;

		while (step <= (1 << 16)) {	/* <--- tweak-me */
			double it_start = spasm_wtime();
			prev_r = r;

			/* random linear combination */
			spasm_vector_zero(x, m);
			for (int i = 0; i < step; i++) {
				int inew = p[npiv + (rand() % (n - npiv))];
				//spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], 1 + (rand() % (prime - 1)), x, prime);
				spasm_lbscatter(A->getColid(), A->getData(), A->getStart(inew), A->getEnd(inew), (rand() %(prime -1)), y, A->field());	
			}
			spasm_eliminate_sparse_pivots(A, npiv, p, x);
			for (int j = 0; j < Sm; j++)	/* gather into y */
				y[j] = x[q[j]];

#pragma omp atomic update
			r += spasm_dense_LU_process(U, y);

			/* this is a barrier */
#pragma omp single
			{
				fprintf(stderr, "\rSchur : %d [%.1fs] -- current rank = %d / step = %d", k, spasm_wtime() - it_start, r, step);
				fflush(stderr);

				k++;
				searched += threads * step;
				gain = r - prev_r;

				if (gain < threads)
					step *= 2;
				else
					step = spasm_max(1, step / 2);
			}
		}

#pragma omp single
		{
			int final_bad = 0;
			k = 0;
			fprintf(stderr, "\n");

			while (final_bad < 3) {
				double it_start = spasm_wtime();
				for (int i = npiv; i < n; i++) {
					int inew = p[i];
					//spasm_scatter(Aj, Ax, Ap[inew], Ap[inew + 1], rand() % prime, x, prime);
					spasm_lbscatter(A->getColid(), A->getData(), A->getStart(inew), A->getEnd(inew), rand()%prime, y, A->field());
				}
				spasm_eliminate_sparse_pivots(A, npiv, p, x);
				for (int j = 0; j < Sm; j++)
					y[j] = x[q[j]];
				int newr = spasm_dense_LU_process(U, y);
				r += newr;
				final_bad += 1 - newr;
				k++;
				fprintf(stderr, "\rSchur : %d [%.1fs] -- current rank = %d / final", k, spasm_wtime() - it_start, r);
				fflush(stderr);
			}
		}
		free(x);
		free(y);
	}
	fprintf(stderr, "\n[schur/rank] Time: %.1fs\n", spasm_wtime() - start);

	free(q);
	spasm_dense_LU_free(U);
	return r;
}


int spasm_is_row_pivotal(const spasm * A, const int *qinv, const int i) {
  //int *Ap, *Aj;

  if(A->getStart(i) == A->getEnd(i)){ // empty row
    return 0;
  }

  //Ap = A->p;
  //Aj = A->j;
  //return (qinv[Aj[Ap[i]]] == i);
  return (qinv[A->getColid(A->getStart(i))] == i);
}

/* make pivot the first entry of the row */
void spasm_prepare_pivot(spasm * A, const int i, const int px) {
  //int *Ap, *Aj;
  //spasm_GFp *Ax;

  //Ap = A->p;
  //Aj = A->j;
  //Ax = A->x;

  //spasm_swap(Aj, Ap[i], px);
  //if (Ax != NULL)
  //spasm_swap(Ax, Ap[i], px);
  //std::cout << "i : " << i << " j : " << A->getStart(px) << std::endl;
  
  spasm_lbswap(A, A->getStart(i), px);
  
  
}


/** Faugère-Lachartre pivot search.
 *
 * The leftmost entry of each row is a candidate pivot. Select the sparsest row
 * with a leftmost entry on the given column. Selected pivots are moved to the
 * front of the row.
 * @param qinv must be initialized to -1
 * @param p can be arbitrary.
 * @return number of pivots found. */
int spasm_find_FL_pivots(spasm * A, int *p, int *qinv) {
  //int n, m, idx_j, *Aj, *Ap, npiv;
  int n, m, idx_j, npiv, px, i;
  double start;

  //n = A->n;
  //m = A->m;
  //Ap = A->p;
  //Aj = A->j;

  n = A->rowdim();
  m = A->coldim();
  
  start = spasm_wtime();

  
  for (i = 0; i < n; i++) {
    int j = -1;
    // for (int px = Ap[i]; px < Ap[i + 1]; px++)
    //   if (j == -1 || Aj[px] < j) {
    // 	j = Aj[px];
    // 	idx_j = px;
    //   }
    for (px = A->getStart(i); px < A->getEnd(i); px++)
      if (j == -1 || A->getColid(px) < j) {
	j = A->getColid(px);
	idx_j = px;
      }

    if (j == -1)	/* Skip empty rows */
      continue;

    
    /* check if it is a sparser pivot */
    if (qinv[j] == -1 || spasm_row_weight(A, i) < spasm_row_weight(A, qinv[j])) {
    qinv[j] = i;
    spasm_prepare_pivot(A, i, idx_j);
    }
  }
  
  /* build p */
  npiv = 0;
  for (int j = 0; j < m; j++)
    if (qinv[j] != -1)
      p[npiv++] = qinv[j];
  
  fprintf(stderr, "[pivots] Faugère-Lachartre: %d pivots found [%.1fs]\n", npiv, spasm_wtime() - start);
  return npiv;
}



/*
 * Leftovers from FL. Column not occuring on previously selected pivot row
 * can be made pivotal, as this will not create alternating cycles.
 * 
 * w[j] = 1 <===> column j does not appear in a pivotal row
 * 
 */
int spasm_find_FL_column_pivots(spasm * A, int *p, int *qinv, int npiv_fl) {
  //int n, m, *Aj, *Ap, npiv, *w;
  int n, m, npiv, *w;
  double start;

  //n = A->n;
  //m = A->m;
  //Ap = A->p;
  //Aj = A->j;
  n = A->rowdim();
  m = A->coldim();
  npiv = npiv_fl;
  w = (int*)spasm_malloc(m * sizeof(int));
  spasm_vector_set(w, 0, m, 1);
  
  start = spasm_wtime();

  /* mark columns on previous pivot rows as obstructed */
  for (int i = 0; i < npiv; i++) {
    int inew = p[i];
    // for (int px = Ap[inew]; px < Ap[inew + 1]; px++)
    //   w[Aj[px]] = 0;

    
    for (int px = A->getStart(inew); px < A->getEnd(inew); px++){
     
      w[A->getColid(px)] = 0;
    }
    
  }

  
  /* find new pivots */
  for (int i = 0; i < n; i++) {
    
    if (spasm_is_row_pivotal(A, qinv, i)){
      continue;
    }
    
    /* does A[i,:] have an entry on an unobstructed column? */
    // for (int px = Ap[i]; px < Ap[i + 1]; px++) {
    //   int j = Aj[px];
    //   if (w[j] == 0)
    // 	continue;	/* this column is closed,
    // 			 * skip this entry */
    
    // 			/* new pivot found! */
    //   if (qinv[j] == -1) {
    // 	p[npiv++] = i;
    // 	qinv[j] = i;
    // 	spasm_prepare_pivot(A, i, px);
    // 	/*
    // 	 * mark the columns occuring on this row as
    // 	 * unavailable
    // 	 */
    // 	for (int px = Ap[i]; px < Ap[i + 1]; px++)
    // 	  w[Aj[px]] = 0;
    
    // 	break;
    //   }
    // }
    
    
    for (int px = A->getStart(i); px < A->getEnd(i); px++) {
      int j = A->getColid(px);
      
      if (w[j] == 0)
	continue;	/* this column is closed,
			 * skip this entry */
      
			/* new pivot found! */
      if (qinv[j] == -1) {
	p[npiv++] = i;
	qinv[j] = i;
	spasm_prepare_pivot(A, i, px);
	/*
	 * mark the columns occuring on this row as
	 * unavailable
	 */
	for (int px = A->getStart(i); px < A->getEnd(i); px++){
	  int j = A->getColid(px);
	  w[j] = 0;
	}
	
	break;
      }
    }
  }
  free(w);
	fprintf(stderr, "[pivots] ``Faugère-Lachartre on columns'': %d pivots found [%.1fs]\n", npiv - npiv_fl, spasm_wtime() - start);
	return npiv;
}



int find_survivor(spasm * A, int i, int *w) {
  //int *Ap = A->p;
  //int *Aj = A->j;

	// for (int px = Ap[i]; px < Ap[i + 1]; px++) {
	// 	int j = Aj[px];
	// 	if (w[j] == 1) {/* potential pivot found */
	// 		spasm_prepare_pivot(A, i, px);
	// 		return j;
	// 	}
	// }
  for (int px = A->getStart(i); px < A->getEnd(i); px++) {
    int j = A->getColid(px);
    if (w[j] == 1) {/* potential pivot found */
      spasm_prepare_pivot(A, i, px);
      return j;
    }
  }
  return -1;
}

/*
 * provide already know pivots, and this looks for more. Updates qinv, but
 * DFS must be performed afterwards
 */
void BFS_enqueue(int *w, int *queue, int *surviving, int *tail, int j) {
	queue[(*tail)++] = j;
	*surviving -= w[j];
	w[j] = -1;
}

void BFS_enqueue_row(int *w, int *queue, int *surviving, int *tail, const spasm::svector_t& Ap, const spasm::svector_t& Aj, int i) {
	for (int px = Ap[i]; px < Ap[i + 1]; px++) {
		/* this is the critical section */
		int j = Aj[px];
		if (w[j] >= 0)
			BFS_enqueue(w, queue, surviving, tail, j);
	}
}

int spasm_find_cycle_free_pivots(spasm * A, int *p, int *qinv, int npiv_start) {
  //int n, m, *Aj, *Ap, processed, v, npiv, retries;
  int n, m, processed, v, npiv, retries;
  double start;

  //n = A->n;
  //m = A->m;
  //Ap = A->p;
  //Aj = A->j;
  n = A->rowdim();
  m = A->coldim();
  
  v = spasm_max(1, spasm_min(1000, n / 100));
  processed = 0;
  retries = 0;
  npiv = npiv_start;
  start = spasm_wtime();
  
#pragma omp parallel
  {
    int *w = (int*)spasm_malloc(m * sizeof(int));
    int *queue = (int*)spasm_malloc(2 * m * sizeof(int));
    int head, tail, npiv_local, surviving, tid;
    
    /* workspace initialization */
    tid = 0;
    spasm_vector_set(w, 0, m, 0);
    
#ifdef USE_OPENMP
    tid = omp_get_thread_num();
    if (tid == 0)
			fprintf(stderr, "[pivots] Greedy pivot search starting on %d threads\n", omp_get_num_threads());
#endif


#pragma omp for schedule(dynamic, 1000)
		for (int i = 0; i < n; i++) {
			/*
			 * for each non-pivotal row, computes the columns
			 * reachable from its entries by alternating paths.
			 * Unreachable entries on the row can be chosen as
			 * pivots. The w[] array is used for marking during
			 * the graph traversal. Before the search: w[j] ==  1
			 * for each non-pivotal entry j on the row w[j] ==  0
			 * otherwise After the search: w[j] ==  1  for each
			 * unreachable non-pivotal entry j on the row w[j] ==
			 * -1  column j is reachable by an alternating path,
			 * or is pivotal (has entered the queue at some
			 * point) w[j] ==  0  column j was absent and is
			 * unreachable
			 */
			if (tid == 0 && (i % v) == 0) {
				fprintf(stderr, "\r[pivots] %d / %d --- found %d new --- %d retries", processed, n - npiv_start, npiv - npiv_start, retries);
				fflush(stderr);
			}
			if (spasm_is_row_pivotal(A, qinv, i))
				continue;

#pragma omp atomic update
			processed++;

			/* we start reading qinv: begining of transaction */
#pragma omp atomic read
			npiv_local = npiv;

			/*
			 * scatters columns of A[i] into w, enqueue pivotal
			 * entries
			 */
			head = 0;
			tail = 0;
			surviving = 0;
			// for (int px = Ap[i]; px < Ap[i + 1]; px++) {
			// 	int j = Aj[px];
			// 	if (qinv[j] < 0) {
			// 		w[j] = 1;
			// 		surviving++;
			// 	} else {
			// 		BFS_enqueue(w, queue, &surviving, &tail, j);
			// 	}
			// }
			for (int px = A->getStart(i); px < A->getEnd(i); px++) {
			  int j = A->getColid(px);
			  if (qinv[j] < 0) {
			    w[j] = 1;
			    surviving++;
			  } else {
			    BFS_enqueue(w, queue, &surviving, &tail, j);
			  }
			}
			
			
			/* BFS. This is where most of the time is spent */
		BFS:
			while (head < tail && surviving > 0) {
			  int j = queue[head++];
			  
			  int I = qinv[j];
			  
			  if (I == -1)
			    continue;	/* j is not pivotal:
					 * nothing to do */
			  //BFS_enqueue_row(w, queue, &surviving, &tail, Ap, Aj, I);

			  BFS_enqueue_row(w, queue, &surviving, &tail, A->getStart(), A->getColid(), I);
			}
			
			/* scan w for surviving entries */
			if (surviving > 0) {
			  int j = find_survivor(A, i, w);
			  int npiv_target = -1;
			  
			  /*
			   * si aucun nouveau pivot n'est arrivé,
			   * ajouter
			   */
#pragma omp critical
			  {
			    if (npiv == npiv_local) {
			      qinv[j] = i;
			      p[npiv] = i;
#pragma omp atomic update
			      npiv++;
			    } else {
#pragma omp atomic read
			      npiv_target = npiv;
			      retries++;
			    }
			  }
			  
			  if (npiv_target < 0)
			    goto cleanup;
			  
			  /*
			   * si on a découvert de nouveaux pivots
			   * aiter... les traiter !
			   */
			  for (; npiv_local < npiv_target; npiv_local++) {
			    int I = p[npiv_local];
			    //int j = Aj[Ap[I]];
			    int j = A->getColid(A->getStart(I));
			    if (w[j] == 0)	/* the new pivot plays
						 * no role here */
			      continue;
			    
			    if (w[j] == 1) {
			      /*
			       * a survivors becomes
			       * pivotal with this pivot
			       */
			      BFS_enqueue(w, queue, &surviving, &tail, j);
			    } else {
			      /* the new pivot has been hit */
			      BFS_enqueue_row(w, queue, &surviving, &tail, A->getStart(), A->getColid(), I);
			    }
			  }
			  goto BFS;
			}
			/* reset w back to zero */
		cleanup:
			// for (int px = Ap[i]; px < Ap[i + 1]; px++)
			//   w[Aj[px]] = 0;
			// for (int px = 0; px < tail; px++)
			//   w[queue[px]] = 0;
			for (int px = A->getStart(i); px < A->getEnd(i); px++)
			  w[A->getColid(px)] = 0;
		        for (int px = 0; px < tail; px++)
			  w[queue[px]] = 0;
		}		/* end for */
		free(w);
		free(queue);
  }			/* end of omp parallel */

  fprintf(stderr, "\r[pivots] greedy alternating cycle-free search: %d pivots found [%.1fs]\n", npiv - npiv_start, spasm_wtime() - start);

  //printf("npiv : %d\n", npiv);
  //printf("p[0] = %d\n", p[0]);
  return npiv;
}

/*
 * return the number of pivots found. @param p : row permutations. Pivotal
 * rows are first. @param qinv : inverse column permutation. q[j] is the row
 * on which the pivot on column j is, or -1 if there is no pivot on column j.
 * both p and qinv must be preallocated
 */
int spasm_find_pivots(spasm * A, int *p, int *qinv) {
  int n, m, k, npiv, top;

  //n = A->n;
  //m = A->m;
  n = A->rowdim();
  m = A->coldim();

  spasm_vector_set(qinv, 0, m, -1);
  npiv = spasm_find_FL_pivots(A, p, qinv);
  npiv = spasm_find_FL_column_pivots(A, p, qinv, npiv);
  npiv = spasm_find_cycle_free_pivots(A, p, qinv, npiv);

	/*
	 * build row permutation. Pivotal rows go first in topological order,
	 * then non-pivotal, non-zero rows, then zero rows
	 */
  int *xj = (int*)spasm_malloc(m * sizeof(int));
  int *marks = (int*)spasm_malloc(m * sizeof(int));
  int *pstack = (int*)spasm_malloc(n * sizeof(int));

  /* topological sort */
  spasm_vector_set(marks, 0, m, 0);
  top = m;
  for (int j = 0; j < m; j++)
    if (qinv[j] != -1 && !marks[j])
      top = spasm_dfs(j, A, top, xj, pstack, marks, qinv);
  k = 0;
  for (int j = top; j < m; j++) {
    int i = qinv[xj[j]];
    if (i != -1)
      p[k++] = i;
  }
  
  for (int i = 0; i < n; i++)
    if (spasm_row_weight(A, i) > 0 && !spasm_is_row_pivotal(A, qinv, i))
      p[k++] = i;
  
  for (int i = 0; i < n; i++)
    if (spasm_row_weight(A, i) == 0)
      p[k++] = i;
  
  free(xj);
  free(pstack);
  free(marks);
  fprintf(stderr, "\r[pivots] %d pivots found\n", npiv);
  return npiv;
}

/*
 * returns a permuted version of A where pivots are pushed to the top-left
 * and form an upper-triangular principal submatrix. qinv is modified.
 */
spasm *spasm_permute_pivots(const spasm * A, const int *p, int *qinv, int npiv) {
  //int k, m, *Ap, *Aj;
  int k, m;

  //m = A->m;
  //Ap = A->p;
  //Aj = A->j;
  m = A->coldim();
  
  /* pivotal columns first */
  k = 0;
  for (int i = 0; i < npiv; i++) {
    //int j = Aj[Ap[p[i]]];   /*the pivot is the first entry of*/
				 /* each row */
    int j = A->getColid(A->getStart(p[i]));     
    qinv[j] = k++;
  }
  
  /* put remaining non-pivotal columns afterwards, in any order */
  for (int j = 0; j < m; j++)
    if (qinv[j] == -1)
      qinv[j] = k++;
  
  return spasm_permute(A, p, qinv, SPASM_WITH_NUMERICAL_VALUES);
}


/** Computes the rank of the input matrix using the hybrid strategy */
int main(int argc, char **argv) {

	/* charge la matrice depuis l'entrée standard */
	int prime, n_times, rank, npiv, n, m, dense_final, gplu_final,
	    allow_transpose, ch;
	double start_time, end_time;
	spasm_triplet *T;
	spasm *A, *B;
	int *p, *qinv;
	double density, sparsity_threshold;
	char nnz[6];

	allow_transpose = 1;	/* transpose ON by default */
	n_times = 3;
	dense_final = 0;
	gplu_final = 0;
	sparsity_threshold = 0.1;
	prime = 42013;

	/* options descriptor */
	struct option longopts[7] = {
		{"sparse-threshold", required_argument, NULL, 's'},
		{"max-recursion", required_argument, NULL, 'm'},
		{"dense-last-step", no_argument, NULL, 'd'},
		{"gplu-last-step", no_argument, NULL, 'g'},
		{"no-transpose", no_argument, NULL, 'a'},
		{"modulus", required_argument, NULL, 'p'},
		{NULL, 0, NULL, 0}
	};

	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 's':
			sparsity_threshold = atof(optarg);
			break;
		case 'm':
			n_times = atoi(optarg);
			break;
		case 'd':
			dense_final = 1;
			break;
		case 'a':
			allow_transpose = 0;
			break;
		case 'g':
			gplu_final = 1;
			break;
		case 'p':
			prime = atoi(optarg);
			break;
		default:
			printf("Unknown option\n");
			exit(1);
		}
	}
	argc -= optind;
	argv += optind;


	T = spasm_load_sms(stdin, prime);
	if (allow_transpose && (T->n < T->m)) {
		fprintf(stderr, "[rank] transposing matrix : ");
		fflush(stderr);
		start_time = spasm_wtime();
		spasm_triplet_transpose(T);
		fprintf(stderr, "%.1f s\n", spasm_wtime() - start_time);
	}
	A = spasm_compress(T, prime);
	spasm_triplet_free(T);
	n = A->rowdim();
	m = A->coldim();
	spasm_human_format(spasm_nnz(A), nnz);
	fprintf(stderr, "start. A is %d x %d (%s nnz)\n", n, m, nnz);

	p = (int*)spasm_malloc(n * sizeof(int));
	qinv = (int*)spasm_malloc(m * sizeof(int));

	start_time = spasm_wtime();
	rank = 0;
	npiv = spasm_find_pivots(A, p, qinv);
	spasm_make_pivots_unitary(A, p, npiv);

	density = spasm_schur_probe_density(A, p, qinv, npiv, 100);	
	
	for (int i = 0; i < n_times; i++) {
		int64_t nnz = (density * (n - npiv)) * (m - npiv);
		char tmp[6];
		spasm_human_format(sizeof(int) * (n - npiv + nnz) + sizeof(spasm_GFp) * nnz, tmp);
		fprintf(stderr, "round %d / %d. Schur complement is %d x %d, estimated density : %.2f (%s byte)\n", i, n_times, n - npiv, m - npiv, density, tmp);

		if (density > sparsity_threshold)
			break;

		/* compute schur complement, update matrix */
		
		B = spasm_schur(A, p, qinv, npiv);
		spasm_csr_free(A);
		
		

		//special case : empty schur.
		if(B==NULL){
		  end_time = spasm_wtime();
		  fprintf(stderr, "done in %.3f s rank = %d\n", end_time - start_time, rank);
		  free(p);
		  free(qinv);
		  return 0;
		}
        

		A = B;
		rank += npiv;
		//n = A->n;
		//m = A->m;
		n = A->rowdim();
		m = A->coldim();

		npiv = spasm_find_pivots(A, p, qinv);
		break;
		
		spasm_make_pivots_unitary(A, p, npiv);
		density = spasm_schur_probe_density(A, p, qinv, npiv, 100);
	}

	exit(1); // Debug end here
	/* ---- final step ---------- */

	/* sparse schur complement : GPLU */
	if (gplu_final || (!dense_final && density < sparsity_threshold)) {
		spasm_lu *LU = spasm_LU(A, p, SPASM_DISCARD_L);
		//rank += LU->U->n;
		rank += LU->U->rowdim();
		spasm_free_LU(LU);
	} else {
		/* dense schur complement */
		int r = spasm_schur_rank(A, p, qinv, npiv);
		fprintf(stderr, "rank = %d + %d\n", npiv, r);
		rank += npiv + r;
	}

	end_time = spasm_wtime();
	fprintf(stderr, "done in %.3f s rank = %d\n", end_time - start_time, rank);
	spasm_csr_free(A);
	free(p);
	free(qinv);
	return 0;
}
