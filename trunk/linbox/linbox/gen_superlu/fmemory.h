

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/util.h"

#include "linbox/gen_superlu/get_perm_c.h"
#include "linbox/gen_superlu/sp_ienv.h"
#include "linbox/gen_superlu/sp_preorder.h"
#include "linbox/gen_superlu/superlu_timer.h"
#include "linbox/gen_superlu/lsame.h"


/* Constants */
#define NO_MEMTYPE  4      /* 0: lusup;
			      1: ucol;
			      2: lsub;
			      3: usub */
#define GluIntArray(n)   (5 * (n) + 5)

/* Internal prototypes 
   void  *Fexpand (int *, MemType,int, int, GlobalLU_t *, Field& F);
   int FLUWorkInit (int, int, int, int **, Field::Element **, LU_space_t, Field& F);
   void  copy_mem_Field (int, void *, void *, Field& F);
   void  FStackCompress (GlobalLU_t *, Field& F);
   void  FSetupSpace (void *, int, LU_space_t *, Field& F);
   void  *Fuser_malloc (int, int, Field& F);
   void  Fuser_free (int, int, Field& F); */

/* External prototypes (in memory.c - prec-indep) */
extern void    copy_mem_int    (int, void *, void *);
extern void    user_bcopy      (char *, char *, int);

/* Headers for 4 types of dynamatically managed memory */
typedef struct e_node {
    int size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

typedef struct {
    int  size;
    int  used;
    int  top1;  /* grow upward, relative to &array[0] */
    int  top2;  /* grow downward */
    void *array;
} LU_stack_t;

/* Variables local to this file */
static ExpHeader *expanders = 0; /* Array of pointers to 4 types of memory */
static LU_stack_t _stack;
static int no_expand;

/* Macros to manipulate stack */
#define StackFull(x)         ( x + _stack.used >= _stack.size )
#define NotDoubleAlign(addr) ( (long int)addr & 7 ) // ??
#define DoubleAlign(addr)    ( ((long int)addr + 7) & ~7L ) // ??
#define TempSpace(m, w)      ( (2*w + 4 + NO_MARKER) * m * sizeof(int) + \
			      (w + 1) * m * sizeof(Field::Element) )
#define Reduce(alpha)        ((alpha + 1) / 2)  /* i.e. (alpha-1)/2 + 1 */

/*
 * Setup the memory model to be used for factorization.
 *    lwork = 0: use system malloc;
 *    lwork > 0: use user-supplied work[] space.
 */
// A.Duran 8/14/2002
template <class Field>
void FSetupSpace(void *work, int lwork, LU_space_t *MemModel, Field& F)
{
    if ( lwork == 0 ) {
	*MemModel = SYSTEM; /* malloc/free */
    } else if ( lwork > 0 ) {
	*MemModel = USER;   /* user provided space */
	_stack.used = 0;
	_stack.top1 = 0;
	_stack.top2 = (lwork/4)*4; /* must be word addressable */
	_stack.size = _stack.top2;
	_stack.array = (void *) work;
    }
}

// A.Duran 8/14/2002
template <class Field>
void *Fuser_malloc(int bytes, int which_end, Field& F)
{
    void *buf;
    
    if ( StackFull(bytes) ) return (NULL);

    if ( which_end == HEAD ) {
	buf = (char*) _stack.array + _stack.top1;
	_stack.top1 += bytes;
    } else {
	_stack.top2 -= bytes;
	buf = (char*) _stack.array + _stack.top2;
    }
    
    _stack.used += bytes;
    return buf;
}

template <class Field>
void Fuser_free(int bytes, int which_end, Field& F)
{
    if ( which_end == HEAD ) {
	_stack.top1 -= bytes;
    } else {
	_stack.top2 += bytes;
    }
    _stack.used -= bytes;
}

/*
 * mem_usage consists of the following fields:
 *    - for_lu (float)
 *      The amount of space used in bytes for the L\U data structures.
 *    - total_needed (float)
 *      The amount of space needed in bytes to perform factorization.
 *    - expansions (int)
 *      Number of memory expansions during the LU factorization.
 */
// A.Duran 8/14/2002
template <class Field>
int FQuerySpace(SuperMatrix<Field> *L, SuperMatrix<Field> *U, int panel_size,
	        mem_usage_t *mem_usage, Field& F)
{
    SCformat<Field> *Lstore;
    NCformat<Field> *Ustore;
    register int n, iword, dword;

    Lstore = (SCformat<Field> *)L->Store;
    Ustore = (NCformat<Field> *)U->Store;
    n = L->ncol;
    iword = sizeof(int);
    dword = sizeof(typename Field::Element);

    /* For LU factors */
    mem_usage->for_lu = (float)( (4*n + 3) * iword + Lstore->nzval_colptr[n] *
				 dword + Lstore->rowind_colptr[n] * iword );
    mem_usage->for_lu += (float)( (n + 1) * iword +
				 Ustore->colptr[n] * (dword + iword) );

    /* Working storage to support factorization */
    mem_usage->total_needed = mem_usage->for_lu +
	(float)( (2 * panel_size + 4 + NO_MARKER) * n * iword +
		(panel_size + 1) * n * dword );

    mem_usage->expansions = --no_expand;

    return 0;
} /* FQuerySpace */


/*
 * Allocate storage for the data structures common to all factor routines.
 * For those unpredictable size, make a guess as FILL * nnz(A).
 * Return value:
 *     If lwork = -1, return the estimated amount of space required, plus n;
 *     otherwise, return the amount of space actually allocated when
 *     memory allocation failure occurred.
 */
// A.Duran 8/14/2002
template <class Field>
int
FLUMemInit(char *refact, void *work, int lwork, int m, int n, int annz,
	  int panel_size, SuperMatrix<Field> *L, SuperMatrix<Field> *U, GlobalLU_t<Field> *Glu,
	  int **iwork, typename Field::Element **dwork, Field& F)
{
    int      info, iword, dword;
    SCformat<Field> *Lstore;
    NCformat<Field> *Ustore;
    int      *xsup, *supno;
    int      *lsub, *xlsub;
    typename Field::Element  *lusup;
    int      *xlusup;
    typename Field::Element  *ucol;
    int      *usub, *xusub;
    int      nzlmax, nzumax, nzlumax;
    int      FILL = sp_ienv(6);
    
    Glu->n    = n;
    no_expand = 0;
    iword     = sizeof(int);
    dword     = sizeof(typename Field::Element);

    if ( !expanders )	
        expanders = (ExpHeader*)SUPERLU_MALLOC(NO_MEMTYPE * sizeof(ExpHeader));
    if ( !expanders ) ABORT("SUPERLU_MALLOC fails for expanders");
    
    if ( lsame_(refact, "N") ) {
	/* Guess for L\U factors */
	nzumax = nzlumax = FILL * annz;
	nzlmax = static_cast<int>( SUPERLU_MAX(1, FILL/4.) * annz ); /* A. Duran */

	if ( lwork == -1 ) {
	    return ( GluIntArray(n) * iword + TempSpace(m, panel_size)
		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword + n );
        } else {
	    FSetupSpace(work, lwork, &Glu->MemModel, F);
	}
	
#ifdef DEBUG		   
	printf("FLUMemInit() called: annz %d, MemModel %d\n", 
		annz, Glu->MemModel);
#endif	
	
	/* Integer pointers for L\U factors */
	if ( Glu->MemModel == SYSTEM ) {
	    xsup   = intMalloc(n+1);
	    supno  = intMalloc(n+1);
	    xlsub  = intMalloc(n+1);
	    xlusup = intMalloc(n+1);
	    xusub  = intMalloc(n+1);
	} else {
	    xsup   = (int *)Fuser_malloc((n+1) * iword, HEAD, F);
	    supno  = (int *)Fuser_malloc((n+1) * iword, HEAD, F);
	    xlsub  = (int *)Fuser_malloc((n+1) * iword, HEAD, F);
	    xlusup = (int *)Fuser_malloc((n+1) * iword, HEAD, F);
	    xusub  = (int *)Fuser_malloc((n+1) * iword, HEAD, F);
	}

	lusup = (typename Field::Element *) Fexpand( &nzlumax, LUSUP, 0, 0, Glu, F );
	ucol  = (typename Field::Element *) Fexpand( &nzumax, UCOL, 0, 0, Glu, F );
	lsub  = (int *)Fexpand( &nzlmax, LSUB, 0, 0, Glu, F );
	usub  = (int *)Fexpand( &nzumax, USUB, 0, 1, Glu, F );

	while ( !lusup || !ucol || !lsub || !usub ) {
	    if ( Glu->MemModel == SYSTEM ) {
		SUPERLU_FREE(lusup); 
		SUPERLU_FREE(ucol); 
		SUPERLU_FREE(lsub); 
		SUPERLU_FREE(usub);
	    } else {
		Fuser_free((nzlumax+nzumax)*dword+(nzlmax+nzumax)*iword, HEAD, F);
	    }
	    nzlumax /= 2;
	    nzumax /= 2;
	    nzlmax /= 2;
	    if ( nzlumax < annz ) {
		printf("Not enough memory to perform factorization.\n");
		return (Fmemory_usage(nzlmax, nzumax, nzlumax, n, F) + n);
	    }
	    lusup = (typename Field::Element *) Fexpand( &nzlumax, LUSUP, 0, 0, Glu, F );
	    ucol  = (typename Field::Element *) Fexpand( &nzumax, UCOL, 0, 0, Glu, F );
	    lsub  = (int *)    Fexpand( &nzlmax, LSUB, 0, 0, Glu, F );
	    usub  = (int *)    Fexpand( &nzumax, USUB, 0, 1, Glu, F );
	}
	
    } else {
	/* refact == 'Y' */
	Lstore   = (SCformat<Field> *)L->Store;
	Ustore   = (NCformat<Field> *)U->Store;
	xsup     = Lstore->sup_to_col;
	supno    = Lstore->col_to_sup;
	xlsub    = Lstore->rowind_colptr;
	xlusup   = Lstore->nzval_colptr;
	xusub    = Ustore->colptr;
	nzlmax   = Glu->nzlmax;    /* max from previous factorization */
	nzumax   = Glu->nzumax;
	nzlumax  = Glu->nzlumax;
	
	if ( lwork == -1 ) {
	    return ( GluIntArray(n) * iword + TempSpace(m, panel_size)
		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword + n );
        } else if ( lwork == 0 ) {
	    Glu->MemModel = SYSTEM;
	} else {
	    Glu->MemModel = USER;
	    _stack.top2 = (lwork/4)*4; /* must be word-addressable */
	    _stack.size = _stack.top2;
	}
	expanders[LSUB].mem  = Lstore->rowind;
	lsub  = (int *)expanders[LSUB].mem;
	// lsub  = (int *)expanders[LSUB].mem  = Lstore->rowind;
	expanders[LUSUP].mem = Lstore->nzval; // A.Duran 8/14/2002
	lusup = (typename Field::Element *)expanders[LUSUP].mem;
	// lusup = (double *)expanders[LUSUP].mem = Lstore->nzval;
	expanders[USUB].mem  = Ustore->rowind; // A.Duran 8/14/2002
	usub  = (int *)expanders[USUB].mem;
	// usub  = (int *)expanders[USUB].mem  = Ustore->rowind;
	expanders[UCOL].mem  = Ustore->nzval;
	ucol  = (typename Field::Element *)expanders[UCOL].mem; // A.Duran 8/14/2002
	// ucol  = (double *)expanders[UCOL].mem  = Ustore->nzval;
	expanders[LSUB].size         = nzlmax;
	expanders[LUSUP].size        = nzlumax;
	expanders[USUB].size         = nzumax;
	expanders[UCOL].size         = nzumax;	
    }

    Glu->xsup    = xsup;
    Glu->supno   = supno;
    Glu->lsub    = lsub;
    Glu->xlsub   = xlsub;
    Glu->lusup   = lusup;
    Glu->xlusup  = xlusup;
    Glu->ucol    = ucol;
    Glu->usub    = usub;
    Glu->xusub   = xusub;
    Glu->nzlmax  = nzlmax;
    Glu->nzumax  = nzumax;
    Glu->nzlumax = nzlumax;
    
    info = FLUWorkInit(m, n, panel_size, iwork, dwork, Glu->MemModel, F);
    if ( info )
	return ( info + Fmemory_usage(nzlmax, nzumax, nzlumax, n, F) + n);
    
    ++no_expand;
    return 0;
    
} /* FLUMemInit */


/* Allocate known working storage. Returns 0 if success, otherwise
   returns the number of bytes allocated so far when failure occurred. */
// A.Duran 8/14/2002
template <class Field>
int
FLUWorkInit(int m, int n, int panel_size, int **iworkptr, 
            typename Field::Element **dworkptr, LU_space_t MemModel, Field& F)
{
    int    isize, dsize, extra;
    typename Field::Element *old_ptr;
    int    maxsuper = sp_ienv(3),
           rowblk   = sp_ienv(4);

    isize = ( (2 * panel_size + 3 + NO_MARKER ) * m + n ) * sizeof(int);
    dsize = (m * panel_size +
	     NUM_TEMPV(m,panel_size,maxsuper,rowblk)) * sizeof(typename Field::Element);
    
    if ( MemModel == SYSTEM ) 
	*iworkptr = (int *) intCalloc(isize/sizeof(int));
    else
	*iworkptr = (int *) Fuser_malloc(isize, TAIL, F);
    if ( ! *iworkptr ) {
	fprintf(stderr, "FLUWorkInit: malloc fails for local iworkptr[]\n");
	return (isize + n);
    }

    if ( MemModel == SYSTEM )
	*dworkptr = (typename Field::Element *) SUPERLU_MALLOC(dsize);
    else {
	*dworkptr = (typename Field::Element *) Fuser_malloc(dsize, TAIL, F);
	if ( NotDoubleAlign(*dworkptr) ) {
	    old_ptr = *dworkptr;
	    *dworkptr = (typename Field::Element*) DoubleAlign(*dworkptr);
	    *dworkptr = (typename Field::Element*) ((double*)*dworkptr - 1);
	    extra = (char*)old_ptr - (char*)*dworkptr;
#ifdef DEBUG	    
	    printf("FLUWorkInit: not aligned, extra %d\n", extra);
#endif	    
	    _stack.top2 -= extra;
	    _stack.used += extra;
	}
    }
    if ( ! *dworkptr ) {
	fprintf(stderr, "malloc fails for local dworkptr[].");
	return (isize + dsize + n);
    }
	
    return 0;
}


/*
 * Set up pointers for real working arrays.
 */
//A.Duran 8/14/2002
template <class Field>
void
FSetRWork(int m, int panel_size, typename Field::Element *dworkptr,
	 typename Field::Element **dense, typename Field::Element **tempv, Field& F)
{
    typename Field::Element zero;
    F.init(zero, 0);

    int maxsuper = sp_ienv(3),
        rowblk   = sp_ienv(4);
    *dense = dworkptr;
    *tempv = *dense + panel_size*m;
    Ffill (*dense, m * panel_size, zero, F);
    Ffill (*tempv, NUM_TEMPV(m,panel_size,maxsuper,rowblk), zero, F);     
}
	
/*
 * Free the working storage used by factor routines.
 */
// A.Duran 8/14/2002
template <class Field>
void FLUWorkFree(int *iwork, typename Field::Element *dwork, GlobalLU_t<Field> *Glu, Field& F)
{
    if ( Glu->MemModel == SYSTEM ) {
	SUPERLU_FREE (iwork);
	SUPERLU_FREE (dwork);
    } else {
	_stack.used -= (_stack.size - _stack.top2);
	_stack.top2 = _stack.size;
/*	dStackCompress(Glu);  */
    }
    
    SUPERLU_FREE (expanders);	
    expanders = 0;
}


/* Expand the data structures for L and U during the factorization.
 * Return value:   0 - successful return
 *               > 0 - number of bytes allocated when run out of space
 */

// A.Duran 8/14/2002
template <class Field>
int
FLUMemXpand(int jcol,
	    int next,          /* number of elements currently in the factors */
	    MemType mem_type,  /* which type of memory to expand  */
	    int *maxlen,       /* modified - maximum length of a data structure */
	    GlobalLU_t<Field> *Glu,    /* modified - global LU data structures */
	    Field& F
	    )
{
  void   *new_mem;
  
#ifdef DEBUG    
  printf("FLUMemXpand(): jcol %d, next %d, maxlen %d, MemType %d\n",
	 jcol, next, *maxlen, mem_type);
#endif    
  
  if (mem_type == USUB) 
    new_mem = Fexpand(maxlen, mem_type, next, 1, Glu, F);
  else
    new_mem = Fexpand(maxlen, mem_type, next, 0, Glu, F);
  
  if ( !new_mem ) {
    int    nzlmax  = Glu->nzlmax;
    int    nzumax  = Glu->nzumax;
    int    nzlumax = Glu->nzlumax;
    fprintf(stderr, "Can't expand MemType %d: jcol %d\n", mem_type, jcol);
    return (Fmemory_usage(nzlmax, nzumax, nzlumax, Glu->n, F) + Glu->n);
  }
  
  switch ( mem_type ) {
  case LUSUP:
    Glu->lusup   = (typename Field::Element *) new_mem;
    Glu->nzlumax = *maxlen;
    break;
  case UCOL:
    Glu->ucol   = (typename Field::Element *) new_mem;
    Glu->nzumax = *maxlen;
    break;
  case LSUB:
    Glu->lsub   = (int *) new_mem;
    Glu->nzlmax = *maxlen;
    break;
  case USUB:
    Glu->usub   = (int *) new_mem;
    Glu->nzumax = *maxlen;
    break;
  }
  
  return 0;
  
}

template <class Field>
void
copy_mem_Field(int howmany, void *old, void *new1, Field& F)
{
    register int i;
    typename Field::Element *dold = (typename Field::Element *)old;
    typename Field::Element *dnew = (typename Field::Element *)new1;
    for (i = 0; i < howmany; i++) dnew[i] = dold[i];
}


/*
 * Expand the existing storage to accommodate more fill-ins.
 */

template <class Field>
void
*Fexpand (
	 int *prev_len,   /* length used from previous call */
	 MemType type,    /* which part of the memory to expand */
	 int len_to_copy, /* size of the memory to be copied to new store */
	 int keep_prev,   /* = 1: use prev_len;
			     = 0: compute new_len to expand */
	 GlobalLU_t<Field> *Glu, /* modified - global LU data structures */
	 Field& F
	)
{
  float EXPAND = 1.5; // ??
  float  alpha; // ??
  void     *new_mem, *old_mem;
  int      new_len, tries, lword, extra, bytes_to_copy;

  alpha = EXPAND;
  
  if ( no_expand == 0 || keep_prev ) /* First time allocate requested */
    new_len = *prev_len;
  else {
    new_len = static_cast<int>( alpha * *prev_len ); /* A. Duran */
  }
  
  if ( type == LSUB || type == USUB ) lword = sizeof(int);
  else lword = sizeof(typename Field::Element);
  
  if ( Glu->MemModel == SYSTEM ) {
    new_mem = (void *) SUPERLU_MALLOC(new_len * lword);
    /*	new_mem = (void *) calloc(new_len, lword); */
    if ( no_expand != 0 ) {
      tries = 0;
      if ( keep_prev ) {
	if ( !new_mem ) return (NULL);
      } else {
	while ( !new_mem ) {
	  if ( ++tries > 10 ) return (NULL);
	  alpha = Reduce(alpha);
	  new_len = static_cast<int>(alpha * *prev_len); /* A. Duran */
	  new_mem = (void *) SUPERLU_MALLOC(new_len * lword); 
	  /*		    new_mem = (void *) calloc(new_len, lword); */
	}
      }
      if ( type == LSUB || type == USUB ) {
	copy_mem_int(len_to_copy, expanders[type].mem, new_mem);
      } else {
	copy_mem_Field(len_to_copy, expanders[type].mem, new_mem, F);
      }
      SUPERLU_FREE (expanders[type].mem);
    }
    expanders[type].mem = (void *) new_mem;
    
  } else { /* MemModel == USER */
    if ( no_expand == 0 ) {
      new_mem = Fuser_malloc(new_len * lword, HEAD, F);
      if ( NotDoubleAlign(new_mem) &&
	   (type == LUSUP || type == UCOL) ) {
	old_mem = new_mem;
	new_mem = (void *)DoubleAlign(new_mem);
	extra = (char*)new_mem - (char*)old_mem;
#ifdef DEBUG		
	printf("expand(): not aligned, extra %d\n", extra);
#endif		
	_stack.top1 += extra;
	_stack.used += extra;
      }
      expanders[type].mem = (void *) new_mem;
    }
    else {
      tries = 0;
      extra = (new_len - *prev_len) * lword;
      if ( keep_prev ) {
	if ( StackFull(extra) ) return (NULL);
      } else {
	while ( StackFull(extra) ) {
	  if ( ++tries > 10 ) return (NULL);
	  alpha = Reduce(alpha);
	  new_len = static_cast<int>(alpha * *prev_len); /* A. Duran */
	  extra = (new_len - *prev_len) * lword;	    
	}
      }
      
      if ( type != USUB ) {
	new_mem = (void*)((char*)expanders[type + 1].mem + extra);
	bytes_to_copy = (char*)_stack.array + _stack.top1
	  - (char*)expanders[type + 1].mem;
	user_bcopy((char *)expanders[type+1].mem, (char *)new_mem, bytes_to_copy);
	
	if ( type < USUB ) {
	  Glu->usub = (int *)expanders[USUB].mem = (void*)((char*)expanders[USUB].mem + extra);
	}
	if ( type < LSUB ) {
	  Glu->lsub = (int *)expanders[LSUB].mem = (void *)((char*)expanders[LSUB].mem + extra);
	}
	if ( type < UCOL ) {
	  expanders[UCOL].mem =
	  Glu->ucol = (typename Field::Element *)expanders[UCOL].mem = (void*)((char*)expanders[UCOL].mem + extra);
	}
	_stack.top1 += extra;
		_stack.used += extra;
		if ( type == UCOL ) {
		  _stack.top1 += extra;   /* Add same amount for USUB */
		  _stack.used += extra;
		}
		
      } /* if ... */
      
    } /* else ... */
  }
  
  expanders[type].size = new_len;
  *prev_len = new_len;
  if ( no_expand ) ++no_expand;
  
  return (void *) expanders[type].mem;
  
} /* Fexpand */


/*
 * Compress the work[] array to remove fragmentation.
 */
template <class Field>
void
FStackCompress(GlobalLU_t<Field> *Glu, Field& F) // ?
{
    register int iword, dword, bytes_to_copy, ndim;
    char    *last, *fragment;
    char     *src, *dest;
    int      *ifrom, *ito;
    typename Field::Element   *dfrom, *dto;
    int      *xlsub, *lsub, *xusub, *usub, *xlusup;
    typename Field::Element  *ucol, *lusup;
    
    iword = sizeof(int);
    dword = sizeof(typename Field::Element);
    ndim = Glu->n;

    xlsub  = Glu->xlsub;
    lsub   = Glu->lsub;
    xusub  = Glu->xusub;
    usub   = Glu->usub;
    xlusup = Glu->xlusup;
    ucol   = Glu->ucol;
    lusup  = Glu->lusup;
    
    dfrom = ucol;
    dto = (typename Field::Element *)((char*)lusup + xlusup[ndim] * dword);
    copy_mem_Field(xusub[ndim], dfrom, dto, F);
    ucol = dto;

    ifrom = lsub;
    ito = (int *) ((char*)ucol + xusub[ndim] * iword);
    copy_mem_int(xlsub[ndim], ifrom, ito);
    lsub = ito;
    
    ifrom = usub;
    ito = (int *) ((char*)lsub + xlsub[ndim] * iword);
    copy_mem_int(xusub[ndim], ifrom, ito);
    usub = ito;
    
    last = (char*)usub + xusub[ndim] * iword;
    fragment = (char*) (((char*)_stack.array + _stack.top1) - last);
    _stack.used -= (long int) fragment;
    _stack.top1 -= (long int) fragment;

    Glu->ucol = ucol;
    Glu->lsub = lsub;
    Glu->usub = usub;
    
#ifdef DEBUG
    printf("FStackCompress: fragment %d\n", fragment);
    /* for (last = 0; last < ndim; ++last)
	print_lu_col("After compress:", last, 0);*/
#endif    
    
}

/*
 * Allocate storage for original matrix A
 */
template <class Field>
void
FallocateA(int n, int nnz, typename Field::Element **a, int **asub, int **xa, Field& F)
{
    *a    = (typename Field::Element *) FieldMalloc(nnz, F);
    *asub = (int *) intMalloc(nnz);
    *xa   = (int *) intMalloc(n+1);
}

template <class Field>
typename Field::Element *FieldMalloc(int n, Field& F)
{
  typename Field::Element *buf;
  buf = (typename Field::Element *) SUPERLU_MALLOC(n * sizeof(typename Field::Element )); 
  if ( !buf ) {
    ABORT("SUPERLU_MALLOC failed for buf in FieldMalloc()\n");
  }
  return (buf);
}


template <class Field>
typename Field::Element *FieldCalloc(int n, Field& F)
{
    typename Field::Element *buf;
    register int i;
    typename Field::Element zero; // A. Duran
    F.init(zero, 0);
    buf = (typename Field::Element *) SUPERLU_MALLOC(n * sizeof(typename Field::Element));
    if ( !buf ) {
	ABORT("SUPERLU_MALLOC failed for buf in FieldCalloc()\n");
    }
    for (i = 0; i < n; ++i) buf[i] = zero;
    return (buf);
}

template <class Field>
int Fmemory_usage(const int nzlmax, const int nzumax, 
		  const int nzlumax, const int n, Field& F)
{
    register int iword, dword;

    iword   = sizeof(int);
    dword   = sizeof(typename Field::Element);
    
    return (10 * n * iword +
	    nzlmax * iword + nzumax * (iword + dword) + nzlumax * dword);

}
