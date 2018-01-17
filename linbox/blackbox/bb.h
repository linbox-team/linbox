#ifndef LB_BB_H
#define LB_BB_H
/* linbox/blackbox/bb.h
 * Copyright (C) 2015 bds for LinBox Team.  See linbox/COPYING.LESSER for License info.
 *
 * Written by bds
 *
 * Blackbox base class for support of separate compilation
 */

/* blackbox base class
 *
 * Algorithms may take BB<Field> parameters and be separately compiled.
 *
 * Non-template functions of BB are pure virtual.  Code bloat is avoided.
 * Template functions may be called by select on the BBType tag.  This 
 * introduces some code bloat.
 */
#include <iostream>
#include "linbox/util/error.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
namespace LinBox {

// for now, only the FIBB tags.
enum BBType {diagonal, permutation, triangular, product, lqup, pluq, other};

template <class Ring>
struct BB 
{
	typedef Ring Field;
	typedef DenseMatrix<Field> ResizableMatrix;
	//typedef DenseMatrix<Field, std::vector<typename Ring::Element> > ResizableMatrix;
	typedef DenseMatrix<Field> Matrix;

	virtual ~BB(){}

	virtual BBType bbTag() const 
	= 0;
	virtual size_t rowdim() const
	= 0;
	virtual size_t coldim() const
	= 0;
	virtual const Field& field() const
	= 0;
	virtual std::ostream& write(std::ostream& os) const
	= 0;
	virtual std::istream& read(std::istream& os) 
	= 0;
	virtual Matrix& applyLeft(Matrix& Y, const Matrix& X) const
	= 0;
	virtual Matrix& applyRight(Matrix& Y, const Matrix& X) const
	= 0;
	template<class OutVector, class InVector>
	OutVector& apply(OutVector& y, const InVector& x) const
	{ switch (bbTag()) 
	  {	//case BBx_tag: static_cast<BBx<Field>*>(this)->apply(y,x); break;
	  	//case permutation: static_cast<Permutation<Field>*>(this)->apply(y, x)
	  	default: throw LinboxError("indirect call to apply not supported for BBType " /* bbTag*/);
	  } 
	  return y;
	}
	template<class OutVector, class InVector>
	OutVector& applyTranspose(OutVector& y, const InVector& x) const
	{ switch (bbTag()) 
	  {	//case BBxtag: static_cast<BBx<Field>*>(this)->applyTranspose(y,x); break;
	  	default: throw LinboxError("indirect call to applyTranspose not supported for BBType " /* bbTag*/);
	  } 
	  return y;
	}

	template<typename BB2>
	void map(BB2& A)
	{ switch (bbTag()) 
	  {	//case bbxtag: static_cast<bbx<Field>*>(this)->map(A); break;
	     // using it's struct rebind;
	  	default: throw LinboxError("indirect call to map not supported for BBType " /* bbTag*/);
	  }
	  return A;
	}

}; // class BB

} // LinBox
#endif // LB_BB_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
