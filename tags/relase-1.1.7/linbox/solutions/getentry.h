/* linbox/solutions/getentry.h
 * Copyright(C) LinBox
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *  -bds
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_getentry_H
#define __LINBOX_getentry_H

#include <vector>

#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

namespace GetEntryTags 
{	struct GenericBB{};	struct Local{};	
	struct SpecialCDB{}; struct SpecialCBD{}; struct SpecialCDD{}; 
}; // namespace GetEntryTags

template<class BB> struct GetEntryCategory { typedef GetEntryTags::GenericBB Tag; };

/** \brief
 * Getting the i,j entry of the blackbox.
 */
template <class BB> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j)
{ 
	typename GetEntryCategory<BB>::Tag t; 
	return getEntry(x, A, i, j, t);
}

/// To ignore methods
template <class BB, class Method> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, Method & m)
{ return getEntry(x, A, i, j);	}

// Generic BBs require use of apply.
template <class BB> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, GetEntryTags::GenericBB t)
{
	typedef typename BB::Field Field;
	typedef typename Field::Element Elt;
	typedef std::vector<Elt> Vector;

	const Field& F = A.field();
	Elt zero; F.init(zero, 0UL);
	Vector v(A.coldim(), zero), w(A.rowdim(), zero);
	for (typename Vector::iterator it = v.begin (); it != v.end (); ++it)
		F.assign (*it, zero);
	F.init(v[j],1UL);
	A.apply (w, v);
	F.assign (x, VectorWrapper::constRef<Field, Vector> (w, i));
	return x;
}

// BBs that offer a local getEntry.
template<class Field> struct GetEntryCategory<DenseMatrix<Field> > { typedef GetEntryTags::Local Tag; };
template<class A, class B> struct GetEntryCategory<SparseMatrix<A,B> > { typedef GetEntryTags::Local Tag; };
template<class A, class B, class C> struct GetEntryCategory<SparseMatrixBase<A,B,C> > { typedef GetEntryTags::Local Tag; };
template<class Field, class Trait> struct GetEntryCategory<Diagonal<Field, Trait> > { typedef GetEntryTags::Local Tag; };
template<class Field> struct GetEntryCategory<ScalarMatrix<Field> > { typedef GetEntryTags::Local Tag; };

template <class BB> 
typename BB::Field::Element& getEntry(typename BB::Field::Element& x, const BB& A, const size_t i, const size_t j, GetEntryTags::Local t )
{ return A.getEntry(x, i, j); }

// Compose< Diagonal, BB > specialization
template <class Field, class Trait, class BB> 
struct GetEntryCategory<Compose<Diagonal<Field, Trait>,BB> > { typedef GetEntryTags::SpecialCDB Tag; };

template <class Field, class Trait, class BB> 
typename Field::Element& getEntry(typename Field::Element& x, const Compose<Diagonal<Field, Trait>, BB>& A, const size_t i, const size_t j, GetEntryTags::SpecialCDB t)
{
    typename Field::Element y;
    getEntry(y, *(A.getLeftPtr()), i, i);
    getEntry(x, *(A.getRightPtr()), i, j);
    return A.field().mulin(x, y);
}

// Compose< BB, Diagonal > specialization
template <class BB, class Field, class Trait> 
struct GetEntryCategory<Compose<BB, Diagonal<Field, Trait> > > { typedef GetEntryTags::SpecialCBD Tag; };

template <class BB, class Field, class Trait> 
typename Field::Element& getEntry(typename Field::Element& x, const Compose<BB, Diagonal<Field, Trait> >& A, const size_t i, const size_t j, GetEntryTags::SpecialCBD t)  
{
    typename Field::Element y;
    getEntry(y, *(A.getLeftPtr()), i, j);
    getEntry(x, *(A.getRightPtr()), j, j);
    return A.field().mulin(x, y);
}

// Compose< Diagonal, Diagonal > specialization
template <class Field, class T1, class T2> 
struct GetEntryCategory<Compose<Diagonal<Field, T1>,Diagonal<Field, T2> > > { typedef GetEntryTags::SpecialCDD Tag; };

template <class Field, class T1, class T2> 
typename Field::Element& getEntry(typename Field::Element& x, const Compose<Diagonal<Field,T1>, Diagonal<Field, T2> >& A, const size_t i, const size_t j, GetEntryTags::SpecialCDD t)
{
    if (i != j) 
        return A.field().init(x, 0UL);
    else {
        typename Field::Element y;
        getEntry(y, *(A.getLeftPtr()), i, i);
        getEntry(x, *(A.getRightPtr()), j, j);
        return A.field().mulin(x, y);
    }
}


}

#endif // __LINBOX_getentry_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
