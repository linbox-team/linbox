/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/const-iterator-type.h
 * Copyright (C) 2002 Zhendong Wan
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>
 * ------------------------------------
 * 2002-08-08  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Rename from constiteratortype.h; reindent
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __CONST_ITERATOR_TYPE_H
#define __CONST_ITERATOR_TYPE_H

#include <vector>

namespace LinBox
{

template<class Iterator>
class Subiterator;

template<class Iterator>
class ConstIteratorType
{
    public:
	typedef Iterator const_iterator;
};
  
template<class Iterator>
class ConstIteratorType<Subiterator<Iterator> >
{
    public:
	typedef Subiterator<typename ConstIteratorType<Iterator>::const_iterator> const_iterator;
};

template<class T>
class ConstIteratorType<T* >
{
    public:
	typedef const T* const_iterator;
};

template<>
class ConstIteratorType<std::vector<bool>::iterator >
{
    public:
	typedef std::vector<bool>::const_iterator const_iterator;
};

}

#include <linbox/vector/subiterator.h>

#endif
