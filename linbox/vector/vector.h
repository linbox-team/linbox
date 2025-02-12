/* linbox/vector/vector.h
 * Copyright (C) 2015 the LinBox group
 *
 * Written by :
 * Jean-Guillaume Dumas
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file vector/vector.h
 * @ingroup vector
 * A \c DenseVector<\c _Field >  is default dense vector of 
 * <code>_Field::Element</code>s.
 *
 */

#ifndef __LINBOX_vector_dense_vector_H
#define __LINBOX_vector_dense_vector_H

#include "linbox/linbox-config.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox { /* BlasVector */

    template<class _Vector>
    class BlasSubvector ;

    template<class _Field, class _Rep=typename RawVector<typename _Field::Element>::Dense>
    class BlasVector ;


    template <typename _Field> struct DenseVectorChooser { typedef BlasVector<_Field> type; }; // to allow specialization of using DenseVector
    template <typename _Field> using DenseVector = typename DenseVectorChooser<_Field>::type;

    //template <typename _Field> struct DenseSubVectorChooser { typedef BlasSubvector<DenseVector<_Field>> type; }; // to allow specialization of using DenseVector
    template <typename _Field> struct DenseSubVectorChooser { typedef typename DenseVector<_Field>::subVectorType type; }; // to allow specialization of using DenseVector
    template <class _Field> using DenseSubvector = typename DenseSubVectorChooser<_Field>::type; 

}


#endif // __LINBOX_vector_dense_vector_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
