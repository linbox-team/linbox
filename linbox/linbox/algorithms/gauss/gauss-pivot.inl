/* linbox/algorithms/gauss-pivot.inl
 * Copyright (C) 2009,2010 The LinBox group
 * Written by JG Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
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
 *.
 *
 * SparseElimination search for pivots
 */

#ifndef __LINBOX_gauss_pivot_INL
#define __LINBOX_gauss_pivot_INL

namespace LinBox
{

	template <class _Field>
	template <class Vector, class D> inline void
	GaussDomain<_Field>::SparseFindPivot (Vector        	&lignepivot,
					      unsigned long 	&indcol,
					      long 		&indpermut,
					      D             	&columns,
					      Element		&determinant) const
	{

		//        std::cerr << "SFP BEG : lignepivot: [";
		//         for(typename Vector::const_iterator refs =  lignepivot.begin();
		//             refs != lignepivot.end() ;
		//             ++refs )
		//             std::cerr << '(' << refs->first << ';' << refs->second << ')';
		//         std::cerr << "]" << std::endl;
		typedef typename Vector::value_type E;

		long nj = (long) lignepivot.size ();


		if (nj > 0) {
			bool pivoting = false;
			indpermut = (long)lignepivot[0].first;

			long ds = (long) --columns[(size_t)indpermut], p = 0;

			for (long j = 1; j < nj; ++j) {
				long dl;
				if ((dl =(long) --columns[lignepivot[(size_t)j].first]) < ds) {
					ds = dl;
					p = j;
				}
			}

			if (p != 0) {
				pivoting = true;
				if (indpermut == static_cast<long>(indcol)) {
					indpermut =(long) lignepivot[(size_t)p].first;
					std::swap( lignepivot[(size_t)p].second, lignepivot[0].second);
				}
				else {
					E ttm = lignepivot[(size_t)p];
					indpermut = (long)ttm.first;

					for (long m = p; m; --m)
						lignepivot[(size_t)m] = lignepivot[(size_t)m-1];

					lignepivot[0] = ttm;
				}
			}

			field().mulin(determinant, lignepivot[0].second);
			if (indpermut != static_cast<long>(indcol)) {
				// std::cerr << "Permuting col: " << indpermut << " <--> " << indcol << std::endl;
				// no need to decrement/increment, already done during the search
				lignepivot[0].first = (unsigned)indcol;
				pivoting = true;
			}

			if (pivoting) field().negin(determinant);
			++indcol;
		}
		else
			indpermut = -1;

		//        std::cerr << "SFP END : lignepivot: [";
		//         for(typename Vector::const_iterator refs =  lignepivot.begin();
		//             refs != lignepivot.end() ;
		//             ++refs )
		//             std::cerr << '(' << refs->first << ';' << refs->second << ')';
		//         std::cerr << "]" << std::endl;
	}


	template <class _Field>
	template <class Vector> inline void
	GaussDomain<_Field>::SparseFindPivot (Vector &lignepivot,
					      unsigned long &indcol,
					      long &indpermut,
					      Element& determinant) const
	{
		long nj = (long)lignepivot.size ();

		if (nj > 0) {
			indpermut = (long) lignepivot[0].first;
			field().mulin(determinant, lignepivot[0].second);
			if (indpermut != static_cast<long>(indcol)){
				// std::cerr << "Permuting col: " << lignepivot[0].first << " <--> " << indcol << std::endl;
				lignepivot[0].first = (unsigned)indcol;
				field().negin(determinant);
			}
			++indcol;
		}
		else
			indpermut = -1;
	}

	template <class _Field>
	template <class Vector> inline void
	GaussDomain<_Field>::FindPivot (Vector &lignepivot,
					unsigned long &k,
					long &indpermut) const
	{
		// Dense lignepivot
		long n = lignepivot.size ();
		long j = (long)k;

		for (; j < n ; ++j )
			if (!field().isZero (lignepivot[(size_t)j])) break ;

		if (j == n )
			indpermut = -1 ;
		else {
			indpermut = j ;
			if (indpermut != (long)k) {
				typename Vector::value_type tmp = lignepivot[(size_t)k] ;
				lignepivot[(size_t)k] = lignepivot[(size_t)j] ;
				lignepivot[(size_t)j] = tmp ;
			}

			++k;
		}
	}

} // namespace LinBox

#endif // __LINBOX_gauss_pivot_INL


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

