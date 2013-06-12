/* linbox/algorithms/gauss-pivot-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <21 Jan 10 15:08:59 Jean-Guillaume.Dumas@imag.fr>
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
 * SparseElimination search for pivots over GF2
 */
#ifndef __LINBOX_gauss_pivot_gf2_INL
#define __LINBOX_gauss_pivot_gf2_INL

namespace LinBox
{
	template <class Vector, class D> inline void
	GaussDomain<GF2>::SparseFindPivotBinary (Vector        	&lignepivot,
						 unsigned long 	&indcol,
						 long 		&indpermut,
						 D             	&columns,
						 bool		&) const //determinant
	{

#if 0
		std::cerr << "SFP BEG : lignepivot: [";
		for(typename Vector::const_iterator refs =  lignepivot.begin();
		    refs != lignepivot.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << "]" << std::endl;
#endif
		typedef typename Vector::value_type E;

		long nj =  (long)lignepivot.size ();

		if (nj > 0) {
			indpermut = (long)lignepivot.front();

			long ds = (long) --columns[(size_t)indpermut], dl, p = 0;

			for (long j = 1; j < nj; ++j) {
				if ((dl = (long) --columns[lignepivot[(size_t)j]]) < ds) {
					ds = dl;
					p = j;
				}
			}

			if (p != 0) {
				if (indpermut == static_cast<long>(indcol)) {
					indpermut = (long)lignepivot[(size_t)p];
				}
				else {
					E ttm = (E)lignepivot[(size_t)p];
					indpermut = (long)ttm;

					for (long m = p; m; --m)
						lignepivot[(size_t)m] = lignepivot[(size_t)m-1];

					lignepivot[0] = ttm;
				}
			}

			if (indpermut != static_cast<long>(indcol)) {
				// std::cerr << "Permuting col: " << indpermut << " <--> " << indcol << std::endl;
				// no need to decrement/increment, already done during the search
				lignepivot[0] = indcol;
			}

			++indcol;
		}
		else
			indpermut = -1;
#if 0
		std::cerr << "SFP END : lignepivot: [";
		for(typename Vector::const_iterator refs =  lignepivot.begin();
		    refs != lignepivot.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << "]" << std::endl;
#endif
	}

	template <class Vector> inline void
	GaussDomain<GF2>::SparseFindPivotBinary (Vector &lignepivot,
						 unsigned long &indcol,
						 long &indpermut,
						 bool& ) const // determinant
	{
		long nj = (long)lignepivot.size ();

		if (nj > 0) {
			indpermut = (long)lignepivot.front();
			if (indpermut != static_cast<long>(indcol)){
				lignepivot.front() = indcol;
			}
			++indcol;
		}
		else
			indpermut = -1;
	}


} // namespace LinBox

#endif // __LINBOX_gauss_pivot_gf2_INL


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

