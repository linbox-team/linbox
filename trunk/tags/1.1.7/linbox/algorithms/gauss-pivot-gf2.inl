/* linbox/algorithms/gauss-pivot-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <21 Jan 10 15:08:59 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
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

	long nj =  lignepivot.size ();

	if (nj > 0) {
            indpermut = lignepivot.front();

            long ds = --columns[indpermut], dl, p = 0;

            for (long j = 1; j < nj; ++j) {
                if ((dl = --columns[lignepivot[j]]) < ds) {
                    ds = dl;
                    p = j;
                }
            }

            if (p != 0) {
                if (indpermut == static_cast<long>(indcol)) {
                    indpermut = lignepivot[p];
                } else {
                    E ttm = lignepivot[p];
                    indpermut = ttm;

                    for (long m = p; m; --m)
                        lignepivot[m] = lignepivot[m-1];

                    lignepivot[0] = ttm;
                }
            }

            if (indpermut != static_cast<long>(indcol)) {
// std::cerr << "Permuting col: " << indpermut << " <--> " << indcol << std::endl;
                    // no need to decrement/increment, already done during the search
                lignepivot[0] = indcol;
            }

            ++indcol;
	} else
            indpermut = -1;

            //        std::cerr << "SFP END : lignepivot: [";
//         for(typename Vector::const_iterator refs =  lignepivot.begin();
//             refs != lignepivot.end() ;
//             ++refs )
//             std::cerr << '(' << refs->first << ';' << refs->second << ')';
//         std::cerr << "]" << std::endl;
    }

    template <class Vector> inline void 
    GaussDomain<GF2>::SparseFindPivotBinary (Vector &lignepivot, 
                                             unsigned long &indcol, 
                                             long &indpermut, 
                                             bool& ) const // determinant
    {
	long nj = lignepivot.size ();

	if (nj > 0) {
            indpermut = lignepivot.front();
            if (indpermut != static_cast<long>(indcol)){
                lignepivot.front() = indcol;
            }
            ++indcol;
	} else
            indpermut = -1;
    }


} // namespace LinBox

#endif // __LINBOX_gauss_pivot_gf2_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
