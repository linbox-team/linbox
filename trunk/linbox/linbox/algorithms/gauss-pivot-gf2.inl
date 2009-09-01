// =================================================================== //
// SparseElimination search for pivots over GF2
// Time-stamp: <01 Sep 09 15:32:35 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef __GAUSS_PIVOT_GF2_INL
#define __GAUSS_PIVOT_GF2_INL

namespace LinBox 
{
    template <class Vector, class D> inline void 
    GaussDomain<GF2>::SparseFindPivotBinary (Vector        	&lignepivot,
                                             unsigned long 	&indcol,
                                             long 		&indpermut,
                                             D             	&columns,
                                             bool		&determinant)
    {
 
//        std::cerr << "SFP BEG : lignepivot: [";
//         for(typename Vector::const_iterator refs =  lignepivot.begin();
//             refs != lignepivot.end() ;
//             ++refs )
//             std::cerr << '(' << refs->first << ';' << refs->second << ')';
//         std::cerr << "]" << std::endl;
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
                                             bool& determinant)
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

#endif // __GAUSS_PIVOT_GF2_INL
