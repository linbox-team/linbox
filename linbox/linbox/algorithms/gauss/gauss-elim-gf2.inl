/* linbox/algorithms/gauss-elim-gf2.inl
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
 * SparseElimination elimination routines over GF2
 */
#ifndef __LINBOX_gauss_elim_gf2_INL
#define __LINBOX_gauss_elim_gf2_INL

namespace LinBox
{
	template <class Vector> inline void
	GaussDomain<GF2>::permuteBinary (Vector              &lignecourante,
					 const unsigned long &indcol,
					 const long &indpermut) const
	{
		const unsigned long k = indcol - 1;

#if 0
		std::cerr << "B PERMUTE: " << indpermut << " <--> " << k << " of  [";
		for(typename Vector::const_iterator refs =  lignecourante.begin();
		    refs != lignecourante.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << "]" << std::endl;
#endif

		// precondition indpermut != k
		if (lignecourante.size () ) {
			typename Vector::iterator kin_it=lignecourante.begin();
			for (; kin_it != lignecourante.end(); ++kin_it)
				if (*kin_it >= k) break;
			if (kin_it !=  lignecourante.end()) {
				typename Vector::iterator pin_it=kin_it;
				for (; pin_it != lignecourante.end(); ++pin_it)
					if (static_cast<long>(*pin_it) >= indpermut) break;
				if ( *kin_it == k) {
					if (pin_it != lignecourante.end()) {
						if ( static_cast<long>(*pin_it) != indpermut) {
							--pin_it;
							// Only k there
							*kin_it = (size_t)indpermut;
							typename Vector::value_type etmp = *kin_it;
							typename Vector::iterator current = kin_it;
							typename Vector::iterator next = kin_it; ++next;
							for( ; current != pin_it; ++current, ++next)
								*current = *next;
							*pin_it = etmp;
						}
					}
					else {
						--pin_it;
						// Only k there
						*kin_it = (size_t)indpermut;
						typename Vector::value_type etmp = *kin_it;
						typename Vector::iterator current = kin_it;
						typename Vector::iterator next = kin_it; ++next;
						for( ; current != pin_it; ++current, ++next)
							*current = *next;
						*pin_it = etmp;
					}
				}
				else {
					if (pin_it != lignecourante.end()) {
						if ( static_cast<long>(*pin_it) == indpermut) {
							// Only indpermut there
							*pin_it = k;
							typename Vector::value_type etmp = *pin_it;
							typename Vector::iterator current = pin_it;
							typename Vector::iterator prev = pin_it; --prev;
							for( ; current != kin_it; --current, --prev)
								*current = *prev;
							*kin_it = etmp;
						} // else Nobody
					} // else Nobody
				}
			} // else rien de supérieur à k dans l
			// donc rien à permuter
		}

#if 0
		std::cerr << "E PERMUTE: " << indpermut << " <--> " << k << " of  [";
		for(typename Vector::const_iterator refs =  lignecourante.begin();
		    refs != lignecourante.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << "]" << std::endl;
#endif
	}


	template <class Vector, class D> inline void
	GaussDomain<GF2>::eliminateBinary (bool             &headpivot,
					   Vector              &lignecourante,
					   const Vector        &lignepivot,
					   const unsigned long indcol,
					   const long indpermut,
					   const unsigned long npiv,
					   D                   &columns) const
	{

		typedef typename Vector::value_type E;

		unsigned long k = indcol - 1;
		unsigned long nj = lignecourante.size () ;
#if 0
		std::cerr << "BEGIN ELIMINATE, k: " << k << ", nj: " << nj << ", indpermut: " << indpermut << ", indcol: " << indcol << std::endl;
		std::cerr << "lignepivot: [";
		for(typename Vector::const_iterator refs =  lignepivot.begin();
		    refs != lignepivot.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << "], lignecour: [";
		for(typename Vector::const_iterator refs =  lignecourante.begin();
		    refs != lignecourante.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << ']' << std::endl;
#endif
		if (nj > 0) {
			unsigned long j_head = 0;

			for (; j_head < nj; ++j_head) {
				if (static_cast<long>(lignecourante[(size_t)j_head]) >= indpermut) break;
#if 0
				std::cerr << "ELIMINATE, j_head: " << j_head << std::endl;
#endif
			}

			if (j_head < nj) {
				if (static_cast<long>(lignecourante[(size_t)j_head]) == indpermut) {
					// -------------------------------------------
					// Permutation
					if ( indpermut != static_cast<long>(k)) {
						if (lignecourante[0] != k) {
							// zero <--> non zero
							E tmp = lignecourante[(size_t)j_head];
							--columns[tmp];
							++columns[k];
							tmp = k;

							for (long l = (long)j_head; l > 0; l--)
								lignecourante[(size_t)l] = lignecourante[(size_t)l-1];

							lignecourante[0] = tmp;
						}
						j_head = 0;
					}
					// -------------------------------------------
					// Elimination
					Vector construit (nj + npiv);

					// construit : <-- j
					// courante  : <-- m
					// pivot     : <-- l
					unsigned long j = 0;
					unsigned long m = j_head + 1;

					// A[i,k] <-- - A[i,k] / A[k,k]
					headpivot = true;
					--columns[lignecourante[(size_t)j_head] ];

					// if A[k,j]=0, then A[i,j] <-- A[i,j]
					while (j < j_head) {
						construit[j] = lignecourante[(size_t)j];
						++j;
					}

					unsigned long j_piv;

					unsigned long l = 0;

					for (; l < npiv; ++l)
						if (lignepivot[l] > k) break;

					// for all j such that (j>k) and A[k,j]!=0
					while (l < npiv) {
						j_piv = lignepivot[l];

						// if A[k,j]=0, then A[i,j] <-- A[i,j]
						while ((m < nj) && (lignecourante[(size_t)m] < j_piv))
							construit[j++] = lignecourante[(size_t)m++];

						// if A[i,j]!=0, then A[i,j] <-- A[i,j] - A[i,k]*A[k,j]
						if ((m < nj) && (lignecourante[(size_t)m] == j_piv)) {
							--columns[lignecourante[(size_t)m++]];
						}
						else {
							++columns[j_piv];
							construit[j++] = E (j_piv);
						}

						++l;
					}

					// if A[k,j]=0, then A[i,j] <-- A[i,j]
					while (m<nj)
						construit[j++] = lignecourante[(size_t)m++];

					construit.resize (j);
					lignecourante = construit;
				}
				else {
					// -------------------------------------------
					// j_head < nj but nothing under the pivot
					// Permutation
#if 0
					std::cerr << "----------------------------------------------------------" << std::endl;
					std::cerr << "j_head < nj" << std::endl;
					std::cerr << "j_head: " << j_head << ", nj: " << nj << ", k:" << k
					// << "lignepivot: " << lignepivot
					// << ", lignecour: " << lignecourante
					<< std::endl;
					std::cerr << "----------------------------------------------------------" << std::endl;
#endif
					if (indpermut != static_cast<long>(k)) {
						if (j_head>0) {
							unsigned long l = 0;

							for (; l < nj; ++l)
								if (lignecourante[(size_t)l] >= k) break;

							if ((l < nj) && (lignecourante[(size_t)l] == k))  {
								// non zero <--> zero
								E tmp ; // = (E)lignecourante[(size_t)l];
								--columns[k];
								++columns[(size_t)indpermut];
								tmp = (E)indpermut;

								unsigned long bjh = j_head-1;
								for (; l < bjh; ++l)
									lignecourante[(size_t)l] = lignecourante[(size_t)l + 1];

								lignecourante[bjh] = tmp;
							} // else // zero <--> zero
						} // else // zero <--> zero
					}
				}
			}
			else {
				// -------------------------------------------
				// j_head >= nj > 0
#if 0
				std::cerr << "----------------------------------------------------------" << std::endl;
				std::cerr << "j_head >= nj > 0" << std::endl;
				std::cerr << "j_head: " << j_head << ", nj: " << nj << ", k:" << k
				// << "lignepivot: " << lignepivot
				// << ", lignecour: " << lignecourante
				<< std::endl;
				std::cerr << "----------------------------------------------------------" << std::endl;
#endif
				if (indpermut != static_cast<long>(k)) {
					unsigned long l = 0;

					for (; l < nj; ++l)
						if (lignecourante[(size_t)l] >= k) break;

					if ((l < nj) && (lignecourante[(size_t)l] == k))  {
						// non zero <--> zero
						E tmp ; // = (E) lignecourante[(size_t)l];
						--columns[k];
						++columns[(size_t)indpermut];
						tmp = (E) indpermut;

						unsigned long bjh = nj - 1;
						for (; l < bjh; ++l)
							lignecourante[(size_t)l] = lignecourante[(size_t)l + 1];

						lignecourante[bjh] = tmp;
					} // else
					// zero <--> zero
				}

			}
		}


#if 0
		std::cerr << "END ELIMINATE, k: " << k << ", nj: " << nj << ", indpermut: " << indpermut << ", indcol: " << indcol << std::endl;
		std::cerr << "lignepivot: [";
		for(typename Vector::const_iterator refs =  lignepivot.begin();
		    refs != lignepivot.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << "], lignecour: [";
		for(typename Vector::const_iterator refs =  lignecourante.begin();
		    refs != lignecourante.end() ;
		    ++refs )
			std::cerr << '(' << refs->first << ';' << refs->second << ')';
		std::cerr << ']' << std::endl;
#endif

	}






} // namespace LinBox

#endif // __LINBOX_gauss_elim_gf2_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
