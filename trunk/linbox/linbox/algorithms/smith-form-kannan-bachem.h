                                                      
/* tests/__LINBOX_smith_form_kannan_bachem.h
 * Copyright (C) 2014 Gavin Harrison,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>,
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
 */

#include "linbox/matrix/dense-matrix.h"

#ifndef __LINBOX_smith_form_kannan_bachem_domain_H
#define __LINBOX_smith_form_kannan_bachem_domain_H

namespace LinBox
{
	template<class MatrixDomain>
	class SmithFormKannanBachemDomain
	{
	public:
		typedef typename MatrixDomain::Field Field;
		typedef typename Field::Element Element;
		typedef BlasMatrix<Field> Rep;

	private:
		MatrixDomain _MD;
		const Field &field() const { return _MD.field(); }

	public:
		SmithFormKannanBachemDomain(const MatrixDomain &MD) : _MD(MD) {}
		SmithFormKannanBachemDomain(const SmithFormKannanBachemDomain &D) : _MD(D._MD) {}

	private:
		void swapRows(Rep &A, int a, int b) const
		{
			for (int i = 0; i < A.coldim(); i++)
			{
				Element tmp1, tmp2;

				A.getEntry(tmp1, a, i);
				A.getEntry(tmp2, b, i);

				A.setEntry(a, i, tmp2);
				A.setEntry(b, i, tmp1);
			}
		}

		void swapCols(Rep &A, int a, int b) const
		{
			for (int i = 0; i < A.rowdim(); i++)
			{
				Element tmp1, tmp2;

				A.getEntry(tmp1, i, a);
				A.getEntry(tmp2, i, b);

				A.setEntry(i, a, tmp2);
				A.setEntry(i, b, tmp1);
			}
		}
		
		// Ensures that if a=b then s=u=v=1 and t=0 to avoid an infinite loop
		void dxgcd(Element &s, Element &t, Element &u, Element &v, const Element &a, const Element &b) const
		{
			if (field().areEqual(a,b))
			{
				field().assign(s, field().one);
				field().assign(t, field().zero);
				field().assign(u, field().one);
				field().assign(v, field().one);
				return;
			}

			Element g;
			field().dxgcd(g,s,t,u,v,a,b);
		}
		
		bool findPivot(Rep &A, int n) const
		{
			for (int i = n; i < A.rowdim(); i++)
			{
				for (int j = n; j < A.coldim(); j++)
				{
					Element tmp;
					A.getEntry(tmp, i, j);

					if (!field().isZero(tmp))
					{
						if (i != n)
							swapRows(A, n, i);

						if (j != n)
							swapCols(A, n, j);
						
						return true;
					}
				}
			}
			
			return false;
		}
		
		void eliminateCol(Rep &A, int p, int o) const
		{
			Element pp, op;
			
			A.getEntry(op, o, p);
			
			if (field().isZero(op))
				return;
			
			A.getEntry(pp, p, p);
			
			Element s,t,u,v;
			dxgcd(s,t,u,v,pp,op);
			
			for (int i = p; i < A.coldim(); i++)
			{
				Element pi,oi;
				
				A.getEntry(pi, p, i);
				A.getEntry(oi, o, i);
				
				Element tmp1, tmp2;
				field().mul(tmp1, s, pi);
				field().mul(tmp2, t, oi);
				field().addin(tmp1, tmp2);
				A.setEntry(p, i, tmp1);
				
				field().mul(tmp1, pi, v);
				field().mul(tmp2, oi, u);
				field().subin(tmp2, tmp1);
				A.setEntry(o, i, tmp2);
			}
		}
		
		void eliminateRow(Rep &A, int p, int o) const
		{
			Element pp, po;
			
			A.getEntry(po, p, o);
			
			if (field().isZero(po))
				return;
			
			A.getEntry(pp, p, p);
			
			Element s,t,u,v;
			dxgcd(s,t,u,v,pp,po);
			
			for (int i = p; i < A.rowdim(); i++)
			{
				Element ip, io;
				
				A.getEntry(ip, i, p);
				A.getEntry(io, i, o);
				
				Element tmp1, tmp2;
				
				field().mul(tmp1, ip, s);
				field().mul(tmp2, io, t);
				field().addin(tmp1, tmp2);
				A.setEntry(i, p, tmp1);
				
				field().mul(tmp1, ip, v);
				field().mul(tmp2, io, u);
				field().subin(tmp2, tmp1);
				A.setEntry(i, o, tmp2);
			}
		}
		
		void reduceOffDiagonal(Rep &A, int s, int e) const
		{
			for (int i = s; i <= e; i++)
			{
				Element nii,ii;
				A.getEntry(ii, i, i);
				field().normalize(nii,ii);
				
				Element tmp;
				field().div(tmp, nii, ii);
				
				if (field().isOne(tmp))
					continue;
				
				// A[i] = A[i,i] * n
				// where n = normalized(A[i,i]) / A[i,i]
				A.setEntry(i, i, nii);
				for (size_t j = i+1; j < A.coldim(); j++)
				{
					Element ij;
					A.getEntry(ij, i, j);
					field().mulin(ij, tmp);
					A.setEntry(i, j, ij);
				}
			}
			
			// Ording of reduction here is an improvement to Kannan/Bachem
			// Introduced by Chou/Collins '82
			// Reduce from bottom to top and left to right
			// * 4 5 6
			// 0 * 2 3
			// 0 0 * 1
			// 0 0 0 *
			for (int i = e-1; i >= s; i--)
			{
				for (int j = i+1; j <= e; j++)
				{
					Element jj, ij, tmp;
					
					A.getEntry(ij, i, j);
					if (field().isZero(ij))
						continue;
					
					A.getEntry(jj, j, j);
					field().quo(tmp, ij, jj);
					
					// A[i] = A[i] - quo(A[i,j], A[j,j]) * A[j]
					for (size_t k = j; k < A.coldim(); k++)
					{
						Element ik, jk;
						
						A.getEntry(ik, i, k);
						A.getEntry(jk, j, k);
						
						field().mulin(jk, tmp);
						field().subin(ik, jk);
						
						A.setEntry(i, k, ik);
					}
				}
			}
		}
		
		// Puts the lower-right n-by-n minor of A into Hermite Normal Form
		void hermite(Rep &A, int n) const
		{
			int dim = A.rowdim();
			
			for (int i = n; i < dim; i++)
			{
				for (int j = n; j < i; j++)
					eliminateCol(A, j, i);
				
				if (!findPivot(A, i))
					return;
				
				reduceOffDiagonal(A, n, i);
			}
		}
		
		bool isRowDiagonalized(const Rep &A, int n) const
		{
			for (int i = n+1; i < A.coldim(); i++)
			{
				Element ni;
				A.getEntry(ni, n, i);
				if (!field().isZero(ni))
					return false;
			}
			return true;
		}
		
		bool pivotDividesRemaining(Rep &A, int n) const
		{
			Element nn;
			A.getEntry(nn, n, n);
			
			for (int i = n+1; i < A.rowdim(); i++)
			{
				for (int j = i; j < A.coldim(); j++)
				{
					Element ij, g;
					A.getEntry(ij, i, j);
					
					if (field().isZero(ij))
						continue;
					
					field().gcd(g, nn, ij);
					
					if (!field().areAssociates(g, nn))
					{
						// Add row i to row n
						for (int k = i; k < A.coldim(); k++)
						{
							Element ik;
							A.getEntry(ik, i, k);
							A.setEntry(n, k, ik);
						}
						return false;
					}
				}
			}
			
			return true;
		}
		
	public:
		template<class Vector>
		Vector &solve(Vector &S, const Rep &A) const
		{
			int dim = A.rowdim();
			linbox_check(A.coldim() == dim && S.size() >= dim);
			
			Rep B(A);
			
			for (int i = 0; i < dim;)
			{
				if (!findPivot(B, i))
					break;
				
				for (int j = i+1; j < dim; j++)
					eliminateRow(B, i, j);
				
				hermite(B, i);
				
				if (!isRowDiagonalized(B, i))
					continue;
				
				if (!pivotDividesRemaining(B, i))
					continue;
				
				i++;
			}
			
			for (int i = 0; i < dim; i++)
			{
				Element ii;
				B.getEntry(ii, i, i);
				S.setEntry(i, ii);
			}
			
			return S;
		}
	};
}

#endif // __LINBOX_smith_form_kannan_bachem_domain_H
