/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain-setup.inl
 * Copyright (C) 2012 Matthew Wezowicz
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_submatrix_adapter_H
#define __LINBOX_submatrix_adapter_H

namespace LinBox{

	/**
	 * Generic submatrix view adapter used internally in the OpenCLMatrixDomain
	 * Feel free to expand its functionality if needed
	 */
	template <class Matrix>
	class SubmatrixAdapter{
	private:
		typedef typename Matrix::Element Element;

		Matrix* _Mat;
		size_t _row;
		size_t _col;
		size_t _r0;
		size_t _c0;

	public:
		SubmatrixAdapter() : _Mat(NULL), _row(0), _col(0), _r0(0), _c0(0){}

		SubmatrixAdapter(const Matrix& M) : _Mat(&(const_cast<Matrix&>(M))), _row(M.rowdim()),
			_col(M.coldim()), _r0(0), _c0(0) {}

		SubmatrixAdapter(const Matrix& M, size_t row, size_t col, size_t Rowdim, size_t Coldim) :
			_Mat(&(const_cast<Matrix&>(M))), _row(Rowdim), _col(Coldim), _r0(row), _c0(col) {}

		size_t rowdim() const{
			return _row;
		}

		size_t coldim() const{
			return _col;
		}

		void setEntry(size_t i, size_t j, const Element& a_ij){
			_Mat->setEntry(_r0 + i, _c0 + j, a_ij);
		}

		const Element& getEntry(size_t i, size_t j) const{
			return _Mat->getEntry(_r0 + i, _c0 + j);
		}
	};

} //end of namespace LinBox

#endif // __LINBOX_submatrix_adapter_H