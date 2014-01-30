/*
 * Copyright (C) the LinBox group
 *
 * Written by BB <bbboyer@ncsu.edu>
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
 */

/** @file linbox/matrix/SparseMatrix/read-write-sparse.h
 * @ingroup sparsematrix
 * @brief
 */

#ifndef __LINBOX_matrix_sparse_matrix_read_write_sparse_H
#define __LINBOX_matrix_sparse_matrix_read_write_sparse_H

#include "linbox/matrix/matrix-traits.h"
#include "linbox/matrix/matrix-category.h"
#include "linbox/util/matrix-stream.h"

namespace LinBox {
	// Small helper classes to make read and write easier
	template <class Matrix>
	class SparseMatrixWriteHelper {
	public:
		typedef typename Matrix::Field          Field;
		typedef typename Field::Element       Element;

	private:
		static std::ostream &writeTriple (const Matrix &A, std::ostream &os, bool oneBased = false);

		static std::ostream &writePretty (const Matrix &A, std::ostream &os
						  , std::string begmat
						  , std::string endmat
						  , std::string begrow
						  , std::string endrow
						  , std::string sepelt
						  , std::string seprow
						 );

		static std::ostream& write(const Matrix &A
				   , std::ostream &os
				   , LINBOX_enum(Tag::FileFormat) format
				   , MatrixCategories::BlackboxTag);

		static std::ostream& write(const Matrix &A
				   , std::ostream &os
				   , LINBOX_enum(Tag::FileFormat) format
				   , MatrixCategories::RowMatrixTag);


	public:
		static std::ostream &write (const Matrix &A
					    , std::ostream &os
					    , LINBOX_enum(Tag::FileFormat) format)
		{
			return write(A,os,format, typename MatrixTraits<Matrix>::MatrixCategory ());
		}
	};

	template <class Matrix>
	class SparseMatrixReadWriteHelper : public SparseMatrixWriteHelper<Matrix> {

	public:
		typedef typename Matrix::Field::Element Element;

	private:

		static std::istream &readTurner    (Matrix &A
						    , std::istream &is
						    , char *buf);

		static std::istream &readGuillaume (Matrix &A
						    , std::istream &is
						    , char *buf);

		static std::istream &readMatlab    (Matrix &A
						    , std::istream &is
						    , char *buf);

		static std::istream &readPretty    (Matrix &A
						    , std::istream &is
						    , char *buf);

		static std::istream &readMagmaCpt  (Matrix &A
						    , std::istream &is
						    , char *buf);

		static std::istream &readMatrixMarket (Matrix &A
						    , std::istream &is
						    , char *buf);

		static std::istream &read (Matrix &A
				    , std::istream &is
				    , LINBOX_enum(Tag::FileFormat) format
				    , MatrixCategories::RowMatrixTag
				   );

		static std::istream &read (Matrix &A
				    , std::istream &is
				    , LINBOX_enum(Tag::FileFormat) format
				    , MatrixCategories::BlackboxTag
				   );


	public:
		static std::istream &read (Matrix &A
					   , std::istream &is
					   , LINBOX_enum(Tag::FileFormat) format)
		{
			return read(A,is,format, typename MatrixTraits<Matrix>::MatrixCategory ());
		}
	};


} // LinBox


#include "linbox/matrix/SparseMatrix/read-write-sparse.inl"
#endif // __LINBOX_matrix_sparse_matrix_read_write_sparse_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
