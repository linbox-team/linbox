
/* Copyright (C) 2012 LinBox
 * Written by bds
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_write_mm_h
#define __LINBOX_write_mm_h

/* write-mm.h
 * Tools to help write fields and matrices in matrix market formats.
 *
 */

#include <iostream>
#include <string>
#include "linbox/integer.h"

namespace LinBox
{

/// Write second line and comment part of matrix market header
template <class Field>
std::ostream& writeMMComment(std::ostream& os, Field& F, std::string name, std::string comment) {
	F.write(os << "% written by LinBox::",std::string("F")) << "; ";
	F.write(os << name << "<", std::string("")) << " >(F)" << std::endl;
	if (comment.size() > 0)
		os << "%" << std::endl << "% " << comment << std::endl << "%" << std::endl;
    return os;
} 

/// Write matrix market header (up to the i,j,val lines) for a sparse or structured matrix. 
template <class BB> 
std::ostream& writeMMCoordHeader(std::ostream& os, BB& A, size_t nnz, std::string name, std::string comment = "") {
	os << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	writeMMComment(os, A.field(), name, comment);
	os << A.rowdim() << " " << A.coldim() << " " << nnz << std::endl;
	return os;
}

/// Write matrix market header (up to the i,j,val lines) for a sparse or structured matrix. 
template <class BB> 
std::ostream& writeMMPatternHeader(std::ostream& os, BB& A, size_t nnz, std::string name, std::string comment = "") {
	os << "%%MatrixMarket matrix coordinate pattern general" << std::endl;
	writeMMComment(os, A.field(), name, comment);
	os << A.rowdim() << " " << A.coldim() << " " << nnz << std::endl;
	return os;
}

/// Write matrix market header (up to the entry lines) for a dense matrix. 
template <class BB> 
std::ostream& writeMMArrayHeader(std::ostream& os, BB& A, std::string name, std::string comment = "") {
	os << "%%MatrixMarket matrix array integer general" << std::endl;
	writeMMComment(os, A.field(), name, comment);
	os << A.rowdim() << " " << A.coldim() << std::endl;
	return os;
}

/// Generic dense matrix writer to matrix market array format.
template <class Mat> 
std::ostream& writeMMArray(std::ostream& os, Mat& A, std::string name, std::string comment = "") {
	writeMMArrayHeader(os, A, name, comment);
	typename Mat::Field::Element x; A.field().init(x, 0);
	for (size_t j = 0; j < A.coldim(); ++j)
		for (size_t i = 0; i < A.rowdim(); ++i)
			os << A.getEntry(x, i, j) << std::endl;
	return os;
}

/// eltype(x) returns a string containing the name of the type of field or ring element x.
std::string eltype(float x) { return "float"; }
std::string eltype(double x) { return "double"; }
std::string eltype(int8_t x) { return "int8_t"; }
std::string eltype(int16_t x) { return "int16_t"; }
std::string eltype(int32_t x) { return "int32"; }
std::string eltype(int64_t x) { return "int64"; }
std::string eltype(integer x) { return "integer"; }
std::string eltype(uint8_t x) { return "uint8_t"; }
std::string eltype(uint16_t x) { return "uint16_t"; }
std::string eltype(uint32_t x) { return "uint32"; }
std::string eltype(uint64_t x) { return "uint64"; }

}  // end of namespace LinBox
#endif // __write_mm_h

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

