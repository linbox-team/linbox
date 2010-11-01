/* Copyright (C) LinBox
 *
 *
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */



#ifndef __LINBOX_kaapi_communicate_H
#define __LINBOX_kaapi_communicate_H

#include <athapascan-1>
#include "linbox/integer.h"
#include "linbox/field/modular-double.h"
#include "linbox/matrix/sparse.h"

/*
 * gmp integers
 */

a1::OStream& operator<<( a1::OStream& out, const Integer& x ) {
	std::ostringstream oss;
	oss << x;
	return out << oss.str();
}

a1::IStream& operator>>( a1::IStream& in, Integer& x ) {
    std::string is;
    in >> is;

	std::istringstream iss(is);
	iss >> x;
    return in;
}

/*
 * modular double
 * we first inherit from them to access protected memebers
 */

namespace kaapi {

template <class T> struct Modular;

/*
 * double specialization
 */
template <>
struct Modular<double> : LinBox::Modular<double>
{
    const double& get_modulus() const { return this->modulus; }
    const unsigned long& get_lmodulus() const { return this->lmodulus; }
    double& get_modulus() { return this->modulus; }
    unsigned long& get_lmodulus() { return this->lmodulus; }
};

} //namespace

a1::OStream& operator<<( a1::OStream& out, const LinBox::Modular<double>& m)
{
    const kaapi::Modular<double>* m_tmp = static_cast<const kaapi::Modular<double>*>(&m);
    return out << m_tmp->get_modulus() << m_tmp->get_lmodulus() ;
}

a1::IStream& operator>>( a1::IStream& in, LinBox::Modular<double>& m)
{
    kaapi::Modular<double>* m_tmp = static_cast<kaapi::Modular<double>*>(&m);
    return in >> m_tmp->get_modulus() >> m_tmp->get_lmodulus() ;
}

#endif //__LINBOX_kaapi_communicate_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
