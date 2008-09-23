#ifndef _KAAPI_COMMUNICATE_H
#define _KAAPI_COMMUNICATE_H

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

#endif
