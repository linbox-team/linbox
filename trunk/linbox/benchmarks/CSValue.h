/* Copyright (C) 2013 LinBox
 * Written by AJS <stachnik@udel.edu>
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

/*! @file   benchmarks/CSValue.h
 * @ingroup benchmarks
 * @brief
 */

#ifndef __LINBOX_CSVALUE_H
#define __LINBOX_CSVALUE_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

namespace LinBox
{

class CSValue {
public:
	virtual ~CSValue() {}

	virtual void print(std::ostream& out) const =0;

	virtual CSValue* clone() const =0;

	virtual int type() const =0;
};

#define CSV_STRING_TYPE 1
#define CSV_INT_TYPE 2
#define CSV_DOUBLE_TYPE 3
#define CSV_DATE_TYPE 4

std::ostream& operator<< (std::ostream& out, const CSValue& v)
{
	v.print(out);
        return out;
}

class CSString : public CSValue {
public:

	CSString() {}

	CSString(const std::string& e) : elt_(e) {}

	std::string getVal() const {return elt_;}

	void print(std::ostream& out) const {out << elt_;}

	CSValue* clone() const {return new CSString(elt_);}

	int type() const {return CSV_STRING_TYPE;}

protected:
	std::string elt_;
};

class CSInt : public CSValue {
public:

	CSInt() : elt_(0) {}

	CSInt(const int e) : elt_(e) {}

	int getVal() const {return elt_;}

	void print(std::ostream& out) const {out << elt_;}

	CSValue* clone() const {return new CSInt(elt_);}

	int type() const {return CSV_INT_TYPE;}

protected:
	int elt_;
};

class CSDouble : public CSValue {
public:

	CSDouble() : elt_(0.0) {}

	CSDouble(const double e) : elt_(e) {}

	double getVal() const {return elt_;}

	void print(std::ostream& out) const {out << elt_;}

	CSValue* clone() const {return new CSDouble(elt_);}

	int type() const {return CSV_DOUBLE_TYPE;}

protected:
	double elt_;
};

class CSDate : public CSValue {
public:

	CSDate() {}

	CSDate(const struct tm& time) : elt_(time) {}

	struct tm getVal() const {return elt_;}

	void print(std::ostream& out) const {out << asctime(&elt_);}

	CSValue* clone() const {return new CSDate(elt_);}

	int type() const {return CSV_DATE_TYPE;}

protected:
	struct tm elt_;
};

}
#endif // __LINBOX_CSVALUE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
