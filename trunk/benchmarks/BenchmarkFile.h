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

/*! @file   benchmarks/BenchmarkFile.h
 * @ingroup benchmarks
 * @brief
 */

#ifndef __LINBOX_BENCHMARKFILE_H
#define __LINBOX_BENCHMARKFILE_H

#include <stdlib.h>
#include <fstream>
#include <map>
#include <vector>

namespace LinBox
{

class BenchmarkFile {
public:
        typedef std::map<std::string,CSValue*>::iterator MetadataIterator;

	BenchmarkFile() : numFields_(0) {}

	~BenchmarkFile();

	void write(std::ostream& out);

	void addMetadata(const std::string& key,const CSValue& val);

        MetadataIterator metadataBegin();

	void addDataField(const std::string& fieldName,const CSValue& val);

        void setType(const std::string& fieldName, const std::string& type);

	void pushBackTest();

        static CSDate getDateStamp();

        static std::string getDateFormat();

protected:
	typedef std::vector<CSValue*> TestLine;
	typedef std::map<std::string,CSValue*> MetadataMap;
	typedef std::map<std::string,int> FieldPosMap;
        typedef std::map<std::string,std::string> TypeMap;

	void printMetadata(std::ostream& out);

	void printFieldTitles(std::ostream& out);

	void printContents(std::ostream& out);

	void printCommaVector(std::ostream& out,const std::vector<CSValue*>& vec);

	void freeTestLine(TestLine& line);


	MetadataMap metadata_;

        TypeMap typeMap_;

	int numFields_;

	FieldPosMap fields_;

	TestLine curTest_;

	std::vector<TestLine> allTests_;
};

}

#include "BenchmarkFile.inl"

#endif // __LINBOX_BENCHMARKFILE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
