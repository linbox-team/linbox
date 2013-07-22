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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file   benchmarks/BenchmarkFile.inl
 * @ingroup benchmarks
 * @brief
 */

#ifndef __LINBOX_BENCHMARKFILE_INL
#define __LINBOX_BENCHMARKFILE_INL

#include "linbox/util/debug.h"

#include <stdlib.h>
#include <fstream>
#include <map>

namespace LinBox
{

void BenchmarkFile::printCommaVector(std::ostream& out,const std::vector<CSValue*>& vec)
{
	bool first=true;
	for (int i=0;i<vec.size();++i) {
		if (!first) {
			out << ", ";
		} else {
			first=false;
		}
		if (vec[i] == NULL) {
			out << "-";
		} else {
			vec[i]->print(out);
		}
	}
	out << std::endl;
}

void BenchmarkFile::printMetadata(std::ostream& out)
{
	typedef MetadataMap::iterator MapIT;

	for (MapIT it=metadata_.begin();it!=metadata_.end();++it) {
		out << it->first << ", ";
		it->second->print(out);
		out << std::endl;
	}
        
        typedef TypeMap::iterator TypeMapIT;

        if (!(typeMap_.empty())) {
                out << "types";
                for (TypeMapIT it=typeMap_.begin();it!=typeMap_.end();++it) {
                        out << ", (" << it->first << "," << it->second << ")";
                }
                out << std::endl;
        }
	out << "end, metadata" << std::endl << std::endl;
}

void BenchmarkFile::printFieldTitles(std::ostream& out)
{
	std::vector<CSValue*> fieldVec(numFields_);

	typedef FieldPosMap::iterator MapIT;

	for (MapIT it=fields_.begin();it!=fields_.end();++it) {
		fieldVec[it->second]=new CSString(it->first);
	}

	printCommaVector(out,fieldVec);

	for (int i=0;i<fieldVec.size();++i) {delete fieldVec[i];}

}

void BenchmarkFile::printContents(std::ostream& out)
{
	for (int i=0;i<allTests_.size();++i) {
		printCommaVector(out,allTests_[i]);
	}
}

void BenchmarkFile::write(std::ostream& out)
{
	printMetadata(out);
	printFieldTitles(out);
	printContents(out);
}

void BenchmarkFile::freeTestLine(TestLine& line)
{
	for (int i=0;i<line.size();++i) {
		delete line[i];
	}
}

BenchmarkFile::~BenchmarkFile() {
	for (int i=0;i<allTests_.size();++i) {
		freeTestLine(allTests_[i]);
	}

	freeTestLine(curTest_);

	typedef MetadataMap::iterator MapIT;
	for (MapIT it=metadata_.begin();it!=metadata_.end();++it) {
		delete it->second;
	}
}

BenchmarkFile::MetadataIterator BenchmarkFile::metadataBegin()
{
        return metadata_.begin();
}

void BenchmarkFile::addMetadata(const std::string& key,const CSValue& val)
{
	metadata_.insert(std::pair<std::string,CSValue*>(key,val.clone()));
}


void BenchmarkFile::setType(const std::string& fieldName, const std::string& type)
{
        typeMap_[fieldName]=type;
}

void BenchmarkFile::addDataField(const std::string& fieldName,const CSValue& val)
{
	typedef FieldPosMap::iterator FieldPosMapIT;

	int fieldPos;

	FieldPosMapIT it = fields_.find(fieldName);
	if (it==fields_.end()) {
		fields_.insert(std::pair<std::string,int>(fieldName,numFields_));
		fieldPos=numFields_;
		++numFields_;
		curTest_.resize(numFields_);
	} else {
		fieldPos=it->second;
	}

	curTest_[fieldPos]=val.clone();
}

void BenchmarkFile::pushBackTest()
{
	allTests_.push_back(curTest_);
	curTest_.clear();
	curTest_.resize(numFields_);
}

}
#endif // __LINBOX_BENCHMARKFILE_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
