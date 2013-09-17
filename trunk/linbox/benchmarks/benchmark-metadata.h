/* Copyright (C) 2013 LinBox
 * Written by BB <bbboyer@ncsu.edu>
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

/*! @file   benchmarks/benchmark-metadata.h
 * @ingroup benchmarks
 * @brief   metadata
  */

#ifndef __LINBOX_benchmarks_benchmark_metadata_H
#define __LINBOX_benchmarks_benchmark_metadata_H

#include "benchmark-utils.h"

//
// Metadata
//


namespace LinBox {

	class MetaData ;

	//! This is the general metadata class
	class MetaData {
	private :
		std::string name ;    //! name tag to refer to.
		// std::string hash ; //! unique id (used to save space)
		svector_t  keys ;     //! keys
		svector_t  vals ;     //! values
		std::vector<MetaData *> metadata ; //! child metadata array

	private :

		//! @internal recursively free memory
		void clean ()
		{
			if (this != NULL) {
				for (size_t i = 0 ; i < metadata.size() ; ++i) {
					metadata[i]->clean();
					delete metadata[i] ;
					metadata[i] = NULL ;
				}
			}
			return;
		}

		//!  @internal copy all, including metadata children
		void deep_copy(const MetaData * md)
		{
			name = md->getName();
			keys = md->getKeys();
			vals = md->getVals();
			size_t md_size = md->getMetaDataSize();
			if (md_size) {
				metadata.resize(md_size,NULL);
				for (size_t i = 0 ; i <md_size ; ++i) {
					metadata[i] = new MetaData ;
					metadata[i]->deep_copy(md->getMetaData(i));
				}
			}

			return;
		}

		//!  @internal adds some metadata child
		void push_back(const MetaData * md)
		{
			size_t md_size = getMetaDataSize();
			metadata.resize(md_size+1,NULL);
			metadata[md_size] = new MetaData ;
			metadata[md_size]->deep_copy(md);
		}

	protected :

		const std::string & getName() const
		{
			return name ;
		}

		const svector_t & getKeys() const
		{
			return keys ;
		}

		const svector_t & getVals() const
		{
			return vals ;
		}

		size_t getMetaDataSize() const
		{
			return metadata.size();
		}

		const MetaData * getMetaData(const size_t & i) const
		{
			return metadata[i];
		}

	public:

		MetaData() :
			name(randomAlNum(8))
			, keys(0)
			, vals(0)
			, metadata(0)
		{}

		~MetaData()
		{
			clean();
		}

		MetaData(const MetaData * md)
		{
			deep_copy(md);
		}

		template<class T>
		void changeValue(const std::string & keyword, T& value)
		{
			index_t i ;
		       	bool ok = findKeyword(i, keys.begin(), keys.end(), keyword);
			if ( !ok )
				// throw LinBoxError("undefined keyword",keyword);
				throw LinBoxError("undefined keyword");

			vals[i] = toString(value);

			return;
		}

		const std::string & getValue(const std::string & keyword)
		{
			index_t i ;
		       	bool ok = findKeyword(i, keys.begin(), keys.end(), keyword);
			if ( !ok )
				// throw LinBoxError("undefined keyword",keyword);
				throw LinBoxError("undefined keyword");

			return vals[i] ;
		}

		void addValue(const std::string & nom, const std::string & val = "N/A")
		{
			keys.push_back(nom);
			vals.push_back(val);
			linbox_check(keys.size() == vals.size());
			return;
		}

		template <class T>
		void addValue(const std::string & nom, const T & val )
		{
			keys.push_back(nom);
			vals.push_back(toString(val));
			linbox_check(keys.size() == vals.size());
			return;
		}

		void addMetaData( const MetaData * md)
		{
			push_back(md);
		}

		void setName (const std::string & nom)
		{
			name = nom ;

			return;
		}

#ifdef __LINBOX_HAVE_TINYXML2
		void writeMetaData(tinyxml2::XMLElement * data, tinyxml2::XMLDocument & doc)
		{
			//! @warning only one name allowed. Todo : matrix1, matrix2,...
			using namespace tinyxml2;
			data = doc.NewElement( name.c_str() );

			for (index_t i = 0 ; i < keys.size() ; ++i  ) {
				data->SetAttribute(keys[i].c_str(),vals[i].c_str());
			}

			for (index_t i = 0 ; i < getMetaDataSize() ; ++i  ) {
				XMLElement * child ;
				writeMetaData(child,doc);
				data->InsertEndChild( child );
			}

			return ;
		}
#endif


	}; // MetaData

} // LinBox

//
// MetaData specialized
//

namespace LinBox {

	class RepresentationMetaData ;
	class MatrixMetaData ;
	class VectorMetaData ;
	class StorageMetaData ;
	class GeneratorMetaData ;
	class FieldMetaData ;
	class SolutionMetaData ;
	class AlgorithmMetaData ;
	class MachineMetaData ;
	class BenchmarkMetaData ;

	//! Field metadata
	class FieldMetaData : public MetaData {
	public :
		FieldMetaData()
		{
			setName("field");
			addValue("name");
			addValue("characteristic");
		}
		// a general field/ring has no exponent.

		template<class Field>
		FieldMetaData( const Field & F)
		{
#if 0
			F.getMetaData(this); // this would also print representation
			also det(A, some_mehtod(), Meta) would do  a dry run and print in Meta.
#endif

			std::ostringstream a ;
			F.write(a);
			addValue("name",a.str());
			addValue("characteristic", F.characteristic());

		}
	}; // FieldMetaData

	//! Matrix metadata
	class MatrixMetaData : public MetaData {
		void initMetadata()
		{
			addValue("rowdim");
			addValue("coldim");
			addValue("nbnz");

			FieldMetaData FMD;
			addMetaData(&FMD);
		}
	public:
		MatrixMetaData()
		{
			setName("matrix");
			initMetadata() ;
		}

		template<class Matrix>
		MatrixMetaData(Matrix & M )
		{
			// M.getMetaData(this);
			addValue("rowdim",M.rowdim());
			addValue("coldim",M.coldim());
			addValue("nbnz",M.size());

			FieldMetaData FMD(M.field());
			addMetaData(&FMD);

		}

	} ; // MatrixMetaData

} // LinBox

//
// typedefs
//

namespace LinBox {
	typedef std::vector<MetaData>  mvector_t ;
} // LinBox



#ifdef LinBoxSrcOnly
#include "benchmarks/benchmark-metadata.C"
#endif

#endif // __LINBOX_benchmarks_benchmark_metadata_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
