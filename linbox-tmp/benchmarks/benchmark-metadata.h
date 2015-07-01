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
		svector_t name ; //!  key/name of the metadata. for instance <matrix id="matrix1"/>
		std::string _hash ; //! unique id (used to save space)
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
			name = md->getIds();
			_hash = md->getHash();
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
			return name[1] ;
		}


		void setHash(const std::string & myhash)
		{
			_hash = myhash ;
		}

		const svector_t & getIds() const
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


		std::string hasher(const std::string & data)
		{
#ifdef HAVE_CXX11
			std::hash<std::string> my_hasher ;
			size_t the_hash = my_hasher(data);
#else
			std::locale loc;
			const std::collate<char>& coll = std::use_facet<std::collate<char> >(loc); // const member ?
			size_t the_hash = coll.hash(data.data(),data.data()+data.length());

#endif
			return toString(the_hash);
		}

	public:

		MetaData() :
			name(2)
			, _hash("")
			, keys(0)
			, vals(0)
			, metadata(0)
		{
			name[0] = "unknown" ;
			name[1] = name[0] + '_' + randomAlNum(8) ;
		}

		~MetaData()
		{
			clean();
		}

		MetaData(const MetaData * md)
		{
			deep_copy(md);
		}

		MetaData(const MetaData & md)
		{
			deep_copy(&md);
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

		void setIds (const std::string & key)
		{
			name[0] = key ;
			name[1] = key + "_" + randomAlNum(8) ;

			return;
		}

#ifdef __LINBOX_HAVE_TINYXML2
		void writeMetaData(tinyxml2::XMLElement ** data, tinyxml2::XMLDocument & doc) const
		{
			//! @warning only one name allowed. Todo : matrix1, matrix2,...
			using namespace tinyxml2;
			*data = doc.NewElement( getIds()[0].c_str() );
			linbox_check(*data);
			(*data)->SetAttribute("id",getIds()[1].c_str());
#ifndef NDEBUG
			(*data)->SetAttribute("hash",getHash().c_str());
#endif

			for (index_t i = 0 ; i < keys.size() ; ++i  ) {
				(*data)->SetAttribute(keys[i].c_str(),vals[i].c_str());
			}

			for (index_t i = 0 ; i < getMetaDataSize() ; ++i  ) {
				XMLElement * child = NULL ;
				getMetaData(i)->writeMetaData(&child,doc);
				(*data)->InsertEndChild( child );
			}

			return ;
		}
#endif

		// should not be too public
		std::string getLocalString() const
	       	{
			std::string res = name[0] + ',' ;
			for (size_t i = 0 ; i < keys.size() ; ++i)
				res += keys[i] + '=' + vals[i] + ';';
			return res ;
		}

		const std::string & getHash() const
		{
			return _hash ;
		}


	}; // MetaData

} // LinBox

//
// MetaData specialized
//

#ifdef HAVE_CXX11
#include <unordered_map>
#endif

namespace LinBox {

	// class RepresentationMetaData ;
	class MatrixMetaData ;
	// class VectorMetaData ;
	class StorageMetaData ;
	class GeneratorMetaData ;
	class FieldMetaData ;
	// class SolutionMetaData ;
	class AlgorithmMetaData ;
	class EnvrironmentMetaData ;
	class BenchmarkMetaData ;

	//! Field metadata
	class FieldMetaData : public MetaData {
	private :

	public :
		FieldMetaData()
		{
			setIds("field");
			addValue("name");
			addValue("characteristic");

			hash();
		}
		// a general field/ring has no exponent.

		template<class Field>
		FieldMetaData( const Field & F)
		{
#if 0
			F.getMetaData(this); // this would also print representation
			also det(A, some_mehtod(), Meta) would do  a dry run and print in Meta.
#endif
			setIds("field");

			std::ostringstream a ;
			F.write(a);
			addValue("name",a.str());
			addValue("characteristic", F.characteristic());

			hash();
		}

		void hash()
		{
			std::string data = getLocalString() ;
			linbox_check(getMetaDataSize() == 0);

			setHash(hasher(data));
		}

	}; // FieldMetaData

	//! Matrix metadata
	// what if MetaData changes something in MatrixMetaData and does not update hash ?
	class MatrixMetaData : public MetaData {
		void initMetadata()
		{
			addValue("rowdim");
			addValue("coldim");
			addValue("nbnz");
			addValue("name");

			FieldMetaData FMD;
			addMetaData(&FMD);
		}
	public:
		MatrixMetaData()
		{
			setIds("matrix");
			initMetadata() ;

			hash();
		}

		template<class Matrix>
		MatrixMetaData(Matrix & M, const std::string nom = "N/A" )
		{
			setIds("matrix");
			// M.getMetaData(this);
			addValue("rowdim",M.rowdim());
			addValue("coldim",M.coldim());
			addValue("nbnz",M.size());
			addValue("name",nom);

			FieldMetaData FMD(M.field());
			addMetaData(&FMD);

			hash();
		}

		void hash()
	       	{
			std::string data = getLocalString() ;
			if (getMetaDataSize() == 1) { // first the field. (@todo search for fields, then search for random generator.)
				data += getMetaData(0)->getLocalString();
			}
			setHash(hasher(data));
		}

	} ; // MatrixMetaData

	//! Environment metadata;
	class EnvironmentMetaData : public MetaData {
		void initMetadata()
		{
			// Machine
			// compiler
		}
	public :
		EnvironmentMetaData()
		{
		}
	}; // EnvironmentMetaData

	//! Benchmark metadata;
	class BenchmarkMetaData : public MetaData {
		void initMetadata()
		{
			// problem
			// machine
		}
	public :
		BenchmarkMetaData()
		{
		}
	}; // BenchmarkMetaData

	//! Algorithm metadata;
	class AlgorithmMetaData : public MetaData {
		void initMetadata()
		{
			// name
			// method
		}
	public :
		AlgorithmMetaData()
		{
		}
	}; // AlgorithmMetaData

	//! Generator metadata;
	class GeneratorMetaData : public MetaData {
		void initMetadata()
		{
		}
	public :
		GeneratorMetaData()
		{
		}
	}; // GeneratorMetaData

	//! Storage metadata;
	class StorageMetaData : public MetaData {
		void initMetadata()
		{
		}
	public :
		StorageMetaData()
		{
		}
	}; // StorageMetaData

} // LinBox

//
// typedefs
//

namespace LinBox {
	typedef std::vector<MetaData>  mvector_t ;
} // LinBox

//
// Metadata Container
//

namespace LinBox {

	struct MetaDataSeries {
		mvector_t MetaDataVec ; // vector of metadatas
		// svector_t PointsIDs ; // vector of points ids
		smatrix_t MetaDataIDs ; // MetaDataIDs[i] is the list of indexes in PointIDs, corresponding to the points associated with metadatas MetaDataVec[i]. Could use some std::map as well.

	public:
		MetaDataSeries() :
			MetaDataVec(0)
			, MetaDataIDs(0)
		{};
		// bool exists(index_t & j, const MetaData & m);
		// uses hashes to search for metadata.
		void push_back(const std::string & pointID, const MetaData & pointMD)
		{

			// size_t j = MetaDataVec.size();
			// re-hash here ?
			std::string hsh = pointMD.getHash();
			size_t i ;
			bool found = false ;
			for (i = 0 ; i < MetaDataVec.size() ; ++i)
				if (hsh == MetaDataVec[i].getHash()) {
					found = true ;
					break;
				}
			linbox_check((!found) && (i == MetaDataVec.size()));

			if (! found) {
				MetaDataVec.push_back(pointMD);
				svector_t used_by(0) ;
				used_by.push_back(pointID);
				MetaDataIDs.push_back(used_by) ;
			}
			else {
				MetaDataIDs[i].push_back(pointID);
			}
		}
	};

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
