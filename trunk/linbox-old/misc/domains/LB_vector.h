#include <vector>
#include <string>
#ifndef _SDVECTOR_H_
#define _SDVECTOR_H_

/* 
 * Sparse_iterators are forward iterators on sequences over an entry type with 
 * a designated special value.   sparse_iterators "skip over" instances 
 * of the designated special value. It is assumed that operator==
 * is defined on the entry type.
 * 
 * SD_vectors are stl vectors with the addition of 
 * the sparse_iterator class and the sparse_iterators
 * begin() and end().
 */ 


class SD_dense_vector_tag {};
class SD_sparse_vector_tag {};


template<class Entry, Entry zip, class vector_tag>
class SD_vector{};

template<class entry>
class SD_vector<entry, SD_dense_vector_tag>: public vector<entry>
{
 public:

  SD_vector(long dim = 0):vector<entry>(dim){}

  typedef vector<entry> super;

  class iterator : public super::iterator
  {
    // //SD_vector& mine;
    //iterator start;
    //iterator me;
    iterator done;
    entry zip;

   public:

    sparse_iterator(){} 

    sparse_iterator(const iterator end, entry& z) 
    { done = end; zip = z; }
    }

    sparse_iterator& nextnz()
    { do ++(*this) while(*this != done && **this == zip);
      return *this;
    }
  }; // class sparse_iterator

  iterator s_begin() const
  { return sparse_iterator(*this, begin());}

  sparse_iterator s_end() const
  { return sparse_iterator(*this, end()); }

  /*
  class const_sparse_iterator : sparse_iterator
  {public:
    const Entry& operator*() const
    { return *me; }
  };

  const_sparse_iterator s_begin() const
  { return const_sparse_iterator(begin(), end());}

  const_sparse_iterator s_end() const
  { return const_sparse_iterator(end(), end()); }
  */

}; // SD_vector< dense .. >

/*************** do sparse case later ***********8
template<class Ring>
class SD_vector<SD_sparse_vector_tag, Ring> : public vector<pair<Ring,long> >
{
  typedef Ring value_type;
  typedef Ring* pointer;
  typedef Ring& reference;
  typedef const Ring& const_reference;
  typedef long size_type;
  typedef long difference_type;
  class iterator
  { vector<pair<Ring,long> >::iterator loc; 
    long index;
   public:
    iterator()
    { loc = begin();// of supertype of enclosing SD_vector type ...
      index = 0;
    }

    const Ring& operator*() const // rvalue
    {
      while( loc < end() && index < (*loc).second() ) loc++;
      if ( loc < end() && index == (*loc).second() ) return (*loc).first;
      else return Ring::zero();
    }
    Ring& operator*() const // lvalue
    {
      while( index < (*loc).second() ) loc++;
      if ( index == (*loc).second() ) return (*loc).first;
      else return v.insert(pair(Ring::zero(),index));
    }
    iterator& operator++() { ++index; return *this; }
    iterator& operator++(int) 
    { iterator p =*this;
      index++; 
      return p; 
    }
    iterator& operator+(const long n) { index += n; return *this; }
   private:
    long index;
  }; // class iterator

  class sparse_iterator
  {
  }; // class sparse_iterator

  iterator begin(){...}
  iterator end(){...}
  sparse_iterator s_begin(){...}
  sparse_iterator s_end(){...}
  ...
}; // SD_vector< sparse .. >


/**************
template<class vector_tag, class Ring>
class Matrix{};

template<class Ring>
class Matrix<SD_dense_vector_tag, Ring>
{
 public:
  Matrix(){}
  typedef Ring entry_domain;
 private:
  typedef SD_vector<SDdense_vector_tag, Ring> row;
  SD_vector<SDdense_vector_tag, row> entries;
}; // dense Matrix class

template<class Ring>
class Matrix<SDsparse_vector_tag, Ring>
{
 public:
  Matrix()
  {}

  // iterators??

  row operator[long i] { return entries[i]; }

 private:
  typedef SD_vector<SDsparse_vector_tag, Ring> row;
  SD_vector<SDsparse_vector_tag, row> entries;
}; // sparse Matrix class

**************/

/***********
uses:

SDvector<SDsparse_vector_tag, myfield> v;
SDvector<SDdense_vector_tag, myfield> w;
v[i] = w[i];

*/
#endif
