/* Copyright (C) 1999,2005 LinBox
 * Written by  JG Dumas
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

#ifndef __LINBOX_sparse_bb_domain_H
#define __LINBOX_sparse_bb_domain_H

// Linbox wrapper for sparse vectors

#include "linbox/vector/sparse.h"

#ifndef _SP_BB_VECTOR_
#ifdef _IBB_VECTOR_
#define _SP_BB_VECTOR_ _IBB_VECTOR_
#else
#include <vector>
#define _SP_BB_VECTOR_ std::vector
#endif // _IBB_VECTOR_
#endif // _SP_BB_VECTOR_

template <class Domain, class I = unsigned long>
class SparseBlackBoxDom 
{
public:
    typedef          Domain                                     Domain_t;
    typedef          Domain                                     Field;
    typedef typename Domain::Element                            Type_t;
    typedef          _SP_BB_VECTOR_<Sparse_Vector<Type_t, I> >  Element;
    typedef          Sparse_Vector<Type_t, I>                   Row_t;
    typedef          Sparse_Vector<Type_t, I>                   value_type;
    typedef          SparseBlackBoxDom< Domain, I >             Self_t;
    typedef          _SP_BB_VECTOR_< Type_t >                   PreferredInMatrix_t;
    typedef          _SP_BB_VECTOR_< Type_t >                   PreferredOutMatrix_t;
protected:
    typedef          Sparse_Vector<Type_t, I>                   SV_t;
    typedef          Element                                    Rep;
    Domain_t _domain;
        //- As a BlackBox is a singleton we can store the only representation
    I _row_dim, _col_dim, _nz_elem;
    double _lognormdet;
    Rep _container;
    
public:
        //--- Default cstors:
    SparseBlackBoxDom() : _domain(),_nz_elem(0) {};
    SparseBlackBoxDom(const Domain& D) : _domain(D),_nz_elem(0) { }
    SparseBlackBoxDom(const Domain& D, char * mat_file) : _domain(D),_nz_elem(0) { read(mat_file) ; }
    
        //--- Cstor of recopy: compiler's generated
    SparseBlackBoxDom(const Self_t& M) : _domain(M._domain),_row_dim(M.n_row()),_col_dim(M.n_col()),_nz_elem(0),_container(M._container) {}

        //--- Usefull to use the same domain to perform other operations
    const Domain_t& getdomain() const { return _domain; }
    const Domain_t& field() const { return _domain; }

        //--- BlackBox size
    size_t n_row(const Rep& a) const { return a.size(); }
    size_t n_col(const Rep& a) const { if (a.size() ) return a[0].actualsize(); else return 0; }
    size_t n_elem(const Rep& a) const { 
        long tot=0;
        for(long s=a.size();s--;)
            tot+=a[s].size();
        return tot;
    }

    size_t  size() const { return _row_dim; }
    size_t n_row() const { return _row_dim; }
    size_t n_col() const { return _col_dim; }
    size_t rowdim() const { return _row_dim; }
    size_t coldim() const { return _col_dim; }
    size_t n_elem() const { return _nz_elem; }
    double lognorm() const { return _lognormdet; }


    template<typename _Tp1> 
    struct rebind 
    { typedef SparseBlackBoxDom<_Tp1> other; };

        //--- initializations
    Rep& init(Rep& a, char * mat_file) const {
        I ni,nj,ne;  double lognorm;
        return read(a,ni,nj,ne,lognorm,mat_file); 
    }

//    Rep& init(Rep& a) const { return a = _container; }
    void init(Rep& a) { _container = a; _row_dim = n_row(a); _col_dim = n_col(a); _nz_elem = n_elem(a); }
    Rep& init(char * mat_file) { return read(mat_file); }
    
  // ***********************************************************
  // Access to the sparse matrix representation of the black-box 
  // Reads to a file in the sparse format if the entries can be 
  // converted to %ld, otherwise ?

    Rep& read (Rep& ca, I& ni, I& nj, I& ne, double& lognorm, char * mat_file) const { 
        char *UT, *File_Name;
        int is_gzipped = 0;
        size_t s = strlen(mat_file);
        if ((mat_file[--s] == 'z') && (mat_file[--s] == 'g') && (mat_file[--s] == '.')) {
            is_gzipped = 1;
            File_Name = tempnam("/tmp","bbx_");
            UT = new char[s+34+strlen(File_Name)];
            sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
            system(UT);
            sprintf(UT,"\\rm %s", File_Name);
        } else
            File_Name = mat_file;

        FILE* FileDes = fopen(File_Name, "r");
        if (FileDes != NULL) {
 	    char * tmp = new char[200]; unsigned long tni, tnj;
            fscanf(FileDes,"%ld %ld %s\n",&tni, &tnj, &tmp) ;
            ni = tni; nj =tnj;
	    // delete [] tmp;
            ca.resize( ni ); ne=0;

            long i,j;
	    long val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
            typename Domain_t::Element cour;
	
            lognorm = 0.0;
            double normrow = 0.0; 

            for(unsigned long ii=0; ii<ni; ++ii) {
                    // No non-zero element yet
                normrow = 0.0;
                ca[ii].resize(0);
		ca[ii].reactualsize(nj);
                while (i == (ii+1)) {
                    _domain.init( cour, val );
                    if (! _domain.isZero( cour )) {
                        ca[ii].push_back( SV_t::value_type(j-1, cour ) );
                        ++ne;
                        normrow += double(cour) * double(cour);
                        
                    }
                    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
                }
                if (normrow > 1.0) lognorm += log( normrow )/2.0;
            }

            fclose(FileDes);
            if (is_gzipped) {
                system(UT);
            }        
        }

        return ca;
    }
     

    void read (I& ni, I& nj, char * mat_file) const { 
        char *UT, *File_Name;
        int is_gzipped = 0;
        size_t s = strlen(mat_file);
        if ((mat_file[--s] == 'z') && (mat_file[--s] == 'g') && (mat_file[--s] == '.')) {
            is_gzipped = 1;
            File_Name = tempnam("/tmp","bbx_");
            UT = new char[s+44+strlen(File_Name)];
            sprintf(UT,"gunzip -c %s | head -1 > %s", mat_file, File_Name);
            system(UT);
            sprintf(UT,"\\rm %s", File_Name);
        } else
            File_Name = mat_file;

        FILE* FileDes = fopen(File_Name, "r");
        if (FileDes != NULL) {
 	    char * tmp = new char[200]; unsigned long tni, tnj;
            fscanf(FileDes,"%ld %ld %s\n",&tni, &tnj, &tmp) ;
            ni = tni; nj =tnj;
        }

        fclose(FileDes);
        if (is_gzipped) {
            system(UT);
        }        
    }
     
    Rep& read (char * mat_file) { return read(_container,_row_dim,_col_dim,_nz_elem, _lognormdet, mat_file); }
    
    Rep& read_transpose (Rep& ca, I& ni, I& nj, I& ne, char * mat_file)  {
        char *UT, *File_Name;
        int is_gzipped = 0;
        size_t s = strlen(mat_file);
        if ((mat_file[--s] == 'z') && (mat_file[--s] == 'g') && (mat_file[--s] == '.')) {
            is_gzipped = 1;
            File_Name = tempnam("/tmp","bbx_");
            UT = new char[s+34+strlen(File_Name)];
            sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
            system(UT);
            sprintf(UT,"\\rm %s", File_Name);
        } else
            File_Name = mat_file;
        
        FILE* FileDes = fopen(File_Name, "r");
        if (FileDes != NULL) {
 	    char * tmp = new char[200]; unsigned long tni, tnj;
            fscanf(FileDes,"%ld %ld %s\n",&tni, &tnj, &tmp) ;
            ni = tnj; nj =tni;
            ca.resize(ni); ne = 0;
            for(long l=0; l<ni; ++l)
            { ca[l].reactualsize(nj); ca[l].resize(0); }
            
            long i,j;
	    long val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
            typename Domain_t::Element cour;
            
            for(long ii=0; ii<nj; ++ii)
                while (i == (ii+1)) {
                    _domain.init( cour, val );
                    if (! _domain.isZero( cour )) {
                        ca[j-1].push_back( SV_t::value_type(I(ii), cour ) );
                        ++ne;
                    }
                    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
                }     
        }

        fclose(FileDes);
        if (is_gzipped) {
            system(UT);
        }
    }

    Rep& read_transpose (char * mat_file) { return read_transpose(_container, _row_dim, _col_dim, _nz_elem, mat_file); }
    

    void write(char * O_File_Name, const Rep& ca ) const {
        FILE* FileDes = fopen(O_File_Name, "w");
        if (FileDes != 0) {
           I nr=n_row(ca),nc=n_col(ca);
           fprintf(FileDes,"%ld %ld M\n",nr,nc);      
           typename SV_t::value_type _entry;
           for (long i=0; i<nr; i++) 
             for (long j=0; j< ca[i].size(); j++) {            
               _entry=ca[i][j];    
               fprintf(FileDes,"%ld %ld %ld\n",i+1,_entry.getindex() +1,_domain.write( _entry.getvalue() ));
	     }
           fprintf(FileDes,"%ld %ld %ld\n",0,0,0);           
        }
        fclose(FileDes);
    }
    
    void write(char * O_File_Name) const { write(O_File_Name, _container) ; }


  // ***********************************************************
  // Access to the sparse matrix representation of the black-box 
  // Reads from a stream in the sparse format if the entries 
  // read from long 
  // GV - PG

    Rep& read (std::istream& is, Rep& ca, I& ni, I& nj, I& ne) const { 

 	    char * tmp = new char[200]; 

            is >> ni;
            is >> nj;
            is >> tmp;

	    // delete [] tmp;
            ca.resize( ni ); ne=0;
//          ca = Rep( ni ); ne=0;

            long i,j;
	    long val;
            is >> i >> j >> val ; 
            typename Domain_t::Element cour;
	
            for(unsigned long ii=0; ii<ni; ++ii) {
                    // No non-zero element yet
//                ca[ii] = SV_t(0,nj);
                ca[ii].resize(0);
		ca[ii].reactualsize(nj);
                while (i == (ii+1)) {
                    _domain.init( cour, val );
                    if (! _domain.isZero( cour )) {
                        ca[ii].push_back( SV_t::value_type(j-1, cour ) );
                        ++ne;
                    }
                    is >> i >> j >> val ; 
                }
            }

        return ca;
    }
     

     Rep& read (std::istream& is)
           { return read(is, _container, _row_dim, _col_dim, _nz_elem); } 



  // ***********************************************************
  // Write a sparse matrix to an output stream 
  // GV - PG
  // JGD : removed -1

    std::ostream& write(std::ostream& os, const Rep& ca ) const {

           I nr=n_row(ca),nc=n_col(ca);
           os << nr << " " << nc << "\n";   
           typename SV_t::value_type _entry;
           for (long i=0; i<nr; i++) 
             for (long j=0; j< ca[i].size(); j++) {            
               _entry=ca[i][j];    
               os << i+1 << " " << _entry.getindex() +1 << " "; 
               _domain.write(os, _entry.getvalue() );
               os << "\n"; 
	     }
	   return os << "\n";        
    }
    
    std::ostream& write(std::ostream& os) const { return write(os, _container) ; }




        //-- BlackBox methods

    template<class OutMatrix, class InMatrix>
    OutMatrix& apply(OutMatrix& res, const InMatrix& vect, const Rep& ca ) const {
        typename SV_t::value_type toto;
        I k,i;
        res.resize(n_row());
        for(k=n_row(); k-- ; )
            for( res[k]=_domain.zero, i=ca[k].size(); i-- ; ) {
                toto = ca[k][i];
                _domain.axpyin(res[k],toto.getvalue(),vect[toto.getindex()]);
            }
        return res;
    }

    template<class OutMatrix, class InMatrix>
    OutMatrix& apply(OutMatrix& res, const InMatrix& vect) const {
        return apply(res, vect, _container);
    }
    
    template<class OutMatrix, class InMatrix>
    OutMatrix& applyTranspose(OutMatrix& res, const InMatrix& vect, const Rep& ca) const {
        typename SV_t::value_type toto;
        I k,i;
        res.resize(n_col());
        for(i=res.size(); i-- ; )
            res[i] = _domain.zero;
        for(k=ca.size(); k-- ; )
            for(i=ca[k].size(); i-- ; ) {
                toto = ca[k][i];
                _domain.axpyin(res[ toto.getindex() ],toto.getvalue(), vect[k]);
            }
        return res;
    }

    template<class OutMatrix, class InMatrix>
    OutMatrix& applyTranspose(OutMatrix& res, const InMatrix& vect) const {
        return applyTranspose(res, vect, _container);
    }

        //-- Special Methods

  const Row_t& operator[] (const I i)  const { return _container[i]; }
  Row_t& operator[] (const I i) { return _container[i]; } ;


   template<class Left, class Right>
    Rep& rank_precondition(const Left& l, const Right& r, Rep& ca) const {
        Type_t tmp;
        for(long ii=ca.size()-1; ii>=0; --ii)
            for(long jj=ca[ii].size()-1; jj>=0; --jj)
                ca[ii][jj].change_value( _domain.mulin( _domain.mul(tmp, l[ii], ca[ii][jj].getvalue()), r[ca[ii][jj].getindex()] ) );
        return ca;
    }
            
    template<class Left, class Right>
    Rep& rank_precondition(const Left& l, const Right& r) {
        return rank_precondition(l, r, _container);
    }
            
    template<class RandGen>
    void precondition(RandGen& g) {
        I ll=0;
        PreferredOutMatrix_t diag_left(_row_dim);
        for(; ll<_row_dim; ++ll)
            _domain.nonzerorandom(g, diag_left[ll]);
 
        PreferredInMatrix_t diag_right(_col_dim);
        for(ll=0;ll<_col_dim; ++ll)
            _domain.nonzerorandom(g, diag_right[ll]);

        rank_precondition(diag_left, diag_right);
    }


    Type_t& trace(Type_t& t, const Rep& ca) {
        t = _domain.zero;
        for(long ii=ca.size()-1; ii>=0; --ii) 
            for(long jj=ca[ii].size()-1; jj>=0; --jj)
               if (ca[ii][jj].getindex() == (ii))
                    _domain.addin(t, ca[ii][jj].getvalue());
        return t;
    }       

    Type_t& trace(Type_t& t) {
        return trace(t, _container);
    }       

    Type_t& trace_ata(Type_t& t, const Rep& ca) {
        t = _domain.zero;
        for(long ii=ca.size()-1; ii>=0; --ii)
            for(long jj=ca[ii].size()-1; jj>=0; --jj)
                _domain.axpyin(t, ca[ii][jj].getvalue(), ca[ii][jj].getvalue());
        return t;
    }       

    Type_t& trace_ata(Type_t& t) {
        return trace_ata(t, _container);
    }       

};

#endif // __LINBOX_sparse_bb_domain_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
