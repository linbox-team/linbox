// =========================================================
// (C) The Linbox Group 1999
// Linbox wrapper for sparse vectors
// file : lin_dom_spv_bb.h
// Time-stamp: <31 Jul 00 21:17:35 Jean-Guillaume.Dumas@imag.fr> 
// =========================================================
#ifndef __SPARSE_ERASE_B_B_DOMAIN_H__
#define __SPARSE_ERASE_B_B_DOMAIN_H__

#include "sparse_vector.h"

#ifndef _SP_BB_VECTOR_
#ifdef _IBB_VECTOR_
#define _SP_BB_VECTOR_ _IBB_VECTOR_
#else
#include <vector.h>
#define _SP_BB_VECTOR_ vector
#endif // _IBB_VECTOR_
#endif // _SP_BB_VECTOR_

template <class Domain>
class SparseBlackBoxEraseDom {
public:
public:
    typedef          Domain                                     Domain_t;
    typedef typename Domain::Rep                                Type_t;

    typedef typename Sparse_Vector<Type_t>::elements            Element_t;

    typedef          _SP_BB_VECTOR_<Sparse_Vector<Type_t> >     element;
    typedef          Sparse_Vector<Type_t>                      Row_t;
    typedef          SparseBlackBoxEraseDom< Domain >           Self_t;
    typedef          _SP_BB_VECTOR_< Type_t >                   PreferredInMatrix_t;
    typedef          _SP_BB_VECTOR_< Type_t >                   PreferredOutMatrix_t;
protected:
    typedef          Sparse_Vector<Type_t>                      SV_t;
    typedef          element                                    Rep;
    typedef          element                                    Storage_t;
    Domain_t _domain;
        /// As a BlackBox is a singleton we can store the only representation
    unsigned long _row_dim, _col_dim, _nz_elem;
    Rep _container;
    
public:
        ///-- Default cstors:
    SparseBlackBoxEraseDom() : _domain(),_nz_elem(0) {};
    SparseBlackBoxEraseDom(const Domain& D) : _domain(D),_nz_elem(0) {}
    SparseBlackBoxEraseDom(const Domain& D, char * mat_file) : _domain(D),_nz_elem(0) { read(mat_file) ; }
    
        ///-- Cstor of recopy: compiler's generated
    SparseBlackBoxEraseDom(const Self_t& M) : _domain(M._domain),_row_dim(M.n_row()),_col_dim(M.n_col()),_nz_elem(0),_container(M._container) {}

        ///-- Usefull to use the same domain to perform other operations
    const Domain_t& getdomain() const { return _domain; }

        ///-- BlackBox size
    long n_row(const Rep& a) const { return a.size(); }
    long n_col(const Rep& a) const { if (a.size() ) return (a[0]).actualsize(); else return 0; }
    long n_elem(const Rep& a) const { 
        long tot=0;
        for(long s=a.size();s--;)
            tot+=a[s].size();
        return tot;
    }

    long size() const { return _row_dim; }
    long n_row() const { return _row_dim; }
    long n_col() const { return _col_dim; }
    long n_elem() const { return _nz_elem; }


        ///-- initializations
    Rep& init(Rep& a, char * mat_file) const {
        unsigned long ni,nj,ne;  
        return read(a,ni,nj,ne,mat_file); 
    }

    unsigned long init_erase(Rep& a, unsigned long& ni, unsigned long& nj, char * mat_file) const {
        unsigned long ne;  
        return read_erase(a,ni,nj,ne,mat_file); 
    }    

    Rep& init(Rep& a) const { return a = _container; }
    Rep& init(char * mat_file) { return read(mat_file); }
    
  // ***********************************************************
  // Access to the sparse matrix representation of the black-box 
  // Reads to a file in the sparse format if the entries can be 
  // converted to %ld, otherwise ?
    Rep& read (Rep& ca, unsigned long& ni, unsigned long& nj, unsigned long& ne, char * mat_file) const { 
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
  
            fscanf(FileDes,"%ld %ld M\n",&ni, &nj) ;
//             ca = Rep( ni ); ne=0;
            ca.resize( ni ); ne=0;

            long i,j, val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
            Type_t cour;
           
            for(long ii=0; ii<ni; ++ii) {
                    // No non-zero element yet
                ca[ii] = SV_t(0,nj);
                while (i == (ii+1)) {
                    _domain.read( cour, val );
                    if (! _domain.iszero( cour )) {
                        ca[ii].push_back( SV_t::value_type(j-1, cour ) );
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

        return ca;
    }
     
  // ***********************************************************
  // Access to the sparse matrix representation of the black-box 
  // Reads to a file in the sparse format if the entries can be 
  // converted to %ld, otherwise ?
    unsigned long read_erase (Rep& ca, unsigned long& ni, unsigned long& nj, unsigned long& ne, char * mat_file) const { 
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
        unsigned long erased = 0, block = 0;
        if (FileDes != NULL) {
  
            fscanf(FileDes,"%ld %ld M\n",&ni, &nj) ;
//             ca = Rep( ni ); ne=0;
            ca.resize( ni ); ne=0;

            long i,j, val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
//             typename Domain_t::element cour;
            Type_t cour;
            vector<long> removecolumns;
           
            for(long ii=0; ii<ni; ++ii) {
                    // No non-zero element yet
//                 ca[ii] = SV_t(0,nj);
                ca[ii].resize(0); ca[ii].reactualsize(nj);
                while (i == (ii+erased+1)) {
//                     _domain.read( cour, val );
                    _domain.assign( cour, val );
                    if (! _domain.iszero( cour )) {
                        if (removecolumns.size()) {
                            vector<long>::const_iterator vvit = removecolumns.begin();
                            for(; vvit != removecolumns.end(); ++vvit) {
                                if (*vvit == j) break;
                                else if (*vvit > j) {
                                    vvit = removecolumns.end();
                                    break;
                                }
                            }
                            if (vvit == removecolumns.end()) {
                                ca[ii].push_back( SV_t::value_type(j-1, cour ) );
                                ++ne;
                            }
                        } else {
                            ca[ii].push_back( SV_t::value_type(j-1, cour ) );
                            ++ne;
                        }
                    }
                    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
                }
                    // Erasing phase
                    // 0 lines are erased
                if (ca[ii].size() == 0) {
                    ++erased;--ii;--ni;
                }
                    // Singleton lines are erased, rank is incremented
                if (ca[ii].size() == 1) {
                    ++block;
                    long rj = ca[ii][0].getindex();
                    if (removecolumns.size()) {
                        vector<long>::iterator rcit = removecolumns.begin();
                        for(; rcit != removecolumns.end(); ++rcit) {
//                             if (*rcit == (rj+1)) {
//                                 --block;
//                                 break;
//                             }
                            if (*rcit > rj) {
                                removecolumns.insert(rcit, rj+1);
                                break;
                            }
                        }
                        if (rcit == removecolumns.end()) removecolumns.push_back(rj+1);
                    } else {
                        removecolumns.push_back(rj+1);
                    }
                    for(long newi=0; newi<ii; ++newi) {
                        for(Row_t::iterator rowit = ca[newi].begin();rowit != ca[newi].end(); ++rowit) {
                            if ((*rowit).getindex() == rj) {
                                ca[newi].erase(rowit);
                                break;
                            } else { 
                                if ((*rowit).getindex() > rj)
                                    break;
                            }
                        }
                    }
                    ++erased;--ii;--ni;
                }
            }
            ca.resize(ni);

            fclose(FileDes);
            if (is_gzipped) {
                system(UT);
            }        
        }

        return block;
    }
     

    void read (unsigned long& ni, unsigned long& nj, char * mat_file) const { 
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
  
            fscanf(FileDes,"%ld %ld M\n",&ni, &nj) ;
        }

        fclose(FileDes);
        if (is_gzipped) {
            system(UT);
        }        
    }
     
    Rep& read (char * mat_file) { return read(_container,_row_dim,_col_dim,_nz_elem, mat_file); }

    unsigned long read_erase (char * mat_file) { return read_erase(_container,_row_dim,_col_dim,_nz_elem, mat_file); }
    
    Rep& read_transpose (Rep& ca, unsigned long& ni, unsigned long& nj, unsigned long& ne, char * mat_file)  {
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
        
            fscanf(FileDes,"%ld %ld M\n",&nj, &ni) ;
//             ca = Rep(ni); ne = 0;
            ca.resize(ni); ne = 0;
            for(long l=0; l<ni; ++l)
                ca[l] = SV_t(0,nj);
            
            long i,j,val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
//             typename Domain_t::element cour;
            Type_t cour;
            
            for(long ii=0; ii<nj; ++ii)
                while (i == (ii+1)) {
//                     _domain.read( cour, val );
                    _domain.assign( cour, val );
                    if (! _domain.iszero( cour )) {
                        ca[j-1].push_back( SV_t::value_type(ii, cour ) );
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
           long nr=n_row(ca),nc=n_col(ca);
           fprintf(FileDes,"%ld %ld M\n",nr,nc);      
           SV_t::value_type _entry;
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

        //-- BlackBox methods

    template<class OutMatrix, class InMatrix>
    OutMatrix& Apply(OutMatrix& res, const InMatrix& vect, const Rep& ca ) const {
        SV_t::value_type toto;
        long k,i;
        res.resize(ca.size());
        for(k=ca.size(); k-- ; )
            for( res[k]=_domain.zero, i=ca[k].size(); i-- ; ) {
                toto = ca[k][i];
                _domain.axpyin(res[k],toto.getvalue(),vect[toto.getindex()]);
            }
        return res;
    }

    template<class OutMatrix, class InMatrix>
    OutMatrix& Apply(OutMatrix& res, const Rep& ca, const InMatrix& vect) const {
//         SV_t::value_type toto;
        Element_t toto;
        long k,i;
        res.resize(ca.size());
        for(k=ca.size(); k-- ; )
            for( res[k]=_domain.zero, i=ca[k].size(); i-- ; ) {
                toto = ca[k][i];
                _domain.axpyin(res[k],toto.getvalue(),vect[toto.getindex()]);
            }
        return res;
    }

    template<class OutMatrix, class InMatrix>
    OutMatrix& Apply(OutMatrix& res, const InMatrix& vect) const {
        return Apply(res, vect, _container);
    }
    
    template<class OutMatrix, class InMatrix>
    OutMatrix& ApplyTrans(OutMatrix& res, const InMatrix& vect, const Rep& ca) const {
        SV_t::value_type toto;
        long k,i;
//         res.resize(n_col());
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
    OutMatrix& ApplyTrans(OutMatrix& res, const Rep& ca, const InMatrix& vect) const {
        SV_t::value_type toto;
        long k,i;
//         res.resize(n_col());
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
    OutMatrix& ApplyTrans(OutMatrix& res, const InMatrix& vect) const {
        return ApplyTrans(res, vect, _container);
    }

        //-- Special Methods

  const Row_t& operator[] (unsigned long i)  const { return _container[i]; }
  Row_t& operator[] (unsigned long i) { return _container[i]; } ;


   template<class Left, class Right>
    Rep& rank_precondition(const Left& l, const Right& r, Rep& ca) const {
        Type_t tmp;
        for(long ii=ca.size()-1; ii>=0; --ii)
            for(long jj=ca[ii].size()-1; jj>=0; --jj)
                ca[ii][jj].change_value( _domain.mulin( _domain.mul(tmp, l[ii], ca[ii][jj].getvalue()), r[ca[ii][jj].getindex()] ) );
        return ca;
    }
            
   template<class Left, class Right>
    Rep& rank_precondition(const Left& l, Rep& ca, const Right& r) const {
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
        long ll=0;
        PreferredOutMatrix_t diag_left(_row_dim);
        for(; ll<_row_dim; ++ll)
            _domain.nonzerorandom(g, diag_left[ll]);
 
        PreferredInMatrix_t diag_right(_col_dim);
        for(ll=0;ll<_col_dim; ++ll)
            _domain.nonzerorandom(g, diag_right[ll]);

        rank_precondition(diag_left, diag_right);
    }


};



#endif // __SPARSE_ERASE_B_B_DOMAIN_H__
