// =========================================================
// (C) The Linbox Group 1999
// Linbox Cassini Bound
// Time-stamp: <07 Feb 02 21:32:09 Jean-Guillaume.Dumas@imag.fr> 
// =========================================================
#ifndef __CASSINI_H__
#define __CASSINI_H__

#ifndef GIVMAX
#define GIVMAX(a,b) ((a)<(b)?(b):(a))
#endif
#ifndef GIVABS
#define GIVABS(a) ((a)<0?-(a):(a))
#endif

#include <vector.h>

template< class TT > class Cassini {
    TT _aat_diag;
    TT _aat_radius;
    TT _aat_radius1;
    unsigned long ni, nj;
public:

    unsigned long n_row() { return ni; }
    unsigned long n_col() { return nj; }
    
    typedef TT Type_t;

        ///-- cstors:
    Cassini(char * mat_file) : _aat_diag(0), _aat_radius(0), _aat_radius1(0) {
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

        vector< TT > diag,w;
        vector<TT>::iterator di;
 
        FILE* FileDes = fopen(File_Name, "r");
        if (FileDes != NULL) {
            char * tmp = new char[80];
            fscanf(FileDes,"%ld %ld %s\n",&ni, &nj, &tmp) ;
	    // delete [] tmp;

            w = vector< TT >(nj); 
            diag = vector< TT > (ni);
            di = diag.begin();
            
            long i,j, val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
           
            for(vector< TT >::iterator wi = w.begin();wi!= w.end();++wi) 
                *wi = TT(0);

            for(long ii=0; di!=diag.end(); ++ii, ++di) {
                *di = TT(0);
                while (i == (ii+1)) {
                    *di += (val * val);          // diag = \sum <C_i,C_i>
                    w[j-1] += GIVABS( val );     // w = |A^t|.One
                    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
                }
                _aat_diag = GIVMAX( _aat_diag, *di );
            }
           
        }
        fclose(FileDes);
        FileDes = fopen(File_Name, "r");
        if (FileDes != NULL) {
            fscanf(FileDes,"%ld %ld M\n",&ni, &nj) ;
            long i,j, val;
            fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
            di = diag.begin();
            for(long ii=0; di!=diag.end(); ++ii, ++di) {
                TT local_radius = TT(0);
                while (i == (ii+1)) {
                    local_radius += GIVABS( val )* w[j-1];
                    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
                }
                local_radius -= *di;
                if ( local_radius > _aat_radius1) {
                    if ( local_radius > _aat_radius) {
                        _aat_radius1 = _aat_radius;
                        _aat_radius = local_radius;
                    } else
                        _aat_radius1 = local_radius;
                }
            }
            
        }
        fclose(FileDes);
        if (is_gzipped) {
            system(UT);
        }

    }

    Type_t& operator() (Type_t& r) {    
        return r = _aat_diag + (TT)sqrt( _aat_radius * _aat_radius1 );
    };   

};



#endif
