// ======================================================================= //
// Linbox project 1999
// Domain for random black-boxes construction, conditioning and statistics
// Author : Gilles.Villard@imag.fr
// Time-stamp: <08 Mar 00 14:05:28 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //


#ifndef _LIN_DOM_RAND_B_B_C_
#define _LIN_DOM_RAND_B_B_C_

#include <stdio.h>
#include <fstream.h>
#include <time.h>
#include <timer.h>
#include <iostream.h>
 
#include "LinBox/lin_spv_bb.h"


template <class Domain>
class RandBB : public SparseBlackBoxDom < Domain > {    
// In fact no need of the init ? 
public:
    typedef          Domain                          Domain_t;
    typedef typename Domain::Rep                     Domain_Rep;
    typedef typename Domain::Residu_t                Domain_Residu_t;
    typedef          RandBB< Domain >                Self_t;
    typedef typename Rep::PreferredInMatrix_t        PreferredInMatrix_t;
    typedef typename Rep::PreferredOutMatrix_t       PreferredOutMatrix_t;
protected:
    typedef typename Rep::Element_t                  Element_t;
    typedef typename Rep::Row_t                      Row_t;
//protected:
//Domain _domain;
//Rep _container;
public:

 
  //-- Default cstors:
    RandBB() : SparseBlackBoxDom < Domain >() {}
    RandBB(const Domain& D) : SparseBlackBoxDom < Domain >(D) {}

 
  //-- Cstor of recopy: compiler's generated
    RandBB(const Self_t& M) : SparseBlackBoxDom < Domain >(M) {}
 
  //-- Usefull to use the same domain to perform other operations
    Domain getdomain() const { return _domain; }
 
  // ***********************************************************
  // Access to the sparse matrix representation of the black-box 
  // Returns the container  

    Rep& access() {return _container;}


  // ***********************************************************
  // Access to the sparse matrix representation of the black-box 
  // Writes to a file in the sparse format if the entries can be 
  // converted to %ld, otherwise ? 

    void access(char* File_Name){     

        FILE* FileDes = fopen(File_Name, "w");
        if (FileDes != 0) {
           Indice nr,nc;
           nr= _container.n_row();
           nc= _container.n_col();          

           fprintf(FileDes,"%ld %ld M\n",nr,nc);      

           Element_t _entry;

           for (Indice i=0; i<nr; i++) 
	     for (Indice j=0; j< _container[i].size(); j++){            
               _entry=_container[i][j];    
               fprintf(FileDes,"%ld %ld %ld\n",i+1,_entry.getindex() +1,_domain.access(_entry.getvalue())); 
           }
              
           fprintf(FileDes,"%ld %ld %ld\n",0,0,0);           
  
        }
        fclose(FileDes);
    }


    // *********************************************************
    // Random diagonal black-box

    Rep& diag(const Indice nr, const Indice nc){
       
         _container = Rep(nr,nc);
         Indice nmin; 
         if (nr <= nc) nmin=nr; else nmin=nc;
         
         Domain_Rep tmp;

         for (Indice i=0; i<nmin; i++){
             _container[i].push_back(Element_t(i,_domain.GIVRANDOM(tmp)));
         }
    return _container;
    } 


    // **************************************************************************
    // Probability function  
    // nc columns over GF(q), probability for an entry in column j to be nonzero
    //
    // From [Wiedemann 86, p56]
    //   
    // Uses the size of the field.

    static double p_col_Wie1 (const Indice j, const Indice nc, const long unsigned dsize){

        double f1=1-pow(dsize,-1.0);
        double f2=2*(log(nc))*pow(nc+1-j,-1.0);
    
        if (f1 < f2) return(f1); else return(f2); 

    }


    // **************************************************************************
    // Probability function  
    // nc columns, probability for an entry in column j to be nonzero
    //
    // Threshold from [Blomer et al. 97]
    //   
    // Does not use the size of the field. Does not use the column index.

    static double p_col_logn (const Indice j, const Indice nc, const long unsigned dsize){

        return (log(nc))*pow(nc,-1.0);
    
    }


    // **************************************************************************
    // Probability function  
    // nc columns, probability for an entry in column j to be nonzero
    //
    // Constant probability cf [Blomer et al. 97, p418]
    //   
    // Uses the size of the field. Does not use the column index.
 

    static double p_col_cst (const Indice j, const Indice nc, const long unsigned dsize){

        return ((dsize-1)*pow(5*dsize,-1.0));
    
    }


    // *************************************************************************
    // Random (sparse) black-box
    // The entries are chosen at random following a column dependent 
    // distribution law : 
    //                      A[i,j]=0 with probability 1-pcol[j,n,q], i,j=1..
    //                      A[i,j] uniform non zero with probability pcol[j,n,q]
    // where the matrix has n columns and the field has q elements


    Rep& rand_sparse_pnz_col(const Indice nr, const Indice nc, 
                                        double (*pcol) (const Indice, const Indice, const Domain_Residu_t)){

         _container = Rep(nr,nc);

         srand48(BaseTimer::seed());
Indice dbg=0;
         double coin,NZ;
         Domain_Rep tmp;
         
         for (Indice j=0; j<nc; j++){

           NZ=pcol(j+1,nc,_domain.size());
	   for (Indice i=0; i<nr; i++){

              coin=drand48();
              if (coin <= NZ){
                   _container[i].push_back(Element_t(j,_domain.nonzerorandom(tmp)));
              dbg+=1;
              }               
	   }
         }
	 cout << " Old gen., " << dbg << " elements non nuls " << endl;     
    return _container;
    }


    // *************************************************************************
    // Random (sparse) black-box
    // The entries are chosen at random following a column independent 
    // distribution law : 
    //                      A[i,j]=0 with probability 1-pcol[1,n,q], i,j=1..
    //                      A[i,j] uniform non zero with probability pcol[1,n,q]
    // where the matrix has n columns and the field has q elements
  
    Rep& rand_huge_pnz(const Indice nr, const Indice nc, 
                                        double (*pcol) (const Indice, const Indice, const Domain_Residu_t)){

         _container = Rep(nr,nc);

         srand48(BaseTimer::seed());

         double invlogpz;
         Domain_Rep tmp;
         double nb_zero;

         invlogpz=pow(log(1.0-pcol(1,nc,_domain.size())),-1.0);

         Indice i=0; 
         Indice j=0;
         Indice quo_ij;

Indice dbg=0;

         nb_zero = log(drand48())*invlogpz;
         if (nb_zero <=1.0) nb_zero=0.0;
         j = (j + ((Indice) nb_zero))%nc;
         i += (j + ((Indice) nb_zero))/nc;

         while ((i<nr) && (j<nc)){
dbg+=1;
              _container[i].push_back(Element_t(j,_domain.nonzerorandom(tmp)));

              nb_zero = log(drand48())*invlogpz;
              if (nb_zero <=1.0) nb_zero=0.0;
              i += (j + 1+ ((Indice) nb_zero))/nc;
              j = (j + 1+ ((Indice) nb_zero))%nc;
 
         }

	 cout << " New gen., " << dbg << " elements non nuls " << endl;         
    return _container;
    }


};


#endif _LIN_DOM_RAND_B_B_C_
