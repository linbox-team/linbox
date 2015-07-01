/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/compose.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __COMPOSE_H
#define __COMPOSE_H


#include "linbox/util/debug.h"
#include "linbox-config.h"
#include "type-chooser.h"
#include "type-checker.h"
#include "linbox/vector/vector-domain.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"
using LinBox::Reader;
using LinBox::Writer;


#include <iostream>
using std::ostream;


#endif

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /* The following enumeration will allow choices for the type of intermediate 
   * vector used by the composed class during the apply function. The user will be 
   * able to choose to use the preferred input vector of A, the preferred output vector of B,
   * a conversion between the two types of vectors, or the default type of vector std::vector<Element>.
   * The choose is there only if the output vector of B and the input vector of A are different. If they
   * are the same, then that vector type will be used for the intermediate vector.
   */

  enum IntermediateVector { DEFAULT, INPUT, OUTPUT, CONVERSION };


  /** @memo Blackbox of a product: C := AB, i.e. Cx := A(Bx).
   * @doc
   * This is a class that multiplies two matrices by implementing an 
   * apply method that calls the apply methods of both of the consituent 
   * matrices, one after the other.
   *
   * This class, like the Black Box archetype from which it is derived, 
   * is templatized by the vector type to which the matrix is applied.  
   * Both constituent matrices must also use this same vector type.
   * For specification of the blackbox members see \Ref{BlackboxArchetype}.
   * 
   * {\bf Template parameter:} must meet the \Ref{Vector} requirement.
   */
  template <class _Blackbox1, class _Blackbox2 = _Blackbox1, IntermediateVector flag = DEFAULT>
    class Compose;

  template <class _Blackbox1,
           class _Blackbox2,
           IntermediateVector flag>
             class Compose 
             {
               public:

                 typedef _Blackbox1 Blackbox1;
                 typedef _Blackbox2 Blackbox2;

                 typedef typename Blackbox1::Field Field;
                 typedef typename Blackbox1::Element Element;
                 typedef typename Blackbox1::PreferredOutputType PreferredOutputType;
                 typedef typename Blackbox2::PreferredInputType PreferredInputType;
                 typedef TypeChooser<typename Blackbox2::PreferredOutputType, typename Blackbox1::PreferredInputType, std::vector<Element> > VectorType; //Choose the intermediate vector type

                 /** Constructor of C := A*B from blackbox matrices A and B.
                  * Build the product A*B of any two black box matrices of compatible dimensions.
                  * Requires A.coldim() equals B.rowdim().
                  */
                 Compose (const Blackbox1 &A, const Blackbox2 &B)
                   : _A_ptr(&A), _B_ptr(&B) 
                   {
                     // Rich Seagraves - "It seems VectorWrapper somehow
                     // became depricated.  Makes the assumption that 
                     // this vector type supports resize"
                     // VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     _z.resize(_A_ptr->coldim());
                   }

                 /** Constructor of C := (*A_ptr)*(*B_ptr).
                  * This constructor creates a matrix that is a product of two black box
                  * matrices: A*B from pointers to them.
                  */
                 Compose (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
                   : _A_ptr(A_ptr), _B_ptr(B_ptr)
                   {
                     linbox_check (A_ptr != (Blackbox1 *) 0);
                     linbox_check (B_ptr != (Blackbox2 *) 0);
                     linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

                     //			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     _z.resize(_A_ptr->coldim());
                   }

                 /** Copy constructor.
                  * Copies the composed matrix (a small handle).  The underlying two matrices
                  * are not copied.
                  */
                 Compose (const Compose<Blackbox1, Blackbox2>& M) 
                   :_A_ptr ( M._A_ptr), _B_ptr ( M._B_ptr)
                   //{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
                   { _z.resize(_A_ptr->coldim());}
                 /*
#ifdef __LINBOX_XMLENABLED
Compose(Reader &R);
#endif
                  */
                 /// Destructor
                 ~Compose () {}

                 /*- Virtual constructor.
                  * Required because constructors cannot be virtual.
                  * Make a copy of the BlackboxArchetype object.
                  * Required by abstract base class.
                  * @return pointer to new blackbox object
                 // 		 */
                 // 		BlackboxArchetype<_Vector> *clone () const
                 // 			{ return new Compose (*this); }

                 /*- Application of BlackBox matrix.
                  * y= (A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */

                 template <class OutVector, class InVector>
                   inline OutVector& apply (OutVector& y, const InVector& x) const
                   {
                     if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _B_ptr->apply (_z, x);
                       _A_ptr->apply (y, _z);
                     }

                     return y;
                   }

                 /*- Application of BlackBox matrix transpose.
                  * y= transpose(A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */
                 template <class OutVector, class InVector>
                   inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                   {
                     if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _A_ptr->applyTranspose (_z, x);
                       _B_ptr->applyTranspose (y, _z);
                     }

                     return y;
                   }

                 /*- Retrieve row dimensions of BlackBox matrix.
                  * This may be needed for applying preconditioners.
                  * Required by abstract base class.
                  * @return integer number of rows of black box matrix.
                  */
                 size_t rowdim (void) const
                 {
                   if (_A_ptr != 0) 
                     return _A_ptr->rowdim ();
                   else 
                     return 0;
                 }

                 /*- Retrieve column dimensions of BlackBox matrix.
                  * Required by abstract base class.
                  * @return integer number of columns of black box matrix.
                  */
                 size_t coldim(void) const 
                 {
                   if (_B_ptr != 0) 
                     return _B_ptr->coldim ();
                   else 
                     return 0;
                 }

                 const Field& field() const {return _A_ptr->field();}

#ifdef __LINBOX_XMLENABLED
                 ostream &write(ostream &os) const
                 {
                   Writer W;
                   if( toTag(W) ) 
                     W.write(os);

                   return os;
                 }

                 bool toTag(Writer &W) const
                 {
                   string s;

                   W.setTagName("MatrixOver");
                   W.setAttribute("rows", Writer::numToString(s, _A_ptr->rowdim()));
                   W.setAttribute("cols", Writer::numToString(s, _B_ptr->coldim()));
                   W.setAttribute("implDetail", "compose");

                   W.addTagChild();
                   W.setTagName("compose");

                   W.addTagChild();
                   _A_ptr->toTag(W);
                   W.upToParent();

                   W.addTagChild();
                   _B_ptr->toTag(W);
                   W.upToParent();

                   W.upToParent();

                   return true;
                 }

#endif


               protected:

                 // Pointers to A and B matrices
                 const Blackbox1 *_A_ptr;
                 const Blackbox2 *_B_ptr;

                 // local intermediate vector
                 mutable typename VectorType::TYPE _z;
                 template <class _OutVector, class _InVector, bool use_list, bool no_conversion> friend class conv_apply_wrap;   
             }; //Default Option


  template <class _Blackbox1,
           class _Blackbox2>
             class Compose<_Blackbox1, _Blackbox2, INPUT> 
             {
               public:

                 typedef _Blackbox1 Blackbox1;
                 typedef _Blackbox2 Blackbox2;

                 typedef typename Blackbox1::Field Field;
                 typedef typename Blackbox1::Element Element;
                 typedef typename Blackbox1::PreferredOutputType PreferredOutputType;
                 typedef typename Blackbox2::PreferredInputType PreferredInputType;

                 /** Constructor of C := A*B from blackbox matrices A and B.
                  * Build the product A*B of any two black box matrices of compatible dimensions.
                  * Requires A.coldim() equals B.rowdim().
                  */
                 Compose (const Blackbox1 &A, const Blackbox2 &B)
                   : _A_ptr(&A), _B_ptr(&B) 
                   {
                     // Rich Seagraves - "It seems VectorWrapper somehow
                     // became depricated.  Makes the assumption that 
                     // this vector type supports resize"
                     // VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     _z.resize(_A_ptr->coldim());
                   }

                 /** Constructor of C := (*A_ptr)*(*B_ptr).
                  * This constructor creates a matrix that is a product of two black box
                  * matrices: A*B from pointers to them.
                  */
                 Compose (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
                   : _A_ptr(A_ptr), _B_ptr(B_ptr)
                   {
                     linbox_check (A_ptr != (Blackbox1 *) 0);
                     linbox_check (B_ptr != (Blackbox2 *) 0);
                     linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

                     //			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     _z.resize(_A_ptr->coldim());
                   }

                 /** Copy constructor.
                  * Copies the composed matrix (a small handle).  The underlying two matrices
                  * are not copied.
                  */
                 Compose (const Compose<Blackbox1, Blackbox2>& M) 
                   :_A_ptr ( M._A_ptr), _B_ptr ( M._B_ptr)
                   //{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
                   { _z.resize(_A_ptr->coldim());}
                 /*
#ifdef __LINBOX_XMLENABLED
Compose(Reader &R);
#endif
                  */
                 /// Destructor
                 ~Compose () {}

                 /*- Virtual constructor.
                  * Required because constructors cannot be virtual.
                  * Make a copy of the BlackboxArchetype object.
                  * Required by abstract base class.
                  * @return pointer to new blackbox object
                 // 		 */
                 // 		BlackboxArchetype<_Vector> *clone () const
                 // 			{ return new Compose (*this); }

                 /*- Application of BlackBox matrix.
                  * y= (A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */

                 template <class OutVector, class InVector>
                   inline OutVector& apply (OutVector& y, const InVector& x) const
                   {
                     if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _B_ptr->apply (_z, x);
                       _A_ptr->apply (y, _z);
                     }

                     return y;
                   }

                 /*- Application of BlackBox matrix transpose.
                  * y= transpose(A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */
                 template <class OutVector, class InVector>
                   inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                   {
                     if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _A_ptr->applyTranspose (_z, x);
                       _B_ptr->applyTranspose (y, _z);
                     }

                     return y;
                   }

                 /*- Retrieve row dimensions of BlackBox matrix.
                  * This may be needed for applying preconditioners.
                  * Required by abstract base class.
                  * @return integer number of rows of black box matrix.
                  */
                 size_t rowdim (void) const
                 {
                   if (_A_ptr != 0) 
                     return _A_ptr->rowdim ();
                   else 
                     return 0;
                 }

                 /*- Retrieve column dimensions of BlackBox matrix.
                  * Required by abstract base class.
                  * @return integer number of columns of black box matrix.
                  */
                 size_t coldim(void) const 
                 {
                   if (_B_ptr != 0) 
                     return _B_ptr->coldim ();
                   else 
                     return 0;
                 }

                 const Field& field() const {return _A_ptr->field();}

#ifdef __LINBOX_XMLENABLED
                 ostream &write(ostream &os) const
                 {
                   Writer W;
                   if( toTag(W) ) 
                     W.write(os);

                   return os;
                 }

                 bool toTag(Writer &W) const
                 {
                   string s;

                   W.setTagName("MatrixOver");
                   W.setAttribute("rows", Writer::numToString(s, _A_ptr->rowdim()));
                   W.setAttribute("cols", Writer::numToString(s, _B_ptr->coldim()));
                   W.setAttribute("implDetail", "compose");

                   W.addTagChild();
                   W.setTagName("compose");

                   W.addTagChild();
                   _A_ptr->toTag(W);
                   W.upToParent();

                   W.addTagChild();
                   _B_ptr->toTag(W);
                   W.upToParent();

                   W.upToParent();

                   return true;
                 }

#endif


               protected:

                 // Pointers to A and B matrices
                 const Blackbox1 *_A_ptr;
                 const Blackbox2 *_B_ptr;

                 // local intermediate vector
                 mutable typename Blackbox1::PreferredInputType _z;
             };// Input option

  template <class _Blackbox1,
           class _Blackbox2>
             class Compose<_Blackbox1, _Blackbox2, OUTPUT> 
             {
               public:

                 typedef _Blackbox1 Blackbox1;
                 typedef _Blackbox2 Blackbox2;

                 typedef typename Blackbox1::Field Field;
                 typedef typename Blackbox1::Element Element;
                 typedef typename Blackbox1::PreferredOutputType PreferredOutputType;
                 typedef typename Blackbox2::PreferredInputType PreferredInputType;

                 /** Constructor of C := A*B from blackbox matrices A and B.
                  * Build the product A*B of any two black box matrices of compatible dimensions.
                  * Requires A.coldim() equals B.rowdim().
                  */
                 Compose (const Blackbox1 &A, const Blackbox2 &B)
                   : _A_ptr(&A), _B_ptr(&B) 
                   {
                     // Rich Seagraves - "It seems VectorWrapper somehow
                     // became depricated.  Makes the assumption that 
                     // this vector type supports resize"
                     // VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     _z.resize(_A_ptr->coldim());
                   }

                 /** Constructor of C := (*A_ptr)*(*B_ptr).
                  * This constructor creates a matrix that is a product of two black box
                  * matrices: A*B from pointers to them.
                  */
                 Compose (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
                   : _A_ptr(A_ptr), _B_ptr(B_ptr)
                   {
                     linbox_check (A_ptr != (Blackbox1 *) 0);
                     linbox_check (B_ptr != (Blackbox2 *) 0);
                     linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

                     //			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     _z.resize(_A_ptr->coldim());
                   }

                 /** Copy constructor.
                  * Copies the composed matrix (a small handle).  The underlying two matrices
                  * are not copied.
                  */
                 Compose (const Compose<Blackbox1, Blackbox2>& M) 
                   :_A_ptr ( M._A_ptr), _B_ptr ( M._B_ptr)
                   //{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
                   { _z.resize(_A_ptr->coldim());}
                 /*
#ifdef __LINBOX_XMLENABLED
Compose(Reader &R);
#endif
                  */
                 /// Destructor
                 ~Compose () {}

                 /*- Virtual constructor.
                  * Required because constructors cannot be virtual.
                  * Make a copy of the BlackboxArchetype object.
                  * Required by abstract base class.
                  * @return pointer to new blackbox object
                 // 		 */
                 // 		BlackboxArchetype<_Vector> *clone () const
                 // 			{ return new Compose (*this); }

                 /*- Application of BlackBox matrix.
                  * y= (A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */

                 template <class OutVector, class InVector>
                   inline OutVector& apply (OutVector& y, const InVector& x) const
                   {
                     if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _B_ptr->apply (_z, x);
                       _A_ptr->apply (y, _z);
                     }

                     return y;
                   }

                 /*- Application of BlackBox matrix transpose.
                  * y= transpose(A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */
                 template <class OutVector, class InVector>
                   inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                   {
                     if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _A_ptr->applyTranspose (_z, x);
                       _B_ptr->applyTranspose (y, _z);
                     }

                     return y;
                   }

                 /*- Retrieve row dimensions of BlackBox matrix.
                  * This may be needed for applying preconditioners.
                  * Required by abstract base class.
                  * @return integer number of rows of black box matrix.
                  */
                 size_t rowdim (void) const
                 {
                   if (_A_ptr != 0) 
                     return _A_ptr->rowdim ();
                   else 
                     return 0;
                 }

                 /*- Retrieve column dimensions of BlackBox matrix.
                  * Required by abstract base class.
                  * @return integer number of columns of black box matrix.
                  */
                 size_t coldim(void) const 
                 {
                   if (_B_ptr != 0) 
                     return _B_ptr->coldim ();
                   else 
                     return 0;
                 }

                 const Field& field() const {return _A_ptr->field();}

#ifdef __LINBOX_XMLENABLED
                 ostream &write(ostream &os) const
                 {
                   Writer W;
                   if( toTag(W) ) 
                     W.write(os);

                   return os;
                 }

                 bool toTag(Writer &W) const
                 {
                   string s;

                   W.setTagName("MatrixOver");
                   W.setAttribute("rows", Writer::numToString(s, _A_ptr->rowdim()));
                   W.setAttribute("cols", Writer::numToString(s, _B_ptr->coldim()));
                   W.setAttribute("implDetail", "compose");

                   W.addTagChild();
                   W.setTagName("compose");

                   W.addTagChild();
                   _A_ptr->toTag(W);
                   W.upToParent();

                   W.addTagChild();
                   _B_ptr->toTag(W);
                   W.upToParent();

                   W.upToParent();

                   return true;
                 }

#endif


               protected:

                 // Pointers to A and B matrices
                 const Blackbox1 *_A_ptr;
                 const Blackbox2 *_B_ptr;

                 // local intermediate vector
                 mutable typename Blackbox2::PreferredOutputType _z;
             }; // Output option 


  template <class _Blackbox1,
           class _Blackbox2>
             class Compose<_Blackbox1, _Blackbox2, CONVERSION> 
             {//Open Conversion Option
               //Nested class to remove dynamic if
               //Class uses meta programming to decide if 
               //a conversion of types is needed.
               //The nested class holds all the private 
               //elements which were held in Compose and 
               //implements all of the member functions.
               //Compose with the conversion option is now a 
               //wrapper of the nested class and calls the 
               //nested class for all of its operations
               template <class _Blackbox1_nest,//Left Blackbox  
               class _Blackbox2_nest,//Right Blackbox 
               class _Inter_Out, //Preferred output vector of Right Blackbox
               class _Inter_In > //Preferred input vector of Left Blackbox
                 class Conversion_nest
                 {   
                   public:
                     typedef _Blackbox1_nest Blackbox1;
                     typedef _Blackbox2_nest Blackbox2;
                     typedef _Inter_Out Inter_Out;
                     typedef _Inter_In Inter_In;
                     /** Constructor of C := A*B from blackbox matrices A and B.
                      * Build the product A*B of any two black box matrices of compatible dimensions.
                      * Requires A.coldim() equals B.rowdim().
                      */
                     Conversion_nest (const Blackbox1 &A, const Blackbox2 &B)
                       : _A_ptr(&A), _B_ptr(&B), VD(_A_ptr->field()) 
                       {
                         // Rich Seagraves - "It seems VectorWrapper somehow
                         // became depricated.  Makes the assumption that 
                         // this vector type supports resize"
                         // VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                         _z.resize(_A_ptr->coldim());
                         _w.resize(_A_ptr->coldim());
                       }

                     Conversion_nest (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
                       : _A_ptr(A_ptr), _B_ptr(B_ptr), VD(_A_ptr->field())
                       {
                         linbox_check (A_ptr != (Blackbox1 *) 0);
                         linbox_check (B_ptr != (Blackbox2 *) 0);
                         linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

                         //			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                         _z.resize(_A_ptr->coldim());
                         _w.resize(_A_ptr->coldim());
                       }
                     /** Copy constructor.
                      * Copies the composed matrix (a small handle).  The underlying two matrices
                      * are not copied.
                      */
                     Conversion_nest (const Compose<Blackbox1, Blackbox2>& M) 
                       :_A_ptr ( M._A_ptr), _B_ptr ( M._B_ptr), VD(_A_ptr->field())
                       //{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
                       { _z.resize(_A_ptr->coldim());
                         -w.resize(_A_ptr->coldim());}
                         /// Destructor
                         ~Conversion_nest () {}

                         template<class OutVector, class InVector>
                           inline OutVector& apply(OutVector& y, const InVector& x) const
                           {
                             if ((_A_ptr != 0) && (_B_ptr != 0)) {
                               _B_ptr->apply (_z, x);
                               VD.copy(_w, _z);
                               _A_ptr->apply (y, _w);

                             }
                             std::cout << "In apply, different interm. types" << std::endl;
                             return y;
                           }
                         template<class OutVector, class InVector>
                           inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                           {
                             if ((_A_ptr != 0) && (_B_ptr != 0)) {
                               _A_ptr->applyTranspose (_z, x);
                               VD.copy(_w, _z);
                               _B_ptr->applyTranspose (y, _w);
                             }

                             return y;
                           }
                         /*- Retrieve row dimensions of BlackBox matrix.
                          * This may be needed for applying preconditioners.
                          * Required by abstract base class.
                          * @return integer number of rows of black box matrix.
                          */
                         size_t rowdim (void) const
                         {
                           if (nested._A_ptr != 0) 
                             return nested._A_ptr->rowdim ();
                           else 
                             return 0;
                         }

                         /*- Retrieve column dimensions of BlackBox matrix.
                          * Required by abstract base class.
                          * @return integer number of columns of black box matrix.
                          */
                         size_t coldim(void) const 
                         {
                           if (nested._B_ptr != 0) 
                             return nested._B_ptr->coldim ();
                           else 
                             return 0;
                         }

                         const typename Blackbox1::Field& field() const {return _A_ptr->field();}
                   protected:
                         // Pointers to A and B matrices
                         const Blackbox1 *_A_ptr;
                         const Blackbox2 *_B_ptr;
                         // Vector Domain to copy output vector into input vector type
                         VectorDomain<typename Blackbox1::Field> VD;
                         // local intermediate vector
                         mutable  Inter_Out _z;
                         mutable  Inter_In _w;
                 };//Perform conversion


               template <class _Blackbox1_nest, //Left Blackbox
                        class _Blackbox2_nest, //Right blackbox
                        class _Inter_Out>      //preferred output of right blackbox = preferred input of left blackbox
                          class Conversion_nest<_Blackbox1_nest,_Blackbox2_nest,_Inter_Out,_Inter_Out>
                          {   
                            public:
                              typedef _Blackbox1_nest Blackbox1;
                              typedef _Blackbox2_nest Blackbox2;
                              typedef _Inter_Out Intermediate;
                              /** Constructor of C := A*B from blackbox matrices A and B.
                               * Build the product A*B of any two black box matrices of compatible dimensions.
                               * Requires A.coldim() equals B.rowdim().
                               */
                              Conversion_nest (const Blackbox1 &A, const Blackbox2 &B)
                                : _A_ptr(&A), _B_ptr(&B) 
                                {
                                  // Rich Seagraves - "It seems VectorWrapper somehow
                                  // became depricated.  Makes the assumption that 
                                  // this vector type supports resize"
                                  // VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                                  _z.resize(_A_ptr->coldim());
                                }

                              Conversion_nest (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
                                : _A_ptr(A_ptr), _B_ptr(B_ptr)
                                {
                                  linbox_check (A_ptr != (Blackbox1 *) 0);
                                  linbox_check (B_ptr != (Blackbox2 *) 0);
                                  linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

                                  //			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                                  _z.resize(_A_ptr->coldim());
                                }
                              /** Copy constructor.
                               * Copies the composed matrix (a small handle).  The underlying two matrices
                               * are not copied.
                               */
                              Conversion_nest (const Compose<Blackbox1, Blackbox2>& M) 
                                :_A_ptr ( M._A_ptr), _B_ptr ( M._B_ptr)
                                //{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
                                { _z.resize(_A_ptr->coldim());
                                }
                              /// Destructor
                              ~Conversion_nest () {}

                              template<class OutVector, class InVector>
                                inline OutVector& apply(OutVector& y, const InVector& x) const
                                {
                                  if ((_A_ptr != 0) && (_B_ptr != 0)) {
                                    _B_ptr->apply (_z, x);
                                    _A_ptr->apply (y, _z);
                                  }

                                  std::cout << "Types are not converted" << std::endl;
                                  return y;
                                }
                              template<class OutVector, class InVector>
                                inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                                {
                                  if ((_A_ptr != 0) && (_B_ptr != 0)) {
                                    _A_ptr->applyTranspose (_z, x);
                                    _B_ptr->applyTranspose (y, _z);
                                  }

                                  return y;
                                }
                              /*- Retrieve row dimensions of BlackBox matrix.
                               * This may be needed for applying preconditioners.
                               * Required by abstract base class.
                               * @return integer number of rows of black box matrix.
                               */
                              size_t rowdim (void) const
                              {
                                if (nested._A_ptr != 0) 
                                  return nested._A_ptr->rowdim ();
                                else 
                                  return 0;
                              }

                              /*- Retrieve column dimensions of BlackBox matrix.
                               * Required by abstract base class.
                               * @return integer number of columns of black box matrix.
                               */
                              size_t coldim(void) const 
                              {
                                if (nested._B_ptr != 0) 
                                  return nested._B_ptr->coldim ();
                                else 
                                  return 0;
                              }

                              const typename Blackbox1::Field& field() const {return _A_ptr->field();}
                            protected:
                              // Pointers to A and B matrices
                              const Blackbox1 *_A_ptr;
                              const Blackbox2 *_B_ptr;
                              // local intermediate vector
                              mutable  Intermediate _z;
                          };// No Conversion
               public:

                 typedef _Blackbox1 Blackbox1;
                 typedef _Blackbox2 Blackbox2;

                 typedef typename Blackbox1::Field Field;
                 typedef typename Blackbox1::Element Element;
                 typedef typename Blackbox1::PreferredOutputType PreferredOutputType;
                 typedef typename Blackbox2::PreferredInputType PreferredInputType;

                 /** Constructor of C := A*B from blackbox matrices A and B.
                  * Build the product A*B of any two black box matrices of compatible dimensions.
                  * Requires A.coldim() equals B.rowdim().
                  */
                 Compose (const Blackbox1 &A, const Blackbox2 &B)
                   : nested(A,B)
                   {
                     // Rich Seagraves - "It seems VectorWrapper somehow
                     // became depricated.  Makes the assumption that 
                     // this vector type supports resize"
                     // VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     //_z.resize(_A_ptr->coldim());
                     //_w.resize(_A_ptr->coldim());
                   }

                 /** Constructor of C := (*A_ptr)*(*B_ptr).
                  * This constructor creates a matrix that is a product of two black box
                  * matrices: A*B from pointers to them.
                  */
                 Compose (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
                   : nested(A_ptr,B_ptr)
                   {
                     //linbox_check (A_ptr != (Blackbox1 *) 0);
                     //linbox_check (B_ptr != (Blackbox2 *) 0);
                     //linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

                     //			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
                     //_z.resize(_A_ptr->coldim());
                     //_w.resize(_A_ptr->coldim());
                   }

                 /** Copy constructor.
                  * Copies the composed matrix (a small handle).  The underlying two matrices
                  * are not copied.
                  */
                 Compose (const Compose<Blackbox1, Blackbox2>& M) 
                   : nested(M.nested)
                   //{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
                   { //_z.resize(_A_ptr->coldim());
                     //-w.resize(_A_ptr->coldim());
                   }
                 /*
#ifdef __LINBOX_XMLENABLED
Compose(Reader &R);
#endif
                  */
                 /// Destructor
                 ~Compose () {}

                 /*- Virtual constructor.
                  * Required because constructors cannot be virtual.
                  * Make a copy of the BlackboxArchetype object.
                  * Required by abstract base class.
                  * @return pointer to new blackbox object
                 // 		 */
                 // 		BlackboxArchetype<_Vector> *clone () const
                 // 			{ return new Compose (*this); }

                 /*- Application of BlackBox matrix.
                  * y= (A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */

                 template <class OutVector, class InVector>
                   inline OutVector& apply (OutVector& y, const InVector& x) const
                   {
                     /*if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _B_ptr->apply (_z, x);
                       if(VectorCheck::same)
                       _A_ptr->apply (y, _z);
                       else{
                       VD.copy(_w, _z);
                       _A_ptr->apply (y, _w);
                       }
                       }
                       return y;*/
                     return nested.apply(y,x);
                   }

                 /*- Application of BlackBox matrix transpose.
                  * y= transpose(A*B)*x.
                  * Requires one vector conforming to the \Ref{LinBox}
                  * vector {@link Archetypes archetype}.
                  * Required by abstract base class.
                  * @return reference to vector y containing output.
                  * @param  x constant reference to vector to contain input
                  */
                 template <class OutVector, class InVector>
                   inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                   {
                     /*if ((_A_ptr != 0) && (_B_ptr != 0)) {
                       _A_ptr->applyTranspose (_z, x);
                       if(VectorCheck::same)
                       _B_ptr->applyTranspose (y, _z);
                       else{
                       VD.copy(_w, _z);
                       _B_ptr->applyTranspose (y, _w);
                       }
                       }

                       return y;*/
                     return nested.applyTranspose(y,x);
                   }

                 /*- Retrieve row dimensions of BlackBox matrix.
                  * This may be needed for applying preconditioners.
                  * Required by abstract base class.
                  * @return integer number of rows of black box matrix.
                  */
                 size_t rowdim (void) const
                 {
                   /*if (nested._A_ptr != 0) 
                     return nested._A_ptr->rowdim ();
                   else 
                     return 0;*/
                     nested.rowdim();
                 }

                 /*- Retrieve column dimensions of BlackBox matrix.
                  * Required by abstract base class.
                  * @return integer number of columns of black box matrix.
                  */
                 size_t coldim(void) const 
                 {
                   /*if (nested._B_ptr != 0) 
                     return nested._B_ptr->coldim ();
                   else 
                     return 0;*/
                     nested.coldim();
                 }

                 const Field& field() const {return nested.field();}

/*#ifdef __LINBOX_XMLENABLED
                 ostream &write(ostream &os) const
                 {
                   Writer W;
                   if( toTag(W) ) 
                     W.write(os);

                   return os;
                 }

                 bool toTag(Writer &W) const
                 {
                   string s;

                   W.setTagName("MatrixOver");
                   W.setAttribute("rows", Writer::numToString(s, nested._A_ptr->rowdim()));
                   W.setAttribute("cols", Writer::numToString(s, nested._B_ptr->coldim()));
                   W.setAttribute("implDetail", "compose");

                   W.addTagChild();
                   W.setTagName("compose");

                   W.addTagChild();
                   nested._A_ptr->toTag(W);
                   W.upToParent();

                   W.addTagChild();
                   nested._B_ptr->toTag(W);
                   W.upToParent();

                   W.upToParent();

                   return true;
                 }

#endif*/


               protected:

                 typedef Conversion_nest<Blackbox1, Blackbox2,typename Blackbox2::PreferredOutputType, typename Blackbox1::PreferredInputType> nested_type;
                 mutable nested_type nested;
             };//Conversion option

  // specialization for _Blackbox1 = _Blackbox2	
  template <class _Blackbox>
    class Compose <_Blackbox, _Blackbox, DEFAULT>
    {
      public:
        typedef _Blackbox Blackbox;

        typedef typename _Blackbox::Field Field;
        typedef typename _Blackbox::Element Element;
        typedef typename _Blackbox::PreferredOutputType PreferredOutputType;
        typedef typename _Blackbox::PreferredInputType PreferredInputType;
        typedef TypeChooser<typename _Blackbox::PreferredOutputType, typename _Blackbox::PreferredInputType, 
                std::vector<typename _Blackbox::Element> > VectorType;


        Compose (const Blackbox& A, const Blackbox& B) {
          _BlackboxL.push_back(&A);
          _BlackboxL.push_back(&B);

          _zl.resize(1);

          _zl.front().resize (A.coldim());
        }

        Compose (const Blackbox* Ap, const Blackbox* Bp) {
          _BlackboxL.push_back(Ap);
          _BlackboxL.push_back(Bp);

          _zl.resize(1);

          _zl.front().resize (Ap ->coldim());
        }

        /** Constructor of C := A*B from blackbox matrices A and B.
         * Build the product A*B of any two black box matrices of compatible dimensions.
         * Requires A.coldim() equals B.rowdim().
         */
        template<class BPVector>
          Compose (const BPVector& v)
          :  _BlackboxL(v.begin(), v.end())
          {

            linbox_check(v.size() > 0);
            _zl.resize(v.size() - 1);
            typename std::vector<const Blackbox*>::iterator b_p;
            typename std::vector<typename VectorType::TYPE>::iterator z_p;
            for ( b_p = _BlackboxL.begin(), z_p = _zl.begin();
                z_p != _zl.end(); ++ b_p, ++ z_p) 
              z_p -> resize((*b_p) -> coldim());
          }

        ~Compose () {}

        template <class OutVector, class InVector>
          inline OutVector& apply (OutVector& y, const InVector& x) const
          {	

            typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
            typename std::vector<typename VectorType::TYPE>::reverse_iterator z_p, pz_p;
            b_p = _BlackboxL.rbegin();
            pz_p = z_p = _zl.rbegin();			

            (*b_p) -> apply(*pz_p, x);
            ++ b_p;  ++ z_p;

            for (; z_p != _zl.rend(); ++ b_p, ++ z_p, ++ pz_p)
              (*b_p) -> apply (*z_p,*pz_p);

            (*b_p) -> apply(y, *pz_p);

            return y;
          }

        /*- Application of BlackBox matrix transpose.
         * y= transpose(A*B)*x.
         * Requires one vector conforming to the \Ref{LinBox}
         * vector {@link Archetypes archetype}.
         * Required by abstract base class.
         * @return reference to vector y containing output.
         * @param  x constant reference to vector to contain input
         */
        template <class OutVector, class InVector>
          inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
          {
            typename std::vector<const Blackbox*>::iterator b_p;
            typename std::vector<typename VectorType::TYPE>::iterator z_p, nz_p;

            b_p = _BlackboxL.begin();
            z_p = nz_p = _zl.begin();

            (*b_p) -> applyTranspose (*z_p, x);

            ++ b_p; ++ nz_p;

            for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p) 
              (*b_p) -> applyTranspose (*nz_p, *z_p);

            (*b_p) -> applyTranspose (y, *z_p);

            return y;
          }

        /*- Retrieve row dimensions of BlackBox matrix.
         * This may be needed for applying preconditioners.
         * Required by abstract base class.
         * @return integer number of rows of black box matrix.
         */
        size_t rowdim (void) const
        {
          return _BlackboxL.front() -> rowdim();
        }

        /*- Retrieve column dimensions of BlackBox matrix.
         * Required by abstract base class.
         * @return integer number of columns of black box matrix.
         */
        size_t coldim(void) const 
        {
          return _BlackboxL[_BlackboxL.size() - 1] -> coldim();
        }

        const Field& field() const {return _BlackboxL.front() -> field();}


      protected:

        // Pointers to A and B matrices
        std::vector<const Blackbox*> _BlackboxL;

        // local intermediate vector
        mutable std::vector<typename VectorType::TYPE > _zl;
    }; //Default option


  template <class _Blackbox>
    class Compose <_Blackbox, _Blackbox, INPUT>
    {
      public:
        typedef _Blackbox Blackbox;

        typedef typename _Blackbox::Field Field;
        typedef typename _Blackbox::Element Element;
        typedef typename _Blackbox::PreferredOutputType PreferredOutputType;
        typedef typename _Blackbox::PreferredInputType PreferredInputType;


        Compose (const Blackbox& A, const Blackbox& B) {
          _BlackboxL.push_back(&A);
          _BlackboxL.push_back(&B);

          _zl.resize(1);

          _zl.front().resize (A.coldim());
        }

        Compose (const Blackbox* Ap, const Blackbox* Bp) {
          _BlackboxL.push_back(Ap);
          _BlackboxL.push_back(Bp);

          _zl.resize(1);

          _zl.front().resize (Ap ->coldim());
        }

        /** Constructor of C := A*B from blackbox matrices A and B.
         * Build the product A*B of any two black box matrices of compatible dimensions.
         * Requires A.coldim() equals B.rowdim().
         */
        template<class BPVector>
          Compose (const BPVector& v)
          :  _BlackboxL(v.begin(), v.end())
          {

            linbox_check(v.size() > 0);
            _zl.resize(v.size() - 1);
            typename std::vector<const Blackbox*>::iterator b_p;
            typename std::vector<PreferredInputType>::iterator z_p;
            for ( b_p = _BlackboxL.begin(), z_p = _zl.begin();
                z_p != _zl.end(); ++ b_p, ++ z_p) 
              z_p -> resize((*b_p) -> coldim());
          }

        ~Compose () {}

        template <class OutVector, class InVector>
          inline OutVector& apply (OutVector& y, const InVector& x) const
          {	

            typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
            typename std::vector<PreferredInputType>::reverse_iterator z_p, pz_p;
            b_p = _BlackboxL.rbegin();
            pz_p = z_p = _zl.rbegin();			

            (*b_p) -> apply(*pz_p, x);
            ++ b_p;  ++ z_p;

            for (; z_p != _zl.rend(); ++ b_p, ++ z_p, ++ pz_p)
              (*b_p) -> apply (*z_p,*pz_p);

            (*b_p) -> apply(y, *pz_p);

            return y;
          }

        /*- Application of BlackBox matrix transpose.
         * y= transpose(A*B)*x.
         * Requires one vector conforming to the \Ref{LinBox}
         * vector {@link Archetypes archetype}.
         * Required by abstract base class.
         * @return reference to vector y containing output.
         * @param  x constant reference to vector to contain input
         */
        template <class OutVector, class InVector>
          inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
          {
            typename std::vector<const Blackbox*>::iterator b_p;
            typename std::vector<PreferredInputType>::iterator z_p, nz_p;

            b_p = _BlackboxL.begin();
            z_p = nz_p = _zl.begin();

            (*b_p) -> applyTranspose (*z_p, x);

            ++ b_p; ++ nz_p;

            for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p) 
              (*b_p) -> applyTranspose (*nz_p, *z_p);

            (*b_p) -> applyTranspose (y, *z_p);

            return y;
          }

        /*- Retrieve row dimensions of BlackBox matrix.
         * This may be needed for applying preconditioners.
         * Required by abstract base class.
         * @return integer number of rows of black box matrix.
         */
        size_t rowdim (void) const
        {
          return _BlackboxL.front() -> rowdim();
        }

        /*- Retrieve column dimensions of BlackBox matrix.
         * Required by abstract base class.
         * @return integer number of columns of black box matrix.
         */
        size_t coldim(void) const 
        {
          return _BlackboxL[_BlackboxL.size() - 1] -> coldim();
        }

        const Field& field() const {return _BlackboxL.front() -> field();}


      protected:

        // Pointers to A and B matrices
        std::vector<const Blackbox*> _BlackboxL;

        // local intermediate vector
        mutable std::vector<typename Blackbox::PreferredInputType > _zl;
    };//Input option

  template <class _Blackbox>
    class Compose <_Blackbox, _Blackbox, OUTPUT>
    {
      public:
        typedef _Blackbox Blackbox;

        typedef typename _Blackbox::Field Field;
        typedef typename _Blackbox::Element Element;
        typedef typename _Blackbox::PreferredOutputType PreferredOutputType;
        typedef typename _Blackbox::PreferredInputType PreferredInputType;


        Compose (const Blackbox& A, const Blackbox& B) {
          _BlackboxL.push_back(&A);
          _BlackboxL.push_back(&B);

          _zl.resize(1);

          _zl.front().resize (A.coldim());
        }

        Compose (const Blackbox* Ap, const Blackbox* Bp) {
          _BlackboxL.push_back(Ap);
          _BlackboxL.push_back(Bp);

          _zl.resize(1);

          _zl.front().resize (Ap ->coldim());
        }

        /** Constructor of C := A*B from blackbox matrices A and B.
         * Build the product A*B of any two black box matrices of compatible dimensions.
         * Requires A.coldim() equals B.rowdim().
         */
        template<class BPVector>
          Compose (const BPVector& v)
          :  _BlackboxL(v.begin(), v.end())
          {

            linbox_check(v.size() > 0);
            _zl.resize(v.size() - 1);
            typename std::vector<const Blackbox*>::iterator b_p;
            typename std::vector<PreferredOutputType>::iterator z_p;
            for ( b_p = _BlackboxL.begin(), z_p = _zl.begin();
                z_p != _zl.end(); ++ b_p, ++ z_p) 
              z_p -> resize((*b_p) -> coldim());
          }

        ~Compose () {}

        template <class OutVector, class InVector>
          inline OutVector& apply (OutVector& y, const InVector& x) const
          {	

            typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
            typename std::vector<PreferredOutputType>::reverse_iterator z_p, pz_p;
            b_p = _BlackboxL.rbegin();
            pz_p = z_p = _zl.rbegin();			

            (*b_p) -> apply(*pz_p, x);
            ++ b_p;  ++ z_p;

            for (; z_p != _zl.rend(); ++ b_p, ++ z_p, ++ pz_p)
              (*b_p) -> apply (*z_p,*pz_p);

            (*b_p) -> apply(y, *pz_p);

            return y;
          }

        /*- Application of BlackBox matrix transpose.
         * y= transpose(A*B)*x.
         * Requires one vector conforming to the \Ref{LinBox}
         * vector {@link Archetypes archetype}.
         * Required by abstract base class.
         * @return reference to vector y containing output.
         * @param  x constant reference to vector to contain input
         */
        template <class OutVector, class InVector>
          inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
          {
            typename std::vector<const Blackbox*>::iterator b_p;
            typename std::vector<PreferredOutputType>::iterator z_p, nz_p;

            b_p = _BlackboxL.begin();
            z_p = nz_p = _zl.begin();

            (*b_p) -> applyTranspose (*z_p, x);

            ++ b_p; ++ nz_p;

            for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p) 
              (*b_p) -> applyTranspose (*nz_p, *z_p);

            (*b_p) -> applyTranspose (y, *z_p);

            return y;
          }

        /*- Retrieve row dimensions of BlackBox matrix.
         * This may be needed for applying preconditioners.
         * Required by abstract base class.
         * @return integer number of rows of black box matrix.
         */
        size_t rowdim (void) const
        {
          return _BlackboxL.front() -> rowdim();
        }

        /*- Retrieve column dimensions of BlackBox matrix.
         * Required by abstract base class.
         * @return integer number of columns of black box matrix.
         */
        size_t coldim(void) const 
        {
          return _BlackboxL[_BlackboxL.size() - 1] -> coldim();
        }

        const Field& field() const {return _BlackboxL.front() -> field();}


      protected:

        // Pointers to A and B matrices
        std::vector<const Blackbox*> _BlackboxL;

        // local intermediate vector
        mutable std::vector<typename Blackbox::PreferredOutputType > _zl;
    };//Output option

  template <class _Blackbox>
    class Compose <_Blackbox, _Blackbox, CONVERSION>
    {
               //Nested class to remove dynamic if
               //Class uses meta programming to decide if 
               //a conversion of types is needed.
               //The nested class holds all the private 
               //elements which were held in Compose and 
               //implements all of the member functions.
               //Compose with the conversion option is now a 
               //wrapper of the nested class and calls the 
               //nested class for all of its operations
        template <class _Blackbox_nest, class _Inter_Out, class _Inter_In> class Conversion_nest
        {   
          public:
            typedef _Blackbox_nest Blackbox;
            typedef _Inter_Out Inter_Out;
            typedef _Inter_In Inter_In;
            Conversion_nest (const Blackbox& A, const Blackbox& B) : VD(A.field())   {
              _BlackboxL.push_back(&A);
              _BlackboxL.push_back(&B);

              _zl.resize(1);
              _wl.resize(1);

              _zl.front().resize (A.coldim());
              _wl.front().resize (A.coldim());
            }

            Conversion_nest (const Blackbox* Ap, const Blackbox* Bp) : VD(Ap->field()){
              _BlackboxL.push_back(Ap);
              _BlackboxL.push_back(Bp);

              _zl.resize(1);
              _wl.resize(1);

              _zl.front().resize (Ap ->coldim());
              _wl.front().resize (Ap ->coldim());
            }

            /** Constructor of C := A*B from blackbox matrices A and B.
             * Build the product A*B of any two black box matrices of compatible dimensions.
             * Requires A.coldim() equals B.rowdim().
             */
            template<class BPVector>
              Conversion_nest (const BPVector& v)
              :  _BlackboxL(v.begin(), v.end()), VD(*(v.begin())->field())
              {

                linbox_check(v.size() > 0);
                _zl.resize(v.size() - 1);
                _wl.resize(v.size() - 1);
                typename std::vector<const Blackbox*>::iterator b_p;
                typename std::vector<Inter_Out >::iterator z_p;
                typename std::vector<Inter_In >::iterator w_p;
                for ( b_p = _BlackboxL.begin(), z_p = _zl.begin(), w_p = _wl.begin();
                    z_p != _zl.end(); ++ b_p, ++ z_p, ++ w_p){
                  z_p -> resize((*b_p) -> coldim());
                  w_p -> resize((*b_p) -> coldim());
                }
              }

            ~Conversion_nest () {}
            /** Constructor of C := A*B from blackbox matrices A and B.
             * Build the product A*B of any two black box matrices of compatible dimensions.
             * Requires A.coldim() equals B.rowdim().
             */

            template<class OutVector, class InVector>
              inline OutVector& apply(OutVector& y, const InVector& x) const
              {
                typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
                typename std::vector<Inter_Out>::reverse_iterator z_p, nz_p;
                typename std::vector<Inter_In>::reverse_iterator w_p, nw_p;
                b_p = _BlackboxL.rbegin();
                z_p = nz_p = _zl.rbegin();			
                w_p = nw_p = _wl.rbegin();			
                (*b_p) -> apply(*z_p, x);
                VD.copy(*w_p, *z_p);
                ++ b_p;  ++ nz_p; ++ nw_p;

                for (; nz_p != _zl.rend(); ++ b_p, ++ nz_p, ++ z_p, ++ w_p, ++ nw_p){
                  (*b_p) -> apply (*nz_p,*w_p);
                  VD.copy(*nw_p, *nz_p);
                }

                (*b_p) -> apply(y, *w_p);

                std::cout << "Blackboxes equal, conversion is performed" << std::endl;
                return y;
              }
            template<class OutVector, class InVector>
              inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
              {
                typename std::vector<const Blackbox*>::const_iterator b_p;
                typename std::vector<Inter_Out>::iterator z_p, nz_p;
                typename std::vector<Inter_In>::iterator w_p, nw_p;

                b_p = _BlackboxL.begin();
                z_p = nz_p = _zl.begin();
                w_p = nw_p = _wl.begin();


                (*b_p) -> applyTranspose (*z_p, x);
                VD.copy(*w_p, *z_p);

                ++ b_p; ++ nz_p, nw_p;

                for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p, ++ w_p, ++ nw_p){ 
                  (*b_p) -> applyTranspose (*nz_p, *w_p);
                  VD.copy(*nw_p, *nz_p);
                }

                (*b_p) -> applyTranspose (y, *w_p);

                return y;
              }
            /*- Retrieve row dimensions of BlackBox matrix.
             * This may be needed for applying preconditioners.
             * Required by abstract base class.
             * @return integer number of rows of black box matrix.
             */
            size_t rowdim (void) const
            {
              return _BlackboxL.front() -> rowdim();
            }

            /*- Retrieve column dimensions of BlackBox matrix.
             * Required by abstract base class.
             * @return integer number of columns of black box matrix.
             */
            size_t coldim(void) const 
            {
              return _BlackboxL[_BlackboxL.size() - 1] -> coldim();
            }

            const typename Blackbox::Field& field() const {return _BlackboxL.front() -> field();}
          protected:
            // Pointers to A and B matrices
            std::vector<const Blackbox*> _BlackboxL;

            //Vector Domain to copy vectors of one type into another
            VectorDomain<typename Blackbox::Field> VD;
            // local intermediate vector
            mutable std::vector<Inter_Out  > _zl;
            mutable std::vector<Inter_In > _wl;
        };//Conversion performed


        template <class _Blackbox_nest, class _Inter_Out> 
          class Conversion_nest<_Blackbox_nest,_Inter_Out,_Inter_Out>
          {   
            public:
              typedef _Blackbox_nest Blackbox;
              typedef _Inter_Out Intermediate;
              Conversion_nest (const Blackbox& A, const Blackbox& B) {
                _BlackboxL.push_back(&A);
                _BlackboxL.push_back(&B);

                _zl.resize(1);

                _zl.front().resize (A.coldim());
                _wl.front().resize (A.coldim());
              }

              Conversion_nest (const Blackbox* Ap, const Blackbox* Bp) {
                _BlackboxL.push_back(Ap);
                _BlackboxL.push_back(Bp);

                _zl.resize(1);

                _zl.front().resize (Ap ->coldim());
              }

              /** Constructor of C := A*B from blackbox matrices A and B.
               * Build the product A*B of any two black box matrices of compatible dimensions.
               * Requires A.coldim() equals B.rowdim().
               */
              template<class BPVector>
                Conversion_nest (const BPVector& v)
                :  _BlackboxL(v.begin(), v.end()) 
                {

                  linbox_check(v.size() > 0);
                  _zl.resize(v.size() - 1);
                  typename std::vector<const Blackbox*>::iterator b_p;
                  typename std::vector<Intermediate >::iterator z_p;
                  for ( b_p = _BlackboxL.begin(), z_p = _zl.begin();
                      z_p != _zl.end(); ++ b_p, ++ z_p){
                    z_p -> resize((*b_p) -> coldim());
                  }
                }

              ~Conversion_nest () {}

              template<class OutVector, class InVector>
                inline OutVector& apply(OutVector& y, const InVector& x) const
                {
                  typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
                  typename std::vector<Intermediate>::reverse_iterator z_p, nz_p;
                  b_p = _BlackboxL.rbegin();
                  z_p = nz_p = _zl.rbegin();			
                  (*b_p) -> apply(*z_p, x);
                  ++ b_p;  ++ nz_p; 

                  for (; nz_p != _zl.rend(); ++ b_p, ++ nz_p, ++ z_p){
                    (*b_p) -> apply (*nz_p,*z_p);
                  }

                  (*b_p) -> apply(y, *z_p);

                  std::cout << "Blackboxes equal, intermediate types equal" << std::endl;
                  return y;
                }
              template<class OutVector, class InVector>
                inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
                {
                  typename std::vector<const Blackbox*>::const_iterator b_p;
                  typename std::vector<Intermediate>::iterator z_p, nz_p;

                  b_p = _BlackboxL.begin();
                  z_p = nz_p = _zl.begin();


                  (*b_p) -> applyTranspose (*z_p, x);

                  ++ b_p; ++ nz_p;

                  for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p){ 
                    (*b_p) -> applyTranspose (*nz_p, *z_p);
                  }

                  (*b_p) -> applyTranspose (y, *z_p);

                  return y;
                }
              /*- Retrieve row dimensions of BlackBox matrix.
               * This may be needed for applying preconditioners.
               * Required by abstract base class.
               * @return integer number of rows of black box matrix.
               */
              size_t rowdim (void) const
              {
                return _BlackboxL.front() -> rowdim();
              }

              /*- Retrieve column dimensions of BlackBox matrix.
               * Required by abstract base class.
               * @return integer number of columns of black box matrix.
               */
              size_t coldim(void) const 
              {
                return _BlackboxL[_BlackboxL.size() - 1] -> coldim();
              }

              const typename Blackbox::Field& field() const {return _BlackboxL.front() -> field();}
            protected:
              // Pointers to A and B matrices
              std::vector<const Blackbox*> _BlackboxL;

              // local intermediate vector
              mutable std::vector<Intermediate > _zl;
          };//No conversion needed
      public:
        typedef _Blackbox Blackbox;

        typedef typename _Blackbox::Field Field;
        typedef typename _Blackbox::Element Element;
        typedef typename _Blackbox::PreferredOutputType PreferredOutputType;
        typedef typename _Blackbox::PreferredInputType PreferredInputType;


        Compose (const Blackbox& A, const Blackbox& B) : nested(A,B)   {
          //_BlackboxL.push_back(&A);
          //_BlackboxL.push_back(&B);

          //_zl.resize(1);
          //_wl.resize(1);

          //_zl.front().resize (A.coldim());
          //_wl.front().resize (A.coldim());
        }

        Compose (const Blackbox* Ap, const Blackbox* Bp) : nested(Ap,Bp){
          //_BlackboxL.push_back(Ap);
          //_BlackboxL.push_back(Bp);

          //_zl.resize(1);
          //_wl.resize(1);

          //_zl.front().resize (Ap ->coldim());
          //_wl.front().resize (Ap ->coldim());
        }

        /** Constructor of C := A*B from blackbox matrices A and B.
         * Build the product A*B of any two black box matrices of compatible dimensions.
         * Requires A.coldim() equals B.rowdim().
         */
        template<class BPVector>
          Compose (const BPVector& v)
          :  nested(v)
          {

            //linbox_check(v.size() > 0);
            //_zl.resize(v.size() - 1);
            //_wl.resize(v.size() - 1);
            //typename std::vector<const Blackbox*>::iterator b_p;
            //typename std::vector<PreferredOutputType >::iterator z_p;
            //typename std::vector<PreferredInputType >::iterator w_p;
            //for ( b_p = _BlackboxL.begin(), z_p = _zl.begin(), w_p = _wl.begin();
            //      z_p != _zl.end(); ++ b_p, ++ z_p, ++ w_p){
            //       z_p -> resize((*b_p) -> coldim());
            //       w_p -> resize((*b_p) -> coldim());
            //}
          }

        ~Compose () {}

        template <class OutVector, class InVector>
          inline OutVector& apply (OutVector& y, const InVector& x) const
          {	

            /*typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
              if(VectorCheck::same){
              typename std::vector<PreferredOutputType >::reverse_iterator z_p, pz_p;
              b_p = _BlackboxL.rbegin();
              pz_p = z_p = _zl.rbegin();			

              (*b_p) -> apply(*pz_p, x);
              ++ b_p;  ++ z_p;

              for (; z_p != _zl.rend(); ++ b_p, ++ z_p, ++ pz_p)
              (*b_p) -> apply (*z_p,*pz_p);

              (*b_p) -> apply(y, *pz_p);
              }
              else{
              typename std::vector<PreferredOutputType>::reverse_iterator z_p, pz_p;
              typename std::vector<PreferredInputType>::reverse_iterator w_p, pw_p;
              b_p = _BlackboxL.rbegin();
              pz_p = z_p = _zl.rbegin();			
              pw_p = w_p = _wl.rbegin();			
              (*b_p) -> apply(*pz_p, x);
              VD.copy(*pw_p, *pz_p);
              ++ b_p;  ++ z_p; ++ w_p;

              for (; z_p != _zl.rend(); ++ b_p, ++ z_p, ++ pz_p, ++ w_p, ++ pw_p){
              (*b_p) -> apply (*z_p,*pw_p);
              VD.copy(*w_p, *z_p);
              }

              (*b_p) -> apply(y, *pw_p);
              }


              return y;*/
            return nested.apply(y,x);
          }

        /*- Application of BlackBox matrix transpose.
         * y= transpose(A*B)*x.
         * Requires one vector conforming to the \Ref{LinBox}
         * vector {@link Archetypes archetype}.
         * Required by abstract base class.
         * @return reference to vector y containing output.
         * @param  x constant reference to vector to contain input
         */
        template <class OutVector, class InVector>
          inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
          {
            /*typename std::vector<const Blackbox*>::iterator b_p;
              if(VectorCheck::same){
              typename std::vector<PreferredOutputType>::iterator z_p, nz_p;

              b_p = _BlackboxL.begin();
              z_p = nz_p = _zl.begin();

              (*b_p) -> applyTranspose (*z_p, x);

              ++ b_p; ++ nz_p;

              for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p) 
              (*b_p) -> applyTranspose (*nz_p, *z_p);

              (*b_p) -> applyTranspose (y, *z_p);
              }
              else{
              typename std::vector<PreferredOutputType>::iterator z_p, nz_p;
              typename std::vector<PreferredInputType>::iterator w_p, nw_p;

              b_p = _BlackboxL.begin();
              z_p = nz_p = _zl.begin();
              w_p = nw_p = _wl.begin();


              (*b_p) -> applyTranspose (*w_p, x);
              VD.copy(*z_p, *w_p);

              ++ b_p; ++ nz_p, nw_p;

              for (; nz_p != _zl.end(); ++ z_p, ++ nz_p, ++ b_p, ++ w_p, ++ nw_p){ 
              (*b_p) -> applyTranspose (*nw_p, *z_p);
              VD.copy(*nz_p, *nw_p);
              }

              (*b_p) -> applyTranspose (y, *z_p);
              }

              return y;*/
            return nested.applyTranspose(y,x);
          }

        /*- Retrieve row dimensions of BlackBox matrix.
         * This may be needed for applying preconditioners.
         * Required by abstract base class.
         * @return integer number of rows of black box matrix.
         */
        size_t rowdim (void) const
        {
          return nested.rowdim();
        }

        /*- Retrieve column dimensions of BlackBox matrix.
         * Required by abstract base class.
         * @return integer number of columns of black box matrix.
         */
        size_t coldim(void) const 
        {
          return nested.coldim();
        }

        const Field& field() const {return nested.field();}


      protected:
        typedef Conversion_nest<Blackbox,
                                typename Blackbox::PreferredOutputType,
                                typename Blackbox::PreferredInputType
                               > nested_type;
        mutable nested_type nested;

    };//Conversion Option


} // namespace LinBox
/*
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/reader-blackbox-factory.h"

namespace LinBox {
template<class _Vector, class _Blackbox1, class _Blackbox2>
Compose<_Vector, _Blackbox1, _Blackbox2>::Compose(Reader &R) 
{
size_t m, n;
ReaderBlackBoxFactory<_Vector> RBFact;

if( !R.expectTagName("MatrixOver")) return;
if( !R.expectAttributeNum("rows", m) || !R.expectAttributeNum("cols", n)) return;

if( !R.expectChildTag()) return;
R.traverseChild();
if(!R.expectTagName("compose") || !R.expectChildTag()) return;
R.traverseChild();
RBFact.reset(R);
if(!RBFact.isBlackBox()) {
R.setErrorString("First Compose Child wasn't a BlackBox");
R.setErrorCode(Reader::OTHER);
return;
}
_A_ptr = static_cast<Blackbox1*>(RBFact.makeBlackBox());

R.upToParent();
if(!R.getNextChild()) {
R.setErrorString("Compose expects two matrices to compose, only got one");
R.setErrorCode(Reader::OTHER);
return;
}
if(!R.expectChildTag()) return;
R.traverseChild();

RBFact.reset(R);
if(!RBFact.isBlackBox()) {
R.setErrorString("Second Compose Child wasn't a BlackBox");
R.setErrorCode(Reader::OTHER);
return;
}

_B_ptr = static_cast<Blackbox2*>(RBFact.makeBlackBox());

R.upToParent();
R.getPrevChild();
R.upToParent();

if(m != _A_ptr->rowdim()) {
R.setErrorString("Given Row dimension of composed matrix does not match actual Row dimensions!");
R.setErrorCode(Reader::OTHER);
return;
}
if(n != _B_ptr->coldim()) {
R.setErrorString("Given Column dimension of composed matrix does not match actual Column dimensions!");
return;
}
_z.resize(_A_ptr->coldim());

return;
}

}
#endif
 */




#endif // __COMPOSE_H
