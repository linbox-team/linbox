/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by:
 *            BB <bbboyer@ncsu.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/** @file linbox/matrix/MatrixDomain/apply-domain.h
 * @brief NO DOC
 */

#ifndef __LINBOX_matrixdomain_apply_domain_H
#define __LINBOX_matrixdomain_apply_domain_H

#include "linbox/integer.h"
#include "linbox/linbox-tags.h"

namespace LinBox {

	struct applyMethod {
		struct Sequential {} ;
		struct OpenMP {} ;
		struct OpenCL {} ;
		struct Cuda {};
	} ;

	template<class _Matrix/*, class Method = applyMethod::Sequential */>
	class applyDomain {
	public:
		typedef typename     _Matrix::Field Field;
		typedef typename Field::Element Element;
	private :
		const Field & _field ;
		MatrixDomain<Field> _MD ;
	private :
		// y = a A x + b y
		template<class _In, class _Out, class Matrix>
		_Out& _apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     , typename ContainerCategories::Vector
			     , typename ContainerCategories::Vector
			    ) const
		{
			if (_field.isZero(a)) {
				_MD.vectorDomain().mulin(y,b);
				return y ;
			}

			if (t == Tag::Transpose::Trans ) {
#if 0
				typename Matrix::ConstColIterator i = A.colBegin ();
				typename _Out::iterator j = y.begin ();
				for (; j != y.end (); ++j, ++i)
					A._VD.dot (*j, x, *i);
#else
				TransposeMatrix<const Matrix> At(A);
				if (_field.isZero(b)) { /* y = a A x */
					_MD.vectorMul(y,At,x); /* y = Ax */

					if (_field.isOne(a)) {
						// return y ;
					}
					else if (_field.isMOne(a) ) { /* y = -y */
						_MD.vectorDomain().negin(y);
					}
					else {/* y = a y */
						_MD.vectorDomain().mulin(y,a);
					}
				}
				else if (_field.isOne(b)) { /* y = a A x + y */
					throw LinBoxError("not implemented yet");
				}
				else if (_field.isMOne(b)) { /* y = a A x - y */
					throw LinBoxError("not implemented yet");
				}
				else { /* y = a A x + b y */
					_MD.vectorDomain().mulin(y,b);
					_In xx(x);
					_MD.vectorDomain().mulin(xx,a);
					_MD.vectorAxpyin(y,At,xx) ;
				}
#endif
			}
			else { /* not Transpose */
				linbox_check( t == Tag::Transpose::NoTrans ) ;
				if (_field.isZero(b)) {/* y = a A x */
					_MD.vectorMul(y,A,x);
					if (_field.isOne(a)) {

					}
					else if (_field.isMOne(a) )
						_MD.vectorDomain().negin(y);
					else  {
						_MD.vectorDomain().mulin(y,a);
					}
				}
				else if (_field.isOne(b)) {/* y = a A x + y */
					throw LinBoxError("not implemented yet");
				}
				else if (_field.isMOne(b)) {/* y = a A x - y */
					throw LinBoxError("not implemented yet");
				}
				else {/* y = a A x + b y */
					_MD.vectorDomain().mulin(y,b);
					_In xx(x);
					_MD.vectorDomain().mulin(xx,a);
					_MD.vectorAxpyin(y,A,xx);
				}

			}
			return y ;
		}

	template<class _In, class _Out, class Matrix>
		_Out& _apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     , typename ContainerCategories::Any
			     , typename ContainerCategories::Any
			     ) const
		{
			// linbox_check(_field.isZero(a) && _field.isOne(a));
			if (t == Tag::Transpose::Trans ) {
				typename Matrix::ConstColIterator i = A.colBegin ();
				typename _Out::iterator j = y.begin ();
				for (; j != y.end (); ++j, ++i)
					A._VD.dot (*j, x, *i);
			}
			else {
				_MD.vectorMul(y,A,x);
			}
			return y;
		}

	public:
		applyDomain(const Field & F) :
			_field(F)
			, _MD(F)
		{}
		/*! Operation \f$y gets a A x + b y\f$ or \f$y gets a A^t x + b y\f$ .
		 * General axpy with x and y vectors or matrices
		 */
		template<class _In, class _Out, class Matrix>
		_Out& apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     // , Method & m
			   ) const
		{
			// std::cout << "called generic" << std::endl;
			return _apply(t,y,a,A,b,x
				      ,typename ContainerTraits<_In> ::ContainerCategory()
				      ,typename ContainerTraits<_Out>::ContainerCategory()
				     );
		}

#if 0
		/*! Operation \f$y gets a A x + b y\f$ or \f$y gets a A^t x + b y\f$ .
		 *(can overwrite A, uses less memory that \c apply)
		 */
		template<class _In, class _Out>
		_Out& applyIn( LINBOX_enum(Tag::Transpose) t
			       , _Out& y
			       , const Element & a
			       , Matrix & A
			       , const Element & b
			       , const _In & x
			       // , Method & m
			     );
#endif

	};

	template<class T, class R>
	class applyDomain<const BlasMatrix<Modular<T>, R > >{
	public:
		typedef Modular<T>                 Field;
		typedef typename Field::Element  Element;
	private:
		const Field & _field ;

		// y = a A x + b y
		template<class _In, class _Out, class Matrix>
		_Out& _apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     , typename ContainerCategories::Vector
			     , typename ContainerCategories::Vector ) const
		{
			// std::cout << "here V" << std::endl;
			// A.write(std::cout,Tag::FileFormat::Maple) << std::endl;
			// std::cout << x << std::endl;


			FFLAS::fgemv((typename Field::Father_t)	_field
				     , (FFLAS::FFLAS_TRANSPOSE)t,
				     A.rowdim(), A.coldim(),
				     _field.one,
				     A.getPointer(), A.getStride(),
				     x.getPointer(), x.getStride(),
				     _field.zero,
				     y.getWritePointer(),y.getStride());

			// FFLAS::fgemm((typename Field::Father_t)	_field
				     // , (FFLAS::FFLAS_TRANSPOSE)t
				     // , FFLAS::FflasNoTrans
				     // , A.rowdim(), 1, A.coldim()
				     // , _field.one,
				     // A.getPointer(), A.getStride(),
				     // x.getPointer(),x.getStride(),
				     // _field.zero,
				     // y.getWritePointer(),y.getStride());

			// std::cout << y << std::endl;

			return y ;
		}

		// y = a A x + b y
		template<class _In, class _Out, class Matrix>
		_Out& _apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     , typename ContainerCategories::Matrix
			     , typename ContainerCategories::Matrix
			     ) const
		{
			// needs to be debuged (trans?). Never Used.
			// std::cout << "here M" << std::endl;
			FFLAS::fgemm((typename Field::Father_t)	_field
				     , (FFLAS::FFLAS_TRANSPOSE)t
				     , FFLAS::FflasNoTrans
				     , A.rowdim(), x.coldim(), A.coldim()
				     , _field.one,
				     A.getPointer(), A.getStride(),
				     x.getPointer(),x.getStride(),
				     _field.zero,
				     y.getWritePointer(),y.getStride());

			return y ;
		}

	public:
		applyDomain(const Field & F) :
			_field(F)
		{}

		/*! Operation \f$y gets a A x + b y\f$ or \f$y gets a A^t x + b y\f$ .
		 * General axpy with x and y vectors or matrices
		 */
		template<class _In, class _Out, class Matrix>
		_Out& apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     // , Method & m
			   ) const
		{
			// std::cout << "called Modular<T>" << std::endl;
			return _apply(t,y,a,A,b,x
				      ,typename ContainerTraits<_In> ::ContainerCategory()
				      ,typename ContainerTraits<_Out>::ContainerCategory() );
		}

	};

	template<class T, class R>
	class applyDomain<const BlasMatrix<ModularBalanced<T>, R > >{
	public:
		typedef ModularBalanced<T>             Field;
		typedef typename Field::Element      Element;
	private:
		const Field & _field ;

		// y = a A x + b y
		template<class _In, class _Out, class Matrix>
		_Out& _apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     , typename ContainerCategories::Vector
			     , typename ContainerCategories::Vector
			     ) const
		{
			// std::cout << "here V" << std::endl;
			FFLAS::fgemv((typename Field::Father_t)	_field, (FFLAS::FFLAS_TRANSPOSE)t,
				     A.rowdim(), A.coldim(),
				     _field.one,
				     A.getPointer(), A.getStride(),
				     x.getPointer(),x.getStride(),
				     // &x[0],ldx,
				     _field.zero,
				     y.getWritePointer(),y.getStride());

			return y ;
		}

		// y = a A x + b y
		template<class _In, class _Out, class Matrix>
		_Out& _apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     , typename ContainerCategories::Matrix
			     , typename ContainerCategories::Matrix
			     ) const
		{
			// needs to be debuged (trans?). Never Used.
			// std::cout << "here M" << std::endl;
			FFLAS::fgemm((typename Field::Father_t)	_field
				     , (FFLAS::FFLAS_TRANSPOSE)t
				     , FFLAS::FflasNoTrans
				     , A.rowdim(), x.coldim(), A.coldim()
				     , _field.one,
				     A.getPointer(), A.getStride(),
				     x.getPointer(),x.getStride(),
				     _field.zero,
				     y.getWritePointer(),y.getStride());

			return y ;
		}

	public:
		applyDomain(const Field & F) :
			_field(F)
		{}

		/*! Operation \f$y gets a A x + b y\f$ or \f$y gets a A^t x + b y\f$ .
		 * General axpy with x and y vectors or matrices
		 */
		template<class _In, class _Out, class Matrix>
		_Out& apply( LINBOX_enum(Tag::Transpose) t
			     , _Out& y
			     , const Element & a
			     , const Matrix & A
			     , const Element & b
			     , const _In & x
			     // , Method & m
			   ) const
		{
			// std::cout << "called ModularBalanced<T>" << std::endl;
			return _apply(t,y,a,A,b,x
				      ,typename ContainerTraits<_In >::ContainerCategory()
				      ,typename ContainerTraits<_Out>::ContainerCategory() );
		}

	};

#if 0
	template<class _Field, class _Rep>
	class applyDomain<BlasMatrix<_Field,_Rep> > {
		// dispatch according to field.
		// if modular and p small -> ModularBalanced<float>
		// if modular and p medium -> ModularBalanced<float>
		// default to triple loop in MatrixDomain
	};

	template<class _Field, class _Rep>
	class applyDomain<BlasMatrix<GivaroExtension<_Field>,_Rep> > {
		// if large enough: Toom Cook
	};

	template<class _Rep>
	class applyDomain<BlasMatrix<PID_integer,_Rep> > {
		// if large : flint
	};

	template<class T, class _Rep>
	class applyDomain<BlasMatrix<Unparametric<T>,_Rep> > {
		// if T == float/double -> BLAS

	};

	template<class _Rep>
	class applyDomain<BlasMatrix<Modular<double>,_Rep> > {
		// Blas
	};

	template<class _Rep>
	class applyDomain<BlasMatrix<Modular<float>,_Rep> > {
		// Blas
	};
#endif


} // LinBox

#endif

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
