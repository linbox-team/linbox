/* Copyright (C) 2007 LinBox
 *
 * Written by Zhendong Wan  <wan@mail.eecis.udel.edu> 
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


#ifndef __LINBOX_rational_solver_adaptive_H
#define __LINBOX_rational_solver_adaptive_H

#include <linbox/field/modular-int32.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/blackbox/dense.h>

namespace LinBox 
{

        // Generic non-numerical solver requires conversion of the vector
        template<class IRing, class OutVector, class InVector>
        struct RationalSolverAdaptiveClass {
            static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, const DenseMatrix<IRing>& M, const InVector& b) {
                linbox_check ((M. rowdim() == M. coldim()) && (b.size() == M.rowdim()) && (num. size() ==M.coldim()));
                typedef Modular<int32> Field;
                RationalSolver<IRing, Field, RandomPrimeIterator, NumericalTraits> numerical_solver;
                SolverReturnStatus ret;
                ret = numerical_solver. solve(num, den, M, b);
                
                if (ret != SS_OK) {
                    RationalSolver<IRing, Field, RandomPrimeIterator> solver;
                    std::vector<typename IRing::Element> Ib; Ib.reserve(b.size());
                    typename IRing::Element tmp;
                    for(typename InVector::const_iterator biter = b.begin();
                        biter != b.end();
                        ++biter) 
                        Ib.push_back( M.field().init(tmp, *biter) );
                    ret = solver. solve(num, den, M, Ib);
                }
                
                return ret;
            }
        };
    
        // Specialization when the vector is already over the ring
        template<class IRing, class OutVector, template<typename T> class Container>
        struct RationalSolverAdaptiveClass<IRing, OutVector, Container<typename IRing::Element> > {
            static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, const DenseMatrix<IRing>& M, const Container<typename IRing::Element> & b) {
                linbox_check ((M. rowdim() == M. coldim()) && (b.size() == M.rowdim()) && (num. size() ==M.coldim()));
                typedef Modular<int32> Field;
                RationalSolver<IRing, Field, RandomPrimeIterator, NumericalTraits> numerical_solver;
                SolverReturnStatus ret;
                ret = numerical_solver. solve(num, den, M, b);
                
                if (ret != SS_OK) {
                    RationalSolver<IRing, Field, RandomPrimeIterator> solver;
                    ret = solver. solve(num, den, M, b);
                }
                
                return ret;
            }
        };


	class RationalSolverAdaptive {
	public:
            template<class IRing, class OutVector, class InVector>
            static SolverReturnStatus solveNonsingular(OutVector& num, typename IRing::Element& den, const DenseMatrix<IRing>& M, const InVector& b) {
                return RationalSolverAdaptiveClass<IRing,OutVector,InVector>::solveNonsingular(num, den, M, b);
            }
	};

}

#endif //__LINBOX_rational_solver_adaptive_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
