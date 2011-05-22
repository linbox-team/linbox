/*
 * Written by Zhendong Wan  <wan@mail.eecis.udel.edu> 
 * Time-stamp: <12 Mar 07 19:45:32 Jean-Guillaume.Dumas@imag.fr> 
 */

#ifndef __LINBOX_RATIONAL_SOLVER_ADAPTIVE_H
#define __LINBOX_RATIONAL_SOLVER_ADAPTIVE_H
#include <linbox/field/modular-int32.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/blackbox/dense.h>

namespace LinBox {

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

#endif
