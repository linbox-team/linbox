/* Copyright (C) 2010 LinBox
 *
 * Time-stamp: <16 Dec 10 16:26:18 Jean-Guillaume.Dumas@imag.fr>
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

#include <linbox/algorithms/cra-domain.h>
#include <linbox/field/modular-double.h>
#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/cra-early-multip.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/cra-full-multip-fixed.h"
#include "linbox/randiter/random-prime.h"

struct Interator {
    std::vector<Integer> _v;
    double maxsize;
    
    Interator(const std::vector<Integer>& v) : _v(v), maxsize(0.0) {
        for(std::vector<Integer>::const_iterator it=_v.begin();
            it != _v.end(); ++it) {
            double ds = naturallog(*it);
            maxsize = (maxsize<ds?ds:maxsize);
        }
    }

    Interator(int n, int s) : _v(n), maxsize(0.0) {
        for(std::vector<Integer>::iterator it=_v.begin();
            it != _v.end(); ++it) {
            Integer::random(*it, s);
            double ds = naturallog(*it);
            maxsize = (maxsize<ds?ds:maxsize);
        }
    }
        
    const std::vector<Integer>& getVector() { return _v; }
    const double getLogSize() { return maxsize; }
    
    template<typename Field>
    std::vector<typename Field::Element>& operator()(std::vector<typename Field::Element>& v, const Field& F) const {
        v.resize(_v.size());
        std::vector<Integer>::const_iterator vit=_v.begin();
        typename std::vector<typename Field::Element>::iterator eit=v.begin();
        for( ; vit != _v.end(); ++vit, ++eit)
            F.init(*eit, *vit);
        return v;
    }
};

struct InteratorIt;
namespace LinBox
{
    template<class Element> struct CRATemporaryVectorTrait<InteratorIt , Element> {
        typedef typename std::vector<double>::iterator Type_t;
    };
}

struct InteratorIt : public Interator {
   
    mutable std::vector<double> _C;

    InteratorIt(const std::vector<Integer>& v) : Interator(v), _C(v.size()) {}
    InteratorIt(int n, int s) : Interator(n,s), _C(n) {}
    
    template<typename Iterator, typename Field>
    Iterator& operator()(Iterator& res, const Field& F) const {
        std::vector<Integer>::const_iterator vit=this->_v.begin();
        std::vector<double>::iterator eit=_C.begin();
        for( ; vit != _v.end(); ++vit, ++eit)
            F.init(*eit, *vit);
        return res=_C.begin();
    }
    
};







bool TestCra(int N, int S) {
    
    std::ostream &report = LinBox::commentator.report
                (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

    report << "TestCra(" << N << ',' << S << ')' << std::endl;

    Interator iteration(N, S);
    InteratorIt iterationIt(iteration.getVector());
    LinBox::RandomPrimeIterator genprime( 24 );
    
    bool pass = true;

    report << "ChineseRemainder<EarlyMultipCRA>(4)" << std::endl;
    {
        LinBox::ChineseRemainder< LinBox::EarlyMultipCRA< LinBox::Modular<double> > > craEM(4UL);
        std::vector<Integer> ResEM(N);
        craEM( ResEM, iteration, genprime);
    
        pass = pass && std::equal( ResEM.begin(), ResEM.end(), iteration.getVector().begin() );
    }

    report << "ChineseRemainder<FullMultipCRA>(" << iteration.getLogSize() << ')' << std::endl;
    {
        LinBox::ChineseRemainder< LinBox::FullMultipCRA< LinBox::Modular<double> > > craFM( iteration.getLogSize() );
        std::vector<Integer> ResFM(N);
        craFM( ResFM, iteration, genprime);
    
        pass = pass && std::equal( ResFM.begin(), ResFM.end(), iteration.getVector().begin() );
    }
    
    report << "ChineseRemainder<FullMultipFixedCRA>(" << N << ',' << iterationIt.getLogSize() << ')' << std::endl;
    {
        LinBox::ChineseRemainder< LinBox::FullMultipFixedCRA< LinBox::Modular<double> > > craFMF( std::pair<size_t,double>(N,iterationIt.getLogSize()) );
        std::vector<Integer> ResFMF(N);
        std::vector<Integer>::iterator ResIT= ResFMF.begin();
        craFMF( ResIT, iterationIt, genprime);
    
        pass = pass && std::equal( iterationIt.getVector().begin(), iterationIt.getVector().end(), ResFMF.begin() );
    }
    
    return pass;
}



#include "test-common.h"

int main (int argc, char **argv)
{
        static size_t n = 1000;
        static size_t s = 30;
        static int iterations = 1;

        static Argument args[] = {
                { 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT
, &n },
                { 's', "-s S", "Set size of test integers.", TYPE_INT
, &s },
                { 'i', "-i I", "Perform each test for I iterations.",     TYPE_INT, &iterations },
                { '\0' }
        };

        parseArguments (argc, argv, args);

        commentator.start("CRA-Domain test suite", "CRADom");
        bool pass = true;

        for(int i=0; i<iterations; ++i)
            pass = TestCra(n,s);
        

        commentator.stop("CRA-Domain test suite");
        return pass ? 0 : -1;
}
