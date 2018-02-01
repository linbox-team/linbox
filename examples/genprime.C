#include "linbox/linbox-config.h"
#include <iostream>
#include "stdlib.h"
#include "linbox/randiter/random-prime.h"

using namespace LinBox;

int main(int argc, char**argv){

    PrimeIterator<IteratorCategories::HeuristicTag> gen(7);

    for (int i=0; i < atoi(argv[1]);i++){
        integer p = *gen;
        std::cout<<p<<std::endl;
        ++gen;
    }
}
