#include "densegfrank.h"

int densegfrank(const DenseMatrix<GF2>& A) 
{
    std::cout << "I am the precompiled implementation calling the generic rank" << std::endl;
    return rankimplementation(A);
}
