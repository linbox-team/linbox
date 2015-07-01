#ifndef __LINBOX_sepraterank_H
#define __LINBOX_sepraterank_H

#include "rank.h"


template<class Mat>
int rank(const Mat& A) {
    std::cout << "I am the separate rank" << std::endl;
    return rankimplementation(A);
}


#ifdef __LINBOX_SEPARATE_COMPILATION
#include "densegfrank.h"

template<>
int rank(const DenseMatrix<GF2>& A) 
{
    std::cout << "I am the specialization of the separate rank" << std::endl;
    return densegfrank(A);
}

#endif

#endif
