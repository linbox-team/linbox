#include <iostream>
#include "gf2.h"
#include "mat.h"
#include "separate.h"


int main() 
{

	DenseMatrix<GF2> A;
        
        int r = rank(A);
        
        std::cout << "the rank is " << r << std::endl;
        

        return 0;
}



