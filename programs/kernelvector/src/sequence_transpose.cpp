#include "linbox/linbox-config.h"
#include <iostream>
#include <givaro/zring.h>
#include <givaro/givinteger.h>
#include "linbox/util/matrix-stream.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"

#include "common_kernel.hpp"

using namespace LinBox;

int main(int argc, char ** argv) 
{
    std::cerr << "command-line:";
    for(int i = 0; i < argc; ++i)
       std::cerr << ' ' << argv[i];
    std::cerr << std::endl;

    
    std::ifstream input (argv[1]);
    Givaro::ZRing<Givaro::Integer> ZZ;
    MatrixStream<Givaro::ZRing<Givaro::Integer>> ms( ZZ, input );

    Givaro::Timer chrono; chrono.start();
    SparseMatrix<Givaro::ZRing<Givaro::Integer>, SparseMatrixFormat::SparseSeq> A ( ms );
    chrono.stop();
    
    std::cerr << "[READ]: " << A.rowdim() << 'x' << A.coldim() << ' ' << chrono << std::endl;

    chrono.start();
    SparseMatrix<Givaro::ZRing<Givaro::Integer>, SparseMatrixFormat::SparseSeq> AT(ZZ, A.coldim(), A.rowdim());
    

    for(auto it = A.IndexedBegin(); it != A.IndexedEnd(); ++it)
        AT.setEntry(it.colIndex(),it.rowIndex(), it.value());

    chrono.stop();
    std::cerr << "[TRSP]: " << AT.rowdim() << 'x' << AT.coldim() << ' ' << chrono << std::endl;



    chrono.start();
    AT.write(std::cout, Tag::FileFormat::Guillaume);
    chrono.stop();
    
    std::cerr << "[WRIT]: " << AT.rowdim() << 'x' << AT.coldim() << ' ' << chrono << std::endl;
    return 0;
}

    
    
