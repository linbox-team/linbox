/* tests/test-vector-read-write.C
 */

#include <string>
#include <fstream>
#include "linbox-config.h"
#include "linbox/vector/vector-domain.h"

std::string home("/home/saunders/lieatlas/");

/*** comment?? */
/** argc[1], model (ops and F)
	argv[2], root and rho
	argv[3], nu
*/
template <class Field, class NV>
int read_facet(NV& nu, std::string facetf, const Field& Q, int rank = 8, int index = 1) {
    // the facet
	//std::cerr << "reading facet "<< std::endl;
	std::ifstream facet, muddle;
	muddle.open(facetf.c_str());

   	std::string path;
    if ( muddle ) {
		muddle.close();
		facet.open(facetf.c_str());
	} else {
		path = home + "facets/" + facetf;
		facet.open(path.c_str());
	    if ( ! facet) std::cerr << 
			"failed on filename " << facetf.c_str() << std::endl <<
			"failed on filename " << path << std::endl;
	}
	if (! facet) std::cerr << "can't find facet " << facetf << std::endl;

	LinBox::VectorDomain<Field> VD(Q);
    nu. resize (rank);
	for(int i = 1; i <= index; ++i) VD.read(facet, nu);
	facet. close();

	/*
	std::cout << "Facet: ";
	VD. write (std::cout, nu);
	std::cout << std::endl;
	*/

}

template <class Field, class LRV, class RV>
int read_rootData(LRV& root, RV& rho, 
			  std::string rootf, const Field& Q, int rank = 8) {

	std::ifstream rootData, muddle;
	muddle.open(rootf.c_str());

   	std::string path;
    if ( muddle) {
		muddle.close();
		rootData.open(rootf.c_str());
	} else {
    	//std::string path("rootData/" + rootf);
		path = home + "rootData/" + rootf;
		rootData.open(path.c_str());
	    if ( ! rootData) std::cerr << 
			"failed on filename " << rootf.c_str() << std::endl <<
			"failed on filename " << path<< std::endl;
	}
	if (! rootData) std::cerr << "can't find rootData " << rootf << std::endl;

	LinBox::VectorDomain<Field> VD(Q);
	root. resize (rank);
	rho. resize (rank);

	for (int i = 0; i < rank; ++i) {
		root[i].resize(rank); VD.read(rootData, root[i]);
	}
	VD.read(rootData, rho);
	rootData.close();
	//
	/*
	std::cout << "Roots:\n";
	for (int i = 0; i < rank; ++i) {
		root[i].resize(rank); VD.write(std::cout, root[i]);
	}
	std::cout << std::endl;
	std::cout << "Rho:\n";
	VD.write (std::cout, rho);
	std::cout << std::endl;
	*/
}


template <class Field, class OV, class OP>
int read_model (OV& op, OP*& F, 
			  std::string modelf, const Field& Q, int rank = 8) {

	std::ifstream model, muddle;
	muddle.open(modelf.c_str());

   	std::string path;
    if ( muddle) {
		muddle.close();
		model.open(modelf.c_str());
	} else {
    	//std::string path("models/" + modelf);
		path = home + "models/" + modelf;
		model.open(path.c_str());
		if ( ! model ) 	std::cerr << 
			"failed on filename " << modelf.c_str() << std::endl <<
			"failed on filename " << path.c_str() << std::endl;
	}
	if (! model) std::cerr << "can't find model " << modelf << std::endl;

	LinBox::VectorDomain<Field> VD(Q);
    // read ops and F from model.
	op. resize (rank); 
	for (int i = 0; i < rank; ++i) op[i] = new OP(Q); 
	F = new OP(Q);
	char c;
	do {model.get(c);}  while ( c != '[' );
	for (int i = 0; i < rank; ++i) { op[i]->read(model, OP::FORMAT_MAGMACPT); }
 	F -> read(model, OP::FORMAT_MAGMACPT);

	// check the reads
	/*
	std::cout << "Operators:\n";
	for (int i = 0; i < rank; ++i) { op[i]->write(cout, OP::FORMAT_PRETTY); }
	std::cout << std::endl;
	std::cout << "F:\n";
	F.write(cout, OP::FORMAT_PRETTY);
	std::cout << std::endl;
	*/
}

template <class Field, class NV, class LRV, class RV, class OV, class OP>
int read_lie (NV& nu, LRV& root, RV& rho, OV& op, OP*& F, 
			  const char* modelf, const char* rootf, const char* facetf, 
			  const Field& Q, int rank = 8) {

    read_rootData(root, rho, std::string(rootf), Q, rank);
	read_facet(nu, facetf, Q, rank);
	read_model(op, F, std::string(modelf), Q, rank);
}
