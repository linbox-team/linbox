#include <iostream.h>
#include <vector.h>
#include "LinBox/param_modular.h" // field
//#include "field_padic.h" // field - needs to be put in linbox.
#include "LinBox/diagonal.h" // matrix
#include "LinBox/hilbert.h" // matrix
//#include "LinBox/compose.h" // matrix

using namespace LinBox;
main()
{
// The field
	cout << " Field: ";

	typedef param_modular Field;
	cout << "param_modular(23)" << endl;
	Field	F(23);

	//typedef Field_padic Field;
	//cout << "Field_padic(23, 2)" << endl;
	//Field	F(23, 2);

// The vector
	cout << " Vector: ";

	typedef vector < Field::element > Vect;
	int	n = 10;
	cout << "vector< Field::element >(" << n << ")" << endl;
	Vect    V(n);
	for (int i = 0; i < n; ++i) { F.init(V[i], i); }
	for (int i = 0; i < n; ++i) cout << V[i] << " ";
	cout << endl;

// The blackbox
	cout << " Blackbox: ";

	////////Diagonal
	//Vect      d(n);
	//for (int i = 0; i < n; ++i) { F.init(d[i], i); }
	//diagonal < Field, Vect > D(F, d);
	//cout << "diagonal(";
	//for (int i = 0; i < n; ++i) cout << d[i] << " ";
	//cout << ")" << endl;
	//Blackbox_archetype<Vect>& A = D;
	//// (I don't see why this is necessary)


	////////Hilbert
	hilbert< Field, Vect > H(F, n);
	Blackbox_archetype<Vect>& A = H;

	////////Composed
	//compose<Vect> A(D, D);

// The application
	cout << " apply:";

	Vect      W(n);
	//W = A.apply(V);  
	//A.apply(W, V);  
	A.apply(W, V, &n);
	//A.apply(W, V, 0);

	cout << " :applied" << endl;

	for (int i = 0; i < n; ++i) cout << W[i] << " ";
	cout << endl;

// The transpose application
	cout << " applyTranspose:";

	A.applyTranspose(V, W);
	cout << " :applied" << endl;
	for (int i = 0; i < n; ++i) cout << V[i] << " ";
	cout << endl;

}
