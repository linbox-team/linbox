#include<cstdlib>
#include "linbox/integer.h"

//#include<iostream>
//using namespace std;
namespace LinBox{

template<class Ring, class Matrix>
Matrix& invhilb(const Ring& R, Matrix& Mat, int n){
	typedef typename Ring::Element Element;
	
	integer p = n;
	integer r;
	Element val;
	int i, j, i_ind, j_ind;
	for(i_ind=0, i = 1; i <= n; i++, i_ind++){
		if(i>1) p = ((n-i+1)*p*(n+i-1))/((i-1)*(i-1));
		r = p*p;	
		R.init(val, (r/(2*i-1)));
		Mat.setEntry(i_ind, i_ind, val);
		for(j_ind=i, j = i+1; j <= n; j++, j_ind++){
			r = -1*((n-j+1)*r*(n+j-1))/((j-1)*(j-1));
			R.init(val, (r/(i+j-1)));
			Mat.setEntry(i_ind, j_ind, val);
			Mat.setEntry(j_ind, i_ind, val);
		}
	}

	return Mat;
}

}
