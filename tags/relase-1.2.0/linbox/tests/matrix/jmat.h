namespace LinBox{

size_t lambda = 1000000000;

template<class Ring, class Matrix>
Matrix& jordanform(const Ring& R, Matrix& Mat, size_t n){
	size_t val = 1;
	typename Ring::Element one, lam;
	R.init(one, val);
	R.init(lam, lambda);

	Mat.setEntry(0, 0, one);

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			if(i == j) Mat.setEntry(i, j, lam);
			if(j+1 == i) Mat.setEntry(i, j, one);
		}
	}

	return Mat;
}

}; //LinBox
