namespace LinBox{

void isPower2(size_t n){
	size_t p = 1;
	while(p < n) p <<= 1;
	if(p != n){} //  THROW ERROR };
}

template<class Ring, class Matrix>
Matrix& hadamard(const Ring& R, Matrix& Mat, size_t n){
	size_t val = 1;
	typename Ring::Element tmp;

	Mat.setEntry(0, 0, R.init(tmp, val));

	for(size_t k = 1; k < n; k*=2){
		for(size_t i = 0; i < k; ++i){
			for(size_t j = 0; j < k; ++j){
				Mat.getEntry(tmp, i, j);
				Mat.setEntry(i+k, j, tmp);
				Mat.setEntry(i, j+k, tmp);
				Mat.setEntry(i+k, j+k, R.negin(tmp));
			}
		}
	}

	return Mat;
}

}
