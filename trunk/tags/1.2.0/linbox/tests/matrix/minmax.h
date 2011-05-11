namespace LinBox{

template<class Ring, class Matrix>
Matrix& minmat(const Ring& R, Matrix& Mat, size_t n){
	size_t val;
	typename Ring::Element tmp;

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			val = i < j ? i : j;
			R.init(tmp, val+1);
			Mat.setEntry(i, j, tmp);
		}
	}

	return Mat;
}

template<class Ring, class Matrix>
Matrix& maxmat(const Ring& R, Matrix& Mat, size_t n){
	size_t val;
	typename Ring::Element tmp;

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			val = i > j ? i : j;
			R.init(tmp, val+1);
			Mat.setEntry(i, j, tmp);
		}
	}

	return Mat;
}

template<class Ring, class Matrix>
Matrix& qlehmer(const Ring& R, Matrix& Mat, size_t n){
	size_t val;
	typename Ring::Element tmp;

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			val = i > j ? i : j;
			val++;
			R.init(tmp, val*val);
			Mat.setEntry(i, j, tmp);
		}
	}

	return Mat;
}

}; //LinBox
