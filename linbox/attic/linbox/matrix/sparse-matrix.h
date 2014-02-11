#if 0 /*  not updated and to be cleaned */
/*! @internal
 * write for CSR format.
 * @bug wrong.
 */
std::ostream & writeSpecialized(std::ostream &os,
		SparseFileFormat::CSR) const
{
	os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
	size_t lig = 0 ;
	size_t i = 0 ;
	while(i < size()) {
		while(lig < _start[i]) {
			os << "-1" << std::endl;
			++lig ;
		}
		while (lig == _start[i]) {
			field().write(_data[i], os << _colid[i] << ' ') << std::endl;
			++i;
		}
		++lig ;
	}
	return os << "0 0 0" << std::endl;
}

/*! @internal
 * write for COO format.
 */
std::ostream & writeSpecialized(std::ostream &os,
		SparseFileFormat::COO) const
{
	os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
	for (size_t i = 0 ; i < rowdim() ; ++i)
		for (size_t j = _start[i] ; j < _start[j+1] ; ++j)
			field().write(_data[j], os << i << ' ' << _colid[j] << ' ') << std::endl;

	return os << "0 0 0" << std::endl;
}



/*! @internal
 * Read for CSR format.
 */
std::istream & readSpecialized(std::istream &is,
		SparseFileFormat::CSR)
{
	size_t nnz = 0 ;
	bool sms = true ;
	std::string firstLine ;
	std::string x ;
	getline(is, firstLine);
	std::istringstream line(firstLine);
	line >> _rownb >> _colnb >> x ;
	size_t mem = 10 ;
	if (!_rownb || _colnb)
		throw LinBoxError("bad input");
	if (x.empty() || x.compare("M")) {  /* SMS */
		// mem = m ;
		_start.reserve(mem);
		_colid.reserve(mem);
		_data.reserve(mem);
	}
	else { /* SMF */
		sms = false ;
		std::istringstream (x) >> nnz ;
		if (!nnz)
			throw LinBoxError("bad input");
		mem = nnz ;
		_start.reserve(nnz);
		_colid.reserve(nnz);
		_data.reserve(nnz);
	}
	Element z ;
	if (sms) { /*  SMS */
		size_t lig = 0 ;
		nnz = 0 ;
		int n ;
		while (is>>n) {
			if (n == 0)
				break;
			while (n == -1) {
				++lig ;
				is >> n ;
			}
			field().read(is,z)  ;
			if (n<0 || lig >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			if (!field().isZero(z)){
				if (mem == nnz) {
					mem+=20 ;
					_start.resize(mem);
					_colid.resize(mem);
					_data.resize (mem);
				}

				_start[nnz]= lig ;
				_colid[nnz]= n ;
				_data[nnz] = z ;
				++nnz ;
			}
		}
		_start.resize(nnz);
		_colid.resize(nnz);
		_data.resize (nnz);

	}
	else { /*  SMF */
		size_t lig = 0 ;
		int n ;
		size_t loc = 0;
		while (is>>n) {
			if (n == 0)
				break;
			while (n == -1) {
				++lig ;
				is >> n ;
			}
			field().read(is,z)  ;
			if (n<0 || lig >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			if (!field().isZero(z)){
				_start[loc]= lig ;
				_colid[loc]= n ;
				_data[loc] = z ;
				++loc ;
			}
		}
		if (loc > nnz)
			throw LinBoxError("bad input");
		_start.resize(loc);
		_colid.resize(loc);
		_data.resize (loc);
	}
	return is ;
}

/*! @internal
 * Read for COO format.
 */
std::istream & readSpecialized(std::istream &is,
		SparseFileFormat::COO)
{
	size_t nnz = 0;
	bool sms = true ;
	std::string firstLine ;
	getline(is, firstLine);
	std::istringstream line(firstLine);
	std::string x ;
	line >> _rownb >> _colnb >> x ;
	size_t mem  = 10 ;
	if (!_rownb || _colnb)
		throw LinBoxError("bad input");
	if (x.empty() || x.compare("M")) {  /* SMS */
		// mem = m ;
		_start.reserve(mem);
		_colid.reserve(mem);
		_data.reserve(mem);
	}
	else { /* SMF */
		sms = false ;
		std::istringstream (x) >> nnz ;
		if (!nnz)
			throw LinBoxError("bad input");
		mem = nnz ;
		_start.reserve(nnz);
		_colid.reserve(nnz);
		_data.reserve(nnz);
	}
	Element z ;
	if (sms) { /*  SMS */
		// size_t lig = 0 ;
		nnz = 0 ;
		int m,n ;
		while (is>>m >> n) {
			if (m == 0 && n == 0)
				break;
			if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			field().read(is,z)  ;
			if (!field().isZero(z)){
				if (mem == nnz) {
					mem+=20 ;
					_start.resize(mem);
					_colid.resize(mem);
					_data.resize (mem);
				}
				_start[nnz]= m ;
				_colid[nnz]= n ;
				_data[nnz] = z ;
				++nnz ;
			}
		}
		_start.resize(nnz);
		_colid.resize(nnz);
		_data.resize (nnz);

	}
	else { /*  SMF */
		size_t loc = 0 ;
		int m,n ;
		while (is>>m >> n) {
			if (m == 0 && n == 0)
				break;
			if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			field().read(is,z)  ;
			if (!field().isZero(z)){
				_start[loc]= m ;
				_colid[loc]= n ;
				_data[loc] = z ;
				++loc ;
			}
		}

		if (loc > nnz)
			throw LinBoxError("bad input");
		_start.resize(loc);
		_colid.resize(loc);
		_data.resize (loc);

	}
	return is ;
}
#endif

#if 0 /*  not updated and to be cleaned */
/*! @internal
 * write for CSR format.
 */
std::ostream & writeSpecialized(std::ostream &os,
		SparseFileFormat::CSR) const
{
	os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
	size_t lig = 0 ;
	size_t i = 0 ;
	while(i < size()) {
		while(lig < _rowid[i]) {
			os << "-1" << std::endl;
			++lig ;
		}
		while (lig == _rowid[i]) {
			field().write(_data[i], os << _colid[i] << ' ') << std::endl;
			++i;
		}
		++lig ;
	}
	return os << "0 0 0" << std::endl;
}

/*! @internal
 * write for COO format.
 */
std::ostream & writeSpecialized(std::ostream &os,
		SparseFileFormat::COO) const
{
	os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
	size_t i = 0 ;
	while(i < size()) {
		field().write(_data[i], os << _rowid[i] << ' ' << _colid[i] << ' ') << std::endl;
		++i;
	}
	return os << "0 0 0" << std::endl;
}



/*! @internal
 * Read for CSR format.
 */
std::istream & readSpecialized(std::istream &is,
		SparseFileFormat::CSR)
{
	size_t nnz = 0 ;
	bool sms = true ;
	std::string firstLine ;
	std::string x ;
	getline(is, firstLine);
	std::istringstream line(firstLine);
	line >> _rownb >> _colnb >> x ;
	size_t mem = 10 ;
	if (!_rownb || _colnb)
		throw LinBoxError("bad input");
	if (x.empty() || x.compare("M")) {  /* SMS */
		// mem = m ;
		_rowid.reserve(mem);
		_colid.reserve(mem);
		_data.reserve(mem);
	}
	else { /* SMF */
		sms = false ;
		std::istringstream (x) >> nnz ;
		if (!nnz)
			throw LinBoxError("bad input");
		mem = nnz ;
		_rowid.reserve(nnz);
		_colid.reserve(nnz);
		_data.reserve(nnz);
	}
	Element z ;
	if (sms) { /*  SMS */
		size_t lig = 0 ;
		nnz = 0 ;
		int n ;
		while (is>>n) {
			if (n == 0)
				break;
			while (n == -1) {
				++lig ;
				is >> n ;
			}
			field().read(is,z)  ;
			if (n<0 || lig >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			if (!field().isZero(z)){
				if (mem == nnz) {
					mem+=20 ;
					_rowid.resize(mem);
					_colid.resize(mem);
					_data.resize (mem);
				}

				_rowid[nnz]= lig ;
				_colid[nnz]= n ;
				_data[nnz] = z ;
				++nnz ;
			}
		}
		_rowid.resize(nnz);
		_colid.resize(nnz);
		_data.resize (nnz);

	}
	else { /*  SMF */
		size_t lig = 0 ;
		int n ;
		size_t loc = 0;
		while (is>>n) {
			if (n == 0)
				break;
			while (n == -1) {
				++lig ;
				is >> n ;
			}
			field().read(is,z)  ;
			if (n<0 || lig >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			if (!field().isZero(z)){
				_rowid[loc]= lig ;
				_colid[loc]= n ;
				_data[loc] = z ;
				++loc ;
			}
		}
		if (loc > nnz)
			throw LinBoxError("bad input");
		_rowid.resize(loc);
		_colid.resize(loc);
		_data.resize (loc);
	}
	return is ;
}

/*! @internal
 * Read for COO format.
 */
std::istream & readSpecialized(std::istream &is,
		SparseFileFormat::COO)
{
	size_t nnz = 0;
	bool sms = true ;
	std::string firstLine ;
	getline(is, firstLine);
	std::istringstream line(firstLine);
	std::string x ;
	line >> _rownb >> _colnb >> x ;
	size_t mem  = 10 ;
	if (!_rownb || _colnb)
		throw LinBoxError("bad input");
	if (x.empty() || x.compare("M")) {  /* SMS */
		// mem = m ;
		_rowid.reserve(mem);
		_colid.reserve(mem);
		_data.reserve(mem);
	}
	else { /* SMF */
		sms = false ;
		std::istringstream (x) >> nnz ;
		if (!nnz)
			throw LinBoxError("bad input");
		mem = nnz ;
		_rowid.reserve(nnz);
		_colid.reserve(nnz);
		_data.reserve(nnz);
	}
	Element z ;
	if (sms) { /*  SMS */
		// size_t lig = 0 ;
		nnz = 0 ;
		int m,n ;
		while (is>>m >> n) {
			if (m == 0 && n == 0)
				break;
			if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			field().read(is,z)  ;
			if (!field().isZero(z)){
				if (mem == nnz) {
					mem+=20 ;
					_rowid.resize(mem);
					_colid.resize(mem);
					_data.resize (mem);
				}
				_rowid[nnz]= m ;
				_colid[nnz]= n ;
				_data[nnz] = z ;
				++nnz ;
			}
		}
		_rowid.resize(nnz);
		_colid.resize(nnz);
		_data.resize (nnz);

	}
	else { /*  SMF */
		size_t loc = 0 ;
		int m,n ;
		while (is>>m >> n) {
			if (m == 0 && n == 0)
				break;
			if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			field().read(is,z)  ;
			if (!field().isZero(z)){
				_rowid[loc]= m ;
				_colid[loc]= n ;
				_data[loc] = z ;
				++loc ;
			}
		}

		if (loc > nnz)
			throw LinBoxError("bad input");
		_rowid.resize(loc);
		_colid.resize(loc);
		_data.resize (loc);

	}
	return is ;
}
#endif

std::ostream & writeSpecialized(std::ostream &os,
		LINBOX_enum(Tag::FileFormat) format) const
{
	switch (format) {
		case (Tag::FileFormat::Maple):
			{

				linbox_check(_colnb > 0);
				os << "[";
				bool firstrow=true;
				size_t idx = 0 ;

				linbox_check(_rownb+1 == _start.size());
				linbox_check(_nbnz == _data.size());
				linbox_check(_nbnz == _colid.size());
				linbox_check(_start[_rownb] == _nbnz);
				for (size_t i = 0 ; i < _rownb ; ++i ) {
					if (firstrow) {
						os << "[";
						firstrow =false;
					}
					else
						os << ", [";


					for (size_t j = 0; j < _colnb ; ++j) {
						if (idx == _nbnz)
							field().write (os, field().zero);
						else if (_colid[idx] == j &&
								_start[i] <= idx && idx < _start[i+1]) {
							field().write (os, _data[idx]);
							++idx;
						}
						else {
							field().write (os, field().zero);
						}

						if (j < _colnb - 1)
							os << ", ";
					}

					os << " ]";
				}

				os << "]";
				linbox_check(idx == _nbnz);

				break;
			}
		default :
			os << "I don't know" << std::endl;

	}
	return os ;

}

std::ostream & writeSpecialized(std::ostream &os,
		LINBOX_enum(Tag::FileFormat) format) const
{
	switch (format) {
		case (Tag::FileFormat::Maple):
			{

				linbox_check(_colnb > 0);
				os << "[";
				bool firstrow=true;
				size_t idx = 0 ;

				linbox_check(_nbnz == _rowid.size());
				linbox_check(_nbnz == _data.size());
				linbox_check(_nbnz == _colid.size());
				for (size_t i = 0 ; i < _rownb ; ++i ) {
					if (firstrow) {
						os << "[";
						firstrow =false;
					}
					else
						os << ", [";


					for (size_t j = 0; j < _colnb ; ++j) {
						if (idx == _nbnz)
							field().write (os, field().zero);
						else if (_colid[idx] == j &&
								_rowid[idx] ==i) {
							field().write (os, _data[idx]);
							++idx;
						}
						else {
							field().write (os, field().zero);
						}

						if (j < _colnb - 1)
							os << ", ";
					}

					os << " ]";
				}

				os << "]";
				linbox_check(idx == _nbnz);

				break;
			}
		default :
			os << "I don't know" << std::endl;

	}
	return os ;

}

std::ostream & writeSpecialized(std::ostream &os,
		LINBOX_enum(Tag::FileFormat) format) const
{
	SparseMatrix<Field,SparseMatrixFormat::CSR> Temp(field());
	this->exporte(Temp);
	Temp.write(os,format);
	return os ;

}

#if 0 /*  not updated and to be cleaned */
/*! @internal
 * write for CSR format.
 * @bug wrong.
 */
std::ostream & writeSpecialized(std::ostream &os,
		SparseFileFormat::CSR) const
{
	os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
	size_t lig = 0 ;
	size_t i = 0 ;
	while(i < size()) {
		while(lig < _start[i]) {
			os << "-1" << std::endl;
			++lig ;
		}
		while (lig == _start[i]) {
			field().write(_data[i], os << _colid[i] << ' ') << std::endl;
			++i;
		}
		++lig ;
	}
	return os << "0 0 0" << std::endl;
}

/*! @internal
 * write for COO format.
 */
std::ostream & writeSpecialized(std::ostream &os,
		SparseFileFormat::COO) const
{
	os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
	for (size_t i = 0 ; i < rowdim() ; ++i)
		for (size_t j = _start[i] ; j < _start[j+1] ; ++j)
			field().write(_data[j], os << i << ' ' << _colid[j] << ' ') << std::endl;

	return os << "0 0 0" << std::endl;
}



/*! @internal
 * Read for CSR format.
 */
std::istream & readSpecialized(std::istream &is,
		SparseFileFormat::CSR)
{
	size_t nnz = 0 ;
	bool sms = true ;
	std::string firstLine ;
	std::string x ;
	getline(is, firstLine);
	std::istringstream line(firstLine);
	line >> _rownb >> _colnb >> x ;
	size_t mem = 10 ;
	if (!_rownb || _colnb)
		throw LinBoxError("bad input");
	if (x.empty() || x.compare("M")) {  /* SMS */
		// mem = m ;
		_start.reserve(mem);
		_colid.reserve(mem);
		_data.reserve(mem);
	}
	else { /* SMF */
		sms = false ;
		std::istringstream (x) >> nnz ;
		if (!nnz)
			throw LinBoxError("bad input");
		mem = nnz ;
		_start.reserve(nnz);
		_colid.reserve(nnz);
		_data.reserve(nnz);
	}
	Element z ;
	if (sms) { /*  SMS */
		size_t lig = 0 ;
		nnz = 0 ;
		int n ;
		while (is>>n) {
			if (n == 0)
				break;
			while (n == -1) {
				++lig ;
				is >> n ;
			}
			field().read(is,z)  ;
			if (n<0 || lig >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			if (!field().isZero(z)){
				if (mem == nnz) {
					mem+=20 ;
					_start.resize(mem);
					_colid.resize(mem);
					_data.resize (mem);
				}

				_start[nnz]= lig ;
				_colid[nnz]= n ;
				_data[nnz] = z ;
				++nnz ;
			}
		}
		_start.resize(nnz);
		_colid.resize(nnz);
		_data.resize (nnz);

	}
	else { /*  SMF */
		size_t lig = 0 ;
		int n ;
		size_t loc = 0;
		while (is>>n) {
			if (n == 0)
				break;
			while (n == -1) {
				++lig ;
				is >> n ;
			}
			field().read(is,z)  ;
			if (n<0 || lig >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			if (!field().isZero(z)){
				_start[loc]= lig ;
				_colid[loc]= n ;
				_data[loc] = z ;
				++loc ;
			}
		}
		if (loc > nnz)
			throw LinBoxError("bad input");
		_start.resize(loc);
		_colid.resize(loc);
		_data.resize (loc);
	}
	return is ;
}

/*! @internal
 * Read for COO format.
 */
std::istream & readSpecialized(std::istream &is,
		SparseFileFormat::COO)
{
	size_t nnz = 0;
	bool sms = true ;
	std::string firstLine ;
	getline(is, firstLine);
	std::istringstream line(firstLine);
	std::string x ;
	line >> _rownb >> _colnb >> x ;
	size_t mem  = 10 ;
	if (!_rownb || _colnb)
		throw LinBoxError("bad input");
	if (x.empty() || x.compare("M")) {  /* SMS */
		// mem = m ;
		_start.reserve(mem);
		_colid.reserve(mem);
		_data.reserve(mem);
	}
	else { /* SMF */
		sms = false ;
		std::istringstream (x) >> nnz ;
		if (!nnz)
			throw LinBoxError("bad input");
		mem = nnz ;
		_start.reserve(nnz);
		_colid.reserve(nnz);
		_data.reserve(nnz);
	}
	Element z ;
	if (sms) { /*  SMS */
		// size_t lig = 0 ;
		nnz = 0 ;
		int m,n ;
		while (is>>m >> n) {
			if (m == 0 && n == 0)
				break;
			if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			field().read(is,z)  ;
			if (!field().isZero(z)){
				if (mem == nnz) {
					mem+=20 ;
					_start.resize(mem);
					_colid.resize(mem);
					_data.resize (mem);
				}
				_start[nnz]= m ;
				_colid[nnz]= n ;
				_data[nnz] = z ;
				++nnz ;
			}
		}
		_start.resize(nnz);
		_colid.resize(nnz);
		_data.resize (nnz);

	}
	else { /*  SMF */
		size_t loc = 0 ;
		int m,n ;
		while (is>>m >> n) {
			if (m == 0 && n == 0)
				break;
			if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
				throw LinBoxError("bad input");
			field().read(is,z)  ;
			if (!field().isZero(z)){
				_start[loc]= m ;
				_colid[loc]= n ;
				_data[loc] = z ;
				++loc ;
			}
		}

		if (loc > nnz)
			throw LinBoxError("bad input");
		_start.resize(loc);
		_colid.resize(loc);
		_data.resize (loc);

	}
	return is ;
}
#endif

