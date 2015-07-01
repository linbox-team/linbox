namespace LinBox { // namespace in which all LinBox code resides.

/*
template< class Field >
void toeplitz_determinant( ZZ_p& res, const Toeplitz<Field>& A );
*/

template< class PField >
typename PField::Coeff& toeplitz_determinant
	( const PField& F, typename PField::Coeff& res, 
          const typename PField::Element& T, size_t n )
{
	short int sign = 1;
	typename PField::Coeff one, temp;
	const typename PField::CoeffField& CField = F.getCoeffField();
	CField.init(one,1);
	CField.init(res,1);
	typename PField::Element f1, f2( T ), fi;
	F.setCoeff( f1, 2*n - 1, one );
	F.init( fi, one );
	
	while( F.deg(f2) >= n ) {
		F.rem( fi, f1, f2 );
		CField.mulin
			( res, CField.powin( F.leadCoeff( temp, f2 ),
				             F.deg(f1) - F.deg(fi) ) );
		if( !((F.deg(f2)-F.deg(f1))%2) && !((F.deg(f1)-n)%2) )
			sign *= -1;
		f1 = f2;
		f2 = fi;
	}

	if( F.deg(f2) == (n-1) ) {
		CField.mulin
			( res, CField.powin( F.leadCoeff( temp, f2 ),
			                     F.deg(f1) - F.deg(f2) ) );
		if( sign == -1 ) {
			typename PField::Coeff negOne;
			CField.init( negOne, -1 );
			CField.mulin( res, negOne );
		}
	}
	else CField.init( res, 0 );

	return res;
}

} // end of namespace LinBox
