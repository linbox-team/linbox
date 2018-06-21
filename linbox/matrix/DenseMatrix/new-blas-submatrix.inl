/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *               2018 revamped by Pascal Giorgi 
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@lirmm.fr
 *               Clément Pernet clement.pernet@imag.fr
 *               Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file matrix/densematrix/new-blas-submatrix.inl
 * @ingroup densematrix
 */


namespace LinBox {

    //////////////////
    // CONSTRUCTORS //
    //////////////////

    /*  constructors */
    
    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (typename BlasSubmatrix<_Matrix>::matrixType &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim) :
        _ptr(M.getPointer()+rowbeg*M.getStride()+colbeg),
        _row (Rowdim),
        _col(Coldim),
		_stride(M.getStride()),
        _field(M.field()),
		,_AD(M.field())
	{
        linbox_check ( rowbeg  <= M.rowdim() ); // allow for NULL matrix
		linbox_check ( colbeg  <= M.coldim() );
        linbox_check ( rowbeg+Rowdim <= M.rowdim() );
        linbox_check ( colbeg+Coldim <= M.coldim() );
    }

    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (BlasSubmatrix<_Matrix>  &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim) :
        _ptr(M.getPointer()+rowbeg*M.getStride()+colbeg),
        _row (Rowdim),
        _col(Coldim),
		_stride(M.getStride()),
        _field(M.field()),
		,_AD(M.field())
	{
        linbox_check ( rowbeg  <= M.rowdim() ); // allow for NULL matrix
		linbox_check ( colbeg  <= M.coldim() );
        linbox_check ( rowbeg+Rowdim <= M.rowdim() );
        linbox_check ( colbeg+Coldim <= M.coldim() );
    }

    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (typename BlasSubmatrix<_Matrix>::matrixType &M):        
        _ptr(M.getPointer()),
        _row (M.rowdim()),
        _col (M.coldim()),
		_stride(M.getStride()),
        _field(M.field()),
		,_AD(M.field())
    {

    }

  

    //////////////////
    //   ELEMENTS   //
    //////////////////

    template < class _Matrix >
    const Element& BlasSubmatrix<_Matrix>::setEntry (size_t i, size_t j, const Element &a_ij) {
        return field().assign(a_ij, _ptr+i*_stride+j);
    }

    template < class _Matrix >
    Element& BlasSubmatrix<_Matrix>::refEntry (size_t i, size_t j) {
        return _ptr+i*_stride+j;
    }

    template < class _Matrix >
    const Element& BlasSubmatrix<_Matrix>::getEntry (size_t i, size_t j) const {
        return _ptr+i*_stride+j;
    }

    template < class _Matrix >
    Element& BlasSubmatrix<_Matrix>::getEntry (Element &x, size_t i, size_t j) const {
        return field().assign(x,_ptr+i*_stride+j);
    }




    ///////////////////
    //      I/O      //
    ///////////////////

    template < class _Matrix >
    std::istream& BlasSubmatrix<_Matrix>::read (std::istream &file)
    {
        // must be improved maybe to allow any matrix format
		Iterator p;
		int m,n;
		char c;
		file>>m>>n>>c;        		
        linbox_check (m == _row);
        linbox_check (m == _col);
		
		if ((c != 'M') && (c != 'm')) {
            for (p = Begin (); p != End (); ++p) {
                field()read (file, *p);
			}
		}
		else { // sparse file format
			int i=0, j=0;
			while (true)
                {
                    file >> i >> j;
                    if (i+j <= 0) break;
                    field()read (file, _ptr+ (i-1)*_stride+(j-1));
                }
		}

		return file;
    }

    template < class _Matrix >
    std::ostream& BlasSubmatrix<_Matrix>::write (std::ostream &os,LINBOX_enum (Tag::FileFormat) f = Tag::FileFormat::MatrixMarket ) const {
        
		ConstRowIterator p;
		switch(f) {
        case (Tag::FileFormat::MatrixMarket ) : /* Matrix Market */
            {
                writeMMArray(os, *this, "BlasSubmatrix");
            }
            break;
        case (Tag::FileFormat::Plain) : /*  raw output */
            {
                integer c;
                int wid;

                field()cardinality (c);

                if (c >0)
                    wid = (int) ceil (log ((double) c) / M_LN10);
                else {
                    wid=1000;
                }
                    
                for (p = rowBegin (); p != rowEnd ();++p) {
                    typename ConstRow::const_iterator pe;

                    os << "  [ ";

                    for (pe = p->begin (); pe != p->end (); ++pe) {
                        os.width (wid);
                        field()write (os, *pe);
                        os << " ";
                    }
                    os << "]" << std::endl;
                }
            }
            break;
        case (Tag::FileFormat::Maple) : /*  maple format */
            {

                os << "Matrix( " << rowdim() << ',' << coldim() << ",\n[" ;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;
                    if (p!=rowBegin()) os << ' ';
                    os << "[ ";
                        
                    for (pe = p->begin (); pe != p->end (); ) {
                        field()write (os, *pe);
                        ++pe ;
                        if (pe != p->end())
                            os << ", ";
                    }

                    os << "]" ;
                    ++p ;
                    if (p != rowEnd() )
                        os << ',' << std::endl;;

                }
                os << "])" ;
            }
            break;
        case (Tag::FileFormat::HTML) : /*  HTML format */
            {

                os << "<table border=\"1\">" ;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;

                    os << "<tr>";

                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os<< "<td>", *pe)<<"</td>";
                        ++pe ;
                    }

                    os << "</tr>" << std::endl;
                    ++p ;

                }
                os << "</table>" ;
            }
            break;
        case (Tag::FileFormat::LaTeX) : /*  LaTex format */
            {
                os << "\\begin{pmatrix} " << std::endl;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;

                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os, *pe);
                        ++pe ;
                        if (pe != p->end())
                            os << "& ";
                    }

                    os << "\\\\" << std::endl;
                    ++p ;

                }
                os << "\\end{pmatrix}" ;
            }
            break;
        default : /*  this is an error */
            {
                throw LinBoxError("unknown format to write matrix in");
            }
		}

		return os;
	}

    

    
}



/* // Local Variables: */
/* // mode: C++ */
/* // tab-width: 4 */
/* // indent-tabs-mode: nil */
/* // c-basic-offset: 4 */
/* // End: */
/* // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s */

