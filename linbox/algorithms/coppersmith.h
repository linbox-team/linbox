/* linbox/algorithms/coppersmith.h
 * evolved from block-wiedemann.h by George Yuhasz
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


#ifndef __LINBOX_coppersmith_H
#define __LINBOX_coppersmith_H

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <givaro/givpoly1crt.h>


#include "linbox/integer.h"
#include "linbox/util/commentator.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/solutions/det.h"

#include "linbox/util/error.h"
#include "linbox/util/debug.h"

namespace LinBox
{

	template <class _Domain>
	class CoppersmithSolver{

	public:
		typedef _Domain 			Domain;
		typedef typename Domain::Field                    Field;
		typedef typename Domain::Element       Element;
		typedef typename Domain::OwnMatrix 	Block;
		typedef typename Domain::Matrix 	Sub;

		inline const Domain & domain() const { return *_MD; }
		inline const Field & field() const { return domain().field(); }
	protected:
		const Domain     *_MD;
		size_t		blocking;

	public:
		CoppersmithSolver(const Domain &MD, size_t blocking_ = 0) :
			 _MD(&MD), blocking(blocking_)
		{}


		template <class Vector, class Blackbox>
		Vector &solveNonSingular (Vector &x, const Blackbox &B, const Vector &y) const
		{
			commentator().start ("Coppersmith solveNonSingular", "solveNonSingular");
#if 1
			std::ostream& report = commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
#endif

			//Set up the projection matrices and their dimensions
			size_t d = B.coldim();
			size_t r,c;
			integer tmp = uint64_t(d);

			//Set the blocking size, Using Pascal Giorgi's convention
			if(blocking==0){
				r=tmp.bitsize()-1;
				c=tmp.bitsize()-1;
			}else
				r=c=blocking;

			//Create the block
			Block U(field(),r,d);
			Block W(field(),d,c-1);
			Block V(field(),d,c);

			//Pick random entries for U and W. W will become the last c-1 columns of V

			U.random();
			W.random();

			//Multiply W by B on the left and place it in the last c-1 columns of V
			Sub V2(V,0,1,d,c-1);
			domain().mul(V2,B,W);

			//Make the first column of V a copy of the right side of the system, y
			for(size_t i=0; i<d; i++)
				V.setEntry(i,0,y[i]);

			//Create the sequence container and its iterator that will compute the projection
			BlackboxBlockContainer<Field, Blackbox > blockseq(&B,field(),U,V);

			//Get the generator of the projection using the Coppersmith algorithm (slightly modified by Yuhasz)
			BlockCoppersmithDomain<Domain, BlackboxBlockContainer<Field, Blackbox> > BCD(domain(), &blockseq,d);
			std::vector<Block> gen;
			std::vector<size_t> deg;
			deg = BCD.right_minpoly(gen);
			report << "Size of blockseq " << blockseq.size() << std::endl;
			report << "Size of gen " << gen.size() << std::endl;
			for(size_t i = 0; i < gen[0].coldim(); i++)
				report << "Column " << i << " has degree " << deg[i] << std::endl;

			//Reconstruct the solution
			//Pick a column of the generator with a nonzero element in the first row of the constant coefficient
			size_t idx = 0;
			if(field().isZero(gen[0].getEntry(0,0))){
				size_t i = 1;
				while(i<c && field().isZero(gen[0].getEntry(0,i)))
					i++;
				if(i==c)
					throw LinboxError(" block minpoly: matrix seems to be singular - abort");
				else
					idx=i;
			}

			//from 1 to the degree of the index column, multiply A^(i-1)V times the idx column of the generator coefficient x^i
			//Accumulate these results in xm
			size_t mu = deg[idx];
			Block BVo(V);
			Block BVe(field(),d,c);
			Block xm(field(),d,1);
			bool odd = true;
			for(size_t i = 1; i < mu+1; i++){
				Sub gencol(gen[i],0,idx,c,1); // BB changed d,1 to c,1
				Block BVgencol(field(),d,1);
				if(odd){
					domain().mul(BVgencol,BVo,gencol);
					domain().addin(xm, BVgencol);
					domain().mul(BVe,B,BVo);
					odd=false;
				}
				else{
					domain().mul(BVgencol,BVe,gencol);
					domain().addin(xm, BVgencol);
					domain().mul(BVo,B,BVe);
					odd=true;
				}

			}

			//For the constant coefficient, loop over the elements in the idx column except the first row
			//Multiply the corresponding column of W (the last c-1 columns of V before application of B) by the generator element
			//Accumulate the results in xm
			for(size_t i = 1; i < c; i++){
				Sub Wcol(W,0,i-1,d,1);
				Block Wcolgen0(field(),d,1);
				domain().mul(Wcolgen0, Wcol, gen[0].getEntry(i,idx));
				domain().addin(xm,Wcolgen0);
			}

			//Multiply xm by -1(move to the correct side of the equation) and divide the the 0,idx entry of the generator constant
			Element gen0inv;
			field().inv(gen0inv,gen[0].getEntry(0,idx));
			field().negin(gen0inv);
			domain().mulin(xm, gen0inv);

#if 0
			//Test to see if the answer works with U
			Block Bxm(field(),d,1), UBxm(field(),r,1), Uycol(field(), r,1);
			Sub ycol(V,0,0,d,1);
			domain().mul(Uycol, U, ycol);
			domain().mul(Bxm, B, xm);
			domain().mul(UBxm, U, Bxm);

			if(domain().areEqual(UBxm, Uycol))
				report << "The solution matches when projected by U" << endl;
			else
				report << "The solution does not match when projected by U" << endl;
#endif


			//Copy xm into x (Change type from 1 column matrix to Vector)
			for(size_t i =0; i<d; i++)
				x[i]=xm.getEntry(i,0);

			commentator().stop ("done", NULL, "solveNonSingular");
			return x;
		}




	}; // end of class CoppersmithSolver

	template <class _Domain>
	class CoppersmithRank{

	public:
		typedef _Domain 			Domain;
		typedef typename Domain::Field                    Field;
		typedef typename Domain::Element       Element;
		typedef typename Domain::Matrix 	Block;
		typedef typename Domain::Submatrix 	Sub;
		typedef typename Field::RandIter	Random;

		inline const Domain & domain() const { return *_MD; }
		inline const Field & field() const { return domain().field(); }
	protected:
		const Domain     *_MD;
		Random		iter;
		size_t		blocking;

		//Compute the determinant of a polynomial matrix at the given set of evaluation points
		//Store the results in the vector dets.
		void EvalPolyMat(std::vector<Element> &dets, std::vector<Element> &values, std::vector<Block> & mat) const {

			size_t deg = mat.size() -1;
			size_t numv = values.size();
			//Compute the determinant of the evaluation at values[i] for each i
			for(size_t i = 0; i<numv; i++){
				//copy the highest matrix coefficient
				Block evalmat(field(), mat[0].rowdim(), mat[0].coldim());
				domain().copy(evalmat,mat[deg]);
				//Evaluate using a horner style evaluation
				typename std::vector<Block>::reverse_iterator addit =  mat.rbegin();
				addit++;
				for(addit; addit != mat.rend(); addit++){
					domain().mulin(evalmat,values[i]);
					domain().addin(evalmat,*addit);
				}//end loop computing horner evaluation
				//Compute the determinant of the evaluation and store it in dets[i]
				dets[i] = det(dets[i],evalmat);
			}//end loop over evaluation points
		}//end evaluation of polynominal matrix determinant

	public:
		CoppersmithRank(const Domain &MD, size_t blocking_ = 0) :
			 _MD(&MD), blocking(blocking_), iter(MD.field())
		{}


		template <class Blackbox>
		size_t rank (const Blackbox &B) const
		{
			commentator().start ("Coppersmith rank", "rank");
#if 1
			std::ostream& report = commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
#endif

			//Set up the projection matrices and their dimensions
			size_t d = B.coldim();
			size_t r,c;
			integer tmp = uint64_t(d);

			//Set the blocking size, Using Pascal Giorgi's convention
			if(blocking==0){
				r=tmp.bitsize()-1;
				c=tmp.bitsize()-1;
			}else
				r=c=blocking;

			//Create the block
			Block U(field(),r,d);
			Block V(field(),d,c);

			//Pick random entries for U and W. W will become the last c-1 columns of V

			U.random();
			V.random();


			BlackboxBlockContainer<Field, Blackbox > blockseq(&B,field(),U,V);

			//Get the generator of the projection using the Coppersmith algorithm (slightly modified by Yuhasz)
			BlockCoppersmithDomain<Domain, BlackboxBlockContainer<Field, Blackbox> > BCD(domain(), &blockseq,d);
			std::vector<Block> gen;
			std::vector<size_t> deg;
			deg = BCD.right_minpoly(gen);
			for(size_t i = 0; i < gen[0].coldim(); i++)
				report << "Column " << i << " has degree " << deg[i] << std::endl;

			//Compute the rank via the determinant of the generator
			//Get the sum of column degrees
			//This is the degree of the determinant via Yuhasz thesis
			//size_t detdeg = std::accumulate(deg.begin(), deg.end(), 0);
			size_t detdeg= 0;
			for(size_t i = 0; i < gen[0].coldim(); i++)
				detdeg+=deg[i];
			//Set up interpolation with one more evaluation point than degree
			size_t numpoints = d+1;
			std::vector<Element> evalpoints(numpoints), evaldets(numpoints);
			for(typename std::vector<Element>::iterator evalit = evalpoints.begin(); evalit != evalpoints.end(); evalit++){
				do{
					//do iter.random(*evalit); while(field().isZero(*evalit));
					iter.random(*evalit);
				}while ((std::find(evalpoints.begin(), evalit, *evalit) != evalit));
			}//end evaluation point construction loop

			//Evaluate the generator determinant at the points
			EvalPolyMat(evaldets, evalpoints, gen);
			for(size_t k = 0; k <numpoints; k++)
				report << evalpoints[k] << "  " << evaldets[k] << std::endl;
			//Construct the polynomial using Givare interpolation
			//Stolen from Pascal Giorgi, linbox/examples/omp-block-rank.C
			typedef Givaro::Poly1CRT< Field >  PolyCRT;
			PolyCRT Interpolator(field(), evalpoints, "x");
			typename PolyCRT::Element Determinant;
			Interpolator.RnsToRing(Determinant,evaldets);
			Givaro::Degree intdetdeg;
			Interpolator.getpolydom().degree(intdetdeg,Determinant);
			Givaro::Degree intdetval;
			Interpolator.getpolydom().val(intdetval,Determinant);
			if(detdeg != (size_t) intdetdeg.value()){
				report << "sum of column degrees " << detdeg << std::endl;
				report << "interpolation degree " << intdetdeg.value() << std::endl;
			}
			report << "sum of column degrees " << detdeg << std::endl;
			report << "interpolation degree " << intdetdeg.value() << std::endl;
			report << "valence (trailing degree) " << intdetval.value() << std::endl;
			for(size_t k = 0; k<gen.size(); k++)
				domain().write(report, gen[k]) << "x^" << k << std::endl;
			Interpolator.write(report << "Interpolated determinant: ", Determinant) << std::endl;
			size_t myrank = size_t(intdetdeg.value() - intdetval.value());
			commentator().stop ("done", NULL, "Coppersmith rank");
			return myrank;
		}




	}; // end of class CoppersmithRank

	//Use the coppersmith block wiedemann to compute the determinant
	template <class _Domain>
	class CoppersmithDeterminant{

	public:
		typedef _Domain 			Domain;
		typedef typename Domain::Field                    Field;
		typedef typename Domain::Element       Element;
		typedef typename Domain::Matrix 	Block;
		typedef typename Domain::Submatrix 	Sub;
		typedef typename Field::RandIter	Random;

		inline const Domain & domain() const { return *_MD; }
		inline const Field & field() const { return domain().field(); }
	protected:
		const Domain     *_MD;
		Random		iter;
		size_t		blocking;

		//Compute the determinant of a polynomial matrix at the given set of evaluation points
		//Store the results in the vector dets.
		void EvalPolyMat(std::vector<Element> &dets, std::vector<Element> &values, std::vector<Block> & mat) const {

			size_t deg = mat.size() -1;
			size_t numv = values.size();
			//Compute the determinant of the evaluation at values[i] for each i
			for(size_t i = 0; i<numv; i++){
				//copy the highest matrix coefficient
				Block evalmat(mat[deg]);
				//Evaluate using a horner style evaluation
				typename std::vector<Block>::reverse_iterator addit =  mat.rbegin();
				addit++;
				for(addit; addit != mat.rend(); addit++){
					domain().mulin(evalmat,values[i]);
					domain().addin(evalmat,*addit);
				}//end loop computing horner evaluation
				//Compute the determinant of the evaluation and store it in dets[i]
				dets[i] = det(dets[i],evalmat);
			}//end loop over evaluation points
		}//end evaluation of polynominal matrix determinant

	public:
		CoppersmithDeterminant(const Domain &MD, size_t blocking_ = 0) :
			 _MD(&MD), blocking(blocking_), iter(MD.field())
		{}


		template <class Blackbox>
		Element det (const Blackbox &B) const
		{
			commentator().start ("Coppersmith determinant", "determinant");
#if 1
			std::ostream& report = commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
#endif

			//Set up the projection matrices and their dimensions
			size_t d = B.coldim();
			size_t r,c;
			integer tmp = uint64_t(d);

			//Use given blocking size, if not given use Pascal Giorgi's convention
			if(blocking==0){
				r=tmp.bitsize()-1;
				c=tmp.bitsize()-1;
			}else
				r=c=blocking;

			//Create the block
			Block U(field(),r,d);
			Block V(field(),d,c);

			//Pick random entries for U and W. W will become the last c-1 columns of V

			U.random();
			V.random();

			//Multiply V by B on the left
			domain().leftMulin(B,V);

			//Create the sequence container and its iterator that will compute the projection
			BlackboxBlockContainer<Field, Blackbox > blockseq(&B,field(),U,V);

			//Get the generator of the projection using the Coppersmith algorithm (slightly modified by Yuhasz)
			BlockCoppersmithDomain<Domain, BlackboxBlockContainer<Field, Blackbox> > BCD(domain(), &blockseq,d);
			std::vector<Block> gen;
			std::vector<size_t> deg;
			deg = BCD.right_minpoly(gen);

			//Compute the determinant via the constant coefficient of the determinant of the generator
			//Get the sum of column degrees
			//This is the degree of the determinant via Yuhasz thesis
			//size_t detdeg = std::accumulate(deg.begin(), deg.end(), 0);
#if 0
			size_t detdeg= 0;
			for(size_t i = 0; i < gen[0].coldim(); i++)
				detdeg+=deg[i];
#endif
			//Set up interpolation with one more evaluation point than degree
			size_t numpoints = 2*d;
			std::vector<Element> evalpoints(numpoints), evaldets(numpoints);
			for(typename std::vector<Element>::iterator evalit = evalpoints.begin(); evalit != evalpoints.end(); evalit++){
				do{
					do iter.random(*evalit); while(field().isZero(*evalit));
				}while ((std::find(evalpoints.begin(), evalit, *evalit) != evalit));
			}//end evaluation point construction loop

			//Evaluate the generator determinant at the points
			EvalPolyMat(evaldets, evalpoints, gen);
			//Construct the polynomial using Givare interpolation
			//Stolen from Pascal Giorgi, linbox/examples/omp-block-rank.C
			typedef Givaro::Poly1CRT<Field>  PolyCRT;
			PolyCRT Interpolator(field(), evalpoints, "x");
			typename PolyCRT::Element Determinant;
			Interpolator.RnsToRing(Determinant,evaldets);
			Givaro::Degree intdetdeg;
			Interpolator.getpolydom().degree(intdetdeg,Determinant);
			Givaro::Degree intdetval(0);
			Interpolator.getpolydom().val(intdetval,Determinant);
			if(d != (size_t)intdetdeg.value()){
				report << "The matrix is singular, determinant is zero" << std::endl;
				return field(0).zero;
			}
			Interpolator.write(report << "Interpolated determinant: ", Determinant) << std::endl;
			Element intdeterminant(field().zero);
			Interpolator.getpolydom().getEntry(intdeterminant,intdetval,Determinant);
			commentator().stop ("done", NULL, "Coppersmith determinant");
			return intdeterminant;
		}

	}; // end of class CoppersmithDeterminant


}// end of namespace LinBox

#endif //__LINBOX_coppersmith_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
