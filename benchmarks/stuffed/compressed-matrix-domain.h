#ifndef __compressed_matrix_domain_h__
#define __compressed_matrix_domain_h__

#include "compressed-field.h"
#include "compressed-matrix.h"
#include "compressed-unit.h"
#include "compressed-word.h"

// for now, we are just handling the cases of full + fully aligned matrices
// the constant of four russians is a compile-time constant

template <typename Field_>
class CompressedMatrixDomain {
	public: // typedefs
		using Field = Field_;
		using Matrix = CompressedMatrix<Field>;
		using Unit = typename Field::Unit;
		using Word = typename Field::Word;
		using Base = typename Field::Base;
		using Element = typename Field::Element;
	
	public: // data
		const Field &f;
	
	public: // constructors
		CompressedMatrixDomain (const Field &f_) : f(f_) {}
	
	public: // multiplication
		void mul_classical  (Matrix &C, const Matrix &A, const Matrix &B) const;
		void mul_compressed (Matrix &C, const Matrix &A, const Matrix &B) const;
		template <uint64_t NT, uint64_t K> void mul_four_russians (Matrix &C, const Matrix &A, const Matrix &B) const;
		void mul_four_russians (Matrix &C, const Matrix &A, const Matrix &B, uint64_t K) const;
	
	public: // binary ops
		void addin  (Matrix &C, const Matrix &A) const;
		void subin  (Matrix &C, const Matrix &A) const;
		void copy   (Matrix &C, const Matrix &A) const;
		void neg    (Matrix &C, const Matrix &A) const;
		void axpyin (Matrix &C, const Matrix &A, const Element a) const;
	
	public: // unary ops
		void negin     (Matrix &C) const;
		void smulin    (Matrix &C, const Element c) const;
		void zero      (Matrix &C) const;
		void normalize (Matrix &C) const;
	
	public: // debug
		void random (Matrix &C) const;
		void print  (Matrix &C, std::ostream &os = std::cerr) const; 
	
	protected: // helpers
		template <uint64_t S> void build_four_russians_table (Unit *dest, const Unit *src, uint64_t dw, uint64_t ss, uint64_t K) const;

};

#include "compressed-matrix-domain.inl"

#endif