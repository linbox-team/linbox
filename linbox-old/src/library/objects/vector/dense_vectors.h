/* File: src/library/objects/vector/sparse_vectors.h
 * Author: Li Chen for the LinBox group
 */

#ifndef _DENSE_Vector_
#define _DENSE_Vector_


#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <cstdlib>


namespace LinBox 
{ 
  /* dense_vector is extended from the vector<element> with constant size
   */
  
  template<class Ring>
    class dense_vector: public vector<typename Ring::element>
    {
      typedef typename Ring::element ring;
      Ring R;
      
    public:
      
      /* construct */
      dense_vector(){}
      
      dense_vector(long col_num){
	ring zero;
	zero = 0;
	resize(col_num, zero);
      }

      /* destruct */
      ~dense_vector(){}
      
      /* To get the memory size of this object */
      long size_byte(){
	return R.size() * size();
      }
      
      /* get the column number of the first non-zero entry */
      long firstcol(){
	vector<ring>::iterator p1;
	for(p1 = begin(); p1 < end(); p1++){
	  if(!R.isZero(*p1)) return p1 - begin();
	}
	return -1; //all entry are 0
      }
      
      /* get the value of the first non-zero entry */
      ring firstvalue(){
	vector<ring>::iterator p1;
	for(p1 = begin(); p1 < end(); p1++){
	  if(!R.isZero(*p1)) return *p1;
	}
	ring zero;
	R.init(zero, 0);
	return zero; ///all entry are 0
      }
      
      /* read a row from a ifstream, it is the corresponding function to writeto */
      dense_vector& readfrom(std::ifstream &from){
	ring m;
	long n;
	
	from >> n; //read the size of the row
	for(int i = 0; i < n; i++){
	  R.read(from, m);
	  push_back(m);
	}
	return *this;
      }
       /* write a row to a ofstream, it is the corresponding function to readfrom*/
      void writeto(std::ofstream &to){
	if(empty()) return;  //nessassary???
	vector<ring>::iterator p1;
	to << size() << " ";  // put the size of the row
	for (p1 = begin(); p1 < end(); p1++){
	  R.write(to, (*p1));
	  to << " ";
	}
	to << endl;
      }
      
     /*  print the row  */
      void print(){
	vector<ring>::iterator p1;
	cout<<"[";
	for (p1 = begin(); p1 < end(); p1++){
	  cout << (*p1) << ", ";
	}
	cout<<"]"<<endl;
      }
      
      
      /*  A.row_operation(B, a) means A = A - A_value_at(a cloumn) * B. because B_value_at(a column) is 1, so after row_operation, A_value_at(a column) = 0. */
      dense_vector& row_operation(dense_vector& pivot_row, long& index){
	vector<ring>::iterator p1, p2;
	if(size() != pivot_row.size()){
	  cout << " the size of the rows of row_operation are diffirent. ERROR!" << endl;
	  exit(1);
	}
	
	ring C, A;
	p1 = begin();
	p2 = pivot_row.begin();
	A = (*(p1 + index));
	if(!R.isZero(A)){
	  R.neg(C, A);
	  ring temp;
	  for(p1; p1 < end(); p1++){
	    R.addin((*p1), R.mul(A, C, (*p2)));
	    p2++;
	  }
	}
	next_non_zero(index); //increase the index to next non_zero entry
	return *this;
      }
      
      /* useless function? but nessesary in sparse_vector */
      long get_nth_col(long n){
	return n; 
      }
      
     /*  get the first non_zero entry and normalize the row  */
      void normalize(){
	vector<ring>::iterator p1;
	ring A, B;
	B = firstvalue();
	if(R.isZero(B) || R.isOne(B)) return;
	R.inv(A, B);
	for(p1 = begin(); p1 < end(); p1++) R.mulin((*p1), A);
      }
      
      /* push back a element in the row */
      void Push_back(const pair<ring, long> x){
	vector<ring>::iterator p1;
	p1 = begin() + x.second - 1;
	(*p1) = x.first;
      }
      
      /* to judge if the row is empty */
      bool Empty(){
	for(vector<ring>::iterator p1 = begin(); p1 < end(); p1++){
	  if(!R.isZero(*p1)) return false;
	}
	return true;
      }
      
      /* to change the index so that index is the next non_zero entry's column */
      void next_non_zero(long& index){
	for(vector<ring>::iterator p1 = begin() + index + 1; p1 < end(); p1++){
	  if(!R.isZero(*p1)){
	    index = p1 - begin();
	    return;
	  }
	}
	index = size();
      }
      
      
    }; // class dense_vector
  
} // namespace LinBox

#endif // _DENSE_Vector_
