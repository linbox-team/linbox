/* File: src/library/objects/vector/sparse_vectors.h
 * Author: Li Chen for the LinBox group
 */

#ifndef _SPARSE_Vector_
#define _SPARSE_Vector_


#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <cstdlib>


namespace LinBox 
{ 
  /* sparse_vector is extended from the vector<element, long> with nonconstant      size 
   */

  template<class Ring>
    class sparse_vector: public vector<pair<typename Ring::element,long> >
    {
      typedef typename Ring::element ring;
      Ring R;
      
    public:
      /* construct */
      sparse_vector(){}
      
      sparse_vector(long col_num){} //do nothing, this is just for the dense vector's case.
      
      /* destruct */
      ~sparse_vector(){}
      
      /* To get the memory size of this object */
      long size_byte(){
	return (R.size() + sizeof(long)) * size();
      }
      
      /* get the column number of the first non-zero entry */
      long firstcol(){ return front().second;}

       /* get the value of the first non-zero entry */
      ring firstvalue(){ return front().first;}
      
       /* read a row from a ifstream, it is the corresponding function to writeto */
      sparse_vector& readfrom(std::ifstream &from){
	pair<ring, long> pair;
	ring m;
	long n;
	
	from >> n;
	while(n != 0){
	  R.read(from, m);
	  pair.first = m;
	  pair.second = n;
	  push_back(pair);
	  from >> n;
	}
	return *this;
      }
      
      /* write a row to a ofstream, it is the corresponding function to readfrom*/
      void writeto(std::ofstream &to){
	if(empty()) return;  //nessassary???
	vector<pair<ring, long> >::iterator p1;
	for (p1 = begin(); p1 < end(); p1++){
	  to << (*p1).second << " ";
	  R.write(to, (*p1).first);
	  to << " ";
	}
	to << 0 << "\n";
      }

      /*  print the row  */
      void print(){
	vector<pair<ring, long> >::iterator p1;
	cout<<"[";
	for (p1 = begin(); p1 < end(); p1++){
	  cout <<"[";
	  cout << (*p1).first;
	  cout <<", "<<(*p1).second;
	  cout <<"]";
	}
	cout<<"]"<<endl;
      }

  
      /*  A.row_operation(B, a) means A = A - A_value_at(a cloumn) * B. because B_value_at(a column) is 1, so after row_operation, A_value_at(a column) = 0. */
      sparse_vector& row_operation(const sparse_vector& pivot_row, const long index){
	pair<ring, long> f;
	sparse_vector new_vec;
	vector<pair<ring,long> >::const_iterator p1, p2, p3, p4;
	
	ring C, A;
	p1 = begin();
	p2 = pivot_row.begin();
	
	p3 = end();
	p4 = pivot_row.end();
	
	A = (*(p1 + index)).first;
	R.neg(C, A);
	
	while(true){
	  if(p1 == p3){
	    for(p2; p2 < p4; p2++){
	      R.mul(f.first, C, (*p2).first);
	      f.second = (*p2).second;
	      new_vec.push_back(f);
	    }
	    *this = new_vec;
	    return *this;
	  }
	  if(p2 == p4){
	    for(p1; p1 < p3; p1++){
	      new_vec.push_back(*p1);
	    }
	    *this = new_vec;
	    return *this;
	  }
	  if((*p2).second > (*p1).second){
	    new_vec.push_back(*p1);
	    p1++;
	  }
	  else{
	    if((*p2).second < (*p1).second){
	      R.mul(f.first, C, (*p2).first);
	      f.second = (*p2).second;
	      new_vec.push_back(f);
	      p2++;
	    }
	    else{
	      R.sub(f.first, (*p1).first, R.mul(f.first, A ,(*p2).first));
	      if (!R.isZero(f.first)){
		f.second = (*p1).second;
		new_vec.push_back(f);
	      }
	      p1++;
	      p2++;
	    }
	  }
	}
      }
      
      /* get the nth non_zero entry's column number */
      long get_nth_col(long n){
	return (*(begin() + n)).second; 
      }
      
       /*  get the first non_zero entry and normalize the row  */
      void normalize(){
	vector<pair<ring,long> >::iterator p1;
	ring A, B;
	B = front().first;
	if(R.isOne(B)) return;
	
	R.inv(A, B);
	for(p1 = begin(); p1 < end(); p1++) R.mulin((*p1).first, A);
      }
      /* push back a element in the row */
      void Push_back(const pair<ring, long> x){
	push_back(x);
      }
      
      /* to judge if the row is empty */
      bool Empty(){
	return empty();
      }
      
      /* to change the index so that index is the next non_zero entry's column */
      void next_non_zero(long& index){ index++; }
      
    }; // class sparse_vector

} // namespace LinBox

#endif // _SPARSE_Vector_
