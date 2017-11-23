with(LinBox);
with(LinearAlgebra):


### CHECKING SOLUTION OVER THE INTEGER

checking_linbox:=proc(n,prec)
local gen,A,b, allpassed, max;
        max:=2^prec;
	gen:=rand(-max..max);
	A:=Matrix(n,gen):
	b:=Vector(n,gen):

allpassed:=true;
if (not Equal(lbSolve(A,b),A^(-1).b)) then
   print("linear system solving failed\n");
   lbSolve(A,B); A^(-1).b;
   allpassed:=false;
end if;
if not(lbDeterminant(A)=Determinant(A)) then
   print("Determinant failed");
   lbDeterminant(A);Determinant(A); 
   allpassed:=false;
end if;
if not (lbRank(A)=Rank(A)) then
   print("Rank failed");
   lbRank(A); Rank(A); 
   allpassed:=false;
end if;
if not(lbMinpoly(A,'x')=MinimalPolynomial(A,'x')) then
   print("Minpoly failed");
   lbMinpoly(A,'x');MinimalPolynomial(A,'x'); 
   allpassed:=false;
end if;
if not (lbCharpoly(A,'x')=CharacteristicPolynomial(A,'x')) then
   print("Minpoly failed");
   lbCharpoly(A,'x');CharacteristicPolynomial(A,'x');
   allpassed:=false;	
end if;
if (not Equal(lbMul(A,A), A.A)) then 
   print("Matrix multiplication failed");	    
   lbMul(A,A);A.A; 
   allpassed:=false;
end if;

if (allpassed) then print("ALL DONE FOLKS") end if;
end proc:



checking_linbox(10,10);
