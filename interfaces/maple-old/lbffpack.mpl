FFPACK:=module()
	export LBMultiply,LBAddMultiply,LSPFactor,LQUPFactor,LBRank,LBDeterminant,LBInverse;
	local fgemm,lsp,lbrank,lbdeterminant,lbinverse,create,libname;
	option package, load=create, unload=gc;
	description "Maple interface to FFPACK package";


	create := proc()
	libname :="/home/pgiorgi/Library/src/linbox/interfaces/maple/.libs/liblbmapleffpack.so"; 
	fgemm:=define_external('fgemm',MAPLE,LIB=libname);
	lsp:=define_external('lsp',MAPLE,LIB=libname);
	lbrank:=define_external('rank',MAPLE,LIB=libname);
	lbdeterminant:=define_external('determinant',MAPLE,LIB=libname);
	lbinverse:=define_external('inverse',MAPLE,LIB=libname);
	
	end proc;

	create();


	LBMultiply:=proc()
		# calls are :
		#		LBMultiply(p,A,B) 		-> return A*B mod p
		#		LBMultiply(p,alpha,A,B)		-> return alpha.A*B mod p
		#		LBMultiply(p,A,B,C)		-> return C:= A*B mod p
		#		LBMultiply(p,alpha,A,B,C)	-> return C:= alpha.A*B mod p

		local m,n,k1,k2,alpha,C,a,b;
	
		if (type(args[2],Matrix)) then
			alpha:=1;
			m,k1 :=  LinearAlgebra:-Dimension(args[2]);
			k2,n :=  LinearAlgebra:-Dimension(args[3]);
			a:=2;b:=3;

		elif (type(args[2],numeric)) then
			alpha:=args[2];
			m,k1 :=  LinearAlgebra:-Dimension(args[3]);
			k2,n :=  LinearAlgebra:-Dimension(args[4]);
			a:=3:b:=4;
		else
			error " 2nd argument should be either a Matrix or a numeric value";

		end if;

		if (k1<> k2) then
			error "Dimensions don't match",k1,k2;
		end if;
	
	
		if (nargs <b+1) then
			C:=LinearAlgebra:-Modular:-Create(args[1],m,n,float[8]);
			fgemm(args[1],m,n,k1,alpha,args[a],args[b],0,C);
		else	
			fgemm(args[1],m,n,k1,alpha,args[a],args[b],0,args[b+1]);
		end if;
				
	end proc;


	LBAddMultiply:=proc()
		# calls are:
		#		 LBAddMultiply(p,alpha,A,B,beta,C) -> return C:= beta.C+alpha.A*B mod p

		local m,n,k1,k2,m1,n1;
	
		m,k1:=  LinearAlgebra:-Dimension(args[3]);
		k2,n:=  LinearAlgebra:-Dimension(args[4]);
		m1,n1:= LinearAlgebra:-Dimension(args[6]);


		
		if ((k1<>k2) or ( m1<> m) or (n1 <> n)) then
			error "Dimensions don't match";
		end if;

		fgemm(args[1],m,n,k1,args[2],args[3],args[4],args[5],args[6]);

				
	end proc;


	LSPFactor:=proc()
		# calls are:
		#		LSPFactor(p,A)		-> return A,P ;   A:=[LSP] mod p (L is lower part of A , S is upper part of A)
		#		LSPFactor(p,A,perm)	-> return A,perm; A:=[LSP],perm:=P

		local m,n,P;
		m,n:= LinearAlgebra:-Dimension(args[2]);
			
		if (nargs <3) then			
			P:=rtable(1..n,subtype=Vector[column],datatype=integer[4],order=C_order);
			#P:=Vector[column](1..n);
			lsp(args[1],m,n,args[2],P,1);
			return args[2],P;
		else
			lsp(args[1],m,n,args[2],args[3],1);
			return args[2],args[3];	
		end if;				

	end proc;


	LQUPFactor:=proc()
		# calls are:
		#		LQUPFactor(p,A)		-> return A,Q,P ;   A:=[LQUP] mod p (L is lower part of A  in compress form, U is upper part of A)
		#		LQUPFactor(p,A,perm1,perm2)	-> return A,perm1,perm2; A:=[LSP],perm:=P

		local m,n,P,Q;
		m,n:= LinearAlgebra:-Dimension(args[2]);
			
		if (nargs <3) then
			Q:=rtable(1..m,subtype=Vector[column],datatype=integer[4],order=C_order);
			P:=rtable(1..n,subtype=Vector[column],datatype=integer[4],order=C_order);
			lsp(args[1],m,n,args[2],P,2,Q);
			return args[2],P,Q;
		else
			lsp(args[1],m,n,args[2],args[3],2,args[4]);	
			return args[2],args[3],args[4];		
		end if;				

	end proc;

	LBRank:=proc()
		# calls are:
		#		LBRank(p,A)	return rank(A mod p)

		local m,n;
		m,n:= LinearAlgebra:-Dimension(args[2]);
		lbrank(args[1],m,n,args[2]);

	end proc;

	LBDeterminant:=proc()
		# calls are:
		#		LBDeterminant(p,A)	return determinant(A) mod p
				local m,n;
		m,n:= LinearAlgebra:-Dimension(args[2]);
		lbdeterminant(args[1],m,n,args[2]);

	end proc;

	
	LBInverse:=proc()
		# calls are:
		#		LBInverse(p,A)	return A^(-1) mod p
		#		LBInverse(p,A,X) return X:=A^(-1) mod p

		local m,n,X;
		m,n:= LinearAlgebra:-Dimension(args[2]);
		if (m<>n) then error "matrix should be square"; end if;

		if (nargs <3) then
			X:=LinearAlgebra:-Modular:-Create(args[1],m,m,float[8]);
			lbinverse(args[1],m,args[2],X);
		else
			lbinverse(args[1],m,args[2],args[3]);
		end if;

	end proc;

	


end module;
