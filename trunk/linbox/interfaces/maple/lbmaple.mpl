module LinBox() 
	local create, destroy, lbXapply, lbXapplyT, lbXinitBB, lbXinitV,
	lbXrank, lbXdet, lbXmp, lbXssolve, lbXstart1, lbXstart2, 
	lbXend, lbXgetMatrix, lbXkillMatrix, lbXgetVector, lbXkillVector, lbXgetEntryVector, lbXgetEntryMatrix, ExToMap, MapToEx, libname, ObjKeyCount;
	export LinBoxBB, LBrank, LBdet, LBminpoly, LBApply, LBApplyTranspose, LinBoxV; # LBDiag, LBSolve;
	option package, load=create, unload=destroy;
	description "Maple interface to LinBox methods & functions";

	 create := proc()
		global `type/LinBoxBB`, `type/LinBoxV`;
		libname := "
/usb/seagrave/linbox/linbox-new/intefaces/maple/
liblbmaple.so";	
		lbXend := define_external('End', MAPLE,LIB=libname);
		lbXapply := define_external('apply', MAPLE,LIB=libname);
		lbXapplyT := define_external('applyT', MAPLE,LIB=libname);
		lbXinitV := define_external('initV', MAPLE,LIB=libname);
#		lbXdiag := define_external('diag', MAPLE,LIB=libname);
		lbXinitBB := define_external('initBB', MAPLE,LIB=libname);
		lbXrank := define_external('rank',MAPLE,LIB=libname);
		lbXdet := define_external('det', MAPLE,LIB=libname);
		lbXmp := define_external('minpoly', MAPLE,LIB=libname);
#		lbXssolve := define_external('ssolve', MAPLE,LIB=libname);
		lbXgetMatrix := define_external('getMatrix', MAPLE,LIB=libname);
		lbXgetVector := define_external('getVector',MAPLE,LIB=libname);
		lbXkillMatrix := define_external('killMatrix', MAPLE,LIB=libname);
		lbXkillVector := define_external('killVector', MAPLE,LIB=libname);
# 		lbXgetEntryVector := define_external('getEntryVector',MAPLE,LIB=libname);
#		lbXgetEntryMatrix := define_external('getEntryMatrix',MAPLE,LIB=libname);
		`type/LinBoxBB` :=  `module`(getMatrix) ;
		`type/LinBoxV` := `module`(getVector) ;
		ObjKeyCount := 0;
	end proc;

	create();

	destroy := proc()
		lbXend();
	end proc;

	ExToMap := proc()
		local b0, s, r, A, i;
		if type(args[1], integer) then
			args[1];
		else
			A := args[1];
			s := A[nops(A)];
			b0 := 2^kernelopts(wordsize);
			for i from nops(A) - 1 to 1 do
				s := s * b0 + A[i];
			od;
			s
		fi;
	end proc;

	MapToEx := proc(num::posint)
		local r;
		if num < kernelopts(maximmediate) then 
			num;
		else 
			r := [disassemble(addressof(num))];
			subsop(1=nops(r),r);
		fi;
	end proc;  

	LinBoxBB := proc( )
		local k, f, D, R, C, rows, cols, s, nnz, M, p; #, check;

		if type(args[1], posint) then
			if nargs = 1 then
				k := args[1];
	#		elif type(args[2], posint) and type(args[3], posint) and type(args[4], proc) then
	#			rows := args[1];
	#			cols := args[2];
	#			k := rand(kernelopts(maximmedite));
	#			check := x -> evalb(rhs(x) <> 0);
	#			s := select(check, {seq(seq( (i,j),apply(args[3],i,j)), i=1..rows),j=1..cols)} );

			elif type(args[2], Matrix) then
				ObjKeyCount := ObjKeyCount + 1;
				k := ObjKeyCount;
				M := args[2]; p := args[1];
				# if p = 0 then
				# ... insert special code for integer ring
				# computation, as opposed to finite field computation
				k := addressof(M);
				if op(3,M)[1] = 'datatype=integer[4]' then
					lbXinitBB(1, k, p, M);
				elif p < kernelopts(maximmediate) then
					s := op(2,M);
					nnz := nops(s);
					rows := op(1,M)[1]; cols := op(1,M)[2];
					D := [seq(rhs(s[i]),i=1..nnz)];
					R := [seq(lhs(s[i])[1],i=1..nnz)];
					C := [seq(lhs(s[i])[2],i=1..nnz)];
					lbXinitBB(2, k, p, D, R, C, rows, cols, nnz);
				else
					s := op(2,M);
					nnz := nops(s);
					rows := op(1,M)[1]; cols := op(1,M)[2];
					D := [seq(MapToEx(rhs(s[i])),i=1..nnz)];
					R := [seq(lhs(s[i])[1],i=1..nnz)];
						C := [seq(lhs(s[i])[2],i=1..nnz)];
					lbXinitBB(3, k, MapToEx(p), D, R, C, rows, cols, nnz);
				fi;
			else
				error "Wrong number or type of parameters!";
			fi;
		else
			error "Wrong number or type of parameters!";
		fi;

		module()
			local destruct;
			export getMatrix, key; # getEntry
			option unload = destruct;
			description "Maple Container for LinBox blackbox";

			key := k;

			getMatrix := proc()
	# Implents difference between Maple 7 & 6(8)
	# Many thanks to MapleSoft for pointing out the need for a workaround
	# for Maple v7
	# If you are running v7, run the "exceptional" function
	# otherwise, run the normal function
				if `Maple 7` = substring(kernelopts(version),1..7) then lbXgetMatrix(k, ExToMap , 1);
				else lbXgetMatrix(k, ExToMap , 2);
				fi;
			end proc;
	
#			getEntry := proc(i::posint, j::posint)
#				lbXgetEntryMatrix(key, i, j);
#			end proc;

			destruct := proc()
				lbXkillMatrix(k):
			end proc;
				
		end module

	end proc;

	LinBoxV := proc()
		local k, index, L, mode, V, p, t, i;
		if type(args[1], posint) then # two possibilities:  Vector from
					      # a existing vector (key), or
					      # build one from a proc
			if nargs = 1 then
				k := args[1];
			elif type(args[2], 'procedure') then


#	Note - There is a slight danger, as this method uses random
#	numbers.  On occaison, one of these numbers might match 
#	a number used a different structure, in which case you will have a
#	copy of a previous blackbox returned to you, which is bad.
#	Looking for a better way to do this.

				ObjKeyCount := ObjKeyCount + 1;
				k := ObjKeyCount;
				mode := 1;
				p := proc(n)
					local r;
					r := MapToEx(apply(args[2], n));
					if type(r, list) then mode := 2;
					r
				end proc; 
				


				L := [seq(p(t), t=1..args[1])];
				lbXinitV(mode, k(), L, args[1]);
			else
				error "Wrong number or type of parameters!";			
			end if;
		elif type(args[1], Vector) then # No key, we're translating from a vector
			ObjKeyCount := ObjKeyCount + 1;
			k := ObjKeyCount;
			V := args[1];
			p := proc(n)
				local r;
				r := MapToEx(rhs(n));
				if type(r,list) then mode := 4; fi;
				r
			end proc;
			
			mode := 3;
			t := sort(op(table(op(2,V))), (x,y)->( lhs(x) < lhs(y) ) );
			index := map(lhs, t);
			L := map(p, t);
			lbXinitV(mode, k, index, L, nops(index) );

		elif type( args[1], list) then 
			ObjKeyCount := ObjKeyCount + 1;

			k := ObjKeyCount;
			lbXinitV(1,k(), args[1], nops(args[1])); 


		else
			error "Wrong number or type of parameters!";
		fi;

		module()
			local destruct;
			export getVector, key; #getEntry
			description "Maple container for LinBox Vector";
			option unload = destruct;
				
			key := k;

			destruct := proc()
				lbXkillVector(key);
			end proc;

#			getEntry := proc(i::posint)
#				lbXgetEntryVector(key, i);
#			end proc;

			getVector := proc()
	# Once again, split b/c of difference between Maple 7 and 6(8)
	# For maple 7, run the special getMatrix function
	# otherwise, run the normal function
				if `Maple 7` = substring(kernelopts(version),1..7) then lbXgetVector(key, ExToMap, 1);
				else lbXgetVector(key, ExToMap, 2);
				fi;
			end proc;
			
		end module
	end proc;


	LBApply := proc()
		local M, V, R;
		if nargs <> 2 then
			error "Wrong number or type of parameters!";
		elif not type(args[1], Matrix) and not type(args[1],LinBoxBB) then
			error "Wrong number or type of parameters!";
		elif not type(args[2],Vector) and not type(args[2], LinBoxV) then
			error "Wrong number or type of parameters!";
		fi;

		# if Matrix is not LinBox BB, make it so
		if type(args[1], Matrix) then 
			M := LinBoxBB(args[1]);
		else
			M := args[1];
		fi;

		# if Vector is not LinBoxV, make it so
		if type(args[2], Vector) then
			V := LinBoxV(args[2]);
		else
			V := args[2];
		fi;
	
		ObjKeyCount := ObjKeyCount + 1;

		lbXapply(M:-key, V:-key, ObjKeyCount );

		if type(args[2], Vector) then # if the user sends the vector
					      # as a maple type, return
					      # as that type
			R := LinBoxV(ObjKeyCount);
			R:-getVector()
		else
			LinBoxV(ObjKeyCount)
		end if;
	end proc;

	LBApplyTranspose := proc()
		local M, V, R;
		if nargs <> 2 then
			error "Wrong number or type of parameters!";
		elif not type(args[1], Matrix) and not type(args[1],LinBoxBB) then
			error "Wrong number or type of parameters!";
		elif not type(args[2],Vector) and not type(args[2], LinBoxV) then
			error "Wrong number or type of parameters!";
		fi;

		# if Matrix is not LinBox BB, make it so
		if type(args[1], Matrix) then 
			M := LinBoxBB(args[1]);
		else
			M := args[1];
		fi;

		# if Vector is not LinBoxV, make it so
		if type(args[2], Vector) then
			V := LinBoxV(args[2]);
		else
			V := args[2];
		fi;
		
		ObjKeyCount := ObjKeyCount + 1;
		lbXapplyT(M:-key, V:-key, ObjKeyCount );

		if type(args[2], Vector) then # if the user sends the vector
					      # as a maple type, return
					      # as that type
			R := LinBoxV(ObjKeyCount);
			R:-getVector();
		else
			LinBoxV(ObjKeyCount);
		fi;
	end proc;

	LBrank := proc()
		local M;
		if type(args[1], LinBoxBB) then
			M := args[1];
		elif type(args[1], Matrix) and type(args[2], posint) then
			M := LinBoxBB(args[1], args[2]);
		else
			error "Wrong number or type of parameters!";
		fi;

		lbXrank(M:-key);
	end proc;

	LBdet := proc()
		local M;
		if type(args[1], LinBoxBB) then
			M := args[1];
		elif type(args[1], Matrix) and type(args[2], posint) then
			M := LinBoxBB(args[1], args[2]);
		else
			error "Wrong number or type of parameters!";
		fi;
		
		lbXdet(M:-key);
	end proc;

	# minpoly function.  If no variable name is given, then a list
	# containg the coefficients in sorted order, from degree lowest to
	# highest is returned (note this is a dense coefficient vector, there
	# are zeros here

	LBminpoly := proc()
		local L, i, x, M;
		if type(args[1], LinBoxBB) then
			M := args[1];
			if nargs = 1 then 
				x := true;
			elif nargs = 2 and type(args[2],name) then
				x := args[2];
			else
				error "Wrong number or type of parameters!";
			fi;
		elif type(args[1],Matrix) and type(args[2], posint) then
			M := LinBoxBB(args[1], args[2]);
			if nargs = 2 then 
				x := true;
			elif nargs = 3 and type(args[3], name) then
				x := args[3];
			else
				error "Wrong number or type of parameters!";
			fi;
		else
			error "Wrong number or type of parameters!";
		fi;
		
		L := lbXmp(M:-key);
		L := [seq(ExToMap(L[i]),i=1..nops(L))];
		if type(x, boolean) then
			L
		else
			add(ExToMap(L[i])*x^(i-1),i=1..nops(L))
		end if;
	end proc;

end module;			
