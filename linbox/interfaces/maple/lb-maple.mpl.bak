LinBox:= module()
        description "Maple interface for LinBox library version 1.0 by Pascal Giorgi (LinBox project - www.linal.org)";
	local lbpath, lbInit, lbEnd, lbStart, lbStop:
	export lbDeterminant, lbRank, lbMinpoly, lbCharpoly, lbSolve, Dev:
	option package, load=lbInit, unload=lbEnd: 
	
	lbpathvalue

	############################################################# 
	# Function to initialize LinBox driver and LinBox interface #
	#############################################################	
	lbStart  := define_external('lbStart' , MAPLE, LIB=lbpath);
	lbStop   := define_external('lbStop'  , MAPLE, LIB=lbpath);


	lbInit := proc() printf("  LinBox/Maple Interface Package beta version 1.0 \n  by Pascal Giorgi (pascal.giorgi@univ-perp.fr)\n");			 
			 lbStart();
		  end proc;
	lbEnd  := proc() end proc;#lbStop();end proc;  


#####################################################################
#############  Definition of the low level submodule ################
#####################################################################

	Dev:= module()
	      description "low level interface for direct LinBox object manipulation";
	        local lbInit, lbEnd;
		export lbDataInfo,lbStart, lbStop,
	        lbElement, lbDomain, lbBlackbox, lbVector, 
	      	lbSetPrimeField, lbSetRationalField, lbSetIntegerRing,
		lbCopyBlackbox, lbBlackboxDimension, lbSetBlackboxAtRandom, lbRebindBlackbox, lbWriteBlackbox, lbSetBlackbox,
		lbCopyVector, lbVectorDimension, lbSetVectorAtRandom, lbRebindVector, lbWriteVector, lbSetVector,
		lbWritePolynomial,
		lbDeterminant, lbRank, lbMinpoly, lbCharpoly, lbSolve,
		lbConvertElement, lbConvertBlackbox, lbConvertVector, lbConvertPolynomial,
		lbCopy, lbWrite, lbDimension, lbRebind, lbRandom:       
		option package, load=lbInit, unload=lbEnd:							

		
        	############################################################# 
		# Function to initialize LinBox driver and LinBox interface #
		#############################################################	
		lbStart  := define_external('lbStart' , MAPLE, LIB=lbpath);
		lbStop   := define_external('lbStop'  , MAPLE, LIB=lbpath);

		######################################################
		# Function to get information from the LinBox driver #
		######################################################
		lbDataInfo := define_external('lbDataInfo', MAPLE, LIB=lbpath);	

		
		########################################
		# Function to set the LinBox data type #
		########################################	
		lbSetPrimeField    := define_external('lbSetPrimeField'     , MAPLE, LIB=lbpath);
		lbSetRationalField := define_external('lbSetRationalField'  , MAPLE, LIB=lbpath);
		lbSetIntegerRing   := define_external('lbSetIntegerRing'    , MAPLE, LIB=lbpath);
		lbSetBlackbox      := define_external('lbSetBlackbox'       , MAPLE, LIB=lbpath);
		lbSetVector        := define_external('lbSetVector'         , MAPLE, LIB=lbpath);


		###############################################
		# Function to create LinBox object from Maple #
		###############################################	
		lbElement   := define_external('lbCreateElement'   , MAPLE, LIB=lbpath);
		lbDomain    := define_external('lbCreateDomain'    , MAPLE, LIB=lbpath);
		lbBlackbox  := define_external('lbCreateBlackbox'  , MAPLE, LIB=lbpath);
		lbVector    := define_external('lbCreateVector'    , MAPLE, LIB=lbpath);


		#######################################
		# Function for LinBox blackbox object #
		#######################################		
		lbCopyBlackbox        := define_external('lbCopyBlackbox'         , MAPLE, LIB=lbpath);
		lbBlackboxDimension   := define_external('lbGetBlackboxDimension' , MAPLE, LIB=lbpath);
		lbSetBlackboxAtRandom := define_external('lbSetBlackboxAtRandom'  , MAPLE, LIB=lbpath);
		lbRebindBlackbox      := define_external('lbRebindBlackbox'       , MAPLE, LIB=lbpath);
		lbWriteBlackbox       := define_external('lbWriteBlackbox'        , MAPLE, LIB=lbpath);
		

		#####################################
		# Function for LinBox vector object #
		#####################################
		lbCopyVector        := define_external('lbCopyVector'         , MAPLE, LIB=lbpath);
		lbVectorDimension   := define_external('lbGetVectorDimension' , MAPLE, LIB=lbpath);
		lbSetVectorAtRandom := define_external('lbSetVectorAtRandom'  , MAPLE, LIB=lbpath);
		lbRebindVector      := define_external('lbRebindVector'       , MAPLE, LIB=lbpath);
		lbWriteVector       := define_external('lbWriteVector'        , MAPLE, LIB=lbpath);
		
		#########################################
		# Function for LinBox polynomial object #
		#########################################
		lbWritePolynomial   := define_external('lbWritePolynomial'    , MAPLE, LIB=lbpath);
	
		##############################
		# Available LinBox solutions #
		##############################
		lbDeterminant := define_external('lbDeterminant'  , MAPLE, LIB=lbpath);
		lbRank        := define_external('lbRank'         , MAPLE, LIB=lbpath);	
		lbMinpoly     := define_external('lbMinpoly'      , MAPLE, LIB=lbpath);
		lbCharpoly    := define_external('lbCharpoly'     , MAPLE, LIB=lbpath);
		lbSolve       := define_external('lbSolve'        , MAPLE, LIB=lbpath);	


		##################################################
		# Conversion from LinBox object to Maple Objects #
		##################################################
		lbConvertElement    := define_external('lbConvertElement'    , MAPLE, LIB=lbpath);
		lbConvertBlackbox   := define_external('lbConvertBlackbox'   , MAPLE, LIB=lbpath);
		lbConvertVector     := define_external('lbConvertVector'     , MAPLE, LIB=lbpath);
		lbConvertPolynomial := define_external('lbConvertPolynomial' , MAPLE, LIB=lbpath);

		####################
		# Higher level API #
		####################			
		lbCopy      := define_external('lbCopy'       , MAPLE, LIB=lbpath);
		lbWrite     := define_external('lbWrite'      , MAPLE, LIB=lbpath);
		lbDimension := define_external('lbDimension'  , MAPLE, LIB=lbpath);
		lbRebind    := define_external('lbRebind'     , MAPLE, LIB=lbpath);
		lbRandom    := define_external('lbRandom'     , MAPLE, LIB=lbpath);	
	
	end module:

##############################################################
#############  end of the low level submodule ################
##############################################################


###################################################################
#############  definition of the high level module ################
###################################################################
		

		################################################
		# Determinant computation from a Maple Matrix  #
		################################################		
		lbDeterminant := proc()
		      local d, A, det, t;			  
		      if (type(args[1], integer)) then	      	      
		      	 if (type(args[2], Matrix)) then	      
			    t := time();
			    A  := Dev:-lbBlackbox(args[1], args[2]);			   					 
			    t:= time()-t;
			 else
			    error("invalid argument",args[2], "must be a Matrix");
			 end if;
		      elif (type(args[1], Matrix)) then
			    t := time();
			    A := Dev:-lbBlackbox(args[1]);			
			    t:= time()-t;					
		      else
			    error("invalid argument",args[1], "must be a Matrix");
		      end if;		 		   
		      d   := Dev:-lbDeterminant(A);		     		     
		      det := Dev:-lbConvertElement(d); 
		      if (printlevel>1) then
			      print("conversion time:",t);
		      end if;
		      return det;			      
		end;

		########################################
		# Rank computation from a Maple matrix #
		########################################		
		lbRank := proc()
		      local A,t;
		      if (type(args[1], integer)) then	      	      
		      	 if (type(args[2], Matrix)) then	      
			    t := time();
			    A:=Dev:-lbBlackbox(args[1], args[2]); 
			    t:=time()-t;
			 else
			    error("invalid argument",args[2], "must be a Matrix");
			 end if;
		      elif (type(args[1], Matrix)) then
			   t := time();
			   A := Dev:-lbBlackbox(args[1]);
			   t:=time()-t;			
		      else
			    error("invalid argument",args[1], "must be a Matrix");
		      end if;
		      if (printlevel>1) then
			      print("conversion time:",t);
		      end if;
		      return Dev:-lbRank(A);			      
		end;

		######################################################
		# Minimal polynomial computation from a Maple matrix #
		######################################################		
		lbMinpoly := proc()
		      local A, p, symb,t;
		      if (type(args[1], integer)) then	      	      
		      	 if (type(args[2], Matrix)) then
			    if (type(args[3], name)) then
			        t:=time(); 	      
			        A:=Dev:-lbBlackbox(args[1], args[2]); 
				t:=time()-t;
				symb:=args[3];
			    else
				 error("invalid argument",args[3], "must be a name");
		            end if;
			 else
			     error("invalid argument",args[2], "must be a Matrix");
			 end if;  
		      elif (type(args[1], Matrix)) then
			   if (type(args[2], name)) then 
			      t:=time(); 
			      A := Dev:-lbBlackbox(args[1]);
			      t:=time()-t;
			      symb:=args[2];			
		           else
			      error("invalid argument",args[2], "must be a name");
			   end if
		      else
			error("invalid argument",args[1], "must be a Matrix");
		      end if;
		      p:= Dev:-lbMinpoly(A); 
		      if (printlevel>1) then
			      print("conversion time:",t);
		      end if;
		      return Dev:-lbConvertPolynomial(p,symb);	 
		end;

		#############################################################
		# Characteristic polynomial computation from a Maple matrix #
		#############################################################		
		lbCharpoly := proc()
		      local A, p, symb,t;		      		      
		      if (type(args[1], integer)) then	      	      
		      	 if (type(args[2], Matrix)) then
			    if (type(args[3], name)) then  
			        t:=time(); 	      
			        A:=Dev:-lbBlackbox(args[1], args[2]); 
				t:=time()-t;
				symb:=args[3];
			    else
				 error("invalid argument",args[3], "must be a name");
		            end if;
			 else
			     error("invalid argument",args[2], "must be a Matrix");
			 end if;  
		      elif (type(args[1], Matrix)) then
			   if (type(args[2], name)) then 
			         t:=time(); 	
				 A := Dev:-lbBlackbox(args[1]);
  			 	 t:=time()-t;
			      symb:=args[2];			
		           else
			      error("invalid argument",args[2], "must be a name");
			   end if
		      else
			error("invalid argument",args[1], "must be a Matrix");
		      end if;
		      p:= Dev:-lbCharpoly(A); 
		      if (printlevel>1) then
			      print("conversion time:",t);
		      end if;
		      return Dev:-lbConvertPolynomial(p,symb);			      
		end;

		#############################################################
		# Linear system solving from Maple matrix and Maple Vector  #
		#############################################################		
		lbSolve := proc()
		      local A,b,x,y,t;
		      if (type(args[1], integer)) then	     	      
			 if (type(args[2], Matrix)) then	      
			    if (type(args[3], Vector)) then 
			        t:=time(); 	
			        A:=Dev:-lbBlackbox(args[1], args[2]);
				b:=Dev:-lbVector(args[1], args[3]); 
				t:=time()-t;							  							 
			    else
				error("invalid argument",args[3], "must be a Vector");				    
			    fi;	
			 else
			    error("invalid argument",args[2], "must be a Matrix");			 
			 fi;   
		      elif (type(args[1], Matrix)) then
	      		   if (type(args[2], Vector)) then 
			        t:=time(); 
				A:=Dev:-lbBlackbox(args[1]);
				b:=Dev:-lbVector(args[2]); 
				t:=time()-t;
			   else
				error("invalid argument",args[2], "must be a Vector");	
			   fi;
		      else
			   error("invalid argument",args[1], "must be a Matrix");
		      fi;		      
		      x:= Dev:-lbSolve(A,b);
		      y:=Dev:-lbConvertVector(x); 
		      if (printlevel>1) then
			      print("conversion time:",t);
		      end if;
		      return y;			      
		end;	
end module:


save_linbox:=proc(repo)
	global savelibname,libname:
	if (FileTools[Exists](cat(repo,"/maple.lib"))) then 
		fremove(cat(repo,"/maple.lib"));
		fremove(cat(repo,"/maple.ind"));
#		fremove(cat(repo,"/maple.hdb"));
		rmdir(repo);
	fi:
	if (not(FileTools[Exists](repo))) then mkdir(repo);fi:
	march( 'create',repo, 1000 ):
	savelibname := repo: libname := savelibname, libname:
	savelib('LinBox'):
end:

save_linbox("LinBoxMaple");


 
