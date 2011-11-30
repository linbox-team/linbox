#!/bin/bash - 

# Copyright(c) 2011 LinBox
# Written by BB <bboyer@imag.fr>
# see COPYING for more details.


# TODO : create .extracted, .configured, .build, .installed and use a switch
# (like --force) when you want to rebuild
# TODO : manage icc/gcc
# TODO : add gmp in givaro and use auto-install in givaro
# TODO : use an optionnal message in die function.

STABLE_GIVARO=3.2.16
GIV_TAR=133
GIV_MD5=134

#switches
DEBUG=""
DEBUG_VAR=""
WARNINGS=""
WARNINGS_VAR=""
OPTIM="--enable-optimization"
OPTIM_VAR=""
CHECK_VAR=""
#options
PREFIX_LOC="/tmp/linbox"
PREFIX_VAR=""
PREFIX="--prefix=$PREFIX_LOC"
BLAS=""
BLAS_VAR=""
NTL="--with-ntl"
NTL_VAR=""
EXTRA=""
EXTRA_VAR=""
IML="--with-iml"
IML_VAR=""
SAGE=""
SAGE_VAR="false"
DRIV=""
DRIV_VAR="false"

MAKEOPT= 
MAKE_VAR=""

DONE="\033[0;36m done !\033[0m"
BEG="\033[1;32m * \033[0m"

#########
#  die  #
#########

die() {
	echo -ne "\n\033[1;31m * \033[0mfailed" ;   
	if [[ -n $1 ]] ; then
		echo " ($1)"
	else 
		echo "."
	fi
	exit -1 ;  
}

cool() {
	echo -e $DONE
}

##############
#   helper   #
##############

help() {
	echo
	echo " script for building and installing linbox the simple way (hopefully)"
	echo
	echo " * usage :"
	echo

	echo " --prefix=MY/PATH      : install all libraries under MY/PATH."
	echo "                         Default : /tmp/"
	echo
	echo " >> Libraries to seach for <<" 
	echo 
	echo " --with-gmp=GMP/PATH   : tell where gmp is."
	echo "                         Default : /usr, /usr/local. No argument is Default"
	echo " --with-blas=BLAS/PATH : same as GMP for BLAS. (will check anyway)"
	echo " --with-ntl=NTL/PATH   : same as GMP for NTL. (default)"
	echo " --with-iml=IML/PATH   : same as GMP for IML. (default)"
	echo " --extra-flags=\"\"      : give extra compiler flags."
	echo "                         Default : empty"
	echo " --make-flags=\"\"      : give extra makefile flags."
	echo "                         Default : empty"
	echo
	echo " >> for the next switches, nothing, Y, y, yes or 1 after \"=\"   <<"
	echo " >> means enabled. Anything else or omission means disabled    <<"
	echo
	echo " --enable-debug        : build in debugging mode."
	echo "                         Default : disabled."
	echo " --enable-check        : run make check."
	echo "                         Default : disabled."
	echo " --enable-warnings     : build with extra compile warnings."
	echo "                         Default : disabled. May be \'full\' "
	echo " --enable-optimization : build with compile-time optimization."
	echo "                         Default : enabled."
	echo " --enable-sage         : build with sage support."
	echo "                         Default : disabled."
	echo " --enable-drivers      : build with drivers support."
	echo "                         Default : disabled."
	echo 
	echo " >> calling helllp <<"
	echo 
	echo " --help, -h, -?        : print help and exit."
}

############
#  parser  # 
############

for i in "$@" ; do
	case "$i" in
		# switches
		"--help"|"-h"|"-?") 
		help
		exit 0
		;;
	"--enable-debug")
		if [ "x$DEBUG_VAR" = "xfalse" ]  ; then  echo "enable-debug or not ?" ;        help ; exit -1 ; fi
		DEBUG="$i";
		DEBUG_VAR="true";
		;;
	"--enable-check")
		if [ "x$CHECK_VAR" = "xfalse" ] ; then   echo "enable-check or not ?";        help ; exit -1; fi
		CHECK_VAR="true";
		;;
	"--enable-warnings")
		if [ "x$WARNINGS_VAR" = "xfalse" ] ; then   echo "enable-warnings or not ?";     help ; exit -1; fi
		WARNINGS="$i";
		WARNINGS_VAR="true";
		;;
	"--enable-optimization")
		if	[ "x$OPTIM_VAR" = "xfalse" ] ; then   echo "enable-optimization or not ?"; help ; exit -1; fi
		OPTIM="$i";
		OPTIM_VAR="true";
		;;
	"--disable-debug")
		if [ "x$DEBUG_VAR" = "xtrue" ]  ; then  echo "enable-debug or not ?" ;        help ; exit -1 ; fi
		DEBUG_VAR="false";
		;;
	"--disable-check")
		if [ "x$CHECK_VAR" = "xtrue" ] ; then   echo "enable-check or not ?";        help ; exit -1; fi
		CHECK_VAR="false";
		;;
	"--disable-warnings")
		if [ "x$WARNINGS_VAR" = "xtrue" ] ; then   echo "enable-warnings or not ?";     help ; exit -1; fi
		WARNINGS_VAR="false";
		;;
	"--disable-optimization")
		if	[ "x$OPTIM_VAR" = "xtrue" ] ; then   echo "enable-optimization or not ?"; help ; exit -1; fi
		OPTIM_VAR="false";
		;;
	"--with-ntl")
		if	[ "x$NTL_VAR" = "xfalse" ] ; then   echo "with-ntl or not ?";            help ; exit -1; fi
		NTL="$i";
		NTL_VAR="true";
		;;
	"--with-gmp")
		if	[ "x$GMP_VAR" = "xfalse" ] ; then   echo "with-gmp or not ?";            help ; exit -1; fi
		GMP="$i"
		GMP_VAR="true"
		;;
	"--with-iml")
		if	[ "x$IML_VAR" = "xfalse" ] ; then   echo "with-iml or not ?";            help ; exit -1; fi
		IML="$i "
		IML_VAR="true"
		;;
	"--enable-sage")
		if	[ "x$SAGE_VAR" = "xfalse" ] ; then  echo "enable-sage or not ?";          help ; exit -1; fi
		SAGE="$i"
		SAGE_VAR="true"
		;;
	"--enable-drivers")
		if	[ "x$DRIV_VAR" = "xfalse" ] ; then  echo "enable-drivers or not ?" ;      help ; exit -1; fi
		DRIV="$i"
		DRIV_VAR="true"
		;;
	"--disable-sage")
		if	[ "x$SAGE_VAR" = "xtrue" ] ; then  echo "enable-sage or not ?";          help ; exit -1; fi
		SAGE=""
		SAGE_VAR="false"
		;;
	"--disable-drivers")
		if	[ "x$DRIV_VAR" = "xtrue" ] ; then  echo "enable-drivers or not ?" ;      help ; exit -1; fi
		DRIV=""
		DRIV_VAR="false"
		;;

	*)
		if [[ ! "$i" =~ --.*=.+ ]] ; then
			echo "bad switch : $i"
			help ;
			exit -1 ;
		fi
		# options (now we can cut)
		QUI="`echo $i | cut -d'=' -f1`"
		QUOI="`echo $i | cut -d'=' -f2`"
		# echo "$QUI = $QUOI"
		case "$QUI" in
			"--prefix")
				if		[ "x$PREFIX_VAR" = "xtrue" ] ; then  echo "prefix already set ?" ;      help ; exit -1; fi
				PREFIX=$i
				PREFIX_LOC=$QUOI
				PREFIX_VAR="true"
				;;
			"--extra-flags")
				if	[ "x$EXTRA_VAR" = "xtrue" ] ; then  echo "extra-flags already set ?" ;      help ; exit -1; fi
				EXTRA="$QUOI"
				EXTRA_VAR="true"
				;;
			"--make-flags")
				if	[ "x$MAKE_VAR" = "xtrue" ] ; then  echo "make-flags already set ?" ;      help ; exit -1; fi
				MAKEOPT="$QUOI"
				MAKE_VAR="true"
				;;
			"--with-gmp")
				if		[ "x$GMP_VAR" = "xtrue" ] ; then  echo "GMP path already set ?" ;      help ; exit -1; fi
				GMP="$i"
				GMP_VAR="true"
				;;
			"--with-blas")
				if	[ "x$BLAS_VAR" = "xtrue" ] ; then  echo "GMP path already set ?" ;      help ; exit -1; fi
				BLAS=$QUI=\"$QUOI\"
				BLAS_VAR="true"
				;;
			"--with-ntl")
				if		[ "x$NTL_VAR" = "xtrue" ] ; then  echo "NTL path already set ?" ;      help ; exit -1; fi
				NTL="$i"
				NTL_VAR="true"
				;;
			"--with-iml")
				if	[ "x$IML_VAR" = "xtrue" ] ; then  echo "IML path already set ?" ;      help ; exit -1; fi
				IML="$i"
				IML_VAR="true"
				;;
			"--enable-optimization")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && OK=1 || OK=0
				if		[ "x$OPTIM_VAR" = "xtrue"  -a "OK" = "0" ] ; then  echo "optim or not optim ?" ;      help ; exit -1; fi
				if		[ "x$OPTIM_VAR" = "xfalse" -a "OK" = "1" ] ; then  echo "optim or not optim ?" ;      help ; exit -1; fi
				if	[[ "x$OK" = "x1" ]] ; then  
					OPTIM=$QUI ; OPTIM_VAR="true" ;
				else
					OPTIM_VAR="false" ;
				fi
				;;
			"--enable-warnings")
				[[ "$QUOI" =~ y|yes|Y|1|full ]] && OK=1 || OK=0
				if [ "x$WARNING_VAR" = "xtrue"  -a "OK" = "0"  ] ; then  echo "warning or not warning ?" ;      help ; exit -1; fi
				if [ "x$WARNING_VAR" = "xfalse" -a "OK" = "1"  ] ; then  echo "warning or not warning ?" ;      help ; exit -1; fi
				if [[ "x$OK" = "x1" ]] ; then
					WARNINGS=$QUI ; WARNING_VAR="true"
				else
					WARNING_VAR="false" 
				fi
				[[ "x$QUOI" = "xfull" ]] && WARNINGS=$i
				;;
			"--enable-debug")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && OK=1 || OK=0
				if		[ "x$DEBUG_VAR" = "xtrue"  -a "OK" = "0"  ] ; then  echo "debug or not debug ?" ;      help ; exit -1; fi
				if		[ "x$DEBUG_VAR" = "xfalse" -a "OK" = "1"  ] ; then  echo "debug or not debug ?" ;      help ; exit -1; fi
				if		[[ "x$OK" = "x1" ]] ; then  
					DEBUG=$QUI ; DEBUG_VAR="true"
				else
					DEBUG_VAR="false" 
				fi
				;;
			"--enable-check")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && OK=1 || OK=0
				if		[ "x$CHECK_VAR" = "xtrue"  -a "OK" = "0"  ] ; then  echo "check or not check ?" ;      help ; exit -1; fi
				if		[ "x$CHECK_VAR" = "xfalse" -a "OK" = "1"  ] ; then  echo "check or not check ?" ;      help ; exit -1; fi
				if		[[ "x$OK" = "x1" ]] ; then 
					CHECK=$QUI ; CHECK_VAR="true" 
				else 
					CHECK_VAR="false"
				fi
				;;
			"--enable-sage")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && OK=1 || OK=0
				if		[ "x$OPTIM_VAR" = "xtrue"  -a "OK" = "0"  ] ; then  echo "sage or not sage ?" ;      help ; exit -1; fi
				if		[ "x$OPTIM_VAR" = "xfalse" -a "OK" = "1"  ] ; then  echo "sage or not sage ?" ;      help ; exit -1; fi
				if		[[ "x$OK" = "x1" ]] ; then 
					OPTIM=$QUI ; OPTIM_VAR="true" 
				else
					OPTIM_VAR="false" 
				fi
				;;
			"--enable-drivers")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && OK=1 || OK=0
				if		[ "x$OPTIM_VAR" = "xtrue"  -a "OK" = "0"  ] ; then  echo "drivers or not drivers ?" ;      help ; exit -1; fi
				if		[ "x$OPTIM_VAR" = "xfalse" -a "OK" = "1"  ] ; then  echo "drivers or not drivers ?" ;      help ; exit -1; fi
				if		[[ "x$OK" = "x1" ]] ; then 
					DRIV=$QUI ; DRIV_VAR="true"
				else
					DRIV_VAR="false" 
				fi
				;;
			*)
				echo "unkown swith option $i" ;
				help ;
				exit -1 ;
				;;
		esac
		;;
esac
done

MAKEPROG="make ${MAKEOPT}"

######################
#  create build dir  #
######################

echo -en "${BEG}Preparing build directory..."
if [ -e build ] ; then
	if [ ! -d build ] ; then
		rm -rf build ;
		mkdir build;
	fi
	# echo -n "emptying build directory..."
	# rm -rf build/
	# mkdir build
else
	# echo -n "creating empty build directory..."
	mkdir build
fi
cool

####################
#  fectch sources  #
####################

cd build ;

### Givaro ###

echo -en "${BEG}fecthing Givaro..."
if [ -f givaro-${STABLE_GIVARO}.tar.gz ] ; then
	echo "already there"
	cool
else
	wget http://ljk.imag.fr/CASYS/LOGICIELS/givaro/Downloads/givaro-${STABLE_GIVARO}.tar.gz  >/dev/null 2>&1 || die
fi


#####################
#  extract sources  #
#####################

### Givaro ###

OK=0
echo -en "${BEG}extracting Givaro..."
tar xzf givaro-$STABLE_GIVARO.tar.gz  && OK=1
[ "$OK" = "1" ] &&  cool   || die 


####################
#  install Givaro  #
####################

cd givaro-$STABLE_GIVARO || die

if [ -f Makefile ] ; then
	echo -e "${BEG}cleaning Givaro..."
	${MAKEPROG} clean || die
	${MAKEPROG} distclean || die 
	# ${MAKEPROG} unistall || die
	cool
fi

echo -e "${BEG}configuring Givaro..."

echo "./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS"
echo "./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS" > configure.givaro.exe
chmod +x configure.givaro.exe
./configure.givaro.exe
rm rf configure.givaro.exe
#./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS || die

echo -e "${BEG}building Givaro..."
echo "${MAKEPROG} CXXFLAGS+=\"$EXTRA\""

if [ -n "$EXTRA" ] ; then
	${MAKEPROG} "CXXFLAGS+=\"$EXTRA\"" || die
else
	${MAKEPROG} || die
fi

if [ "$CHECK_VAR" = "true" ] ; then
	echo -e "${BEG}checking Fflas-Ffpack..."
	${MAKEPROG} check || die
fi

echo -e "${BEG}installing Givaro..."
${MAKEPROG} install || die

#return in build
cd ..

cool

#return in linbox
cd ..

#####################
#  cleaning LinBox  #
#####################

if [ -f Makefile ] ; then
	echo -e "${BEG}cleaning LinBox..."
	${MAKEPROG} clean || die
	${MAKEPROG} distclean || die 
	# ${MAKEPROG} unistall || die
	cool
fi

echo -e "${BEG}configuring LinBox..."


GIVARO="--with-givaro=$PREFIX_LOC"

if [ -x autogen.sh ] ;  then 
	echo "./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $WARNINGS $IML $SAGE $DRIV"
	./autogen.sh "$PREFIX" "$DEBUG" "$OPTIM" "$GMP" "$BLAS" "$NTL" "$GIVARO" "$WARNINGS" "$IML" "$SAGE" "$DRIV" || die
else
	echo "./configure $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $WARNINGS $IML $SAGE $DRIV"
	# ./configure $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $WARNINGS  $IML $SAGE $DRIV || die
	./configure "$PREFIX" "$DEBUG" "$OPTIM" "$GMP" "$BLAS" "$NTL" "$GIVARO" "$WARNINGS" "$IML" "$SAGE" "$DRIV" || die
fi

echo -e "${BEG}building LinBox..."
echo "${MAKEPROG} CXXFLAGS+=\"$EXTRA\""

if [ -n "$EXTRA" ] ; then
	${MAKEPROG} "CXXFLAGS+=\"$EXTRA\"" || die
else
	${MAKEPROG} || die
fi

if [ "$CHECK_VAR" = "true" ] ; then
	echo -e "${BEG}checking LinBox..."
	${MAKEPROG} check || die
fi

echo -e "${BEG}installing LinBox..."
${MAKEPROG} install || die

cool

echo
echo -e "${BEG}Don't forget to run something like"
echo -e " *   'export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PREFIX_LOC/lib'"
echo -e " * to ensure you don't get undefined symbols !"
echo 
echo -e " * Happy LinBoxing ! (installed in $PREFIX_LOC)"
echo
cool
