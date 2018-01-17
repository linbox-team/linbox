#!/bin/bash - 

# Copyright(c) 2011 LinBox
# Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
#         Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
# see COPYING for more details.


# TODO : create .extracted, .configured, .build, .installed and use a switch
# (like --force) when you want to rebuild
# TODO : manage icc/gcc
# TODO : add gmp in givaro and use auto-install in givaro
# TODO : use an optionnal message in die function.
# TODO : put openblas and use auto-install in fflas-ffpack
# TODO : rpath instead of LD_LIBRARY_PATH ?

#############################
## Only for stable fetching 
#############################
#### stable .tar.gz 
STABLE_LB=1.5.2
STABLE_FFLAS=2.3.1
STABLE_GIVARO=4.0.4
STABLE_OPENBLAS=0.2.20
MD5SUFF=md5
#############################

function decompress {
#tar xf $1
    gunzip -c $1 | tar xf -
}



#switches
STABLE_VAR="false"
DEBUG=""
DEBUG_VAR=""
WARNINGS=""
WARNINGS_VAR=""
OPTIM="--enable-optimization"
OPTIM_VAR=""
CHECK_VAR=""
#options
PREFIX_LOC="/tmp"
PREFIX_VAR=""
PREFIX="--prefix=$PREFIX_LOC"
BLAS=""
BLAS_VAR="false"
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
OPENBLAS=""
OPENBLAS_VAR="false"
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
    echo " --stable=[yes,no]     : install latest stable versions or latest git versions."
    echo "                         Default : no, even if switch ommitted. No argument means no"

    echo " --prefix=MY/PATH      : install all libraries under MY/PATH."
    echo "                         Default : /tmp"
    echo
    echo " >> Libraries to search for <<" 
    echo 
    echo " If some library cannot be linked, don't forget to export LD_LIBRARY_PATH ! "
    echo 
    echo " --with-gmp=GMP/PATH   : tell where gmp is."
    echo "                         Default : /usr, /usr/local. No argument is Default"
    echo " --with-blas-libs=BLAS/PATH : same as GMP for BLAS. (will check anyway)"
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
    echo " --enable-openblas     : fetch openblas."
    echo "                         Default : disabled (use local blas)."
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

# Recover command line, with double-quotes
CMDLINE=""
for arg in "$@"
  do
  WHO="`echo $arg | cut -d'=' -f1`"
  WHAT="`echo $arg | cut -s -d'=' -f2`"
  if test "x$WHAT" = "x"; then
      CMDLINE="$CMDLINE $WHO"
  else
      CMDLINE="$CMDLINE $WHO=\"$WHAT\""
  fi
done
echo  "$0 $CMDLINE" | tee -a auto-install.log

# Parse command line
for i in "$@" ; do
    case "$i" in
		# switches
	"--help"|"-h"|"-?") 
	help
	exit 0
	;;
	"--stable")
	if [ "x$STABLE_VAR" = "xfalse" ] ; then echo "stable or not ?";              help ; exit -1 ;  fi
	STABLE_VAR=true;
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
	"--enable-openblas")
	if [ "x$OPENBLAS_VAR" = "xfalse" ]  ; then  echo "enable-openblas or not ?" ;        help ; exit -1 ; fi
	OPENBLAS="$i";
	OPENBLAS_VAR="true";
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
	if [ "x$OPTIM_VAR" = "xtrue" ] ; then   echo "enable-optimization or not ?"; help ; exit -1; fi
	OPTIM_VAR="false";
	;;
	"--disable-openblas")
	if [ "x$OPENBLAS_VAR" = "xtrue" ]  ; then  echo "enable-openblas or not ?" ;        help ; exit -1 ; fi
	OPENBLAS_VAR="false";
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
	# echo "QUIQUOI: $QUI = $QUOI"
	case "$QUI" in
	    "--stable")
	    OK=2
	    [[ "x$QUOI" = "xyes" ]]  &&  OK=1  
	    [[ "x$QUOI" = "xno"  ]]  &&  OK=0   
	    if [[ "$OK" = "2" ]] ; then 
		echo "stable=[yes/no] !" ; help ; exit -1 ; 
	    fi
	    [[ "$OK" = "1" ]] && STABLE_VAR="true" || STABLE_VAR="false"
	    ;;
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
	    "--with-blas-libs")
	    if	[ "x$BLAS_VAR" = "xtrue" ] ; then  echo "BLAS path already set ?" ;      help ; exit -1; fi
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
	    "--enable-openblas")
	    [[ "$QUOI" =~ y|yes|Y|1 ]] && OK=1 || OK=0
	    if		[ "x$OPENBLAS_VAR" = "xtrue"  -a "OK" = "0" ] ; then  echo "openblas or not openblas ?" ;      help ; exit -1; fi
	    if		[ "x$OPENBLAS_VAR" = "xfalse" -a "OK" = "1" ] ; then  echo "openblas or not openblas ?" ;      help ; exit -1; fi
	    if	[[ "x$OK" = "x1" ]] ; then  
		OPENBLAS=$QUI ; OPENBLAS_VAR="true" ;
	    else
		OPENBLAS_VAR="false" ;
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
	    if		[ "x$SAGE_VAR" = "xtrue"  -a "OK" = "0"  ] ; then  echo "sage or not sage ?" ;      help ; exit -1; fi
	    if		[ "x$SAGE_VAR" = "xfalse" -a "OK" = "1"  ] ; then  echo "sage or not sage ?" ;      help ; exit -1; fi
	    if		[[ "x$OK" = "x1" ]] ; then 
		SAGE=$QUI ; SAGE_VAR="true" 
	    else
		SAGE_VAR="false" 
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
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${PREFIX_LOC}/lib/pkgconfig
echo "PKG_CONFIG_PATH=$PKG_CONFIG_PATH"| tee -a auto-install.log
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PREFIX_LOC/lib
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"| tee -a auto-install.log



##################################
#  check within LinBox or fetch  #
##################################

if [ ! \( -x autogen.sh -o -x configure \) ] ; then

### Extract LinBox sources ###
    echo -en "${BEG}fetching LinBox..."| tee -a auto-install.log
    if [ "$STABLE_VAR" = "true" ]; then
	if [ -f linbox-${STABLE_LB}.tar.gz ] ; then
	    echo -ne " already there!\n"| tee -a auto-install.log
	    echo -ne "${BEG}fetching md5sum" | tee -a auto-install.log;
	    [ -f fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} ] && rm fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} ;
	    wget --no-check-certificate https://github.com/linbox-team/linbox/releases/download/v${STABLE_LB}/linbox-${STABLE_LB}.tar.gz.${MD5SUFF} >/dev/null 2>&1 || die
	    [ -f linbox-${STABLE_LB}.tar.gz.${MD5SUFF} ] || die
	    cool| tee -a auto-install.log
	    echo -ne "${BEG}"
	    md5sum -c linbox-${STABLE_LB}.tar.gz.${MD5SUFF} || die
	else
	    wget https://github.com/linbox-team/linbox/releases/download/v${STABLE_LB}/linbox-${STABLE_LB}.tar.gz >/dev/null 2>&1 || die
	    [ -f linbox-${STABLE_LB}.tar.gz ] &&  cool || die
	    echo -ne "${BEG}fetching md5sum" | tee -a auto-install.log; 
	    wget --no-check-certificate https://github.com/linbox-team/linbox/releases/download/v${STABLE_LB}/linbox-${STABLE_LB}.tar.gz.${MD5SUFF} >/dev/null 2>&1 || die
	    cool| tee -a auto-install.log
	    echo -ne "${BEG}"
	    md5sum -c linbox-${STABLE_LB}.tar.gz.${MD5SUFF} || die
	fi
	OK=0
	echo -en "${BEG}extracting LinBox..."
	decompress linbox-${STABLE_LB}.tar.gz  && OK=1
	[ "$OK" = "1" ] &&  cool | tee -a auto-install.log  || die 
	cd linbox-${STABLE_LB} &&  cool   || die 
    else
	OK=0 ;
	git clone --depth 1 https://github.com/linbox-team/linbox.git 2>&1 >/dev/null && OK=1
	[ "$OK" = "1" ] &&  cool | tee -a auto-install.log || die
	cd linbox &&  cool   || die 
    fi
    mv ../auto-install.log . || die
fi

######################
#  create build dir  #
######################

#first tee creates a new log.
echo -en "${BEG}Preparing build directory..."| tee -a auto-install.log
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
cool| tee -a auto-install.log

####################
#  fetch sources  #
####################

cd build ;

### Givaro ###

echo -en "${BEG}fetching Givaro..."| tee -a ../auto-install.log
if [ "$STABLE_VAR" = "true" ]; then
    if [ -f givaro-${STABLE_GIVARO}.tar.gz ] ; then
	echo -ne " already there!\n"
	echo -ne "${BEG}fetching md5sum" ; 
	[ -f givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} ] && rm givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} ;
	wget --no-check-certificate https://github.com/linbox-team/givaro/releases/download/v${STABLE_GIVARO}/givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} >/dev/null 2>&1 || die
	[ -f givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} ] || die
	cool| tee -a ../auto-install.log
	echo -ne "${BEG}"
	md5sum -c givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} || die
    else
	wget --no-check-certificate https://github.com/linbox-team/givaro/releases/download/v${STABLE_GIVARO}/givaro-${STABLE_GIVARO}.tar.gz >/dev/null 2>&1 || die
	[ -f givaro-${STABLE_GIVARO}.tar.gz ] &&  cool || die
	echo -ne "${BEG}fetching md5sum" ; 
	wget --no-check-certificate https://github.com/linbox-team/givaro/releases/download/v${STABLE_GIVARO}/givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} >/dev/null 2>&1 || die
	cool| tee -a ../auto-install.log
	echo -ne "${BEG}"
	md5sum -c givaro-${STABLE_GIVARO}.tar.gz.${MD5SUFF} || die
    fi
else
    OK=0 ;
    git clone --depth 1 https://github.com/linbox-team/givaro.git 2>&1 >/dev/null && OK=1
    [ "$OK" = "1" ] &&  cool | tee -a ../auto-install.log || die 
fi

### Fflas-ffpack ###

echo -en "${BEG}fetching Fflas-Ffpack..."| tee -a ../auto-install.log
if [ "$STABLE_VAR" = "true" ]; then
    if [ -f fflas-ffpack-${STABLE_FFLAS}.tar.gz ] ; then
	echo -ne " already there!\n"
	echo -ne "${BEG}fetching md5sum" ; 
	[ -f fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} ] && rm fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} ;
	wget --no-check-certificate https://github.com/linbox-team/fflas-ffpack/releases/download/v${STABLE_FFLAS}/fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} >/dev/null 2>&1 || die
	[ -f fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} ] &&  cool || die
	cool
	echo -ne "${BEG}"
	md5sum -c fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} || die
    else
	wget --no-check-certificate https://github.com/linbox-team/fflas-ffpack/releases/download/v${STABLE_FFLAS}/fflas-ffpack-${STABLE_FFLAS}.tar.gz >/dev/null 2>&1 || die
	[ -f fflas-ffpack-${STABLE_FFLAS}.tar.gz ] &&  cool || die
	echo -ne "${BEG}fetching md5sum" ; 
	wget --no-check-certificate https://github.com/linbox-team/fflas-ffpack/releases/download/v${STABLE_FFLAS}/fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} >/dev/null 2>&1 || die
	cool
	echo -ne "${BEG}"
	md5sum -c fflas-ffpack-${STABLE_FFLAS}.tar.gz.${MD5SUFF} || die
    fi
else
    OK=0 ;
    git clone --depth=1 https://github.com/linbox-team/fflas-ffpack.git 2>&1 >/dev/null && OK=1
    [ "$OK" = "1" ] &&  cool | tee -a ../auto-install.log || die
fi


### OpenBlas ###
if [ "$OPENBLAS_VAR" = "true" ]; then
    echo -en "${BEG}fetching OpenBlas..."| tee -a ../auto-install.log
    if [ "$STABLE_VAR" = "true" ]; then
	if [ -f v${STABLE_OPENBLAS}.tar.gz ] ; then
	    echo -ne " already there!\n"
	else
	    wget --no-check-certificate http://github.com/xianyi/OpenBLAS/archive/v${STABLE_OPENBLAS}.tar.gz >/dev/null 2>&1 || die
	    [ -f v${STABLE_OPENBLAS}.tar.gz ] &&  cool || die
	fi
    else
	OK=0 ;
	git clone --depth=1 https://github.com/xianyi/OpenBLAS.git 2>&1 >/dev/null && OK=1
	[ "$OK" = "1" ] &&  cool | tee -a ../auto-install.log || die
    fi
fi

#####################
#  extract sources  #
#####################

### Givaro ###

OK=0
if [ "$STABLE_VAR" = "true" ]; then
    echo -en "${BEG}extracting Givaro..."| tee -a ../auto-install.log
    decompress givaro-${STABLE_GIVARO}.tar.gz  && OK=1
    [ "$OK" = "1" ] &&  cool | tee -a ../auto-install.log  || die 
fi

### Fflas-ffpack ###

OK=0
if [ "$STABLE_VAR" = "true" ]; then
    echo -en "${BEG}extracting Fflas-Ffpack..."| tee -a ../auto-install.log
    decompress fflas-ffpack-${STABLE_FFLAS}.tar.gz  && OK=1
    [ "$OK" = "1" ] &&  cool  | tee -a ../auto-install.log || die
fi

### OpenBlas ###

if [ "$OPENBLAS_VAR" = "true" ]; then
    OK=0
    if [ "$STABLE_VAR" = "true" ]; then
	echo -en "${BEG}extracting OpenBlas..."| tee -a ../auto-install.log
	decompress v${STABLE_OPENBLAS}.tar.gz  && OK=1
	[ "$OK" = "1" ] &&  cool | tee -a ../auto-install.log  || die
    fi
fi

####################
#  install Givaro  #
####################

if [ "$STABLE_VAR" = "true" ]; then
    cd givaro-${STABLE_GIVARO} || die
else
    cd givaro/ || die
fi

if [ -f Makefile ] ; then
    echo -e "${BEG}cleaning Givaro..."| tee -a ../../auto-install.log
    ${MAKEPROG} clean | tee -a ../../auto-install.log|| die
    ${MAKEPROG} distclean | tee -a ../../auto-install.log|| die 
	# ${MAKEPROG} unistall || die
    cool
fi

echo -e "${BEG}configuring Givaro..."

if [ "$STABLE_VAR" = "true" ]; then
    echo "./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS "
    echo "./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS " > configure.givaro.exe
    chmod +x configure.givaro.exe
    ./configure.givaro.exe | tee -a ../../auto-install.log
    rm -rf configure.givaro.exe
	#./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS || die
else
    echo "./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $WARNINGS"
    echo "./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $WARNINGS" > autogen.givaro.exe
    chmod +x autogen.givaro.exe
    ./autogen.givaro.exe| tee -a ../../auto-install.log
    rm -rf autogen.givaro.exe
	#./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $WARNINGS || die
fi

echo -e "${BEG}building Givaro..."| tee -a ../../auto-install.log
echo "${MAKEPROG} CXXFLAGS+=\"$EXTRA\" LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\""| tee -a ../../auto-install.log

if [ -n "$EXTRA" ] ; then
    ${MAKEPROG} "CXXFLAGS+=\"$EXTRA\" LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\"" | tee -a ../../auto-install.log|| die
else
    ${MAKEPROG} | tee -a ../../auto-install.log|| die
fi

if [ "$CHECK_VAR" = "true" ] ; then
    echo -e "${BEG}checking Fflas-Ffpack..."| tee -a ../../auto-install.log
    ${MAKEPROG} check | tee -a ../../auto-install.log|| die
fi

echo -e "${BEG}installing Givaro..."| tee -a ../../auto-install.log
${MAKEPROG} install | tee -a ../../auto-install.log|| die

#return in build
cd ..

cool| tee -a ../auto-install.log

######################
#  install OpenBlas  #
######################

if [ "$OPENBLAS_VAR" = "true" ]; then

    if [ "$STABLE_VAR" = "true" ]; then
	cd OpenBLAS-${STABLE_OPENBLAS} || die
    else
	cd OpenBLAS/ || die
    fi
    
    if [ -f Makefile ] ; then
	echo -e "${BEG}cleaning OpenBLAS..."| tee -a ../../auto-install.log
	${MAKEPROG} clean | tee -a ../../auto-install.log|| die
	${MAKEPROG} distclean | tee -a ../../auto-install.log|| die 
	# ${MAKEPROG} unistall || die
	cool
    fi
    
    OPENBLAS_FLAGS="USE_THREADS=0 USE_THREAD=0 CC=gcc FC=gfortran PREFIX=$PREFIX_LOC"
    
    echo -e "${BEG}building OpenBLAS..."| tee -a ../../auto-install.log
    echo "${MAKEPROG} ${OPENBLAS_FLAGS} CXXFLAGS+=\"$EXTRA\" LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\""| tee -a ../../auto-install.log
    
    if [ -n "$EXTRA" ] ; then
	${MAKEPROG} ${OPENBLAS_FLAGS} "CXXFLAGS+=\"$EXTRA\" LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\"" | tee -a ../../auto-install.log|| die
    else
	${MAKEPROG} ${OPENBLAS_FLAGS} | tee -a ../../auto-install.log|| die
    fi
    
    echo -e "${BEG}installing OpenBLAS..."| tee -a ../../auto-install.log
    ${MAKEPROG} ${OPENBLAS_FLAGS} install | tee -a ../../auto-install.log|| die
    
#return in build
    cd ..
    
    cool| tee -a ../auto-install.log

    if [ "$BLAS_VAR" = "false" ]; then
	BLAS="--with-blas-libs="\""-L${PREFIX_LOC}/lib -lopenblas -lpthread -lgfortran"\"
	BLAS_VAR=true
    fi
fi


##########################
#  install fflas-ffpack  #
##########################

if [ "$STABLE_VAR" = "true" ]; then
    cd fflas-ffpack-${STABLE_FFLAS}/ || die
else
    cd fflas-ffpack/ || die
fi


if [ -f Makefile ] ; then
    echo -e "${BEG}cleaning Fflas-Ffpack..."| tee -a ../../auto-install.log
    ${MAKEPROG} clean | tee -a ../../auto-install.log|| die
    ${MAKEPROG} distclean | tee -a ../../auto-install.log|| die 
	# ${MAKEPROG} unistall || die
    cool| tee -a ../../auto-install.log
fi

echo -e "${BEG}configuring Fflas-Ffpack..."| tee -a ../../auto-install.log

if [ "$STABLE_VAR" = "true" ]; then
    echo "./configure  $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS"| tee -a ../../auto-install.log
    echo "./configure  $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS" > configure.fflas.exe
    chmod +x configure.fflas.exe
    ./configure.fflas.exe| tee -a ../../auto-install.log
    rm -rf configure.fflas.exe
	#./configure  "$PREFIX" "$DEBUG" "$OPTIM" "$BLAS"  "$WARNINGS" || die
else
    echo "./autogen.sh $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS"| tee -a ../../auto-install.log
    echo "./autogen.sh $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS" > configure.fflas.exe
    chmod +x configure.fflas.exe
    ./configure.fflas.exe| tee -a ../../auto-install.log
    rm -rf configure.fflas.exe
	#./autogen.sh "$PREFIX" "$DEBUG" "$OPTIM" "$BLAS"  "$WARNINGS" || die
fi

echo -e "${BEG}building Fflas-Ffpack..."| tee -a ../../auto-install.log
echo "${MAKEPROG} CXXFLAGS+=\"$EXTRA\""
if [ -n "$EXTRA" ] ; then
    ${MAKEPROG} "CXXFLAGS+=\"$EXTRA\"" | tee -a ../../auto-install.log|| die
else
    ${MAKEPROG} | tee -a ../../auto-install.log|| die
fi

if [ "$CHECK_VAR" = "true" ] ; then
    echo -e "${BEG}checking Fflas-Ffpack..."| tee -a ../../auto-install.log
    ${MAKEPROG} check | tee -a ../../auto-install.log|| die
fi


echo -e "${BEG}installing Fflas-Ffpack..."
${MAKEPROG} install | tee -a ../../auto-install.log|| die

cool| tee -a ../../auto-install.log
#return in build
cd ..

#return in linbox
cd ..


#####################
#  cleaning LinBox  #
#####################

if [ -f Makefile ] ; then
    echo -e "${BEG}cleaning LinBox..."| tee -a ./auto-install.log
    ${MAKEPROG} clean | tee -a ./auto-install.log|| die
    ${MAKEPROG} distclean | tee -a ./auto-install.log|| die 
	# ${MAKEPROG} unistall || die
    cool| tee -a ./auto-install.log
fi

echo -e "${BEG}configuring LinBox..."| tee -a ./auto-install.log

echo ""| tee -a ./auto-install.log
echo -e "${BEG}Don't forget to run something like"| tee -a ./auto-install.log
echo -e " *   'export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PREFIX_LOC/lib'"| tee -a ./auto-install.log
echo -e " * to ensure you don't get undefined symbols !"| tee -a ./auto-install.log
echo  ""| tee -a ./auto-install.log

if [ -x autogen.sh ] ;  then 
    echo "./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $WARNINGS $IML $SAGE $DRIV"| tee -a ./auto-install.log
    ./autogen.sh "$PREFIX" "$DEBUG" "$OPTIM" "$GMP" "$BLAS" "$NTL" "$WARNINGS" "$IML" "$SAGE" "$DRIV" | tee -a ./auto-install.log|| die
else
    echo "./configure $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $WARNINGS $IML $SAGE $DRIV"| tee -a ./auto-install.log
	# ./configure $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $WARNINGS  $IML $SAGE $DRIV || die
    ./configure "$PREFIX" "$DEBUG" "$OPTIM" "$GMP" "$BLAS" "$NTL" "$WARNINGS" "$IML" "$SAGE" "$DRIV" | tee -a ./auto-install.log|| die
fi

echo -e "${BEG}building LinBox..."| tee -a ./auto-install.log
echo "${MAKEPROG} CXXFLAGS+=\"$EXTRA\" LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\""| tee -a ./auto-install.log

if [ -n "$EXTRA" ] ; then
    ${MAKEPROG} "CXXFLAGS+=\"$EXTRA\" LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\"" | tee -a ./auto-install.log|| die
else
    ${MAKEPROG} "LDFLAGS+=\"-Wl,-rpath,$PREFIX_LOC\""| tee -a ./auto-install.log|| die
fi

if [ "$CHECK_VAR" = "true" ] ; then
    echo -e "${BEG}checking LinBox..."| tee -a auto-install.log
    ${MAKEPROG} check | tee -a auto-install.log|| die
fi

echo -e "${BEG}installing LinBox..."| tee -a auto-install.log
${MAKEPROG} install | tee -a auto-install.log|| die

cool| tee -a auto-install.log

echo    " " | tee -a auto-install.log
echo -e "${BEG}Don't forget to run something like"| tee -a auto-install.log
echo -e " *   'export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PREFIX_LOC/lib'"| tee -a auto-install.log
echo -e " * to ensure you don't get undefined symbols !"| tee -a auto-install.log
echo  "" | tee -a auto-install.log
echo -e " * Happy LinBoxing ! (installed in $PREFIX_LOC)"| tee -a auto-install.log
echo " "| tee -a auto-install.log
cool| tee -a auto-install.log
