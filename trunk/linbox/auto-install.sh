#!/bin/bash - 

# Copyright(c) 2011 LinBox
# Written by BB <bboyer@imag.fr>
# see COPYING for more details.


# TODO : create .extracted, .configured, .build, .installed and use a switch
# (like --force) when you want to rebuild
# TODO : manage icc/gcc
# TODO : add gmp in givaro and use auto-install in givaro
# TODO : use an optionnal message in die function.
# TODO : add options to make like '-j'

STABLE_FFLAS=1.4.0
STABLE_GIVARO=3.4.1
GIV_TAR=123
GIV_MD5=124

#switches
STABLE=true
DEBUG=""
WARNINGS=""
OPTIM="--enable-optimization"
CHECK=false
#options
PREFIX_LOC="/tmp"
PREFIX="--prefix=$PREFIX_LOC"
BLAS=""
NTL=""
EXTRA=""

DONE="\033[0;36m done !\033[0m"
BEG="\033[1;32m * \033[0m"

#########
#  die  #
#########

die() {
	echo -ne "\n\033[1;31m * \033[0mfailed" ;   
	if [ -n $1 ] ; then
		echo " ($1)"
	else 
		echo ""
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
	echo " --stable              : install latest stable versions or latest svn versions."
	echo "                          Default : yes."
	echo " --prefix=MY/PATH      : install all libraries under MY/PATH."
	echo "                         Default : /tmp/"
	echo " --with-gmp=GMP/PATH   : tell where gmp is."
	echo "                         Default : /usr, /usr/local"
	echo " --with-blas=BLAS/PATH : same for Blas installation."
	echo " --with-ntl=NTL/PATH   : same for NTL."
	echo " --extra-flags=\"\"      : give extra compiler flags."
	echo "                         Default : empty"
	echo " --enable-debug        : build in debugging mode."
	echo "                         Default : no."
	echo " --enable-check        : run make check."
	echo "                         Default : no."
	echo " --enable-warnings     : build with extra compile warnings."
	echo "                         Default : no."
	echo " --enable-optimization : build with compile-time optimization."
	echo "                         Default : yes."
	echo " --help, -h, -?        : print help and exit."
}

############
#  parser  # 
############

for i in $* ; do
	case "$i" in
		# switches
	"--help"|"-h"|"-?") 
		help
		exit 0
		;;
	"--no-stable")
		STABLE=false;
		;;
	"--stable")
		#default
		;;
	"--enable-debug")
		DEBUG=$i;
		;;
	"--no-enable-debug")
		DEBUG="";
		;;
	"--enable-check")
		CHECK=true;
		;;
	"--no-enable-check")
		CHECK=false;
		;;
	"--enable-warnings")
		WARNINGS=$i;
		;;
	"--no-enable-warnings")
		WARNINGS="";
		;;
	"--enable-optimization")
		OPTIM=$i;
		;;
	"--no-enable-optimization")
		OPTIM="";
		;;
	"--with-ntl")
		NTL=$i
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
				PREFIX=$i
				PREFIX_LOC=$QUOI
				;;
			"--extra-flags")
				EXTRA=$QUOI
				;;
			"--with-gmp")
				GMP=$i
				;;
			"--with-blas")
				BLAS=$i
				;;
			"--with-ntl")
				NTL=$i
				;;
			"--enable-optimization")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && OPTIM="--enable-optimization" || OPTIM=""
				;;
			"--enable-warnings")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && WARNINGS="--enable-warnings" || WARNINGS=""
				;;
			"--enable-debug")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && DEBUG="--enable-debug" || DEBUG=""
				;;
			"--enable-check")
				[[ "$QUOI" =~ y|yes|Y|1 ]] && CHECK="true" || CHECK="false"
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

### Fflas-ffpack ###

echo -en "${BEG}fecthing Fflas-Ffpack..."
if [ "$STABLE" = "true" ]; then
	if [ -f fflas-ffpack-$STABLE_FFLAS.tar.gz ] ; then
		echo "already there"
	   	echo -ne "${BEG}fetching md5sum" ; 
		[ -f fflas-ffpack-$STABLE_FFLAS.tar.gz.md5sum ] && rm fflas-ffpack-${STABLE_FFLAS}.tar.gz.md5sum ;
		wget http://linalg.org/fflas-ffpack-$STABLE_FFLAS.tar.gz.md5sum >/dev/null 2>&1 || die
		[ -f fflas-ffpack-$STABLE_FFLAS.tar.gz.md5sum ] || die
		cool
		echo -ne "${BEG}"
		md5sum -c fflas-ffpack-$STABLE_FFLAS.tar.gz.md5sum || die
	else
		wget http://linalg.org/fflas-ffpack-$STABLE_FFLAS.tar.gz  >/dev/null 2>&1 || die
		[ -f fflas-ffpack-$STABLE_FFLAS.tar.gz ] || die
	   	echo -ne "${BEG}fetching md5sum" ; 
		wget http://linalg.org/fflas-ffpack-$STABLE_FFLAS.tar.gz.md5sum >/dev/null 2>&1 || die
		cool
		echo -ne "${BEG}"
		md5sum -c fflas-ffpack-$STABLE_FFLAS.tar.gz.md5sum || die
	fi
else
	OK=0 ;
	svn co svn://linalg.org/fflas-ffpack 2>&1 >/dev/null && OK=1 
	[ "$OK" = "1" ] &&  cool  || die
fi

### Givaro ###

echo -en "${BEG}fecthing Givaro..."
if [ "$STABLE" = "true" ]; then
	if [ -f givaro-$STABLE_GIVARO.tar.gz ] ; then
		echo "already there"
	   	echo -ne "${BEG}fetching md5sum" ; 
		[ -f givaro-$STABLE_GIVARO.tar.gz.md5sum ] && rm givaro-${STABLE_GIVARO}.tar.gz.md5sum ;
		wget --no-check-certificate https://forge.imag.fr/frs/download.php/$GIV_MD5/givaro-$STABLE_GIVARO.tar.gz.md5sum >/dev/null 2>&1 || die
		[ -f givaro-$STABLE_GIVARO.tar.gz.md5sum ] || die
		cool
		echo -ne "${BEG}"
		md5sum -c givaro-$STABLE_GIVARO.tar.gz.md5sum || die
	else
		wget --no-check-certificate https://forge.imag.fr/frs/download.php/$GIV_TAR/givaro-$STABLE_GIVARO.tar.gz >/dev/null 2>&1 || die
		[ -f givaro-$STABLE_GIVARO.tar.gz ] || die
	   	echo -ne "${BEG}fetching md5sum" ; 
		wget --no-check-certificate https://forge.imag.fr/frs/download.php/$GIV_MD5/givaro-$STABLE_GIVARO.tar.gz.md5sum >/dev/null 2>&1 || die
		cool
		echo -ne "${BEG}"
		md5sum -c givaro-$STABLE_GIVARO.tar.gz.md5sum || die
	fi
else
	OK=0 ;
	svn co svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/givaro/trunk 2>&1 >/dev/null && OK=1 
	[ "$OK" = "1" ] &&  cool  || die 
fi


#####################
#  extract sources  #
#####################

### Fflas-ffpack ###

OK=0
if [ "$STABLE" = "true" ]; then
	echo -en "${BEG}extracting Fflas-Ffpack..."
	tar xzf fflas-ffpack-$STABLE_FFLAS.tar.gz  && OK=1
	[ "$OK" = "1" ] &&  cool   || die
fi

### Givaro ###

OK=0
if [ "$STABLE" = "true" ]; then
	echo -en "${BEG}extracting Givaro..."
	tar xzf givaro-$STABLE_GIVARO.tar.gz  && OK=1
	[ "$OK" = "1" ] &&  cool   || die 
fi

##########################
#  install fflas-ffpack  #
##########################

if [ "$STABLE" = "true" ]; then
	cd fflas-ffpack-$STABLE_FFLAS/ || die
else
	cd fflas-ffpack/ || die
fi


if [ -f Makefile ] ; then
	echo -e "${BEG}cleaning Fflas-Ffpack..."
	make clean || die
	make distclean || die 
	# make unistall || die
	cool
fi

echo -e "${BEG}configuring Fflas-Ffpack..."

if [ "$STABLE" = "true" ]; then
	echo "./configure  $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS"
	./configure  $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS || die
else
	echo "./autogen.sh $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS"
	./autogen.sh $PREFIX $DEBUG $OPTIM $BLAS $WARNINGS || die
fi

echo -e "${BEG}building Fflas-Ffpack..."
echo "make CXXFLAGS+=\"$EXTRA\""
if [ -n "$EXTRA" ] ; then
	make "CXXFLAGS+=\"$EXTRA\"" || die
else
	make || die
fi

if [ "$CHECK" = "true" ] ; then
	echo -e "${BEG}checking Fflas-Ffpack..."
	make check || die
fi


echo -e "${BEG}installing Fflas-Ffpack..."
make install || die

cool
#return in build
cd ..

####################
#  install Givaro  #
####################

if [ "$STABLE" = "true" ]; then
	cd givaro-$STABLE_GIVARO || die
else
	cd trunk/ || die
fi

if [ -f Makefile ] ; then
	echo -e "${BEG}cleaning Givaro..."
	make clean || die
	make distclean || die 
	# make unistall || die
	cool
fi

echo -e "${BEG}configuring Givaro..."

if [ "$STABLE" = "true" ]; then
	echo "./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS"
	./configure  $PREFIX $DEBUG $OPTIM $GMP $WARNINGS || die
else
	echo "./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $WARNINGS"
	./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $WARNINGS || die
fi

echo -e "${BEG}building Givaro..."
echo "make CXXFLAGS+=\"$EXTRA\""

if [ -n "$EXTRA" ] ; then
	make "CXXFLAGS+=\"$EXTRA\"" || die
else
	make || die
fi

if [ "$CHECK" = "true" ] ; then
	echo -e "${BEG}checking Fflas-Ffpack..."
	make check || die
fi

echo -e "${BEG}installing Givaro..."
make install || die

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
	make clean || die
	make distclean || die 
	# make unistall || die
	cool
fi

echo -e "${BEG}configuring LinBox..."


GIVARO="--with-givaro=$PREFIX_LOC"
FFLAFLAS="--with-fflas-ffpack=$PREFIX_LOC"

if [ -x autogen.sh ] ;  then 
	echo "./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $FFLAFLAS $WARNINGS"
	./autogen.sh $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $FFLAFLAS $WARNINGS || die
else
	echo "./configure $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $FFLAFLAS $WARNINGS"
	./configure $PREFIX $DEBUG $OPTIM $GMP $BLAS $NTL $GIVARO $FFLAFLAS $WARNINGS || die
fi

echo -e "${BEG}building LinBox..."
echo "make CXXFLAGS+=\"$EXTRA\""

if [ -n "$EXTRA" ] ; then
	make "CXXFLAGS+=\"$EXTRA\"" || die
else
	make || die
fi

if [ "$CHECK" = "true" ] ; then
	echo -e "${BEG}checking LinBox..."
	make check || die
fi

echo -e "${BEG}installing LinBox..."
make install || die

cool

echo
echo -e "${BEG}Don't forget to run something like"
echo -e " *   'export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$PREFIX_LOC/lib'"
echo -e " * to ensure you don't get undefined symobols !"
echo 
echo -e " * Happy LinBoxing ! (installed in $PREFIX_LOC)"
echo
cool
