#!/bin/bash
# This file is part of the LinBox library.
# It is distributed under the terms of the LGPL licence version 2.1 or later 
# (see COPYING)
# Created by AB - 2014/12/03
# Modified by AC - 2016/06/20
# Modified by CP - 2016/06/22

# Some influential environment variables:
#	CXX			C++ compiler command
#	CXXFLAGS	C++ compiler flags

# Note: This script is intended to be launched
# by the Jenkins web interface whenever it needs
# to compile the project.
# It is launched from the svn:trunk root directory.
# But should be stored in /<slave_jenkins_path>/makers/

SOURCE_DIRECTORY=$( cd "$( dirname "$0" )" && pwd )

#=============================#
# Change only these variables #
#=============================#
ARCH=`pwd | awk -F/ '{print $(NF-4)}'`
CXX=`pwd | awk -F/ '{print $(NF-2)}'`
NTL=`pwd | awk -F/ '{print $NF}'`
JENKINS_DIR=${SOURCE_DIRECTORY%%/workspace/*}
LOCAL_DIR="$JENKINS_DIR"/local
# Add path to compilers (if needed)
export PATH=$PATH:/usr/local/bin:"$LOCAL_DIR/$CXX/bin"
echo $PATH
# Add specific locations (if needed)
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":/usr/local/lib:"$LOCAL_DIR/$CXX/lib":"$LOCAL_DIR/$CXX/lib"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="$LOCAL_DIR/$CXX/lib/pkgconfig"
echo "PKG_CONFIG_PATH = $PKG_CONFIG_PATH"
VERSION=`grep VERSION Makefile.am  | sed s/VERSION=// `
echo "VERSION = $VERSION"
distdir=linbox-$VERSION
distarchive=$distdir.tar.gz
echo "distarchive = $distarchive"

# Where to install linbox binaries
# Keep default for local installation.
PREFIX_INSTALL="$LOCAL_DIR/$CXX/$NTL"

# Job Linbox with Ntl option flag
if [ "$NTL" == "withNTL" ]; then
  LINBOX_NTLFLAG="--with-ntl=$LOCAL_DIR/$CXX"
fi

# /!\ Warning /!\ This could be an issue if you changed
# the local installation directory
rm -rf "$PREFIX_INSTALL"/bin/linbox* "$PREFIX_INSTALL"/include/linbox* "$PREFIX_INSTALL"/lib/liblinbox*

#================#
# Setup Variables#
#================#

if [ "$CXX" == "icpc" ]; then
     distribution=`uname -m`
     if [ "$distribution" == "i686" ]; then 	
	source /usr/local/bin/compilervars.sh ia32
     else
	source /usr/local/bin/compilervars.sh intel64
     fi
fi

#==================================#
# Automated installation and tests #
#==================================#

echo "|=== JENKINS AUTOMATED SCRIPT ===| ./autogen.sh CXX=$CXX CXXFLAGS=$CXXFLAGS"
./autogen.sh CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix="$PREFIX_INSTALL" "$LINBOX_NTLFLAG"
V="$?"; if test "x$V" != "x0";then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make dist"
make dist
V="$?"; if test "x$V" != "x0";then exit "$V"; fi

rm -rf work
mkdir -p work/builddir
mv $distarchive work
cd work
tar zxvf $distarchive
cd builddir

echo "|=== JENKINS AUTOMATED SCRIPT ===| ../$distdir/configure CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix=$PREFIX_INSTALL $LINBOX_NTLFLAG"
../$distdir/configure CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix="$PREFIX_INSTALL" "$LINBOX_NTLFLAG"

echo "|=== JENKINS AUTOMATED SCRIPT ===| make install"
make install
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make perfpublisher"
make perfpublisher

echo "|=== JENKINS AUTOMATED SCRIPT ===| make examples"
make examples
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi
(cd examples && make clean)

