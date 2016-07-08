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
CXX=`pwd | awk -F/ '{print $(NF-2)}'`
NTL=`pwd | awk -F/ '{print $NF}'`
JENKINS_DIR=${SOURCE_DIRECTORY%%/workspace/*}
LOCAL_DIR="$JENKINS_DIR"/local
# Add path to compilers (if needed)
export PATH=$PATH:/usr/local/bin:"$LOCAL_DIR/$CXX/bin"
echo $PATH
# Add specific locations (if needed)
LD_LIBRARY_PATH="$LD_LIBRARY_PATH":/usr/local/lib:"$LOCAL_DIR/$CXX/lib"
PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:"$LOCAL_DIR/$CXX/lib/pkgconfig":"$LOCAL_DIR/$CXX"/withSSE

# Where to install linbox binaries
# Keep default for local installation.
PREFIX_INSTALL="$LOCAL_DIR/$CXX"

# Job Linbox with Ntl option flag
if [ "$NTL" == "withNTL" ]; then
  LINBOX_NTLFLAG="--with-ntl=$PREFIX_INSTALL"
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

# Particular case for Fedora23: g++=g++-5.3
#vm_name=`uname -n | cut -d"-" -f1`
#if [[ "$vm_name" == "fedora" && "$CXX" == "g++-5.3" ]]; then
#   CXX="g++"
#fi

#==================================#
# Automated installation and tests #
#==================================#

echo "|=== JENKINS AUTOMATED SCRIPT ===| ./autogen.sh CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix=$PREFIX_INSTALL --with-gmp=$GMP_PATH --with-givaro=$GIVARO_PATH $LINBOX_NTLFLAG"
./autogen.sh CXX=$CXX CXXFLAGS=$CXXFLAGS --prefix="$PREFIX_INSTALL" --with-gmp=$GMP_PATH --with-givaro=$GIVARO_PATH "$LINBOX_NTLFLAG"
V="$?"; if test "x$V" != "x0";then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make install"
make install
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

echo "|=== JENKINS AUTOMATED SCRIPT ===| make perfpublisher"
make perfpublisher

echo "|=== JENKINS AUTOMATED SCRIPT ===| make examples"
make examples
V="$?"; if test "x$V" != "x0"; then exit "$V"; fi

