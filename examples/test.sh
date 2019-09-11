#!/bin/bash 

# written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
# part of LinBox, see COPYING
# SED="sed"
# case "`uname`" in
#   Darwin*) SED="gsed" ;;
# esac

set -o nounset                              # Treat unset variables as an error

fail() {
	echo "fail"
}

success() {
	echo "ok"
}


pass="true"

echo -n "check rank ... "
result=`./rank  data/test.matrix 7`
[ "$result" -eq "9" ] && success || { fail ;  pass="false" ; }

echo -n "check rank ... "
result=`./rank  data/test.matrix`
[ "$result" -eq "10" ] && success || { fail ; pass="false" ; }


[ "$pass" == "true" ] && { success ; exit 0  ;} || { fail ; exit 1 ;};

