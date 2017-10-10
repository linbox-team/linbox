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
rank_cmd="Rank\sis\s"
./rank  data/test.matrix 7 > linbox-tmp.data
result=`cat linbox-tmp.data | grep ${rank_cmd} | $SED 's/'"$rank_cmd"'\([0-9]*\).*/\1/'`
[ "$result" -eq "9" ] && success || { fail ;  pass="false" ; }

echo -n "check rank ... "
./rank  data/test.matrix > linbox-tmp.data
result=`cat linbox-tmp.data | grep ${rank_cmd} | $SED 's/'"$rank_cmd"'\([0-9]*\).*/\1/'`
[ "$result" -eq "10" ] && success || { fail ; pass="false" ; }


[ "$pass" == "true" ] && { success ; exit 0  ;} || { fail ; exit 1 ;};

