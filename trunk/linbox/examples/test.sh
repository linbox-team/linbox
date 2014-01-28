#!/bin/bash - 

# written by bb <bbboyer@ncsu.edu> , part of LinBox, see COPYING

set -o nounset                              # Treat unset variables as an error

fail() {
	echo "fail"
}

success() {
	echo "ok"
}


pass=true

echo -n "check rank ... "
rank_cmd="Rank\sis\s"
./rank  data/test.matrix 7 > linbox-tmp.data
result=`cat linbox-tmp.data | grep ${rank_cmd} | sed 's/'"$rank_cmd"'\([0-9]*\).*/\1/'`
[ "${result}" -eq "9" ] && success || { fail ;  pass=false ; }

echo -n "check rank ... "
./rank  data/test.matrix > linbox-tmp.data
result=`cat linbox-tmp.data | grep ${rank_cmd} | sed 's/'"$rank_cmd"'\([0-9]*\).*/\1/'`
[ "${result}" -eq "10" ] && success || { fail ; pass=false ; }


[ "${pass}" -eq "true" ] && { echo -e "\nSUCCESS" ; exit 0  ;} || { echo -e "\nFAIL" ; exit 1 ;};

