# usage 1mesure_LUdivine_dense n p nb, to compute nb  LUP on a random n*n dense matrix
# over Z/pZ
# Remember to set the COMPUT_L variable to 1 in Test_LUdivine.C
mkfifo /tmp/Ad
./dense_generator $1 $1 $2> /tmp/Ad &
./Test_LUdivine $2 /tmp/Ad  $3
rm -f /tmp/Ad
