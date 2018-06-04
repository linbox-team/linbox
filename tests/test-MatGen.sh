#!/bin/bash
OPT_A=2
OPT_B=3
OPT_C=1
OPT_D=1
while getopts ":d:n:i:p:" OPTION
do
	case $OPTION in
	d ) OPT_A=$OPTARG ;;
	n ) OPT_B=$OPTARG ;;
	i ) OPT_C=$OPTARG ;;
	p ) OPT_D=$OPTARG ;;
	esac
done
shift $(($OPTIND - 1))

if ((OPT_D<2)); then
echo "Running a set of tests for matrix " $OPT_B "by" $OPT_B " with" $OPT_A "digits without parallelization"
    eval "./test-MatGen -d $OPT_A -n $OPT_B"
elif type $mpirun; then
         COUNTER=0
	 DIM=250
echo "Running a set of tests for matrix from 250 by 250 to" $OPT_B "by" $OPT_B " with" $OPT_A "digits with parallelization" 
echo "Each test will be done for" $OPT_C "iteration(s) with" $OPT_D "process(es) and a step size of 50"
         until [ $COUNTER -gt $OPT_C ]; do
		until [ $DIM -gt $OPT_B ]; do
	if ! (eval "mpirun -np $OPT_D test-MatGen -d $OPT_A -n $DIM"); then 
 echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
 echo "         mpirun -np $OPT_D test-MatGen -d $OPT_A -n $DIM failed            "; 
 echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	  exit 1; 
	else
 	  eval "mpirun -np $OPT_D test-MatGen -d $OPT_A -n $DIM"
	fi
        let COUNTER+=1
        let DIM+=50
		done
         done   
else
"Running a set of tests for matrix " $OPT_B "by" $OPT_B " with" $OPT_A "digits without parallelization as mpirun is not available"
	eval "./test-MatGen -d $OPT_A -n $OPT_B"
fi
