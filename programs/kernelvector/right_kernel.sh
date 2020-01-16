#!/bin/bash

stopping() {
        RESULT=$(cat $1 | grep "ERROR")
    if [ "x$RESULT" != "x" ]; then
	echo "[STOP]: $2 failed, stopping everything." 
	exit 1
    fi
}

SECONDS=0
export OMP_NUM_THREADS=$7
# Run all the kernel evaluation on a matrix file

# Check and create folder if needed
if [ "$#" -lt 7 ]; then
    echo "Usage: $0 <matrixFile> <outputFolder> <fieldCharacteristic> <leftBlockSize> <rightBlockSize> <nodes> <maxThread> [step]"
    exit 1
fi
echo "command-line: $0 $*"




##### prime is bounded by 2^63 -> to go up to 128bit needs to define FLAG128=-D__USE_128bits -D__USE_SIMD=0 in Makefile.am
prime="$3"
maxprime=4294967296
if [ "$prime" -gt "$maxprime" ]; then
    echo "[WARNING]: prime is (>2^32), it might not work, if failure try to recompile with -D__USE_128bits -D__USE_SIMD=0."
    exit 1
fi

mkdir -p $2

startingstep=0
if [ "$#" -eq 8 ]; then
    startingstep=$8
    echo "[INFO]: Starting at step $startingstep"
fi

matrix_file=$1
original_matrix_file=$1
matrix_path=`dirname $1`
matrix_name=`basename $1`

rowcol=(`head -1 $matrix_file | awk '{print $1,$2}'`)
row=${rowcol[0]}
col=${rowcol[1]}
echo "[INFO]: Matrix ($row x $col) '$matrix_name' found"

##### 1 Sequence ######
if [ "$startingstep" -lt 1 ]; then
    # Sequence
    echo "[INFO]: Checking if the given matrix has not more rows than columns."   
    if [ $row -gt $col ]; then
	echo "[STOP]: The given matrix has more row than column. You should use left_kernel.sh script"
	exit 1
    fi    
    if [ $row -ne $col ]; then
	echo -n "[INFO] ($2/sequence_log.txt): extracting a square submatrix into"
	./bin/sequence_decomp_kernel $matrix_file $2 2> $2/sequence_log.txt
	#matrix_file=`echo $matrix_file | sed 's/\.[^\.]*$/\.matB\.sms/'`
	matrix_file="$2/`echo $matrix_name | sed 's/\.[^\.]*$/\.matB\.sms/'`"
	echo " $matrix_file ... DONE"
    fi	    
    echo -n "[KERN] ($2/sequence_log.txt): Generating sequence and checkpoints files ... "
    ./bin/sequence_sequence -f $matrix_file -d $2 -q $3 -s $4 -t $5 -B 2>> $2/sequence_log.txt
    echo "DONE"
fi


##### 2 Verif sequence ######
if [ "$startingstep" -lt 2 ]; then
# Verif krylov and checkpoints
# This part could be done in parallel, while sequence and checkpoints files are generated
    echo -n "[KERN] ($2/verif_log.txt): Verifying sequence and checkpoints ... "
    ./bin/verif_precompute -d $2 2> $2/verif_log_prec.txt
    NCHECKPOINTS=$(cat $2/chk_header.sdmp | cut -d' ' -f3)
    ./bin/verif_par_verif -d $2 -i $NCHECKPOINTS 2> $2/verif_log_par.txt
    cat $2/verif_log_prec.txt $2/verif_log_par.txt > $2/verif_log.txt
    echo "DONE"

    stopping $2/verif_log.txt "Sequence verification"
fi


##### 3 Polynomial ######
if [ "$startingstep" -lt 3 ]; then

# Combine files - But keeping disc space OK
     echo "[INFO]: Merging sequence files."
    cat $2/seq_header.sdmp > $2/seq.sdmp
    for f in $2/seq_0*.sdmp
      do
      cat $f >> $2/seq.sdmp
      rm $f
    done
    
# Convert from and to SDMP to compute polynomial
    echo -n "[KERN]  ($2/minpoly_log.txt) : Generating polynomial ... "
    OMP_NUM_THREADS=1 ./bin/minpoly_minpoly -r -p $3 -i $2/seq.sdmp -o $2/poly.sdmp 2>&1 2> $2/minpoly_log.txt
    ./bin/verif_verif_minpoly -d $2 2>> $2/minpoly_log.txt
    echo "DONE"

    stopping $2/minpoly_log.txt "Minpoly verification"

fi

##### 4 Evaluation : gathering files ######
if [ "$startingstep" -lt 4 ]; then
    
# Combine files - But keeping disc space OK
    echo "[INFO]: Merging checkpoints files."
    cat $2/chk_header.sdmp > $2/chk.sdmp
    for f in $2/chk_0*.sdmp
      do
      cat $f >> $2/chk.sdmp
      rm $f
    done
    
fi

##### 5 Evaluation : constructing the kernel vector  ######
if [ "$startingstep" -lt 5 ]; then
    
    # Get kernel vector
    echo -n "[KERN] ($2/evaluation_log.txt): Prepare distribution for evaluation ... "
    ./bin/evaluation_g0kernel -d $2 2> $2/evaluation_log.txt
    ./bin/verif_g0kernel -d $2 2>> $2/evaluation_log.txt
    echo "DONE"
    stopping $2/evaluation_log.txt "building g0 kernel"
    
    # Eval on each node
    echo -n "[KERN] ($2/evaluation_log.txt): Building the kernel vector of the square submatrix ... "
   ./bin/evaluation_evaluation -d $2  2>> $2/evaluation_log.txt
   echo "DONE"
  

    
    # Final checks
    echo -n "[KERN] ($2/evaluation_log.txt): Final checks ... "
    ./bin/evaluation_finalchecks -d $2 -r $4 2>> $2/evaluation_log.txt
    echo "DONE"

    stopping $2/evaluation_log.txt "building kernel vector"


    # Building the kernel of the initial matrix (not only its main square submatrix)
    echo -n "[KERN] ($2/evaluation_log.txt): Building the kernel vector from the original matrix ... "
    ./bin/evaluation_build_kernel -d $2 -f $original_matrix_file -p $3  2>> $2/evaluation_log.txt
    echo "DONE"

    stopping $2/evaluation_log.txt "building original kernel vector"
    
    RESULT=$(cat $2/evaluation_log.txt | grep "/!")
    
    echo -e "\n--> $RESULT\n"
    
    echo "Done! Saved result in $2/kernel.dv."
    
fi


duration=$SECONDS
echo -e "\nComputation time: $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."




