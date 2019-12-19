#!/bin/bash

# Run all the kernel evaluation on a matrix file

# Check and create folder if needed
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <matrixFile> <outputFolder> <fieldCharacteristic> <leftBlockSize> <rightBlockSize> <nodes> [step]"
    exit 1
fi
echo "command-line: $0 $*"

mkdir -p $2

startingstep=0
if [ "$#" -eq 7 ]; then
    startingstep=$7
    echo "Starting at step $startingstep"
fi

matrix_file=$1
rowcol=(`head -1 $matrix_file | awk '{print $1,$2}'`)
row=${rowcol[0]}
col=${rowcol[1]}
echo "Matrix ($row x $col) found"

new_matrix_file=$1
##### 1 Sequence ######
if [ "$startingstep" -lt 1 ]; then
    # Sequence
    echo "Checking if the given matrix has not more columns rows than row..."   
    if [ $row -lt $col ]; then
	echo "The given matrix has more column  than rows. You should use right_kernel.sh script"
	exit 1
    fi

    # compute the transpose matrix and call right kernel script
    new_matrix_file=`echo $matrix_file | sed 's/\.[^\.]*$/-trans\.sms/'`
    echo "Transposing the matrix into file $new_matrix_file ($col x $row)"
    ./bin/sequence_transpose $matrix_file > $new_matrix_file 2> $2/transpose_log.txt
fi

./right_kernel.sh $new_matrix_file $2 $3 $4 $5 $6 $startingstep

