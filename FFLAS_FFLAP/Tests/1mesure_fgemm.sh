./dense_generator $1 $1 >/tmp/A &
./dense_generator $1 $1 >/tmp/B &
./Test_fgemm 101 /tmp/A /tmp/B 100 1