#!/usr/bin/python

from random import *;

# Test matrix 0: A random dense 4x4 matrix

f = open ('test-0.matrix', 'w');

f.write ("4 4 M\n");

for i in range (4):
    for j in range (4):
        f.write (`i + 1` + " " + `j + 1` + " " + `randint (-10, 10)` + "\n");

f.write ("0 0 0\n");

f.close ();

# Test vector 1: A random dense vector of length 200

f = open ('test-1.vector', 'w');

for i in range(200):
    f.write (`randint(-10,10)` + "\n");

f.close ();

# Test vector 2: A random sparse vector of length 200

data=[];

f = open ('test-2.vector', 'w');

for i in range (200):
    data.append (0);

for i in range (40):
    data[randint(0,199)] = randint(-10,10);

for i in data:
    f.write (`i` + "\n");

f.close ();

# Test vectors 3-1, 3-2: Random sparse vectors differing by only one nonzero entry

data=[];

f1 = open ('test-3-1.vector', 'w');
f2 = open ('test-3-2.vector', 'w');

for i in range (200):
    data.append (0);

for i in range (40):
    data[randint(0,199)] = randint(-10,10);

for i in data:
    f1.write (`i` + "\n");

newdata = randint (-10,10);
col = randint (0, 199);
while data[col] == newdata or data[col] == 0:
    col = randint (0, 199);

data[col] = newdata;

for i in data:
    f2.write (`i` + "\n");

f1.close ();
f2.close ();

# Test vectors 4-1, 4-2: Random sparse vectors; 4-2 is the same as 4-1 except that is it missing an entry

data=[];

f1 = open ('test-4-1.vector', 'w');
f2 = open ('test-4-2.vector', 'w');

for i in range (200):
    data.append (0);

for i in range (39):
    data[randint(0,199)] = randint(-10,10);

for i in data:
    f2.write (`i` + "\n");

col = randint (0, 199);
while data[col] != 0:
    col = randint (0, 199);

data[col] = randint (-10, 10);

for i in data:
    f1.write (`i` + "\n");

f1.close ();
f2.close ();

# Test vectors 5-1, 5-2: Random orthogonal sparse vectors

data1=[];
data2=[];

f1 = open ('test-5-1.vector', 'w');
f2 = open ('test-5-2.vector', 'w');

for i in range (200):
    data1.append (0);
    data2.append (0);

for i in range (49):
    data1[randint(0,198)] = randint(-10,10);
    data2[randint(0,198)] = randint(-10,10);

sum = 0;
for i in range (200):
    sum = sum + data1[i] * data2[i];

while data1[199] == 0:
    data1[199] = randint(-10,10);

data2[199] = -sum / data1[199];

for i in range (200):
    f1.write (`data1[i]` + "\n");
    f2.write (`data2[i]` + "\n");

f1.close ();
f2.close ();

# Test vector 5: Random sparse vector whose entries are all 1, 0, and -1

data=[];

f = open ('test-5.vector', 'w');

for i in range (200):
    data.append (0);

for i in range (40):
    data[randint(0,199)] = randint(-1,1);

for i in data:
    f.write (`i` + "\n");

f.close ();

# Test matrix 1: A random dense 200x200 matrix

f = open ('test-1.matrix', 'w');

f.write ("200 200 M\n");

for i in range (200):
    for j in range (200):
        f.write (`i + 1` + " " + `j + 1` + " " + `randint (-10, 10)` + "\n");

f.write ("0 0 0\n");

f.close ();

# Test matrix 2: A random sparse 200x200 matrix

f = open ('test-2.matrix', 'w');

f.write ("200 200 M\n");

for i in range (200):
    data = [];
    
    for j in range (200):
        data.append (0);

    for j in range (40):
        data[randint(0,199)] = randint(-10,10);

    for j in range (200):
        if data[j] != 0:
            f.write (`i` + " " + `j` + " " + `data[j]` + "\n");

f.write ("0 0 0\n");

f.close ();

# Test matrix 3: A random sparse 200x200 1,0,-1 -matrix

f = open ('test-3.matrix', 'w');

f.write ("200 200 M\n");

for i in range (200):
    data = [];
    
    for j in range (200):
        data.append (0);

    for j in range (40):
        data[randint(0,199)] = randint(-1,1);

    for j in range (200):
        if data[j] != 0:
            f.write (`i` + " " + `j` + " " + `data[j]` + "\n");

f.write ("0 0 0\n");

f.close ();

