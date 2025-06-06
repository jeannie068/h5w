# HW5 Placement Legalization
Implement an existing algorithm, published in the ISPD-08 paper entitled “Abacus: fast legalization of standard cell circuits with minimal movement” by Spindler, Schlichtmann and Johannes, to legalize a given global placement result with minimal total displacement (measured by Euclidean distance).

## How to Compile
In `HW5/src/` directory, enter the following command:
```
$ make
```
An executable file `hw5` will be generated in `HW5/bin/`.

If you want to remove it, please enter the following command:
```
$ make clean
```

## How to Run
In `HW5/bin/` directory, enter the following command:
Format: 
```
$ ./hw4 <input file> <output file>
```

E.g.,
```
$ ./hw5 ../testcase/public1.txt ../output/public1.out
```


