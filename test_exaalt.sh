#!/bin/bash

mpirun -np 2 /home/tkurc/codar/Example-Heat_Transfer/stage_write/stage_write  output.bp staged.bp FLEXPATH "" MPI "" &
mpirun -np 4 ./pt_reader_global states_list.txt output.bp FLEXPATH "" output.txt 1
