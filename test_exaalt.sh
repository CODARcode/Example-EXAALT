#!/bin/bash

mpirun -np 8 /home/tkurc/codar/Example-Heat_Transfer/stage_write/stage_write  output.bp staged.bp FLEXPATH "" MPI "" &
mpirun -np 8 ./pt_producer_global states_list.txt 128 output.bp FLEXPATH "" 
