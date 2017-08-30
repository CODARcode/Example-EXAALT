#!/bin/bash

mpirun -np 8 ../stage_write/stage_write output.bp staged.bp FLEXPATH "" MPI "" &
mpirun -np 8 ../pt_producer_global states_local.txt 1024 output.bp FLEXPATH "" 
