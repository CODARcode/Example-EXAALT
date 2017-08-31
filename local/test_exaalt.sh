#!/bin/bash

mpirun -np 8 ../stage_write/stage_write output.bp staged.bp FLEXPATH "" MPI "" &
mpirun -np 8 ../pt_producer_global states_local.txt 128 output.bp FLEXPATH "" 
