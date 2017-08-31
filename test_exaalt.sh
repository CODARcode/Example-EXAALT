#!/bin/bash

aprun -n 8 ./stage_write/stage_write  output.bp staged.bp FLEXPATH "" MPI "" &
aprun -n 16 ./pt_producer_global states_list.txt 6000 output.bp FLEXPATH "" 
