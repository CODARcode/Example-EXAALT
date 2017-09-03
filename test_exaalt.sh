#!/bin/bash

aprun -n 16 ./stage_write/stage_write  output.bp staged.bp FLEXPATH "" POSIX "have_metadata_file=0" 16 1 &
aprun -n 16 ./pt_producer_global states_list.txt 32 output.bp FLEXPATH "" 
