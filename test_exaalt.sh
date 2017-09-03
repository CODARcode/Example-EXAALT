#!/bin/bash

aprun -n 64 ./stage_write/stage_write  output.bp staged.bp FLEXPATH "" POSIX "have_metadata_file=0" &
aprun -n 64 ./pt_producer_global states_list.txt 8192 output.bp FLEXPATH "" "none"
