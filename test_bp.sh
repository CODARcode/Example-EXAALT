#!/bin/bash

aprun -n 64 ./pt_producer_global states_list.txt 8192 output.bp POSIX "have_metadata_file=0" "zlib:9"
