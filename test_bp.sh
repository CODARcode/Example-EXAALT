#!/bin/bash

aprun -n 64 ./pt_producer_global states_list.txt 1024000 output.bp POSIX "have_metadata_file=0" "zlib:9"
