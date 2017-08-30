# EXECUTABLES=pt_producer_global compute_stats pt_test
EXECUTABLES=pt_producer_global

ADIOS_DIR=`adios_config -d` 
ADIOS_INCS=`adios_config -c`
ADIOS_LIBS=`adios_config -l` 

C_FLAGS=-O3 -Wall

CXX=mpicxx
OPENMP=-DNOOMP

# Makefile rules
all: $(EXECUTABLES)

pt_producer_global: pt_producer_global.c 
	$(CXX) $(OPENMP) $(C_FLAGS) $(ADIOS_INCS) -o pt_producer_global pt_producer_global.c $(ADIOS_LIBS) 

pt_test: pt_test.c 
	$(CXX) $(OPENMP) $(C_FLAGS) $(ADIOS_INCS) -o pt_test pt_test.c $(ADIOS_LIBS) 

compute_stats: compute_stats.c 
	$(CXX) $(OPENMP) $(C_FLAGS) $(ADIOS_INCS) -o compute_stats compute_stats.c $(ADIOS_LIBS) -lm

clean:
	rm -f $(EXECUTABLES)

