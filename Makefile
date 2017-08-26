EXECUTABLES=pt_reader pt_reader_global pt_producer_global

ADIOS_DIR=`adios_config -d` 
ADIOS_INCS=`adios_config -c`
ADIOS_LIBS=`adios_config -l` 

C_FLAGS=-O3

CXX=mpicxx
OPENMP=-DNOOMP

# Makefile rules
all: $(EXECUTABLES)

pt_reader: pt_reader.c 
	$(CXX) $(OPENMP) $(C_FLAGS) $(ADIOS_INCS) -o pt_reader pt_reader.c $(ADIOS_LIBS) 

pt_producer_global: pt_producer_global.c 
	$(CXX) $(OPENMP) $(C_FLAGS) $(ADIOS_INCS) -o pt_producer_global pt_producer_global.c $(ADIOS_LIBS) 

pt_reader_global: pt_reader_global.c 
	$(CXX) $(OPENMP) $(C_FLAGS) $(ADIOS_INCS) -o pt_reader_global pt_reader_global.c $(ADIOS_LIBS) 

clean:
	rm -f $(EXECUTABLES)
