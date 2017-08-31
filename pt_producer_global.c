#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"
#include "pt_structs.h"

int text_read_state(FILE *fp, pt_atoms *atoms_array, char allocate_atoms_array) 
{
	int  i, len;
	char *str;
	int  num_atoms, num_types;

	len = MAX_LINE_LEN;
	if ((str=(char*)malloc(sizeof(char)*len))==NULL) return 1;

	/* Skip the header string */
	if (fgets(atoms_array->header_str,len,fp)==NULL) return 1;

	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 

	/* Read number of atoms */
	if (fgets(str,len,fp)==NULL) return 1;
	sscanf(str,"%d",&num_atoms);
	/* Read number of types */
	if (fgets(str,len,fp)==NULL) return 1;
	sscanf(str,"%d",&num_types);

	// Allocate space for atoms
	if (allocate_atoms_array) {
		if (alloc_atoms_array(atoms_array,num_atoms,num_types)!=0) return 1;
	}

	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 
	
	/* Read simulation cell dimensions */
	if (fgets(str,len,fp)==NULL) return 1;
	sscanf(str,"%lg %lg",&(atoms_array->cell_dims[0][0]),&(atoms_array->cell_dims[0][1]));
	if (fgets(str,len,fp)==NULL) return 1;
	sscanf(str,"%lg %lg",&(atoms_array->cell_dims[1][0]),&(atoms_array->cell_dims[1][1]));
	if (fgets(str,len,fp)==NULL) return 1;
	sscanf(str,"%lg %lg",&(atoms_array->cell_dims[2][0]),&(atoms_array->cell_dims[2][1]));

	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 
	/* Skip "Masses" text line */
	if (fgets(str,len,fp)==NULL) return 1; 
	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 

	/* Read masses of atom types */
	for (i=0;i<atoms_array->num_types;i++) {
		if (fgets(str,len,fp)==NULL) return 1;
		sscanf(str,"%d %lg",&(atoms_array->type_id[i]),&(atoms_array->masses[i]));
	}

	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 
	/* Skip "Atoms # atomic" text line */
	if (fgets(str,len,fp)==NULL) return 1; 
	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 

	// Read atoms
	for (i=0;i<atoms_array->num_atoms;i++) {
		if (fgets(str,len,fp)==NULL) return 1;
		sscanf(str,"%d %d %lg %lg %lg %d %d %d",
				&(atoms_array->atom_id[i]),
				&(atoms_array->atom_type[i]),
				&(atoms_array->px[i]),
				&(atoms_array->py[i]),
				&(atoms_array->pz[i]),
				&(atoms_array->imx[i]),
				&(atoms_array->imy[i]),
				&(atoms_array->imz[i]));
	}

	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 
	/* Skip "Velocities" text line */
	if (fgets(str,len,fp)==NULL) return 1; 
	/* Skip empty line */
	if (fgets(str,len,fp)==NULL) return 1; 

	// Read atom velocities
	for (i=0;i<atoms_array->num_atoms;i++) {
		if (fgets(str,len,fp)==NULL) return 1;
		sscanf(str,"%d %lg %lg %lg",
				&(atoms_array->atom_vid[i]),
				&(atoms_array->vx[i]),
				&(atoms_array->vy[i]),
				&(atoms_array->vz[i]));
	}

	free(str);
	
	return 0;
}

#define GROUP_NAME "pt_exaalt_global"
int adios_declare(int64_t *gh, const char *transport_method, const char* opts)
{
    	adios_declare_group (gh, GROUP_NAME, NULL, adios_stat_default);
    	adios_select_method (*gh, transport_method, opts, "");

	adios_define_var (*gh, "num_procs", "", adios_integer, "", "", "");
	adios_define_var (*gh, "proc_no", "", adios_integer, "", "", "");
    	adios_define_var (*gh, "num_atoms", "", adios_integer, "", "", "");
    	adios_define_var (*gh, "num_types", "", adios_integer, "", "", "");

    	adios_define_var (*gh, "state_id", "", adios_byte, "1,64", "num_procs,64", "proc_no,0");
	adios_define_var (*gh, "cell_dims", "", adios_double, "1,6", "num_procs,6", "proc_no,0");

	adios_define_var (*gh, "type_id", "", adios_integer, "1,num_types", "num_procs,num_types", "proc_no,0");
	adios_define_var (*gh, "masses", "", adios_double, "1,num_types", "num_procs,num_types", "proc_no,0");

	adios_define_var (*gh, "atom_id", "", adios_integer, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "atom_type", "", adios_integer, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "px", "", adios_double, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "py", "", adios_double, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "pz", "", adios_double, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "imx", "", adios_integer, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "imy", "", adios_integer, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "imz", "", adios_integer, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	
	adios_define_var (*gh, "atom_vid", "", adios_integer, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "vx", "", adios_double, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "vy", "", adios_double, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");
	adios_define_var (*gh, "vz", "", adios_double, "1,num_atoms", "num_procs,num_atoms", "proc_no,0");

	return 0;	
}

int adios_write_state(char *adios_file, char *fmode, MPI_Comm comm, int num_procs, int proc_no, pt_atoms *atoms_array)
{
	int64_t fh;

    	int err = adios_open (&fh, GROUP_NAME, adios_file, fmode, comm);
	if (err != MPI_SUCCESS) {
		fprintf(stderr,"Error opening file: %s\n",adios_file);
		return 1;
	}

	uint64_t groupsize;
	uint64_t groupTotalSize;

	groupsize  = 4*sizeof(int);      /* num_procs, proc_no, num_atoms, num_types */
	groupsize += 64*sizeof(char);    /* state_id */
	groupsize += 3*2*sizeof(double); /* cell_dims */
	groupsize += atoms_array->num_types*sizeof(int);    /* type_id */
	groupsize += atoms_array->num_types*sizeof(double); /* masses */

	groupsize += atoms_array->num_atoms*sizeof(int)*2;    /* atom_id, atom_type */
	groupsize += atoms_array->num_atoms*sizeof(double)*3; /* px,py,pz */
	groupsize += atoms_array->num_atoms*sizeof(int)*3;    /* imx, imy, imz */
	groupsize += atoms_array->num_atoms*sizeof(int);      /* atom_vid */
	groupsize += atoms_array->num_atoms*sizeof(double)*3; /* vx, vy, vz */

    	adios_group_size (fh, groupsize, &groupTotalSize);

	adios_write (fh, "num_procs", &num_procs);
	adios_write (fh, "proc_no",   &proc_no);
        adios_write (fh, "num_atoms", &atoms_array->num_atoms);
        adios_write (fh, "num_types", &atoms_array->num_types);

        adios_write (fh, "state_id",  atoms_array->state_id);
	adios_write (fh, "cell_dims", atoms_array->cell_dims);
	
	adios_write (fh, "type_id", atoms_array->type_id);
	adios_write (fh, "masses",  atoms_array->masses);

	adios_write (fh, "atom_id", atoms_array->atom_id);
	adios_write (fh, "atom_type", atoms_array->atom_type);
	adios_write (fh, "px", atoms_array->px);
	adios_write (fh, "py", atoms_array->py);
	adios_write (fh, "pz", atoms_array->pz);
	adios_write (fh, "imx", atoms_array->imx);
	adios_write (fh, "imy", atoms_array->imy);
	adios_write (fh, "imz", atoms_array->imz);

	adios_write (fh, "atom_vid", atoms_array->atom_vid);
	adios_write (fh, "vx", atoms_array->vx);
	adios_write (fh, "vy", atoms_array->vy);
	adios_write (fh, "vz", atoms_array->vz);
	
    	adios_close(fh);

	return 0;
}

char **create_filename_array(FILE *fpi)
{
	int i;

	// Read in the entire state list file
	fseek(fpi, 0L, SEEK_END);
	long fsize = ftell(fpi);
	fseek(fpi, 0L, SEEK_SET);
	char *file_buffer = (char*)malloc(sizeof(char)*fsize); 
	if (file_buffer==NULL) { 
		fprintf(stderr, "Cannot allocate array to read file.\n");
		return NULL;
	}
	fread(file_buffer, fsize, 1, fpi);

	// find and count all the lines in the file
	int line_cnt   = 0;
	char *file_str = file_buffer;
	for (i=0;i<fsize;i++) 
		if (file_str[i]=='\n') line_cnt++;
	char **filename_ptr = (char**)malloc(sizeof(char*)*(line_cnt+1));
	if (filename_ptr==NULL) {
		fprintf(stderr, "Cannot allocate array to store file names.\n");
		return NULL;
	} 

	// Create a list of filenames
	file_str = file_buffer;
	filename_ptr[0] = file_buffer;
	line_cnt = 1;
	for (i=0;i<(fsize-1);i++) {
		if (file_str[i]=='\n') {
			file_str[i] = '\0';
			filename_ptr[line_cnt] = &file_str[i+1];
			line_cnt++;
		}
	}
	return filename_ptr;
}

void free_filename_array(char **filename_ptr)
{
	free(filename_ptr[0]);
	free(filename_ptr);
}

#define MAX_NUM_STATES 10240
int main(int argc, char *argv[])
{
	pt_atoms atoms_array;
	MPI_Comm comm = MPI_COMM_WORLD;
	int64_t gh;
	int  comm_rank, comm_size;
	int  i;
	// char filename[256];
	FILE *fpi, *fpp;

    	MPI_Init (&argc, &argv);
    	MPI_Comm_rank (comm, &comm_rank);
	MPI_Comm_size (comm, &comm_size);

	if (argc!=6) {
		printf("new Usage: %s <input file list> <num of states> <bp file> <transport method> <transport opts>\n",
			argv[0]);
		return 1;
	}

	char *input_file_list = argv[1];
	int  num_states  = (int)atoi(argv[2]);
	char *bp_file = argv[3];
	char *transport_method = argv[4];
	char *transport_opts = argv[5];

	/* Create adios structure */
	adios_init_noxml(comm);
	adios_declare(&gh,transport_method,transport_opts);

	if ((fpi=fopen(input_file_list,"r"))==NULL) return 1;
	char **filename_array;
	filename_array = create_filename_array(fpi);
	if (filename_array==NULL) return 1;	

	MPI_Barrier(comm);
	int first_time = 1;
	int state_cnt  = 0;
	double io_time = 0.0, io_time_start = 0.0, io_time_end = 0.0;
	for (i=0;i<num_states;i++) {
		if (comm_rank==(i%comm_size)) {
			if ((fpp=fopen(filename_array[state_cnt],"r"))==NULL) { 
				printf("Cannot open file... %s\n",filename_array[i]);
				return 1;
			}
			if (first_time) {	
				if (text_read_state(fpp,&atoms_array,1)!=0) {
					printf("Read Error.....\n");
					return 1;
				}

				io_time_start = MPI_Wtime();
				adios_write_state(bp_file,(char*)"w",comm,comm_size,comm_rank,&atoms_array);
				io_time_end   = MPI_Wtime();
				first_time =0;
			} else {
				if (text_read_state(fpp,&atoms_array,0)!=0) {
					printf("Read Error.....\n");
					return 1;
				}
				io_time_start = MPI_Wtime();
				adios_write_state(bp_file,(char*)"a",comm,comm_size,comm_rank,&atoms_array);
				io_time_end   = MPI_Wtime();
			}
			fclose(fpp);
			io_time += (io_time_end-io_time_start);	
		} 
		state_cnt++;
<<<<<<< HEAD
		if (state_cnt>=MAX_NUM_STATES) state_cnt = 0;
=======
		if (state_cnt>=MAX_NUM_STATES) {
			state_cnt = 0;
			fseek(fpi,0,SEEK_SET);
		}
>>>>>>> 4bc0151c339284bcf5fd5d33c0af4ddf4f2be3cc
	}
	fclose(fpi);
	free_atoms_array(&atoms_array); 
	free_filename_array(filename_array);

	MPI_Barrier(comm);
	printf("Rank: %d io_time: %lf\n",comm_rank,io_time); fflush(stdout);
	MPI_Barrier(comm);
    	adios_finalize (comm_rank);
    	MPI_Finalize ();

	return 0;
}
