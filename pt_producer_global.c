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

#define HEADER_STRING "LAMMPS data file via write_data, version 30 Nov 2016, timestep = 0\n"
#define MAX_LINE_LEN  1024
#define SC_NUM_ATOMS  147
#define SC_NUM_TYPES  1

typedef struct _pt_atoms {
	char    header_str[128];     /* header string for txt files     */
	char    state_id[64];        /* state/configuration id          */
	double  cell_dims[3][2];     /* minx,maxx, miny,maxy, minz,maxz */

	/* atom type and mass */
	int    num_types;
	int    *type_id;
	double *masses;	

	/* len: num_atoms */	
	int    num_atoms;
	int    *atom_id;
	int    *atom_type; 
	double *px, *py, *pz;
	int    *imx, *imy, *imz;
	int    *atom_vid;
	double *vx, *vy, *vz;
} pt_atoms;

typedef struct _atomic_stats {
	int    cnt;
	double mean_px;
	double mean_py;
	double mean_pz;
	double mean_vx;
	double mean_vy;
	double mean_vz;
	double std_px;
	double std_py;
	double std_pz;
	double std_vx;
	double std_vy;
	double std_vz;
} atomic_stats;

void set_stats(atomic_stats *stats_atoms)
{
	stats_atoms->cnt = 1470000;
	stats_atoms->mean_px = 15.003000; stats_atoms->std_px = 8.931469;
	stats_atoms->mean_py = 14.996714; stats_atoms->std_py = 8.953612;
	stats_atoms->mean_pz = 14.995469; stats_atoms->std_pz = 8.948379;
	stats_atoms->mean_vx = -0.000738; stats_atoms->std_vx = 1.600756;
	stats_atoms->mean_vy = -0.000068; stats_atoms->std_vy = 1.599595;
	stats_atoms->mean_vz =  0.000286; stats_atoms->std_vz = 1.600011;
}

#define RANDOM_VALUE (((double)rand())/RAND_MAX)
int generate_random_atoms(pt_atoms *atoms_array, atomic_stats *stats_atoms)
{
	int i;
	for (i=0;i<atoms_array->num_atoms;i++) {
		atoms_array->px[i] = (RANDOM_VALUE*(2*stats_atoms->std_px))+(stats_atoms->mean_px-stats_atoms->std_px);
		atoms_array->py[i] = (RANDOM_VALUE*(2*stats_atoms->std_py))+(stats_atoms->mean_py-stats_atoms->std_py);
		atoms_array->pz[i] = (RANDOM_VALUE*(2*stats_atoms->std_pz))+(stats_atoms->mean_pz-stats_atoms->std_pz);
		atoms_array->vx[i] = (RANDOM_VALUE*(2*stats_atoms->std_vx))+(stats_atoms->mean_vx-stats_atoms->std_vx);
		atoms_array->vy[i] = (RANDOM_VALUE*(2*stats_atoms->std_vy))+(stats_atoms->mean_vy-stats_atoms->std_vy);
		atoms_array->vz[i] = (RANDOM_VALUE*(2*stats_atoms->std_vz))+(stats_atoms->mean_vz-stats_atoms->std_vz);

		atoms_array->atom_type[i] = 1;
		atoms_array->atom_id[i]   = i+1;
		atoms_array->atom_vid[i]  = i+1;
	}
	return 0;
}

void init_atoms_array(pt_atoms *atoms_array)
{
	atoms_array->num_atoms = 0;
	atoms_array->num_types = 0;

	atoms_array->type_id = NULL;
	atoms_array->masses = NULL;	

	atoms_array->atom_id = NULL;
	atoms_array->atom_type = NULL; 
	atoms_array->px = NULL;
	atoms_array->py = NULL;
	atoms_array->pz = NULL;
	atoms_array->imx = NULL;
	atoms_array->imy = NULL;
	atoms_array->imz = NULL;
	atoms_array->atom_vid = NULL;
	atoms_array->vx = NULL;
	atoms_array->vy = NULL;
	atoms_array->vz = NULL;
}

void free_atoms_array(pt_atoms *atoms_array)
{
	if (atoms_array->type_id) free(atoms_array->type_id);
	if (atoms_array->masses) free(atoms_array->masses);	

	if (atoms_array->atom_id) free(atoms_array->atom_id);
	if (atoms_array->atom_type) free(atoms_array->atom_type); 
	if (atoms_array->px) free(atoms_array->px);
	if (atoms_array->py) free(atoms_array->py);
	if (atoms_array->pz) free(atoms_array->pz);
	if (atoms_array->imx) free(atoms_array->imx);
	if (atoms_array->imy) free(atoms_array->imy);
	if (atoms_array->imz) free(atoms_array->imz);
	if (atoms_array->atom_vid) free(atoms_array->atom_vid);
	if (atoms_array->vx) free(atoms_array->vx);
	if (atoms_array->vy) free(atoms_array->vy);
	if (atoms_array->vz) free(atoms_array->vz);

	init_atoms_array(atoms_array);
}

int alloc_atoms_array(pt_atoms *atoms_array, int num_atoms, int num_types) 
{
	atoms_array->num_atoms = num_atoms;
	atoms_array->num_types = num_types;

	// Allocate space for atoms
	atoms_array->atom_id = (int*)malloc(sizeof(int)*atoms_array->num_atoms);
	if (atoms_array->atom_id==NULL) return 1;
	atoms_array->atom_type = (int*)malloc(sizeof(int)*atoms_array->num_atoms);
	if (atoms_array->atom_type==NULL) return 1;
	atoms_array->px = (double*)malloc(sizeof(double)*atoms_array->num_atoms);
	if (atoms_array->px==NULL) return 1;
	atoms_array->py = (double*)malloc(sizeof(double)*atoms_array->num_atoms);
	if (atoms_array->py==NULL) return 1;
	atoms_array->pz = (double*)malloc(sizeof(double)*atoms_array->num_atoms);
	if (atoms_array->pz==NULL) return 1;
	atoms_array->imx = (int*)malloc(sizeof(int)*atoms_array->num_atoms);
	if (atoms_array->imx==NULL) return 1;
	atoms_array->imy = (int*)malloc(sizeof(int)*atoms_array->num_atoms);
	if (atoms_array->imy==NULL) return 1;
	atoms_array->imz = (int*)malloc(sizeof(int)*atoms_array->num_atoms);
	if (atoms_array->imz==NULL) return 1;
	atoms_array->atom_vid = (int*)malloc(sizeof(int)*atoms_array->num_atoms);
	if (atoms_array->atom_vid==NULL) return 1;
	atoms_array->vx = (double*)malloc(sizeof(double)*atoms_array->num_atoms);
	if (atoms_array->vx==NULL) return 1;
	atoms_array->vy = (double*)malloc(sizeof(double)*atoms_array->num_atoms);
	if (atoms_array->vy==NULL) return 1;
	atoms_array->vz = (double*)malloc(sizeof(double)*atoms_array->num_atoms);
	if (atoms_array->vz==NULL) return 1;

	atoms_array->masses  = (double*)malloc(sizeof(double)*atoms_array->num_types);
	if (atoms_array->masses==NULL) return 1;
	atoms_array->type_id = (int*)malloc(sizeof(int)*atoms_array->num_types);
	if (atoms_array->type_id==NULL) return 1;
	
	return 0;	
}

int text_read_state(FILE *fp, pt_atoms *atoms_array, char allocate_atoms_array) 
{
	int  i, len, str_len;
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

int main(int argc, char *argv[])
{
	pt_atoms atoms_array;
	pt_atoms atoms_out;
	atomic_stats stats_atoms;
	MPI_Comm comm = MPI_COMM_WORLD;
	int64_t gh;
	int  comm_rank, comm_size;
	int  i;
	char filename[256];
	FILE *fpi, *fpp;

    	MPI_Init (&argc, &argv);
    	MPI_Comm_rank (comm, &comm_rank);
	MPI_Comm_size (comm, &comm_size);

	if (argc!=7) {
		printf("Usage: %s <input file list> <num of states> <num of randoms> <bp file> <transport method> <transport opts>\n",
			argv[0]);
		return 1;
	}

	char *input_file_list = argv[1];
	int  num_states  = (int)atoi(argv[2]);
	int  num_randoms = (int)atoi(argv[3]);
	char *bp_file = argv[4];
	char *transport_method = argv[5];
	char *transport_opts = argv[6];

	/* Create adios structure */
	adios_init_noxml(comm);
	adios_declare(&gh,transport_method,transport_opts);

	if ((fpi=fopen(argv[1],"r"))==NULL) return 1;

	int first_time = 1;
	for (i=0;i<num_states;i++) {
		fscanf(fpi,"%s",filename);
		if (comm_rank==(i%comm_size)) {
			if ((fpp=fopen(filename,"r"))==NULL) return 1;

			if (first_time) {	
				if (text_read_state(fpp,&atoms_array,1)!=0) {
					printf("Read Error.....\n");
					return 1;
				}
				adios_write_state(bp_file,(char*)"w",comm,comm_size,comm_rank,&atoms_array);
				first_time =0;
			} else {
				if (text_read_state(fpp,&atoms_array,0)!=0) {
					printf("Read Error.....\n");
					return 1;
				}
				adios_write_state(bp_file,(char*)"a",comm,comm_size,comm_rank,&atoms_array);
			}
		
			fclose(fpp);	
			/* free_atoms_array(&atoms_array); */
		} 
	}
	fclose(fpi);

	srand(0);
	set_stats(&stats_atoms);
	for (i=0;i<num_randoms;i++) {
		if (comm_rank==(i%comm_size)) {
			if (first_time) {	
				if (alloc_atoms_array(&atoms_array, SC_NUM_ATOMS, SC_NUM_TYPES)!=0) { 
					printf("Cannot allocate array.\n");
					return 1;
				}
				atoms_array.masses[0]  = 191.0;
				atoms_array.type_id[0] = 1;
				if (generate_random_atoms(&atoms_array,&stats_atoms)!=0) {
					printf("Random atom generation Error.....\n");
					return 1;
				}
				adios_write_state(bp_file,(char*)"w",comm,comm_size,comm_rank,&atoms_array);
				first_time =0;
			} else {
				if (generate_random_atoms(&atoms_array,&stats_atoms)!=0) {
					printf("Random atom generation Error.....\n");
					return 1;
				}
				adios_write_state(bp_file,(char*)"a",comm,comm_size,comm_rank,&atoms_array);
			}
			/* free_atoms_array(&atoms_array); */
		} 
	}

	MPI_Barrier(comm);
    	adios_finalize (comm_rank);
    	MPI_Finalize ();

	return 0;
}
