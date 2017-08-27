#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"

#define HEADER_STRING "LAMMPS data file via write_data, version 30 Nov 2016, timestep = 0\n"

#define MAX_LINE_LEN  1024
typedef struct _pt_atoms {
	char    *header_str;
	int     state_id[64];
	int     num_atoms;
	int     num_types;

	double cell_dims[3][2]; /* minx,maxx, miny,maxy, minz,maxz */
	int    *type_id;
	double *masses;	

	/* len: num_atoms */	
	int    *atom_id;
	int    *atom_type; 
	double *px, *py, *pz;
	int    *imx, *imy, *imz;
	int    *atom_vid;
	double *vx, *vy, *vz;
} pt_atoms;

void free_atoms_array(pt_atoms *atoms_array)
{
	free(atoms_array->header_str);
	free(atoms_array->type_id);
	free(atoms_array->masses);	

	free(atoms_array->atom_id);
	free(atoms_array->atom_type); 
	free(atoms_array->px);
	free(atoms_array->py);
	free(atoms_array->pz);
	free(atoms_array->imx);
	free(atoms_array->imy);
	free(atoms_array->imz);
	free(atoms_array->atom_vid);
	free(atoms_array->vx);
	free(atoms_array->vy);
	free(atoms_array->vz);
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

    	adios_define_var (*gh, "state_id", "", adios_integer, "1,64", "num_procs,64", "proc_no,0");
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

	groupsize  = 4*sizeof(int);   /* num_procs, proc_no, num_atoms, num_types */
	groupsize += 64*sizeof(char); /* state_id */
	groupsize += 3*2*sizeof(int); /* cell_dims */
	groupsize += atoms_array->num_types*sizeof(int);    /* type_id */
	groupsize += atoms_array->num_types*sizeof(double); /* masses */

	groupsize += atoms_array->num_atoms*sizeof(int)*2;    /* atom_id, atom_type */
	groupsize += atoms_array->num_atoms*sizeof(double)*3; /* px,py,pz */
	groupsize += atoms_array->num_atoms*sizeof(int)*3;    /* imx, imy, imz */
	groupsize += atoms_array->num_atoms*sizeof(int);      /* atom_vid */
	groupsize += atoms_array->num_atoms*sizeof(double)*3; /* vx, vy, vz */

    	adios_group_size (fh, groupsize, &groupTotalSize);

	adios_write (fh, "num_procs", &num_procs);
	adios_write (fh, "proc_no", &proc_no);
        adios_write (fh, "num_atoms", &atoms_array->num_atoms);
        adios_write (fh, "num_types", &atoms_array->num_types);

        adios_write (fh, "state_id",  &atoms_array->state_id);

	adios_write (fh, "cell_dims", atoms_array->cell_dims);
	
	adios_write (fh, "type_id", atoms_array->type_id);
	adios_write (fh, "masses", atoms_array->masses);

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

ADIOS_FILE *adios_open_state(char *adios_file, enum ADIOS_READ_METHOD read_method, char *opts, MPI_Comm comm)
{
    	adios_read_init_method(read_method, comm, opts);

	printf("I am here...\n");

	ADIOS_FILE* afile = adios_read_open_file(adios_file, read_method, comm);
	
#if 0
	ADIOS_FILE* afile = adios_read_open(adios_file, read_method, comm, ADIOS_LOCKMODE_ALL, 30.0);
	printf("I am here...\n");
	int num_iter = 0;
	while (adios_errno == err_file_not_found && num_iter<20) {
		printf("I am here before open stream....\n"); fflush(stdout);
		afile = adios_read_open(adios_file, read_method, comm, ADIOS_LOCKMODE_ALL, 30.0);
		usleep((unsigned int)100000); // sleep for 100 milliseconds
		num_iter++;
	}
	if (adios_errno == err_end_of_stream) {
		fprintf (stderr, "Stream terminated before open: %s\n", adios_errmsg());
		return NULL;
	} else if (afile == NULL) {
		fprintf (stderr, "Error at opening: %s\n", adios_errmsg());
		return NULL;
	}
#endif 
	printf("I am here before return...\n");

	return afile;
}

int adios_read_state(ADIOS_FILE *afile, pt_atoms *atoms_array)
{

#if 0
	while(adios_errno != err_end_of_stream && adios_errno!=err_step_notready){ 
#endif 
		printf("I am waiting to read...\n"); fflush(stdout);
		adios_schedule_read (afile, NULL, "state_id", 0, 1, &atoms_array->state_id);
		adios_schedule_read (afile, NULL, "num_atoms", 0, 1, &atoms_array->num_atoms);
		adios_schedule_read (afile, NULL, "num_types", 0, 1, &atoms_array->num_types);
		adios_perform_reads (afile, 1);
#if 0
		while (adios_errno==err_step_notready) {
			adios_schedule_read (afile, NULL, "state_id", 0, 1, &atoms_array->state_id);
			adios_schedule_read (afile, NULL, "num_atoms", 0, 1, &atoms_array->num_atoms);
			adios_schedule_read (afile, NULL, "num_types", 0, 1, &atoms_array->num_types);
			adios_perform_reads (afile, 1);
		}
#endif 

		printf("I am here %d %d %d....\n",atoms_array->state_id,atoms_array->num_atoms,atoms_array->num_types);

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

		printf("Next set of reads.\n");
	
		uint64_t start[2];
		uint64_t count[2];	
		start[0] = 0; start[1] = 0;
		count[0] = 3; count[1] = 2;
            	ADIOS_SELECTION *sel  = adios_selection_boundingbox (2,start,count);

		uint64_t start1[1], count1[1];
		start1[0] = 0; count1[0] = atoms_array->num_atoms;
		ADIOS_SELECTION *sel1 = adios_selection_boundingbox (1,start1,count1);

		start1[0] = 0; count1[0] = atoms_array->num_types;
		ADIOS_SELECTION *sel2 = adios_selection_boundingbox (1,start1,count1);
		
		adios_schedule_read (afile, sel, "cell_dims", 0,1,atoms_array->cell_dims);
		adios_schedule_read (afile, sel2, "type_id", 0,1,atoms_array->type_id);
		adios_schedule_read (afile, sel2, "masses", 0,1,atoms_array->masses);
		adios_schedule_read (afile, sel1, "atom_id", 0,1,atoms_array->atom_id);
		adios_schedule_read (afile, sel1, "atom_type", 0,1,atoms_array->atom_type);
		adios_schedule_read (afile, sel1, "px", 0,1,atoms_array->px);
		adios_schedule_read (afile, sel1, "py", 0,1,atoms_array->py);
		adios_schedule_read (afile, sel1, "pz", 0,1,atoms_array->pz);
		adios_schedule_read (afile, sel1, "imx", 0,1,atoms_array->imx);
		adios_schedule_read (afile, sel1, "imy", 0,1,atoms_array->imy);
		adios_schedule_read (afile, sel1, "imz", 0,1,atoms_array->imz);
		adios_schedule_read (afile, sel1, "atom_vid", 0,1,atoms_array->atom_vid);
		adios_schedule_read (afile, sel1, "vx", 0,1,atoms_array->vx);
		adios_schedule_read (afile, sel1, "vy", 0,1,atoms_array->vy);
		adios_schedule_read (afile, sel1, "vz", 0,1,atoms_array->vz);
		adios_perform_reads (afile, 1);
#if 0
		while (adios_errno==err_step_notready) {
			printf("I am waiting...\n");
			adios_schedule_read (afile, sel, "cell_dims", 0,1,atoms_array->cell_dims);
			adios_schedule_read (afile, sel2, "type_id", 0,1,atoms_array->type_id);
			adios_schedule_read (afile, sel2, "masses", 0,1,atoms_array->masses);
			adios_schedule_read (afile, sel1, "atom_id", 0,1,atoms_array->atom_id);
			adios_schedule_read (afile, sel1, "atom_type", 0,1,atoms_array->atom_type);
			adios_schedule_read (afile, sel1, "px", 0,1,atoms_array->px);
			adios_schedule_read (afile, sel1, "py", 0,1,atoms_array->py);
			adios_schedule_read (afile, sel1, "pz", 0,1,atoms_array->pz);
			adios_schedule_read (afile, sel1, "imx", 0,1,atoms_array->imx);
			adios_schedule_read (afile, sel1, "imy", 0,1,atoms_array->imy);
			adios_schedule_read (afile, sel1, "imz", 0,1,atoms_array->imz);
			adios_schedule_read (afile, sel1, "atom_vid", 0,1,atoms_array->atom_vid);
			adios_schedule_read (afile, sel1, "vx", 0,1,atoms_array->vx);
			adios_schedule_read (afile, sel1, "vy", 0,1,atoms_array->vy);
			adios_schedule_read (afile, sel1, "vz", 0,1,atoms_array->vz);
			adios_perform_reads (afile, 1);
		}
#endif 

#if 0	
		adios_release_step(afile);
    		adios_advance_step(afile, 0, 3.0);
		printf("ADIOS advance step....\n"); fflush(stdout);
    	}
#endif 
	return 0;
}


int read_state(FILE *fp, pt_atoms *atoms_array) 
{
	int len = MAX_LINE_LEN;
	char *str; 
	char tmp_str[128];
	int i;

	str = (char*)malloc(sizeof(char)*len);
	if (str==NULL) return 1;

	// Read and copy the header string
	str = fgets(str,len,fp);
	if (str==NULL) return 1;
	int str_len = strlen(str);	
	printf("READING: %s\n",str);
	atoms_array->header_str = (char*)malloc(sizeof(char)*(str_len+1));
	strcpy(atoms_array->header_str,str);

	// Read number of atoms
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	str = fgets(str,len,fp); 
	if (str==NULL) return 1;
	sscanf(str,"%d",&(atoms_array->num_atoms));
	// Read number of types
	str = fgets(str,len,fp); 
	if (str==NULL) return 1;
	sscanf(str,"%d",&(atoms_array->num_types));

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
	
	// Read simulation cell dimensions
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	str = fgets(str,len,fp); 
	if (str==NULL) return 1;
	sscanf(str,"%lg %lg",&(atoms_array->cell_dims[0][0]),&(atoms_array->cell_dims[0][1]));
	str = fgets(str,len,fp); 
	if (str==NULL) return 1;
	sscanf(str,"%lg %lg",&(atoms_array->cell_dims[1][0]),&(atoms_array->cell_dims[1][1]));
	str = fgets(str,len,fp); 
	if (str==NULL) return 1;
	sscanf(str,"%lg %lg",&(atoms_array->cell_dims[2][0]),&(atoms_array->cell_dims[2][1]));

	// Read masses
	atoms_array->masses  = (double*)malloc(sizeof(double)*atoms_array->num_types);
	atoms_array->type_id        = (int*)malloc(sizeof(int)*atoms_array->num_types);
	if (atoms_array->masses==NULL) return 1;
	if (atoms_array->type_id==NULL) return 1;
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	str = fgets(str,len,fp); // Skip Masses string 
	if (str==NULL) return 1;
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	for (i=0;i<atoms_array->num_types;i++) {
		str = fgets(str,len,fp); 
		if (str==NULL) return 1;
		sscanf(str,"%d %lg",&(atoms_array->type_id[i]),&(atoms_array->masses[i]));
	}

	// Read atoms
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	str = fgets(str,len,fp); // Skip atoms string  
	if (str==NULL) return 1;
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	
	for (i=0;i<atoms_array->num_atoms;i++) {
		str = fgets(str,len,fp); // Skip empty line
		if (str==NULL) return 1;
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

	// Read atom velocities
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	str = fgets(str,len,fp); // Skip Velocities string  
	if (str==NULL) return 1;
	str = fgets(str,len,fp); // Skip empty line
	if (str==NULL) return 1;
	for (i=0;i<atoms_array->num_atoms;i++) {
		str = fgets(str,len,fp); // Skip empty line
		if (str==NULL) return 1;
		sscanf(str,"%d %lg %lg %lg",
				&(atoms_array->atom_vid[i]),
				&(atoms_array->vx[i]),
				&(atoms_array->vy[i]),
				&(atoms_array->vz[i]));
	}
	
	return 0;
}

int write_state(FILE *fp, pt_atoms *atoms_array)
{
	int i;

	fprintf(fp,"%s",atoms_array->header_str);	
	fprintf(fp,"\n");
	fprintf(fp,"%d atoms\n",atoms_array->num_atoms);	
	fprintf(fp,"%d atom types\n",atoms_array->num_types);	
	fprintf(fp,"\n");

	fprintf(fp,"%.16e %.16e xlo xhi\n",atoms_array->cell_dims[0][0],atoms_array->cell_dims[0][1]);
	fprintf(fp,"%.16e %.16e ylo yhi\n",atoms_array->cell_dims[1][0],atoms_array->cell_dims[1][1]);
	fprintf(fp,"%.16e %.16e zlo zhi\n",atoms_array->cell_dims[2][0],atoms_array->cell_dims[2][1]);
	fprintf(fp,"\n");
	fprintf(fp,"Masses\n");
	fprintf(fp,"\n");
	for (i=0;i<atoms_array->num_types;i++) 
		fprintf(fp,"%d %lg\n",atoms_array->type_id[i],atoms_array->masses[i]);
	fprintf(fp,"\n");	

	fprintf(fp,"Atoms # atomic\n");
	fprintf(fp,"\n");
	for (i=0;i<atoms_array->num_atoms;i++) {
		fprintf(fp,"%d %d %.16le %.16le %.16le %d %d %d\n",
				atoms_array->atom_id[i],
				atoms_array->atom_type[i],
				atoms_array->px[i],
				atoms_array->py[i],
				atoms_array->pz[i],
				atoms_array->imx[i],
				atoms_array->imy[i],
				atoms_array->imz[i]);
	}
	fprintf(fp,"\n");

	fprintf(fp,"Velocities\n");
	fprintf(fp,"\n");
	for (i=0;i<atoms_array->num_atoms;i++) {
		fprintf(fp,"%d %.16le %.16le %.16le\n",
				atoms_array->atom_vid[i],
				atoms_array->vx[i],
				atoms_array->vy[i],
				atoms_array->vz[i]);
	}

	return 0;
}

int main(int argc, char *argv[])
{
	pt_atoms atoms_array;
	pt_atoms atoms_out;
	MPI_Comm comm = MPI_COMM_WORLD;
	int64_t gh;
	int comm_rank, comm_size;
	enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_BP;
	int i;
	char filename[256];

    	MPI_Init (&argc, &argv);
    	MPI_Comm_rank (comm, &comm_rank);
	MPI_Comm_size (comm, &comm_size);

	if (argc!=7) {
		printf("Usage: %s <input state file> <bp file> <transport method> <opts> <output file> <1:write/0:read>\n",argv[0]);
		return 1;
	}

	if (!strcmp(argv[6],"1")) {

		/* Create adios structure */
		adios_init_noxml(comm);
		adios_declare(&gh,argv[3],argv[4]);

		FILE *fpi = fopen(argv[1],"r");
		if (fpi==NULL) return 1;

		int first_time = 1;
		for (i=0;i<256;i++) {
			fscanf(fpi,"%s",filename);
			if (comm_rank==(i%comm_size)) {
				FILE *fpp = fopen(filename,"r");
				if (fpp==NULL) return 1;

				if (read_state(fpp,&atoms_array)!=0) {
					printf("Read Error.....\n");
					return 1;
				}
				if (first_time) {	
					adios_write_state(argv[2],(char*)"w",comm,comm_size,comm_rank,&atoms_array);
					first_time =0;
				} else {
					adios_write_state(argv[2],(char*)"a",comm,comm_size,comm_rank,&atoms_array);
				}
			
				fclose(fpp);	
				free_atoms_array(&atoms_array);
			} else {
			}
		}
		fclose(fpi);
		MPI_Barrier(comm);
    		adios_finalize (comm_rank);
		printf("Done writing the bp file...\n"); fflush(stdout);
	} else if (!strcmp(argv[6],"0")) {
		printf("Reading the file....\n"); fflush(stdout);

		adios_init_noxml(comm);

		FILE *fpo = fopen(argv[5],"w");
		if (fpo==NULL) return 1;

		int len = strlen(HEADER_STRING)+1;
		atoms_out.header_str = (char*)malloc(sizeof(char)*len);
		if (atoms_out.header_str==NULL) return 1;
		strcpy(atoms_out.header_str,HEADER_STRING);

		ADIOS_FILE *afile  = adios_open_state(argv[2], read_method, (char*)"", comm);
		if (afile==NULL) return 1;
		adios_read_state(afile,&atoms_out);
    		adios_read_close(afile);	
    		MPI_Barrier(comm);
    		adios_read_finalize_method(read_method);

		if (write_state(fpo,&atoms_out)!=0) {
			printf("Write error....\n");
			return 1;
		}
		fclose(fpo);
	}

    	MPI_Finalize ();

	return 0;
}


#if 0
	FILE *fpo = fopen(argv[2],"w");
	if (fpo==NULL) return 1;

	if (read_state(fpi,&atoms_array)!=0) {
		printf("Read Error.....\n");
		return 1;
	} 

	if (write_state(fpo,&atoms_array)!=0) {
		printf("Write error....\n");
		return 1;
	}

	free_atoms_array(&atoms_array);

	return 0;
}
#endif 
