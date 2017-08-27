#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"
#include "pt_structs.h"

ADIOS_FILE *adios_open_state(char *adios_file, enum ADIOS_READ_METHOD read_method, char *opts, MPI_Comm comm)
{
    	adios_read_init_method(read_method, comm, opts);
	ADIOS_FILE* afile = adios_read_open_file(adios_file, read_method, comm);
	return afile;
}

int adios_read_state(ADIOS_FILE *afile, pt_atoms *atoms_array, int proc_no, int time_step)
{
	int num_atoms, num_types;

	adios_schedule_read (afile, NULL, "num_atoms", 0, 1, &num_atoms);
	adios_schedule_read (afile, NULL, "num_types", 0, 1, &num_types);
	adios_perform_reads (afile, 1);

	printf("I am here %d %d....\n",num_atoms,num_types);
	if (alloc_atoms_array(atoms_array,num_atoms,num_types)!=0) return 1;  

	// printf("Next set of reads.: %d %d\n",proc_no,time_step);
	
	uint64_t start[2];
	uint64_t count[2];	
	start[0] = proc_no; start[1] = 0;
	count[0] = 1; count[1] = 6;
       	ADIOS_SELECTION *sel_dims  = adios_selection_boundingbox (2,start,count);

	adios_schedule_read (afile, sel_dims, "cell_dims", time_step,1,atoms_array->cell_dims);
	adios_perform_reads (afile, 1);
	
	printf("I am after cell_dims...\n"); fflush(stdout);

	uint64_t start1[2], count1[2];
	start1[0] = proc_no; start1[0] = 0;
	count1[0] = 1; count1[1] = atoms_array->num_types;
	ADIOS_SELECTION *sel_types = adios_selection_boundingbox (2,start1,count1);
		
	adios_schedule_read (afile, sel_types, "type_id", time_step,1,atoms_array->type_id);
	adios_schedule_read (afile, sel_types, "masses", time_step,1,atoms_array->masses);
	adios_perform_reads (afile, 1);	

	printf("I am here...\n");

	uint64_t start2[2], count2[2];
	start2[0] = proc_no; start2[1] = 0;
	count2[0] = 1; count2[1] = atoms_array->num_atoms;
	ADIOS_SELECTION *sel_atoms = adios_selection_boundingbox (2,start2,count2);

	adios_schedule_read (afile, sel_atoms, "atom_id", time_step,1,atoms_array->atom_id);
	adios_schedule_read (afile, sel_atoms, "atom_type", time_step,1,atoms_array->atom_type);
	adios_schedule_read (afile, sel_atoms, "px", time_step,1,atoms_array->px);
	adios_schedule_read (afile, sel_atoms, "py", time_step,1,atoms_array->py);
	adios_schedule_read (afile, sel_atoms, "pz", time_step,1,atoms_array->pz);
	adios_schedule_read (afile, sel_atoms, "imx", time_step,1,atoms_array->imx);
	adios_schedule_read (afile, sel_atoms, "imy", time_step,1,atoms_array->imy);
	adios_schedule_read (afile, sel_atoms, "imz", time_step,1,atoms_array->imz);
	adios_schedule_read (afile, sel_atoms, "atom_vid", time_step,1,atoms_array->atom_vid);
	adios_schedule_read (afile, sel_atoms, "vx", time_step,1,atoms_array->vx);
	adios_schedule_read (afile, sel_atoms, "vy", time_step,1,atoms_array->vy);
	adios_schedule_read (afile, sel_atoms, "vz", time_step,1,atoms_array->vz);
	adios_perform_reads (afile, 1);

	return 0;
}

int text_write_state(FILE *fp, pt_atoms *atoms_array)
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
	pt_atoms atoms_out;
	MPI_Comm comm = MPI_COMM_WORLD;
	int comm_rank, comm_size;
	enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_BP;
	int i;
	char filename[256];

    	MPI_Init (&argc, &argv);
    	MPI_Comm_rank (comm, &comm_rank);
	MPI_Comm_size (comm, &comm_size);

	if (argc!=5) {
		printf("Usage: %s <bp file> <output file> <select_id> <time step>\n",argv[0]);
		return 1;
	}

	char *bp_file  = argv[1];
	char *out_file = argv[2];
	int  select_id = atoi(argv[3]);
	int  time_step = atoi(argv[4]);

	printf("Reading the file....\n"); fflush(stdout);

	adios_init_noxml(comm);

	FILE *fpo = fopen(out_file,"w");
	if (fpo==NULL) return 1;

	int len = strlen(HEADER_STRING)+1;
	strcpy(atoms_out.header_str,HEADER_STRING);

	ADIOS_FILE *afile  = adios_open_state(bp_file, read_method, (char*)"", comm);
	if (afile==NULL) return 1;
	adios_read_state(afile,&atoms_out,select_id,time_step);
    	adios_read_close(afile);	
    	MPI_Barrier(comm);
    	adios_read_finalize_method(read_method);

	if (text_write_state(fpo,&atoms_out)!=0) {
		printf("Write error....\n");
		return 1;
	}
	fclose(fpo);

    	MPI_Finalize ();

	return 0;
}

