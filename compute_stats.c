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

int compute_mean(FILE *fp, pt_atoms *atoms_array, atomic_stats *stats_atoms) 
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

		stats_atoms->mean_px += atoms_array->px[i]; 
		stats_atoms->mean_py += atoms_array->py[i]; 
		stats_atoms->mean_pz += atoms_array->pz[i]; 
		stats_atoms->cnt += 1;
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

		stats_atoms->mean_vx += atoms_array->vx[i]; 
		stats_atoms->mean_vy += atoms_array->vy[i]; 
		stats_atoms->mean_vz += atoms_array->vz[i]; 
	}
	
	return 0;
}

int compute_std(FILE *fp, pt_atoms *atoms_array, atomic_stats *stats_atoms) 
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

		stats_atoms->std_px  += (atoms_array->px[i]-stats_atoms->mean_px)*(atoms_array->px[i]-stats_atoms->mean_px);
		stats_atoms->std_py  += (atoms_array->py[i]-stats_atoms->mean_py)*(atoms_array->py[i]-stats_atoms->mean_py);
		stats_atoms->std_pz  += (atoms_array->pz[i]-stats_atoms->mean_pz)*(atoms_array->pz[i]-stats_atoms->mean_pz);
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

		stats_atoms->std_vx  += (atoms_array->vx[i]-stats_atoms->mean_vx)*(atoms_array->vx[i]-stats_atoms->mean_vx);
		stats_atoms->std_vy  += (atoms_array->vy[i]-stats_atoms->mean_vy)*(atoms_array->vy[i]-stats_atoms->mean_vy);
		stats_atoms->std_vz  += (atoms_array->vz[i]-stats_atoms->mean_vz)*(atoms_array->vz[i]-stats_atoms->mean_vz);
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
	atomic_stats stats_atoms;
	

    	MPI_Init (&argc, &argv);
    	MPI_Comm_rank (comm, &comm_rank);
	MPI_Comm_size (comm, &comm_size);

	if (argc!=3) {
		printf("Usage: %s <input state file> <number of states>\n",argv[0]);
		return 1;
	}

	char *state_file = argv[1];
	int  num_states  = (int)atoi(argv[2]);

	stats_atoms.cnt = 0;
	stats_atoms.mean_px = 0.0;
	stats_atoms.mean_py = 0.0;
	stats_atoms.mean_pz = 0.0;
	stats_atoms.std_px = 0.0;
	stats_atoms.std_py = 0.0;
	stats_atoms.std_pz = 0.0;
	stats_atoms.mean_vx = 0.0;
	stats_atoms.mean_vy = 0.0;
	stats_atoms.mean_vz = 0.0;
	stats_atoms.std_vx = 0.0;
	stats_atoms.std_vy = 0.0;
	stats_atoms.std_vz = 0.0;
	
	FILE *fpi = fopen(state_file,"r");
	if (fpi==NULL) return 1;

	for (i=0;i<num_states;i++) {
		fscanf(fpi,"%s",filename);
		FILE *fpp = fopen(filename,"r");
		if (fpp==NULL) return 1;
		compute_mean(fpp,&atoms_array,&stats_atoms);
		fclose(fpp);	
		free_atoms_array(&atoms_array);
	}
	fclose(fpi);

	printf("ATOMS: %d %d\n",stats_atoms.cnt,num_states);

	stats_atoms.mean_px = stats_atoms.mean_px/stats_atoms.cnt;
	stats_atoms.mean_py = stats_atoms.mean_py/stats_atoms.cnt;
	stats_atoms.mean_pz = stats_atoms.mean_pz/stats_atoms.cnt;
	stats_atoms.mean_vx = stats_atoms.mean_vx/stats_atoms.cnt;
	stats_atoms.mean_vy = stats_atoms.mean_vy/stats_atoms.cnt;
	stats_atoms.mean_vz = stats_atoms.mean_vz/stats_atoms.cnt;

	fpi = fopen(state_file,"r");
	if (fpi==NULL) return 1;

	for (i=0;i<num_states;i++) {
		fscanf(fpi,"%s",filename);
		FILE *fpp = fopen(filename,"r");
		if (fpp==NULL) return 1;
		compute_std(fpp,&atoms_array,&stats_atoms);
		fclose(fpp);	
		free_atoms_array(&atoms_array);
	}
	fclose(fpi);

	stats_atoms.std_px = sqrt(stats_atoms.std_px/stats_atoms.cnt);	
	stats_atoms.std_py = sqrt(stats_atoms.std_py/stats_atoms.cnt);	
	stats_atoms.std_pz = sqrt(stats_atoms.std_pz/stats_atoms.cnt);	
	stats_atoms.std_vx = sqrt(stats_atoms.std_vx/stats_atoms.cnt);	
	stats_atoms.std_vy = sqrt(stats_atoms.std_vy/stats_atoms.cnt);	
	stats_atoms.std_vz = sqrt(stats_atoms.std_vz/stats_atoms.cnt);	

	printf("stats_atoms->cnt = %d;\n",stats_atoms.cnt);
	printf("stats_atoms->mean_px = %lf; stats_atoms->std_px = %lf;\n",stats_atoms.mean_px,stats_atoms.std_px);
	printf("stats_atoms->mean_py = %lf; stats_atoms->std_py = %lf;\n",stats_atoms.mean_py,stats_atoms.std_py);
	printf("stats_atoms->mean_pz = %lf; stats_atoms->std_pz = %lf;\n",stats_atoms.mean_pz,stats_atoms.std_pz);
	printf("stats_atoms->mean_vx = %lf; stats_atoms->std_vx = %lf;\n",stats_atoms.mean_vx,stats_atoms.std_vx);
	printf("stats_atoms->mean_vy = %lf; stats_atoms->std_vy = %lf;\n",stats_atoms.mean_vy,stats_atoms.std_vy);
	printf("stats_atoms->mean_vz = %lf; stats_atoms->std_vz = %lf;\n",stats_atoms.mean_vz,stats_atoms.std_vz);


	MPI_Barrier(comm);
    	MPI_Finalize ();

	return 0;
}
