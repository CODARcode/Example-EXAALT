#ifndef PT_STRUCTS_H
#define PT_STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

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

#endif /* PT_STRUCTS_H */
