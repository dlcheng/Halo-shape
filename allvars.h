#ifndef ALLVAR_H
#define ALLVAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"

typedef struct gadget_head GADGET_HEAD;
typedef struct gadget_head *PPT_GADGET_HEAD;
typedef struct part_snapshot PART_SNAPSHOT;
typedef struct part_snapshot *PPT_PART_SNAPSHOT;
typedef struct neig NEIG;
typedef struct neig *PPT_NEIG;
typedef struct part_halo PART_HALO;
typedef struct part_halo *PPT_PART_HALO;
typedef struct link_list LINK_LIST;
typedef struct link_list *PPT_LINK_LIST;

struct gadget_head
{  
  unsigned int Npart[6];
  double Massarr[6];
  double Time;
  double Redshift;
  int FlagSfr;
  int FlagFeedback;
  int Nall[6];
  int  FlagCooling;
  int NumFiles;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int FlagAge;
  int FlagMetals;
  int NallHW[6];
  int Flag_entr_ics;
  char unused[60];          
};

struct part_snapshot
{
  unsigned int id;
  float pos[3];
#ifdef CAL_ANGULAR_MOMENTUM
  float vel[3];
#endif  
  double mass;
};

struct neig
{
  PPT_PART_HALO p_ngb;                /* pointing to its neighbor */
  double d;                            /* the distance to the neighbor */
};

struct part_halo
{
  unsigned int id;
  float pos[3];
#ifdef CAL_ANGULAR_MOMENTUM  
  float vel[3];
#endif
  double mass;
  double density;

  NEIG ngb[NEIGHBORS_NUM];

  int flag[4];              /* 
                             flag[0] marks wheather it has been visited or is under visiting or will be visited in the future
                             flag[1] marks the it is belonged to the boundary, 1 for YES, 0 for NO
                             flag[2] marks the shell number possiblely assigned, start from 1
                             flag[3] marks the whether the particle is inside the isodensity shell, 1 for YES, 0 for NO
                            */
};
 
struct link_list
{
  double d;
  int place;            
};


/* SNAPSHOT IN THE MEMORY */
extern PPT_PART_SNAPSHOT p_snapshot;                  /* pointing to the whole snaphot */     
extern double rho_crit;                                /* the critical density */
extern GADGET_HEAD g_head;

/*HALO IN THE MEMORY */
extern PPT_PART_HALO p_halo;                          /* global pointer to the halo */

/*GLOBAL NUMBERS */
extern unsigned int all_part_num;                    /* all partilce number in the snapshot file */
extern int total_halo_num_each_file;                  /* the total halo number as specified in the *.particle file */
extern int prcessed_halo_num_each_file;
extern int prcessed_halo_num;
extern int omitted_halo_num;
extern int min_halo_part_num;
extern int file_id;

/*SINGLE HALO INFORMATION */
extern int halo_id;                                   /* the ID of a halo, start from 1 */
extern int halo_part_num;                             /* the particle number inside the halo */
extern PART_HALO halo_cen;                           /* the center of the halo */
extern double halo_red;                              /* redshift of the simulation */
extern int halo_num_inside;                           /* particle number in the mass shell */
extern int halo_num_shell;                            /* particle number on the mass shell */
extern float shell_mass_cen[3];                       /* mass center from particles on the shell */
extern float grid_body_mass_cen[3];                   /* mass center from particles inside the shell */

/*RELATED TO THE INERTIA TENSOR */
extern double inertia_tensor[9];
extern double eigen_val[3];                           /* eigen values */
extern double eigen_vec[3][3];                        /* eigen vector in row instead of column */

/*RELATED TO HALO SHAPE PARAMETERS */
extern double ellip[3];
#ifdef CAL_ANGULAR_MOMENTUM
extern double halo_am[3];                              /* the angular momentum of the halo */
#endif

/* THE FILES */
extern FILE *f_ahf_1;
extern FILE *f_ahf_2;
extern FILE *f_gadget;
#ifdef CAL_ANGULAR_MOMENTUM
extern FILE *f_gadget_vel;
#endif 
extern FILE *f_output;
#ifdef OUTPUT_NEIGHBORS_RELATION
extern FILE *f_my_neigh;
#endif
extern FILE *f_sorted_gadget;


/* THE LINK_NUM AND MASS_SHELL */
extern int min_link_num;
extern int link_num_step;
extern int link_num_step_num;
extern int link_num;                                  /* link_num value now */

extern double max_over_den;
extern double over_dec_fac;
extern int mass_shell_num;
extern double over_den;                               /* over density value now */

/*OTHER FEATURES */

#ifdef SPECIFY_HALO_MASS_RANGE_CASE
extern int halo_part_left;
extern int halo_part_right;
extern double mass_bin_data[MASS_BIN_NUM][47];       /* mass_bin_data[][0] -> number of halos inside the bin
                                                        * mass_bin_data[][1] -> a/c
                                                        * mass_bin_data[][2] -> b/c
                                                        * mass_bin_data[][3] -> a/b
                                                        * mass_bin_data[][4-23] -> distribution of a/c
                                                        * mass_bin_data[][24-43] -> distribution of b/c
                                                        * mass_bin_data[][44] -> error of a/c
                                                        * mass_bin_data[][45] -> error of b/c
                                                        * mass_bin_data[][46] -> error of a/b
                                                        */
                               
#endif

#endif
