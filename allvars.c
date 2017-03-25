#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allvars.h"

/* SNAPSHOT IN THE MEMORY */
 PPT_PART_SNAPSHOT p_snapshot;                                     /* pointing to the whole snaphot */     
 double rho_crit;                                                 /* the critical density */
 GADGET_HEAD g_head;

/*HALO IN THE MEMORY */
 PPT_PART_HALO p_halo;                                             /* global pointer to the halo */

/*GLOBAL NUMBERS */
 unsigned int all_part_num;                                       /* all partilce number in the snapshot file */
 int total_halo_num_each_file;                                    /* the total halo number as specified in the *.particle file */
 int prcessed_halo_num_each_file;
 int prcessed_halo_num;
 int omitted_halo_num;
 int min_halo_part_num;
 int file_id;
 
/*SINGLE HALO INFORMATION */
 int halo_id;                                                      /* the ID of a halo, start from 1 */
 int halo_part_num;                                                /* the particle number inside the halo */
 PART_HALO halo_cen;                                                /* the center of the halo */
 double halo_red;                                                  /* redshift of the simulation */
 int halo_num_inside;                                               /* particle number in the mass shell */
 int halo_num_shell;                                                /* particle number on the mass shell */
 float shell_mass_cen[3];                                           /* mass center from particles on the shell */
 float grid_body_mass_cen[3];                                       /* mass center from particles inside the shell */

/*RELATED TO THE INERTIA TENSOR */
 double inertia_tensor[9];
 double eigen_val[3];                                               /* eigen values */
 double eigen_vec[3][3];                                            /* eigen vector in row instead of column */

/*RELATED TO HALO SHAPE PARAMETERS */
 double ellip[3];
#ifdef CAL_ANGULAR_MOMENTUM
 double halo_am[3];                              /* the angular momentum of the halo */
#endif

/* THE FILES */
 FILE *f_ahf_1;
 FILE *f_ahf_2;
 FILE *f_gadget;
#ifdef CAL_ANGULAR_MOMENTUM
 FILE *f_gadget_vel;
#endif 
 FILE *f_output;
#ifdef OUTPUT_NEIGHBORS_RELATION
 FILE *f_my_neigh;
#endif
 FILE *f_sorted_gadget;


/* THE LINK_NUM AND MASS_SHELL */
 int min_link_num;
 int link_num_step;
 int link_num_step_num;
 int link_num;                                                        /* link_num value now */

 double max_over_den;
 double over_dec_fac;
 int mass_shell_num;
 double over_den;                                                    /* over density value now */
/*OTHER FEATURES */
#ifdef SPECIFY_HALO_MASS_RANGE_CASE
 int halo_part_left;
 int halo_part_right;
 double mass_bin_data[MASS_BIN_NUM][47];
#endif
