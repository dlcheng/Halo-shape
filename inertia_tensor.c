#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define WORK_SIZE 90

int mass_shell_info()
{
 int i;
 int num_inside = 0;
 int num_shell = 0;
 float shell_x = 0.0;
 float shell_y = 0.0;
 float shell_z = 0.0;
 float grid_x = 0.0;
 float grid_y = 0.0;
 float grid_z = 0.0;

 for(i = 0 ; i < halo_part_num; i++)
  {
    if(p_halo[i].flag[3] == 1)
     {
       num_inside++;
       grid_x += p_halo[i].pos[0];
       grid_y += p_halo[i].pos[1];
       grid_z += p_halo[i].pos[2];

     if(p_halo[i].flag[2] != 0)
       {
         num_shell++;
         shell_x += p_halo[i].pos[0];
         shell_y += p_halo[i].pos[1];
         shell_z += p_halo[i].pos[2];
       }
     }   
  }

  if(num_inside == 0)
    return 0;
    
#ifdef SHELL_INERTIA_TENSOR
  if(num_shell == 0)
    return 0;
#endif  
  
  if(num_inside > 0)
    {
    grid_x = grid_x /(float) num_inside;
    grid_y = grid_y /(float) num_inside;
    grid_z = grid_z /(float) num_inside;
    }

  if(num_shell > 0)
    {
    shell_x = shell_x /(float) num_shell;
    shell_y = shell_y /(float) num_shell;
    shell_z = shell_z /(float) num_shell;
    } 
  
  halo_num_inside = num_inside;
  halo_num_shell = num_shell;
  shell_mass_cen[0] = shell_x;
  shell_mass_cen[1] = shell_y;
  shell_mass_cen[2] = shell_z;  
  grid_body_mass_cen[0] = grid_x;
  grid_body_mass_cen[1] = grid_y;
  grid_body_mass_cen[2] = grid_z;
  
  return 1;

}                                                          /* void mass_shell_info */

void get_inertia_tensor()
{
 float origin[3];                                         /* the origin of the coordinate */
 float relate[3];                                         /* the position related to the center */
 int i;

#ifdef HALO_CENTER_AS_ORIGIN
 copy_center(origin,halo_cen.pos);
#endif

#ifdef SHELL_MASS_CENTER_AS_ORIGIN
 copy_center(origin,shell_mass_cen);
#endif

#ifdef GRID_BODY_MASS_CENTER_AS_ORIGIN
 copy_center(origin,grid_body_mass_cen);
#endif

 for(i=0 ; i<9 ; i++)
  inertia_tensor[i] = 0.0;                                 /* all set to zero before calculation */

 for(i=0 ; i<halo_part_num ; i++)
  {
  relate_center(p_halo[i].pos,origin,relate);
#ifdef GRID_BODY_INERTIA_TENSOR
  if(p_halo[i].flag[3] == 1)                               /* particle labels as inside the contour */
   {
   inertia_tensor[0] += p_halo[i].mass * (relate[1] * relate[1] + relate[2] * relate[2]);
   inertia_tensor[1] += -1.0 * p_halo[i].mass * relate[0] * relate[1];
   inertia_tensor[2] += -1.0 * p_halo[i].mass * relate[0] * relate[2];
   inertia_tensor[3] = inertia_tensor[1];
   inertia_tensor[4] += p_halo[i].mass * (relate[0] * relate[0] + relate[2] * relate[2]);
   inertia_tensor[5] += -1.0 * p_halo[i].mass * relate[1] * relate[2];
   inertia_tensor[6] = inertia_tensor[2];
   inertia_tensor[7] = inertia_tensor[5];
   inertia_tensor[8] += p_halo[i].mass * (relate[0] * relate[0] + relate[1] * relate[1]);
   }
#endif
#ifdef SHELL_INERTIA_TENSOR
  if(p_halo[i].flag[2] != 0)
   {
   inertia_tensor[0] += p_halo[i].mass * (relate[1] * relate[1] + relate[2] * relate[2]);
   inertia_tensor[1] += -1.0 * p_halo[i].mass * relate[0] * relate[1];
   inertia_tensor[2] += -1.0 * p_halo[i].mass * relate[0] * relate[2];
   inertia_tensor[3] = inertia_tensor[1];
   inertia_tensor[4] += p_halo[i].mass * (relate[0] * relate[0] + relate[2] * relate[2]);
   inertia_tensor[5] += -1.0 * p_halo[i].mass * relate[1] * relate[2];
   inertia_tensor[6] = inertia_tensor[2];
   inertia_tensor[7] = inertia_tensor[5];
   inertia_tensor[8] += p_halo[i].mass * (relate[0] * relate[0] + relate[1] * relate[1]);
   }
#endif
  }
}                                                          /* end get_inertia_tensor */


void eigen_inertia_tensor()
{
 gsl_matrix_view m = gsl_matrix_view_array(inertia_tensor, 3, 3);
 gsl_vector *eval = gsl_vector_alloc(3);
 gsl_matrix *evec = gsl_matrix_alloc(3,3);

 gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(WORK_SIZE);

 gsl_eigen_symmv(&m.matrix,eval,evec,w);

 gsl_eigen_symmv_free(w);

 gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);  /* sort the eigen value as increasing way */

 int i,j;
 for(i=0 ; i<3 ; i++)
  {
  eigen_val[i] = gsl_vector_get(eval,i);
  for(j=0 ; j<3 ; j++)
  eigen_vec[i][j] = gsl_matrix_get(evec,j,i);            /* eigen vector in row */
  }
}                                                        /* end eigen_inertia_tensor */

void get_halo_shape()
{ 
 double m;                                              /* mass of total particles under consideration */
#ifdef GRID_BODY_INERTIA_TENSOR
 m = g_head.Massarr[1] * halo_num_inside;
 ellip[0] = sqrt(2.5 * (eigen_val[0] + eigen_val[1] - eigen_val[2]) / m);
 ellip[1] = sqrt(2.5 * (eigen_val[0] + eigen_val[2] - eigen_val[1]) / m);
 ellip[2] = sqrt(2.5 * (eigen_val[1] + eigen_val[2] - eigen_val[0]) / m);
#endif

#ifdef SHELL_INERTIA_TENSOR
 m = g_head.Massarr[1] * halo_num_shell;
 ellip[0] = sqrt(0.75 * (eigen_val[0] + eigen_val[1] - eigen_val[2]) / m);
 ellip[1] = sqrt(0.75 * (eigen_val[0] + eigen_val[2] - eigen_val[1]) / m);
 ellip[2] = sqrt(0.75 * (eigen_val[1] + eigen_val[2] - eigen_val[0]) / m);
#endif 
}                                                        /* end get_halo_shape */

void copy_center(float * a, float *b )                   /* copy a from b */
{
 int i;
 for(i=0 ; i<3 ; i++)
   a[i] = b[i];
}                                                        /* end copy_center */

void relate_center(float * a, float *b , float *c)       /* c = a - b */
{
 int i;
 for(i=0 ; i<3 ; i++)
  c[i] = a[i] - b[i];
}                                                       /* end relate_center */

#undef WORK_SIZE
