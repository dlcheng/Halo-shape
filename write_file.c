#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void write_shape_file()
{
#ifndef CAL_ANGULAR_MOMENTUM	
  fprintf(f_output,
"%d	%d	%.6e	%.6e	%.6e	%.6e	%d	\
%d	%.6e	%.6e	%.6e	%d	%.6e	%.6e	\
%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	\
%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e\n",
halo_id, halo_part_num, halo_cen.pos[0], halo_cen.pos[1], halo_cen.pos[2], 
over_den, link_num,halo_num_shell, shell_mass_cen[0], shell_mass_cen[1], 
shell_mass_cen[2], halo_num_inside, grid_body_mass_cen[0], grid_body_mass_cen[1],
grid_body_mass_cen[2], ellip[0], ellip[1], ellip[2],
eigen_vec[0][0], eigen_vec[0][1], eigen_vec[0][2], eigen_vec[1][0],
eigen_vec[1][1],eigen_vec[1][2],eigen_vec[2][0],eigen_vec[2][1],eigen_vec[2][2]);     /* write halo shape results */
#else
  fprintf(f_output,
"%d	%d	%.6e	%.6e	%.6e	%.6e	%d	\
%d	%.6e	%.6e	%.6e	%d	%.6e	%.6e	\
%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	\
%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e\n",
halo_id, halo_part_num, halo_cen.pos[0], halo_cen.pos[1], halo_cen.pos[2], 
over_den, link_num,halo_num_shell, shell_mass_cen[0], shell_mass_cen[1], 
shell_mass_cen[2], halo_num_inside, grid_body_mass_cen[0], grid_body_mass_cen[1],
grid_body_mass_cen[2], ellip[0], ellip[1], ellip[2],
eigen_vec[0][0], eigen_vec[0][1], eigen_vec[0][2], eigen_vec[1][0],
eigen_vec[1][1],eigen_vec[1][2],eigen_vec[2][0],eigen_vec[2][1],eigen_vec[2][2], halo_am[0], halo_am[1], halo_am[2]);     /* write halo shape results */
#endif
}                                                            /* write file */

void write_particle_file()
{
  char part_inside[500];
  char part_shell[500];
  int i;
  FILE *fp_inside,*fp_shell;
  
  sprintf(part_inside,"%s%s%d%s%.0f%s%d%s",OUTPUT_FOLDER,"halo_id=",halo_id,"_O=",over_den,"_N=",link_num,"_inside.txt");
  sprintf(part_shell,"%s%s%d%s%.0f%s%d%s",OUTPUT_FOLDER,"halo_id=",halo_id,"_O=",over_den,"_N=",link_num,"_shell.txt");

 if((fp_inside = fopen(part_inside,"w+")) == NULL)
   {
   printf("Can't create file %s, we have to stop!\n",part_inside);
   exit(0);
   } 

 if((fp_shell = fopen(part_shell,"w+")) == NULL)
   {
   printf("Can't create file %s, we have to stop!\n",part_shell);
   exit(0);
   } 

 for(i=0 ; i<halo_part_num ; i++)
   {
    if(p_halo[i].flag[3] == 1)
      {
      fprintf(fp_inside,"%.3e	%.3e	%.3e\n",p_halo[i].pos[0],p_halo[i].pos[1],p_halo[i].pos[2]);
      if(p_halo[i].flag[2] != 0)
        fprintf(fp_shell,"%.3e	%.3e	%.3e\n",p_halo[i].pos[0],p_halo[i].pos[1],p_halo[i].pos[2]);
      }
   }

 fclose(fp_inside);
 fclose(fp_shell);
}                                                          /* end write_particle_file */

void write_whole_halo_particle()
{
  char file_string[500];
  int i;
  FILE *fp;
  sprintf(file_string,"%s%s%d%s",OUTPUT_FOLDER,"halo_id=",halo_id,"_whole.txt");
 
  if((fp = fopen(file_string,"w+")) == NULL)
   {
   printf("Can't create file %s, we have to stop!\n",file_string);
   exit(0);
   } 

  for(i=0 ; i<halo_part_num ; i++)
    fprintf(fp,"%.3e	%.3e	%.3e\n",p_halo[i].pos[0],p_halo[i].pos[1],p_halo[i].pos[2]);

  fclose(fp);
}                                                          /* end write_all_particle */
