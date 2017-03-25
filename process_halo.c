#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


void process_exist_halo()
{
  int i,j;
#ifdef OUTPUT_WHOLE_HALO
  write_whole_halo_particle();
#endif  
  make_link_list();
  printf("Nearby neighbors done.\n");
  cal_halo_part_density();
  printf("Density construction done.\n");

  for(i=0 ; i<mass_shell_num ; i++)
    for(j=0 ; j<=link_num_step_num ; j++)
      {       
      update_link_num_mass_shell(j,i);                  /* update link_num and over density */
      det_mass_shell(&halo_cen);                        /* with mass shell number i, start from 1, and Link number = min_link_num+j */

      if(mass_shell_info() == 1) 
        {
	    get_inertia_tensor();
        printf("Inertia tensor for mass shell %d with link_num %d done.\n",i+1,link_num);
        eigen_inertia_tensor();
        get_halo_shape();
        
#ifdef CAL_ANGULAR_MOMENTUM 
        cal_angular_momentum();
#endif        

#if defined(SPECIFY_HALO_MASS_RANGE_CASE) && defined(FIX_LINK_NUM)
        collect_mass_bin_data();
#endif

        write_shape_file();
        
        prcessed_halo_num_each_file++;
        prcessed_halo_num++;        
        
#if defined(OUTPUT_PARTICLES_INFO) && defined(FIX_LINK_NUM)     /* only output if link_num is fixed */
        write_particle_file();
#endif 
	    }
	  else
	    {
		printf("There is no particle has density higher than the threshold in this halo.\n"); 	
	    }  
	    
      refresh_flags();
      }
  printf("Halo shape done.\n\n");
  free(p_halo);

}                                                       /* end process_exit_halo */

void refresh_flags()
{
  int i,j;
  for(i=0 ; i<halo_part_num ; i++)
    for(j=0 ; j<4 ; j++)
      p_halo[i].flag[j] = 0;
}                                                       /* end refresh_flags */

void update_link_num_mass_shell(int a, int b)
{
  link_num = min_link_num + a * link_num_step;
  over_den = max_over_den / pow(over_dec_fac , b);
}                                                       /*end update_link_num_mass_shell */

#ifdef CAL_ANGULAR_MOMENTUM
void cal_angular_momentum()
{
  float vel_mean[3];
  double total_halo_mass;
  float relate_pos[3];
  float relate_vel[3];
  int i;
  
  vel_mean[0] = vel_mean[1] = vel_mean[2] = 0.0;
  total_halo_mass = 0;
  
  for(i=0 ; i<halo_part_num; i++)
    {
	 vel_mean[0] += p_halo[i].vel[0] * p_halo[i].mass;	
	 vel_mean[1] += p_halo[i].vel[1] * p_halo[i].mass;		 
	 vel_mean[2] += p_halo[i].vel[2] * p_halo[i].mass;	  	
	 total_halo_mass += p_halo[i].mass;
    }
	
  for(i=0 ; i<3; i++)
    vel_mean[i] = vel_mean[i] / total_halo_mass;
    
  halo_am[0] = halo_am[1] = halo_am[2] = 0;  	
  
  for(i=0 ; i<halo_part_num; i++)
   {
    relate_center(p_halo[i].pos, halo_cen.pos, relate_pos);
    relate_center(p_halo[i].vel, vel_mean, relate_vel);
    
    halo_am[0] += p_halo[i].mass * (relate_pos[1] * relate_vel[2] - relate_pos[2] * relate_vel[1]);   /* x direction */
    halo_am[1] += p_halo[i].mass * (relate_pos[2] * relate_vel[0] - relate_pos[0] * relate_vel[2]);   /* y direction */
    halo_am[2] += p_halo[i].mass * (relate_pos[0] * relate_vel[1] - relate_pos[1] * relate_vel[0]);   /* z direction */	   
   }
   
  for(i=0; i<3; i++)
    halo_am[i] = halo_am[i] * pow(1 + halo_red, -1.5);   /* transfer to physical definition of AM */
}                                                        /* end cal_angular_momentum */	
#endif
