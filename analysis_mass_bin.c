#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SPECIFY_HALO_MASS_RANGE_CASE 

void collect_mass_bin_data()
{
  int i;
  int j;
  int k;

  i = (int) (log10(halo_part_num * 1e10 * g_head.Massarr[1] / 
              SPECIFY_HALO_MASS_RANGE_LEFT) / log10(2.0) );
  if(i >= MASS_BIN_NUM)
    i = MASS_BIN_NUM - 1;
  if(i < 0)
    i = 0;

  if(ellip[0] > 0 && ellip[1] > 0 && ellip[2] > 0)
    {
      mass_bin_data[i][0]++;
      mass_bin_data[i][1] += ellip[0] / ellip[2];                               /* a/c */
      mass_bin_data[i][2] += ellip[1] / ellip[2];                               /* b/c */
      mass_bin_data[i][3] += ellip[0] / ellip[1];                               /* a/b */
      j =  (int) (ellip[0] / ellip[2] / 0.05);                                  /* a/c distribution */
      k =  (int) (ellip[1] / ellip[2] / 0.05);                                  /* b/c distribution */
      if(j >= 20)
        j = 19;	   
      if(k >= 20)
        k = 19;	        
      mass_bin_data[i][j+4]++;
      mass_bin_data[i][k+24]++;
     
    }
}                                                                               /* end collect_mass_bin_data */


void mass_bin_data_out_put()
{
  int temp_halo_id;
  int temp_halo_part_num;
  double temp_halo_cen[3];
  double temp_over_den;
  int temp_link_num;
  int temp_halo_num_shell;
  double temp_shell_mass_cen[3];
  int temp_halo_num_inside;
  double temp_grid_body_mass_cen[3];
  double a;
  double b;
  double c;
  double vec_a[3];
  double vec_b[3];
  double vec_c[3]; 	          /* to read the shape output file */
#ifdef CAL_ANGULAR_MOMENTUM 
  double am[3];
#endif  
			
  int i, j, n_halo;
  int n_halo_temp[MASS_BIN_NUM];
  double d_factor;
  
  char file_name_mass[500];
  char file_name_a_c[500];
  char file_name_b_c[500];
  FILE *fp, *fp_a_c, *fp_b_c;
  

  for(i = 0; i < MASS_BIN_NUM; i++)
    {
    for(j = 1; j < 4; j++)
      if(mass_bin_data[i][0] > 0)
        mass_bin_data[i][j] = mass_bin_data[i][j] / mass_bin_data[i][0];  /* gives the average */
    
    n_halo_temp[i] = 0;          /* consistency check */
    }
            
  rewind(f_output);              /* file pointer rewind */
  
#ifndef CAL_ANGULAR_MOMENTUM  
  for(n_halo = 0; n_halo < prcessed_halo_num; n_halo++)
    {
    fscanf(f_output,
"%d	%d	%lf	%lf	%lf	%lf	%d	\
%d	%lf	%lf	%lf	%d	%lf	%lf	\
%lf	%lf	%lf	%lf	%lf	%lf	\
%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",
    &temp_halo_id, &temp_halo_part_num, &temp_halo_cen[0], &temp_halo_cen[1],
    &temp_halo_cen[2], &temp_over_den, &temp_link_num, &temp_halo_num_shell, 
    &temp_shell_mass_cen[0], &temp_shell_mass_cen[1], &temp_shell_mass_cen[2],
    &temp_halo_num_inside, &temp_grid_body_mass_cen[0], &temp_grid_body_mass_cen[1],
    &temp_grid_body_mass_cen[2], &a, &b, &c, &vec_a[0], &vec_a[1], &vec_a[2],
    &vec_b[0], &vec_b[1], &vec_b[2], &vec_c[0], &vec_c[1], &vec_c[2]);
#else
  for(n_halo = 0; n_halo < prcessed_halo_num; n_halo++)
    {
    fscanf(f_output,
"%d	%d	%lf	%lf	%lf	%lf	%d	\
%d	%lf	%lf	%lf	%d	%lf	%lf	\
%lf	%lf	%lf	%lf	%lf	%lf	\
%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n",
    &temp_halo_id, &temp_halo_part_num, &temp_halo_cen[0], &temp_halo_cen[1],
    &temp_halo_cen[2], &temp_over_den, &temp_link_num, &temp_halo_num_shell, 
    &temp_shell_mass_cen[0], &temp_shell_mass_cen[1], &temp_shell_mass_cen[2],
    &temp_halo_num_inside, &temp_grid_body_mass_cen[0], &temp_grid_body_mass_cen[1],
    &temp_grid_body_mass_cen[2], &a, &b, &c, &vec_a[0], &vec_a[1], &vec_a[2],
    &vec_b[0], &vec_b[1], &vec_b[2], &vec_c[0], &vec_c[1], &vec_c[2], &am[0], &am[1], &am[2]);
#endif    
    i = (int) (log10(temp_halo_part_num * 1e10 * g_head.Massarr[1] / 
              SPECIFY_HALO_MASS_RANGE_LEFT) / log10(2.0) );
    if(i >= MASS_BIN_NUM)
      i = MASS_BIN_NUM - 1;
    if(i < 0)
      i = 0;
      
    if(a > 0 && b > 0 && c > 0)
      {
	   n_halo_temp[i]++; 
	   mass_bin_data[i][44] += pow(a/c - mass_bin_data[i][1], 2);	  
	   mass_bin_data[i][45] += pow(b/c - mass_bin_data[i][2], 2);	  
	   mass_bin_data[i][46] += pow(a/b - mass_bin_data[i][3], 2);	    
      }               
    }
   
   for(i = 0; i < MASS_BIN_NUM; i++)
     {
	  if(n_halo_temp[i] != mass_bin_data[i][0])
	    {
	    printf("Error detected in calculate the error bar!\n");
	    exit(0); 
	    }
	  else
	    {
		 if(mass_bin_data[i][0] > 1)
		  {
		  d_factor = mass_bin_data[i][0] * (mass_bin_data[i][0] - 1);  
		  mass_bin_data[i][44] = sqrt(mass_bin_data[i][44] / d_factor);  
		  mass_bin_data[i][45] = sqrt(mass_bin_data[i][45] / d_factor); 
		  mass_bin_data[i][46] = sqrt(mass_bin_data[i][46] / d_factor); 			  
		  }                      
	    }
     }
      
  sprintf(
          file_name_mass,
          "%s%s%.2f%s",
          OUTPUT_FOLDER,
          "mass_a_c_z_",
          g_head.Redshift,
          ".txt"
         );
 fp = fopen(file_name_mass, "w+");
	
  for(i = 0; i < MASS_BIN_NUM; i++)
    if(mass_bin_data[i][0] > 1)
     fprintf(
            fp,
            "%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.6e	%.0f\n", 
            (pow(2.0,i) + pow(2.0,i+1)) / 2.0 * SPECIFY_HALO_MASS_RANGE_LEFT, 
            mass_bin_data[i][1],
            mass_bin_data[i][44],
            mass_bin_data[i][2],
            mass_bin_data[i][45],         
            mass_bin_data[i][3], 
            mass_bin_data[i][46],
            mass_bin_data[i][0]
            );

 fclose(fp);
 
  for(i = 0; i < MASS_BIN_NUM; i++)
    {
     if(mass_bin_data[i][0] > 0)  /* only output the good bins */
       {
       sprintf(
              file_name_a_c,
              "%s%s%.2f%s%s%d%s%s%d%s%.3e%s",
              OUTPUT_FOLDER,
              "a_c_distribution_z_",
              g_head.Redshift,
              "_(","2^",i,"-","2^",i+1,")_",
              SPECIFY_HALO_MASS_RANGE_LEFT,
              ".txt"
              );
             
      sprintf(
              file_name_b_c,
              "%s%s%.2f%s%s%d%s%s%d%s%.3e%s",
              OUTPUT_FOLDER,
              "b_c_distribution_z_",
              g_head.Redshift,
              "_(","2^",i,"-","2^",i+1,")_",
              SPECIFY_HALO_MASS_RANGE_LEFT,
              ".txt"
              );
     fp_a_c = fopen(file_name_a_c, "w+");
     fp_b_c = fopen(file_name_b_c, "w+");
   
     for(j=4; j<24; j++)
        {
        fprintf(fp_a_c,
               "%lf	%.6e	%.0f\n", 
               (j - 4  + 0.5) * 0.05, 
               mass_bin_data[i][j] / mass_bin_data[i][0] / 0.05,
               mass_bin_data[i][j]);  
        fprintf(fp_b_c,
               "%lf	%.6e	%.0f\n", 
               (j - 4  + 0.5) * 0.05, 
               mass_bin_data[i][j+20] / mass_bin_data[i][0] / 0.05,
               mass_bin_data[i][j+20]);              
          
        }
       
     fclose(fp_a_c);
     fclose(fp_b_c);
      }
   }

}                                                                      /* end mass_bin_data_out_put */
                                                    
#endif
