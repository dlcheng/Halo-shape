#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int construct_ahf_halo()
{ 
  int halo_part_from_file_1;
  int i;
#ifndef AHF_V_1_0  
  float f_temp;
#endif
  char a_temp;
#ifdef AHF_V_1_0
  int temp_ID;
  int temp_hostHalo;
  int temp_numSubStruct;
  float temp_Mvir;
  int temp_type;
#endif 

  halo_id++;                                                                    /* halo ID increase */

/*
#ifdef SPECIFY_HALO_ID_RANGE_CASE
  if(halo_id > SPECIFY_HALO_ID_RANGE_RIGHT)
    {
     printf("\n%d halos have been processed.\n",prcessed_halo_num_each_file);
     printf("\nReach the expected halo number for analyzing.\n");
     printf("Job finished.\n");
     return 0;
    }
#endif
*/

  if(halo_id > total_halo_num_each_file)
    {
    printf("\n%d halos have been processed.\n",prcessed_halo_num_each_file);
    printf("%d halos have mass larger than the specified range.\n", omitted_halo_num);
    printf("Finish all %d halos in file %d\n",total_halo_num_each_file, file_id);
    return 0;
    }
#ifndef AHF_V_1_0
  fscanf(f_ahf_1,"%d %f %f %f %f",&halo_part_from_file_1,
                   &f_temp,&halo_cen.pos[0],&halo_cen.pos[1],&halo_cen.pos[2]);
#else
  fscanf(f_ahf_1,"%d %d %d %f %d %f %f %f",
         &temp_ID, &temp_hostHalo, &temp_numSubStruct, &temp_Mvir,
         &halo_part_from_file_1, &halo_cen.pos[0],&halo_cen.pos[1],&halo_cen.pos[2]);
#endif
  do{
    a_temp = fgetc(f_ahf_1);
    }while(a_temp != 10);

  re_scale_center();                                                             /* rescale halo center since the output 
	                                                                                of AHF is in Mpc/h,
                                                                                    while Gadget is Kpc/h */
  fscanf(f_ahf_2, "%d\n", &halo_part_num);

#ifdef SPECIFY_HALO_MASS_RANGE_CASE
  if(halo_part_num < halo_part_left)
    {
     printf("\n%d halos within the mass range.\n", prcessed_halo_num_each_file);
     printf("%d halos have mass larger than the higher end.\n",omitted_halo_num);
     printf("Reach the expected halo mass low end.\n");
     printf("Job finished in file %d.\n\n", file_id);
     return 0;
    }
#endif

  if(halo_part_num != halo_part_from_file_1)
    {
    printf("In the construction of AHF halo, particle number from two AHF files does not equal, something is wrong here!\n");
    exit(0);
    }

   if(halo_part_num <  min_halo_part_num )
     {
      printf("\nHalo ID: %d at (%.3e,%.3e,%.3e)[kpc/h] has particle number less than %d, we stop doing further analysics!\n",
             halo_id,halo_cen.pos[0],halo_cen.pos[1],halo_cen.pos[2],min_halo_part_num);
      return 0; 
     }

/*
#ifdef SPECIFY_HALO_ID_RANGE_CASE
   if(halo_id < SPECIFY_HALO_ID_RANGE_LEFT)
     { 
     int id;
     for(i=0 ; i<halo_part_num ; i++)
       fscanf(f_ahf_2,"%d\n",&id);
     printf("Skip halo ID: %d at (%.3e,%.3e,%.3e)[kpc/h]\n",halo_id,halo_cen.pos[0],halo_cen.pos[1],halo_cen.pos[2]);
     }
   else
     init_halo();

#endif
*/

#ifdef SPECIFY_HALO_MASS_RANGE_CASE
   if(halo_part_num > halo_part_right)
     { 
     int id;
     for(i=0 ; i<halo_part_num ; i++)
       {
#ifndef AHF_V_1_0	        
       fscanf(f_ahf_2,"%d\n",&id);
#else
       fscanf(f_ahf_2,"%d %d\n",&id, &temp_type);
#endif       
       }
     printf("Skip halo ID: %d at (%.3e,%.3e,%.3e)[kpc/h] with %d particles in file %d \n",
            halo_id,halo_cen.pos[0],halo_cen.pos[1],halo_cen.pos[2], halo_part_num, file_id);
     }
   else
     init_halo();
#endif

#ifdef NORMAL_MODE
   init_halo();
#endif

   return 1;
}                                                                 /* end construct_ahf_halo */

void re_scale_center()
{
 int i;
 for(i=0 ; i<3 ; i++)
  {
#ifndef AHF_V_1_0
   halo_cen.pos[i] = halo_cen.pos[i] * MPCTOKPC;
#else
   halo_cen.pos[i] = halo_cen.pos[i];
#endif   
  }
}                                                                 /* void re_scale_center */

void init_halo()
{  
   int i;
#ifdef AHF_V_1_0
   int temp_type;
#endif   
   if((p_halo = (PPT_PART_HALO) calloc(halo_part_num,sizeof(PART_HALO))) == NULL)
     {
     printf("Fail to allocate memory for halo, ID:%d at (%.3e,%.3e,%.3e)[kpc/h]\n",
             halo_id,halo_cen.pos[0],halo_cen.pos[1],halo_cen.pos[2]);
     exit(0);
     }

   for(i=0 ; i<halo_part_num ; i++)
     {
#ifndef AHF_V_1_0		 
     fscanf(f_ahf_2,"%d\n",&p_halo[i].id);
#else     
     fscanf(f_ahf_2,"%d	%d\n",&p_halo[i].id, &temp_type);     
#endif     
     init_single_halo_part(&p_halo[i]);                           /* fill the cell with particle information from Gadget file and set
                                                                     Flags to zero */
     }
   
   printf("Processing halo ID: %d at (%.3e,%.3e,%.3e)[kpc/h] with %d particles in file %d\n",
           halo_id, halo_cen.pos[0], halo_cen.pos[1], halo_cen.pos[2], halo_part_num, file_id);

}                                                                 /* end init_halo */

void init_single_halo_part(PPT_PART_HALO a)
{
  int i;
  int j;

#if (ID_START == 0)
  j = a->id;
#endif

#if (ID_START == 1)
  j = a->id - 1;
#endif

  if(a->id != p_snapshot[j].id)
    {
    printf("Error detected in finding particle from Gadget file, ID can not match!\n");
    exit(0);
    }

  for(i=0 ; i<3 ; i++)
    { 
     a->pos[i] = p_snapshot[j].pos[i];
#ifdef CAL_ANGULAR_MOMENTUM
     a->vel[i] = p_snapshot[j].vel[i];
#endif    
    }
  a->mass = p_snapshot[j].mass;
  a->density = 0.0;

  for(i=0 ; i<4 ; i++)
    a->flag[i] = 0;                                                 /* set all flags to 0 at first */
}                                                                   /* end init_single_part */
