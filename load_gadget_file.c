#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void load_sort_snapshot()
{
  int block_size;
  int i;
  float pos[3];
#ifdef CAL_ANGULAR_MOMENTUM
  float vel[3];
#endif  
  unsigned int id;
  
#ifdef STORE_SORTED_GADGET_FILE  
  char filename_sort[500];
  
  sprintf(filename_sort,"%s%s%s",INPUT_FOLDER,GADGET_FILE,"_gadget_sorted.dat");
  
  if((f_sorted_gadget = fopen(filename_sort,"rb")) != NULL)
    {
     printf("\nSorted gadget file %s%s detected.\n",GADGET_FILE,"_gadget_sorted.dat");
     load_sorted_gadget_file();
    }
  else
    {
    printf("\nSorted gadget file is not detected, will generate new one\n");
    file_open_first();                                               /* open gadget file first */

    fread(&block_size,sizeof(int),1,f_gadget);
    fread(&g_head,256,1,f_gadget);
    fread(&block_size,sizeof(int),1,f_gadget);
    if(block_size != 256)
     {
      printf("Error detected in reading the head block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
     }
#ifdef CAL_ANGULAR_MOMENTUM
    fread(&block_size,sizeof(int),1,f_gadget_vel);
    fread(&g_head,256,1,f_gadget);
    fread(&block_size,sizeof(int),1,f_gadget_vel);
    if(block_size != 256)
     {
      printf("Error detected in reading the head block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
     }
     for(i=0 ; i<g_head.Npart[1] ; i++)
       fread(pos, 3*sizeof(float), 1, f_gadget_vel);
       
    fread(&block_size,sizeof(int),1,f_gadget_vel);
    fread(&block_size,sizeof(int),1,f_gadget_vel);           
#endif     
   
    halo_red = g_head.Redshift;
    all_part_num = g_head.Npart[1];                            /* initialize these two values */
    rho_crit = crit_density();                                 /* calculate the critical density of the universe with Gadget unit */  

    if((p_snapshot = (PPT_PART_SNAPSHOT) calloc(all_part_num,sizeof(PART_SNAPSHOT))) == NULL)
      {
      printf("Failed to allocate memory for the Gadget file %s!\n",GADGET_FILE);
      exit(0);
      }

    fread(&block_size,sizeof(int),1,f_gadget);
    for(i=0 ; i<all_part_num ; i++)
      {
      fread(pos,3*sizeof(float),1,f_gadget);      
      p_snapshot[i].pos[0] = pos[0];
      p_snapshot[i].pos[1] = pos[1];
      p_snapshot[i].pos[2] = pos[2];
      p_snapshot[i].mass = g_head.Massarr[1];
#ifdef CAL_ANGULAR_MOMENTUM
      fread(vel, 3*sizeof(float), 1, f_gadget_vel);
      p_snapshot[i].vel[0] = vel[0];
      p_snapshot[i].vel[1] = vel[1];
      p_snapshot[i].vel[2] = vel[2];
#endif      
      }
    fread(&block_size,sizeof(int),1,f_gadget);
    if(block_size != 3*sizeof(float)*all_part_num)
      {
      printf("Error detected in copying the position block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
      }

    fread(&block_size,sizeof(int),1,f_gadget);
    for(i=0 ; i<all_part_num ; i++)
      fread(pos,3*sizeof(float),1,f_gadget);
    fread(&block_size,sizeof(int),1,f_gadget);

    fread(&block_size,sizeof(int),1,f_gadget);
     for(i=0 ; i<all_part_num ; i++)
      {
       fread(&id,sizeof(int),1,f_gadget);
       p_snapshot[i].id = id;
      }
    fread(&block_size,sizeof(int),1,f_gadget);
  
    if(block_size != sizeof(int)*all_part_num)
      {
      printf("Error detected in reading the ID block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
      }  

    sort_gadget(0,all_part_num-1);
    make_sorted_gadget_file(filename_sort);
  }                                                                  /* end first else */
#else
    file_open_first();                                               /* open gadget file first */

    fread(&block_size,sizeof(int),1,f_gadget);
    fread(&g_head,256,1,f_gadget);
    fread(&block_size,sizeof(int),1,f_gadget);
    if(block_size != 256)
     {
      printf("Error detected in reading the head block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
     }
#ifdef CAL_ANGULAR_MOMENTUM
    fread(&block_size,sizeof(int),1,f_gadget_vel);
    fread(&g_head,256,1,f_gadget);
    fread(&block_size,sizeof(int),1,f_gadget_vel);
    if(block_size != 256)
     {
      printf("Error detected in reading the head block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
     }
     for(i=0 ; i<g_head.Npart[1] ; i++)
       fread(pos, 3*sizeof(float), 1, f_gadget_vel);
       
    fread(&block_size,sizeof(int),1,f_gadget_vel);
    fread(&block_size,sizeof(int),1,f_gadget_vel);           
#endif   
    halo_red = g_head.Redshift;
    all_part_num = g_head.Npart[1];                            /* initialize these two values */
    rho_crit = crit_density();                                 /* calculate the critical density of the universe with Gadget unit */  

    if((p_snapshot = (PPT_PART_SNAPSHOT) calloc(all_part_num,sizeof(PART_SNAPSHOT))) == NULL)
      {
      printf("Failed to allocate memory for the Gadget file %s!\n",GADGET_FILE);
      exit(0);
      }

    fread(&block_size,sizeof(int),1,f_gadget);
    for(i=0 ; i<all_part_num ; i++)
      {
      fread(pos,3*sizeof(float),1,f_gadget);
      p_snapshot[i].pos[0] = pos[0];
      p_snapshot[i].pos[1] = pos[1];
      p_snapshot[i].pos[2] = pos[2];
      p_snapshot[i].mass = g_head.Massarr[1];
#ifdef CAL_ANGULAR_MOMENTUM
      fread(vel, 3*sizeof(float), 1, f_gadget_vel);
      p_snapshot[i].vel[0] = vel[0];
      p_snapshot[i].vel[1] = vel[1];
      p_snapshot[i].vel[2] = vel[2];
#endif 
      }
    fread(&block_size,sizeof(int),1,f_gadget);
    if(block_size != 3*sizeof(float)*all_part_num)
      {
      printf("Error detected in copying the position block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
      }

    fread(&block_size,sizeof(int),1,f_gadget);
    for(i=0 ; i<all_part_num ; i++)
      fread(pos,3*sizeof(float),1,f_gadget);
    fread(&block_size,sizeof(int),1,f_gadget);

    fread(&block_size,sizeof(int),1,f_gadget);
     for(i=0 ; i<all_part_num ; i++)
      {
       fread(&id,sizeof(int),1,f_gadget);
       p_snapshot[i].id = id;
      }
    fread(&block_size,sizeof(int),1,f_gadget);
  
    if(block_size != sizeof(int)*all_part_num)
      {
      printf("Error detected in reading the ID block of Gadget file %s!\n",GADGET_FILE);
      exit(0);
      }  

    sort_gadget(0,all_part_num-1);
#endif  
  
  printf("Gadget file loaded.\n\n");
}                                                                    /* end load_snapshot */

double crit_density()
{
 double a;
 a = 1.0 / (1.0 + halo_red);

 return a * a * a * (g_head.Omega0 / a / a /a + g_head.OmegaLambda) * 2.7754e-8; 
                                                                     /* critical comving density in Gadget unit */

}                                                                    /* end crit_density */

void sort_gadget(int l, int r)
{
  int i = l;
  int j = r;
  unsigned int key = p_snapshot[(i+j)/2].id;
  while(i < j)
    {
      for(;(i<r)&&(p_snapshot[i].id < key);i++);
      for(;(j>l)&&(p_snapshot[j].id > key);j--);
      if(i <= j)
        {
         swap_gadget(&p_snapshot[i],&p_snapshot[j]);
         i++;
         j--;
         }
      }
   if(i<r)
      sort_gadget(i,r);
   if(j>l)
      sort_gadget(l,j);

}                                                                       /* end sort_gadget */

void swap_gadget(PPT_PART_SNAPSHOT a , PPT_PART_SNAPSHOT b)
{
 int i;
 PART_SNAPSHOT temp;

 temp.id = b->id;
 for(i=0; i<3; i++)
   {
   temp.pos[i] = b->pos[i];
#ifdef CAL_ANGULAR_MOMENTUM
   temp.vel[i] = b->vel[i]; 
#endif     
   }
 temp.mass = b->mass;                                                     /* store b */
  
 b->id = a->id;
 for(i=0 ; i<3; i++)
   {
   b->pos[i] = a->pos[i];
#ifdef CAL_ANGULAR_MOMENTUM
   b->vel[i] = a->vel[i];  
#endif    
   }
 b->mass = a->mass;                                                       /* b = a*/

 a->id = temp.id;
 for(i=0 ; i<3; i++)
   {
   a->pos[i] = temp.pos[i];
#ifdef CAL_ANGULAR_MOMENTUM
   a->vel[i] = temp.vel[i]; 
#endif      
   }
 a->mass = temp.mass;                                                     /* a = temp */
 
}                                                                         /* end swap_gadget */

void load_sorted_gadget_file()
{ 
  int i;

  fread(&g_head,256,1,f_sorted_gadget);

  halo_red = g_head.Redshift;
  all_part_num = g_head.Npart[1];                            /* initialize these two values */
  rho_crit = crit_density();                                 /* calculate the critical density of the universe with Gadget unit */ 

  if((p_snapshot = (PPT_PART_SNAPSHOT) calloc(all_part_num,sizeof(PART_SNAPSHOT))) == NULL)
    {
    printf("Failed to allocate memory for the Gadget file %s!\n",GADGET_FILE);
    exit(0);
    }

  for(i=0 ; i<all_part_num; i++)
   fread(&p_snapshot[i], sizeof(PART_SNAPSHOT), 1, f_sorted_gadget);

  fclose(f_sorted_gadget);
}                                                                         /* load_sorted_gadget_file */

void make_sorted_gadget_file(char *filename_sort)
{
 int i;
 f_sorted_gadget = fopen(filename_sort,"wb");

 fwrite(&g_head,256,1,f_sorted_gadget);
 
 for(i=0; i<all_part_num; i++)
  fwrite(&p_snapshot[i], sizeof(PART_SNAPSHOT), 1, f_sorted_gadget);

 fclose(f_sorted_gadget);
}                                                                          /* make_sorted_gadget_file */
