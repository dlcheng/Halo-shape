#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef ENABLE_OMP
#include <omp.h>
#endif
#include "allvars.h"
#include "proto.h"


void make_link_list()
{
  link_center();
  link_halo_particles();
}                                                                      /* end link_list */


void link_center()
{
  int i;
  PPT_LINK_LIST lt;
  
  if((lt = (PPT_LINK_LIST) calloc(halo_part_num, sizeof(LINK_LIST))) == NULL)
   {
   printf("Fail to allocate memory to construct the LINK_LIST of halo center for halo ID:%d at (%.3e,%.3e,%.3e)[kpc/h]\n",
           halo_id,halo_cen.pos[0],halo_cen.pos[1],halo_cen.pos[2]);
   exit(0);
   }
  
  for(i=0 ; i<halo_part_num ; i++)
   {
   lt[i].d = distance(halo_cen.pos,p_halo[i].pos);
   lt[i].place = i;
   }
   
  LINK_LIST temp;
  int j;

  for(i=0 ; i<NEIGHBORS_NUM ; i++)
    for(j = i+1 ; j< halo_part_num ; j++)
      if(lt[j].d < lt[i].d)
         {
          temp.d = lt[i].d;
          temp.place = lt[i].place;
          lt[i].d = lt[j].d;
          lt[i].place = lt[j].place;
          lt[j].d = temp.d;
          lt[j].place = temp.place;
         }
 
  for(i=0 ; i<NEIGHBORS_NUM ; i++)
    {
    halo_cen.ngb[i].d = lt[i].d;
    halo_cen.ngb[i].p_ngb = &p_halo[lt[i].place];
    }

  free(lt);
}                                                                               /* end link_center */


void link_halo_particles()
{
  int i;
  PPT_LINK_LIST lt;
  FILE *link_table;
  char table_name[500];

  sprintf(table_name,"%s%s%d%s%.3e%s%.3e%s%.3e%s%s",
          NEIGHBORS_FOLDER,"halo_",halo_part_num,"_(",halo_cen.pos[0],",",
          halo_cen.pos[1],",",halo_cen.pos[2],")","_neighbors.dat");
  if((link_table = fopen(table_name,"rb")) != NULL)
    {
    make_link_from_table(link_table);
    printf("Neigbors link relation recovered from %s%d%s%.3e%s%.3e%s%.3e%s%s.\n",
           "halo_",halo_part_num,"_(",halo_cen.pos[0],",",halo_cen.pos[1],",",
           halo_cen.pos[2],")","_neighbors.dat");
    fclose(link_table);
    }
  else
    {
    printf("Find no the neigbors link relation file.\n");
#ifdef OUTPUT_NEIGHBORS_RELATION
    char file_string[500];

    if(halo_part_num >= 2000)  /* only write a neighbors file if halo with particles >=2000  */
      {
      printf("Will generate new one.\n");    
      sprintf(file_string,"%s%s%d%s%.3e%s%.3e%s%.3e%s%s",
              NEIGHBORS_FOLDER,"halo_",halo_part_num,"_(",
              halo_cen.pos[0],",",halo_cen.pos[1],",",halo_cen.pos[2],")","_neighbors.dat");  
      if((f_my_neigh = fopen(file_string,"wb")) == NULL)
        {
        printf("Can't write file %s, we have to stop!\n",file_string);
        exit(0);
        }
      }
#endif
  
    if((lt = (PPT_LINK_LIST) calloc(halo_part_num -1 , sizeof(LINK_LIST))) == NULL)
     {
     printf("Fail to allocate memory to construct the LINK_LIST for halo ID:%d at (%.3e,%.3e,%.3e)[kpc/h]\n",
             halo_id,halo_cen.pos[0],halo_cen.pos[1],halo_cen.pos[2]);
     exit(0);
     }
  
    for(i=0 ; i<halo_part_num ; i++)
      {
      fill_link_list(i,lt);
      sort_link_list(lt);
      get_neig_info(i,lt);
#ifdef OUTPUT_NEIGHBORS_RELATION
      if(halo_part_num >= 2000)
        output_neigh_file(lt);
#endif
      }
      
#ifdef OUTPUT_NEIGHBORS_RELATION
    if(halo_part_num >= 2000)
      fclose(f_my_neigh);
#endif
    free(lt);
    }                                                                           /* end first else */

}                                                                               /* end link_halo_particles */

void make_link_from_table(FILE * link_table)
{
 int i,j;
 int ngb_id;
 for(i=0 ; i<halo_part_num ; i++)
   for(j=0 ; j<NEIGHBORS_NUM ; j++)
     {
     fread(&ngb_id,sizeof(int),1,link_table);
     p_halo[i].ngb[j].p_ngb = &p_halo[ngb_id];
     p_halo[i].ngb[j].d = distance(p_halo[i].pos,p_halo[ngb_id].pos);
     }
}                                                                               /* end make_link_from_table */

void fill_link_list(int a, PPT_LINK_LIST lt)
{
#ifdef ENABLE_OMP
  int i;
  int thread_count[1000];          /* large enough so thart no node can have 1000 cores */
  int thread_id;
  PPT_LINK_LIST lt_temp;
  
#pragma omp parallel shared(p_halo, a, thread_count, lt, halo_part_num) private(thread_id, lt_temp, i)
{ 
	  
  int place_start = 0;
  int num;
  thread_id = omp_get_thread_num();
  lt_temp = (PPT_LINK_LIST) calloc(halo_part_num -1 , sizeof(LINK_LIST));
  thread_count[thread_id] = 0;
	
  #pragma omp for 	
  for(i=0; i<halo_part_num ; i++)
    if(i != a)
    {
     lt_temp[thread_count[thread_id]].d = distance(p_halo[a].pos,p_halo[i].pos);
     lt_temp[thread_count[thread_id]].place = i;
     thread_count[thread_id]++;
    }
      
  for(num=0; num < thread_id; num++)
    place_start += thread_count[num];
    
  for(num=0; num < thread_count[thread_id]; num++)
    {
	lt[place_start + num].d = lt_temp[num].d;
	lt[place_start + num].place = lt_temp[num].place;	
    }  
    
  free(lt_temp);
}

#else
  int i;
  int j = 0;
  for(i=0; i<halo_part_num ; i++)
    if(a != i)
    {
     lt[j].d = distance(p_halo[a].pos,p_halo[i].pos);
     lt[j].place = i;
     j++;
    }
#endif
    
}                                                                               /* end fill_link_list */

void sort_link_list(PPT_LINK_LIST lt)
{
 int i,j;
 LINK_LIST temp;

 for(i = 0 ; i<NEIGHBORS_NUM ; i++)
   for(j = i+1 ; j< (halo_part_num - 1) ; j++)
     if(lt[j].d < lt[i].d)
        {
         temp.d = lt[i].d;
         temp.place = lt[i].place;
         lt[i].d = lt[j].d;
         lt[i].place = lt[j].place;
         lt[j].d = temp.d;
         lt[j].place = temp.place;
        }

}                                                                               /* end sort_link_list */


void get_neig_info(int a,PPT_LINK_LIST lt)
{
  int i;
  for(i=0 ; i<NEIGHBORS_NUM ; i++)
    {
    p_halo[a].ngb[i].d = lt[i].d;
    p_halo[a].ngb[i].p_ngb = &p_halo[lt[i].place];
    }

}                                                                               /* end get_neig_info */

double distance(float * a, float *b)
{
  int i;
  double r = 0.0;
  for(i=0 ; i<3 ; i++)
    r += (a[i] - b[i]) * (a[i] - b[i]);

  return sqrt(r);
}                                                                               /*end distance */

void output_neigh_file(PPT_LINK_LIST lt)
{
 int i;
#ifdef OUTPUT_NEIGHBORS_RELATION
 for(i=0 ; i<NEIGHBORS_NUM ; i++)
   fwrite(&lt[i].place,sizeof(int),1,f_my_neigh);
#endif
}                                                                               /* end output_neigh_file */
