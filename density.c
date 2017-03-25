#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void cal_halo_part_density()
{
  int i,j;
  double h;                                                        /* the distance to the 32th particle */
  
  for(i=0 ; i<halo_part_num ; i++)
   {
    h = p_halo[i].ngb[NEIGHBORS_NUM-1].d;
    p_halo[i].density = p_halo[i].mass * sph_smooth_kernel(0.0, h); /* density caused by the particle itself */
    for(j=0 ; j<NEIGHBORS_NUM ; j++)
      p_halo[i].density += p_halo[i].ngb[j].p_ngb->mass * sph_smooth_kernel(p_halo[i].ngb[j].d,h);
   }

}                                                                   /* end cal_halo_part_density */

double sph_smooth_kernel(double r, double h)
{
 if(r >= 0.0 && r <= 0.5 * h)
   return (1.0 - 6.0 * (r/h) * (r/h) + 6.0 * (r/h) * (r/h) * (r/h)) * 8.0 / PI / h / h / h;

 if(r > 0.5 * h && r <=h)
   return (2.0 * ( 1.0 - r/h) * ( 1.0 - r/h))* 8.0 / PI / h / h / h;
 
   return 0.0;                                                      /*else return 0 */
}                                                                    /* end sph_smooth_kernel */
