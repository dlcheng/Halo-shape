#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void det_mass_shell(PPT_PART_HALO a)
{
 int i;

 if(a->id == -1)                                                     /* is this the central virtual particle */
   for(i=0 ; i < link_num ; i++)
     det_mass_shell(a->ngb[i].p_ngb);                                /* visiting all the branchs from the center particle */

 else
 {
  if(a->flag[0] == 0)                                                /* open for visiting */
    {
    a->flag[0] = 1;                                                  /* marked as not acessable again */

     if(a->density >= (1.0 - MASS_SHELL_APT) * over_den * rho_crit)
       {
       a->flag[3] = 1;                                                /* this particle is inside the mass shell */
       if(a->density <= (1.0 + MASS_SHELL_APT) * over_den * rho_crit)
         a->flag[2] = 1;                                              /* particles are identified in the mass shell */
       for(i=0 ; i < link_num ; i++)
         det_mass_shell(a->ngb[i].p_ngb);                             /* visiting all the branchs */
       }
     else
      a->flag[1] = 1;                                                 /* these particles are defined as boundary particle, but actually
                                                                         not used */
    }

 }
}                                                                     /* end det_mass_shell */
