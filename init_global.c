#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "allvars.h"
#include "proto.h"

void init_global()
{ 
 min_link_num = MIN_LINK_NUM;
 link_num_step = LINK_NUM_STEP;
 link_num_step_num = LINK_NUM_STEP_NUM;
 max_over_den = MAX_OVER_DEN;
 over_dec_fac = OVER_DEC_FAC;
 mass_shell_num = MASS_SHELL_NUM;
 min_halo_part_num = MIN_HALO_PART_NUM;
 file_id = 0;
 prcessed_halo_num = 0;

#ifdef FIX_LINK_NUM
 link_num_step = 0;
 link_num_step_num = 0;
#endif
 
 increase_stack_size();                                           /* increase the stack size to 64M */

 if((min_link_num + link_num_step * link_num_step_num)> NEIGHBORS_NUM)
   {
   printf("The maximum Link number (LINK_NUM) has exceed the maximum, try some value <%d instead!\n",NEIGHBORS_NUM);
   exit(0);
   }
   
/*
#ifdef SPECIFY_HALO_ID_RANGE_CASE
 if(SPECIFY_HALO_ID_RANGE_LEFT > SPECIFY_HALO_ID_RANGE_RIGHT)
   {
   printf("Wrong ID range, we have to stop!\n");
   exit(0);
   }
#endif
*/

#ifdef SPECIFY_HALO_MASS_RANGE_CASE
 if(SPECIFY_HALO_MASS_RANGE_LEFT > SPECIFY_HALO_MASS_RANGE_RIGHT)
   {
   printf("Wrong MASS range, we have to stop!\n");
   exit(0);
   }
#endif

 load_sort_snapshot();                                           /* load gadget file to memory and sort the ID */
 file_open_output();

#ifdef SPECIFY_HALO_MASS_RANGE_CASE
 halo_part_left = (int) (SPECIFY_HALO_MASS_RANGE_LEFT / 1e10 / g_head.Massarr[1]);
 halo_part_right = (int) (SPECIFY_HALO_MASS_RANGE_RIGHT / 1e10 / g_head.Massarr[1]);

 if(halo_part_left < min_halo_part_num)
   {
   printf("The lower end of mass range is smaller than MIN_HALO_PART_NUM\n");
   printf("Reset MIN_HALO_PART_NUM to 0.\n");
   min_halo_part_num = 0;
   }
#endif

 halo_cen.id = -1;                                                /* the center of the halo is marked "-1", while usual particles will
                                                                     start from 0 */
 halo_cen.mass = 0.0;                                             /* no mass for this virtual particle */
 

#ifdef SPECIFY_HALO_MASS_RANGE_CASE 
 int i,j;
 for(i=0; i<MASS_BIN_NUM; i++)
   for(j=0; j<47; j++)
     mass_bin_data[i][j] = 0.0;
#endif

}                                                                 /* end init_global */


void increase_stack_size()
{
  const rlim_t kStackSize = 64L * 1024L * 1024L;                  /* min stack size = 64 Mb */
  struct rlimit rl;
  int result;
  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
    {
     if (rl.rlim_cur < kStackSize)
      {
      rl.rlim_cur = kStackSize;
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)
       fprintf(stderr, "setrlimit returned result = %d\n", result);
      }
    }
}                                                                 /* end increase_stack_size */
