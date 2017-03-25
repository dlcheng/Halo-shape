#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int main()
{
 init_global();
 
 while(init_each_input_file())
  {
   while(construct_ahf_halo())
     {
/*		 
#ifdef SPECIFY_HALO_ID_RANGE_CASE
      if(halo_id < SPECIFY_HALO_ID_RANGE_LEFT)
        printf("Do nothing.\n\n");
      else
        process_exist_halo();
#endif
*/
#ifdef SPECIFY_HALO_MASS_RANGE_CASE
      if(halo_part_num > halo_part_right)
        {
        printf("Do nothing.\n\n");
        omitted_halo_num++;
        }
      else
        process_exist_halo();
#endif

#ifdef NORMAL_MODE
      process_exist_halo();
#endif
     }
  file_close_input();
  file_id++;
  }         /* end first while */
  
#ifdef SPECIFY_HALO_MASS_RANGE_CASE 
 mass_bin_data_out_put();
#endif

 free(p_snapshot);
 file_close_output();

 return 0;
}                                                   /* end main */


