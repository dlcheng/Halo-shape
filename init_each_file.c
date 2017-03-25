#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int init_each_input_file()                            /* success return 1 else return 0 */
{
 char file_string[500];
 char a_temp;
 
 sprintf(file_string, "%s%s%s%.4d%s%.3f%s",
         INPUT_FOLDER, AHF_FILE_BASE, ".", file_id, ".z", halo_red, ".AHF_halos");
 if((f_ahf_1 = fopen(file_string,"r")) == NULL)
   {
   if(file_id == 0)
     {
     printf("Can't open file %s\nPls check the setting of input folder.\n",file_string);
     exit(0);
     }
   else
     {
     printf("There are no further more files, total file number %d\n", file_id);
     return 0;
     }
   }

 sprintf(file_string, "%s%s%s%.4d%s%.3f%s",
         INPUT_FOLDER, AHF_FILE_BASE, ".", file_id, ".z", halo_red, ".AHF_particles");
 if((f_ahf_2 = fopen(file_string,"r")) == NULL)
   {
   if(file_id == 0)
     {
     printf("Can't open file %s\nPls check the setting of input files path.\n",file_string);
     exit(0);
     }
   else 
     {
     printf("\n\nThere are no further more files, total file number %d\n", file_id);
     return 0;
     }
   }

 fscanf(f_ahf_2,"%d\n", &total_halo_num_each_file);            /* read the total halo number from the particle file */
 
 a_temp = fgetc(f_ahf_1);
 if(a_temp == 35)
   do{
     a_temp = fgetc(f_ahf_1);
     }while(a_temp != 10);                                    /* read the useless first line if any*/
 else
   rewind(f_ahf_1);
 
 halo_id = 0;                                                 /* halo ID start from 1 */
 prcessed_halo_num_each_file = 0;
 omitted_halo_num = 0;   
 
 printf("Start file %d\n\n", file_id);
      
 return 1; 
}                             /* end init_input_file */

