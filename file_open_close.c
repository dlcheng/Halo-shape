#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void file_open_first()
{
 char file_string[500];
	
 sprintf(file_string,"%s%s",INPUT_GADGET_FOLDER,GADGET_FILE);
 if((f_gadget = fopen(file_string,"rb")) == NULL)
   {
   printf("Can't open file %s, we have to stop!\n",file_string);
   exit(0);
   }

#ifdef CAL_ANGULAR_MOMENTUM
   f_gadget_vel = f_gadget;
#endif  
}                                                    /* end file_open_first, fist open the gadget file */


void file_open_output()
{
 char file_string[500];
 sprintf(file_string,"%s%s%.3f%s",OUTPUT_FOLDER,"halo_shape_z=",halo_red,".txt");
 if((f_output = fopen(file_string,"w+")) == NULL)
   {
   printf("Can't open file %s, we have to stop!\n",file_string);
   exit(0);
   }

}                                                    /* end file_open_then */

void file_close_input()
{
 fclose(f_ahf_1);
 fclose(f_ahf_2);
}                                                    /* end file_close_input */

void file_close_output()
{
 fclose(f_output);
}                                                    /* end file_close */
