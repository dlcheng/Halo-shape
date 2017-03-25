#------------------------------------------Always recommended 
OPT += -DFIX_LINK_NUM                      # fix the link number
OPT += -DID_START=0                        # ID starts from 0
#OPT += -DENABLE_OMP
OPT += -DCAL_ANGULAR_MOMENTUM              # calculate the angular momentum of the halo
#------------------------------------------The version of the AHF
OPT += -DAHF_V_1_0                         # if AHF version >= 1.0

#------------------------------------------Multi mass particles
#OPT += -DMULTI_MASS                       #coming soon

#------------------------------------------Origin choice
#OPT += -DHALO_CENTER_AS_ORIGN
#OPT += -DSHELL_MASS_CENTER_AS_ORIGIN
OPT += -DGRID_BODY_MASS_CENTER_AS_ORIGIN

#-----------------------------------------Which inertia tensor
#OPT += -DSHELL_INERTIA_TENSOR
OPT += -DGRID_BODY_INERTIA_TENSOR

#-----------------------------------------Analysize sequence, choose only one of them
#OPT += -DNORMAL_MODE                     #from beginning to the end
OPT += -DSPECIFY_HALO_MASS_RANGE_CASE     #only consider halos within the mass range, recommended

#-----------------------------------------Additional output choice
#OPT += -DOUTPUT_WHOLE_HALO               #write the whole halo to a file
#OPT += -DOUTPUT_PARTICLES_INFO           #write additional files containing particles' info on the mass shell and inside 
#OPT += -DSTORE_SORTED_GADGET_FILE        #store the sorted gadget file
OPT += -DOUTPUT_NEIGHBORS_RELATION        #write a file contain the neighbor relations

#--------------------------------------- Select target computer

SYSTYPE="dlcheng"
#SYSTYPE="ITSC"

#--------------------------------------- Adjust settings for target computer

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =    -O3 -Wall
GSL_INCL =  -I/home/dalong/Install/gsl/include
GSL_LIBS =  -L/home/dalong/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/gsl-1.14/lib  -Wl,"-R /usr/local/gsl-1.14/lib" 
GSL_INCL =  -I/usr/local/gsl-1.14/include
endif



OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   = AHF_Halo_Shape+

OBJS   = allvars.o analysis_mass_bin.o construct_ahf_halo.o density.o determine_mass_shell.o \
         file_open_close.o init_each_file.o init_global.o inertia_tensor.o \
         link_list.o load_gadget_file.o main.o process_halo.o write_file.o

INCL   = allvars.h proto.h define.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) -fopenmp


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm -fopenmp

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)  -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)


