/* Constants */
#define MPCTOKPC                1000               /* rescale halo center since the output of AHF is in Mpc/h, while Gadget is Kpc/h */
#define PI                      3.14159265359

/* Parameters */
#define NEIGHBORS_NUM           32                 /* number of neighbors used in calculating the local density */

#define MIN_LINK_NUM            10                 /* number of nearby neighbors are considered in the linking process  */
#define LINK_NUM_STEP_NUM       0                  /* number of steps, if no steps set this to 0 */
#define LINK_NUM_STEP           1                  /* the link number increased in each step */ 

#define MAX_OVER_DEN            1000               /* maximum over density */
#define MASS_SHELL_NUM          1                  /* number of mass shell, starts from 1 */
#define OVER_DEC_FAC            2                  /* over density decrease factor per mass shell */

#define MIN_SHELL_PART_NUM      500                /* the minimum allowed number of particles in the most inside shell, not used yet */
#define MIN_HALO_PART_NUM       100                /* the minimum allowed number of particles in the halo */
#define MASS_SHELL_APT          0.03               /* 0.97~1.03 overdensity region defined as mass shell */

/* Folder paths */
#define INPUT_FOLDER            "./Halo-files/fnl--1000/z-2/"
#define INPUT_GADGET_FOLDER     "                            "
#define OUTPUT_FOLDER           "./Result/fnl--1000/z-2/"
#define NEIGHBORS_FOLDER        "./Halo-files/fnl--1000/z-2/Neighbors-files/"

#define GADGET_FILE             "Gadget-files/snapshot_z_2"  /* the gadget file, in the future multi files support would be added */
#define AHF_FILE_BASE           "B100-512"


/* For some special behavirs */
#ifdef SPECIFY_HALO_MASS_RANGE_CASE 
#define SPECIFY_HALO_MASS_RANGE_LEFT   1e13                                               /* in unit of M_sun / h */
#define MASS_BIN_NUM                   6                                                  /* number of mass bin */
#define SPECIFY_HALO_MASS_RANGE_RIGHT  SPECIFY_HALO_MASS_RANGE_LEFT * pow(2,MASS_BIN_NUM) /* in unit of M_sun / h */
#endif


#ifdef SPECIFY_HALO_ID_RANGE_CASE                  /* only consider halo within the ID range, halo ID start form 1*/
#define SPECIFY_HALO_ID_RANGE_LEFT     1700   
#define SPECIFY_HALO_ID_RANGE_RIGHT    1800
#endif


