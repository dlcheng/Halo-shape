#ifndef ALLVAR_H
 #include "allvars.h"
#endif

void increase_stack_size();

void init_global();
void process_exist_halo();

void file_open_first();
void file_open_output();
void file_close_input();
void file_close_output();

int init_each_input_file();

void load_sort_snapshot();
double crit_density();
void sort_gadget(int l, int r);
void swap_gadget(PPT_PART_SNAPSHOT a, PPT_PART_SNAPSHOT b);
void load_sorted_gadget_file();
void make_sorted_gadget_file(char *filename);

int construct_ahf_halo();
void re_scale_center();
void init_halo();
void init_single_halo_part(PPT_PART_HALO a);


void make_link_list();
void link_center();
void link_halo_particles();
void make_link_from_table(FILE *link_table);
void fill_link_list(int a, PPT_LINK_LIST lt);
void sort_link_list(PPT_LINK_LIST lt);
void get_neig_info(int a,PPT_LINK_LIST lt);
double distance(float * a, float *b);
void output_neigh_file(PPT_LINK_LIST lt);


void cal_halo_part_density();
double sph_smooth_kernel(double r, double h);


void det_mass_shell(PPT_PART_HALO a);


int mass_shell_info();
void get_inertia_tensor();
void eigen_inertia_tensor();
void get_halo_shape();
void copy_center(float * a, float *b );
void relate_center(float * a, float *b , float *c);

void write_shape_file();
void write_particle_file();
void write_whole_halo_particle();

void process_exist_halo();
void refresh_flags();
void update_link_num_mass_shell(int link_step, int mass_shell);

#ifdef SPECIFY_HALO_MASS_RANGE_CASE 
void collect_mass_bin_data();
void mass_bin_data_out_put();
#endif
