#ifndef GEMINI3D_H
#define GEMINI3D_H

#ifdef __cplusplus
extern "C" {
#endif

enum { LMAX = 1000 };

struct params {
  // order and lengths must match in Fortran and C
  // see gemini_main.f90 "cparams"
  bool fortran_nml;
  bool fortran_cli;
  bool debug;
  bool dryrun;
  char out_dir[LMAX];
  // .ini [base]
  int ymd[3];
  float UTsec0;
  float tdur;
  float dtout;
  float activ[3];
  float tcfl;
  float Teinf;
  // .ini
};

// cfgC for now is just passed around as void**, to modify it would have to be struct *


/* This is now housed in the main C program */
extern int gemini_main(struct params *, int*, int*);

extern void help_gemini_bin();

/* interfaces to libgemini */
extern void cli_config_gridsize_C(struct params*, int*, int*, void**);
extern void get_fullgrid_size_C(int*, int*, int*);
extern void get_config_vars_C(void**, bool*, int*, double*, double*);
extern void get_subgrid_size_C(int*, int*, int*);
extern void get_species_size_C(int*);
extern void gemini_alloc_C(void**,double**, double**, double**,void**);
//extern void memblock_from_C(double**, double**, double**);
extern void set_start_values_C(int*, double*, double*, double*, double*);
extern void msisinit_C();
extern void init_neutralBG_C(double*, double*, int*, double*, double*, double*);
extern void set_update_cadence_C(int*);
extern void neutral_atmos_winds_C(int*, double*);
extern void check_finite_output_C(double*);
extern void gemini_dealloc_C(double**, double**, double**);
extern void get_solar_indices_C(void**, double*, double*);
extern void v12rhov1_C();
extern void T2rhoe_C();
extern void interface_vels_allspec_C(int*);
extern void sweep3_allparams_C(double*);
extern void sweep1_allparams_C(double*);
extern void sweep2_allparams_C(double*);
extern void rhov12v1_C();
extern void clean_param_C(int*, int*, void**, double**);
extern void VNRicht_artvisc_C();
extern void compression_C(double*);
extern void rhoe2T_C();
extern void energy_diffusion_C(double*);
extern void source_loss_allparams_C(double*, double*, int*, double*, double*, double*, bool*, double*, double*);
extern void check_dryrun_C(void**);

extern void dateinc_C(double*, int*, double*);

/* interfaces for libgemini_mpi */
// some of these will very likely need to be rewritten when used with forestclaw
extern void mpisetup_C();
extern void mpiparms_C(int*, int*);
extern void init_procgrid_C(int*, int*, int*, int*);
extern void read_grid_C(void**, int*, void**);
extern void outdir_fullgridvaralloc_C(void**, void**, int*, int*, int*);
extern void get_initial_state_C(void**, double**, double**, void**, int*, void**, double*, int*, double*);
extern void init_Efieldinput_C(double*, double*, int*, double*);
extern void pot2perpfield_C();
extern void BGfield_Lagrangian_C(void**, int*, void**, double**, void**, double*, double*);
extern void init_precipinput_C(double*, double*, int*, double*);
extern void init_neutralperturb_C(double*, int*, double*);
extern void get_initial_drifts_C(void**, int*, void**, double**, double**, double**, void**);
extern void neutral_atmos_wind_update_C(double*, double*);
extern void neutral_perturb_C(double*, double*, int*, double*, double*, double*);
extern void electrodynamics_C(int*, double*, double*, int*, double*);
extern void get_gavg_Tinf_C(void**, double*, double*);
extern void clear_dneu_C(void**);
extern void halo_interface_vels_allspec_C(int*);
extern void set_global_boundaries_allspec_C(int*);
extern void halo_allparams_C(int*, void**, double**, double**);
extern void RK2_prep_mpi_allspec_C(int*, void**, double**);
extern void dt_select_C(int*, double*, double*, double*, double*);
extern void check_fileoutput_C(double*, double*, double*, double*, int*, int*, double*);

#ifdef __cplusplus
}
#endif

#endif
