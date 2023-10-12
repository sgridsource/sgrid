/* rename_symbols.h */
/* Wolfgang Tichy, July 2022 */


/* In this header file we can put defines that rename certain symbols.
   When sgrid is compiled into a library this is necessary to avoid
   conflicts. */

/**********************************************/
/* symbols in sgrid core and output functions */
/**********************************************/
/* main/main */
#define free_everything SGRID_free_everything
#define initialize_libraries SGRID_initialize_libraries
#define read_command_line SGRID_read_command_line
#define make_output_directory SGRID_make_output_directory
#define parse_command_line_options SGRID_parse_command_line_options
#define findparameterindex SGRID_findparameterindex
#define findparameter SGRID_findparameter
#define translatevalue SGRID_translatevalue
#define set_numericalvalue_byIndex SGRID_set_numericalvalue_byIndex
#define set_booleanvalue_byIndex SGRID_set_booleanvalue_byIndex
#define makeparameter SGRID_makeparameter
#define npdbmax SGRID_npdbmax
#define setparameter SGRID_setparameter
#define parse_parameter_file SGRID_parse_parameter_file
#define printparameter SGRID_printparameter
#define printparameters SGRID_printparameters
#define iterate_parameters SGRID_iterate_parameters
#define print_pdb_i1_i2 SGRID_print_pdb_i1_i2
#define print_parameter_database SGRID_print_parameter_database
#define create_copy_of_pdb1_in_pdb2 SGRID_create_copy_of_pdb1_in_pdb2
#define make_empty_pdb SGRID_make_empty_pdb
#define copy_pdb SGRID_copy_pdb
#define free_pdb_contents SGRID_free_pdb_contents
#define free_pdb SGRID_free_pdb
#define Yo SGRID_Yo
#define prdivider SGRID_prdivider
#define initTimeIn_s SGRID_initTimeIn_s
#define getTimeIn_s SGRID_getTimeIn_s
#define prTimeIn_s SGRID_prTimeIn_s
#define getClockTimeIn_s SGRID_getClockTimeIn_s
#define prClockTimeIn_s SGRID_prClockTimeIn_s
#define min2 SGRID_min2
#define min3 SGRID_min3
#define max2 SGRID_max2
#define max3 SGRID_max3
#define min_in_1d_array SGRID_min_in_1d_array
#define max_in_1d_array SGRID_max_in_1d_array
#define min2_in_1d_array SGRID_min2_in_1d_array
#define max2_in_1d_array SGRID_max2_in_1d_array
#define min3_in_1d_array SGRID_min3_in_1d_array
#define max3_in_1d_array SGRID_max3_in_1d_array
#define finit SGRID_finit
#define trim_whitespace SGRID_trim_whitespace
#define get_par_from_str SGRID_get_par_from_str
#define remove_dir SGRID_remove_dir
#define system_emu SGRID_system_emu
#define unlock_curr_til_EOF SGRID_unlock_curr_til_EOF
#define finalexit SGRID_finalexit
#define print_errno SGRID_print_errno
#define copy_file SGRID_copy_file
#define copy_file_into_dir SGRID_copy_file_into_dir
#define dmalloc SGRID_dmalloc
#define imalloc SGRID_imalloc
#define cmalloc SGRID_cmalloc
#define pmalloc SGRID_pmalloc
#define prvarlist SGRID_prvarlist
#define vlalloc SGRID_vlalloc
#define vlfree SGRID_vlfree
#define vlpushone SGRID_vlpushone
#define vlpush SGRID_vlpush
#define vlpushvl SGRID_vlpushvl
#define vldropone SGRID_vldropone
#define vldrop SGRID_vldrop
#define vldropn SGRID_vldropn
#define vlduplicate SGRID_vlduplicate
#define vlenable SGRID_vlenable
#define vldisable SGRID_vldisable
#define VLPtrEnable1 SGRID_VLPtrEnable1
#define VLDisableFree SGRID_VLDisableFree
#define AddDuplicate SGRID_AddDuplicate
#define AddDuplicateEnable SGRID_AddDuplicateEnable
#define vlsetconstant SGRID_vlsetconstant
#define vlcopy SGRID_vlcopy
#define varcopy SGRID_varcopy
#define vlswap SGRID_vlswap
#define varswap SGRID_varswap
#define vlaverage SGRID_vlaverage
#define vlsubtract SGRID_vlsubtract
#define vladd SGRID_vladd
#define varadd SGRID_varadd
#define vladdto SGRID_vladdto
#define printCI SGRID_printCI
#define printbface SGRID_printbface
#define printbfaces SGRID_printbfaces

/* main/main/utilities.c */
#define construct_argv SGRID_construct_argv
#define system2 SGRID_system2
#define system3 SGRID_system3
#define lock_curr_til_EOF SGRID_lock_curr_til_EOF

/* main/MemoryMan */
#define make_grid SGRID_make_grid
#define add_empty_bface SGRID_add_empty_bface
#define remove_bface SGRID_remove_bface
#define enablevar SGRID_enablevar
#define disablevar SGRID_disablevar
#define enablevarlist SGRID_enablevarlist
#define disablevarlist SGRID_disablevarlist

/* utility/output */
#define output0d_value SGRID_output0d_value
#define fopen_vtk SGRID_fopen_vtk
#define write_raw_vtk_data SGRID_write_raw_vtk_data

/* utility/checkpoint/wolfio.c */
#define fgotonext SGRID_fgotonext
#define fgetparameter SGRID_fgetparameter
#define extract_after_EQ SGRID_extract_after_EQ
#define extrstr_before_after_EQ SGRID_extrstr_before_after_EQ
#define fscanline SGRID_fscanline
#define fscan_str_using_getc SGRID_fscan_str_using_getc
#define fscanf1 SGRID_fscanf1

/* utility/Coordinates */
#define set_consistent_flags_in_all_bfaces SGRID_set_consistent_flags_in_all_bfaces
#define CubedSphere_sigma_AB SGRID_CubedSphere_sigma_AB
#define CubedSphere_sigma SGRID_CubedSphere_sigma
#define xyz_of_lamAB_CubSph SGRID_xyz_of_lamAB_CubSph
#define lamAB_of_xyz_CubSph SGRID_lamAB_of_xyz_CubSph
#define dlamAB_dxyz_CubSph SGRID_dlamAB_dxyz_CubSph
#define r_of_lam_sig0sig1 SGRID_r_of_lam_sig0sig1
#define lam_of_r_sig0sig1 SGRID_lam_of_r_sig0sig1
#define dr_dlam_of_lam_sig0sig1 SGRID_dr_dlam_of_lam_sig0sig1
#define dlam_dr_of_lam_sig0sig1 SGRID_dlam_dr_of_lam_sig0sig1
#define r_dr_dlam_of_lamAB_CubSph SGRID_r_dr_dlam_of_lamAB_CubSph
#define ThetaPhi_of_AB_CubSph SGRID_ThetaPhi_of_AB_CubSph
#define ThetaPhi_dThetaPhidAB_of_AB_CubSph SGRID_ThetaPhi_dThetaPhidAB_of_AB_CubSph
#define rho_of_lam_sig0sig1 SGRID_rho_of_lam_sig0sig1
#define lam_of_rho_sig0sig1 SGRID_lam_of_rho_sig0sig1
#define drho_dlam_of_rho_sig0sig1 SGRID_drho_dlam_of_rho_sig0sig1
#define xyz_of_rhoAB_CubSph SGRID_xyz_of_rhoAB_CubSph
#define rhoAB_of_xyz_CubSph SGRID_rhoAB_of_xyz_CubSph
#define drhoAB_dxyz_CubSph SGRID_drhoAB_dxyz_CubSph
#define XYZ_on_face SGRID_XYZ_on_face
#define set_AB_min_max_from_Din SGRID_set_AB_min_max_from_Din
#define arrange_12CubSph_into_empty_cube SGRID_arrange_12CubSph_into_empty_cube
#define two_full_cubes_touching_at_x0 SGRID_two_full_cubes_touching_at_x0
#define sphere_around_two_full_cubes_touching_at_x0 SGRID_sphere_around_two_full_cubes_touching_at_x0
#define two_spheres_around_two_full_cubes SGRID_two_spheres_around_two_full_cubes
#define cart_partials SGRID_cart_partials

/* utility/evolve */
#define evolve_test_analyze SGRID_evolve_test_analyze

/***********************************/
/* symbols in sgrid projects       */
/***********************************/
/* physics/ADMvars */
#define computeADMconstraints SGRID_computeADMconstraints

/* Projects/EoS_T0/EoS_T0.c */
#define EoS_T0_rho0_P_rhoE_from_hm1 SGRID_EoS_T0_rho0_P_rhoE_from_hm1
#define epsl_of_rho0_rhoE SGRID_epsl_of_rho0_rhoE

/* more from Projects/EoS_T0 */
#define hm1_of_rho0_epsl_P SGRID_hm1_of_rho0_epsl_P
#define h_of_rho0_epsl_P SGRID_h_of_rho0_epsl_P
#define alloc_PwP_globals SGRID_alloc_PwP_globals
#define PwP_rho0 SGRID_PwP_rho0
#define PwP_kappa SGRID_PwP_kappa
#define PwP_n SGRID_PwP_n
#define PwP_k SGRID_PwP_k
#define PwP_q SGRID_PwP_q
#define PwP_P SGRID_PwP_P
#define free_PwP_globals SGRID_free_PwP_globals
#define PwP_select_polytrope_n_kappa_k_of_hm1 SGRID_PwP_select_polytrope_n_kappa_k_of_hm1
#define PwP_n_pieces SGRID_PwP_n_pieces
#define PwP_select_polytrope_n_kappa_k_of_P SGRID_PwP_select_polytrope_n_kappa_k_of_P
#define PwP_polytrope_of_hm1 SGRID_PwP_polytrope_of_hm1
#define PwP_polytrope_rho0_of_hm1 SGRID_PwP_polytrope_rho0_of_hm1
#define PwP_polytrope_P_of_hm1 SGRID_PwP_polytrope_P_of_hm1
#define PwP_polytrope_hm1_of_P SGRID_PwP_polytrope_hm1_of_P
#define PwP_polytrope_rho0_rhoE_of_P SGRID_PwP_polytrope_rho0_rhoE_of_P
#define PwP_test_piecewise_n_kappa_k SGRID_PwP_test_piecewise_n_kappa_k
#define PwP_finalize SGRID_PwP_finalize
