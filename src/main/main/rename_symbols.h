/* rename_symbols.h */
/* Wolfgang Tichy, July 2022 */


/* In this header file we can put defines that rename certain symbols.
   When sgrid is compiled into a library this is necessary to avoid
   conflicts. */

/* sgrid core and output functions */
#define make_grid SGRID_make_grid
#define sgrid_memory_persists SGRID_memory_persists
#define free_everything SGRID_free_everything
#define write_raw_vtk_data SGRID_write_raw_vtk_data

/* main/main/utilities.c */
#define construct_argv SGRID_construct_argv
#define system2 SGRID_system2
#define system3 SGRID_system3
#define lock_curr_til_EOF SGRID_lock_curr_til_EOF

/* utility/checkpoint/wolfio.c */
#define fgotonext SGRID_fgotonext
#define fgetparameter SGRID_fgetparameter
#define extract_after_EQ SGRID_extract_after_EQ
#define extrstr_before_after_EQ SGRID_extrstr_before_after_EQ
#define fscanline SGRID_fscanline
#define fscan_str_using_getc SGRID_fscan_str_using_getc
#define fscanf1 SGRID_fscanf1

/* Projects/EoS_T0/EoS_T0.c */
#define EoS_T0 SGRID_EoS_T0
