/* main.h */
/* Wolfgang Tichy, April 2005 */

void parse_parameter_file(char *parfile);
int iterate_parameters(int next);

int get_next_argument_line(int argc, char **argv, char *buffer,
                           int *largc, char ***largv, int *lnum);
int read_command_line(int argc, char **argv);
int parse_command_line_options(void);
int make_output_directory(void);
int initialize_grid(tGrid *g);
int evolve_grid(tGrid *grid);
int finalize_grid(tGrid *g);


/* skeleton.c */
void remove_all_Funs(void);

/* variables.c */
void free_global_vdb(void);

/* parameters.c */
void free_global_parameter_database(void);
