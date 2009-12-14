/* main.h */
/* Wolfgang Tichy, April 2005 */

void parse_parameter_file(char *parfile);
int iterate_parameters(void);

int read_command_line(int argc, char **argv);
int parse_command_line_options();
int make_output_directory();
int initialize_grid(tGrid *g);
int evolve_grid(tGrid *grid);
int finalize_grid(tGrid *g);


