/* checkpoint.h */
/* Wolfgang Tichy 11/2008 */


/* checkpoint.c */
char *checkpoint_filename(char *suffix);
int checkpoint_checkforfiles(char *suffix);
int checkpoint_checkifrestart(tGrid* grid);
FILE *checkpoint_openfiles(char *suffix);
tVarList *checkpoint_varlist(tGrid *grid);
void checkpoint_copy_output();

void checkpoint_read(tGrid *grid);
void checkpoint_read_ParsAndIterations_local(tGrid *grid, FILE *fp);
void checkpoint_read_Vars_local(tGrid *grid, FILE *fp);

void checkpoint_write(tGrid *grid);
void checkpoint_write_local(tGrid *grid, FILE *fp);

/* fake the mpi calls this module makes */
int sgrid_mpi_size(void);
int sgrid_mpi_rank(void);
int sgrid_mpi_allreduce_sum_int(int n);
void sgrid_mpi_barrier(void);
#define processor0 1
