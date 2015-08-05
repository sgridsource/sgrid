/* sgrid_MPI.c */
/* Wolfgang Tichy, August 2015 */
/* a few wrappers for MPI */

#include "sgrid.h"


/* include files for MPI */
#ifdef USEMPI
#include <mpi.h>
#endif


/* Wrappers for MPI_Init and MPI_Finalize */
int sgrid_MPI_Init(int *pargc, char ***pargv)
{
#ifdef USEMPI
  return MPI_Init(pargc, pargv);
#endif
  return 0;
}
int sgrid_MPI_Finalize(void)
{
#ifdef USEMPI
  return MPI_Finalize();
#endif
  return 0;
}

/* return MPI rank, if MPI is not compiled in return 0 */
int sgrid_MPI_rank(void)
{
  int rank=0;
#ifdef USEMPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  return rank;
}

/* number of processes we are running with, or 1 if no MPI */
int sgrid_MPI_size(void)
{
  int size=1;
#ifdef USEMPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  return size;
}

/* barrier */
int sgrid_MPI_barrier(void)
{
  int ret=0;
#ifdef USEMPI
  ret = MPI_Barrier(MPI_COMM_WORLD);
#endif
  return ret;
}
