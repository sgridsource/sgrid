/* sgrid_GridIterators.c */
/* Wolfgang Tichy 8/2008 */

#include "sgrid.h"
#include "GridIterators.h"



int sgrid_GridIterators() 
{
  printf("Adding GridIterators\n");

  /* functions */
  AddFun(PRE_GRID, Init_Newton, "initialize Newton");

  /* global parameter */
  AddPar("GridIterators_verbose", "yes", "talk about it");
  AddPar("GridIterators_GMRES_restart", "max", 
         "Restart parameter for GMRES[some #,max");
  AddPar("GridIterators_setABStozero_below", "0", 
         "some numbers will be set equal to 0 if their magnitude "
         "is below this value [any positive number]");
  AddPar("GridIterators_Newtonstep", "full",
         "how we take Newton steps [full,backtrack]");
  AddPar("GridIterators_templates_RESID_mode", "tol/norm(b)",
         "how we compute input RESID from tol [tol/norm(b),tol]");
  AddPar("GridIterators_UMFPACK_version", "di",
         "umfpack version to be called [di,dl]");
  AddPar("GridIterators_templates_as_Preconditioner", "GMRES",
         "which templates wrapper we use as precon [GMRES,BICGSTAB,CGS]");
  AddPar("GridIterators_Preconditioner_type", "fd",
         "whether we use spectral or finite diff [spectral,fd]");
  AddPar("GridIterators_Preconditioner_itmax", "1000",
         "max number of iterations in Preconditioner");
  AddPar("GridIterators_Preconditioner_reltol", "0.001",
         "relative tol in Preconditioner");
  AddPar("GridIterators_SOR_omega", "1.5",
         "omega par in SOR algorithm, we like: 0<omega<2, omega=1 is Gauss-Seidel");
  AddPar("GridIterators_SOR_matrix", "SetMatrixLines_slowly",
         "how we set the matrix for SOR "
         "[SetMatrixLines_slowly,SetMatrixLines_forSortedVars_slowly]");

  /* check whether there is more to do */
  if (!Getv("physics", "GridIterators")) return 0;

  /* functions */

  /* variables */
  //AddVar("u", "", "scalar field for Poisson equation");

/* temporary variables, should be inside if statement to avoid
     cluttering the variable name space 
  AddVar("GridIterators_p",  "", "temporary variable");
  AddVar("GridIterators_ph", "", "temporary variable");
  AddVar("GridIterators_rt", "", "temporary variable");
  AddVar("GridIterators_s",  "", "temporary variable");
  AddVar("GridIterators_sh", "", "temporary variable");
  AddVar("GridIterators_t",  "", "temporary variable");
  AddVar("GridIterators_v",  "", "temporary variable");
*/
  return 0;
}
