/* bicgstab.c */
/* Wolfgang Tichy 8/2007 */

/* BiCGSTAB: preconditioned biconjugate gradient stabilized method 
   cmp. "Templates for the solution of linear systems", Barett et al
        http://www.netlib.org/templates					*/

#include "sgrid.h"
#include "GridIterators.h"



int bicgstab(tVarList *x, tVarList *b, tVarList *r, tVarList *c,
	     int itmax, double tol, double *normres,
	     void (*lop)(tVarList *, tVarList *, tVarList *), 
	     void (*precon)(tVarList *, tVarList *, tVarList *))
{
  tGrid *grid = b->grid;
  int boxi;
  double alpha = 0, beta = 0;
  double rho = 0, rho1 = 1, rhotol = 1e-50;
  double omega = 0, omegatol = 1e-50;
  int i, ii, j;
  int pr = Getv("GridIterators_verbose", "yes");

  /* temporary storage */
  tVarList *p  = AddDuplicateEnable(x, "_GridIterators_p");
  tVarList *ph = AddDuplicateEnable(x, "_GridIterators_ph");
  tVarList *rt = AddDuplicateEnable(x, "_GridIterators_rt");
  tVarList *s  = AddDuplicateEnable(x, "_GridIterators_s");
  tVarList *sh = AddDuplicateEnable(x, "_GridIterators_sh");
  tVarList *t  = AddDuplicateEnable(x, "_GridIterators_t");
  tVarList *v  = AddDuplicateEnable(x, "_GridIterators_v");
     
#if 0
  /* hack for scalar only: variables passed to may need 
     lop need boundary info */
  if (VarFallOff(x->index[0])) {
    VarNameSetBoundaryInfo("GridIterators_ph", 0, 1, 0);
    VarNameSetBoundaryInfo("GridIterators_sh", 0, 1, 0);
  }
#endif

  /* check */
  if (pr) printf("bicgstab:  itmax %d, tol %e\n", itmax, tol);

  /* compute initial residual rt = r = b - A x */
  lop(r, x, c);
  for (j = 0; j < r->n; j++)
    forallboxes(grid,boxi)
    {
      tBox *box = grid->box[boxi];      
      double *prt = vlldataptr(rt, box, j);
      double *pr  = vlldataptr(r,  box, j);
      double *pb  = vlldataptr(b,  box, j);
      forallpoints(box,i) prt[i] = pr[i] = pb[i] - pr[i];
    }
  *normres = norm2(r);
  if (pr) printf("bicgstab: %5d  %10.3e\n", 0, *normres);
  if (*normres <= tol) return 0;

  /* cgs iteration */
  for (ii = 0; ii < itmax; ii++) {
    rho = dot(rt, r);
    if (fabs(rho) < rhotol) break;

    /* compute direction vector p */
    if (ii == 0)
      vlcopy(p, r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      for (j = 0; j < r->n; j++)
        forallboxes(grid,boxi)
        {
          tBox *box = grid->box[boxi];
	  double *pp = vlldataptr(p, box, j);
	  double *rp = vlldataptr(r, box, j);
	  double *vp = vlldataptr(v, box, j);
	  forallpoints(box,i) 
	    pp[i] = rp[i] + beta * (pp[i] - omega * vp[i]);
        }
    }

    /* compute direction adjusting vector ph and scalar alpha */
    precon(ph, p, c);
    lop(v, ph, c);
    alpha = rho/dot(rt, v);
    for (j = 0; j < r->n; j++)
      forallboxes(grid,boxi)
      {
        tBox *box = grid->box[boxi];
        double *sp = vlldataptr(s, box, j);
        double *rp = vlldataptr(r, box, j);
        double *vp = vlldataptr(v, box, j);
        forallpoints(box,i)
	  sp[i] = rp[i] - alpha * vp[i];
      }

    /* early check of tolerance */
    *normres = norm2(s);
    if (*normres <= tol) {
      for (j = 0; j < r->n; j++)
        forallboxes(grid,boxi)
        {
          tBox *box = grid->box[boxi];
	  double *xp  = vlldataptr(x, box,   j);
	  double *php = vlldataptr(ph, box,  j);
	  forallpoints(box,i) 
	    xp[i] += alpha * php[i];
        }
      if (pr) printf("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n", 
		     ii+1, *normres, alpha, beta, omega);
      break;
    }

    /* compute stabilizer vector sh and scalar omega */
    precon(sh, s, c);
    lop(t, sh, c);
    omega = dot(t, s) / dot (t, t);

    /* compute new solution approximation */
    for (j = 0; j < r->n; j++)
      forallboxes(grid,boxi)
      {
        tBox *box = grid->box[boxi];
        double *rp = vlldataptr(r, box,  j);
        double *sp = vlldataptr(s, box,  j);
        double *tp = vlldataptr(t, box,  j);
        double *xp = vlldataptr(x, box,  j);
        double *php = vlldataptr(ph, box,  j);
        double *shp = vlldataptr(sh, box,  j);
        forallpoints(box,i) {
	  xp[i] += alpha * php[i] + omega * shp[i];
	  rp[i] = sp[i] - omega * tp[i];
        }
      }

    /* are we done? */
    *normres = norm2(r);
    if (0) printf("bicgstab: %5d  %10.3e\n", ii+1, *normres);
    if (pr) printf("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n", 
		  ii+1, *normres, alpha, beta, omega);
    if (*normres <= tol) break;
    rho1 = rho;
    if (fabs(omega) < omegatol) break;
  }

  /* free temporary storage */
  VLDisableFree(p);
  VLDisableFree(ph);
  VLDisableFree(rt);
  VLDisableFree(s);
  VLDisableFree(sh);
  VLDisableFree(t);
  VLDisableFree(v);

  /* iteration failed */
  if (ii > itmax) return -1;

  /* breakdown */
  if (fabs(rho) < rhotol) return -10;
  if (fabs(omega) < omegatol) return -11;

  /* success! */
  return ii+1;
}



