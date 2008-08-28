/* BSSN.c */
/* Bernd Bruegmann 6/02 */
/* Wolfgang Tichy  4/2004 */

#include "sgrid.h"
#include "BSSN.h"

tVarList *BSSNvars;



/* evolve in the interior and 
   for those boundary points set by special evolution
*/ 
void BSSN_evolve(tVarList *unew, tVarList *upre, double dt, tVarList *ucur)
{
  BSSN_rhs(unew, upre, dt, ucur);

  /* wether addDissipation is called after each ICN (or RK) step: */
  if(Getv("evolve_Dissipation", "yes")) 
     addDissipation(unew, upre, dt, ucur);

  if(Getv("BSSN_filter_time", "afterRHS"))
    BSSN_ChooseAndApplyFilters(unew);

  set_boundary(unew, upre, dt, ucur);

  if(Getv("BSSN_filter_time", "afterBC"))
    BSSN_ChooseAndApplyFilters(unew);

  if(Getv("BSSN_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(unew);
}





/* derive ADM variables from BSSN variables
   g_ij = exp(4 phi) gtilde_ij
   A_ij = exp(4 phi) Atilde_ij
   K_ij = A_ij + 1/3 K g_ij
*/
int BSSNtoADM(tGrid *grid) 
{
  int bi;

  for(bi = 0; bi < grid->nboxes; bi++)
  {
    tBox *box = grid->box[bi];
    int n, i;
    int ipsi    = Ind("psi");
    int ig      = Ind("gxx");
    int iK      = Ind("Kxx");
    int iphi    = Ind("BSSN_phi");
    int itrK    = Ind("BSSN_K");
    int igtilde = Ind("BSSN_gxx");
    int iAtilde = Ind("BSSN_Axx");
    double *psi = box->v[ipsi];
    double *phi = box->v[iphi];
    double *trK = box->v[itrK];
    double *g, *K, *gtilde, *Atilde;
    double fg, fK;

    /* for now always use conformal factor 
       int usepsi = 1; */
    if (0) printf("  deriving ADM variables from BSSN variables\n");

    forallpoints(box, i) {
      fg = exp(4 * phi[i]);
      fK = fg * pow(psi[i], 4.0);   // use conformal factor
      for (n = 0; n < 6; n++) {
        g = box->v[ig+n];
        K = box->v[iK+n];
        gtilde = box->v[igtilde+n];
        Atilde = box->v[iAtilde+n];
        
        g[i] = fg * gtilde[i];     
        K[i] = fK * (Atilde[i] + trK[i] * gtilde[i] / 3.0);
      }
    }
  }
  return 0;
}




/* derive BSSN variables from ADM variables
   p  = detgbar^(-1/3) 
   p0 = psi^(-4)

   gtilde_ij = p gbar_ij
   Ktilde_ij = p p0 K_ij

   phi = - log(p) / 4
   K   = gtildeinv^ij Ktilde_ij
   Atilde_ij = Ktilde_ij - gtilde_ij K / 3

   G^i = - del_j gtildeinv^ji
*/
void ADMtoBSSN(tGrid *grid) 
{
  printf("  deriving BSSN variables from ADM variables\n");

  /* compute */
  BSSN_init(grid, Ind("gxx"), Ind("Kxx"), Ind("psi"), 
	    Ind("BSSN_gxx"), Ind("BSSN_Axx"), Ind("BSSN_Gx"), 
	    Ind("BSSN_K"), Ind("BSSN_phi"),
	    Ind("alpha"), Ind("BSSN_alphaDensity"));

  /* synchronize and apply boundary conditions since Gamma was obtained 
     through finite differencing */
  //bampi_synchronize(level, Ind("BSSN_Gx"));
  // admc_setbound(GH);
}






/* initialize BSSN after initial data has been computed in POST_INITIALDATA
*/
int BSSN_startup(tGrid *grid)
{
  double vgauge;

  printf("Initializing BSSN: ");

  /* set gauge speed for lapse and related quantities */
  if (Getv("BSSN_lapse", "1+log"))
    vgauge = sqrt(Getd("BSSN_lapseharmonicf"));
  else if (Getv("BSSN_lapse", "1+log2"))
    vgauge = sqrt(4.0/3.0);
  else
    vgauge = 1;

  /* set boundary information for BSSN evolution: 
     farlimit, falloff, propagation speed 
  */
  VarNameSetBoundaryInfo("BSSN_gxx", 1, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_gxy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_gxz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_gyy", 1, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_gyz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_gzz", 1, 1, 1.0);

  VarNameSetBoundaryInfo("BSSN_Axx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Axy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Axz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Ayy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Ayz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Azz", 0, 1, 1.0);

  VarNameSetBoundaryInfo("BSSN_Gx",  0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Gy",  0, 1, 1.0);
  VarNameSetBoundaryInfo("BSSN_Gz",  0, 1, 1.0);

  VarNameSetBoundaryInfo("BSSN_K",   0, 1, vgauge);
  VarNameSetBoundaryInfo("BSSN_phi", 0, 1, vgauge);

  VarNameSetBoundaryInfo("alpha",    1, 1, vgauge);
  VarNameSetBoundaryInfo("betax",    0, 1, 1.0); 
  VarNameSetBoundaryInfo("betay",    0, 1, 1.0);
  VarNameSetBoundaryInfo("betaz",    0, 1, 1.0);
  VarNameSetBoundaryInfo("betadotx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("betadoty", 0, 1, 1.0);
  VarNameSetBoundaryInfo("betadotz", 0, 1, 1.0);

  VarNameSetBoundaryInfo("BSSN_alphaDensity", 1, 1, vgauge);

  /* create a variable list for BSSN evolutions 
     note that we include lapse and shift directly
  */
  BSSNvars = vlalloc(grid);
  vlpush(BSSNvars, Ind("BSSN_gxx"));
  vlpush(BSSNvars, Ind("BSSN_phi"));
  vlpush(BSSNvars, Ind("BSSN_Axx"));
  vlpush(BSSNvars, Ind("BSSN_K"));
  vlpush(BSSNvars, Ind("BSSN_Gx"));
  vlpush(BSSNvars, Ind("alpha"));
  vlpush(BSSNvars, Ind("betax"));
  vlpush(BSSNvars, Ind("betadotx"));
  vlpush(BSSNvars, Ind("BSSN_alphaDensity"));
  if (0) prvarlist(BSSNvars);
  enablevarlist(BSSNvars);

  /* register evolved variables */
  evolve_vlregister(BSSNvars);
  
  /* register evolution routine */
  evolve_rhsregister(BSSN_evolve);

//  /* filter all newly computed vars */
//  if(Getv("BSSN_filter_unew", "simple"))
//    evolve_algebraicConditionsregister(BSSN_filter_unew);
//  else if(Getv("BSSN_filter_unew", "naive_Ylm"))
//    evolve_algebraicConditionsregister(BSSN_naive_Ylm_filter);

  /* set K identically to zero only if we are doing maximal slicing */
  if (Getv("BSSN_forceKzero", "yes")) {
    if (!GetvLax("Gauge", "maximal"))
      errorexit("\"BSSN_forceKzero = yes\" only makes sense in conjunction "
                "with maximal slicing");
  }

  /* translate initial data in ADM variables to BSSN variables */
  ADMtoBSSN(grid);
  //set_boundary_symmetry(level, BSSNvars);


  if(Getv("BSSN_filter_type", "coordinateDependentFilter"))
    BSSN_ChooseAndApplyFilters(BSSNvars);

  if(Getv("BSSN_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(BSSNvars);
  
  /* enable all derivative vars */
  enablevar(grid, Ind("BSSN_dphix"));
  enablevar(grid, Ind("BSSN_ddphixx"));
  enablevar(grid, Ind("BSSN_dKx"));
  enablevar(grid, Ind("BSSN_dGxx"));
  enablevar(grid, Ind("BSSN_dalpx"));
  enablevar(grid, Ind("BSSN_ddalpxx"));
  enablevar(grid, Ind("BSSN_dbetaxx"));
  enablevar(grid, Ind("BSSN_ddbetaxxx"));

  return 0;
}


/* Filter BSSN vars */
int filter_VarList(tVarList *vl)
{
  tGrid *grid = vl->grid;
  int bi, j;
  
  forallboxes(grid,bi)
  {
    tBox *box = grid->box[bi];
//    double *BP = calloc(boxBoundaryPointList->npoints[bi], sizeof( double ) );
    
    /* for all variables */
    for (j = 0; j < vl->n; j++)
    {
      int vi = vl->index[j];
      double *u = box->v[vi];
      int pi,ijk;

      //printf("filter_VarList: VarName[vi]=%s\n", VarName(vi));

      /* save boundary values */
      //if(boxBoundaryPointList->npoints[bi]>0)
      //  forPointList_inbox(boxBoundaryPointList, box, pi , ijk)
      //    BP[pi] = u[ijk];
      
      /* use filters */
      spec_filter1(box, 1, u);
      spec_filter1(box, 2, u);
      spec_filter1(box, 3, u);

      /* restore boundary values */
      //if(boxBoundaryPointList->npoints[bi]>0)
      //  forPointList_inbox(boxBoundaryPointList, box, pi , ijk)
      //    u[ijk]=BP[pi];

      /* true simpleExcision hack */
/*
      {
        int n1=box->n1;
        int n2=box->n2;
        int n3=box->n3;
        int i,j,k;

        forplane1(i,j,k, n1,n2,n3, 0)
        {
          ijk=Index(i,j,k);
          u[ijk] = u[ijk+1];
        }
      }
*/
    }

    /* Hack: set deriv of ui on ExcisionBoundary */
/*
    {
      int ui = Ind("alpha");
      set_boundary_normalderiv_leftBound(selectedBoundaryPointList, 
                                         1, ui, 0.0);
    }
*/
//    free(BP);
  }
  return 0;
}


/* Filter BSSN vars */
int BSSN_filter(tGrid *grid)
{
  //printf("BSSN_filter\n");
  BSSN_ChooseAndApplyFilters(BSSNvars);

  return 0;
}


/* filter all newly computed vars */
void BSSN_filter_unew(tVarList *unew, tVarList *upre)
{
  tGrid *grid = unew->grid;
  int b;
  double n2frac = Getd("BSSN_filter_n2frac");
  double n3frac = Getd("BSSN_filter_n3frac");
  int shift2 = Geti("BSSN_filter_shift2");
  int shift3 = Geti("BSSN_filter_shift3");
  int ellipse= Getv("BSSN_filter_YZregion", "ellipse");

  coordinateDependentFilter(unew);

  /* filter high freq. angular modes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, jf,kf, vi;
    double jmax, kmax;

    jf=0.5*n2frac*n2; jf*=2; jf+=shift2-2;
    jmax=jf;
        
    kf=0.5*n3frac*n3; kf*=2; kf+=shift3-2;
    kmax=kf;

    /* filter all vars */
    for(vi=0; vi<unew->n; vi++)
    {
      double *var = vlldataptr(unew, box, vi);
      double *temp1 = box->v[Ind("temp1")];

      spec_Coeffs(box, var, temp1);
      forallijk(i,j,k)
      {
        int tj = 2*((j+1)/2);
        int tk = 2*((k+1)/2);

        if(ellipse) { if((tj/jmax)*(tj/jmax) + (tk/kmax)*(tk/kmax) > 1.0)
                        temp1[Index(i,j,k)]=0.0; }
        else        if(k>kf || j>jf) temp1[Index(i,j,k)]=0.0;
          
      }
      spec_Eval(box, var, temp1);
    }
  } /* end forallboxes */
}

/* filter unew with Ylm */
void BSSN_naive_Ylm_filter(tVarList *unew, tVarList *upre)
{
  Naive_YlmFilter(unew);
}

/* filter unew with 2/3 rule in X-direction */
void filter_with2o3rule_inX(tVarList *unew, tVarList *upre)
{
  tGrid *grid = unew->grid;
  int b;

  /* loop over boxes */
  forallboxes(grid,b)
  {
    tBox *box = grid->box[b];
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    int i,j,k, vi;

    /* filter all vars */
    for(vi=0; vi<unew->n; vi++)
    {
      double *var = vlldataptr(unew, box, vi);
      double *vc  = box->v[Ind("temp1")];

      /* compute coeffs of var in X-dir */
      spec_analysis1(box, 1, box->Mcoeffs1, var, vc);

      /* remove all coeffs with i>=2*(n1/3) */
      for(k=0; k<n3; k++)
        for(j=0; j<n2; j++)
          for(i=(n1/3)*2; i<n1; i++)
              vc[Index(i,j,k)]=0.0;

      /* use modified coeffs to change var */
      spec_synthesis1(box, 1, box->Meval1, var, vc);
    }
  }
}    

void BSSN_ChooseAndApplyFilters(tVarList *vl)
{
  if(!Getv("BSSN_filter_vars", "no"))
  {
    if(Getv("BSSN_filter_type", "old_spec_filters"))
      filter_VarList(vl);
    if(Getv("BSSN_filter_type", "coordinateDependentFilter"))
      coordinateDependentFilter(vl);
    if(Getv("BSSN_filter_type", "simple"))
      BSSN_filter_unew(BSSNvars, NULL);
    if(Getv("BSSN_filter_type", "naive_Ylm"))
      Naive_YlmFilter_lmshift(vl, 0);
    if(Getv("BSSN_filter_type", "Ylm_lmshift"))
    {
      tGrid *grid = vl->grid;
      tVarList *vl_flt;

      vl_flt = vlalloc(grid);
      vlpush(vl_flt, vl->index[6]);   // BSSN_phi
      vlpush(vl_flt, vl->index[13]);  // BSSN_K
      vlpush(vl_flt, vl->index[17]);  // alpha
      vlpush(vl_flt, vl->index[24]);  // BSSN_alphaDensity
      Naive_YlmFilter_lmshift(vl, -2);
      vlfree(vl_flt);

      vl_flt = vlalloc(grid);
      vlpush(vl_flt, vl->index[14]);  // BSSN_Gx
      vlpush(vl_flt, vl->index[18]);  // betax
      vlpush(vl_flt, vl->index[21]);  // betadotx
      Naive_YlmFilter_lmshift(vl, -1);
      vlfree(vl_flt);

      vl_flt = vlalloc(grid);
      vlpush(vl_flt, vl->index[0]);  // BSSN_gxx
      vlpush(vl_flt, vl->index[7]);  // BSSN_Axx
      Naive_YlmFilter_lmshift(vl, 0);
      vlfree(vl_flt);
    }
    if(Getv("BSSN_filter_type", "X2/3"))
      filter_with2o3rule_inX(vl, NULL);
  }
}
