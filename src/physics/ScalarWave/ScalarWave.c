/* ScalarWave.c */
/* Wolfgang Tichy  8/2005 */

#include "sgrid.h"
#include "ScalarWave.h"




/* evolve in the interior and 
   for those boundary points set by special evolution
*/ 
void ScalarWave_evolve(tVarList *unew, tVarList *upre, double dt, 
                       tVarList *ucur)
{
  tGrid *grid = ucur->grid;
  int b;
  double nonlin = Getd("ScalarWave_nonlinearity");
  
  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *cpsi = vlldataptr(ucur, box, 0);
    double *ppsi = vlldataptr(upre, box, 0);
    double *npsi = vlldataptr(unew, box, 0);
    int ipsi = (ucur)->index[0];
    double *cpsidot = vlldataptr(ucur, box, 1);
    double *ppsidot = vlldataptr(upre, box, 1);
    double *npsidot = vlldataptr(unew, box, 1);
    int i;
    //double *psix = box->v[Ind("ScalarWave_dpsix")];
    //double *psiy = box->v[Ind("ScalarWave_dpsix")+1];
    //double *psiz = box->v[Ind("ScalarWave_dpsix")+2];
    double *psixx = box->v[Ind("ScalarWave_ddpsixx")];
    //double *psixy = box->v[Ind("ScalarWave_ddpsixx")+1];
    //double *psixz = box->v[Ind("ScalarWave_ddpsixx")+2];
    double *psiyy = box->v[Ind("ScalarWave_ddpsixx")+3];
    //double *psiyz = box->v[Ind("ScalarWave_ddpsixx")+4];
    double *psizz = box->v[Ind("ScalarWave_ddpsixx")+5];

    allDerivsOf_S(box, ipsi,
                  Ind("ScalarWave_dpsix"), Ind("ScalarWave_ddpsixx"));

    forallpoints(box, i)
    {
      double rpsidot, rpsi;
      
      rpsidot = psixx[i] + psiyy[i] + psizz[i]
                - nonlin * cpsi[i]*cpsi[i]*cpsi[i]/(1.0+cpsi[i]*cpsi[i]);
      rpsi = cpsidot[i];
      
      if(dt!=0.0)
      {
        npsidot[i] = ppsidot[i] + dt * rpsidot;
        npsi[i] = ppsi[i] + dt * rpsi;
      }
      else
      {
        npsidot[i] = rpsidot;
        npsi[i] = rpsi;
      }
    }
  }

  /* wether addDissipation is called after each ICN (or RK) step: */
  if(Getv("evolve_Dissipation", "yes")) 
     addDissipation(unew, upre, dt, ucur);

  set_boundary(unew, upre, dt, ucur);

  if(Getv("ScalarWave_coordinateDependentFilter", "yes"))
    coordinateDependentFilter(unew);
  
  if(Getv("ScalarWave_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(unew);
}





/* initialize ScalarWave */
int ScalarWave_startup(tGrid *grid)
{
  tVarList *ScalarWavevars;
  int b;
  double A         = Getd("ScalarWave_A");
  double sigma     = Getd("ScalarWave_sigma");
  double r0        = Getd("ScalarWave_r0");
  double sigmax    = Getd("ScalarWave_sigmax");
  double sigmay    = Getd("ScalarWave_sigmay");
  double sigmaz    = Getd("ScalarWave_sigmaz");
  double x0        = Getd("ScalarWave_x0");
  double y0        = Getd("ScalarWave_y0");
  double z0        = Getd("ScalarWave_z0");
  int spherical    = Getv("ScalarWave_waveform", "spherical");
  
  printf("Initializing ScalarWave:\n");

  /* set boundary information for ScalarWave evolution: 
     farlimit, falloff, propagation speed 
  */
  VarNameSetBoundaryInfo("ScalarWave_psi",     0, 1, 1.0);

  VarNameSetBoundaryInfo("ScalarWave_dpsix",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_dpsiy",   0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_dpsiz",   0, 1, 1.0);

  VarNameSetBoundaryInfo("ScalarWave_ddpsixx", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_ddpsixy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_ddpsixz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_ddpsiyy", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_ddpsiyz", 0, 1, 1.0);
  VarNameSetBoundaryInfo("ScalarWave_ddpsizz", 0, 1, 1.0);

  VarNameSetBoundaryInfo("ScalarWave_psidot",  0, 1, 1.0);

  /* create a variable list for ScalarWave evolutions 
     note that we include lapse and shift directly
  */
  ScalarWavevars = vlalloc(grid);
  vlpush(ScalarWavevars, Ind("ScalarWave_psi"));
  vlpush(ScalarWavevars, Ind("ScalarWave_psidot"));
  if (0) prvarlist(ScalarWavevars);
  enablevarlist(ScalarWavevars);

  /* register evolved variables */
  evolve_vlregister(ScalarWavevars);
  
  /* register evolution routine */
  evolve_rhsregister(ScalarWave_evolve);

  /* set initial data in boxes */
  forallboxes(grid,b)
  {  
    tBox *box = grid->box[b];
    int i,j,k;
    int n1=box->n1;
    int n2=box->n2;
    int n3=box->n3;
    double *pX = box->v[Ind("X")];
    double *pY = box->v[Ind("Y")];
    double *pZ = box->v[Ind("Z")];
    double *px = box->v[Ind("x")];
    double *py = box->v[Ind("y")];
    double *pz = box->v[Ind("z")];

    double *psi = box->v[Ind("ScalarWave_psi")];
    double *psidot= box->v[Ind("ScalarWave_psidot")];
//test
/*
forallpoints(box,i)
{
psi[i]=1.0/(2.0*PI);
//psi[i]=(pX[i]*pX[i])*fabs(sin(acos(pz[i]/pX[i]))); //correct weight, but not analytic
//printf("acos(pz[i]/pX[i])=%f pY[i]=%f\n", acos(pz[i]/pX[i]), pY[i]);
psi[i]=(pX[i]*pX[i]);
psi[i]=(pX[i]*pX[i])*fabs(sin(pZ[i]));
psi[i]=( sin(1*(pY[i]+PI/n2)) + cos(4*(pY[i]+PI/n2)) )*1.0/(2.0*PI);
psi[i]=(pX[i]*pX[i])*sin(pZ[i]-1)*cos(pZ[i])*( sin(1*(pY[i]+PI/n2)) + cos(4*(pY[i]+PI/n2)) )/(1.4374629934615632542716*(-2.6435590640814554674876));
psi[i]=(pX[i]*pX[i]);
psi[i]=1;
psi[i]=sin(pZ[i]-1)*cos(pZ[i])*( sin(1*(pY[i]+PI/n2)) + cos(4*(pY[i]+PI/n2)) )/(1.4374629934615632542716*(-2.6435590640814554674876));
//psi[i] *= fabs(sin(acos(pz[i]/pX[i])));
//psi[i]=(pX[i]*pX[i])*fabs(sin(acos(pz[i]/pX[i])))/2; 
}
spec_Integral1(box,1, psi, psidot);
//spec_2dIntegral(box,2, psi, psidot);
//spec_3dIntegral(box, psi, psidot);

//spec_Integral1(box, 1, psi, psidot);
//spec_Integral1(box, 2, psidot, psidot);
//spec_Integral1(box, 3, psidot, psidot);

spec_3dIntegral(box, psi, psidot);

spec_sphericalDF2dIntegral(box, psi, psidot);
spec_Integral1(box, 1, psidot, psidot);

//spec_sphericalDF3dIntegral(box, psi, psidot);

i=0; j=0; k=0;
printf("psidot[Index(%d,%d,%d)]=%.22f\n", i,j,k, psidot[Index(i,j,k)]);  

forallpoints(box,i) 
  if(psidot[i] != psidot[0]) printf("i=%d ", i);
*/
//end test
    forallpoints(box,i)
    {
      double x = pX[i];
      double y = pY[i];
      double z = pZ[i];
      double r, f, df;

      if(px!=NULL) 
      {
        x = px[i];
        y = py[i];
        z = pz[i];
      }
      
      r = sqrt(x*x + y*y + z*z);

      if(spherical)
      {
        f = A*exp( -(r-r0)*(r-r0) / (2.0*sigma*sigma) );
        df= ( -1.0/(sigma*sigma) )*(r-r0)*f;
      }
      else
      {
        f = A*exp( -0.5*( 
               (x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay)
             + (z-z0)*(z-z0)/(sigmaz*sigmaz)   ) );
        df= -f*(
                 (x-x0)*(x/r)/(sigmax*sigmax) + (y-y0)*(y/r)/(sigmay*sigmay)
                +(z-z0)*(z/r)/(sigmaz*sigmaz)  );
      }
      psi[i]    =  f/r;
      psidot[i] = df/r;
    }
    /* set psi=0 on boundaries */
    for(k=0; k<n3; k++)
      for(j=0; j<n2; j++)
      {
        psi[Index(0,j,k)]    = psi[Index(n1-1,j,k)] = 0;
        psidot[Index(0,j,k)] = psidot[Index(n1-1,j,k)] = 0;
      }
  }

  //set_boundary_symmetry(level, ScalarWavevars);

  if(Getv("ScalarWave_coordinateDependentFilter", "yes"))
    coordinateDependentFilter(ScalarWavevars);

  if(Getv("ScalarWave_reset_doubleCoveredPoints", "yes"))
    reset_doubleCoveredPoints(ScalarWavevars);
  
  /* enable all derivative vars */
  enablevar(grid, Ind("ScalarWave_dpsix"));
  enablevar(grid, Ind("ScalarWave_ddpsixx"));

  /* enable all rho vars */
  enablevar(grid, Ind("ScalarWave_rho"));
  enablevar(grid, Ind("ScalarWave_2dInt_rho"));
  /* enablevar(grid, Ind("ScalarWave_temp")); */

  return 0;
}


/* compute rho and integrals of rho */
int ScalarWave_analyze(tGrid *grid)
{
  int b;
  double nonlin;

  if( ! timeforoutput_index(grid, Ind("ScalarWave_rho")) ) return 0;

  nonlin = Getd("ScalarWave_nonlinearity");

  for (b = 0; b < grid->nboxes; b++)
  {
    tBox *box = grid->box[b];
    double *psi    = box->v[Ind("ScalarWave_psi")];
    double *psidot = box->v[Ind("ScalarWave_psidot")];
    double *psix   = box->v[Ind("ScalarWave_dpsix")];
    double *psiy   = box->v[Ind("ScalarWave_dpsix")+1];
    double *psiz   = box->v[Ind("ScalarWave_dpsix")+2];
    double *rho    = box->v[Ind("ScalarWave_rho")];
    double *I2rho  = box->v[Ind("ScalarWave_2dInt_rho")];
    /* double *I3rho  = box->v[Ind("ScalarWave_temp")]; */
    int i;

    FirstDerivsOf_S(box, Ind("ScalarWave_psi"), Ind("ScalarWave_dpsix"));

    /* compute rho */
    forallpoints(box, i)
      rho[i] = 0.5*(  psidot[i]*psidot[i]  // is this correct???
                 + psix[i]*psix[i] + psiy[i]*psiy[i] + psiz[i]*psiz[i] )
               + nonlin*( psi[i]*psi[i] - log(1+psi[i]*psi[i]) )*0.5;

    spec_sphericalDF2dIntegral(box, rho, I2rho);
    /* spec_Integral1(box, 1, I2rho, I3rho); */
  }
  return 0;
}
