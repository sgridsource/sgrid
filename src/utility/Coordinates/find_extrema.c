/* find_extrema.c */
/* Wolfgang Tichy 2/2018 */

#include "sgrid.h"
#include "Coordinates.h"



/* par type used for box_extremum_of_F */
typedef struct T_box_c1c2c3_struct {
  tBox *box;   /* box */
  double *c1;  /* coeffs of fvec1 */
  double *c2;  /* coeffs of fvec2 */
  double *c3;  /* coeffs of fvec3 */
  double C[4]; /* optional coordinate values, e.g. C[1]=X */
  int dir;     /* direction 1,2, or 3 */
  int status;  /* status code */
} t_box_c1c2c3_struct;



/* get fvec from coeffs, e.g. coeffs could come from dF/dX of F */
void f3vec_from_c1c2c3_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  t_box_c1c2c3_struct *par = (t_box_c1c2c3_struct *) p;
  tBox *box  = par->box;
  double *cX = par->c1;
  double *cY = par->c2;
  double *cZ = par->c3;
  double X = vec[1];
  double Y = vec[2];
  double Z = vec[3];
  double X0 = box->bbox[0];
  double X1 = box->bbox[1];
  double Y0 = box->bbox[2];
  double Y1 = box->bbox[3];
  double Z0 = box->bbox[4];
  double Z1 = box->bbox[5];
  double f1, f2, f3;
  int ox,oy,oz, ux,uy,uz; /* overrun & underrun flags */

  /* set overrun & underrun flags */
  ox=oy=oz=ux=uy=uz=0;
  if(X<X0) { X=X0; ux=1; }
  if(Y<Y0) { Y=Y0; uy=2; }
  if(Z<Z0) { Z=Z0; uz=4; }
  if(X>X1) { X=X1; ox=8; }
  if(Y>Y1) { Y=Y1; oy=16; }
  if(Z>Z1) { Z=Z1; oz=32; }

  //printf("P: %d %g %g %g\n", box->b, X,Y,Z);
  //quick_Array_output(box, cX, "cX", X,X*1000, 0);

  /* interpolate */
  fvec[1] = f1 = spec_interpolate(box, cX, X,Y,Z);
  fvec[2] = f2 = spec_interpolate(box, cY, X,Y,Z);
  fvec[3] = f3 = spec_interpolate(box, cZ, X,Y,Z);

  /* put zero on boundary and continue linearly if over- or underrun */
  if(ux) fvec[1] =  f1 * (vec[1]-X0);
  if(uy) fvec[2] =  f2 * (vec[2]-Y0);
  if(uz) fvec[3] =  f3 * (vec[3]-Z0);
  if(ox) fvec[1] = -f1 * (vec[1]-X1);
  if(oy) fvec[2] = -f2 * (vec[2]-Y1);
  if(oz) fvec[3] = -f3 * (vec[3]-Z1);
  par->status = ox | oy | oz | ux | uy | uz;
//printf("ox oy oz ux uy uz = %d %d %d %d %d %d", ox, oy, oz, ux, uy, uz);
//printf("vec[1],vec[2],vec[3]=%.15g %.15g %.15g\n", vec[1],vec[2],vec[3]);
//printf("X,Y,Z=%.15g %.15g %.15g\n", X,Y,Z);
//printf("f1,f2,f3 = %.15g %.15g %.15g\n", f1,f2,f3);
//printf("fvec[1],fvec[2],fvec[3] = %.15g %.15g %.15g\n", fvec[1],fvec[2],fvec[3]);
}

/* Find min/max of function F on the grid in a box */
/* Note:  
   dF/dx = dF/dX dX/dx = 0  | assume Jacobian matrix dX/dx is invertible
     <==> dF/dX = 0  */
int box_extremum_of_F(tBox *box, int Fi,
                      double *X, double *Y, double *Z, double *Fextr)
{
  double X0 = box->bbox[0];
  double X1 = box->bbox[1];
  double Y0 = box->bbox[2];
  double Y1 = box->bbox[3];
  double Z0 = box->bbox[4];
  double Z1 = box->bbox[5];
  int ox,oy,oz, ux,uy,uz; /* overrun & underrun flags */
  int N = box->nnodes;
  double *cx = dmalloc(N);
  double *cy = dmalloc(N);
  double *cz = dmalloc(N);
  double *F  = box->v[Fi];
  double Xvec[4];
  t_box_c1c2c3_struct par[1];
  int check, stat, status;

  /* get all 3 derivs */
  spec_Deriv1(box, 1, F, cx);
  spec_Deriv1(box, 2, F, cy);
  spec_Deriv1(box, 3, F, cz);
  /* transform derivs to coeffs, can be done in place */
  spec_Coeffs(box, cx, cx);
  spec_Coeffs(box, cy, cy);
  spec_Coeffs(box, cz, cz);

  //printf("A: %d %g %g %g %s:%p\n", box->b, *X,*Y,*Z, VarName(Fi), F);
  //quick_Array_output(box, cx, "cx", *X,*X*1000, 0);

  /* set Xvec */
  Xvec[1] = *X;
  Xvec[2] = *Y;
  Xvec[3] = *Z;

  /* use newton_linesrch_itsP to find max, use par for parameters */
  par->box = box;
  par->c1  = cx;
  par->c2  = cy;
  par->c3  = cz;
  stat = newton_linesrch_itsP(Xvec, 3, &check, f3vec_from_c1c2c3_VectorFuncP,
                              (void *) par, 1000, dequaleps);
  if(check || stat<0) printf("  --> check=%d stat=%d\n", check, stat);
  if(check || stat<0) printf("  --> Xvec=(%.15g,%.15g,%.15g)\n",
                               Xvec[1],Xvec[2],Xvec[3]);
  *X = Xvec[1];
  *Y = Xvec[2];
  *Z = Xvec[3];

  /* check if we are outside box, set overrun & underrun flags */
  ox=oy=oz=ux=uy=uz=0;
  if(*X<X0) { *X=X0; ux=1; }
  if(*Y<Y0) { *Y=Y0; uy=2; }
  if(*Z<Z0) { *Z=Z0; uz=4; }
  if(*X>X1) { *X=X1; ox=8; }
  if(*Y>Y1) { *Y=Y1; oy=16; }
  if(*Z>Z1) { *Z=Z1; oz=32; }
  status = -(ox | oy | oz | ux | uy | uz);
  if(status==0) status = stat;

  /* set coeffs of F in cx, and set Fextr at rhis point */
  spec_Coeffs(box, F, cx);
  *Fextr = spec_interpolate(box, cx, *X,*Y,*Z);
  //printf("B: %d %g %g %g\n", box->b, *X,*Y,*Z);
  free(cz);
  free(cy);
  free(cx);
  return status;
}


/* for zbrent_itsP: */
/* get value from coeffs, e.g. coeffs could come from dF/dX of F */
double f1_from_c1_dir_ZP(double vec, void *p)
{
  t_box_c1c2c3_struct *par = (t_box_c1c2c3_struct *) p;
  tBox *box  = par->box;
  double *co = par->c1;
  double C1  = par->C[1];
  double C2  = par->C[2];
  int dir    = par->dir;
  double X,Y,Z, f1;
  /* set X,Y,Z */
  switch(dir)
  {
    case 1:  X=vec; Y=C1; Z=C2;  break;
    case 2 : Y=vec; X=C1; Z=C2;  break;
    case 3:  Z=vec; X=C1; Y=C2;  break;
    default: errorexit("dir must be 1,2,3");
  }
  //printf("Pd: %d %g %g %g\n", box->b, X,Y,Z);
  //quick_Array_output(box, co, "co", vec,vec*1000, 0);

  /* interpolate */
  f1 = spec_interpolate(box, co, X,Y,Z);
  //printf("X,Y,Z=%g %g %g  f1=%g\n", X,Y,Z, f1);
  return f1;
}
/* get fvec from coeffs, e.g. coeffs could come from dF/dX of F */
void f1vec_from_c1_dir_VectorFuncP(int n, double *vec, double *fvec, void *p)
{
  fvec[1] = f1_from_c1_dir_ZP(vec[1], p);
}


/* find extremum along coord line, line is picked by direction dir.
   e.g. if dir=1, C1=Y, C2=Z, and it will look along X and return X in C */
int box_extremum_of_F_in_dir(tBox *box, int Fi, int dir, double C1, double C2,
                             double *C, double *Fextr)
{
  int N = box->nnodes;
  double *c = dmalloc(N);
  double *F = box->v[Fi];
  double Cl = box->bbox[(dir-1)*2];   /* bounds of C: C \in [Cl,Ch] */
  double Ch = box->bbox[(dir-1)*2+1];
  double X,Y,Z;
  t_box_c1c2c3_struct par[1];
  int stat, status;

  /* get deriv in c */
  spec_Deriv1(box, dir, F, c);
  /* transform deriv c to coeffs c, can be done in place */
  spec_Coeffs(box, c, c);

  //printf("Ad: %d %g %g %g %s:%p\n", box->b, X,Y,Z, VarName(Fi), F);
  //quick_Array_output(box, c, "c", *C,*C*1000, 0);

  /* set pars */
  par->box = box;
  par->c1  = c;
  par->dir = dir;
  par->C[1] = C1;
  par->C[2] = C2;

  //printf("Cl=%g Ch=%g\n", Cl,Ch);
  /* bracket root */
  stat = zbrac_P(f1_from_c1_dir_ZP, &Cl,&Ch, (void *) par);
  if(stat<0) return stat;
  //printf("Cl=%g Ch=%g\n", Cl,Ch);

  //printf("C=%g\n", *C);
  /* look for root in [Cl,Ch] */
  status = zbrent_itsP(C, f1_from_c1_dir_ZP, Cl,Ch, (void *) par,
                       1000, dequaleps);
  //printf("C=%g\n", *C);

  /* set coeffs of F in cx, and set Fextr at rhis point */
  spec_Coeffs(box, F, c);
  switch(dir)
  {
    case 1:  X=*C; Y=C1; Z=C2;  break;
    case 2 : Y=*C; X=C1; Z=C2;  break;
    case 3:  Z=*C; X=C1; Y=C2;  break;
    default: errorexit("dir must be 1,2,3");
  }
  *Fextr = spec_interpolate(box, c, X,Y,Z);
  //printf("Bd: %d %g %g %g\n", box->b, X,Y,Z);
  free(c);
  return status;
}
