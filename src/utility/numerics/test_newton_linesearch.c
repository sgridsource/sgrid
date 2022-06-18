#define finit isfinite

#include "lu.c"
#include "newton_linesearch.c"

void func(int n, double x[], double f[], void *par)
{
  double *r2 = par;
  f[0] = x[0]*x[0] + x[1]*x[1] - *r2;
  f[1] = x[1] - exp(x[0]);
  printf("x=(%.11g, %.11g)  f(x)=(%g, %g)\n", x[0],x[1], f[0],f[1]);
}

int main()
{
  int its, check;
  double rad2 = 4;
  double x[] = {2., 0.};

  its = newton_linesearch(2, x, &check, func,&rad2,0, 50, 1e-8);
  printf("root is at: x=(%.15g, %.15g)   "
         "its=%d check=%d\n", x[0],x[1], its,check);
}
