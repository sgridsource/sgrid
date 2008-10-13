/* DV_CircSchwSource.c */
/* Wolfgang Tichy  8/2007 */

#include "sgrid.h"
#include "DV_CircSchwSource.h"
#include "Constants.h"
#include "Source.h"


/* initialize DV_CircSchwSource */
int DV_CircSchwSource_startup(tGrid *grid)
{
  double x,y;
  double M=1.0;
  int i,j;
  printf("Initializing DV_CircSchwSource:\n");
  
  //set_mass_radius(M,10.0*M);
  set_parameters(M, 1.0, 10.0*M, Getd("DV_CircSchwSource_DVWindow_n"), 
                 Getd("DV_CircSchwSource_DVWindow_width"));

//  for(i=0;i<40;i++){
//        x=9.9*M+i*0.005*M;
//        for(j=0;j<40;j++){
//            y=-0.1*M+j*0.005*M;
//            printf ("%f %f %3.10e \n", x, y,SourceInKerrSchild(4.0*log(2.0)*M,x,y,0.0));
//        }
//        printf("\n");
//  }
  return 0;
}
