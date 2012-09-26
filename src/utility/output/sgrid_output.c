/* sgrid_output.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 12/99 */

#include "sgrid.h"
#include "output.h"



int sgrid_output() 
{
  int nboxes = Geti("nboxes");
  int b;

  printf("Adding output\n");

  /* functions */
  AddFun(OUTPUT, write_grid, "write data to files");
  /* variables */

  /* parameters */
  AddPar("0doutiter", "-1", "when to output based on iterations");
  AddPar("1doutiter", "-1", "when to output based on iterations");
  AddPar("2doutiter", "-1", "when to output based on iterations");
  AddPar("3doutiter", "-1", "when to output based on iterations");

  AddPar("0douttime", "-1", "when to output based on time");
  AddPar("1douttime", "-1", "when to output based on time");
  AddPar("2douttime", "-1", "when to output based on time");
  AddPar("3douttime", "-1", "when to output based on time");

  AddPar("1doutinterpolate", "no", 
	 "whether to interpolate onto standard coord");
  AddPar("2doutinterpolate", "no", 
	 "whether to interpolate onto standard coord");
  AddPar("3doutinterpolate", "no", 
	 "whether to interpolate onto standard coord");
  AddPar("0doutput_normtype", "integral", 
	 "how we compute norms such as rms [integral,L2norm]");
  AddPar("0doutput_VolumeIntegralJacobian", "fromCoordinates", 
	 "which Jacobian we use for vol. integrals [fromCoordinates,one]");

  AddPar("0doutput", "", "variables to output along axes");
  AddPar("1doutput", "", "variables to output along axes");
  AddPar("2doutput", "", "variables to output on coordinate planes");
  AddPar("3doutput", "", "variables to output in 3d");

  AddPar("outputall",   "no", "whether to output all components");
  AddPar("0doutputall", "no", "whether to output all components");
  AddPar("1doutputall", "no", "whether to output all components");
  AddPar("2doutputall", "no", "whether to output all components");
  AddPar("3doutputall", "no", "whether to output all components");

  AddPar("0doutputmax",  "", "variables with output of maximum");
  AddPar("0doutputmin",  "", "variables with output of minimum");
  AddPar("0doutputrms", "", "variables with output of l2-norm");
  AddPar("0doutputmeanAbs", "", "variables with output of l_infinity-norm");
  AddPar("0doutputmean",    "", "variables with output of mean");
  AddPar("0doutputVolInt",  "", "variables with output of Volume Integral");

  AddPar("1doutputX", "", "variables to output along X axis");
  AddPar("1doutputY", "", "variables to output along Y axis");
  AddPar("1doutputZ", "", "variables to output along Z axis");
  AddPar("1doutputD", "", "variables to output along Diagonal");

  AddPar("2doutputXY", "", "variables to output in X-Y-plane");
  AddPar("2doutputXZ", "", "variables to output in X-Z-plane");
  AddPar("2doutputYZ", "", "variables to output in Y-Z-plane");
  AddPar("2doutputXD", "", "variables to output in X-Y=Z-plane");
  AddPar("2doutputYD", "", "variables to output in Y-X=Z-plane");
  AddPar("2doutputZD", "", "variables to output in Z-X=Y-plane");

  AddPar("2dformat", "opendx binary float",
	 "format for 2d output (opendx,vtk,text,binary,float,double)"); 

  AddPar("3dformat", "opendx binary float", "format for 3d output "
	 "(opendx,vtk,hdf5,amrunion,text,binary,float,double,"
	 "add_rotant_points)"); 
  AddPar("hdf5_timesteps_in_one_file", "yes", 
         "whether all timesteps are put into one hdf5 file");

  for(b=0; b<nboxes; b++)
  {
    char str[1000];

    /* origin for output */
    snprintf(str, 999, "outputX0_box%d", b);
    AddPar(str,   "0", "origin for output in X");
    snprintf(str, 999, "outputY0_box%d", b);
    AddPar(str,   "0", "origin for output in Y");
    snprintf(str, 999, "outputZ0_box%d", b);
    AddPar(str,   "0", "origin for output in Z");
    /* get X, Y, Z to be used for output */
    snprintf(str, 999, "outputReplaceXby_box%d", b);
    AddPar(str,   "X", "origin for output in X");
    snprintf(str, 999, "outputReplaceYby_box%d", b);
    AddPar(str,   "Y", "origin for output in Y");
    snprintf(str, 999, "outputReplaceZby_box%d", b);
    AddPar(str,   "Z", "origin for output in Z");
  }
  return 0;
}
