/* sgrid_Coordinates.c */
/* Wolfgang Tichy 2003 */

#include "sgrid.h"
#include "Coordinates.h"


int sgrid_Coordinates(void) 
{
  int nboxes = Geti("nboxes");
  int b;
  
  printf("Adding Coordinates\n");

  /* functions */
  AddFun(POST_GRID, set_box_CI_struct_from_pars,
         "read box_CI pars and set box->CI struct");
  AddFun(COORDINATES, init_CoordTransform_And_Derivs, 
         "initialize coords and coord transforms");

  /* variables */
  AddVar("x", "", "cartesian x coordinate");
  AddVar("y", "", "cartesian y coordinate");
  AddVar("z", "", "cartesian z coordinate");
  AddVar("dXd", "i", "coordinate derivative dX/dx^i");
  AddVar("dYd", "i", "coordinate derivative dY/dx^i");
  AddVar("dZd", "i", "coordinate derivative dZ/dx^i");
  AddVar("ddXdd", "(ij)", "2nd coordinate derivative d^2 X/dx^i dx^j");
  AddVar("ddYdd", "(ij)", "2nd coordinate derivative d^2 Y/dx^i dx^j");
  AddVar("ddZdd", "(ij)", "2nd coordinate derivative d^2 Z/dx^i dx^j");
  AddVar("oX", "", "X-coord from other boxes for each point in box-bfaces");
  AddVar("oY", "", "Y-coord from other boxes for each point in box-bfaces");
  AddVar("oZ", "", "Z-coord from other boxes for each point in box-bfaces");
  AddVar("Temp1", "",     "temporary variable(e.g. to store coeffs)");
  /* create vars that contain cub. sph. sigma_{0/1} and their derivs */
  AddVar("Coordinates_CubedSphere_sigma01",     "", "sigma_{0/1}");
  AddVar("Coordinates_CubedSphere_dsigma01_dA", "", "d/dA sigma_{0/1}");
  AddVar("Coordinates_CubedSphere_dsigma01_dB", "", "d/dB sigma_{0/1}");
  
  /* parameters */
  for(b=0; b<nboxes; b++)
  {
    char str[1000];

    snprintf(str, 999, "box%d_Coordinates", b);
    AddPar(str, "Cartesian", 
           "coordinates used in box [Cartesian, Polar, ...]");

    snprintf(str, 999, "box%d_CI_s", b);
    AddPar(str, "", "box->CI->s part of tCoordInfo");
    snprintf(str, 999, "box%d_CI_xc", b);
    AddPar(str, "", "box->CI->xc part of tCoordInfo");
    snprintf(str, 999, "box%d_CI_dom", b);
    AddPar(str, "", "box->CI->dom part of tCoordInfo");
    snprintf(str, 999, "box%d_CI_type", b);
    AddPar(str, "", "box->CI->type part of tCoordInfo");

    snprintf(str, 999, "box%d_CoordinateTransforms_generic", b);
    AddPar(str, "no", 
           "select dXdx or ddXdxdx, to compute them using spectral derivs "
           "[no,dXdx,ddXdxdx]");
  }

  AddPar("Coordinates_verbose", "yes", "verbose [yes,no]");
  AddPar("CoordinateTransforms_stored", "yes",
         "whether we store Coordinate Transforms in dXdx,... ddXddxx,...");
  AddPar("Coordinates_newtTOLF", "1e-10", "newton tolerence");
  AddPar("Coordinates_newtMAXITS", "100000", "max. newton iterations");
  AddPar("Coordinates_useDD", "no", "if we use DD in D_and_DD_Of_S..");
  AddPar("compactSphericalDF_r0", "-1", "radius r at xi=0");
  AddPar("tan_stretch_s", "0", "how much we stretch [0,Xmax]");
  AddPar("Coordinates_Spherical3_c",  "5",
         "constant c in Spherical3 coord trafo");
  AddPar("Coordinates_AnsorgNS_sigma_pm_vars", "no",
         "create vars that contain sigma_{+-} and their derivs [yes,no]");
  AddPar("Coordinates_AnsorgNS_b", "1", "value of x if A=1 in AnsorgNS0/3");
  AddPar("Coordinates_AnsorgNS_set_sigma_pm_pointers", "yes",
         "initialize sigma_{+-} in init_CoordTransform_And_Derivs [yes,no]");
  AddPar("Coordinates_AnsorgNS_Bshift", "no",
         "shift B by 0.5/((1+N mod 2)*N) [no,yes]");
  AddPar("Coordinates_AnsorgNS_version", "AnsorgNS", "version we use" 
         " AnsorgNS -> old version in Coordinates.c, "
         "NAnsorgNS -> new version in Coordinates_AnsorgNS.c and coordtrans_AnsorgNS?.m, "
         "DDAnsorgNS -> like AnsorgNS but add some 2nd derivs from NAnsorgNS");
  if(Getv("Coordinates_AnsorgNS_sigma_pm_vars", "yes"))
  {
    AddVar("Coordinates_AnsorgNS_sigma_pm",       "", "sigma_{+-}");
    AddVar("Coordinates_AnsorgNS_dsigma_pm_dB",   "", "d/dB sigma_{+-}");
    AddVar("Coordinates_AnsorgNS_dsigma_pm_dphi", "", "d/dphi sigma_{+-}");
    AddVar("Coordinates_AnsorgNS_ddsigma_pm_dBdB",    "", "(d/dB)^2 sigma_{+-}");
    AddVar("Coordinates_AnsorgNS_ddsigma_pm_dBdphi",  "", "d/dB d/dphi sigma_{+-}");
    AddVar("Coordinates_AnsorgNS_ddsigma_pm_dphidphi","", "(d/dphi)^2 sigma_{+-}");
    AddPar("Coordinates_AnsorgNS_dsigma_pm_dphi_ZeroOnAxis", "no", "whether "
           "dsigma_pm_dphi and ddsigma_pm_dphidphi are 0 on x-axis [no,yes]");
  }
  AddPar("Coordinates_set_bfaces", "no", "whether we set bfaces in the "
         "function init_CoordTransform_And_Derivs [no,yes]");

  return 0;
}
