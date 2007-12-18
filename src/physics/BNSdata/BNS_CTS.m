(* BNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP ,FSigma,
              dPsi[a],   dB[a,b],   dalphaP[a],    dSigma[a],
             ddPsi[a,b],ddB[a,b,c],ddalphaP[a,b], ddSigma[a,b],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
	     g[a,b], alpha, beta[a], K[a,b], q, vRS[a], x, y}

constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  Cif == nonlin,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_Psi, \ 
			Ind(\"BNSdata_Psix\"), Ind(\"BNSdata_Psixx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_B1, \
			Ind(\"BNSdata_Bxx\"), Ind(\"BNSdata_Bxxx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_alphaP, \
			Ind(\"BNSdata_alphaPx\"), Ind(\"BNSdata_alphaPxx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_Sigma, \
			Ind(\"BNSdata_Sigmax\"), Ind(\"BNSdata_Sigmaxx\"));",
  Cif == else,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lPsi, \ 
					index_dlPsi1, index_ddlPsi11);",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_lB1, \
					index_dlB11, index_ddlB111);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lalphaP, \
					index_dlalphaP1, index_ddlalphaP11);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lSigma, \
					index_dlSigma1, index_ddlSigma11);",
  Cif == end,

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  Cinstruction == "
        double xmax1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,0.0,0.0);
        double xmin1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,1.0,0.0);
        double xmax2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,1.0,0.0);
        double xmin2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,0.0,0.0);
        double R1  = 0.5*(xmax1-xmin1);
        double R2  = 0.5*(xmax2-xmin2);
	double rh1 = 0.0;
        double rh2 = 0.0;

        if(bi==0 || bi==5)  rh1 = -3.0/(R1*R1*R1);
        if(bi==3 || bi==4)  rh2 = -6.0/(R2*R2*R2); ",

n ==2,
kappa==1,
Omega ==0,
q == 1,

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega x,
  OmegaCrossR3 == 0,

  (* shift *)
  beta[a] == B[a] + OmegaCrossR[a],

  (* get 1st derivs of B *)
  (* Note: if B^i = beta^i - (Omega \times r)^i 
	=> vecLapB = vecLapbeta , LB = Lbeta, 
	since the L of any Killingvec is zero *)
  gdB == delta[a,b] dB[a,b],
  LB[a,b] == dB[a,b] + dB[b,a] -(2/3) delta[a,b] gdB,
  LBdo[a,b] == delta[a,c] delta[b,d] LB[c,d], 
  LBLB == LB[a,b] LBdo[a,b],

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,
  Psi5   == Psi*Psi4,

  (* irrotational part of 3-vel in rotating frame*)
  vRI[a] == dSigma[a],

  (* vR[a] is 3-vel. in rotating frame *)
  vR[a] == vRS[a] + vRI[a],

  (* vI[a] is vel in inertial frame *)
  vI[a] == vR[a] + OmegaCrossR[a],

  (* compute square of u^0 *)
  uzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vI[b]) (beta[c] + vI[c]),

  (* rest mass density, pressure, and total energy density *)
  rho0 == Power[q/kappa, n],
  P    == q rho0,
  rhoE == rho0 (1 + n q),

  (* fluid vars in 3+1 *)
  rho  == alpha2 (rhoE + P) uzerosqr - P,
  j[a] == alpha (rhoE + P) uzerosqr (vI[a]+beta[a]),
  S    == 3P - rhoE + rho,

  (* dLnalphaPsim6[i] = \partial_i ln(alpha Psi^{-6}) 
			= \partial_i ln(alphaP Psi^{-7})
			= Psi^7 alphaP^{-1} \partial_i(alphaP Psi^{-7}) *)
  dLnalphaPsim6[a] == dalphaP[a]/alphaP - 7 dPsi[a]/Psi,


  (* decide if use non-linear ot linear equations *)
  Cif == nonlin, (* non-linear case *)

    vecLapB[a] == delta[b,c] (ddB[a,b,c] + (1/3) ddB[b,c,a]),

    (* equations for Psi, B[a], alphaP, Sigma *)
    FPsi    == delta[b,c] ddPsi[b,c] + Psi5 LBLB/(32 alpha2) +
               2Pi Psi5 rho, 
    FB[a]   == vecLapB[a] - LB[a,b] dLnalphaPsim6[b] -
               16Pi alpha Psi4 j[a],
    FalphaP == delta[b,c] ddalphaP[b,c] - alphaP (
               (7/8) Psi4 LBLB/(4 alpha2) + 2Pi Psi4 (rho+2S) ),

    Cif == (bi==0 || bi==3), (* ell. eqn. inside stars *)
      FSigma  == delta[b,c] ddSigma[b,c] - 0,
              (*   -(vRS[a] + dSigma[a])(...), *)
    Cif == else,
      FSigma  == Sigma,  (* set Sigma=0 outside stars *)
    Cif == end,

  Cif == else, (* linear case *)

    alphaP2 == alphaP*alphaP,
    alphaP3 == alphaP2*alphaP,
    Psi3    == Psi*Psi2,
    Psi6    == Psi4*Psi2,
    Psi7    == Psi4*Psi3,

    gdlB == delta[a,b] dlB[a,b],
    LlB[a,b] == dlB[a,b] + dlB[b,a] -(2/3) delta[a,b] gdlB,
    LlBdo[a,b] == delta[a,c] delta[b,d] LlB[c,d], 
    LlBLlB == LlB[a,b] LlBdo[a,b],
    vecLaplB[a] == delta[b,c] (ddlB[a,b,c] + (1/3) ddlB[b,c,a]),

    (* linearized alpha == alphaP/Psi and vI[a] *)
    lvI[a] == dlSigma[a], 
    lalpha == lalphaP/Psi - alphaP lPsi/Psi2,
    (* linearized
      uzerosqr = alpha2 - Psi4 delta[b,c] (beta[b] + vI[b]) (beta[c] + vI[c])*)
    luzerosqr == 2alpha lalpha - 
                 4 Psi3 lPsi delta[b,c] (beta[b] + vI[b]) (beta[c] + vI[c]) -
                 Psi4 delta[b,c] 2 (lB[b] + lvI[b]) (beta[c] + vI[c]),

    (* rho  == alpha2 (rhoE + P) uzerosqr - P,
       j[a] == alpha (rhoE + P) uzerosqr (vI[a]+beta[a]),
       S    == 3P - rhoE + rho, *)
    (* linearized fluid vars in 3+1 *)
    lrho  == 2 alpha lalpha (rhoE + P) uzerosqr + 
             alpha2 (rhoE + P) luzerosqr,
    lj[a] == lalpha (rhoE + P) uzerosqr (vI[a]+beta[a]) + 
             alpha (rhoE + P) luzerosqr (vI[a]+beta[a]) +
             alpha (rhoE + P) uzerosqr (lvI[a]+lB[a]),
    lS    == lrho,

    ldLnalphaPsim6[a] == dlalphaP[a]/alphaP - dalphaP[a] lalphaP/alphaP2 -
                         7 dlPsi[a]/Psi + 7 dPsi[a] lPsi/Psi2,

    (* linearized equations for Psi, B[a], alphaP, Sigma *)
    FlPsi    == delta[b,c] ddlPsi[b,c] + 7 Psi6 lPsi LBLB/(32 alphaP2) -
                (Psi7 LBLB/(16 alphaP3)) lalphaP + 
                (Psi5/(32 alpha2)) 2 LBdo[a,b] LlB[a,b] +
                2Pi 5 Psi4 lPsi rho + 2Pi Psi5 lrho, 
    FlB[a]   == vecLaplB[a] - LlB[a,b] dLnalphaPsim6[b] - 
                LB[a,b] ldLnalphaPsim6[b] - 16Pi lalphaP Psi3 j[a] -
                16Pi alphaP 3 Psi2 lPsi j[a] - 16Pi alpha Psi4 lj[a],
    FlalphaP == delta[b,c] ddlalphaP[b,c] - lalphaP (
                (-7/32) Psi6 LBLB/(alphaP2) + 2Pi Psi4 (rho+2S) ) - alphaP (
                (21/16)Psi5 lPsi LBLB/(alphaP2) +
                (7/16) (Psi4/alpha2)LBdo[a,b] LlB[a,b] + 2Pi (
                 4 Psi3 lPsi (rho+2S) + Psi4 (lrho+2lS) ) ), 

    Cif == (bi==0 || bi==3), (* ell. eqn. inside stars *)
      FlSigma  == delta[b,c] ddlSigma[b,c] - 0,
    Cif == else,
      FlSigma  == lSigma,  (* set Sigma=0 outside stars *)
    Cif == end,

  Cif == end,


  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] :=  K[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Kup[a_,b_]  := Kup[b,a]  /; !OrderedQ[{a,b}]

ddPsi[a_,b_]     := ddPsi[b,a]    /; !OrderedQ[{a,b}]
ddB[a_,b_,c_]    := ddB[a,c,b]    /; !OrderedQ[{b,c}]
LB[a_,b_]        := LB[b,a]       /; !OrderedQ[{a,b}]
ddalphaP[a_,b_]  := ddalphaP[b,a] /; !OrderedQ[{a,b}]
ddSigma[a_,b_]   := ddSigma[b,a]  /; !OrderedQ[{a,b}]

ddlPsi[a_,b_]     := ddlPsi[b,a]    /; !OrderedQ[{a,b}]
ddlB[a_,b_,c_]    := ddlB[a,c,b]    /; !OrderedQ[{b,c}]
LlB[a_,b_]        := LlB[b,a]       /; !OrderedQ[{a,b}]
ddlalphaP[a_,b_]  := ddlalphaP[b,a] /; !OrderedQ[{a,b}]
ddlSigma[a_,b_]   := ddlSigma[b,a]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_CTS.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_CTS(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];
  pr["tGrid *grid = vlu->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{FPsi, FB[a], FalphaP ,FSigma}, "vlFu"];
  prdecvl[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];

  prdecvl[{FlPsi, FlB[a], FlalphaP ,FlSigma}, "vlJdu"];
  prdecvl[{ lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvl[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvlindices[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];
  prdecvlindices[{lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvlindices[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvarname[{dPsi[a]},       "BNSdata_Psix"];
  prdecvarname[{ddPsi[a,b]},    "BNSdata_Psixx"];
  prdecvarname[{dB[a,b]}, 	"BNSdata_Bxx"];
  prdecvarname[{ddB[a,b,c]},    "BNSdata_Bxxx"];
  prdecvarname[{dalphaP[a]},    "BNSdata_alphaPx"];
  prdecvarname[{ddalphaP[a,b]}, "BNSdata_alphaPxx"];
  prdecvarname[{dSigma[a]},     "BNSdata_Sigmax"];
  prdecvarname[{ddSigma[a,b]},  "BNSdata_Sigmaxx"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  prdecvarname[{g[a,b]}, "gxx"];
  prdecvarname[{alpha},  "alpha"];
  prdecvarname[{beta[a]},"betax"];
  prdecvarname[{K[a,b]}, "Kxx"];
  prdecvarname[{q},      "BNSdata_q"];
  prdecvarname[{vRS[a]}, "BNSdata_vRSx"];

  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["} /* end of boxes */\n"];
  pr["\n\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"
