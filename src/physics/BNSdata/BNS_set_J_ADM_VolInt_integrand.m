(* BNS_set_J_ADM_VolInt_integrand.m 
   Wolfgang Tichy  11/2008       *)

(* set integrand for inner volume integrals of J_ADM *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dSigma[a], q, vRS[a], x, y, Integ}


constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_S(box, Ind(\"BNSdata_Sigma\"),
					Ind(\"BNSdata_Sigmax\"));",

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* set lapse and Psi4 *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xCM),
  OmegaCrossR3 == 0,

  (* shift in rotating frame *)
  beta[a] == B[a] + OmegaCrossR[a],

  (* irrotational part of 3-vel in rotating frame*)
  vRI[a] == dSigma[a],

  (* vR[a] is 3-vel. in rotating frame *)
  vR[a] == vRS[a] + vRI[a],

  (* compute square of u^0 in rotating frame *)
  oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
  uzerosqr == 1.0/oouzerosqr,

  (* rest mass density, pressure, and total energy density *)
  rho0 == Power[q/kappa, n],
  P    == q rho0,
  rhoE == rho0 (1 + n q),

  (* if q=0 all matter vars are zero, which can be enforced by uzerosqr=0 *)
  Cif == (q==0),
    uzerosqr == 0,
  Cif == end,

  (* get fluid var j in 3+1 *)
  jup[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),

 (* Integrand for rest mass with 3-volume element factor Psi^6 *)
  Integ == Psi4 ( (x-xCM)*jup2 - y*jup1 ) * Psi4*Psi2,

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_set_J_ADM_VolInt_integrand.c"

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

  pr["void BNS_set_J_ADM_VolInt_integrand(tGrid *grid, int iInteg)\n"];
  pr["{\n"];

  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"BNSdata_x_CM\");\n"];
  pr["\n"];

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

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  prdecvarname[{Psi},     "BNSdata_Psi"];
  prdecvarname[{B[a]},    "BNSdata_Bx"];
  prdecvarname[{alphaP},  "BNSdata_alphaP"];
  prdecvarname[{Sigma},   "BNSdata_Sigma"];
  prdecvarname[{q},      "BNSdata_q"];
  prdecvarname[{vRS[a]}, "BNSdata_vRSx"];
  prdecvarname[{dSigma[a]},     "BNSdata_Sigmax"];

  prdecvar[{Integ}, "iInteg"];

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
