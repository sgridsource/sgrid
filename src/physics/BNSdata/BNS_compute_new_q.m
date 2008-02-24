(* BNS_compute_new_q.m 
   Wolfgang Tichy  2/2008       *)

(* compute the updated q *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dSigma[a], q, vRS[a], x, y}

constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_S(box,index_BNSdata_Sigma,
                                   Ind(\"BNSdata_Sigmax\"));",

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xCM),
  OmegaCrossR3 == 0,

  (* shift in rotating frame (in the inertial frame beta^i = B^i) *)
  beta[a] == B[a] + OmegaCrossR[a],

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,

  (* irrotational part of 3-vel in rotating frame*)
  vRI[a] == dSigma[a],

  (* vR[a] is 3-vel. in rotating frame *)
  vR[a] == vRS[a] + vRI[a],

  (* vI[a] is vel in inertial frame *)
  (* vI[a] == vR[a] + OmegaCrossR[a], *)

  (* compute u^0 in rotating frame *)
  oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
  Cif == (oouzerosqr==0),
    oouzerosqr == -1,
  Cif == end,
  Cif == (oouzerosqr<0),
    uzero == -1, (* -Sqrt[-1/oouzerosqr], *)
  Cif == else,
    uzero == Sqrt[1/oouzerosqr],
  Cif == end,

  (* killing vec xi^i in rotating frame, xi^a = (xi^0, xi^i) *)
  xi[a] == 0,

  (* constant F *)
  F == -(alpha2 - Psi4 delta[b,c] (beta[b] + xi[b]) (beta[c] + vR[c]) ) *
	uzero,

  (* set q *)
  Cif == ((bi==0) || (bi==1) || (bi==5)),  
    q == (C1/F - 1.0)/(n+1.0),
  Cif == else,
    q == (C2/F - 1.0)/(n+1.0),
  Cif == end,

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_compute_new_q.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_compute_new_q(tGrid *grid)\n"];
  pr["{\n"];

  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double C1 = Getd(\"BNSdata_C1\");\n"];
  pr["double C2 = Getd(\"BNSdata_C2\");\n"];
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

  (* prdecvarname[{g[a,b]},    "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{Psi},       "BNSdata_Psi"];
  prdecvarname[{alphaP},    "BNSdata_alphaP"];
  prdecvarname[{B[a]},      "BNSdata_Bx"];
  prdecvarname[{q},         "BNSdata_q"];
  prdecvarname[{vRS[a]},    "BNSdata_vRSx"];
  prdecvarname[{Sigma},     "BNSdata_Sigma"];
  prdecvarname[{dSigma[a]}, "BNSdata_Sigmax"];

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
